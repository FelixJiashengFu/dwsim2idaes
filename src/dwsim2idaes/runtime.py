"""Helpers used by the model scripts that dwsim2idaes generates.

Nothing here knows about DWSIM files: it only works with an IDAES model and the
reference stream table (the converged DWSIM results) embedded in the script.
"""
from __future__ import annotations

import contextlib
import logging
import math
from typing import Dict, List, Optional

import pyomo.environ as pyo
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.network import Port
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core.base.property_base import StateBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import solve_indexed_blocks
from idaes.core.util import model_statistics as _stats
from idaes.core.util.model_statistics import degrees_of_freedom

_log = logging.getLogger("dwsim2idaes")

MIN_MOLE_FRAC = 1e-8

# IPOPT's defaults push the starting point away from variable bounds and start
# with a large barrier parameter. That is right for a cold start, but here the
# starting point is the converged DWSIM solution, where phase flows sit exactly
# on their bounds (a vanished phase). These settings keep IPOPT near that point.
WARM_START = {"bound_push": 1e-8, "bound_frac": 1e-8, "mu_init": 1e-6}


# ---------------------------------------------------------------------------
# small utilities
# ---------------------------------------------------------------------------
def get(m, path: str):
    """Resolve ``"u_heater.outlet"`` relative to ``m.fs``."""
    obj = m.fs
    for part in path.split("."):
        obj = getattr(obj, part)
    return obj


def _t0(m):
    return m.fs.time.first()


def _member(port, name: str, t, idx=None):
    comp = getattr(port, name)
    return comp[t] if idx is None else comp[t, idx]


def _is_phase_outlet(port) -> bool:
    return port.local_name in ("vap_outlet", "liq_outlet")


def state_block(port, t):
    """State block behind a port, or None (e.g. the phase outlets of a Flash)."""
    if _is_phase_outlet(port):
        return None
    try:
        v = port.temperature[t]
        blk = v.parent_block()
        return blk if hasattr(blk, "phase_frac") else None
    except Exception:  # port members are expressions
        return None


def clean_composition(x: Dict[str, float], components) -> Dict[str, float]:
    """Floor zero mole fractions (log-form equilibrium needs x > 0) and renormalise."""
    xs = {j: max(float(x.get(j, 0.0)), MIN_MOLE_FRAC) for j in components}
    tot = sum(xs.values())
    return {j: v / tot for j, v in xs.items()}


# ---------------------------------------------------------------------------
# specifications
# ---------------------------------------------------------------------------
def fix_feed(port, *, flow_mol, temperature, pressure, mole_frac, vapor_fraction=None):
    """Fix a feed stream (F, T, P, x) at a unit inlet port.

    With ``vapor_fraction`` (a two-phase feed) the vapour fraction is fixed
    instead of the temperature: close to a bubble or dew point, a small
    difference between the DWSIM and IDAES property models would otherwise
    move a lot of material from one phase to the other at the given T.
    """
    blk = port.parent_block()
    t = blk.flowsheet().time.first()
    comps = [i[1] for i in port.mole_frac_comp.keys() if i[0] == t]
    xs = clean_composition(mole_frac, comps)
    port.flow_mol[t].fix(flow_mol)
    port.temperature[t].fix(temperature)
    port.pressure[t].fix(pressure)
    for j, v in xs.items():
        port.mole_frac_comp[t, j].fix(v)
    if vapor_fraction is not None:
        sb = state_block(port, t)
        if sb is None or "Vap" not in sb.phase_list or "Liq" not in sb.phase_list:
            return
        port.temperature[t].unfix()
        name = "dwsim_feed_vapor_fraction"
        while hasattr(sb, name):
            name += "_"
        sb.add_component(name, pyo.Constraint(expr=sb.phase_frac["Vap"] == float(vapor_fraction)))
        if not any(c.local_name.startswith("eps_t") for c in sb.component_objects(pyo.Param, descend_into=False)):
            # smooth VLE: the phase split stops depending on T outside the two-phase
            # window, so keep T between the bubble and dew points explicitly
            pair = ("Vap", "Liq")
            sb.add_component(name + "_tbub", pyo.Constraint(expr=sb.temperature >= sb.temperature_bubble[pair]))
            sb.add_component(name + "_tdew", pyo.Constraint(expr=sb.temperature <= sb.temperature_dew[pair]))


def close_loop(arc, *, flow_mol, mole_frac):
    """Fix amount and composition where a closed loop is cut open.

    The flow and mole-fraction equalities of ``arc`` are deactivated and the
    destination port gets the given values; temperature and pressure stay
    connected, so the loop is still closed energetically.
    """
    blk = arc.expanded_block
    for name in ("flow_mol_equality", "mole_frac_comp_equality"):
        c = getattr(blk, name, None)
        if c is None:
            raise ValueError(f"{arc.name}: no {name}; is the arc expanded?")
        c.deactivate()
    port = arc.destination
    t = arc.flowsheet().time.first() if hasattr(arc, "flowsheet") else 0.0
    for t in port.flow_mol.index_set() if hasattr(port.flow_mol, "index_set") else [t]:
        port.flow_mol[t].fix(flow_mol)
        comps = [j for (tt, j) in port.mole_frac_comp.keys() if tt == t]
        for j, v in clean_composition(mole_frac, comps).items():
            port.mole_frac_comp[t, j].fix(v)


def fix_outlet_vapor_fraction(unit, vapor_fraction: float):
    """Specify the molar vapour fraction at the outlet of a heater/cooler."""
    t = unit.flowsheet().time.first()
    sb = unit.control_volume.properties_out[t]
    vf = float(vapor_fraction)
    # a saturated outlet is specified through the dew/bubble temperature: the
    # smoothed VLE formulations only approach a vapour fraction of exactly 0 or 1
    if vf >= 1.0 - 1e-6:
        expr = sb.temperature == sb.temperature_dew[("Vap", "Liq")]
    elif vf <= 1e-6:
        expr = sb.temperature == sb.temperature_bubble[("Vap", "Liq")]
    else:
        expr = sb.phase_frac["Vap"] == vf
    unit.dwsim_vapor_fraction_spec = pyo.Constraint(expr=expr)


def fix_temperature_change(unit, delta_t: float):
    """Specify T_out - T_in for a heater/cooler."""
    t = unit.flowsheet().time.first()
    cv = unit.control_volume
    unit.dwsim_delta_t_spec = pyo.Constraint(
        expr=cv.properties_out[t].temperature - cv.properties_in[t].temperature
        == delta_t * pyo.units.K)


def fix_outlet_mass_flow(unit, outlet: str, flow_mass: float):
    """Specify the mass flow (kg/s) of one outlet of a splitter."""
    t = unit.flowsheet().time.first()
    sb = getattr(unit, outlet + "_state")[t]
    setattr(unit, f"dwsim_{outlet}_mass_flow_spec",
            pyo.Constraint(expr=sb.flow_mass == flow_mass * pyo.units.kg / pyo.units.s))


def set_mixer_pressure(unit, mode: str):
    """Outlet pressure rule for a Mixer built with ``MomentumMixingType.none``.

    ``mode`` is ``"average"`` or ``"maximum"`` (DWSIM's other option, the
    minimum, is IDAES' built-in behaviour).
    """
    t = unit.flowsheet().time.first()
    ps = [getattr(unit, n + "_state")[t].pressure for n in unit.config.inlet_list]
    out = unit.mixed_state[t].pressure
    if mode == "average":
        expr = out == sum(ps) / len(ps)
    elif mode == "maximum":
        eps = 1e-2 * pyo.units.Pa**2
        acc = ps[0]
        for p in ps[1:]:
            acc = 0.5 * (acc + p + ((acc - p) ** 2 + eps) ** 0.5)
        expr = out == acc
    else:
        raise ValueError(mode)
    unit.dwsim_pressure_rule = pyo.Constraint(expr=expr)


@contextlib.contextmanager
def fast_statistics():
    """Speed up IDAES' degrees-of-freedom counting.

    IDAES walks every constraint separately and re-expands shared named
    expressions each time. For cubic equations of state with many compounds
    that costs tens of seconds per call, and unit initialisation calls it
    repeatedly. Sharing one named-expression cache across the walk returns the
    same set of variables orders of magnitude faster.
    """
    original = _stats.variables_in_activated_equalities_set

    def variables_in_activated_equalities_set(block):
        cache: dict = {}
        var_set = ComponentSet()
        for c in _stats.activated_equalities_generator(block):
            var_set.update(identify_variables(c.body, named_expression_cache=cache))
        greybox = getattr(_stats, "greybox_variables", None)
        if greybox is not None:
            var_set.update(greybox(block))
        return var_set

    _stats.variables_in_activated_equalities_set = variables_in_activated_equalities_set
    try:
        yield
    finally:
        _stats.variables_in_activated_equalities_set = original


def dof(m) -> int:
    with fast_statistics():
        return degrees_of_freedom(m)


def check_dof(m):
    n = dof(m)
    if n != 0:
        raise RuntimeError(
            f"the converted model has {n} degrees of freedom (expected 0); "
            "the DWSIM specifications were not mapped one-to-one")
    return n


# ---------------------------------------------------------------------------
# initialisation from the DWSIM solution
# ---------------------------------------------------------------------------
def _set_port_guess(port, t, ref: dict):
    """Write the DWSIM stream state into the (free) variables behind a port."""
    def put(member, value, idx=None):
        try:
            v = _member(port, member, t, idx)
        except (AttributeError, KeyError):
            return
        if getattr(v, "is_variable_type", lambda: False)():
            if not v.fixed and value is not None and math.isfinite(value):
                v.set_value(value)

    put("flow_mol", max(ref["flow_mol"], 0.0))
    put("temperature", ref["temperature"])
    put("pressure", ref["pressure"])
    try:
        comps = [i[1] for i in port.mole_frac_comp.keys() if i[0] == t]
    except AttributeError:
        return
    xs = clean_composition(ref["x"], comps)
    for j, v in xs.items():
        put("mole_frac_comp", v, j)


def _set_state_guess(sb, ref: dict):
    """Also seed phase split and phase compositions when DWSIM reports them."""
    try:
        vf = ref.get("vapor_fraction")
        if vf is None or not math.isfinite(vf):
            return
        vf = min(max(vf, 0.0), 1.0)
        F = max(ref["flow_mol"], 0.0)
        if hasattr(sb, "flow_mol_phase"):
            for p, frac in (("Vap", vf), ("Liq", 1.0 - vf)):
                if p in sb.phase_list and not sb.flow_mol_phase[p].fixed:
                    sb.flow_mol_phase[p].set_value(max(F * frac, 1e-10 * max(F, 1.0)))
        for p, key in (("Vap", "y"), ("Liq", "xl")):
            comp = ref.get(key) or ref["x"]
            if p in sb.phase_list and hasattr(sb, "mole_frac_phase_comp"):
                xs = clean_composition(comp, [j for (pp, j) in sb.mole_frac_phase_comp if pp == p])
                for j, v in xs.items():
                    if not sb.mole_frac_phase_comp[p, j].fixed:
                        sb.mole_frac_phase_comp[p, j].set_value(v)
    except Exception as err:  # guesses are best effort
        _log.debug("could not seed phase guess: %s", err)


def set_initial_guess(m, streams: Dict[str, dict]):
    """Load the DWSIM solution into every stream port as an initial guess.

    Entries with a ``"block"`` key address a state block directly (used for
    states inside a unit, e.g. the two-phase mixture in a flash drum).
    """
    t = _t0(m)
    for tag, ref in streams.items():
        if "block" in ref:
            sb = get(m, ref["block"])[t]
            for name, val in (("flow_mol", ref["flow_mol"]), ("temperature", ref["temperature"]),
                              ("pressure", ref["pressure"])):
                v = getattr(sb, name)
                if not v.fixed and math.isfinite(val):
                    v.set_value(max(val, 0.0))
            xs = clean_composition(ref["x"], list(sb.mole_frac_comp.keys()))
            for j, v in xs.items():
                if not sb.mole_frac_comp[j].fixed:
                    sb.mole_frac_comp[j].set_value(v)
            _set_state_guess(sb, ref)
            continue
        for path in ref.get("ports", []):
            port = get(m, path)
            _set_port_guess(port, t, ref)
            sb = state_block(port, t)
            if sb is not None:
                _set_state_guess(sb, ref)


def _data_objects(block, ctype):
    for b in (block.values() if block.is_indexed() else [block]):
        yield from b.component_data_objects(ctype, descend_into=True)


def _snapshot(unit):
    """Record values, fixed flags and activity so a failed init can be undone."""
    vs = [(v, v.value, v.fixed) for v in _data_objects(unit, pyo.Var)]
    cs = [(c, c.active) for c in _data_objects(unit, pyo.Constraint)]
    bs = [(b, b.active) for b in _data_objects(unit, pyo.Block)]
    return vs, cs, bs


def _restore(snapshot):
    vs, cs, bs = snapshot
    for v, val, fixed in vs:
        v.set_value(val, skip_validation=True)
        v.fixed = fixed
    for comp, active in cs + bs:
        if comp.active != active:
            comp.activate() if active else comp.deactivate()


def _ok(results) -> bool:
    try:
        return pyo.check_optimal_termination(results)
    except Exception:
        return False


def _quiet_solve(block, options: Optional[dict] = None):
    """Solve a sub-block, returning True on success and never raising."""
    try:
        pyomo_log = logging.getLogger("pyomo")
        level = pyomo_log.level
        pyomo_log.setLevel(logging.CRITICAL)
        try:
            opt = get_solver(options=dict(WARM_START if options is None else options))
            if block.is_indexed():
                return _ok(solve_indexed_blocks(opt, [block], tee=False))
            return _ok(opt.solve(block, tee=False))
        finally:
            pyomo_log.setLevel(level)
    except Exception:
        return False


def _seed_isentropic_state(unit, t):
    """Start the isentropic state of a compressor/expander at the outlet guess."""
    iso = getattr(unit, "properties_isentropic", None)
    cv = getattr(unit, "control_volume", None)
    if iso is None or cv is None:
        return
    src, dst = cv.properties_out[t], iso[t]
    for v in src.component_data_objects(pyo.Var, descend_into=False):
        d = dst.find_component(v.local_name)
        if d is not None and not d.fixed and v.value is not None:
            d.set_value(v.value, skip_validation=True)


def _seed_energy_terms(unit, t):
    """Give a free duty or work variable the value the energy balance implies."""
    cv = getattr(unit, "control_volume", None)
    bal = getattr(cv, "enthalpy_balances", None) if cv is not None else None
    if bal is None:
        return
    for name in ("heat_duty", "work_mechanical"):
        var = getattr(unit, name, None)
        if var is not None and not var[t].fixed:
            try:
                calculate_variable_from_constraint(var[t], bal[t])
            except Exception as err:
                _log.debug("could not seed %s of %s: %s", name, unit.name, err)
            return


def _direct_unit_init(unit, outlvl) -> bool:
    """Fallback initialisation: start from the DWSIM state instead of from scratch.

    1. every state block with a fully defined state is solved on its own at the
       DWSIM temperature, pressure, flow and composition already loaded in it;
    2. the unit is then solved as a square problem with its inlets fixed.
    """
    t = unit.flowsheet().time.first()
    for sb in unit.component_objects(pyo.Block, descend_into=True):
        if not isinstance(sb, StateBlock):
            continue
        snapshot = _snapshot(sb)
        try:
            # fixes the state at its current (DWSIM) value, solves, releases it
            sb.initialize(outlvl=outlvl)
        except Exception as err:
            _restore(snapshot)
            _log.debug("state block %s not initialised: %s", sb.name, str(err)[:120])
    _seed_isentropic_state(unit, t)
    _seed_energy_terms(unit, t)
    freed = []
    for port in unit.component_objects(Port, descend_into=False):
        if "inlet" not in port.local_name:
            continue
        for v in port.iter_vars(fixed=False):
            v.fix()
            freed.append(v)
    try:
        return degrees_of_freedom(unit) == 0 and _quiet_solve(unit)
    finally:
        for v in freed:
            v.unfix()


def initialize(m, streams: Dict[str, dict], units: List[str], outlvl=logging.WARNING) -> List[str]:
    """Initialise the flowsheet without sequential decomposition.

    Because DWSIM already converged the flowsheet, every unit's inlet state is
    known. Each unit is therefore initialised on its own, at the DWSIM inlet
    conditions, in any order; recycle loops need no tearing. If a unit's own
    IDAES routine fails, the unit is solved directly from the DWSIM state.
    Returns the units for which both attempts failed (the full solve may still
    succeed from the DWSIM guess).
    """
    set_initial_guess(m, streams)
    failed = []
    with fast_statistics():
        for name in units:
            unit = get(m, name)
            snapshot = _snapshot(unit)
            try:
                unit.initialize(outlvl=outlvl)
                continue
            except Exception as err:
                _restore(snapshot)
                _log.info("IDAES initialisation of %s failed (%s); solving it from the DWSIM state",
                          name, str(err)[:200])
            if not _direct_unit_init(unit, outlvl):
                _restore(snapshot)
                failed.append(name)
                _log.warning("could not initialise %s; continuing from the DWSIM guess", name)
    return failed


def set_vle_smoothing(m, rel: float = 1e-4, floor: float = 1e-8):
    """Scale the smoothing parameters of the complementarity VLE formulation.

    IDAES uses ``eps = 1 mol/s`` by default, which on streams of a few mol/s
    blurs the phase boundary (a superheated vapour shows a few tenths of a
    percent of liquid). Here every state block gets ``eps = rel * flow``.
    Returns the number of parameters changed.
    """
    n = 0
    for p in m.component_data_objects(pyo.Param, descend_into=True):
        name = p.parent_component().local_name
        if not (name.startswith("eps_t") or name.startswith("eps_z")):
            continue
        sb = p.parent_block()
        try:
            flow = abs(pyo.value(sb.flow_mol))
        except Exception:
            continue
        p.set_value(max(rel * flow, floor))
        n += 1
    return n


def solve(m, tee: bool = False, solver: Optional[str] = None, options: Optional[dict] = None,
          smoothing: Optional[float] = 1e-4):
    """Solve the whole flowsheet simultaneously with IPOPT.

    The VLE smoothing parameters are first scaled to ``smoothing`` times each
    stream's flow (see :func:`set_vle_smoothing`) and IPOPT is warm-started
    (:data:`WARM_START`). If that does not converge, the starting point is
    restored and IPOPT's cold-start defaults are tried, and finally IDAES'
    default smoothing (1 mol/s). ``results.dwsim2idaes_smoothing`` records the
    smoothing that was used. Pass ``options`` to give IPOPT options yourself.
    """
    if options is not None:
        attempts = [(options, smoothing)]
    else:
        attempts = [(WARM_START, smoothing), ({}, smoothing)]
        if smoothing is not None:
            attempts.append((WARM_START, None))
    snapshot = _snapshot(m) if len(attempts) > 1 else None
    results = None
    for i, (opts, eps) in enumerate(attempts):
        if i > 0:
            _restore(snapshot)
        if eps is None:
            _reset_vle_smoothing(m)
        else:
            set_vle_smoothing(m, eps)
        opt = get_solver(solver, options=dict(opts)) if solver else get_solver(options=dict(opts))
        pyomo_log = logging.getLogger("pyomo")
        level = pyomo_log.level
        if i < len(attempts) - 1:
            pyomo_log.setLevel(logging.CRITICAL)      # a failed first attempt is not news
        try:
            results = opt.solve(m, tee=tee, load_solutions=False)
        except Exception as err:
            if i == len(attempts) - 1:
                raise
            _log.info("solve attempt %d raised %s; trying the next settings", i + 1, str(err)[:120])
            continue
        finally:
            pyomo_log.setLevel(level)
        if len(results.solution) > 0:
            m.solutions.load_from(results)
        results.dwsim2idaes_smoothing = eps
        if _ok(results):
            if eps is None and smoothing is not None:
                _log.warning("converged only with the default VLE smoothing of IDAES (1 mol/s); "
                             "phase fractions close to a phase boundary are approximate")
            break
    return results


def _reset_vle_smoothing(m):
    for p in m.component_data_objects(pyo.Param, descend_into=True):
        name = p.parent_component().local_name
        if name.startswith("eps_t") or name.startswith("eps_z"):
            p.set_value(1.0)


# ---------------------------------------------------------------------------
# results and comparison with DWSIM
# ---------------------------------------------------------------------------
def _val(x):
    try:
        return float(pyo.value(x))
    except Exception:
        return float("nan")


def stream_results(m, streams: Dict[str, dict]) -> Dict[str, dict]:
    """IDAES stream table in the same shape as the embedded DWSIM table."""
    t = _t0(m)
    out = {}
    for tag, ref in streams.items():
        if not ref.get("ports"):
            continue                      # internal state, not a DWSIM stream
        path = ref["ports"][0]
        port = get(m, path)
        row = {
            "flow_mol": _val(port.flow_mol[t]),
            "temperature": _val(port.temperature[t]),
            "pressure": _val(port.pressure[t]),
            "x": {},
        }
        for (tt, j) in port.mole_frac_comp.keys():
            if tt == t:
                row["x"][j] = _val(port.mole_frac_comp[tt, j])
        sb = state_block(port, t)
        if path.endswith("vap_outlet"):
            row["vapor_fraction"] = 1.0
        elif path.endswith("liq_outlet"):
            row["vapor_fraction"] = 0.0
        elif sb is not None:
            row["vapor_fraction"] = _val(sb.phase_frac["Vap"]) if "Vap" in sb.phase_list else 0.0
        else:
            row["vapor_fraction"] = float("nan")
        out[tag] = row
    return out


def compare(m, streams: Dict[str, dict], flow_floor: float = 1e-6) -> List[dict]:
    """Per-stream deviations between the solved IDAES model and DWSIM."""
    res = stream_results(m, streams)
    fmax = max([abs(r["flow_mol"]) for r in streams.values() if r.get("ports")] + [1e-30])
    rows = []
    for tag, ida in res.items():
        ref = streams[tag]
        negligible = abs(ref["flow_mol"]) <= flow_floor * fmax
        dx = 0.0
        if not negligible:
            for j, v in ida["x"].items():
                dx = max(dx, abs(v - ref["x"].get(j, 0.0)))
        dvf = ida["vapor_fraction"] - ref.get("vapor_fraction", float("nan"))
        if ref.get("supercritical"):
            dvf = 0.0           # "vapour" vs "liquid" is only a label above the critical point
        rows.append({
            "stream": tag,
            "is_feed": bool(ref.get("is_feed")),
            "negligible_flow": negligible,
            "T_dwsim": ref["temperature"], "T_idaes": ida["temperature"],
            "dT": ida["temperature"] - ref["temperature"],
            "P_dwsim": ref["pressure"], "P_idaes": ida["pressure"],
            "dP_rel": (ida["pressure"] - ref["pressure"]) / ref["pressure"] if ref["pressure"] else float("nan"),
            "F_dwsim": ref["flow_mol"], "F_idaes": ida["flow_mol"],
            "dF_rel": (ida["flow_mol"] - ref["flow_mol"]) / fmax,
            "vf_dwsim": ref.get("vapor_fraction", float("nan")), "vf_idaes": ida["vapor_fraction"],
            "dvf": dvf,
            "dx_max": dx,
        })
    return rows


def summarize(rows: List[dict]) -> dict:
    """Worst deviations over the streams the model actually computes."""
    calc = [r for r in rows if not r["is_feed"] and not r["negligible_flow"]]

    def worst(key):
        vals = [abs(r[key]) for r in calc if math.isfinite(r[key])]
        return max(vals) if vals else 0.0

    return {
        "streams_compared": len(calc),
        "max_abs_dT_K": worst("dT"),
        "max_rel_dT": max([abs(r["dT"]) / r["T_dwsim"] for r in calc if r["T_dwsim"]] or [0.0]),
        "max_rel_dP": worst("dP_rel"),
        "max_rel_dF": worst("dF_rel"),
        "max_abs_dvf": worst("dvf"),
        "max_abs_dx": worst("dx_max"),
    }


def format_comparison(rows: List[dict]) -> str:
    head = (f"{'stream':<24}{'T DWSIM':>10}{'T IDAES':>10}{'dT':>8}"
            f"{'P DWSIM':>12}{'P IDAES':>12}{'F DWSIM':>12}{'F IDAES':>12}"
            f"{'vf DW':>8}{'vf ID':>8}{'max dx':>9}")
    lines = [head, "-" * len(head)]
    for r in rows:
        mark = " (feed)" if r["is_feed"] else ""
        lines.append(
            f"{(r['stream'] + mark)[:23]:<24}{r['T_dwsim']:>10.2f}{r['T_idaes']:>10.2f}{r['dT']:>8.2f}"
            f"{r['P_dwsim']:>12.4g}{r['P_idaes']:>12.4g}{r['F_dwsim']:>12.5g}{r['F_idaes']:>12.5g}"
            f"{r['vf_dwsim']:>8.4f}{r['vf_idaes']:>8.4f}{r['dx_max']:>9.1e}")
    s = summarize(rows)
    lines.append("")
    lines.append(
        f"{s['streams_compared']} calculated streams compared: "
        f"max |dT| = {s['max_abs_dT_K']:.3g} K ({s['max_rel_dT']:.2%}), "
        f"max |dP| = {s['max_rel_dP']:.2%}, max |dF| = {s['max_rel_dF']:.2%}, "
        f"max |d vapour fraction| = {s['max_abs_dvf']:.3g}, max |dx| = {s['max_abs_dx']:.2g}")
    return "\n".join(lines)


def unit_duties(m, duties: Dict[str, dict]) -> List[dict]:
    """Compare unit duties/powers (W) with the values stored by DWSIM."""
    t = _t0(m)
    rows = []
    for tag, d in duties.items():
        unit = get(m, d["unit"])
        var = getattr(unit, d["var"])
        val = _val(var[t]) * d.get("sign", 1.0) / d.get("efficiency", 1.0)
        ref = d["dwsim"]
        rows.append({"unit": tag, "kind": d["var"], "dwsim_W": ref, "idaes_W": val,
                     "rel": (val - ref) / ref if ref else float("nan")})
    return rows


def format_duties(rows: List[dict]) -> str:
    if not rows:
        return ""
    lines = [f"{'unit':<24}{'quantity':<18}{'DWSIM [W]':>14}{'IDAES [W]':>14}{'rel. dev.':>11}",
             "-" * 81]
    for r in rows:
        lines.append(f"{r['unit'][:23]:<24}{r['kind']:<18}{r['dwsim_W']:>14.6g}{r['idaes_W']:>14.6g}{r['rel']:>11.2%}")
    return "\n".join(lines)
