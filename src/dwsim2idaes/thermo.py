"""Build IDAES generic (modular) property packages from DWSIM property data.

The input is a small JSON-friendly ``thermo spec`` (see :func:`thermo_spec`)
so that generated model scripts can carry their own property data.
"""
from __future__ import annotations

import math
from typing import Dict, Iterable, Optional

from pyomo.environ import log, value, units as pyunits

from idaes.core import Component, LiquidPhase, VaporPhase
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import (
    CubicComplementarityVLE,
    SmoothVLE,
)
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
    LogBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity, log_fugacity
from idaes.models.properties.modular_properties.state_definitions import FTPx

from . import correlations as corr
from .ir import Flowsheet, PropertyPackage

SUPPORTED_PACKAGES = {"PengRobinson": "PR", "SRK": "SRK", "Raoult": "Ideal"}


class ThermoError(Exception):
    pass


# ---------------------------------------------------------------------------
# Pure-component methods that reproduce the DWSIM correlations
# ---------------------------------------------------------------------------
def _T(b, T):
    return pyunits.convert(T, to_units=pyunits.K) / pyunits.K      # dimensionless K


def _Tref(b) -> float:
    return value(pyunits.convert(b.params.temperature_ref, to_units=pyunits.K))


def _dw(cobj) -> dict:
    return cobj.config.parameter_data["dwsim"]


class DwsimPure:
    """DWSIM/ChemSep correlations for ideal gas Cp, Psat, Hvap and liquid density.

    Reference state is the one DWSIM uses: ideal gas at 298.15 K (and 101325 Pa
    for entropy) has zero enthalpy/entropy; heats of formation are not included.
    """

    class cp_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            pass

        @staticmethod
        def return_expression(b, cobj, T):
            cp = corr.cp_ig_value(_dw(cobj)["cp_ig"], _T(b, T)) * pyunits.J / pyunits.mol / pyunits.K
            return pyunits.convert(cp, b.params.get_metadata().derived_units.HEAT_CAPACITY_MOLE)

    class enth_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            pass

        @staticmethod
        def return_expression(b, cobj, T):
            h = corr.enth_ig_value(_dw(cobj)["cp_ig"], _T(b, T), _Tref(b)) * pyunits.J / pyunits.mol
            return pyunits.convert(h, b.params.get_metadata().derived_units.ENERGY_MOLE)

    class entr_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            pass

        @staticmethod
        def return_expression(b, cobj, T):
            s = corr.entr_ig_value(_dw(cobj)["cp_ig"], _T(b, T), _Tref(b)) * pyunits.J / pyunits.mol / pyunits.K
            return pyunits.convert(s, b.params.get_metadata().derived_units.ENTROPY_MOLE)

    class pressure_sat_comp:
        @staticmethod
        def build_parameters(cobj):
            pass

        @staticmethod
        def return_expression(b, cobj, T, dT=False):
            d = _dw(cobj)
            du = b.params.get_metadata().derived_units
            Tk = _T(b, T)
            if dT:
                # central difference: exact enough for every DWSIM form and
                # keeps one code path for all equation numbers
                h = 1e-2
                dp = (corr.psat_value(d["psat"], Tk + h, d["tc"], d["pc"], d["omega"])
                      - corr.psat_value(d["psat"], Tk - h, d["tc"], d["pc"], d["omega"])) / (2 * h)
                return pyunits.convert(dp * pyunits.Pa / pyunits.K, du.PRESSURE / du.TEMPERATURE)
            p = corr.psat_value(d["psat"], Tk, d["tc"], d["pc"], d["omega"]) * pyunits.Pa
            return pyunits.convert(p, du.PRESSURE)

    class dens_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            pass

        @staticmethod
        def return_expression(b, cobj, T):
            d = _dw(cobj)
            rho = corr.dens_liq_value(d["dens_liq"], _T(b, T), d["tc"], d["pc"]) * pyunits.mol / pyunits.m**3
            return pyunits.convert(rho, b.params.get_metadata().derived_units.DENSITY_MOLE)

    class enth_mol_liq_comp:
        """DWSIM "Ideal" mode: H_liq = H_ig(T) - Hvap(T)."""

        @staticmethod
        def build_parameters(cobj):
            pass

        @staticmethod
        def return_expression(b, cobj, T):
            d = _dw(cobj)
            hv = corr.hvap_value(d["hvap"], _T(b, T), d["tc"]) * pyunits.J / pyunits.mol
            h = DwsimPure.enth_mol_ig_comp.return_expression(b, cobj, T)
            return h - pyunits.convert(hv, b.params.get_metadata().derived_units.ENERGY_MOLE)

    class entr_mol_liq_comp:
        """S_liq = S_ig(T, Psat) - Hvap/T (thermodynamically consistent form)."""

        @staticmethod
        def build_parameters(cobj):
            pass

        @staticmethod
        def return_expression(b, cobj, T):
            d = _dw(cobj)
            Tk = _T(b, T)
            hv = corr.hvap_value(d["hvap"], Tk, d["tc"])
            ps = corr.psat_value(d["psat"], Tk, d["tc"], d["pc"], d["omega"])
            pref = pyunits.convert(b.params.pressure_ref, to_units=pyunits.Pa) / pyunits.Pa
            ds = (-hv / Tk - corr.R_GAS * log(ps / pref)) * pyunits.J / pyunits.mol / pyunits.K
            s = DwsimPure.entr_mol_ig_comp.return_expression(b, cobj, T)
            return s + pyunits.convert(ds, b.params.get_metadata().derived_units.ENTROPY_MOLE)


# ---------------------------------------------------------------------------
# thermo spec  (Flowsheet + package  ->  JSON-friendly dict)
# ---------------------------------------------------------------------------
def thermo_spec(fs: Flowsheet, pkg: PropertyPackage, compounds: Optional[Iterable[str]] = None,
                t_range=None, p_range=None, flow_max=None) -> dict:
    """Collect everything the IDAES property package needs, as plain data."""
    if pkg.kind not in SUPPORTED_PACKAGES:
        raise ThermoError(
            f"property package '{pkg.kind}' is not supported yet "
            f"(supported: {', '.join(sorted(SUPPORTED_PACKAGES))})")
    model = SUPPORTED_PACKAGES[pkg.kind]
    names = list(compounds) if compounds is not None else list(fs.compounds)
    ts = [s.temperature for s in fs.streams.values() if math.isfinite(s.temperature)]
    ps = [s.pressure for s in fs.streams.values() if math.isfinite(s.pressure) and s.pressure > 0]
    fl = [s.flow_mol for s in fs.streams.values() if math.isfinite(s.flow_mol)]
    t_min = min(ts + [298.15])
    t_lo, t_hi = t_range or (min(0.5 * t_min, max(50.0, t_min - 150.0)), max(ts + [298.15]) + 300.0)
    p_lo, p_hi = p_range or (min(ps + [101325.0]) / 50.0, max(ps + [101325.0]) * 20.0)
    f_hi = flow_max or max(fl + [1.0]) * 100.0
    f_nom = max(fl + [1.0])

    comps, notes = {}, []
    for n in names:
        c = fs.compounds[n]
        if c.is_petroleum_fraction:
            raise ThermoError(f"compound '{n}' is a petroleum fraction / hypothetical; not supported yet")
        if not (c.tc > 0 and c.pc > 0 and c.mw > 0):
            raise ThermoError(f"compound '{n}' has no critical properties or molar weight")
        try:
            cp = corr.cp_ig_spec(c.db, c.cp_ig.eq, c.cp_ig.coeffs)
        except corr.CorrelationError as e:
            raise ThermoError(f"compound '{n}': {e}") from e
        entry = {
            "name": n, "mw": c.mw, "tc": c.tc, "pc": c.pc, "omega": c.omega,
            "cp_ig": cp,
            "psat": corr.psat_spec(c.db, c.psat.eq, c.psat.coeffs),
            "dens_liq": corr.dens_liq_spec(c.db, c.dens_liq.eq, c.dens_liq.coeffs, c.tc, c.pc, c.omega, c.zra),
        }
        if model == "Ideal":
            try:
                entry["hvap"] = corr.hvap_spec(c.db, c.hvap.eq, c.hvap.coeffs)
            except corr.CorrelationError as e:
                raise ThermoError(f"compound '{n}': {e}") from e
        comps[c.key] = entry

    kij = {}
    if model in ("PR", "SRK"):
        for pair, v in pkg.kij.items():
            a, b = pair.split("|")
            if a in names and b in names and v != 0.0:
                kij[f"{fs.compounds[a].key}|{fs.compounds[b].key}"] = v

    mode = pkg.options.get("EnthalpyEntropyCpCvCalculationMode", "")
    if model in ("PR", "SRK") and mode == "LeeKesler":
        notes.append("the DWSIM file asks for Lee-Kesler enthalpy/entropy; the IDAES model uses "
                     f"{model} departure functions, so duties can differ slightly")
    if model == "Ideal" and mode in ("Excess", "ExpData"):
        notes.append(f"DWSIM enthalpy mode '{mode}' is mapped to the ideal form H_liq = H_ig - Hvap")
    return {
        "model": model,
        "package": pkg.kind,
        "components": comps,
        "kij": kij,
        "bounds": {
            "flow_mol": [0.0, f_nom, f_hi],
            "temperature": [t_lo, 298.15, t_hi],
            "pressure": [p_lo, 101325.0, p_hi],
        },
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# thermo spec  ->  IDAES GenericParameterBlock configuration
# ---------------------------------------------------------------------------
def build_property_config(spec: Dict, phases: str = "VLE") -> dict:
    """Return the ``configuration`` dict for ``GenericParameterBlock(**cfg)``.

    ``phases`` selects the phase set: ``"VLE"`` (vapour-liquid equilibrium),
    ``"Vap"`` (vapour only) or ``"Liq"`` (liquid only). The single-phase
    variants share the state variables of the VLE package, so units built on
    different variants can be connected directly with Arcs.
    """
    if phases not in ("VLE", "Vap", "Liq"):
        raise ValueError(f"phases must be 'VLE', 'Vap' or 'Liq', not {phases!r}")
    model = spec["model"]
    cubic = model in ("PR", "SRK")
    form = log_fugacity if cubic else fugacity
    components = {}
    for key, c in spec["components"].items():
        pdata = {
            "mw": (c["mw"] * 1e-3, pyunits.kg / pyunits.mol),
            "pressure_crit": (c["pc"], pyunits.Pa),
            "temperature_crit": (c["tc"], pyunits.K),
            "omega": c["omega"],
            "dwsim": {k: c.get(k) for k in ("tc", "pc", "omega", "cp_ig", "psat", "hvap", "dens_liq")},
        }
        comp = {
            "type": Component,
            "cp_mol_ig_comp": DwsimPure.cp_mol_ig_comp,
            "enth_mol_ig_comp": DwsimPure.enth_mol_ig_comp,
            "entr_mol_ig_comp": DwsimPure.entr_mol_ig_comp,
            "pressure_sat_comp": DwsimPure.pressure_sat_comp,
            "dens_mol_liq_comp": DwsimPure.dens_mol_liq_comp,
            "parameter_data": pdata,
        }
        if phases == "VLE":
            comp["phase_equilibrium_form"] = {("Vap", "Liq"): form}
        if not cubic:
            comp["enth_mol_liq_comp"] = DwsimPure.enth_mol_liq_comp
            comp["entr_mol_liq_comp"] = DwsimPure.entr_mol_liq_comp
        components[key] = comp

    if cubic:
        ctype = CubicType.PR if model == "PR" else CubicType.SRK
        eos = {"equation_of_state": Cubic, "equation_of_state_options": {"type": ctype}}
        # The complementarity formulation needs no bubble/dew points and copes
        # with supercritical compounds; the IDAES binaries only ship it for PR.
        vle_state = CubicComplementarityVLE if model == "PR" else SmoothVLE
        bubble_dew = LogBubbleDew
    else:
        eos = {"equation_of_state": Ideal}
        vle_state, bubble_dew = SmoothVLE, IdealBubbleDew

    b = spec["bounds"]
    cfg = {
        "components": components,
        "phases": {k: v for k, v in (
            ("Liq", dict(type=LiquidPhase, **eos)),
            ("Vap", dict(type=VaporPhase, **eos)),
        ) if phases in ("VLE", k)},
        "base_units": {
            "time": pyunits.s, "length": pyunits.m, "mass": pyunits.kg,
            "amount": pyunits.mol, "temperature": pyunits.K,
        },
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (b["flow_mol"][0], b["flow_mol"][1], b["flow_mol"][2], pyunits.mol / pyunits.s),
            "temperature": (b["temperature"][0], b["temperature"][1], b["temperature"][2], pyunits.K),
            "pressure": (b["pressure"][0], b["pressure"][1], b["pressure"][2], pyunits.Pa),
        },
        "pressure_ref": (101325.0, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
        "include_enthalpy_of_formation": False,
    }
    if phases == "VLE":
        cfg["phases_in_equilibrium"] = [("Vap", "Liq")]
        cfg["phase_equilibrium_state"] = {("Vap", "Liq"): vle_state}
        cfg["bubble_dew_method"] = bubble_dew
    if cubic:
        kappa = {}
        keys = list(spec["components"])
        for i in keys:
            for j in keys:
                kappa[(i, j)] = spec["kij"].get(f"{i}|{j}", spec["kij"].get(f"{j}|{i}", 0.0))
        cfg["parameter_data"] = {f"{model}_kappa": kappa}
    return cfg
