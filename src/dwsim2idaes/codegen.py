"""Turn a DWSIM :class:`~dwsim2idaes.ir.Flowsheet` into an IDAES model script."""
from __future__ import annotations

import json
import math
import pprint
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

from . import __version__
from .ir import Flowsheet, Stream, Unit
from .thermo import SUPPORTED_PACKAGES, ThermoError, thermo_spec

# units that only pass a material stream through (tear/recycle blocks)
PASSTHROUGH = {"Recycle"}
# objects that carry no material and are ignored
LOGICAL = {
    "Adjust", "Spec", "EnergyRecycle", "PIDController", "PythonController", "LevelGauge",
    "AnalogGauge", "DigitalGauge", "Input", "Switch", "ControllerBlock",
}


# A stream counts as single-phase when the other phase holds less than this molar
# fraction (flash outlets in DWSIM are saturated only to the flash tolerance).
PHASE_TOL = 1e-4


class ConversionError(Exception):
    """The flowsheet cannot be converted; ``problems`` lists the reasons."""

    def __init__(self, problems: List[str]):
        self.problems = list(problems)
        super().__init__("; ".join(self.problems))


@dataclass
class UnitPlan:
    key: str
    tag: str
    kind: str
    idaes: str = ""
    build: List[str] = field(default_factory=list)
    specs: List[str] = field(default_factory=list)
    in_ports: Dict[int, str] = field(default_factory=dict)
    out_ports: Dict[int, str] = field(default_factory=dict)
    init: List[str] = field(default_factory=list)
    internal_arcs: List[Tuple[str, str, str]] = field(default_factory=list)
    duty: Optional[dict] = None
    # guesses for states inside the unit: block path -> (streams in, T, P, vapour, liquid)
    internal: Dict[str, dict] = field(default_factory=dict)


@dataclass
class Conversion:
    code: str
    flowsheet: Flowsheet
    thermo: dict
    streams: Dict[str, dict]
    units: List[UnitPlan]
    notes: List[str]
    cut_units: List[str]
    ignored_units: List[str]

    def summary(self) -> str:
        lines = [f"{len(self.units)} unit operations, {len(self.streams)} material streams, "
                 f"{len(self.thermo['components'])} compounds, property model {self.thermo['model']}"]
        for u in self.units:
            lines.append(f"  {u.tag:<24} {u.kind:<16} -> {u.idaes}")
        for n in self.notes:
            lines.append(f"  note: {n}")
        return "\n".join(lines)


def _num(x: float) -> str:
    return repr(float(f"{float(x):.12g}"))


class _Ctx:
    def __init__(self, fs: Flowsheet, phase_mode: str = "auto"):
        self.fs = fs
        self.notes: List[str] = []
        self.imports: set = set()
        self.phase_mode = phase_mode
        self.phase_sets: set = set()
        flows = [s.flow_mol for s in fs.streams.values() if math.isfinite(s.flow_mol)]
        self.fmax = max(flows + [1e-30])

    # -- phase set selection ------------------------------------------------
    def label(self, s: Stream) -> str:
        """'V', 'L', '2' (two-phase) or '' (negligible flow) from the DWSIM result."""
        if not (s.flow_mol > 1e-9 * self.fmax):
            return ""
        if s.vapor_fraction >= 1.0 - PHASE_TOL:
            return "V"
        if s.vapor_fraction <= PHASE_TOL:
            return "L"
        return "2"

    def pseudo_critical(self, s: Stream) -> Tuple[float, float]:
        tc = sum(x * self.fs.compounds[n].tc for n, x in s.z.items() if n in self.fs.compounds)
        pc = sum(x * self.fs.compounds[n].pc for n, x in s.z.items() if n in self.fs.compounds)
        return tc, pc

    def supercritical(self, s: Stream) -> bool:
        tc, pc = self.pseudo_critical(s)
        return s.pressure > pc or s.temperature > tc

    def props(self, streams: List[Optional[Stream]], force_vle: bool = False) -> str:
        """Property block for a unit, chosen from the phases DWSIM found."""
        ss = [s for s in streams if s is not None]
        labels = {self.label(s) for s in ss} - {""}
        if force_vle or self.phase_mode == "vle" or "2" in labels or not labels:
            mode = "VLE"
        elif labels == {"V"}:
            mode = "Vap"
        elif labels == {"L"}:
            mode = "Liq"
        elif all(s.pressure > self.pseudo_critical(s)[1] for s in ss if self.label(s)):
            mode = "Vap"          # above the critical pressure there is one fluid root
        else:
            mode = "VLE"
        self.phase_sets.add(mode)
        return {"VLE": "m.fs.props", "Vap": "m.fs.props_vap", "Liq": "m.fs.props_liq"}[mode]

    def unit_props(self, u: Unit, force_vle: bool = False) -> str:
        ss = [self.fs.streams[i] for i in list(u.inlets.values()) + list(u.outlets.values())]
        return self.props(ss, force_vle)

    def note(self, msg: str):
        if msg not in self.notes:
            self.notes.append(msg)

    def stream(self, sid: str) -> Stream:
        return self.fs.streams[sid]

    def f(self, unit: Unit, name: str, default: float = float("nan")) -> float:
        v = unit.params.get(name)
        if v is None or isinstance(v, list):
            return default
        try:
            x = float(v)
        except ValueError:
            return default
        return default if math.isnan(x) else x

    def s(self, unit: Unit, name: str, default: str = "") -> str:
        v = unit.params.get(name)
        return default if v is None or isinstance(v, list) else str(v)

    def energy_power(self, unit: Unit) -> Optional[float]:
        for eid in unit.energy:
            p = self.fs.energy_streams[eid].power
            if math.isfinite(p):
                return p
        return None

    def out_stream(self, unit: Unit, idx: int = 0) -> Optional[Stream]:
        sid = unit.outlets.get(idx)
        return self.fs.streams[sid] if sid else None

    def in_stream(self, unit: Unit, idx: int = 0) -> Optional[Stream]:
        sid = unit.inlets.get(idx)
        return self.fs.streams[sid] if sid else None


def _require(ctx: _Ctx, unit: Unit, n_in: Optional[int] = None, n_out: Optional[int] = None):
    if n_in is not None and len(unit.inlets) != n_in:
        raise ConversionError([f"{unit.kind} '{unit.tag}' has {len(unit.inlets)} connected inlet(s), expected {n_in}"])
    if n_out is not None and len(unit.outlets) != n_out:
        raise ConversionError([f"{unit.kind} '{unit.tag}' has {len(unit.outlets)} connected outlet(s), expected {n_out}"])


# ---------------------------------------------------------------------------
# unit handlers
# ---------------------------------------------------------------------------
def _heater(ctx: _Ctx, u: Unit, cooler: bool = False) -> UnitPlan:
    _require(ctx, u, 1, 1)
    ctx.imports.add("Heater")
    out = ctx.out_stream(u)
    mode = ctx.s(u, "CalcMode", "OutletTemperature")
    vf_spec = ctx.f(u, "OutletVaporFraction", out.vapor_fraction)
    vf_ok = 0.0 <= vf_spec <= 1.0 and abs(vf_spec - out.vapor_fraction) < 1e-3
    p = UnitPlan(u.key, u.tag, u.kind, idaes="Heater")
    props = ctx.unit_props(u, force_vle=(mode == "OutletVaporFraction" and vf_ok))
    p.build.append(f"m.fs.{u.key} = Heater(property_package={props}, has_pressure_change=True)")
    p.in_ports, p.out_ports, p.init = {0: f"{u.key}.inlet"}, {0: f"{u.key}.outlet"}, [u.key]
    eff = ctx.f(u, "Eficiencia", 100.0) / 100.0 or 1.0
    sign = -1.0 if cooler else 1.0
    dp = ctx.f(u, "DeltaP", 0.0)
    sp = p.specs
    sp.append(f"u.deltaP.fix({_num(-dp)})  # Pa (DWSIM pressure drop {dp:g} Pa)")
    if mode == "OutletTemperature":
        sp.append(f"u.outlet.temperature.fix({_num(ctx.f(u, 'OutletTemperature', out.temperature))})  # K")
    elif mode in ("HeatAdded", "HeatRemoved", "HeatAddedRemoved"):
        q = ctx.f(u, "DeltaQ", 0.0) * 1000.0 * eff * sign
        sp.append(f"u.heat_duty.fix({_num(q)})  # W into the fluid")
    elif mode == "EnergyStream":
        e = ctx.energy_power(u)
        if e is None:
            e = ctx.f(u, "DeltaQ", 0.0) * 1000.0
        sp.append(f"u.heat_duty.fix({_num(e * eff)})  # W, from the attached energy stream")
        if cooler:
            ctx.note(f"cooler '{u.tag}' takes its duty from an energy stream; DWSIM adds that energy "
                     "to the fluid, and so does the converted model")
    elif mode == "TemperatureChange":
        dt = ctx.f(u, "DeltaT", 0.0)
        sp.append(f"rt.fix_temperature_change(u, {_num(dt)})  # K, T_out - T_in")
    elif mode == "OutletVaporFraction" and not vf_ok:
        sp.append(f"u.outlet.temperature.fix({_num(out.temperature)})  # K, DWSIM result")
        ctx.note(f"'{u.tag}' specifies an outlet vapour fraction of "
                 f"{ctx.f(u, 'OutletVaporFraction', float('nan')):g} but the DWSIM result has "
                 f"{out.vapor_fraction:.4g}; the outlet temperature is fixed at the DWSIM result instead")
    elif mode == "OutletVaporFraction":
        vf = ctx.f(u, "OutletVaporFraction", out.vapor_fraction)
        sp.append(f"rt.fix_outlet_vapor_fraction(u, {_num(vf)})")
        if vf <= 0.0 or vf >= 1.0:
            ctx.note(f"'{u.tag}' specifies a saturated outlet (vapour fraction {vf:g}); the converted "
                     f"model sets the outlet temperature to the {'dew' if vf >= 1 else 'bubble'} point")
    else:
        raise ConversionError([f"{u.kind} '{u.tag}': calculation mode '{mode}' is not supported"])
    dq = ctx.f(u, "DeltaQ")
    if math.isfinite(dq):
        p.duty = {"unit": u.key, "var": "heat_duty", "dwsim": dq * 1000.0, "sign": sign, "efficiency": eff}
        if mode == "EnergyStream" and cooler:
            p.duty["sign"] = 1.0
    return p


def _cooler(ctx: _Ctx, u: Unit) -> UnitPlan:
    return _heater(ctx, u, cooler=True)


def _valve(ctx: _Ctx, u: Unit) -> UnitPlan:
    _require(ctx, u, 1, 1)
    ctx.imports.update({"PressureChanger", "ThermodynamicAssumption"})
    p = UnitPlan(u.key, u.tag, u.kind, idaes="PressureChanger (isenthalpic)")
    p.build.append(
        f"m.fs.{u.key} = PressureChanger(property_package={ctx.unit_props(u)}, compressor=False,\n"
        f"        thermodynamic_assumption=ThermodynamicAssumption.adiabatic)")
    p.in_ports, p.out_ports, p.init = {0: f"{u.key}.inlet"}, {0: f"{u.key}.outlet"}, [u.key]
    mode = ctx.s(u, "CalcMode", "OutletPressure")
    out = ctx.out_stream(u)
    if mode == "DeltaP":
        p.specs.append(f"u.deltaP.fix({_num(-ctx.f(u, 'DeltaP', 0.0))})  # Pa")
    else:
        pout = ctx.f(u, "OutletPressure", out.pressure)
        if mode != "OutletPressure":
            pout = out.pressure
            ctx.note(f"valve '{u.tag}' uses a Kv/Cv model ({mode}); the outlet pressure is fixed at the "
                     "DWSIM result instead")
        p.specs.append(f"u.outlet.pressure.fix({_num(pout)})  # Pa")
    return p


def _pump(ctx: _Ctx, u: Unit) -> UnitPlan:
    _require(ctx, u, 1, 1)
    ctx.imports.add("Pump")
    p = UnitPlan(u.key, u.tag, u.kind, idaes="Pump")
    p.build.append(f"m.fs.{u.key} = Pump(property_package={ctx.unit_props(u)})")
    p.in_ports, p.out_ports, p.init = {0: f"{u.key}.inlet"}, {0: f"{u.key}.outlet"}, [u.key]
    eff = ctx.f(u, "Eficiencia", 75.0) / 100.0 or 0.75
    mode = ctx.s(u, "CalcMode", "Delta_P")
    out, inl = ctx.out_stream(u), ctx.in_stream(u)
    p.specs.append(f"u.efficiency_pump.fix({_num(eff)})")
    if mode == "Delta_P":
        p.specs.append(f"u.deltaP.fix({_num(ctx.f(u, 'DeltaP', out.pressure - inl.pressure))})  # Pa")
    elif mode == "OutletPressure":
        p.specs.append(f"u.outlet.pressure.fix({_num(ctx.f(u, 'Pout', out.pressure))})  # Pa")
    else:
        p.specs.append(f"u.outlet.pressure.fix({_num(out.pressure)})  # Pa")
        ctx.note(f"pump '{u.tag}' is specified by {mode}; the outlet pressure is fixed at the DWSIM "
                 "result instead")
    dq = ctx.f(u, "DeltaQ")
    if math.isfinite(dq):
        p.duty = {"unit": u.key, "var": "work_mechanical", "dwsim": dq * 1000.0}
    return p


def _compressor(ctx: _Ctx, u: Unit, expander: bool = False) -> UnitPlan:
    _require(ctx, u, 1, 1)
    ctx.imports.update({"PressureChanger", "ThermodynamicAssumption"})
    what = "expander" if expander else "compressor"
    p = UnitPlan(u.key, u.tag, u.kind, idaes=f"PressureChanger (isentropic {what})")
    p.build.append(
        f"m.fs.{u.key} = PressureChanger(property_package={ctx.unit_props(u)}, compressor={not expander},\n"
        f"        thermodynamic_assumption=ThermodynamicAssumption.isentropic)")
    p.in_ports, p.out_ports, p.init = {0: f"{u.key}.inlet"}, {0: f"{u.key}.outlet"}, [u.key]
    eff = ctx.f(u, "AdiabaticEfficiency", 75.0) / 100.0 or 0.75
    if ctx.s(u, "ProcessPath", "Adiabatic") != "Adiabatic":
        ctx.note(f"{what} '{u.tag}' follows a polytropic path in DWSIM; the converted model uses the "
                 "equivalent adiabatic (isentropic) efficiency stored in the file")
    mode = ctx.s(u, "CalcMode", "OutletPressure")
    out = ctx.out_stream(u)
    sign = -1.0 if expander else 1.0
    p.specs.append(f"u.efficiency_isentropic.fix({_num(eff)})")
    if mode == "OutletPressure":
        p.specs.append(f"u.outlet.pressure.fix({_num(ctx.f(u, 'POut', out.pressure))})  # Pa")
    elif mode == "Delta_P":
        p.specs.append(f"u.deltaP.fix({_num(sign * ctx.f(u, 'DeltaP', 0.0))})  # Pa")
    elif mode in ("PowerRequired", "PowerGenerated", "EnergyStream"):
        w = ctx.f(u, "DeltaQ", 0.0) * 1000.0
        if mode == "EnergyStream" and ctx.energy_power(u) is not None:
            w = ctx.energy_power(u)
        p.specs.append(f"u.work_mechanical.fix({_num(sign * w)})  # W")
    else:
        p.specs.append(f"u.outlet.pressure.fix({_num(out.pressure)})  # Pa")
        ctx.note(f"{what} '{u.tag}' is specified by {mode}; the outlet pressure is fixed at the DWSIM "
                 "result instead")
    dq = ctx.f(u, "DeltaQ")
    if math.isfinite(dq):
        # DWSIM stores the power an expander generates, IDAES the work done on the fluid
        p.duty = {"unit": u.key, "var": "work_mechanical", "dwsim": dq * 1000.0, "sign": sign}
    return p


def _expander(ctx: _Ctx, u: Unit) -> UnitPlan:
    return _compressor(ctx, u, expander=True)


def _mixer_build(ctx: _Ctx, key: str, n_in: int, behaviour: str, props: str) -> Tuple[List[str], List[str], List[str]]:
    """Return (build lines, spec lines, inlet names) for a Mixer."""
    ctx.imports.update({"Mixer", "MomentumMixingType"})
    names = [f"inlet_{i + 1}" for i in range(n_in)]
    specs: List[str] = []
    if n_in == 1:
        mom = "equality"
    elif behaviour == "Minimum":
        mom = "minimize"
    else:
        mom = "none"
        specs.append(f"rt.set_mixer_pressure(m.fs.{key}, {behaviour.lower()!r})")
    build = [f"m.fs.{key} = Mixer(property_package={props}, inlet_list={names!r},\n"
             f"        momentum_mixing_type=MomentumMixingType.{mom})"]
    return build, specs, names


def _mixer(ctx: _Ctx, u: Unit) -> UnitPlan:
    _require(ctx, u, None, 1)
    if not u.inlets:
        raise ConversionError([f"Mixer '{u.tag}' has no inlet"])
    p = UnitPlan(u.key, u.tag, u.kind, idaes="Mixer")
    behaviour = ctx.s(u, "PressureCalculation", "Minimum")
    idxs = sorted(u.inlets)
    build, specs, names = _mixer_build(ctx, u.key, len(idxs), behaviour, ctx.unit_props(u))
    p.build, p.specs = build, [s.replace(f"m.fs.{u.key}", "u") for s in specs]
    p.in_ports = {i: f"{u.key}.{n}" for i, n in zip(idxs, names)}
    p.out_ports, p.init = {0: f"{u.key}.outlet"}, [u.key]
    return p


def _splitter(ctx: _Ctx, u: Unit) -> UnitPlan:
    _require(ctx, u, 1, None)
    idxs = sorted(u.outlets)
    if not idxs:
        raise ConversionError([f"Splitter '{u.tag}' has no outlet"])
    ctx.imports.update({"Separator", "SplittingType"})
    p = UnitPlan(u.key, u.tag, u.kind, idaes="Separator (total flow split)")
    names = [f"outlet_{i + 1}" for i in idxs]
    p.build.append(
        f"m.fs.{u.key} = Separator(property_package={ctx.unit_props(u)}, outlet_list={names!r},\n"
        f"        split_basis=SplittingType.totalFlow, ideal_separation=False)")
    p.in_ports, p.init = {0: f"{u.key}.inlet"}, [u.key]
    p.out_ports = {i: f"{u.key}.{n}" for i, n in zip(idxs, names)}
    mode = ctx.s(u, "OperationMode", "SplitRatios")
    inl = ctx.in_stream(u)
    ratios = []
    raw = u.params.get("Ratios")
    if isinstance(raw, list):
        for r in raw:
            try:
                ratios.append(float(r))
            except (TypeError, ValueError):
                ratios.append(float("nan"))
    for k, (i, n) in enumerate(zip(idxs[:-1], names[:-1])):
        out = ctx.stream(u.outlets[i])
        if mode == "SplitRatios" and i < len(ratios) and math.isfinite(ratios[i]):
            p.specs.append(f"u.split_fraction[0, {n!r}].fix({_num(ratios[i])})")
        elif mode == "StreamMassFlowSpec":
            w = ctx.f(u, "StreamFlowSpec" if k == 0 else "Stream2FlowSpec", out.flow_mass)
            p.specs.append(f"rt.fix_outlet_mass_flow(u, {n!r}, {_num(w)})  # kg/s")
        elif mode == "StreamMoleFlowSpec":
            w = ctx.f(u, "StreamFlowSpec" if k == 0 else "Stream2FlowSpec", out.flow_mol)
            p.specs.append(f"u.{n}.flow_mol.fix({_num(w)})  # mol/s")
        else:
            frac = out.flow_mol / inl.flow_mol if inl.flow_mol else 0.0
            p.specs.append(f"u.split_fraction[0, {n!r}].fix({_num(frac)})  # from the DWSIM result")
    return p


def _vessel(ctx: _Ctx, u: Unit) -> UnitPlan:
    ctx.imports.add("Flash")
    idxs = sorted(u.inlets)
    if not idxs:
        raise ConversionError([f"Vessel '{u.tag}' has no inlet"])
    liq2 = ctx.out_stream(u, 2)
    if liq2 is not None and liq2.flow_mol > 1e-9 * max(ctx.stream(u.inlets[idxs[0]]).flow_mol, 1e-30):
        raise ConversionError([f"Vessel '{u.tag}' separates a second liquid phase (three-phase "
                               "separation is not supported yet)"])
    p = UnitPlan(u.key, u.tag, u.kind, idaes="Flash")
    ctx.unit_props(u, force_vle=True)
    p.build.append(f"m.fs.{u.key} = Flash(property_package=m.fs.props, has_heat_transfer=True,\n"
                   f"        has_pressure_change=True)")
    behaviour = ctx.s(u, "PressureCalculation", "Minimum")
    if len(idxs) > 1:
        mk = u.key + "_mix"
        build, specs, names = _mixer_build(ctx, mk, len(idxs), behaviour, "m.fs.props")
        p.build = build + p.build
        p.specs += specs
        p.in_ports = {i: f"{mk}.{n}" for i, n in zip(idxs, names)}
        p.internal_arcs.append((u.key + "_feed", f"{mk}.outlet", f"{u.key}.inlet"))
        p.init = [mk, u.key]
        p.idaes = "Mixer + Flash"
    else:
        p.in_ports = {idxs[0]: f"{u.key}.inlet"}
        p.init = [u.key]
    p.out_ports = {0: f"{u.key}.vap_outlet", 1: f"{u.key}.liq_outlet"}
    for i in u.outlets:
        if i not in (0, 1) and not (i == 2 and liq2 is not None):
            raise ConversionError([f"Vessel '{u.tag}': unexpected outlet connector {i}"])
    vap, liq = ctx.out_stream(u, 0), ctx.out_stream(u, 1)
    p.internal[f"{u.key}.control_volume.properties_out"] = {
        "inlets": [u.inlets[i] for i in idxs], "vapor": vap.id if vap else None,
        "liquid": liq.id if liq else None}
    over_t = ctx.s(u, "OverrideT", "false").lower() == "true"
    over_p = ctx.s(u, "OverrideP", "false").lower() == "true"
    mode = ctx.s(u, "CalculationMode", "Legacy")
    ref = ctx.out_stream(u, 0) or ctx.out_stream(u, 1)
    if mode in ("HeatingCoolingIsothermic", "HeatingCoolingIsobaric"):
        ctx.note(f"vessel '{u.tag}' uses mode {mode}; the flash temperature is fixed at the DWSIM result")
        over_t = True
    if over_t:
        t = ctx.f(u, "FlashTemperature", ref.temperature) if mode == "Legacy" else ref.temperature
        p.specs.append(f"u.control_volume.properties_out[0].temperature.fix({_num(t)})  # K, specified flash temperature")
    else:
        p.specs.append("u.heat_duty.fix(0.0)  # adiabatic flash")
    if over_p:
        p.specs.append(f"u.control_volume.properties_out[0].pressure.fix({_num(ctx.f(u, 'FlashPressure', ref.pressure))})  # Pa")
    else:
        p.specs.append("u.deltaP.fix(0.0)")
    return p


def _heat_exchanger(ctx: _Ctx, u: Unit) -> UnitPlan:
    _require(ctx, u, 2, 2)
    if set(u.inlets) != {0, 1} or set(u.outlets) != {0, 1}:
        raise ConversionError([f"HeatExchanger '{u.tag}': both sides must be connected"])
    ctx.imports.update({"HeatExchanger", "HeatExchangerFlowPattern", "delta_temperature_underwood_callback"})
    s_in = {i: ctx.stream(u.inlets[i]) for i in (0, 1)}
    s_out = {i: ctx.stream(u.outlets[i]) for i in (0, 1)}
    hot = 0 if s_in[0].temperature >= s_in[1].temperature else 1
    cold = 1 - hot
    co = ctx.s(u, "FlowDir", "CounterCurrent") == "CoCurrent"
    p = UnitPlan(u.key, u.tag, u.kind, idaes="HeatExchanger (0D, LMTD)")
    p.build.append(
        f"m.fs.{u.key} = HeatExchanger(\n"
        f"        hot_side={{\"property_package\": {ctx.props([s_in[hot], s_out[hot]])}, \"has_pressure_change\": True}},\n"
        f"        cold_side={{\"property_package\": {ctx.props([s_in[cold], s_out[cold]])}, \"has_pressure_change\": True}},\n"
        f"        delta_temperature_callback=delta_temperature_underwood_callback,\n"
        f"        flow_pattern=HeatExchangerFlowPattern.{'cocurrent' if co else 'countercurrent'})")
    p.in_ports = {hot: f"{u.key}.hot_side_inlet", cold: f"{u.key}.cold_side_inlet"}
    p.out_ports = {hot: f"{u.key}.hot_side_outlet", cold: f"{u.key}.cold_side_outlet"}
    p.init = [u.key]
    th1, th2, tc1, tc2 = s_in[hot].temperature, s_out[hot].temperature, s_in[cold].temperature, s_out[cold].temperature
    d1, d2 = (th1 - tc1, th2 - tc2) if co else (th1 - tc2, th2 - tc1)
    if not (d1 > 0 and d2 > 0):
        raise ConversionError([f"HeatExchanger '{u.tag}': the DWSIM result has a temperature cross "
                               "(no valid LMTD), which the IDAES 0D exchanger cannot represent"])
    dph, dpc = ctx.f(u, "HotSidePressureDrop", 0.0), ctx.f(u, "ColdSidePressureDrop", 0.0)
    sp = p.specs
    sp.append(f"u.hot_side.deltaP.fix({_num(-dph)})  # Pa")
    sp.append(f"u.cold_side.deltaP.fix({_num(-dpc)})  # Pa")
    area = ctx.f(u, "Area", 1.0)
    f_lmtd = ctx.f(u, "LMTD_F", 1.0) or 1.0
    ucoef = ctx.f(u, "OverallCoefficient", 0.0) * f_lmtd
    q = ctx.f(u, "Q", 0.0) * 1000.0
    mode = ctx.s(u, "CalculationMode", "CalcTempHotOut")
    defined = ctx.s(u, "DefinedTemperature", "Cold_Fluid")
    fix_tc = f"u.cold_side_outlet.temperature.fix({_num(ctx.f(u, 'ColdSideOutletTemperature', tc2))})  # K"
    fix_th = f"u.hot_side_outlet.temperature.fix({_num(ctx.f(u, 'HotSideOutletTemperature', th2))})  # K"
    fix_a = f"u.area.fix({_num(area)})  # m2"
    fix_u = f"u.overall_heat_transfer_coefficient.fix({_num(ucoef)})  # W/m2/K"
    if mode == "CalcTempHotOut":
        sp += [fix_tc, fix_a]
    elif mode == "CalcTempColdOut":
        sp += [fix_th, fix_a]
    elif mode == "CalcBothTemp":
        sp += [f"u.heat_duty.fix({_num(q)})  # W", fix_a]
    elif mode == "CalcBothTemp_UA":
        sp += [fix_u, fix_a]
    elif mode == "CalcArea":
        sp += [fix_u, fix_th if defined == "Hot_Fluid" else fix_tc]
    else:
        sp += [f"u.cold_side_outlet.temperature.fix({_num(tc2)})  # K, DWSIM result", fix_a]
        ctx.note(f"heat exchanger '{u.tag}' is specified by {mode}; the cold outlet temperature is "
                 "fixed at the DWSIM result instead")
    if f_lmtd != 1.0 and "overall_heat_transfer_coefficient" in "\n".join(sp):
        ctx.note(f"heat exchanger '{u.tag}': the LMTD correction factor {f_lmtd:g} is folded into U")
    loss = ctx.f(u, "HeatLoss", 0.0) * 1000.0
    if not math.isfinite(loss):
        loss = 0.0
    if loss != 0.0:
        # DWSIM: the energy lost to the surroundings leaves the exchanger, so
        # hot-side heat + cold-side heat = -loss instead of 0
        p.build += [
            f"m.fs.{u.key}.heat_loss = pyo.Var(m.fs.time, initialize={_num(loss)}, units=pyo.units.W)",
            f"m.fs.{u.key}.heat_loss.fix()",
            f"m.fs.{u.key}.unit_heat_balance.deactivate()",
            f"m.fs.{u.key}.unit_heat_balance_with_loss = pyo.Constraint(m.fs.time, rule=lambda b, t:\n"
            f"        0 == b.hot_side.heat[t] + b.cold_side.heat[t] + b.heat_loss[t])"]
        ctx.note(f"heat exchanger '{u.tag}': the heat loss of {loss / 1000:g} kW is kept as a fixed "
                 "loss to the surroundings")
    # IDAES heat_duty is the heat the cold side receives
    q_cold = q - loss if mode == "CalcTempColdOut" else q
    p.duty = {"unit": u.key, "var": "heat_duty", "dwsim": q_cold}
    return p


HANDLERS: Dict[str, Callable[[_Ctx, Unit], UnitPlan]] = {
    "Heater": _heater,
    "Cooler": _cooler,
    "Valve": _valve,
    "Pump": _pump,
    "Compressor": _compressor,
    "Expander": _expander,
    "Mixer": _mixer,
    "Splitter": _splitter,
    "Vessel": _vessel,
    "HeatExchanger": _heat_exchanger,
}

SUPPORTED_UNITS = sorted(HANDLERS)


# ---------------------------------------------------------------------------
# conversion
# ---------------------------------------------------------------------------
def _blocker(fs: Flowsheet, u: Unit) -> Optional[str]:
    """Why this particular unit cannot be converted although its kind is supported."""
    if str(u.params.get("OverrideCalculationRoutine", "")).lower() == "true":
        return "a custom script replaces its calculation routine in DWSIM"
    ctx = _Ctx(fs)
    ins = [fs.streams[i] for i in u.inlets.values() if i in fs.streams]
    outs = [fs.streams[i] for i in u.outlets.values() if i in fs.streams]
    labels = {ctx.label(s) for s in ins + outs} - {""}
    if u.kind == "Pump" and ins and ctx.label(ins[0]) in ("V", "2"):
        return ("its inlet contains vapour; DWSIM then applies the incompressible-liquid formula "
                "to a gas, which the IDAES Pump does not reproduce")
    present = {n for s in ins + outs if ctx.label(s) for n, x in s.z.items() if x > 1e-12}
    if len(present) == 1 and ("2" in labels or {"V", "L"} <= labels):
        return ("a pure fluid changes phase across it; temperature and pressure do not fix the phase "
                "split of a pure fluid in the IDAES vapour-liquid formulation")
    if u.kind in ("Compressor", "Expander") and ins and outs:
        s_in = ins[0]
        if (ctx.label(s_in) in ("V", "L") and ctx.label(outs[0]) == "2"
                and s_in.pressure > ctx.pseudo_critical(s_in)[1]):
            return ("a dense single-phase inlet above the critical pressure expands into two phases; "
                    "the inlet entropy is not unique in the IDAES vapour-liquid formulation, so the "
                    "isentropic calculation does not converge")
    return None


def analyze(fs: Flowsheet) -> dict:
    """Classify a flowsheet without converting it (used by ``inspect``)."""
    kinds: Dict[str, int] = {}
    for u in fs.units.values():
        kinds[u.kind] = kinds.get(u.kind, 0) + 1
    unsupported = sorted(k for k in kinds if k not in HANDLERS and k not in PASSTHROUGH and k not in LOGICAL)
    pk_ids = {u.package_id for u in fs.units.values() if u.kind in HANDLERS and u.package_id}
    pk_ids |= {s.package_id for s in fs.streams.values() if s.package_id}
    pk_kinds = sorted({fs.packages[i].kind for i in pk_ids if i in fs.packages}) or \
        sorted({p.kind for p in list(fs.packages.values())[:1]})
    bad_pk = [k for k in pk_kinds if k not in SUPPORTED_PACKAGES]
    return {
        "unit_kinds": kinds,
        "unsupported_units": unsupported,
        "package_kinds": pk_kinds,
        "unsupported_packages": bad_pk,
        "multiple_packages": len(pk_kinds) > 1,
        "convertible": not unsupported and not bad_pk and len(pk_kinds) == 1,
    }


def convert(fs: Flowsheet, cut_unsupported: bool = False, phases: str = "auto") -> Conversion:
    """Convert ``fs``; raise :class:`ConversionError` if that is not possible.

    With ``cut_unsupported=True`` unit operations that have no IDAES mapping
    are removed and the streams they produce become fixed feeds at their DWSIM
    values, so the rest of the flowsheet can still be converted.

    ``phases="auto"`` (default) gives every unit the smallest phase set that
    covers what DWSIM found around it (vapour only, liquid only or VLE);
    ``phases="vle"`` uses the full vapour-liquid package everywhere.
    """
    if phases not in ("auto", "vle"):
        raise ValueError("phases must be 'auto' or 'vle'")
    ctx = _Ctx(fs, phases)
    problems: List[str] = []

    blocked = {uid: _blocker(fs, u) for uid, u in fs.units.items() if u.kind in HANDLERS}
    converted = {uid: u for uid, u in fs.units.items() if u.kind in HANDLERS and not blocked[uid]}
    ignored = [u for u in fs.units.values() if u.kind in LOGICAL]
    unsupported = [u for u in fs.units.values()
                   if u.kind not in HANDLERS and u.kind not in PASSTHROUGH and u.kind not in LOGICAL]
    if unsupported and not cut_unsupported:
        kinds = sorted({u.kind for u in unsupported})
        problems.append("unsupported unit operations: " + ", ".join(
            f"{k} ({', '.join(u.tag for u in unsupported if u.kind == k)})" for k in kinds))
    for u in unsupported:
        ctx.note(f"{u.kind} '{u.tag}' is not supported and was cut out: its outlet streams are "
                 "feeds fixed at the DWSIM values")
    for uid, reason in blocked.items():
        if not reason:
            continue
        u = fs.units[uid]
        if not cut_unsupported:
            problems.append(f"{u.kind} '{u.tag}' cannot be converted: {reason}")
        ctx.note(f"{u.kind} '{u.tag}' was cut out ({reason}): its outlet streams are feeds fixed at "
                 "the DWSIM values")
        unsupported.append(u)
    for u in ignored:
        if u.kind in ("Adjust", "Spec"):
            ctx.note(f"{u.kind} block '{u.tag}' is not translated; the variable it manipulates keeps "
                     "its converged DWSIM value")
    if not converted:
        problems.append("the flowsheet contains no supported unit operation")
    if problems:
        raise ConversionError(problems)

    def resolve_source(s: Stream, depth: int = 0) -> Optional[Tuple[Unit, int]]:
        if s.source is None or depth > 50:
            return None
        src = fs.units[s.source]
        if src.kind in PASSTHROUGH:
            up = src.inlets.get(0)
            return resolve_source(fs.streams[up], depth + 1) if up else None
        if s.source in converted:
            return src, s.source_port
        return None

    # ---- property package -------------------------------------------------
    active_streams: Dict[str, Stream] = {}
    for u in converted.values():
        for sid in list(u.inlets.values()) + list(u.outlets.values()):
            active_streams[sid] = fs.streams[sid]
    for s in fs.streams.values():
        if resolve_source(s) is not None:
            active_streams[s.id] = s
    pk_ids = [u.package_id for u in converted.values() if u.package_id]
    pk_ids += [s.package_id for s in active_streams.values() if s.package_id]
    pk_ids = [i for i in dict.fromkeys(pk_ids) if i in fs.packages]
    if not pk_ids:
        pk_ids = list(fs.packages)[:1]
    if not pk_ids:
        raise ConversionError(["the file defines no property package"])
    kinds = sorted({fs.packages[i].kind for i in pk_ids})
    if len(kinds) > 1:
        raise ConversionError(["the flowsheet uses several property packages ("
                               + ", ".join(kinds) + "); only one per flowsheet is supported"])
    pkg = fs.packages[pk_ids[0]]
    if len(pk_ids) > 1:
        ctx.note(f"{len(pk_ids)} property package instances of type {pkg.kind} are merged into one")

    for s in active_streams.values():
        if not (math.isfinite(s.temperature) and math.isfinite(s.pressure) and math.isfinite(s.flow_mol)):
            problems.append(f"stream '{s.tag}' has no calculated state in the file (was the flowsheet solved?)")
        if s.liquid2_fraction > 1e-6:
            ctx.note(f"DWSIM found two liquid phases in stream '{s.tag}'; the IDAES model has a single "
                     "liquid phase, so results around it will deviate")
        if s.solid_fraction > 1e-6:
            ctx.note(f"DWSIM found a solid phase in stream '{s.tag}'; solids are not modelled")
    if problems:
        raise ConversionError(problems)

    present = [n for n in fs.compounds
               if any(s.z.get(n, 0.0) > 1e-12 for s in active_streams.values())]
    dropped = [n for n in fs.compounds if n not in present]
    if dropped:
        ctx.note("compounds absent from every converted stream were left out: " + ", ".join(dropped))
    if not present:
        raise ConversionError(["no compound has a non-zero mole fraction in any stream"])
    if len(present) == 1:
        boiling = [s.tag for s in active_streams.values() if ctx.label(s) == "2"]
        if boiling:
            raise ConversionError([
                "single-compound flowsheet with two-phase streams (" + ", ".join(boiling[:4])
                + ("..." if len(boiling) > 4 else "") + "): temperature and pressure cannot define the "
                "state of a boiling pure fluid, and the enthalpy-based state definition this needs is "
                "not implemented yet"])
    try:
        thermo = thermo_spec(fs, pkg, present)
    except ThermoError as e:
        raise ConversionError([str(e)]) from e
    for n in thermo.pop("notes"):
        ctx.note(n)
    if thermo["model"] == "SRK":
        ctx.note("Soave-Redlich-Kwong support is experimental: IPOPT has crashed on some SRK models "
                 "(the same flowsheets solve with Peng-Robinson)")
    key_of = {n: fs.compounds[n].key for n in present}

    # ---- units ------------------------------------------------------------
    plans: Dict[str, UnitPlan] = {}
    for uid, u in list(converted.items()):
        try:
            plans[uid] = HANDLERS[u.kind](ctx, u)
        except ConversionError as e:
            if not cut_unsupported:
                problems.extend(e.problems)
                continue
            del converted[uid]
            unsupported.append(u)
            ctx.note(f"{u.kind} '{u.tag}' was cut out ({'; '.join(e.problems)}): its outlet streams "
                     "are feeds fixed at the DWSIM values")
    if problems:
        raise ConversionError(problems)
    if not converted:
        raise ConversionError(["no unit operation is left to convert"])

    # ---- connectivity -----------------------------------------------------
    arcs: List[Tuple[str, str, str, str]] = []          # (name, src port, dst port, tag)
    feeds: List[Tuple[str, Stream]] = []                # (dst port, stream)
    stream_ports: Dict[str, List[str]] = {}
    feed_ids = set()
    used_arc = set()
    for uid, u in converted.items():
        plan = plans[uid]
        for idx, sid in sorted(u.inlets.items()):
            s = fs.streams[sid]
            dst = plan.in_ports.get(idx)
            if dst is None:
                raise ConversionError([f"{u.kind} '{u.tag}': inlet connector {idx} is not mapped"])
            src = resolve_source(s)
            if src is None:
                feeds.append((dst, s))
                feed_ids.add(sid)
                stream_ports.setdefault(sid, []).append(dst)
                if s.spec_type != "Temperature_and_Pressure" and s.source is None:
                    ctx.note(f"feed '{s.tag}' is specified by {s.spec_type} in DWSIM; the converted "
                             "model fixes its temperature and pressure at the DWSIM values")
            else:
                su, sidx = src
                sp = plans[su.id].out_ports.get(sidx)
                if sp is None:
                    raise ConversionError([f"{su.kind} '{su.tag}': outlet connector {sidx} is not mapped"])
                name = s.key
                while name in used_arc:
                    name += "_x"
                used_arc.add(name)
                arcs.append((name, sp, dst, s.tag))
                stream_ports.setdefault(sid, []).append(dst)
    for s in fs.streams.values():
        src = resolve_source(s)
        if src is not None:
            su, sidx = src
            sp = plans[su.id].out_ports.get(sidx)
            if sp is not None:
                stream_ports.setdefault(s.id, []).insert(0, sp)

    closures = _close_loops(ctx, converted, plans, arcs, feeds)

    streams: Dict[str, dict] = {}
    tags_seen = set()
    for sid, ports in stream_ports.items():
        s = fs.streams[sid]
        tag = s.tag
        while tag in tags_seen:
            tag += "'"
        tags_seen.add(tag)
        entry = {
            "ports": ports,
            "is_feed": sid in feed_ids,
            "temperature": s.temperature, "pressure": s.pressure, "flow_mol": s.flow_mol,
            "vapor_fraction": s.vapor_fraction,
            "x": {key_of[n]: v for n, v in s.z.items() if n in key_of},
        }
        if ctx.supercritical(s):
            entry["supercritical"] = True
        if PHASE_TOL < s.vapor_fraction < 1 - PHASE_TOL:
            entry["y"] = {key_of[n]: v for n, v in s.vapor.x.items() if n in key_of}
            entry["xl"] = {key_of[n]: v for n, v in s.liquid.x.items() if n in key_of}
        streams[tag] = entry

    # states inside units (the two-phase mixture in a flash drum)
    for plan in plans.values():
        for path, d in plan.internal.items():
            ins = [fs.streams[i] for i in d["inlets"]]
            vap = fs.streams[d["vapor"]] if d["vapor"] else None
            liq = fs.streams[d["liquid"]] if d["liquid"] else None
            ftot = sum(s.flow_mol for s in ins)
            if ftot <= 0:
                continue
            z = {n: sum(s.flow_mol * s.z.get(n, 0.0) for s in ins) / ftot for n in key_of}
            ref = vap if (vap is not None and vap.flow_mol > 0) else liq
            fv = vap.flow_mol if vap is not None else 0.0
            entry = {
                "block": path, "temperature": ref.temperature, "pressure": ref.pressure,
                "flow_mol": ftot, "vapor_fraction": min(max(fv / ftot, 0.0), 1.0),
                "x": {key_of[n]: v for n, v in z.items()},
            }
            if vap is not None and liq is not None and vap.flow_mol > 0 and liq.flow_mol > 0:
                entry["y"] = {key_of[n]: v for n, v in vap.z.items() if n in key_of}
                entry["xl"] = {key_of[n]: v for n, v in liq.z.items() if n in key_of}
            streams["(inside) " + plan.tag] = entry

    duties = {p.tag: p.duty for p in plans.values() if p.duty}
    code = _emit(ctx, fs, thermo, list(plans.values()), arcs, feeds, streams, duties, key_of, closures)
    return Conversion(code=code, flowsheet=fs, thermo=thermo, streams=streams,
                      units=list(plans.values()), notes=ctx.notes,
                      cut_units=[u.tag for u in unsupported], ignored_units=[u.tag for u in ignored])


def _close_loops(ctx, converted, plans, arcs, feeds) -> List[Tuple[str, Stream]]:
    """Find groups of units that no feed reaches (closed loops, e.g. a refrigerant
    cycle) and pick one connection in each where the amount and composition are fixed."""
    port_unit = {}
    for uid, plan in plans.items():
        for port in list(plan.in_ports.values()) + list(plan.out_ports.values()):
            port_unit[port] = uid
    succ: Dict[str, set] = {uid: set() for uid in converted}
    for _, sp, dp, _ in arcs:
        su, du = port_unit.get(sp), port_unit.get(dp)
        if su and du:
            succ[su].add(du)
    reached = set()
    todo = [port_unit[d] for d, _ in feeds if d in port_unit]
    while todo:
        u = todo.pop()
        if u in reached:
            continue
        reached.add(u)
        todo.extend(succ[u] - reached)
    out: List[Tuple[str, Stream]] = []
    closed = set()
    by_tag = {s.tag: s for s in ctx.fs.streams.values()}
    for name, sp, dp, tag in arcs:
        su, du = port_unit.get(sp), port_unit.get(dp)
        if su in reached or su in closed or du in closed or su is None:
            continue
        # close the loop at this arc; everything downstream of it is now reached
        out.append((name, by_tag[tag]))
        stack = [du]
        while stack:
            u = stack.pop()
            if u in closed:
                continue
            closed.add(u)
            stack.extend(succ[u] - closed)
        ctx.note(f"the units around stream '{tag}' form a closed loop without a feed; its flow and "
                 "composition are fixed at the DWSIM values")
    return out


# ---------------------------------------------------------------------------
# code emission
# ---------------------------------------------------------------------------
_IMPORT_LINES = {
    "Heater": "from idaes.models.unit_models import Heater",
    "Mixer": "from idaes.models.unit_models import Mixer",
    "MomentumMixingType": "from idaes.models.unit_models import MomentumMixingType",
    "Separator": "from idaes.models.unit_models import Separator",
    "SplittingType": "from idaes.models.unit_models.separator import SplittingType",
    "Flash": "from idaes.models.unit_models import Flash",
    "Pump": "from idaes.models.unit_models import Pump",
    "PressureChanger": "from idaes.models.unit_models import PressureChanger",
    "ThermodynamicAssumption": "from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption",
    "HeatExchanger": "from idaes.models.unit_models import HeatExchanger",
    "HeatExchangerFlowPattern": "from idaes.models.unit_models.heat_exchanger import HeatExchangerFlowPattern",
    "delta_temperature_underwood_callback":
        "from idaes.models.unit_models.heat_exchanger import delta_temperature_underwood_callback",
}


def _literal(name: str, obj) -> str:
    return f"{name} = " + pprint.pformat(obj, width=100, sort_dicts=False, compact=True)


def _emit(ctx, fs, thermo, plans, arcs, feeds, streams, duties, key_of, closures=()) -> str:
    L: List[str] = []
    a = L.append
    a('"""IDAES model of the DWSIM flowsheet %r.' % fs.source_file)
    a("")
    a(f"Generated by dwsim2idaes {__version__}. Units are SI (K, Pa, mol/s, W).")
    a("Run this file to build, initialise and solve the model and to compare the")
    a("result with the DWSIM solution stored in the original file; import it and")
    a("call ``build_model()`` to use the flowsheet in your own optimisation study.")
    if ctx.notes:
        a("")
        a("Conversion notes:")
        for n in ctx.notes:
            a("  * " + n)
    a('"""')
    a("import logging")
    a("")
    a("import pyomo.environ as pyo")
    a("from pyomo.network import Arc")
    a("")
    a("from idaes.core import FlowsheetBlock")
    a("from idaes.models.properties.modular_properties import GenericParameterBlock")
    for k in sorted(ctx.imports):
        a(_IMPORT_LINES[k])
    a("")
    a("from dwsim2idaes import runtime as rt")
    a("from dwsim2idaes.thermo import build_property_config")
    a("")
    a("# Property data read from the DWSIM file (compound constants, correlations, kij).")
    a(_literal("THERMO", thermo))
    a("")
    a("# Converged DWSIM stream results: used as the initial guess and as the reference")
    a("# for validation. 'ports' is where the stream lives in the IDAES flowsheet")
    a("# (first the port that produces it, then the ports that consume it).")
    a(_literal("STREAMS", streams))
    a("")
    a("# Duties and powers reported by DWSIM (W), for validation.")
    a(_literal("DUTIES", duties))
    a("")
    a(_literal("UNITS", [n for p in plans for n in p.init]))
    a("")
    a("")
    a("def build_model():")
    a('    """Build the flowsheet: property package, unit models and connections."""')
    a(f"    m = pyo.ConcreteModel(name={fs.name!r})")
    a("    m.fs = FlowsheetBlock(dynamic=False)")
    sets = ctx.phase_sets or {"VLE"}
    if "VLE" in sets:
        a("    m.fs.props = GenericParameterBlock(**build_property_config(THERMO))")
    if "Vap" in sets:
        a("    # vapour-only variant for units where DWSIM found no liquid")
        a('    m.fs.props_vap = GenericParameterBlock(**build_property_config(THERMO, phases="Vap"))')
    if "Liq" in sets:
        a("    # liquid-only variant for units where DWSIM found no vapour")
        a('    m.fs.props_liq = GenericParameterBlock(**build_property_config(THERMO, phases="Liq"))')
    a("")
    for p in plans:
        a(f"    # {p.kind} {p.tag!r}")
        for line in p.build:
            a("    " + line)
    a("")
    a("    # connections (one Arc per DWSIM material stream; recycle blocks vanish)")
    for p in plans:
        for name, src, dst in p.internal_arcs:
            a(f"    m.fs.{name} = Arc(source=m.fs.{src}, destination=m.fs.{dst})")
    for name, src, dst, tag in arcs:
        a(f"    m.fs.{name} = Arc(source=m.fs.{src}, destination=m.fs.{dst})  # {tag!r}")
    a('    pyo.TransformationFactory("network.expand_arcs").apply_to(m)')
    a("    return m")
    a("")
    a("")
    a("def set_specifications(m):")
    a('    """Fix the same degrees of freedom that are specified in DWSIM."""')
    a("    # feed streams")
    for dst, s in feeds:
        a(f"    rt.fix_feed(m.fs.{dst},  # {s.tag!r}")
        a(f"                flow_mol={_num(s.flow_mol)}, temperature={_num(s.temperature)}, "
          f"pressure={_num(s.pressure)},")
        x = {key_of[n]: float(f"{v:.12g}") for n, v in s.z.items() if n in key_of}
        a(f"                mole_frac={x!r})")
    if closures:
        a("")
        a("    # closed loops without a feed: the circulating amount is not fixed by the")
        a("    # steady-state balances, so it is set to the DWSIM value at one stream")
    for arc, s in closures:
        x = {key_of[n]: float(f"{v:.12g}") for n, v in s.z.items() if n in key_of}
        a(f"    rt.close_loop(m.fs.{arc}, flow_mol={_num(s.flow_mol)}, mole_frac={x!r})  # {s.tag!r}")
    for p in plans:
        if not p.specs:
            continue
        a("")
        a(f"    # {p.kind} {p.tag!r}")
        a(f"    u = m.fs.{p.key}")
        for line in p.specs:
            a("    " + line)
    a("")
    a("")
    a("def initialize(m, outlvl=logging.WARNING):")
    a('    """Initialise every unit at its DWSIM inlet state (no tear streams needed)."""')
    a("    return rt.initialize(m, STREAMS, UNITS, outlvl=outlvl)")
    a("")
    a("")
    a("def build_and_solve(tee=False):")
    a("    m = build_model()")
    a("    set_specifications(m)")
    a("    rt.check_dof(m)")
    a("    initialize(m)")
    a("    results = rt.solve(m, tee=tee)")
    a("    return m, results")
    a("")
    a("")
    a('if __name__ == "__main__":')
    a("    m, results = build_and_solve()")
    a('    print("solver:", results.solver.termination_condition)')
    a("    print(rt.format_comparison(rt.compare(m, STREAMS)))")
    a("    print()")
    a("    print(rt.format_duties(rt.unit_duties(m, DUTIES)))")
    return "\n".join(L) + "\n"
