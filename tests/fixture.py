"""Write small, synthetic DWSIM flowsheet files for the tests.

The files contain only the XML elements dwsim2idaes reads. Compound constants
are the ChemSep values DWSIM ships for methane and propane. The stream results
are not produced by DWSIM: they are the converged values of the converted
IDAES model (rounded), so the solve test is a regression test of the whole
pipeline, while the comparisons against real DWSIM results live in
``validation/`` and need the DWSIM sample files.
"""
from __future__ import annotations

import zipfile
import xml.etree.ElementTree as ET

COMPOUNDS = {
    "Methane": {
        "CAS_Number": "74-82-8", "Formula": "CH4", "OriginalDB": "ChemSep", "CurrentDB": "ChemSep",
        "Molar_Weight": "16.04246", "Critical_Temperature": "190.56", "Critical_Pressure": "4599000",
        "Acentric_Factor": "0.011", "Z_Rackett": "0.2876", "Normal_Boiling_Point": "111.66",
        "IsPF": "0", "IsHYPO": "0",
        "IdealgasCpEquation": "16", "VaporPressureEquation": "101",
        "VaporizationEnthalpyEquation": "106", "LiquidDensityEquation": "105",
        "Vapor_Pressure_TMIN": "83.65", "Vapor_Pressure_TMAX": "191.03",
        "HVap_TMIN": "90.67", "HVap_TMAX": "190.55",
        "Liquid_Density_Tmin": "90.68", "Liquid_Density_Tmax": "191.04",
        "Ideal_Gas_Heat_Capacity_Const_A": "33151.9", "Ideal_Gas_Heat_Capacity_Const_B": "-1220.001",
        "Ideal_Gas_Heat_Capacity_Const_C": "12.0907", "Ideal_Gas_Heat_Capacity_Const_D": "-0.000384791",
        "Ideal_Gas_Heat_Capacity_Const_E": "9.896403E-08",
        "Vapor_Pressure_Constant_A": "39.98844", "Vapor_Pressure_Constant_B": "-1337.308",
        "Vapor_Pressure_Constant_C": "-3.580049", "Vapor_Pressure_Constant_D": "3.20698E-05",
        "Vapor_Pressure_Constant_E": "2",
        "HVap_A": "14418000", "HVap_B": "2.3055", "HVap_C": "-5.4199", "HVap_D": "5.658",
        "HVap_E": "-2.1286",
        "Liquid_Density_Const_A": "1.894", "Liquid_Density_Const_B": "0.23603",
        "Liquid_Density_Const_C": "191.05", "Liquid_Density_Const_D": "0.21974",
        "Liquid_Density_Const_E": "0",
    },
    "Propane": {
        "CAS_Number": "74-98-6", "Formula": "CH3CH2CH3", "OriginalDB": "ChemSep", "CurrentDB": "ChemSep",
        "Molar_Weight": "44.09562", "Critical_Temperature": "369.83", "Critical_Pressure": "4248000",
        "Acentric_Factor": "0.152", "Z_Rackett": "0.2763", "Normal_Boiling_Point": "231.02",
        "IsPF": "0", "IsHYPO": "0",
        "IdealgasCpEquation": "16", "VaporPressureEquation": "101",
        "VaporizationEnthalpyEquation": "106", "LiquidDensityEquation": "105",
        "Vapor_Pressure_TMIN": "85.47", "Vapor_Pressure_TMAX": "371.25",
        "HVap_TMIN": "85.44", "HVap_TMAX": "369.85",
        "Liquid_Density_Tmin": "85.47", "Liquid_Density_Tmax": "369.85",
        "Ideal_Gas_Heat_Capacity_Const_A": "37840.4", "Ideal_Gas_Heat_Capacity_Const_B": "-445.5789",
        "Ideal_Gas_Heat_Capacity_Const_C": "11.83871", "Ideal_Gas_Heat_Capacity_Const_D": "0.000653764",
        "Ideal_Gas_Heat_Capacity_Const_E": "-2.200137E-07",
        "Vapor_Pressure_Constant_A": "55.2725", "Vapor_Pressure_Constant_B": "-3396.946",
        "Vapor_Pressure_Constant_C": "-5.423393", "Vapor_Pressure_Constant_D": "8.953731E-06",
        "Vapor_Pressure_Constant_E": "2",
        "HVap_A": "30459000", "HVap_B": "1.2001", "HVap_C": "-2.1107", "HVap_D": "1.9732",
        "HVap_E": "-0.65316",
        "Liquid_Density_Const_A": "1.3186", "Liquid_Density_Const_B": "0.27005",
        "Liquid_Density_Const_C": "369.86", "Liquid_Density_Const_D": "0.27852",
        "Liquid_Density_Const_E": "0",
    },
}

PACKAGE_ID = "PP-test"


def _sub(parent, tag, text=None, **attrib):
    e = ET.SubElement(parent, tag, {k: str(v) for k, v in attrib.items()})
    if text is not None:
        e.text = str(text)
    return e


def _phase(parent, pid, props, x):
    ph = _sub(parent, "Phase")
    _sub(ph, "ID", pid)
    pr = _sub(ph, "Properties")
    for k, v in props.items():
        _sub(pr, k, repr(float(v)))
    comps = _sub(ph, "Compounds")
    for n, v in x.items():
        c = _sub(comps, "Compound")
        _sub(c, "Name", n)
        _sub(c, "MoleFraction", repr(float(v)))


def stream(sid, tag, T, P, F, z, vf, h=0.0, y=None, x=None, spec="Temperature_and_Pressure"):
    """A material stream record (T in K, P in Pa, F in mol/s, h in kJ/kg)."""
    return dict(id=sid, tag=tag, T=T, P=P, F=F, z=z, vf=vf, h=h, y=y or z, x=x or z, spec=spec)


def unit(uid, tag, cls, inlets, outlets, **params):
    """A unit operation; ``cls`` is the DWSIM class name, e.g. 'Heater'."""
    return dict(id=uid, tag=tag, cls=cls, inlets=inlets, outlets=outlets, params=params)


def build_xml(streams, units, compounds=("Methane", "Propane"), kij=None, package="PengRobinson"):
    root = ET.Element("DWSIM_Simulation_Data")
    gi = _sub(root, "GeneralInfo")
    _sub(gi, "BuildVersion", "8.8.3.0")
    sim = _sub(root, "SimulationObjects")
    gfx = _sub(root, "GraphicObjects")
    pps = _sub(root, "PropertyPackages")
    pp = _sub(pps, "PropertyPackage")
    _sub(pp, "ID", PACKAGE_ID)
    _sub(pp, "Type", f"DWSIM.Thermodynamics.PropertyPackages.{package}PropertyPackage")
    _sub(pp, "Tag", package)
    ips = _sub(pp, "InteractionParameters")
    for (a, b), v in (kij or {}).items():
        _sub(ips, "InteractionParameter", Compound1=a, Compound2=b, Value=v)
    cs = _sub(root, "Compounds")
    for n in compounds:
        c = _sub(cs, "Compound")
        _sub(c, "Name", n)
        for k, v in COMPOUNDS[n].items():
            _sub(c, k, v)

    mw = {n: float(COMPOUNDS[n]["Molar_Weight"]) for n in compounds}
    for s in streams:
        o = _sub(sim, "SimulationObject")
        _sub(o, "Type", "DWSIM.Thermodynamics.Streams.MaterialStream")
        _sub(o, "Name", s["id"])
        _sub(o, "SpecType", s["spec"])
        _sub(o, "PropertyPackage", PACKAGE_ID)
        _sub(o, "Calculated", "true")
        phs = _sub(o, "Phases")
        mmw = sum(mw[n] * v for n, v in s["z"].items())
        _phase(phs, "0", {"temperature": s["T"], "pressure": s["P"], "molarflow": s["F"],
                          "massflow": s["F"] * mmw / 1000.0, "enthalpy": s["h"]}, s["z"])
        _phase(phs, "1", {"molarfraction": 1.0 - s["vf"]}, s["x"])
        _phase(phs, "2", {"molarfraction": s["vf"]}, s["y"])
    for u in units:
        o = _sub(sim, "SimulationObject")
        _sub(o, "Type", f"DWSIM.UnitOperations.UnitOperations.{u['cls']}")
        _sub(o, "Name", u["id"])
        _sub(o, "PropertyPackage", PACKAGE_ID)
        _sub(o, "Calculated", "true")
        for k, v in u["params"].items():
            if isinstance(v, (list, tuple)):
                e = _sub(o, k)
                for item in v:
                    _sub(e, "double", item)
            elif v is not None:
                _sub(o, k, v)

    feeds_to = {}
    made_by = {}
    for u in units:
        for i, sid in enumerate(u["inlets"]):
            feeds_to[sid] = (u["id"], i)
        for i, sid in enumerate(u["outlets"]):
            made_by[sid] = (u["id"], i)
    for s in streams:
        g = _sub(gfx, "GraphicObject")
        _sub(g, "Name", s["id"])
        _sub(g, "Tag", s["tag"])
        _sub(g, "ObjectType", "MaterialStream")
        ic, oc = _sub(g, "InputConnectors"), _sub(g, "OutputConnectors")
        if s["id"] in made_by:
            _sub(ic, "Connector", IsAttached="True", AttachedFromObjID=made_by[s["id"]][0])
        if s["id"] in feeds_to:
            _sub(oc, "Connector", IsAttached="True", AttachedToObjID=feeds_to[s["id"]][0])
    for u in units:
        g = _sub(gfx, "GraphicObject")
        _sub(g, "Name", u["id"])
        _sub(g, "Tag", u["tag"])
        _sub(g, "ObjectType", u["cls"])
        ic, oc = _sub(g, "InputConnectors"), _sub(g, "OutputConnectors")
        for sid in u["inlets"]:
            if sid is None:
                _sub(ic, "Connector", IsAttached="False")
            else:
                _sub(ic, "Connector", IsAttached="True", AttachedFromObjID=sid)
        for sid in u["outlets"]:
            if sid is None:
                _sub(oc, "Connector", IsAttached="False")
            else:
                _sub(oc, "Connector", IsAttached="True", AttachedToObjID=sid)
    return ET.ElementTree(root)


def write(path, streams, units, zipped=None, **kw):
    """Write a .dwxml (plain) or .dwxmz (zipped) file and return its path."""
    tree = build_xml(streams, units, **kw)
    path = str(path)
    if zipped or (zipped is None and path.endswith(".dwxmz")):
        import io
        buf = io.BytesIO()
        tree.write(buf, encoding="utf-8", xml_declaration=True)
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("flowsheet.xml", buf.getvalue())
    else:
        tree.write(path, encoding="utf-8", xml_declaration=True)
    return path


# ---------------------------------------------------------------------------
# a cooler followed by a flash drum (methane/propane, Peng-Robinson)
# ---------------------------------------------------------------------------
Z = {"Methane": 0.6, "Propane": 0.4}


def cooler_flash(path, outlet_T=230.0, guesses=None):
    """FEED -> COOLER (outlet temperature spec, 50 kPa drop) -> S2 -> FLASH -> VAP / LIQ.

    ``guesses`` overrides the stored "DWSIM" results of the calculated streams.
    """
    g = {
        # converged values of the converted model at outlet_T = 230 K
        "S2": dict(T=outlet_T, P=1.95e6, F=10.0, vf=0.544359),
        "VAP": dict(T=outlet_T, P=1.95e6, F=5.44359, vf=1.0,
                    z={"Methane": 0.934189, "Propane": 0.065811}),
        "LIQ": dict(T=outlet_T, P=1.95e6, F=4.55641, vf=0.0,
                    z={"Methane": 0.200741, "Propane": 0.799259}),
    }
    for k, v in (guesses or {}).items():
        g[k].update(v)
    streams = [
        stream("MAT-feed", "FEED", 300.0, 2.0e6, 10.0, Z, 1.0),
        stream("MAT-s2", "S2", g["S2"]["T"], g["S2"]["P"], g["S2"]["F"], Z, g["S2"]["vf"]),
        stream("MAT-vap", "VAP", g["VAP"]["T"], g["VAP"]["P"], g["VAP"]["F"], g["VAP"]["z"], 1.0),
        stream("MAT-liq", "LIQ", g["LIQ"]["T"], g["LIQ"]["P"], g["LIQ"]["F"], g["LIQ"]["z"], 0.0),
    ]
    units = [
        unit("COOL-1", "COOLER", "Cooler", ["MAT-feed"], ["MAT-s2"],
             CalcMode="OutletTemperature", OutletTemperature=outlet_T, DeltaP=50000, Eficiencia=100),
        unit("FLASH-1", "FLASH", "Vessel", ["MAT-s2"], ["MAT-vap", "MAT-liq"]),
    ]
    return write(path, streams, units, kij={("Methane", "Propane"): 0.0119})


# ---------------------------------------------------------------------------
# a vapour-only network: valve, compressor, mixer, splitter, heat exchanger
# ---------------------------------------------------------------------------
def gas_network(path, heat_loss_kw=0.0):
    """Methane-rich gas, all streams above their dew point (vapour-only package).

    A -> VALVE -> B -> COMP -> C --\\
                                    MIX -> E -> SPLIT -> F (hot side of HX) -> H
    D (second feed) ---------------/              \\-> G
    COLD (feed, cold side of HX) -> HX -> COLD_OUT

    Stream values are rough guesses: the test checks balances and specs only.
    """
    zg = {"Methane": 0.95, "Propane": 0.05}
    streams = [
        stream("MAT-a", "A", 300.0, 3.0e6, 5.0, zg, 1.0),
        stream("MAT-b", "B", 290.0, 1.0e6, 5.0, zg, 1.0),
        stream("MAT-c", "C", 380.0, 2.5e6, 5.0, zg, 1.0),
        stream("MAT-d", "D", 350.0, 2.5e6, 3.0, zg, 1.0),
        stream("MAT-e", "E", 370.0, 2.5e6, 8.0, zg, 1.0),
        stream("MAT-f", "F", 370.0, 2.5e6, 4.8, zg, 1.0),
        stream("MAT-g", "G", 370.0, 2.5e6, 3.2, zg, 1.0),
        stream("MAT-h", "H", 330.0, 2.45e6, 4.8, zg, 1.0),
        stream("MAT-cold", "COLD", 280.0, 1.0e6, 6.0, zg, 1.0),
        stream("MAT-coldout", "COLD_OUT", 320.0, 0.98e6, 6.0, zg, 1.0),
    ]
    units = [
        unit("VALV-1", "VALVE", "Valve", ["MAT-a"], ["MAT-b"], CalcMode="OutletPressure",
             OutletPressure=1.0e6),
        unit("COMP-1", "COMP", "Compressor", ["MAT-b"], ["MAT-c"], CalcMode="OutletPressure",
             POut=2.5e6, AdiabaticEfficiency=75, ProcessPath="Adiabatic"),
        unit("MIX-1", "MIX", "Mixer", ["MAT-c", "MAT-d"], ["MAT-e"], PressureCalculation="Minimum"),
        unit("SPLT-1", "SPLIT", "Splitter", ["MAT-e"], ["MAT-f", "MAT-g"], OperationMode="SplitRatios",
             Ratios=[0.6, 0.4]),
        unit("HX-1", "HX", "HeatExchanger", ["MAT-f", "MAT-cold"], ["MAT-h", "MAT-coldout"],
             CalculationMode="CalcTempHotOut", ColdSideOutletTemperature=320.0, Area=10.0,
             HotSidePressureDrop=50000, ColdSidePressureDrop=20000, FlowDir="CounterCurrent",
             HeatLoss=heat_loss_kw),
    ]
    return write(path, streams, units, kij={("Methane", "Propane"): 0.0119})
