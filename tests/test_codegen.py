import pytest

import fixture
from dwsim2idaes import ConversionError, convert, read_flowsheet
from dwsim2idaes.codegen import analyze


def _fs(tmp_path, maker=fixture.cooler_flash, **kw):
    return read_flowsheet(maker(tmp_path / "f.dwxmz", **kw))


def test_generated_script_is_valid_python(tmp_path):
    c = convert(_fs(tmp_path))
    compile(c.code, "generated.py", "exec")
    assert "Heater(" in c.code and "Flash(" in c.code
    assert "def build_model" in c.code and "def set_specifications" in c.code
    assert c.thermo["model"] == "PR"
    assert c.thermo["kij"] == {"Methane|Propane": 0.0119}
    assert [u.idaes for u in c.units] == ["Heater", "Flash"]


def test_cooler_spec_and_pressure_drop(tmp_path):
    c = convert(_fs(tmp_path))
    assert "u.outlet.temperature.fix(230.0)" in c.code
    assert "u.deltaP.fix(-50000.0)" in c.code


def test_phase_sets_are_chosen_from_dwsim_results(tmp_path):
    c = convert(_fs(tmp_path, fixture.gas_network))
    assert "props_vap" in c.code and "m.fs.props =" not in c.code
    c = convert(_fs(tmp_path, fixture.gas_network), phases="vle")
    assert "props_vap" not in c.code


def test_heat_loss_is_kept(tmp_path):
    c = convert(_fs(tmp_path, fixture.gas_network, heat_loss_kw=5.0))
    assert "unit_heat_balance_with_loss" in c.code
    assert any("heat loss of 5 kW" in n for n in c.notes)


def _with_unit(tmp_path, cls, **params):
    """The cooler/flash fixture with one extra (possibly unsupported) unit on the vapour."""
    import xml.etree.ElementTree as ET
    path = fixture.cooler_flash(tmp_path / "f.dwxml")
    tree = ET.parse(path)
    root = tree.getroot()
    sim, gfx = root.find("SimulationObjects"), root.find("GraphicObjects")
    o = ET.SubElement(sim, "SimulationObject")
    ET.SubElement(o, "Type").text = f"DWSIM.UnitOperations.UnitOperations.{cls}"
    ET.SubElement(o, "Name").text = "X-1"
    for k, v in params.items():
        ET.SubElement(o, k).text = str(v)
    import copy
    src = next(el for el in sim if el.findtext("Name") == "MAT-vap")
    s = copy.deepcopy(src)
    s.find("Name").text = "MAT-out"
    sim.append(s)
    for name, tag, ins, outs in (("X-1", "X", ["MAT-vap"], ["MAT-out"]), ("MAT-out", "OUT", [], [])):
        g = ET.SubElement(gfx, "GraphicObject")
        ET.SubElement(g, "Name").text = name
        ET.SubElement(g, "Tag").text = tag
        ic, oc = ET.SubElement(g, "InputConnectors"), ET.SubElement(g, "OutputConnectors")
        for i in ins:
            ET.SubElement(ic, "Connector", IsAttached="True", AttachedFromObjID=i)
        for o_ in outs:
            ET.SubElement(oc, "Connector", IsAttached="True", AttachedToObjID=o_)
    tree.write(path)
    return read_flowsheet(path)


def test_unsupported_unit_is_reported_and_can_be_cut(tmp_path):
    fs = _with_unit(tmp_path, "DistillationColumn")
    info = analyze(fs)
    assert info["unsupported_units"] == ["DistillationColumn"] and not info["convertible"]
    with pytest.raises(ConversionError) as e:
        convert(fs)
    assert "DistillationColumn (X)" in str(e.value)
    c = convert(fs, cut_unsupported=True)
    assert c.cut_units and "X" in " ".join(c.notes)
    compile(c.code, "generated.py", "exec")


def test_custom_script_units_are_not_translated(tmp_path):
    fs = _with_unit(tmp_path, "Heater", OverrideCalculationRoutine="true", CalcMode="OutletTemperature",
                    OutletTemperature=230)
    with pytest.raises(ConversionError, match="custom script"):
        convert(fs)
    assert convert(fs, cut_unsupported=True).cut_units


def test_pump_with_vapour_inlet_is_refused(tmp_path):
    fs = _with_unit(tmp_path, "Pump", CalcMode="Delta_P", DeltaP=1e5)
    with pytest.raises(ConversionError, match="inlet contains vapour"):
        convert(fs)


def test_expander_power_sign(tmp_path):
    fs = _with_unit(tmp_path, "Expander", CalcMode="Delta_P", DeltaP=5e5, AdiabaticEfficiency=80,
                    DeltaQ=12.5)
    c = convert(fs)
    assert "u.deltaP.fix(-500000.0)" in c.code
    plan = [u for u in c.units if u.kind == "Expander"][0]
    assert plan.duty["dwsim"] == 12500.0 and plan.duty["sign"] == -1.0


def test_logical_blocks_are_ignored(tmp_path):
    fs = _with_unit(tmp_path, "Recycle")
    c = convert(fs)
    assert "Arc(" in c.code
