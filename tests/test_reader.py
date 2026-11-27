import pytest

import fixture
from dwsim2idaes.reader import DwsimReadError, _f, read_flowsheet


@pytest.mark.parametrize("ext", ["dwxmz", "dwxml"])
def test_reads_structure(tmp_path, ext):
    fs = read_flowsheet(fixture.cooler_flash(tmp_path / f"cf.{ext}"))
    assert list(fs.compounds) == ["Methane", "Propane"]
    assert fs.dwsim_version == "8.8.3.0"
    assert {u.tag: u.kind for u in fs.units.values()} == {"COOLER": "Cooler", "FLASH": "Vessel"}
    flash = fs.unit_by_tag("FLASH")
    assert fs.streams[flash.inlets[0]].tag == "S2"
    assert [fs.streams[flash.outlets[i]].tag for i in (0, 1)] == ["VAP", "LIQ"]
    s2 = fs.stream_by_tag("S2")
    assert s2.source == fs.unit_by_tag("COOLER").id and s2.sink == flash.id


def test_stream_values_and_units(tmp_path):
    fs = read_flowsheet(fixture.cooler_flash(tmp_path / "cf.dwxmz"))
    feed = fs.stream_by_tag("FEED")
    assert feed.temperature == 300.0 and feed.pressure == 2.0e6 and feed.flow_mol == 10.0
    assert feed.z == pytest.approx({"Methane": 0.6, "Propane": 0.4})
    assert feed.vapor_fraction == 1.0
    assert feed.flow_mass == pytest.approx(10.0 * (0.6 * 16.04246 + 0.4 * 44.09562) / 1000)


def test_compound_and_package_data(tmp_path):
    fs = read_flowsheet(fixture.cooler_flash(tmp_path / "cf.dwxmz"))
    ch4 = fs.compounds["Methane"]
    assert ch4.tc == 190.56 and ch4.pc == 4599000 and ch4.db == "ChemSep"
    assert ch4.cp_ig.eq == "16" and ch4.cp_ig.A == 33151.9
    assert ch4.psat.tmin == 83.65
    pkg = next(iter(fs.packages.values()))
    assert pkg.kind == "PengRobinson"
    assert pkg.kij == {"Methane|Propane": 0.0119}


def test_lists_in_unit_parameters(tmp_path):
    fs = read_flowsheet(fixture.gas_network(tmp_path / "g.dwxmz"))
    assert fs.unit_by_tag("SPLIT").params["Ratios"] == ["0.6", "0.4"]


def test_number_parsing():
    assert _f("1.5") == 1.5
    assert _f("1,5") == 1.5           # a few old files use a decimal comma
    assert _f("Infinity") == float("inf")
    assert _f("") != _f("")           # NaN default


def test_errors(tmp_path):
    with pytest.raises(DwsimReadError):
        read_flowsheet(tmp_path / "missing.dwxmz")
    bad = tmp_path / "bad.dwxml"
    bad.write_text("<root><Something/></root>")
    with pytest.raises(DwsimReadError):
        read_flowsheet(str(bad))
