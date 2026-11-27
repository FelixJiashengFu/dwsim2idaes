import fixture
from dwsim2idaes.cli import main


def test_inspect_and_convert(tmp_path, capsys):
    path = fixture.cooler_flash(tmp_path / "cf.dwxmz")
    assert main(["inspect", path]) == 0
    out = capsys.readouterr().out
    assert "Methane" in out and "convertible" in out
    target = tmp_path / "model.py"
    assert main(["convert", path, "-o", str(target)]) == 0
    compile(target.read_text(), str(target), "exec")


def test_survey(tmp_path, capsys):
    fixture.cooler_flash(tmp_path / "a.dwxmz")
    fixture.gas_network(tmp_path / "b.dwxml")
    assert main(["survey", str(tmp_path)]) == 0
    assert "2 of 2 files convertible" in capsys.readouterr().out


def test_unreadable_file(tmp_path, capsys):
    bad = tmp_path / "bad.dwxml"
    bad.write_text("<x/>")
    assert main(["inspect", str(bad)]) == 2
