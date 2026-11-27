"""Build, initialise and solve converted models (needs IPOPT)."""
import pyomo.environ as pyo
import pytest

import fixture
from dwsim2idaes import convert_file, validate

pytestmark = pytest.mark.solve


@pytest.fixture(autouse=True, scope="module")
def _need_ipopt():
    from idaes.core.solvers import get_solver
    if not get_solver().available(exception_flag=False):
        pytest.skip("ipopt is not available")


def test_cooler_flash_reproduces_reference(tmp_path):
    v = validate(convert_file(fixture.cooler_flash(tmp_path / "cf.dwxmz")), keep_model=True)
    assert v.optimal and v.degrees_of_freedom == 0 and not v.init_failed
    # the fixture stores the converged values of this model: deviations must vanish
    assert v.summary["max_rel_dF"] < 1e-5
    assert v.summary["max_abs_dx"] < 1e-5
    assert v.summary["max_abs_dvf"] < 1e-5
    m = v.model
    vap, liq = m.fs.u_FLASH.vap_outlet, m.fs.u_FLASH.liq_outlet
    total = pyo.value(vap.flow_mol[0] + liq.flow_mol[0])
    assert total == pytest.approx(10.0, rel=1e-8)
    assert pyo.value(m.fs.u_COOLER.outlet.temperature[0]) == pytest.approx(230.0)
    assert pyo.value(m.fs.u_COOLER.heat_duty[0]) < 0


def test_gas_network_balances(tmp_path):
    v = validate(convert_file(fixture.gas_network(tmp_path / "g.dwxmz", heat_loss_kw=5.0)), keep_model=True)
    assert v.optimal and v.degrees_of_freedom == 0
    m = v.model
    split = m.fs.u_SPLIT
    assert pyo.value(split.outlet_1.flow_mol[0]) == pytest.approx(0.6 * 8.0, rel=1e-6)
    assert pyo.value(m.fs.u_MIX.outlet.pressure[0]) == pytest.approx(2.5e6)       # minimum of inlets
    assert pyo.value(m.fs.u_COMP.outlet.pressure[0]) == pytest.approx(2.5e6)
    hx = m.fs.u_HX
    q_hot, q_cold = pyo.value(hx.hot_side.heat[0]), pyo.value(hx.cold_side.heat[0])
    assert q_hot + q_cold == pytest.approx(-5000.0, abs=1e-3)                      # 5 kW lost
    assert pyo.value(hx.cold_side_outlet.temperature[0]) == pytest.approx(320.0)
    # valve is isenthalpic: Joule-Thomson cooling of the gas
    assert pyo.value(m.fs.u_VALVE.outlet.temperature[0]) < 300.0
