import math

import pytest

from dwsim2idaes import correlations as c

CH4_CP = [33151.9, -1220.001, 12.0907, -0.000384791, 9.896403e-08]          # ChemSep eq 16, J/kmol/K
CH4_PSAT = [39.98844, -1337.308, -3.580049, 3.20698e-05, 2.0]               # eq 101, Pa
CH4_HVAP = [14418000, 2.3055, -5.4199, 5.658, -2.1286]                      # eq 106, J/kmol
CH4_RHO = [1.894, 0.23603, 191.05, 0.21974, 0.0]                            # eq 105, kmol/m3


def test_equation_forms():
    A, B, C, D, E = CH4_CP
    T = 400.0
    assert c.eval_eq(16, A, B, C, D, E, T) == pytest.approx(A + math.exp(B / T + C + D * T + E * T * T))
    A, B, C, D, E = CH4_PSAT
    assert c.eval_eq(101, A, B, C, D, E, 150.0) == pytest.approx(
        math.exp(A + B / 150.0 + C * math.log(150.0) + D * 150.0 ** E))


def test_methane_vapour_pressure_at_normal_boiling_point():
    spec = c.psat_spec("ChemSep", "101", CH4_PSAT)
    assert c.psat_value(spec, 111.66, 190.56, 4.599e6, 0.011) == pytest.approx(101325, rel=0.01)


def test_lee_kesler_fallback_when_no_correlation():
    spec = c.psat_spec("ChemSep", "0", [0, 0, 0, 0, 0])
    assert spec["form"] == "lee_kesler"
    assert c.psat_value(spec, 111.66, 190.56, 4.599e6, 0.011) == pytest.approx(101325, rel=0.05)


def test_methane_cp_and_enthalpy_integral():
    spec = c.cp_ig_spec("ChemSep", "16", CH4_CP)
    cp298 = c.cp_ig_value(spec, 298.15)
    assert cp298 == pytest.approx(35.7, rel=0.01)          # J/mol/K, textbook value
    # enthalpy = integral of cp (trapezoid check)
    n, t0, t1 = 4000, 298.15, 800.0
    h = sum(c.cp_ig_value(spec, t0 + (t1 - t0) * (i + 0.5) / n) for i in range(n)) * (t1 - t0) / n
    assert c.enth_ig_value(spec, t1, t0) == pytest.approx(h, rel=1e-6)
    s = sum(c.cp_ig_value(spec, T) / T for T in [t0 + (t1 - t0) * (i + 0.5) / n for i in range(n)]) * (t1 - t0) / n
    assert c.entr_ig_value(spec, t1, t0) == pytest.approx(s, rel=1e-6)


def test_methane_hvap_and_liquid_density():
    h = c.hvap_spec("ChemSep", "106", CH4_HVAP)
    assert c.hvap_value(h, 111.66, 190.56) == pytest.approx(8.2e3, rel=0.03)   # J/mol
    assert c.hvap_value(h, 250.0, 190.56) >= 0.0                               # capped above Tc
    d = c.dens_liq_spec("ChemSep", "105", CH4_RHO, tc=190.56, pc=4.599e6, omega=0.011, zra=0.2876)
    assert c.dens_liq_value(d, 111.66, 190.56, 4.599e6) == pytest.approx(26.3e3, rel=0.02)  # mol/m3
    r = c.dens_liq_spec("DWSIM", "", [0] * 5, tc=190.56, pc=4.599e6, omega=0.011, zra=0.2876)
    assert r["form"] == "rackett"
    assert c.dens_liq_value(r, 111.66, 190.56, 4.599e6) == pytest.approx(26.3e3, rel=0.05)


def test_unsupported_equation_raises():
    with pytest.raises(c.CorrelationError):
        c.cp_ig_spec("ChemSep", "999", CH4_CP)
