"""DWSIM / ChemSep temperature-dependent property correlations.

:func:`eval_eq` mirrors ``PropertyPackage.CalcCSTDepProp`` in the DWSIM source.
It is written with Pyomo's ``exp``/``log`` so the same code evaluates plain
floats (for checks and fitting) and builds Pyomo expressions (for the model).
"""
from __future__ import annotations

import math
from typing import Sequence, Tuple

import numpy as np
from pyomo.environ import exp, log

R_GAS = 8.314462618  # J/mol/K

SUPPORTED_EQS = {
    1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16, 17, 100, 101, 102, 103, 104, 105, 106, 115, 119,
    207, 208,
}


class CorrelationError(Exception):
    pass


def is_int_eq(eq) -> bool:
    try:
        int(str(eq).strip())
        return True
    except (TypeError, ValueError):
        return False


def eval_eq(eq, A, B, C, D, E, T, Tc=None):
    """Evaluate DWSIM equation number ``eq`` at ``T`` (float or Pyomo expr)."""
    n = int(str(eq).strip())
    if n == 1:
        return A
    if n == 2:
        return A + B * T
    if n == 3:
        return A + B * T + C * T**2
    if n == 4:
        return A + B * T + C * T**2 + D * T**3
    if n in (5, 100):
        return A + B * T + C * T**2 + D * T**3 + E * T**4
    if n == 6:
        return A + B * T + C * T**2 + D * T**3 + E / T**2
    if n in (10, 207):
        return exp(A - B / (T + C))
    if n == 11:
        return math.exp(A)
    if n == 12:
        return exp(A + B * T)
    if n == 13:
        return exp(A + B * T + C * T**2)
    if n == 14:
        return exp(A + B * T + C * T**2 + D * T**3)
    if n == 15:
        return exp(A + B * T + C * T**2 + D * T**3 + E * T**4)
    if n == 16:
        return A + exp(B / T + C + D * T + E * T**2)
    if n == 17:
        return A + exp(B + C * T + D * T**2 + E * T**3)
    if n == 101:
        return exp(A + B / T + C * log(T) + D * T**E)
    if n == 102:
        return A * T**B / (1 + C / T + D / T**2)
    if n == 103:
        return A + B * exp(-C / (T**D))
    if n == 104:
        return A + B / T + C / T**2 + D / T**8 + E / T**9
    if n == 105:
        return A / (B ** (1 + (1 - T / C) ** D))
    if n == 106:
        if not Tc:
            raise CorrelationError("equation 106 needs the critical temperature")
        Tr = T / Tc
        return A * (1 - Tr) ** (B + C * Tr + D * Tr**2 + E * Tr**3)
    if n == 115:
        return exp(A + B / T + C * log(T) + D * T**2 + E / T**2)
    if n == 119:
        return exp(A / T + B + C * T + D * T**2 + E * log(T))
    if n == 208:
        return 10 ** (A - B / (T + C))
    raise CorrelationError(f"DWSIM temperature equation {n} is not supported")


def smooth_min(a, b, eps=1e-2):
    """Smooth ``min(a, b)``; works for floats and Pyomo expressions."""
    return 0.5 * (a + b - ((a - b) ** 2 + eps) ** 0.5)


# --------------------------------------------------------------------------
# Ideal-gas heat capacity, J/mol/K
# --------------------------------------------------------------------------
def cp_ig(db: str, eq: str, coeffs: Sequence[float], T: float) -> float:
    A, B, C, D, E = coeffs
    if db in ("DWSIM", ""):
        return A + B * T + C * T**2 + D * T**3 + E * T**4          # kJ/kmol/K
    if db == "CheResources":
        return (A + B * T + C * T**2 + D * T**3) * 4.1868           # cal/mol/K
    if db in ("ChemSep", "User"):
        if not is_int_eq(eq):
            raise CorrelationError(f"ideal gas Cp given as free-form equation '{eq}'")
        v = eval_eq(eq, A, B, C, D, E, T) / 1000.0                  # J/kmol/K
        return v if v != 0.0 else 3.5 * R_GAS
    raise CorrelationError(f"compound database '{db}' is not supported for ideal gas Cp")


def cp_ig_spec(db: str, eq: str, coeffs: Sequence[float]) -> dict:
    """JSON-friendly description of the ideal gas Cp correlation in J/mol/K.

    ``{"form": "poly", "c": [a0..a4]}`` for polynomial sources (integrated
    analytically) or ``{"form": "eq", "eq": n, "c": [...], "scale": s}`` for the
    exponential ChemSep forms (integrated by Gauss-Legendre quadrature).
    """
    A, B, C, D, E = coeffs
    if db in ("DWSIM", ""):
        if not any(coeffs):
            raise CorrelationError("no ideal gas Cp coefficients")
        return {"form": "poly", "c": [A, B, C, D, E]}
    if db == "CheResources":
        return {"form": "poly", "c": [A * 4.1868, B * 4.1868, C * 4.1868, D * 4.1868, 0.0]}
    if db in ("ChemSep", "User"):
        if not is_int_eq(eq):
            raise CorrelationError(f"ideal gas Cp given as free-form equation '{eq}'")
        n = int(eq)
        if not any(coeffs) or n == 0:
            return {"form": "poly", "c": [3.5 * R_GAS, 0.0, 0.0, 0.0, 0.0]}
        if n in (1, 2, 3, 4, 5, 100):
            keep = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 100: 5}[n]
            return {"form": "poly", "c": [(c / 1000.0 if i < keep else 0.0) for i, c in enumerate(coeffs)]}
        if n not in SUPPORTED_EQS:
            raise CorrelationError(f"ideal gas Cp equation {n} is not supported")
        return {"form": "eq", "eq": n, "c": [A, B, C, D, E], "scale": 1e-3}
    raise CorrelationError(f"compound database '{db}' is not supported for ideal gas Cp")


_GL_X, _GL_W = np.polynomial.legendre.leggauss(10)
_GL_X = [float(0.5 * (x + 1.0)) for x in _GL_X]     # nodes on [0, 1]
_GL_W = [float(0.5 * w) for w in _GL_W]


def cp_ig_value(spec: dict, T):
    """Ideal gas Cp, J/mol/K."""
    c = spec["c"]
    if spec["form"] == "poly":
        return c[0] + c[1] * T + c[2] * T**2 + c[3] * T**3 + c[4] * T**4
    return spec["scale"] * eval_eq(spec["eq"], c[0], c[1], c[2], c[3], c[4], T)


def enth_ig_value(spec: dict, T, Tref: float):
    """Integral of Cp dT from Tref to T, J/mol."""
    c = spec["c"]
    if spec["form"] == "poly":
        return sum(c[i] / (i + 1) * (T ** (i + 1) - Tref ** (i + 1)) for i in range(5) if c[i] != 0.0) \
            if any(c) else 0.0 * T
    dT = T - Tref
    return dT * sum(w * cp_ig_value(spec, Tref + x * dT) for x, w in zip(_GL_X, _GL_W))


def entr_ig_value(spec: dict, T, Tref: float):
    """Integral of Cp/T dT from Tref to T, J/mol/K."""
    c = spec["c"]
    if spec["form"] == "poly":
        s = c[0] * log(T / Tref)
        for i in range(1, 5):
            if c[i] != 0.0:
                s = s + c[i] / i * (T**i - Tref**i)
        return s
    # substitute u = ln T so the integrand Cp(e^u) is smooth and well resolved
    lr = log(T / Tref)
    return lr * sum(w * cp_ig_value(spec, Tref * exp(x * lr)) for x, w in zip(_GL_X, _GL_W))


# --------------------------------------------------------------------------
# Vapour pressure (Pa), heat of vaporisation (J/mol), liquid density (mol/m3)
# as (kind, parameters) descriptions that are JSON friendly.
# --------------------------------------------------------------------------
def psat_spec(db: str, eq: str, coeffs: Sequence[float]) -> dict:
    A, B, C, D, E = coeffs
    lk = {"form": "lee_kesler"}
    if db in ("DWSIM", ""):
        if A == 0 and B == 0:
            return lk
        return {"form": "eq", "eq": 101, "c": [A, B, C, D, E], "scale": 1.0}
    if db == "CheResources":
        return {"form": "eq", "eq": 10, "c": [A, B, C, 0.0, 0.0], "scale": 133.322368}
    if db in ("ChemSep", "CoolProp", "User", "ChEDL Thermo", "KDB"):
        if not is_int_eq(eq) or int(eq) == 0 or int(eq) not in SUPPORTED_EQS:
            return lk
        return {"form": "eq", "eq": int(eq), "c": [A, B, C, D, E], "scale": 1.0}
    if db == "Biodiesel" and is_int_eq(eq) and int(eq) in SUPPORTED_EQS:
        return {"form": "eq", "eq": int(eq), "c": [A, B, C, D, E], "scale": 1000.0}
    return lk


def psat_value(spec: dict, T, tc: float, pc: float, omega: float):
    if spec["form"] == "eq":
        A, B, C, D, E = spec["c"]
        return spec["scale"] * eval_eq(spec["eq"], A, B, C, D, E, T, tc)
    Tr = T / tc
    f0 = 5.92714 - 6.09648 / Tr - 1.28862 * log(Tr) + 0.169347 * Tr**6
    f1 = 15.2518 - 15.6875 / Tr - 13.4721 * log(Tr) + 0.43577 * Tr**6
    return pc * exp(f0 + omega * f1)


def hvap_spec(db: str, eq: str, coeffs: Sequence[float]) -> dict:
    A, B, C, D, E = coeffs
    if db in ("DWSIM", ""):
        return {"form": "eq106", "c": [A / 1000.0, B, C, D, 0.0]}      # J/kmol -> J/mol
    if db in ("ChemSep", "User") and is_int_eq(eq) and int(eq) == 106 and A != 0:
        return {"form": "eq106", "c": [A / 1000.0, B, C, D, E]}
    raise CorrelationError(f"no usable heat of vaporisation correlation (database '{db}', equation '{eq}')")


def hvap_value(spec: dict, T, tc: float):
    A, B, C, D, E = spec["c"]
    Te = smooth_min(T, 0.9999 * tc)
    Tr = Te / tc
    return A * (1 - Tr) ** (B + C * Tr + D * Tr**2 + E * Tr**3)


def dens_liq_spec(db: str, eq: str, coeffs: Sequence[float], tc: float, pc: float,
                  omega: float, zra: float) -> dict:
    A, B, C, D, E = coeffs
    if db in ("ChemSep", "User", "ChEDL Thermo") and is_int_eq(eq) and int(eq) in (105, 106) and A != 0:
        return {"form": "eq", "eq": int(eq), "c": [A, B, C, D, E]}     # kmol/m3
    z = zra if zra and zra > 0 else 0.29056 - 0.08775 * omega
    return {"form": "rackett", "zra": z}


def dens_liq_value(spec: dict, T, tc: float, pc: float):
    Te = smooth_min(T, 0.99 * tc)
    if spec["form"] == "eq":
        A, B, C, D, E = spec["c"]
        if spec["eq"] == 105:
            Te = smooth_min(Te, 0.999 * C)
        return 1000.0 * eval_eq(spec["eq"], A, B, C, D, E, Te, tc)
    Tr = Te / tc
    v = R_GAS * tc / pc * spec["zra"] ** (1 + (1 - Tr) ** (2.0 / 7.0))   # m3/mol
    return 1.0 / v
