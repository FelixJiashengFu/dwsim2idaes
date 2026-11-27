"""Intermediate representation of a DWSIM flowsheet.

Everything in here is plain data in SI units (K, Pa, mol/s, kg/s, W, J/mol),
independent of both DWSIM and IDAES, and can be dumped to JSON.
"""
from __future__ import annotations

import dataclasses
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


def sanitize(name: str) -> str:
    """Turn a DWSIM display name into a valid Python/Pyomo identifier."""
    s = re.sub(r"[^0-9a-zA-Z_]", "_", name.strip())
    s = re.sub(r"_+", "_", s).strip("_")
    if not s:
        s = "x"
    if s[0].isdigit():
        s = "_" + s
    return s


@dataclass
class Correlation:
    """A DWSIM/ChemSep temperature-dependent correlation ``y = f(eq; A..E; T)``."""

    eq: str = ""
    A: float = 0.0
    B: float = 0.0
    C: float = 0.0
    D: float = 0.0
    E: float = 0.0
    tmin: Optional[float] = None
    tmax: Optional[float] = None

    @property
    def coeffs(self):
        return (self.A, self.B, self.C, self.D, self.E)


@dataclass
class Compound:
    name: str                     # DWSIM name, e.g. "N-decane"
    key: str                      # identifier used in the IDAES model
    cas: str = ""
    formula: str = ""
    db: str = ""                  # OriginalDB: ChemSep, DWSIM, User, ...
    mw: float = 0.0               # g/mol
    tc: float = 0.0               # K
    pc: float = 0.0               # Pa
    omega: float = 0.0
    zra: float = 0.0              # Rackett compressibility
    nbp: float = 0.0              # K
    is_petroleum_fraction: bool = False
    cp_ig: Correlation = field(default_factory=Correlation)
    psat: Correlation = field(default_factory=Correlation)
    hvap: Correlation = field(default_factory=Correlation)
    dens_liq: Correlation = field(default_factory=Correlation)


@dataclass
class PropertyPackage:
    id: str
    kind: str                     # e.g. "PengRobinson", "SRK", "Raoult", "NRTL"
    tag: str = ""
    kij: Dict[str, float] = field(default_factory=dict)   # "A|B" -> kij
    options: Dict[str, str] = field(default_factory=dict)


@dataclass
class PhaseState:
    fraction: float = 0.0                       # molar phase fraction
    x: Dict[str, float] = field(default_factory=dict)


@dataclass
class Stream:
    id: str
    tag: str
    key: str
    temperature: float = float("nan")           # K
    pressure: float = float("nan")              # Pa
    flow_mol: float = float("nan")              # mol/s
    flow_mass: float = float("nan")             # kg/s
    enth_mass: float = float("nan")             # J/kg (DWSIM reference state)
    vapor_fraction: float = float("nan")        # molar
    z: Dict[str, float] = field(default_factory=dict)
    vapor: PhaseState = field(default_factory=PhaseState)
    liquid: PhaseState = field(default_factory=PhaseState)
    liquid2_fraction: float = 0.0               # second liquid phase, if any
    solid_fraction: float = 0.0
    spec_type: str = "Temperature_and_Pressure"
    package_id: str = ""
    calculated: bool = True
    source: Optional[str] = None                # unit id feeding this stream
    source_port: int = -1
    sink: Optional[str] = None                  # unit id consuming this stream
    sink_port: int = -1


@dataclass
class EnergyStream:
    id: str
    tag: str
    key: str
    power: float = float("nan")                 # W


@dataclass
class Unit:
    id: str
    tag: str
    key: str
    kind: str                                   # DWSIM class: Heater, Vessel, ...
    params: Dict[str, Any] = field(default_factory=dict)
    inlets: Dict[int, str] = field(default_factory=dict)    # port -> stream id
    outlets: Dict[int, str] = field(default_factory=dict)
    energy: List[str] = field(default_factory=list)         # energy stream ids
    package_id: str = ""
    calculated: bool = True


@dataclass
class Flowsheet:
    name: str = ""
    source_file: str = ""
    dwsim_version: str = ""
    compounds: Dict[str, Compound] = field(default_factory=dict)    # by DWSIM name
    packages: Dict[str, PropertyPackage] = field(default_factory=dict)
    streams: Dict[str, Stream] = field(default_factory=dict)        # by id
    energy_streams: Dict[str, EnergyStream] = field(default_factory=dict)
    units: Dict[str, Unit] = field(default_factory=dict)
    reactions: int = 0
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    def stream_by_tag(self, tag: str) -> Stream:
        for s in self.streams.values():
            if s.tag == tag or s.key == tag:
                return s
        raise KeyError(tag)

    def unit_by_tag(self, tag: str) -> Unit:
        for u in self.units.values():
            if u.tag == tag or u.key == tag:
                return u
        raise KeyError(tag)
