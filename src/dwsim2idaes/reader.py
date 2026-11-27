"""Read DWSIM flowsheet files (``.dwxmz`` zipped XML or ``.dwxml`` plain XML)."""
from __future__ import annotations

import io
import math
import os
import zipfile
import xml.etree.ElementTree as ET
from typing import Dict, Optional

from .ir import (
    Compound,
    Correlation,
    EnergyStream,
    Flowsheet,
    PhaseState,
    PropertyPackage,
    Stream,
    Unit,
    sanitize,
)

# DWSIM phase indices inside a material stream
_PH_MIX, _PH_LIQ, _PH_VAP, _PH_LIQ1, _PH_LIQ2, _PH_SOLID = "0", "1", "2", "3", "4", "7"

_PACKAGE_KINDS = {
    "PengRobinsonPropertyPackage": "PengRobinson",
    "SRKPropertyPackage": "SRK",
    "RaoultPropertyPackage": "Raoult",
}

_MAX_XML_BYTES = 400 * 1024 * 1024


class DwsimReadError(Exception):
    pass


def _f(text: Optional[str], default: float = float("nan")) -> float:
    """Parse a DWSIM number. Files are written with invariant culture, but a
    few old ones use a decimal comma."""
    if text is None:
        return default
    t = text.strip()
    if not t:
        return default
    try:
        return float(t)
    except ValueError:
        pass
    if t.count(",") == 1 and "." not in t:
        try:
            return float(t.replace(",", "."))
        except ValueError:
            pass
    low = t.lower()
    if low in ("nan", "-nan"):
        return float("nan")
    if low in ("infinity", "+infinity", "∞"):
        return float("inf")
    if low in ("-infinity", "-∞"):
        return float("-inf")
    return default


def _b(text: Optional[str]) -> bool:
    return (text or "").strip().lower() == "true"


def load_xml(path: str) -> ET.Element:
    """Return the XML root of a DWSIM file. No external entities are resolved
    (ElementTree never does) and the uncompressed size is capped."""
    if not os.path.isfile(path):
        raise DwsimReadError(f"file not found: {path}")
    with open(path, "rb") as fh:
        head = fh.read(4)
    if head[:2] == b"PK":
        with zipfile.ZipFile(path) as z:
            xmls = [i for i in z.infolist() if i.filename.lower().endswith(".xml")]
            if not xmls:
                raise DwsimReadError("archive contains no XML flowsheet")
            info = max(xmls, key=lambda i: i.file_size)
            if info.file_size > _MAX_XML_BYTES:
                raise DwsimReadError("flowsheet XML is unreasonably large")
            data = z.read(info)
        return ET.parse(io.BytesIO(data)).getroot()
    if os.path.getsize(path) > _MAX_XML_BYTES:
        raise DwsimReadError("flowsheet XML is unreasonably large")
    return ET.parse(path).getroot()


def _unique(key: str, used: set) -> str:
    base, n = key, 2
    while key.lower() in used:
        key = f"{base}_{n}"
        n += 1
    used.add(key.lower())
    return key


def _corr(c: ET.Element, eq_tag: str, prefix: str, tmin: str = "", tmax: str = "") -> Correlation:
    return Correlation(
        eq=(c.findtext(eq_tag) or "").strip(),
        A=_f(c.findtext(prefix + "A"), 0.0),
        B=_f(c.findtext(prefix + "B"), 0.0),
        C=_f(c.findtext(prefix + "C"), 0.0),
        D=_f(c.findtext(prefix + "D"), 0.0),
        E=_f(c.findtext(prefix + "E"), 0.0),
        tmin=_f(c.findtext(tmin), None) if tmin else None,
        tmax=_f(c.findtext(tmax), None) if tmax else None,
    )


def _read_compounds(root: ET.Element, fs: Flowsheet) -> None:
    used: set = set()
    for c in root.findall("Compounds/Compound"):
        name = (c.findtext("Name") or "").strip()
        if not name:
            continue
        fs.compounds[name] = Compound(
            name=name,
            key=_unique(sanitize(name), used),
            cas=(c.findtext("CAS_Number") or "").strip(),
            formula=(c.findtext("Formula") or "").strip(),
            db=(c.findtext("OriginalDB") or c.findtext("CurrentDB") or "").strip(),
            mw=_f(c.findtext("Molar_Weight"), 0.0),
            tc=_f(c.findtext("Critical_Temperature"), 0.0),
            pc=_f(c.findtext("Critical_Pressure"), 0.0),
            omega=_f(c.findtext("Acentric_Factor"), 0.0),
            zra=_f(c.findtext("Z_Rackett"), 0.0),
            nbp=_f(c.findtext("Normal_Boiling_Point"), 0.0),
            is_petroleum_fraction=(c.findtext("IsPF") or "0").strip() == "1"
            or (c.findtext("IsHYPO") or "0").strip() == "1",
            cp_ig=_corr(c, "IdealgasCpEquation", "Ideal_Gas_Heat_Capacity_Const_"),
            psat=_corr(c, "VaporPressureEquation", "Vapor_Pressure_Constant_",
                       "Vapor_Pressure_TMIN", "Vapor_Pressure_TMAX"),
            hvap=_corr(c, "VaporizationEnthalpyEquation", "HVap_", "HVap_TMIN", "HVap_TMAX"),
            dens_liq=_corr(c, "LiquidDensityEquation", "Liquid_Density_Const_",
                           "Liquid_Density_Tmin", "Liquid_Density_Tmax"),
        )


def _read_packages(root: ET.Element, fs: Flowsheet) -> None:
    for p in root.findall("PropertyPackages/PropertyPackage"):
        pid = (p.findtext("ID") or "").strip()
        cls = (p.findtext("Type") or "").split(".")[-1]
        kind = _PACKAGE_KINDS.get(cls, cls.replace("PropertyPackage", "") or cls)
        pkg = PropertyPackage(id=pid, kind=kind, tag=(p.findtext("Tag") or "").strip())
        for child in p:
            if len(child) == 0 and child.text and child.tag not in ("ID", "Type", "Tag", "ComponentDescription"):
                pkg.options[child.tag] = child.text.strip()
        for ip in p.findall("InteractionParameters/InteractionParameter"):
            a, b = ip.get("Compound1"), ip.get("Compound2")
            # cubic packages store a single "Value"; activity models store
            # several named attributes and are not read here.
            if a in fs.compounds and b in fs.compounds and ip.get("Value") is not None:
                v = _f(ip.get("Value"), 0.0)
                if not math.isnan(v):
                    pkg.kij[f"{a}|{b}"] = v
        fs.packages[pid] = pkg


def _phase(obj: ET.Element, pid: str) -> Optional[ET.Element]:
    for ph in obj.findall("Phases/Phase"):
        if (ph.findtext("ID") or "").strip() == pid:
            return ph
    return None


def _phase_props(ph: Optional[ET.Element]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    if ph is None:
        return out
    for props in ph.findall("Properties"):
        if len(props):
            for c in props:
                if c.tag != "Type":
                    out[c.tag] = _f(c.text)
    return out


def _phase_x(ph: Optional[ET.Element]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    if ph is None:
        return out
    for c in ph.findall("Compounds/Compound"):
        n = (c.findtext("Name") or "").strip()
        x = _f(c.findtext("MoleFraction"), 0.0)
        out[n] = 0.0 if math.isnan(x) else x
    return out


def _read_stream(obj: ET.Element, sid: str) -> Stream:
    mix = _phase_props(_phase(obj, _PH_MIX))
    vap_ph, liq_ph = _phase(obj, _PH_VAP), _phase(obj, _PH_LIQ)
    vap, liq = _phase_props(vap_ph), _phase_props(liq_ph)
    liq2, sol = _phase_props(_phase(obj, _PH_LIQ2)), _phase_props(_phase(obj, _PH_SOLID))

    def frac(d):
        v = d.get("molarfraction", 0.0)
        return 0.0 if (v is None or math.isnan(v)) else v

    s = Stream(id=sid, tag=sid, key=sid)
    s.temperature = mix.get("temperature", float("nan"))
    s.pressure = mix.get("pressure", float("nan"))
    s.flow_mol = mix.get("molarflow", float("nan"))
    s.flow_mass = mix.get("massflow", float("nan"))
    h = mix.get("enthalpy", float("nan"))           # kJ/kg
    s.enth_mass = h * 1000.0
    s.z = _phase_x(_phase(obj, _PH_MIX))
    s.vapor = PhaseState(fraction=frac(vap), x=_phase_x(vap_ph))
    s.liquid = PhaseState(fraction=frac(liq), x=_phase_x(liq_ph))
    s.vapor_fraction = s.vapor.fraction
    s.liquid2_fraction = frac(liq2)
    s.solid_fraction = frac(sol)
    s.spec_type = (obj.findtext("SpecType") or "Temperature_and_Pressure").strip()
    s.package_id = (obj.findtext("PropertyPackage") or "").strip()
    s.calculated = _b(obj.findtext("Calculated"))
    return s


_SKIP_PARAM_TAGS = {
    "Type", "ComponentDescription", "ComponentName", "Name", "Name2", "Annotation",
    "ErrorMessage", "LastUpdated", "DynamicProperties", "DynamicPropertiesDescriptions",
    "DynamicPropertiesUnitTypes", "AttachedUtilities", "PropertyPackage", "EditorState",
}


def _unit_params(obj: ET.Element) -> dict:
    """Scalar fields plus simple lists (e.g. splitter ratios) of a unit."""
    params: dict = {}
    for c in obj:
        if c.tag in _SKIP_PARAM_TAGS:
            continue
        if len(c) == 0:
            if c.text is not None and c.text.strip() != "":
                params[c.tag] = c.text.strip()
        else:
            kids = list(c)
            if all(len(k) == 0 for k in kids) and len(kids) <= 200:
                if all(k.attrib for k in kids) and not any((k.text or "").strip() for k in kids):
                    params[c.tag] = [dict(k.attrib) for k in kids]
                else:
                    params[c.tag] = [(k.text or "").strip() for k in kids]
    return params


def read_flowsheet(path: str) -> Flowsheet:
    """Parse a DWSIM file into a :class:`Flowsheet`."""
    root = load_xml(path)
    if root.find("SimulationObjects") is None or root.find("GraphicObjects") is None:
        raise DwsimReadError("not a DWSIM flowsheet (no SimulationObjects/GraphicObjects)")

    fs = Flowsheet(name=os.path.splitext(os.path.basename(path))[0], source_file=os.path.basename(path))
    fs.dwsim_version = (root.findtext("GeneralInfo/BuildVersion") or "").strip()
    fs.reactions = len(root.findall("Reactions/Reaction"))
    _read_compounds(root, fs)
    _read_packages(root, fs)

    # --- graphic objects: tags and connectivity --------------------------------
    tags: Dict[str, str] = {}
    gtypes: Dict[str, str] = {}
    gobjs: Dict[str, ET.Element] = {}
    for g in root.find("GraphicObjects"):
        gid = (g.findtext("Name") or "").strip()
        if not gid:
            continue
        tags[gid] = (g.findtext("Tag") or "").strip() or gid
        gtypes[gid] = (g.findtext("ObjectType") or "").strip()
        gobjs[gid] = g

    used: set = set()
    for obj in root.find("SimulationObjects"):
        oid = (obj.findtext("Name") or "").strip()
        cls = (obj.findtext("Type") or "").split(".")[-1]
        tag = tags.get(oid, oid)
        if cls == "MaterialStream":
            s = _read_stream(obj, oid)
            s.tag, s.key = tag, _unique("s_" + sanitize(tag), used)
            fs.streams[oid] = s
        elif cls == "EnergyStream":
            fs.energy_streams[oid] = EnergyStream(
                id=oid, tag=tag, key=_unique("e_" + sanitize(tag), used),
                power=_f(obj.findtext("EnergyFlow")) * 1000.0)        # kW -> W
        else:
            fs.units[oid] = Unit(
                id=oid, tag=tag, key=_unique("u_" + sanitize(tag), used), kind=cls,
                params=_unit_params(obj),
                package_id=(obj.findtext("PropertyPackage") or "").strip(),
                calculated=_b(obj.findtext("Calculated")))

    # --- connectivity -----------------------------------------------------------
    for uid, unit in fs.units.items():
        g = gobjs.get(uid)
        if g is None:
            continue
        for i, c in enumerate(g.find("InputConnectors") if g.find("InputConnectors") is not None else []):
            if (c.get("IsAttached") or "").lower() != "true":
                continue
            src = c.get("AttachedFromObjID", "")
            if src in fs.streams:
                unit.inlets[i] = src
                fs.streams[src].sink, fs.streams[src].sink_port = uid, i
            elif src in fs.energy_streams and src not in unit.energy:
                unit.energy.append(src)
        for i, c in enumerate(g.find("OutputConnectors") if g.find("OutputConnectors") is not None else []):
            if (c.get("IsAttached") or "").lower() != "true":
                continue
            dst = c.get("AttachedToObjID", "")
            if dst in fs.streams:
                unit.outlets[i] = dst
                fs.streams[dst].source, fs.streams[dst].source_port = uid, i
            elif dst in fs.energy_streams and dst not in unit.energy:
                unit.energy.append(dst)
        ec = g.find("EnergyConnector")
        for c in (ec if ec is not None else []):
            if (c.get("IsAttached") or "").lower() == "true":
                dst = c.get("AttachedToObjID", "")
                if dst in fs.energy_streams and dst not in unit.energy:
                    unit.energy.append(dst)
    return fs
