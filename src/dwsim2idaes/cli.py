"""Command line interface: ``dwsim2idaes inspect|convert|validate|survey``."""
from __future__ import annotations

import argparse
import glob
import json
import logging
import math
import os
import sys
from typing import List, Optional

from . import __version__
from .codegen import ConversionError, SUPPORTED_UNITS, analyze
from .reader import DwsimReadError, read_flowsheet
from .thermo import SUPPORTED_PACKAGES


def _quiet():
    for name in ("idaes", "pyomo"):
        logging.getLogger(name).setLevel(logging.ERROR)


def _conversion_args(p):
    p.add_argument("--cut-unsupported", action="store_true",
                   help="cut out unit operations that cannot be converted; the streams they produce "
                        "become feeds fixed at the DWSIM values")
    p.add_argument("--phases", choices=("auto", "vle"), default="auto",
                   help="'auto' (default): give each unit the phase set DWSIM found around it; "
                        "'vle': vapour-liquid equilibrium everywhere")


def _cmd_inspect(a) -> int:
    fs = read_flowsheet(a.file)
    info = analyze(fs)
    print(f"{fs.source_file}  (DWSIM {fs.dwsim_version or 'unknown version'})")
    print(f"  compounds ({len(fs.compounds)}): " + ", ".join(fs.compounds))
    for p in fs.packages.values():
        ok = "supported" if p.kind in SUPPORTED_PACKAGES else "not supported"
        print(f"  property package: {p.kind} ({p.tag}) - {ok}")
    print(f"  material streams: {len(fs.streams)}, energy streams: {len(fs.energy_streams)}")
    print("  unit operations:")
    for kind, n in sorted(info["unit_kinds"].items()):
        state = ("converted" if kind in SUPPORTED_UNITS else
                 "unsupported" if kind in info["unsupported_units"] else "passed through / ignored")
        print(f"    {n:3d} x {kind:<28} {state}")
    try:
        from .codegen import convert
        c = convert(fs, cut_unsupported=a.cut_unsupported, phases=a.phases)
        print("\nconvertible" + (" (with cuts)" if c.cut_units else "") + ":")
        print(c.summary())
        return 0
    except ConversionError as e:
        print("\nnot convertible:")
        for p in e.problems:
            print("  - " + p)
        if not a.cut_unsupported:
            print("  (try --cut-unsupported to convert the rest of the flowsheet)")
        return 2


def _cmd_convert(a) -> int:
    from .api import convert_file
    c = convert_file(a.file, cut_unsupported=a.cut_unsupported, phases=a.phases)
    out = a.output or os.path.splitext(os.path.basename(a.file))[0].replace(" ", "_") + "_idaes.py"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(c.code)
    print(c.summary())
    print(f"\nwrote {out}  (run it with: python {out})")
    return 0


def _clean(obj):
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean(v) for v in obj]
    return obj


def _cmd_validate(a) -> int:
    _quiet()
    from .api import convert_file, validate
    c = convert_file(a.file, cut_unsupported=a.cut_unsupported, phases=a.phases)
    if not a.json:
        print(c.summary())
        print()
    v = validate(c, tee=a.tee)
    if a.json:
        print(json.dumps(_clean({
            "file": a.file, "status": v.status, "optimal": v.optimal, "seconds": v.seconds,
            "degrees_of_freedom": v.degrees_of_freedom, "init_failed": v.init_failed,
            "notes": c.notes, "cut_units": c.cut_units, "summary": v.summary,
            "streams": v.streams, "duties": v.duties}), indent=1))
    else:
        print(v.report())
    return 0 if v.optimal else 1


def _files(patterns: List[str]) -> List[str]:
    out: List[str] = []
    for pat in patterns:
        if os.path.isdir(pat):
            for ext in ("*.dwxmz", "*.dwxml"):
                out += glob.glob(os.path.join(pat, "**", ext), recursive=True)
        else:
            out += glob.glob(pat, recursive=True) or ([pat] if os.path.exists(pat) else [])
    return sorted(dict.fromkeys(out))


def _cmd_survey(a) -> int:
    from .codegen import convert
    files = _files(a.paths)
    if not files:
        print("no DWSIM files found", file=sys.stderr)
        return 2
    ok = 0
    reasons: dict = {}
    for f in files:
        name = os.path.basename(f)
        try:
            fs = read_flowsheet(f)
            c = convert(fs, cut_unsupported=a.cut_unsupported, phases=a.phases)
            ok += 1
            cut = f", {len(c.cut_units)} cut" if c.cut_units else ""
            print(f"OK    {name}  ({len(c.units)} units{cut}, {c.thermo['model']}, "
                  f"{len(c.thermo['components'])} compounds)")
        except ConversionError as e:
            first = e.problems[0]
            print(f"NO    {name}  - {first[:150]}")
            for p in e.problems:
                key = p.split(":")[0].split("(")[0].strip()[:70]
                reasons[key] = reasons.get(key, 0) + 1
        except DwsimReadError as e:
            print(f"ERROR {name}  - {e}")
    print(f"\n{ok} of {len(files)} files convertible"
          + (" (with --cut-unsupported)" if a.cut_unsupported else ""))
    for k, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"  {n:4d}  {k}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        prog="dwsim2idaes",
        description="Convert DWSIM flowsheets (.dwxmz/.dwxml) into equation-oriented IDAES/Pyomo models.")
    ap.add_argument("--version", action="version", version=f"dwsim2idaes {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("inspect", help="list what a file contains and whether it can be converted")
    p.add_argument("file")
    _conversion_args(p)
    p.set_defaults(func=_cmd_inspect)

    p = sub.add_parser("convert", help="write a Python script with the IDAES model")
    p.add_argument("file")
    p.add_argument("-o", "--output", help="output .py file (default: <name>_idaes.py)")
    _conversion_args(p)
    p.set_defaults(func=_cmd_convert)

    p = sub.add_parser("validate", help="convert, solve with IPOPT and compare with the DWSIM results")
    p.add_argument("file")
    p.add_argument("--tee", action="store_true", help="show the IPOPT log")
    p.add_argument("--json", action="store_true", help="print a machine-readable report")
    _conversion_args(p)
    p.set_defaults(func=_cmd_validate)

    p = sub.add_parser("survey", help="check which files in a folder can be converted")
    p.add_argument("paths", nargs="+", help="files, folders or glob patterns")
    _conversion_args(p)
    p.set_defaults(func=_cmd_survey)

    a = ap.parse_args(argv)
    try:
        return a.func(a)
    except ConversionError as e:
        print("cannot convert this flowsheet:", file=sys.stderr)
        for prob in e.problems:
            print("  - " + prob, file=sys.stderr)
        if not getattr(a, "cut_unsupported", False):
            print("  (try --cut-unsupported to convert the rest of the flowsheet)", file=sys.stderr)
        return 2
    except DwsimReadError as e:
        print(f"cannot read {a.file if hasattr(a, 'file') else 'file'}: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
