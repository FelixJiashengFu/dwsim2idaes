"""High-level entry points: convert a DWSIM file, run and validate the result."""
from __future__ import annotations

import logging
import time
import types
from dataclasses import dataclass, field
from typing import List, Optional

from .codegen import Conversion, convert
from .reader import read_flowsheet


def convert_file(path: str, cut_unsupported: bool = False, phases: str = "auto") -> Conversion:
    """Read a ``.dwxmz``/``.dwxml`` file and convert it (see :func:`codegen.convert`)."""
    return convert(read_flowsheet(path), cut_unsupported=cut_unsupported, phases=phases)


def load_module(code: str, name: str = "dwsim2idaes_generated") -> types.ModuleType:
    """Execute generated model code and return it as a module object."""
    mod = types.ModuleType(name)
    mod.__file__ = f"<{name}>"
    exec(compile(code, mod.__file__, "exec"), mod.__dict__)
    return mod


@dataclass
class Validation:
    """Outcome of building, initialising and solving a converted flowsheet."""
    status: str
    optimal: bool
    degrees_of_freedom: int
    init_failed: List[str]
    seconds: float
    streams: List[dict] = field(default_factory=list)
    duties: List[dict] = field(default_factory=list)
    summary: dict = field(default_factory=dict)
    model: Optional[object] = None

    def report(self) -> str:
        from . import runtime as rt
        lines = [f"solver status: {self.status} ({self.seconds:.1f} s)"]
        if self.init_failed:
            lines.append("units that kept the DWSIM values as their starting point: "
                         + ", ".join(self.init_failed))
        lines += ["", rt.format_comparison(self.streams)]
        if self.duties:
            lines += ["", rt.format_duties(self.duties)]
        return "\n".join(lines)


def validate(conversion: Conversion, tee: bool = False, keep_model: bool = False,
             outlvl=logging.ERROR) -> Validation:
    """Build, initialise and solve the converted model; compare it with DWSIM."""
    from . import runtime as rt
    mod = load_module(conversion.code)
    t0 = time.time()
    m = mod.build_model()
    mod.set_specifications(m)
    dof = rt.dof(m)
    failed = mod.initialize(m, outlvl=outlvl)
    try:
        results = rt.solve(m, tee=tee)
        status, optimal = str(results.solver.termination_condition), rt._ok(results)
    except Exception as err:            # e.g. IPOPT terminated by a signal
        status, optimal = f"solver failed: {str(err)[:120]}", False
    rows = rt.compare(m, mod.STREAMS)
    return Validation(
        status=status,
        optimal=optimal,
        degrees_of_freedom=dof,
        init_failed=list(failed),
        seconds=time.time() - t0,
        streams=rows,
        duties=rt.unit_duties(m, mod.DUTIES),
        summary=rt.summarize(rows),
        model=m if keep_model else None,
    )
