# dwsim2idaes

**Convert DWSIM flowsheets into equation-oriented IDAES / Pyomo models.**

[中文说明](README_zh.md)

[DWSIM](https://dwsim.org) is a free, sequential-modular process simulator: it
calculates one unit operation after another and iterates around recycles.
[IDAES](https://idaes.org) builds the same kind of flowsheet as one large system
of equations in [Pyomo](https://www.pyomo.org), which a Newton-type solver
(IPOPT) solves all at once. The equation-oriented form is what you need for
optimisation, parameter estimation, sensitivity analysis or design studies with
many degrees of freedom, but writing an IDAES model by hand for a flowsheet
that already exists in DWSIM is tedious and error-prone.

`dwsim2idaes` reads a DWSIM file (`.dwxmz` or `.dwxml`) and writes a plain,
readable Python script that builds the equivalent IDAES flowsheet, fixes the
same specifications, initialises it from the converged DWSIM solution, solves
it, and compares every stream with the values DWSIM stored in the file.

We searched GitHub, the IDAES community pages and recent papers and did not
find another tool that does this conversion.

```text
$ dwsim2idaes validate "Cavett's Problem.dwxmz"
12 unit operations, 24 material streams, 15 compounds, property model PR
...
solver status: optimal (33 s)
19 calculated streams compared: max |dT| = 0.966 K (0.31%), max |dP| = 0.00%,
max |dF| = 0.01%, max |d vapour fraction| = 0.006, max |dx| = 3.1e-05
COOL-001  work_mechanical   7.29781e+06   7.30077e+06    0.04%
...
```

## Install

```bash
pip install git+https://github.com/FelixJiashengFu/dwsim2idaes
idaes get-extensions        # installs IPOPT and the IDAES helper libraries
```

Python 3.10 or newer and IDAES 2.8 or newer are needed (developed with IDAES
2.13 and Pyomo 6.11). DWSIM itself is **not** needed: the converter only reads
the file.

## Use

```bash
dwsim2idaes inspect  plant.dwxmz              # what is in the file, can it be converted?
dwsim2idaes convert  plant.dwxmz -o plant.py  # write the IDAES model
python plant.py                               # build, initialise, solve, compare with DWSIM
dwsim2idaes validate plant.dwxmz              # the last two steps in one go (--json for a report)
dwsim2idaes survey   samples/                 # which files in a folder are convertible
```

Add `--cut-unsupported` to convert the part of a flowsheet that can be
converted: unit operations without an IDAES counterpart (columns, reactors,
CAPE-OPEN units, ...) are cut out and the streams they produce become feeds
fixed at their DWSIM values. That is often enough to optimise the heat
exchanger network or compression train around a column.

From Python:

```python
from dwsim2idaes import convert_file, load_module, validate

conv = convert_file("plant.dwxmz")
print(conv.summary())                 # unit mapping and conversion notes
mod = load_module(conv.code)          # or: open("plant.py", "w").write(conv.code)
m = mod.build_model()                 # a normal Pyomo ConcreteModel with m.fs
mod.set_specifications(m)
mod.initialize(m)
# ... unfix a specification, add an objective, solve with IPOPT ...
```

The generated script contains everything as literals at the top: `THERMO`
(compound constants, correlations and interaction parameters read from the
file), `STREAMS` (the DWSIM results, used as initial guess and reference) and
`DUTIES`. The functions `build_model()`, `set_specifications()` and
`initialize()` are short and meant to be edited. See
[`examples/`](examples/) for two generated scripts.

## How it works

1. **Read** the XML inside the file: compounds with their constants and
   temperature-dependent correlations, the property package and its binary
   interaction parameters, unit operations with their settings, material and
   energy streams with their converged results, and the connectors.
2. **Map** thermodynamics to the IDAES modular property framework. The ideal
   gas heat capacity, vapour pressure, heat of vaporisation and liquid density
   use the same equations and coefficients as DWSIM (the ChemSep/DIPPR
   equation numbers DWSIM uses, e.g. 16, 101, 105, 106, with DWSIM's unit
   conventions for each database); ideal-gas
   enthalpy and entropy are integrated exactly by Gauss-Legendre quadrature.
   The reference state is DWSIM's (ideal gas at 298.15 K, no enthalpy of
   formation), so enthalpies can be compared directly.
3. **Map** each unit operation to an IDAES model and its DWSIM calculation mode
   to the matching fixed variables (see the table below).
4. **Pick phases per unit.** A unit where DWSIM found only vapour (or only
   liquid) gets a vapour-only (liquid-only) copy of the property package.
   This avoids the ill-posed phase split of a single-phase fluid in a
   vapour-liquid formulation. All copies share the state variables, so they
   connect with ordinary Arcs.
5. **Initialise without tearing.** DWSIM already converged the flowsheet, so
   the inlet of every unit is known. Each unit is initialised on its own at
   its DWSIM inlet state, in any order; recycle blocks simply disappear (the
   recycle becomes an Arc) and no tear streams are needed. The full model is
   then solved simultaneously, with IPOPT warm-started at that point.
6. **Compare** every calculated stream (T, P, flow, vapour fraction,
   composition) and every duty and power with DWSIM.

## What is supported

| DWSIM | IDAES | DWSIM calculation modes |
|---|---|---|
| Heater, Cooler | `Heater` | outlet temperature, heat added/removed, energy stream, temperature change, outlet vapour fraction (saturated outlets via the dew/bubble point); efficiency, pressure drop |
| Valve | `PressureChanger` (isenthalpic) | outlet pressure, pressure drop (Kv/Cv modes use the DWSIM outlet pressure) |
| Pump | `Pump` | pressure increase, outlet pressure, efficiency (curves use the DWSIM outlet pressure) |
| Compressor, Expander | `PressureChanger` (isentropic) | outlet pressure, pressure change, power, energy stream; adiabatic efficiency |
| Mixer | `Mixer` | outlet pressure = minimum, maximum or average of the inlets |
| Splitter | `Separator` | split ratios, mass or molar flow of an outlet |
| Vessel (flash drum) | `Flash` (+ `Mixer` for several inlets) | adiabatic at inlet conditions, fixed temperature and/or pressure |
| Heat exchanger | `HeatExchanger` (0D, LMTD) | hot/cold outlet temperature, duty, U and A, area; heat loss; pressure drops |
| Recycle | (vanishes, becomes an Arc) | |
| Adjust, Spec, energy recycle, controllers | ignored; the manipulated variable keeps its DWSIM value | |

Property packages: **Peng-Robinson** (cubic EoS with DWSIM's kij;
complementarity VLE), **Raoult's law / ideal** (ideal gas + ideal liquid with
DWSIM's vapour pressure and heat of vaporisation correlations) and
**Soave-Redlich-Kwong** (experimental, see below).

## Validation

Every one of the 175 sample flowsheets in the DWSIM repository was read
without error. With `--cut-unsupported`, 60 of them contain something to
convert; the others consist only of unsupported units, use an unsupported
property package (NRTL, UNIQUAC, UNIFAC, steam tables, ...) or several
packages at once. The results of solving those 60 conversions and comparing
them with the DWSIM results stored in the files are in
[`docs/validation.md`](docs/validation.md).

Highlights for complete flowsheets (nothing cut out):

- **Cavett's problem** (12 units, 15 compounds, PR, nested recycles): solved in
  about 30 s; flows within 0.01%, the three compressor powers within 0.04%,
  temperatures within 0.02 K except one flash liquid where DWSIM's own outlet
  is 1 K below the vessel temperature it reports.
- **Compression and expansion** (8 compressors/expanders, PR): powers within
  0.2%, outlet temperatures within 0.15 K (two dense-phase expanders cut, see
  limitations).
- **Collins helium liquefier** (22 units, ideal): duties and powers within
  0.5% except two units where pure helium condenses (cut, see limitations).
- **Humid air cooled to its dew point** (PR): dew point exact, duty 0.1%.

Two things to keep in mind when reading the comparison: DWSIM's own results
are only as converged as its recycle tolerances (some sample files close
their mass balance only to 0.3%), and where DWSIM found two liquid phases the
IDAES model, which has one liquid phase, must deviate.

## Limitations and roadmap

- **Not yet converted:** distillation and absorption columns, reactors,
  component separators, CAPE-OPEN and scripted units, three-phase vessels,
  solids, petroleum fractions. Use `--cut-unsupported` around them.
- **Property packages:** activity-coefficient models (NRTL, UNIQUAC, UNIFAC),
  Lee-Kesler-Plöcker, steam tables and CoolProp are not mapped yet.
  Peng-Robinson is mapped with the standard PR departure functions even
  when the DWSIM file asks for Lee-Kesler enthalpies.
- **SRK is experimental:** IPOPT terminated with a segmentation fault on 3 of
  the 5 SRK samples, while the same flowsheets solve with Peng-Robinson. The
  cause has not been found yet.
- **Pure fluids that change phase** inside a unit cannot be described with
  temperature and pressure as state variables (the phase split of a boiling
  pure fluid is not fixed by T and P). Such units are refused or cut. A
  pressure-enthalpy state definition would lift this.
- **Dense-phase expansion:** an expander whose single-phase inlet is above
  the critical pressure and whose outlet is two-phase does not converge in
  the vapour-liquid formulation and is refused or cut.
- **Nearly pure two-phase mixtures with the ideal (Raoult) package:** the
  smooth VLE formulation is badly conditioned when all K-values are close to
  one, and the phase split of such streams can drift.
- Adjust and Spec blocks are not translated into equality constraints; the
  variable they manipulate stays at its converged DWSIM value.
- Only DWSIM's stored results were used for validation; DWSIM itself was not
  run, and the comparison is against whatever state the file was saved in.

Contributions of further unit mappings and property packages are welcome;
`tests/fixture.py` shows how to build small synthetic DWSIM files for tests.

## Development

```bash
pip install -e ".[dev]"
idaes get-extensions
pytest                      # 27 tests; tests marked "solve" need IPOPT
```

## License

BSD 3-Clause, see [LICENSE](LICENSE). DWSIM is a trademark of its authors and
is distributed under the GPL; this package contains no DWSIM code and reads
its file format only. The compound constants used in the tests are ChemSep
values as distributed with DWSIM.
