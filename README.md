# DARTESS

**D**art-leader **A**cceleration & **R**etardation via **T**hermo-**E**lectric **S**imulation
**S**cheme — a one-dimensional, time-domain model of a horizontal dart leader propagating along a
partially decayed intracloud lightning channel.

**Current version: v1.0** (recorded in `dartess.m` and saved in every result as `drs.version`).

The channel is a lossy, nonlinear transmission line solved with the Telegrapher's Equations
(FDTD). Its resistance per unit length is not prescribed. It is derived from a temperature- and
field-dependent conductivity that evolves in time, so the speed of the reactivation front is an
output of the model, not an input. Optional features include a quasi-electrostatic solver
(`iQES = 1`), kept for validation, and an ambient electric field `Eamb(x) = E0 + s·x`.

DARTESS was adapted from two open-source lightning models: the TERS return-stroke model
([open-lightning-models/ters1](https://github.com/open-lightning-models/ters1)), which provides
the transmission-line framework, and the TELLS upward-leader model
([open-lightning-models/tells1](https://github.com/open-lightning-models/tells1)), which
provides the thermodynamic conductivity closure and the quasi-electrostatic kernel.

![Scenario I versus Scenario II](movies/scenarioI_vs_II.gif)

*Two dart leaders, overlaid. They differ only in the channel's initial conductivity beyond
x = -1000 m. Panels show charge, current, electric field, sheath radius, temperature, and
conductivity along the channel. This animation and a second one (Telegrapher's Equations versus
the quasi-electrostatic solver) are described in `documentation.md`, Sec. 9.*

## Requirements

MATLAB, developed and tested on R2016a. No toolboxes are needed.

## Files

| File | Purpose |
|---|---|
| `dartess.m` | the model: input parameters, solver, and saving of results |
| `plotWaveFront.m` | front position, speed, and field at the front, for one or more saved runs |
| `gen1Dmov.m` | animation (optionally an MP4) of the charge, current, field, sheath radius, temperature, and conductivity profiles |
| `src/getTopResDir.m` | where results are stored (edit `topresdir` there to change it) |
| `src/getWaveSpeed.m` | front-speed binning and statistics, used by `plotWaveFront.m` |
| `src/getEambE0s.m` | reads `E0` and `s` of a saved run's ambient field |
| `reproducibility.csv` | input parameters of every simulation in the main article, one row per run |
| `documentation.md` | full description of the equations, numerics, inputs and outputs |

## How to run

Run every script from the folder that contains it, so that it finds `src/`.

1. Edit the input block near the top of `dartess.m`. Every parameter has a one-line comment.
   The values shipped are the reference case of the main article: Scenario I, `Ilead` = 100 A,
   `sigma_coff` = 1e-3 S/m, `dx` = 0.5 m, 800 µs.
2. Run `dartess.m`. Batch mode is convenient for long runs:

   ```bash
   matlab -nodisplay -nosplash -nodesktop -r "try, run('dartess.m'); catch e, disp(getReport(e)); end; exit"
   ```

3. Results go to `runs/<caseName>/drs.mat`, next to a copy of the script that produced them.
   A run with an existing `caseName` overwrites it.
4. Plot with `plotWaveFront.m`, or animate with `gen1Dmov.m`. Both read the runs listed in their
   `simCases` variable (default: `referenceCase`, the shipped `caseName`).

## Output: the `drs` structure

| Field | Contents |
|---|---|
| `drs.version` | code version that produced the run |
| `drs.t` | time of every step (the step is adaptive; `diff(drs.t)` is its history) |
| `drs.tsave` | snapshot times; snapshot 1 is the initial condition |
| `drs.xI`, `drs.xV` | the two staggered spatial grids (current/field nodes and charge/potential nodes) |
| `drs.Qsave`, `drs.Vsave` | charge per unit length and **total** potential (ambient + self-consistent) snapshots, on `xV` |
| `drs.Vamb` | ambient potential on `xV` (zero unless `Eamb` is set) |
| `drs.Isave`, `drs.Rsave` | current and resistance per unit length snapshots, on `xI` |
| `drs.rQsave` | corona sheath (charge) radius snapshots, on `xV` |
| `drs.Tsave`, `drs.Ssave`, `drs.Eksave` | temperature, conductivity, and breakdown-field snapshots, on `xI` |
| `drs.T`, `drs.S` | final temperature and conductivity |
| `drs.Igt` | current at the leader tip (`x = xmax`), every step, in the code's internal sign |
| `drs.xIwavpos`, `drs.Ewavpos` | front position and electric field magnitude at the front, every step |
| `drs.params.*` | input and derived parameters, including `Eamb`, `EambE0`, `EambS` |
| `drs.flags.*` | `iQES`, `iupdateR`, `iSavVars`, `caseName` |

The electric field is not saved. It is `E = -dV/dx`, which on the decreasing-`x` grid is
`E = +diff(Vsave,1,2)/dx` (on the grid `xI(2:end-1)`). `documentation.md` describes the
equations, numerics, and every input in detail.

## Citation

da Silva, C. L., Pantuso, J. G., et al., "Dart leader propagation in the absence of driving
fields," in preparation.

## License

GNU Affero General Public License v3.0 (`LICENSE`), the same as TERS.

Author: Caitano L. da Silva, New Mexico Institute of Mining and Technology.
