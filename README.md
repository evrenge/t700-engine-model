# t700-engine-model

A Python replication of the real-time turboshaft engine model in:

> Ballin, M. G., *A High Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
> NASA TM-100991, July 1988.

The subject is the **General Electric T700-GE-700** — the engine on the UH-60A Black Hawk.
Ballin built a component-level thermodynamic model that had to compute in under
**10 milliseconds per frame** on 1988 hardware while a pilot flew a helicopter around it,
and the report is largely the record of how he made an implicit engine model run in real
time and what it cost him.

**[Read the engineering report →](https://evrenge.github.io/t700-engine-model/)** — what was
built, how the source data was recovered from a scanned document, which modelling decisions
depart from the report and what each costs, and the validation against every number the
report prints, with all twelve figure sets.

The source document is included as `docs/ballin-tm100991.pdf` — NTRS accession N88-26378,
not subject to US copyright under 17 U.S.C. §105. Everything here is derived from it.

## Results

Every figure below is recomputed by `validation/report.py`, which re-runs the model rather
than quoting a stored number.

| | |
|---|---|
| **Table B.1, the printed engine state** | **rms 0.22 %<!--rj:table_b1.rms-->** over 21<!--rj:table_b1.n--> numbers — worst shaft power −0.81 %<!--rj:table_b1.worst--> at descent |
| **Table 1, printed eigenvalues** | all **27**<!--rj:eigenvalues.notes.table_1_stats.n--> modes across three model variants: **rms 1.9 %<!--rj:eigenvalues.notes.table_1_stats.rms-->**, worst 4.7 %<!--rj:eigenvalues.notes.table_1_stats.worst--> (ours more damped) on the 2-DOF NG mode at hover, modes paired by dominant state, with the report's own ±2 % extraction step |
| **Appendix B, 297 printed elements** | zero structure exact; **all 206 non-zero elements** carry a numeric per-element comparison. `A` rms 6.3 %<!--rj:appendix_b.A.rms--> over 158<!--rj:appendix_b.A.n-->, `b` rms 6.4 %<!--rj:appendix_b.b.rms--> over 42<!--rj:appendix_b.b.n--> — the heat-sink rows carry the report's own TC_T41 inconsistency (open question #31) |
| **Closed loop vs Table B.1** | **worst 0.25 %<!--rj:closed_loop.worst-->** over NG, NP, Wf, Ps3 and shaft power at three trims — **with fuel flow as an output**, averaged over the governor's own limit cycle |
| Figure 6, NG vs fuel flow | mean 0.00 %NG<!--rj:steady.figure_6.mean-->, rms 0.23<!--rj:steady.figure_6.rms-->, worst −1.01<!--rj:steady.figure_6.worst--> over all 29<!--rj:steady.figure_6.n--> printed points |
| Figure 7, shaft power vs fuel flow | mean +0.18 %<!--rj:steady.figure_7.mean-->, rms 0.68 %<!--rj:steady.figure_7.rms-->, worst +2.15<!--rj:steady.figure_7.worst--> over 14<!--rj:steady.figure_7.n--> of 15<!--rj:steady.figure_7.printed_points--> |
| Figure 8, Ps3 vs NG | mean +0.06 %<!--rj:steady.figure_8.mean-->, rms 0.27 %<!--rj:steady.figure_8.rms-->, worst +0.65<!--rj:steady.figure_8.worst--> over 26<!--rj:steady.figure_8.n--> of 27<!--rj:steady.figure_8.printed_points--> |
| Figures 9 and 10, fuel transients | whole-curve rms **0.62<!--rj:whole_curve_summary.min-->–3.80 %<!--rj:whole_curve_summary.max-->** of each panel's excursion, mean 1.64 %<!--rj:whole_curve_summary.mean--> over ten panels |
| The report's own claim about Figure 6 | reproduced: we overestimate fuel by **+6.31 %<!--rj:fuel_overestimate.ours_peak_pct-->** in the 81–86 %NG band it names, against its own line's +5.99 %<!--rj:fuel_overestimate.ballin_peak_pct--> and its printed "as much as five percent" |
| Frame cost | **31 µs** at a held trim, against the report's 10 ms budget |
| Tests | `pytest` inside the `t700-dev` container runs them all; `ruff` keeps lint and formatting clean |

## What is in here

```
src/t700/        the model. NumPy and the standard library, nothing else
  engine.py      the gas path, Eqs. 1-49, in the report's own order
  realtime.py    the real-time frame, Eqs. 69-80 -- the model Ballin ran
  trim.py        equilibrium solver with continuation
  maps.py        the eleven function tables, loaded and interpolated
  linear.py      the Jacobian, extracted by the report's own method
  appendix_b.py  the published linear models, for element-by-element comparison
  constants.py   Table A.1, every value cited to a page
  corrections.py theta, delta, corrected speed and flow
  units.py       the report's US customary units; nothing converts inline
  control/       the fuel control: ECU, HMU, the closed loop, 73 constants, 8 schedules
  thermo/        gas properties behind one interface
data/maps/       11 component maps, digitized from printed figures
data/schedules/  all 8 of Appendix C's scheduling functions
data/linear/     the 297 printed Appendix B matrix elements
data/reference/  transient and sweep traces to validate against
tools/           the digitizers; each reproduces its CSV byte for byte
tools/figure/    the one library they share: raster, frame, ticks, axis, glyph,
                 provenance -- each stage of reading a printed figure, once
tests/           fast unit tests, and the architecture rules as executable checks
validation/      comparisons against the report, and the figure set
docs/notes/      5,800 lines: the report transcribed -- equations, symbols, appendix
                 inventories -- plus the open-questions ledger
```

## The rule this was built under

`CLAUDE.md` carries the full set. One matters more than the rest:

> **Every number carries a citation to a page of the report.** Never a value from general
> turbomachinery knowledge, from another engine model, or from a plausible-looking guess.
> If the report does not give it, that goes in the open-questions ledger and the work stops
> until it is decided.

A model that runs on invented numbers is worse than one that does not run, because it looks
like it works.

A second rule did most of the day-to-day work: **never tune a constant to make a comparison
pass.** Successive rounds of recalibration took shaft power at the hover and level trims to
−0.05 %, and not one data value was adjusted to fit — what was wrong each time was our own
pixels-to-numbers map, not the report's data.

Both rules, and the architecture constraints that go with them, are enforced by tests rather
than by intention where that is possible at all: `tests/test_architecture.py` walks every core
module for banned dependencies, clocks and RNGs; `tests/test_scope_tolerances.py` fails if a
tolerance printed in `SCOPE.md` differs from the one the code asserts; and
`bash tools/reproduce_all.sh` reruns every digitizer and diffs its output against the
committed CSV.

## Status

Phases 0–6 of `SCOPE.md` are complete. The engine trims and runs transients in **either of
the two configurations Ballin published** — with or without the station 4.1 heat-sink model
— and the Appendix C fuel control is closed around it, so fuel flow is an **output** of the
model rather than an input to it.

What is not claimed, and why, is set out in `SCOPE.md`. In short: Figures 11–15 need the
Gen Hel UH-60A blade-element simulation, which this report consumes and does not contain;
Table 3 was produced with functions from a different engine and is evidence of method, never
pass/fail; and nothing above 100 %NG is claimed, because Ballin claims nothing there either.

Seven questions remain open and four are partly closed, out of 66 logged in
`docs/notes/open-questions.md`, each row stating what was asked and what was decided.
Every open one is a place where the report is silent and a
choice was made, not a place where something is unexplained.

How those numbers were reached — what was wrong, how it was found, and what each
correction moved — is in git history rather than in the tree. `git log -p` over
`docs/notes/` and `SCOPE.md` is the record.

## Running it

The model needs NumPy. Analysis and digitizing also want SciPy, Matplotlib and Pillow —
and **poppler** (`pdftoppm`, `pdfimages`), which is not a Python package and so cannot be
declared as a dependency; the digitizers shell out to it to raster the source scan.

```bash
pip install -e '.[tools,dev]'
pytest
PYTHONPATH=src python validation/report.py   # every comparison, as figures and report.json
```

`.[dev]` pulls in `.[tools]`, because `testpaths` includes `validation/` and those modules
import SciPy at module level. A plain `pip install .` also works: the wheel carries
`data/` inside the package.

`Containerfile` builds the environment the published numbers were produced under, and
`docs/notes/environment.md` records its versions.

## The report site

`site/` is generated. To rebuild it after a model change:

```bash
PYTHONPATH=src python validation/report.py   # 24 figures + assets/report.json
python tools/build_site.py                   # index.html
```

`validation/report.py` recomputes every number and redraws every figure from the model, in
both light and dark, so a rebuild is a fresh measurement rather than a re-render.
`.github/workflows/pages.yml` publishes it on any push that touches `site/`. See
`site/README.md`.

## Licence and provenance

This implementation is MIT — see `LICENSE`.

The source document is a NASA technical memorandum, not subject to copyright in the
United States under 17 U.S.C. §105. Two qualifications are written out in `LICENSE`:
§105 is a US-only provision, and the report reproduces at least one figure it does not
own — Figure 2 is captioned "(From ref. 6.)", a General Electric training guide. No
number in this repository comes from that figure.

The scan itself is pinned by checksum in `docs/notes/environment.md`, along with its NTRS
accession (N88-26378) and the exact software versions every number above was produced
under. That matters more than it sounds: every citation in the codebase is `pdf p.NN`
against *this* digitization, whose page numbering comes from its front matter and would
differ in another copy.
