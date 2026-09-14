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
| **Table B.1, the printed engine state** | **rms 0.24 %** over 21 numbers — worst shaft power −0.90 % at descent |
| **Table 1, printed eigenvalues** | all **27** modes across three model variants: rms 5.5 %, worst +14.9 % on the 5-DOF P45 mode at hover |
| **Appendix B, 297 printed elements** | zero structure exact; **all 206 non-zero elements** carry a numeric per-element comparison. `A` rms 7.8 % over 158, `b` rms 5.2 % over 42 |
| **Closed loop vs Table B.1** | **worst 0.34 %** over NG, NP, Wf, Ps3 and shaft power at three trims — **with fuel flow as an output**, averaged over the governor's own limit cycle |
| Figure 6, NG vs fuel flow | mean +0.00 %NG, rms 0.23, worst −1.01 over all 29 printed points |
| Figure 7, shaft power vs fuel flow | mean −0.14 %, rms 0.61 %, worst −1.74 over all 15 |
| Figure 8, Ps3 vs NG | mean −0.32 %, rms 0.42 %, worst −0.94 over 26 of 27 |
| Figures 9 and 10, fuel transients | whole-curve rms **0.67–3.82 %** of each panel's excursion, mean 1.68 % over ten panels |
| The report's own claim about Figure 6 | reproduced: we overestimate fuel by **+6.17 %** in the 81–86 %NG band it names, against its own line's +5.82 % and its printed "as much as five percent" |
| Frame cost | **31 µs** at a held trim, against the report's 10 ms budget |
| Tests | 1104 passing, 3 skipped, lint and formatting clean; the suite runs in 34 s |

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
tests/           fast unit tests, and the architecture rules as executable checks
validation/      comparisons against the report, and the figure set
docs/notes/      ~7,500 lines: equations, symbols, inventories, open questions
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
than by intention where that is possible at all: `tests/test_imports.py` walks every core
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

Six questions remain open and four are partly closed, out of 64 logged in
`docs/notes/open-questions.md`. Every open one is a place where the report is silent and a
choice was made, not a place where something is unexplained.

`docs/notes/engineering-log.md` is the narrative record of how the numbers above were
reached — what was wrong, how it was found, and what each correction moved.

## Running it

The model needs NumPy. Analysis and digitizing also want SciPy, Matplotlib and Pillow —
and **poppler** (`pdftoppm`, `pdfimages`), which is not a Python package and so cannot be
declared as a dependency; the digitizers shell out to it to raster the source scan.

```bash
pip install -e '.[tools,dev]'
pytest
python validation/plot_validation.py     # writes the figure set
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
