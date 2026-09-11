# t700-engine-model

A Python replication of the real-time turboshaft engine model in:

> Ballin, M. G., *A High Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
> NASA TM-100991, July 1988.

The subject is the **General Electric T700-GE-700** — the engine on the UH-60A Black
Hawk. Ballin built a component-level thermodynamic model that had to compute in under
**10 milliseconds per frame** on 1988 hardware while a pilot flew a helicopter around it,
and the report is largely the record of how he made an implicit engine model run in real
time and what it cost him.

The source document is included as `docs/ballin-tm100991.pdf` (NASA technical memorandum,
public domain). Everything here is derived from it.

## Where it stands

| | |
|---|---|
| Gas generator speed vs Table B.1 | **±0.13 %** across three trim conditions |
| Shaft power vs Table B.1 | **±0.05 %** at hover and level, −0.73 % at descent |
| Jacobian eigenvalues vs Table 1 | 7 of 12 modes within 4 %, 9 within 8 %; worst −22.6 % |
| Fuel-step transient, Figure 9 | endpoints within 0.9–3.9 %; T41 overshoot 1.67× Ballin's, down from 3.41× |
| Tests | 463, with lint and formatting clean |

Phases 0–4 of `SCOPE.md` are complete: the report is ingested, the data captured, and the
engine trims and runs transients in **either of the two configurations Ballin published**
— with or without the station 4.1 heat-sink model. The fuel control system (Phase 5) is
not built.

## What is actually in here

```
src/t700/        the model. NumPy and the standard library, nothing else
  engine.py      the gas path, Eqs. 1-49, in the report's own order
  realtime.py    the real-time frame, Eqs. 69-80 -- the model Ballin ran
  trim.py        equilibrium solver with continuation
  maps.py        the eleven function tables, loaded and interpolated
  thermo/        gas properties behind one interface
data/maps/       11 component maps, digitized from printed figures
data/reference/  transient and sweep traces to validate against
tools/           the digitizers; each reproduces its CSV byte for byte
docs/notes/      ~6,000 lines: equations, symbols, inventories, open questions
validation/      comparisons against the report, and the figure set
```

## The rules this was built under

They are in `CLAUDE.md`, and one matters more than the rest:

> **Every number carries a citation to a page of the report.** Never a value from general
> turbomachinery knowledge, from another engine model, or from a plausible-looking guess.
> If the report does not give it, that goes in the open-questions ledger and the work
> stops until it is decided.

A model that runs on invented numbers is worse than one that does not run, because it
looks like it works.

A second rule did most of the day-to-day work: **never tune a constant to make a
comparison pass.** Shaft power went from −0.75 % to −0.05 % over five rounds of
recalibration and not one data value was adjusted to fit. What was wrong every time was
our own pixels-to-numbers map, not the report's data.

## Two models, not one

Ballin ran the engine with and without a station 4.1 heat-sink model and published results
for both, which makes it a switch rather than a feature. He states which is which for the
linear results and leaves it to inference for the transients:

| Result | Heat sink | Basis |
|---|---|---|
| Table 1 eigenvalues, Appendix B B1–B6 | off | printed, pdf pp.27, 67 |
| Appendix B B7–B12 | on | printed, pdf p.67 |
| Figures 9, 10 | on | inferred — pdf pp.20, 38, 47 |

Eq. 50 has unit DC gain, so the switch cannot move a trim. Every steady-state result here
is valid in both configurations at once, and that is measured rather than argued: the two
agree to ~1e-15 relative. `docs/notes/heat-sink-configuration.md` carries the evidence.

## Things the report gets wrong, and one it does not

Reading a 1988 scan carefully turns up genuine defects. They are reproduced as printed and
recorded, never silently corrected:

- **Figure A3's `0.11` axis label** is misplaced by 19.5 px — a 6σ outlier against every
  other label on the page. The bleed plateau is 0.1090, not 0.110.
- **`TC_T41`'s units** are printed `sec^(9/5)` in two places; Eq. 51 requires `sec^(1/5)`.
- **Figure A10's axis label is inverted** relative to its own values.
- **Eq. 74 contradicts its own prose.** Settled by experiment rather than argument:
  refining the time step, the prose reading converges first-order and the printed reading
  converges to nothing.

And two that looked like report defects and were **ours**: a torque threshold that seemed
absurd until we noticed we had read rotor-hub torque as engine torque, and an exhaust
pressure ratio we nearly inverted.

## Known gaps

- **The transient is still too sharp**, though much less so. T41 now overshoots by 1.67×
  Ballin's figure, down from 3.41×. Two thirds of that gap was the station 4.1 heat-sink
  model (Eqs. 48–53) — and part of it was our own error: Figures 9 and 10 were generated
  with the heat sink *on*, so comparing our 5-DOF model against them was a category
  mistake rather than an incomplete model. The remaining 1.67× is unexplained.
- **Table B.1 and Figures 6–7 disagree with each other**, by up to 5.5 % on shaft power at
  low power. We track Table B.1, which is printed numbers rather than a plot.
- Figures 11–15 are **not reproducible** — they need the Gen Hel UH-60A blade-element
  simulation, which this report consumes and does not contain. That was Ballin's boundary
  too.

Open questions are tracked in `docs/notes/open-questions.md` — 43 logged, 19 closed.

## Running it

The model needs NumPy. Analysis and digitizing also want SciPy, Matplotlib and Pillow.

```bash
pip install -e '.[tools,dev]'
pytest
python validation/plot_validation.py     # writes the figure set
```

## Licence and provenance

The source document is a NASA technical memorandum and is in the public domain. This
implementation is offered in the same spirit — see `LICENSE`.
