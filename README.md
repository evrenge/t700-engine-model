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
| **Table B.1's full printed state** | **rms 0.21 % over 21 numbers** — Ps3 0.31 %, P41 0.30 %, P45 0.25 %, T45 0.16 %, NG 0.13 %, T41 0.08 %, shp 0.77 % |
| Jacobian eigenvalues vs Table 1 | 7 of 12 modes within 4 %, 9 within 8 %; worst −22.6 % |
| Fuel-step transients, Figures 9 and 10 | **whole-curve RMS 2.1-7.5 % of each panel's excursion, mean 3.8 %** |
| Appendix B, 297 printed elements | zero structure exact; P3/P41 block <1.5 %; `b` 0.1 % |
| Tests | 793 passing, with lint and formatting clean |

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

All five of his linear models extract from one call — `extract(trim, wf, dof="6dof")` —
and comparing them element by element against Appendix B rather than by eigenvalue turned
up three things the report never printed: the UH-60A load inertia (`9.2503 × J_PT`,
consistent to 0.186 % across three flight conditions), `dQreq/dNP`, and the heat-sink lead
to lag ratio. What it also showed is that our remaining derivative error is not physics but
**interpolation** — see `docs/notes/derivative-ambiguity.md`.

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

And three that looked like report defects and were **ours**: a torque threshold that seemed
absurd until we noticed we had read rotor-hub torque as engine torque; an exhaust pressure
ratio we nearly inverted; and a 3–6 % "inconsistency" between Figures 8 and 10 that came
from comparing a decelerating engine against an equilibrium locus. A chop *must* run below
that locus — fuel is cut, T41 falls, and the choked station 4.1 nozzle then passes the same
flow at a lower P41. Ballin's own Figure 10 sits 2.3–7.0 % below his own Figure 8, and his
Figure 9 accel sits up to +3.9 % above it. Both signs are required, and he ties the two
figures together himself on pdf p.39.

## Known gaps

- **The transient error was a realization artifact, and it is fixed.** Eq. 50 writes the
  station 4.1 heat sink as a lead-lag, but it is a *collapse* of the two printed
  heat-transfer equations (48) and (49) — and that collapse is only valid when the
  coefficients are constant, which Eqs. 51 and 53 make sure they are not. We carried Eq. 50's
  lead-lag memory `x`, and `x = (τb/τa)·T_m`: the gain is baked into the stored state, so
  every time the coefficients moved, `x` referred to the old gain. On a fuel step the gain
  swings 0.62 → 0.67 within a few frames, and the stale memory amplified the excursion.
  Integrating the metal temperature instead — Eqs. 48–49 as printed — has no such artifact,
  and is identical whenever the coefficients are constant, so Appendix B, Table 1 and every
  trim are untouched. Whole-curve RMS, mean over nine panels: **9.8 % → 3.6 %**; worst panel
  15.8 % → 7.0 %; the Figure 9 T41 overshoot 1.67× Ballin's → **0.87×**.
- **What looked like the report disagreeing with itself at 775 lbm/hr was an unsettled
  transient, and we match both sides of it.** Figure 6 gives 99.69 %NG against Figure 9's
  98.88, Figure 8 gives Ps3 244.2 against 235.3, Figure 7 gives 1724.8 shp against 1648.2 —
  spreads of 0.82, 3.77 and 4.64 %. But **Ballin's Figure 9 trace is still climbing at
  +0.142 %NG/s when its record ends** at t = 4.46 s, so the two figures never described the
  same instant. Our own model reproduces 51 / 39 / 62 % of each spread purely as settling,
  and matches the sweep to −0.26 / −1.92 / −0.91 % *and* the transient to +0.15 / +0.30 /
  +0.81 % at the same time — which a real contradiction would forbid.
- **The transient error is localised to the torque balance, and the torque balance is at the
  noise floor.** Comparing Ps3 against NG instead of against time discards the unprinted step
  time and the rate, leaving the thermodynamic path through the state space: **ours matches
  Ballin's to 0.59 % on Figure 9 and to 1.07 % over 78–90 %NG on Figure 10**, and the
  off-equilibrium excursion tracks in sign and shape. So nothing is left in the compressor
  map, the pressure solution, the mass balance, or the heat sink. `dNG/dt` at matched NG then
  agrees to 3.3 % over 82–88 %NG and falls to 0.75× at 76 % — where the net torque is
  **7.5 % of the turbine torque it is the difference of**, so 1 % on either term moves the
  rate by 13.3 %. The 25 % rate deficit is ~1.9 % on a torque, under the read error below.
- **These figures' read error is measured, not estimated.** Figures 9–10 carry one panel whose
  true value is printed — `WFPH` is the input and the caption states both levels — so
  digitizing it measures the error directly: **+1.83 % and +1.67 % at the 400 lbm/hr level**,
  0.26 and 0.38 % at the far levels. Deviations under ~1.8 % against these figures are under
  the reference's own error and are not chased.
- **Seven hypotheses were tested and eliminated on the way there**, and the list is worth
  keeping because each cost real work: the heat-sink time constants (a lead-lag with unit DC
  gain cannot create an overshoot, and Eqs. 50–53 were all re-read off the raster); the volume
  dynamics (integrating Eqs. 42–47 agrees to 0.11 % at settle and makes the spike *worse*); the
  frame size; the input shape; the function-table clamping (≤7.1 °R of 112); host-frame sampling
  of his plot (≤9.5 °R); and both `f1` interpolation schemes — a shape-preserving cubic across
  speed lines, and the beta lines that Figure A1 actually prints.
- **Two fixes, both replacing an invention with something printed.** The pressure iterations ran
  to `tol=1e-10` where the report states 0.1 percent in ten and eight passes (pdf p.37); and the
  heat sink is now Eqs. 48–49 rather than the collapsed Eq. 50. Neither was aimed at a figure.
- **A false equilibrium below flight idle.** Run the Figure 10 chop past ~25 s — twenty times
  the 4.5 s of record — and the frame settles at 75.7 %NG where the differential model trims at
  67.0 %. It is the P45 *tolerance*: 1e-6 reaches the right root, the printed 1e-3 does not,
  because at 125 lbm/hr `f9` is extrapolated 136,704 times in a 200 s run and a loose tolerance
  on a flat iteration function invites a false fixed point. The printed criterion is kept —
  below-flight-idle fuel control is a feature the report says was eliminated (pdf p.38) — and
  the behaviour is pinned by a test so it stays known.
- **Table B.1 and Figures 6–7 disagree with each other**, by up to 5.5 % on shaft power at
  low power. We track Table B.1, which is printed numbers rather than a plot.
- Figures 11–15 are **not reproducible** — they need the Gen Hel UH-60A blade-element
  simulation, which this report consumes and does not contain. That was Ballin's boundary
  too.

Open questions are tracked in `docs/notes/open-questions.md` — **49 logged, 33 closed, 5 partly closed, 11 open**. **Nine of the eleven are things the report simply does not print**: the initialization rule for the opened iteration, what the 0.1 % time-step criterion is measured on, the iteration counts, the integration algorithm, the relaxation parameter, the 10 ms against 14 ms conflict, and the power turbine speed and step time behind Figures 6-10. Those cannot be closed by working harder. Of the remaining three, one waits on Phase 5, one records a contradiction between two of the report's own datasets, and three are work we have not done: linearizing the discrete real-time map, Figure C30's seven crossing curves, and cleaning the off-curve samples out of the Figure 9/10 traces so transient *shape* can be compared at all.

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
