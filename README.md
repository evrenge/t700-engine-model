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
| **Table B.1's full printed state** | **rms 0.20 % over 21 numbers** — worst shp 0.72 %, then Ps3 0.31 %, P41 0.30 %, P45 0.25 %, T45 0.16 %, NG 0.13 %, T41 0.08 % |
| Jacobian eigenvalues vs Table 1 | 7 of 12 modes within 4 %, 9 within 8 %; worst −22.6 % |
| Fuel-step transients, Figures 9 and 10 | **whole-curve RMS 2.0–7.6 % of each panel's excursion, mean 3.8 %** over nine panels |
| Appendix B, 297 printed elements | all loaded, zero structure exact across the four DOF variants the report prints; ~150 numerically compared element by element; `b` worst 0.15 % |
| **Closed loop vs Table B.1** | **worst 0.74 % over NG, NP, Wf, Ps3, shp at three trims — with fuel flow as an *output*** |
| Tests | 865 passing, 3 skipped, lint and formatting clean |

Phases 0–4 and 6 of `SCOPE.md` are complete for the engine: the report is ingested, the data
captured, and the engine trims and runs transients in **either of the two configurations
Ballin published** — with or without the station 4.1 heat-sink model.

**Phase 5 is complete.** `control/hmu.py` implements Figures C9–C22, `control/ecu.py`
Figures C1–C8, and `control/loop.py` closes them on the engine. So fuel flow is now an
**output** of the model rather than an input to it.

## What is actually in here

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
docs/notes/      ~6,900 lines: equations, symbols, inventories, open questions
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
comparison pass.** Successive rounds of recalibration took shaft power at the hover and
level trims to −0.05 %, and not one data value was adjusted to fit; what was wrong every
time was our own pixels-to-numbers map, not the report's data. (The descent trim still
sits at −0.72 %, and that is the honest headline number — see the table above.)

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
agree to better than 1e-12 relative, against a 1e-9 test bound. `docs/notes/heat-sink-configuration.md` carries the evidence.

All five of his linear models extract from one call — `extract(trim, wf, dof="6dof")` —
and comparing them element by element against Appendix B rather than by eigenvalue turned
up three things the report never printed: the UH-60A load inertia (`9.2503 × J_PT`,
consistent to 0.186 % across three flight conditions), `dQreq/dNP`, and the heat-sink lead
to lag ratio. What it also showed is that our remaining derivative error is not physics but
**interpolation** — see `docs/notes/derivative-ambiguity.md`.

## Things the report gets wrong, and one it does not

Reading a 1988 scan carefully turns up genuine defects. They are reproduced as printed and
recorded, never silently corrected:

- **Figure A3's `0.11` axis label** is misplaced by 17.0 px, where every other label on the
  page sits within 5 px. The bleed plateau is 0.1091, not 0.110.
- **`TC_T41`'s units** are printed `sec^(9/5)` in two places; Eq. 51 requires `sec^(1/5)`.
- **Figure A10's axis label is inverted** relative to its own values.
- **Eq. 29 does not conserve energy**, by construction and by the report's own words: the
  station 4.5 mixed enthalpy is "proportional to" the station 4.4 value, `H45 = K_H45·H44`
  with `K_H45 = 0.9623`, a single multiplicative fraction rather than a mass-weighted mix.
  Reproduced as printed.

And four that looked like report defects and were **ours**: a torque threshold that seemed
absurd until we noticed we had read rotor-hub torque as engine torque; an exhaust pressure
ratio we nearly inverted; **Eq. 74 "contradicting its own prose"**, which was a carry bug in
our own alternative branch — it held `WA31` where Eq. 74 holds `WA3_bl`, making the
recurrence an involution that oscillated instead of converging. With the bleed carried, both
readings converge first-order to the same limit and differ by at most 0.11 % on NG across a
whole fuel step at the report's own 7 ms frame, which is inside the 0.1 %-per-frame criterion
the report works to. And a 3–6 % "inconsistency" between Figures 8 and 10 that came
from comparing a decelerating engine against an equilibrium locus. A chop *must* run below
that locus — fuel is cut, T41 falls, and the choked station 4.1 nozzle then passes the same
flow at a lower P41. Over its full record Ballin's own Figure 10 runs **−1.68 % to −8.73 %**
against his own Figure 8, and his Figure 9 accel runs **−1.43 % to +7.35 %**. Both signs are
required, and he ties the two figures together himself on pdf p.39.

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
  trim are untouched. Whole-curve RMS, mean over nine panels: **9.8 % → 3.8 %**; worst panel
  15.8 % → 7.6 %; the Figure 9 T41 overshoot 1.67× Ballin's → **0.87×**.
- **What looked like the report disagreeing with itself at 775 lbm/hr was an unsettled
  transient, and we match both sides of it.** Figure 6 gives 99.69 %NG against Figure 9's
  98.88, Figure 8 gives Ps3 244.2 against 235.3, Figure 7 gives 1724.8 shp against 1648.2 —
  spreads of 0.82, 3.77 and 4.64 %. But **Ballin's Figure 9 trace is still climbing when its
  record ends** at t = 4.46 s — at +0.24 to +0.28 %NG/s measured over its last 1.5–3 s — so
  the two figures never described the same instant. (The trace is quantised at 0.030 %NG per
  pixel row, so the *instantaneous* end slope is not resolved; an earlier "+0.142 %NG/s" was
  one fit window's answer.) Our own model reproduces 51 / 39 / 62 % of each spread purely as
  settling, and matches the sweep to −0.26 / −1.92 / −0.91 % *and* the transient to +0.15 /
  +0.30 / +0.81 % at the same time — which a real contradiction would forbid.
- **The state trajectory is right; what is left is the rate along it.** Comparing Ps3 against
  NG instead of against time discards the unprinted step time, leaving the thermodynamic path
  through the state space. On **Figure 10** ours matches Ballin's to **1.07 % over 78–90 %NG
  against a 7.24 % chord baseline** — the chord being the straight line between the two
  endpoint trims, which the steady-state tests already pin, so being 6.8× closer than it is
  real shape information. Figure 10's agreement also *improves* as the frame shrinks
  (1.25 / 1.07 / 0.94 / 0.86 % at dt = 10 / 7 / 5 / 3.5 ms). `dNG/dt` at matched NG runs
  0.75× Ballin's at 76 %NG and approaches 1 by 88 %, and at 76 % the net torque is **7.5 % of
  the turbine torque it is the difference of**, so 1 % on either term moves the rate 13.3 %.
  **Three things this does not show**, all measured rather than conceded: Figure 9's panel is
  nearly vacuous (his Ps3(NG) there is a straight line to R² = 0.9994, chord error 0.30 %, so
  its 0.59 % agreement adds almost nothing to the endpoint trims); the metric is *blind* to
  the volume constants (±30 % on `K_V3`/`K_V41` is bit-identical, since the real-time
  formulation solves the pressures algebraically); and it cannot exonerate the inertia — ±20 %
  on `J_GT` moves it by about as much as the residual itself.
- **These figures' read error is measured, not estimated — in percent of full scale.** Figures
  9–10 carry one panel whose true value is printed (`WFPH` is the input and the caption states
  both levels), so digitizing it measures the error directly: **+0.97 % and +1.34 % of full
  scale**. The currency matters more than the number. The digitizer maps pixel rows to values
  by a straight line between the frame rows, so the error is a *pixel offset*; the same offset
  reads +1.83 % of value at 400 on a 250–1000 axis and **+0.21 % on the 80–100 PCNG axis**.
  This README briefly carried the 1.83 % as a universal floor, which is nine times too
  permissive on PCNG — the channel most transient agreements are quoted on. The evidence for
  the pixel model is a prediction: the two figures read the same 400 lbm/hr trim on
  different spans, and calibrating on WFPH alone predicts their PCNG disagreement as
  **+0.340 %NG against a measured +0.340**.
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
- **A false equilibrium below flight idle.** Run the Figure 10 chop far past its 4.5 s of
  record — the test integrates 120 s — and the frame settles at 75.7 %NG where the differential model trims at
  67.0 %. It is the P45 *tolerance*: 1e-6 reaches the right root, the printed 1e-3 does not,
  because at 125 lbm/hr `f9` is extrapolated 136,704 times in a 200 s run and a loose tolerance
  on a flat iteration function invites a false fixed point. The printed criterion is kept —
  below-flight-idle fuel control is a feature the report says was eliminated (pdf p.38) — and
  the behaviour is pinned by a test so it stays known.
- **Table B.1 and Figures 6–7 disagree with each other**, by up to 5.5 % on shaft power at
  low power. We track Table B.1, which is printed numbers rather than a plot.
- Figures 11–15 are **not reproducible** — they need the Gen Hel UH-60A blade-element
  simulation, which this report consumes and does not contain. That was Ballin's boundary
  too. **So Phase 5 has no published trace to overlay** — which was known before it
  started, not discovered at the end. The check that *is* available turns out to be a
  strong one, because the engine and the control were digitized from different halves of
  the report and neither knows the other exists: close the loop, hold the load at Table
  B.1's printed shaft torque, and the model settles on the printed trims **with fuel flow
  as an output** — worst 0.74 % across NG, NP, Wf, Ps3 and shp at three conditions, and
  that one is inherited (the descent shaft power is 0.72 % out open-loop too). The
  governor also absorbs a 35–70 % collective sweep while the settled state moves under
  0.2 % and `SPDG` travels from −0.40 to +1.62, which is what a governor is *for* and is
  the thing a mis-wired loop fails first. Building it caught three errors that would
  otherwise have sat undetected: `XCPC` is percent of maximum and not inches; the torque
  motor's null is 0.49496 and not the 0.44 its summing junction suggests; and `W45` must
  be lbm/sec, not the nomenclature's lb/hr, or the thermocouple lookup runs 3600× off its
  printed axis.
- **Figure C30 took a third method after two failed, and the failures are the useful part.**
  `F_HM7`'s seven curves cross, so the rank-in-y rule that identifies every other multi-curve
  figure does not apply. Reading the printed digits fails: the curve runs through each glyph
  as a bar that correlates about equally with all seven templates (30 of 97 confident), and
  stripping the bar off first destroys the glyph, because any annulus wide enough to fit the
  local curve direction also reaches the neighbouring curves. Tracing through the crossings
  fails too: a one-to-one column tracer follows the geometry well but at a crossing the
  minimum-displacement assignment is not necessarily the true one. What works is noticing
  that **the curves are exact polylines between their markers** — a straight line from
  (96.76, 3.78) to (104, 1.97) predicts curve 4's ink to ±0.02 WFPAC at four intermediate
  columns — so the figure is described by its vertices and the job is to find those. One
  systematic error fell out of it: a run centre measured *at* a marker is the glyph's
  centroid, not the line's, and digits are not vertically symmetric — curve 7 read **0.03
  WFPAC high at every marker**. Fitting each segment from glyph-free columns and
  intersecting adjacent fits removed it. Segment-on-ink coverage **0.9850** over 185
  segments is the acceptance test, and the tool fails below 0.98.

Open questions are tracked in `docs/notes/open-questions.md` — **51 logged, 34 closed, 5 partly
closed, 12 open**.

**Seven of the eleven are things the report simply does not print**, and no amount of work
closes them: the initialization rule for the opened iteration (#23), the integration algorithm
(#26), the relaxation parameter and Lipschitz constant (#27), the 10 ms against 14 ms conflict
(#33), the power turbine speed behind Figures 6–8 (#35) and behind Figures 9–10 (#36), and the
fuel-step time (#37).

The other three are work, not gaps in the source: **#45** the model running outside its own
digitized envelope on the 775 lbm/hr step; **#46** the Table B.1 against Figure 7 low-power
disagreement; and **#49** linearizing the discrete real-time map. (**#50**, Figure C30's
seven crossing curves, closed on 2026-09-13 — see below.)

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
