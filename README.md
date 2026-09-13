# t700-engine-model

A Python replication of the real-time turboshaft engine model in:

> Ballin, M. G., *A High Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
> NASA TM-100991, July 1988.

The subject is the **General Electric T700-GE-700** — the engine on the UH-60A Black
Hawk. Ballin built a component-level thermodynamic model that had to compute in under
**10 milliseconds per frame** on 1988 hardware while a pilot flew a helicopter around it,
and the report is largely the record of how he made an implicit engine model run in real
time and what it cost him.

The source document is included as `docs/ballin-tm100991.pdf` — NASA TM-100991, NTRS
accession N88-26378, not subject to US copyright under 17 U.S.C. §105 (see `LICENSE` for
the two qualifications on that). Everything here is derived from it.

## Where it stands

| | |
|---|---|
| **Table B.1's full printed state** | **rms 0.24 % over 21 numbers** — worst shp 0.91 % at descent |
| Jacobian eigenvalues vs Table 1 | worst **+7.6 %** on the 2-DOF NG mode, down from −22.6 % |
| Fuel-step transients, Figures 9 and 10 | **whole-curve RMS 1.3–7.7 % of each panel's excursion, mean 3.5 %** over nine panels — Figure 9 1.3–5.2 %, Figure 10 1.9–7.7 % |
| Appendix B, 297 printed elements | all loaded, zero structure exact across the four DOF variants the report prints; ~150 numerically compared element by element; `b` worst 0.15 % |
| **Closed loop vs Table B.1** | **worst 0.62 % over NG, NP, Wf, Ps3, shp at three trims — with fuel flow as an *output*** |
| Tests | 997 passing, 3 skipped, lint and formatting clean |

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
docs/notes/      ~7,400 lines: equations, symbols, inventories, open questions
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

## What the real-time approximation actually costs

Table 1 prints two eigenvalue columns for what should be the same system — a 2-DOF model
and an order-reduced 5-DOF — and Ballin's differ, −2.69 against −2.81 at hover. **Ours are
identical**, and that is a theorem rather than a bug: linearizing an exactly solved
quasi-steady system and taking the Schur complement of the full Jacobian are the same
operation, so any model that converges its pressures must print two identical columns.

His differ because his 2-DOF is a *separately coded* nonlinear program (pdf p.27) carrying
the real-time numerics. Linearizing **our** `realtime.step` as a six-state discrete map —
including Eq. 74's carried mass flow, which the continuous model does not have — and
converting with `logm(A_d)/dt` puts a number on that:

| hover NG mode | value | vs Ballin's −2.69 |
|---|---|---|
| continuous 2-DOF / reduced-5 | −2.444 | −9.1 % |
| **discrete map at the report's own 7 ms frame** | **−2.728** | **+1.4 %** |

So most of the gap is the frame, not the physics, and Table 1 column 4 becomes an
independent check we pass at two trims of three. Refining the frame converges it back to
the continuous value first-order. Two by-products: **two** of the frame map's six modes are
discretization artifacts rather than dynamics — their *discrete* eigenvalues stay a fixed
fraction per **frame** as the frame shrinks, so neither has a continuous limit. P45's
one-frame lag runs 0.675 → 0.590 → 0.538 → 0.507 at dt = 14 / 7 / 3.5 / 1.75 ms, and
Eq. 74's carried compressor flow is a complex pair going the other way, 0.154 → 0.191.
(This README attributed the first series to Eq. 74 until 2026-09-13. By participation
factor that mode is P45's — P[P45] = 0.68 against P[WA31] = 0.07 at the 7 ms frame — and it
is a mode because Eq. 26 takes the *entering* P45, so it feeds forward a frame whatever the
iteration does. Raising the P45 pass cap from 8 to 400 leaves it at −75.330 to the digit.) and the printed 0.1 % iteration tolerance costs a further 5.2–8.7 % on the slow mode
(10–19 % before its stopping rule was corrected).

## Known gaps

- **The pressure loops were stopping on the wrong quantity, and correcting it split the
  two transients apart.** [pdf p.37] prints no convergence tolerance. What it prints is a
  cost — "eleven arithmetic operations required for each pass" — a pass count, and an
  outcome: "up to ten iterations may be required for convergence with **less than 0.1
  percent error**, resulting in a total of 110 arithmetic operations." `TOL_PRESSURE =
  1e-3` was that 0.1 percent applied to the iterate **step**. For a linearly convergent
  iteration the two differ by `ρ/(1−ρ)`, and ρ measures **0.887–0.902** here, so the test
  delivered about eight times the error it was named for: a median true P3 error of
  **0.25 %** on the Figure 9 accel and a worst frame of 4.21 %. Stopping on the error
  itself — estimated from the iterates, so no new constant enters — gives a median of
  **0.075 %**, and makes the loop take a mean of 4.1 passes with the ten-pass cap reached
  on 15.7 % of frames, which is the report's own budget rather than the *one* pass the
  step test was exiting on.

  **Figure 9 — the accel, and the report's own stated test case for this iteration —
  improved on four of five panels**: pcng 2.57 → 1.70 %, Ps3 2.14 → 1.29 %, torq45
  2.38 → 2.06 %, T41 and T45 unmoved. Figure 10's chop initially degraded badly, to
  T45 13.15 % with a floor of NGc 69.90 %, and that was attributed to `f1`'s data hole.
  **It was not `f1`.** See the next entry.

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
  15.8 % → 7.6 %; the Figure 9 T41 overshoot 1.67× Ballin's → **0.83×**. (Those are the
  numbers at that change. The stopping-rule correction above then took the mean to 4.4 %
  and the worst panel to 13.2 %, all of it on Figure 10 — and solving Eq. 80 took them back
  to **3.5 %** and **7.7 %**, better than either.)
- **What looked like the report disagreeing with itself at 775 lbm/hr was an unsettled
  transient, and we match both sides of it.** Figure 6 gives 99.69 %NG against Figure 9's
  98.88, Figure 8 gives Ps3 244.2 against 235.3, Figure 7 gives 1724.8 shp against 1648.2 —
  spreads of 0.82, 3.77 and 4.64 %. But **Ballin's Figure 9 trace is still climbing when its
  record ends** at t = 4.46 s — at +0.24 to +0.28 %NG/s measured over its last 1.5–3 s — so
  the two figures never described the same instant. (The trace is quantised at 0.030 %NG per
  pixel row, so the *instantaneous* end slope is not resolved; an earlier "+0.142 %NG/s" was
  one fit window's answer.) Our own model reproduces 50 / 38 / 61 % of each spread purely as
  settling, and matches the sweep to −0.26 / −1.92 / −0.91 % *and* the transient to +0.15 /
  +0.30 / +0.81 % at the same time — which a real contradiction would forbid.
- **The state trajectory is right; what is left is the rate along it.** Comparing Ps3 against
  NG instead of against time discards the unprinted step time, leaving the thermodynamic path
  through the state space. On **Figure 10** ours matches Ballin's to **1.30 % over 78–90 %NG
  against a 7.24 % chord baseline** — the chord being the straight line between the two
  endpoint trims, which the steady-state tests already pin, so being about seven times
  closer than it is real shape information. Figure 10's agreement also *improves* as the
  frame shrinks. `dNG/dt` at matched NG runs **0.92 / 1.01 / 0.97×** Ballin's at
  80 / 84 / 88 %NG and **0.76×** at 76 %. (It read 2.09× at 76 % for part of 2026-09-13,
  when Eq. 80's non-convergence was letting the chop run away downward; solving Eq. 80
  restored it.) At 80 %NG the net torque is **6.6 % of the turbine torque it is the
  difference of**, so 1 % on either term moves the rate 15 %.
  **Three things this does not show**, all measured rather than conceded: Figure 9's panel is
  worse than vacuous (his Ps3(NG) there is a straight line to R² = 0.9994, chord error 0.30 %,
  while our deviation is 1.66 % — larger than the curvature it would have to explain); the
  metric is *blind* to
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
- **Three fixes, each replacing an invention with something printed.** The pressure iterations
  ran to `tol=1e-10` where the report states 0.1 percent in ten and eight passes (pdf p.37);
  that 0.1 percent was then applied to the iterate step rather than the error the report
  describes, which under-converged the loop eightfold; and the heat sink is now Eqs. 48–49
  rather than the collapsed Eq. 50. None was aimed at a figure, and the middle one cost us
  Figure 10.
- **Eq. 80's iteration was the whole of the sub-idle problem, and it was also most of
  Figure 10's.** Below roughly 160 lbm/hr, f9's elasticity passes −1 and Eq. 80's printed
  fixed-point iteration stops converging: its fixed point becomes repelling, and eight
  passes land wherever they land. Everything that had been recorded as separate defects was
  that one fact.

  - Held at its own 125 lbm/hr trim the model settled at **75.743 %NG** against a 67.039 %
    differential trim, and which root it found depended on nothing but the **parity** of
    `MAX_ITER_P45` — 75.743 for even, 67.610 for odd, magnitude irrelevant out to 21 passes.
  - The frame map had its own sub-idle equilibria, up to **−7.07 %NG** at 150 lbm/hr,
    identical at tol 1e-3 and 1e-9 and dead still at the end.
  - Figure 10's chop ran away to a floor of NGc **69.90 %**, which was attributed to `f1`'s
    65–80 % data hole and cost a raise of the whole-curve ratchet from 8 to 14 %.

  Established constructively: solve the **same** printed equation, with the same f9 data,
  by a bracketed root-find instead of successive substitution, and all three vanish. The
  iteration now runs exactly as printed and falls back to bisection only when it fails to
  meet the report's own criterion — bit-identical where it converges (Figure 9 does not
  move; 400 lbm/hr exits on pass one). Results:

  | | before | after |
  |---|---|---|
  | sub-idle trims, 125–175 lbm/hr | −7.07 to −0.21 % | **≤0.003 %** |
  | dependence on `MAX_ITER_P45` | parity-selected, two roots | **none, caps 2–21 identical** |
  | Figure 10 floor | 69.90 %NGc | **74.08** (Ballin: 74.18) |
  | Figure 10 worst panel | 13.15 % | **7.72 %** |
  | nine-panel mean | 4.40 % | **3.51 %** |

  This is a **deliberate departure** from the method the report describes — [pdf p.37]
  accounts four arithmetic operations and one table lookup per pass, which is successive
  substitution specifically — taken because the printed method does not converge in a
  regime the report never ran: below-flight-idle fuel control was one of the features
  eliminated from the real-time model [pdf p.38], and the printed "eight iterations
  resulted in an error equal to less than 0.1 percent" was measured on a step from flight
  idle to full power. Where Ballin measured it, it holds and we reproduce it. The cost is
  the pass count: **~23 a frame** below 160 lbm/hr against the report's budget of eight,
  and 1–2 above flight idle. Open question #57.

  **And `f1`'s "data hole" turned out not to exist — see below.**

- **`f1` is interpolated on a beta grid now, and the "data hole" was never the problem.**
  Beta is the fraction of a speed line's own pressure-ratio span — 0 at the choked end, 1 at
  surge — so one beta names corresponding points on lines whose spans differ by a factor of
  six. Interpolating at constant *pressure ratio*, which is what this model did until
  2026-09-14, asks the shorter line for a ratio it cannot physically reach: the 65 % line
  stops at Ps3/P2 = 3.753 where the 80 % line reaches 6.671. Over the Figure 10 chop that
  blended a real number with a clamped one on **389 frames**, and it was recorded as open
  question #58, "f1's 15-point data hole".

  There is no hole. Blending at equal beta blends the surge limits too, so a query inside
  both lines is inside the blend: **zero of 715 frames leave the map.** The extrapolation was
  an artifact of the evaluation.

  **Ballin's own beta lines were tried and rejected on measurement.** Figure A1 prints six
  dotted construction lines joining the k-th marker of all eleven speed lines — those *are*
  beta lines, and they are his. But the marker positions carry digitizing noise: as a
  fraction of each line's span the k=1 markers run 0.785, 0.817, 0.831, 0.832, **0.814,
  0.811**, 0.855 … — non-monotone — while the spans themselves are smooth. Using the markers
  as the correspondence propagates that jitter into every interpolated value. An analytic
  beta scores better (Table B.1 0.2439 vs 0.2471 %, worst Table 1 mode 7.56 vs 8.18 %) and is
  far less grid-sensitive: 0.2438–0.2439 % across every grid from 32 to 96 points, against
  0.2471–0.2488 % for the printed markers.

  **This is a departure from replication and is recorded as one** (open question #60). The
  speed lines are *drawn* as polylines, so the seven-point table is Ballin's and linear
  interpolation between those points is his model — `maps.f1_as_printed()` still provides it.
  The characteristic they sample is smooth; the corners are an artifact of plotting a coarse
  table. The grid is produced by `tools/digitize_a1.py` (`--raw-only` skips it) so it carries
  its own provenance header and the reproducibility gate checks it like any other data file.

  | | constant-x | printed beta | **analytic beta** |
  |---|---|---|---|
  | frames reading `f1` off the map | 389 | 0 | **0** |
  | worst Table 1 NG mode | −22.6 % | +11.8 % | **+7.56 %** |
  | Table B.1 rms | **0.2018 %** | 0.2355 % | 0.2439 % |
  | Figure 10 worst panel | 13.15 % | 7.50 % | **7.50 %** |
  | stored points per speed line | 7 | 7 | 56 |

  The Table B.1 cost is the honest price: a 21 % degradation of the project's primary
  evidence, taken because the report disagrees with *itself* by 5.5 % on that quantity and a
  scheme that manufactures extrapolation is wrong whatever it scores. **What it did not fix
  is #43** — the derivative ambiguity moved from descent to level rather than going away.

- **Table B.1 and Figures 6–7 disagree with each other**, by up to 5.5 % on shaft power at
  low power. We track Table B.1, which is printed numbers rather than a plot.
- Figures 11–15 are **not reproducible** — they need the Gen Hel UH-60A blade-element
  simulation, which this report consumes and does not contain. That was Ballin's boundary
  too. **So Phase 5 has no published trace to overlay** — which was known before it
  started, not discovered at the end. The check that *is* available turns out to be a
  strong one, because the engine and the control were digitized from different halves of
  the report and neither knows the other exists: close the loop, hold the load at Table
  B.1's printed shaft torque, and the model settles on the printed trims **with fuel flow
  as an output** — worst **0.62 %** across NG, NP, Wf, Ps3 and shp at three conditions, on
  descent fuel flow. **The governor does not inherit the open-loop shaft-power error, it
  converts it**: descent shp is −0.724 % open-loop and **+0.101 %** closed, because fuel
  flow is an output and the loop trims it until the torque matches, absorbing the error
  into Wf instead. This paragraph said "worst 0.74 %, and that one is inherited (the
  descent shaft power is 0.72 % out open-loop too)" until 2026-09-13; the 2026-09-13
  accuracy audit rebuilt the tree at the commit that wrote it and measured 0.696 % with
  descent shp +0.115 %, so both the number and the attribution were wrong when written. The
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

Open questions are tracked in `docs/notes/open-questions.md` — **56 logged, 40 closed, 5 partly
closed, 11 open**. The 2026-09-13 audits closed #25 (the pressure loops' stopping test and
pass count) and #45 (Eq. 80's repelling fixed point), and opened #56 (the Jacobian's
perturbation step, which `SCOPE.md` had claimed was the report's ±2 % and never was).

**Seven of the eleven are things the report simply does not print**, and no amount of work
closes them: the initialization rule for the opened iteration (#23), the integration algorithm
(#26), the relaxation parameter and Lipschitz constant (#27), the 10 ms against 14 ms conflict
(#33), the power turbine speed behind Figures 6–8 (#35) and behind Figures 9–10 (#36), and the
fuel-step time (#37).

Of the three that remained work rather than gaps in the source, two are now closed. **#50**
(Figure C30's seven crossing curves) was extracted on 2026-09-13. **#45** — the model running
outside its own digitized envelope on the 775 lbm/hr step — turned out to be an artifact of
*prescribing* a fuel step the control would never command: with the loop closed, a harder
demand leaves `f1`, `f7` and `f9` entirely unclamped, because that is what the acceleration
limit is for. **#49**, linearizing the discrete real-time map, is now done — see below.

Phase 5 added **#54** (is a hysteresis constant the whole band or the half?) and **#55** (the
report states no initial condition for any of the 22 block diagrams). Neither can be closed
from the source, so both are **bounded by measurement instead**: #54 is worth ≤0.09 % on NG
and at worst 1.08 % on descent shaft power, and #55 ≤0.61 % on the settled state, because
the governor absorbs it. No conclusion here turns on either.

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
