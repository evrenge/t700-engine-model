# Scope — T700 Replication

What this project set out to reproduce, what it deliberately does not, and the tolerances
every comparison is judged against. Facts carry PDF page citations into NASA TM-100991;
`docs/notes/index.md` is the section map.

## Objective

Reproduce, in Python, the engine and control-system model of NASA TM-100991, and demonstrate
that it matches the results published in that report to a stated tolerance.

"Done" means: the report's result figures regenerate from this code, the deviations are
quantified, and every deviation is either within tolerance or explained.

The governing rule is **replication, not improvement**. Where the report is dated,
simplified or arguably wrong, it is reproduced as written and the objection recorded in
`docs/notes/open-questions.md`.

## What the report contains

This shaped the scope, so it is recorded rather than left implicit.

- A **5-state nonlinear model** — NG, NP, P3, P41, P45: two spool speeds plus three
  control-volume pressures [pdf p.27]. The fuel control adds **22 more states plus 4 memory
  elements**, so the assembled model is ~31 states.
- A **fixed-step real-time formulation** with a 10 ms maximum step, set by a 0.1 % allowable
  inter-step error, whose response is time-step dependent because the compressor mass-flow
  iteration is *opened* across the frame boundary [pdf p.38].
- **Appendix A** — engine model constants (Table A.1) and eleven function tables, printed
  only as plots [pdf pp.55–66].
- **Appendix B** — small-perturbation linear models: three trim conditions × four DOF
  variants, 12 figures, **297 printed matrix elements** [pdf pp.67–76]. The 2- and 3-DOF
  models are separately extracted, not order-reduced from the larger ones.
- **Appendix C** — the complete fuel control system: one constants table, 22 block diagrams
  and 8 function plots. **No equations at all** [pdf pp.77–101].
- The helicopter load enters as an input, `Qreq = Qmr + Qtr + Qacc + Qdamp` [pdf p.28], with
  rotor torques supplied by the **Gen Hel UH-60A blade-element simulation** — a separate
  program, not modelled in this report.

## In scope

- Gas path, component by component: inlet, compressor, combustor, gas generator turbine,
  free power turbine, exhaust.
- The three control volumes and their pressure states.
- Rotor dynamics: gas generator spool and free power turbine spool.
- The Appendix C fuel control — hydromechanical and electrical units, governing, scheduling,
  limiting.
- The report's fixed-step integration scheme, including the opened compressor iteration and
  its time-step sensitivity.
- The load interface: `Qreq` as an input, plus the load-torque derivatives the report uses
  for its linear models [pdf p.28].

## Out of scope

- **The Gen Hel UH-60A rotor and airframe model.** The report does not contain it; it
  consumes it. `Qreq` is a boundary condition to us exactly as the rotor model is an external
  supplier to the report. This is the report's own boundary, not one we imposed.
- Engine start, windmill and shutdown regimes — the real-time model eliminates fuel control
  below flight idle [pdf p.38].
- Other T700 variants (-701, -701C, -401) and later FADEC control units.
- Hardware I/O, real-time OS integration, pilot-in-the-loop display.

## Deferred

Designed for, not built: a C++ lookup-table thermodynamics backend (`t700.thermo` is an
interface from day one so it drops in without touching component code — build it when
profiling shows property evaluation dominating, not before); any JIT; and a wall-clock-paced
real-time execution harness.

## Status

All seven phases are complete. Each had a gate, and what each gate actually delivered is
below rather than what it was hoped to.

| Phase | Gate | Outcome |
|---|---|---|
| 0 — Ingest the report | We know precisely what data the report gives and what it withholds | Per-page text, a verified 106-entry symbol table, and a figure-by-figure inventory of all three appendices |
| 1 — Capture the data | Every function the model needs exists as data, or is logged as unavailable | All **11 engine maps** and all **8 control schedules** digitized, plus five result figures as ~8,100 reference points |
| 2 — Units and thermodynamics | Property values reproduce a worked example in the report | `t700.units`, and `t700.thermo` behind one interface |
| 3 — Steady-state engine | Gas path matches the report's trim conditions | Table B.1's full printed state, **rms 0.22 %** over 21 numbers |
| 4 — Dynamics | Open-loop response has the report's shape and time constants | Both heat-sink configurations Ballin published, as a switch; whole-curve rms **0.61–3.79 %** of panel excursion |
| 5 — Control system | Closed-loop transients match the report's figures | **Not achievable, and known before the phase started** — Figures 11–15 all need Gen Hel. The gate actually met is the one available: the loop settles on Table B.1's printed trims **with fuel flow as an output**, worst **0.27 %** |
| 6 — Validation | Deviations quantified against the report | All four printed surfaces compared; **all 206 non-zero Appendix B elements** carry a numeric per-element comparison, and all **27** printed eigenvalues |

The published report at `site/index.html` sets out what each of those means.

## Which source wins when the report disagrees with itself

**Table B.1 is primary.** It is printed *numbers*, so it carries no digitizing error of ours
at all; it underpins Appendix B's linear models, whose eigenvalues our own Jacobian
reproduces; and Table 1's modes come from the same trims. It is the most internally
corroborated data in the report.

**Figures 6–8 are secondary** — evidence of *trend across the operating range*, not of
absolute value. They are plots we digitized, so they carry our read error on top of whatever
they are.

**They no longer disagree by much, and this paragraph used to say otherwise.** It stated
+5.49 / +1.32 / −0.14 % on shaft power at the three trims and drew a rule from it: do not
tune below the report's own internal spread, because at descent that spread is 5.5 % and we
sit 0.90 % inside it. Those numbers predate the 2026-09-14 digitizer corrections (open
questions #46, #63, #64), which found the disagreement was ours and not the report's. They
were prose that nothing computed, so nothing caught them going stale.

`validation/report.py::internal_spread` computes it now, into `report.json`, from the same
digitized figures every other comparison uses:

| trim | Wf, lbm/hr | Table B.1 vs the figure | our deviation from Table B.1 | nearest printed marker |
|---|---|---|---|---|
| hover | 476.3 | −0.034 % shp, +0.002 %NG | −0.058 % shp, +0.001 %NG | 5.0 lbm/hr |
| level 80 kt | 349.3 | +0.133 % shp, −0.011 %NG | −0.077 % shp, +0.120 %NG | 0.6 lbm/hr |
| descent 80 kt | 267.7 | **+1.306 %** shp, +0.020 %NG | **−0.806 %** shp, −0.048 %NG | 2.8 lbm/hr |

**The consequence changes with them, and it has now changed twice.** The report's widest
self-disagreement is 1.306 % on shaft power at descent, not 5.5 % and not the 0.734 % this
table carried between 2026-09-14 and 2026-09-18. Our own 0.806 % at that trim has not moved
at all across any of it — so the same number was outside the report's spread under the old
reading and is **inside** it under the current one, without anything about our model
changing. Three of the six comparisons sit inside; descent shaft power is one of them.

Take the lesson rather than the number. This row has now been wrong in two opposite
directions, both times because the figures it rests on were re-digitized underneath it —
2026-09-14's corrections, then the 2026-09-16/18 recalibrations. `internal_spread` computes
it into `report.json` so the *value* cannot go stale, but the paragraph around it is still
prose, and prose is what went stale both times. The rule it exists to state is unaffected
and still holds: do not tune below the source's own spread. Whether this particular residual
shelters under that rule is a fact about the digitization, not about the model, and it will
move again.

Read the last column before quoting a row. Figure 7 prints 15 markers over its whole range
and the value at a trim is interpolated between them; on its steep low-power end the curve
climbs ~3 shp per lbm/hr, so hover's 4.9 lbm/hr of interpolation is itself worth 1.6 % of
hover power. Only the level row is close to a direct read.

## Defects in the source, and four that were ours

Reading a 1988 scan carefully turns up genuine errors in it. Under the replication rule they
are reproduced as printed and recorded, never silently corrected. Each cites its ledger row.

| Defect | pdf p. | What we do | Row |
|---|---|---|---|
| Figure A3's `0.11` axis label is misplaced by 17.0 px, where every other label on the page sits within 5 px | 58 | read the plateau as 0.1091 | #42 |
| `TC_T41`'s units are printed `sec^(9/5)` in two places; Eq. 51 requires `sec^(1/5)` | 55 | units wrong, value 0.29 right and used | #4 |
| Figure A10's ordinate label `PS9/P49` is inverted relative to its own values | 65 | use the plotted values as `f10` directly | #5 |
| Eq. 29 does not conserve energy, by construction and by the report's own words — `H45 = K_H45·H44` is one multiplicative fraction, not a mass-weighted mix | 23 | reproduced as printed; costs −1.52 % on the overall energy balance | #52 |

**Four more looked like the report's and were ours**, which is the more useful list: a torque
threshold that seemed absurd until we noticed we had read rotor-hub torque as engine torque
(#11); an exhaust pressure ratio we nearly inverted (#5); Eq. 74 "contradicting its own
prose", which was a carry bug in our own alternative branch (#22); and a 3–6 % inconsistency
between Figures 8 and 10 that came from comparing a decelerating engine against an
equilibrium locus (#46, and "Comparing a transient to a steady-state figure" below).

## What the real-time approximation costs

Table 1 prints two eigenvalue columns for what should be the same system — a 2-DOF model and
an order-reduced 5-DOF — and Ballin's differ, −2.69 against −2.81 at hover. **Ours are
identical**, and that is a theorem rather than a bug: linearizing an exactly solved
quasi-steady system and taking the Schur complement of the full Jacobian are the same
operation, so any model that converges its pressures must print two identical columns.

His differ because his 2-DOF is a separately coded nonlinear program [pdf p.27] carrying the
real-time numerics. Linearizing our `realtime.step` as a six-state discrete map — including
Eq. 74's carried mass flow, which the continuous model does not have — puts a number on it:

| NG mode, per second | hover | level | descent |
|---|---|---|---|
| Table 1, printed | −2.69 | −2.23 | −1.82 |
| continuous 2-DOF / reduced-5 | −2.539 (−5.6 %) | −2.399 (+7.6 %) | −1.897 (+4.2 %) |
| **discrete map at the report's own 7 ms frame** | **−2.605 (−3.2 %)** | **−2.321 (+4.1 %)** | **−1.802 (−1.0 %)** |

So part of the gap is the frame, not the physics — the discrete map is closer at all three
trims — which is what makes Table 1's fourth column an independent check rather than a
puzzle. Two of the frame map's six modes are discretization artifacts rather than dynamics:
their discrete eigenvalues stay a fixed fraction per *frame* as the frame shrinks, so
neither has a continuous limit. Bounded by `test_discrete_map.DISCRETE_TOL_PCT` above.

## Tolerances

**The report states almost none.** Six numbers exist in the whole document, and most are
deviations Ballin *accepted* rather than accuracy he *met*:

| Claim | Value | Scope | pdf p. |
|---|---|---|---|
| Max inter-step error | 0.1 % | sets the 10 ms step cap | 38 |
| Fuel consumption overestimate | ≤5 % | only 81–86 % NG | 39 |
| NG overestimate, hardware case | 4 % | lowest-power case only | 39 |
| NG overestimate, open-loop step | 1–2 % | Figs. 9–10 | 39 |
| Rotor speed agreement | 0.2 % | Fig. 12 only | 50 |
| Eigenvalue mismatch called "good" | 4 % | Table 1 comparison | 29 |
| NG validity ceiling | 100 % | above this the model is not claimed | 39 |

Everything else he says is qualitative. So our tolerances are **declared by us and labelled
as ours**. Provisionally, for quantities the report shows in a reproducible figure:

**Every number in this table is now one the suite asserts, and names the symbol that holds
it.** Until 2026-09-19 it named none, so the four bounds a reader is most likely to take as
our declared accuracy were enforced by nothing — and two of them were not the operative
number at all: the gas path is held an order tighter than the ±3 % printed here, and the
±5 % on fuel flow was asserted nowhere.

| Quantity | Steady state | Transient | Basis |
|---|---|---|---|
| NG | ±2 % | ±2 % | `test_trim_points.NG_TOL_PCT`, asserted at all three trims. Ballin's own accepted 1–4 % overestimates |
| Torque / SHP | ±3 % | ±3 % | `test_trim_points.SHP_TOL_PCT` and `test_trim_points.TORQUE_TOL_PCT`. Ours, undeclared by the report |
| Gas path P, T | ±0.5 % | ±1 % | `test_trim_points.STATE_TOL_PCT` and `test_fuel_step.FINAL_TOL_PCT`. This row asked for ±3 % until 2026-09-19 and nothing asserted it; the printed stations are held to 0.5 % and the settled transient to 1 %, both far tighter, so the looser band was never the operative number |
| Fuel flow | see the 81–86 %NG row below | — | The ±5 % this row printed until 2026-09-19 was asserted nowhere. What *is* asserted is the report's own printed claim about itself, in the per-comparison table |
| Appendix B elements | see the per-block bounds below, **not** a single number | — | his own 4 % "good agreement" |

### Per-comparison bounds

Every number in the Tolerance column is checked against the constant the code asserts by
`tests/test_scope_tolerances.py`. Widening a bound is allowed; doing it in one place only is
not.

| Comparison | Tolerance | Basis |
|---|---|---|
| Phase-plane state trajectory, Ps3 vs NG | +/-2.5 % | `test_figure_consistency.PHASE_PLANE_TOL_PCT`. Measured **0.81 %** worst on Fig. 9 and **1.78 %** on Fig. 10, rms 0.67 and 1.20, over the grids the test asserts. These have moved four times: the … |
| Settled state vs the Fig. 6-8 sweep | +/-0.5 % | `test_figure_consistency.BOTH_SIDES_TOL_PCT`. **Tightened 3 -> 0.5 on 2026-09-14**, never widened: the two digitizer corrections of that day (#63, #64) took the worst of the six comparisons from 1.92 … |
| Fig. 9/10 initial trim vs the figure | +/-1.5 % | `test_fuel_step.INITIAL_TOL_PCT`; this is the same physics the steady tests already check |
| Steady sweep vs Figure 6's real-time series | rms 0.35 %NG, worst 1.25 %NG | `test_steady_sweeps.FIG6_RMS_PCT_NG` and `test_steady_sweeps.FIG6_WORST_PCT_NG`. **New 2026-09-14.** Figure 6 was the one sweep the model did not track, at a whole-curve -2.0 to -0.4 %NG. That was the digitizer reading … |
| Fig. 9/10 settled state vs the figure | +/-1 %, both figures | `test_fuel_step.FINAL_TOL_PCT`. Tightened 20 -> 5/8 when the printed convergence criterion replaced an invented tolerance, and **5/8 -> 1 on 2026-09-14**, never widened. The 8 was covering a … |
| Whole-curve RMS, normalised by each panel's excursion | **4 % ratchet** | `test_whole_curve.WHOLE_CURVE_RMS_CEILING_PCT`. Not a tolerance -- set just above the worst measured panel so a regression fails and an improvement is free. Current over ten panels, **riser … |
| Step-edge jump across the riser hole | +/-12 % of panel excursion | `test_whole_curve.STEP_JUMP_TOL_PCT`. **New 2026-09-14.** The reference traces have a 64-165 ms sampling hole at the fuel step -- a line tracer cannot follow a vertical edge -- and `np.interp` drew a … |
| Table B.1's printed station state | +/-0.5 % | `test_trim_points.STATE_TOL_PCT`. Much tighter than the +/-3 % gas-path row because these are printed *numbers*, with no read error of ours in them at all |
| Figures 9/10 read error | +/-0.6 % of full scale | `test_figure_consistency.READ_ERROR_CEILING_PCT_FS` -- see the section above, and note the currency. **Tightened 1.6 -> 0.6 on 2026-09-14**, never widened: 1.6 was covering an uncorrected page skew … |
| HMU collective at a Table B.1 trim | 5-95 % of maximum | `test_hmu_trim.COLLECTIVE_RANGE_PCT`. The report prints no collective for these trims, so this is a believability band, not a comparison. What it tests is that Appendix C and Table B.1 -- digitized … |
| **Closed loop vs Table B.1** | **+/-0.5 %** | `test_closed_loop.CLOSED_LOOP_TOL_PCT`. Measured worst **0.245 %**, descent fuel flow, as a **mean over the limit cycle** -- the loop settles to a cycle, not a point (#61), and the terminal sample … |
| **Table 1 eigenvalues** | **+/-8 %** | `test_trim_points.EIGENVALUE_TOL_PCT`. **A characterization, not an accuracy claim** — a regression detector set above the measurement, currently **4.70 %** worst over all 27 printed modes (rms 1.75 %) with the report's ±2 % extraction step [pdf p.28], adopted 2026-10-12, when this tightened 30 -> 8. On the tangent step used before, the worst was 14.87 %, carrying the #31/#43 knot residual. Ballin calls a 4 % mismatch "good agreement" for his own comparison [pdf p.29]; ours is not that, and this row exists so the difference is visible rather than implied. **Undocumented until 2026-09-19**, while its own docstring claimed SCOPE.md declared it — the loosest bound in the project, on the one surface that is a restatement of another |
| Fuel overestimate, 81-86 %NG, vs the report's own claim | +/-8 % | `test_steady_sweeps.FUEL_OVERESTIMATE_CEILING_PCT`. The report's sharpest printed steady-state claim about itself — "overestimate fuel consumption by as much as five percent" [pdf p.39] — and nothing checked it until 2026-09-14. A ceiling above a measured 6.31 %, against Ballin's own 5.99 %; ours exceeds his claim, which is the finding. Documented here 2026-09-19 |
| Discrete real-time map vs Table 1 col. 4 | +/-8 % (hover), +/-15 % (level), +/-7 % (descent) | `test_discrete_map.DISCRETE_TOL_PCT`. The loose one carries the #31/#43 interpolation residual, which the continuous model has too and which this comparison is not measuring. **It was descent and is … |


**A third category exists beyond accuracy claims and structural pins: characterizations.**
A bound fitted around the current measurement is a regression detector wearing a tolerance's
clothes. They are legitimate, and they are *not* evidence of agreement with the report — the
distinction matters when a number from one is quoted as a result.

### Appendix B, element by element

There is no single per-element tolerance. Twenty-four distinct bounds exist across
`test_appendix_b_elements.py`, `test_all_dof_models.py` and `test_heat_sink_taus.py`, and
they differ by two orders of magnitude because the printed elements do:

| block | bound | test |
|---|---|---|
| `b`, the fuel column, 5-DOF | +/-0.2 % | `test_appendix_b_elements.FUEL_COLUMN_TOL_PCT` |
| 6-DOF T41 column | +/-1.0 % | `test_appendix_b_elements.SIX_DOF_T41_COLUMN_TOL_PCT` |
| P3 / P41 pressure block, 5-DOF | +/-3.0 % | `test_appendix_b_elements.PRESSURE_BLOCK_TOL_PCT` |
| P3 / P41 pressure block, 6-DOF | ratio 1.00-1.02 | `test_appendix_b_elements.SIX_DOF_PRESSURE` |
| P3 column outside that block, 5-DOF | +/-8.0 % | `test_appendix_b_elements.P3_COLUMN_TOL_PCT` |
| T41 row | ratio 0.90-0.96, spread under 1.0 % | `test_appendix_b_elements.T41_ROW` and `test_appendix_b_elements.T41_ROW_SPREAD_PCT` |
| 3-DOF T41 column | ratio 0.88-1.12 | `test_all_dof_models.THREE_DOF_T41_COLUMN` |
| 3-DOF T41 feedthrough `d` | +/-7.0 % | `test_all_dof_models.THREE_DOF_FEEDTHROUGH_TOL_PCT` |
| 3-DOF NG column, less descent `A(T41,NG)` | ratio 0.76-1.12 | `test_appendix_b_elements.THREE_DOF_NG_COLUMN` |
| descent `A(T41,NG)`, 3-DOF | ratio 0.6-1.0 | `test_all_dof_models.DESCENT_T41_NG` |
| 3-DOF `b`, NG and NP | ratio 0.89-0.99 | `test_appendix_b_elements.THREE_DOF_B_GAS` |
| 3-DOF `b`, T41 | ratio 0.74-0.88 | `test_appendix_b_elements.THREE_DOF_B_T41` |
| heat-sink lead and lag taus, ours over those B7/B9/B11 imply | ratio 0.84-0.88 and 0.89-0.93 | `test_heat_sink_taus.TAU1_RATIO` and `test_heat_sink_taus.TAU2_RATIO` |
| 6-DOF fuel column with those taus | +/-0.5 % | `test_heat_sink_taus.SIX_DOF_SWAP_TOL_PCT` |
| 3-DOF T41 row with those taus, hover and descent | +/-2.0 % | `test_heat_sink_taus.THREE_DOF_SWAP_TOL_PCT` |
| 3-DOF T41 row with those taus at the report's step, all three trims | +/-8.0 % | `test_heat_sink_taus.SHIPPED_THREE_DOF_SWAP_TOL_PCT` |
| TC_T41 those taus imply, against the printed 0.29 | 0.315-0.325, and f_hs within +/-1.0 % of its plot | `test_heat_sink_taus.IMPLIED_TC_T41` and `test_heat_sink_taus.IMPLIED_F_HS_TOL_PCT` |
| 6-DOF `b`, gas path | ratio 0.93-0.96 | `test_appendix_b_elements.SIX_DOF_B_GAS` |
| 6-DOF `b`, T41 | ratio 0.93-0.97 | `test_appendix_b_elements.SIX_DOF_B_T41` |
| rest of the 6-DOF `A` | ratio 0.8-1.26 | `test_appendix_b_elements.SIX_DOF_REST` |
| `A(NP,NP)` without dQreq/dNP | -56 to -50 % | `test_appendix_b_elements.NP_DIAGONAL_RESIDUAL_PCT` |
| f7 group, by trim | +/-8 / 12 / 5 % | `test_appendix_b_elements.F7_GROUP_TOL_PCT` |
| 5-DOF NP row, with the derived load inertia | +/-17.0 % | `test_appendix_b_elements.NP_ROW_TOL_PCT` |
| 2-DOF elements | ratio 0.70-1.20 `A`, +/-6.0 % `b` | `test_all_dof_models.TWO_DOF_A` and `test_all_dof_models.TWO_DOF_B_TOL_PCT` |


**Coverage: 297 printed elements, 206 of them non-zero, and all 206 carry a numeric
per-element comparison.** Two elements are excluded rather than bounded, each with a test
recording why — the 6-DOF `A(P45,P3)`, which is a near-cancellation printed to four figures,
and the `A(NP,NP)` diagonals, which are Gen Hel's.

### The figures' own read error

Figures 9 and 10 each contain one panel whose true value is printed — the fuel flow is the
*input*, and both its levels are stated in the caption — so digitizing it measures the read
error of those figures directly. It is a **pixel offset**: a fixed fraction of each panel's
span, not of the value it lands on.

| panel | Fig. 9 floor | Fig. 10 floor |
|---|---|---|
| WFPH | 1.83 % | 1.67 % |
| PS3 | 1.24 % | 1.70 % |
| **PCNG** | **0.21 %** | **0.59 %** |
| T41 | 0.44 % | 0.61 % |
| T45 | 0.62 % | 0.85 % |
| TORQ45 | 2.22 % | 2.91 % |

The currency matters more than the number: **on PCNG the floor is 0.21 %**, nine times
tighter than a single per-page figure would suggest, and PCNG is what most transient
agreements are quoted on. The evidence for the pixel model is a prediction — the two figures
read the *same* 400 lbm/hr trim on axes of different span, and calibrating only on the fuel
panels predicts **+0.340 %NG** on PCNG against a measured **+0.340**. Pinned by
`validation/test_figure_consistency.py`.

### Comparing a transient to a steady-state figure

**A transient trace must not be compared against a steady-state locus at matched speed.**
Figures 6–8 are equilibrium sweeps; Figures 9–10 are transients. During a chop the choked
station 4.1 nozzle passes the same flow at a lower pressure once the temperature falls, so
Ps3 sits *below* its equilibrium value at that NG by construction — Ballin's own Figure 10
runs −1.68 % to −8.73 % against his own Figure 8. The effect is quasi-steady, not a pressure
lag: the pressures are fast and the spool is slow.

The valid comparisons are trace against trace, settled state against the sweep, and the
**phase plane** — Ps3 against NG, which discards the time axis and with it the unprinted step
time.

## Validation surfaces available, in order of what they are worth

The ordering is **evidence that is printed as numbers and localizes a disagreement outranks
evidence that is digitized or aggregated.** Two consequences do the work: a printed number
carries none of our own read error, and an element-wise comparison says *which* partial
derivative moved where a mode only says that something did.

| | Surface | pdf p. | Scale | our read error | what only this surface constrains |
|---|---|---|---|---|---|
| 1 | **Table B.1** | 67 | 21 numbers | none | the equilibrium every other comparison is taken about. If this is wrong, every derivative below is evaluated at the wrong point |
| 2 | **Appendix B** | 68–76 | 297 elements, 206 non-zero | none | every partial derivative, **localizably**. The three error families in `test_appendix_b_elements.py` — the missing load inertia, the NP diagonal residual, the `f7` group — were each found here and none could have been found from a mode |
| 3 | **Figures 9–10** | 45–46 | ~8,000 points, 6 panels each | ±0.6 % of full scale | the integrated time response: volume dynamics and heat-sink lags over a whole trajectory, which no printed table touches |
| 4 | **Closed loop vs Table B.1** | — | 3 trims | none | the fuel control. No printed surface covers it; this is internal consistency, not agreement with the report |
| 5 | **Figures 6–8** | 40–42 | ~70 markers | ours | steady behaviour across the **whole operating range** rather than three points, so a map error shows as a trend in a region instead of one unattributable number |
| 6 | **Table 1 eigenvalues** | 31 | 12 comparable of 27 | none | *nothing the rows above do not already carry.* The modes are an aggregation of Appendix B's matrices at Table B.1's trims, and they cover only the 2-DOF, 5-DOF and reduced-5 models — the 3-DOF and 6-DOF print no eigenvalues at all |
| — | Figures 11–15 | 48–53 | — | — | **not reproducible**: they need Gen Hel |
| — | Table 3 | 44 | 6 cases | — | **not a target**; see below |

Two things this ordering is not. It is **not** a ranking of how hard each comparison is to
pass — the loosest bound in the project sits on row 6 and one of the tightest on row 1. And
it is **not** an instruction to stop running the lower rows: row 6 is cheap, and a mode that
moves when no element did would be a finding about our own linearization.

Where it bites is when two rows disagree, which has happened twice on the digitisation
branch: a calibration that improved rows 1, 3 and 4 while degrading row 2 was taken, and the
per-figure frame-model choice that improved row 1 while costing one element of row 2 was also
taken, with the cost recorded in open question #73 rather than tuned away.

### Two things not to do

- **Table 3 (pdf p.44) is not a validation target.** It looks like the best steady-state
  benchmark in the report — six hardware test cases — but its model column was produced with
  compressor and turbine functions from the *NASA-Lewis prototype engine*, supplied in place
  of the standard ones and printed nowhere. Its own caption says it does not represent
  specification T700 performance. Use it as evidence of method, never as pass/fail.
- **Do not fix the T41/T45 trim bias.** Ballin's model underestimates turbine temperature at
  trim in every validation case and he says so [pdf p.54]. Our replication should reproduce
  that bias with the same sign. A version that matches the hardware better than his did is a
  worse replication.

## Admitted limits that become ours

From the Conclusions [pdf p.54]: the 100 %NG ceiling; invalid at low NP; T41/T45 trim
underestimated throughout; direction-independent heat-sink constants; no T45 heat sink; no
inlet guide vane dynamics; damping slightly high; load-demand compensation weaker than the
real aircraft.

One positive result worth keeping in view: omitting high-speed inter-volume mass-flow
dynamics "was found to be unnecessary", which is what justifies the 5-state formulation in
the first place.

## Open questions

Tracked in `docs/notes/open-questions.md` — 66 logged, 55 closed, 4 partly closed, 7 open,
with the tally checked by `tests/test_open_questions_ledger.py` rather than maintained by
hand. Numbering is deliberately non-contiguous; never renumber a row. Each row states what
was asked and what was decided; the investigations that reached those decisions are in git
history, and that test fails if a closed row states a verdict without a finding.

Anything the report does not answer goes there rather than into the code as a guess.
