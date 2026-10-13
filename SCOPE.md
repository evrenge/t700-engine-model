# Scope — T700 Replication

The objective: reproduce, in Python, the engine and control-system model of NASA TM-100991,
and show that it matches the results published there to a stated tolerance. "Done" means the
report's result figures regenerate from this code, the deviations are quantified, and every
deviation is within tolerance or explained. The rule is **replication, not improvement**:
where the report is dated, simplified or arguably wrong, it is reproduced as written and the
objection recorded in `docs/notes/open-questions.md`. `docs/notes/index.md` is the section map.

## What the report contains

- A **5-state nonlinear model** (NG, NP, P3, P41, P45: two spool speeds, three
  control-volume pressures) [pdf p.27]. The fuel control adds **22 more states plus 4 memory
  elements**, so the assembled model is ~31 states.
- A **fixed-step real-time formulation**: a 10 ms maximum step, set by a 0.1 % allowable
  inter-step error, and a response that depends on the step because the compressor mass-flow
  iteration is *opened* across the frame boundary [pdf p.38].
- **Appendix A**: Table A.1's constants and eleven function tables, printed only as plots
  [pdf pp.55–66]. **Appendix B**: linear models at three trims × four DOF variants, **297
  printed matrix elements** [pdf pp.67–76]; the 2- and 3-DOF models are separately extracted,
  not order-reduced. **Appendix C**: the fuel control as one constants table, 22 block
  diagrams and 8 function plots, with **no equations at all** [pdf pp.77–101].
- The load enters as `Qreq = Qmr + Qtr + Qacc + Qdamp` [pdf p.28], with rotor torques from
  the **Gen Hel UH-60A blade-element simulation**, a separate program the report does not
  contain.

## Scope

**In:** the gas path component by component (inlet, compressor, combustor, gas generator
turbine, free power turbine, exhaust) and its three control volumes; both spools; the
Appendix C fuel control (hydromechanical and electrical units, governing, scheduling,
limiting); the report's fixed-step scheme with its opened iteration and step sensitivity; and
`Qreq` as an input, with the load-torque derivatives the linear models use [pdf p.28].

**Out:** the Gen Hel rotor and airframe model (the report consumes it, so `Qreq` is a
boundary condition to us exactly as it is to the report); start, windmill and shutdown, since
the real-time model eliminates fuel control below flight idle [pdf p.38]; other T700 variants
and FADEC units; hardware I/O, real-time OS integration, pilot display.

**Deferred:** a C++ lookup-table thermodynamics backend (`t700.thermo` is an interface so it
can drop in; build it when profiling shows property evaluation dominating), any JIT, and a
wall-clock-paced real-time harness.

## Status

All seven phases are complete. Each had a gate; the last column is what each delivered.

| Phase | Gate | Outcome |
|---|---|---|
| 0 — Ingest the report | We know precisely what data the report gives and what it withholds | Per-page text, a verified symbol table, and a figure-by-figure inventory of all three appendices |
| 1 — Capture the data | Every function the model needs exists as data, or is logged as unavailable | All **11 engine maps** and all **8 control schedules** digitized, plus the five result figures 6–10 as reference traces |
| 2 — Units and thermodynamics | Property values reproduce a worked example in the report | `t700.units`, and `t700.thermo` behind one interface |
| 3 — Steady-state engine | Gas path matches the report's trim conditions | Table B.1's full printed state, **rms 0.22 %<!--rj:table_b1.rms-->** over 21<!--rj:table_b1.n--> numbers |
| 4 — Dynamics | Open-loop response has the report's shape and time constants | Both heat-sink configurations Ballin published, as a switch; whole-curve rms **0.62<!--rj:whole_curve_summary.min-->–3.80 %<!--rj:whole_curve_summary.max-->** of panel excursion |
| 5 — Control system | Closed-loop transients match the report's figures | **Not achievable, and known before the phase started** — Figures 11–15 all need Gen Hel. The gate met instead: the loop settles on Table B.1's printed trims **with fuel flow as an output**, worst **0.25 %<!--rj:closed_loop.worst-->** |
| 6 — Validation | Deviations quantified against the report | All four printed surfaces compared; **all 206 non-zero Appendix B elements** carry a numeric per-element comparison, and all **27**<!--rj:eigenvalues.notes.table_1_stats.n--> printed eigenvalues |

## Which source wins when the report disagrees with itself

**Table B.1 is primary**: printed *numbers*, free of our read error, and the trims that
Appendix B's linear models and Table 1's modes are taken at, so it is the most corroborated
data in the report. **Figures 6–8 are secondary**: evidence of *trend across the operating
range*, not of absolute value, carrying our read error on top of their own. Where the two
disagree, do not tune below the report's own spread, which
`validation/report.py::internal_spread` computes into `report.json`:

| trim | Wf, lbm/hr | Table B.1 vs the figure | our deviation from Table B.1 | nearest printed marker |
|---|---|---|---|---|
| hover | 476.3 | −0.028 %<!--rj:internal_spread.rows.1.theirs_pct--> shp, +0.000 %NG<!--rj:internal_spread.rows.0.theirs_pct--> | −0.060 %<!--rj:internal_spread.rows.1.ours_pct--> shp, +0.000 %NG<!--rj:internal_spread.rows.0.ours_pct--> | 5.0 lbm/hr<!--rj:internal_spread.rows.1.nearest_marker_pph--> |
| level 80 kt | 349.3 | +0.148 %<!--rj:internal_spread.rows.3.theirs_pct--> shp, −0.015 %NG<!--rj:internal_spread.rows.2.theirs_pct--> | −0.077 %<!--rj:internal_spread.rows.3.ours_pct--> shp, +0.120 %NG<!--rj:internal_spread.rows.2.ours_pct--> | 0.6 lbm/hr<!--rj:internal_spread.rows.3.nearest_marker_pph--> |
| descent 80 kt | 267.7 | **+1.329 %<!--rj:internal_spread.rows.5.theirs_pct-->** shp, +0.015 %NG<!--rj:internal_spread.rows.4.theirs_pct--> | **−0.806 %<!--rj:internal_spread.rows.5.ours_pct-->** shp, −0.048 %NG<!--rj:internal_spread.rows.4.ours_pct--> | 2.7 lbm/hr<!--rj:internal_spread.rows.5.nearest_marker_pph--> |

The report's widest self-disagreement is 1.329 %<!--rj:internal_spread.widest_theirs_pct-->
on shaft power at descent, and our shaft power there lies between the table and the figure
(#69). Whether a residual shelters under this rule is a fact about the digitization, not the
model: re-digitizing the figures moves the spread. Read the last column before quoting a row:
a trim's value is interpolated between printed markers, so only the level row is close to a
direct read.

## Defects in the source, and four that were ours

Under the replication rule the report's own errors are reproduced as printed and recorded,
never silently corrected. Each cites its ledger row.

| Defect | pdf p. | What we do | Row |
|---|---|---|---|
| Figure A3's `0.11` axis label sits well below its tick, the only misplaced label on the page | 58 | read the plateau as 0.1091 | #42 |
| `TC_T41`'s units are printed `sec^(9/5)` in two places; Eq. 51 requires `sec^(1/5)` | 55 | units wrong, value 0.29 right and used | #4 |
| Figure A10's ordinate label `PS9/P49` is inverted relative to its own values | 65 | use the plotted values as `f10` directly | #5 |
| Eq. 29 does not conserve energy, by construction and by the report's own words — `H45 = K_H45·H44` is one multiplicative fraction, not a mass-weighted mix | 24 | reproduced as printed, so the overall energy balance does not close (#52 measures by how much) | #52 |

**Four more looked like the report's and were ours**: a torque threshold that seemed absurd
because we had read rotor-hub torque as engine torque (#11); an exhaust pressure ratio we
nearly inverted (#5); Eq. 74 seeming to contradict its own prose, which was a bug in our own
alternative branch (#22); and an inconsistency between Figures 8 and 10 that came from
comparing a decelerating engine against an equilibrium locus (see "Comparing a transient to
a steady-state figure" below).

## What the real-time approximation costs

Table 1 [pdf p.31] prints three linear models, made three ways [pdf pp.27-28]: the 5-DOF by
±2 % central differences of the complete simulation, the 2-DOF by the same on "a
reduced-order nonlinear simulation which contained the volume dynamics approximation", and
the reduced-order 5-DOF by setting the 5-DOF's pressure derivatives to zero. `t700.linear`
makes each of ours the same way, and `test_all_dof_models.NG_MODE_TOL_PCT` holds each NG
mode to its own column. With exact derivatives the 2-DOF and the reduction would be one
model (`tests/test_linear.py::test_two_dof_equals_reduced_five_at_a_tangent_step`); at the
±2 % step they part by a few percent, Ballin's and ours alike (#49).

`validation/test_discrete_map.py` measures what our own 7 ms frame [pdf p.47] costs. It
linearizes `realtime.step` as a six-state map, including Eq. 74's carried flow, and compares
it with our continuous model. Refining the frame converges it on the continuous model at
first order; two of its six modes are frame artifacts with no continuous limit; the frame
makes the NG mode faster by under `test_discrete_map.FRAME_COST_TOL_PCT`. It is not a check
on Table 1: Ballin read his derivatives "at the input of each integrator" [pdf p.28], not
across a frame.

## Tolerances

**The report states almost none.** These are the numbers it gives, and most are deviations
Ballin *accepted* rather than accuracy he *met*:

| Claim | Value | Scope | pdf p. |
|---|---|---|---|
| Max inter-step error | 0.1 % | sets the 10 ms step cap | 38 |
| Fuel consumption overestimate | ≤5 % | only 81–86 % NG | 39 |
| NG overestimate, hardware case | 4 % | lowest-power case only | 39 |
| NG overestimate, open-loop step | 1–2 % | Figs. 9–10 | 39 |
| Rotor speed agreement | 0.2 % | Fig. 12 only | 50 |
| Eigenvalue mismatch called "good agreement" | 4 % | Table 1 comparison | 29 |
| NG validity ceiling | 100 % | above this the model is not claimed | 39 |

Everything else he says is qualitative, so our tolerances are **declared by us and labelled
as ours**. Every bound below names the symbol the suite asserts, and
`tests/test_scope_tolerances.py` fails if the two differ. Widening a bound is allowed;
doing it in one place only is not. Some bounds are **characterizations**, set just above
today's measurement to catch a regression: they are *not* evidence of agreement with the
report, so never quote one as a result.

| Quantity | Steady state | Transient | Basis |
|---|---|---|---|
| NG | ±2 % | ±2 % | `test_trim_points.NG_TOL_PCT`, asserted at all three trims; Ballin's own accepted 1–4 % overestimates. The transient column is backed by that trim symbol only: no transient test asserts it |
| Torque / SHP | ±3 % | ±3 % | `test_trim_points.SHP_TOL_PCT` and `test_trim_points.TORQUE_TOL_PCT`, at the trims; ours, undeclared by the report. As for NG, only the trim symbols back the transient column |
| Gas path P, T | ±0.5 % | ±1 % | `test_trim_points.STATE_TOL_PCT` on the printed stations and `test_fuel_step.FINAL_TOL_PCT` on the settled transient |
| Fuel flow | see the 81–86 %NG row below | — | What is asserted is the report's own printed claim about itself, in the per-comparison table |
| Appendix B elements | see the per-block bounds below, **not** a single number | — | his own 4 % "good agreement" [pdf p.29] |

### Per-comparison bounds

| Comparison | Tolerance | Basis |
|---|---|---|
| Phase-plane state trajectory, Ps3 vs NG | +/-2.5 % | `test_figure_consistency.PHASE_PLANE_TOL_PCT`. Measured **0.81 %<!--rj:phase_plane.figure_9.worst-->** worst on Fig. 9 and **1.78 %<!--rj:phase_plane.figure_10.worst\|abs-->** on Fig. 10, rms 0.67<!--rj:phase_plane.figure_9.rms--> and 1.27<!--rj:phase_plane.figure_10.rms-->, at the measured step times (#37), over the grids the test asserts |
| Settled state vs the Fig. 6-8 sweep | +/-0.5 % | `test_figure_consistency.BOTH_SIDES_TOL_PCT`. At 775 lbm/hr the sweeps and Figure 9 disagree because Figure 9 has not settled when its record ends; the model matches both at once, its settled state against the sweep and its state at the record's end against Figure 9 |
| Fig. 9/10 initial trim vs the figure | +/-1.5 % | `test_fuel_step.INITIAL_TOL_PCT`; this is the same physics the steady tests already check |
| Steady sweep vs Figure 6's real-time series | rms 0.35 %NG, worst 1.25 %NG | `test_steady_sweeps.FIG6_RMS_PCT_NG` and `test_steady_sweeps.FIG6_WORST_PCT_NG`. Both sit above the measured values (README, Figure 6 row), so they catch a regression rather than certify agreement |
| Fig. 9/10 settled state vs the figure | +/-1 %, both figures | `test_fuel_step.FINAL_TOL_PCT`. Compared over the reference's own last samples, up to where Ballin's record stops at the panel frame, not to the end of our run: on Figure 10 nothing has settled by then (#68) |
| Whole-curve RMS, normalised by each panel's excursion | **4 % ratchet** | `test_whole_curve.WHOLE_CURVE_RMS_CEILING_PCT`. Not a tolerance: set just above the worst measured panel, so a regression fails and an improvement is free. Today, over ten panels with the riser excluded, **0.62<!--rj:whole_curve_summary.min-->–3.80 %<!--rj:whole_curve_summary.max-->**, mean 1.64 %<!--rj:whole_curve_summary.mean--> |
| Step-edge jump across the riser hole | +/-12 % of panel excursion | `test_whole_curve.STEP_JUMP_TOL_PCT`. A line tracer cannot follow a vertical edge, so each reference trace has a sampling hole at the fuel step; the chord `np.interp` would draw across it is excluded from the whole-curve metric (#66), and the edge's size is asserted here instead, since the hole's two endpoints are Ballin's |
| Table B.1's printed station state | +/-0.5 % | `test_trim_points.STATE_TOL_PCT`. Printed *numbers*, with no read error of ours in them |
| Figures 9/10 read error | +/-0.6 % of full scale | `test_figure_consistency.READ_ERROR_CEILING_PCT_FS`, on the fuel panels whose levels the captions print, and on the two figures' shared 400 lbm/hr trim (see "The figures' own read error" below) |
| HMU collective at a Table B.1 trim | 5-95 % of maximum | `test_hmu_trim.COLLECTIVE_RANGE_PCT`. The report prints no collective for these trims, so this is a believability band, not a comparison. It tests that Appendix C and Table B.1, digitized independently, agree at all |
| **Closed loop vs Table B.1** | **+/-0.5 %** | `test_closed_loop.CLOSED_LOOP_TOL_PCT`, the Phase 5 gate. Measured worst **0.246 %<!--rj:closed_loop.worst-->**, descent fuel flow, as a **mean over the limit cycle**: the loop settles to a cycle, not a point (#61), so a terminal sample depends on where the run stops. Fuel flow is an output, so the open-loop descent shaft-power error shows up in fuel flow instead |
| **Table 1 eigenvalues** | **+/-6 %** | `test_trim_points.EIGENVALUE_TOL_PCT`. **A characterization, not an accuracy claim**: a regression detector set a quarter above today's **4.70 %<!--rj:eigenvalues.notes.table_1_stats.worst-->** worst over all 27<!--rj:eigenvalues.notes.table_1_stats.n--> printed modes (rms 1.94 %<!--rj:eigenvalues.notes.table_1_stats.rms-->), each paired by dominant state and asserted by `test_all_dof_models.test_every_printed_table_1_mode` at Table B.1's 20895 rpm, with the report's ±2 % extraction step [pdf p.28]. Ballin calls a 4 % mismatch "good agreement" for his own comparison [pdf p.29]; ours is not that, and this row keeps the difference visible |
| Fuel overestimate, 81-86 %NG, vs the report's own claim | +/-8 % | `test_steady_sweeps.FUEL_OVERESTIMATE_CEILING_PCT`. The report's sharpest printed steady-state claim about itself: the real-time model tends to "overestimate fuel consumption by as much as five percent" [pdf p.39]. A ceiling above our measured 6.34 %<!--rj:fuel_overestimate.ours_peak_pct-->, against Ballin's own line's 5.98 %<!--rj:fuel_overestimate.ballin_peak_pct-->: ours exceeds his claim, and that is the finding |
| Table 1 NG mode, each model against its own column | +/-2 % (full), +/-6 % (quasi-steady), +/-3 % (reduced) | `test_all_dof_models.NG_MODE_TOL_PCT`. A characterization set a quarter above today's measurement: each of our three linear models, made the way Ballin made his [pdf pp.27-28], against its own column (#49) |
| What the 7 ms frame costs the NG mode | under 16 % (hover), 5 % (level), 2 % (descent) | `test_discrete_map.FRAME_COST_TOL_PCT`. Our real-time map against our continuous model, both with the pressure iterations converged; set a quarter above today's measurement. A characterization of our frame, not a comparison with Table 1 (see "What the real-time approximation costs") |

### Appendix B, element by element

There is no single per-element tolerance. The bounds below live in
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
per-element comparison.** Two are excluded rather than bounded, each with a test recording
why: the 6-DOF `A(P45,P3)`, a near-cancellation printed to four figures, and the `A(NP,NP)`
diagonals, which are Gen Hel's.

### The figures' own read error

In each of Figures 9 and 10 the fuel-flow panel is the *input*, and the caption prints both
its levels, so digitizing it measures the figure's read error directly. The error is a
**pixel offset**, a fixed fraction of each panel's span rather than of the value, so it reads
very differently on a wide axis and a narrow one. The two figures also start from the same
400 lbm/hr trim on axes of different span, so their pre-step plateaus must agree within the
read error on every panel. Both checks are in `validation/test_figure_consistency.py`.

### Comparing a transient to a steady-state figure

**Never compare a transient trace against a steady-state locus at matched speed.** During a
chop the choked station 4.1 nozzle passes the same flow at a lower pressure once the
temperature falls, so Ps3 sits *below* its equilibrium value at that NG by construction, and
Ballin's own Figure 10 runs below his own Figure 8. The effect is quasi-steady, not a
pressure lag. The valid comparisons are trace against trace, settled state against the
sweep, and the **phase plane** (Ps3 against NG), which discards the time axis and with it
the unprinted step time.

## Validation surfaces available, in order of what they are worth

**Evidence that is printed as numbers and localizes a disagreement outranks evidence that is
digitized or aggregated**: a printed number carries none of our read error, and an
element-wise comparison says *which* partial derivative moved where a mode only says that
something did.

| | Surface | pdf p. | Scale | our read error | what only this surface constrains |
|---|---|---|---|---|---|
| 1 | **Table B.1** | 67 | 21 numbers | none | the equilibrium every other comparison is taken about. If this is wrong, every derivative below is evaluated at the wrong point |
| 2 | **Appendix B** | 68–76 | 297 elements, 206 non-zero | none | every partial derivative, **localizably**. The three error families in `test_appendix_b_elements.py` — the missing load inertia, the NP diagonal residual, the `f7` group — were each found here and none could have been found from a mode |
| 3 | **Figures 9–10** | 45–46 | 6 traced panels each | ±0.6 % of full scale | the integrated time response: volume dynamics and heat-sink lags over a whole trajectory, which no printed table touches |
| 4 | **Closed loop vs Table B.1** | — | 3 trims | none | the fuel control. No printed surface covers it; this is internal consistency, not agreement with the report |
| 5 | **Figures 6–8** | 40–42 | ~70 markers | ours | steady behaviour across the **whole operating range** rather than three points, so a map error shows as a trend in a region instead of one unattributable number |
| 6 | **Table 1 eigenvalues** | 31 | 27 modes, all compared | none | *nothing the rows above do not already carry.* The modes are an aggregation of Appendix B's matrices at Table B.1's trims, and they cover only the 2-DOF, 5-DOF and reduced-5 models — the 3-DOF and 6-DOF print no eigenvalues at all |
| — | Figures 11–15 | 48–53 | — | — | **not reproducible**: they need Gen Hel |
| — | Table 3 | 44 | 6 cases | — | **not a target**; see below |

The ordering is **not** a ranking of how hard each comparison is to pass — the loosest bound
in the project sits on row 6 and one of the tightest on row 1 — and **not** an instruction to
stop running the lower rows: a mode that moves when no element did would be a finding about
our own linearization. Where two rows disagree, the higher one wins and the cost to the
lower one is recorded rather than tuned away (as in #73).

### Two things not to do

- **Table 3 (pdf p.44) is not a validation target.** It looks like the best steady-state
  benchmark in the report, six hardware test cases, but its model column was produced with
  compressor and turbine functions from the *NASA-Lewis prototype engine*, used in place of
  the standard ones and printed nowhere [pdf p.39]. Use it as evidence of method, never as
  pass/fail.
- **Do not fix the T41/T45 trim bias.** Ballin's model underestimates turbine temperature at
  trim, "a characteristic of the real-time model which was found for all validation cases"
  [pdf pp.39, 47]. A replication reproduces that bias with the same sign; one that matches
  the hardware better than his did is a worse replication.

## Admitted limits that become ours

Ballin's own [pdf pp.39, 47, 50]: the 100 %NG ceiling; invalid at low NP; T41/T45 trim
underestimated throughout; direction-independent heat-sink constants; no T45 heat sink; no
inlet guide vane dynamics; damping slightly high; load-demand compensation weaker than the
real aircraft. And one positive result: modelling the high-speed dynamics of mass flow
between the control volumes "was found to be unnecessary" [pdf p.54], which is what
justifies the 5-state formulation.

## Open questions

Tracked in `docs/notes/open-questions.md`, whose header carries the tally that
`tests/test_open_questions_ledger.py` checks. Numbering is deliberately non-contiguous; never
renumber a row. Each row states what was asked and what was decided, and that test fails if
a closed row states a verdict without a finding. Anything the report does not answer goes
there rather than into the code as a guess.
