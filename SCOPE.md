# Scope — T700 Replication

> **This is a proposal, not an agreement.** You chose "full model + controls" and
> "digitize figures + tables"; those are settled. Everything else here — the seven
> phases, the gates, and the entire tolerance table — is my draft and was agreed by
> nobody. Cut, reorder, or delete any of it.

Status: **first revision**, written after a structural pass over the report. Facts below
carry PDF page citations; anything without one is still a hypothesis. See
`docs/notes/index.md` for the section map and `docs/notes/open-questions.md` for gaps.

## Objective

Reproduce, in Python, the engine and control-system model of NASA TM-100991, and
demonstrate that it matches the results published in that report to a stated tolerance.

"Done" means: the report's result figures regenerate from this code, the deviations are
quantified, and every deviation is either within tolerance or explained.

## What the report actually contains

This shaped the scope, so it is recorded here rather than left implicit:

- A **5-state nonlinear model**: NG, NP, P3, P41, P45 — two spool speeds plus three
  control-volume pressures, printed as a vector at [pdf p.27]. The fuel control adds
  **22 more states plus 4 memory elements**, so the assembled model is ~31 states.
- The Appendix B **2- and 3-DOF models are separately extracted, not order-reduced** from
  the larger ones [pdf p.27, p.33] — so Phase 6 has two distinct tests, not one.
- A **fixed-step real-time formulation** with a maximum time step of 10 ms, set by a
  0.1% allowable inter-step error, and with response that is time-step dependent because
  of an *opened* compressor mass-flow iteration [pdf p.38].
- **Appendix A** — engine model constants (Table A.1) and function tables, plus plots of
  mass-flow and energy functions [pdf pp.55–66].
- **Appendix B** — small-perturbation linear models: three trim conditions x four DOF
  variants (2, 3, 5, 6), 12 figures, 297 matrix elements [pdf pp.67-76].
- **Appendix C** — the complete T700 fuel control system model [pdf pp.77–101].
- The helicopter load enters as Qreq = Qmr + Qtr + Qacc + Qdamp [pdf p.28, Eq. 56], with
  rotor torques supplied by the **Gen Hel UH-60A blade-element simulation** — a separate
  program (ref. 2), not modelled in this report.

## In scope

The T700-GE-700 as modelled in the report:

- Gas path, component by component: inlet, compressor, combustor, gas generator turbine,
  free power turbine, exhaust.
- The three control volumes and their pressure states (P3, P41, P45).
- Rotor dynamics: gas generator spool (NG) and power turbine spool (NP).
- The fuel control system of Appendix C — hydromechanical and electrical units, governing,
  scheduling, limiting.
- The report's fixed-step integration scheme, including the opened compressor mass-flow
  iteration and its time-step sensitivity.
- The load interface: Qreq as an input, plus the load-torque derivatives ∂Qreq/∂NP and
  ∂Qreq/∂ṄP that the report uses for its linear models [pdf p.28, Eq. 55].

## Out of scope

- **The Gen Hel UH-60A rotor and airframe model.** The report does not contain it; it
  consumes it. Qreq is an input boundary condition to us, exactly as the rotor model is an
  external supplier to the report. This is the report's own boundary, not one we imposed.
- Engine start, windmill, and shutdown regimes, unless Appendix C turns out to model them.
- Other T700 variants (-701, -701C, -401) and later FADEC control units.
- Hardware I/O, real-time OS integration, pilot-in-the-loop display.

## Deferred (designed for, not built now)

- **C++ lookup-table thermodynamics backend.** `t700.thermo` is an interface from day one
  so this drops in behind it without touching component code. Build it when profiling
  shows property evaluation dominating frame time — not before.
- Numba or any other JIT.
- Wall-clock-paced real-time execution harness.

## Phases

Each phase has a gate. Do not start the next until the gate passes or the failure is
explicitly accepted.

**Phase 0 — Ingest the report.** **COMPLETE.**
Per-page text extracted (done). Remaining: verify the nomenclature against rendered page
images into `docs/notes/symbols.md`; inventory every figure and table in
`docs/notes/inventory-appendix-{a,b,c}.md` with whether it is transcribable or must be
digitized. (This line said `data-inventory.md` until 2026-09-12; no such file was ever
written -- the inventories were split per appendix instead.)
*Gate:* we know precisely what data the report gives and what it withholds.

**Phase 1 — Capture the data.** **COMPLETE.** All eleven engine maps and **all eight**
Appendix C scheduling functions are captured; Fig. C30 (`F_HM7`), the last and hardest,
landed 2026-09-13 (open question #50).
Transcribe Appendix A constants and function tables; digitize the mass-flow and energy
function plots and any component maps into `data/` as CSV with provenance headers.
Digitize the report's result figures into `data/reference/` as validation traces, with a
recorded read-error estimate.
*Gate:* every function the model needs exists as data, or is logged as unavailable.

**Phase 2 — Units and thermodynamics.** **COMPLETE.**
`t700.units`; `t700.thermo` interface plus NumPy backend, matching the report's gas
property treatment — the constants named `KH2`, `KH3a`, `KH32`, `KH411`, `KH412`, `KH45`
in Table A.1 suggest linear enthalpy fits rather than full property tables [pdf p.55].
*Gate:* property values reproduce any worked example in the report.

**Phase 3 — Steady-state engine.** **COMPLETE and validated.**
Components and control volumes assembled; equilibrium solved across the operating range
with fuel flow prescribed.
*Gate:* steady-state gas path variables match the report's trim conditions within
tolerance.

**Phase 4 — Dynamics.** **COMPLETE**, including the station 4.1 heat-sink model as a
switch (Eqs. 48–53), so the engine runs in either configuration Ballin published. Shape
now compared, not just endpoints: the Figure 9 T41 overshoot is **0.87× his**, having gone
3.41× -> 1.67× when the heat sink landed as Eq. 50's lead-lag and 1.67× -> 0.87× when it was
rebuilt on the printed Eqs. 48-49. (This line reported 1.67× as the present state until
2026-09-12.)
Five-state integration at the report's frame rate, including the opened compressor
mass-flow iteration. Reproduce the report's time-step sensitivity (0.1% at 10 ms) as a
test — matching the *error behaviour*, not just the answer.
*Gate:* open-loop response has the report's shape and time constants.

**Phase 5 — Control system.** **IN PROGRESS. The HMU is built; the ECU is not.**
`src/t700/control/hmu.py` implements Figures C9-C22 -- the droop line, the torque motor,
the Ps3 and NG sensors, the collective rigging, and the idle / acceleration / deceleration
cam cascade -- against `src/t700/control/schedules.py`, which loads all eight digitized
scheduling functions. It runs open-loop with the ECU trim held at its null. What remains is
the **ECU** (Figs. C1-C8) and then closing the loop around the engine.

Its data is complete.
`src/t700/control/constants.py` holds **73** constants (57 from Table C.1 plus 16 read off
the figures; the other 3 of the 76 read from the report live in `t700.constants` and
`t700.units`). **All 8 scheduling functions are digitized** and committed under
`data/schedules/` -- F_EC1 and F_HM1 through F_HM7. **What remains is the 22 block
diagrams**, which are the whole of the phase.

*Gate warning, worth deciding before starting rather than after:* the report's only
closed-loop figures are 11-15, and those need the Gen Hel UH-60A simulation, which is out
of scope by construction. So **Phase 5 will have no figure to validate against.** Its
checks have to be structural -- does each block diagram reproduce as drawn -- plus the
report's one printed closed-loop number, the 0.2 % rotor-speed agreement on Fig. 12
[pdf p.50].
Appendix C, implemented and closed around the engine.
*Gate:* closed-loop transients match the report's figures.

**Phase 6 — Validation.** **SUBSTANTIALLY COMPLETE**, and the remaining gaps are in the
source rather than the work. Done: Table B.1's full printed state, 21 numbers at rms
0.20 % (worst: shp −0.72 % at the descent trim, NG +0.13 % at level); Table 1
eigenvalues; **all 297 Appendix B elements loaded and structurally checked across the
four DOF variants Appendix B prints**, of which roughly 150 carry a numeric per-element
comparison -- the 2-DOF and 5-DOF models fully, the 3-DOF and 6-DOF through their T41
row/column, `d`, and the structural zeros. The 3-DOF and 6-DOF `b` vectors (27 elements)
are not compared to anything. `DOF.REDUCED_FIVE` has **no** Appendix B figure, so "all
five linear models" was wrong: the printed elements come from four; the steady-state sweeps of Figures 6–8 across their full range; Figures 9–10
transients in the configuration they were actually generated in.

Not done, with reasons: Figures 11–15 need Gen Hel and are out of scope by construction;
and Phase 5 below. (The GE reference series on Figures 6-8 was listed here as
"digitized but uncompared" until 2026-09-12; `validation/test_ge_reference.py` compares
all three series and is explicit that agreeing with GE is not evidence of anything.)
Two independent lines:
1. *Figure reproduction* — regenerate the report's result figures; `/validate` reports
   every deviation.
2. *Jacobian comparison* — extract stability derivatives from our nonlinear model using
   the report's own central-difference method (±2% perturbation, integrations suppressed
   after trim [pdf p.28]) and compare element by element against the published A and B
   matrices of Appendix B, at all three trims, across the 2/3/5/6-DOF variants.

   Two limits found during extraction, both real:
   - **The NP row is not independently reproducible.** Row 2 and b(2) are built from
     dQreq/dNP and dQreq/dNPdot (Eqs. 57-62), which come from the external Gen Hel
     UH-60A simulation. Exclude that row, or treat those derivatives as an input.
   - **Compare the 6-DOF matrices element-wise only, never by eigenvalue.** As printed,
     6-DOF trim 1 is open-loop unstable (+3.07 /sec) through precision loss, not a typo:
     row 6 differences ~2e5 terms to give ~1e3, so +/-50 of rounding moves the slow mode
     from -0.99 to +7.3. Do not chase that instability.

The second line is the stronger test and it is why Appendix B matters: it validates the
entire linearized model, term by term, against printed numbers — no figure digitization,
no read error, no eyeballing. Any structural error in a component equation shows up as a
wrong matrix element that points straight at the offending partial derivative.

## Which source wins when the report disagrees with itself

Established 2026-09-11, and it governs every comparison from here:

**Table B.1 is primary.** It is printed *numbers*, so it carries no digitizing error of
ours at all; it underpins Appendix B's linear models, whose eigenvalues our own Jacobian
reproduces; and Table 1's modes come from the same trims. It is the most internally
corroborated data in the report.

**Figures 6-8 are secondary** -- evidence of *trend across the operating range*, not of
absolute value. They are plots we digitized, so they carry our read error on top of
whatever they are.

They genuinely disagree, and systematically: Table B.1 against Figure 7 differs by
**+5.49 %, +1.32 %, -0.14 %** on shaft power at Wf = 267.7, 349.3, 476.3 lbm/hr -- growing
as power falls. See open question #46 for what was ruled out.

The practical consequence: **do not tune the model below the report's own internal
spread.** At the descent condition that spread is 5.5 % on power, and we sit 0.73 % from
Table B.1. Chasing further is fitting to one of two sources that contradict each other.

## Tolerances

**Rewritten from the report, 2026-09-10.** The guessed table that stood here is gone. The
body was read (pdf pp.38-54) and it turns out **Ballin states no percentage tolerance for
pressures, temperatures, torque, or fuel flow anywhere.** Six numbers exist in the whole
report, and most of them are deviations he *accepted*, not accuracy he *met*:

| Claim | Value | Scope | pdf p. |
|---|---|---|---|
| Max inter-step error | 0.1% | sets the 10 ms step cap | 38 |
| Fuel consumption overestimate | =<5% | only 81-86% NG | 39 |
| NG overestimate, hardware case | 4% | lowest-power case only | 39 |
| NG overestimate, open-loop step | 1-2% | Figs. 9-10 | 39 |
| Rotor speed agreement | 0.2% | Fig. 12 only | 50 |
| Eigenvalue mismatch called "good" | 4% | Table 1 comparison | 29 |
| NG validity ceiling | 100% | above this the model is not claimed | 38 |

Everything else Ballin says is qualitative -- "fair agreement", "agree very well". So our
tolerances cannot be derived from the report; they have to be **declared by us and labelled
as ours**. Provisionally, and only for quantities the report shows in a reproducible figure:

| Quantity | Steady state | Transient | Basis |
|---|---|---|---|
| NG | +/-2% | +/-2% | Ballin's own accepted 1-4% overestimates |
| Gas path P, T | +/-3% | +/-3% | ours, undeclared by the report |
| Fuel flow | +/-5% | +/-5% | his 5% band, widened to all NG |
| Torque / SHP | +/-3% | +/-3% | ours, undeclared by the report |
| Appendix B matrix elements | +/-5% per element | -- | his own 4% "good agreement" |

**Three of those page citations were wrong until 2026-09-12** and pointed at the figures
being discussed rather than at the sentence: the 4 % and the 1-2 % were cited to pp.44 and
45-46, and the 0.2 % to p.47. All three of the first two live in one paragraph of **p.39**,
read off the raster; the 0.2 % is on **p.50**. `docs/notes/body-validation.md` had them
right, which is the argument for citing the notes rather than re-deriving a page number.

Transient tolerances are floored by the digitizing error of the source figure, which
Phase 1 records per figure.

### The figures' own read error, measured rather than estimated -- and in percent of FULL SCALE

**Added 2026-09-12, corrected the same day.** Figures 9 and 10 contain one panel whose true
value is printed: `WFPH` is the *input*, and both of its levels are stated in the caption
("from 400 to 775" and "from 400 to 125 lbm per hour"). Digitizing a panel whose answer is
known measures the read error of those figures directly:

| Figure | level | printed | digitized | % of value | **% of full scale** |
|---|---|---|---|---|---|
| 9 | pre-step | 400.0 | 407.30 | +1.83 % | **+0.97 %** |
| 9 | post-step | 775.0 | 777.04 | +0.26 % | +0.27 % |
| 10 | pre-step | 400.0 | 406.69 | +1.67 % | **+1.34 %** |
| 10 | post-step | 125.0 | 125.48 | +0.38 % | +0.10 % |

**The right-hand column is the one to use.** `calibrate()` in `tools/digitize_fig910.py`
maps pixel rows to values by a straight line between the two frame rows, so a read error
is a *pixel offset* -- a fixed fraction of the panel's span, not of whatever value it lands
on. The 1.83 % is large only because 400 sits low on a 250-1000 axis. Converted through
each panel's own span, the floor is:

| panel | Fig. 9 span | floor | Fig. 10 span | floor |
|---|---|---|---|---|
| WFPH | 750 | +1.83 % | 500 | +1.67 % |
| PS3 | 200 | +1.24 % | 200 | +1.70 % |
| **PCNG** | 20 | **+0.21 %** | 40 | **+0.59 %** |
| T41 | 1000 | +0.44 % | 1000 | +0.61 % |
| T45 | 1000 | +0.62 % | 1000 | +0.85 % |
| TORQ45 | 400 | +2.22 % | 400 | +2.91 % |

This file briefly said "any deviation below about 1.8 % against Figures 9-10 is at or under
the read error". That is **wrong, and wrong in the permissive direction on the channel that
matters most**: on PCNG the floor is 0.21 %, nine times tighter, and PCNG is what most of
the project's headline transient agreements are quoted on. On TORQ45 it is too tight.

The evidence for the pixel-offset model is a prediction, not an argument. Figures 9 and 10
read the *same* 400 lbm/hr trim on axes of different span, so their disagreement is pure
read error: calibrating only on the two WFPH panels predicts **+0.340 %NG** on PCNG against
a measured **+0.340**, and +0.73 psia on PS3 against a measured +0.59. A percent-of-value
model cannot produce that. Pinned by
`validation/test_figure_consistency.py::test_the_read_error_is_a_pixel_offset_not_a_fraction_of_value`.

Two limits of the measurement, both real. It is a *pre-step* offset: Figure 10's WFPH
post-step level, whose truth is exactly 125, drifts from +3.7 % at t = 0.9 s to -7.4 % at
t = 4.4 s -- a panel skew of about -0.8 %FS/s that the digitizer does not deskew, so the
mean lands at +0.38 % only by cancellation. And the model predicts PCNG, PS3 and T41 but
not T45 or TORQ45, the latter being independently known-bad reference data.

### Comparing a transient to a steady-state figure

**A transient trace must not be compared against a steady-state locus at matched speed.**
Figures 6-8 are equilibrium sweeps; Figures 9-10 are transients. During a chop, fuel is
cut, T41 falls, and the choked station 4.1 nozzle then passes the same flow at a lower
P41 -- so Ps3 sits *below* its equilibrium value at that NG, by construction. Over its full
record Ballin's own Figure 10 line runs **-1.68 % to -8.73 %** against his own Figure 8, and
his Figure 9 accel runs **-1.43 % to +7.35 %**. Those are the physics, not a defect in
either figure. (Ranges corrected 2026-09-12: "2.3-7.0 %" and "up to 3.9 %" were values at
four sampled instants, quoted as if they were the extremes.)

The effect is **quasi-steady, not a pressure lag** -- the pressures are fast and the spool
is slow. Pinning NG on the steady locus and cutting fuel to 125 lbm/hr, the algebraic
pressure solution alone gives -5.9 / -9.3 / -12.0 / -15.5 % at 78 / 82 / 86 / 90 %NG, and
the choked-nozzle scaling sqrt(theta41) predicts -6.3 / -10.2 / -13.8 / -17.1 %. No other
path (T3 = T2*f2(PR), the bleed schedules) comes close, and the sign cannot reverse.

The valid comparisons are (a) trace against trace, (b) settled state against the sweep, and
(c) **the phase plane** -- Ps3 against NG, which discards the time axis and with it the
unprinted step time (open question #37). Tolerances for the last two, ours:

| Comparison | Tolerance | Basis |
|---|---|---|
| Phase-plane state trajectory, Ps3 vs NG | +/-2.5 % | measured 0.59 % worst on Fig. 9, 1.07 % on Fig. 10 over the 78-90 %NG grid the test asserts; 2.15 % at 76 %NG, the last point of his record, which the test does not check |
| Settled state vs the Fig. 6-8 sweep | +/-3 % | the gas-path P,T row above |
| Figures 9 and 10 at their shared 400 lbm/hr trim | +/-1.5 % | five of six panels agree to 0.54 % worst (t45); the sixth, torq45, disagrees by 4.84 % and is excluded as bad reference data |
| Fig. 9/10 initial trim vs the figure | +/-1.5 % | `test_fuel_step.INITIAL_TOL_PCT`; this is the same physics the steady tests already check |
| Fig. 9/10 settled state vs the figure | +/-5 % (fig 9), +/-8 % (fig 10) | `test_fuel_step.FINAL_TOL_PCT`; tightened from 20 % when the printed convergence criterion replaced an invented tolerance, never widened |
| Whole-curve RMS, normalised by each panel's excursion | **8 % ratchet** | `test_whole_curve.WHOLE_CURVE_RMS_CEILING_PCT`. Not a tolerance -- set just above the worst measured panel so a regression fails and an improvement is free. Lower it when the model improves; never raise it. Current: mean 3.79 %, range 1.96-7.62 % over nine panels |
| Table B.1's printed station state | +/-0.5 % | `test_trim_points.STATE_TOL_PCT`. Much tighter than the +/-3 % gas-path row because these are printed *numbers*, with no read error of ours in them at all |
| Figures 9/10 read error | +/-1.6 % of full scale | `test_figure_consistency.READ_ERROR_CEILING_PCT_FS`, measured on the WFPH panel -- see the section above, and note the currency |
| HMU collective at a Table B.1 trim | 5-95 % of maximum | `test_hmu_trim.COLLECTIVE_RANGE_PCT`. The report prints no collective for these trims, so this is a believability band, not a comparison. What it tests is that Appendix C and Table B.1 -- digitized independently -- agree at all |

**Every tolerance any test asserts is in this file.** Five were inline in test files until
2026-09-12, two of them carrying their own docstring note saying they belonged here. A
tolerance that lives only beside the assertion it governs is one nobody compares against
its siblings, which is how two tests come to bound the same quantity differently.

### Two things not to do

- **Table 3 (pdf p.44) is not a validation target.** It looks like the best steady-state
  benchmark in the report -- six hardware test cases with NG, NP, WA2, P3, T3, T45 -- but
  its model column was produced with compressor-mass-flow and turbine-energy functions
  derived from the *NASA-Lewis prototype engine*, supplied by NASA Lewis "in place of the
  standard functions" [pdf p.39]. Those functions appear nowhere in the report. Its own
  caption says "This does not represent specification T700 performance." Use it as
  evidence of method, never as pass/fail.
- **Do not fix the T41/T45 trim bias.** Ballin's model underestimates turbine temperature
  at trim in *every* validation case and he says so [pdf p.54]. Our replication should
  reproduce that bias with the same sign. A version that matches the hardware better than
  Ballin's did is a worse replication.

## Validation figures actually available

Five of the report's ten result figures can become automated reference traces; the other
five cannot, because they run the engine inside the Gen Hel UH-60A blade-element
helicopter simulation, which this report does not contain.

| Figure | pdf p. | What | Points |
|---|---|---|---|
| 6 | 40 | NG vs Wf, 3 models | ~30 |
| 7 | 41 | SHP vs Wf, 3 models | ~30 |
| 8 | 42 | Ps3 vs NG, 3 models | ~30 |
| 9 | 45 | fuel step up, 400 -> 775 lbm/hr, **6** panels | ~4,100 |
| 10 | 46 | fuel step down, 400 -> 125 lbm/hr, **6** panels | ~3,870 |
| 11-15 | 48-53 | closed-loop UH-60A vs flight test | not reproducible |

**Corrected 2026-09-12.** This table was a pre-Phase-1 estimate and had never been
updated against what was actually captured. Two things it got wrong beyond the counts:
Figures 9 and 10 have **six** panels each, not five (the sixth, `WFPH`, is the input --
and it is the panel that later gave the only *measured* read error in the project, above);
and they are **not** marker-based. Figures 6-8 and the `+` reference series of 9-10 are
marker series that digitize by point extraction, but Ballin's own model output on 9-10 is
a **continuous line, curve-traced** at 600-740 samples per panel, which is where the
whole-curve shape comparison comes from. Actual total: **~8,100** reference points, not
~550. See `docs/notes/digitizing-recipe.md`.

## Admitted limits that become ours

From the Conclusions [pdf p.54] and the validation section. These are Ballin's stated
limitations, and a replication inherits them rather than repairing them: the 100% NG
ceiling; invalid at low NP; T41/T45 trim underestimated throughout; direction-independent
heat-sink constants; no T45 heat sink; no inlet guide vane dynamics; damping slightly
high; load-demand compensation weaker than the real aircraft.

One positive result worth keeping in view: omitting high-speed inter-volume mass-flow
dynamics "was found to be unnecessary", which is what justifies the 5-state formulation
in the first place.

## Open questions

Tracked in `docs/notes/open-questions.md` -- 49 rows, numbering deliberately
non-contiguous (see the note at the top of that file). Anything the report does not answer
goes there rather than into the code as a guess.
