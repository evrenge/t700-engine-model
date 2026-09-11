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
`docs/notes/data-inventory.md` with whether it is transcribable or must be digitized.
*Gate:* we know precisely what data the report gives and what it withholds.

**Phase 1 — Capture the data.** **COMPLETE for the engine; Appendix C's 8 control schedules remain.**
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

**Phase 4 — Dynamics.** **COMPLETE; transition *shape* not yet compared, only endpoints.**
Five-state integration at the report's frame rate, including the opened compressor
mass-flow iteration. Reproduce the report's time-step sensitivity (0.1% at 10 ms) as a
test — matching the *error behaviour*, not just the answer.
*Gate:* open-loop response has the report's shape and time constants.

**Phase 5 — Control system.** **NOT STARTED.** 76 constants transcribed; 8 schedules and 22 block diagrams remain.
Appendix C, implemented and closed around the engine.
*Gate:* closed-loop transients match the report's figures.

**Phase 6 — Validation.** **PARTIAL.** Trim, transient endpoints and 5 eigenvalues done; the 297-element Appendix B comparison remains.
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
| Fuel consumption overestimate | =<5% | only 81-86% NG | 38-39 |
| NG overestimate, hardware case | 4% | lowest-power case only | 44 |
| NG overestimate, open-loop step | 1-2% | Figs. 9-10 | 45-46 |
| Rotor speed agreement | 0.2% | Fig. 12 only | 47 |
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

Transient tolerances are floored by the digitizing error of the source figure, which
Phase 1 records per figure.

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
| 9 | 45 | fuel step up, 400 -> 775 lbm/hr, 5 panels | ~230 |
| 10 | 46 | fuel step down, 400 -> 125 lbm/hr, 5 panels | ~230 |
| 11-15 | 48-53 | closed-loop UH-60A vs flight test | not reproducible |

All five are marker-based, so they digitize by point extraction rather than curve tracing.
Budget ~550 reference points. See `docs/notes/digitizing-recipe.md`.

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

Tracked in `docs/notes/open-questions.md` -- 37 rows, numbering deliberately
non-contiguous (see the note at the top of that file). Anything the report does not answer
goes there rather than into the code as a guess.
