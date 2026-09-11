# Report index — NASA TM-100991 (Ballin, 1988)

File: `docs/ballin-tm100991.pdf`, 104 PDF pages.

## Citation convention

**Cite PDF page numbers**, written `[TM-100991 pdf p.NN]`. Printed page numbers *do*
exist on the page images (the text layer just loses them): **printed = pdf − 14** across
the appendices — pdf p.56 is printed p.42. Cite the PDF page; mention the printed one only
when quoting a caption that names it. The scan carries no usable
printed page numbers in its text layer, so the PDF page is the only unambiguous handle.
The cover states "95 p", so printed and PDF pages differ by front matter — do not mix the
two.

## Section map (PDF pages, from a text-layer scan — refine as sections are read)

| PDF pages | Content |
|---|---|
| 1–6 | Cover, report documentation, **errata sheet** |
| 7–14 | Draft title page, nomenclature |
| 9-13 | Nomenclature, one alphabetical list, 106 entries -- fully transcribed to `symbols.md` |
| 15 | Summary, Introduction |
| 17 | Engine and fuel control system — physical description |
| 19–25 | Engine model — component equations |
| 26–37 | Real-time implementation |
| 38–44 | Fuel control system model; validation |
| 45-53 | Real-time model results (open-loop steps, closed-loop UH-60A flight comparisons) |
| 54 | Conclusions |
| 55–66 | **Appendix A — T700 engine model constants and function tables** |
| 67-76 | **Appendix B — Linear models**: 3 trim conditions x **4 DOF variants (2/3/5/6)**, 12 figures |
| 77–101 | **Appendix C — T700 fuel control system model** |
| 102 | References |

## Errata (PDF p.6)

Official errata sheet from Ames (Shelley J. Scarich, M/S 241-13):

> On p. 2, figure 1, a correction should be noted. In the upper middle portion of the
> figure, GAS INDUCED FLAP REGRESSIVE should be CONTROL-SYSTEM-INDUCED FLAP REGRESSIVE.

Figure 1 only — no numerical erratum. Recorded so nobody rediscovers it.

## Facts established so far (each with its page)

- **Appendix C prints no numbered equations at all** [pdf pp.77-101]. It is nomenclature,
  one constants table (Table C.1, 57 rows, pp.83-84) and **30 figures** — 22 block
  diagrams and 8 plots. The model is specified *graphically*. A further **19 constants
  appear only inside block diagrams** and are absent from Table C.1.
- **Total model state count is ~31, not 5.** The engine contributes 5 (below); the fuel
  control adds **22 continuous states + 4 memory elements** [Appendix C] — 19 first-order
  lags/lead-lags, 3 limited integrators, 3 hysteresis memories and one 0.015 s transport
  delay. One lag (TAU45, the thermocouple harness) has a **time constant that is itself a
  lookup on the state**.
- **Engine states, 5-DOF nonlinear model**: NG, NP, P3, P41, P45 — two rotor speeds plus
  three control-volume pressures. **Confirmed** — the vector is printed explicitly at [pdf p.27, below Eq. 54] and the
  modes are listed in Table 1 [pdf p.31].
- **A 6th state, T41, exists** in the station-4.1 heat-sink variant; the 6-DOF state vector
  {NG, NP, P3, P41, P45, T41} is printed explicitly [pdf p.32, below Eq. 65] and is the
  only linear-model state vector written out anywhere in the report. The 3-DOF variant is
  {NG, NP, T41} — **inferred, not printed** (from the structure of G2 and d̄, pdf p.33).
- **The Appendix B 2- and 3-DOF models are NOT order-reduced** from the 5- and 6-DOF
  ones. They are separately *extracted* from a quasi-steady nonlinear simulation
  [pdf p.27, p.33]. Only the third column of Table 1 ("Reduced-Order 5-DOF") is an order
  reduction, and that one is not in Appendix B. Two different tests, both available.
- **Control input is a single scalar**: Wf, fuel flow, lbm/sec [nomenclature pdf p.13].
  Table B.1 quotes trim fuel flow in lbm/**hr** — different unit, same quantity.
- **Appendix B's three trims** [Table B.1, pdf p.67], all sea level, UH-60A at 16825 lbm,
  NP = 20895 rpm: hover (911.1 hp/engine, NG 41638, Wf 476.3 lbm/hr); 80 kt level
  (552.6 hp, NG 39768, Wf 349.3); 80 kt / 1000 ft-min descent, γ = -7.08 deg (302.6 hp,
  NG 38072, Wf 267.7).
- **The NP row of any linear model is not independently reproducible.** Row 2 and b(2) are
  not plain Jacobian entries — Eqs. 57-62 build them from ∂Qreq/∂NP and ∂Qreq/∂ṄP, which
  come from the external Gen Hel UH-60A simulation. Validation of that row needs those
  load derivatives, which this report does not supply.
- **Time step limit**: max 10 ms, set by a 0.1% maximum allowable error between time
  steps. Response is time-step dependent because of an *opened* (broken-open across
  frames) compressor mass-flow iteration [pdf p.38].
- **Load torque** Qreq = Qmr + Qtr + Qacc + Qdamp — main rotor, tail rotor, accessories,
  gearbox damping [pdf p.28, Eq. 56]. (Corrected: an earlier draft of this file listed
  `Qgb` and `Qgg`, which are OCR inventions. The report prints **Q_eng**, torque to the
  helicopter gearbox from both engines, and **Q_GT**, gas generator torque output —
  see `symbols.md`.) The rotor torques come from the **Gen Hel UH-60A
  blade-element simulation** (an external program, NASA CR ref. 2), *not* from this
  report. Qreq is therefore an input boundary condition to our model.
- **Load dynamics linearization**: ∂Qreq/∂NP (steady-state torque vs rotor speed) and
  ∂Qreq/∂ṄP (rotor inertia appearing as a lag-hinge shear transient) [pdf p.28, Eq. 55].
- **The report contradicts itself on the units of `TC_T41`.** Printed as
  `lb_m^(4/5)·sec^(9/5)/deg R^(1/2)` in both the nomenclature [pdf p.12] and Table A.1
  [pdf p.55]; Eq. 51 [pdf p.26] requires `sec^(1/5)`. The glyph is unambiguously a 9 at
  600 dpi, so this is the report's own typesetting error, printed twice — not a scan
  artefact. The value 0.29 is unaffected.
- **Station 2 conditions**: stagnation effects roughly offset aircraft inlet losses across
  most of the flight envelope, so P2 and T2 are taken as ambient [pdf p.22].
- **Appendix A holds one numeric table and eleven plots** [pdf pp.55-66]. Table A.1
  (p.55) gives 33 scalar constants and is the appendix's *only* numeric data. Figures
  A1-A11 (pp.56-66) carry every one of the model's functional relationships f₁…f₁₀ and
  f_hs — as printed plots, with no breakpoint lists. The appendix title "...AND FUNCTION
  TABLES" is misleading; the body text at p.22 ("plots of mass-flow and energy
  functions") is accurate. **All eleven must be digitized; none can be transcribed.**

## OCR health

The text layer is a 1988 scan and is unreliable for symbols and numbers. Observed damage:

- `K` renders as `g` or `i_` (`gdamp` = K_damp, `gb1` = **K_bl** — b-ell, bleed, *not*
  b-one, confirmed at 600 dpi against the nomenclature; `i_H2` = K_H2, `Kqc_` = K_QC₁)
- subscripts flatten into the line or vanish (`KT41_`, `Q_w`, `P_3` → `/63`)
- `lbf` → `Ibf`, `lb!`, `lbl`; `deg R` → `de9 1:l`
- minus signs and decimal points drop

**Never transcribe a number from the text layer alone.** Render the page and read it:

```bash
pdftoppm -f 55 -l 55 -r 200 -png docs/ballin-tm100991.pdf /tmp/p
```

Use the text layer for locating things, page images for reading them.

`symbols.md` documents a further trap the text layer cannot show: the report sets **station
indices full size** (`T41` is T-4-1 inline) but **true subscripts small and lowered**, and
sixteen constants carry a subscript *on* a subscript (`K_H41₁` vs `K_H41₂`) where the two
members are slope and intercept of one linear fit and carry **different units**.

## Why the real-time model cannot be integrating the pressures (established 2026-09-11)

Computing our own 5-DOF Jacobian and comparing with Table 1 [pdf p.31]. Central
differences at a relative step of 1e-5; the eigenvalues are converged to 5 significant
figures over rel = 1e-4 .. 1e-8 and agree between one-sided and central schemes, so these
are properties of the model and not of the differencing.

**Re-measured 2026-09-11.** An earlier version of this table reported hover modes of
-4881.0 / -3787.0 / -54.19 / -2.91. Those numbers cannot be reproduced by this code at any
step size or scheme, and `src/t700/engine.py` has not changed since it was written, so they
were never output of the committed model. Replaced rather than explained.

| trim | mode | ours (1/sec) | Table 1 | deviation |
|---|---|---|---|---|
| hover | P41/P45 fast | -4870.39 | -4900. | **-0.6 %** |
| hover | P41/P45 slow | -2603.70 | -3060. | **-14.9 %** |
| hover | P3 | -55.42 | -51.6 | +7.4 % |
| hover | NG | -2.29 | -2.66 | **-14.0 %** |
| hover | NP | -2.79 | -0.565 | not comparable |
| level 80 kt | P41/P45 fast | -4617.32 | -4640. | **-0.5 %** |
| level 80 kt | P41/P45 slow | -3932.36 | -4040. | -2.7 % |
| level 80 kt | P3 | -53.38 | -52.2 | +2.3 % |
| level 80 kt | NG | -2.13 | -2.08 | +2.3 % |
| level 80 kt | NP | -2.13 | -0.446 | not comparable |
| descent 80 kt | P41/P45 fast | -4562.35 | -4530. | +0.7 % |
| descent 80 kt | P41/P45 slow | -4510.95 | -4430. | +1.8 % |
| descent 80 kt | P3 | -55.50 | -52.6 | +5.5 % |
| descent 80 kt | NG | -1.35 | -1.75 | **-22.6 %** |
| descent 80 kt | NP | -1.66 | -0.357 | not comparable |

**Pairing is structural, not positional.** Column NP of `A` is exactly zero off the
diagonal at all three trims -- the same structure Appendix B prints -- so `A[1,1]` *is*
the NP eigenvalue and is identified that way. Sorting the eigenvalues and pairing them
against Table 1's rows in order instead **swaps the NG and NP modes** at hover and
descent, because our NP mode is more damped than our NG mode while Ballin's is far less.
An earlier revision of this table did exactly that and reported the two slow modes as
+5.0 % and -5.0 %; they are in fact -14.0 % and -22.6 %.

The **NP mode is not comparable at any trim**: its eigenvalue depends on dQreq/dNP and on
the load inertia added to J_PT [Eq. 46], both from the Gen Hel UH-60A simulation, which
this report consumes and does not contain. We hold Qreq constant and default `j_load` to
zero. Note the direction: ours is 4-5x *more* damped than Ballin's at every trim, which is
what a missing rotor inertia would do -- open question #6.

**Seven of twelve comparable modes are inside the 4 % Ballin calls "good agreement"
[pdf p.29], nine inside 8 %, and three are outliers**: hover's slow pressure mode
(-14.9 %) and the NG mode at hover (-14.0 %) and descent (-22.6 %). The NG mode is the one
that matters for handling qualities, and it is the one with no clear pattern -- good at
level, poor either side of it.

A structural check that it is not a general defect: Ballin's two fast modes progressively
coalesce as power falls -- the ratio between them is 1.60 at hover, 1.15 at level, 1.02 at
descent. Ours reproduces that trend closely, at 1.87, 1.17 and 1.01. We match the
coalescence and misplace only the separation, at the single condition where the two modes
are furthest apart and therefore most distinguishable.

**The structural consequence.** The fastest mode has a time constant of **0.205 ms**, so
explicit integration of the five states is stable only for dt < 0.41 ms. The report runs
the engine at **7 ms** -- seventeen times over that limit. So the real-time model
**cannot** be integrating P3, P41 and P45 explicitly, and this is not a matter of
preference: it would diverge on the first frame.

That is exactly why Eqs. 74-80 exist. Under the quasi-steady approximation the three
pressures become **algebraic** and are solved each frame by the opened compressor loop plus
the inner P3/P41 fixed-point sweep and the independent P45 iteration. What remains
integrated is NG (tau = 18 ms) and NP (tau ~ 340 ms), both comfortably
stable at 7 and 14 ms.

It also explains the Conclusions' remark [pdf p.54] that omitting the high-speed
inter-volume mass-flow dynamics "was found to be unnecessary": those dynamics are the two
fast modes, they are 20-100x faster than anything the pilot or the rotor can feel, and
carrying them would have cost a 20x smaller time step.

**So the transient model must implement Eqs. 74-80, not Euler on Eqs. 42-44.**
