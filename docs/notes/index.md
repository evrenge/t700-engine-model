# Report index — NASA TM-100991 (Ballin, 1988)

File: `docs/ballin-tm100991.pdf`, 104 PDF pages.

## Citation convention

**Cite PDF page numbers**, written `[TM-100991 pdf p.NN]`. The text layer loses the printed
page numbers, though the page images show them: **printed = pdf − 14** across the
appendices (pdf p.56 is printed p.42). The cover states "95 p", so the two numberings
differ by front matter; never mix them. Mention a printed page only when quoting a caption
that names one.

## Section map (PDF pages)

The `body-*.md` and `inventory-appendix-*.md` notes map their ranges page by page.

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

## Errata (pdf p.5)

Official errata sheet from Ames (Shelley J. Scarich, M/S 241-13):

> On p. 2, figure 1, a correction should be noted. In the upper middle portion of the
> figure, GAS INDUCED FLAP REGRESSIVE should be CONTROL-SYSTEM-INDUCED FLAP REGRESSIVE.

Figure 1 only; there is no numerical erratum.

## Facts established (each with its page)

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
  three control-volume pressures. The vector is printed below Eq. 54 [pdf p.27] and the
  modes are listed in Table 1 [pdf p.31].
- **A 6th state, T41, exists** in the station-4.1 heat-sink variant; the 6-DOF state vector
  {NG, NP, P3, P41, P45, T41} is printed below Eq. 65 [pdf p.32]. These two are the only
  state vectors the report prints. The 3-DOF variant is {NG, NP, T41} — **inferred, not
  printed** (from the structure of G2 and d̄, pdf p.33).
- **The Appendix B 2- and 3-DOF models are NOT order-reduced** from the 5- and 6-DOF
  ones. They are separately *extracted* from a quasi-steady nonlinear simulation
  [pdf p.27, p.33]. Only Table 1's last column ("Reduced-Order 5-DOF Model") is an order
  reduction, and that one is not in Appendix B. Two different tests, both available.
- **Control input is a single scalar**: Wf, fuel flow, lbm/sec [nomenclature pdf p.13].
  Table B.1 quotes trim fuel flow in lbm/**hr** — different unit, same quantity.
- **Appendix B's three trims** [Table B.1, pdf p.67], all sea level, UH-60A at 16825 lbm,
  NP = 20895 rpm: hover (911.1 hp/engine, NG 41638, Wf 476.3 lbm/hr); 80 kt level
  (552.6 hp, NG 39768, Wf 349.3); 80 kt / 1000 ft-min descent, γ = -7.08 deg (302.6 hp,
  NG 38072, Wf 267.7).
- **The NP row of a linear model cannot be checked independently.** Row 2 and b(2) are not
  plain Jacobian entries: Eqs. 57-62 build them from ∂Qreq/∂NP and ∂Qreq/∂ṄP, which come
  from the external Gen Hel UH-60A simulation and are not printed. They are recovered from
  Appendix B itself (ledger #6), so that row is fitted to Appendix B rather than tested
  against it.
- **Time step limit**: max 10 ms, set by a 0.1 % maximum allowable error between time
  steps. Response is time-step dependent because of an *opened* (broken-open across
  frames) compressor mass-flow iteration [pdf p.38].
- **Load torque** Qreq = Qmr + Qtr + Qacc + Qdamp — main rotor, tail rotor, accessories,
  gearbox damping [pdf p.28, Eq. 56]. The report also prints **Q_eng**, torque to the
  helicopter gearbox from both engines, and **Q_GT**, gas generator torque output (see
  `symbols.md`); `Qgb` and `Qgg` in the text layer are OCR inventions. The rotor torques
  come from the **Gen Hel UH-60A blade-element simulation** (an external program, NASA CR
  ref. 2), *not* from this report. Qreq is therefore an input boundary condition to our
  model.
- **Load dynamics linearization**: ∂Qreq/∂NP (steady-state torque vs rotor speed) and
  ∂Qreq/∂ṄP (rotor inertia appearing as a lag-hinge shear transient) [pdf p.28, Eq. 55].
- **The report contradicts itself on the units of `TC_T41`.** Printed as
  `lb_m^(4/5)·sec^(9/5)/deg R^(1/2)` in both the nomenclature [pdf p.12] and Table A.1
  [pdf p.55]; Eq. 51 [pdf p.26] requires `sec^(1/5)`. The glyph is unambiguously a 9 at
  600 dpi, so this is the report's own typesetting error, printed twice — not a scan
  artefact (ledger #4). The printed value is 0.29; Appendix B was generated with about
  0.32 (ledger #31).
- **Station 2 conditions**: stagnation effects roughly offset aircraft inlet losses across
  most of the flight envelope, so P2 and T2 are taken as ambient [pdf p.22].
- **Appendix A holds one numeric table and eleven plots** [pdf pp.55-66]. Table A.1
  (p.55) gives 33 scalar constants and is the appendix's *only* numeric data. Figures
  A1-A11 (pp.56-66) carry every one of the model's functional relationships f₁…f₁₀ and
  f_hs — as printed plots, with no breakpoint lists. The appendix title "...AND FUNCTION
  TABLES" is misleading; the body text at p.22 ("plots of mass-flow and energy
  functions") is accurate. **All eleven must be digitized; none can be transcribed.**
- **The real-time model cannot integrate the three pressures explicitly.** Table 1's
  fastest modes are −4900, −4640 and −4530 /sec at the three trims [pdf p.31], time
  constants near 0.2 ms, so explicit integration of P3, P41 and P45 is stable only below
  about 0.4 ms; the engine runs at 7 ms [pdf p.47]. That is why Eqs. 74-80 exist: the
  pressures become algebraic, solved each frame by the opened compressor loop, the inner
  P3/P41 fixed-point sweep and the separate P45 iteration, and only NG and NP are
  integrated. It is also what the Conclusions mean by "The modeling of high-speed dynamics
  which represent changes of mass flow between the internal control volumes was found to
  be unnecessary" [pdf p.54]. Our eigenvalues against Table 1: `SCOPE.md` and ledger #56.

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
