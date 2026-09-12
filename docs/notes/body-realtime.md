# Real-Time Implementation — PDF pp.26–37 (printed pp.12–23)

Source: Ballin, *A High-Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
NASA TM-100991, July 1988. File `docs/ballin-tm100991.pdf`.

**Method.** Every equation, table cell and quotation below was read from a rendered page
image. Pages 26–37 were rendered at 300 dpi (`docs/extracted/p-0NN.png`); pp.25, 26, 28,
29, 31, 32, 33, 35, 36, 37 were re-rendered at 600 dpi and cropped equation-by-equation
for every sign, subscript and overdot. The OCR text layer was used only to locate things
and, once, as a corroborating second reading of Eq. 74. Printed page = PDF page − 14
throughout this range.

**Scope note.** The section titled **REAL-TIME IMPLEMENTATION** begins mid-page on
pdf p.26 and runs to the top of pdf p.38. Everything above the heading on p.26
(Eqs. 50–53) is the tail of the preceding *Heat-Sink Model* subsection; it is transcribed
here because it falls in the page range, and because Eqs. 50–53 are what Eqs. 63–65
linearize. The first four lines of pdf p.38 are also quoted, because they are the sentence
that open-question #2 (the opened compressor iteration) was raised against.

---

## 1. Page-by-page map

| PDF p. | Printed p. | Content |
|---|---|---|
| 26 | 12 | End of *Heat-Sink Model*: **Eqs. 50, 51, 52, 53**. **REAL-TIME IMPLEMENTATION** section heading. Subsection ***Volume Dynamics Approximation*** opens: why the volume dynamics must be approximated as quasi-steady, and why that needs checking for a small turboshaft coupled to a blade-element rotor sim (refs. 9–12) |
| 27 | 13 | The extracted 5-DOF linear model: **Eq. 54** (`ẋ̄ = A x̄ + b̄ W_f`), the **5-DOF state vector printed explicitly**, and `A` (5×5) and `b̄` (5×1) written out as symbolic partial derivatives. Closes: the model with the volume-dynamics approximation has the two turbine-speed DOF |
| 28 | 14 | Order reduction of the 5-DOF to two states (Ṗ3 = Ṗ41 = Ṗ45 = 0). **The derivative-extraction method**: integrations suppressed, central difference, ±2 %, rotor re-equilibrated, integer rotor revolutions. Load linearization: **Eqs. 55, 56, 57, 58, 59** |
| 29 | 15 | **Eqs. 60, 61, 62** — the `A(2,j)`, `A(2,2)` and `b(2)` elements. Discussion of figure 4 and of Appendix B; linearity vs power level. Discussion of Table 1 (five stable nonoscillatory modes, NP decoupled, three pressure modes out of band). Subsection ***Effects of Heat-Sink Dynamics*** opens |
| 30 | 16 | **Figure 4** — three stacked panels (fuel flow % of trim; gas turbine speed %; power turbine speed %) vs time 0–16 s, comparing Gen Hel UH-60A simulation, 5-DOF, 2-DOF and reduced-order linear models |
| 31 | 17 | **Table 1** — eigenvalues of the three linear models at three trims. **Eq. 63** (heat-sink lead–lag), **Eq. 64**, two unnumbered relations for `ΔT41_ns` and `ΔṪ41_ns`, **Eq. 65** (`F₁ẋ̄ = F₂x̄ + G₁W_f + G₂Ẇ_f`) |
| 32 | 18 | The **6-DOF state vector printed explicitly**, then `F₁` (6×6), `F₂` (6×6), `G₁` (6×1) written out symbolically |
| 33 | 19 | `G₂` (6×1). **Eqs. 66, 67** plus four unnumbered definitions (`A`, `b̄`, `C`, `d̄`) — "as derived by Chen". Discussion of figure 5; station 4.5 heat sink found small and **not** represented; pointer to Figs. B7–B12 |
| 34 | 20 | **Figure 5** — same three panels, 6-DOF vs 5-DOF linear models (with and without station 4.1 heat sink) |
| 35 | 21 | Subsection ***Real-Time Modeling Considerations***: real-time constraints at Ames, asynchronous host, time step set at run time. Why iteration is needed and why it must be bounded. The two previously-used one-pass methods; **Eq. 68** (the NASA Lewis z-domain convergence lag) |
| 36 | 22 | Rejection of both previous methods. Subsection ***Real-Time Iteration Solution***: CDC 7600, Gen Hel UH-60A. **Eqs. 69, 70, 71** (fixed point, Lipschitz, successive overrelaxation), the two nested algebraic loops, **Eqs. 72, 73** (quasi-steady continuity at stations 3 and 4.1) and **Eq. 74 — the opened compressor mass-flow iteration** |
| 37 | 23 | Outer-loop convergence and its time-step-dependent error. The inner loop: **Eqs. 75, 76, 77, 78**; the P45 loop: **Eqs. 79, 80**. Accuracy statements (10 iterations / <0.1 %; 8 iterations / <0.1 % of steady state). Advantages of the method |
| *38* | *24* | *(outside range, quoted here)* "…may be obtained within a specified error tolerance. **Because of the opened compressor mass-flow iteration, however, response is dependent on time step. For the T700 implementation, a restriction on time step to a maximum value of 10 msec was established based on a maximum allowable error of 0.1 percent between time steps.**" Then subsection *Pressure Function Table Solution* (a future alternative, not implemented) |

Two figures (4, 5) and one table (Table 1) live in this range. No other tables.

---

## 2. THE OPENED COMPRESSOR MASS-FLOW ITERATION — answer to open question #2

### 2.1 The closed iteration: what the loops are

[TM-100991 pdf p.36, verbatim]

> "Two nested algebraic loops are present in the representation of pressures `P3` and
> `P41`. The outer loop arises from the calculation of compressor mass flow which is
> indirectly a function of pressure `P3`. The inner loop is a result of the
> interdependence of pressures across the combustor."

Under the quasi-steady (volume-dynamics) approximation the two pressure states are no
longer integrated; their derivatives are set to zero and the resulting **coupled algebraic
equations** must be solved every frame:

```
Ṗ3  = 0 = K_V3  · T3  · (WA3 − WA3_bl − WA31)                    (72)  [p.36]
Ṗ41 = 0 = K_V41 · T41 · (WA31 + W_f − W41)                       (73)  [p.36]
```

**The outer loop** closes through the compressor map. Tracing it with the body equations
of pdf pp.22–23:

```
P3 → Ps3 = K_Ps3·P3              (Eq. 4,  p.22)
   → WA2_c = f₁(Ps3/P2, NG_c)    (Eq. 7,  p.22)   ← the expensive step
   → WA2                          (Eq. 9,  p.22)
   → WA3    = WA2 − WA24_bl       (Eq. 17, p.23)
     WA3_bl = WA2 · (B3 + K_b3)   (Eq. 16, p.23)
   → WA31   = WA3 − WA3_bl        (Eq. 72 with Ṗ3 = 0)
   → P3                           (Eq. 76, p.37, the combustor pressure-drop relation)
```

**The inner loop** closes between `P3` and `P41` across the combustor: Eq. 76 gives `P3`
from `P41` and `WA31`; Eq. 78 gives `P41` from `P3`.

A **closed** scheme would iterate the outer loop — re-entering the compressor map with
each new `P3` — with the inner `P3`/`P41` iteration nested inside it, until both loops are
mutually converged inside a single frame.

### 2.2 What "opening" it means

[TM-100991 pdf p.36, verbatim]

> "Because the calculation of compressor flows involves a significant amount of
> computation, the outer loop was separated from the inner loop by assuming that the mass
> flow entering the combustor for a given interval is approximately equal to the mass flow
> leaving the compressor in the previous interval."

```
WA31_(n) ≈ WA3_(n) − WA3_bl(n−1)                                  (74)  [p.36]
```

The outer loop is **cut, not iterated**. It executes exactly once per time step, and the
value of `WA31` that the inner loop then works with is not the value a converged outer
loop would have produced — it is built from compressor-side mass flows that are one frame
stale. The inner loop *is* iterated (§2.5), so the only thing that is "opened" is the
compressor→`WA31`→`P3` path. That is precisely what pdf p.38 names the "opened compressor
mass-flow iteration".

> ⚠ **Index inconsistency, verified, not a scan artefact.** As printed, Eq. 74 puts `WA3`
> at index `(n)` and only `WA3_bl` at `(n−1)`. Read three times from independent 600-dpi
> crops of the glyph, and the OCR text layer independently agrees (`WA31(,_) ,_ WA3(,_) -
> WA3b_(,___)`). But the prose one line above says the combustor inlet flow equals **the
> mass flow leaving the compressor in the previous interval**, which is `WA3(n−1) −
> WA3_bl(n−1)`; and lagging only the small diffuser bleed would *not* open the loop at
> all, because `WA3(n)` still depends on `P3(n)` through Eqs. 4/7/9/17. Reproduce as
> printed, but see open question #22 (Eq. 74 prose vs equation) — this is the one place in the section where the
> equation and its own sentence disagree.

### 2.3 What state carries across frames, and how it is initialized

**State carried:** one mass flow, `lbm/sec` — the previous frame's diffuser bleed
discharge `WA3_bl(n−1)` as printed, or the previous frame's whole compressor exit flow
`(WA3 − WA3_bl)(n−1)` as the prose has it. Either way this is **memory that is not one of
the five (or six) integrator states**: an extra saved scalar in the frame-to-frame state.
It is the only such quantity introduced by the real-time scheme.

**Initialization: the report does not say.** Searched pp.26–38 and p.47. The nearest
statement is [p.37–38] "An equilibrium solution may be obtained within a specified error
tolerance", which implies a trim solve exists that could seed it, but no initialization
rule is printed. Open question #23 (Eq. 74 initialization at t = 0).

### 2.4 Why the response is time-step dependent

[TM-100991 pdf p.37, verbatim, opening paragraph]

> "This outer loop proved to be strongly convergent. This was verified by examining the
> value of the Lipschitz constant for extreme transient operating cases. With a
> sufficiently small time step, one iteration of the outer loop produces near-equilibrium
> results. However, the approximation results in a value of iteration error which is a
> function of the time step. The time step must therefore be kept as small as possible to
> achieve minimum error and acceptable transient response."

[TM-100991 pdf p.38, verbatim]

> "Because of the opened compressor mass-flow iteration, however, response is dependent on
> time step."

The mechanism, stated plainly: a converged closed iteration would give the same answer for
any Δt, because the algebraic loop has no memory. The opened loop replaces convergence
with a one-frame delay, so the residual it leaves is the amount the compressor mass flow
moves in one frame — of order `(d WA3/dt)·Δt`. That residual is a *modelling* error, not a
truncation error, and it feeds straight into `WA31`, hence into `P3` and `P41`, hence into
`NG`. **Halving Δt halves it. Changing Δt changes the transient.** This is why the report
insists (p.35) that "Simulations must therefore be designed to be independent of time step
size" and then has to concede on p.38 that this one is not.

Note this is the *opposite* of the two rejected methods' failure mode — for those, the
complaint (p.36) was that "the user has no control over the amount of error which is
produced" and "a change of time step will directly affect the magnitude of the error".
Ballin's method keeps the second property but at least makes the first one explicit and
tunable through Δt.

### 2.5 The 0.1 percent criterion — three distinct statements, do not conflate them

| # | Where | Applies to | Exact printed wording |
|---|---|---|---|
| A | p.37 | inner `P3`/`P41` fixed-point loop, **within** a frame | "Under the most extreme conditions, up to ten iterations may be required for convergence with less than 0.1 percent error, resulting in a total of 110 arithmetic operations." |
| B | p.37 | `P45` fixed-point loop, **within** a frame | "Stability and speed of execution were tested under extreme conditions of an instantaneous step in fuel flow from flight idle to full power. Eight iterations resulted in an error equal to less than 0.1 percent of the steady-state value." |
| C | p.38 | the **opened outer loop**, i.e. frame-to-frame | "For the T700 implementation, a restriction on time step to a maximum value of 10 msec was established based on a maximum allowable error of 0.1 percent between time steps." |

Only **C** is the criterion that sets the 10 ms cap, and it is the one that answers the
question as asked. What the report states about it:

- it is an **error between time steps** — i.e. a frame-to-frame discrepancy, consistent
  with the lag in Eq. 74, not an intra-frame convergence tolerance;
- it is **0.1 percent**, so relative;
- the limit it produced for this engine is **Δt ≤ 10 msec**.

What the report does **not** state about C, and which we therefore cannot reproduce
exactly: *which variable* the 0.1 percent is measured on (`WA31`? `P3`? `NG`?), whether it
is relative to the instantaneous value or to a steady-state value (B says "of the
steady-state value" explicitly; C says nothing), and under what transient it was evaluated
(A and B both name "the most extreme conditions" / "flight idle to full power"; C does
not). Open question #24 -- **closed 2026-09-12**: the referent is the *response*, frame to
frame, and the mechanism is named in the same sentence. This paragraph is kept because the
three sub-questions above it are still the right ones to have asked, but it should no longer
be read as "we cannot reproduce this".

A and B are *measurements of worst-case iteration counts*, not runtime policy. Note the
arithmetic in A: 11 operations per pass ("a total of eleven arithmetic operations required
for each pass"), 10 passes, 110 operations.

### 2.6 Answer in one paragraph

The T700 real-time model solves its three quasi-steady pressures by fixed-point iteration.
Two algebraic loops exist: an inner one coupling `P3` and `P41` across the combustor, and
an outer one through the compressor map (`P3 → Ps3 → WA2_c → WA3, WA3_bl → WA31 → P3`).
Iterating the outer loop would mean re-evaluating the compressor performance map several
times per frame, which is the single most expensive computation in the model. Ballin
therefore *opens* it: the outer loop is executed exactly once per frame and the combustor
inlet flow is formed from a compressor mass flow carried over from the previous frame
(Eq. 74), so the loop is never closed within a frame at all — the "iteration" happens once
per time step instead of once per pass. The inner `P3`/`P41` loop is then iterated
normally (up to ten passes, <0.1 % error) with that `WA31` held fixed. The price is that
the outer-loop residual is proportional to how far the compressor flow moves in one frame,
so the engine's transient response depends on Δt; a maximum allowable 0.1 percent error
between time steps caps Δt at 10 ms.

---

## 3. The integration scheme

### 3.1 The integrator itself is never named — a firm negative finding

**The report does not name a numerical integration algorithm anywhere.** Verified by
grepping the text layer of all 104 pages for `euler`, `runge`, `trapezoid`, `adams`,
`predictor`, `corrector`, `solver`, `second-order`: the only hits for "integrat*" in the
body are "suppressing all integrations" (p.28), "suppressing the NP integration" (p.39),
"integrated" in the introduction (p.15), and the ECU/HMU *integrators* of Appendix C
(pp.81–88). Nothing describes how `ẋ` becomes `x`.

What the report gives instead is the **form**: the five states are printed as explicit
integrals [pdf p.25],

```
P3  = K_V3  ∫ T3 (WA3 − WA3_bl − WA31) dt            (42)
P41 = K_V41 ∫ T41(WA31 − W_f − W41) dt               (43)   ← sign, see §7 item 5
P45 = K_V45 ∫ T45(W41 − W45 + B₃ K_bl WA2) dt        (44)
NG  = (60/2π) ∫ (Q_GT − Q_C)/J_GT dt                 (45)
NP  = (60/2π) ∫ (Q_PT − Q_req)/J dt                  (47)
```

and three of these (42, 43, 44) are then *not* integrated at all in the real-time model —
their derivatives are set to zero (Eqs. 72, 73, 79) and they are solved algebraically. So
in the real-time implementation **only `NG` and `NP` are integrated** (plus `T41` if the
heat-sink lag is realized as a state). Open question #26 (no integration algorithm named).

### 3.2 Frame times and rates — stated, and multirate

Three separate numbers, from three separate pages. They must be kept apart.

| Quantity | Value | Citation |
|---|---|---|
| Maximum allowable engine time step | **10 msec** | pdf p.38 |
| Host (Gen Hel UH-60A / rotor) time step actually used | **14 msec** | pdf p.47 |
| Engine model update interval actually used | **7 msec** (twice per rotor cycle) | pdf p.47 |
| `NP` degree-of-freedom update interval | **14 msec** (once per rotor cycle) | pdf p.47 |

[TM-100991 pdf p.47, verbatim, read from the page image]

> "A time step of 14 msec was used; this is a typical value for real-time execution of the
> blade-element rotor. In order to meet the cycling requirements of the iteration model,
> the T700 program was updated twice for each rotor routine cycle, or once every 7 msec.
> The power-turbine-speed degree of freedom corresponds to that of the drive train and
> rotor hub, and was therefore updated every 14 msec."

So the shipped configuration is **multirate, 2:1**: everything in the engine model runs at
7 ms (satisfying the 10 ms cap), except the `NP` state, which runs at 14 ms because it is
physically the drive-train/rotor-hub degree of freedom and belongs to the host's rotor
frame. Note that 14 ms > the 10 ms cap — the cap applies to the compressor mass-flow
iteration, not to `NP`, but the tension is worth recording (open question #33 — 10 ms cap vs 14 ms NP update).

### 3.3 Order of operations within a frame — what is stated and what is not

**Stated, and load-bearing:**

1. **The outer compressor loop runs once, first, on stale data.** Eq. 74 uses a mass flow
   from interval `(n−1)`. §2 above.
2. **The inner `P3`/`P41` loop is a two-equation fixed-point sweep that uses
   partially-updated neighbours.** From p.37: `P3` is obtained from Eq. 76 (which takes
   the current `P41` and the frame's fixed `WA31`); `P41` is then obtained from Eq. 78
   (which takes the `P3` just computed). This is a Gauss–Seidel-style alternation, not a
   simultaneous solve — "A fixed-point iteration was then used to solve the two equations
   for `P3` and `P41`. The iteration converges linearly, with a total of eleven arithmetic
   operations required for each pass. No relaxation algorithm is required for
   convergence." (p.37).
3. **`P45` is solved independently of the other two pressures** — "Pressure upstream of
   the power turbine, `P45`, may be solved independently of the other pressures, but
   because it is a nonlinear function of its own value, iterative techniques must be used"
   (p.37). Its iteration is Eq. 80, one function-table lookup and four arithmetic
   operations per pass, with a numerator that "is a constant over the iteration".
4. **Successive overrelaxation is available but not used for the inner loop.** Eq. 71
   defines it with a relaxation parameter `R < 1` for slower, more stable convergence, and
   p.36 says "A value of `R` was chosen to allow monotonic convergence over the operating
   range of the engine" — but p.37 says of the `P3`/`P41` loop "No relaxation algorithm is
   required for convergence". Where `R` *is* applied is not stated.
5. **`NP` is updated at half the rate of everything else** (§3.2) — so within a 14 ms
   host frame the engine is stepped twice against an `NP` that only moves on the second.

**Not stated, and needed:**

- Whether `NG` and `NP` are integrated before or after the pressure solve, and whether the
  pressure solve sees the old or the new speeds.
- How many iterations are actually run at runtime, and against what stopping test. "Up to
  ten" (inner) and "eight" (`P45`) are *worst-case measurements under extreme step inputs*,
  not a stated policy; p.35 records that Ames "require[s] a minimum amount of multiple-pass
  coding or iteration", which is a pressure to fix the count, but no count is printed.
  Open question #25 (runtime iteration counts and stopping test).
- The numerical value of the relaxation parameter `R` and of the Lipschitz constant `L`.
  Open question #27 (numeric values of R and L).

### 3.4 The real-time constraints Ballin was working under [pdf p.35]

Recorded because they explain every choice above.

> "Typical restrictions that have been used in previous real-time models are the use of
> third-order polynomial curve fits for nonlinear functions, the use of integer exponents
> and square roots only, and allowing no iteration between time steps (ref. 12)."

> "Present implementation restrictions at Ames require a minimum amount of multiple-pass
> coding or iteration and the use of the Ames simulation-standard function table processor
> routines as much as possible. Also, Ames host computers are operated asynchronous to
> peripheral hardware such as the image-generation computers. The time-step size is
> determined for a given simulation depending on computer speed and the amount of
> peripheral computer system loading. **Simulations must therefore be designed to be
> independent of time step size, the only specifiable restrictions being maximum or
> minimum allowed values.**"

Hardware: "It has been executed successfully in conjunction with the Gen Hel UH-60A
simulation using a **CDC 7600** computer and an Ames-developed real-time operating system."
[p.36]

Why not Newton: "Because multiple passes are required to determine elements of the
Jacobian for quadratic convergence methods, they were found to involve more computation
than linearly convergent methods." [p.36]

---

## 4. Every numbered equation, pp.26–37

Transcribed exactly as printed. Sub-subscripted constants are distinguished per
`symbols.md`. Overdots were checked individually at 600 dpi; where an overdot collides
with a fraction bar it is called out.

### 4.1 pdf p.26 — tail of the Heat-Sink Model

```
 T_go     T41     ( Mc_pm/(hA_m) − Mc_pm/(W_g c_pg) ) s + 1
 ──── = ────── = ────────────────────────────────────────── .          (50)
 T_gi    T41_ns            Mc_pm/(hA_m) s + 1
```

```
 Mc_pm              √T41
 ───── = TC_T41 · ────────                                             (51)
 hA_m              W41^(4/5)
```

> `W41` carries the stacked exponent 4/5 (both digits unambiguous at 600 dpi). This is the
> equation that contradicts the printed units of `TC_T41` — see open question #4 (sec exponent 9/5 vs 1/5), not
> re-opened here.

```
 T41_sgn = f_hs(NG_c)                                                  (52)
```

```
 Mc_pm      T41_sgn
 ────── = ─────────                                                    (53)
 W_g c_pg     W41
```

### 4.2 pdf p.27 — the extracted 5-DOF linear model

```
 ẋ̄ = A x̄ + b̄ W_f                                                      (54)
```

with, printed immediately below (**this is the 5-DOF state vector, written out — see §6**):

```
        ⎧ NG  ⎫
        ⎪ NP  ⎪
 x̄ =    ⎨ P3  ⎬
        ⎪ P41 ⎪
        ⎩ P45 ⎭
```

```
       ⎡ ∂ṄG/∂NG   ∂ṄG/∂NP   ∂ṄG/∂P3   ∂ṄG/∂P41   ∂ṄG/∂P45  ⎤
       ⎢ ∂ṄP/∂NG   ∂ṄP/∂NP   ∂ṄP/∂P3   ∂ṄP/∂P41   ∂ṄP/∂P45  ⎥
 A =   ⎢ ∂Ṗ3/∂NG   ∂Ṗ3/∂NP   ∂Ṗ3/∂P3   ∂Ṗ3/∂P41   ∂Ṗ3/∂P45  ⎥
       ⎢ ∂Ṗ41/∂NG  ∂Ṗ41/∂NP  ∂Ṗ41/∂P3  ∂Ṗ41/∂P41  ∂Ṗ41/∂P45 ⎥
       ⎣ ∂Ṗ45/∂NG  ∂Ṗ45/∂NP  ∂Ṗ45/∂P3  ∂Ṗ45/∂P41  ∂Ṗ45/∂P45 ⎦
```

```
       ⎡ ∂ṄG/∂W_f  ⎤
       ⎢ ∂ṄP/∂W_f  ⎥
 b̄ =   ⎢ ∂Ṗ3/∂W_f  ⎥
       ⎢ ∂Ṗ41/∂W_f ⎥
       ⎣ ∂Ṗ45/∂W_f ⎦
```

(The overdot sits on the leading letter of each station-indexed symbol: `Ṗ41` is printed
`P` with a dot, then `41`.)

### 4.3 pdf p.28 — load linearization

```
             ∂Q_req         ∂Q_req
 ΔQ_req = ─────────· ΔNP + ─────── · ΔṄP                               (55)
              ∂NP           ∂ṄP
```

> Second denominator: `∂ṄP`, dot present. In the printed second term both the numerator
> `∂Q_req` and the denominator `∂ṄP` are as shown; the *first* term's denominator is plain
> `∂NP` with no dot. Verified at 600 dpi.

```
 Q_req = Q_mr + Q_tr + Q_acc + Q_damp                                  (56)
```

```
              ∂Q_eng          ∂Q_eng          ∂Q_eng
 ΔQ_eng = Σ ────────· Δx_j + ──────── · ΔNP + ──────── · ΔW_f          (57)
                ∂x_j            ∂NP             ∂W_f
```

```
 ΔṄP = (ΔQ_eng − ΔQ_req)/J                                             (58)
```

```
                 1              ⎧ 1     ∂Q_eng          1 ⎛ ∂Q_eng   ∂Q_req ⎞          1  ∂Q_eng       ⎫
 ΔṄP = ───────────────────────  ⎨ ─ Σ ────────· Δx_j + ─ ⎜ ────── − ────── ⎟ ΔNP  +  ─ ──────── ΔW_f ⎬  (59)
        1 + (1/J)(∂Q_req/∂ṄP)   ⎩ J     ∂x_j            J ⎝  ∂NP      ∂NP   ⎠          J   ∂W_f        ⎭
```

> Denominator of the leading fraction: `1 + (1/J)(∂Q_req/∂ṄP)` — the overdot on the `N` is
> fused into the fraction bar above it, visible as a downward tick at 600 dpi. The two
> `∂NP` inside the round bracket carry **no** dot; those are steady-state derivatives.

### 4.4 pdf p.29 — the NP row of the state and control matrices

```
             ∂Q_eng/∂x_j
 A_{2,j} = ───────────────────── ,   j = 1, 3 − 5                      (60)
            J + ∂Q_req/∂ṄP
```

```
             ∂Q_eng/∂NP − ∂Q_req/∂NP
 A_{2,2} = ───────────────────────────                                 (61)
                 J + ∂Q_req/∂ṄP
```

```
             ∂Q_eng/∂W_f
 b_2 = ─────────────────────                                           (62)
         J + ∂Q_req/∂ṄP
```

> All three denominators are `J + ∂Q_req/∂ṄP` (dot present, fused with the bar). Eq. 61's
> two numerator terms are `∂NP` with **no** dot. Checked side by side at 600 dpi: the
> denominator bars carry the tick, the numerator bars do not. Dimensionally consistent:
> `J` is `ft·lb_f·sec²`, and `∂Q_req/∂ṄP` is `ft·lb_f / (rad/sec²)`.
>
> Note the index range is printed `j = 1, 3 − 5` — i.e. j ∈ {1, 3, 4, 5}, all states
> except NP itself.

### 4.5 pdf p.31 — heat-sink linear representation

```
  T41       τ₁ s + 1
 ────── = ───────────                                                  (63)
 T41_ns     τ₂ s + 1
```

```
            1              τ₁                1
 ΔṪ41 = − ── ΔT41  +  ── ΔṪ41_ns  +  ── ΔT41_ns                        (64)
            τ₂             τ₂               τ₂
```

> Term 2 carries an overdot on `T41_ns`; term 3 does not. Verified at 600 dpi. The overdot
> in the report sits over the `4` of `T41` rather than the `T` in this equation — a
> typesetting habit, not a different symbol.

Then two **unnumbered** displayed relations, "where":

```
              ∂T41_ns            ∂T41_ns
 ΔT41_ns = Σ ───────── Δx_i  +  ───────── ΔW_f
                ∂x_i               ∂W_f
```

```
              ∂T41_ns            ∂T41_ns
 ΔṪ41_ns = Σ ───────── Δẋ_i  +  ───────── ΔẆ_f
                ∂x_i               ∂W_f
```

> Note the *numerators* of the second relation are `∂T41_ns` (not `∂Ṫ41_ns`); only the
> perturbations `Δẋ_i` and `ΔẆ_f` carry dots. That is what makes `F₁`/`G₂` below come out
> the way they do.

```
 F₁ ẋ̄ = F₂ x̄ + G₁ W_f + G₂ Ẇ_f                                        (65)
```

### 4.6 pdf p.32 — the 6-DOF arrays

Printed immediately below Eq. 65:

```
        ⎧ NG  ⎫
        ⎪ NP  ⎪
        ⎪ P3  ⎪
 x̄ =    ⎨ P41 ⎬
        ⎪ P45 ⎪
        ⎩ T41 ⎭
```

```
        ⎡ 1   0   0   0   0   0 ⎤
        ⎢ 0   1   0   0   0   0 ⎥
        ⎢ 0   0   1   0   0   0 ⎥
 F₁ =   ⎢ 0   0   0   1   0   0 ⎥
        ⎢ 0   0   0   0   1   0 ⎥
        ⎣ r₆₁ r₆₂ r₆₃ r₆₄ r₆₅  1 ⎦
```
with row 6 printed as
```
 r₆₁ = −(τ₁/τ₂)(∂T41_ns/∂NG)      r₆₄ = −(τ₁/τ₂)(∂T41_ns/∂P41)
 r₆₂ = −(τ₁/τ₂)(∂T41_ns/∂NP)      r₆₅ = −(τ₁/τ₂)(∂T41_ns/∂P45)
 r₆₃ = −(τ₁/τ₂)(∂T41_ns/∂P3)
```

```
        ⎡ ∂ṄG/∂NG   ∂ṄG/∂NP   ∂ṄG/∂P3   ∂ṄG/∂P41   ∂ṄG/∂P45   ∂ṄG/∂T41_ns  ⎤
        ⎢ ∂ṄP/∂NG   ∂ṄP/∂NP   ∂ṄP/∂P3   ∂ṄP/∂P41   ∂ṄP/∂P45   ∂ṄP/∂T41_ns  ⎥
        ⎢ ∂Ṗ3/∂NG   ∂Ṗ3/∂NP   ∂Ṗ3/∂P3   ∂Ṗ3/∂P41   ∂Ṗ3/∂P45   ∂Ṗ3/∂T41_ns  ⎥
 F₂ =   ⎢ ∂Ṗ41/∂NG  ∂Ṗ41/∂NP  ∂Ṗ41/∂P3  ∂Ṗ41/∂P41  ∂Ṗ41/∂P45  ∂Ṗ41/∂T41_ns ⎥
        ⎢ ∂Ṗ45/∂NG  ∂Ṗ45/∂NP  ∂Ṗ45/∂P3  ∂Ṗ45/∂P41  ∂Ṗ45/∂P45  ∂Ṗ45/∂T41_ns ⎥
        ⎣ s₆₁       s₆₂       s₆₃       s₆₄       s₆₅        −1/τ₂        ⎦
```
with row 6 printed as
```
 s₆₁ = (1/τ₂)(∂T41_ns/∂NG)        s₆₄ = (1/τ₂)(∂T41_ns/∂P41)
 s₆₂ = (1/τ₂)(∂T41_ns/∂NP)        s₆₅ = (1/τ₂)(∂T41_ns/∂P45)
 s₆₃ = (1/τ₂)(∂T41_ns/∂P3)
```

> Column 6 of `F₂` differentiates with respect to `T41_ns`, **not** `T41` — the `ns`
> subscript is present on all six entries and was re-read at 600 dpi. `F₂(6,6) = −1/τ₂`.

```
        ⎡ ∂ṄG/∂W_f            ⎤
        ⎢ ∂ṄP/∂W_f            ⎥
        ⎢ ∂Ṗ3/∂W_f            ⎥
 G₁ =   ⎢ ∂Ṗ41/∂W_f           ⎥
        ⎢ ∂Ṗ45/∂W_f           ⎥
        ⎣ (1/τ₂)(∂T41_ns/∂W_f)⎦
```

### 4.7 pdf p.33 — the Chen transformation

```
        ⎡ 0                    ⎤
        ⎢ 0                    ⎥
        ⎢ 0                    ⎥
 G₂ =   ⎢ 0                    ⎥
        ⎢ 0                    ⎥
        ⎣ (τ₁/τ₂)(∂T41_ns/∂W_f)⎦
```

```
 ż̄ = A z̄ + b̄ W_f                                                      (66)
```

```
 x̄ = C z̄ + d̄ W_f                                                      (67)
```

then four **unnumbered** displayed definitions, "where":

```
 A = F₁⁻¹ F₂
```
```
 b̄ = F₁⁻¹ G₁ + F F₁⁻¹ G₂
```
```
 C = I
```
```
 d̄ = F₁⁻¹ G₂
```

> The bare `F` in the `b̄` definition is printed exactly as shown — an upright italic `F`
> with no subscript, re-read at 600 dpi. It is defined nowhere in the report and is not in
> the nomenclature. Open question #32 (the bare F in b̄).

### 4.8 pdf p.35 — the rejected NASA Lewis convergence lag

```
 x_{n+1}        K z
 ────── = ─────────────                                                (68)
   x_n      z − e^(−KΔt)
```

> Background only — this is refs. 9/10's method, which Ballin evaluated and **did not
> use**: "the method was found to give poor results under large power-transient operation
> as experienced by small turboshaft engines. Stability is also lowered under low-power
> conditions." Note Eq. 68 uses true subscripts `n+1`, `n` (time indices), unlike
> Eqs. 69–71 and 74 which use parenthesized `(n)`.

### 4.9 pdf p.36 — the real-time iteration solution

```
 x_(n) = f(x_(n−1))                                                    (69)
```

```
 |f(x_(n−1)) − f(x_(n))| ≤ L |x_(n−1) − x_(n)|                         (70)
```

```
 x_(n) = x_(n−1) + R · [ f(x_(n)) − x_(n−1) ]                          (71)
```

> **Printed exactly as shown.** `f(x_(n))` appears on the right-hand side while `x_(n)` is
> on the left, which makes Eq. 71 implicit as printed; the standard successive-
> overrelaxation form would be `f(x_(n−1))`. Re-read at 600 dpi — the subscript is `(n)`,
> not `(n−1)`. Reproduce as printed and record the objection: open question #28 (Eq. 71 implicit as printed).

```
 Ṗ3 = 0 = K_V3 · T3 · (WA3 − WA3_bl − WA31)                            (72)
```

```
 Ṗ41 = 0 = K_V41 · T41 · (WA31 + W_f − W41)                            (73)
```

> `+ W_f`, verified at 600 dpi. This **contradicts Eq. 43 on pdf p.25**, which prints
> `WA31 − W_f − W41`. Eq. 73 is the physically correct one (fuel adds mass to the
> combustor) and is the one consistent with Eq. 77 below. Open question #29 (Eq. 73 +W_f vs Eq. 43 −W_f).

```
 WA31_(n) ≈ WA3_(n) − WA3_bl(n−1)                                      (74)
```

> **The opened compressor mass-flow iteration.** See §2, including the flagged index
> inconsistency with the prose.

### 4.10 pdf p.37 — the inner loop and the P45 loop

```
 P3² − P3 · P41 − K_dpb · T3 · WA31² = 0                               (75)
```

```
        1 ⎛                                            ⎞
 P3 =  ── ⎜ P41 + √( P41² + 4 · K_dpb · T3 · WA31² )   ⎟               (76)
        2 ⎝                                            ⎠
```

```
 WA31 = W41 − W_f                                                      (77)
```

```
                   T3 · K_dpb · ( K_WGT·P41/√θ₄₁ − W_f )²
 P41 = P3 − ──────────────────────────────────────────────             (78)
                                  P3
```

```
 Ṗ45 = 0 = K_V45 · T45 · (W41 − W45 + B₃ B₄ WA2)                       (79)
```

> `B₄` is printed unambiguously (an italic capital B with a subscript 4, re-read at 600
> dpi). **`B₄` is not in the nomenclature** — only `B₁`, `B₂`, `B₃` are — and Eq. 44 on
> pdf p.25 writes the identical physical term as `B₃ K_bl WA2`. Open question #30 (B₄ absent from the nomenclature).

```
        (W41 + B₃ B₄ WA2) √θ₄₅
 P45 = ────────────────────────                                        (80)
             f₉( P_s9 / P45 )
```

> Consistent with Eqs. 33/34 on pdf p.24 (`W45_c = f₉(P_s9/P45)`, `W45 = W45_c · P45/√θ₄₅`)
> with `W45` eliminated using Eq. 79. The `f₉` argument is `P_s9/P45` — small subscript `s`
> on the `P`, station-9 *static* pressure.

**Internal consistency check of Eqs. 75–78, performed here as a transcription test.**
Eq. 76 is the positive root of Eq. 75 by the quadratic formula:
`P3 = ½(P41 + √(P41² + 4·K_dpb·T3·WA31²))` — exact match. Eq. 78 follows from Eq. 75
rearranged as `P41 = P3 − K_dpb·T3·WA31²/P3`, with `WA31` replaced by Eq. 77
(`W41 − W_f`) and `W41` by the choked-nozzle relation `K_WGT·P41/√θ₄₁` — exact match. A
dropped sign, a lost factor of 4, or a mis-read `K_dpb`/`K_WGT` anywhere in this group
would break the chain. It does not break: **Eqs. 75–78 are verified end to end, not merely
read twice.** (Eq. 74 sits outside this chain and cannot be checked this way — which is
part of why flag §8 item 1 matters.)

### 4.11 Count

**31 numbered equations** (50 through 80 inclusive), **6 unnumbered displayed relations**
(2 on p.31, 4 on p.33), and **8 displayed symbolic arrays** (`x̄`, `A`, `b̄` on p.27;
`x̄`, `F₁`, `F₂`, `G₁` on p.32; `G₂` on p.33).

---

## 5. The linear-model derivation

### 5.1 Extraction method [pdf p.28, verbatim]

> "Stability derivatives were extracted from the nonlinear simulations by suppressing all
> integrations after trimming the simulation at a desired operating condition. Each state
> or control was then individually perturbed and the effects of its change were determined
> in the state derivatives at the input of each integrator. The central-difference
> extraction method that was used perturbs each state or control variable in both positive
> and negative directions. A perturbation step-size of plus or minus two percent of
> equilibruim [*sic*] conditions was considered to be adequate based on the rotor speed
> deviations experienced by current-generation helicopters. During the derivative
> extraction process, the rotor was allowed to reach a new equilibrium after each
> perturbation to eliminate dynamics of the rotor blade lag degree-of-freedom. The rotor
> was allowed to rotate an exact-integer number of revolutions between extractions,
> resulting in rotor blade azimuths which were unchanging for each extraction. In this way,
> first harmonic changes of required torque which occur in blade-element rotor simulations
> were eliminated from the linear model."

Reproducible procedure, stated concretely:

1. Trim the nonlinear model at the operating condition.
2. **Suppress all integrations** — freeze the states.
3. For each state `x_j` and for the control `W_f` in turn, perturb by **±2 % of the
   equilibrium value**, evaluate the state derivatives at the input of each integrator, and
   form the **central difference** `(ẋ⁺ − ẋ⁻)/(2·0.02·x_j⁰)`.
4. Between extractions, let the rotor re-equilibrate and rotate a whole number of
   revolutions so the blade azimuths are identical every time.

Steps 3–4 are ours to reproduce; step 4 only matters when the load comes from the
blade-element rotor sim, which for us it does not (`Q_req` is a boundary condition).

### 5.2 The load-torque derivatives, and why the NP row is special

The load model is **two linear elements** [p.28]:

- `∂Q_req/∂NP` — "used to model steady-state change of torque with a change of rotor
  speed";
- `∂Q_req/∂ṄP` — "needed to model the rotor inertia which manifests itself as a shear
  force transient at the rotor hub lag-hinge".

Both come from applying the same small-perturbation extraction to the **Gen Hel UH-60A
blade-element helicopter simulation** — an external program. Eq. 55 combines them.

`J` is defined narrowly [p.28]: "the sum of inertias of the power turbine, drive train,
and rotor hub. **Rotor blade inertias are not included in `J`**; effects of these inertias
are reflected in the `∂Q_req/∂ṄP` term."

The consequence for us: Eqs. 60–62 show that `∂Q_req/∂ṄP` divides **every element of the
NP row** and `b(2)`, and `∂Q_req/∂NP` additionally shifts `A(2,2)`. Neither derivative is
printed anywhere in the report. Row 2 and `b(2)` of every Appendix B matrix are therefore
not independently reproducible from this report — confirming what
`inventory-appendix-b.md` §6.5 already records.

`A(2,2)` can be back-solved, though, and cheaply: at each trim the 2-, 3-, 5- and 6-DOF
models all print the *same* `A(2,2)` (−0.5650, −0.4461, −0.3567), and Table 1 gives it as
the NP eigenvalue. That is one printed number per trim constraining the two unknown load
derivatives — enough to check a hypothesis, not enough to determine both.

### 5.3 The three linear models of Table 1 — and how they are actually formed

**This is where the task brief's framing needs correcting.** The Appendix B 2-DOF and
3-DOF models are **not** order-reduced from the 5-DOF and 6-DOF models. There are three
distinct constructions, and only one of them is an order reduction:

| Model | How it is built | Where it appears |
|---|---|---|
| **5-DOF** | Extracted from *the complete nonlinear digital simulation* (full volume dynamics) | Table 1 col. 3; Figs. B2, B4, B6 |
| **2-DOF** | **Extracted** from *a reduced-order nonlinear simulation which contains the volume dynamics approximation* [p.27] — a different nonlinear model, not a reduction of the linear one | Table 1 col. 4; Figs. B1, B3, B5 |
| **Reduced-order 5-DOF** | **Order-reduced** from the 5-DOF linear model: "Order reduction was performed by setting the state derivatives `Ṗ3`, `P̈41`, and `Ṗ45` to zero and solving the equations for the remaining two states." [p.28] | Table 1 col. 5 **only** — it is *not* in Appendix B |

[p.27, closing line] "The model containing the nonlinear volume dynamics approximation is
represented by the two turbine-speed degrees-of-freedom."
[p.28, opening line] "For a further comparison, a third linear model was created by
eliminating the states corresponding to the volume dynamics in the five-degree-of-freedom
model."

The same pattern holds with the heat sink [p.33, verbatim]: "Figures B7 through B12
present the six-degree-of-freedom **extracted** linear models and three-state **extracted**
models which contain the volume dynamics approximation for the three flight conditions."
So the 3-DOF models are extracted from a reduced-order nonlinear simulation too — there is
no printed order reduction of the 6-DOF model anywhere.

**Consequence for validation.** To reproduce Figs. B1/B3/B5 (2-DOF) and B7/B9/B11 (3-DOF)
we must extract from *our quasi-steady nonlinear model* (pressures solved algebraically),
not from a reduction of our full-dynamics Jacobian. To reproduce Table 1's third column we
must do the algebraic reduction of our 5-DOF `A`. These are two different tests and both
are available.

The order reduction itself, written out: partition `ẋ̄ = A x̄ + b̄ W_f` into speeds
`x̄₁ = {NG, NP}` and pressures `x̄₂ = {P3, P41, P45}`; set `ẋ̄₂ = 0`, giving
`x̄₂ = −A₂₂⁻¹(A₂₁ x̄₁ + b̄₂ W_f)`; substitute to get
`ẋ̄₁ = (A₁₁ − A₁₂A₂₂⁻¹A₂₁) x̄₁ + (b̄₁ − A₁₂A₂₂⁻¹b̄₂) W_f`. **The report does not print this
algebra** — only the sentence quoted above. Recorded here as the obvious reading; it is
checkable against Table 1 column 5.

### 5.4 The heat-sink linear model (Eqs. 63–67)

The station 4.1 heat sink is a lead–lag on `T41` (Eq. 63). Because `T41_ns` is a function
of the other five states *and* of `W_f`, differentiating the lead term introduces a `Ẇ_f`
input — hence the descriptor form Eq. 65 rather than plain state-space, and hence the Chen
transformation of Eqs. 66–67 to recover a standard form. The physical state is recovered
as `x̄ = z̄ + d̄ W_f` (since `C = I`).

The extraction step for the extra state [p.31, verbatim]:

> "During derivative extraction, `T41` was held fixed at the trim value and the change in
> `T41` was determined for each state perturbation. The time constants `τ₁` and `τ₂` are
> those values corresponding to the trim operating condition."

> The second `T41` is printed with no `ns` subscript, verified at 600 dpi. Read
> literally the sentence is circular (hold `T41` fixed, measure the change in `T41`); it
> is almost certainly `T41_ns`, which is what Eqs. 64/F₁/F₂/G₁/G₂ actually require.
> Transcribed as printed. Open question #31 (heat-sink τ₁, τ₂ values) covers the related fact that no numeric `τ₁`,
> `τ₂` appear anywhere in the report.

Effect of including the heat sink [p.33, verbatim]:

> "Figure 5 illustrates the importance of the use of a station 4.1 heat-sink model.
> Response is significantly slowed compared to the five-degree-of-freedom model. The
> overall effect of the heat-sink representation is a larger variation of rotor shaft
> speed and decreased closed-loop system stability. Analysis-oriented simulations such as
> the performance standard component-model program developed by the engine manufacturer
> also model a heat sink effect at station 4.5. Because its influence was found to be
> small, the station 4.5 heat-sink effect is not represented in the real-time model."

### 5.5 Table 1 — full transcription [TM-100991 pdf p.31, printed p.17]

**"Table 1: Eigenvalues of two extracted perturbation models and the reduced-order
model."** Read cell by cell from the 300-dpi page image; blank cells are blank in the
original (the 2-DOF and reduced-order models have only two modes each). Units are not
printed; they are 1/sec.

| Trim Condition | Mode | 5-DOF Model | 2-DOF Model | Reduced-Order 5-DOF Model |
|---|---|---|---|---|
| 1 | NG  | -2.66  | -2.69  | -2.81  |
| 1 | NP  | -0.565 | -0.565 | -0.565 |
| 1 | P3  | -51.6  | | |
| 1 | P41 | -4900. | | |
| 1 | P45 | -3060. | | |
| 2 | NG  | -2.08  | -2.23  | -2.16  |
| 2 | NP  | -0.446 | -0.446 | -0.446 |
| 2 | P3  | -52.2  | | |
| 2 | P41 | -4640. | | |
| 2 | P45 | -4040. | | |
| 3 | NG  | -1.75  | -1.82  | -1.83  |
| 3 | NP  | -0.357 | -0.357 | -0.357 |
| 3 | P3  | -52.6  | | |
| 3 | P41 | -4430. | | |
| 3 | P45 | -4530. | | |

Reading notes: the trailing decimal points on the four-digit values (`-4900.`, `-3060.`,
`-4640.`, `-4040.`, `-4430.`, `-4530.`) are printed. Trim numbers 1/2/3 are set once,
vertically centred over each five-row block. Mode labels are italic (`NG`, `NP`, `P3`,
`P41`, `P45`).

This table **verifies the Appendix B transcription end to end** — the eigenvalues of the
three transcribed 5×5 `A` matrices reproduce column 3 to three significant figures, and
the 2-DOF diagonals reproduce column 4 exactly. See `inventory-appendix-b.md` §6.3. Use it
the same way on our own extracted matrices.

What the report says about the table [p.29, verbatim]:

> "The five-degree-of-freedom extracted model consists of five stable, nonoscillatory
> modes. The power turbine state is completely decoupled from the other states. The
> neglected dynamics of the simplified load representation therefore have no effect on the
> gas-generator turbine or volume dynamics. Three modes, corresponding to the pressure
> states, are outside the bandwidth of interest for the modeling of rotor and propulsion
> system dynamics. These modes are well above the second drive-train torsional mode and
> the rotor blade first in-plane elastic mode of a typical helicopter. Also, the
> reduced-order two-state model eigenvalues are seen to be in good agreement with those of
> the extracted two-state model. **The volume dynamics approximation is therefore
> considered to be valid for real-time and nonreal-time engine dynamic modeling.**"

---

## 6. Corrections to existing notes

1. **`inventory-appendix-b.md` §2.1 says the 6-DOF state vector on pdf p.32 "is the only
   place in the report where a linear-model state vector is written out." That is wrong.**
   The **5-DOF state vector `{NG, NP, P3, P41, P45}` is printed explicitly on pdf p.27**,
   directly below Eq. 54, in the same brace notation. The p.32 vector is the second such
   display, not the only one. (The 2-DOF and 3-DOF orderings remain uninferred/unprinted,
   so §2.4's `# UNVERIFIED-ORDERING` flag still stands for the 3-DOF case.)

2. **`inventory-appendix-b.md` §6.5 and `SCOPE.md` Phase 6 describe the Appendix B 2-/3-DOF
   models as order-reduced.** They are not; they are separately *extracted* from a
   quasi-steady nonlinear simulation. Only Table 1's third column is an order reduction,
   and it is not in Appendix B. See §5.3 above. This changes how Phase 6 validates them.

3. **`SCOPE.md` "What the report actually contains"** lists the 10 ms cap without the
   actual run configuration. Add: the report's own results were produced with a **14 ms
   host step, the engine updated every 7 ms, and `NP` updated every 14 ms** [pdf p.47].

---

## 7. The report's own accuracy claims (pp.26–38)

These are the statements that can replace guessed tolerances. Quoted exactly, with page.

| # | Claim | Page | What it bounds |
|---|---|---|---|
| 1 | "Under the most extreme conditions, up to ten iterations may be required for convergence with **less than 0.1 percent error**, resulting in a total of 110 arithmetic operations." | p.37 | The inner `P3`/`P41` fixed-point iteration, within a frame. Worst case, extreme transient |
| 2 | "Stability and speed of execution were tested under extreme conditions of an instantaneous step in fuel flow from flight idle to full power. **Eight iterations resulted in an error equal to less than 0.1 percent of the steady-state value.**" | p.37 | The `P45` fixed-point iteration, within a frame. Explicitly relative to steady state |
| 3 | "For the T700 implementation, a restriction on time step to a maximum value of **10 msec** was established based on a **maximum allowable error of 0.1 percent between time steps**." | p.38 | The opened outer loop. Frame-to-frame. This is the one that sets Δt |
| 4 | "This outer loop proved to be **strongly convergent**. This was verified by examining the value of the Lipschitz constant for extreme transient operating cases. With a sufficiently small time step, one iteration of the outer loop produces **near-equilibrium results**." | p.37 | Qualitative; no number given for `L` |
| 5 | "The reduced-order linear model response is **nearly identical** to that of the extracted linear models." | p.29 | Qualitative — Fig. 4 |
| 6 | "Both the complete five-degree-of-freedom and the two-degree-of-freedom perturbation models are **good representations** of engine response. **Errors associated with changes in trim are small**, and transient responses of internal engine states are **reproduced well**." | p.29 | Qualitative — Fig. 4 |
| 7 | "the reduced-order two-state model eigenvalues are seen to be **in good agreement** with those of the extracted two-state model" | p.29 | Qualitative — Table 1. Numerically this is 4 %, 3 % and 0.5 % on the NG mode at trims 1/2/3, and exact on the NP mode |
| 8 | (about **ref. 11's** method, not Ballin's) "When applied to a large turbofan engine, the difference from a full iteration method was **less than one percent** for severe ramp inputs" | p.36 | Not a claim about this model. Do not adopt |

**Recommended use.** Only claims 1–3 are numeric and only 3 is a *model-level* tolerance.
None of them bounds `NG`, `NP`, pressures, temperatures or torque against reference data —
those bounds, if the report states them at all, are in the VALIDATION section (pp.38–44)
and in the results (pp.45–53), which are another agent's range. Concretely usable now:

- **Claim 3 is a test, not a tolerance**, and it is the one `SCOPE.md` Phase 4 already
  asks for: run the model at Δt and Δt/2 and assert the frame-to-frame difference metric
  stays under 0.1 % at Δt = 10 ms. **Still not implemented** -- `validation/test_transient.py`
  checks the convergence *order* at 4/2/1/0.5 ms, not the 0.1 % at 10 ms that `SCOPE.md`
  Phase 4 asks for. (#24 is closed, so the referent is no longer the obstacle: it is the
  frame-to-frame response.)
- Claim 7 quantifies, indirectly, how much of a mismatch Ballin was willing to call "good
  agreement" on a *mode*: **up to about 4 %**. That is a defensible tolerance for the
  eigenvalue comparison, and it is looser than `SCOPE.md`'s guessed ±5 % per matrix
  element is tight.
- Nothing here justified changing the state tolerances in `SCOPE.md` at the time. **Both
  halves of that sentence are now spent**: pp.38-54 were read (`body-validation.md`), and
  `SCOPE.md`'s tolerance table was rewritten from the report on 2026-09-10 -- it has no
  ±1 % row any more. The ±5 %/element figure is no longer "guessed" either; SCOPE now
  attributes it to Ballin's own 4 % "good agreement".

---

## 8. Ambiguities and flags raised by this section

| # | Item | Page | Status |
|---|---|---|---|
| 1 | Eq. 74 prints `WA3_(n)` but its own prose says the previous interval | p.36 | Verified at 600 dpi *and* against the OCR layer; both give `(n)`. Prose/equation genuinely disagree. Open question #22 (Eq. 74 prose vs equation) |
| 2 | Eq. 71 prints `f(x_(n))` on the RHS — implicit as printed | p.36 | Verified at 600 dpi. Standard SOR would be `f(x_(n−1))`. Open question #28 (Eq. 71 implicit as printed) |
| 3 | `B₄` in Eqs. 79, 80 vs `K_bl` in Eq. 44 | pp.37, 25 | `B₄` verified at 600 dpi; not in the nomenclature. Open question #30 (B₄ absent from the nomenclature) |
| 4 | `+W_f` (Eq. 73) vs `−W_f` (Eq. 43) in the station 4.1 continuity equation | pp.36, 25 | Both verified at 600 dpi. Eq. 73 is physically right and consistent with Eq. 77. Open question #29 (Eq. 73 +W_f vs Eq. 43 −W_f) |
| 5 | Bare `F` in `b̄ = F₁⁻¹G₁ + F F₁⁻¹G₂` | p.33 | Verified at 600 dpi; undefined in the report. Open question #32 (the bare F in b̄) |
| 6 | "the change in `T41` was determined" should read `T41_ns` | p.31 | Verified at 600 dpi; no subscript printed. Recorded, does not block |
| 7 | Overdots on `∂ṄP` in Eqs. 55, 59, 60, 61, 62 collide with the fraction bar | pp.28, 29 | Resolved: the tick is present in every denominator and absent in every numerator. Checked side by side at 600 dpi. Dimensional check confirms |
| 8 | No integration algorithm is named anywhere in the report | all | Verified by grep over all 104 text pages. Open question #26 (no integration algorithm named) |
| 9 | 10 ms cap (p.38) vs 14 ms `NP` update (p.47) | pp.38, 47 | Both verified from page images. Open question #33 (10 ms cap vs 14 ms NP update) |
| 10 | `R` (relaxation) and `L` (Lipschitz) have no printed values | pp.36, 37 | Open question #27 (numeric values of R and L) |
