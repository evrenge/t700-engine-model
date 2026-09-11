# Report body — Engine model equations, PDF pp.15–25

Source: Ballin, *A High-Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
NASA TM-100991, 1988 — `docs/ballin-tm100991.pdf`.

**Status: complete for pp.15–25.** Every numbered equation on these pages (Eqs. 1–49) is
transcribed below. Everything here was read from page images rendered at 300 dpi with
600 dpi re-reads and pixel-level zooms on every sign, sub-subscript and broken glyph.
Nothing in this file comes from the OCR text layer.

Printed page = PDF page − 14 across this range (pdf p.15 is printed p.1). Cite the PDF
page: `[TM-100991 pdf p.NN, Eq. N]`.

Symbols are cross-referenced to `symbols.md` (the 106-entry nomenclature). Typesetting
conventions from that file are preserved here: **station indices are set full size**
(`T41`, `P45`, `WA31`), **true subscripts are set small** (`H41_ns`, `P_s3`, `W_f`), and
**sub-subscripted constants are two distinct constants** (`K_H41₁` ≠ `K_H41₂`, different
units).

---

## 1. Page-by-page map, pdf pp.15–25

| PDF | Printed | Content |
|---|---|---|
| 15 | 1 | Title block (*A HIGH FIDELITY REAL-TIME SIMULATION OF A SMALL TURBOSHAFT ENGINE*, Mark G. Ballin, Ames Research Center). **SUMMARY** (one paragraph). **INTRODUCTION** begins — motivation, rotorcraft handling qualities, why partial-derivative engine models are inadequate. No equations. |
| 16 | 2 | **Figure 1**, full page: *Modal frequencies of interest in engine and fuel-control design and modeling (ref. 3).* Frequency axis 0.1–50 Hz, log. Bars above the axis = airframe/rotor modes, hatched bars below = drive-system modes. Spans labelled HELICOPTER MANUFACTURER, PERFORMANCE, HANDLING QUALITIES, AUTOMATIC FLIGHT CONTROLS, DYNAMICS, ENGINE MANUFACTURER. **This is the figure the errata sheet (pdf p.6) corrects**: printed `GAS INDUCED FLAP REGRESSIVE` should read `CONTROL-SYSTEM-INDUCED FLAP REGRESSIVE`. No equations. |
| 17 | 3 | Introduction continues (component-type modelling as the alternative; Gen Hel / UH-60A / GE T700 selection). Acknowledgements paragraph (R. T. N. Chen, R. McFarland, D. Gilmore of GE). **ENGINE AND FUEL CONTROL SYSTEM** section begins — physical description of the T700-GE-700. No equations. |
| 18 | 4 | **Figure 2**, full page: *General Electric T700-GE-700 engine. (From ref. 6.)* Two cutaway drawings; lower one labels COMPRESSOR STATOR, COMPRESSOR, COMBUSTOR, POWER TURBINE, GAS GENERATOR, POWER TURBINE SHAFT, COMPRESSOR ANTI-ICING AND START BLEED MANIFOLD. Carries the NASA "ORIGINAL PAGE IS OF POOR QUALITY" stamp. No equations. |
| 19 | 5 | Fuel control system description — HMU (vane pump + mechanical cams: acceleration, deceleration, topping, idle schedules vs inlet temperature and NG; collective feed-forward; variable geometry) and ECU (isochronous NP governing, overtemperature protection, load sharing; torque motor adjusts HMU demand **downward** so electrical failure gives maximum power). Drive train and overrunning clutch. **ENGINE MODEL** section begins — general free-turbine turboshaft description, control volumes, inputs. No equations. |
| 20 | 6 | The **Figure 3 station-definition paragraph** (transcribed in §3 below). Then **Development History** (NASA Lewis hybrid model ref. 7 as the basis; simplifications; CSMP; then the Ames FORTRAN real-time version, which *added* a station-4.9 total-pressure function, a station-4.1 nonadiabatic (heat-sink) model, and a power-turbine damping factor). Then **Model Component Equations** subsection begins. No equations. |
| 21 | 7 | **Figure 3**, full page: *Primary components of a small turboshaft engine.* Described in full in §3. No equations. |
| 22 | 8 | **Eqs. (1)–(11)** — inlet/ambient, station 2 enthalpy, station 3 static pressure, corrected speed, compressor mass flow and temperature rise. |
| 23 | 9 | **Eqs. (12)–(24)** — bleed fractions and bleed flows, diffuser flow, combustor flow / FAR / efficiency / enthalpy rise, station 4.1 temperature and enthalpy including the heat-sink call. |
| 24 | 10 | **Eqs. (25)–(38)** — gas generator turbine, station 4.4 → 4.5 mixing, power turbine, station 4.9, exhaust. |
| 25 | 11 | **Eqs. (39)–(49)** — torques, the three control-volume pressure integrals, the two rotor-speed integrals, then **Heat-Sink Model** subsection begins with Eqs. (48)–(49). |

**Equation numbering is consecutive with no gaps**: 1–11 on p.22, 12–24 on p.23, 25–38 on
p.24, 39–49 on p.25. **49 numbered equations in this range.** The engine model continues
past this range: Eqs. (50)–(53) complete the heat-sink model on pdf p.26, and the
REAL-TIME IMPLEMENTATION section starts on that same page.

---

## 2. Every numbered equation, pp.22–25

Notation used below: `·` is the report's printed multiplication dot; `√` over a single
symbol is the report's radical; `∫ … dt` is the report's integral sign with `dt` at the
right. All station digits shown full-size are full-size on the page.

### 2.1 Inlet and station 2 — pdf p.22

Preceding text: *"Pressure and temperature stagnation effects have been found to roughly
offset the aircraft inlet losses over most of the aircraft flight envelope. Total pressure
P2 and temperature T2 are therefore estimated to be equal to the static ambient
conditions."*

**(1)** `[pdf p.22]`

```
P2 = P1 = P_amb
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| P2 | station 2 total pressure, lb_f/in² | yes |
| P1 | station 1 total pressure, lb_f/in² | yes |
| P_amb | ambient pressure, lb_f/in² | yes (model input) |

**(2)** `[pdf p.22]`

```
T2 = T1 = T_amb
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| T2 | station 2 temperature, deg R | yes |
| T1 | station 1 temperature, deg R | yes |
| T_amb | ambient temperature, deg R | yes (model input) |

Preceding text for (3): *"Station 2 enthalpy is calculated from station 2 temperature."*

**(3)** `[pdf p.22]`

```
H2 = K_H2 · T2
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| H2 | station 2 enthalpy, Btu/lb_m | yes |
| K_H2 | station 2 enthalpy coefficient, Btu/lb_m·deg R (0.239, Table A.1 row 8) | yes |

Note: single-coefficient fit, **no intercept** — unlike H3 (Eq. 11) and H41 (Eq. 24),
which are two-constant fits. Do not add an intercept to H2.

Preceding text for (4): *"Static pressure at station 3 is represented as a linear function
of the total pressure."*

**(4)** `[pdf p.22]`

```
P_s3 = K_Ps3 · P3
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| P_s3 | station 3 static pressure (a fuel-control-system parameter), lb_f/in² | yes |
| K_Ps3 | fraction of station 3 total pressure used to represent static pressure (0.956) | yes |
| P3 | station 3 total pressure, lb_f/in² — **a model state** | yes |

Typesetting: the report prints `P` with a small lowered `s` then a full-size `3` —
`P_s3`, "static pressure at station 3". Same pattern as `P_s9`.

Preceding text for (5)–(6): *"The gas-generator turbine speed is corrected by the square
root of a nondimensional temperature parameter."*

**(5)** `[pdf p.22]`

```
        T2
θ_2 = ------
       T_std
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| θ_2 | ratio of inlet temperature to standard-day temperature | yes |
| T_std | standard atmosphere temperature at sea level, deg R | yes — **value not in Table A.1** (see open questions) |

**(6)** `[pdf p.22]`

```
        NG
NG_c = -----
       √θ_2
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| NG_c | corrected gas generator speed, rpm | yes |
| NG | gas generator / compressor speed, rpm — **a model state** | yes |

### 2.2 Compressor — pdf p.22

Preceding text: *"…corrected compressor airflow is determined from a performance map
which is a function of corrected station-3 static pressure and corrected gas generator
speed only. This assumes that for any instant in time, the variable-geometry stator vanes
are in the nominal design position for a given operating point of the compressor."*

**(7)** `[pdf p.22]`

```
             ( P_s3        )
WA2_c = f₁ ( ------ , NG_c )
             (  P2         )
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| WA2_c | corrected station 2 mass flow of atmospheric gas, lb_m/sec | yes |
| f₁ | compressor mass-flow map, Fig. A1, pdf p.56 | **NO — see §7** |

**(8)** `[pdf p.22]`

```
        P2
δ_2 = ------
       P_std
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| δ_2 | ratio of inlet pressure to sea-level pressure | yes |
| P_std | standard-day pressure, lb_f/in² | yes — **value not in Table A.1** (see open questions) |

**(9)** `[pdf p.22]`

```
                   δ_2
WA2 = WA2_c · -----------
                  √θ_2
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| WA2 | station 2 mass flow of atmospheric gas, lb_m/sec | yes |

Preceding text for (10)–(11): *"Temperature and enthalpy changes from stations 2 to 3 are
determined as a function of the change of pressure across the compressor component."*

**(10)** `[pdf p.22]`

```
                ( P_s3 )
T3 = T2 · f₂ ( ------ )
                (  P2  )
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| T3 | station 3 temperature, deg R | yes |
| f₂ | compressor temperature map, Fig. A2, pdf p.57 | **NO — see §7** |

Note `f₂` is a **temperature ratio multiplier** on T2, not an additive rise.

**(11)** `[pdf p.22]`

```
H3 = K_H3₁ T3 + K_H3₂
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| H3 | station 3 enthalpy, Btu/lb_m | yes |
| K_H3₁ | station 3 enthalpy coefficient, Btu/lb_m·deg R (0.2496) | yes |
| K_H3₂ | station 3 enthalpy coefficient, Btu/lb_m (**−8.4**, negative) | yes |

**Sub-subscript trap.** `K_H3₁` and `K_H3₂` are two different constants with different
units, slope and intercept of one linear fit. There is no plain `K_H3`.

### 2.3 Bleeds — pdf p.23

Preceding text: *"The nonlinear fractions of bleed flow are also determined from function
maps. Actual bleed flows are then determined from the bleed fractions. Flow extracted at
station 2.4 acts to pressurize the seals, preventing oil loss and keeping hot gases, dust,
and moisture out of the oil sumps. In addition, the flow is used to pressurize the
power-turbine balance piston, which provides a forward force on the power turbine to
alleviate some of the thrust load on the turbine bearing. The diffuser discharge air,
WA3_bl, is used to cool the combustion liner and the gas-generator-turbine blades and
shrouds. The fraction of the bleed flow used to cool the gas generator is reintroduced to
the gas path downstream of the gas generator turbine. The turbine blades are cooled
internally and the bleed gas is then released through a shower-head-type series of
nozzles."*

**(12)** `[pdf p.23]`

```
B_1 = f₃ (NG_c)
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| B_1 | seal-pressurization bleed fraction | yes |
| f₃ | Fig. A3, pdf p.58 | **NO — see §7** |

**(13)** `[pdf p.23]`

```
B_2 = f₄ (WA2_c)
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| B_2 | power-turbine-balance bleed fraction | yes |
| f₄ | Fig. A4, pdf p.59 | **NO — see §7** |

**(14)** `[pdf p.23]`

```
B_3 = f₅ (WA2_c)
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| B_3 | compressor-diffuser bleed fraction | yes |
| f₅ | Fig. A5, pdf p.60 | **NO — see §7** |

**(15)** `[pdf p.23]`

```
WA24_bl = WA2(B_1 + B_2)
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| WA24_bl | station 2.4 bleed flow of atmospheric gas, lb_m/sec | yes |

**(16)** `[pdf p.23]`

```
WA3_bl = WA2(B_3 + K_b3)
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| WA3_bl | diffuser bleed discharge flow, lb_m/sec | yes |
| K_b3 | station 3 bleed-flow coefficient (0.0025, constant) | yes |

Note both bleed flows are scaled by **WA2**, the compressor-inlet flow — not by the local
station flow. `K_b3` is an additive *constant* bleed fraction on top of the mapped `B_3`.

Preceding text for (17): *"The fraction of mass flow at the diffuser is determined by
subtracting the compressor stage-4 bleed flow from flow entering the compressor."*

**(17)** `[pdf p.23]`

```
WA3 = WA2 − WA24_bl
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| WA3 | station 3 mass flow of atmospheric gas, lb_m/sec | yes |

### 2.4 Combustor and station 4.1 — pdf p.23

Preceding text: *"The combustor mass flow and efficiency are nonlinear functions of inlet
and outlet pressures, inlet temperature, and the fuel-to-air ratio. Change of enthalpy
across the combustor is also a function of the heating value of the fuel."*

**(18)** `[pdf p.23]`

```
           ______________
          / P3(P3 − P41)
WA31 =   / ---------------
       \/    K_dpb · T3
```

The radical covers the whole fraction.

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| WA31 | station 3.1 (combustor) mass flow of atmospheric gas, lb_m/sec | yes |
| P41 | station 4.1 total pressure, lb_f/in² — **a model state** | yes |
| K_dpb | combustor pressure drop coefficient, lb_f²·sec²/lb_m²·in⁴·deg R (0.03045) | yes |

**(19)** `[pdf p.23]`

```
        W_f
FAR = ------
       WA31
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| FAR | ratio of fuel to atmospheric gas in combustor | yes |
| W_f | fuel flow, lb_m/sec — **the single control input** | yes |

**(20)** `[pdf p.23]`

```
η = f₆(FAR)
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| η | combustor efficiency | yes |
| f₆ | Fig. A6, pdf p.61 | **NO — see §7** |

**(21)** `[pdf p.23]`

```
          H3 + η · FAR · HVF
H41_ns = --------------------
               1 + FAR
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| H41_ns | station 4.1 enthalpy **not** including heat-sink effects, Btu/lb_m | yes |
| HVF | heating value of fuel, Btu/lb_m (18300.0) | yes |

The `1 + FAR` denominator is the mass-weighting of the fuel addition. `_ns` is a small
lowered true subscript.

Preceding text for (22): *"The value of the temperature at station 4.1 without heat-sink
dynamics is a function of enthalpy."*

**(22)** `[pdf p.23]`

```
T41_ns = K_T41₁ H41_ns + K_T41₂
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| T41_ns | station 4.1 temperature **not** including heat-sink effects, deg R | yes |
| K_T41₁ | station 4.1 temperature coefficient, lb_m·deg R/Btu (3.322) | yes |
| K_T41₂ | station 4.1 temperature coefficient, deg R (288.7) | yes |

Preceding text for (23): *"The temperature at station 4.1 is expressed as a transfer
function with slowly varying coefficients. (See the following subsection.) If no heat-sink
representation is used, T41 = T41_ns."*

**(23)** `[pdf p.23]`

```
T41 = T41_ns · f_s(T41, T41_ns, W41, NG_c)
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| T41 | station 4.1 temperature, deg R — **the 6th state in the heat-sink variant** | yes |
| W41 | station 4.1 mass flow of combustion gases, lb_m/sec | yes |
| f_s | the heat-sink transfer function | **NO — see §7** |

**Read at 600 dpi and pixel-zoomed: the subscript is a single lowercase italic `s`, not
`hs`.** `f_s` (Eq. 23) and `f_hs` (Eq. 52, pdf p.26, Fig. A11) are **two different
objects** and must not be conflated. `f_s` is not an Appendix A figure — it is the
transfer function built analytically from Eqs. (48)–(53) on pp.25–26; `f_hs` is one lookup
*inside* it.

**Eq. 23 is implicit**: `T41` appears on both sides. It is also coupled to Eq. 28 through
`W41` (which needs θ_41 from Eq. 25, which needs T41). See open questions.

**(24)** `[pdf p.23]`

```
H41 = K_H41₁ T41 + K_H41₂
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| H41 | station 4.1 enthalpy, Btu/lb_m | yes |
| K_H41₁ | station 4.1 enthalpy coefficient, Btu/lb_m·deg R (0.3010) | yes |
| K_H41₂ | station 4.1 enthalpy coefficient, Btu/lb_m (**−86.905**, negative) | yes |

### 2.5 Gas generator turbine — pdf p.24

Preceding text: *"As described in reference 7, a critical velocity parameter is used to
calculate an enthalpy change parameter for the gas generator turbine as a function of
pressure ratio only. Actual enthalpy drop is calculated by multiplying this parameter by
the squared critical velocity ratio. Exit enthalpy is then determined."*

**(25)** `[pdf p.24]`

```
θ_41 = K_TH41₁ T41 + K_TH41₂
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| θ_41 | station 4.1 squared-critical-velocity ratio | yes |
| K_TH41₁ | station 4.1 velocity ratio coefficient, 1/deg R (0.0018326) | yes |
| K_TH41₂ | station 4.1 velocity ratio coefficient, nondimensional (0.0856) | yes |

**(26)** `[pdf p.24]`

```
                    ( P45 )
ΔH_GT = θ_41 · f₇ ( ----- )
                    ( P41 )
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| ΔH_GT | gas-generator-turbine enthalpy drop, Btu/lb_m | yes |
| P45 | station 4.5 total pressure, lb_f/in² — **a model state** | yes |
| f₇ | Fig. A7, pdf p.62 | **NO — see §7** |

The argument is the **expansion** ratio P45/P41 (< 1), not P41/P45.

**(27)** `[pdf p.24]`

```
H44 = H41 − ΔH_GT
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| H44 | station 4.4 enthalpy, Btu/lb_m | yes |

Preceding text for (28): *"Over the normal operating range of the engine, a choked nozzle
equation is adequate to calculate the mass flow."*

**(28)** `[pdf p.24]`

```
                 P41
W41 = K_WGT · ---------
                √θ_41
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| K_WGT | station 4.1 flow coefficient, lb_m·in²/lb_f·sec (0.0876) | yes |

### 2.6 Station 4.5 mixing — pdf p.24

Preceding text: *"At station 4.5, gases from the station 4.4 and cooling-bleed flow from
the compressor are mixed before passing through the power turbine. The enthalpy of the
mixed gases is proportional to enthalpy at station 4.4. Temperature is then determined
from enthalpy."*

**(29)** `[pdf p.24]`

```
H45 = K_H45 H44
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| H45 | station 4.5 enthalpy, Btu/lb_m | yes |
| K_H45 | fraction of station 4.4 enthalpy used to represent station 4.5 enthalpy (0.9623) | yes |

The cooling-bleed dilution enters the **energy** balance only as this single constant
fraction. It enters the **mass** balance separately, in Eq. 44.

**(30)** `[pdf p.24]`

```
T45 = K_T45₁ H45 + K_T45₂
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| T45 | station 4.5 temperature (power turbine **inlet**), deg R | yes |
| K_T45₁ | lb_m·deg R/Btu (3.519) | yes |
| K_T45₂ | deg R (179.1) | yes |

### 2.7 Power turbine — pdf p.24

Preceding text: *"A power-turbine enthalpy drop parameter is determined as a function of
pressure ratio across the turbine, as described for the gas-generator turbine. In
addition, mass flow is determined as a function of pressure ratio only."*

**(31)** `[pdf p.24]`

```
θ_45 = K_TH45₁ T45 + K_TH45₂
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| θ_45 | station 4.5 squared-critical-velocity ratio | yes |
| K_TH45₁ | 1/deg R (0.0018326) | yes |
| K_TH45₂ | nondimensional (0.0856) | yes |

Table A.1 gives `K_TH45₁ = K_TH41₁` and `K_TH45₂ = K_TH41₂` numerically. They are still
four separate table rows; keep them separate.

**(32)** `[pdf p.24]`

```
                    ( P49 )
ΔH_PT = θ_45 · f₈ ( ----- )
                    ( P45 )
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| ΔH_PT | power turbine enthalpy drop, Btu/lb_m | yes |
| P49 | station 4.9 total pressure, lb_f/in² | yes |
| f₈ | Fig. A8, pdf p.63 | **NO — see §7** |

**(33)** `[pdf p.24]`

```
            ( P_s9 )
W45_c = f₉ ( ------ )
            ( P45  )
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| W45_c | corrected station 4.5 mass flow of combustion gases, lb_m/sec | yes |
| P_s9 | station 9 static pressure, lb_f/in² | yes |
| f₉ | Fig. A9, pdf p.64 | **NO — see §7** |

The argument is `P_s9/P45` — the **exhaust static** over the PT inlet total, *not*
`P49/P45`. Two different pressure ratios appear on this page (Eq. 32 uses P49/P45); do not
interchange them.

**(34)** `[pdf p.24]`

```
                   P45
W45 = W45_c · ----------
                 √θ_45
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| W45 | station 4.5 mass flow of combustion gases, lb_m/sec | yes |

Note the de-correction here uses **P45**, not δ — this is an absolute-pressure form, unlike
Eq. 9 which uses δ_2.

Preceding text for (35)–(36): *"Station 4.9 enthalpy and temperature are determined from
the change in enthalpy across the power turbine."*

**(35)** `[pdf p.24]`

```
H49 = H45 − ΔH_PT
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| H49 | station 4.9 enthalpy, Btu/lb_m | yes |

Scan note: the `9` of `H49` on the left-hand side is broken in the scan (prints as an open
`Ɔ`-like glyph at 600 dpi). It is unambiguously `H49` — Eq. 36 on the next line reads
`T49 = K_T49₁ H49 + K_T49₂` with a clean `9`, and station 4.9 is the only candidate.

**(36)** `[pdf p.24]`

```
T49 = K_T49₁ H49 + K_T49₂
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| T49 | station 4.9 temperature, deg R | yes |
| K_T49₁ | lb_m·deg R/Btu (3.516) | yes |
| K_T49₂ | deg R (172.3) | yes |

Preceding text for (37)–(38): *"Station 4.9 pressure is determined from a nonlinear
function map of gas generator speed (provided by GE)."*

**(37)** `[pdf p.24]`

```
P_s9 = P_amb
```

**(38)** `[pdf p.24]`

```
P49 = P_s9 · f₁₀(NG_c)
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| f₁₀ | Fig. A10, pdf p.65 | **NO — see §7** |

`f₁₀` is a **multiplier on P_s9 to obtain P49**, so `f₁₀ = P49/P_s9 > 1`. Figure A10's
ordinate is printed `PS9/P49` — inverted with respect to this equation. Resolved by
decision (open question #5): digitize the figure as printed and take the reciprocal at the
point of use.

### 2.8 Torques — pdf p.25

Preceding text: *"Values of torque which are output or required by the engine are
determined from the changes in energy across the compressor and turbines. Effects of the
compressor interstage bleed flows are accounted for with an empirically determined
function. A damping factor based on change of power turbine speed is an additional term in
the power turbine torque equation."*

**(39)** `[pdf p.25]`

```
                   60     1
Q_C = 778.12 · ------ · ---- {WA2(K_QC₁ H3 − H2) + WA3 · K_QC₂ H3}
                   2π     NG
```

Braces `{ }` are the report's own. Verified at 600 dpi: the sign inside the first
parenthesis is a **minus**; `K_QC₁` multiplies `H3` only, and `H2` is unscaled.

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| Q_C | torque required by compressor, ft·lb_f | yes |
| K_QC₁ | empirically-determined compressor torque coefficient (0.71) | yes |
| K_QC₂ | empirically-determined compressor torque coefficient (0.29) | yes |
| 778.12 | (no symbol; printed as a bare number) — ft·lb_f per Btu | **NO — see §7** |

`K_QC₁ + K_QC₂ = 1.00` exactly, in the printed Table A.1 values. That is a property of the
numbers, not something the report states.

**(40)** `[pdf p.25]`

```
                    60     1
Q_GT = 778.12 · ------ · ---- · W41 ΔH_GT
                    2π     NG
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| Q_GT | torque output of gas generator, ft·lb_f | yes |

**(41)** `[pdf p.25]`

```
                    60     1                             2π
Q_PT = 778.12 · ------ · ---- · W45 ΔH_PT − K_damp · ----- · (NP − NP_des)
                    2π     NP                            60
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| Q_PT | torque output of power turbine, ft·lb_f | yes |
| K_damp | power turbine speed damping coefficient, ft·lb_f·sec/rad (0.06854) | yes |
| NP | power turbine / output shaft speed, rpm — **a model state** | yes |
| NP_des | design power turbine speed, rpm (20900.0) | yes |

Note the two conversion factors run **opposite ways**: `60/2π` converts rpm to rad/sec in
the denominator of the power term; `2π/60` converts the rpm error to rad/sec in the damping
term. Both are printed. The damping term is the item the Development History paragraph
(p.20) says was **added at Ames** to match low-power transient response — it is not in the
NASA Lewis hybrid model.

### 2.9 Control-volume pressure dynamics — pdf p.25

Preceding text: *"Conservation of mass is applied to determine intervolume pressure
dynamics for stations 3, 4.1, and 4.5."*

**(42)** `[pdf p.25]`

```
P3 = K_V3 ∫ T3(WA3 − WA3_bl − WA31) dt
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| K_V3 | station 3 volume coefficient, lb_f/in²·lb_m·deg R (0.97) | yes |

**(43)** `[pdf p.25]`

```
P41 = K_V41 ∫ T41(WA31 − W_f − W41) dt
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| K_V41 | station 4.1 volume coefficient, lb_f/in²·lb_m·deg R (6.17) | yes |

**⚠ SIGN, verified and confirmed as printed.** The fuel flow `W_f` enters with a **minus**.
Read at 600 dpi and pixel-zoomed, and compared against a known `+` glyph elsewhere on the
same page (Eq. 44 `+ B_3 K_bl WA2`, Eq. 46 `J_PT + J_load`), which show a clear vertical
stroke that these two signs do not have. Both operators in Eq. 43 are minus. A mass balance
on the station-4.1 volume needs `WA31 + W_f − W41`. **Transcribed as printed; logged as an
open question — do not silently correct it in code.**

**(44)** `[pdf p.25]`

```
P45 = K_V45 ∫ T45(W41 − W45 + B_3 K_bl WA2) dt
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| K_V45 | station 4.5 volume coefficient, lb_f/in²·lb_m·deg R (13.63) | yes |
| K_bl | fraction of diffuser bleed gas used to cool the gas generator turbine blades (0.7826) | yes |

**The returned cooling flow is `B_3 · K_bl · WA2`, not `K_bl · WA3_bl`.** Since
`WA3_bl = WA2(B_3 + K_b3)` (Eq. 16), the `K_b3` component of the station-3 bleed never
comes back into the gas path. That asymmetry is deliberate in the printed equations and
must be reproduced.

### 2.10 Rotor dynamics — pdf p.25

Preceding text: *"Turbine speeds are determined as a function of the externally-applied
torques by assuming conservation of angular momentum."*

**(45)** `[pdf p.25]`

```
       60    ⌠ Q_GT − Q_C
NG = ----- · │ ------------ dt
       2π    ⌡    J_GT
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| J_GT | moment of inertia of the rigid mass representing the gas generator, compressor and associated shafting, ft·lb_f·sec² (0.0445) | yes |

**(46)** `[pdf p.25]`

```
J = J_PT + J_load
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| J | moment of inertia of all mass rigidly attached to the engine output shaft, ft·lb_f·sec² | yes |
| J_PT | moment of inertia of the power turbine and output shaft, ft·lb_f·sec² (0.062) | yes |
| J_load | moment of inertia of the load mass, ft·lb_f·sec² | yes — **value not in Table A.1** (open question #6) |

**(47)** `[pdf p.25]`

```
       60    ⌠ Q_PT − Q_req
NP = ----- · │ -------------- dt
       2π    ⌡       J
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| Q_req | torque required by external load with respect to power turbine speed, ft·lb_f | yes — **an external input from Gen Hel, not computed here** |

### 2.11 Heat-Sink Model, first two equations — pdf p.25

Section heading **Heat-Sink Model** appears on p.25 between Eq. 47 and Eq. 48.

Preceding text: *"Early in the validation effort, the modeling of nonadiabatic processes
occurring in areas of large temperature gradients was found to be necessary for correct
transient response. There is a transient energy transfer which is caused by the metal mass
of the turbine absorbing heat from the hot gas. This is especially significant at station
4.1, immediately downstream of the combustor. Known as the heat-sink effect, this
phenomenon can be modeled as a lumped parameter system with the heat transfer
equations."*

**(48)** `[pdf p.25]`

```
        dT_m
c_pm M ------ = h A_m (T_gi − T_m)
         dt
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| c_pm | specific heat of metal (units not given in nomenclature) | yes |
| M | mass of metal that absorbs heat energy from gas | yes (Python: `m_metal`) |
| T_m | temperature of exposed metal | yes |
| h | convective heat transfer coefficient | yes (Python: `h_conv`) |
| A_m | area of exposed metal | yes |
| T_gi | temperature of gas entering the heat-sink control volume | yes |

`T_gi` is printed as `T` with a small lowered `g` carrying a further sub-subscript `i` —
the same two-level subscript typesetting as the `K…₁/₂` family. Likewise `T_go`.

**(49)** `[pdf p.25]`

```
                                   dT_m
c_pg W_g (T_gi − T_go) = c_pm M ------
                                    dt
```

| Symbol | Meaning | In `symbols.md`? |
|---|---|---|
| c_pg | specific heat of gas | yes |
| W_g | mass flow rate of gas into the heat-sink control volume | yes |
| T_go | temperature of gas leaving the heat-sink control volume | yes |

Closing text on p.25: *"The first equation relates the change in metal temperature to the
difference between metal and gas temperatures, and the second equation represents the
effect of this change of metal temperature on the gas."*

The heat-sink model **continues onto pdf p.26** with Eqs. (50)–(53), which cast (48)–(49)
as the lead/lag transfer function `T_go/T_gi = T41/T41_ns` that Eq. 23 calls `f_s`, define
the time constant `M c_pm / h A_m` via `TC_T41` and `W41`, and supply `T41_sgn = f_hs(NG_c)`
from Fig. A11. Those four equations are outside this transcription's page range and belong
to whoever covers pp.26–37.

**None of Eqs. (48)–(49) has a dot-over-variable derivative.** Both use explicit Leibniz
`dT_m/dt`. No overdot notation appears anywhere in pp.15–25; it first appears in the
real-time / linear-model sections beyond p.25.

---

## 3. Figure 3 — the station diagram, pdf p.21

**Caption as printed:** `Figure 3: Primary components of a small turboshaft engine.`

### What kind of figure it is

Figure 3 is a **line-art axial cross-section of the engine** (upper and lower halves about
a horizontal centreline, with rotating hardware filled solid black and stationary
flowpath/casing drawn in outline). It is **not** a block or flow-network diagram. It has
**no flow arrows other than WF**, **no bleed arrows, no bleed lines, and no control-volume
boxes**. The topology it defines is the *axial position* of each station along the gas
path, and nothing more. Everything about where bleeds are taken and returned comes from
the surrounding text (pdf p.20 and p.23), not from the figure.

Read at 300 dpi and re-read at 600 dpi in two overlapping crops.

### Station labels — the complete set, left to right

A row of numerals runs across the top of the figure, the leftmost prefixed by the word
`STATION`. Each numeral drops a straight vertical leader line into the gas path. In
printed order:

`STATION 1` · `2` · `2.4` · `3` · `3.1` · `4.1` · `4.4` · `4.5` · `4.9` · `9`

**Ten stations. There is no station 4 and no station 5 anywhere in the figure or the
report.** The gas-generator turbine exit is 4.4, the power turbine inlet is 4.5, the power
turbine exit is 4.9, the exhaust is 9. `T45` is therefore the power-turbine *inlet*
temperature and must never be mapped onto another engine's "T4.5".

### Where each leader lands in the gas path

| Station | Leader terminates at | Report's own words (pdf p.20 unless noted) |
|---|---|---|
| 1 | The engine inlet, at the front face of the inlet duct / particle-separator housing, upstream of everything | Figure 3 labels it; the report gives no verbal definition. Eqs. 1–2 set P1 = P_amb, T1 = T_amb |
| 2 | The compressor face, immediately aft of the curved inlet duct and the variable inlet guide vanes, at the first axial rotor stage | *"The six compressor stages and variable-geometry flow vanes are represented with the compressor component, between stations 2 and 3."* |
| 2.4 | Within the axial compressor blading, at the 4th stage | *"Bleed flow which is used for seal pressurization and power turbine balance is extracted from stage 4 of the compressor, which is represented by station 2.4."* |
| 3 | The diffuser passage immediately aft of the centrifugal impeller | *"Station 3 is the compressor diffuser, or outlet. Flow is bled at this point to cool the combustor and gas generator turbine."* |
| 3.1 | The front of the annular combustor liner (the combustor dome) | *"Station 3.1 is the combustor…"* |
| 4.1 | The combustor exit / gas-generator-turbine first-stage nozzle inlet | *"…and station 4.1 is the mixing volume representing the combustor outlet and gas-generator-turbine inlet."* |
| 4.4 | Immediately aft of the second gas-generator-turbine stage | *"Station 4.4 represents the thermodynamic state of output gases passing through this turbine, not including the effects of the cooling bleed flow."* |
| 4.5 | Between the gas-generator turbine and the power turbine, at the power turbine inlet | *"These effects are added at station 4.5, which represents the power turbine inlet."* |
| 4.9 | Immediately aft of the second power-turbine stage | *"Station 4.9 is the power turbine outlet…"* |
| 9 | The exhaust duct exit, at the rear of the engine | *"…and station 9 is the engine exhaust."* |

### Component labels

A second row of labels runs across the **bottom**, each with a leader up into the drawing:

`COMPRESSOR` (leader into the axial compressor blading) · `COMBUSTOR` (leader into the
combustor liner) · `GAS GENERATOR` (leader into the gas-generator turbine blading, between
the 4.1 and 4.4 leaders) · `POWER TURBINE` (leader into the power turbine blading, between
the 4.5 and 4.9 leaders).

Four components, exactly the four the text names: *"A typical free-turbine turboshaft
engine consists of four major components: the compressor, the combustor, the gas
generator, and the power turbine"* [pdf p.19].

### The one flow annotation

`WF` labels a short **downward arrow** entering the combustor liner from above, positioned
between the `3.1` and `4.1` leaders. This is the only mass-flow arrow in the figure and
the only external input drawn. It is fuel flow, `W_f` in the equations (the figure sets
`WF` in full-size capitals; the nomenclature and equations use `W` with a small lowered
`f`).

### Shaft arrangement as drawn

Two concentric shafts are drawn solid black along the centreline. The gas-generator shaft
carries the compressor blading and the gas-generator turbine blading. The power-turbine
shaft carries the power turbine blading and runs **forward, coaxially, through the whole
engine and out the front** — visible as the inner solid bar extending left past the
compressor to the output flange at the extreme left of the drawing. This matches the text:
*"The power turbine, which is uncooled, has a coaxial driveshaft which extends forward
through the front of the engine where it is connected to the output shaft assembly"*
[pdf p.17].

### Flow topology the model actually uses (assembled from the equations, not the figure)

The figure shows only the axial layout. The complete in/out accounting, from Eqs. 15–17,
42–44 and the p.23 paragraph:

```
  ambient ──► [1] ──► [2]  WA2  ──► compressor
                              │
            station 2.4 bleed │  WA24_bl = WA2(B_1 + B_2)      Eq.15
            LEAVES THE MODEL  ▼   (seal pressurization + power turbine balance piston)
                              ·
                       WA3 = WA2 − WA24_bl                     Eq.17
                              │
                              ▼
                            [3]  ── control volume, state P3 ── Eq.42
                              │
            station 3 bleed   │  WA3_bl = WA2(B_3 + K_b3)      Eq.16
            ──────────────────┤   (cools combustion liner + GG turbine blades/shrouds)
                              │        of which B_3·K_bl·WA2 is REINTRODUCED at [4.5]
                              │        and the K_b3 part is not
                              ▼
                            [3.1] WA31 = f(P3, P3−P41, T3)     Eq.18
                              │   + W_f (fuel)                 Eqs.19-21
                              ▼
                            [4.1] ── control volume, state P41 ── Eq.43
                              │   W41 = K_WGT · P41/√θ_41      Eq.28
                              ▼
                        gas generator turbine  (ΔH_GT, Q_GT)
                              │
                              ▼
                            [4.4] (enthalpy only, H44)
                              │
                              ├──◄── returned cooling bleed B_3·K_bl·WA2
                              ▼
                            [4.5] ── control volume, state P45 ── Eq.44
                              │   W45 = W45_c · P45/√θ_45      Eq.34
                              ▼
                        power turbine  (ΔH_PT, Q_PT)
                              │
                              ▼
                            [4.9] (H49, T49, P49 = P_s9·f₁₀)
                              │
                              ▼
                            [9]   P_s9 = P_amb
```

Stations that carry **no full state**: 3.1 is a mass flow only (`WA31`); 4.4 is an enthalpy
only (`H44`); 9 is a static pressure only (`P_s9`); 1 is `P1`, `T1` set to ambient.

---

## 4. The control-volume formulation — how P3, P41, P45 are formed

This is the mechanism that makes the model a **5-state real-time model** rather than an
iterative component-matching model, and it is contained entirely in Eqs. (42)–(44) plus the
three constants `K_V3`, `K_V41`, `K_V45`.

### What the report says

One sentence, pdf p.25, immediately above Eq. 42:

> *"Conservation of mass is applied to determine intervolume pressure dynamics for stations
> 3, 4.1, and 4.5."*

And the framing on pdf p.19:

> *"The four major components are separated by fluid mixing volumes, each of which is
> associated with flow passages within the engine where thermodynamic states of the gas are
> quantifiable. States of the gas in each control volume are expressed in terms of pressure,
> temperature, and mass flow. … Conservation of mass is used to determine the values of mass
> flows into and out of each control volume."*

That is the report's entire prose justification. It does **not** derive the equations and
does **not** state the ideal-gas relation behind them.

### The three equations, side by side

| State | Equation | Volume temperature | Net mass flow into the volume |
|---|---|---|---|
| P3 | (42) `P3 = K_V3 ∫ T3(WA3 − WA3_bl − WA31) dt` | T3 | in: `WA3`; out: `WA3_bl` (diffuser bleed), `WA31` (to combustor) |
| P41 | (43) `P41 = K_V41 ∫ T41(WA31 − W_f − W41) dt` | T41 | in: `WA31`; out: `W_f` (**sic — printed as an outflow**), `W41` (to GG turbine) |
| P45 | (44) `P45 = K_V45 ∫ T45(W41 − W45 + B_3 K_bl WA2) dt` | T45 | in: `W41`, returned cooling bleed `B_3·K_bl·WA2`; out: `W45` (to power turbine) |

All three have the identical shape

```
P<n> = K_V<n> ∫ T<n> · (Σ flows) dt
```

— the volume's own temperature multiplies the *net* mass flow imbalance, and the integral
of that product, scaled by one constant, **is** the total pressure. Pressure is a pure
integrator state: nothing solves for it, nothing iterates on it.

### What `K_V3`, `K_V41`, `K_V45` do

The nomenclature calls each of them a *"station N volume coefficient"* [pdf p.11] and gives
units `lb_f/in²·lb_m·deg R`. Table A.1 [pdf p.55] gives:

| Constant | Value | Station |
|---|---|---|
| `K_V3` | 0.97 | compressor diffuser |
| `K_V41` | 6.17 | combustor outlet / GG turbine inlet mixing volume |
| `K_V45` | 13.63 | power turbine inlet |

**Units notation.** The printed string `lb_f/in²·lb_m·deg R` is ambiguous about what the
solidus covers. Eq. 42 settles it: `[psia] = K_V · [deg R]·[lb_m/sec]·[sec]`, so
`K_V` must be **lb_f/in² per (lb_m · deg R)** — psia per pound-mass per degree Rankine.
Read the printed units as `lb_f / (in² · lb_m · deg R)`. This is a reading of the
notation, forced by the equation; the report never writes the parentheses.

Each `K_V` is therefore inversely proportional to its control volume's physical size — a
larger volume gives a smaller pressure change per unit of mass imbalance. The ordering of
the printed values is consistent with that (`K_V3` smallest at the diffuser, `K_V45`
largest at the smallest volume), but **the report does not state a volume anywhere**, does
not give the gas constant used, and does not say how the three numbers were obtained. They
are single empirical constants, to be used as printed. Do not attempt to back out a
physical volume from them.

### Why this makes the model real-time

Because P3, P41 and P45 are integrator outputs rather than the result of an inner
iteration, one pass through Eqs. (1)–(49) evaluates the whole engine with **no inner loop**:
the states {NG, NP, P3, P41, P45} come in, the algebra runs forward station by station, and
five derivatives come out. That is the property the whole Volume Dynamics Approximation
section on pp.26–37 exists to defend — the report's argument there is that the volume
dynamics can be *approximated* (quasi-steady) without losing the bandwidth that a
blade-element rotor simulation needs. The nonlinear model of Eqs. 42–44 as printed is the
5-DOF baseline against which that approximation is checked.

**One inner iteration does survive in the printed equations, and it is not in the
pressures**: Eq. 23 is implicit in T41 and coupled to Eq. 28 through `W41` and `θ_41`. See
§7.

---

## 5. Rotor dynamics — NG and NP

Both rotors are single-degree-of-freedom rigid inertias driven by a torque difference,
written as integrals with the rpm↔rad/sec conversion `60/2π` outside the integral.

### 5.1 Gas generator spool, NG — Eq. (45)

```
       60    ⌠ Q_GT − Q_C
NG = ----- · │ ------------ dt
       2π    ⌡    J_GT
```

- **Driving torque**: `Q_GT`, the gas-generator turbine output, Eq. 40 —
  `778.12 · (60/2π) · (1/NG) · W41 ΔH_GT`.
- **Resisting torque**: `Q_C`, the compressor, Eq. 39 —
  `778.12 · (60/2π) · (1/NG) {WA2(K_QC₁ H3 − H2) + WA3 · K_QC₂ H3}`.
  The two-term form with `K_QC₁`/`K_QC₂` is the report's way of accounting for the
  interstage bleed: the `WA2` term carries the full inlet flow at a reduced enthalpy
  weighting, the `WA3` term adds the post-bleed flow at the complementary weighting.
  Prose confirms: *"Effects of the compressor interstage bleed flows are accounted for with
  an empirically determined function."*
- **Inertia**: `J_GT` = 0.0445 ft·lb_f·sec² — *"the rigid mass which represents the gas
  generator, compressor, and the associated shafting"*.
- **No external load, no damping term.** The gas generator spool sees only turbine minus
  compressor.

### 5.2 Power turbine / output shaft, NP — Eqs. (46)–(47)

```
J = J_PT + J_load

       60    ⌠ Q_PT − Q_req
NP = ----- · │ -------------- dt
       2π    ⌡       J
```

- **Driving torque**: `Q_PT`, Eq. 41, which is the power-turbine work term **minus a speed
  damping term** `K_damp · (2π/60) · (NP − NP_des)`. The damping term is proportional to
  the *displacement of NP from design speed*, not to `dNP/dt`, despite the prose calling it
  *"a damping factor based on change of power turbine speed"*. Reproduce the equation, not
  the prose.
- **Resisting torque**: `Q_req`, *"torque required by external load with respect to power
  turbine speed"*. Per pdf p.28, Eq. 56, `Q_req = Q_mr + Q_tr + Q_acc + Q_damp` — main
  rotor, tail rotor, accessories, gearbox damping — all of which come from the **Gen Hel
  UH-60A blade-element simulation**, an external program. `Q_req` is an input boundary
  condition to this model, not something it computes.
- **Inertia**: `J = J_PT + J_load`. `J_PT` = 0.062 ft·lb_f·sec² (Table A.1). **`J_load` has
  no value anywhere in the report** — like `Q_req` it is presumed to be a Gen Hel quantity
  (open question #6).
- The nomenclature also defines `J` on its own as *"moment of inertia of all mass rigidly
  attached to the engine output shaft"*, which is exactly what Eq. 46 constructs.

### 5.3 What is asymmetric between the two

| | NG (Eq. 45) | NP (Eq. 47) |
|---|---|---|
| Torque in | Q_GT (gas generator turbine) | Q_PT (power turbine) |
| Torque out | Q_C (compressor) | Q_req (external — Gen Hel) |
| Damping term | none | inside Q_PT, `−K_damp·(2π/60)·(NP − NP_des)` |
| Inertia | J_GT, a single constant | J = J_PT + J_load, the second term external |
| Fully self-contained? | yes | **no** — both Q_req and J_load are external |

The clutch described on pdf p.19 (*"power from each engine must also pass through a clutch
which disengages the engine from the gearbox if transmission speed exceeds that of the
engine shaft"*) is **not** represented in Eqs. 45–47 or anywhere else in pp.15–25.

---

## 6. Where each of `f₁`…`f₁₀`, `f_s`, `f_hs` is used

Complete wiring map. Every one of `f₁`–`f₁₀` and `f_hs` is called from exactly one
equation, with exactly one argument list.

| Function | Called by | Argument(s) | Produces | Appendix A figure |
|---|---|---|---|---|
| `f₁` | **Eq. 7**, p.22 | `P_s3/P2` **and** `NG_c` (two arguments — the only 2-D map) | `WA2_c`, corrected compressor mass flow, lb_m/sec | Fig. A1, pdf p.56 (11 speed lines) |
| `f₂` | **Eq. 10**, p.22 | `P_s3/P2` | a **multiplier on T2** giving T3 | Fig. A2, pdf p.57 |
| `f₃` | **Eq. 12**, p.23 | `NG_c` | `B_1`, seal-pressurization bleed fraction | Fig. A3, pdf p.58 |
| `f₄` | **Eq. 13**, p.23 | `WA2_c` | `B_2`, power-turbine-balance bleed fraction | Fig. A4, pdf p.59 |
| `f₅` | **Eq. 14**, p.23 | `WA2_c` | `B_3`, compressor-diffuser bleed fraction | Fig. A5, pdf p.60 |
| `f₆` | **Eq. 20**, p.23 | `FAR` | `η`, combustor efficiency | Fig. A6, pdf p.61 |
| `f₇` | **Eq. 26**, p.24 | `P45/P41` (expansion ratio, < 1) | enthalpy-drop **parameter**, multiplied by `θ_41` to give `ΔH_GT` | Fig. A7, pdf p.62 |
| `f₈` | **Eq. 32**, p.24 | `P49/P45` (expansion ratio, < 1) | enthalpy-drop **parameter**, multiplied by `θ_45` to give `ΔH_PT` | Fig. A8, pdf p.63 |
| `f₉` | **Eq. 33**, p.24 | `P_s9/P45` — **static exhaust over PT-inlet total**, not `P49/P45` | `W45_c`, corrected PT mass flow, lb_m/sec | Fig. A9, pdf p.64 |
| `f₁₀` | **Eq. 38**, p.24 | `NG_c` | a **multiplier on P_s9** giving `P49` (so `f₁₀ = P49/P_s9 > 1`) | Fig. A10, pdf p.65 — ordinate printed inverted, see open question #5 |
| `f_hs` | **Eq. 52, pdf p.26** — outside this page range | `NG_c` | `T41_sgn`, the station 4.1 heat-sink constant | Fig. A11, pdf p.66 |
| `f_s` | **Eq. 23**, p.23 | `T41`, `T41_ns`, `W41`, `NG_c` (four arguments) | a **multiplier on T41_ns** giving `T41` | **not a figure** — the transfer function assembled from Eqs. 48–53, pp.25–26 |

### Things to get right when wiring these in

1. **`f₂` and `f₁₀` are multipliers, not values.** `T3 = T2 · f₂(·)` and
   `P49 = P_s9 · f₁₀(·)`. Neither returns the station quantity directly.
2. **`f₇` and `f₈` return a parameter, not the enthalpy drop.** The actual drop is
   `θ · f(·)`, where θ is the squared-critical-velocity ratio from Eqs. 25 / 31.
3. **`f₄` and `f₅` share the same abscissa (`WA2_c`); `f₃` and `f₁₀` and `f_hs` share
   `NG_c`.** Three different quantities keyed on `NG_c` — do not merge their tables.
4. **`f₉`'s argument is `P_s9/P45`, `f₈`'s is `P49/P45`.** Two different numerators on the
   same page, one static and one total. The distinction is load-bearing.
5. **`f₁` is the only two-dimensional map** and the only one taking a speed *and* a
   pressure ratio.
6. **`f_s` and `f_hs` are different objects.** Eq. 23 prints a single lowercase `s`
   (confirmed by pixel-level zoom at 600 dpi). `f_hs` is one lookup nested inside `f_s`.
7. **None of `f₁`…`f₁₀`, `f_s`, `f_hs` appears in the nomenclature.** The 106-entry symbol
   table has no `f` entries at all. Their only definitions are the equations that call them
   and the Appendix A captions.

---

## 7. Symbols used in Eqs. 1–49 that are NOT in `symbols.md`

The nomenclature (pdf pp.9–13, 106 entries) covers every state, station quantity, constant
and Greek symbol appearing in Eqs. 1–49 — with exactly these exceptions:

| Not in the symbol table | Where it appears | What it is | Consequence |
|---|---|---|---|
| `f₁`, `f₂`, `f₃`, `f₄`, `f₅`, `f₆`, `f₇`, `f₈`, `f₉`, `f₁₀` | Eqs. 7, 10, 12, 13, 14, 20, 26, 32, 33, 38 | the ten Appendix A function maps | Their meaning is defined only by the calling equation and the Appendix A caption. No unit is stated for any of them. |
| `f_s` | Eq. 23 | the station 4.1 heat-sink transfer function | Not an Appendix A figure. Defined by Eqs. 48–53, pp.25–26. |
| `f_hs` | Eq. 52 (pdf p.26, outside this range) | the `T41_sgn` lookup, Fig. A11 | Listed here for completeness. |
| `778.12` | Eqs. 39, 40, 41 — printed as a bare number, three times | ft·lb_f per Btu (mechanical equivalent of heat), converting the Btu/lb_m enthalpies to ft·lb_f torque | Not a nomenclature entry and **not a Table A.1 row**. Cite the equation page, not a handbook, and use the report's 778.12 exactly. |
| `60`, `2π`, `1`, `1 + FAR` etc. | throughout | unit conversions and structural constants | Not symbols; no action. |

Two symbols **are** in the table but have **no value anywhere in the report**:
`T_std` (Eq. 5) and `P_std` (Eq. 8). Neither is in Table A.1. See open questions.

---

## 8. Cross-checks and observations recorded, not acted on

These are properties of the printed equations. They are recorded so nobody rediscovers
them; none of them licenses changing a transcription.

- **Eq. 43's `− W_f` is the one real anomaly** in pp.22–25. Verified as printed at pixel
  level against a known `+` on the same page. Logged as an open question.
- **`K_QC₁ + K_QC₂ = 0.71 + 0.29 = 1.00`** exactly, from Table A.1. Suggestive of a
  weighted split, but the report never says so.
- **`K_TH45₁ = K_TH41₁ = 0.0018326` and `K_TH45₂ = K_TH41₂ = 0.0856`** numerically. Still
  four separate constants in four separate rows; keep them separate in code so a later
  correction to one does not silently move the other.
- **`H2` (Eq. 3) is a one-constant fit; `H3` (Eq. 11) and `H41` (Eq. 24) are two-constant
  fits with negative intercepts** (−8.4 and −86.905 Btu/lb_m). The three enthalpies do not
  share a form.
- **The `K_b3` fraction of the station-3 bleed never returns to the gas path.** Eq. 16
  takes `WA2(B_3 + K_b3)` out; Eq. 44 puts `B_3 K_bl WA2` back.
- **The whole bleed schedule is referenced to `WA2`**, the compressor inlet flow, never to
  the local station flow.
- **`Q_C`, `Q_GT`, `Q_PT` all divide by the *instantaneous* rotor speed** (`1/NG`, `1/NP`),
  so all three are singular at zero speed. The report does not discuss starting or
  zero-speed behaviour anywhere in pp.15–25.
- **No overdot (ẋ) notation appears anywhere in pp.15–25.** Rate quantities are written
  either as explicit integrals (Eqs. 42–47) or as Leibniz derivatives (Eqs. 48–49). Overdot
  notation begins later, in the linear-model material.
- **The station 4.9 pressure function `f₁₀` and the station 4.1 heat-sink model are both
  Ames additions**, not part of the NASA Lewis hybrid model they started from — stated
  explicitly on pdf p.20, along with the power-turbine damping factor of Eq. 41.

---

## 9. Provenance

- Pages rendered with `pdftoppm -f 15 -l 25 -r 300 -png docs/ballin-tm100991.pdf
  docs/extracted/p`, then read as images.
- pp.21, 22, 23, 24, 25 re-rendered at **600 dpi** and cropped/zoomed for: the Figure 3
  leader lines (both halves), Eq. 23's `f` subscript (pixel-level), Eq. 39's brace
  contents and the `K_QC₁ H3 − H2` sign, Eqs. 42–44 in full, Eq. 43's two operators against
  the `+` glyphs of Eqs. 44 and 46, Eqs. 25–36 and 7–11 for the sub-subscripted constants.
- No value or symbol in this file was taken from `docs/extracted/p0NN.txt`.
- Constant *values* quoted in the symbol tables above are cross-references to Table A.1 as
  already transcribed in `inventory-appendix-a.md` [TM-100991 pdf p.55]; they are given for
  convenience and are not part of this page range.

**Ambiguities marked**: one — the broken `9` glyph in the left-hand side of Eq. 35 (`H49`),
resolved unambiguously by Eq. 36 and by station context. Nothing else in pp.15–25 was
ambiguous at 600 dpi.
