# Body — Fuel Control System Model, Validation, Results, Conclusions

Source: Ballin, *A High-Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
NASA TM-100991, 1988. **PDF pages 38–54** (printed pages 24–40; offset = 14).

Method: every page in this range was rendered at 300 dpi and **read from the image**.
Tables 2 and 3 and the accuracy sentences were re-rendered and re-read at 600 dpi
(and one glyph at 1200 dpi). Nothing here comes from the OCR text layer; the text layer
was used only to locate strings.

---

## 0. Headline findings

1. **There are no numbered equations in pp.38–54.** Not one. The range is prose, two
   tables (2, 3) and ten figures (6–15). The body therefore adds **no equations for the
   fuel control system** — Appendix C's block diagrams remain the *only* specification of
   the control system, exactly as `inventory-appendix-c.md` reports.
2. **The body's "FUEL CONTROL SYSTEM MODEL" section is one paragraph** [pdf p.38]. It is
   a scoping statement (what was simplified away), not a model description.
3. **Of the ten result figures, five are reproducible from this report alone**
   (Figs. 6, 7, 8, 9, 10) and **five are not** (Figs. 11–15 all require the external Gen
   Hel UH-60A blade-element helicopter simulation).
4. **Table 3 is a trap.** Its "Real-time model" column was produced with
   *compressor-mass-flow and turbine-energy functions derived from the NASA-Lewis test
   engine, supplied by NASA Lewis* — **not** the Appendix A functions [pdf p.39]. A model
   built from Appendix A cannot reproduce Table 3 and must not be tuned toward it.
5. **This range contains no linear-vs-nonlinear comparison.** The only such figure in the
   report is **Figure 5, pdf p.34** (linear model responses with and without the station
   4.1 heat sink), which is outside this range.

---

## 1. Page-by-page map, pp.38–54

| PDF p. | Printed p. | Content |
|---|---|---|
| 38 | 24 | End of "Pressure Function Table Solution" (time-step / 0.1 % error statement carries over from p.37); **§ Pressure Function Table Solution**; **§ FUEL CONTROL SYSTEM MODEL** (one paragraph); **§ VALIDATION** (one paragraph); **§ Static Trim** begins |
| 39 | 25 | Static Trim continued — discussion of Figs. 6–8 and of Tables 2–3; **§ Dynamic Response** begins; discussion of Figs. 9–10 begins |
| 40 | 26 | **Figure 6** — NG vs fuel flow, three simulation models |
| 41 | 27 | **Figure 7** — shaft horsepower vs fuel flow, three simulation models |
| 42 | 28 | **Figure 8** — station 3 static pressure vs NG, three simulation models |
| 43 | 29 | **Table 2** — test conditions for the NASA-Lewis experimental test engine (6 cases) |
| 44 | 30 | **Table 3** — real-time model steady state vs NASA-Lewis test article (6 cases × 2 rows) |
| 45 | 31 | **Figure 9** — response to a fuel-flow step increase 400 → 775 lbm/hr (6 stacked panels) |
| 46 | 32 | **Figure 10** — response to a fuel-flow step decrease 400 → 125 lbm/hr (6 stacked panels) |
| 47 | 33 | Text: end of Fig. 9 discussion; Fig. 10 discussion; closed-loop validation setup (Gen Hel UH-60A, time steps, flight-test data source); Fig. 11 discussion; Fig. 12 discussion begins |
| 48 | 34 | **Figure 11** — closed-loop response to a collective input at high forward speed (5 panels) |
| 49 | 35 | **Figure 12** — closed-loop response to a lateral cyclic input (5 panels) |
| 50 | 36 | Text: end of Fig. 12 discussion; Figs. 13, 14, 15 discussion; closing assessment of dynamic response and recommended refinements. **Page ends mid-page — no figure.** |
| 51 | 37 | **Figure 13** — closed-loop response to a collective input at hover (5 panels) |
| 52 | 38 | **Figure 14** — closed-loop response to a pedal input at hover (5 panels) |
| 53 | 39 | **Figure 15** — closed-loop response to a large collective input (4 panels) |
| 54 | 40 | **§ CONCLUSIONS** (two paragraphs). Page otherwise blank. |

---

## 2. THE RESULT-FIGURE CATALOGUE

Everything the report measures itself against. Read the "reproducible?" column first.

### 2.0 Summary table

| Fig. | pdf p. | What it compares | Reproducible by us? | Digitizing difficulty |
|---|---|---|---|---|
| 6 | 40 | steady-state NG vs Wf: real-time model vs 2 GE models | **yes** (needs trim NP — not printed) | easy–medium |
| 7 | 41 | steady-state SHP vs Wf: real-time model vs 2 GE models | **yes** (same caveat) | easy–medium |
| 8 | 42 | steady-state Ps3 vs NG: real-time model vs 2 GE models | **yes** (same caveat) | easy–medium |
| 9 | 45 | open-loop transient, fuel step **up**: real-time model vs GE status-81 | **yes** (NP held; value not printed) | medium |
| 10 | 46 | open-loop transient, fuel step **down**: real-time model vs GE status-81 | **yes** (same caveat) | medium |
| 11 | 48 | closed-loop, collective at 90 kt: Gen Hel + T700 vs UH-60A flight test | **no** — needs Gen Hel | medium–hard |
| 12 | 49 | closed-loop, lateral cyclic at 55 kt: same | **no** — needs Gen Hel | medium–hard |
| 13 | 51 | closed-loop, collective at hover: same | **no** — needs Gen Hel | medium |
| 14 | 52 | closed-loop, pedal at hover: same | **no** — needs Gen Hel | medium |
| 15 | 53 | closed-loop, large collective, 87 kt: Gen Hel + T700 vs Sikorsky test | **no** — needs Gen Hel | hard (4/rev oscillation) |

Common to **all ten**: linear–linear axes, **no gridlines**, boxed frame with minor tick
marks on all four sides, black-and-white line art. All are legible at 300 dpi; 600 dpi
resolves every overprinted marker cluster I checked.

---

### 2.1 Figure 6 — [pdf p.40, printed p.26]

**Caption, exactly as printed:**
> Figure 6: T700 gas generator speed vs. fuel flow for three simulation models at sea-level standard atmosphere conditions.

- **Comparison:** the report's real-time model against **two** GE analysis models. Not
  against hardware. Three-way model-to-model comparison.
- **x axis:** `FUEL FLOW, lb/hr` — **100 to 900**, labelled every 100.
- **y axis:** `GAS GENERATOR SPEED, percent` — **65.0 to 105.0**, labelled every 5.0.
- **Traces:** 3 discrete-marker series, keyed in an in-plot legend:
  - `●  REAL-TIME MODEL`
  - `△  G.E. UNBALANCED TORQUE MODEL`
  - `×  G.E. STATUS-81 MODEL`
- **Test condition:** steady-state trim sweep, sea-level standard atmosphere. Fuel flow is
  the independent variable. **The power turbine speed / load at which the sweep was run is
  not printed.**
- **Digitizing assessment: reproducible as a reference trace — easy to medium.**
  Discrete markers, no curve tracing, no gridlines. Counted from the image: **≈26–30 ●,
  ≈20–24 ×, ≈16–19 △** (counts carry ±2 because the three series overprint where the
  models agree — clusters at ≈81 %, 83 %, 85 %, 87 %, 96 %, 99 % NG need 600 dpi to
  separate, and one ● / △ pair at ≈79 % is a single blob even at 600 dpi).
  Marker diameter ≈0.28 % NG and ≈9 lb/hr, so centroid read error is about
  **±0.1 % NG and ±3 lb/hr**. Only the `●` series is our reference; the other two are
  context. **Budget: ~30 points.**
  *Artefact warning:* a small speck sits at roughly (315 lb/hr, 90.5 %) with no legend
  shape — it is not one of the three marker glyphs. Treat it as scan dirt, not data.

---

### 2.2 Figure 7 — [pdf p.41, printed p.27]

**Caption, exactly as printed:**
> Figure 7: T700 horsepower vs. fuel flow for three simulation models at sea-level standard atmosphere conditions.

- **Comparison:** same three models, same sweep, different output.
- **x axis:** `FUEL FLOW, lb/hr` — **100 to 900**, labelled every 100.
- **y axis:** `SHAFT HORSEPOWER` — **−200 to 2000**, labelled every 200. **No unit
  printed on the axis** beyond "HORSEPOWER"; the caption confirms horsepower.
- **Traces:** 3, same legend and same marker glyphs as Fig. 6.
- **Test condition:** as Fig. 6 — steady-state trim sweep at sea-level standard.
  Shaft horsepower depends on power turbine speed; **the NP at which it was evaluated is
  not printed.**
- **Digitizing assessment: easy to medium.** Discrete markers; the three series lie almost
  on top of one another below ~600 lb/hr, where the ● and × glyphs merge repeatedly, so
  the low-power end needs 600 dpi. **≈27 ●, ≈24 ×, ≈16 △**, each ±2. Centroid read error
  ≈**±10 hp, ±3 lb/hr**. **Budget: ~30 points** for the `●` series.

---

### 2.3 Figure 8 — [pdf p.42, printed p.28]

**Caption, exactly as printed:**
> Figure 8: T700 station 3 static-pressure vs. gas-generator speed for three simulation models at sea-level standard atmosphere conditions.

- **Comparison:** same three models. This is the plot of the *fuel control system's own
  input signal*, Ps3 — the report says so: "Static pressure at station 3, an input to the
  control system, is shown as a function of gas generator speed in figure 8" [p.39].
- **x axis:** `GAS GENERATOR SPEED, percent` — **65 to 105**, labelled every 5.
- **y axis:** `STATION 3 STATIC PRESSURE, psia` — **40 to 280**, labelled every 40.
- **Traces:** 3, same legend and glyphs.
- **Test condition:** as Figs. 6–7.
- **Digitizing assessment: easy.** The best-separated of the three scatter plots — the
  three series diverge visibly above 85 % NG and only overlap below ~80 %. **≈29 ●,
  ≈22 ×, ≈17 △**, each ±2. Centroid read error ≈**±1 psia, ±0.3 % NG**.
  **Budget: ~30 points.**

---

### 2.4 Figure 9 — [pdf p.45, printed p.31] — the primary open-loop transient reference

**Caption, exactly as printed:**
> Figure 9: Response to a step increase in fuel flow from 400 to 775 *lb<sub>m</sub>* per hour.

- **Comparison:** the real-time model (solid line) against the **GE performance-standard
  status-81 simulation**, plotted as discrete `+` markers under the legend name
  `REFERENCE STANDARD MODEL DATA`. (The identification of "reference standard" with the
  status-81 model comes from the body text on p.39 — "Open loop response was validated by
  comparison with the GE performance-standard status-81 simulation. Time history
  simulation data were provided by GE for large-step fuel flow inputs." The legend itself
  does not name status-81.)
- **Layout:** six stacked panels sharing one x axis.
- **x axis (shared):** `TIME, sec` — **0 to 5.0**, labelled every 1.0.
- **y axes, top to bottom** (none of the six carries a printed unit — see the unit note
  below):

  | Panel | Label | Range | Ticks labelled |
  |---|---|---|---|
  | 1 | `WFPH` | 250 to 1000 | 250, 500, 750, 1000 |
  | 2 | `PS3` | 100 to 300 | 100, 150, 200, 250, 300 |
  | 3 | `PCNG` | 80 to 100 | 80, 85, 90, 95, 100 |
  | 4 | `T41` | 2000 to 3000 | 2000, 2200, 2400, 2600, 2800, 3000 |
  | 5 | `T45` | 1250 to 2250 | 1250, 1450, 1650, 1850, 2050, 2250 |
  | 6 | `TORQ45` | 100 to 500 | 100, 200, 300, 400, 500 |

- **Traces:** 2 per panel — one solid line (real-time model) and one `+` marker series
  (reference standard) — **except panel 1 (`WFPH`, the input), which carries the solid
  line only.** So 11 traces on the page.
- **Test condition / input:** step increase in fuel flow **400 → 775 lbm/hr**, applied at
  **t ≈ 0.5 s** (read off the figure; the step time is not stated in text). Per the body
  [p.39]: the load is "a simple model representing the dynamometer used for testing of the
  NASA-Lewis experimental engine … variable, based on a simulated collective-pitch control
  input which was used to trim the power turbine at the design speed for a specified fuel
  flow", and **"the power turbine speeds were held constant by suppressing the NP
  integration. Output torque was therefore used as a measure of engine power."**
  So: NP integration off, NP held at (presumably) design speed — **the held value is not
  printed on the figure or in the text.**
- **Digitizing assessment: reproducible as a reference trace — medium.**
  The `+` markers are the reference to digitize; they are evenly spaced at
  **Δt ≈ 0.1 s** and run from t ≈ 0 to t ≈ 4.55 s → **≈46 markers per panel**, ×5 data
  panels = **≈230 points** for the reference series. Markers are individually resolved at
  600 dpi across every panel; the only crowding is in panels 2 and 6 late in the run,
  where consecutive markers touch but do not merge. The solid line is a continuous trace
  and would need conventional curve tracing (~60–100 samples/panel) if we want the
  report's own model output too — but for validation we only need the `+` reference and
  our own model's line. No gridlines; ticks are frequent enough on both axes to fix the
  affine transform precisely.
  **Unit note:** the panel labels are bare names. From the Appendix C nomenclature and the
  caption: `WFPH` = fuel flow, lbm/hr (**inferred** — WFPH is *not* in the Appendix C
  nomenclature; the caption's "400 to 775 lbm per hour" and the panel's 400→775 step fix
  it); `PS3` = station 3 static pressure, psia; `PCNG` = "rotational speed of compressor
  and gas generator in percent" [pdf p.78]; `T41`, `T45` = deg R; `TORQ45` = "power
  turbine torque (identical to QPT), ft·lbf" [pdf p.80].

---

### 2.5 Figure 10 — [pdf p.46, printed p.32]

**Caption, exactly as printed:**
> Figure 10: Response to a step decrease in fuel flow from 400 to 125 *lb<sub>m</sub>* per hour.

- **Comparison:** identical setup to Fig. 9 — real-time model (solid) vs reference
  standard model data (`+`).
- **x axis (shared):** `TIME, sec` — **0 to 5.0**, labelled every 1.0.
- **y axes, top to bottom:**

  | Panel | Label | Range | Ticks labelled |
  |---|---|---|---|
  | 1 | `WFPH` | 0 to 500 | 0, 250, 500 |
  | 2 | `PS3` | 0 to 200 | 0, 50, 100, 150, 200 |
  | 3 | `PCNG` | 60 to 100 | 60, 70, 80, 90, 100 |
  | 4 | `T41` | 1500 to 2500 | 1500, 1700, 1900, 2100, 2300, 2500 |
  | 5 | `T45` | 1000 to 2000 | 1000, 1200, 1400, 1600, 1800, 2000 |
  | 6 | `TORQ45` | 0 to 400 | 0, 100, 200, 300, 400 |

- **Traces:** 2 per panel, 1 in the input panel — 11 total, as Fig. 9.
- **Test condition / input:** step **decrease** 400 → 125 lbm/hr at **t ≈ 0.55 s** (read
  off the figure). Same suppressed-NP dynamometer setup as Fig. 9. The body calls this
  "a decrease in fuel flow to below-idle power" [p.47].
- **Digitizing assessment: medium**, same as Fig. 9. **≈46 `+` markers per panel at
  Δt ≈ 0.1 s**, ×5 data panels = **≈230 points**. Slightly harder than Fig. 9 in the
  `T45` panel, where the post-step markers sit on a nearly flat line and the `+` glyphs
  abut for the whole 1.0–4.5 s stretch; and in the `TORQ45` panel, where the model and
  reference both go to zero and the markers pile onto the axis. Both resolve at 600 dpi.

---

### 2.6 Figure 11 — [pdf p.48, printed p.34]

**Caption, exactly as printed:**
> Figure 11: Closed-loop response of engine, fuel control system, and blade-element helicopter simulation to a collective input at high forward speed.

- **Comparison:** **simulation vs flight test.** Legend, exactly as printed:
  - `———  SIKORSKY GEN HEL WITH T700 AND GEARBOX (AMES VERSION)`
  - `– – –  1982 USAAEFA UH-60A TEST DATA`
- **x axis (shared):** `TIME, sec` — **2.0 to 8.0**, labelled every 1.0.
- **y axes, top to bottom:**

  | Panel | Label | Range | Ticks labelled |
  |---|---|---|---|
  | 1 | `COLLECTIVE, in.` | 0 to 10.0 | 0, 2.0, 4.0, 6.0, 8.0, 10.0 |
  | 2 | `TOT. TORQUE, ft-lb × 10⁻³` | 18.0 to 28.0 | 18.0, 20.0, 22.0, 24.0, 26.0, 28.0 |
  | 3 | `GAS GEN. TURBINE, %` | 88.0 to 93.0 | 88.0, 89.0, 90.0, 91.0, 92.0, 93.0 |
  | 4 | `VEH. FUEL, lb/hr` | 500.0 to 900.0 | 500.0, 600.0, 700.0, 800.0, 900.0 |
  | 5 | `MAIN ROTOR SPEED, %` | 99.0 to 102.0 | 99.0, 99.5, 100.0, 100.5, 101.0, 101.5, 102.0 |

- **Traces:** 2 per panel (solid = simulation, dashed = test), 10 total.
- **Test condition / input:** single-axis **collective** input; aircraft at *light gross
  weight*, trimmed at **90 kt equivalent airspeed**; all stability augmentation disabled;
  input amplitude "typically not more than 15 percent of the control travel"; helicopter
  time step 14 ms with the T700 updated every 7 ms and the NP DOF every 14 ms [pdf p.47].
  Panel 2 is "the output torque of both engines" [p.47] — i.e. total, not per-engine.
- **Digitizing assessment: NOT usable as a validation reference for this project.**
  Both traces are products of a helicopter model we do not have. Mechanically it is a
  medium–hard continuous-curve trace: solid and dashed lines cross repeatedly in panels 2,
  3 and 5, and the `MAIN ROTOR SPEED` panel carries a visible high-frequency ripple on the
  solid trace that would alias unless sampled at ≲0.02 s. ~150–250 samples per panel.
  Digitize only if the Gen Hel boundary condition is ever supplied.

---

### 2.7 Figure 12 — [pdf p.49, printed p.35]

**Caption, exactly as printed:**
> Figure 12: Closed-loop response of engine, fuel control system, and blade-element helicopter simulation to a lateral cyclic input at high forward speed.

- **Comparison:** same two sources as Fig. 11 (identical legend text).
- **x axis (shared):** `TIME, sec` — **0.0 to 6.0**, labelled every 1.0.
- **y axes, top to bottom:**

  | Panel | Label | Range | Ticks labelled |
  |---|---|---|---|
  | 1 | `LATERAL CYCLIC, in.` | −5.0 to 5.0 | −5.0, −3.0, −1.0, 1.0, 3.0, 5.0 |
  | 2 | `TOT. TORQUE, ft-lb × 10⁻³` | 18.0 to 23.0 | 18.0, 19.0, 20.0, 21.0, 22.0, 23.0 |
  | 3 | `GAS GEN. TURBINE, %` | 88.0 to 91.0 | 88.0, 88.5, 89.0, 89.5, 90.0, 90.5, 91.0 |
  | 4 | `VEH. FUEL, lb/hr` | 500.0 to 700.0 | 500.0, 550.0, 600.0, 650.0, 700.0 |
  | 5 | `MAIN ROTOR SPEED, %` | 97.0 to 99.0 | 97.0, 97.5, 98.0, 98.5, 99.0 |

- **Traces:** 2 per panel, 10 total.
- **Test condition / input:** **lateral cyclic** input at **55 knots** [body, pdf p.47:
  "Figure 12 shows system response to a lateral cyclic control input at an airspeed of 55
  knots."]. **The caption says "high forward speed"** — see Ambiguity A-1 below.
- **Digitizing assessment: NOT usable as a validation reference.** Mechanically
  medium–hard, as Fig. 11; the `MAIN ROTOR SPEED` panel again carries a dense ripple on
  the solid trace.

---

### 2.8 Figure 13 — [pdf p.51, printed p.37]

**Caption, exactly as printed:**
> Figure 13: Closed-loop response of engine, fuel control system, and blade-element helicopter simulation to a collective input at hover.

- **Comparison:** same two sources (identical legend text as Fig. 11).
- **x axis (shared):** `TIME, sec` — **1.0 to 7.0**, labelled every 1.0.
- **y axes, top to bottom:**

  | Panel | Label | Range | Ticks labelled |
  |---|---|---|---|
  | 1 | `COLLECTIVE, in.` | 0 to 10.0 | 0, 2.0, 4.0, 6.0, 8.0, 10.0 |
  | 2 | `TOT. TORQUE, ft-lb × 10⁻³` | 20.0 to 50.0 | 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0 |
  | 3 | `GAS GEN. TURBINE, %` | 92.0 to 96.0 | 92.0, 93.0, 94.0, 95.0, 96.0 |
  | 4 | `VEH. FUEL, lb/hr` | 600.0 to 1100.0 | 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0 |
  | 5 | `MAIN ROTOR SPEED, %` | 98.0 to 101.0 | 98.0, 98.5, 99.0, 99.5, 100.0, 100.5, 101.0 |

- **Traces:** 2 per panel, 10 total.
- **Test condition / input:** **1-inch-down collective** input at **hover** [pdf p.50].
  This is the report's high-power case. The collective *trims* deliberately do not match
  — see the accuracy quotes.
- **Digitizing assessment: NOT usable as a validation reference.** Mechanically the
  easiest of Figs. 11–15: the two traces are well separated in every panel and never
  cross except in panels 3 and 5. ~120 samples per panel.

---

### 2.9 Figure 14 — [pdf p.52, printed p.38]

**Caption, exactly as printed:**
> Figure 14: Closed-loop response of engine, fuel control system, and blade-element helicopter simulation to a pedal input at hover.

- **Comparison:** same two sources (identical legend text as Fig. 11).
- **x axis (shared):** `TIME, sec` — **2.0 to 10.0**, labelled every 2.0.
- **y axes, top to bottom:**

  | Panel | Label | Range | Ticks labelled |
  |---|---|---|---|
  | 1 | `PEDALS, in.` | −3 to 3 | −3, −2, −1, 0, 1, 2, 3 |
  | 2 | `TOT. TORQUE, ft-lb × 10⁻³` | 20 to 50 | 20, 25, 30, 35, 40, 45, 50 |
  | 3 | `GAS GEN. TURBINE, %` | 92 to 97 | 92, 93, 94, 95, 96, 97 |
  | 4 | `VEH. FUEL, lb/hr` | 600 to 1200 | 600, 700, 800, 900, 1000, 1100, 1200 |
  | 5 | `MAIN ROTOR SPEED, %` | 98 to 102 | 98, 99, 100, 101, 102 |

- **Traces:** 2 per panel, 10 total.
- **Test condition / input:** **pedal** input at **hover**. The report's stated point:
  "Another high power case is shown **without the influences of the load demand
  compensation**" [pdf p.50] — a pedal input does not move the load-demand spindle, so
  this figure isolates the governor from the collective feed-forward.
- **Digitizing assessment: NOT usable as a validation reference.** Mechanically medium;
  traces are parallel and offset in panels 2–4, converging in panel 5.

---

### 2.10 Figure 15 — [pdf p.53, printed p.39]

**Caption, exactly as printed:**
> Figure 15: Closed-loop response of engine, fuel control system, and blade-element helicopter simulation to a large collective input.

- **Comparison:** simulation vs **a different test source**. Legend, exactly as printed:
  - `———  SIKORSKY GEN HEL WITH T700 AND GEARBOX (AMES VERSION)`
  - `– – –  SIKORSKY UH-60A TEST DATA`
  (note: **"SIKORSKY UH-60A TEST DATA"**, not the "1982 USAAEFA" data of Figs. 11–14;
  the body confirms — "The test data for this case were generated by Sikorsky" [p.50]).
- **x axis (shared):** `TIME, sec` — **0 to 6.0**, labelled every 1.0.
- **y axes, top to bottom — only 4 panels, not 5:**

  | Panel | Label | Range | Ticks labelled |
  |---|---|---|---|
  | 1 | `COLLECTIVE, in.` | 0 to 10.0 | 0, 2.0, 4.0, 6.0, 8.0, 10.0 |
  | 2 | `NORMAL ACCELERATION, g` | −1.0 to 2.0 | −1.0, −0.5, 0.0, 0.5, 1.0, 1.5, 2.0 |
  | 3 | `MAIN ROTOR TORQUE, ft-lb×10⁻³` | 0 to 30.0 | 0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0 |
  | 4 | `MAIN ROTOR SPEED, %` | 96.0 to 104.0 | 96.0, 98.0, 100.0, 102.0, 104.0 |

- **Traces:** 2 per panel, 8 total.
- **Test condition / input:** a **large transient** — "Control inputs were significant on
  all four axes, although the primary input was a lowering of collective from the trim
  position at **87 knots** to the full-down position in **one second**. This results in
  nearly zero-G flight immediately after the input" [pdf p.50]. Vehicle stabilator in an
  off-nominal position, producing large-amplitude first-harmonic rotor-hub moments.
- **Digitizing assessment: NOT usable as a validation reference, and the hardest of the
  ten anyway.** The `MAIN ROTOR TORQUE` panel carries a 4/rev oscillation on the solid
  trace at ≈0.06 s period, riding a large mean excursion; capturing it needs ≳500 samples
  and the two traces interpenetrate through the whole 2–5 s window. The report itself
  says the *test* trace's frequency is wrong ("an aliased frequency caused by a large
  test-data time step" [p.50]), so even a perfect digitization of the dashed trace is
  digitizing an artefact.

---

### 2.11 Tables 2 and 3 — the steady-state hardware comparison

Both transcribed **cell by cell from 600 dpi page images**. Every cell below was read
individually; none is OCR-derived.

**Table 2 caption, exactly as printed:**
> Table 2: Test conditions for NASA-Lewis experimental test engine.

| Test Case | W<sub>f</sub><br>lb<sub>m</sub>/hr | P2<br>PSIA | T2<br>°R | P49<br>PSIA | Load Torque<br>ft − lb<sub>f</sub> |
|---|---|---|---|---|---|
| 1 | 140.1 | 14.37 | 516.7 | 14.37 | 30.1 |
| 2 | 297.2 | 14.17 | 515.6 | 14.43 | 90.1 |
| 3 | 372.0 | 14.16 | 508.3 | 14.46 | 148.3 |
| 4 | 458.4 | 14.09 | 508.0 | 14.60 | 206.5 |
| 5 | 560.6 | 14.02 | 507.2 | 14.63 | 274.3 |
| 6 | 694.4 | 13.92 | 507.2 | 14.72 | 360.8 |

Column headings are set in italic maths type: `W_f` (true subscript f), `P2`, `T2`, `P49`
(full-size station indices, per the report's convention). Units are printed in a third
header row. Note **P49 is the *exhaust* boundary condition** — the report says "the inlet
and exhaust conditions, fuel flow, and load conditions were specified so that internal
states could be compared" [p.39].

**Table 3 caption, exactly as printed:**
> Table 3: Comparison of real-time model steady-state operation with NASA-Lewis test article. This does not represent specification T700 performance.

| Test Case | Data Type | NG<br>% | NP<br>% | WA2<br>lb<sub>m</sub>/sec | P3<br>PSIA | T3<br>°R | T45<br>°R |
|---|---|---|---|---|---|---|---|
| 1 | Experimental engine | 65.9 | 52.6 | 3.20 | 58.0 | 832.0 | 1413. |
| 1 | Real-time model | 70.0 | 83.8 | 3.59 | 63.4 | 855.9 | 1372. |
| 2 | Experimental engine | 84.7 | 95.7 | 5.16 | 113.1 | 1026. | 1577. |
| 2 | Real-time model | 86.4 | 100.6 | 5.52 | 119.0 | 1043. | 1555. |
| 3 | Experimental engine | 87.7 | 95.7 | 6.16 | 139.0 | 1081. | 1626. |
| 3 | Real-time model | 88.1 | 97.8 | 6.27 | 141.7 | 1084. | 1629. |
| 4 | Experimental engine | 90.4 | 95.7 | 6.92 | 161.1 | 1127. | 1731. |
| 4 | Real-time model | 90.7 | 97.3 | 7.00 | 163.3 | 1130. | 1734. |
| 5 | Experimental engine | 92.6 | 95.7 | 7.66 | 184.8 | 1173. | 1838. |
| 5 | Real-time model | 92.8 | 98.2 | 7.73 | 186.4 | 1178. | 1855. |
| 6 | Experimental engine | 95.9 | 95.7 | 8.50 | 211.9 | 1228. | 1974. |
| 6 | Real-time model | 96.1 | 99.4 | 8.53 | 213.8 | 1233. | 2009. |

Transcription notes, all confirmed at 600 dpi:
- Case 1 prints `832.0` / `855.9` for T3 (one decimal); cases 2–6 print T3 and T45 as
  four digits with a **trailing period and no decimal digit** (`1026.`, `1577.`). That is
  the report's own formatting, not a lost digit.
- NP is 95.7 % for the experimental engine in **all** of cases 2–6 — the dynamometer load
  torque was adjusted to hold a specified power turbine speed, while "the model power
  turbine speed was allowed to vary in order to achieve trim with the test load torque"
  [p.39]. So the NP column is *not* a like-for-like comparison by construction.
- Columns are `NG`, `NP`, `WA2`, `P3`, `T3`, `T45` — full-size station indices in italic
  maths type; note `WA2` is printed with a space, `W A2`, but is the single symbol WA2.

**Reproducibility of Table 3 — read this before using it:** the "Real-time model" rows
were **not** produced with the Appendix A functions. [pdf p.39]:

> "The experimental engine is a prototype model which does not reproduce specification engine performance. Because of this, compressor-mass-flow and turbine-energy functions derived from this test engine were used in place of the standard functions which represent specification performance. These functions were supplied by NASA Lewis."

Those substitute functions are **nowhere in the report**. Table 3 is therefore evidence
about the *modelling method*, not a reproducible benchmark for a model built from
Appendix A. It should not appear in `validation/` as a pass/fail target, and no constant
may ever be tuned toward it. (The caption's own second sentence — "This does not
represent specification T700 performance" — is saying the same thing.)

---

## 3. Every accuracy and fidelity statement, quoted exactly

Ordered by page. These replace the provisional tolerances in `SCOPE.md`.

### Integration / numerics

> "Because of the opened compressor mass-flow iteration, however, response is dependent on time step. For the T700 implementation, a restriction on time step to a maximum value of 10 msec was established based on a maximum allowable error of 0.1 percent between time steps." — [pdf p.38]

> "A time step of 14 msec was used; this is a typical value for real-time execution of the blade-element rotor. In order to meet the cycling requirements of the iteration model, the T700 program was updated twice for each rotor routine cycle, or once every 7 msec. The power-turbine-speed degree of freedom corresponds to that of the drive train and rotor hub, and was therefore updated every 14 msec." — [pdf p.47]

### Validation basis

> "The model was validated by comparison of static trim performance and dynamic response with available information, which was provided by GF and NASA Lewis." — [pdf p.38]
> *(The page really does print "GF"; verified at 1200 dpi — the glyph has no lower bar and no lower serif. It is a typo for GE. See Ambiguity A-2.)*

> "Steady-state simulation performance was verified to be within normal limits of operation by comparison with data supplied from two computer simulations developed by the engine manufacturer. The GE status-81 model is a comprehensive analysis-oriented model which is used for detailed representation of the engine thermodynamic cycle. The GE unbalanced torque model is a simplified model which is optimized to reproduce engine dynamics for control system design purposes." — [pdf pp.38–39]

### Steady state — Figures 6–8

> "In all cases, the real-time model displays an acceptable steady-state performance." — [pdf p.39]

> "The steady-state operation of the real-time model is bounded by the two analysis models over the range of operation except between 81 and 86 percent of gas generator speed. In this area, the real-time model tends to overestimate fuel consumption by as much as five percent." — [pdf p.39]

> "At higher gas generator speeds, the real-time model displays the characteristic of the GE unbalanced torque model, requiring less fuel for a given fuel flow than the GE status-81 model." — [pdf p.39]
> *(printed as written — "requiring less fuel for a given fuel flow" is the report's own wording and is self-contradictory; see Ambiguity A-3.)*

> "Because of the limits of the data supplied in the real-time model functional relations, the maximum gas generator speed for which the model is valid is 100 percent. This is adequate for the intended use of the model because the fuel control system prevents steady-state operation outside of this range." — [pdf p.39]

> "The real-time model displays agreement with the status-81 model over the entire operating range of the engine. The unbalanced torque model shows less available horsepower for a given fuel flow at high power levels." — [pdf p.39, on Fig. 7]

> "The same trends shown in figure 6 are seen, with the real-time model showing closer agreement with the unbalanced torque model over most of the range of operation." — [pdf p.39, on Fig. 8]

### Steady state — against hardware, Tables 2–3

> "Fair agreement is seen in the medium- and high-power test cases (cases 3 through 6). The real-time model tends to overestimate power output by a small percentage in all cases. Internal temperatures and pressures agree very well. At lower power settings, agreement is poorer. The model overestimates gas generator speed by 4 percent in case 1. Although the model is not valid for such low power turbine speeds, the major difference seen in these speeds is caused by the difference in gas generator speeds and not in the power turbine model." — [pdf p.39]

### Open-loop dynamics — Figures 9–10

> "As shown in figure 9, the two simulations are in close agreement for a step increase from midpower to high power. Gas generator speed is overestimated by 1 to 2 percent; this is reflected in the trim differences between the real-time model and the status-81 model in figure 6." — [pdf p.39]

> "Trim values of station 4.1 and 4.5 temperature are slightly underestimated by the real-time program; this is a characteristic of the real-time model which was found for all validation cases. Output torque, station 3 static pressure, and temperature responses are in good agreement." — [pdf pp.39–47]

> "Low power engine performance is shown in figure 10, which represents a decrease in fuel flow to below-idle power. Gas generator dynamics are accurately represented. Other real-time model outputs display a slightly different dynamic characteristic, although zero torque is reached for both simulations at approximately the same time. The real-time model is initially less responsive. Under a simplifying assumption, the real-time heat-sink model constants were made independent of the direction of power change. Additional sophistication of the heat-sink representation may be warranted." — [pdf p.47]

### Closed-loop dynamics — Figures 11–15

> "The correlation of secondary effects such as propulsion system response depends on adequate correlation of the vehicle and rotor system responses with the test data. Although the blade-element helicopter simulation is considered to be a high-fidelity model, differences do exist which are reflected in the engine and fuel-control-system responses." — [pdf p.47]

> "The amplitudes of the inputs are small, typically not more than 15 percent of the control travel." — [pdf p.47]

> "Collective trim is in good agreement with that of the flight test. Power required by the aircraft is therefore correctly represented by the vehicle simulation. The output torque of both engines is shown by the second plot in the figure; this is in good agreement with test data for the initial trim. The trends displayed by the rotor speed time histories are correct. The data suggest that the test aircraft contains greater load demand compensation than the simulation, however, resulting in less rotor speed droop after the initial input. Gas generator speed and output torque time histories are also more responsive to the initial input." — [pdf p.47, Fig. 11]

> "After the input, the test vehicle and simulation time histories diverged in pitch, with the vehicle reaching a 6-degree nose-up attitude at 6 seconds, while the simulation achieved a 2-degree nose-up attitude." — [pdf p.47, Fig. 11]

> "The lateral cyclic trim position is in good agreement with the test, while the collective trim position is slightly higher in the simulation. The engine torque is therefore higher." — [pdf p.47, Fig. 12]

> **"Rotor speed is shown to agree to within 0.2 percent over the duration of the run."** — [pdf p.50, Fig. 12] *(the only numerical closed-loop tolerance the report states)*

> "In this case, the collective trims do not match well, a result of the simplified vehicle simulation of rotor downwash impingement on the fuselage and stabilator when in hover (ref. 13). Because the simulation requires less power, its torque, fuel flow, and gas generator trim values are lower than their test data counterparts. Trends in the data are reproduced well, however. As in all cases, the simulation fuel flow appears to be more oscillatory than the test data. This is attributed to the location of the sensor used in the test vehicle; it was mounted upstream of both engines and therefore did not correctly represent the fuel flow transients." — [pdf p.50, Fig. 13]

> "As in all hover cases, the engine operates at a slightly lower power output because of the incorrect collective-trim position. As shown by the figure, the pedal trim is correct. Rotor speed changes in the test data appear to be be greater, despite greater gas-generator speed changes." — [pdf p.50, Fig. 14] *("to be be" is the report's own duplication.)*

> "After 7.5 seconds, the test vehicle and the simulation had diverged in pitch, resulting in a poor rotor speed match after this point." — [pdf p.50, Fig. 14]

> "The frequency shown by the simulation time history is correct; it is equal to the number of blades times the rotor rotational frequency. The frequency shown by the test data is not correct. It is an aliased frequency caused by a large test-data time step. Rotor speed trends are shown to be in good agreement." — [pdf p.50, Fig. 15]

### The report's own closing assessment of dynamic fidelity

> "Dynamic response of the propulsion system is at a level of fidelity comparable to that of the blade-element helicopter simulation. Propulsion system damping is slightly greater than that indicated by the test data, and the load demand compensation is greater in the test data. For cases which are not influenced by the load demand compensation, rotor speed variations appear to be slightly larger in the test data, although there is evidence in some instances that changes in torque required by the vehicle are greater than those required by the vehicle simulation. All mechanical actuator and sensor nonlinearities were modeled with lags, transport delays, and hysteresis loops as provided by GE. Better correlation may possibly be attained by modification of these simple models. Greater model sophistication may also be necessary. Additions may include an explicit variable-geometry guide vane model with dynamics and heat-sink model constants which are a function of increasing and decreasing power. A small effective lag may be added with the inclusion of a T45 heat-sink model." — [pdf p.50]

### What this gives us as tolerances

Everything the report quantifies about itself, gathered:

| Claim | Value | Scope | Page |
|---|---|---|---|
| max inter-step integration error | **0.1 %** | sets the 10 ms step cap | 38 |
| fuel consumption overestimate vs GE models | **≤ 5 %** | only in 81–86 % NG | 39 |
| NG upper validity limit | **100 % NG** | steady state | 39 |
| NG overestimate, lowest-power hardware case | **4 %** | Table 3 case 1 | 39 |
| NG overestimate, open-loop step | **1 to 2 %** | Fig. 9 | 39 |
| rotor speed agreement, closed loop | **within 0.2 %** | Fig. 12 only | 50 |

Everything else is qualitative: "acceptable", "fair agreement", "agree very well",
"good agreement", "accurately represented", "close agreement", "slightly underestimated",
"trends … reproduced well". **The report states no percentage tolerance for pressures,
temperatures, torque, or fuel flow anywhere.** So the provisional table in `SCOPE.md`
cannot be wholly replaced — only these six rows are the report's own numbers, and they
are mostly *deviations the report accepted*, not tolerances it met.

---

## 4. The fuel control system section in the body [pdf p.38]

The whole of the section, transcribed exactly:

> **FUEL CONTROL SYSTEM MODEL**
>
> The real-time digital implementation of the fuel control system consists of simplified versions of the ECU and the HMU and models of pressure, temperature, torque, and speed sensors. Collective pitch rigging to the load demand spindle for the UH-60A implementation is also provided. Each of these is modeled as an explicit entity such that all interfaces between components represent interfaces in the actual control system. Simplifications were made by eliminating models and control logic which are beyond the scope of real-time simulation requirements or are not needed because of simplifications to the engine model. Eliminated model features include automatic engine start-up capability, fuel control below flight-idle power, position control of the variable-geometry inlet-guide-vanes, redundancy models, and redundancy model logic. The complete fuel control system model is presented in Appendix C.

That is the entire section. **One paragraph, no equations, no constants, no figure of its
own.**

### How it relates to Appendix C

| Body statement (p.38) | Where it lands in Appendix C |
|---|---|
| "simplified versions of the ECU and the HMU" | Fig. C1 (ECU top level, p.85) and Fig. C9 (HMU top level, p.89) |
| "models of pressure, temperature, torque, and speed sensors" | Fig. C11 Ps3 sensor (p.90); Fig. C5 thermocouple harness (p.87); Fig. C2 load-share torque lags (p.85); Fig. C15 NG spool sensor (p.91) |
| "Collective pitch rigging to the load demand spindle for the UH-60A implementation" | **Fig. C13** — collective pitch to load demand spindle rigging (UH-60A), p.90 |
| "all interfaces between components represent interfaces in the actual control system" | why Appendix C is 22 separate block diagrams rather than one lumped transfer function |
| **eliminated:** automatic engine start-up | nothing in Appendix C models start — consistent |
| **eliminated:** fuel control below flight-idle power | explains why the idle schedule (Figs. C20, C28, C29) is a *floor*, not a start-up path |
| **eliminated:** position control of the variable-geometry inlet-guide-vanes | no IGV block anywhere in Appendix C. The physical description [pdf p.19] says the real engine controls "compressor variable geometry … as a function of inlet temperature and gas generator speed"; the model does not. The Conclusions [p.54] list "an explicit variable-geometry guide vane model with dynamics" as future work |
| **eliminated:** redundancy models and redundancy model logic | explains the **"ONE ENGINE IMPLEMENTATION"** switch drawn *open* in Fig. C2 (p.85) and the two-engine ZLOLIM logic of Fig. C8 (p.88) being effectively single-valued |

The **physical** (not model) description of the fuel control system — HMU vane pump,
mechanical cams, ECU torque motor, fail-to-max-power behaviour, load sharing — is on
**pdf p.19** (printed p.5), *not* in this range. Worth reading alongside Appendix C:

> "It consists of a hydromechanical control unit (HMU) for fuel metering as a function of schedules of gas generator speed and power demand, and an electrical control unit (ECU) which performs isochronous power-turbine speed governing and overtemperature protection (ref. 6). … The torque motor adjusts the HMU fuel demand downward, so that an electrical system failure results in maximum power." — [pdf p.19]

That last sentence is a useful sign check on Fig. C9: TMRU enters the WFPDM summing
junction with a **negative** sign, and `inventory-appendix-c.md` reads exactly that.

---

## 5. Every numbered equation in pp.38–54

**None.** Zero numbered equations appear on any page from 38 to 54 inclusive. Confirmed by
reading all seventeen page images; pp.40–46 and 48–49 and 51–53 are figures and tables,
pp.38, 39, 47, 50 and 54 are running text with no display maths, and no equation number
appears anywhere in the range. The report's numbered equations stop before p.38 (the last
one seen in the preceding range is Eq. 65 area, pdf p.32) and do not resume.

---

## 6. Conclusions [pdf p.54, printed p.40] — transcribed in full

> **CONCLUSIONS**
>
> As the maneuvering envelope of helicopters is widened for increasingly demanding mission tasks, the associated large-amplitude transients in aircraft power require high-fidelity modeling of the propulsion system. The real-time digital simulation of a small turboshaft engine fills this need in pilot-in-the-loop handling qualities investigations involving such power transients. Applications include real-time studies of the effect of rotor speed variation on handling qualities, investigations of new fuel control and flight control methodologies, and simulations of rotorcraft engine degradation and failure. The model adequately reproduces trim performance over the complete flight-power operating range as well as dynamics associated with changing load conditions. Engine degradation is easily modeled by modifying compressor or turbine flow and energy functions. The digital fuel control system model is separate and may be modified or replaced depending on user requirements.
>
> Validation results suggest that the static and dynamic fidelity of the model is within the limits of the fidelity of current rotorcraft simulations. The modeling of high-speed dynamics which represent changes of mass flow between the internal control volumes was found to be unnecessary. Several refinements were found to be necessary to obtain correct propulsion system response, however. These include estimates of heat transfer to the engine components downstream of the combustor, estimates of losses between the power turbine outlet and the engine exhaust, power-turbine-speed damping, and sensor and actuator dynamics and nonlinearities.

### What Ballin claims the model does well

- Reproduces **trim performance over the complete flight-power operating range**.
- Reproduces **dynamics associated with changing load conditions**.
- Static and dynamic fidelity **"within the limits of the fidelity of current rotorcraft
  simulations"** — i.e. he claims the engine model is not the weak link; the helicopter
  model is.
- **Engine degradation** is modelled simply by modifying the compressor/turbine flow and
  energy functions (this is a claim about the *architecture*, and it is why Appendix A's
  functions must stay data, not be baked into code).
- The fuel control model is **separable and replaceable**.

### The limitations he admits — these are ours

1. **NG validity ceiling of 100 %** [p.39] — set by the extent of the supplied functional
   relations, not by physics. Our digitized Appendix A curves inherit exactly this range.
2. **Not valid at low power turbine speeds** [p.39, discussing Table 3 case 1].
3. **Systematic bias:** T41 and T45 trim values are **slightly underestimated in every
   validation case** [p.47]. Expect our model to show the same signed bias and do not
   "fix" it.
4. **Heat-sink constants are direction-independent** — the same constants are used for
   increasing and decreasing power, which is a *simplifying assumption* he flags twice
   (p.47 and p.50) and which he identifies as the cause of the sluggish initial response
   in Fig. 10. Replicate as written; the Conclusions' "heat-sink model constants which are
   a function of increasing and decreasing power" is future work, not this model.
5. **No T45 heat-sink model** — "A small effective lag may be added with the inclusion of
   a T45 heat-sink model" [p.50]. Only station 4.1 has one.
6. **No variable-geometry inlet-guide-vane dynamics** [p.38, p.50].
7. **Propulsion system damping is slightly greater than the test data indicates** [p.50].
8. **Load demand compensation is weaker than the real aircraft's** [p.47, p.50] — the
   report sees this in three separate figures (11, 13, 15) as too much rotor droop / a
   missing rotor overspeed.
9. **Response is time-step dependent** because of the opened compressor mass-flow
   iteration; 10 ms is the cap for 0.1 % inter-step error [p.38].
10. **High-speed inter-volume mass-flow dynamics were deliberately omitted** and found
    unnecessary [p.54]. This is a *positive* finding about the 5-state formulation and
    should be recorded as such — it justifies the state count.

---

## 7. What this means for Phase 6

- **Digitize five figures** into `data/reference/`: Figs. 6, 7, 8 (steady state; ~30
  points each, `●` series only) and Figs. 9, 10 (transient; ~230 `+` points each across
  five panels). Total ≈ 550 reference points. All are marker-based, so digitizing is
  **point extraction, not curve tracing** — the same happy situation as Appendix C's
  schedule plots.
- **Transcribe, do not digitize,** Tables 2 and 3 — already done above. But Table 3 is
  **not** a pass/fail target (§2.11).
- **Figures 11–15 are out of reach** without the Gen Hel UH-60A simulation, which is out
  of scope by `SCOPE.md`'s own boundary. Five of the report's ten result figures are
  therefore unreproducible *by construction*, and that should be stated plainly in
  `SCOPE.md` rather than discovered in Phase 6.
- **The strongest validation line in the report remains Appendix B's matrices**, exactly
  as `SCOPE.md` already argues — they are printed numbers with no read error. Figs. 6–10
  are the second line.

---

## 8. Ambiguities and marked items

**A-1. Figure 12's caption contradicts the body.** The caption reads "to a lateral cyclic
input at **high forward speed**" [p.49]; the body says "figure 12 shows system response to
a lateral cyclic control input at an **airspeed of 55 knots**" [p.47]. 55 kt is not high
forward speed for a UH-60A, and Fig. 11 (the genuinely high-speed case) is 90 kt. Read at
600 dpi; both are as printed. Most likely the caption was copied from Fig. 11. Use 55 kt.

**A-2. p.38 prints "GF" for "GE".** "…which was provided by **GF** and NASA Lewis."
Verified at 1200 dpi: the second glyph has neither a lower horizontal bar nor a lower
serif — it is an F, not a damaged E. A report typo, recorded so nobody re-reads the page.

**A-3. p.39, "requiring less fuel for a given fuel flow".** Printed exactly so. It is
internally contradictory; from Fig. 6 the intended sense is *"requiring less fuel for a
given gas generator speed"* — but that is inference, and the sentence is transcribed above
as printed. Do not silently correct it in a quotation.

**A-4. p.50, "appear to be be greater".** The report duplicates "be". Transcribed as
printed.

**A-5. `WFPH` is undefined.** It labels panel 1 of Figs. 9 and 10 and appears nowhere in
the main nomenclature (pdf pp.9–13) or the Appendix C nomenclature (pdf pp.77–82). From
the captions ("from 400 to 775 lbm per hour") and the panel values it is fuel flow in
lbm/hr — i.e. `W_f` expressed per hour. **Inferred, not printed.**

**A-6. "REFERENCE STANDARD MODEL DATA" is not named in the legend.** The legend of
Figs. 9 and 10 says only that. The identification with the GE status-81 model rests on the
body sentence "Open loop response was validated by comparison with the GE
performance-standard status-81 simulation" [p.39]. High confidence, but it is a
cross-reference, not a legend entry.

**A-7. A stray speck in Figure 6** at approximately (315 lb/hr, 90.5 %), matching none of
the three legend glyphs. Treated as scan dirt. Whoever digitizes Fig. 6 must not pick it
up as a fourth series.

**A-8. Figure 15 uses a different test data source** from Figs. 11–14 — "SIKORSKY UH-60A
TEST DATA" rather than "1982 USAAEFA UH-60A TEST DATA". Confirmed from both the legend
and the body ("The test data for this case were generated by Sikorsky", p.50).
