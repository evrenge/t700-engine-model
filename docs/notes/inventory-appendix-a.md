# Inventory — Appendix A: T700 Engine Model Constants and Function Tables

Source: `docs/ballin-tm100991.pdf`, **PDF pages 55–66** (printed pages 41–52).
Every number below was read from the rendered page image, never from the text layer.
Page images: `docs/extracted/p-055.png` … `p-066.png` (200 dpi); ambiguous glyphs re-rendered
at 300 dpi and 600 dpi with `pdftoppm -r 300|600`.

---

## THE HEADLINE ANSWER

> **Appendix A contains exactly one numeric table and eleven plots.**
>
> The only numbers printed anywhere in Appendix A are the **33 scalar constants of
> Table A.1** (pdf p.55). **Every one of the model's eleven functional relationships —
> `f₁` through `f₁₀` and `f_hs` — exists only as a printed x–y plot.** There is no
> numeric function table, no breakpoint list, no tabulated map, and no digital appendix,
> anywhere in pp.55–66.
>
> The appendix title says "…AND FUNCTION TABLES". **That title is misleading.** What is
> printed are function *plots*. The report body is the accurate description: p.22 calls
> them "plots of mass-flow and energy functions".

**Consequence for the next phase:** all eleven engine map/schedule relationships must be
**digitized from the figures**. Nothing about them can be transcribed. This resolves
open question #1.

The digitizing job is smaller than "eleven maps" sounds, because the figures are plotted
as **discrete marked data points joined by straight or lightly-smoothed segments**, not as
freehand curves. Counting markers: five figures (A4, A5, A6, A7, A11) carry **six or fewer
points each**; only Figure A1 is a genuine two-dimensional map. Estimated total point
count across all eleven figures: **≈ 190 points**, of which ≈ 77 belong to Figure A1.

---

## 1. Summary table — every functional relationship Appendix A supplies

**Status column updated 2026-09-12.** It read **must-digitize** against all eleven
relationships long after all eleven were digitized; `data/maps/` holds 183 points across 11
files, and `tools/reproduce_all.sh` reproduces every one byte for byte.

| Relationship | Model use (report equation) | Where in App. A | Status | Points to recover |
|---|---|---|---|---|
| 33 scalar engine constants | throughout | pdf p.55, Table A.1 | **transcribed** (below) | 33 values, done |
| `f₁`: WA2c = f₁(PS3/P2, NGc) | compressor mass flow, Eq. 7 (p.22) | pdf p.56, Fig. A1 | **digitized** | 11 curves × ~6 pts + 11 flat extensions ≈ 77 |
| `f₂`: T3 = T2·f₂(PS3/P2) | compressor temperature, Eq. 10 (p.22) | pdf p.57, Fig. A2 | **digitized** | ~20 |
| `f₃`: B1 = f₃(NGc) | seal-pressurization bleed, Eq. 12 (p.23) | pdf p.58, Fig. A3 | **digitized** | ~17 |
| `f₄`: B2 = f₄(WA2c) | power-turbine-balance bleed, Eq. 13 (p.23) | pdf p.59, Fig. A4 | **digitized** | 3 |
| `f₅`: B3 = f₅(WA2c) | impeller-tip-leakage + turbine cooling bleed, Eq. 14 (p.23) | pdf p.60, Fig. A5 | **digitized** | 3 |
| `f₆`: η = f₆(FAR) | combustor efficiency, Eq. 20 (p.23) | pdf p.61, Fig. A6 | **digitized** | 2 (a constant) |
| `f₇`: ΔH_GT = θ₄₁·f₇(P45/P41) | gas-generator turbine energy, Eq. 26 (p.24) | pdf p.62, Fig. A7 | **digitized** | 6 |
| `f₈`: ΔH_PT = θ₄₅·f₈(P49/P45) | power turbine energy, Eq. 32 (p.24) | pdf p.63, Fig. A8 | **digitized** | 12 |
| `f₉`: W45c = f₉(Ps9/P45) | power turbine mass flow, Eq. 33 (p.24) | pdf p.64, Fig. A9 | **digitized** | ~23 |
| `f₁₀`: P49 = Ps9·f₁₀(NGc) | exhaust pressure loss, Eq. 38 (p.24) | pdf p.65, Fig. A10 | **digitized** | ~20 |
| `f_hs`: T41sgn = f_hs(NGc) | station 4.1 heat-sink constant, Eq. 52 (p.26) | pdf p.66, Fig. A11 | **digitized** | 6 |

`f_s` (Eq. 23, p.23 — `T41 = T41_ns · f_s(T41, T41_ns, W41, NGc)`) is **not** an Appendix A
figure. It is the heat-sink transfer function, defined analytically in the body by
Eqs. 48–53 (pp.25–26) from `TC_T41` (Table A.1) and `f_hs` (Fig. A11). No gap.

---

## 2. Page-by-page map, pdf pp.55–66

| PDF p. | Printed p. | Content |
|---|---|---|
| 55 | 41 | Appendix A title page. One sentence of body text: *"Model constants and functional relationships specific to the T700-GE-700 engine are given below."* Then **Table A.1, "T700 Engine Constants"** — 33 rows, three columns (Constant / Value / Units). Nothing else on the page. |
| 56 | 42 | **Figure A1** — `f₁`, compressor mass flow. The only 2-D map in the appendix (11 speed lines). |
| 57 | 43 | **Figure A2** — `f₂`, compressor temperature. Single curve. |
| 58 | 44 | **Figure A3** — `f₃`, seal-pressurization bleed fraction. Single curve. |
| 59 | 45 | **Figure A4** — `f₄`, power-turbine-balance bleed fraction. Single 3-point piecewise line. |
| 60 | 46 | **Figure A5** — `f₅`, impeller tip leakage and turbine cooling-bleed fraction. Single 3-point piecewise line. |
| 61 | 47 | **Figure A6** — `f₆`, combustor efficiency. Single horizontal line, 2 endpoints. |
| 62 | 48 | **Figure A7** — `f₇`, gas-generator turbine energy. Single 6-point curve. |
| 63 | 49 | **Figure A8** — `f₈`, power turbine energy. Single 12-point curve. |
| 64 | 50 | **Figure A9** — `f₉`, power turbine mass flow. Single ~23-point curve. |
| 65 | 51 | **Figure A10** — `f₁₀`, exhaust pressure loss. Single ~20-point curve. |
| 66 | 52 | **Figure A11** — `f_hs`, station 4.1 heat-sink constant. Single 6-point piecewise line. |

Each figure occupies its whole page: one plot, one caption, printed page number at the
foot. No text, no tabulated values, no secondary panels on any of pp.56–66.

---

## 3. Table A.1 — full transcription

**Caption exactly as printed** (two centred lines above the rule):

```
Table A.1
T700 Engine Constants
```

**Column headings exactly as printed:** `Constant` | `Value` | `Units`

Transcribed cell by cell from `docs/extracted/p-055.png` at 200 dpi, then re-read at
300 dpi in two halves and at 600 dpi for the fractional-exponent cell. Values are given
with the **printed precision preserved exactly** — no rounding, no normalising.
Citation for every row: `[TM-100991 pdf p.55, Table A.1]`.

| # | Constant | Value | Units |
|---|---|---|---|
| 1 | HVF | 18300.0 | Btu/lb_m |
| 2 | J_GT | 0.0445 | ft · lb_f · sec² |
| 3 | J_PT | 0.062 | ft · lb_f · sec² |
| 4 | K_bl | 0.7826 | nondimensional |
| 5 | K_b3 | 0.0025 | nondimensional |
| 6 | K_damp | 0.06854 | ft · lb_f · sec/rad |
| 7 | K_dpb | 0.03045 | lb_f² · sec² / lb_m² · in⁴ · deg R |
| 8 | K_H2 | 0.239 | Btu/lb_m · deg R |
| 9 | K_H3₁ | 0.2496 | Btu/lb_m · deg R |
| 10 | K_H3₂ | -8.4 | Btu/lb_m |
| 11 | K_H41₁ | 0.3010 | Btu/lb_m · deg R |
| 12 | K_H41₂ | -86.905 | Btu/lb_m |
| 13 | K_H45 | 0.9623 | nondimensional |
| 14 | K_Ps3 | 0.956 | nondimensional |
| 15 | K_QC₁ | 0.71 | nondimensional |
| 16 | K_QC₂ | 0.29 | nondimensional |
| 17 | K_T41₁ | 3.322 | lb_m · deg R/Btu |
| 18 | K_T41₂ | 288.7 | deg R |
| 19 | K_T45₁ | 3.519 | lb_m · deg R/Btu |
| 20 | K_T45₂ | 179.1 | deg R |
| 21 | K_T49₁ | 3.516 | lb_m · deg R/Btu |
| 22 | K_T49₂ | 172.3 | deg R |
| 23 | K_TH41₁ | 0.0018326 | 1/deg R |
| 24 | K_TH41₂ | 0.0856 | nondimensional |
| 25 | K_TH45₁ | 0.0018326 | 1/deg R |
| 26 | K_TH45₂ | 0.0856 | nondimensional |
| 27 | K_V3 | 0.97 | lb_f/in² · lb_m · deg R |
| 28 | K_V41 | 6.17 | lb_f/in² · lb_m · deg R |
| 29 | K_V45 | 13.63 | lb_f/in² · lb_m · deg R |
| 30 | K_WGT | 0.0876 | lb_m · in²/lb_f · sec |
| 31 | NG_des | 44700.0 | rpm |
| 32 | NP_des | 20900.0 | rpm |
| 33 | TC_T41 | 0.29 | lb_m^(4/5) · sec^(9/5) / deg R^(1/2)  ⚠ see A-1 below |

**Row count: 33.** Nothing in the table is elided or continued on another page.

### Notes on specific cells

- **Subscript typography.** The report sets station indices full size (`K_T41₁` is
  K-sub-(T-4-1)-sub-1) and true subscripts small and lowered. Sixteen constants carry a
  subscript *on* a subscript — the `₁`/`₂` pairs are the slope and intercept of one linear
  fit and **carry different units** (e.g. rows 11/12, 17/18). Do not merge them.
- **`K_bl` is b-ell, not b-one.** Read at 600 dpi the subscript is an italic lowercase `l`,
  visibly different from the digit `3` in `K_b3` directly below it. Confirmed against the
  nomenclature, pdf p.10: *"K_bl — fraction of diffuser bleed gas which is used to cool the
  gas generator turbine blades"*. (`bl` = bleed.) The text layer renders it `gb1`, which
  is where the earlier `Kb1` reading in `index.md` came from — that reading is wrong.
- **Negative values are real.** `K_H3₂ = -8.4` and `K_H41₂ = -86.905` both carry a clearly
  printed minus sign at 300 dpi. The text layer drops minus signs, so these two were
  specifically re-checked.
- **`K_TH41₁` and `K_TH45₁` are identical** (0.0018326), as are `K_TH41₂` and `K_TH45₂`
  (0.0856). This is what the page prints; it is not a transcription slip.
- **`K_dpb` units**: the superscripts are `lb_f` squared, `sec` squared, `lb_m` squared,
  `in` to the fourth. Read at 300 dpi.

### A-1. The one ambiguous / problematic cell

**`TC_T41` units — the `sec` exponent numerator.**

Printed as stacked fractional exponents at a size the 1988 phototypesetter and the 2009
scan barely resolve. Read at 300 dpi and again at 600 dpi with pixel-level magnification:

- `lb_m` exponent numerator **4**, denominator **5** — unambiguous.
- `deg R` exponent **1/2** — unambiguous.
- `sec` exponent denominator **5** — same glyph as the `lb_m` denominator, unambiguous.
- `sec` exponent numerator: **`9`**. At 600 dpi (12× magnification) the glyph is a closed
  bowl with a descending right-hand stem. It is a 9. Candidates in order: **9**, then 2.
  It is *not* a 1 — a 1 in this italic face is a bare stem with no bowl.

**But it is dimensionally wrong.** Eq. 51 (pdf p.26) reads

    M·c_pm / (h·A_m)  =  TC_T41 · √T41 / W41^(4/5)

The left side is a time constant, in sec. With `W41` in lb_m/sec and `T41` in deg R, this
requires

    TC_T41  =  lb_m^(4/5) · sec^(1/5) / deg R^(1/2)

i.e. **sec^(1/5), not sec^(9/5)**. So the report prints a `9` where its own Eq. 51
demands a `1`. The same `9/5` is printed in the nomenclature (pdf p.12), so it is a
consistent typesetting error in the report, not a scan artefact.

*This does not affect the number.* `TC_T41 = 0.29` is what the model uses; the units cell
is a label. But it is logged as an open question so nobody later "fixes" the value to make
the units work.

**No other cell in Table A.1 is ambiguous.** Every digit of every value resolves cleanly.

---

## 4. Figures A1–A11 — descriptions for the digitizing job

Common to all eleven figures:

- **No gridlines.** None of the eleven plots has a single gridline.
- Plot area is a **closed box frame with minor tick marks on all four sides**, tick labels
  on the left and bottom only. The box corners are exact and give clean registration
  points for a digitizer.
- Data are drawn as **discrete plot markers joined by segments**, not as continuous
  freehand curves. In every figure except A1 the marker is a small `×`; in A1 the marker
  is the curve's own index digit `1`–`11`.
- All eleven are line-art at high contrast; the 200 dpi renders are clean, with no
  bleed-through, no skew worth correcting, and no broken axes. Where a figure needs more
  than 200 dpi it is noted below.
- Axis titles are set in a sans-serif face and are fully legible at 200 dpi.

---

### Figure A1 — pdf p.56 (printed 42)

**Caption as printed:** `Figure A1.- T700 real-time model function f₁—compressor mass flow.`

**Relationship:** corrected compressor (station 2) mass flow as a function of compressor
static pressure ratio and corrected gas generator speed. Model use: Eq. 7, pdf p.22 —
`WA2c = f₁(Ps3/P2, NGc)`. This is the compressor map, and it is the only two-dimensional
function in the appendix.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `COMPRESSOR STATIC PRESSURE RATIO, PS3/P2` | 0 → 20 | 0, 10, 20 | ≈ every 1.0 |
| y | `STATION 2 CORRECTED MASS FLOW, WA2C, LBM/SEC` | 2 → 12 | 2, 4, 6, 8, 10, 12 | every 0.2 |

Units: y is lbm/sec; x is a dimensionless ratio.

**Parameter lines: 11.** Parameter swept is **corrected gas generator speed, `NGc`, in
percent.** An in-plot legend block (lower right quadrant, clear of all curves) gives:

```
SYMBOL  1:   65% NGC        SYMBOL  7:   92% NGC
SYMBOL  2:   80% NGC        SYMBOL  8:   94% NGC
SYMBOL  3:   82% NGC        SYMBOL  9:   96% NGC
SYMBOL  4:   85% NGC        SYMBOL 10:   98% NGC
SYMBOL  5:   87% NGC        SYMBOL 11:  100% NGC
SYMBOL  6:   89% NGC
```

Note the uneven speed spacing — a 15-point gap from 65 % to 80 %, then 2–3 points apart.

**Structure of each speed line** (this is what a digitizer has to reproduce):

1. A **long horizontal solid segment** at constant WA2c, running from Ps3/P2 ≈ 1 rightward
   to the first plotted symbol. This is the low-pressure-ratio extension, not decoration.
2. A **short sloped solid curve** carrying **6 plotted symbols**, rising to a rounded knee
   and then turning down to the right (the choke/stall side).
3. A **vertical dotted line at Ps3/P2 ≈ 1** joining the left-hand ends of all eleven
   horizontal segments.
4. **Six dotted lines** run diagonally across the whole family, each connecting the k-th
   symbol of every speed line — cross-plot/constant-index lines. They are construction
   lines, not data curves; do not digitize them as curves, but they are useful because
   they confirm the map is on a regular 11 × 6 grid.

**Data volume:** 11 × 6 = 66 plotted symbols, plus 11 horizontal-extension left endpoints
≈ **77 points**.

**Legibility at 200 dpi:** curves 7–11 (92–100 % NGc) are cleanly separated and each symbol
is individually readable. Curves 1–4 (65–85 % NGc) crowd into the lower-left corner: at
200 dpi, symbols on curves 2 and 3 touch and the `2`/`3` glyphs partially overlap.
**Digitize this figure at 600 dpi**, and expect to disambiguate curves 2–4 by following the
dotted cross-lines rather than by reading the glyphs. Curve 1 (65 %) is isolated at the
bottom and is easy.

---

### Figure A2 — pdf p.57 (printed 43)

**Caption as printed:** `Figure A2.- T700 real-time model function f₂—compressor temperature.`

**Relationship:** compressor temperature ratio vs compressor static pressure ratio. Model
use: Eq. 10, pdf p.22 — `T3 = T2 · f₂(Ps3/P2)`. Note it is a **function of pressure ratio
only** — no speed parameter, unlike `f₁`.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `COMPRESSOR STATIC PRESSURE RATIO, PS3/P2` | 0 → 20 | 0, 10, 20 | ≈ every 1.0 |
| y | `COMPRESSOR TEMPERATURE RATIO, T3/T2` | 1.1 → 3.0 | every 0.1 (1.1, 1.2 … 3.0) | 5 per 0.1 interval |

Both axes dimensionless.

**Parameter lines: 1** (no parameter swept).

**Markers:** ≈ 20 `×` markers, one at roughly each integer value of Ps3/P2 from 1 to 20.
Curve is smooth and monotonically increasing, 1.31 at the left end rising to 2.81 at x=20.

**Legibility at 200 dpi: excellent.** Every marker is separate and sits on a labelled
y-decade band. This figure could be digitized at 200 dpi. The dense y tick labelling
(every 0.1) makes vertical registration easy.

---

### Figure A3 — pdf p.58 (printed 44)

**Caption as printed:** `Figure A3.- T700 real-time model function f₃—seal-pressurization bleed fraction.`

**Relationship:** seal-pressurization bleed fraction vs corrected gas generator speed.
Model use: Eq. 12, pdf p.23 — `B1 = f₃(NGc)`.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `CORRECTED GAS GENERATOR SPEED, NGC, %` | 65 → 100 | 70, 80, 90, 100 | every 1 % |
| y | `BLEED FRACTION, B₁` | -0.02 → 0.12 | every 0.01 | 5 per 0.01 interval |

x in percent; y dimensionless. The x axis left edge is at 65 % and is **unlabelled** — the
first tick label is 70. The y axis extends to -0.02 although no data goes below 0.

**Parameter lines: 1.**

**Markers:** 13 `×` markers (this said "≈ 17"; the file carries 13). Shape: flat at ≈0.11 from 65 % to ≈78 %, a smooth shoulder,
a steep near-linear fall between ≈83 % and ≈88 %, a knee, then flat at 0.00 from ≈89 % to
100 %. Markers are concentrated in the 78–89 % transition; the two flat regions carry only
their endpoints.

**Legibility at 200 dpi: excellent.** Even on the steep segment the markers are several
pixels apart. 200–300 dpi is enough.

---

### Figure A4 — pdf p.59 (printed 45)

**Caption as printed:** `Figure A4.- T700 real-time model function f₄—power-turbine-balance bleed fraction.`

**Relationship:** power-turbine-balance-piston bleed fraction vs station 2 corrected mass
flow. Model use: Eq. 13, pdf p.23 — `B2 = f₄(WA2c)`.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `STATION 2 CORRECTED MASS FLOW, WA2_C, lbm/sec` | 3.0 → 12.0 | 4.0, 6.0, 8.0, 10.0, 12.0 | every 0.2 |
| y | `BLEED FRACTION, B₂` | 0.0088 → 0.0108 | every 0.0002 | 4 per interval |

x in lbm/sec; y dimensionless. Note the very expanded y scale — the whole variation is
0.0016.

**Parameter lines: 1.**

**Markers: exactly 3** `×` markers — one on the left axis (WA2c ≈ 3.0), one at WA2c ≈ 5.0,
one on the right axis at WA2c = 12.0. Two straight segments: a steep fall from ≈0.01057 to
0.0090, then dead flat to 12.0.

**Legibility at 200 dpi: excellent.** This is effectively a 3-breakpoint table drawn as a
graph. Digitizing risk is low; the main care needed is the y scale, where one pixel is
worth roughly 1.5×10⁻⁶.

---

### Figure A5 — pdf p.60 (printed 46)

**Caption as printed:** `Figure A5.- T700 real-time model function f₅—impeller tip leakage and turbine cooling-bleed fraction.`

**Relationship:** combined impeller-tip-leakage and gas-generator-turbine cooling bleed
fraction vs station 2 corrected mass flow. Model use: Eq. 14, pdf p.23 — `B3 = f₅(WA2c)`.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `STATION 2 CORRECTED MASS FLOW, WA2_C, lbm/sec` | 3 → 12 | 4, 6, 8, 10, 12 | every 0.2 |
| y | `BLEED FRACTION, B₃` | 0.0770 → 0.0860 | every 0.0010 | 4 per interval |

**Parameter lines: 1.**

**Markers: exactly 3** `×` markers — left axis (WA2c ≈ 3), a knee at WA2c ≈ 8.4, and the
right axis at 12. Two straight segments: a shallow fall from ≈0.0851 to ≈0.0846, then a
steep fall to ≈0.0779.

**Legibility at 200 dpi: excellent.** Another effective 3-breakpoint table. The only
judgement call is the exact x of the knee (between the 8 and 9 ticks).

---

### Figure A6 — pdf p.61 (printed 47)

**Caption as printed:** `Figure A6.- T700 real-time model function f₆—combustor efficiency.`

**Relationship:** combustor efficiency vs fuel-to-air ratio. Model use: Eq. 20, pdf p.23 —
`η = f₆(FAR)`.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `FUEL-TO-RATIO, FAR` *(printed exactly so — the word "AIR" is missing; a typo in the report)* | 0.010 → 0.020 | 0.010, 0.012, 0.014, 0.016, 0.018, 0.020 | every 0.0002 |
| y | `COMBUSTOR EFFICIENCY, ETA` | 0.88 → 1.10 | every 0.02 | 4 per interval |

Both dimensionless.

**Parameter lines: 1.**

**Markers: exactly 2** `×` markers, one on each vertical axis. The line between them is
**perfectly horizontal** — combustor efficiency is modelled as a constant, independent of
FAR, over the whole plotted range. The level sits between the 0.98 and 1.00 tick labels,
nearer 0.98.

**Legibility at 200 dpi: excellent.** Digitizing this reduces to reading one y value to
three or four significant figures; it is worth doing at 600 dpi purely because the whole
function is that single number, and the y scale (0.02 per labelled interval, ~90 px at
200 dpi) means a one-pixel error is ~0.0002.

---

### Figure A7 — pdf p.62 (printed 48)

**Caption as printed:** `Figure A7.- T700 real-time model function f₇—gas-generator turbine energy.`

**Relationship:** gas-generator turbine enthalpy-drop parameter vs gas-generator turbine
pressure ratio. Model use: Eq. 26, pdf p.24 — `ΔH_GT = θ₄₁ · f₇(P45/P41)`.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `GAS GENERATOR TURBINE PRESSURE RATIO, P45/P41` | ≈0.20 → 0.35 | 0.21, 0.23, 0.25, 0.27, 0.29, 0.31, 0.33, 0.35 | every 0.005 |
| y | `GAS GENERATOR ENTHALPY DROP PARAMETER, BTU/LBM` | ≈20 → 52.5 | 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5 | 4 per interval |

y in Btu/lbm; x dimensionless. **The bottom of the y frame is below the 22.5 label and is
unlabelled** (≈20) — do not assume the frame bottom is 22.5 when setting up the
transformation. Use the labelled ticks, not the frame, for registration.

**Parameter lines: 1.**

**Markers: exactly 6** `×` markers. Five of them are bunched between x ≈ 0.200 and
x ≈ 0.232 on a sharply-curving knee (values falling ≈48 → ≈37.5); the sixth sits alone at
x = 0.35, y ≈ 24. The segment from the fifth marker to the sixth is a single long straight
line covering three-quarters of the plot width.

**Legibility at 200 dpi: adequate for the shape, poor for the knee.** The five left-hand
markers span roughly 90 px at 200 dpi and two of them nearly touch.
**Digitize this figure at 600 dpi.** It matters: with only six points, every one carries
one-sixth of the function.

---

### Figure A8 — pdf p.63 (printed 49)

**Caption as printed:** `Figure A8.- T700 real-time model function f₈—power turbine energy.`

**Relationship:** power turbine enthalpy-drop parameter vs power turbine pressure ratio.
Model use: Eq. 32, pdf p.24 — `ΔH_PT = θ₄₅ · f₈(P49/P45)`.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `POWER TURBINE PRESSURE RATIO, P49/P45` | 0.30 → 0.85 | 0.35, 0.45, 0.55, 0.65, 0.75, 0.85 | every 0.01 |
| y | `POWER TURBINE ENTHALPY DROP PARAMETER, BTU/LBM` | -5.0 → 40.0 | every 2.5 | 4 per interval |

y in Btu/lbm; x dimensionless. Note the x axis starts at 0.30 and the leftmost label is
0.35 — the left frame edge is unlabelled.

**Parameter lines: 1.**

**Markers: 12** `×` markers, evenly spaced at Δ(P49/P45) = 0.05 from 0.30 to 0.85. Curve is
near-linear with a slight convexity, falling from ≈34.2 at 0.30 to ≈-1.2 at 0.85. It
crosses zero at ≈0.83.

**Legibility at 200 dpi: excellent.** Well-separated markers, dense y tick labelling.
200–300 dpi is enough.

**Note for the model:** `f₈` goes **negative** above P49/P45 ≈ 0.83. That is real, not a
plotting artefact — the y axis is deliberately extended to -5.0 to show it.

---

### Figure A9 — pdf p.64 (printed 50)

**Caption as printed:** `Figure A9.- T700 real-time model function f₉—power turbine mass flow.`

**Relationship:** power turbine corrected mass flow vs power turbine pressure ratio.
Model use: Eq. 33, pdf p.24 — `W45c = f₉(Ps9/P45)`.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `POWER TURBINE PRESSURE RATIO, PS9/P45` | 0.30 → 0.85 | 0.35, 0.45, 0.55, 0.65, 0.75, 0.85 | every 0.01 |
| y | `POWER TURBINE CORRECTED MASS FLOW, W45C, LBM/SEC` | 0.25 → 0.39 | every 0.02 (0.25, 0.27 … 0.39) | 4 per interval |

y in lbm/sec; x dimensionless.

**Careful:** the independent variable here is `Ps9/P45` — station 9 **static** over station
4.5 total — whereas Figure A8 on the facing page uses `P49/P45`, total over total. Both
were checked against the equations at 300 dpi: Eq. 32 prints `P49/P45`, Eq. 33 prints
`Ps9/P45`. The two figures are **not** on the same abscissa despite identical numeric
ranges and near-identical axis titles. Do not cross-index them.

**Parameter lines: 1.**

**Markers: ≈ 23** `×` markers, evenly spaced at Δ ≈ 0.025 from 0.30 to 0.85. Shape is a
flat top (choked) at ≈0.372 out to ≈0.45, then an accelerating fall to ≈0.259 at 0.85.

**Legibility at 200 dpi: good.** On the flat left-hand portion consecutive markers nearly
touch (they are ~25 px apart with a marker ~10 px across), but they never merge.
**300 dpi recommended**, 600 dpi if the flat region turns out to matter.

---

### Figure A10 — pdf p.65 (printed 51)

**Caption as printed:** `Figure A10.- T700 real-time model function f₁₀—exhaust pressure loss.`

**Relationship:** exhaust pressure ratio vs corrected gas generator speed. Model use:
Eq. 38, pdf p.24 — `P49 = Ps9 · f₁₀(NGc)`.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `CORRECTED GAS GENERATOR SPEED, NGC, %` | 65 → 100 | 65, 70, 75, 80, 85, 90, 95, 100 | every 1 % |
| y | `EXHAUST PRESSURE RATIO, PS9/P49` | ≈1.012 → 1.14 | 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14 | 4 per interval |

x in percent; y dimensionless. **The bottom of the y frame is below the 1.02 label and is
unlabelled** (≈1.012) — register on labelled ticks, not the frame.

**Parameter lines: 1.**

**Markers: ≈ 20** `×` markers, **unevenly spaced** — sparse on the long fall from 65 % to
85 %, then closely packed through the minimum at 90–92 %, then sparse again on the rise to
100 %. Shape: falls from ≈1.123 at 65 % to a shallow minimum of ≈1.020 at ≈90–91 %, then
rises to ≈1.065 at 100 %.

**Legibility at 200 dpi: good, except at the minimum.** Three or four markers between 89 %
and 92 % sit within a few pixels of each other and of the curve's own line weight.
**Digitize at 600 dpi**, and expect the marker count in the trough to be a judgement call.

**⚠ Axis-label / equation inconsistency.** The y axis is labelled `PS9/P49`, but the
plotted values are all > 1 and Eq. 38 uses `f₁₀` as a **multiplier on Ps9 to obtain P49**,
which requires `f₁₀ = P49/Ps9`. Taken literally the axis label is inverted with respect to
the equation. Both the label (600 dpi) and the equation (300 dpi) were re-read and both are
as stated. Logged as an open question — this must be settled before `f₁₀` is wired in,
because using the label as printed inverts the exhaust pressure loss.

---

### Figure A11 — pdf p.66 (printed 52)

**Caption as printed:** `Figure A11.- T700 real-time model function f_hs—station 4.1 heat-sink constant.`

**Relationship:** the station 4.1 heat-sink constant vs corrected gas generator speed.
Model use: Eq. 52, pdf p.26 — `T41sgn = f_hs(NGc)`, which feeds the heat-sink transfer
function of Eq. 50 via Eq. 53.

**Axes**

| | Label as printed | Range | Tick labels | Minor ticks |
|---|---|---|---|---|
| x | `CORRECTED GAS GENERATOR SPEED, NGC, %` | 65 → 100 | 65, 70, 75, 80, 85, 90, 95, 100 | every 1 % |
| y | `HEAT-SINK CONSTANT, T41SGN` | 3.50 → 8.00 | every 0.50 | 4 per interval |

x in percent. y units are not stated on the figure and `T41_sgn` is given no units in the
nomenclature either (pdf p.12 lists it as "heat-sink function" with no unit).

**Parameter lines: 1.**

**Markers: exactly 6** `×` markers, at NGc = **65, 70, 75, 80, 85 and 100 %** — five on
the 5 %-grid then one long flat run from 85 % to 100 %. Shape: flat at 4.00 from 65 to
70 %, a steep rise to ≈5.75 at 75 %, ≈6.20 at 80 %, ≈7.49 at 85 %, then dead flat to 100 %.
Markers land exactly on labelled x ticks, which makes the abscissa exact and reduces the
job to reading six ordinates.

**Legibility at 200 dpi: excellent.** Effectively a 6-breakpoint table. This is the
easiest figure in the appendix to digitize.

---

## 5. What Appendix A does not contain

Appendix A's own body text is a single sentence and makes no cross-reference, so there is
nothing it cites and fails to deliver. What follows is the complementary list: things the
engine model needs and Appendix A does **not** supply. All were checked against the
component equations on pdf pp.22–26.

- **Numeric values for `f₁`–`f₁₀` and `f_hs`.** Not printed anywhere. This is the finding
  above, not a side note: the appendix supplies these eleven relationships *only* as ink.
- **Breakpoint (independent-variable) values for every function.** Even the x locations of
  the plotted markers have to be recovered by digitizing. The one exception is Figure A11,
  whose markers land on labelled ticks.
- **`J_load`** — Eq. 46, pdf p.25: `J = J_PT + J_load`. Table A.1 gives `J_PT` but not
  `J_load`, the helicopter drivetrain inertia reflected to the power turbine. It is an
  external (Gen Hel) quantity, consistent with `Q_req` being an input boundary condition.
- **Ambient / atmosphere model** (`P_amb`, `T_amb`) — inputs, not model constants.
- **Fuel control system constants** — those are in Appendix C (pdf pp.77–101), not here.
- **`778.12`** (ft·lbf per Btu) appears inline in Eqs. 39–41 rather than in Table A.1. It
  is a unit conversion, so it belongs in `t700.units`, not in the constants module — but
  note that the report's value is 778.12, and cite the equation page, not a handbook.

Every other symbol appearing in the component equations of pp.22–26 resolves either to a
Table A.1 constant or to an Appendix A figure. Specifically checked: `K_QC₁`/`K_QC₂`
(Eq. 39), `K_damp`/`NP_des` (Eq. 41), `K_V3`/`K_V41`/`K_V45` (Eqs. 42–44), `K_WGT`
(Eq. 28), `K_dpb` (Eq. 18), `TC_T41` (Eq. 51) — all present.

---

## 6. Provenance summary

- Table A.1: **transcribed** from `docs/extracted/p-055.png` (200 dpi), verified against
  300 dpi and 600 dpi re-renders. 33 rows, 1 flagged cell (`TC_T41` units, item A-1).
- Figures A1–A11: **not digitized here.** Only captions, axis labels, ranges, tick
  structure, parameter lines and marker counts were read. No data point was read off any
  curve — that is a separate, deliberate step with its own tooling.
- Citation form for anything taken from this appendix:
  `[TM-100991 pdf p.55, Table A.1]` for constants,
  `[TM-100991 pdf p.NN, Fig. ANN]` plus `digitized` for map data.
