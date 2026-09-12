# Symbol table

Transcribed from the report's NOMENCLATURE section, **PDF pp.9–13** (printed pages iii–vii).
PDF pp.8 and 14 are blank; p.7 is the draft title page and p.15 begins the Summary.

**Status: complete.** All 106 entries transcribed. Every row was read from a page image
rendered at 300 dpi (600 dpi and pixel-level zoom for the ambiguous cells); nothing in this
file comes from the OCR text layer.

Rules:
- Every row cites the PDF page it came from.
- Every row is read from a rendered page image, not the OCR text layer.
- Python name is the report symbol lowercased, with a unit suffix where the unit is not
  obvious. `NP` → `np_` (trailing underscore: `np` is NumPy).

## How the report typesets symbols (read this before using the table)

The report distinguishes two kinds of trailing digits, and the distinction is *not*
recoverable from the OCR text layer:

- **Station indices are set full size, inline with the letter**: `H41`, `P45`, `T41`,
  `W41`, `WA2`, `WA31`. There is no subscript there — `T41` means "station 4.1
  temperature", written as T followed by 4 and 1.
- **True subscripts are set small and lowered**: `A_m`, `B_1`, `H41_ns`, `P_amb`, `Q_req`,
  `δ_2`, `θ_41`.
- **A few constants carry a subscript on a subscript** — K with subscript `H3` and a
  further sub-subscript `1`. Written below with a Unicode subscript digit for the second
  level: `K_H3₁`, `K_H3₂`, `K_H41₁`, `K_H41₂`, `K_QC₁`, `K_QC₂`, `K_T41₁`, `K_T41₂`,
  `K_T45₁`, `K_T45₂`, `K_T49₁`, `K_T49₂`, `K_TH41₁`, `K_TH41₂`, `K_TH45₁`, `K_TH45₂`.
  **These are 16 distinct constants, not 8.** `K_H41₁` and `K_H41₂` have different units
  (Btu/lb_m·deg R and Btu/lb_m) and are the slope and intercept of one linear fit; the same
  pattern holds for the other pairs. Flattening `K_H41₁` to `KH41` loses information and
  will silently corrupt the enthalpy fits.

Unit-suffix policy applied in the Python-name column:
- `deg R` → `_degR`; `lb_m/sec` → `_pps`; `lb_f/in²` → `_psia`; `ft·lb_f` → `_ftlbf`;
  `ft·lb_f·sec²` → `_ftlbfs2`; `Btu/lb_m` → `_btulbm`; `rpm` → `_rpm`; `sec` → `_sec`.
- Dimensionless quantities (fractions, ratios, efficiencies) get no suffix.
- Empirical constants with compound units (the `K…` family, `TC_T41`) get **no** suffix —
  a suffix for `lb_f²·sec²/lb_m²·in⁴·deg R` would be unreadable. Their units live in the
  Units column here and in the Appendix A citation at the point of definition.
- Where the nomenclature prints **no** unit, the Python name gets no unit suffix and the
  Units cell says *not given*. Do not invent one; confirm it from the defining equation
  before adding a suffix.

## Grouping

The report's nomenclature is **one single alphabetical list** — there is no separate
subscript, acronym, or superscript section. Latin symbols run A→W, then Greek symbols
δ, ΔH, η, τ, θ. The tables below preserve the report's own printed order, page by page.
Two order oddities are the report's, not transcription errors:
`WA24_bl` is printed *before* `WA2_c`, and `τ` is printed *before* `θ`.

## Latin symbols

| Report symbol | Python name | Meaning | Units | Source |
|---|---|---|---|---|
| A_m | `am` | area of exposed metal; used in heat-sink representation | *not given* | pdf p.9 |
| B_1 | `b1` | seal-pressurization bleed fraction | – (fraction) | pdf p.9 |
| B_2 | `b2` | power-turbine-balance bleed fraction | – (fraction) | pdf p.9 |
| B_3 | `b3` | compressor-diffuser bleed fraction | – (fraction) | pdf p.9 |
| c_pg | `cpg` | specific heat of gas; used in heat-sink representation | *not given* | pdf p.9 |
| c_pm | `cpm` | specific heat of metal; used in heat-sink representation | *not given* | pdf p.9 |
| FAR | `far` | ratio of fuel to atmospheric gas in combustor | – (ratio) | pdf p.9 |
| h | `h_conv` | convective heat transfer coefficient; used in heat-sink representation | *not given* | pdf p.9 |
| H2 | `h2_btulbm` | station 2 enthalpy | Btu/lb_m | pdf p.9 |
| H3 | `h3_btulbm` | station 3 enthalpy | Btu/lb_m | pdf p.9 |
| H41 | `h41_btulbm` | station 4.1 enthalpy | Btu/lb_m | pdf p.9 |
| H41_ns | `h41ns_btulbm` | station 4.1 enthalpy **not** including heat-sink effects | Btu/lb_m | pdf p.9 |
| H44 | `h44_btulbm` | station 4.4 enthalpy | Btu/lb_m | pdf p.9 |
| H45 | `h45_btulbm` | station 4.5 enthalpy | Btu/lb_m | pdf p.9 |
| H49 | `h49_btulbm` | station 4.9 enthalpy | Btu/lb_m | pdf p.9 |
| HVF | `hvf_btulbm` | heating value of fuel | Btu/lb_m | pdf p.9 |
| J | `j_ftlbfs2` | moment of inertia of all mass rigidly attached to the engine output shaft | ft·lb_f·sec² | pdf p.9 |
| J_GT | `jgt_ftlbfs2` | moment of inertia of the rigid mass which represents the gas generator, compressor, and the associated shafting | ft·lb_f·sec² | pdf p.9 |
| J_load | `jload_ftlbfs2` | moment of inertia of the load mass | ft·lb_f·sec² | pdf p.9 |
| J_PT | `jpt_ftlbfs2` | moment of inertia of the power turbine and output shaft | ft·lb_f·sec² | pdf p.9 |
| K_b3 | `kb3` | station 3 bleed-flow coefficient | *not given* | pdf p.9 |
| K_bl | `kbl` | fraction of diffuser bleed gas which is used to cool the gas generator turbine blades | – (fraction) | pdf p.10 |
| K_damp | `kdamp` | power turbine speed damping coefficient | ft·lb_f·sec/rad | pdf p.10 |
| K_dpb | `kdpb` | combustor pressure drop coefficient | lb_f²·sec²/lb_m²·in⁴·deg R | pdf p.10 |
| K_H2 | `kh2` | station 2 enthalpy coefficient | Btu/lb_m·deg R | pdf p.10 |
| K_H3₁ | `kh31` | station 3 enthalpy coefficient | Btu/lb_m·deg R | pdf p.10 |
| K_H3₂ | `kh32` | station 3 enthalpy coefficient | Btu/lb_m | pdf p.10 |
| K_H41₁ | `kh411` | station 4.1 enthalpy coefficient | Btu/lb_m·deg R | pdf p.10 |
| K_H41₂ | `kh412` | station 4.1 enthalpy coefficient | Btu/lb_m | pdf p.10 |
| K_H45 | `kh45` | fraction of station 4.4 enthalpy used to represent station 4.5 enthalpy | – (fraction) | pdf p.10 |
| K_Ps3 | `kps3` | fraction of total pressure at station 3 used to represent static pressure at station 3 | – (fraction) | pdf p.10 |
| K_QC₁ | `kqc1` | empirically-determined compressor torque coefficient | *not given* | pdf p.10 |
| K_QC₂ | `kqc2` | empirically-determined compressor torque coefficient | *not given* | pdf p.10 |
| K_T41₁ | `kt411` | station 4.1 temperature coefficient | lb_m·deg R/Btu | pdf p.10 |
| K_T41₂ | `kt412` | station 4.1 temperature coefficient | deg R | pdf p.10 |
| K_T45₁ | `kt451` | station 4.5 temperature coefficient | lb_m·deg R/Btu | pdf p.10 |
| K_T45₂ | `kt452` | station 4.5 temperature coefficient | deg R | pdf p.10 |
| K_T49₁ | `kt491` | station 4.9 temperature coefficient | lb_m·deg R/Btu | pdf p.10 |
| K_T49₂ | `kt492` | station 4.9 temperature coefficient | deg R | pdf p.10 |
| K_TH41₁ | `kth411` | station 4.1 velocity ratio coefficient | 1/deg R | pdf p.10 |
| K_TH41₂ | `kth412` | station 4.1 velocity ratio coefficient | – (none printed) | pdf p.10 |
| K_TH45₁ | `kth451` | station 4.5 velocity ratio coefficient | 1/deg R | pdf p.10 |
| K_TH45₂ | `kth452` | station 4.5 velocity ratio coefficient | – (none printed) | pdf p.11 |
| K_V3 | `kv3` | station 3 volume coefficient | lb_f/in²·lb_m·deg R | pdf p.11 |
| K_V41 | `kv41` | station 4.1 volume coefficient | lb_f/in²·lb_m·deg R | pdf p.11 |
| K_V45 | `kv45` | station 4.5 volume coefficient | lb_f/in²·lb_m·deg R | pdf p.11 |
| K_WGT | `kwgt` | station 4.1 flow coefficient | lb_m·in²/lb_f·sec | pdf p.11 |
| M | `m_metal` | mass of metal that absorbs heat energy from gas; used in heat-sink representation | *not given* | pdf p.11 |
| NG | `ng_rpm` | rotational speed of compressor and gas generator | rpm | pdf p.11 |
| NG_c | `ngc_rpm` | corrected compressor and gas generator speed; a nonphysical value which is independent of inlet conditions | rpm | pdf p.11 |
| NG_des | `ngdes_rpm` | design rotational speed of the gas generator and compressor | rpm | pdf p.11 |
| NP | `np_rpm` (bare: `np_`) | rotational speed of power turbine and output shaft | rpm | pdf p.11 |
| NP_des | `npdes_rpm` | design rotational speed of power turbine and output shaft | rpm | pdf p.11 |
| P1 | `p1_psia` | station 1 total pressure | lb_f/in² | pdf p.11 |
| P2 | `p2_psia` | station 2 total pressure | lb_f/in² | pdf p.11 |
| P3 | `p3_psia` | station 3 total pressure | lb_f/in² | pdf p.11 |
| P41 | `p41_psia` | station 4.1 total pressure | lb_f/in² | pdf p.11 |
| P45 | `p45_psia` | station 4.5 total pressure | lb_f/in² | pdf p.11 |
| P49 | `p49_psia` | station 4.9 total pressure | lb_f/in² | pdf p.11 |
| P_amb | `pamb_psia` | ambient pressure | lb_f/in² | pdf p.11 |
| P_s3 | `ps3_psia` | station 3 static pressure; a fuel-control-system parameter | lb_f/in² | pdf p.11 |
| P_s9 | `ps9_psia` | station 9 static pressure | lb_f/in² | pdf p.11 |
| P_std | `pstd_psia` | standard-day pressure | lb_f/in² | pdf p.11 |
| Q_acc | `qacc_ftlbf` | torque required for helicopter accessory power | ft·lb_f | pdf p.11 |
| Q_C | `qc_ftlbf` | torque required by compressor | ft·lb_f | pdf p.12 |
| Q_damp | `qdamp_ftlbf` | torque required by helicopter gearbox (gearbox damping) | ft·lb_f | pdf p.12 |
| Q_eng | `qeng_ftlbf` | torque transmitted to helicopter gearbox from both engines | ft·lb_f | pdf p.12 |
| Q_GT | `qgt_ftlbf` | torque output of gas generator | ft·lb_f | pdf p.12 |
| Q_mr | `qmr_ftlbf` | torque required by helicopter main rotor | ft·lb_f | pdf p.12 |
| Q_PT | `qpt_ftlbf` | torque output of power turbine | ft·lb_f | pdf p.12 |
| Q_req | `qreq_ftlbf` | torque required by external load with respect to power turbine speed | ft·lb_f | pdf p.12 |
| Q_tr | `qtr_ftlbf` | torque required by helicopter tail rotor | ft·lb_f | pdf p.12 |
| T1 | `t1_degR` | station 1 temperature | deg R | pdf p.12 |
| T2 | `t2_degR` | station 2 temperature | deg R | pdf p.12 |
| T3 | `t3_degR` | station 3 temperature | deg R | pdf p.12 |
| T41 | `t41_degR` | station 4.1 temperature | deg R | pdf p.12 |
| T41_ns | `t41ns_degR` | station 4.1 temperature **not** including heat-sink effects | deg R | pdf p.12 |
| T41_sgn | `t41sgn` | heat-sink function | *not given* | pdf p.12 |
| T45 | `t45_degR` | station 4.5 temperature | deg R | pdf p.12 |
| T49 | `t49_degR` | station 4.9 temperature | deg R | pdf p.12 |
| T_amb | `tamb_degR` | ambient temperature | deg R | pdf p.12 |
| TC_T41 | `tct41` | empirically determined constant used in station 4.1 heat-sink representation | lb_m^(4/5)·sec^(9/5)/deg R^(1/2) — see note **[A]** | pdf p.12 |
| T_gi | `tgi` | temperature of gas entering the heat-sink representation control volume | *not given* | pdf p.12 |
| T_go | `tgo` | temperature of gas leaving the heat-sink representation control volume | *not given* | pdf p.12 |
| T_m | `tm` | temperature of exposed metal as used in heat-sink representation | *not given* | pdf p.12 |
| T_std | `tstd_degR` | standard atmosphere temperature at sea-level | deg R | pdf p.13 |
| W41 | `w41_pps` | station 4.1 mass flow rate of combustion gases | lb_m/sec | pdf p.13 |
| W45 | `w45_pps` | station 4.5 mass flow rate of combustion gases | lb_m/sec | pdf p.13 |
| W45_c | `w45c_pps` | corrected station 4.5 mass flow rate of combustion gases | lb_m/sec | pdf p.13 |
| WA2 | `wa2_pps` | station 2 mass-flow rate of atmospheric gas | lb_m/sec | pdf p.13 |
| WA24_bl | `wa24bl_pps` | station 2.4 bleed flow of atmospheric gas | lb_m/sec | pdf p.13 |
| WA2_c | `wa2c_pps` | corrected station 2 mass-flow rate of atmospheric gas | lb_m/sec | pdf p.13 |
| WA3 | `wa3_pps` | station 3 mass flow rate of atmospheric gas | lb_m/sec | pdf p.13 |
| WA31 | `wa31_pps` | station 3.1 mass flow rate of atmospheric gas | lb_m/sec | pdf p.13 |
| WA3_bl | `wa3bl_pps` | diffuser bleed discharge flow | lb_m/sec | pdf p.13 |
| W_f | `wf_pps` | fuel flow | lb_m/sec | pdf p.13 |
| W_g | `wg` | mass flow rate of gas into heat-sink representation control volume | *not given* | pdf p.13 |

## Greek symbols

Printed after the Latin list, in the report's own order (note τ precedes θ).

| Report symbol | Python name | Meaning | Units | Source |
|---|---|---|---|---|
| δ_2 | `delta2` | ratio of inlet pressure to sea-level pressure | – (ratio) | pdf p.13 |
| ΔH_GT | `dhgt_btulbm` | gas-generator-turbine enthalpy drop | Btu/lb_m | pdf p.13 |
| ΔH_PT | `dhpt_btulbm` | power turbine enthalpy drop | Btu/lb_m | pdf p.13 |
| η | `eta` | combustor efficiency | – (fraction) | pdf p.13 |
| τ_1 | `tau1_sec` | heat-sink **lead** time constant for a trim operating condition; used in small-perturbation representation | sec | pdf p.13 |
| τ_2 | `tau2_sec` | heat-sink **lag** time constant for a trim operating condition; used in small-perturbation representation | sec | pdf p.13 |
| θ_2 | `theta2` | ratio of inlet temperature to standard-day temperature | – (ratio) | pdf p.13 |
| θ_41 | `theta41` | station 4.1 squared-critical-velocity ratio | – (ratio) | pdf p.13 |
| θ_45 | `theta45` | station 4.5 squared-critical-velocity ratio | – (ratio) | pdf p.13 |

**Total: 106 entries** (p.9: 21, p.10: 21, p.11: 22, p.12: 21, p.13: 21).

## Marked cells

**[A] `TC_T41` units — partially ambiguous.** Printed as stacked fractional exponents at a
size the 1988 phototypesetting and the 2009 scan barely resolve. Read at 300 dpi, re-read at
600 dpi, and inspected pixel-by-pixel:

- `lb_m` exponent: **4/5** — confident. The 4 and the 5 are both unambiguous at pixel level.
- `deg R` exponent: **1/2** — confident.
- `sec` exponent numerator: **9** — read as a closed bowl with a descending right stem, so
  9 rather than 2 or 0. Confidence: good but not certain. Candidates, in order: **9**, then 2.
  The denominator matches the `lb_m` denominator glyph exactly, so it is **5**.

So: `lb_m^(4/5) · sec^(9/5) / deg R^(1/2)`. **Settled 2026-09-12 (open question #4): the
glyph is unambiguously a 9 at 600 dpi, and the printed units are therefore the report's own
defect, not a scan artefact.** Eq. 51 [pdf p.26] reads `M c_pm/(h A_m) = TC_T41·sqrt(T41) /
W41^(4/5)`, whose left side is a time, so `TC_T41` must carry `sec^(1/5)`; the printed
`sec^(9/5)` cannot be made consistent with the equation that uses it at any value of the
constant. The heat-sink equations are pdf **pp.25-26** (Eqs. 48-49 on p.25, Eqs. 50-53 on
p.26) -- this note cited pp.23-24. The value 0.29 is unaffected. No other cell in the
nomenclature was ambiguous.

## Corrections to earlier placeholder rows in this file

The previous OCR-derived draft of this table contained one wrong symbol and one wrong
Python name. Both are fixed above:

- **`Qgb` does not exist.** The symbol for "torque transmitted to helicopter gearbox from
  both engines" is **`Q_eng`** [pdf p.12]. Anything already written against `qgb` must be
  renamed `qeng_ftlbf`.
- **`Qgg` does not exist.** "Torque output of gas generator" is **`Q_GT`** [pdf p.12] →
  `qgt_ftlbf`.
- `Qc` is printed **`Q_C`** (capital C subscript) [pdf p.12].

## Name-collision notes

**Decided, and CLAUDE.md's naming table now carries all four.** This section was written
before Phase 1 as a recommendation; it is kept for the reasoning, and the table rows above
use the decided names. Transcribing the full list turned up three single-letter names beyond
`NP` → `np_`:

- `h` (convective heat transfer coefficient) sits one character from `h2`, `h3`, `h41`
  (enthalpies). A local `h` in a thermo routine will read as an enthalpy to anyone skimming.
- `m` (mass of heat-sink metal) is the conventional loop/mass variable name.
- `j` (moment of inertia) is the conventional loop index.

Resolved as `h_conv`, `m_metal`, and inertias always suffixed so bare `j` never appears
(the code uses `J_GT`, `J_PT`, `J_LOAD_UH60A`). `h_conv` and `m_metal` do not appear in
`src/` because the quantities themselves never do: Eqs. 51 and 53 give their *groupings*
(`M c_pm/(h A_m)` and `M c_pm/(W_g c_pg)`) directly as time constants, so the individual
factors are never needed.

## Station numbering

Ballin's station scheme is **not** the conventional 0/2/25/3/4/45/5 turbofan numbering, and
it must not be silently mapped onto it. The complete set of stations, from Figure 3
[pdf p.21] and the paragraph that describes it [pdf p.20]:

| Station | What the report says it is | Source | Symbols carrying this index |
|---|---|---|---|
| 1 | Engine inlet. Fig. 3 labels station 1 at the engine inlet; the report gives no verbal definition of it. Eqs. (1)–(2) set P2 = P1 = P_amb and T2 = T1 = T_amb, so station 1 carries ambient conditions. | pdf p.21 (Fig. 3), pdf p.22 (Eqs. 1–2) | P1, T1 |
| 2 | Compressor inlet. The six compressor stages and variable-geometry flow vanes are represented "between stations 2 and 3". P2 and T2 are estimated equal to static ambient conditions, because stagnation effects roughly offset aircraft inlet losses over most of the flight envelope. | pdf p.20, pdf p.22 | P2, T2, H2, WA2, WA2_c, K_H2, δ_2, θ_2 |
| 2.4 | Stage 4 of the compressor. Bleed flow used for seal pressurization and power turbine balance is extracted here. | pdf p.20 | WA24_bl |
| 3 | The compressor diffuser, or outlet. Flow is bled at this point to cool the combustor and gas generator turbine. | pdf p.20 | P3, P_s3, T3, H3, WA3, WA3_bl, B_3, K_b3, K_V3, K_H3₁, K_H3₂, K_Ps3 |
| 3.1 | The combustor. | pdf p.20 | WA31 |
| 4.1 | The mixing volume representing the combustor outlet and gas-generator-turbine inlet. | pdf p.20 | P41, T41, T41_ns, T41_sgn, H41, H41_ns, W41, θ_41, K_V41, K_WGT, K_H41₁, K_H41₂, K_T41₁, K_T41₂, K_TH41₁, K_TH41₂, TC_T41 |
| 4.4 | The thermodynamic state of output gases passing through the gas generator turbine, **not** including the effects of the cooling bleed flow. | pdf p.20 | H44 (and K_H45, the fraction of station 4.4 enthalpy used to represent station 4.5) |
| 4.5 | The power turbine inlet. The cooling-bleed effects excluded at 4.4 are added here. | pdf p.20 | P45, T45, H45, W45, W45_c, θ_45, K_V45, K_H45, K_T45₁, K_T45₂, K_TH45₁, K_TH45₂ |
| 4.9 | The power turbine outlet. | pdf p.20 | P49, T49, H49, K_T49₁, K_T49₂ |
| 9 | The engine exhaust. | pdf p.20 | P_s9 |

Points worth holding on to:

- There is **no station 4** and **no station 5**. The gas-generator turbine exit is 4.4,
  the power turbine inlet is 4.5, the power turbine exit is 4.9, and the exhaust is 9.
  `T45` is therefore the power turbine *inlet* temperature, not a "T4.5" in any other
  engine's sense.
- Stations 4.4 and 4.5 differ **only** by the cooling bleed mixing; H45 is derived from H44
  through the fraction K_H45.
- Not every station carries a full state. Station 3.1 (the combustor) appears only as a mass
  flow, WA31. Station 4.4 appears only as an enthalpy, H44. Station 9 appears only as a
  static pressure, P_s9. Station 1 appears only as P1 and T1, both set to ambient.
- Figure 3 [pdf p.21] is the authority for the layout: it labels, left to right,
  STATION 1, 2, 2.4, 3, 3.1, 4.1, 4.4, 4.5, 4.9, 9, with WF (fuel flow) entering at the
  combustor, and names the four components COMPRESSOR, COMBUSTOR, GAS GENERATOR,
  POWER TURBINE.
