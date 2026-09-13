# Inventory — Appendix C, T700 Fuel Control System Model

Source: Ballin, *A High-Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
NASA TM-100991, 1988. **PDF pages 77–101** (printed pages 63–87; offset = 14).

Method: every value below was read from the rendered page image (200 dpi, re-rendered at
300 or 600 dpi where digits or block-diagram signs were tight). Nothing here comes from
the text layer. Where the text layer disagrees with the image, the image wins and the
discrepancy is noted.

**Headline finding: Appendix C contains no numbered equations and no printed function
tables.** It is 6 pages of nomenclature, one 2-page constants table (Table C.1), and
**30 figures** (C1–C30): 22 block diagrams and 8 function plots. The control system's
"equations" exist only as block-diagram transfer functions and nonlinearity icons, and
the eight scheduling functions (F_EC1, F_HM1…F_HM7) exist only as **plots that must be
digitized** — no breakpoint tables are printed anywhere in the appendix.

---

## 1. Page-by-page map

| PDF p. | Printed p. | Content |
|---|---|---|
| 77 | 63 | Appendix C title; one-line intro; "T700 Fuel Control System Variables and Constants" begins — AWFP … CT16 |
| 78 | 64 | Nomenclature CTPL … PS3HYS |
| 79 | 65 | Nomenclature PS3L … TL2 |
| 80 | 66 | Nomenclature TLGE … WFPAC |
| 81 | 67 | Nomenclature WFPDC … ZK3 |
| 82 | 68 | Nomenclature ZK5 … ZLOLIM (page otherwise blank) |
| 83 | 69 | **Table C.1** — T700 Fuel Control System Constants (AWFP … T17), 28 rows |
| 84 | 70 | **Table C.1 Concluded** (T45COR … ZLOLIM), 29 rows |
| 85 | 71 | Fig. C1 — T700 electrical control unit (ECU), top level; Fig. C2 — ECU load-share speed trim |
| 86 | 72 | Fig. C3 — ECU governor rate compensation; Fig. C4 — ECU governor dynamics |
| 87 | 73 | Fig. C5 — ECU thermocouple harness; Fig. C6 — T45 compensation |
| 88 | 74 | Fig. C7 — ECU proportional-plus-integral (P+I) compensation; Fig. C8 — P+I integrator limit logic |
| 89 | 75 | Fig. C9 — T700 hydromechanical control unit (HMU), top level; Fig. C10 — torque motor dynamics and compensation |
| 90 | 76 | Fig. C11 — compressor static discharge pressure sensor; Fig. C12 — load demand (fuel flow command) dynamics; Fig. C13 — collective pitch to load demand spindle rigging (UH-60A) |
| 91 | 77 | Fig. C14 — topping line schedule; Fig. C15 — NG spool sensor dynamics; Fig. C16 — power available spindle input schedule; Fig. C17 — load demand spindle input schedule |
| 92 | 78 | Fig. C18 — limit selector logic; Fig. C19 — fuel flow metering valve |
| 93 | 79 | Fig. C20 — idle schedule; Fig. C21 — acceleration schedule; Fig. C22 — deceleration schedule |
| 94 | 80 | Fig. C23 — function F_EC1, ECU thermocouple sensor time constant (plot, 4 curves) |
| 95 | 81 | Fig. C24 — function F_HM1, HMU topping line schedule (plot) |
| 96 | 82 | Fig. C25 — function F_HM2, HMU power-available input schedule (plot) |
| 97 | 83 | Fig. C26 — function F_HM3, HMU load-demand-compensation schedule function 1 (plot) |
| 98 | 84 | Fig. C27 — function F_HM4, HMU load-demand-compensation schedule function 2 (plot) |
| 99 | 85 | Fig. C28 — function F_HM5, HMU idle-schedule function 1 (plot) |
| 100 | 86 | Fig. C29 — function F_HM6, HMU idle-schedule function 2 (plot) |
| 101 | 87 | Fig. C30 — function F_HM7, HMU maximum fuel-parameter limit during acceleration (plot, 7 curves) |

---

## 2. Control architecture as the report presents it

The report splits the fuel control into **two units**, named exactly so:

- **ECU — "T700 electrical control unit"** [Fig. C1, p.85]. Governs *power turbine speed*
  NP and limits *T4.5*. Output is a single voltage **SPDG**, the trim-demand signal.
- **HMU — "T700 hydromechanical control unit"** [Fig. C9, p.89]. Governs *gas generator
  speed* NG on a droop line, applies the cam limits (idle / acceleration / deceleration),
  and meters fuel. It takes SPDG from the ECU through the torque motor and produces WF.

The ECU **trims** the HMU; the HMU is the primary loop and works alone if the ECU trim is
zero. The ECU's whole authority is exercised through TMRU, the Wf/Ps3 trim from the
torque motor.

### 2.1 ECU signal path [Fig. C1, p.85]

```
PCNP ──►[ 1/(CTPL·s+1) ]──►(+)
                              ├──► SPDER ──►[ GOVERNOR RATE ]──► SPDS1 ──►[ GOVERNOR  ]──► SPDSF ──┐
XQLO ────────────────────►(−)│              [ COMPENSATION  ]              [ DYNAMICS ]           │
PCPRF ───────────────────►(−)┘                    ▲                                               ▼
                                                  │ TRQL                          ┌─────────────────────┐
TORQ45 ──►[ LOAD SHARE  ]──► XQLO ────────────────┘                               │  MAXIMUM  ERROR     │
          [ SPEED TRIM  ]──► TRQL                                                 │     SELECTOR        │
                                                                                  └──────────┬──────────┘
T45 ──►[ THERMOCOUPLE ]──► T45EL ──►(+)                                                      │ SPDSS
       [   HARNESS    ]              ├──► ET45 ──►[ TEMPERATURE  ]──► TSIG ───────────────►(▲)│
       [  DYNAMICS    ]  T45REF ──►(−)             [ COMPENSATION ]                           ▼
                                                                              [ PROPORTIONAL PLUS INTEGRAL ]──► SPDG
                                                                              [        COMPENSATION        ]
```

Signs read at 600 dpi from the summing junctions:

- `SPDER = PCNP/(CTPL·s+1) [+] − XQLO [−] − PCPRF [−]`  (%NP)
- `ET45 = T45EL [+] − T45REF [−]`  (°R)

So SPDER is **positive when sensed NP exceeds the cockpit reference**, i.e. the ECU is a
*speed-topping/droop-trim* path, not a classical error = reference − measurement.

The **MAXIMUM ERROR SELECTOR** takes SPDSF (speed path) and TSIG (T4.5 path) and passes
the larger to the P+I compensation. That is the ECU's limit-arbitration point: T4.5
limiting takes over from NP governing when the temperature error signal is the larger.

### 2.2 ECU sub-blocks

- **Load share speed trim** [Fig. C2, p.85] — TORQ45 (engine 1) through two lags
  1/(TL1·s+1), 1/(TL2·s+1) gives TRQL. TRQL(engine 1) [−] is differenced against
  TRQL(engine 2) [+] to form TRQER. A switch labelled **"ONE ENGINE IMPLEMENTATION"** is
  drawn *open* in this path — with it open the load-share input is zero. TRQER then goes
  through gain-lag 1.1/(CT7·s+1), a saturation with upper limit **CE** and lower limit
  **DBIAS**, then DBIAS is subtracted, then ZK8/(T17·s+1) → **XQLO** (engine 1).
  Self-consistency check: with the switch open, clamp(0, DBIAS, CE) − DBIAS = 0, so
  XQLO = 0 and the load-share path is inert in the single-engine model.
- **Governor rate compensation** [Fig. C3, p.86] — see §4 for the full block list. Forms
  SPDS1 from SPDER, TRQL and the nonlinear NP-loop-gain switch B4.
- **Governor dynamics** [Fig. C4, p.86] — `SPDSF = SPDS1 · (T11·s+1)/(CT2·s+1) · 1/(CT12·s+1)`.
- **Thermocouple harness** [Fig. C5, p.87] — T45 + T45COR → lag 1/(TLGE·s+1) → T45L; the
  power-turbine flow parameter W45R = W45 · √T45L / P45 is formed; F_EC1(T45L, W45R)
  gives the *variable* time constant TAU45; T45EL = T45L/(TAU45·s+1).
- **T45 compensation** [Fig. C6, p.87] — `TSIG = ET45 · ZK3/(CT9·s+1) · (T8·s+1)/(T10·s+1)`.
- **P+I compensation** [Fig. C7, p.88] — parallel XKPROP path and XKINTG→1/s→limiter
  [ZLOLIM, ZHILIM] path, summed to SPDSP, then 1/(CT16·s+1) → SPDG.
- **P+I integrator limit logic** [Fig. C8, p.88] — ZLOLIM is *not* a constant: it is a
  two-level function of engine-2 sensed torque, −1.0 below a threshold of 180 and −0.3
  above it. Table C.1's ZLOLIM = −1.0 is therefore the single-engine / low-torque value.

### 2.3 HMU signal path [Fig. C9, p.89]

```
T2 ───►[ TOPPING LINE SCHEDULE ]──► WFPTP ─────────────────────────────►(+)
                                                                          │
NG ───►[ NG SPEED SENSOR ]──► PCNGHL ──►(−)                               │
                              NGREF ──►(+)  ├──►[ KNDRP ]──────────────►(+)├──► WFPDM ──►┌──────────┐
                                                                          │             │ SELECTOR │──► HMUSEL
PAS ──►[ POWER AVAILABLE SPINDLE ]──► WFPRF ───────────────────────────►(−)│             │  LOGIC   │
       [    INPUT SCHEDULE       ]                                        │             └──▲──▲──▲──┘
                                                                          │                │  │  │
XLDSA►[ LOAD DEMAND SPINDLE ]►[ LOAD DEMAND ]► DWFPL ►(+)                 │          WFPDC │  │  │ WFIDM
      [   INPUT SCHEDULE    ] [  DYNAMICS   ]           ├►[ max(·,0) ]►(−)┘        (decel) │  │  │ (idle)
SPDG ►[ TORQUE MOTOR DYNAMICS ]────────────► TMRU ────►(+)                          WFPAC ─┘  │  │
      [  AND COMPENSATION     ]  (FROM ECU)                                         (accel)   │  │

HMUSEL ──►(×)──► WFMV ──►[ FUEL METERING VALVE ]──► WF        PS3 ──►[ PS3 SENSOR ]──► PS3L ──►(×)
```

Signs read at 600 dpi from the WFPDM summing junction:

- `WFPDM = WFPTP [+] + KNDRP·(NGREF − PCNGHL) [+] − WFPRF [−] − max(DWFPL + TMRU, 0) [−]`

This is the HMU **droop line**: WFPTP (topping, a function of T2) is the ceiling, and
WFPRF (power-available spindle) and the load-demand/torque-motor term are *decrements*
from it. Fig. C25 confirms the sign convention — WFPRF falls monotonically from 10.4 at
PAS = 28° to −0.1 at PAS = 120°, so **more PAS ⇒ less subtraction ⇒ more fuel**. Likewise
Fig. C17 forms DWFP so that more load demand ⇒ smaller (even negative) DWFP ⇒ more fuel.

The unlabelled saturation block between the DWFPL/TMRU summer and the WFPDM summer has a
flat segment on the negative side and a unity-slope segment on the positive side with no
limit labels — i.e. **max(·, 0)**, a lower limit at zero and no upper limit.

**Limit arbitration** [Fig. C18, p.92] is a fixed cascade, in this order:

`HMUSEL = MAX( MIN( MAX(WFPDM, WFIDM), WFPAC ), WFPDC )`

i.e. idle floor first (select max), then the acceleration ceiling (select min), then the
deceleration floor (select max). HMUSEL is multiplied by PS3L to give WFMV in lbm/hr.

---

## 3. Constants and schedules — full transcription

### 3.1 Table C.1 — T700 Fuel Control System Constants

Transcribed cell by cell from the page images re-rendered at 300 dpi
[TM-100991 pdf pp.83–84, Table C.1 and Table C.1 Concluded]. Printed precision preserved
exactly (`1.0` is not `1`, `0.010` is not `0.01`).

**Part 1 — pdf p.83 (28 rows)**

| Constant | Value | Units |
|---|---|---|
| AWFP | 0.05909 | lb_m · in² / lb_f · hr · % |
| B6 | 1.0 | nondimensional |
| BWFP | -3.927 | lb_m · in² / lb_f · hr |
| CB | 0.75 | % NP |
| CE | 13.21 | % NP |
| CH | 0.05 | % NG |
| CLMV | 0.03 | sec |
| CLLDS | 0.2 | sec |
| CNTL | 0.025 | sec |
| CORR | 20.0 | nondimensional |
| CR | 20.0 | nondimensional |
| CT2 | 0.088 | sec |
| CT7 | 1.0 | sec |
| CT9 | 0.010 | sec |
| CT12 | 0.088 | sec |
| CT13 | 1.0 | sec |
| CT14 | 0.10 | sec |
| CT16 | 1.0 | sec |
| CTPL | 0.010 | sec |
| CTPS3 | 0.010 | sec |
| DBIAS | 3.25 | % NP |
| KNDRP | 0.25 | lb_m · in² / lb_f · hr · % |
| NGREF | 101.0 | % NG |
| PS3HYS | 0.375 | lb_f / in² |
| T8 | 0.40 | sec |
| T10 | 0.060 | sec |
| T11 | 0.87 | sec |
| T17 | 0.010 | sec |

**Part 2 — pdf p.84 (29 rows)**

| Constant | Value | Units |
|---|---|---|
| T45COR | 11.0 | °R |
| T45REF | 2004.0 | °R |
| TL1 | 0.010 | sec |
| TL2 | 0.010 | sec |
| TLGE | 0.077 | sec |
| TMDB | 2.0 | ma |
| TMGN | 0.0159 | in / ma · sec |
| TMLG | 84.0 | lb_m · in² / lb_f · hr · in |
| TMLVG | 55.0 | volts / in |
| WFMAX | 785.0 | lb_m / hr |
| WFMIN | 65.0 | lb_m / hr |
| WFPDCH | 2.10 | lb_m · in² / lb_f · hr |
| WFPDCL | 1.45 | lb_m · in² / lb_f · hr |
| XHILIM | 0.0787 | in / sec |
| XKINTG | 0.18 | nondimensional |
| XKPROP | 0.20 | nondimensional |
| XLDHYS | 2.5 | deg |
| XLOLIM | -0.0427 | in / sec |
| YHILIM | 40.0 | nondimensional |
| YLOLIM | 0.0 | nondimensional |
| ZHILIM | 3.5 | nondimensional |
| ZK1 | 1.7 | nondimensional |
| ZK3 | 0.045 | nondimensional |
| ZK5 | 0.40 | nondimensional |
| ZK7 | 0.30 | nondimensional |
| ZK8 | 0.231 | nondimensional |
| ZK9 | 0.625 | nondimensional |
| ZK10 | 0.375 | nondimensional |
| ZLOLIM | -1.0 | nondimensional |

**57 rows total.** Confidence: high — all read from 300 dpi crops with clear glyphs; no
cell was ambiguous. Two OCR traps worth recording, both caught by reading the image:

- text layer gives `CORR 2010 ....` → **image reads 20.0**;
- text layer gives `T10 0.050` → **image reads 0.060**.

Note also that the text layer's `T11` renders as `Tll`, and `PS3 HYS` / `T45 COFI` /
`X KINTG` / `Z K 1` are letter-spacing artefacts, not different symbols.

### 3.2 Constants that appear only inside figures (not in Table C.1)

These are printed numerically on the block diagrams and are just as load-bearing.

| Value | Where | What it is |
|---|---|---|
| 1.1 | Fig. C2, p.85 | numerator gain of the load-share lag `1.1/(CT7·s+1)` |
| −1 | Fig. C3, p.86 | lower input threshold of the nonlinear NP-loop-gain relay (%NP) |
| 4 | Fig. C3, p.86 | upper input threshold of the same relay (%NP) |
| 1 | Fig. C3, p.86 | output level of that relay (and of the CORR relay) |
| 1000 | Fig. C3, p.86 | gain between that relay and the B6 gain |
| −0.44 | Fig. C10, p.89 | bias summed with SPDG at the torque motor input (volts) |
| 0.04 | Fig. C10, p.89 | lead time constant of the torque motor compensation, sec |
| 0.2 | Fig. C10, p.89 | lag time constant of the torque motor compensation, sec |
| 564.0 | Fig. C10, p.89 | torque motor forward gain |
| −31.0 | Fig. C10, p.89 | bias summed after the 564.0 gain (ma) |
| 0.914 | Fig. C13, p.90 | collective-pitch-to-spindle rigging gain, deg/% (UH-60A) |
| 5.34 | Fig. C13, p.90 | collective-pitch-to-spindle rigging bias, deg (UH-60A) |
| 100 / NGDES | Fig. C15, p.91 | NG → PCNG conversion; **NGDES is not defined in Appendix C** |
| 3.26 | Fig. C17, p.91 | bias from which WFQPS3 is subtracted in the load-demand path |
| 1/3600 | Fig. C19, p.92 | lbm/hr → lbm/sec conversion in the metering valve |
| 0.015 | Fig. C19, p.92 | fuel transport delay, sec (`e^(−0.015 s)`) |
| 180 | Fig. C8, p.88 | engine-2 torque threshold in the P+I integrator limit logic |
| −0.3 | Fig. C8, p.88 | ZLOLIM above that threshold |
| −1.0 | Fig. C8, p.88 | ZLOLIM below that threshold (matches Table C.1) |

### 3.3 Schedule breakpoint tables

**There are none.** The eight functions F_EC1 and F_HM1…F_HM7 are given only as plots
(Figs. C23–C30). See §5.

---

## 4. Relation catalogue (in place of an equation catalogue)

Appendix C prints **no numbered equations**. What follows catalogues each block-diagram
relation: figure, page, what it computes, the symbols on each side, and whether it carries
state. "State" means an integrator, a lag/lead-lag, a hysteresis memory, or a transport
delay.

### ECU

| # | Fig. / p. | Computes | Relation (symbols as printed) | State? |
|---|---|---|---|---|
| C-1 | C1 / 85 | sensed NP | PCNP → 1/(CTPL·s+1) | **lag** |
| C-2 | C1 / 85 | speed error | SPDER = PCNP_lagged − XQLO − PCPRF | no |
| C-3 | C1 / 85 | T4.5 error | ET45 = T45EL − T45REF | no |
| C-4 | C1 / 85 | error arbitration | SPDSS = MAXIMUM ERROR SELECTOR(SPDSF, TSIG) | no |
| C-5 | C2 / 85 | sensed PT torque | TRQL = TORQ45 · 1/(TL1·s+1) · 1/(TL2·s+1) | **2 lags** |
| C-6 | C2 / 85 | torque error | TRQER = TRQL(eng 2) − TRQL(eng 1); zero when the "ONE ENGINE IMPLEMENTATION" switch is open | no |
| C-7 | C2 / 85 | load-share trim | XQLO = [clamp(TRQER·1.1/(CT7·s+1), DBIAS, CE) − DBIAS] · ZK8/(T17·s+1) | **2 lags** |
| C-8 | C3 / 86 | NP deadband | deadband on SPDER with breakpoints ±CB | no |
| C-9 | C3 / 86 | rate-comp lags | deadband output → 1/(CT13·s+1) → 1/(CT14·s+1) | **2 lags** |
| C-10 | C3 / 86 | pseudo-derivative | (x_CT13 − x_CT14) · (ZK5/CT14) — constant gain, ZK5 ÷ CT14 = 4.0 | no |
| C-11 | C3 / 86 | compensated error | U = ZK9·SPDER + ZK10·deadband(SPDER) + (ZK5/CT14)·(x_CT13 − x_CT14) | no |
| C-12 | C3 / 86 | speed-error relay | g(SPDER) = 1 for SPDER < −1 or SPDER > 4, else 0 | no |
| C-13 | C3 / 86 | torque integrator | Y = ∫(TRQL − CR) dt, limited to [YLOLIM, YHILIM] | **integrator (limited)** |
| C-14 | C3 / 86 | nonlinear gain switch | B4 = 1 if (1000·B6·g(SPDER) + Y) > CORR, else 0 | no |
| C-15 | C3 / 86 | rate-comp output | SPDS1 = ZK7·U + B4·(ZK1·U) | no |
| C-16 | C4 / 86 | governor dynamics | SPDSF = SPDS1 · (T11·s+1)/(CT2·s+1) · 1/(CT12·s+1) | **lead-lag + lag** |
| C-17 | C5 / 87 | harness bias | T45E = T45 + T45COR | no |
| C-18 | C5 / 87 | harness lag | T45L = T45E / (TLGE·s+1) | **lag** |
| C-19 | C5 / 87 | PT flow parameter | W45R = W45 · √(T45L) · (1/P45) | no |
| C-20 | C5 / 87 | variable time constant | TAU45 = F_EC1(T45L, W45R) | no (table lookup) |
| C-21 | C5 / 87 | sensed T4.5 | T45EL = T45L / (TAU45·s+1) — **time constant varies with the state** | **lag (variable τ)** |
| C-22 | C6 / 87 | T45 compensation | TSIG = ET45 · ZK3/(CT9·s+1) · (T8·s+1)/(T10·s+1) | **lag + lead-lag** |
| C-23 | C7 / 88 | P+I | SPDSP = XKPROP·SPDSS + clamp(∫XKINTG·SPDSS dt, ZLOLIM, ZHILIM) | **integrator (limited)** |
| C-24 | C7 / 88 | P+I output lag | SPDG = SPDSP / (CT16·s+1) | **lag** |
| C-25 | C8 / 88 | integrator lower limit | ZLOLIM = −1.0 if TRQL(eng 2) < 180, else −0.3 | no |

### HMU

| # | Fig. / p. | Computes | Relation (symbols as printed) | State? |
|---|---|---|---|---|
| C-26 | C9 / 89 | fuel-flow demand | WFPDM = WFPTP + KNDRP·(NGREF − PCNGHL) − WFPRF − max(DWFPL + TMRU, 0) | no |
| C-27 | C9 / 89 | metered demand | WFMV = HMUSEL · PS3L | no |
| C-28 | C10 / 89 | torque motor fwd path | e = (SPDG − 0.44) − TMLVG·x_lim; then e·(0.04·s+1)/(0.2·s+1)·564.0 − 31.0 | **lead-lag** |
| C-29 | C10 / 89 | torque motor deadband | deadband ±TMDB on that current signal | no |
| C-30 | C10 / 89 | torque motor integrator | x = ∫TMGN·(deadband output) dt, limited to [XLOLIM, XHILIM] | **integrator (limited)** |
| C-31 | C10 / 89 | LVDT feedback / output | feedback = TMLVG·x (volts); TMRU = TMLG·x | no |
| C-32 | C11 / 90 | PS3 sensor | PS3L = hysteresis(PS3/(CTPS3·s+1), width PS3HYS) | **lag + hysteresis** |
| C-33 | C12 / 90 | load demand dynamics | DWFPL = DWFP / (CLLDS·s+1) | **lag** |
| C-34 | C13 / 90 | collective rigging | XLDSA = 0.914·XCPC + 5.34; XLDSH = hysteresis(XLDSA, width XLDHYS) | **hysteresis** |
| C-35 | C14 / 91 | topping schedule | WFPTP = F_HM1(T2) | no (table lookup) |
| C-36 | C15 / 91 | NG sensor | PCNG = NG·100/NGDES; PCNGHL = hysteresis(PCNG, CH)/(CNTL·s+1) | **hysteresis + lag** |
| C-37 | C16 / 91 | power-available schedule | WFPRF = F_HM2(PAS) | no (table lookup) |
| C-38 | C17 / 91 | load demand schedule | PNG = F_HM3(XLDSH); WFQPS3 = F_HM4(XLDSH); DWFP = KNDRP·(NGREF − PNG) + (3.26 − WFQPS3) | no |
| C-39 | C18 / 92 | limit selector | HMUSEL = MAX( MIN( MAX(WFPDM, WFIDM), WFPAC ), WFPDC ) | no |
| C-40 | C19 / 92 | metering valve | WF = clamp(WFMV, WFMIN, WFMAX) · 1/(CLMV·s+1) · (1/3600) · e^(−0.015 s) | **lag + transport delay** |
| C-41 | C20 / 93 | idle schedule | WFIRF = F_HM5(T2); PCNGI = F_HM6(T2); WFIDM = KNDRP·(PCNGI − PCNGHL) − WFIRF | no |
| C-42 | C21 / 93 | acceleration schedule | WFPAC = F_HM7(T2, PCNGHL) | no (2-D table lookup) |
| C-43 | C22 / 93 | deceleration schedule | WFPDC = clamp(AWFP·PCNGHL + BWFP, WFPDCL, WFPDCH) | no |

Sanity check on C-43 with the Table C.1 values: at PCNGHL = 100 %NG,
0.05909·100 − 3.927 = 1.982, inside [1.45, 2.10] — consistent scale for a Wf/Ps3 floor.

---

## 5. Figure list and digitizing assessment

### 5.1 Block diagrams (describe, do not digitize) — 22 figures

C1, C2 (p.85); C3, C4 (p.86); C5, C6 (p.87); C7, C8 (p.88); C9, C10 (p.89);
C11, C12, C13 (p.90); C14, C15, C16, C17 (p.91); C18, C19 (p.92); C20, C21, C22 (p.93).

All are line art and fully legible at 200 dpi. The only elements that needed 600 dpi were
the ± signs on summing junctions (Figs. C1, C9, C17, C20) and the threshold braces on
nonlinearity icons (Figs. C2, C3, C8). Those are all resolved above; none remain
ambiguous.

Icon conventions used throughout, worked out by cross-checking against the nomenclature:

- **Saturation icon** (flat–ramp–flat): the top brace is the *upper output limit*, the
  bottom brace the *lower output limit*. Examples: ZHILIM/ZLOLIM (C7), XHILIM/XLOLIM
  (C10), WFMAX/WFMIN (C19), WFPDCH/WFPDCL (C22), CE/DBIAS (C2).
- **Deadband icon** (ramp–flat–ramp through the origin): braces mark *input* thresholds.
  Examples: ±CB (C3), ±TMDB (C10).
- **Relay / step icon** (flat levels with dashed verticals): braces mark *input*
  thresholds and the number inside the box is the *output level*. Examples: −1 / 4 with
  output 1 (C3), CORR with output 1 → B4 (C3), 180 with outputs −1.0 / −0.3 (C8).
- **Parallelogram icon**: hysteresis, brace = width. Examples: PS3HYS (C11), XLDHYS (C13),
  CH (C15).
- **Boxed ⊠**: multiplier. Examples: W45R formation (C5), B4·ZK1·U (C3), HMUSEL·PS3L (C9).

### 5.2 Function plots (must be digitized) — 8 figures

All eight are **linear-linear, no gridlines**, with a boxed frame and minor tick marks on
all four sides. Crucially, **every plot marks its own breakpoints**: single-parameter
curves use `×` markers at the knots, multi-curve plots use the digit of the SYMBOL key as
the marker. Digitizing these is therefore breakpoint extraction, not curve tracing.

| Fig. | p. | Function | x axis (range, units) | y axis (range, units) | Curves | Markers | Legibility at 200 dpi | Difficulty |
|---|---|---|---|---|---|---|---|---|
| C23 | 94 | F_EC1 — ECU thermocouple sensor time constant | POWER TURBINE FLOW PARAMETER, W45R — 0 to 16, nondimensional | SENSOR TIME CONSTANT, TAU45, sec — 1 to 6 | 4, parameter T45L = 1260.0 / 1660.0 / 2060.0 / 2460.0 deg R | digits 1–4 at each knot; ~5 knots per curve (W45R ≈ 1, 2.5, 6.5, 10, 15) | good; curves 2/3/4 converge near W45R = 15 | **easy–medium** (converging tail) |
| C24 | 95 | F_HM1 — HMU topping line schedule | INLET TEMPERATURE, T2, deg R — 390 to 640 | FUEL FLOW TOPPING LINE PARAMETER, WFPTP — −0.50 to 4.00 | 1 | `×`, 9 knots | excellent | **easy** |
| C25 | 96 | F_HM2 — HMU power-available input schedule | POWER AVAILABLE SPINDLE ANGLE, deg — 20 to 120 | FUEL FLOW POWER-AVAILABLE PARAMETER, WFPRF — −2 to 12 | 1, monotonically decreasing | `×`, **11** knots (this row said 13 until 2026-09-12; `data/schedules/fhm2_power_available.csv` carries 11) | excellent | **easy** |
| C26 | 97 | F_HM3 — HMU load-demand-compensation schedule function 1 | LOAD DEMAND SPINDLE ANGLE, XLDSA, deg — 0 to 100 | GAS GENERATOR SPEED DEMAND PARAMETER, PNG — 75.0 to 112.5 | 1, rising then flat above 80° | `×`, 11 knots at XLDSA = 0,10,…,100 | excellent | **easy** |
| C27 | 98 | F_HM4 — HMU load-demand-compensation schedule function 2 | LOAD DEMAND SPINDLE ANGLE, XLDSA, deg — 0 to 100 | FUEL FLOW DELTA DEMAND PARAMETER, WFQPS3 — 2.1 to 3.7 | 1, rising then flat above 80° | `×`, 11 knots at XLDSA = 0,10,…,100 | excellent | **easy** |
| C28 | 99 | F_HM5 — HMU idle-schedule function 1 | INLET TEMPERATURE, T2, deg R — 390 to 640 | FUEL FLOW IDLE SCHEDULE LIMIT PARAMETER, WFIRF — 2.05 to 2.65 | 1, monotone rising | `×`, 13 knots | excellent | **easy** |
| C29 | 100 | F_HM6 — HMU idle-schedule function 2 | INLET TEMPERATURE, T2, deg R — 390 to 640 | NG DROOP LINE IDLE REFERENCE PARAMETER, PCNGI — 62 to ~75 (top tick unlabelled) | 1, a straight line | `×`, only **2** knots, at the two ends | excellent | **trivial** (two endpoints define it) |
| C30 | 101 | F_HM7 — HMU maximum fuel-parameter limit during acceleration | SENSED GAS TURBINE SPEED, PCNGHL, % — 50 to 110 | FUEL FLOW ACCELERATION LIMIT PARAMETER, WFPAC — 0 to 5 | 7, parameter T2 = 395.0 / 430.0 / 460.0 / 480.0 / 519.0 / 555.0 / 590.0 deg R | digits 1–7 at each knot | **the hard one** — all seven curves cross in a tight bundle around PCNGHL 84–95, WFPAC 3.3–3.8, where the digit markers overlap and overprint | **hard** — expect to need 600 dpi and careful per-curve tracing through the crossing region |

No data values were read off any curve; the ranges above are axis labels only.

---

## 6. Control system state variables

These add to the engine's 5 states (NG, NP, P3, P41, P45). **22 continuous states**, plus
3 hysteresis memories and 1 transport-delay buffer.

### Continuous states — first-order lags and lead-lags (19)

| # | State of | Time constant | Value | Fig. / p. |
|---|---|---|---|---|
| 1 | NP speed sensor | CTPL | 0.010 s | C1 / 85 |
| 2 | Load-share torque sensor, stage 1 | TL1 | 0.010 s | C2 / 85 |
| 3 | Load-share torque sensor, stage 2 | TL2 | 0.010 s | C2 / 85 |
| 4 | Load-share circuit gain-lag (numerator 1.1) | CT7 | 1.0 s | C2 / 85 |
| 5 | Load-share path output | T17 | 0.010 s | C2 / 85 |
| 6 | Governor rate compensation, first lag | CT13 | 1.0 s | C3 / 86 |
| 7 | Governor rate compensation, second lag | CT14 | 0.10 s | C3 / 86 |
| 8 | Governor dynamics lead-lag (T11·s+1)/(CT2·s+1) | CT2 (lead T11) | 0.088 s (0.87 s) | C4 / 86 |
| 9 | Governor dynamics lag | CT12 | 0.088 s | C4 / 86 |
| 10 | Thermocouple harness lag | TLGE | 0.077 s | C5 / 87 |
| 11 | Thermocouple sensor lag — **τ is state-dependent**, TAU45 = F_EC1(T45L, W45R) | TAU45 | ≈ 1.7 to 5.0 s from Fig. C23 | C5 / 87 |
| 12 | T4.5 compensation lag (numerator ZK3) | CT9 | 0.010 s | C6 / 87 |
| 13 | T4.5 compensation lead-lag (T8·s+1)/(T10·s+1) | T10 (lead T8) | 0.060 s (0.40 s) | C6 / 87 |
| 14 | P+I output lag | CT16 | 1.0 s | C7 / 88 |
| 15 | Torque motor compensation lead-lag (0.04·s+1)/(0.2·s+1) | 0.2 s (lead 0.04 s) | — | C10 / 89 |
| 16 | PS3 sensor lag | CTPS3 | 0.010 s | C11 / 90 |
| 17 | Load demand (fuel flow command) dynamics | CLLDS | 0.2 s | C12 / 90 |
| 18 | NG spool sensor lag | CNTL | 0.025 s | C15 / 91 |
| 19 | Fuel metering valve lag | CLMV | 0.03 s | C19 / 92 |

### Continuous states — limited integrators (3)

| # | Integrator | Input | Output limits | Fig. / p. |
|---|---|---|---|---|
| 20 | Engine torque integrator, `1/s`, for the nonlinear NP-loop-gain circuit | TRQL − CR | [YLOLIM, YHILIM] = [0.0, 40.0] | C3 / 86 |
| 21 | ECU P+I integrator, `1/s` | XKINTG · SPDSS | [ZLOLIM, ZHILIM] = [−1.0 or −0.3, 3.5] — lower limit is itself switched by Fig. C8 | C7 / 88 |
| 22 | Torque motor integrator, `TMGN/s` | deadbanded torque-motor current | [XLOLIM, XHILIM] = [−0.0427, 0.0787] in/sec — **note the units: these are rate limits on a position state** | C10 / 89 |

### Discrete / memory elements (4)

| # | Element | Parameter | Value | Fig. / p. |
|---|---|---|---|---|
| 23 | NG speed sensor hysteresis | CH | 0.05 % NG | C15 / 91 |
| 24 | PS3 sensor hysteresis | PS3HYS | 0.375 lb_f/in² | C11 / 90 |
| 25 | Load demand spindle hysteresis | XLDHYS | 2.5 deg | C13 / 90 |
| 26 | Fuel transport delay | `e^(−0.015 s)` | 0.015 s | C19 / 92 |

### Memoryless nonlinearities (no state, listed for completeness)

Deadbands ±CB (C3) and ±TMDB (C10); relays with thresholds −1/4 and CORR producing B4
(C3); the ZLOLIM switch at 180 (C8); the saturations CE/DBIAS (C2), max(·,0) (C9),
WFMIN/WFMAX (C19), WFPDCL/WFPDCH (C22); the three-stage MAX/MIN/MAX selector (C18).

---

## 7. Ambiguities and things marked

Nothing was left as a digit-level ambiguity — all 57 Table C.1 cells and all figure-embedded
numbers resolved cleanly at 300–600 dpi. The open items are *interpretive*, and are also
appended to `open-questions.md`:

1. **WFIDM sign (Fig. C20, p.93) — RESOLVED, see open-questions.md #10.** The summing
   junction is unambiguous at 600 dpi: KNDRP·(PCNGI − PCNGHL) enters `+`, WFIRF enters `−`,
   giving `WFIDM = KNDRP·(PCNGI − PCNGHL) − WFIRF`. **The signs are right and an earlier
   draft of this file drew the wrong conclusion from them.**

   The `≈ −2.4` is `−WFIRF` evaluated *at the idle reference*, where `PCNGI − PCNGHL ≈ 0`.
   That is exactly what a droop-line floor should do — it must not override the governor at
   or above idle. With KNDRP = 0.25 and WFIRF ≈ 2.05–2.65, WFIDM turns positive once sensed
   NG droops `WFIRF/KNDRP` ≈ 8.2–10.6 percentage points below the idle reference, and
   outbids a 1.5–4 demand at ≈ 14–27 points below. That is deep sub-idle — and the report
   [pdf p.38] explicitly eliminates *"fuel control below flight-idle power"*. The floor is
   inactive in every modelled condition **because the model stops above its regime**, not
   because of a sign error. Implement as printed.
2. **CR / TRQL units (Fig. C3, p.86) — RESOLVED, see open-questions.md #11.** TRQL is
   ft·lbf and is *not* normalised anywhere; the path TORQ45 → 1/(TL1·s+1) → 1/(TL2·s+1) →
   TRQL is unity DC gain, with no gain, divide or lookup before Figs. C3 and C8. The scale
   is small because this is *shaft*-referenced torque: Table B.1 (p.67) prints "Engine
   torque (ref. shaft)" = 229.0 / 138.9 / 76.06 ft·lbf at the three trims (= 5252·hp/NP at
   NP = 20895 rpm), against "Required torque (ref. hub)" = 32865. / 20747. / 10792., which
   is the rotor-hub number and never reaches the ECU. So CR = 20.0 ft·lbf ≡ 79.6 hp (a
   sub-idle "engine is loaded" threshold — p.77 calls CR "engine torque-level threshold for
   nonlinear NP loop gain circuit") and the Fig. C8 threshold 180 ft·lbf ≡ 716 hp (between
   the 80-kt and hover trims). B4 is **not** permanently 1: Y = ∫(TRQL − CR)dt clamped to
   [0.0, 40.0] crosses CORR = 20.0 in 0.10–0.36 s at flight power and unwinds to 0 in ≥ 2 s
   once TRQL < 20 ft·lbf, which is what ZK1 ("used during high power operation", p.81)
   requires. What remains is only a labelling defect: Table C.1's "nondimensional" for CR,
   CORR, YHILIM and YLOLIM is the report's own word (verified at 600 dpi) but is wrong —
   CR is ft·lbf and YLOLIM/YHILIM are ft·lbf·sec. Use the printed *numbers*; ignore that
   units column for these four.
3. **B4 has no value in Table C.1** — correctly so, since it is computed (Fig. C3), but
   worth stating: the nomenclature calls it a "discrete switch" and the figure makes it
   {0, 1}.
4. **NGDES is used but never defined in Appendix C** (Fig. C15, the `100/NGDES` gain).
   Presumably an Appendix A constant.
5. **F_HM3 / F_HM4 argument** — Fig. C17 feeds them XLDSH (post-hysteresis), but the
   x-axes of Figs. C26 and C27 are labelled XLDSA (pre-hysteresis). **Implemented as the
   block diagram draws it: the wire carries XLDSH.** The axis label is the more likely
   slip, since a schedule is naturally drawn against the physical spindle angle while the
   wire carries whatever the block upstream produced. The difference is bounded by
   XLDHYS = 2.5 deg.
6. **The unlabelled saturation in Fig. C9** — read as `max(·, 0)` from the icon; it carries
   no limit labels at all, unlike every other saturation in the appendix.
7. **No breakpoint tables for F_EC1 and F_HM1…F_HM7.** The eight schedule functions exist
   only as Figs. C23–C30 and had to be digitized before any of the fuel control could run.
   **All eight now are** (Fig. C30 last, 2026-09-13, open question #50 — closed), and they
   load through `src/t700/control/schedules.py`.
8. **`XCPC` is percent, not inches** [nomenclature, pdf p.81, read off the raster:
   *"helicopter collective pitch position in percent of maximum, percent"*]. Worth stating
   because the Fig. C13 rigging `XLDSA = 0.914·XCPC + 5.34` only reaches the load demand
   spindle's printed 0–100 deg range if XCPC is a percentage; read as inches it spans just
   5.3–14.5 deg and the whole control sits against its deceleration floor.
9. **The torque motor null is not 0.44** [Fig. C10, pdf p.89]. The first summing junction
   forms `SPDG − 0.44`, but the forward path is `e·(0.04s+1)/(0.2s+1)·564.0 − 31.0`, so at
   `e = 0` it delivers −31.0, far outside the ±TMDB = ±2.0 deadband, and the integrator
   ramps to `XLOLIM`. The signal sits still at `0.44 + 31.0/564.0 = 0.49496`. See #55.
