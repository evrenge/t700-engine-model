# Appendix B — Linear Models: full inventory and transcription

Source: Ballin, *A High-Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
NASA TM-100991, July 1988. File `docs/ballin-tm100991.pdf`.

**Scope of this note:** PDF pages 67–76 (printed pages 53–62; printed = pdf − 14).
Appendix B begins on pdf p.67 and ends on pdf p.76 — pdf p.77 opens Appendix C, and
pdf p.66 is the last page of Appendix A.

**Method:** every number below was read from a page image, never from the text layer.
Pages were re-rendered at 300 dpi (`pdftoppm -r 300`) and cropped column-by-column;
each matrix element was read at least twice from independent crops. The text layer of
these pages is badly corrupted (it drops minus signs, reorders vector elements, and on
pdf p.69 it omitted an entire `b` element) and was used only to locate content.

---

## 1. Page-by-page map

| PDF p. | Printed p. | Content |
|---|---|---|
| 67 | 53 | Appendix B title, one-paragraph preamble, **Table B.1** (trim conditions) |
| 68 | 54 | **Fig. B1** 2-DOF, trim 1 · **Fig. B2** 5-DOF, trim 1 |
| 69 | 55 | **Fig. B3** 2-DOF, trim 2 · **Fig. B4** 5-DOF, trim 2 |
| 70 | 56 | **Fig. B5** 2-DOF, trim 3 · **Fig. B6** 5-DOF, trim 3 |
| 71 | 57 | **Fig. B7** 3-DOF, trim 1 |
| 72 | 58 | **Fig. B8** 6-DOF, trim 1 |
| 73 | 59 | **Fig. B9** 3-DOF, trim 2 |
| 74 | 60 | **Fig. B10** 6-DOF, trim 2 |
| 75 | 61 | **Fig. B11** 3-DOF, trim 3 |
| 76 | 62 | **Fig. B12** 6-DOF, trim 3 |

Twelve figures, four DOF variants (2, 3, 5, 6), three trim conditions.

> **Correction to the task brief.** Appendix B contains **2-, 3-, 5- and 6-DOF** models,
> not 2/3/5. The 3-DOF and 6-DOF variants are the heat-sink models. Preamble, pdf p.67,
> verbatim:
>
> > "Small-perturbation linear models for three trim conditions are given below. The two-
> > and three-degree-of-freedom models approximate the dynamics between the control volumes
> > to be instantaneous; the five- and six-degree-of-freedom models contain complete
> > dynamics. The heat-sink model is contained in the three- and six-degree-of-freedom
> > models."

---

## 2. State and control vectors

**Appendix B itself prints no state-vector definition, no control-vector definition, and
no units on any matrix.** The figures show bare `A`, `b̄`, `C`, `d̄` arrays. The
definitions below come from the main text and the nomenclature, cited per item.

### 2.1 The 6-DOF state vector — printed explicitly

[TM-100991 pdf p.32], displayed equation immediately below Eq. 65, read from the page image:

```
      ⎧ NG  ⎫
      ⎪ NP  ⎪
 x̄ =  ⎨ P3  ⎬
      ⎪ P41 ⎪
      ⎪ P45 ⎪
      ⎩ T41 ⎭
```

This is the only place in the report where a linear-model state vector is written out.

### 2.2 The 5-DOF state vector — confirmed

The brief's inference (NG, NP, P3, P41, P45, from pdf p.28) is **correct**, and the report
states it directly in **Table 1 [pdf p.31]**, whose "Mode" column lists, for the 5-DOF
model, the five rows `NG`, `NP`, `P3`, `P41`, `P45` in that order. It is corroborated by
pdf p.32 (above), whose 6-DOF vector is the 5-DOF vector with `T41` appended, and by
pdf p.28: "Order reduction was performed by setting the state derivatives Ṗ3, Ṗ41, and
Ṗ45 to zero and solving the equations for the remaining two states."

Independent numerical confirmation: the eigenvalues of the three 5×5 `A` matrices as
transcribed below reproduce Table 1 (pdf p.31) to three significant figures — see §6.

### 2.3 The 2-DOF state vector — confirmed

`{NG, NP}`. Table 1 [pdf p.31] lists only the `NG` and `NP` rows under the "2-DOF Model"
column. Confirmed numerically: the printed 2-DOF `A` is lower-triangular, so its
eigenvalues are its diagonal, and they match Table 1's 2-DOF column exactly.

### 2.4 The 3-DOF state vector — inferred, flagged

**`{NG, NP, T41}`. This ordering is not printed anywhere in the report.** It is inferred
from three pieces of printed evidence:

1. pdf p.67: the heat-sink model is in the 3- and 6-DOF models; pdf p.29 says the
   heat sink "requires the addition of an extra state. This state, gas temperature at
   station 4.1 …". So 3-DOF = 2-DOF states + T41.
2. `d̄` (Eq. 67 feedthrough) is `F1⁻¹G2`, and `G2` [pdf p.33] is zero in every row except
   the T41 row. In every printed 3-DOF `d̄` (Figs. B7, B9, B11) rows 1 and 2 are
   `0.0000E+0` and only **row 3** is nonzero — so row 3 is the T41 row.
3. The 3-DOF `A(2,2)` equals the 2-DOF and 5-DOF `A(2,2)` exactly at every trim
   (−0.5650E+0 / −0.4461E+0 / −0.3567E+0), pinning row 2 = NP; column 2 is otherwise all
   zeros, as it is in the 5-/6-DOF models.

Confidence: high, but mark it `# UNVERIFIED-ORDERING` at any code site that depends on it.

### 2.5 Control (input) vector

**Single scalar input: `Wf`, fuel flow, `lbm/sec`** [nomenclature, pdf p.13, read from
image: "W_f  fuel flow, lb_m/sec"]. `b̄` is a column vector in every variant, and Eq. 66
[pdf p.33] is `ż̄ = A z̄ + b̄ Wf`. There is no second input column anywhere in Appendix B.

Note the collision with Table B.1, which reports trim fuel flow in **lbm/hr**. The
matrices' input is per-second (see §6 for the consistency check).

### 2.6 Units (derived, not printed)

Nomenclature units [pdf pp.11–13, read from images]:

| Symbol | Meaning | Unit |
|---|---|---|
| NG | rotational speed of compressor and gas generator | rpm |
| NP | rotational speed of power turbine and output shaft | rpm |
| P3 | station 3 **total** pressure | lbf/in² |
| P41 | station 4.1 total pressure | lbf/in² |
| P45 | station 4.5 total pressure | lbf/in² |
| T41 | station 4.1 temperature | deg R |
| Wf | fuel flow | lbm/sec |
| τ1, τ2 | heat-sink lead / lag time constants at trim | sec |

Therefore, by Eqs. 66–67:

- `A(i,j) = ∂ẋi/∂xj` → units `[xi] · s⁻¹ / [xj]`
- `b̄(i)  = ∂ẋi/∂Wf` → units `[xi] · s⁻¹ / (lbm/s)`
- `d̄(i)` → units `[xi] / (lbm/s)` (an output offset, not a rate)

### 2.7 What `C` and `d̄` mean

Only the 3- and 6-DOF (heat-sink) figures carry `C` and `d̄`. From pdf p.31–33:

- Eq. 65: `F1 ẋ̄ = F2 x̄ + G1 Wf + G2 Ẇf` — the heat-sink lead/lag (Eq. 63,
  `T41/T41ns = (τ1 s + 1)/(τ2 s + 1)`) introduces a **Ẇf** term.
- Eq. 66: `ż̄ = A z̄ + b̄ Wf`
- Eq. 67: `x̄ = C z̄ + d̄ Wf`
- with `A = F1⁻¹F2`, `b̄ = F1⁻¹G1 + F F1⁻¹G2`, `C = I`, `d̄ = F1⁻¹G2`
  ("as derived by Chen", pdf p.33).

So `z̄` is a transformed state, and the **physical** state vector is recovered as
`x̄ = z̄ + d̄·Wf`. `C = [I]` is printed literally in Figs. B7–B12 — an identity of order 3
or 6 as appropriate; no numeric identity matrix is spelled out.

> **Flag — printed symbol `F`.** The printed `b̄ = F1⁻¹G1 + F F1⁻¹G2` [pdf p.33] uses a
> symbol `F` that is defined nowhere in the report (`F1`, `F2`, `G1`, `G2` are all
> defined; bare `F` is not, and it is not in the nomenclature). The standard Chen
> transformation for `ẋ = Ax + B0u + B1u̇` gives `b̄ = B0 + A·B1`, which would make `F` ≡
> `A` = `F1⁻¹F2`. **Transcribed as printed; do not assume.** Add to open-questions.

---

## 3. Table B.1 — the three trim conditions

**TABLE B.1 — SMALL-PERTURBATION MODEL TRIM CONDITIONS** [TM-100991 pdf p.67, printed
p.53]. Whole table transcribed from the 300-dpi page image, cell by cell.

| Trim Condition | 1ᵃ | 2ᵇ | 3ᶜ |
|---|---|---|---|
| Aircraft weight, *lbₘ* | 16825. | 16825. | 16825. |
| CG station, *in* | 355.0 | 355.0 | 355.0 |
| CG waterline, *in* | 248.2 | 248.2 | 248.2 |
| CG buttline, *in* | 0.0 | 0.0 | 0.0 |
| Equivalent airspeed, *kts* | 0.0 | 80.0 | 80.0 |
| Flightpath angle, *deg* | 0.0 | 0.0 | **-7.08** |
| Altitude, *ft* | 0.0 | 0.0 | 0.0 |
| Required torque (ref. hub), *ft·lb_f* | 32865. | 20747. | 10792. |
| Engine torque (ref. shaft), *ft·lb_f* | 229.0 | 138.9 | 76.06 |
| Horsepower per engine | 911.1 | 552.6 | 302.6 |
| Power turbine speed, *rpm* | 20895. | 20895. | 20895. |
| Gas Generator Speed, *rpm M* [sic] | 41638. | 39768. | 38072. |
| Fuel flow per engine, *lbₘ/hr* | 476.3 | 349.3 | 267.7 |
| P2, *PSIA* | 14.696 | 14.696 | 14.696 |
| T2, *deg R* | 518.67 | 518.67 | 518.67 |
| P_s3, *PSIA* | 176.34 | 142.13 | 114.27 |
| P41, *PSIA* | 174.28 | 140.26 | 112.77 |
| T41, *deg R* | 2292. | 2102. | 1982. |
| P45, *PSIA* | 37.42 | 30.66 | 25.54 |
| T45, *deg R* | 1632. | 1501. | 1424. |

Footnotes, printed below the table:

- ᵃ Trim condition 1: **hover**
- ᵇ Trim condition 2: **level flight at 80 knots**
- ᶜ Trim condition 3: **1000 ft/min descending flight at 80 knots**

Ambient conditions are given only as `P2 = 14.696 PSIA`, `T2 = 518.67 deg R`,
`Altitude = 0.0 ft` — i.e. sea-level standard day at all three trims. Note pdf p.22: P2
and T2 are taken as ambient (inlet losses assumed offset by stagnation effects), so these
are also the ambient values.

Reading notes on Table B.1:

- **`rpm M`** — the gas-generator-speed row really is printed `rpm M`, with a stray
  italic `M` after the unit; verified at 300 dpi and at 2× zoom. Almost certainly a
  typesetting artifact; the unit is rpm (cf. nomenclature p.11). Recorded as printed.
- **`P_s3`** — the row label carries a small subscript glyph between `P` and `3`. It reads
  as `s`, giving `Ps3` = *station 3 static pressure, a fuel-control-system parameter*
  [nomenclature pdf p.11]. **This is not the same symbol as the state `P3`** (station 3
  *total* pressure). Either the table means the static pressure (in which case Appendix B
  gives no trim value for the state P3) or the subscript is a typo for total. The glyph is
  blotted at 300 dpi; candidates `s` (strongly favoured, matches the nomenclature entry)
  or `a`. **Flagged as ambiguous — see §7.**
- Trim-3 engine torque is `76.06` (the text layer's `76106` is OCR damage).
- Flightpath angle for trim 3 is **negative**: `-7.08 deg` (descent). Sign verified.

---

## 4. The matrices — 2-DOF and 5-DOF (no heat sink)

Model form for these two variants: `ẋ̄ = A x̄ + b̄ Wf`. No `C`, no `d̄` is printed.

Row *i* = the derivative of state *i*; column *j* = the state being perturbed.

### 4.1 Figure B1 — 2 DOF, trim condition 1 [pdf p.68]

States (rows and columns), in order: **NG, NP**

```
        NG            NP
NĠ [ -0.2693E+1    0.0000E+0 ]
NṖ [  0.3865E+0   -0.5650E+0 ]

b̄ = [ 0.1296E+6 ]   NG row
    [ 0.2105E+5 ]   NP row
```

### 4.2 Figure B2 — 5 DOF, trim condition 1 [pdf p.68]

States, in order: **NG, NP, P3, P41, P45**

```
          NG            NP            P3            P41           P45
NĠ  [ -0.4474E+1    0.0000E+0   -0.6928E+3    0.1462E+4   -0.2942E+4 ]
NṖ  [ -0.2444E-1   -0.5650E+0   -0.4128E+2    0.3302E+2    0.2321E+3 ]
Ṗ3  [  0.7524E+0    0.0000E+0   -0.4230E+3    0.4071E+3    0.0000E+0 ]
Ṗ41 [  0.0000E+0    0.0000E+0    0.4115E+4   -0.4526E+4    0.0000E+0 ]
Ṗ45 [  0.1063E+1    0.0000E+0   -0.2695E+3    0.8101E+3   -0.3063E+4 ]

b̄ = [ 0.8238E+5 ]  NG
    [ 0.6174E+4 ]  NP
    [ 0.0000E+0 ]  P3
    [ 0.1891E+6 ]  P41
    [ 0.4021E+5 ]  P45
```

Note `A(2,1) = -0.2444E-1` — **negative mantissa AND negative exponent**. This element is
positive at trims 2 and 3. Verified at 300 dpi and at 2× zoom.

### 4.3 Figure B3 — 2 DOF, trim condition 2 [pdf p.69]

```
        NG            NP
NĠ [ -0.2233E+1    0.0000E+0 ]
NṖ [  0.3128E+0   -0.4461E+0 ]

b̄ = [ 0.1494E+6 ]   NG
    [ 0.1663E+5 ]   NP
```

(The text layer prints only `0.1663E+5` and drops `0.1494E+6` entirely — the image has
both.)

### 4.4 Figure B4 — 5 DOF, trim condition 2 [pdf p.69]

```
          NG            NP            P3            P41           P45
NĠ  [ -0.3659E+1    0.0000E+0   -0.6173E+3    0.1107E+4   -0.1548E+4 ]
NṖ  [  0.1283E-1   -0.4461E+0   -0.2898E+2    0.2551E+2    0.1825E+3 ]
Ṗ3  [  0.6340E+0    0.0000E+0   -0.4090E+3    0.3899E+3    0.0000E+0 ]
Ṗ41 [  0.0000E+0    0.0000E+0    0.3898E+4   -0.4287E+4    0.0000E+0 ]
Ṗ45 [  0.8461E+0    0.0000E+0   -0.2551E+3    0.9497E+3   -0.4038E+4 ]

b̄ = [ 0.8411E+5 ]  NG
    [ 0.4888E+4 ]  NP
    [ 0.0000E+0 ]  P3
    [ 0.1882E+6 ]  P41
    [ 0.4221E+5 ]  P45
```

`A(2,1) = +0.1283E-1` — **positive**, unlike trim 1. Deliberately re-checked at 2× zoom;
there is no minus sign. The repeated mantissa `2551` in `A(2,4) = 0.2551E+2` and
`A(5,3) = -0.2551E+3` is genuine (both re-read at 1.6× zoom).

### 4.5 Figure B5 — 2 DOF, trim condition 3 [pdf p.70]

```
        NG            NP
NĠ [ -0.1816E+1    0.0000E+0 ]
NṖ [  0.3615E+0   -0.3567E+0 ]

b̄ = [ 0.1531E+6 ]   NG
    [ 0.1435E+5 ]   NP
```

### 4.6 Figure B6 — 5 DOF, trim condition 3 [pdf p.70]

```
          NG            NP            P3            P41           P45
NĠ  [ -0.3696E+1    0.0000E+0   -0.5959E+3    0.9809E+3   -0.1034E+4 ]
NṖ  [  0.3555E-1   -0.3567E+0   -0.1940E+2    0.1774E+2    0.1690E+3 ]
Ṗ3  [  0.7624E+0    0.0000E+0   -0.3989E+3    0.3783E+3    0.0000E+0 ]
Ṗ41 [  0.0000E+0    0.0000E+0    0.3803E+4   -0.4179E+4    0.0000E+0 ]
Ṗ45 [  0.6997E+0    0.0000E+0   -0.2619E+3    0.1000E+4   -0.4431E+4 ]

b̄ = [ 0.8525E+5 ]  NG
    [ 0.3436E+4 ]  NP
    [ 0.0000E+0 ]  P3
    [ 0.1876E+6 ]  P41
    [ 0.4496E+5 ]  P45
```

---

## 5. The matrices — 3-DOF and 6-DOF (heat-sink models)

Model form: `ż̄ = A z̄ + b̄ Wf`; `x̄ = C z̄ + d̄ Wf` with `C = [I]` (Eqs. 66–67, pdf p.33).
`C = [ I ]` is printed literally in all six figures; no numeric entries.

### 5.1 Figure B7 — 3 DOF, trim condition 1 [pdf p.71]

States, in order: **NG, NP, T41** (ordering inferred — see §2.4)

```
          NG            NP           T41
NĠ  [ -0.1618E+1    0.0000E+0    0.1472E+2 ]
NṖ  [  0.5569E+0   -0.5650E+0    0.2432E+1 ]
Ṫ41 [  0.5148E-1    0.0000E+0   -0.9891E+0 ]

b̄ = [  0.8500E+5 ]  NG
    [  0.1362E+5 ]  NP
    [ -0.3011E+4 ]  T41      <-- negative

C = [ I ]

d̄ = [  0.0000E+0 ]  NG
    [  0.0000E+0 ]  NP
    [  0.5294E+4 ]  T41
```

### 5.2 Figure B8 — 6 DOF, trim condition 1 [pdf p.72]

States, in order: **NG, NP, P3, P41, P45, T41** [explicit, pdf p.32]

```
          NG            NP            P3            P41           P45           T41
NĠ  [ -0.4474E+1    0.0000E+0   -0.1419E+3    0.9129E+3   -0.2942E+4    0.1049E+2 ]
NṖ  [ -0.2444E-1   -0.5650E+0    0.0000E+0   -0.7748E+1    0.2321E+3    0.7867E+0 ]
Ṗ3  [  0.7524E+0    0.0000E+0   -0.4230E+3    0.4071E+3    0.0000E+0    0.0000E+0 ]
Ṗ41 [  0.0000E+0    0.0000E+0    0.5252E+4   -0.5653E+4    0.0000E+0    0.2229E+2 ]
Ṗ45 [  0.1063E+1    0.0000E+0    0.0000E+0    0.5445E+3   -0.3063E+4    0.5123E+1 ]
Ṫ41 [ -0.2659E+2    0.0000E+0    0.1984E+6   -0.2118E+6    0.0000E+0    0.7785E+3 ]

b̄ = [ 0.5533E+5 ]  NG
    [ 0.4147E+4 ]  NP
    [ 0.0000E+0 ]  P3
    [ 0.1316E+6 ]  P41
    [ 0.2701E+5 ]  P45
    [ 0.4600E+7 ]  T41

C = [ I ]

d̄ = [ 0.0000E+0 ]  NG
    [ 0.0000E+0 ]  NP
    [ 0.0000E+0 ]  P3
    [ 0.0000E+0 ]  P41
    [ 0.0000E+0 ]  P45
    [ 0.5271E+4 ]  T41
```

`A(6,6) = +0.7785E+3` is **positive** and carries no minus sign — read three times from
three independent 300-dpi crops. See §6.3 before drawing any conclusion from it.

### 5.3 Figure B9 — 3 DOF, trim condition 2 [pdf p.73]

```
          NG            NP           T41
NĠ  [ -0.1305E+1    0.0000E+0    0.1439E+2 ]
NṖ  [  0.4215E+0   -0.4461E+0    0.1596E+1 ]
Ṫ41 [  0.3364E-1    0.0000E+0   -0.8622E+0 ]

b̄ = [  0.9239E+5 ]
    [  0.1030E+5 ]
    [ -0.2673E+4 ]           <-- negative

C = [ I ]

d̄ = [  0.0000E+0 ]
    [  0.0000E+0 ]
    [  0.6076E+4 ]
```

### 5.4 Figure B10 — 6 DOF, trim condition 2 [pdf p.74]

```
          NG            NP            P3            P41           P45           T41
NĠ  [ -0.3659E+1    0.0000E+0   -0.1185E+3    0.6127E+3   -0.1548E+4    0.8968E+1 ]
NṖ  [  0.1283E-1   -0.4461E+0    0.0000E+0   -0.3099E+1    0.1825E+3    0.5212E+0 ]
Ṗ3  [  0.6340E+0    0.0000E+0   -0.4090E+3    0.3899E+3    0.0000E+0    0.0000E+0 ]
Ṗ41 [  0.0000E+0    0.0000E+0    0.4909E+4   -0.5287E+4    0.0000E+0    0.1868E+2 ]
Ṗ45 [  0.8461E+0    0.0000E+0   -0.4315E+1    0.7010E+3   -0.4038E+4    0.4501E+1 ]
Ṫ41 [ -0.2279E+2    0.0000E+0    0.1889E+6   -0.2017E+6    0.0000E+0    0.6630E+3 ]

b̄ = [ 0.5425E+5 ]  NG
    [ 0.3153E+4 ]  NP
    [ 0.0000E+0 ]  P3
    [ 0.1260E+6 ]  P41
    [ 0.2722E+5 ]  P45
    [ 0.4473E+7 ]  T41

C = [ I ]

d̄ = [ 0.0000E+0, 0.0000E+0, 0.0000E+0, 0.0000E+0, 0.0000E+0, 0.6048E+4 ]ᵀ
```

### 5.5 Figure B11 — 3 DOF, trim condition 3 [pdf p.75]

```
          NG            NP           T41
NĠ  [ -0.7428E+0    0.0000E+0    0.1217E+2 ]
NṖ  [  0.4617E+0   -0.3567E+0    0.1130E+1 ]
Ṫ41 [  0.1659E-1    0.0000E+0   -0.8819E+0 ]

b̄ = [  0.9056E+5 ]
    [  0.8540E+4 ]
    [ -0.3507E+4 ]           <-- negative

C = [ I ]

d̄ = [  0.0000E+0 ]
    [  0.0000E+0 ]
    [  0.7076E+4 ]
```

Note `A(1,1) = -0.7428E+0` here — exponent `+0`, not `+1` as at trims 1 and 2. The OCR
text layer misassigns this element to row 2; the image is unambiguous.

### 5.6 Figure B12 — 6 DOF, trim condition 3 [pdf p.76]

```
          NG            NP            P3            P41           P45           T41
NĠ  [ -0.3696E+1    0.0000E+0   -0.1143E+3    0.5044E+3   -0.1034E+4    0.7507E+1 ]
NṖ  [  0.3555E-1   -0.3567E+0    0.0000E+0   -0.1431E+1    0.1690E+3    0.3025E+0 ]
Ṗ3  [  0.7624E+0    0.0000E+0   -0.3989E+3    0.3783E+3    0.0000E+0    0.0000E+0 ]
Ṗ41 [  0.0000E+0    0.0000E+0    0.4767E+4   -0.5133E+4    0.0000E+0    0.1544E+2 ]
Ṗ45 [  0.6997E+0    0.0000E+0   -0.7397E+1    0.7474E+3   -0.4431E+4    0.3959E+1 ]
Ṫ41 [ -0.3040E+2    0.0000E+0    0.2036E+6   -0.2172E+6    0.0000E+0    0.6080E+3 ]

b̄ = [ 0.5287E+5 ]  NG
    [ 0.2131E+4 ]  NP
    [ 0.0000E+0 ]  P3
    [ 0.1210E+6 ]  P41
    [ 0.2789E+5 ]  P45
    [ 0.4767E+7 ]  T41

C = [ I ]

d̄ = [ 0.0000E+0, 0.0000E+0, 0.0000E+0, 0.0000E+0, 0.0000E+0, 0.7044E+4 ]ᵀ
```

---

## 6. Eigenvalues, time constants, transfer functions, mode descriptions

### 6.1 Appendix B prints none

There are **no** eigenvalues, time constants, transfer functions or mode descriptions
anywhere in pdf pp.67–76. The appendix is Table B.1 plus twelve bare matrix figures.
Everything in this section comes from the main text and is cited as such.

### 6.2 Table 1 [pdf p.31] — eigenvalues (outside Appendix B, but the matching data)

**Table 1: Eigenvalues of two extracted perturbation models and the reduced-order model.**
Transcribed from the 200-dpi page image; blank cells are blank in the original.

| Trim | Mode | 5-DOF Model | 2-DOF Model | Reduced-Order 5-DOF Model |
|---|---|---|---|---|
| 1 | NG  | -2.66   | -2.69  | -2.81  |
| 1 | NP  | -0.565  | -0.565 | -0.565 |
| 1 | P3  | -51.6   | | |
| 1 | P41 | -4900.  | | |
| 1 | P45 | -3060.  | | |
| 2 | NG  | -2.08   | -2.23  | -2.16  |
| 2 | NP  | -0.446  | -0.446 | -0.446 |
| 2 | P3  | -52.2   | | |
| 2 | P41 | -4640.  | | |
| 2 | P45 | -4040.  | | |
| 3 | NG  | -1.75   | -1.82  | -1.83  |
| 3 | NP  | -0.357  | -0.357 | -0.357 |
| 3 | P3  | -52.6   | | |
| 3 | P41 | -4430.  | | |
| 3 | P45 | -4530.  | | |

Units are not printed; they are 1/sec (the modes of `ẋ = Ax`).

Mode descriptions [pdf p.29]: "The five-degree-of-freedom extracted model consists of five
stable, nonoscillatory modes. The power turbine state is completely decoupled from the
other states. … Three modes, corresponding to the pressure states, are outside the
bandwidth of interest for the modeling of rotor and propulsion system dynamics."

The "power turbine state is completely decoupled" claim is visible in every matrix
transcribed above: **column 2 (NP) is identically zero except for the diagonal element**,
in the 5-DOF, the 6-DOF and the 3-DOF alike.

### 6.3 Independent check of this transcription (performed here, not from the report)

The eigenvalues of the three transcribed 5×5 `A` matrices were computed and compared to
Table 1. A single dropped minus sign anywhere would destroy this agreement:

| Trim | computed from transcription | Table 1 (pdf p.31) |
|---|---|---|
| 1 | -4899.80, -3062.36, -51.65, -2.668, -0.5650 | -4900., -3060., -51.6, -2.66, -0.565 |
| 2 | -4645.31, -4038.04, -52.23, -2.087, -0.4461 | -4640., -4040., -52.2, -2.08, -0.446 |
| 3 | -4525.72, -4432.49, -52.63, -1.753, -0.3567 | -4530., -4430., -52.6, -1.75, -0.357 |

Agreement to 3 significant figures at all three trims. **The three 5-DOF matrices are
therefore verified end-to-end, not merely read twice.** The 2-DOF matrices are
lower-triangular, so their diagonals are their eigenvalues, and those match Table 1's
2-DOF column exactly (-2.693/-2.233/-1.816 vs -2.69/-2.23/-1.82).

**Caution on the 6-DOF matrices.** As transcribed, the 6-DOF `A` at trim 1 has an
*unstable* eigenvalue, +3.07 /sec (trims 2 and 3 are stable, with slow modes -0.095 and
-0.391). Every element of the trim-1 6×6 was re-read three times; the transcription is
what is printed. The explanation is loss of precision in the printed values, not a typo:
`A` row 6 is `F1⁻¹F2` row 6, a difference of terms of order 2×10⁵ that produces a result
of order 10³, so the printed 4-significant-digit values carry ±50 absolute uncertainty in
`A(6,3)` and `A(6,4)`. Perturbing those two elements **within their printed rounding**
moves the slow mode anywhere from -0.99 to +7.3 /sec.

> **Consequence for validation:** compare our extracted 6-DOF derivatives to Figs. B8/B10/B12
> **element by element only**. Do not compare eigenvalues of the 6-DOF models, and do not
> treat the printed trim-1 instability as a claim by the report. The 5-DOF and 2-DOF
> matrices, by contrast, are well conditioned and their eigenvalues are a valid check.

### 6.4 Time constants and the heat-sink transfer function [pdf p.31, Eq. 63]

Appendix B prints no time constants. The main text gives the heat-sink transfer function:

```
  T41 / T41_ns = (τ1 s + 1) / (τ2 s + 1)                            (Eq. 63, pdf p.31)
```

with τ1 the heat-sink lead and τ2 the heat-sink lag time constant at a trim condition,
both in sec [nomenclature pdf p.13]. **Numerical values of τ1 and τ2 are not printed
anywhere in the report** (searched pp.28–34 and Appendix B). Note that `A(3,3)` of the
3-DOF models is *not* simply `-1/τ2`; it is the (3,3) element of `F1⁻¹F2`.

### 6.5 Extraction method, for reproducing these matrices [pdf p.28]

Confirming the brief: "Stability derivatives were extracted from the nonlinear simulations
by suppressing all integrations after trimming the simulation at a desired operating
condition. Each state or control was then individually perturbed … The central-difference
extraction method that was used perturbs each state or control variable in both positive
and negative directions. A perturbation step-size of plus or minus two percent of
equilibruim [sic] conditions was considered to be adequate…"

Also relevant to reproducing row 2 (NP): the NP row is *not* a plain Jacobian row. Eqs.
57–62 [pdf pp.28–29] build it from the load-torque derivatives ∂Qreq/∂NP and ∂Qreq/∂ṄP,
which come from the **Gen Hel UH-60A blade-element simulation** (an external program), and
give
`A(2,j) = (∂Qeng/∂xj)/(J + ∂Qreq/∂ṄP)` for j = 1, 3–5, and
`A(2,2) = (∂Qeng/∂NP − ∂Qreq/∂NP)/(J + ∂Qreq/∂ṄP)`, `b(2) = (∂Qeng/∂Wf)/(J + ∂Qreq/∂ṄP)`.
Our model will not reproduce row 2 or `b(2)` without those external load derivatives.

---

## 7. Ambiguous / flagged items

Every matrix element in §4 and §5 was resolved unambiguously at 300 dpi. **No matrix cell
is marked ambiguous.** The flags are all in the surrounding text:

| # | Item | Page | Status |
|---|---|---|---|
| 1 | Table B.1 row label `P_s3, PSIA` — subscript glyph | pdf p.67 | **Ambiguous glyph.** Candidates: `s` (strongly favoured; `Ps3` = station 3 *static* pressure is a real nomenclature entry, pdf p.11) or `a`. If `s` is right, this row is *not* the trim value of the state `P3` (station 3 total pressure) and Appendix B gives no trim value for that state. |
| 2 | Table B.1 unit `rpm M` for gas generator speed | pdf p.67 | Printed as shown, verified at 2× zoom. Stray italic `M`; unit is rpm. Recorded as printed. |
| 3 | Symbol `F` in `b̄ = F1⁻¹G1 + F F1⁻¹G2` | pdf p.33 | Undefined in the report. Probably `A = F1⁻¹F2` (standard Chen form). Transcribed as printed; do not assume. |
| 4 | 3-DOF state ordering `{NG, NP, T41}` | pdf pp.71/73/75 | Not printed. Inferred from `d̄` structure + `G2` (pdf p.33) + the NP diagonal. High confidence, but not a printed fact. |
| 5 | Units of the `A`, `b̄`, `d̄` entries | pdf pp.68–76 | Not printed in Appendix B. Derived from Eqs. 66–67 and the nomenclature (§2.6). |
| 6 | 6-DOF trim-1 `A` is open-loop unstable as printed (+3.07 /sec) | pdf p.72 | Transcription verified 3×. Rounding artifact of the printed 4-digit precision, not a report claim — see §6.3. |
| 7 | 6-DOF `A(5,3)` is exactly `0.0000E+0` at trim 1 but `-0.4315E+1` / `-0.7397E+1` at trims 2/3 | pdf pp.72/74/76 | Verified; recorded as printed. Likewise `A(2,3) = 0.0000E+0` in all three 6-DOF models while the 5-DOF `A(2,3)` is large and nonzero. |
| 8 | 5-DOF `A(2,1)` sign flips with trim: `-0.2444E-1` (T1), `+0.1283E-1` (T2), `+0.3555E-1` (T3) | pdf pp.68–70 | All three verified at 2× zoom. Genuine sign change, not a scan artifact. |

---

## 8. Summary of what was transcribed

- 1 table (Table B.1): 20 rows × 3 columns + 3 footnotes.
- 3 × 2-DOF: `A` 2×2, `b̄` 2×1.
- 3 × 5-DOF: `A` 5×5, `b̄` 5×1.
- 3 × 3-DOF: `A` 3×3, `b̄` 3×1, `C = [I]`, `d̄` 3×1.
- 3 × 6-DOF: `A` 6×6, `b̄` 6×1, `C = [I]`, `d̄` 6×1.

**30 numeric arrays across 12 figures** (plus 6 symbolic `C = [I]`), and Table 1
(pdf p.31) transcribed as supporting data.
Total numeric matrix elements: 3×(4+2) + 3×(25+5) + 3×(9+3+3) + 3×(36+6+6)
= 18 + 90 + 45 + 144 = **297**.
