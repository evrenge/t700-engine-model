# Every figure in the report, and what we did with it

Written 2026-09-16, after two corrections that had been applied unevenly:

* **A — straightness measured in data space.** `tools/check_straightness.py`. Whether the
  printed curve bends between its markers decides whether `np.interp` loses anything. The
  first version of that check measured the bow in *pixels*, which is wrong on every page
  here: the frames are tilted quadrilaterals (opposite edges differ by up to 0.356 deg on
  Figure A8) and the axis maps are polynomials (0.20-0.25 % of span). Pushing the traced
  centreline through the figure's own calibration first removes both.
* **B — the two ticked edges blended, not pooled.** Each axis has two edges carrying
  ticks. Fitting them together averages away a real difference; `digitize_native.blend_axis`
  has always fitted each alone and blended across the plot, and now everything does.

The table records, for every figure the report prints: is it a plot, do we use it, did A
and B reach it, and if not, why not. **The "why not" column is the point of this file** --
a correction applied to most of a set and silently not to the rest is worse than one not
applied at all, because the set looks uniform.

## Main body -- Figures 1 to 15

| Fig | pdf p. | Plot? | Used? | A | B | Note |
|---|---|---|---|---|---|---|
| 1 | 20 | chart | no | n/a | n/a | modal frequencies from ref. 3, not T700 data. The errata on pdf p.6 corrects a label here and nothing numeric |
| 2 | 21 | drawing | no | n/a | n/a | engine cutaway |
| 3 | 22 | drawing | no | n/a | n/a | turboshaft component schematic |
| 4 | 30 | plot | **no** | n/a | n/a | Gen Hel UH-60A vs linear models. Needs the rotor simulation the report consumes and does not contain -- out of scope by the report's own boundary |
| 5 | 34 | plot | **no** | n/a | n/a | 6-DOF vs 5-DOF linear models, same reason |
| **6** | 40 | plot | **yes** | n/a | n/a | NG vs fuel flow. `digitize_fig678`. Comparison data, not a function table: there is no interpolation between markers, so A does not apply. Its axis map is LINEAR anchored on the printed frame values, with the interior ticks held out as the check -- no axis polynomial exists to blend, so B does not apply. Tilt and shear ARE corrected: all four frames fitted, `x_of`/`y_of` blend each opposite pair |
| **7** | 41 | plot | **yes** | n/a | n/a | shaft horsepower vs fuel flow, as Figure 6 |
| **8** | 42 | plot | **yes** | n/a | n/a | Ps3 vs NG, as Figure 6 |
| **9** | 45 | plot | **yes** | n/a | n/a | fuel step 400->775, six stacked panels. `digitize_fig910`. Per-panel frames, both fitted with slopes, so tilt is corrected; linear axis maps anchored on frames. Held-out tick check: worst interior major 0.0154 s on a 0-5 s axis |
| **10** | 46 | plot | **yes** | n/a | n/a | fuel step 400->125, as Figure 9. Held-out tick check 0.0013 s |
| 11 | 48 | plot | **no** | n/a | n/a | closed-loop collective at high forward speed. Gen Hel rotor input |
| 12 | 49 | plot | **no** | n/a | n/a | closed-loop lateral cyclic. Gen Hel rotor input |
| 13 | 51 | plot | **no** | n/a | n/a | closed-loop collective at hover. Gen Hel rotor input |
| 14 | 52 | plot | **no** | n/a | n/a | closed-loop pedal at hover. Gen Hel rotor input |
| 15 | 53 | plot | **no** | n/a | n/a | closed-loop large collective. Different test data source from Figs 11-14 |
| Table 1 | 31 | table | **yes** | n/a | n/a | 27 printed eigenvalues, transcribed |

## Appendix A -- Table A.1 and Figures A1 to A11

Every one is a plot and every one is used: these are the engine's function tables.

| Fig | pdf p. | A | B | Note |
|---|---|---|---|---|
| Table A.1 | 55 | n/a | n/a | 33 scalar constants, transcribed |
| **A1** `f1` | 56 | **no** | **yes** | The beta map is derived from the raw table, not traced from ink, so `check_straightness` has nothing to trace; `digitize_a1` reports its own per-speed-line fit residual instead (0.23-0.68 px). B applied 2026-09-16 -- and A1 settles pooled-vs-blended by itself: its eleven flat extensions are horizontal lines spanning the plot, and their end-to-end difference falls from 0.573 px mean to 0.409 px |
| **A2** `f2` | 57 | **no** | **already** | Header does not publish the frame arcs, so the check cannot rebuild the map. Validated instead by a stronger free check: the 20 abscissae land on the integer lattice to 0.0155 max, and the fit was never told to expect integers. Blended since it was written -- `digitize_native` |
| **A3** `f3` | 58 | yes | yes | worst bow 0.79 sigma |
| **A4** `f4` | 59 | yes | yes | worst bow 0.91 sigma |
| **A5** `f5` | 60 | yes | yes | worst bow 0.61 sigma |
| **A6** `f6` | 61 | yes | yes | worst bow 0.02 sigma. B collapsed its endpoint difference from 0.18 sigma to 0.01 -- independent support for #62's reading that the drawn slope is page skew |
| **A7** `f7` | 62 | yes | **already** | worst bow 1.04 sigma. Blended since it was written |
| **A8** `f8` | 63 | yes | yes | worst bow 0.16 sigma |
| **A9** `f9` | 64 | **no** | **already** | Header does not publish the frame arcs. Validated by the tool's own verification overlay. Blended since written |
| **A10** `f10` | 65 | yes | yes | worst bow 0.24 sigma. The printed axis LABEL is inverted; the stored values are `f10 = P49/Ps9` directly -- open question #5 |
| **A11** `f_hs` | 66 | yes | yes | worst bow 0.03 sigma |

## Appendix B -- Figures B1 to B12

| Fig | pdf p. | Plot? | Used? | A | B | Note |
|---|---|---|---|---|---|---|
| B1-B12 | 68-76 | **no** | **yes** | n/a | n/a | Printed state-space matrices, not plots. 297 elements transcribed into `data/linear/b01..b12.csv`; 206 non-zero elements carry a numeric per-element comparison. Nothing is digitized, so neither correction applies |
| Table B.1 | 67 | table | **yes** | n/a | n/a | the three trim conditions, transcribed |

## Appendix C -- Table C.1 and Figures C1 to C30

| Fig | pdf p. | Plot? | Used? | A | B | Note |
|---|---|---|---|---|---|---|
| C1-C22 | 79-93 | **no** | **yes** | n/a | n/a | Block diagrams of the fuel control. Read to build the control model; nothing to digitize. Appendix C prints no numbered equations at all |
| Table C.1 | 83-84 | table | **yes** | n/a | n/a | 57 control constants, transcribed |
| **C23** `F_EC1` | 94 | plot | **yes** | **no** | yes | Traced multi-curve, identified by rank in y; not a marker-to-marker table, so the straightness question does not arise. Validated by the gap-threshold sweep (<=0.05 in y). B applied 2026-09-16 |
| **C24** `F_HM1` | 95 | plot | yes | yes | yes | worst bow 0.34 sigma |
| **C25** `F_HM2` | 96 | plot | yes | yes | yes | worst bow 0.12 sigma |
| **C26** `F_HM3` | 97 | plot | yes | yes | yes | worst bow 0.06 sigma |
| **C27** `F_HM4` | 98 | plot | yes | yes | yes | worst bow 0.22 sigma |
| **C28** `F_HM5` | 99 | plot | yes | yes | yes | worst bow 0.18 sigma |
| **C29** `F_HM6` | 100 | plot | yes | yes | yes | worst bow 4.18 sigma -- the only figure above its own read error. Two printed markers, so Ballin's table is two points and the bow is the draughtsman's pen; there is nothing in the report to represent it with |
| **C30** `F_HM7` | 101 | plot | **yes** | **no** | yes | Seven curves, 192 vertices, extracted as vertices rather than traced, so there is no centreline to test. Validated by requiring every segment between consecutive knots to lie on ink. B applied 2026-09-16; coverage 0.9850 -> 0.9879. The audit then found the column selector still weighting the blend at a fixed mid-plot 0.5, worth up to 4.96 columns, and fixed it -> 0.9881 |

## Summary

* **19 digitized figures.** B now reaches all 19 -- 13 on 2026-09-16 (`8048c60`), C23 and
  C30 the same day (`67ed8fa`), f1 the same day (`c448f51`), and f2, f7, f9 always had it.
* **A covers 14 of them.** The five it does not are not oversights: f1's beta map is
  derived rather than traced, f2 and f9 do not publish their arc coefficients in the
  header, and C23 and C30 are not marker-to-marker tables. Each carries a stronger check
  of its own, named above.
* **Figures 6-10 of the main body are comparison data, not function tables.** Neither
  correction applies to them in the form applied to the maps, and both their tools already
  fit both frames of each pair, so tilt and shear are corrected. What they do NOT have is a
  polynomial axis map: they are linear between printed frame values, with the interior
  ticks held out as the check. Figure 9's own check reports a worst interior major tick
  0.0154 s off a 0-5 s axis -- 0.31 %, the same order as the bow the Appendix A figures
  needed a cubic to remove. **That is the one place in this table where a correction is
  arguably owed and has not been made**, and it matters because these five figures are the
  validation targets.

## What the audit itself turned up

Writing the table was not bookkeeping; grepping every producer for a blend call found two
places where the correction had been applied to the values and not to the coordinate that
selects them:

* `digitize_c30.py` built its column lookup once down the middle of the plot
  (`blend_value(CX, uu, 0.5)`). Every reported coordinate was already recomputed honestly
  from the chosen column, so no stored number was wrong -- but the *choice* of column was
  off by up to 4.96 columns for a knot near either edge. Fixing it needed the blend weight
  without a row, since finding the row is what the lookup is for; inverting the y
  calibration gives it directly, because the homography has already normalized the
  vertical coordinate. Segment-on-ink 0.9879 -> 0.9881, worst seed miss unchanged at 0.066.
  An intermediate attempt that estimated the row with the tool's straight pixel map made
  things worse (0.9848, seed miss 0.106): that map ignores the column, so on a frame this
  tilted it is wrong by the full tilt at the right-hand end.
* `digitize_a1.py` reported the printed construction line as a single number. It is a
  vertical line under a blended x map, so it is a band: 0.9913 at the bottom scale, 0.9895
  at the top. The header now prints both. The spread is 0.0018, which is the useful part --
  it says the 0.9904-against-0.974 comparison does not rest on which edge you read.

Two other fixed weights survive and are correct: `_vgrid` in `digitize_c30` inverts the y
map to produce a blend weight, so its own weight is second-order, and `_GX`/`_GY` in
`digitize_a1` place the tracing guide, which is documented at the site as approximate and
never produces a stored value.
