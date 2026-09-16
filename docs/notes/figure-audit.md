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
| **6** | 40 | plot | **yes** | n/a | **yes** | NG vs fuel flow. `digitize_fig678`. Comparison data, not a function table: no interpolation happens between markers, so A does not apply. Tilt was always corrected -- all four frames fitted. The axis maps were **linear** until 2026-09-16, when the audit tested the held-out ticks for shape rather than size and found the ordinates bow; they are now anchored polynomials, per edge, blended across |
| **7** | 41 | plot | **yes** | n/a | **yes** | shaft horsepower vs fuel flow. `digitize_fig678`. Comparison data, not a function table: no interpolation happens between markers, so A does not apply. Tilt was always corrected -- all four frames fitted. The axis maps were **linear** until 2026-09-16, when the audit tested the held-out ticks for shape rather than size and found the ordinates bow; they are now anchored polynomials, per edge, blended across |
| **8** | 42 | plot | **yes** | n/a | **yes** | Ps3 vs NG. `digitize_fig678`. Comparison data, not a function table: no interpolation happens between markers, so A does not apply. Tilt was always corrected -- all four frames fitted. The axis maps were **linear** until 2026-09-16, when the audit tested the held-out ticks for shape rather than size and found the ordinates bow; they are now anchored polynomials, per edge, blended across |
| **9** | 45 | plot | **yes** | n/a | **no** | fuel step 400->775, six stacked panels. `digitize_fig910`. Per-panel frames, both fitted with slopes, so tilt is corrected. Tested for a bow the same way and does **not** have one -- the residuals are largest at t=1 and t=4 and smallest mid-panel, which is the opposite of a bow, so the linear map stands. The tick detector was fixed here instead |
| **10** | 46 | plot | **yes** | n/a | **no** | fuel step 400->125, six stacked panels. `digitize_fig910`. Per-panel frames, both fitted with slopes, so tilt is corrected. Tested for a bow the same way and does **not** have one -- the residuals are largest at t=1 and t=4 and smallest mid-panel, which is the opposite of a bow, so the linear map stands. The tick detector was fixed here instead |
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

* **19 digitized figures.** B now reaches all 19 -- 13 on 2026-09-16, C23 and C30 the same
  day, f1 the same day, and f2, f7, f9 always had it.
* **A covers 14 of them.** The five it does not are not oversights: f1's beta map is
  derived rather than traced, f2 and f9 do not publish their arc coefficients in the
  header, and C23 and C30 are not marker-to-marker tables. Each carries a stronger check
  of its own, named above.
* **Figures 6-10 are comparison data, not function tables**, so A does not apply to them in
  the form it takes on the maps. Both tools always fitted every frame, so tilt and shear
  were never the problem. What they did have was a **linear** axis map, and the audit is
  what caught it.

## What the audit itself turned up

Writing the table was not bookkeeping. Four things came out of it.

**1. C30 selected its columns with an unblended coordinate.** `digitize_c30` built its
column lookup once down the middle of the plot (`blend_value(CX, uu, 0.5)`). Every reported
coordinate was already recomputed honestly from the chosen column, so no stored number was
wrong -- but the *choice* of column was off by up to 4.96 columns for a knot near either
edge. The blend weight has to arrive without a row, since finding the row is what the
lookup is for; inverting the y calibration gives it directly, because the homography has
already normalized the vertical coordinate. Segment-on-ink 0.9879 -> 0.9881. An
intermediate attempt that estimated the row from the tool's straight pixel map made things
worse (0.9848, seed miss 0.066 -> 0.106): that map ignores the column, so on a frame this
tilted it is wrong by the full tilt at the right-hand end.

**2. A1 reported a band as a number.** The printed construction line is vertical under a
blended x map, so it reads 0.9913 at the bottom scale and 0.9895 at the top. The header now
prints both. The spread, 0.0018, is the useful part: it says the 0.9904-against-0.974
comparison does not depend on which edge you read.

**3. Figure 9 and 10's tick check was looking in the wrong place.** Its band ran rows 4 to
16 above the frame. Counting inked columns row by row on the p.45 PCNG panel: the frame's
own antialiased shoulder puts 124-344 of 1672 columns in rows 2-5, the ticks put 35-45 in
rows 6-11, and rows 12 and up hold only the trace. The band swallowed the shoulder at one
end and the trace at the other, so it found 62 groups across the twelve panels of which
only 14 were near an interior major, and reported a contaminated 0.1079 s. Rows 6-11 find
35 of 48, worst 0.0154 s. **No digitized value changed** -- the tick check only ever fed
the header -- but the number in that header is now real.

**4. Figures 6, 7 and 8 were read off a straight axis that is not straight.** This is the
one that moved data. The old check asked only how *big* the held-out tick residuals were;
asking what *shape* they had settles it in one number. A printed scale that genuinely bows
gives residuals pinned near zero at both frames and one-signed in between -- correlated
with t(1-t). Figure 8's ordinate correlates -0.82 (left frame) and -0.79 (right), 21 of 22
residuals share a sign, peaking at 0.36 % of span mid-axis and tapering to 0.06 % at both
frames. That is the same 0.2-0.4 % of span that made every Appendix A axis a cubic.

The fix keeps the frames pinning the map and lets the ticks say only how it bends between
them: `value = linear(t) + t(1-t)P(t)`, where the `t(1-t)` factor vanishes at both frames,
so the printed frame values stay exact by construction. The degree of `P` is chosen by
leave-one-out -- refit without each tick, predict it -- which buys back the check that
fitting the ticks would otherwise spend. Held-out rms, the honest number:

| edge | straight | after the bow | |
|---|---|---|---|
| Fig 6 left / right (%NG) | 0.0453 / 0.0547 | 0.0247 / 0.0238 | -45 % / -56 % |
| Fig 7 left / right (SHP) | 1.951 / 3.447 | 1.178 / 1.377 | -40 % / -60 % |
| Fig 8 left / right (psia) | 0.441 / 0.471 | 0.106 / 0.145 | -76 % / -69 % |

The abscissae bow far less (strongest -0.35 against the ordinates' -0.82) and leave-one-out
still prefers a small bend on each; that is the right outcome to accept rather than
overrule, since the test is whether a bend predicts ticks it was not shown.

Against the model, Figure 8 is the vindication: rms 0.45509 -> 0.26557 (-42 %) and its
**bias all but vanishes**, -0.34615 -> +0.05793. The model had been reading systematically
low against Figure 8 and it was the figure that was misread. Figure 6 is unchanged in rms
(0.23163 -> 0.23125) with its bias likewise collapsing to zero.

**Figure 7's percentage rms gets worse, 0.59357 -> 0.67353, and that is worth stating
plainly rather than burying.** Scored in horsepower instead, the same points are
essentially unchanged: rms 1.881 -> 1.931 hp over a 147-1789 hp range. What actually
happened is that a -0.94 hp mean bias became +0.96 hp, and Figure 7's deviation is measured
*relative to SHP*, which the figure carries down to 147 hp -- so a one-horsepower shift is
1.0 % there and 0.02 % at the top. The per-point change is monotone in exactly that way:
+1.030 at the lowest fuel flow, +0.018 at the highest. And the bias of either sign is
smaller than the figure resolves: Figure 7's own held-out read error is 1.38 hp. The CSV
headers now publish that number, so a deviation can be judged against what the figure can
actually say.

Two fixed blend weights survive and are correct: `_vgrid` in `digitize_c30` inverts the y
map to produce a blend weight, so its own weight is second-order, and `_GX`/`_GY` in
`digitize_a1` place the tracing guide, documented at the site as approximate and never
producing a stored value.

## Branch against main, as it stands

143 of the 214 numbers in `site/assets/report.json` differ. Main's `report.json` and its
`data/` were committed together in `ab96c34`, so the comparison is like for like. The
branch carries both the digitizing work in this file and the solver work that preceded it
(the degenerate trim root), and the table separates them where it can.

| Measure | main | branch | |
|---|---|---|---|
| **Table B.1 trims** rms | 0.24389 | **0.21927** | -10.1 % |
| Table B.1 worst | -0.90042 | **-0.79721** | -11.5 % |
| **Table 1 eigenvalues** rms | 5.49377 | **5.37007** | -2.3 % |
| **Figure 6** rms | 0.23369 | **0.23125** | -1.0 % |
| Figure 6 mean bias | 0.00030 | **-0.00054** | to zero |
| **Figure 7** rms | 0.60617 | 0.67353 | +11.1 %, see below |
| **Figure 8** rms | 0.42468 | **0.26557** | **-37.5 %** |
| Figure 8 mean bias | -0.31676 | **0.05793** | -82 % |
| **Closed loop** rms | 0.09813 | **0.09423** | -4.0 % |
| **Phase plane** fig 9 rms | 0.78108 | **0.67873** | -13.1 % |
| Phase plane fig 10 rms | 1.24526 | 1.25046 | +0.4 % |
| **Whole curve** fig 9 (5 channels) | 1.47388 / 0.76013 / 0.70575 / 0.83566 / 0.67254 | **1.41089 / 0.71655 / 0.66291 / 0.82549** / 0.69556 | 4 better, torq45 +3.4 % |
| Whole curve fig 10 (5 channels) | 1.44656 / 1.69841 / 3.81958 / 2.93741 / 2.42067 | **all five lower** | -0.2 to -1.9 % |
| **Appendix B** *A* rms | 7.76001 | 7.82576 | +0.9 % |
| Appendix B *b* rms | 5.22201 | 5.37290 | +2.9 % |
| 2-DOF *A* / *b* rms | 4.92973 / 2.32987 | **4.88433 / 2.24455** | -0.9 % / -3.7 % |
| 3-DOF *A* / *b* rms | 6.03312 / 6.27824 | 7.48801 / 6.81401 | +24.1 % / +8.5 %, open question #73 |
| 5-DOF *A* / *b* rms | 6.90878 / 0.08863 | **6.46452 / 0.07516** | -6.4 % / -15.2 % |
| 6-DOF *A* / *b* rms | 8.96979 / 7.10823 | 9.02510 / 7.13810 | +0.6 % / +0.4 % |
| **Steady sweep** points solved | 159 | **170** | +6.9 % |

### The three that moved the wrong way, and what each is

**Figure 7, +11.1 %.** Two separate things are folded into that number and neither is a
model regression.

*One:* main scored 15 points and the branch scores 14. The point main dropped is at 806.5
lb/hr, and main **trimmed it successfully** -- off the top of the compressor map. The
branch refuses it, which is the degenerate-root fix doing exactly what it was written for.
Main's 0.60617 includes a solve that should never have converged.

*Two:* on the 14 points both versions score, the bow correction moved rms 0.59357 ->
0.67353. In horsepower those same points are unchanged: 1.881 -> 1.931 hp over a 147-1789
hp range. A -0.94 hp bias became +0.96 hp, and Figure 7 scores *relative* to an SHP the
figure carries down to 147, so one horsepower is 1.0 % at the bottom and 0.02 % at the top.
Both biases are below the figure's own held-out read error of 1.38 hp.

**3-DOF blocks, +24.1 % on *A*.** This is open question #73, already known and already
carrying a deliberately red test. Not introduced here.

**Appendix B overall, +0.9 % and +2.9 %.** Driven by the 3-DOF block above; 2-DOF and
5-DOF both improve, 6-DOF is flat to within half a percent.

### The independent check nobody tuned

`internal_spread` measures **the report against itself**: Table B.1 prints numbers for the
three trims, Figures 6 and 7 plot the same three points, and they disagree. The Table B.1
side is transcribed and carries no digitizing error of ours at all, so moving the figure
closer to it is evidence the figure is being read better -- and this number was never a
target.

| trim, quantity | main | branch |
|---|---|---|
| hover NG | 0.0431 | **0.0024** |
| hover SHP | -0.2892 | **-0.0342** |
| level 80 kt NG | 0.0236 | **-0.0111** |
| level 80 kt SHP | -0.3068 | **0.1325** |
| descent 80 kt NG | 0.0331 | **0.0199** |
| descent 80 kt SHP | 0.7345 | 1.3062 |

Five of six shrink; excluding the descent SHP row the mean falls 0.139 -> 0.040, **-71 %**.
The one that grows is the descent trim's shaft power -- the lowest-power point on the
figure, which is both where the relative metric is worst conditioned and the subject of
#73.

## Was each correction implementable, and what did it give

Every plot the replication reads. **A** is data-space straightness, **B** the per-edge
blend. "n/a" is a reasoned no, not an omission; the reason is in the row.

| Plot | A possible | A result | B possible | B result (half the two edges' disagreement, x / y) |
|---|---|---|---|---|
| **A1** f1 | **no** -- beta map is derived, not traced | tool reports its own per-speed-line fit, 0.23-0.68 px | yes | applied. Free check: 11 flat extensions, end-to-end 0.573 -> **0.409** px mean, 1.005 -> 0.921 worst |
| **A2** f2 | **no** -- header omits the frame arcs | free check: 20 abscissae on the integer lattice, 0.0155 max | already | blended since written |
| **A3** f3 | yes | bow **0.79 sigma** | yes | 0.0117 / 4.18e-05 |
| **A4** f4 | yes | bow **0.91 sigma** | yes | 0.0045 / 4.16e-07 |
| **A5** f5 | yes | bow **0.61 sigma** | yes | 0.00475 / 1.98e-06 |
| **A6** f6 | yes | bow **0.02 sigma** | yes | 5.4e-06 / 3.91e-05. Endpoint difference 0.18 -> **0.01 sigma** |
| **A7** f7 | yes | bow **1.04 sigma** -- above its read error | already | blended since written |
| **A8** f8 | yes | bow **0.16 sigma** | yes | 0.000126 / 0.00914 |
| **A9** f9 | **no** -- header omits the frame arcs | tool's own verification overlay | already | blended since written |
| **A10** f10 | yes | bow **0.24 sigma** | yes | 0.0166 / 2.15e-05 |
| **A11** f_hs | yes | bow **0.03 sigma** | yes | 0.0206 / 0.00478 |
| **C23** F_EC1 | **no** -- traced multi-curve, not marker-to-marker | gap-threshold sweep, <=0.05 in y | yes | applied 2026-09-16 |
| **C24** F_HM1 | yes | bow **0.34 sigma** | yes | 0.0966 / 0.00171 |
| **C25** F_HM2 | yes | bow **0.12 sigma** | yes | 0.0534 / 0.00419 |
| **C26** F_HM3 | yes | bow **0.06 sigma** | yes | 0.0672 / 0.0142 |
| **C27** F_HM4 | yes | bow **0.22 sigma** | yes | 0.0505 / 0.000482 |
| **C28** F_HM5 | yes | bow **0.18 sigma** | yes | 0.0817 / 0.000127 |
| **C29** F_HM6 | yes | bow **4.18 sigma** -- see below | yes | 0.0162 / 0.00103 |
| **C30** F_HM7 | **no** -- vertices, not a trace | segment ink coverage **0.9881** | yes | applied; coverage 0.9850 -> 0.9879 -> **0.9881** |
| **Fig 6** | **n/a** -- comparison data, nothing interpolates between markers | -- | **yes, and it was owed** | held-out tick rms 0.0453/0.0547 -> **0.0247/0.0238** (-45 %/-56 %) |
| **Fig 7** | n/a, same | -- | **yes, and it was owed** | 1.951/3.447 -> **1.178/1.377** (-40 %/-60 %) |
| **Fig 8** | n/a, same | -- | **yes, and it was owed** | 0.441/0.471 -> **0.106/0.145** (**-76 %/-69 %**) |
| **Fig 9** | n/a, same | -- | **no, and tested** | residuals largest at t=1 and t=4, smallest mid-panel: the opposite of a bow. Linear map stands; tick detector fixed instead, 14 -> 35 of 48 ticks found |
| **Fig 10** | n/a, same | -- | **no, and tested** | as Figure 9 |

**Two figures sit above their own read error and neither is fixable from the report.**
C29 (F_HM6) at 4.18 sigma prints exactly two markers, so Ballin's table *is* two points and
the bow between them is the draughtsman's pen -- there is nothing in the report to
represent it with. A7 (f7) at 1.04 sigma is a genuine printed curvature between sparse
markers; representing it would mean inventing intermediate points the report does not give.
Both are recorded rather than smoothed, which is what the replication rule requires.

**Coverage, finally:** A reaches 14 of the 19 digitized figures and is reasoned-out on the
other 5. B reaches **all 19**, plus Figures 6-8 of the main body, and is reasoned-out on
Figures 9-10 by the same test that justified it elsewhere.
