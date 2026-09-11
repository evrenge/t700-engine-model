# Digitizing recipe

`tools/digitize.py` turns the report's printed plots into CSV. Nineteen figures in the
appendices carry model data (11 engine functions, 8 control schedules) and five more in
the body carry validation references. None of them exists as numbers anywhere.

Worked example, Figure A2 (`f2`, compressor temperature), which is the simplest case and
the one to copy.

## Four steps

```bash
# 0. render the page (300 dpi is enough; go to 600 where markers crowd)
pdftoppm -f 57 -l 57 -r 300 -png docs/ballin-tm100991.pdf validation/out/digitize/a2

# 1. detect the plot frame and get a ready-made calibration
#    the four values are what the printed axis labels say at each frame edge
python3 tools/digitize.py axes --page validation/out/digitize/a2-057.png \
    --x0 0 --x1 20 --y0 1.1 --y1 3.0
#    -> --calib "591:0.0,2255:20.0,2636:1.1,293:3.0"
#    -> --roi 603,305,2243,2624

# 2. sweep the threshold until the count is stable, and check it against the
#    independent visual count in inventory-appendix-a.md
python3 tools/digitize.py marks --page ... --roi ... --min-fill 0.40

# 3. extract
python3 tools/digitize.py extract --page ... --roi ... --min-fill 0.40 --calib ... \
    --figure A2 --pdf-page 57 --quantity "f2 -- compressor temperature" \
    --x-label "..." --y-label "..." --out data/maps/f2_compressor_temperature.csv

# 4. VERIFY, and look at the image
python3 tools/digitize.py verify --page ... --csv ... --calib ... --roi ... \
    --out validation/out/digitize/a2_verify.png
```

Step 4 is not optional. Both of the failures below were invisible in the CSV and obvious
in the overlay.

## Two failure modes, both hit on the first figure

**1. Markers sit ON the joining line, so blob detection is useless.** The curve and all
twenty of its crosses form a single connected component; `find_marks` returns one
meaningless centroid for the lot. Use `--mode density` (the default): a marker is locally
*denser* than the line running through it, so box-filtering the ink mask and taking local
maxima finds them. `find_marks` is still there for genuinely isolated markers.

**2. Calibrating by eye off the grid overlay put every point out by 0.3.** Reading tick
pixel coordinates from a downscaled image cost ~25 px, which on Figure A2 is 0.3 in
PS3/P2 — large enough to move every value, small enough to look entirely plausible in the
CSV. **Always calibrate with `axes`, never by eye.** That is what the subcommand is for.

A third, smaller one: a marker sitting exactly on the frame gets clipped by the `axes`
default 12 px inset. Figure A2's last point, at PS3/P2 = 20, needed the ROI extended past
the right frame. If the count comes up one short, look at the corners first.

## Checking a result without trusting the tool

Figure A2's markers were drawn at each integer of PS3/P2. The extraction recovers
`1.02, 2.03, 3.00, 4.01, 5.00, ... 19.00, 19.99` — integers to within 0.04, and the tool
was never told to expect them. When a figure has a regular structure like that, it is a
free check on the calibration. Use it where it exists.

Otherwise the checks are: the point count against the independent visual count recorded
in `inventory-appendix-a.md`, a stable plateau across several `--min-fill` values, and
the verify overlay.

## Provenance

Every CSV carries a header naming the figure, the PDF page, the method, and the settings.
`digitized` is weaker provenance than `transcribed` and the header says so. Per CLAUDE.md,
nothing enters `data/` without it.

## Status

| Figure | Quantity | Points | State |
|---|---|---|---|
| A2, pdf p.57 | `f2` compressor temperature | 20 | done, verified |
| A1, A3–A11 | `f1`, `f3`–`f10`, `f_hs` | ~170 | not started |
| C23–C30 | `F_EC1`, `F_HM1`–`F_HM7` | — | not started |
| 6–10, pdf pp.40–46 | validation references | ~550 | not started |

Figure A1 (11 crossing speed lines) and Figure C30 (7 curves bundling where markers
overprint) are the two that will need per-curve ROIs and probably 600 dpi.

### Update — A3, A4, A5, A6, A11 done (2026-09-10)

| Figure | Quantity | Points | Inventory said | min-fill | State |
|---|---|---|---|---|---|
| A3, pdf p.58 | `f3` seal-pressurization bleed | 13 | ~17 | 0.45 (stable 0.40–0.50) | done, verified |
| A4, pdf p.59 | `f4` power-turbine-balance bleed | 3 | 3 ✓ | 0.55 (stable 0.50–0.65) | done, verified |
| A5, pdf p.60 | `f5` tip leakage + turbine cooling bleed | 3 | 3 ✓ | 0.55 (stable 0.50–0.55) | done, verified |
| A6, pdf p.61 | `f6` combustor efficiency | 2 | 2 ✓ | 0.45 (stable 0.35–0.60) | done, verified |
| A11, pdf p.66 | `f_hs` station 4.1 heat-sink constant | 6 | 6 ✓ | 0.40 (stable 0.35–0.45) | done, verified |

All five at 300 dpi. Outputs in `data/maps/`; overlays in `validation/out/digitize/`.

**A third failure mode: `axes` mis-locates the frame when the page carries a long flat
data line.** `find_frame` averages every strong row in the top half and every strong row
in the bottom half. On Figure A3 the flat 0.11 run and the flat 0.00 run are strong rows,
and the bottom frame line is diluted by scan skew, so `axes` returned **top=362,
bottom=2826 against a measured 313, 2672** — 0.009 in B1, five times the error that
calibrating by eye caused on A2. It was invisible in the CSV. Every figure whose function
goes flat (A3, A4, A5, A6, A11 — five of the eleven) is exposed. **Check the frame `axes`
reports against the printed ticks before you use it.**

**Fourth: these pages are keystoned, not square.** The printed frame is not axis-parallel
in the render — the bottom frame line drops 10 px across the plot width on p.58, 12 px on
p.59, 16 px on p.60, 15 px on p.66, while the top drops a different amount (2 px on p.59).
No single row/column pair describes the frame, so any two-point calibration is wrong
somewhere on the page. Fix it in the image before digitizing: fit the four frame lines,
intersect for the corners, and map that quadrilateral onto its enclosing rectangle with a
PIL `PERSPECTIVE` transform. Residual drift afterwards is under 1 px on all five pages,
and `axes` then finds left/right/top exactly (the flat-line problem above still bites the
bottom). The rectified renders are saved beside the originals as `a3r-058.png` etc.

**Calibration used here, and why it beats the frame corners.** Least squares over *every*
printed major tick on the axis, with the two frame lines included as the extreme ticks —
8 to 19 points per axis instead of 2. Residuals came out ≤2 px in x and ≤5.7 px in y on
all five figures. The y residuals are not noise: they trace a smooth ±5 px bow (the scan's
paper curl), which is the accuracy floor of a two-point linear map on these pages. Worth
about 0.2 % of full scale, so it matters only where the figure *is* one number — on A6 the
line was calibrated against the two ticks that bracket it (1.00 and 0.98) rather than the
global fit, which moves the answer 0.9846 → 0.9851.

**Marker counts.** Four of the five matched the inventory. A3 did not: **13 markers, not
~17.** Recounted marker by marker at 4× zoom, and the extraction agrees — the inventory's
figure is an estimate (it is written "≈ 17") and is too high. The 13 x values come out on
integer NGc to within 0.09 (65, 78–85, 87, 88, 89, 100 — note 86 is genuinely absent),
which is the free check the section above asks for. A11 gives the same kind of check and
passes it harder: x on 65/70/75/80/85/100 and y on 4.00, 4.00, 5.75, 6.20, 7.50, 7.50.

**The `--roi` corner gotcha bit on all five.** Every one of these figures prints a marker
*on* a vertical frame line — nine markers across the five — so the ROI has to run past the
frame, not inside it. Going too far past picks up the y tick labels; the working ROI puts
its left edge about one marker half-width (17 px at 300 dpi) outside the frame, and the
label ink that leaks in is then rejected by the min-fill sweep rather than by the ROI.

**One thing to know before reading B1 off A3 by eye.** The printed `0.11` y tick *label*
on p.58 sits ~19 px (0.0011) below its own tick mark — the only label on the page that
does; all the others are within 5 px. The flat top of `f3` lines up with that displaced
label, which makes 0.110 look right. Measured against the tick it is **0.1089–0.1091**.
The digitized value follows the tick.

### Update — A7, A8, A9, A10 done (2026-09-10)

| Figure | Quantity | Points | Inventory said | dpi | min-fill | State |
|---|---|---|---|---|---|---|
| A7, pdf p.62 | `f7` gas-generator turbine energy | 6 | 6 ✓ | 300 (rectified at 600, halved) | 0.40 (stable 0.38–0.44) | done, verified |
| A8, pdf p.63 | `f8` power turbine energy | 12 | 12 ✓ | 300 | 0.42 (stable 0.40–0.46) | done, verified |
| A9, pdf p.64 | `f9` power turbine mass flow | 23 | ~23 ✓ | 300 | 0.40 (stable 0.38–0.45) | done, verified |
| A10, pdf p.65 | `f10` exhaust pressure loss | 18 | ~20 ✗ | 600 | 0.42 (stable 0.40–0.46) | done, verified |

All four went through `tools/rectify_page.py rectify` first, then density mode with
`--window 21 --min-sep 40` at 300 dpi (`42`/`40` at 600). All four needed the ROI run 14 px
*past* the left and right frame: every one of these figures prints its first and last
marker **on** a vertical frame line.

**These pages are skewed as well as keystoned.** Rotation and keystone are separate
defects and `rectify` fixes both; it is worth knowing how big the rotation alone is,
because it is what breaks frame *detection* rather than just calibration: **+0.77° on
p.62**, +0.21° on p.63, +0.37° on p.64, +0.39° on p.65. On p.62 that alone drifts the frame
22 px across the plot width, which is why `find_frame` there returns the x-axis *title*:
a line drifting across 22 rows puts less ink in any one row than a block of text does.
Measured residual keystone after `rectify`: 0.0 px on all four.

**`rectify` cannot do p.62 at 300 dpi.** `edge_fits` finds horizontal lines with
`long_lines`, which wants 80 % continuous ink, and p.62's bottom frame line is broken
enough at 300 dpi that only the top line is ever seen — so `edge_fits` returns the same
line twice, the four "corners" are degenerate and `np.linalg.solve` raises
`LinAlgError: Singular matrix`. At 600 dpi both lines are found. **Workaround: render
p.62 at 600, rectify there, and halve with Lanczos.** The geometry is then the 600 dpi
fit and the detection behaviour is 300 dpi's, which is what A7 needs (see below).

**Calibrating an axis whose frame edge carries no printed label.** `ticks` gives a
least-squares fit over the printed majors, but A7 (x-left *and* y-bottom), A8 (x-left),
A9 (x-left) and A10 (y-bottom) print no number at that edge, so the value there has to be
extrapolated. Two independent computed routes, and they agree:

* *Tick lattice.* Take the median spacing of the long ticks and count frame widths in it.
  A7's x lattice is 110.61 px per 0.01 and the frame spans **14.980** of them from the
  printed 0.35 → left edge = 0.20. A8's is 151.0 px per 0.05, span **11.020** → 0.30.
  A9's, 150.98 px, span **10.988** → 0.30. A10's x spans **6.99** steps of 5 % against
  printed 65 and 100 at *both* edges — a pure check, and it passes at 34.96 % against 35.
* *Printed tick labels.* Take the text band nearest the frame (the band further out is the
  axis title) and fit value against label centroid. A10's 7-label fit reproduces the
  printed 1.14 at the top frame **to 1.4×10⁻⁴** and puts the unlabelled bottom frame at
  1.01023. A9's 8-label fit reproduces *both* printed edges, 0.39 to 1×10⁻⁴ and 0.25 to
  2×10⁻⁴.

**Label centroids carry a ~4–5 px systematic** — glyphs like `52.5` put the period's ink
low, which drags the centroid down; the A3 note above hit the same thing from the other
side. So use the label fit for the **scale** and a printed frame-edge value for the
**offset**, never a label centroid as an anchor. On A8, where both y frame edges *are*
printed, the 19-label fit lands 0.21 BTU/LBM off at the top and 0.078 off at the bottom —
0.17 % of range — which is exactly this bias, and is why the printed frame values are used
and the label fit only checks them.

**Free structural checks, both passed.** A8's 12 recovered abscissae come out evenly
spaced at 0.0500 ± 0.0008 — exactly 0.30 to 0.85 in steps of 0.05. A9's 23 come out at
0.0250 ± 0.0006 — 0.30 to 0.85 in steps of 0.025. The tool was told neither.

**How much rectifying was worth.** All four were first done on pages corrected for
*rotation only*, then redone on rectified pages. Counts were identical; values moved by
0.11 % (`f8`), 0.17 % (`f9`), 0.28 % (`f7`) and 0.42 % (`f10`) of the y range. Small, real,
and invisible in the CSV — which is the point.

**A10's marker count disagrees with the inventory: 18, not ~20.** The inventory wrote
"≈ 20" and flagged the trough count as a judgement call, so this is the flagged case
resolving downward, not a detection failure. Verified by reading the whole curve at 600 dpi
in three overlapping crops — 7 markers on the 65–85.5 % descent, 8 from 85.5 to 94 %
through the trough, 5 from 94 to 100 % on the rise, sharing endpoints → 18. Nothing is
merged and nothing is missed; an independent 300 dpi render finds the same 18 and agrees
on y to 6×10⁻⁵.

**A7 was flagged for 600 dpi, and 300 dpi is the right answer for the *markers*.** 600 dpi
resolves the five knee markers no better than 300 (they are 55 px apart at 300 with a
~22 px marker) and it introduces a trap: where the curve leaves the first marker it runs
16 px from the left frame, and that wedge of two near-parallel lines is *denser* than the
marker itself at 600 dpi. At `--window 46 --min-fill 0.40` the 600 dpi render returns a
stable, plausible, **wrong** six — the wedge in, the real first marker out, moving that
point by 0.66 BTU/LBM. At 300 dpi the marker outscores the wedge, and `--min-sep 100`
(half-window 50, under the 55 px marker spacing) merges the pair onto the marker. So:
rectify at 600 because `rectify` needs it, digitize at 300 because the detector needs it.
**A count that is stable is not the same as a count that is right.** Read the overlay.

**A8 and A9 are not on the same abscissa.** Facing pages, identical 0.30–0.85 range,
near-identical axis titles, and different independent variables: A8 is `P49/P45`
(total/total, Eq. 32), A9 is `Ps9/P45` (**static** station 9 over total, Eq. 33). Both CSVs
say so in their `x:` header, at length, because conflating them yields a plausible and
wrong model.

**`f8` goes negative** above `P49/P45` ≈ 0.83, reaching −1.32 at 0.85. That is data; the
printed y axis is extended to −5.0 to show it.

**A10 is stored exactly as printed — `PS9/P49`, 1.02–1.14 — and is *not* inverted.** This is
open question #5, closed by decision: the CSV stays faithful to the page and the consumer
takes the reciprocal, `f10 = P49/PS9 = 1/y`, at the point of use with that row cited beside
the `1/x`. The `y:` header says so in as many words. Using the stored numbers directly as
`f10` inverts the exhaust pressure loss.

### Update — A1 done (2026-09-10), the 2-D map

| Figure | Quantity | Points | Inventory said | dpi | State |
|---|---|---|---|---|---|
| A1, pdf p.56 | `f1` compressor mass flow, `WA2c = f1(Ps3/P2, NGc)` | 77 | ≈77 ✓ | 600 (rectified) | done, verified |

Output `data/maps/f1_compressor_mass_flow.csv`: **five columns**, `ngc_pct, ps3_p2,
wa2c_lbm_per_s, k, symbol`. 11 speed lines × 7 points. The speed-line parameter values are
printed in an in-plot legend and are `65, 80, 82, 85, 87, 89, 92, 94, 96, 98, 100 % NGC`
for symbols 1–11. Extraction is in **`tools/digitize_a1.py`** — this figure needed enough
bespoke logic that a shell recipe could not reproduce it; the script reruns and reproduces
the CSV byte for byte.

**The `--roi` per-curve advice above is wrong for this figure, and the right answer is
better.** At 600 dpi each speed line — its flat extension, its joining polyline and all
seven digit markers — is a **single 8-connected ink component**, and the eleven components
are disjoint. `ndimage.label` separates the family exactly, including curves 2–4 in the
crowded lower left where no ROI band can. Sort the eleven largest components by mean `y`
and you have symbol 11 down to symbol 1. The dotted construction lines fall out as small
compact components and are masked off before any marker is measured. **On a figure with
crossing curves, try component labelling before reaching for an ROI.**

**The figure carries its own attribution check, and on a 2-D map that is the check that
matters.** Six dotted lines each join the k-th marker of all eleven speed lines (a seventh,
vertical, joins the left ends at Ps3/P2 ≈ 1). They are printed ink. For every k and every
adjacent pair of speed lines, sample the straight link between the two extracted markers
and take the median distance to the nearest printed dot: **all 60 links score 6.1–9.8 px**
against a ~7 px dot half-spacing, while deliberately wrong pairings (k of one line to k+1
of the next) score 16–23 px. That is what licenses saying every point is on the right
speed line — not the fact that the extraction produced 6 markers per curve.

**The marker is the curve's own index digit, so one estimator does not fit all.** `1`–`11`,
which is one glyph on nine curves and a two-glyph *string* on two, and a digit's ink
centroid is not its string centre — the density peak of `10` lands on the `0`, 21 px right
of the data point. Per curve: density peak (window 33, min-sep 27, min-fill 0.30) for the
single digits; the enclosed bowl of the `0` for symbol 10, offset by the −21.9 px measured
at that curve's own left endpoint; the pair of `1` strokes from a 33×1 vertical opening for
symbol 11; vertical strokes alone for symbol 1.

**A bounding-box refinement was built, validated, and then rejected — worth recording as a
failure mode.** Isolating "columns whose ink extent exceeds the joining line's" and taking
the bbox centre is unbiased on the eleven left-endpoint glyphs (sd 3.0 px against a known
truth) and looks like a strict improvement. It is not: where two markers sit within a glyph
width it merges or clips their column runs and moves them up to 10 px, turning the regular
83/67/66/49/32 px marker spacing shared by curves 7–11 into 85/58/66/49/40. **An estimator
validated only where the data is easy can be worse than the seed where the data is hard.**
The link check caught it; the CSV would not have.

**Two internal ground truths this figure hands you free, both used to measure and neither
to adjust.** (a) The eleven flat-extension left ends sit on one printed vertical dotted
line, x = 1075.7 px ± 2.2 over 97 dots — so the estimator's x error is directly measurable
(+4.3 px bias, sd 3.0). (b) Each knee is the right end of its own flat segment, so its `y`
*must* be that segment's level; the segments fit to rms 0.4–1.2 px and the knee glyphs land
−8.7…+8.4 px from them.

**600 dpi does not add detail on this scan and the inventory's advice to use it is only
half right.** `pdfimages -list` shows p.56 is a single **300 dpi 1-bit CCITT image**
(2544×3300), so a 600 dpi render is a 2× upsample. It is still the right render — the
interpolated grey edges give better sub-pixel centroids and stable morphology — but it
cannot resolve what the scan did not record. **Check `pdfimages -list` before assuming more
dpi buys resolution.**

**Curve 1 (65 % NGc) is the honest limit.** Its six sloped markers span 0.59 in Ps3/P2 with
a glyph ~0.12 wide, so they overprint into one ink mass; at native 300 dpi they are 7–13 px
apart with a ~10 px glyph. They are *not* unattributable — a vertical opening resolves six
distinct `1` strokes, and all six sit −4.3…+2.4 px from their construction lines, the
tightest residuals of the family. But the bbox refinement cannot run there and the marker
spacing is the one place the family pattern is not monotone, so those five points carry
≈ ±0.08 in Ps3/P2 against ±0.05 elsewhere. The CSV header says so per point.

**Calibration now costs more than reading the markers.** Three defensible calibrations —
all minor ticks + frame lines, major ticks + frame lines, frame lines alone — move the same
pixel by up to **0.026 in Ps3/P2 and 0.015 in WA2c**, which is larger than the estimator
scatter. The all-ticks fit (used) puts the left construction line at Ps3/P2 = 0.974; frame
lines alone put it at 0.9995. The tick lattice and the frame lines genuinely disagree: the
x = 10 major tick sits 3.8 px right of the frame midpoint, and the four y majors imply a
scale 0.33 % shorter than the frame span. Nothing on the page resolves it, so the CSV
states both and tells the consumer not to assume the 1.0.

**`axes` on this page returns `top=1655 bottom=5358`** — the eleven flat speed lines
averaged into the "top", the x-axis title text row taken as the "bottom". Measured frame is
`814 / 5010`. In 300 dpi terms that is a 420 px error against A3's 49 px — same bug,
nearly an order of magnitude worse, because eleven flat runs out-vote one frame line.


### Update — A2 redone, and `find_frame` fixed (2026-09-10)

`digitize.py find_frame` was rewritten after two agents independently found it broken. The
first version scored a row by **total ink** and averaged every strong row in each half of
the page. It was wrong nearly everywhere, and always invisibly:

| Page | What it returned | Measured | Error |
|---|---|---|---|
| A1, p.56 | top=1655 (mean of eleven flat speed lines) | 814 | 420 px at 300 dpi |
| A3, p.58 | top=362 (pulled by the flat plateaux) | 313 | 49 px = 0.009 in B1 |
| A7, p.62 | bottom = the x-axis **title** | — | text out-inks a thin line |

**What actually discriminates a frame line** is neither total ink nor longest run — ticks
break the lines, so on A3 the top frame's longest run is 994 px and the bottom's only 351.
It is the **fraction of the plot width that is inked on that row**. The frame spans it;
text and data do not. The vertical frame lines are solid enough for a longest-run scan to
find directly, and they supply the width to measure against.

**The new `find_frame` refuses rather than guesses.** On a rotated page no row reaches the
threshold, because a tilted line has no single row — so it raises and names
`rectify_page.py` instead of returning a plausible wrong answer. Verified: it returns
(467, 313, 2134, 2672) on the rectified A3 against a measured 313/2672, and refuses on the
un-rectified A2.

**A2 was redone** on a rectified page with the corrected frame. The values moved only
0.12 % of range in x and 0.17 % in y — A2 has no flat run, which is the bug's trigger, so
it was largely spared. The free integer check tightened from 0.0457 to 0.0300 maximum
deviation, which is the real confirmation.

**A fifth gotcha, found during the redo:** `--min-sep` must exceed its default of 14 where
a marker sits on a frame line. A2's marker at PS3/P2 = 20 splits into two density peaks
otherwise, giving 21 points with a spurious 19.92 beside the real 20.00. The count is
stable at 20 over `min-sep` 25–60. **A count that will not plateau usually means a split
marker, not a missing one** — check the extracted values for a near-duplicate before
lowering the threshold.

### Update — A7 re-digitized (2026-09-11)

`f7` is the accuracy bottleneck of the engine model — +1 % on it moves shaft power by
+3.1 % (hover), +2.7 % (level), +5.3 % (descent), and all three published trim conditions
sit at P45/P41 = 0.2151/0.2188/0.2271, inside the crowded knee. So A7 was redone from
scratch. New extractor: **`tools/digitize_a7.py`**, which reruns and reproduces the six
rows of `data/maps/f7_gg_turbine_energy.csv` exactly.

| | old (2026-09-10) | new | move |
|---|---|---|---|
| 1 | 0.200045, 48.2414 | 0.199937, 48.2815 | −1.1e-4, +0.040 (+0.08 %) |
| 2 | 0.215072, 40.3295 | 0.214923, 40.3985 | −1.5e-4, +0.069 (+0.17 %) |
| 3 | 0.220051, 39.0561 | 0.219917, 39.1268 | −1.3e-4, +0.071 (+0.18 %) |
| 4 | 0.225030, 38.2906 | 0.224969, 38.3672 | −0.6e-4, +0.077 (+0.20 %) |
| 5 | 0.230100, 37.5252 | 0.229970, 37.5898 | −1.3e-4, +0.065 (+0.17 %) |
| 6 | 0.350045, 24.5962 | 0.350096, 24.6064 | +0.5e-4, +0.010 (+0.04 %) |

Per-point uncertainty, 1σ: **x ±0.0001, y ±0.03 BTU/LBM**. So the x moves are 0.5–1.5σ
(inside) and the four knee y moves are 2.2–2.6σ (**outside**). Marker 1 and marker 6 moved
1.3σ and 0.3σ in y — inside.

**The old marker detection was already right; the old *calibration* was not.** This is the
finding, and it is the one that generalises. Reversing the old header's stated calibration
recovers old pixel centres that agree with the new ones to about 0.5 px, and applying the
old y map to the new centres reproduces the old CSV to 0.007 BTU/LBM. Every bit of the
0.07 BTU/LBM move is the y axis, not the ink.

**A sixth failure mode: a linear axis map cannot fit these scans, and the tick lattice will
tell you so.** The page bows along the scan direction — A7's y major-tick spacing runs
181 → 177 → 180.5 px from top to bottom, a smooth 2 % swing, worth ±4 px against a straight
line. The A3–A11 note above called this "the accuracy floor of a two-point linear map";
it is not a floor, it is a **fittable, and checkable, quadratic**:

* linear over the 24 printed y majors → top frame = **52.5624** against a **printed 52.5**;
* quadratic over the same 24, never told about 52.5 → top frame = **52.4947**;
* quadratic on the left-edge ticks alone predicts the measured top frame row to **0.22 px**,
  on the right-edge ticks alone to **0.48 px**.

Residual falls from 0.033 to 0.019 BTU/LBM rms. **Where a figure prints a value at a frame
edge that the tick fit is not told, that value is a free test of the calibration model.**
Every other axis on every other figure in this appendix should be re-checked this way; A7's
0.18 % was invisible because the linear fit's residuals looked like noise and were not —
they were an arch, negative at both ends and positive in the middle.

**Don't be fooled by the printed label at the frame.** A7's `52.5` label sits 3.8 px *below*
the top frame line while every other y label sits within 0.9 px of its own tick. That
offset is what makes the linear fit look defensible. The lattice, not the label, fixes the
anchor — same lesson as the displaced `0.11` on A3, from the other direction.

**Do not resample the page.** `pdfimages -list` shows p.62 is one 300 dpi 1-bit CCITT
image. `rectify_page.py` warps pixels, which needs interpolation and (on this page) a
600 dpi render just to find the frame. Instead: pull the scan out with `pdfimages`, fit the
four frame lines sub-pixel *in the scan*, and carry the rotation and keystone as a
**homography applied to the measured points**. Nothing interpolates the data, and the four
edge slopes come out −0.552/−0.649/+0.766/+1.018°, i.e. a genuine quadrilateral, not one
rotation. A full independent re-run on a 900 dpi `pdftoppm` render agrees on x to 2e-5 and
on the marker rows to 0.13 px; its y is up to 0.06 BTU/LBM different only because its
greyscale-threshold frame and ticks are worse — its unconstrained lattice misses the
printed 52.5 by 0.09 where the native bitmap misses it by 0.005. **300 dpi native wins, and
"go to 600 where markers crowd" is the wrong instinct on a 300 dpi scan.**

**What replaced density mode: a generative model of the local ink.** Model the window as
`stroke1 ∪ stroke2 ∪ incoming half-segment ∪ outgoing half-segment`, all four radiating
from one point; render it by 5×5 supersampling; score `Σ(coverage − ink)²` over an unmasked
disc. The glyph template is global — the plotter drew the same mark six times — and comes
out as strokes at −44.05° (half-length 14.88 px) and +46.09° (15.22 px), 3.46 px wide.
This scores **background as well as ink**, which is exactly why the wedge cannot win: two
near-parallel lines 16 px from the frame do not have four arms. Frame lines and ticks are
**masked, never modelled** (m1's mask also covers a minor tick hidden underneath it at
y = 608.6; m6's covers one at y ≈ 2323). Independent check — mask the joining line out
altogether and fit the glyph alone: agrees to **≤0.62 px** at all six.

**Rejected: reading the interior points as polyline vertices.** The curve *is* five straight
segments (fitted rms 1.0–1.15 px; the 1624 px segment 5 is straight to a 0.18 px sagitta),
so each interior data point is also a kink, and intersecting the fitted segments should give
it. It does not give it *accurately*: the kinks are 15.2°, 14.0°, **1.1°** and 10.5° on
segments only 78–107 px long, so the m4 intersection is indeterminate to tens of pixels and
lands 5.7 px from the truth. Kept only as a consistency check — the fitted centres sit
within 1.2 px of the independently fitted adjacent segment lines, which is what licenses
saying the vertex is the data point.

**Two free checks this figure hands you.** (a) The six abscissae land on a 0.005 lattice —
0.200, 0.215, 0.220, 0.225, 0.230, 0.350 — to a maximum of **1.1 px** (1.0e-4), sd 6e-5,
with nothing in the extraction told to expect it. That is where the x uncertainty number
comes from. (b) **A marker count that does not depend on a detector:** walk the ink density
along each fitted segment in 10 px bins; every excess over the joining-line baseline sits at
a segment endpoint and there is no interior excess anywhere, including on the 593 px segment
1 and the 1624 px segment 5. Six markers, no seventh.

**A tick filter that returns the right count and the wrong ticks.** m1 and m6 sit *on* a
vertical frame line and out-protrude a major tick, so the protrusion scan reports a false
major at each. A greedy "drop anything closer than 0.6 of the median spacing to the last
keeper" rule keeps the *marker* at y = 616 and drops the true 47.5 tick at y = 664 — and
still returns 12 majors, so the assert passes and the y calibration is quietly wrong by
0.06 BTU/LBM. The fix in `digitize_a7.py` is to keep the largest subset lying on **one
arithmetic lattice**. **Assert on structure, not on count.**

**Does this close the 0.70–2.30 % shaft-power gap? Only part of it.** At the three trim
conditions `f7` rises by 0.077 %, 0.089 % and 0.153 %, which through the measured
sensitivities is **+0.24 %, +0.24 % and +0.81 % of SHP**. What remains of `f7`'s
uncertainty is ±0.26 %, ±0.23 % and ±0.38 % of SHP at 1σ — so `f7` could plausibly move
SHP by at most another ~0.5–0.8 % even at 2σ, and only at the descent condition does that
reach the top of the observed gap. **A7 is now measured to the limit of the ink; a residual
of order 1 % at the worst condition has to come from somewhere else** (`f9` is the other
high-leverage map, and the calibration-model finding above applies to every figure in the
appendix).

### Update — A1 recalibrated (2026-09-11)

`f1` is the largest map in the engine (11 speed lines × 7 points, the only 2-D one) and
everything downstream of the compressor rides on it, so A7's calibration finding was
applied to it. **The extraction survived; the calibration did not.** No marker pixel moved.
`tools/digitize_a1.py` still reruns and reproduces `data/maps/f1_compressor_mass_flow.csv`
byte for byte.

**The held-out frame-edge test, run on both axes.** A1 is unusually generous here: *both*
axes print their end values *at* the frame lines — 0 and 20 in x, 2 and 12 in y — so a fit
made on the interior ticks alone (19 minor ticks per edge in x, 49 per edge in y, all four
edges measured on the native scan) has **four** numbers left to predict.

| axis | model | tick-fit rms | predicts the two printed frame values to | |
|---|---|---|---|---|
| x | linear | 0.80 px | −2.38 px and −0.95 px | **fail** |
| x | quadratic | 0.63 px | −1.09 px and +0.32 px | pass |
| x | **cubic** | 0.46 px | +0.48 px and −1.24 px | **used** |
| x | quintic | 0.38 px | −1.26 px and +0.82 px | pass |
| y | linear | 1.94 px | −4.00 px and +3.71 px | **fail** |
| y | quadratic | 1.94 px | −4.36 px and +3.35 px | **fail** |
| y | **cubic** | 0.91 px | +0.79 px and −1.79 px | **used** |
| y | quintic | 0.74 px | −1.05 px and +1.10 px | pass |

**On y the quadratic is worthless — identical to the line — because this page's bow is odd,
not even.** A7's was a symmetric arch and a quadratic caught it. A1's y residual to a line
runs −2.7 px at the top, +3.1 px at *v* ≈ 0.68, −3.0 px at *v* ≈ 0.40 and +3.1 px at the
bottom: two interior turning points, so degree 3 is the first that can fit it, and it halves
the tick residual when it does. The minor-tick spacing is 42.4 px near the top frame, 41.2
in the middle, 42.1 near the bottom. **"Try a quadratic" is not the lesson from A7; "fit the
lowest order that passes the held-out test" is, and on a different page that is a different
order.** Degrees 1 and 2 came out indistinguishable, as did 3 and 4, and 5 and 6 — the
signature of an odd distortion, and a cheap way to read the parity off the fit table.

**The frame-vs-tick disagreement was ours, and it is resolved.** The 2026-09-10 header said
"the tick lattice and the frame lines genuinely disagree … nothing on the page resolves it",
on the evidence that a linear tick fit implies a y scale 0.33 % short of the frame span. That
reproduces on the un-warped scan (0.37 %) — and it is entirely an artefact of forcing a
straight line through a bowed lattice. **The model-free version of the test settles it: fit
only the five ticks nearest each frame line and extrapolate.** The lattice then predicts all
four printed frame values to **1.41 px or better** (x: −0.06, −1.41, −0.07, −0.38; y: +0.61,
+0.55, −0.64, +1.16). Locally the ticks and the frame agree everywhere; only a global straight
line makes them fight. Neither old calibration was "right": the frame lines were right, and
the all-ticks fit was right about the data and wrong about the model. **A local extrapolation
from the nearest few ticks is the cheapest possible diagnostic and should be run before
anyone writes "the page disagrees with itself".** (The other half of the old complaint — the
x = 10 major tick sitting 1.9 px right of the frame midpoint at 300 dpi — is *real*, and
reproduces exactly in the native scan. That one was a correct measurement, not a warp
artefact; it is simply small, and the cubic absorbs it.)

**A third, unprinted check the figure hands you, now passed better.** The eleven flat
extensions end on one printed vertical dotted construction line, which on a compressor map
ought to be Ps3/P2 = 1. Old linear calibration: **0.974**. Frame-lines-only: 0.998. New
cubic: **0.9907**. The report never prints that 1.0, so it cannot be used as a datum — but a
model that was told nothing about it and moves most of the way there is telling you
something.

**Do not warp pixels — but if you already have, measure what it cost instead of redoing it.**
The markers are still read on the rectified 600 dpi render, because the 2× upsample's
interpolated grey edges are the only thing that resolves SYMBOL 1's six overprinting markers
(7–13 px apart under a ~10 px glyph at native 300 dpi; a morphological opening cannot split
them). The geometry is now measured on the native scan and the two are joined by (u,v), the
unit square on the frame corners. That bridge is **measured, not assumed**: twelve structures
visible in both rasters — the printed vertical construction line and the eleven flat
extensions — were located independently in each and agree to **0.25 px** (rms 0.12). Running
the density-peak estimator directly on the native bitmap likewise reproduces the 600 dpi
marker centres to ≤1.0 px in x and ≤0.7 px in y. **The warp was worth a quarter pixel; the
linear axis map was worth four.** Frame slopes on this page are −0.011°, +0.050°, +0.150°,
+0.428° — a genuine quadrilateral, and it is carried as a homography on the points.

**What moved.** All of it is calibration.

| | move | in σ |
|---|---|---|
| `ps3_p2`, the eleven k=0 rows | +0.017 (0.974 → 0.991) | 2.4 |
| `ps3_p2`, the 66 sloped rows | −0.011 … +0.007 | ≤0.35 |
| `wa2c`, the 22 k=0/k=1 rows | −0.0098 … +0.0104 | ≤2.8 |
| `wa2c`, the 55 glyph rows | −0.0098 … +0.0104 | ≤0.89 |

In percent, `wa2c` moved −0.215 % … +0.115 %. **The sign follows the bow, not the data:**
rows above WA2c ≈ 7.0 moved up by 0.003–0.010, rows between 4.5 and 6.8 moved down by
0.001–0.010, and SYMBOL 1 near 3.0 barely moved at all. Nobody would have found that pattern
by staring at the extracted numbers — it is only visible once the tick lattice is fitted as
a curve. The largest σ-moves land on the *most precisely located* quantities (the flat-segment
levels and the construction line), which is the signature of a calibration error rather than
a reading error, exactly as on A7.

**Per-point 1σ, three measured terms in quadrature:**

| rows | 1σ `ps3_p2` | 1σ `wa2c` |
|---|---|---|
| k=0, all speed lines | 0.007 | 0.0037 |
| k=1, SYMBOL 2–11 | 0.031 | 0.0037 |
| k=2…6, SYMBOL 2–11 | 0.031 | 0.0117 |
| k=1…6, SYMBOL 1 | 0.079 | 0.0037 / 0.0117 |

The calibration term is now the *spread of the models that passed the test* (x: anchored
quadratic, quartic, quintic; y: quartic and quintic — the failed y quadratic does not count),
evaluated at the 77 points: 0.0061 in Ps3/P2 and 0.0032 in WA2c. **That is a third of what
the old header had to quote** (0.026 and 0.015 across three "defensible" calibrations),
because the held-out test disqualifies two of the three.

**The attribution check is untouched and was re-run: all 60 links score 6.1–9.8 px, none
fail.** Marker pixels did not move, so this could not have changed — but a recalibration is
exactly the moment someone quietly breaks it, so it is re-run every time.

**A generative glyph model was considered and NOT built, and the figure says why.** A7's
marker is one two-stroke cross drawn six times, so a single template serves and the model
beats density mode outright. A1's markers are eleven different digit strings, 12–18 px wide
at native 300 dpi, with the joining line running through them — eleven templates, fitting
shape to three or four ink pixels per stroke. Before building it, the estimator was tested
against printed ink: for each sloped marker, take the construction-line dots 15–55 px away
on both sides (**the construction lines are curved, so they must be fitted locally — a global
straight or polynomial fit misses them by 7–20 px and is worthless**), fit total least
squares with the marker excluded, and measure the perpendicular offset. 47 of the 66 have
enough dots: mean **+0.19 px**, sd 2.54 px, and solving all 47 for a *common* x displacement
of the whole family returns **−0.09 px**. There was no bias for a glyph model to remove.
**Test the estimator you have against printed ink before writing the one you imagine.**

**One measured offset that is real and deliberately not corrected.** The density peak on the
eight isolated single-digit *left-endpoint* glyphs sits +2.42 px (sd 1.44) right of the
printed vertical construction line. The obvious explanation — one-sided line ink, since the
flat extension leaves rightward only — is **falsified**: deleting that line from the window
changes the answer by 0.00 px. Whatever it is (glyph anchor convention, or the dotted line
drawn fractionally left), it does not appear on the sloped markers, where the independent
test above finds −0.09 px. So it is recorded and not applied. **A bias measured at one
special point is not a bias until you find it where the data actually is.**

### Update — A6, A8, A10 recalibrated (2026-09-11)

The three functions on the power-turbine and combustor path, redone from scratch with A7's
method. New extractor **`tools/digitize_a6810.py`**, which reruns and reproduces all three
CSVs — `f6_combustor_efficiency.csv`, `f8_pt_energy.csv`,
`f10_exhaust_pressure_loss.csv` — byte for byte, headers included. Overlays in
`validation/out/digitize/a6810/`.

Why these three: in a sweep of the trimmed engine `f8` and `f10` move shaft power **without
moving gas generator speed at all** (ΔNG = 0.000 to five figures), so the 0.08–0.11 % speed
agreement cannot constrain them; `f6` adds +1.44 / +1.79 / +3.01 % of SHP per 1 % at the
three Table B.1 trims.

| | order chosen (x, y) | markers | largest x move | largest y move |
|---|---|---|---|---|
| A6, p.61 | linear, **quartic** | 2 (inventory 2 ✓) | +1.5e-5 = 1.2σ | −3.8e-5 = **0.1σ** |
| A8, p.63 | linear, **quartic** | 12 (12 ✓) | −6.6e-4 = 2.0σ | +0.131 BTU/LBM = **6.3σ** |
| A10, p.65 | linear, **quartic** | 18 (≈20 ✗) | +0.106 % = 1.6σ | −3.2e-4 = **4.1σ** |

(Largest in σ; the largest in absolute terms is a different row in two cases — A8 x −1.05e-3
at 1.8σ, A10 y −3.9e-4 at 2.4σ.)

**The held-out frame test is decisive on all six axes, and on none of the three y axes does
it pick A7's quadratic.** Fit the interior majors alone, then predict the printed frame value
the fit was never given. Errors in per-mille of the axis range (bold = chosen; on A6 x the
quadratic scores marginally better but is inside the 1.25× parsimony tie-break, and its
leave-one-out is worse, so the linear stands; on A10 x the quartic scores best and is
rejected — see below):

| axis | linear | quadratic | cubic | quartic | used |
|---|---|---|---|---|---|
| A6 x (0.010 and 0.020 printed) | **0.44** | 0.39 | 1.44 | 1.76 | linear |
| A6 y (0.88 and 1.10 printed) | 1.79 | 3.18 | 2.73 | **0.50** | quartic |
| A8 x (0.85 printed) | **0.26** | 0.53 | 0.55 | 0.65 | linear |
| A8 y (−5.0 and 40.0 printed) | 2.51 | 2.57 | 0.88 | **0.65** | quartic |
| A10 x (65 and 100 printed) | 1.94 | 2.19 | 2.60 | *1.11* | linear |
| A10 y (1.14 printed) | 1.81 | 0.35 | 1.70 | **0.36** | quartic |

**A6's quadratic is worse than its linear.** That kills "A7 found a quadratic, so try a
quadratic" as a rule outright — the A1 note above reached the same conclusion from the other
direction. The parity trick that note describes reads straight off the leave-one-tick-out
column: A8's y pairs as (1,2)(3,4), gain at **odd** order 3 → an S-shaped distortion; A6's y
pairs as (2,3) with separate gains at 2 and 4 → an **even** distortion with a real quartic
term; A10's y improves steadily with no pairing at all.

**The bow is per-page, not a scanner signature — do not carry a model between figures.**
Deviation of the y tick lattice from a straight line through it, per-mille of the axis
range, sampled at *v* = 0.05 … 0.95:

| | 0.05 | 0.20 | 0.35 | 0.50 | 0.65 | 0.80 | 0.95 |
|---|---|---|---|---|---|---|---|
| A6 | −1.19 | −1.06 | +0.85 | +1.82 | +0.76 | −0.52 | −2.08 |
| A7 | −0.10 | −1.00 | +0.30 | +0.93 | +1.10 | −0.24 | −1.73 |
| A8 | +0.55 | −0.96 | −0.48 | +0.44 | +1.68 | −0.20 | −1.64 |
| A10 | −0.06 | +0.66 | +0.55 | −0.90 | −1.12 | +0.40 | +0.84 |

A6 and A7 share an arch; A8 is that arch with an extra wiggle; **A10 is inverted**. Same
amplitude everywhere, ±1–2 per-mille ≈ ±3–5 px, always far above the 0.2–0.5 per-mille
tick-fit noise. Four pages, four shapes.

**A seventh failure mode: the held-out frame test alone will over-fit, and leave-one-out is
the guard.** On A10's x axis the quartic wins the frame test (1.11 per-mille against the
linear's 1.94) while its leave-one-tick-out error is **40 % worse** (1.94 against 1.38) —
because that axis has only 6 interior majors per edge, no curvature at all in its residuals,
and a genuine ~2 px offset *between* the two edges that no polynomial in one variable can
absorb. The quartic was buying the frame with noise. Rule used here and worth copying:
**choose the order that minimises the held-out frame error, but only among orders whose
leave-one-tick-out error is within 1.2× of the best; ties within 1.25× go to the lower
order.** Both criteria are honest — neither is a fit statistic on the data being fitted.

**A new free check, and it fixes two unlabelled frames.** Measure how many major steps the
frame spans, by the ticks alone. On the four axes here where *both* frame values are printed
it lands within 0.08 of an integer every time — A6 x 10.008, A6 y 11.015, A8 y 18.071,
A10 x 7.008 — so on the two where one frame is unlabelled the integer fixes that frame:

* A8 x spans **11.017** steps of 0.05 → the unlabelled left frame is **0.30 exactly**;
* A10 y spans **12.974** steps of 0.01 → the unlabelled bottom frame is **1.01 exactly**,
  which retires the inventory's "≈ 1.012" and the old header's 1.010233.

(Read the span as a rounding test, not as a precision measurement: it comes from a straight
line through the lattice, so on a bowed axis the bow's odd component leaks into it. The
margin to the wrong integer is ≥ 0.9 steps everywhere here, so the rounding is never close.)

**And it explains why extrapolating the ticks to an unlabelled frame lands short.** On the
three x axes — the three that turn out to be genuinely linear, so the span means what it
says — the printed tick lattice is **0.08–0.16 % tighter** than (frame span)/N: A6 x
+0.076 %, A10 x +0.114 %, A8 x +0.155 %. A tick-only fit extrapolated one step past the last
tick therefore undershoots by ~2 px. On A8 that is the difference between 0.29932 and
0.30000, and it was enough to tilt all twelve abscissae: before anchoring the left frame the
0.05-lattice check came out max 6.3e-4 with a visible slope, after it max **5.0e-4, sd
3.1e-4**, against a mean σ_x of 3.8e-4. **Anchor on the frame; never extrapolate the lattice
to it.**

**Marker seeding by morphological opening, not density.** The joining line is ~4 px wide and
a 7×7 opening erases it and every tick; where two strokes cross, ~5 px of ink survives in
every direction. One core per marker, no threshold sweep, nothing to tune, and it is a
*shape* discriminator like the generative fit rather than a density one. It returned 12, 19
and 2 cores inside the frame with no help. **A10's ~70 % marker is the one interesting
case:** at 7×7 it splits into two cores 9 px apart, at 6×6 they are a single 138 px core —
the scan broke the glyph's waist. Merging cores closer than 12 px gives 18, and the
independent interior-ink walk finds no excess on any of the 17 joining segments. **A10 has
18 markers, not the inventory's "≈ 20"** — same conclusion as 2026-09-10, now from a
detector that shares nothing with the one that reached it then.

**An eighth failure mode, and it is a SciPy default. `minimize(..., method="Nelder-Mead")`
builds its initial simplex by perturbing each coordinate by 5 % of its own value.** At an
image column of x ≈ 2269 that is a **113 px** step — far outside a ±25 px scoring window,
where the cost is flat because no model ink reaches it. On A6's right-hand marker, whose
disc is mostly masked (it sits on the frame line, among minor ticks every 21 px), the fit
took that step and never came back: the reported centre was **125 px** off the page's data
line, the cost was low and constant, and every printed diagnostic looked converged. Fix:
optimise the **offset** from the current centre with an explicit 2 px initial simplex, and
assert that no marker moves more than 10 px from its seed. **Never hand Nelder-Mead an
absolute pixel coordinate.**

**Mask the frame symmetrically or the abscissa is biased.** Ticks protrude inward only. A
marker sitting on a frame line has the inward half of its disc stripped by the tick mask and
the outward half left intact, which pulls the fitted centre outward. All three figures print
a marker on a vertical frame (A6 both, A8 first and last, A10 first and last), so the mask
here covers |d| < k+5 on *both* sides of each frame line. With that in place the six
frame-sitting markers land 0.06–1.96 px from their frame lines, and A6's two abscissae come
out 0.009990 and 0.020003 against printed frame values of 0.010 and 0.020 — a free check the
extraction was never told about, passed at **1.0σ and 0.25σ**.

**A local two-tick interpolation is bow-free, costs nothing, and should be run on every
figure.** Interpolate linearly between the two majors that bracket the point, on each ticked
edge, and average. It cannot see a global bow because it never fits one. Chosen model minus
two-tick, worst point: A6 **5.9e-6** in y, A8 0.0225 BTU/LBM, A10 8.7e-5. Where they agree,
the polynomial is not wandering between ticks.

**A6 was already right, and that is the most useful single result here.** Its 2026-09-10
value came from calibrating against the two ticks that bracket the line (1.00 and 0.98)
rather than a global linear fit — i.e. the bow-free estimator above, chosen for the right
reason at the time. New η = **0.985038** and **0.984967** against old 0.985076 and 0.984982:
**0.1σ and 0.04σ**. The whole function is one number, 0.98500 ± 0.00041, and the two markers
agree with each other to 7.1e-5 (0.1σ), which is the figure's own statement that it is a
constant.

**A8's move is calibration, not detection — except at the marker on the frame.** Push the
new pixel centres through the old two-point linear y map and the old CSV comes back to
0.0657 worst, 0.0169 rms. Per point the decomposition is: calibration +0.097…−0.066
BTU/LBM, detection +0.004…+0.015 everywhere *except* m1 at P49/P45 = 0.30, where detection
alone contributes **+0.066** — that is the marker sitting on the left frame, which the old
density detector read 3.4 px low. The linear map's own error explains the rest: it puts the
top frame at 40.113 against a printed 40.0, i.e. 5.9 px, and m1 is the point nearest the top.

**A10's y moved 1–4σ downward at every one of its 18 points, and that is the label-centroid
bias biting.** The old y scale came from a 7-point least-squares fit over the printed tick
*labels*; the A7–A10 note above already measured that label centroids carry a 4–5 px
systematic. Fitting the tick *marks* instead moves every value down by 1.3e-4…3.9e-4. The x
values moved +0.106 % at worst (1.6σ), which is the 6-major-per-edge x calibration, the
weakest of the six axes here (σ_x = 0.047…0.065 % NGc, against σ_y = 7.8e-5…1.65e-4).

**Per-point 1σ, quadrature of four measured terms** — tick-fit residual, spread over
calibration orders 1–4, half the disagreement between the two ticked edges calibrated
separately, and the glyph centre (worst case over window radius ±3 px, frame-mask width 3/5
px, tick-mask width 3/6 px, and masking the joining line out entirely and fitting the two
strokes alone; that worst case is **0.18 px on A6, 0.35 on A10, 0.45 on A8**, and on all
three it is the line-masked-out variant that sets it). The CSVs carry it per row as
`sigma_x` and `sigma_y`:

| | 1σ x | 1σ y |
|---|---|---|
| A6 | 1.0e-5 … 1.2e-5 FAR | 4.07e-4 η |
| A8 | 3.2e-4 … 5.8e-4 | 0.021 … 0.057 BTU/LBM |
| A10 | 0.047 … 0.065 % NGc | 7.8e-5 … 1.65e-4 |

The order term is deliberately conservative — it includes the linear map that the held-out
test rejects. Over the orders that pass, it is 2–3× smaller.

**What it bought, and what it cannot buy.** Re-trimming at the three Table B.1 conditions
with the old maps injected and then the new ones, everything else held:

| | NG, before → after | SHP, before → after | ΔSHP |
|---|---|---|---|
| hover | −0.109 % → −0.109 % | −0.643 % → **−0.406 %** | +0.239 |
| 80 kt level | +0.082 % → +0.082 % | −0.466 % → **−0.230 %** | +0.237 |
| 80 kt descent | −0.081 % → −0.082 % | −1.319 % → **−1.120 %** | +0.201 |

NG does not move, as the sensitivity sweep said it would not. Almost all of the gain is `f8`,
which rose +0.21…+0.23 % at the trim arguments (P49/P45 = 0.404 / 0.491 / 0.600); `f10` fell
0.03 %…0.00 % and `f6` fell 0.002 %, worth +0.03 % and −0.005 % of SHP.

**Remaining headroom at 1σ, shifting every knot of a map coherently:**

| | hover | 80 kt level | 80 kt descent |
|---|---|---|---|
| `f6` | ±0.059 % | ±0.074 % | ±0.125 % |
| `f8` (x and y in quadrature) | ±0.227 % | ±0.278 % | ±0.201 % |
| `f10` | ±0.021 % | ±0.021 % | ±0.078 % |
| **all three** | **±0.235 %** | **±0.288 %** | **±0.249 %** |

So the residual SHP deficit is **inside 1σ at the 80 kt level trim, 1.7σ at hover, and 4.5σ
at descent**. These three maps are finished: even at 2σ they cannot supply the descent
condition's remaining 1.1 %. `f7` (+5.45 % of SHP per 1 % at descent, and +2.90 % at hover) and
`f9` (+2.08 % at descent, **+7.28 % at hover** — the largest single lever in the engine) are
where the rest has to come from. `f7` is already re-digitized; `f9` was in progress in
`tools/digitize_a9.py` while this was written, so read its own update section rather than
this sentence.

### Update — A2, A9 recalibrated (2026-09-11)

`f2` and `f9` are the two highest-leverage maps left in the engine. A sensitivity sweep of
the trimmed model, measured here by perturbing each CSV and re-solving, gives per **+1 %**
on the map:

| map | ΔSHP hover | ΔSHP level | ΔSHP descent |
|---|---|---|---|
| `f9` | +7.50 % | +2.03 % | +2.01 % |
| `f2` | −2.09 % | −3.69 % | −8.05 % |

So a 0.2 % calibration error in either is worth up to 1.5 % of shaft power. Both were
redone from the native scan with A7's method plus two upgrades that A7's own numbers
implied but did not use. New extractors: **`tools/digitize_a2.py`** and
**`tools/digitize_a9.py`** over the shared **`tools/digitize_native.py`**; each reruns and
reproduces its CSV rows byte for byte.

| | old (2026-09-10) | new | move | 1σ |
|---|---|---|---|---|
| A2, worst y | 2.805290 at x=20 | 2.808680 | **+0.00339 (+0.121 %)** | ±0.00055 |
| A2, worst x | 13.976 at k=14 | 14.00065 | +0.0247 (+2.5σ) | ±0.0097 |
| A9, worst y | 0.285285 at x=0.80 | 0.284948 | **−0.000337 (−0.118 %)** | ±0.000027 |
| A9, worst x | 0.42598 at k=6 | 0.42523 | −0.00075 (−4.7σ) | ±0.00016 |

Every `f2` ordinate moved, by −0.0028 at the low end through +0.0034 at the high end
(−0.22 % … +0.12 %), which is **2.5–6.2σ** at fourteen of the twenty points. Every `f9`
ordinate moved too: **+0.00024 at x ≤ 0.60, −0.00034 at x ≥ 0.75**, a rotation of the whole
curve, 3–12σ at nineteen of twenty-three. The direction reverses at x ≈ 0.70. Neither move
is the linear-vs-polynomial correction on its own (that is at most +0.0018 on `f2` and
+0.00016 on `f9`); most of it is the old rectify-then-calibrate pipeline's frame location.

**Net effect on the trims:** SHP against Table B.1 goes −0.41 % → **−0.26 %** (hover),
−0.23 % → −0.18 % (level), −1.12 % → **−0.91 %** (descent). `f9` supplies about two thirds
of that and `f2` the rest.

#### The seventh failure mode: the frame edges are arcs, so a homography is wrong too

A7 fixed the *axis map*. It left the *geometry* as a homography between four corners —
i.e. four straight edges. **They are not straight.** Fit each printed frame edge as a
quadratic in the coordinate along it and read off the sagitta:

| page | left | right | top | bottom |
|---|---|---|---|---|
| A2, p.57 | **−1.16 px** | **−1.42 px** | −0.98 px | +0.52 px |
| A9, p.64 | −0.51 px | −0.88 px | −0.72 px | +0.43 px |
| A7, p.62 | −0.44 px | −0.75 px | −0.38 px | −0.14 px |

These are *printed straight lines*, so the sagitta measures the page distortion directly,
with no reliance on any tick. Replace the homography with a curvilinear blend between the
four fitted arcs (`digitize_native.make_sv`: `s` is 0 on the left arc and 1 on the right,
`v` 0 on the bottom and 1 on the top, linear in between) and A2's free integer check —
its 20 abscissae against 1…20, which nothing was told to expect — tightens from **sd
0.0115 (0.95 px) to sd 0.0060 (0.50 px)**. Under the straight-edge geometry the residual
carries a −0.02 plateau across markers 9–14, sitting exactly where mid-page is; the arcs
remove it. On A7's flatter page the correction is worth little, which is why it did not
show up there.

#### The bow is odd, not even: on these two pages the y axis wants a CUBIC

A7's y major-tick spacing ran 181 → 177 → 180.5 px and a quadratic fixed it. A2 and A9 run
the same shape but symmetrically, and a symmetric spacing "smile" integrates to an
**S-shaped**, i.e. cubic, value function. Deviation of the printed y long ticks from a
straight lattice, in pixels, top to bottom:

* A2 left: `0 +1.4 +1.3 +2.2 +4.0 +3.9 +4.3 +3.8 +2.6 +1.0 −0.6 −1.3 −1.9 −2.4 −3.1 −2.7 −2.8 −1.4 0`
* A9 left: `0 +0.4 +3.2 +5.2 +5.6 +4.5 +2.9 +0.3 −1.2 −3.7 −4.4 −2.5 −1.6 −0.2 0`
* A9 right: `0 +1.1 +4.0 +5.3 +4.7 +4.0 +3.5 +1.2 +0.1 −3.0 −4.2 −2.3 −1.5 −0.2 0`

Left and right agree to about 1 px, so this is the page, not the detector. A quadratic
cannot represent an odd function and buys almost nothing: on A2's left edge the interior
residual goes 1.436 → 1.212 → **0.524** milli-units rms for orders 1, 2, 3; on A9's left
edge 15.1 → 15.0 → **6.2**e-5. **Order 2 is not "the" answer; the tick lattice tells you
which order, per page.**

#### The held-out frame test, done on four edges

| axis, edge | printed | linear predicts | quadratic | cubic |
|---|---|---|---|---|
| A2 y, left | 1.1 / 3.0 | 1.09809 / 3.00303 | 1.09605 / 3.00099 | **1.10012 / 2.99694** |
| A2 y, right | 1.1 / 3.0 | 1.09607 / 3.00221 | 1.09489 / 3.00104 | **1.09965 / 2.99633** |
| A2 x, bottom | 0.0 / 20.0 | 0.0205 / 20.0290 | 0.0101 / 20.0185 | **0.0335 / 19.9946** |
| A2 x, top | 0.0 / 20.0 | −0.0053 / 20.0180 | −0.0120 / 20.0113 | **0.0096 / 19.9897** |
| A9 y, left | 0.25 / 0.39 | 0.249846 / 0.390282 | 0.249816 / 0.390252 | **0.250194 / 0.389873** |
| A9 y, right | 0.25 / 0.39 | 0.249794 / 0.390196 | 0.249765 / 0.390168 | **0.250097 / 0.389846** |
| A9 x, bottom | — / 0.85 | 0.30038 / 0.850425 | 0.30028 / 0.850331 | **0.30076 / 0.849844** |
| A9 x, top | — / 0.85 | 0.29977 / 0.850355 | 0.30000 / 0.850582 | **0.30068 / 0.849901** |

On A9's x the test is decisive — only the cubic lands on the printed 0.85 to better than
2e-4 — and it hands back a **measurement of the unprinted left frame: 0.3008 / 0.3007**,
the first independent confirmation of the 0.30 that the 2026-09-10 header inferred from a
tick-count argument. On A2's x, linear and cubic both land within 2 px and the free integer
check is what settles it. **Where the frame test is a wash, look for a structural check;
where there is no structural check, the frame test is all you have — and one printed value
is not enough (see A7 below).**

A second test, which needs no printed edge value at all, is in `digitize_native.tail_cv`:
fit the middle 64 % of the tick lattice and predict the two ends. It agrees on A9's x
(order 3 halves the held-out rms on the top edge) and is too noisy to decide anything on
an 18-tick edge like A2's y. Report it; do not lean on it.

#### Calibrate on EVERY tick, and label the fine lattice through the coarse fit

A9 prints minor ticks at 0.001 in y (141 per edge) and 0.01 in x (56 per edge) — ten times
the data the long ticks give. They are usable, but only if they are labelled correctly, and
the straight line between the two frame corners is several pixels out mid-page, which is
the very defect being measured: it silently discarded a **contiguous block of 29 good minor
ticks** on A9's left edge. `label_by_map` labels the fine lattice through a cubic fitted to
the long ticks instead, and then all 141 survive. Long-tick and all-tick calibrations then
agree to 4.2e-4 in x and 1.2e-5 in y, and the difference is carried as a systematic.

A2's minor ticks are the counter-example: 0.01 spacing is only 12.3 px there and the
per-tick noise is ~3.7 px, so slot assignment fails outright (7 duplicate slots, 2 empty)
and a fit on them predicts the bottom frame 0.03 out. **Try the fine lattice; check that
every slot is filled exactly once; drop it if not.**

#### `label_ticks` replaces counting, and it caught a live error

`label_ticks` attaches a printed value to each tick from the two frame-edge values — a
*discrete* choice, so it cannot leak into the continuous fit that is later tested against
those same values — and drops anything more than 0.3 slots off the lattice. That is what
removes the marker arms A2's marker 20 leaves on the right frame (false ticks at y = 529.3
and 541.8, which also destroy the true 2.8 tick, so the right edge calibrates on 17 majors
and the left on 18) and the false long tick A9's marker 1 leaves at y = 612.6.

It also caught this, which is the reason to prefer it over any count-based filter: a run of
the A2 y detector returned **exactly 189 minor ticks, the exactly correct number**, while
containing one spurious entry and missing one real one. The off-by-one labelling moved the
axis by 0.03 in T3/T2 and every check based on the count passed. **A count that matches is
not a check.** (A7's note says "assert on structure, not on count"; this is that lesson
arriving from the other direction — the count was right and the ticks were wrong.)

#### Markers: same generative glyph model, and two independent checks it passes

Template `stroke1 ∪ stroke2 ∪ incoming half-segment ∪ outgoing half-segment`, 5×5
supersampled, scored `Σ(coverage − ink)²` on ink and background over a 20 px disc, with the
template global and only the centres per-point. A2 comes out as −45.00°/15.43 px and
+44.87°/15.28 px, both strokes 4.27 px wide; A9 as −44.96°/15.00 and +45.00°/15.00, glyph
3.60 px wide but the **joining line only 2.40 px** — on A9 the line is much thinner than the
glyph and must be a separate parameter.

* *Mask the joining polyline out of the score entirely and refit the glyph alone*: centres
  agree to **0.49 px max, 0.11 px rms** (A2) and **0.48 px / 0.11 px** (A9). A7's was
  ≤0.62 px.
* *Distance from each centre to the independently fitted joining segments on either side*
  (fitted from ink 20–24 px clear of every marker, so they know nothing about the glyph):
  A2 **mean +0.064 px, sd 0.220, max 0.52**; A9 **mean −0.023 px, sd 0.401, max 1.95**.
  A7's was 1.2 px. The near-zero *means* are the useful part: they exclude a common-mode
  glyph bias perpendicular to the curve at about 8σ, which is what lets the free-lattice
  offset below be attributed to the axis map rather than to the ink.

#### Free structural checks, and what their residual offset means

| figure | lattice | max | sd | mean |
|---|---|---|---|---|
| A2 | integers 1…20 in PS3/P2 | 0.0155 (1.29 px) | 0.0060 (0.50 px) | **+0.0054 (0.45 px)** |
| A9 | 0.025 steps, 0.300…0.850 | 0.00050 (1.51 px) | 0.00019 (0.56 px) | **+0.00013 (0.40 px)** |

Both pass, and both leave the *same* residual: every abscissa about **+0.4 px** to the
right of its lattice slot, on two different figures. The segment check above rules out a
marker-position bias (a +0.45 px displacement along x would put A2's centres 0.34 px off
their segments perpendicular, and the measured mean is +0.06 ± 0.04). So it is an offset in
the x axis map, inside its own ±1–2 px frame-test spread. **It is left in the CSVs, not
corrected out** — the files record what is on the page and the script reproduces them — but
a consumer that wants the structural values should read A2's abscissae as the integers.
Snapping A2's x to integers would move SHP by −0.03 % (hover) to −0.17 % (descent).

#### Read at magnification: what was actually seen

* **A2 marker 20** sits *on* the right frame line *and* on top of the printed 2.8 major
  tick; its lower-left arm is inseparable from the incoming joining segment. Frame band and
  tick both masked; the 2.8 tick dropped from the right-edge fit.
* **A2 marker 1**'s upper-right arm merges with the outgoing segment, which is exactly why
  the segment is *modelled* rather than avoided.
* **A9 marker 1** sits on the left frame with 0.001 minor ticks running under it;
  **A9 marker 23** sits on the right frame, 25 px above the 0.26 long tick. Both masked.
* **A9's joining polyline is 1 px wide over the flat run and the scan has broken it in a
  dozen places.** "Largest 8-connected component" therefore returns a *fragment* — it cut
  the first two markers off, and the lattice seeding then put six seeds on the wrong part
  of the curve. `curve_component(..., close=3)` closes the mask to decide membership and
  returns the original ink. **Check the component's bounding box against the frame before
  trusting it.**
* A9's glyphs show **six** arms, not four: the four stroke ends plus the two joining
  segments leaving left and right along a nearly flat curve.
* Neither figure has a split or merged marker. 20 and 23, matching `inventory-appendix-a.md`.

#### Uncertainty, and whether either map can still explain 1 % of SHP

Per-point 1σ by quadrature of three named terms, all printed by the tools: *calibration*
(the polynomial's prediction standard error at that point, from the tick scatter), *glyph*
(the polyline-masked disagreement × local scale), and *model* (the spread over 8–16
defensible calibrations — frame arcs order 1 or 2, each axis order 1 or 3 or 4, edges
blended or pooled, and on A9 long ticks vs all ticks).

| | 1σ x | 1σ y | 1σ y as % of value |
|---|---|---|---|
| `f2` | 0.0066–0.0107 in PS3/P2 | 0.00044–0.00063 in T3/T2 | 0.019–0.042 % |
| `f9` | 0.00015–0.00044 in Ps9/P45 | 0.000024–0.000037 LBM/SEC | 0.007–0.011 % |

Pushed through the model as a rigid shift of the whole curve and re-trimmed, that is:

| perturbation | ΔSHP hover | level | descent |
|---|---|---|---|
| `f9` +1σ in y | +0.069 % | +0.020 % | +0.023 % |
| `f9` +1σ in x | +0.033 % | +0.024 % | +0.040 % |
| `f2` +1σ in y | −0.072 % | −0.103 % | −0.277 % |
| `f2` +1σ in x | +0.067 % | +0.105 % | +0.334 % |

Combined in quadrature: **±0.12 % (hover), ±0.15 % (level), ±0.44 % (descent)** of SHP at
1σ. **Neither map can account for 1 % of shaft power any more at hover or in level flight,
even at 2σ.** At the descent trim `f2` alone still could — 2σ is ±0.88 % — and it is `f2`'s
sensitivity there (−8.05 % SHP per +1 %) rather than any weakness in the reading that makes
it so. The `f2` x term is the largest single contributor and is the one the integer lattice
above would remove.

#### Flag: A7's own y model is not settled, and `f7` may still move ~0.25 %

Applying this method's held-out test to p.62 turns up something the A7 redo could not see,
because it tested against **one** printed value. A7's y has a second known value: the bottom
frame is exactly 20.0, fixed by an integer tick count (13 steps of 2.5 from the printed
52.5), not by extrapolation. Fitted on the 12 interior majors per edge alone:

| order | interior rms | predicts bottom (20.0) | predicts top (52.5) |
|---|---|---|---|
| 1 | 0.030 | 19.965 | 52.568 |
| 2 (in use) | 0.017 | **19.892** | 52.496 |
| 3 | 0.011 | 19.952 | 52.437 |
| 4 | 0.0058 | **20.016** | **52.500** |

The quadratic passes the 52.5 test and **fails the 20.0 test by 0.108 BTU/LBM**; the quartic
passes both and cuts the interior residual threefold. Swapping order 2 for order 4 moves
`f7` by −0.017, +0.079, +0.097, +0.071, +0.061, −0.089 BTU/LBM at the six points —
**+0.20 to +0.25 % at the three trim knees**, three times A7's stated ±0.03 σ, and worth
about +0.6 % of SHP through `f7`'s +3.1 %/+2.7 %/+5.3 % sensitivities. This is **not** acted
on here (A7 is outside this task) and it is not a claim that order 4 is right — 5 parameters
on 12 ticks is thin. It is a claim that **A7's y order was chosen on one held-out number and
one number cannot separate an even from an odd correction.** p.62 also prints minor ticks;
redoing A7 on the fine lattice, with the second frame value in the test, would settle it.

#### Not yet applied elsewhere

`tools/digitize_a6810.py` already picks its axis order by the held-out frame test (with a
leave-one-out guard against a high order winning on noise), so A6, A8 and A10 have the
*axis* half of this. None of A1, A3–A8, A10, A11 has the *geometry* half: they all use a
homography, i.e. four straight edges. The sagitta table above says that is worth 0.4–1.4 px
mid-page, and on A2 it was the difference between an 0.95 px and an 0.50 px free-lattice
scatter. Cheap to retrofit — `frame_curves(order=2)` and `make_sv` replace `frame` and
`homography` — and worth doing the next time any of those figures is touched. Nothing about
this changes their published values by more than their stated read error, so it is not a
reason to redo them on its own.

#### Files

* `tools/digitize_native.py` — the shared machinery (native bitmap, frame arcs, `make_sv`,
  `ticks`/`label_ticks`/`label_by_map`, `held_out`/`tail_cv`, the glyph model).
* `tools/digitize_a2.py`, `tools/digitize_a9.py` — one per figure; both rerun and reproduce
  their CSVs exactly. Run from the repo root with `PYTHONPATH=tools`.
* Overlays: `validation/out/digitize/a2redo/a2_verify*.png`,
  `validation/out/digitize/a9redo/a9_verify*.png`.
* Both CSVs gained **`sx` and `sy` columns** carrying the per-point 1σ in data units.
  `maps.load_curve` selects columns by name, so the extra columns are inert to the model.

### Update — A7 calibration settled (2026-09-11)

The A2/A9 note above flagged A7's y model as unsettled: its order had been chosen on **one**
printed frame value, and one value cannot separate an even axis correction from an odd one.
Settled here, on p.62's **minor-tick lattice**, with **both** y frame values in the test.
Extractor rewritten as `tools/digitize_a7.py` over `tools/digitize_native.py`; it reruns and
reproduces the six rows of `data/maps/f7_gg_turbine_energy.csv`, which now carries `sx`/`sy`.

**Answer: y is order 5, x is order 3. The quadratic in use was wrong — but the published
values move by at most 1.4σ, and the flag's "~0.25 %" was 6–8× too large.** Both halves of
that sentence are the finding.

#### The unlabelled frame, fixed by counting, and then used as a second test

The A6/A8/A10 free check ports directly and it is what makes the whole exercise possible.
Measured by a **straight line through the ticks alone**, A7's frame spans

| edge | in majors | nearest | in minors | nearest |
|---|---|---|---|---|
| left (y) | 13.029 of 2.5 | **13** | 130.29 of 0.25 | 130 |
| right (y) | 13.033 | **13** | 130.33 | 130 |
| bottom (x) | 14.996 of 0.01 | **15** | 74.98 of 0.002 | 75 |
| top (x) | 15.022 | **15** | 75.11 | 75 |

so both unlabelled frames are lattice points: **bottom = 52.5 − 13×2.5 = 20.0 exactly**,
**left = 0.35 − 15×0.01 = 0.20 exactly**, with 0.97 of a major step of margin to the wrong
integer. Read as a rounding test only — the straight line leaks the bow's odd part into the
span — but the margin is nowhere near close. That integer is the second known value the
order test needs.

#### The held-out table, and what the minor ticks buy

Fit the interior ticks alone; predict both frame values; guard with leave-one-tick-out.
Left edge, **129 minor ticks at 0.25** (right edge, 125, in brackets):

| order | interior rms | LOO | predicts 20.0 | predicts 52.5 |
|---|---|---|---|---|
| 1 | 0.0333 [0.0402] | 0.0338 [0.0408] | −0.026 [−0.022] | +0.050 [+0.065] |
| 2 *(was in use)* | 0.0233 [0.0263] | 0.0240 [0.0271] | **−0.081 [−0.091]** | −0.004 [−0.004] |
| 3 | 0.0144 [0.0153] | 0.0150 [0.0159] | −0.029 [−0.032] | −0.055 [−0.062] |
| 4 | 0.0086 [0.0089] | 0.0090 [0.0093] | +0.009 [+0.008] | −0.018 [−0.023] |
| **5** | 0.0073 [0.0078] | 0.0077 [0.0081] | **−0.009 [−0.008]** | **−0.001 [−0.007]** |
| 6 | 0.0073 [0.0078] | 0.0076 [0.0082] | −0.012 [−0.010] | −0.004 [−0.008] |

(BTU/LBM.) Orders 4, 5, 6 pass the LOO guard on both edges; order 5 has the smallest
worst-frame error on both. **The quadratic hits the value it was chosen on and misses the
one it was not**, by 0.081–0.091, with 3× order 5's interior residual — exactly the failure
the flag predicted. On the **long ticks alone** the same procedure stops at order 4: with 12
ticks per edge the LOO guard rejects 5 and 6 outright. **The fine lattice is not a
refinement here, it is the thing that lets the test reach the right order at all.** x, on 74
minor ticks per edge at 0.002: LOO rejects orders 1 and 2 (1.8–2.0× the best), and 3–6 tie on
the frame test within 1.26×, so parsimony takes **order 3**.

#### Did the minor ticks separate even from odd? Yes — and the answer is "both"

Legendre coefficients of the y value function, in page pixels, left edge [right]:

| | P2 (even) | P3 (odd) | P4 (even) | P5 (odd) | P6 (even) |
|---|---|---|---|---|---|
| coefficient | −3.77 [−4.83] | −3.52 [−4.14] | +2.56 [+2.71] | +1.13 [+1.08] | −0.22 [−0.09] |
| significance | 35σ [41σ] | 27σ [30σ] | 17σ [17σ] | 6.9σ [6.0σ] | 1.2σ [0.4σ] |

**The bow is even *and* odd at comparable amplitude, through degree 5.** That is why no
quadratic and no cubic can pass both frame tests, and why one printed frame value could
never have decided between them — a pure-even and a pure-odd correction can hit the same
single value from opposite sides. It also retires "A7 found a quadratic" and "A2/A9 found a
cubic" as rules: this page needs a quintic, and 129 ticks can tell you so where 12 cannot.
On x the same decomposition gives P3 = −1.6 [−1.8] px at 14σ with everything else at 1–4σ:
a linear map plus a cubic wiggle, which is the order the rule picked.

#### An eighth failure mode of its own: anchor the frame value HARD

`blend_axis` carries a printed frame value as **two extra data points** among the ticks —
2.6 % of the weight on a 76-tick edge. The frame value is *known*, and its pixel position
comes from a ~1900-row fit of the frame line, so it is the best-determined datum on the
axis. With the soft anchor A7's order-3 x map misses the left frame by **2.0 px** — and
**marker 1 sits on that frame**, so the miss lands straight in the CSV (x = 0.200229 instead
of 0.200044, 2 px off a lattice the extraction was never told about). Fixed by fitting
`lo + (hi−lo)t + t(1−t)q(t)`, exact at both ends by construction. On A7's y at order 5 hard
and soft anchoring differ by ≤0.0006 BTU/LBM, so this is an x-axis, marker-on-the-frame
problem, not a general one — but it will bite on any figure whose first or last marker sits
on a frame line, which is most of them.

#### A ninth: a seed ten pixels out, and a converged-looking wrong answer

Marker seeds come from a 6×6 morphological opening (six cores, no threshold sweep, and an
independent marker count). For a marker sitting **on** a frame line the core is the union of
glyph, frame and tick, so its centroid is up to **10 px** off. Handed that, the offset-simplex
Nelder-Mead fit slid marker 1 **18 px down the joining line** into a clean-looking local
minimum: cost converged, template plausible, every printed diagnostic normal, and f7 low by
**0.27 BTU/LBM** at that point. A 1 px grid search of the hard-edged glyph model over ±16 px
before the simplex fixes it — a grid search cannot slide, it evaluates everywhere. With it,
all six centres land within 0.66 px of their seeds and **within 0.25 px of the first redo's
independently fitted centres**. Rule: *seed a local optimiser with a global search, and
assert on the distance it then moves.*

`digitize_native.fit_glyphs` gained an optional `simplex=` argument that optimises the
centre as an **offset** with an explicit 2 px initial simplex (the A6/A8/A10 fix, now in the
shared machinery). Default `None` keeps the old path bit for bit, so A2 and A9 are untouched.

#### What actually moved, and why so little

| | first redo (2026-09-11) | settled | move | 1σ |
|---|---|---|---|---|
| 1 | 0.199937, 48.2815 | 0.200044, 48.2596 | +1.1e-4 (+1.4σ), −0.022 (−1.3σ) | 7.6e-5, 0.017 |
| 2 | 0.214923, 40.3985 | 0.214956, 40.4095 | +0.3e-4 (+0.4σ), +0.011 (+0.8σ) | 7.4e-5, 0.014 |
| 3 | 0.219917, 39.1268 | 0.219942, 39.1397 | +0.3e-4 (+0.5σ), +0.013 (+1.1σ) | 4.7e-5, 0.011 |
| 4 | 0.224969, 38.3672 | 0.224989, 38.3785 | +0.2e-4 (+0.4σ), +0.011 (+1.0σ) | 4.5e-5, 0.012 |
| 5 | 0.229970, 37.5898 | 0.229971, 37.5973 | +0.0e-4 (+0.0σ), +0.008 (+0.7σ) | 4.5e-5, 0.012 |
| 6 | 0.350096, 24.6064 | 0.350109, 24.6074 | +0.1e-4 (+0.2σ), +0.001 (+0.1σ) | 8.4e-5, 0.015 |

**Everything is inside 1.4σ, and the flag's estimate of "+0.20 to +0.25 % at the knee" was
6–8× too large.** The flag compared *unanchored long-tick* fits of order 2 and order 4. The
deployed maps are anchored, and the corrections partly cancel. Step by step, at the knee:

| step | m2 | m3 | m4 | m5 |
|---|---|---|---|---|
| old calibration, old centres → published CSV | −0.0006 | −0.0006 | −0.0006 | −0.0004 |
| new centres instead of old (i.e. **detection**) | −0.002 | +0.001 | +0.001 | −0.001 |
| anchor the derived 20.0 at the bottom frame | −0.011 | −0.011 | −0.010 | −0.010 |
| fine lattice instead of long ticks | +0.004 | +0.005 | +0.004 | +0.004 |
| frame arcs instead of the homography | −0.002 | −0.002 | −0.002 | −0.003 |
| order 2 → 3 | +0.018 | +0.013 | +0.010 | +0.006 |
| order 3 → 4 | +0.006 | +0.010 | +0.012 | +0.014 |
| order 4 → 5 | +0.005 | +0.004 | +0.003 | +0.002 |
| per-edge fit, blended, hard anchors | −0.008 | −0.007 | −0.006 | −0.005 |
| **net** | **+0.011** | **+0.013** | **+0.011** | **+0.008** |

The first row is the check that licenses the comparison: **rebuilding the superseded
calibration in full and applying it to the superseded centres reproduces the published CSV
to 0.0006 BTU/LBM**, so this is like for like. The second row is the whole detection
contribution and it is ≤0.002 — **the ink was right both times; every one of these updates
has been the axis map.** The order alone is worth **+0.029 to +0.022** at the knee (0.06 %),
which is real and is the flag's point vindicated in direction if not in size; anchoring the
bottom frame takes about a third of it back, and the per-edge blend most of the rest.

Free checks on the settled values: the six abscissae land on the 0.005 lattice to **max
1.20 px, sd 0.64 px** (the 1.20 is marker 6, which genuinely sits that far outside the right
frame line; the four interior markers are within 0.65 px); the local **two-tick
interpolation** — bow-free, since it never fits a bow — agrees with the quintic to **0.013
BTU/LBM** at every point, which is the statement that the order-5 map is not wandering
between ticks; masking the joining polyline out entirely instead of modelling it moves the
centres ≤0.71 px; and the detector-free interior-ink walk still finds six markers.

#### Effect on the trims, and where f7 now stands

At the three Table B.1 trim arguments f7 rises **+0.070 % (hover), +0.049 % (level),
+0.032 % (descent)** — note it rises *most* at hover, where the trim sits on the steep first
segment, and least at descent, contrary to the previous redo. Re-trimming:

| | NG before → after | SHP before → after |
|---|---|---|
| hover | −0.050 % → **+0.014 %** | −0.257 % → **−0.048 %** |
| 80 kt level | +0.105 % → +0.132 % | −0.180 % → **−0.049 %** |
| 80 kt descent | −0.046 % → −0.026 % | −0.910 % → **−0.730 %** |

The implied sensitivities (3.0, 2.7, 5.6 % of SHP per 1 % of f7) match the sweep's
+3.1/+2.7/+5.3, so nothing surprising is happening in the loop.

Remaining 1σ headroom in f7, shifting the whole curve coherently and re-trimming: y ±0.100 /
±0.081 / ±0.171 and x ±0.289 / ±0.095 / ±0.104 percentage points of SHP, i.e. **±0.31 %
(hover), ±0.13 % (level), ±0.20 % (descent)** in quadrature. The x term dominates at hover
because that trim sits on the −526 BTU/LBM segment. So hover and level are now **inside 1σ
of f7 alone**, and the descent condition's remaining −0.73 % is **3.7σ away from anything f7
can supply**. A7 is finished; the descent deficit is somewhere else.

#### Not yet applied elsewhere

The two things this note adds — **hard frame anchoring** and **grid-searched seeds before the
local fit** — are not in `digitize_a2.py`, `digitize_a9.py` or `digitize_a6810.py`. The
anchoring matters wherever a marker sits on a frame line (A2's marker 20, A9's 1 and 23,
A6's both, A8's and A10's first and last); the seeding matters wherever a seed can be far
off, which is the same list. Neither is expected to move a published value by more than its
stated read error, so neither is on its own a reason to redo those figures — but the next
time any of them is touched, both are cheap to carry across. A1, A3, A4, A5, A11 still use
the straight-edge homography as well.

## Attempted 2026-09-11: porting A3, A4, A5 and A11 to the `digitize_a6810` machinery

`f3`, `f4`, `f5` and `f_hs` are still the 2026-09-10 extraction by the generic
`tools/digitize.py` (300 dpi render, `rectify_page` homography, least-squares tick fit)
while the other seven maps were redone the next day on the native bitmap with frame arcs,
held-out frame calibration and leave-one-tick-out order selection. The four were measured
first and found to be worth at most 0.21 % on any NG trim and 0.37 % on `tau1/tau2` per
1 % of map movement, so this was a consistency exercise, not a correction.

**It did not complete, and the reasons are worth recording so the next attempt starts
further along.**

### What worked: the axis calibration on A11

Frames located on the native bitmap at left 475.0, right 2133.5, top 305.5, bottom 2836.7
(plot 1658 x 2531 px). Both axes chose **order 1** on the held-out frame test, and the
lattice closed on the printed span: x spanned 7.0435 major steps against a printed 7, y
spanned 9.0112 against 9.

That is an independent cross-check of the existing extraction, by a different method on a
different raster:

| | printed | new calibration recovers |
|---|---|---|
| x low frame | 65 | 64.806 |
| x high frame | 100 | 99.894 |
| y low frame | 3.50 | **3.49728** |
| y high frame | 8.00 | **8.00012** |

The y axis -- the one that sets `f_hs` and therefore `tau_b` -- comes back to within
**0.0044 in T41SGN**, or 0.06 % of the 7.497 plateau, against the older extraction's stated
0.0104. Two independent calibrations of the same figure agree to better than a tenth of a
per cent, which is worth more than either alone.

### What did not: the marker seeder does not transfer to A11

`seeds()` uses a 7x7 morphological opening, which works where markers are isolated 'x'
glyphs. A11's six markers sit **on** the joining polyline, which is why the original used a
density detector with a window and a fill threshold. Sweeping the opening size gives 15
components at 7x7, 8 at 8x8, 2 at 9x9 and nothing at 10x10 -- it never passes through 6.
The glyph strokes are too thin relative to the line for an opening to separate them.

Also noted: A11's **bottom frame is faint**. Tracing it yields 105 points against
1442-1558 for the other three edges. The frame fit there rests on far less ink than
elsewhere.

### And A3, A4, A5 need their tick lattices read off the page

All three fail in `label_ticks` with "no major near slot ...", which means the `x`/`y`
config -- low and high frame values, major step, count of interior majors, which frame
values are actually printed -- does not match what the page prints. Those were inferred
from `inventory-appendix-a.md`'s axis table rather than read off the rendered page, and
the inventory records the *labelled* majors without saying whether unlabelled majors sit
between them, which is exactly what the config needs.

**Next attempt should start by rendering each page and reading its tick lattice
directly**, the way `digitize_a7.py`'s header shows was done for A7. Until then the
existing extraction stands: it is sound, its uncertainty is stated, and the measured
effect of any plausible revision is below 0.05 % on every result in this repository.

## Appendix C scheduling functions — state as of 2026-09-11

The eight functions `F_EC1` and `F_HM1`..`F_HM7` exist only as plots (Figs. C23-C30,
pdf pp.94-101). All eight pages are single 300 dpi 1-bit CCITT rasters, so
`digitize_a6810`'s native-bitmap machinery applies directly rather than needing a render.

**`F_HM1` (Fig. C24) is done** -- `data/schedules/fhm1_topping_line.csv`, nine points with
per-point uncertainties. Both axes calibrated to better than 1.3 per-mille of range on the
held-out frame test, the tick lattices closing at 5.0100 against a printed 5 in x and
9.0118 against 9 in y. The six interior abscissae land on a 20 deg R lattice to within
0.19 against a mean sigma_x of 0.22, which the extraction was never told to expect.

Two changes were needed and both are now parameters with defaults that leave A6, A8 and
A10 reproducing byte for byte:

* **`open_sz`** -- the Appendix C plots are drawn with a lighter pen. On C24 the default
  7x7 morphological opening *erases* the ninth marker outright; it is not filtered out, it
  is never found. 6x6 keeps it.
* **`rail`** -- 6x6 also keeps a rail of fragments along the frame lines. A6/A8/A10 need
  the default, which keeps components up to 14 px *outside* a frame because their end
  markers sit on it; C24 needs the opposite, discarding anything within 20 px *inside* it,
  because its markers all stand at least 33 px clear.

### What blocks the remaining seven

**A bug in the frame locator, not in the tool.** The helper that finds the printed frame
by taking the two strongest horizontal lines picks the wrong bottom line on at least
C25 and C26: the detected frame spans fifteen major intervals where the printed labels say
fourteen, putting every predicted tick slot about 188 px out and failing
`_pick_majors` with "no major near slot ...". The tick candidates themselves are clean --
C25 gives fourteen evenly spaced at 167.5 px, C26 sixteen at 157 px -- so once the frame
is right the rest should follow.

Fix the frame locator before adding more configs. Candidate approach: find the frame from
the *tick lattice* rather than from ink alone -- the correct bottom frame is the one that
makes the interior tick count match the printed label count.

### Axis lattices already read off the pages

| fig | function | x | y | markers |
|---|---|---|---|---|
| C24 | `F_HM1` | T2 390-640, step 50, both frames printed | WFPTP -0.50-4.00, step 0.50, both | 9 **done** |
| C25 | `F_HM2` | PAS 20-120, step 20, both | WFPRF -2-12, step 1, both | 11 |
| C26 | `F_HM3` | XLDSA 0-100, step 20, both | PNG 75.0-112.5, step 2.5, both | 11 |
| C23 | `F_EC1` | W45R 0-16 | TAU45 1-6 | 4 curves, digit markers |
| C27 | `F_HM4` | XLDSA 0-100 | WFQPS3 2.1-3.7 | 11 |
| C28 | `F_HM5` | T2 390-640 | WFIRF 2.05-2.65 | 13 |
| C29 | `F_HM6` | T2 390-640 | PCNGI 62-~75 (top tick unlabelled) | 2 |
| C30 | `F_HM7` | PCNGHL 50-110 | WFPAC 0-5 | 7 curves, crossing bundle |

C23 and C30 will need more than a frame fix: their markers are printed **digits**
identifying which curve a point belongs to, not `x` glyphs, and C30's seven curves cross
in a tight bundle around PCNGHL 84-95.

Frame windows for all eight pages are in the session notes and are re-derivable in a few
seconds; they are deliberately not pasted here, because the ones for C25 and C26 are
wrong and the next attempt should recompute them with a fixed locator rather than inherit
a known-bad table.
