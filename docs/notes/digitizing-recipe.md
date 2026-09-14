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

Step 4 is not optional. Every trap below was invisible in the CSV and obvious in the overlay.

## Nine traps, all of them hit at least once

They are numbered as they were found, not by severity. Each cost a wrong CSV that looked
right. The numbering is not the report's and means nothing outside this file.

**1. Markers sit ON the joining line, so blob detection is useless.** The curve and its
crosses form one connected component and `find_marks` returns a single meaningless centroid.
Use `--mode density` (the default): a marker is locally *denser* than the line through it,
so box-filtering the ink mask and taking local maxima finds them.

**2. Calibrating by eye off a grid overlay puts every point out.** Reading tick coordinates
from a downscaled image cost ~25 px, which on Figure A2 is 0.3 in PS3/P2 — enough to move
every value, small enough to look plausible in the CSV. **Always calibrate with `axes`.**

**3. A marker on the frame gets clipped** by the `axes` default 12 px inset. Figure A2's
last point, at PS3/P2 = 20, needed the ROI extended past the right frame. If the count comes
up one short, look at the corners first.

**4. `axes` mis-locates the frame when the page carries a long flat data line.**
`find_frame` averages every strong row in each half; on Figure A3 the flat 0.11 and 0.00 runs
are strong rows and the skew-diluted bottom frame is not, so it returned top/bottom
362/2826 against a measured 313/2672 — 0.009 in B1, five times the cost of trap 2. Five of
the eleven maps go flat somewhere. **Check the frame `axes` reports against the printed
ticks before using it.**

**5. A linear axis map cannot fit these scans, and the tick lattice says so.** The page bows
along the scan direction: A7's y major spacing runs 181 → 177 → 180.5 px top to bottom, a
smooth 2 % swing worth ±4 px against a straight line. It is a fittable, checkable quadratic
— linear over A7's 24 printed y majors predicts the top frame at 52.5624 against a printed
52.5; quadratic, never told the answer, gives 52.4947.

**6. The held-out frame test alone over-fits; leave-one-out is the guard.** On A10's x the
quartic wins the frame test (1.11 per-mille against linear's 1.94) while its
leave-one-tick-out error is 40 % worse. **Choose the order minimising held-out frame error,
but only among orders whose leave-one-tick-out error is within 1.2x of the best; ties within
1.25x go to the lower order.**

**7. The frame edges are arcs, so a homography is wrong too** — it assumes four straight
edges. Fit each printed edge as a quadratic along its own coordinate and read the sagitta.
Related: **do not resample the page.** Warping pixels needs interpolation; fit the
distortion into the axis map instead. If you already warped, measure what it cost rather
than redoing the work.

**8. `Nelder-Mead` builds its initial simplex by perturbing each coordinate by 5 % of its
own value.** At an image column of x ~ 2269 that is a 113 px step, outside any sane scoring
window, where the cost is flat because no model ink reaches it. On A6's right-hand marker
the fit took that step and never came back: the centre was 125 px off the data line and
every printed diagnostic looked converged. **Never hand Nelder-Mead an absolute pixel
coordinate** — optimise the offset from the current centre with an explicit 2 px simplex,
and assert no marker moves more than 10 px from its seed.

**9. Mask the frame symmetrically or the abscissa is biased.** Ticks protrude inward only,
so a marker on a frame line has the inward half of its disc stripped and the outward half
left, pulling the fitted centre outward. Mask |d| < k+5 on *both* sides of each frame line.

## Checks the figure gives you for free

The tool is never told to expect any of these, which is what makes them checks.

- **Regular structure.** Figure A2's markers sit at each integer of PS3/P2; the extraction
  recovers 1.02, 2.03, 3.00 ... 19.99 — integers to within 0.04.
- **A printed frame value held out of the fit.** Fit the ticks without it and predict it.
  On A9's x only the cubic lands on the printed 0.85 to better than 2e-4, and it hands back
  a measurement of the *unprinted* left frame, 0.3008, confirming a 0.30 that had been
  inferred from a tick count.
- **A local two-tick interpolation.** Interpolate between the two majors bracketing the
  point on each ticked edge and average. It cannot see a global bow because it never fits
  one; where it agrees with the chosen polynomial, that polynomial is not wandering between
  ticks. Worst disagreement: A6 5.9e-6, A8 0.0225 BTU/LBM, A10 8.7e-5.

**Two rules that outrank any of them: assert on structure, not on count — and a count that
matches is not a check.** The A2 y detector once returned exactly 189 minor ticks, the
exactly correct number, while holding one spurious tick and missing one real one; the
off-by-one labelling moved the axis by 0.03 in T3/T2 and every count-based check passed.

## Read error carried by the maps

What the digitization contributes to a trim, measured by pushing +/-1 sigma of each map's
fitted axis uncertainty through the model:

| map | hover | level | descent |
|---|---|---|---|
| `f9` +1 sigma in y | +0.069 % | +0.020 % | +0.023 % |
| `f9` +1 sigma in x | +0.033 % | +0.024 % | +0.040 % |
| `f2` +1 sigma in y | -0.072 % | -0.103 % | -0.277 % |
| `f2` +1 sigma in x | +0.067 % | +0.105 % | +0.334 % |

Per-map axis uncertainty, worst over the three trims: `f6` +/-0.125 %, `f8` +/-0.278 %,
`f10` +/-0.078 %. `f6` is one number, 0.98500 +/- 0.00041, and its two markers agree with
each other to 7.1e-5 — which is the figure's own statement that it is a constant.

SCOPE.md carries the *figures'* read error separately, measured a different way: Figures 9
and 10 each print one panel whose true value is stated, so digitizing it measures that
page's error directly.

## Which tool owns which figure

`tools/reproduce_all.sh` runs all of them and diffs the result against `data/`. Verified
against each tool's own output paths rather than from memory.

| tool | figures | writes |
|---|---|---|
| `digitize_a1.py` | A1 (p.56) | `f1` compressor mass flow — the one 2-D map, plus its beta-line form |
| `digitize_a2.py` | A2 (p.57) | `f2` |
| `digitize_single_curve.py` | A3–A6, A8, A10, A11 (pp.58–66); C24–C29 (pp.95–100) | `f3`, `f4`, `f5`, `f6`, `f8`, `f10`, `f_hs`, and `F_HM1`–`F_HM6` — thirteen figures |
| `digitize_a7.py` | A7 (p.62) | `f7` |
| `digitize_a9.py` | A9 (p.64) | `f9` |
| `digitize_multi_curve.py` | C23 | `F_EC1`. Refuses C30 by design — its curves cross, so rank in y is not identity |
| `digitize_c30.py` | C30 (p.101) | `F_HM7`, the acceleration fuel limit — the one figure needing per-curve tracking |
| `digitize_fig678.py` | 6, 7, 8 | `data/reference/fig06–08_*.csv`, three series each |
| `digitize_fig910.py` | 9, 10 | `data/reference/fig09–10_*.csv`, six panels each |
| `digitize.py`, `digitize_native.py`, `rectify_page.py` | — | shared libraries; they write nothing |

`f2`, `f7` and `f9` are the three whose CSVs are hand-maintained rather than generated —
their headers carry provenance no generator writes — so those tools *check* their file
instead of rewriting it. CLAUDE.md records why that distinction matters.

---

The per-figure working notes that used to follow — nine dated "Update: A7 recalibrated"
sections, 1,350 lines of them — were cut on 2026-09-14 once every map reproduced byte for
byte under `tools/reproduce_all.sh`. They are in git history:
`git log -p docs/notes/digitizing-recipe.md`.
