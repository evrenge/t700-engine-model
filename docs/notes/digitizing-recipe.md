# Digitizing recipe

Nineteen figures in the appendices carry model data — 11 engine functions, 8 control
schedules — and five more in the body carry validation references. **None of them exists as
numbers anywhere.** This file is what has been learned turning them into CSV: how the
pipeline works, what it refuses to do, and the traps that each cost a wrong file that
looked right.

Rewritten 2026-09-18. The version before it opened with a worked example built on
`tools/digitize.py`, a CLI that was deleted that week, and taught a two-point linear axis
calibration that traps 5, 6 and 7 further down the same page had already demolished. If
you are looking for it, it is in `git log -p docs/notes/digitizing-recipe.md`.

## The shape of the thing

One library, `tools/figure/`, and a script per figure that supplies the figure's own facts
and nothing else:

    raster.py       the page bitmap: the scan's own 300 dpi 1-bit image, never resampled
    frame.py        the four printed frame edges, and the pixel -> plot-coordinate map
    ticks.py        the printed tick lattice, and each tick's printed value
    axis.py         the axis map: order chosen by held-out test, anchored, blended
    glyph.py        the marker: a generative model of the local ink, fitted
    provenance.py   the CSV and its header -- one writer, and the parser that reads it back

Until 2026-09-16 that machinery existed five times over, and the cost was not the
duplication. It was that an improvement made for one figure did not reach the others and
nothing said so — the frame model, the Nelder-Mead guards and the anchoring rule were each
fixed in one copy and not the rest, and it took reading two files side by side to see it.
`tests/test_tool_duplication.py` now fails if the same function body appears in two files
under `tools/`.

## What the tool decides, and what you decide

You supply, per figure: the page, the printed value at each frame edge, which of those are
actually printed, the major-tick step, how many markers to expect, and the marker radius.
Everything below that is measured or chosen by a test:

| question | how it is answered |
|---|---|
| where the plot frame is | `find_frame_box`, by shearing the page and taking the extreme rows that span it |
| straight edges or arcs | both fitted, both calibrated; keep whichever predicts the printed frame values better |
| how many minors per major | major spacing ÷ median tick gap, which must come out an integer |
| majors alone or the full lattice | whichever predicts the printed frame values better, per axis |
| polynomial order of the axis map | lowest held-out frame error among orders whose leave-one-tick-out error is within 1.2× of best |
| where each marker is | a generative model of two strokes plus the joining half-segments, fitted sub-pixel |

**The recurring shape: the figure prints a number, you fit without it, you predict it.**
That is the only test available that is not a fit statistic, and every model choice here is
made on it.

## The traps

Numbered as they were found, not by severity. Each cost a wrong CSV that looked right.

**1. Markers sit ON the joining line, so blob detection is useless.** The curve and its
markers are one connected component. Model the marker instead: two strokes plus the two
joining half-segments, all radiating from one point, scored against ink **and background**
over an unmasked disc. Scoring the background is what stops a near-parallel pair of lines
from outscoring a four-armed glyph.

**2. Calibrating by eye off a grid overlay puts every point out.** Reading tick coordinates
off a downscaled image cost ~25 px, which on Figure A2 is 0.3 in PS3/P2 — enough to move
every value, small enough to look plausible in the CSV.

**3. A marker on a frame gets clipped** by an inset, and **4. the frame finder is fooled by
a long flat data line** — `find_frame` averaged every strong row, and a flat run of data is
a strong row. On Figure A3 it returned top=362 against a measured 313. Both are why
`find_frame_box` thresholds on coverage of the plot *width* after shearing, rather than on
row strength: a footer rule is too short to qualify, a data line qualifies but lies inside,
and a frame — however broken — still spans the plot once the shear is undone.

**5. A linear axis map cannot fit these scans, and the tick lattice says so.** The page bows
along the scan direction: A7's y major spacing runs 181 → 177 → 180.5 px top to bottom, a
smooth 2 % swing worth ±4 px against a straight line. It is a fittable, checkable
polynomial.

**6. The held-out frame test alone over-fits; leave-one-out is the guard.** On A10's x the
quartic won the frame test while its leave-one-tick-out error was 40 % worse. Choose the
order minimising held-out frame error, but only among orders whose LOO error is within 1.2×
of the best; ties within 1.25× go to the lower order. **And that guard goes inert on a fine
lattice** — with 200 ticks the LOO error is set by individual tick scatter and barely moves
with order, so every order clears it. What protects you then is that the frame test is held
out: an order that bought the interior with noise predicts the frame worse and loses.

**7. The frame edges image as arcs — but fitting them as arcs is usually WORSE.** This trap
used to read "the frame edges are arcs, so a homography is wrong too", inferred from the
measured sagittae, which run 0.3 to 3.5 px on a printed straight line. Measured on 2026-09-18
across all thirteen single-curve figures, that does not generalise: **eleven of them predict
their own printed frame values better under the straight model.** Two do better with arcs
— A10 by 11×, C26 by 1.2× — and the tool now fits both and lets each figure choose, with
arcs required to win by 10 % because straight is the incumbent. (A11 was first measured as
a third, by 4.3×: a minor tick had stood in for a major hidden under a marker in its
straight trial. Fixed 2026-10-12.)

Two lessons, and the second is the larger one. A measured physical effect is not by itself
a reason to model it: the sagitta is real on every page and modelling it helps on three.
And **this question could not be answered on the coarse lattice.** On the majors alone A8
preferred arcs by 17.7×; with ten times the ticks it prefers straight by 7.9×. The tick
lattice had to be fixed first or the answer would have been confidently backwards.

Related, and unchanged: **do not resample the page.** Warping pixels needs interpolation,
and interpolation invents ink. Fit the distortion into the map instead.

**8. `Nelder-Mead` builds its initial simplex by perturbing each coordinate by 5 % of its
own value.** At an image column of ~2269 that is a 113 px first step, outside any scoring
window, where the cost is flat because no model ink reaches it. On A6's right-hand marker
the fit took that step and never came back: 125 px off the data line, every printed
diagnostic converged-looking. **Never hand Nelder-Mead an absolute pixel coordinate** —
optimise the offset from the current centre with an explicit 2 px simplex.

**9. Assert that no marker runs away from its seed.** Trap 8's cure is not complete without
it: 10 px is the bound in use, and A7 asserts 4.

**10. Mask the frame symmetrically or the abscissa is biased.** Ticks protrude inward only,
so a marker on a frame line has the inward half of its disc stripped and the outward half
left, pulling the fitted centre outward. Mask |d| < k+5 on *both* sides.

**11. The page prints five to twenty-five times more ticks than the majors, and they are
free.** Every one of these figures prints a minor lattice; until 2026-09-16 the tools
detected it, used it to build their masks and threw it away before calibrating. Using it
took the summed held-out frame error over 26 axes from 36.99 to 22.82 per-mille. Two things
make it work, and it fails without either: **measure the minor step** from the ticks
themselves rather than reading it off the page — major spacing ÷ median gap, which must come
out an integer, and does to within 1.7 % on all 52 edges — and **label the fine ticks
through a fit made on the majors alone**, never through the straight line between the frame
corners, which is several pixels out in mid-page and discards a contiguous block of good
ticks. A9 lost 29 of its 141 that way.

A11 is the cautionary one for reading a lattice by eye: it prints **three** tick lengths —
a minor, a longer one every ten minors at the half-major, and the major — so counted by eye
its half-major ticks pass for majors and give 10 where the answer is 20.

**12. A shared function with two plausible input conventions is a trap.** `protrusion` takes
a polynomial, highest power first; A1 holds its edges as `(intercept, slope)`. Passing the
pair unreversed probes a line with its slope and intercept exchanged, and it walked the
probe a million rows off the page. Loud, immediately — but the next one might not be.

## Checks the figure gives you for free

The tool is never told to expect any of these, which is what makes them checks.

- **Regular structure.** A2's markers sit at each integer of PS3/P2 and the extraction
  recovers them to within 0.016; A9's abscissae fall on a 0.025 lattice and land within
  0.00048.
- **A printed frame value held out of the fit.** The backbone of every model choice above.
  On A9's x it also handed back a measurement of the *unprinted* left frame, 0.3008,
  confirming a 0.30 that had been inferred from a tick count.
- **A local two-tick interpolation.** Interpolate between the two ticks bracketing each
  point on each ticked edge and average. It cannot see a global bow because it never fits
  one, so where it agrees with the chosen polynomial, that polynomial is not wandering
  between ticks. On a fine lattice the brackets are 5–25× tighter, which makes this much
  sharper than it was.
- **Two markers that should be one number.** A6 prints a horizontal line and two markers;
  they agree to 0.016 σ, which is the figure's own statement that `f6` is a constant.

**Two rules that outrank any of them: assert on structure, not on count — and a count that
matches is not a check.** The A2 y detector once returned exactly 189 minor ticks, the
exactly correct number, while holding one spurious tick and missing one real one; the
off-by-one labelling moved the axis by 0.03 in T3/T2 and every count-based check passed.

## Read error, and what it is worth

Each CSV carries a per-point 1 sigma in its own units, built by quadrature from the tick-fit
residual, the spread over defensible polynomial orders, the disagreement between the two
ticked edges, and the glyph-centre robustness.

**What that uncertainty costs the model, remeasured 2026-10-12.** Each map is shifted by
exactly +1 sigma of its own per-point error, on one axis at a time, and all three trims plus
both open-loop transients are resolved. Steady column is the worst shaft-power change over
the three trims; transient is the largest change in any whole-curve panel rms, in percentage
points of that panel's excursion. Every map that carries a per-point sigma is here; `f1`
carries none:

| map, +1 sigma | SHP, worst of 3 trims | worst transient panel |
|---|---|---|
| `f7` x | **0.290 %** | fig9 t45 **+0.346** |
| `f2` x | **0.240 %** | fig10 %NG **-0.221** |
| `f8` y | **0.222 %** | fig9 torq45 +0.095 |
| `f2` y | 0.177 % | fig10 %NG +0.131 |
| `f7` y | 0.154 % | fig9 %NG -0.110 |
| `f8` x | 0.132 % | fig9 torq45 +0.095 |
| `f3` x | 0.129 % | fig10 torq45 +0.039 |
| `f6` y | 0.113 % | fig9 torq45 +0.053 |
| `f10` x | 0.064 % | fig10 torq45 +0.019 |
| `f9` y | 0.057 % | fig9 %NG -0.091 |
| `f10` y, `f3` y | <= 0.028 % | <= 0.014 |
| `f9` x, `f5` x/y | <= 0.014 % | <= 0.028 |
| `f_hs` x/y | **0.000 %** | fig10 %NG +0.019 / -0.014 |
| `f4` x/y, `f6` x | 0.000 % | <= 0.001 |

Three things in that table are worth more than the numbers.

**The ranking has changed, and the old table had the leader wrong.** It reported `f9` as the
most expensive map, at +0.069 % on hover shaft power. `f9` is now tenth at 0.057 %, and
`f7` leads at 0.290 % -- absent from the 2026-09-19 version of this table, which the sweep
had skipped. That is not a regression: `f9` was re-digitized and its read error
fell. The old numbers were measured on 2026-09-11 against maps that have since moved, and
they sat in this file being quoted for eight days after they stopped being true.

**`f_hs` costs exactly nothing at trim and is only visible in a transient.** It sets the
station-4.1 heat-sink time constant [Eq. 52-53], which does no work at equilibrium, so a
trim-only sensitivity table reports zero for it however wrong it is. It appears at 0.019
points of a Figure 10 panel. Any table like this one that stops at the trims will
systematically understate every dynamic-only quantity -- which is why the transient column
is here.

**A perturbation that improves the report comparison is not a licence to apply it.** Shifting
`f2` by +1 sigma in x moves Table B.1's rms by -0.063, i.e. better. That is a measurement of
how close the calibration sits to the target, not an instruction: tuning a digitized map
toward a printed table is the move CLAUDE.md forbids, and the sigma it would be tuned within
is the reason the comparison has any headroom at all.

Reproduce with `tools/sensitivity_sweep.py`, which perturbs a throwaway git worktree so the
real tree is never touched and runs one subprocess per case, because the loaded `Curve`
objects are cached and their arrays are frozen.

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
| `tools/figure/` | — | the shared library; writes nothing |
| `rectify_page.py` | — | squares up a page for A1's component labelling. **Resamples**, which nothing else does |
| `check_straightness.py` | — | audit: is each printed curve straight between its markers, measured in data space |
| `compare_models.py` | — | audit: what the frame model and the anchoring rule cost, figure by figure |
| `minor_lattice.py` | — | audit: what minor lattice each edge resolves, and how cleanly |
| `sensitivity_sweep.py` | — | audit: what each map's read error is worth, pushed through the model |

`f2`, `f7` and `f9` are the three whose CSVs are hand-maintained rather than generated —
their headers carry provenance no generator writes — so those tools *check* their file
instead of rewriting it, and return non-zero on any difference. CLAUDE.md records why that
distinction matters.

`digitize_fig678.py` and `digitize_fig910.py` still carry their own frame and axis code:
the body figures are stacked panels with tilted frames rather than single framed plots, and
they were not converted when the appendix tools were. They share the raster and the run
helpers and nothing else. That is the largest piece of duplication left.
