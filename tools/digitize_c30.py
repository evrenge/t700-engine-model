#!/usr/bin/env python3
"""Figure C30 -- `F_HM7`, the HMU maximum fuel-parameter limit during acceleration.

The last of Appendix C's eight scheduling functions, and the only one that needed its own
tool. Open question #50.

## Why `digitize_appc_multi` cannot do this one

That tool identifies a curve by its marker's **rank in y** within a column, which works
whenever the curves are ordered by the parameter and never cross. C30's seven curves
cross, twice: they converge into a bundle around PCNGHL 85-95 and come out reordered --
after it the order top to bottom is 4, 3, then 5, 6, 7, with 1 and 2 already descending.
Rank is not identity there, and the tool says so and refuses the figure.

## What this tool does instead

**It does not read the digits, and it does not trace through the crossings.** Both were
tried and both fail; the record is in open question #50. Template matching on the printed
digits is too weak because the curve runs through every glyph as a bar that correlates
about equally with all seven templates (30 of 97 glyphs confident), and stripping the bar
off first destroys the glyph, because any annulus wide enough to fit the local curve
direction also reaches the neighbouring curves. A column tracer with one-to-one
curve-to-run assignment follows the geometry well -- it reproduces the fan to the
digitized markers -- but at a crossing the minimum-displacement assignment is not
necessarily the true one, so the labels scramble.

What works is the observation that makes both unnecessary: **the curves are exact
polylines between their markers.** Measured on curve 4, whose flat section and descent are
unambiguous: a straight line from (96.76, 3.78) to (104, 1.97) predicts 3.22 / 2.72 /
2.47 / 2.22 at PCNGHL 99 / 101 / 102 / 103, against ink measured at 3.21-3.24 / 2.70-2.73
/ 2.45-2.49 / 2.20-2.23. So the figure is fully described by its vertices, and the job is
to find each curve's vertices rather than to follow its line.

`KNOTS` below declares, for each curve, the columns at which it is identifiable and a
**seed** for where it sits. The seed is the human-read part -- rank in y where the fan is
open and monotone in T2, the printed digit where it is not. The tool then **measures** the
value: at each declared column it takes the ink run nearest the seed and reports that
run's centre through the calibrated axes. So every number in the output is a measurement
of the scan, the seeds only say which run to measure, and a mis-seeded knot fails loudly
rather than silently shifting a curve.

## How the identity was established, knot by knot

* **PCNGHL 55-82**: the fan is open and ordered by T2, so rank is identity. Verified
  independently rather than assumed -- the seven glyphs at PCNGHL 68.1 all sit at one
  abscissa and their WFPAC values order 7:4.03 > 6:3.89 > 5:3.73 > 4:3.51 > 3:3.39 >
  2:3.22 > 1:3.01, exactly monotone in T2.
* **The bundle and the drops**: from the printed digits, read off the 300 dpi raster at
  600 % magnification. Curve 1's "1" at (87.2, 3.41) and (94.5, 0.99); curve 2's "2" at
  (90.4, 3.54), (92.5, 3.02) and (96.7, 1.99); curve 3's "3" at (94.7, 3.69); curve 4's
  "4" flat at 3.78 across (92.5), (94.8) and (96.8); and the terminal group, where the
  three surviving curves are labelled 5, 6, 7 top to bottom at PCNGHL 103.8.
* **There is no inversion at the left-hand end.** The curves start at PCNGHL 54 and curve
  1 is the lowest throughout; columns 61, 62 and 63 each resolve seven clean runs and are
  monotone in T2.

## What the output is checked against

Every segment between consecutive knots is sampled and required to lie on ink. Over 185
segments the mean coverage is **0.9850**, and the two weakest curves are 5 (T2 = 519.0,
0.970 mean, 0.64 worst) and 7 (T2 = 590.0, 0.959 mean, 0.48 worst). That is the acceptance
test, and `main` fails below the 0.98 floor.

(This said "mean coverage 0.9965, every curve between 0.990 and 1.000, worst single
segment 0.83" until 2026-09-13 -- flattering, on the figure this project calls its hardest
extraction. The CSV header the tool emits, and README, both carried the correct 0.9850;
only the docstring was wrong. Found by the 2026-09-13 accuracy audit.)

    python3 tools/digitize_c30.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))

from digitize_a6810 import (  # noqa: E402
    _apply,
    auto_frame_spec,
    calibrate_axis,
    corner,
    frame,
    homography,
    major_values,
    native_bitmap,
    ticks,
)

PAGE = 101
OUT = Path("data/schedules/fhm7_accel_limit.csv")

CFG = dict(
    page=PAGE,
    fig="C30",
    func="F_HM7",
    title=(
        "HMU maximum fuel-parameter limit during acceleration"
        " (Fig. C21, pdf p.93: WFPAC = F_HM7(T2, PCNGHL))"
    ),
    xcol="pcnghl_pct",
    ycol="wfpac",
    xdesc="PCNGHL sensed gas turbine speed, percent",
    ydesc="WFPAC fuel flow acceleration limit parameter, nondimensional",
    param="T2",
    param_unit="deg R",
    x=dict(lo=50.0, hi=110.0, step=10.0, n=5, printed=("lo", "hi")),
    y=dict(lo=0.0, hi=5.0, step=1.0, n=4, printed=("lo", "hi")),
)

# [pdf p.101, in-plot legend] SYMBOL 1..7 -> T2, deg R
T2_OF = {1: 395.0, 2: 430.0, 3: 460.0, 4: 480.0, 5: 519.0, 6: 555.0, 7: 590.0}

# Columns where the fan is open and rank in y IS identity, listed top (curve 7) to
# bottom (curve 1). At 55, 57 and 59 the three hottest curves share one ink run and are
# seeded to it deliberately; they are separable from 61 on.
FAN = {
    55.0: [3.79, 3.79, 3.79, 3.77, 3.74, 3.72, 3.59],
    57.0: [3.79, 3.79, 3.79, 3.77, 3.75, 3.64, 3.47],
    59.0: [3.84, 3.84, 3.84, 3.75, 3.71, 3.56, 3.34],
    61.0: [3.90, 3.89, 3.84, 3.71, 3.65, 3.48, 3.23],
    62.0: [3.93, 3.90, 3.83, 3.68, 3.61, 3.43, 3.19],
    63.0: [3.95, 3.91, 3.81, 3.65, 3.57, 3.39, 3.15],
    64.1: [3.98, 3.91, 3.80, 3.62, 3.53, 3.34, 3.10],
    68.1: [4.03, 3.89, 3.73, 3.51, 3.39, 3.22, 3.01],
    72.1: [4.01, 3.83, 3.64, 3.41, 3.29, 3.15, 2.99],
    76.1: [3.93, 3.73, 3.53, 3.34, 3.26, 3.15, 3.11],
    80.1: [3.83, 3.65, 3.47, 3.34, 3.32, 3.28, 3.23],
    82.0: [3.78, 3.62, 3.49, 3.40, 3.37, 3.36, 3.33],
}

# Knots through the bundle and the drops, where rank is not identity. Seeded from the
# printed digits; see the module docstring. Curve 1 carries its own dense set because its
# marker at 87.2 sits exactly on its knee, so the glyph centroid is not the line position.
BUNDLE = {
    1: [
        (83.5, 3.39),
        (85.0, 3.404),
        (86.5, 3.402),
        (87.0, 3.310),
        (88.0, 3.055),
        (88.5, 2.927),
        (89.0, 2.809),
        (90.0, 2.543),
        (90.5, 2.417),
        (91.5, 2.155),
        (92.0, 2.030),
        (93.0, 1.663),
        (93.5, 1.419),
        (94.0, 1.179),
        (94.5, 0.986),
        (103.8, 0.985),
    ],
    2: [
        (83.5, 3.42),
        (85.0, 3.47),
        (89.0, 3.56),
        (90.0, 3.55),
        (90.44, 3.54),
        (91.0, 3.39),
        (93.0, 2.90),
        (95.0, 2.39),
        (96.0, 2.15),
        (98.0, 1.69),
        (99.0, 1.47),
        (100.0, 1.25),
        (101.0, 1.03),
        (102.0, 0.98),
        (103.8, 0.98),
    ],
    3: [
        (83.5, 3.42),
        (85.0, 3.49),
        (87.2, 3.57),
        (89.0, 3.64),
        (91.0, 3.68),
        (93.0, 3.69),
        (94.69, 3.69),
        (96.0, 3.36),
        (98.0, 2.87),
        (99.0, 2.61),
        (100.0, 2.35),
        (101.0, 2.08),
        (102.0, 1.74),
        (103.0, 1.32),
        (103.8, 1.00),
    ],
    4: [
        (83.5, 3.44),
        (85.0, 3.49),
        (87.2, 3.58),
        (89.0, 3.68),
        (91.0, 3.72),
        (92.0, 3.76),
        (93.0, 3.78),
        (94.0, 3.78),
        (96.0, 3.78),
        (97.0, 3.77),
        (99.0, 3.22),
        (100.0, 2.98),
        (101.0, 2.72),
        (102.0, 2.47),
        (103.0, 2.22),
        (103.8, 1.99),
    ],
    5: [
        (83.5, 3.50),
        (85.0, 3.53),
        (87.2, 3.61),
        (89.0, 3.65),
        (91.0, 3.69),
        (93.0, 3.67),
        (94.0, 3.64),
        (96.0, 3.60),
        (98.0, 3.56),
        (99.0, 3.53),
        (100.0, 3.51),
        (101.0, 3.49),
        (102.0, 3.35),
        (103.0, 3.11),
        (103.8, 2.92),
    ],
    6: [
        (83.5, 3.61),
        (85.0, 3.61),
        (86.0, 3.62),
        (87.2, 3.64),
        (89.0, 3.68),
        (91.0, 3.72),
        (93.0, 3.63),
        (94.0, 3.60),
        (96.0, 3.56),
        (98.0, 3.51),
        (99.0, 3.49),
        (100.0, 3.46),
        (101.0, 3.44),
        (102.0, 3.30),
        (103.0, 3.04),
        (103.8, 2.84),
    ],
    7: [
        (83.5, 3.75),
        (85.0, 3.74),
        (87.2, 3.74),
        (89.0, 3.70),
        (91.0, 3.73),
        (93.0, 3.59),
        (94.0, 3.56),
        (96.0, 3.52),
        (98.0, 3.47),
        (99.0, 3.44),
        (100.0, 3.41),
        (101.0, 3.39),
        (102.0, 3.23),
        (103.0, 2.94),
        (103.8, 2.71),
    ],
}

SEED_TOL_WFPAC = 0.06
"""How far a declared seed may sit, in WFPAC, from the run it names. Wide enough for a
seed read off a plot by eye, tight enough that naming the wrong run fails. The curves are
never closer together than this except where they genuinely merge into one run."""

SEED_MARGIN_WFPAC = 0.10
"""How much better the named run must be than the next one, when the seed itself is not
within SEED_TOL_WFPAC. This is the test that actually matters: a seed must identify one
run unambiguously, even when its own value is imprecise."""

COVERAGE_FLOOR = 0.98
"""Mean fraction of segment samples that must lie within 6 px of ink. Measured 0.9850."""


def plot_mask(ink, box):
    """The plot interior with the in-plot legend block removed."""
    lo_x, hi_x, lo_y, hi_y = box
    work = ink.copy()
    work[int(lo_y) : int(hi_y), int(lo_x) : int(hi_x)] = False
    return work


def runs_in_column(work, col, top, bot):
    d = np.diff(np.concatenate(([0], work[top:bot, col].view(np.int8), [0])))
    a = np.flatnonzero(d == 1) + top
    b = np.flatnonzero(d == -1) + top - 1
    return list(zip(a.astype(float), b.astype(float), strict=True))


def main(argv=None) -> int:
    ink = native_bitmap(PAGE)
    fr = frame(ink, auto_frame_spec(PAGE))
    print(f"\n================ {CFG['fig']} {CFG['func']} (pdf p.{PAGE})")
    maj, _allt, _corners = ticks(ink, fr, CFG)

    TL = corner(fr["top"], fr["left"])
    TR = corner(fr["top"], fr["right"])
    BL = corner(fr["bottom"], fr["left"])
    BR = corner(fr["bottom"], fr["right"])
    H = homography([TL, TR, BR, BL], [(0, 1), (1, 1), (1, 0), (0, 0)])

    def uv_of(nm, pos):
        a, b = fr[nm]
        if nm in ("left", "right"):
            return _apply(H, [a + b * q for q in pos], pos)
        return _apply(H, pos, [a + b * q for q in pos])

    def paired(nm, vals, descending):
        pos = np.asarray(maj[nm], float)
        vv = np.asarray(vals, float)[::-1] if descending else np.asarray(vals, float)
        good = ~np.isnan(pos)
        pos, vv = pos[good], vv[good]
        o = np.argsort(pos)
        return pos[o], vv[o]

    xv, yv = major_values(CFG["x"]), major_values(CFG["y"])
    pb, xb = paired("bottom", xv, False)
    pt, xt = paired("top", xv, False)
    pl, yl = paired("left", yv, True)
    pr, yr = paired("right", yv, True)
    print("calibration:")
    CX = calibrate_axis(
        (uv_of("bottom", list(pb))[0], xb), (uv_of("top", list(pt))[0], xt), CFG["x"], "x"
    )
    CY = calibrate_axis(
        (uv_of("left", list(pl))[1], yl), (uv_of("right", list(pr))[1], yr), CFG["y"], "y"
    )

    # a straight pixel map, used only to turn a declared seed into a column and a row
    lef, rig = float(np.ravel(TL[0])[0]), float(np.ravel(TR[0])[0])
    top, bot = float(np.ravel(TL[1])[0]), float(np.ravel(BL[1])[0])
    to_col = lambda x: lef + (x - 50.0) / 60.0 * (rig - lef)  # noqa: E731
    to_row = lambda v: top + (5.0 - v) / 5.0 * (bot - top)  # noqa: E731

    work = plot_mask(ink, (to_col(51.0), to_col(82.0), to_row(1.7), to_row(0.2)))
    lo_r, hi_r = int(top) + 12, int(bot) - 12

    # column index <-> PCNGHL, through the calibrated axis rather than the raw frame
    cols = np.arange(int(lef) + 2, int(rig) - 1)
    uu, _ = _apply(H, cols.astype(float), np.full(cols.size, (top + bot) / 2.0))
    xcal = np.polyval(CX["c"], uu)

    def column_for(x):
        return int(cols[int(np.argmin(np.abs(xcal - x)))])

    def wfpac_of(col, row):
        u, v = _apply(H, [float(col)], [float(row)])
        return float(np.polyval(CY["c"], v[0]))

    seeds = {k: list(BUNDLE[k]) for k in T2_OF}
    for x, col_vals in FAN.items():
        for k in T2_OF:
            seeds[k].append((x, col_vals[7 - k]))
    for k in seeds:
        seeds[k] = sorted(set(seeds[k]))

    def clean_runs(c):
        """Runs in a column, keeping only those narrow enough to be line and not glyph."""
        return [(a, b) for a, b in runs_in_column(work, c, lo_r, hi_r) if (b - a) <= 8]

    def fit_segment(k, x0, v0, x1, v1):
        """Least-squares line through the glyph-free columns strictly inside a segment.

        The printed markers are the polyline's vertices, but a run centre measured *at* a
        marker is the glyph's centroid, not the line's position -- and the digits are not
        vertically symmetric, so that bias has a sign. Curve 7 reads 0.03 WFPAC high at
        every one of its markers. Fitting the straight part and intersecting adjacent fits
        recovers the vertex without ever measuring inside a glyph.
        """
        xs, ys = [], []
        lo, hi = column_for(x0 + 0.7), column_for(x1 - 0.7)
        for c in range(lo, hi + 1):
            rs = clean_runs(c)
            if not rs:
                continue
            t = float(np.polyval(CX["c"], _apply(H, [float(c)], [0.0])[0][0])) - x0
            t = 0.0 if x1 == x0 else t / (x1 - x0)
            want = v0 + max(0.0, min(1.0, t)) * (v1 - v0)
            cand = [((a + b) / 2.0, wfpac_of(c, (a + b) / 2.0)) for a, b in rs]
            centre, got = min(cand, key=lambda q: abs(q[1] - want))
            if abs(got - want) < 0.10:
                xs.append(float(c))
                ys.append(centre)
        if len(xs) < 4:
            return None
        return np.polyfit(np.array(xs), np.array(ys), 1)

    measured, worst_seed = {}, 0.0
    for k in sorted(seeds):
        pts = []
        for x, v in seeds[k]:
            c = column_for(x)
            rs = runs_in_column(work, c, lo_r, hi_r)
            if not rs:
                raise AssertionError(f"curve {k}: no ink in the column at PCNGHL {x}")
            cands = sorted(
                (((a + b) / 2.0, wfpac_of(c, (a + b) / 2.0)) for a, b in rs),
                key=lambda q: abs(q[1] - v),
            )
            centre, got = cands[0]
            miss = abs(got - v)
            worst_seed = max(worst_seed, miss)
            # The test is not "is the seed accurate" -- a seed read off a plot by eye is
            # not -- but "does it name one run and not another". On a steep segment a
            # 4 px column offset moves WFPAC by 0.07, so a fixed tolerance would reject a
            # perfectly good seed; what must hold is that the runner-up is clearly worse.
            runner_up = abs(cands[1][1] - v) if len(cands) > 1 else float("inf")
            if miss > SEED_TOL_WFPAC and runner_up - miss < SEED_MARGIN_WFPAC:
                raise AssertionError(
                    f"curve {k} at PCNGHL {x}: seed WFPAC {v} names a run at {got:.3f} "
                    f"({miss:.3f} away) with the next at {cands[1][1]:.3f} -- ambiguous. "
                    f"The seed is wrong or the figure has been re-rendered."
                )
            u, _vv = _apply(H, [float(c)], [centre])
            pts.append((float(np.polyval(CX["c"], u[0])), got, float(c), centre))
        pts = sorted(pts)
        # replace each interior vertex by the intersection of its two fitted segments
        sx = [q[0] for q in pts]
        sv = [q[1] for q in pts]
        fits = [fit_segment(k, sx[i], sv[i], sx[i + 1], sv[i + 1]) for i in range(len(pts) - 1)]
        out = list(pts)
        for i in range(1, len(pts) - 1):
            f0, f1 = fits[i - 1], fits[i]
            if f0 is None or f1 is None or abs(f0[0] - f1[0]) < 1e-9:
                continue
            cx_i = (f1[1] - f0[1]) / (f0[0] - f1[0])
            if not (pts[i][2] - 30.0 <= cx_i <= pts[i][2] + 30.0):
                continue
            cy_i = f0[0] * cx_i + f0[1]
            u, _ = _apply(H, [cx_i], [cy_i])
            out[i] = (float(np.polyval(CX["c"], u[0])), wfpac_of(cx_i, cy_i), cx_i, cy_i)
        measured[k] = out
    print(
        f"measured {sum(len(v) for v in measured.values())} knots; "
        f"worst seed-to-run miss {worst_seed:.3f} WFPAC of {SEED_TOL_WFPAC} allowed"
    )

    dist = ndimage.distance_transform_edt(~ink)
    cover = []
    for k in sorted(measured):
        p = measured[k]
        fr_k = []
        for i in range(len(p) - 1):
            # in raw pixels, from the measured run centres -- never round-tripped back
            # through a coordinate map, which is its own source of error
            xa, ya = p[i][2], p[i][3]
            xb, yb = p[i + 1][2], p[i + 1][3]
            n = max(12, int(abs(xb - xa)))
            t = np.linspace(0.0, 1.0, n)
            d = dist[
                np.round(ya + t * (yb - ya)).astype(int),
                np.round(xa + t * (xb - xa)).astype(int),
            ]
            fr_k.append(float((d < 6).mean()))
        cover.extend(fr_k)
        print(
            f"   curve {k} (T2={T2_OF[k]:5.1f}): {len(p):2d} knots, "
            f"segment-on-ink {np.mean(fr_k):.3f} mean, {min(fr_k):.2f} worst"
        )
    mean_cov = float(np.mean(cover))
    print(f"   overall {mean_cov:.4f} over {len(cover)} segments")
    if mean_cov < COVERAGE_FLOOR:
        raise AssertionError(
            f"segment-on-ink coverage {mean_cov:.4f} is below the {COVERAGE_FLOOR} floor; "
            f"the knot table no longer describes the printed curves"
        )

    rows = [(T2_OF[k], q[0], q[1]) for k in sorted(measured) for q in measured[k]]
    write_csv(rows, CX, CY, mean_cov, worst_seed)
    return 0


def write_csv(rows, CX, CY, cover, worst_seed):
    lines = [
        f"# source: TM-100991 pdf p.{PAGE}, Figure {CFG['fig']}",
        f"# quantity: {CFG['func']} -- {CFG['title']}",
        "# method: digitized -- the seven curves are PIECEWISE LINEAR between their printed",
        "#   markers, so this extracts their vertices rather than tracing their lines. Each",
        "#   value below is the centre of an ink run measured in the named column; the",
        "#   identity of the curve that run belongs to comes from rank in y where the fan is",
        "#   open and ordered by T2 (PCNGHL 55-82), and from the printed digit markers where",
        "#   it is not (the bundle at 85-95 and the drops beyond it). Rank is NOT identity",
        "#   through the bundle: the curves emerge from it ordered 4, 3, 5, 6, 7 with 1 and 2",
        "#   already descending, which is why tools/digitize_appc_multi.py refuses this page.",
        "# CHECK. Every segment between consecutive knots is sampled and required to lie on",
        f"#   ink: mean coverage {cover:.4f}; the worst seed named a run"
        f" {worst_seed:.3f} WFPAC away.",
        f"# {CFG['param'].lower()}_degR: {CFG['param']} curve parameter, {CFG['param_unit']}",
        f"# {CFG['xcol']}: {CFG['xdesc']}",
        f"# {CFG['ycol']}: {CFG['ydesc']}",
        f"# points: {len(rows)}",
        "# digitized: 2026-09-13 by tools/digitize_c30.py -- reruns and reproduces this file",
        "#   exactly",
        f"# RASTER. pdf p.{PAGE} is a single 300 dpi 1-bit CCITT image (2544x3300); this works",
        "#   on that bitmap directly and never resamples it.",
        "# CALIBRATION. Frame, tick lattice and axis polynomials come from",
        "#   tools/digitize_a6810.py unchanged, so this figure is calibrated exactly as the",
        "#   other fourteen are.",
        f"#   x: order {CX['order']}, interior-tick residual {CX['resid']:.4g}"
        f" over {CX['n']} ticks",
        f"#   y: order {CY['order']}, interior-tick residual {CY['resid']:.4g}"
        f" over {CY['n']} ticks",
        "# NOT transcribed. Values carry read error.",
        f"{CFG['param'].lower()}_degR,{CFG['xcol']},{CFG['ycol']}",
    ]
    for t2, x, v in rows:
        lines.append(f"{t2:.1f},{x:.4f},{v:.4f}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT} ({len(rows)} points)")


if __name__ == "__main__":
    raise SystemExit(main())
