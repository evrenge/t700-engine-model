#!/usr/bin/env python3
"""Figure C30 (pdf p.101) -- `F_HM7`, the HMU maximum fuel-parameter limit during acceleration.

The one Appendix C scheduling function with its own tool (open question #50).
`digitize_multi_curve` identifies a curve by rank in y, but C30's seven curves cross: they
converge into a bundle around PCNGHL 85-95 and come out reordered (4, 3, then 5, 6, 7, with
1 and 2 already descending). Reading the digits and tracing through the crossings both
fail; the record is in #50.

What works: the curves are exact polylines between their markers (on curve 4 a straight
line between markers predicts the ink at PCNGHL 99-103 to 0.02 WFPAC), so the job is to find
each curve's vertices. `FAN` and `BUNDLE` declare, per curve, the columns where it is
identifiable and a seed for its value. The tool then measures: at each declared column it
takes the ink run nearest the seed and reports its centre through the calibrated axes. The
seeds only name which run; a seed that does not name one run unambiguously fails loudly.

Where the seeds come from [pdf p.101]:

* PCNGHL 55-82: the fan is open and ordered by T2, so rank is identity -- checked at
  PCNGHL 68.1, where the seven glyphs share one abscissa and their WFPAC values are exactly
  monotone in T2. Columns 61-63 each resolve seven clean runs: no inversion at the left end.
* The bundle and the drops: the printed digits, read off the 300 dpi raster at 600 %.
  Curve 1's "1" at (87.2, 3.41) and (94.5, 0.99); curve 2's "2" at (90.4, 3.54), (92.5,
  3.02) and (96.7, 1.99); curve 3's "3" at (94.7, 3.69); curve 4's "4" flat at 3.78 across
  92.5, 94.8 and 96.8; the three surviving curves labelled 5, 6, 7 top to bottom at 103.8.

Check: every segment between consecutive knots is sampled and must lie on ink; `main`
fails below `COVERAGE_FLOOR`.

    python3 tools/digitize_c30.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))

from figure.axis import blend_value, calibrate_under
from figure.frame import auto_frame_spec, runs
from figure.raster import native_bitmap

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
        # confirmed independently: 1.488 against a 1.486 midpoint of its two neighbours
        (93.5, 1.488),
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
        # confirmed independently: 3.005 against a 3.023 midpoint of its two neighbours
        (103.0, 3.005),
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
"""Mean fraction of segment samples that must lie within 6 px of ink.

The fine-lattice calibration, which several-fold improves the held-out frame error, puts
the measured coverage exactly on this floor (0.9800). The frame test is the stronger
evidence, but if this check ever fails, look here first."""


def plot_mask(ink, box):
    """The plot interior with the in-plot legend block removed."""
    lo_x, hi_x, lo_y, hi_y = box
    work = ink.copy()
    work[int(lo_y) : int(hi_y), int(lo_x) : int(hi_x)] = False
    return work


def runs_in_column(work, col, top, bot):
    """`(first, last)` row of every ink run in one column between `top` and `bot`."""
    return [(float(a + top), float(b + top - 1)) for a, b in runs(work[top:bot, col])]


def main(argv=None) -> int:
    ink = native_bitmap(PAGE)
    print(f"\n================ {CFG['fig']} {CFG['func']} (pdf p.{PAGE})")
    # the single-curve calibration path, with straight frame edges
    cal = calibrate_under(ink, auto_frame_spec(PAGE), CFG, 1)
    TL, TR, _BR, BL = cal["corners"]
    to_uv, CX, CY = cal["to_uv"], cal["CX"], cal["CY"]

    # a straight pixel map, used only to turn a declared seed into a column and a row
    lef, rig = float(np.ravel(TL[0])[0]), float(np.ravel(TR[0])[0])
    top, bot = float(np.ravel(TL[1])[0]), float(np.ravel(BL[1])[0])
    to_col = lambda x: lef + (x - 50.0) / 60.0 * (rig - lef)  # noqa: E731
    to_row = lambda v: top + (5.0 - v) / 5.0 * (bot - top)  # noqa: E731

    work = plot_mask(ink, (to_col(51.0), to_col(82.0), to_row(1.7), to_row(0.2)))
    lo_r, hi_r = int(top) + 12, int(bot) - 12

    # Column index <-> PCNGHL, through the calibrated axis. The bottom and top x scales
    # disagree by about five columns, so the blend weight comes from the wanted WFPAC
    # (inverting the y calibration), never from a row: on this tilted frame a row guessed by
    # `to_row` is wrong by the full tilt at the right-hand end.
    cols = np.arange(int(lef) + 2, int(rig) - 1)
    _tgrid = np.linspace(0.0, 1.0, 2001)
    _vgrid = np.asarray([blend_value(CY["edgec"], t, 0.5) for t in _tgrid], float)
    _uu, _ = to_uv(cols.astype(float), np.full(cols.size, (top + bot) / 2.0))

    def column_for(x, v):
        across = float(np.interp(v, _vgrid, _tgrid))
        return int(cols[int(np.argmin(np.abs(blend_value(CX["edgec"], _uu, across) - x)))])

    def wfpac_of(col, row):
        u, v = to_uv([float(col)], [float(row)])
        return float(blend_value(CY["edgec"], v[0], u[0]))

    seeds = {k: list(BUNDLE[k]) for k in T2_OF}
    for x, col_vals in FAN.items():
        for k in T2_OF:
            seeds[k].append((x, col_vals[7 - k]))
    for k in seeds:
        seeds[k] = sorted(set(seeds[k]))

    def clean_runs(c):
        """Runs in a column, keeping only those narrow enough to be line and not glyph."""
        return [(a, b) for a, b in runs_in_column(work, c, lo_r, hi_r) if (b - a) <= 8]

    def fit_segment(x0, v0, x1, v1):
        """Least-squares line through the glyph-free columns strictly inside a segment.

        The printed markers are the polyline's vertices, but a run centre measured *at* a
        marker is the glyph's centroid, not the line's position -- and the digits are not
        vertically symmetric, so that bias has a sign. Curve 7 reads 0.03 WFPAC high at
        every one of its markers. Fitting the straight part and intersecting adjacent fits
        recovers the vertex without ever measuring inside a glyph.
        """
        xs, ys = [], []
        lo, hi = column_for(x0 + 0.7, v0), column_for(x1 - 0.7, v1)
        for c in range(lo, hi + 1):
            rs = clean_runs(c)
            if not rs:
                continue
            uc, vc = to_uv([float(c)], [float(sum(rs[0]) / 2.0)])
            t = float(blend_value(CX["edgec"], uc[0], vc[0])) - x0
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
            c = column_for(x, v)
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
            # The test is whether the seed names one run, not whether it is accurate: on a
            # steep segment a 4 px column offset moves WFPAC by 0.07, so instead the
            # runner-up must be clearly worse.
            runner_up = abs(cands[1][1] - v) if len(cands) > 1 else float("inf")
            if miss > SEED_TOL_WFPAC and runner_up - miss < SEED_MARGIN_WFPAC:
                raise AssertionError(
                    f"curve {k} at PCNGHL {x}: seed WFPAC {v} names a run at {got:.3f} "
                    f"({miss:.3f} away) with the next at {cands[1][1]:.3f} -- ambiguous. "
                    f"The seed is wrong or the figure has been re-rendered."
                )
            u, vpt = to_uv([float(c)], [centre])
            pts.append((float(blend_value(CX["edgec"], u[0], vpt[0])), got, float(c), centre))
        pts = sorted(pts)
        # replace each interior vertex by the intersection of its two fitted segments
        sx = [q[0] for q in pts]
        sv = [q[1] for q in pts]
        fits = [fit_segment(sx[i], sv[i], sx[i + 1], sv[i + 1]) for i in range(len(pts) - 1)]
        out = list(pts)
        for i in range(1, len(pts) - 1):
            f0, f1 = fits[i - 1], fits[i]
            if f0 is None or f1 is None or abs(f0[0] - f1[0]) < 1e-9:
                continue
            cx_i = (f1[1] - f0[1]) / (f0[0] - f1[0])
            if not (pts[i][2] - 30.0 <= cx_i <= pts[i][2] + 30.0):
                continue
            cy_i = f0[0] * cx_i + f0[1]
            u, vpt = to_uv([cx_i], [cy_i])
            out[i] = (
                float(blend_value(CX["edgec"], u[0], vpt[0])),
                wfpac_of(cx_i, cy_i),
                cx_i,
                cy_i,
            )
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
        "#   already descending, which is why tools/digitize_multi_curve.py refuses this page.",
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
        "#   tools/digitize_single_curve.py unchanged, so this figure is calibrated exactly as the",
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
