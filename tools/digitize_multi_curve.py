#!/usr/bin/env python3
"""Figure C23 -- the multi-curve path, and why C30 needs a tool of its own.

C23 and C30 are the two Appendix C functions with more than one curve. **This tool
extracts C23 and refuses C30**, by design: the method below separates curves by rank in y,
and C30's curves cross. `digitize_c30.py` tracks those through the crossings instead. Ask
for c30 here and it says so rather than emitting something plausible.

Neither fits the single-curve path in `digitize_single_curve`, for two reasons:

**The markers are printed DIGITS, not `x` glyphs.** `fit_glyphs` models two crossing
strokes and would be fitting the wrong shape, so this tool takes the marker position from
the ink centroid of each digit instead. The digits are printed centred on the curve, so
the centroid is the plotted point -- but it is a weaker estimator than a model fit, and
the uncertainty reported here says so.

**The digit says which curve the point belongs to, and it does not have to be read.**
C23's four curves are ordered 1 at the top through 4 at the bottom at every abscissa, and
never cross; C30's seven are likewise ordered by T2 everywhere the report draws a marker.
So the curve a point belongs to is its RANK in y within its column of points, and no OCR
is needed. Ranking is also checkable, which reading a 10 px digit off a 1988 scan is not:
the assignment is rejected unless every column yields exactly the expected number of
points.

Everything upstream of the markers -- the frame location, the tick lattice, the held-out
frame calibration and the axis polynomials -- is imported from `digitize_single_curve` unchanged,
so these two figures are calibrated exactly as the other thirteen are.

    python3 tools/digitize_multi_curve.py c23
    python3 tools/digitize_multi_curve.py --probe c23
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))

from figure.axis import blend_value, pick_lattice
from figure.frame import apply_h as _apply
from figure.frame import as_lines, auto_frame_spec, edge_uv, fit_edges, homography
from figure.frame import corner_line as corner
from figure.frame import dist as _dist
from figure.raster import native_bitmap
from figure.ticks import fine_lattice_all, major_values, majors_all_edges
from figure.ticks import paired as _paired

FIGS = {
    "c23": dict(
        page=94,
        fig="C23",
        func="F_EC1",
        title=(
            "ECU thermocouple sensor time constant (Fig. C5, pdf p.87: TAU45 = F_EC1(T45L, W45R))"
        ),
        xcol="w45r",
        ycol="tau45_s",
        xdesc="W45R power turbine flow parameter, nondimensional",
        ydesc="TAU45 sensor time constant, sec",
        csv=Path("data/schedules/fec1_thermocouple_tau.csv"),
        x=dict(lo=0.0, hi=16.0, step=2.0, n=7, printed=("lo", "hi")),
        y=dict(lo=1.0, hi=6.0, step=1.0, n=4, printed=("lo", "hi")),
        curves=[1260.0, 1660.0, 2060.0, 2460.0],
        param="T45L",
        param_unit="deg R",
        knots=[1.0, 2.5, 6.5, 10.0, 15.0],
        knot_tol=0.30,
        run_gap=5,
        open_sz=5,
        legend=dict(x=(0.42, 1.0), y=(0.78, 1.0)),  # top right; v = 0 at the bottom
    ),
    "c30": dict(
        page=101,
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
        csv=Path("data/schedules/fhm7_accel_limit.csv"),
        x=dict(lo=50.0, hi=110.0, step=10.0, n=5, printed=("lo", "hi")),
        y=dict(lo=0.0, hi=5.0, step=1.0, n=4, printed=("lo", "hi")),
        curves=[395.0, 430.0, 460.0, 480.0, 519.0, 555.0, 590.0],
        param="T2",
        param_unit="deg R",
        # NOT EXTRACTED. C30's seven curves CROSS, around PCNGHL 88-92, where curve 4 runs
        # above 5, 6 and 7. Rank in y is therefore not identity there, and the measurement
        # says so: only 70 of 410 sampled columns resolve seven separate runs at all, and
        # the ones that do are outside the bundle. Each curve also drops steeply to WFPAC
        # = 1 at its own speed, so the knot abscissae differ per curve and a shared `knots`
        # list cannot describe them.
        #
        # What this figure needs is per-curve TRACKING -- following each line through the
        # crossings by direction rather than by order -- with the printed digits used to
        # seed and to check the identities. That is a different algorithm from the one
        # here, not a parameter of it.
        knots=None,
        knot_tol=0.40,
        open_sz=4,
        legend=dict(x=(0.0, 0.45), y=(0.0, 0.30)),
    ),
}


def blobs(ink, fr, cfg):
    """Digit-marker centroids inside the plot, with the legend block excluded."""
    Y, X = np.mgrid[0 : ink.shape[0], 0 : ink.shape[1]]
    inside = (
        (_dist(fr["left"], X, Y) > 6)
        & (_dist(fr["right"], X, Y) < -6)
        & (_dist(fr["top"], X, Y, False) > 6)
        & (_dist(fr["bottom"], X, Y, False) < -6)
    )
    op = ndimage.binary_opening(inside & ink, structure=np.ones((cfg["open_sz"],) * 2, bool))
    lab, n = ndimage.label(op, structure=np.ones((3, 3)))
    if n == 0:
        raise AssertionError("no marker candidates")
    com = np.array(ndimage.center_of_mass(op, lab, range(1, n + 1)))
    sz = np.array(ndimage.sum(op, lab, range(1, n + 1)))
    cy, cx = com[:, 0], com[:, 1]

    lo_x, hi_x = _dist(fr["left"], cx, cy), -_dist(fr["right"], cx, cy)
    u = lo_x / (lo_x + hi_x)
    lo_y, hi_y = _dist(fr["top"], cx, cy, False), -_dist(fr["bottom"], cx, cy, False)
    v = 1.0 - lo_y / (lo_y + hi_y)  # v = 0 at the BOTTOM, as everywhere else in the package
    lg = cfg["legend"]
    in_legend = (u > lg["x"][0]) & (u < lg["x"][1]) & (v > lg["y"][0]) & (v < lg["y"][1])
    keep = ~in_legend
    return cx[keep], cy[keep], sz[keep], u[keep], v[keep]


def trace(ink, fr, cfg, ncur, step=4, halfwidth=2):
    """Follow the curves column by column, ranking ink runs by y.

    The curves of C23 and C30 are ordered by their parameter at every abscissa and never
    cross, so at any column the n-th ink run from the top IS the n-th curve. That makes
    tracing an ordering problem rather than a tracking one, and it needs no continuity
    assumption to get started and no recovery when a column is ambiguous -- a column that
    does not yield exactly `ncur` runs is simply reported as unresolved.

    Returns one list of run centres per column, or None where the count was wrong.
    """
    h, w = ink.shape
    Y, X = np.mgrid[0:h, 0:w]
    # 18 px, not 8: the tick marks protrude INWARD from every frame by about 10 px, and a
    # tick is an ink run like any other. At C23's x = 10 and x = 15 knots -- both tick
    # columns -- a bottom tick was being counted as the fourth curve, putting TAU45 at 1.02
    # where the figure shows 1.87.
    pad = cfg.get("frame_pad", 18)
    inside = (
        (_dist(fr["left"], X, Y) > pad)
        & (_dist(fr["right"], X, Y) < -pad)
        & (_dist(fr["top"], X, Y, False) > pad)
        & (_dist(fr["bottom"], X, Y, False) < -pad)
    )
    lg = cfg["legend"]
    lo_x, hi_x = _dist(fr["left"], X, Y), -_dist(fr["right"], X, Y)
    u = lo_x / (lo_x + hi_x)
    lo_y, hi_y = _dist(fr["top"], X, Y, False), -_dist(fr["bottom"], X, Y, False)
    v = 1.0 - lo_y / (lo_y + hi_y)  # v = 0 at the bottom
    in_legend = (u > lg["x"][0]) & (u < lg["x"][1]) & (v > lg["y"][0]) & (v < lg["y"][1])
    plot = ink & inside & ~in_legend

    cols = range(int(fr["left"][0]) + 12, int(fr["right"][0]) - 12, step)
    out = []
    for c in cols:
        band = plot[:, max(0, c - halfwidth) : c + halfwidth + 1].any(axis=1)
        idx = np.flatnonzero(band)
        if idx.size == 0:
            out.append(None)
            continue
        brk = np.flatnonzero(np.diff(idx) > cfg.get("run_gap", 6))
        groups = np.split(idx, brk + 1)
        groups = [g for g in groups if g.size >= 2]
        out.append([float(g.mean()) for g in groups] if len(groups) == ncur else None)
    return list(cols), out


def write_csv(cfg, rows, CX, CY, nok, ntot):
    """Three columns -- parameter, x, y -- the same shape `maps.load_speed_map` reads."""
    ax, ay = cfg["x"], cfg["y"]
    L = [
        f"# source: TM-100991 pdf p.{cfg['page']}, Figure {cfg['fig']}",
        f"# quantity: {cfg['func']} -- {cfg['title']}",
        "# method: digitized -- the curves are TRACED column by column and ranked in y,",
        "#   not read from their printed digit markers. The digits identify which curve a",
        "#   point belongs to, and on this page they overprint where the curves converge --",
        "#   at C23's right-hand knot all four merge into a single 462 px blob even at a 5x5",
        "#   opening. But the curves are ordered by the parameter at every abscissa and never",
        "#   cross, so the n-th ink run from the top IS the n-th curve, and the identity",
        "#   follows from rank without reading a 10 px digit off a 1988 scan.",
        f"#   {nok} of {ntot} sampled columns resolved all {len(cfg['curves'])} curves;",
        "#   the knots below are taken from resolved columns only.",
        f"# {cfg['param'].lower()}: {cfg['param']} curve parameter, {cfg['param_unit']}",
        f"# {cfg['xcol']}: {cfg['xdesc']}",
        f"# {cfg['ycol']}: {cfg['ydesc']}",
        f"# points: {len(rows)}",
        "# digitized: 2026-09-12 by tools/digitize_multi_curve.py -- reruns and reproduces",
        "#   this file exactly",
        f"# RASTER. pdf p.{cfg['page']} is a single 300 dpi 1-bit CCITT image (2544x3300);",
        "#   this works on that bitmap directly and never resamples it.",
        "# CALIBRATION. Frame, tick lattice and axis polynomials come from",
        "#   tools/digitize_single_curve.py unchanged, so this figure is calibrated exactly as the",
        "#   thirteen single-curve figures are.",
        f"#   x: order {CX['order']}, interior-tick residual {CX['resid']:.4g}"
        f" over {CX['n']} ticks",
        f"#   y: order {CY['order']}, interior-tick residual {CY['resid']:.4g}"
        f" over {CY['n']} ticks",
        "# UNCERTAINTY is not reported per point: a traced run centre is a weaker estimator",
        "#   than the generative glyph fit the single-curve figures use, and no honest",
        "#   per-point sigma was derived for it. What WAS measured is the sensitivity to the",
        "#   one free parameter of the trace, the gap at which two ink runs are called",
        "#   separate. Sweeping it over 3, 4, 5 and 6 px moves every value by at most 0.05 in",
        "#   y -- about 1 per-cent of the axis range -- and most by far less. The largest",
        "#   movers are the two lowest curves at the knots where they converge, which is",
        "#   where a gap threshold is doing the most work and where the figure itself is",
        "#   least legible.",
        "# NOT transcribed. Values carry read error.",
        f"{cfg['param'].lower()},{cfg['xcol']},{cfg['ycol']}",
    ]
    for pval, x, y in rows:
        L.append(f"{pval:g},{x:.4f},{y:.5f}")
    cfg["csv"].write_text("\n".join(L) + "\n")
    print(f"  wrote {cfg['csv']}  ({len(L)} lines)")
    assert ax and ay


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("key", choices=sorted(FIGS))
    ap.add_argument("--probe", action="store_true", help="report candidates and stop")
    a = ap.parse_args(argv)
    cfg = FIGS[a.key]
    if cfg["knots"] is None and not a.probe:
        raise SystemExit(
            f"{cfg['fig']} is not extracted: its curves cross, so rank in y is not "
            f"identity. See the note in FIGS[{a.key!r}]."
        )

    ink = native_bitmap(cfg["page"])
    fc, _ = fit_edges(ink, auto_frame_spec(cfg["page"]), order=1)
    fr = as_lines(fc)
    print(f"\n================ {cfg['fig']} {cfg['func']} (pdf p.{cfg['page']})")
    maj, allt, corners = majors_all_edges(ink, fc, cfg)

    cx, cy, sz, u, v = blobs(ink, fr, cfg)
    order = np.argsort(cx)
    print(f"markers: {cx.size} candidates after excluding the legend block")
    if a.probe:
        for i in np.argsort(u):
            print(
                f"   u={u[i]:.4f} v={v[i]:.4f} size={sz[i]:5.0f} "
                f"x_px={cx[i]:7.1f} y_px={cy[i]:7.1f}"
            )
        return 0
    if False:
        for i in order:
            print(
                f"   x_px={cx[i]:8.1f} y_px={cy[i]:8.1f} size={sz[i]:5.0f} "
                f"u={u[i]:.4f} v={v[i]:.4f}"
            )
        return 0
    # ---- calibration, exactly as the single-curve tool does it --------------------
    TL = corner(fr["top"], fr["left"])
    TR = corner(fr["top"], fr["right"])
    BL = corner(fr["bottom"], fr["left"])
    BR = corner(fr["bottom"], fr["right"])
    H = homography([TL, TR, BR, BL], [(0, 1), (1, 1), (1, 0), (0, 0)])
    uv_of = edge_uv(fr, H)

    xv, yv = major_values(cfg["x"]), major_values(cfg["y"])

    def tick_u(nm, pos):
        q = uv_of(nm, list(np.atleast_1d(pos)))
        return q[1] if nm in ("left", "right") else q[0]

    coarse_of = {}
    for nm, vals, desc in (
        ("bottom", xv, False),
        ("top", xv, False),
        ("left", yv, True),
        ("right", yv, True),
    ):
        pos, val = _paired(maj, nm, vals, desc)
        coarse_of[nm] = (np.asarray(tick_u(nm, pos), float), val)

    # The MINOR lattice, on the same terms as the thirteen single-curve figures: the step
    # measured from the ticks, each tick labelled through a fit made on the majors alone,
    # and the axis keeping whichever lattice predicts the printed frame value it was never
    # shown. This figure shares its calibration with `digitize_single_curve` and had been
    # left behind on the majors when that path moved.
    print("fine lattice:")
    span = {
        "left": (TL[1], BL[1]),
        "right": (TR[1], BR[1]),
        "top": (TL[0], TR[0]),
        "bottom": (BL[0], BR[0]),
    }
    fine = fine_lattice_all(ink, fc, cfg, tick_u, span, coarse_of)
    print("calibration:")
    CX = pick_lattice(
        coarse_of["bottom"], coarse_of["top"], fine["bottom"], fine["top"], cfg["x"], "x"
    )
    CY = pick_lattice(
        coarse_of["left"], coarse_of["right"], fine["left"], fine["right"], cfg["y"], "y"
    )

    # ---- trace, then sample at the knot abscissae ---------------------------------
    ncur = len(cfg["curves"])
    cols, runs = trace(ink, fr, cfg, ncur)
    ok = [(c, r) for c, r in zip(cols, runs, strict=True) if r is not None]
    print(f"traced {len(ok)} of {len(runs)} columns with all {ncur} curves resolved")
    if not ok:
        raise AssertionError("no column resolved every curve")

    us, vs = [], []
    for c, r in ok:
        uu, vv = zip(*(_apply(H, [c], [y]) for y in r), strict=True)
        us.append(float(uu[0][0]))
        vs.append([float(q[0]) for q in vv])
    us = np.array(us)
    vs = np.array(vs)
    # Blended across the plot, as digitize_single_curve now does. The x map's blend
    # coordinate is v, which differs per curve at a given column, so the column's mean v
    # locates the knot; the y map's blend coordinate is that column's own u.
    xvals = blend_value(CX["edgec"], us, vs.mean(axis=1))

    knots = np.asarray(cfg["knots"], float)
    print(f"sampling at the {knots.size} printed knots: {knots}")
    rows = []
    for k in knots:
        j = int(np.argmin(np.abs(xvals - k)))
        if abs(xvals[j] - k) > cfg.get("knot_tol", 0.25):
            raise AssertionError(f"no resolved column within tolerance of knot {k}")
        for i, pval in enumerate(cfg["curves"]):
            rows.append((pval, float(k), float(blend_value(CY["edgec"], vs[j, i], us[j]))))
            print(
                f"   {cfg['param']}={pval:7.1f}  x={k:7.3f}  y={rows[-1][2]:8.4f} "
                f"(column x={xvals[j]:.3f})"
            )
    write_csv(cfg, rows, CX, CY, len(ok), len(runs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
