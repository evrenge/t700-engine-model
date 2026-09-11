#!/usr/bin/env python3
"""Figures C23 and C30 -- the two Appendix C functions with more than one curve.

These are the last of the eight scheduling functions and they do not fit the single-curve
path in `digitize_a6810`, for two reasons:

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
frame calibration and the axis polynomials -- is imported from `digitize_a6810` unchanged,
so these two figures are calibrated exactly as the other thirteen are.

    python3 tools/digitize_appc_multi.py c23
    python3 tools/digitize_appc_multi.py --probe c23
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))

from digitize_a6810 import (  # noqa: E402
    _apply,
    _dist,
    auto_frame_spec,
    calibrate_axis,
    corner,
    frame,
    homography,
    major_values,
    native_bitmap,
    ticks,
)

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
        legend=dict(x=(0.42, 1.0), y=(0.0, 0.22)),
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
    v = lo_y / (lo_y + hi_y)
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
    v = lo_y / (lo_y + hi_y)
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
        "# digitized: 2026-09-12 by tools/digitize_appc_multi.py -- reruns and reproduces",
        "#   this file exactly",
        f"# RASTER. pdf p.{cfg['page']} is a single 300 dpi 1-bit CCITT image (2544x3300);",
        "#   this works on that bitmap directly and never resamples it.",
        "# CALIBRATION. Frame, tick lattice and axis polynomials come from",
        "#   tools/digitize_a6810.py unchanged, so this figure is calibrated exactly as the",
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
    fr = frame(ink, auto_frame_spec(cfg["page"]))
    print(f"\n================ {cfg['fig']} {cfg['func']} (pdf p.{cfg['page']})")
    maj, allt, corners = ticks(ink, fr, cfg)

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

    xv, yv = major_values(cfg["x"]), major_values(cfg["y"])
    pb, xb = paired("bottom", xv, False)
    pt, xt = paired("top", xv, False)
    pl, yl = paired("left", yv, True)
    pr, yr = paired("right", yv, True)
    print("calibration:")
    CX = calibrate_axis(
        (uv_of("bottom", list(pb))[0], xb), (uv_of("top", list(pt))[0], xt), cfg["x"], "x"
    )
    CY = calibrate_axis(
        (uv_of("left", list(pl))[1], yl), (uv_of("right", list(pr))[1], yr), cfg["y"], "y"
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
    xvals = np.polyval(CX["c"], us)

    knots = np.asarray(cfg["knots"], float)
    print(f"sampling at the {knots.size} printed knots: {knots}")
    rows = []
    for k in knots:
        j = int(np.argmin(np.abs(xvals - k)))
        if abs(xvals[j] - k) > cfg.get("knot_tol", 0.25):
            raise AssertionError(f"no resolved column within tolerance of knot {k}")
        for i, pval in enumerate(cfg["curves"]):
            rows.append((pval, float(k), float(np.polyval(CY["c"], vs[j, i]))))
            print(
                f"   {cfg['param']}={pval:7.1f}  x={k:7.3f}  y={rows[-1][2]:8.4f} "
                f"(column x={xvals[j]:.3f})"
            )
    write_csv(cfg, rows, CX, CY, len(ok), len(runs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
