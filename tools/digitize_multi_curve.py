#!/usr/bin/env python3
"""Figure C23 (`F_EC1`, pdf p.94) -- the multi-curve path, and why C30 has its own tool.

C23's four curves are ordered 1 at the top through 4 at the bottom at every abscissa and
never cross, so the curve a point belongs to is its RANK in y within its column. The curves
are traced column by column and sampled at the printed knots; the printed digit markers
are never read (they overprint where the curves converge). A column is used only if it
yields exactly one ink run per curve.

C30's curves cross, so rank is not identity there: `digitize_c30.py` handles it, and c30
is refused here as an unknown figure rather than emitting something plausible.

Calibration (frame, tick lattice, held-out axis polynomials) is
`figure.axis.calibrate_under`, the path `digitize_single_curve` uses.

    python3 tools/digitize_multi_curve.py c23
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from figure.axis import blend_value, calibrate_under
from figure.frame import auto_frame_spec
from figure.frame import dist as _dist
from figure.raster import native_bitmap

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
        legend=dict(x=(0.42, 1.0), y=(0.78, 1.0)),  # top right; v = 0 at the bottom
    ),
}


def trace(ink, fr, cfg, ncur, step=4, halfwidth=2):
    """Follow the curves column by column, ranking ink runs by y.

    C23's curves never cross, so at any column the n-th ink run from the top IS the n-th
    curve: an ordering problem, not a tracking one, with no continuity assumption.

    Returns the sampled columns and, per column, the run centres, or None where the count
    was not `ncur`.
    """
    h, w = ink.shape
    Y, X = np.mgrid[0:h, 0:w]
    # 18 px clears the ticks, which protrude ~10 px inward and are ink runs like any other
    # (at the x = 10 and 15 knots a bottom tick would count as the fourth curve).
    pad = 18
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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("key", choices=sorted(FIGS))
    a = ap.parse_args(argv)
    cfg = FIGS[a.key]

    ink = native_bitmap(cfg["page"])
    print(f"\n================ {cfg['fig']} {cfg['func']} (pdf p.{cfg['page']})")
    # ---- calibration: the single-curve path, with straight frame edges ------------
    cal = calibrate_under(ink, auto_frame_spec(cfg["page"]), cfg, 1)
    fr, to_uv, CX, CY = cal["fr"], cal["to_uv"], cal["CX"], cal["CY"]

    # ---- trace, then sample at the knot abscissae ---------------------------------
    ncur = len(cfg["curves"])
    cols, runs = trace(ink, fr, cfg, ncur)
    ok = [(c, r) for c, r in zip(cols, runs, strict=True) if r is not None]
    print(f"traced {len(ok)} of {len(runs)} columns with all {ncur} curves resolved")
    if not ok:
        raise AssertionError("no column resolved every curve")

    us, vs = [], []
    for c, r in ok:
        uu, vv = zip(*(to_uv([c], [y]) for y in r), strict=True)
        us.append(float(uu[0][0]))
        vs.append([float(q[0]) for q in vv])
    us = np.array(us)
    vs = np.array(vs)
    # Blended across the plot. The x map's blend coordinate is v, which differs per curve
    # at a given column, so the column's mean v locates the knot; the y map's blend
    # coordinate is that column's own u.
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
