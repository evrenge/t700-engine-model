"""Digitize Figure A1 (pdf p.56) -- f1, compressor mass flow. The one 2-D map in the report.

Figure A1 is the only figure in TM-100991 the single-curve path cannot do, so it gets its
own script. Four things make it different:

1. **Eleven speed lines that cross and crowd.** A per-curve `--roi` band cannot separate
   curves 1-4 in the lower left. It turns out not to be needed: at 600 dpi each speed line
   -- flat extension, joining polyline and all seven digit markers -- is a single
   8-connected ink component, and the eleven components are disjoint. Connected-component
   labelling separates the family exactly.

2. **The marker is the curve's own index digit, not a cross.** `1`..`11`, so the marker is
   one glyph on nine curves and a two-glyph string on two, and a digit's ink centroid is
   not its string centre. Each curve gets the estimator its digit shape allows: density
   peak for the single digits, the bowl of the `0` for SYMBOL 10, the pair of `1` strokes
   for SYMBOL 11, and vertical strokes alone for SYMBOL 1, whose markers overprint.

3. **The figure carries its own check.** Six dotted construction lines each join the k-th
   marker of all eleven speed lines, and a seventh vertical one joins the left ends. They
   are printed ink, they are separable from the data (small isolated components), and they
   settle the attribution question that matters on a 2-D map: which speed line and which
   index does this marker belong to. Every extracted point is tested against them.

4. **The axes are calibrated on the SCAN, not on a render, and not linearly.** `pdfimages
   -list` shows p.56 is one 300 dpi 1-bit CCITT image, so the 600 dpi page is a 2x upsample
   -- fine for sub-pixel centroids, wrong for geometry. The four frame lines are fitted
   sub-pixel in the native bitmap and their rotation and keystone are carried as a
   HOMOGRAPHY applied to the measured points; nothing interpolates the geometry. The tick
   lattice then bows along the page, so the value map is a CUBIC in the homography
   coordinate, chosen by a held-out test the figure gives away free: both axes print a
   value at both frame lines (0/20 in x, 2/12 in y), so a fit made on the interior ticks
   alone can be asked to predict four numbers it was never told. See `calibrate`.

Run from the repo root; writes data/maps/f1_compressor_mass_flow.csv and its beta-gridded
derivative. Inputs:

    pdfimages -f 56 -l 56 -png docs/ballin-tm100991.pdf validation/out/digitize/a1redo/native
    pdftoppm -f 56 -l 56 -r 600 -png docs/ballin-tm100991.pdf validation/out/digitize/a1-600dpi
    rectify_page.rectify(<that render>, validation/out/digitize/a1r_600-056.png)   # no CLI

All three are pulled automatically if missing -- see `native_bitmap` and `ensure_page` --
so this script runs from a fresh clone with only the PDF and poppler present. Until
2026-09-13 only the first was, and the other two were gitignored intermediates no tool
produced, which made `tools/reproduce_all.sh` unrunnable outside the one machine that
happened to still have them.

This module lives in tools/ and may use SciPy and PIL. Nothing here is imported by
src/t700/ -- and nothing imports *it*: it is a script that works at module level, and
importing it raises.
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
from scipy.spatial import cKDTree

import rectify_page
from figure.axis import blend_value as blend
from figure.axis import fit_soft_anchored as _anchored
from figure.frame import apply_h as _apply
from figure.frame import corner_line as corner
from figure.frame import homography, robust_keep
from figure.frame import trace as _trace
from figure.raster import load_gray, native_bitmap, render
from figure.ticks import groups as _tick_groups
from figure.ticks import protrusion as _protrusion

# This file is a script whose work happens at module level, so *importing* it used to run
# the whole digitization and overwrite data/maps/f1_compressor_mass_flow.csv as a side
# effect -- `python3 -c "import digitize_a1"` silently rewrote a committed data file. It
# was the only tool in tools/ with no `if __name__ == "__main__"` guard. Refusing the
# import outright is the honest form of that guard here: there is no library in this
# module to import, every name below is an intermediate of one figure's extraction, and a
# partial run would leave the CSV half written.
if __name__ != "__main__":
    raise ImportError(
        "tools/digitize_a1.py is a script, not a module: everything it does happens at "
        "import time, including rewriting data/maps/f1_compressor_mass_flow.csv. Run it "
        "instead -- `python3 tools/digitize_a1.py` -- from the repository root."
    )

NAT = Path("validation/out/digitize/a1redo")
PAGE_NO = 56
PAGE = Path("validation/out/digitize/a1r_600-056.png")


def ensure_page() -> None:
    """Produce the rectified 600 dpi render this script reads, if it is not already there.

    `PAGE` lives under `validation/out/`, which is gitignored, and until 2026-09-13 no
    tool in `tools/reproduce_all.sh` produced it -- so the reproducibility gate that
    CLAUDE.md names as one of two `live` mechanisms could not run from a fresh clone at
    all. It failed loudly rather than falsely passing, which is the better of the two
    failures, but a gate nobody can run is not enforcing anything.

    Both steps are the two commands this module's docstring already told a reader to run
    by hand, in the same order and with the same arguments. Neither introduces a number:
    `pdftoppm` rasterizes the page and `rectify_page.rectify` fits the four printed frame
    lines and maps their quadrilateral onto its enclosing rectangle.
    """
    if not PAGE.exists():
        print(f"rectifying -> {PAGE}")
        rendered = render(PAGE_NO, 600, PAGE.parent / "a1")
        left, top, right, bottom = rectify_page.rectify(rendered, PAGE)
        print(f"  frame {left},{top},{right},{bottom}")


ensure_page()

img = load_gray(PAGE)
RAW = img < 0.6
X0, X1, Y0, Y1 = 914, 4221, 818, 5006  # plot interior, just inside the rectified frame
RFRAME = (910.0, 814.0, 4225.0, 5010.0)  # rectified frame: left, top, right, bottom (px)

# ======================================================================================
# 1. CALIBRATION -- measured on the native 300 dpi scan, never on a render
# ======================================================================================
# The old calibration for this figure was a linear least squares over every printed tick
# on a pixel-warped ("rectified") 600 dpi render, and it left a disagreement nothing on
# the page could settle: the tick lattice implied a y scale 0.37 % shorter than the frame
# span, and the x = 10 major tick sat 1.9 px right of the frame midpoint.  Both facts are
# real -- they reproduce here, in the un-warped scan -- but the first is not a conflict in
# the data, it is a WRONG MODEL.  The page bows: the y minor-tick spacing runs 42.4 px
# near the top frame, 41.2 px in the middle and 42.1 px near the bottom, a smooth 3 %
# swing worth +-3 px against a straight line.  Force a straight line through that and the
# ends go 4 px wrong in opposite directions, which is exactly the "0.37 % disagreement".


NINK = native_bitmap(PAGE_NO)


def _robust(t, v, tol=2.0, it=8):
    """Straight fit with outlier rejection, as (intercept, slope, kept residuals)."""
    keep = robust_keep(t, v, tol, it)
    s, i = np.polyfit(t[keep], v[keep], 1)
    return float(i), float(s), (v - (i + s * t))[keep]


def nframe(ink) -> dict[str, tuple[float, float]]:
    """The four frame lines as (intercept, slope): verticals x(y), horizontals y(x)."""
    out = {}
    print("native frame lines (300 dpi scan):")
    for nm, along, args in (
        ("left", 0, (432, 482, 420, 2490)),
        ("right", 0, (2092, 2140, 420, 2490)),
        ("top", 1, (388, 428, 480, 2090)),
        ("bottom", 1, (2484, 2528, 480, 2090)),
    ):
        a, b, r = _robust(*_trace(ink, along, *args, maxw=5))
        out[nm] = (a, b)
        print(
            f"  {nm:6s} {a:10.3f} + {b:+.6f}*t   rms {r.std():.2f} px  "
            f"({np.degrees(np.arctan(b)):+.3f} deg)"
        )
    return out


FR = nframe(NINK)
TL = corner(FR["top"], FR["left"])
TR = corner(FR["top"], FR["right"])
BL = corner(FR["bottom"], FR["left"])
BR = corner(FR["bottom"], FR["right"])
HM = homography([TL, TR, BR, BL], [(0, 1), (1, 1), (1, 0), (0, 0)])
HI = homography([(0, 1), (1, 1), (1, 0), (0, 0)], [TL, TR, BR, BL])
NW = float(np.hypot(*(np.array(BR) - np.array(BL))))  # frame width  in native px
NH = float(np.hypot(*(np.array(TL) - np.array(BL))))  # frame height in native px
print(f"  corners TL{np.round(TL, 2)} TR{np.round(TR, 2)} BR{np.round(BR, 2)} BL{np.round(BL, 2)}")


def to_uv(px, py):
    """Native scan pixel -> the unit square whose corners are the four frame corners."""
    return _apply(HM, px, py)


def to_native(u, v):
    return _apply(HI, u, v)


def calibrate(ink, fr):
    """Every printed tick on all four edges, and the held-out frame-edge model test.

    x runs 0..20 with a minor tick every 1.0 and one long major at 10; y runs 2..12 with a
    minor every 0.2 and majors at 4/6/8/10.  Both axes print their end values AT the frame
    lines, so fitting the interior ticks alone leaves four numbers to predict that the fit
    was never given.  That is the whole model-selection argument, and it is free."""

    def T(nm, along, sign, lo, hi):
        # maxlen 40, not the shared default 60: A1's ticks are short and its flat
        # speed-line extensions are long, so a deeper probe starts finding data.
        # `fr` holds (intercept, slope); the shared probe follows a POLYNOMIAL, highest
        # power first, so the pair is reversed here rather than guessed at inside it.
        a, b = fr[nm]
        pr = _protrusion(ink, (b, a), along, sign, lo, hi, maxlen=40)
        # a peak of 20 or more is a frame corner or a marker sitting on the line
        return np.array([c for c, pk, _ in _tick_groups(pr, lo, 2) if pk < 20])

    bot = T("bottom", 1, -1, 462, 2110)
    top = T("top", 1, +1, 462, 2110)
    lef = T("left", 0, +1, 412, 2492)
    rig = T("right", 0, -1, 412, 2508)
    n = (len(bot), len(top), len(lef), len(rig))
    assert n == (19, 19, 49, 49), f"tick counts {n}, wanted (19,19,49,49)"
    # Assert on STRUCTURE, not on count (A7's lesson: a tick filter can return the right
    # number of ticks and the wrong ticks).  A minor-tick lattice has no gaps: every
    # consecutive spacing must be within 10 % of that edge's median, or a tick was missed
    # or invented and the count assertion above happened to survive it.
    for nm, arr in (("bottom", bot), ("top", top), ("left", lef), ("right", rig)):
        g = np.diff(np.sort(arr))
        assert np.abs(g / np.median(g) - 1.0).max() < 0.10, (
            f"{nm} tick lattice is irregular: spacings {np.round(g, 1)}"
        )
    XV = np.arange(1.0, 20.0, 1.0)  # x minors, 1..19
    YV = np.arange(11.8, 2.0, -0.2)  # y minors, 11.8 down to 2.2
    # a tick's own edge supplies the second coordinate, so it lands on the frame line
    ub, _ = to_uv(bot, [fr["bottom"][0] + fr["bottom"][1] * x for x in bot])
    ut, _ = to_uv(top, [fr["top"][0] + fr["top"][1] * x for x in top])
    _, vl = to_uv([fr["left"][0] + fr["left"][1] * y for y in lef], lef)
    _, vr = to_uv([fr["right"][0] + fr["right"][1] * y for y in rig], rig)
    U, XX = np.concatenate([ub, ut]), np.concatenate([XV, XV])
    V, YY = np.concatenate([vl, vr]), np.concatenate([YV, YV])

    print("calibration -- held-out frame-edge test (the fit is never told these four):")
    for nm, edges, vals, e0, e1, scale, unit in (
        ("x", (ub, ut), XV, 0.0, 20.0, NW / 20.0, "PS3/P2"),
        ("y", (vl, vr), YV, 2.0, 12.0, NH / 10.0, "lbm/s"),
    ):
        uu = np.concatenate(edges)
        yy = np.concatenate([vals, vals])
        for d in (1, 2, 3, 5):
            c = np.polyfit(uu, yy, d)
            r = yy - np.polyval(c, uu)
            p0, p1 = np.polyval(c, 0.0), np.polyval(c, 1.0)
            print(
                f"  {nm} deg{d}: ticks-only fit rms {r.std() * scale:.2f} px"
                f"   frame lines -> {p0:+.5f} / {p1:.5f} {unit}"
                f"   err {(p0 - e0) * scale:+.2f} / {(p1 - e1) * scale:+.2f} px"
            )
        # the model-free version: extrapolate only the five ticks nearest each frame line
        loc = []
        for arr in edges:
            o = np.argsort(arr)
            a, b = arr[o], np.asarray(vals)[o]
            loc.append((np.polyval(np.polyfit(a[:5], b[:5], 1), 0.0) - e0) * scale)
            loc.append((np.polyval(np.polyfit(a[-5:], b[-5:], 1), 1.0) - e1) * scale)
        print(
            f"  {nm} local: the five ticks nearest each frame line predict it to "
            f"{max(abs(q) for q in loc):.2f} px "
            f"({', '.join(f'{q:+.2f}' for q in loc)}) -- lattice and frame agree LOCALLY"
        )
    cx = _anchored(U, XX, 3, 0.0, 20.0)
    cy = _anchored(V, YY, 3, 2.0, 12.0)
    # ... and the same cubic fitted to each ticked edge ALONE, so the map used at a point
    # is the one that point's own end of the plot measured. Pooling the two, which this
    # file did until 2026-09-16, averages away a real difference -- it is the same defect
    # `figure.axis.blend_axis` exists to avoid on A2, A7 and A9, and the one
    # `digitize_single_curve` was corrected for the same day.
    cx_edges = (_anchored(ub, XV, 3, 0.0, 20.0), _anchored(ut, XV, 3, 0.0, 20.0))
    cy_edges = (_anchored(vl, YV, 3, 2.0, 12.0), _anchored(vr, YV, 3, 2.0, 12.0))
    for nm, (e0, e1), tt in (("x", cx_edges, U), ("y", cy_edges, V)):
        d = np.abs(np.polyval(e0, tt) - np.polyval(e1, tt))
        sc = (NW / 20.0) if nm == "x" else (NH / 10.0)
        print(
            f"  {nm} edge disagreement: max {d.max():.5f} ({d.max() * sc:.2f} px), "
            f"mean {d.mean():.5f} -- blended across the plot, not pooled"
        )
    rx = XX - np.polyval(cx, U)
    ry = YY - np.polyval(cy, V)
    print("  chosen: CUBIC on both axes, refitted with the printed frame values anchored")
    print(
        f"    x(u) = {np.array2string(cx, precision=7)}   tick resid rms "
        f"{rx.std() * NW / 20.0:.2f} px, frame err "
        f"{np.polyval(cx, 0.0):+.5f} / {np.polyval(cx, 1.0) - 20:+.5f}"
    )
    print(
        f"    y(v) = {np.array2string(cy, precision=7)}   tick resid rms "
        f"{ry.std() * NH / 10.0:.2f} px, frame err "
        f"{np.polyval(cy, 0.0) - 2:+.5f} / {np.polyval(cy, 1.0) - 12:+.5f}"
    )
    # spread over the alternatives, at the points the map is actually read at
    alt = {}
    for d in (2, 4, 5):
        alt[("x", d)] = _anchored(U, XX, d, 0.0, 20.0)
        alt[("y", d)] = _anchored(V, YY, d, 2.0, 12.0)
    return cx, cy, alt, cx_edges, cy_edges


CX, CY, CALT, CXE, CYE = calibrate(NINK, FR)


def to_data(px, py):
    """A RECTIFIED 600 dpi pixel -> data.  The rectified frame is a rectangle by
    construction, so (px,py) -> (u,v) is exact; (u,v) is the same square the native
    homography maps onto, which is what lets the markers be measured on the render and
    calibrated on the scan.  `bridge_check` below measures the error of that step."""
    u = (np.asarray(px, float) - RFRAME[0]) / (RFRAME[2] - RFRAME[0])
    v = (RFRAME[3] - np.asarray(py, float)) / (RFRAME[3] - RFRAME[1])
    return float(blend(CXE, u, v)), float(blend(CYE, v, u))


def to_uv_rect(px, py):
    return (
        (np.asarray(px, float) - RFRAME[0]) / (RFRAME[2] - RFRAME[0]),
        (RFRAME[3] - np.asarray(py, float)) / (RFRAME[3] - RFRAME[1]),
    )


# ---------------------------------------------------------------- 2. ink layers
lab, n = ndimage.label(RAW[Y0:Y1, X0:X1], structure=np.ones((3, 3)))
sizes = ndimage.sum(RAW[Y0:Y1, X0:X1], lab, range(1, n + 1))
objs = ndimage.find_objects(lab)
DOTS = np.zeros_like(RAW)
dotc = []
for i in range(n):
    if not (8 <= sizes[i] <= 220):
        continue
    sl = objs[i]
    if (sl[0].stop - sl[0].start) > 24 or (sl[1].stop - sl[1].start) > 24:
        continue
    m = lab[sl] == i + 1
    DOTS[Y0 + sl[0].start : Y0 + sl[0].stop, X0 + sl[1].start : X0 + sl[1].stop] |= m
    cy, cx = ndimage.center_of_mass(m)
    dotc.append((cx + sl[1].start + X0, cy + sl[0].start + Y0))
INK = RAW & ~DOTS
dotc = np.array(dotc)
lg = (dotc[:, 0] > 2930) & (dotc[:, 1] > 3330) & (dotc[:, 1] < 4720)
vt = (dotc[:, 0] > 1035) & (dotc[:, 0] < 1115)
CROSSDOTS = dotc[~lg & ~vt]
XV = dotc[vt][:, 0].mean()
print(f"vertical construction line x = {XV:.2f} px  (n={vt.sum()} dots)")

order = np.argsort(sizes)[::-1][:11] + 1
cur = []
for cid in order:
    yy, xx = np.where(lab == cid)
    cur.append((yy.mean(), cid, xx + X0, yy + Y0))
cur.sort()
LABELS = list(range(11, 0, -1))
CURVE = {lbl: (c[1], c[2], c[3]) for lbl, c in zip(LABELS, cur, strict=True)}

# ---------------------------------------------------------------- 3. horizontal extensions
HFIT = {}
for lbl in LABELS:
    cid, cx, cy = CURVE[lbl]
    m = np.zeros((Y1 - Y0, X1 - X0), bool)
    m[cy - Y0, cx - X0] = True
    op = ndimage.binary_opening(m, structure=np.ones((1, 201)))
    yy, xx = np.where(op)
    xa, xb = xx.min() + X0, xx.max() + X0
    y0 = int(yy.mean()) + Y0
    pts = []
    for x in range(xa + 40, xb - 40, 15):
        col = img[y0 - 14 : y0 + 15, x]
        if col.min() > 0.5:
            continue
        mm = col < 0.5
        if mm.sum() > 14:
            continue  # a glyph, not the bare line
        w = 1.0 - col[mm]
        rr = np.arange(y0 - 14, y0 + 15)[mm]
        pts.append((x, (rr * w).sum() / w.sum()))
    pts = np.array(pts)
    p = np.polyfit(pts[:, 0], pts[:, 1], 1)
    HFIT[lbl] = (p, xa, xb, np.std(pts[:, 1] - np.polyval(p, pts[:, 0])))
    print(
        f"curve {lbl:2d}: horizontal fit slope {p[0]:+.5f}  y@1075={np.polyval(p, XV):8.2f} "
        f" x=[{xa},{xb}]  rms {HFIT[lbl][3]:.2f} px  n={len(pts)}"
    )


# ---------------------------------------------------------------- 4. glyph seeds
def density_peaks(lbl, win=33, sep=27, thr=0.30, xcut=None):
    cid, cx, cy = CURVE[lbl]
    m = np.zeros((Y1 - Y0, X1 - X0), float)
    m[cy - Y0, cx - X0] = 1.0
    d = ndimage.uniform_filter(m, size=win, mode="constant")
    d[m == 0] = 0
    if xcut:
        d[:, : xcut - X0] = 0
    pk = (d == ndimage.maximum_filter(d, size=sep, mode="constant")) & (d > thr)
    pl_, pn = ndimage.label(pk)
    c = ndimage.center_of_mass(pk, pl_, range(1, pn + 1))
    return sorted((p[1] + X0, p[0] + Y0) for p in c)


def vstrokes(lbl, L=33):
    cid, cx, cy = CURVE[lbl]
    m = np.zeros((Y1 - Y0, X1 - X0), bool)
    m[cy - Y0, cx - X0] = True
    op = ndimage.binary_opening(m, structure=np.ones((L, 1)))
    pl_, pn = ndimage.label(op, structure=np.ones((3, 3)))
    o = []
    for i in range(1, pn + 1):
        yy, xx = np.where(pl_ == i)
        if len(yy) < 60:
            continue
        o.append((xx.mean() + X0, yy.mean() + Y0, yy.max() - yy.min() + 1))
    return sorted(o)


def holes(lbl):
    cid, cx, cy = CURVE[lbl]
    m = np.zeros((Y1 - Y0, X1 - X0), bool)
    m[cy - Y0, cx - X0] = True
    hl = ndimage.binary_fill_holes(m) & ~m
    pl_, pn = ndimage.label(hl)
    o = []
    for i in range(1, pn + 1):
        yy, xx = np.where(pl_ == i)
        if len(yy) < 55:
            continue
        o.append((xx.mean() + X0, yy.mean() + Y0))
    # merge holes split by the joining line
    o.sort()
    mg = []
    for x, y in o:
        if mg and x - mg[-1][-1][0] < 12:
            mg[-1].append((x, y))
        else:
            mg.append([(x, y)])
    return [(np.mean([p[0] for p in g]), np.mean([p[1] for p in g])) for g in mg]


def merge_close(pts, tol=20.0):
    out = []
    for x, y in pts:
        if out and abs(x - out[-1][-1][0]) < tol and abs(y - out[-1][-1][1]) < tol:
            out[-1].append((x, y))
        else:
            out.append([(x, y)])
    return [(np.mean([q[0] for q in g]), np.mean([q[1] for q in g])) for g in out]


SEEDS = {}
for lbl in range(9, 1, -1):  # curves 2..9: single-digit glyphs
    SEEDS[lbl] = merge_close(density_peaks(lbl, xcut=int(HFIT[lbl][2]) - 140))
SEEDS[1] = [(x, y) for x, y, h in vstrokes(1) if x > 1300]  # six '1' strokes
h10 = [p for p in holes(10) if p[0] > 3400]  # six '0' bowls
OFF10 = np.mean([p[0] for p in holes(10) if p[0] < 1200]) - XV
SEEDS[10] = [(x - OFF10, y) for x, y in h10]
print(f"\ncurve 10 '0'-bowl offset from the string centre: {OFF10:.1f} px")


# ---- curve 11: split the stroke blob where two '1's merged vertically, then pair by y
def vstrokes11():
    m = np.zeros((Y1 - Y0, X1 - X0), bool)
    m[1560 - Y0 : 1820 - Y0, 3600 - X0 : 4050 - X0] = INK[1560:1820, 3600:4050]
    op = ndimage.binary_opening(m, structure=np.ones((33, 1)))
    pl_, pn = ndimage.label(op, structure=np.ones((3, 3)))
    out = []
    for i in range(1, pn + 1):
        yy, xx = np.where(pl_ == i)
        if len(yy) < 60:
            continue
        if yy.max() - yy.min() + 1 > 60:  # two strokes stacked
            cut = (yy.max() + yy.min()) // 2
            for s in (yy <= cut, yy > cut):
                out.append((xx[s].mean() + X0, yy[s].mean() + Y0))
        else:
            out.append((xx.mean() + X0, yy.mean() + Y0))
    return sorted(out)


st = vstrokes11()
grp, g = [], [st[0]]
for x, y in st[1:]:
    if abs(y - g[-1][1]) <= 8 and x - g[-1][0] <= 40:
        g.append((x, y))
    else:
        grp.append(g)
        g = [(x, y)]
grp.append(g)
SEEDS[11] = [(np.mean([p[0] for p in q]), np.mean([p[1] for p in q])) for q in grp]
print("\ncurve 11 strokes after split:", [(round(x), round(y)) for x, y in st])
print("curve 11 glyph centres:", [(round(x, 1), round(y, 1)) for x, y in SEEDS[11]])
for lbl in LABELS:
    print(
        f"seeds curve {lbl:2d}: n={len(SEEDS[lbl])} "
        + " ".join(f"({x:.0f},{y:.0f})" for x, y in SEEDS[lbl])
    )


# ---------------------------------------------------------------- 5. refine + assemble
# The bounding-box refinement was validated on the eleven left endpoints (isolated glyphs)
# but does NOT generalise to the sloped runs: where two markers are within a glyph width it
# merges or clips their column runs and moves them by up to 10 px, wrecking the spacing.
# Structural seeds are used unrefined.
PX = {}
for lbl in LABELS:
    pts = list(SEEDS[lbl])
    p, xa, xb, rms = HFIT[lbl]
    left = (XV, float(np.polyval(p, XV)))
    knee = (pts[0][0], float(np.polyval(p, pts[0][0])))  # knee y == the flat level
    PX[lbl] = [left, knee] + pts[1:]
    print(
        f"curve {lbl:2d}: knee glyph y {pts[0][1]:.1f} vs flat level {knee[1]:.1f} "
        f"(delta {pts[0][1] - knee[1]:+.1f} px)"
    )

# ---- estimator check on the eleven left endpoints (true x = the construction line)
lp = []
for lbl in range(9, 1, -1):
    pk = density_peaks(lbl, xcut=None)
    lp.append(min(pk, key=lambda q: q[0])[0])
lp = np.array(lp)
print("\nleft-endpoint estimator check (single-digit curves 2-9, density peak):")
print(
    f"   x mean {lp.mean():.2f} vs construction line {XV:.2f}"
    f"  ->  bias {lp.mean() - XV:+.2f} px, sd {lp.std():.2f} px"
)
kd = []
for lbl in LABELS:
    pfit = HFIT[lbl][0]
    kd.append(SEEDS[lbl][0][1] - float(np.polyval(pfit, SEEDS[lbl][0][0])))
kd = np.array(kd)
print(
    f"knee-glyph y minus fitted flat level: {kd.min():.1f}..{kd.max():+.1f} px,"
    f" sd {kd.std():.2f} px"
)

# ---------------------------------------------------------------- 6. link validation
tree = cKDTree(CROSSDOTS)


def link(a, b, skip=45):
    a = np.array(a)
    b = np.array(b)
    L = np.hypot(*(b - a))
    if L < 2 * skip + 20:
        return None
    t = np.arange(skip, L - skip, 8.0) / L
    d, _ = tree.query(a + np.outer(t, b - a))
    return float(np.median(d))


print("\nlink check (median px to nearest construction dot; <11 == the dotted line is there)")
bad = 0
for j in range(1, 11):
    row = []
    for k in range(1, 7):
        s = link(PX[j][k], PX[j + 1][k])
        row.append("  --  " if s is None else f"{s:6.1f}")
        if s is not None and s > 11:
            bad += 1
    print(f"  {j:2d}->{j + 1:2d} " + " ".join(row))
print("links failing the 11 px test:", bad)


# ---------------------------------------------------------------- 3b. bridge validation
# The markers are measured on the rectified 600 dpi render; the axes are calibrated on the
# native scan.  What joins them is the assumption that `rectify_page` produced a true
# homography, so that (px - left)/(right - left) on the rectified page is the same `u` the
# native homography produces.  That assumption is not asserted, it is MEASURED: twelve
# structures visible on both rasters -- the printed vertical construction line and the
# eleven flat extensions -- are located independently in each and compared in (u,v).
def native_structures():
    nx0, nx1, ny0, ny1 = 458, 2112, 410, 2494
    sub = NINK[ny0:ny1, nx0:nx1]
    lab, n = ndimage.label(sub, structure=np.ones((3, 3)))
    sz = ndimage.sum(sub, lab, range(1, n + 1))
    objs = ndimage.find_objects(lab)
    dots = []
    for i in range(n):
        if not (2 <= sz[i] <= 60):
            continue
        sl = objs[i]
        if (sl[0].stop - sl[0].start) > 12 or (sl[1].stop - sl[1].start) > 12:
            continue
        cy, cx = ndimage.center_of_mass(lab[sl] == i + 1)
        dots.append((cx + sl[1].start + nx0, cy + sl[0].start + ny0))
    dots = np.array(dots)
    vt = (dots[:, 0] > 520) & (dots[:, 0] < 558)
    cl = np.polyfit(dots[vt][:, 1], dots[vt][:, 0], 1)  # x = cl[1] + cl[0]*y
    cur = []
    for cid in np.argsort(sz)[::-1][:11] + 1:
        yy, xx = np.where(lab == cid)
        cur.append((yy.mean(), xx + nx0, yy + ny0))
    cur.sort()
    flats = {}
    for lbl, (_, cx, cy) in zip(range(11, 0, -1), cur, strict=True):
        m = np.zeros((ny1 - ny0, nx1 - nx0), bool)
        m[cy - ny0, cx - nx0] = True
        op = ndimage.binary_opening(m, structure=np.ones((1, 101)))
        yy, xx = np.where(op)
        xa, xb, y0 = xx.min() + nx0, xx.max() + nx0, int(yy.mean()) + ny0
        pts = []
        for x in range(xa + 20, xb - 20, 4):
            rr = np.flatnonzero(NINK[y0 - 8 : y0 + 9, x])
            if rr.size == 0 or rr.size > 8 or rr.max() - rr.min() + 1 != rr.size:
                continue  # a glyph or a gap, not the bare line
            pts.append((x, (y0 - 8) + (rr.min() + rr.max()) / 2.0))
        pts = np.array(pts)
        p = np.polyfit(pts[:, 0], pts[:, 1], 1)
        flats[lbl] = (p, float(np.std(pts[:, 1] - np.polyval(p, pts[:, 0]))), len(pts))
    return cl, flats, int(vt.sum())


def bridge_check(xv_rect, hfit_rect):
    ncl, nflat, ndot = native_structures()
    uu, _ = to_uv(np.polyval(ncl, np.linspace(500, 2400, 200)), np.linspace(500, 2400, 200))
    ur, _ = to_uv_rect(xv_rect, 0.0)
    print("\nbridge check -- rectified 600 dpi render vs native 300 dpi scan, in (u,v):")
    print(
        f"  vertical construction line: native u {uu.mean():.6f} (sd {uu.std():.6f}, "
        f"{ndot} dots)  rectified u {float(ur):.6f}   diff {(uu.mean() - float(ur)) * NW:+.2f} px"
    )
    d = []
    for lbl in range(11, 0, -1):
        p, rms, npt = nflat[lbl]
        x = float((ncl[1] + ncl[0] * p[1]) / (1 - ncl[0] * p[0]))  # flat meets the dotted line
        vn = float(to_uv(x, float(np.polyval(p, x)))[1][0])
        vr = float(to_uv_rect(0.0, float(np.polyval(hfit_rect[lbl][0], xv_rect)))[1])
        d.append((vn - vr) * NH)
        print(
            f"  flat level SYMBOL {lbl:2d}: native fit rms {rms:.2f} px over {npt:3d} samples"
            f"   native v {vn:.6f}  rectified v {vr:.6f}   diff {d[-1]:+.2f} px"
        )
    d = np.array(d)
    print(
        f"  => the two rasters agree to {np.abs(d).max():.2f} px (rms {d.std():.2f}); the "
        f"pixel warp moves nothing in this map by more than a quarter pixel."
    )
    return ncl, nflat, d


# ---------------------------------------------------------------- 7. data coordinates
NGC = {1: 65, 2: 80, 3: 82, 4: 85, 5: 87, 6: 89, 7: 92, 8: 94, 9: 96, 10: 98, 11: 100}
NCL, NFLAT, BRIDGE = bridge_check(XV, HFIT)

# ---- 1-sigma uncertainty, per point.  Three measured terms, added in quadrature:
#   estimator  the marker position estimator, from the figure's own ground truths
#   calib      the spread of the defensible calibration models AT the points read
#   raster     the rectified-render vs native-scan disagreement measured above
# All are quoted in native 300 dpi pixels first, then converted through the local scale.
SPX = NW / 20.0  # native px per unit PS3/P2
SPY = NH / 10.0  # native px per unit WA2c
UVX = np.array([to_uv_rect(PX[lbl][k][0], PX[lbl][k][1])[0] for lbl in LABELS for k in range(7)])
UVY = np.array([to_uv_rect(PX[lbl][k][0], PX[lbl][k][1])[1] for lbl in LABELS for k in range(7)])
# only models that PASSED the held-out frame-edge test count as defensible alternatives:
# on x that is degrees 2, 4 and 5; on y the quadratic fails, so only 4 and 5.
SPREAD_X = max(
    float(np.abs(np.polyval(CX, UVX) - np.polyval(CALT[("x", d)], UVX)).max()) for d in (2, 4, 5)
)
SPREAD_Y = max(
    float(np.abs(np.polyval(CY, UVY) - np.polyval(CALT[("y", d)], UVY)).max()) for d in (4, 5)
)
RASTER = float(np.abs(BRIDGE).max())
EST_X = 2.5  # px: perp scatter to the printed construction curves, resolved onto x
EST_X1 = 6.5  # px: SYMBOL 1, whose six sloped markers overprint into one ink mass
EST_Y = 2.35  # px: knee glyph minus its own fitted flat level, sd over eight curves
EST_Y0 = 0.30  # px: the flat-segment fit itself (rms 0.20-0.57 px over 38-318 samples)
EST_X0 = 0.20  # px: the printed vertical construction line (96 dots)


def sigma(lbl, k):
    ex = (EST_X0 if k == 0 else (EST_X1 if lbl == 1 else EST_X)) / SPX
    ey = (EST_Y0 if k <= 1 else EST_Y) / SPY
    cx = np.hypot(SPREAD_X, RASTER / SPX)
    cy = np.hypot(SPREAD_Y, RASTER / SPY)
    return float(np.hypot(ex, cx)), float(np.hypot(ey, cy))


print("\n  NGc  k    Ps3/P2      WA2c     1sig_x   1sig_y")
DATA, SIG = {}, {}
for lbl in LABELS:
    DATA[lbl] = [to_data(*q) for q in PX[lbl]]
    SIG[lbl] = [sigma(lbl, k) for k in range(7)]
    for k, (a, b) in enumerate(DATA[lbl]):
        sx, sy = SIG[lbl][k]
        print(f"  {NGC[lbl]:3d} {k:2d}  {a:8.3f} {b:9.4f}   {sx:6.3f} {sy:7.4f}")


# ---- verify overlay: the calibration grid drawn back onto the scan, markers ringed
def overlay():
    im = Image.fromarray((~NINK * 255).astype(np.uint8)).convert("RGB")
    d = ImageDraw.Draw(im)

    def at(u, v):
        a, b = to_native(u, v)
        return (float(a[0]), float(b[0]))

    g = np.linspace(0.0, 1.0, 20001)
    for val in range(21):  # every printed x minor tick
        u = float(g[np.argmin(np.abs(np.polyval(CX, g) - val))])
        d.line([at(u, 0.0), at(u, 1.0)], fill=(0, 120, 255) if val % 10 == 0 else (255, 60, 60))
    for val in np.arange(2.0, 12.001, 0.5):  # every fifth printed y minor tick
        v = float(g[np.argmin(np.abs(np.polyval(CY, g) - val))])
        d.line([at(0.0, v), at(1.0, v)], fill=(0, 120, 255) if val % 2 == 0 else (255, 60, 60))
    for lbl in LABELS:
        for k in range(7):
            u, v = to_uv_rect(*PX[lbl][k])
            x, y = at(float(u), float(v))
            d.ellipse([x - 7, y - 7, x + 7, y + 7], outline=(0, 170, 0), width=2)
    im.crop((380, 360, 2180, 2560)).save(NAT / "a1_verify.png")
    im.crop((380, 360, 900, 900)).resize((1040, 1080), Image.LANCZOS).save(
        NAT / "a1_verify_topleft.png"
    )
    print(f"wrote {NAT / 'a1_verify.png'} and a1_verify_topleft.png")


overlay()

# ---------------------------------------------------------------- 8. output
os.makedirs("data/maps", exist_ok=True)
out = "data/maps/f1_compressor_mass_flow.csv"
# The date the calibration was last established -- NOT the date this script runs.
# Stamping `date.today()` here made the output a function of the clock, so
# `tools/reproduce_all.sh` reported a diff on f1 on every day after the file was
# committed. A gate that always cries wolf gets ignored, which is worse than no gate:
# this is the most load-bearing map in the project (the only 2-D one, and Ps3 has ~13x
# leverage on dNG/dt). Bump this literal when the calibration actually changes, and say
# in the commit message what moved.
RECALIBRATED = "2026-09-11"
# The construction line is vertical, and the x map is blended across the plot, so it does
# not read as one number: it reads as a band. Quoting only the mid-plot value would hide
# exactly the spread this calibration exists to measure, so both ends are reported.
_ULU = float(to_uv_rect(XV, 0.0)[0])
UL_LO = float(blend(CXE, _ULU, 0.0))
UL_HI = float(blend(CXE, _ULU, 1.0))
UL = 0.5 * (UL_LO + UL_HI)
HDR = [
    "# source: TM-100991 pdf p.56, Figure A1",
    "# quantity: f1 -- compressor mass flow.  WA2c = f1(Ps3/P2, NGc), Eq. 7, pdf p.22",
    "# method: digitized -- plot-symbol centroids, calibrated on the native 300 dpi scan",
    "#",
    "# THIS IS A 2-D MAP.  Columns:",
    "#   ngc_pct        parameter: corrected gas generator speed NGC, percent.  INDEPENDENT.",
    "#                  Value as printed in the in-plot SYMBOL legend for that speed line.",
    "#   ps3_p2         x axis, as printed: 'COMPRESSOR STATIC PRESSURE RATIO, PS3/P2'.",
    "#                  Compressor static pressure ratio, nondimensional.  INDEPENDENT.",
    "#   wa2c_lbm_per_s y axis, as printed: 'STATION 2 CORRECTED MASS FLOW, WA2C, LBM/SEC'.",
    "#                  Station 2 corrected mass flow, lbm/sec.  DEPENDENT = f1.",
    "#   k              index of the symbol along its own speed line, 0..6 (see below).",
    "#   symbol         the digit printed as the plot marker for that speed line, 1..11.",
    "# Rows are grouped by speed line and ordered by increasing ps3_p2 within it.",
    "#",
    "# Structure of a speed line (11 lines x 7 points = 77):",
    "#   k=0  left end of the FLAT LEFT EXTENSION -- a horizontal solid segment at constant",
    "#        WA2c running from PS3/P2 ~= 1 rightward to the first sloped symbol.  It is part",
    "#        of the function, not decoration: f1 is constant in PS3/P2 below the knee, so a",
    "#        consumer must keep k=0 and k=1 as two breakpoints of one flat segment.",
    "#        All eleven k=0 points sit on one printed vertical dotted construction line.",
    "#   k=1  the knee: right end of the flat extension and first symbol of the sloped run.",
    "#        Its wa2c is taken from the fitted flat segment (rms 0.4-1.2 px at 600 dpi over",
    "#        21-171 samples; 0.20-0.57 px over 38-318 samples refitted on the native scan),",
    "#        not from the glyph, because the segment locates it an order of magnitude",
    "#        better.  k=0 and k=1 are the two ends of that one fitted segment, so their",
    "#        wa2c differ only by its residual slope: under 0.005 lbm/sec.  Treat the pair",
    "#        as one flat segment.",
    "#   k=2..6  the remaining five symbols on the sloped run, WA2c falling with PS3/P2.",
    "#",
    "# Speed lines, exactly as printed in the in-plot legend:",
    "#   SYMBOL  1:  65% NGC     SYMBOL  7:  92% NGC",
    "#   SYMBOL  2:  80% NGC     SYMBOL  8:  94% NGC",
    "#   SYMBOL  3:  82% NGC     SYMBOL  9:  96% NGC",
    "#   SYMBOL  4:  85% NGC     SYMBOL 10:  98% NGC",
    "#   SYMBOL  5:  87% NGC     SYMBOL 11: 100% NGC",
    "#   SYMBOL  6:  89% NGC",
    "#",
    "# CALIBRATION  (recalibrated 2026-09-11; the marker pixels are unchanged from the",
    "# 2026-09-10 extraction and are documented further down)",
    "#   Geometry is measured on the page's OWN raster and never on a render.  `pdfimages",
    "#   -list` shows p.56 is a single 300 dpi 1-bit CCITT image (2544x3300), so any higher",
    "#   dpi is an upsample.  The four frame lines are fitted sub-pixel in that bitmap:",
    *[
        f"#     {nm:6s} {'x' if nm in ('left', 'right') else 'y'} = {FR[nm][0]:9.3f}"
        f" {FR[nm][1]:+.6f} {'y' if nm in ('left', 'right') else 'x'}"
        f"   ({np.degrees(np.arctan(FR[nm][1])):+.3f} deg)"
        for nm in ("left", "right", "top", "bottom")
    ],
    "#   Four different slopes, i.e. a genuine quadrilateral, not one rotation.  It is",
    "#   carried as a HOMOGRAPHY applied to the measured points -- the four frame corners",
    "#   map to the corners of a unit square (u,v) -- so no interpolation ever touches the",
    "#   data.  The previous run warped the pixels instead; see RASTER BRIDGE below.",
    "#",
    "#   THE AXIS MAP IS A CUBIC IN (u,v), NOT A LINE.  The scan bows: the y minor-tick",
    "#   spacing runs 42.4 px near the top frame, 41.2 px in the middle and 42.1 px near",
    "#   the bottom, a smooth swing worth +-3 px against a straight line.  The model is",
    "#   chosen by a test the figure gives away free and that the fit is never told:",
    "#   BOTH axes print their end values AT the frame lines (0 and 20 in x; 2 and 12 in",
    "#   y), so a fit made on the interior ticks alone -- 19 minor ticks per edge in x,",
    "#   49 per edge in y, all four edges measured -- has four numbers left to predict.",
    "#",
    "#     axis  model      tick-fit rms   predicts the two printed frame values to",
    "#     x     linear        0.80 px      -2.38 px  and  -0.95 px      FAIL",
    "#     x     quadratic     0.63 px      -1.09 px  and  +0.32 px      pass",
    "#     x     CUBIC         0.46 px      +0.48 px  and  -1.24 px      pass  <- used",
    "#     x     quintic       0.38 px      -1.26 px  and  +0.82 px      pass",
    "#     y     linear        1.94 px      -4.00 px  and  +3.71 px      FAIL",
    "#     y     quadratic     1.94 px      -4.36 px  and  +3.35 px      FAIL",
    "#     y     CUBIC         0.91 px      +0.79 px  and  -1.79 px      pass  <- used",
    "#     y     quintic       0.74 px      -1.05 px  and  +1.10 px      pass",
    "#",
    "#   On y the QUADRATIC IS WORTHLESS -- identical to the line, because the bow is odd,",
    "#   not even.  Cubic is the lowest order that works, and it halves the tick residual.",
    "#   Having passed, the chosen model is refitted with the four printed frame values",
    "#   included as anchors:",
    "#   THE TWO TICKED EDGES ARE FITTED SEPARATELY AND BLENDED, not pooled (2026-09-16).",
    "#   Pooling averages away a real difference between the two lattices. Figure A1 gives",
    "#   the check free: its eleven flat extensions are horizontal lines spanning the plot,",
    "#   and nothing in the calibration is told so. Their end-to-end difference falls from",
    "#   0.573 px mean / 1.005 px worst when pooled to 0.409 / 0.921 when blended.",
    f"#     PS3/P2, bottom edge (u) = {np.array2string(CXE[0], precision=8, separator=' ')}",
    f"#     PS3/P2, top edge    (u) = {np.array2string(CXE[1], precision=8, separator=' ')}",
    f"#     WA2C,   left edge   (v) = {np.array2string(CYE[0], precision=8, separator=' ')}",
    f"#     WA2C,   right edge  (v) = {np.array2string(CYE[1], precision=8, separator=' ')}",
    "#   blended by the other coordinate: x by v, y by u. The pooled maps they replace were",
    f"#     PS3/P2(u) = {np.array2string(CX, precision=8, separator=' ')}",
    f"#     WA2C(v)   = {np.array2string(CY, precision=8, separator=' ')}",
    "#   (numpy polynomial order, highest power first; u = 0 at the left frame line and 1",
    "#   at the right, v = 0 at the bottom frame line and 1 at the top.)",
    "#",
    "#   THE OLD FRAME-vs-TICK DISAGREEMENT IS RESOLVED, AND IT WAS OUR MODEL, NOT THE",
    "#   PAGE.  The 2026-09-10 header said 'the tick lattice and the frame lines genuinely",
    "#   disagree ... nothing on the page resolves it', on the evidence that a linear tick",
    "#   fit implies a y scale 0.33 % short of the frame span.  It reproduces here (0.37 %",
    "#   on the un-warped scan) and it is an artefact of forcing a straight line through a",
    "#   bowed lattice: fit only the FIVE TICKS NEAREST each frame line and extrapolate,",
    "#   and the lattice predicts all four printed frame values to 1.41 px or better",
    "#   (x: -0.06, -1.41, -0.07, -0.38; y: +0.61, +0.55, -0.64, +1.16).  Locally the",
    "#   ticks and the frame agree everywhere.  Neither of the two old calibrations was",
    "#   right: the frame lines were right, and the all-ticks fit was right about the data",
    "#   and wrong about the model.  The x = 10 major tick really does sit 1.9 px right of",
    "#   the frame midpoint (that part was a correct measurement, not a warp artefact).",
    f"#   The printed vertical construction line now reads PS3/P2 = {UL:.4f} "
    f"({UL_LO:.4f} at the bottom scale, {UL_HI:.4f} at the top), against 0.974",
    "#   under the old linear fit and 0.998 under a frame-lines-only calibration.  It is",
    "#   still not asserted to be exactly 1.0 -- the report does not print that -- but the",
    "#   cubic moves it most of the way there without ever being told to.",
    "#",
    "# RASTER BRIDGE  (why marker pixels from a 600 dpi render may be calibrated on the",
    "# 300 dpi scan)",
    "#   The markers are still measured on validation/out/digitize/a1r_600-056.png, because",
    "#   the 2x upsample's interpolated grey edges are what resolve SYMBOL 1's six",
    "#   overprinting markers -- at native 300 dpi they are 7-13 px apart under a ~10 px",
    "#   glyph and a morphological opening cannot separate them.  The rectified frame is a",
    "#   rectangle by construction, so a rectified pixel maps to (u,v) exactly.  That the",
    "#   two rasters share one (u,v) is MEASURED, not assumed: twelve structures visible in",
    "#   both -- the printed vertical construction line and the eleven flat extensions --",
    f"#   were located independently in each and agree to {RASTER:.2f} px"
    f" (rms {BRIDGE.std():.2f} px)",
    "#   in native pixels.  Re-running the density-peak estimator directly on the native",
    "#   bitmap likewise reproduces the 600 dpi marker centres to <=1.0 px in x and <=0.7 px",
    "#   in y on the 48 markers of the single-digit curves.  The pixel warp is therefore",
    "#   worth a quarter pixel here, and the geometry it would have cost is now avoided.",
    "#",
    "# HOW EACH SPEED LINE WAS SEPARATED AND DIGITIZED  (unchanged from 2026-09-10)",
    "#   Per-curve ROIs were not needed and were not used: at 600 dpi each speed line -- its",
    "#   flat extension, its joining polyline and all seven of its digit markers -- is a",
    "#   single 8-connected ink component, and the eleven components are disjoint (so are",
    "#   they at native 300 dpi).  Labelling separates the family exactly, including curves",
    "#   2-4 in the crowded lower left, which an ROI band cannot do.  The six dotted",
    "#   cross-plot lines and the vertical dotted line are separate small components and",
    "#   were masked out before measuring any symbol.  Speed lines identified by y order,",
    "#   top (SYMBOL 11) to bottom (SYMBOL 1).",
    "#   Symbol position estimators, one per curve, chosen by the shape of that curve's digit:",
    "#     SYMBOL 2-9   local ink-density peak on the component mask: window 33 px, min-sep 27,",
    "#                  min-fill 0.30.  Six peaks on the sloped run for every one of the eight.",
    "#     SYMBOL 10    the enclosed bowl of the '0' (fill-holes), shifted by the measured",
    "#                  -21.9 px from bowl centre to string centre.  That offset is read off",
    "#                  this curve's own left endpoint, where the true x is known.",
    "#     SYMBOL 11    the two '1' strokes, found by vertical morphological opening (33x1) and",
    "#                  paired by y; string centre = pair midpoint.  One pair overprints the",
    "#                  next and was split by its stroke height (86 px = two 44 px strokes).",
    "#     SYMBOL 1     vertical morphological opening only.  See ACCURACY -- this is the one",
    "#                  curve whose markers overprint each other.",
    "#   A bounding-box refinement of the density peaks was tried and REJECTED: it is unbiased",
    "#   on the isolated left-endpoint glyphs but where two markers are within a glyph width it",
    "#   merges or clips their column runs and moves them up to 10 px, destroying the spacing",
    "#   regularity and worsening the construction-line check.  The seeds are used unrefined.",
    "#   A generative glyph model (the method that replaced density mode on Figure A7) was",
    "#   considered and NOT built: A7's marker is one two-stroke cross drawn six times, so one",
    "#   template serves; A1's markers are eleven different digit strings 12-18 px wide at",
    "#   native 300 dpi with the joining line running through them, so it would need eleven",
    "#   templates and would be fitting shape to 3-4 ink pixels per stroke.  The check below",
    "#   shows the existing estimator is already unbiased against printed ink, so there was",
    "#   nothing for it to win.",
    "#",
    "# ATTRIBUTION CHECK (this is a 2-D map: a point on the wrong speed line is worse than a",
    "# missing point).  The figure prints six dotted construction lines, each joining the k-th",
    "# symbol of all eleven speed lines.  For every k and every adjacent pair of speed lines,",
    "# the straight link between the two extracted symbols was tested against the printed dots:",
    "# median distance to the nearest dot, sampled along the link.  All 60 links score 6.1-9.8",
    "# px (dot half-spacing is ~7 px).  Deliberately wrong pairings (k of one line to k+1 of the",
    "# next) score 16-23 px on the separated curves.  No link fails.  Every one of the 77 points",
    "# is attributed to a speed line and an index confirmed by printed ink, including SYMBOL 1.",
    "# The marker pixels did not move in this recalibration, so this check is unchanged and was",
    "# re-run to confirm it.",
    "#",
    "# ACCURACY  (measured against the figure's own internal ground truths, never used to",
    "# adjust the estimator; all pixel figures are native 300 dpi)",
    "#   1. The eleven k=0 markers sit on one printed vertical dotted line (96 dots, its own",
    "#      position good to 0.20 px).  The density estimator on the eight single-digit",
    "#      left-endpoint glyphs returns +2.42 px right of it, sd 1.44 px.  That offset does",
    "#      NOT propagate to the sloped markers and is not corrected: cutting the outgoing",
    "#      flat line out of the window changes it by 0.00 px, so it is not one-sided line",
    "#      ink, and the sloped-marker test in 3 finds no such displacement.  The output does",
    "#      not use the estimator for k=0 anyway; it uses the construction line, which is why",
    "#      every k=0 row carries the same PS3/P2.",
    "#   2. Each k=1 knee must sit at its own flat-extension level.  Knee glyph y minus that",
    "#      fitted level: -3.80..+3.49 px, sd 2.35 px = 0.011 in WA2c.  (k=1's WA2c in the",
    "#      output is the fitted level, not the glyph, so this measures the estimator that IS",
    "#      used for k=2..6.)",
    "#   3. Every sloped marker should lie on its printed construction line.  Those lines are",
    "#      CURVED, so they were fitted locally: for each marker, the dots 15-55 px away along",
    "#      the line on both sides, marker excluded, total least squares.  47 of the 66 sloped",
    "#      markers have enough dots.  Perpendicular offset: mean +0.19 px, sd 2.54 px",
    "#      (2.15 px after removing the per-k chord sagitta), range -5.14..+6.71.  Solving the",
    "#      47 offsets for a COMMON x displacement of the whole family returns -0.09 px: the",
    "#      density-peak estimator has no measurable x bias where it is actually used.",
    "#",
    "# UNCERTAINTY, 1 sigma, per point.  Three measured terms in quadrature:",
    "#   estimator  from checks 1-3 above",
    "#   calibration  the spread of the OTHER models that passed the held-out test, around",
    "#                the chosen cubic, evaluated AT the 77 points (x: anchored quadratic,",
    "#                quartic, quintic; y: quartic and quintic -- the y quadratic failed",
    f"#                and does not count): {SPREAD_X:.4f} in PS3/P2, {SPREAD_Y:.4f} in WA2c",
    f"#   raster       the rectified-vs-native bridge, {RASTER:.2f} px",
    "#",
    "#     rows                     1sig ps3_p2   1sig wa2c",
    f"#     k=0, all speed lines        {sigma(2, 0)[0]:.3f}       {sigma(2, 0)[1]:.4f}",
    f"#     k=1, SYMBOL 2..11           {sigma(2, 1)[0]:.3f}       {sigma(2, 1)[1]:.4f}",
    f"#     k=2..6, SYMBOL 2..11        {sigma(2, 2)[0]:.3f}       {sigma(2, 2)[1]:.4f}",
    f"#     k=1, SYMBOL 1               {sigma(1, 1)[0]:.3f}       {sigma(1, 1)[1]:.4f}",
    f"#     k=2..6, SYMBOL 1            {sigma(1, 2)[0]:.3f}       {sigma(1, 2)[1]:.4f}",
    "#",
    "#   SYMBOL 1 (65% NGC) is the weakest: its six sloped markers span only 0.59 in PS3/P2",
    "#   for a glyph ~0.12 wide, so they overprint into one ink mass.  The six '1' strokes",
    "#   are still individually resolved and every one is on a printed construction line, but",
    "#   its x carries 6.5 px rather than 2.5.",
    "#",
    "# WHAT MOVED FROM THE 2026-09-10 VALUES (all of it calibration; no marker pixel moved)",
    "#   ps3_p2:  the eleven k=0 rows  +0.017 (0.974 -> 0.991), which is 2.4 sigma;",
    "#            the 66 sloped rows   -0.011 .. +0.007, at most 0.35 sigma.",
    "#   wa2c:    -0.0098 .. +0.0104, i.e. -0.215 % .. +0.115 %.  On the 22 k=0/k=1 rows,",
    "#            whose y comes from a fitted line and is therefore the most precisely",
    "#            located quantity on the figure, that is up to 2.8 sigma; on the 55 glyph",
    "#            rows it is at most 0.89 sigma.  The sign of the move follows the bow, not",
    "#            the data: rows above WA2c ~ 7.0 moved up by 0.003-0.010, rows between 4.5",
    "#            and 6.8 moved down by 0.001-0.010, and SYMBOL 1 near 3.0 barely moved.",
    "#            Nobody would have found that pattern by looking at the extracted numbers.",
    "#   NOT transcribed. Values carry read error; see validation/out/digitize/a1_zoom_c*.png.",
    "# points: 77 (11 speed lines x 7)",
    f"# digitized: 2026-09-10, recalibrated {RECALIBRATED}, by tools/digitize_a1.py; rerunnable",
]
with open(out, "w", newline="") as fh:
    for line in HDR:
        fh.write(line + "\n")
    w = csv.writer(fh)
    w.writerow(["ngc_pct", "ps3_p2", "wa2c_lbm_per_s", "k", "symbol"])
    for lbl in LABELS:
        for k, (a, b) in enumerate(DATA[lbl]):
            w.writerow([NGC[lbl], f"{a:.3f}", f"{b:.4f}", k, lbl])
print("wrote", out)

# ======================================================================================
# THE BETA-GRIDDED MAP -- traced from the printed line, not interpolated between markers
# ======================================================================================
# Beta is the standard compressor-map coordinate: 0 at the choked end of a speed line, 1 at
# surge, so that one beta names corresponding points on lines whose pressure-ratio ranges
# differ by a factor of six. Interpolating between speed lines at constant beta is what
# performance models do, and it is what `t700.maps.SpeedMap` does with this file.
#
# WHY NOT THE PRINTED CONSTRUCTION LINES. Figure A1 draws six dotted lines joining the k-th
# marker of all eleven speed lines, and those ARE beta lines -- Ballin's own. Using them is
# the obvious choice and it is not the best one, because the marker positions carry
# digitizing noise and using them as the correspondence propagates it. Expressed as a
# fraction of each line's own span, the k=1 markers sit at
#
#     0.7846 0.8174 0.8306 0.8318 0.8135 0.8113 0.8549 0.8681 0.8818 0.8914 0.8973
#
# -- up, up, up, DOWN, DOWN, up -- which is jitter, not physics (rms second difference
# 0.0216). The span they are fractions of is smooth. So beta here is defined analytically:
#
#     beta = (Ps3/P2 - X0) / (Ps3/P2 at surge - X0),   X0 = 0.9907, the plot's left frame
#
# WHAT THE CURVE COMES FROM, AND WHY IT CHANGED ON 2026-09-15. Until this revision y came
# from `scipy.interpolate.PchipInterpolator` through the seven printed markers -- an
# INTERPOLANT, whose residual at every marker is zero by construction. Six noisy points per
# line is very little to build a map from, and an exact interpolant through them turns the
# glyph-centroid read error straight into curvature: five of the eleven lines came out with
# a non-monotone secant sequence (94 %: -0.015 -0.170 -0.117 -0.396 -0.398), and
# Fritsch-Carlson then amplified it, because it sets each node's derivative to a weighted
# harmonic mean of the adjacent secants and so pins a secant that is flat within noise to a
# dead-flat shelf followed by a cliff. Local dWA2c/d(Ps3/P2) differed from any smooth fit by
# 14-135 % on ten of the eleven lines.
#
# The markers were never the best measurement on this page. The printed LINE is: it is
# continuous ink, and it is drawn as straight segments between Ballin's own table points, so
# it carries far more samples of the same function. This section traces it.
#
#   1. TRACE. Cut perpendicular to a guide curve, take the ink run nearest the guide, keep
#      it when the run is thin enough to be line rather than glyph. The first guide is the
#      seven-point extraction above -- good to a couple of pixels, far better than needed to
#      aim a cross-section -- and each fit becomes the next guide. Three passes.
#
#      A sequential walker was tried first and abandoned: it has to coast blind through
#      every glyph, and tuning the coast to rescue SYMBOL 2 broke SYMBOL 9, retightening it
#      broke SYMBOL 5 and 6. Guide sampling carries no state along the curve, so it cannot
#      drift; the guide supplies the direction everywhere and the perpendicular is always
#      right. All eleven lines converge, 1104-1343 samples each.
#
#   2. FIT. Values at 26 knots, linear between, with the shape the physics gives free
#      imposed as hard constraints -- first differences <= 0 (a compressor passes less
#      corrected flow against a higher pressure ratio) and second differences <= 0 (flat at
#      choke, steepening into surge) -- plus a second-difference roughness penalty. Smooth
#      and monotone by construction, so no abrupt change can be invented. The knee is found,
#      not assumed: it is where the fitted curve leaves the flat level.
#
#   3. BOTH MEASUREMENTS, EACH AT ITS OWN VARIANCE. The traced samples (sigma 0.45 px, the
#      residual a well-traced line reaches) and the seven printed markers (sigma 2.35 px,
#      the estimator figure from ACCURACY check 2 above) are fitted together, weighted by
#      1/sigma^2. Where the trace is rich the markers are swamped and change nothing. Where
#      the glyphs overprint into one ink mass they are the only measurement left: SYMBOL 1's
#      six '1' glyphs leave FOUR OF ITS FIVE SEGMENTS with no traced sample at all, and
#      without the markers the fit there is a smoothness prior filling a vacuum -- it put
#      surge 21 px high and drew a corner. This is one rule for all eleven lines, not a
#      special case for SYMBOL 1.
#
# MEASURED. Fit residual 0.23-0.60 px rms per line. Surge endpoints land within -7.2..+8.6 px
# of the seven-point extraction, SYMBOL 1 at 0.00. The stored map moves from the PCHIP
# version by 1.76 px rms, 8.58 px worst. Against the report: Table B.1 rms 0.2439 -> 0.2211 %,
# worst -0.900 -> -0.806 %, Table 1 eigenvalue rms 5.494 -> 5.379, whole-curve Figure 9 %NG
# 1.474 -> 1.400 and Ps3 0.760 -> 0.717, closed loop 0.0981 -> 0.0945; against Figure 8 rms
# 0.4247 -> 0.4564 and Appendix B b rms 5.222 -> 5.348. See open questions #60 and #69.
#
# STILL A DEPARTURE FROM REPLICATION -- open question #60. Figure A1's speed lines are DRAWN
# as polylines, so the seven-point file is Ballin's table and linear interpolation between
# those points is his model; `maps.f1_as_printed()` still provides it.
#
# GRID. 8 values in the choked part and 48 in the unchoked, split at beta = 0.70. The knee
# sits below beta 0.90 on every line, so the split is below every one of them and the whole
# unchoked run gets the fine spacing. The choked run is flat and needs nothing.
#
# `--raw-only` skips this file.
BETA_CHOKED = 8
BETA_UNCHOKED = 48
BETA_SPLIT = 0.70
TRACE_KNOTS = 26
TRACE_PASSES = 3
SIG_TRACE_PX = 0.45
SIG_MARKER_PX = 2.35
CROSS_OFFSETS = np.arange(-9.0, 9.0001, 0.25)

if "--raw-only" not in sys.argv:
    from scipy.optimize import minimize

    dense_out = "data/maps/f1_compressor_mass_flow_beta.csv"
    x_left = min(DATA[lbl][0][0] for lbl in LABELS)

    _UG = np.linspace(0.0, 1.0, 400001)
    # The forward map is blended, so its inverse depends on BOTH coordinates. This inverse
    # is used only to place the tracing guide, where a mid-plot blend is far better than
    # the couple of pixels the guide needs; `native_to_data`, which produces every stored
    # value, blends exactly.
    _GX, _GY = blend(CXE, _UG, 0.5), blend(CYE, _UG, 0.5)

    def data_to_native(x, y):
        """Data coordinates -> native 300 dpi scan pixel. Approximate inverse, guide only."""
        u = np.interp(x, _GX, _UG)
        v = np.interp(y, _GY, _UG)
        a, b = to_native(u, v)
        return float(np.atleast_1d(a)[0]), float(np.atleast_1d(b)[0])

    def native_to_data(px, py):
        u, v = to_uv(px, py)
        uu, vv = float(np.atleast_1d(u)[0]), float(np.atleast_1d(v)[0])
        return float(blend(CXE, uu, vv)), float(blend(CYE, vv, uu))

    NLAB, _ = ndimage.label(NINK, structure=np.ones((3, 3)))

    def speed_line_mask(lbl):
        """The 8-connected component carrying one speed line: its flat extension, its
        joining polyline and all seven digit markers are one piece of ink, and the eleven
        components are disjoint. Seeded at the knee, whose position is already known."""
        px, py = data_to_native(*DATA[lbl][1])
        win = NLAB[int(py) - 10 : int(py) + 11, int(px) - 10 : int(px) + 11]
        ids, counts = np.unique(win[win > 0], return_counts=True)
        return NLAB == ids[int(np.argmax(counts))]

    def line_width(mask, lbl):
        """Perpendicular ink thickness, measured on the flat extension where the slope is
        zero and nothing else is printed."""
        x0 = int(data_to_native(*DATA[lbl][0])[0]) + 45
        x1 = int(data_to_native(*DATA[lbl][1])[0]) - 25
        seen = []
        for x in range(x0, max(x1, x0 + 1)):
            idx = np.nonzero(mask[:, x])[0]
            if idx.size and idx.max() - idx.min() + 1 <= 6:
                seen.append(idx.max() - idx.min() + 1)
        return float(np.median(seen)) if seen else 3.0

    def _ink_at(mask, xs, ys):
        xi = np.clip(np.round(xs).astype(int), 0, mask.shape[1] - 1)
        yi = np.clip(np.round(ys).astype(int), 0, mask.shape[0] - 1)
        return mask[yi, xi]

    def sample_along(mask, gx, gy, width):
        """One cross-section per guide point. Returns the accepted (x, y) native pixels."""
        out = []
        zero = int(np.argmin(np.abs(CROSS_OFFSETS)))
        for i in range(1, len(gx) - 1):
            dx, dy = gx[i + 1] - gx[i - 1], gy[i + 1] - gy[i - 1]
            norm = float(np.hypot(dx, dy))
            if norm < 1e-9:
                continue
            dx, dy = dx / norm, dy / norm
            nx, ny = -dy, dx
            prof = _ink_at(mask, gx[i] + CROSS_OFFSETS * nx, gy[i] + CROSS_OFFSETS * ny)
            if not prof.any():
                continue
            idx = np.nonzero(prof)[0]
            groups = np.split(idx, np.nonzero(np.diff(idx) != 1)[0] + 1)

            def gap(g):
                return 0 if g[0] <= zero <= g[-1] else min(abs(g[0] - zero), abs(g[-1] - zero))

            g = min(groups, key=gap)
            if gap(g) * 0.25 > 5.0:
                continue
            lo = CROSS_OFFSETS[g[0]] - 0.125
            hi = CROSS_OFFSETS[g[-1]] + 0.125
            if hi - lo > width + 1.6:  # a glyph is sitting on the line here
                continue
            c = 0.5 * (lo + hi)
            out.append((gx[i] + c * nx, gy[i] + c * ny))
        return out

    def shape_fit(knots, X, Y, W, lam):
        """Least squares at the knots under first- and second-difference sign constraints."""
        m = knots.size
        A = np.zeros((X.size, m))
        idx = np.clip(np.searchsorted(knots, X) - 1, 0, m - 2)
        frac = (X - knots[idx]) / (knots[idx + 1] - knots[idx])
        A[np.arange(X.size), idx] = 1.0 - frac
        A[np.arange(X.size), idx + 1] = frac
        d1 = np.zeros((m - 1, m))
        d2 = np.zeros((m - 2, m))
        for i in range(m - 1):
            d1[i, i], d1[i, i + 1] = -1.0, 1.0
        for i in range(m - 2):
            d2[i, i], d2[i, i + 1], d2[i, i + 2] = 1.0, -2.0, 1.0
        sw = np.sqrt(W)
        Aw, Yw = A * sw[:, None], Y * sw
        res = minimize(
            lambda v: float(((Aw @ v - Yw) ** 2).sum() + lam * ((d2 @ v) ** 2).sum()),
            np.interp(knots, X, Y),
            jac=lambda v: 2.0 * (Aw.T @ (Aw @ v - Yw)) + 2.0 * lam * (d2.T @ (d2 @ v)),
            constraints=[
                {"type": "ineq", "fun": lambda v: -(d1 @ v), "jac": lambda v: -d1},
                {"type": "ineq", "fun": lambda v: -(d2 @ v), "jac": lambda v: -d2},
            ],
            method="SLSQP",
            options={"ftol": 1e-16, "maxiter": 900},
        )
        resid = (A @ res.x - Y) * SPY
        return res.x, float(np.sqrt((resid**2).mean())), float(np.abs(resid).max())

    def traced_curve(lbl):
        """Trace one speed line and fit it. Returns (knots, values, rms px, max px, n)."""
        mask = speed_line_mask(lbl)
        width = line_width(mask, lbl)
        xs_end = DATA[lbl][6][0]
        knots = np.linspace(x_left, xs_end, TRACE_KNOTS)
        mx = np.array([q[0] for q in DATA[lbl]], float)
        my = np.array([q[1] for q in DATA[lbl]], float)
        values = np.interp(knots, mx, my)
        rms = mx_px = float("nan")
        n_used = 0
        for _ in range(TRACE_PASSES):
            q = np.linspace(knots[0], knots[-1], 1400)
            gy = np.interp(q, knots, values)
            guide = np.array([data_to_native(a, b) for a, b in zip(q, gy, strict=True)])
            pts = sample_along(mask, guide[:, 0], guide[:, 1], width)
            if len(pts) < 12:
                break
            D = np.array([native_to_data(a, b) for a, b in pts])
            keep = (D[:, 0] >= x_left) & (D[:, 0] <= xs_end)
            X, Y = D[keep, 0], D[keep, 1]
            W = np.full(X.size, 1.0 / SIG_TRACE_PX**2)
            X = np.concatenate([X, mx[1:]])
            Y = np.concatenate([Y, my[1:]])
            W = np.concatenate([W, np.full(mx.size - 1, 1.0 / SIG_MARKER_PX**2)])
            order = np.argsort(X, kind="stable")
            X, Y, W = X[order], Y[order], W[order]
            n_used = int(X.size)
            best = None
            for lam in (3e-5, 1e-4, 3e-4):
                v, r, m_ = shape_fit(knots, X, Y, W, lam)
                if best is None or r < best[1] * 0.98:
                    best = (v, r, m_)
            values, rms, mx_px = best
        return knots, np.asarray(values, float), rms, mx_px, n_used

    beta = np.concatenate(
        [
            np.linspace(0.0, BETA_SPLIT, BETA_CHOKED, endpoint=False),
            np.linspace(BETA_SPLIT, 1.0, BETA_UNCHOKED),
        ]
    )
    print("\ntracing the printed joining line of each speed line:")
    CURVES, chord, fit_rms, fit_max = {}, 0.0, [], []
    for lbl in LABELS:
        knots, values, rms, mx_px, n_used = traced_curve(lbl)
        CURVES[lbl] = (knots, values)
        fit_rms.append(rms)
        fit_max.append(mx_px)
        xs_end = DATA[lbl][6][0]
        xg = x_left + beta * (xs_end - x_left)
        fine = np.linspace(xg[0], xg[-1], 4000)
        chord = max(
            chord,
            float(
                np.abs(
                    np.interp(fine, xg, np.interp(xg, knots, values))
                    - np.interp(fine, knots, values)
                ).max()
            ),
        )
        knee_idx = np.nonzero((values[0] - values) * SPY > 1.0)[0]
        knee = float(knots[knee_idx[0]]) if knee_idx.size else float("nan")
        print(
            f"  SYMBOL {lbl:2d} ({NGC[lbl]:3d}%): {n_used:5d} samples, fit rms {rms:4.2f} px, "
            f"max {mx_px:5.2f} px, knee Ps3/P2 {knee:6.3f}"
        )
    POOLED_RMS = float(np.sqrt(np.mean(np.square(fit_rms))))
    WORST_FIT = float(max(fit_max))

    DHDR = [
        "# source: TM-100991 pdf p.56, Figure A1 -- DERIVED from the printed line",
        "# quantity: f1 -- compressor corrected mass flow, WA2c = f1(Ps3/P2, NGc), Eq. 7",
        "# method: digitized -- the printed joining LINE of each speed line, traced by",
        "#         perpendicular cross-section against an iterated guide, then fitted with a",
        "#         monotone, concave smoothing curve and resampled onto a common beta grid.",
        "#         The seven printed markers enter the same fit at their own variance.",
        "#",
        "# BETA is the fraction of a speed line's own pressure-ratio span:",
        f"#     beta = (Ps3/P2 - {x_left:.4f}) / (Ps3/P2 at surge - {x_left:.4f})",
        "# so beta = 0 at the choked left end of every line and 1 at its surge point, and one",
        "# beta names corresponding points on lines whose spans differ by a factor of six.",
        "#",
        "# NOT the six dotted construction lines Figure A1 prints, although those are beta",
        "# lines and are Ballin's own. Their marker positions carry digitizing noise: as a",
        "# fraction of each line's span the k=1 markers run 0.785 0.817 0.831 0.832 0.814",
        "# 0.811 0.855 0.868 0.882 0.891 0.897 -- non-monotone, rms 2nd difference 0.0216 --",
        "# while the spans themselves are smooth.",
        "#",
        f"# GRID: {BETA_CHOKED} values below beta = {BETA_SPLIT} and {BETA_UNCHOKED} above,",
        f"# {len(beta)} per line. The knee sits below beta 0.90 on every line, so the split is",
        "# below every one of them and the whole unchoked run gets the fine spacing.",
        "#",
        "# COLUMNS:",
        "#   ngc_pct         corrected gas generator speed NGc of the speed line, percent",
        "#   beta            position along the line, nondimensional: 0 at the choked end,",
        "#                   1 at surge. The interpolation coordinate.",
        "#   ps3_p2          compressor static pressure ratio Ps3/P2, nondimensional",
        "#   wa2c_lbm_per_s  station 2 corrected mass flow WA2c, lbm/sec -- the value of f1",
        "#   symbol          the plot symbol of that speed line in Figure A1, 1..11",
        "#",
        "# THE CURVE IS TRACED FROM INK, NOT INTERPOLATED BETWEEN MARKERS. Until 2026-09-15",
        "# y came from a Fritsch-Carlson interpolant through the seven printed markers, whose",
        "# residual at every marker is zero by construction; six noisy points per line is very",
        "# little to build a map from, and an exact interpolant turns the glyph read error",
        "# straight into curvature. Five of eleven lines carried a non-monotone secant",
        "# sequence and local dWA2c/d(Ps3/P2) was off by 14-135 %. The printed LINE is the",
        "# better measurement -- continuous ink, ~1200 samples per line instead of six points.",
        "#",
        "# TRACE: cut perpendicular to a guide curve, take the ink run nearest the guide, keep",
        f"# it when thinner than the local line width + 1.6 px. {TRACE_PASSES} passes, each fit",
        "# becoming the next guide, so nothing is carried along the curve and nothing drifts.",
        f"# FIT: values at {TRACE_KNOTS} knots, first differences <= 0 and second differences",
        "# <= 0 imposed as hard constraints, plus a second-difference roughness penalty. The",
        "# knee is found, not assumed. BOTH measurements are used, each weighted by its own",
        f"# variance: traced samples at {SIG_TRACE_PX} px and the printed markers at",
        f"# {SIG_MARKER_PX} px. Where the trace is rich the markers are swamped; where the",
        "# glyphs overprint into one mass they are all that is left. SYMBOL 1 has NO traced",
        "# sample in four of its five segments, and without the markers its fit put surge",
        "# 21 px high and drew a corner.",
        "#",
        "#   NOT transcribed. Values carry read error: the traced centreline, the axis",
        "#   calibration it is carried through, and the printed markers that share the fit",
        "#   all contribute. See the raw seven-point CSV's UNCERTAINTY block for the",
        "#   marker term and validation/out/digitize/a1_zoom_c*.png for the ink.",
        "#",
        f"# Fit residual, per line: {POOLED_RMS:.2f} px pooled rms, {WORST_FIT:.2f} px worst.",
        f"# The stored polyline approximates the fitted curve to {chord:.2e} lbm/sec.",
        "#",
        "# USING THIS FILE IS A DEPARTURE FROM REPLICATION -- open question #60. The raw",
        "# seven-point CSV is Ballin's table, and `maps.f1_as_printed()` loads it.",
        f"# points: {11 * len(beta)} ({len(beta)} per line x 11 speed lines)",
        "# derived: by tools/digitize_a1.py from the page raster; rerunnable, deterministic",
    ]
    with open(dense_out, "w", newline="") as fh:
        for line in DHDR:
            fh.write(line + "\n")
        w = csv.writer(fh)
        w.writerow(["ngc_pct", "beta", "ps3_p2", "wa2c_lbm_per_s", "symbol"])
        for lbl in LABELS:
            knots, values = CURVES[lbl]
            xs_end = DATA[lbl][6][0]
            for b in beta:
                x = x_left + b * (xs_end - x_left)
                y = float(np.interp(x, knots, values))
                w.writerow([NGC[lbl], f"{b:.6f}", f"{x:.4f}", f"{y:.5f}", lbl])
    print(f"wrote {dense_out}  (fit rms {POOLED_RMS:.2f} px, polyline chord {chord:.2e} lbm/sec)")
