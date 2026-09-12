"""Digitize Figure A1 (pdf p.56) -- f1, compressor mass flow. The one 2-D map in the report.

Figure A1 is the only figure in TM-100991 that `tools/digitize.py`'s CLI cannot do, so it
gets its own script rather than a shell recipe. Four things make it different:

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

Run from the repo root; writes data/maps/f1_compressor_mass_flow.csv and the per-curve
CSVs that `digitize.py verify` consumes. Inputs:

    pdfimages -f 56 -l 56 -png docs/ballin-tm100991.pdf validation/out/digitize/a1redo/native
    pdftoppm -f 56 -l 56 -r 600 -png docs/ballin-tm100991.pdf validation/out/digitize/a1_600
    python3 tools/rectify_page.py rectify --page validation/out/digitize/a1_600-056.png \
        --out validation/out/digitize/a1r_600-056.png

The first is pulled automatically if missing. This module lives in tools/ and may use
SciPy and PIL. Nothing here is imported by src/t700/.
"""

from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
from scipy.spatial import cKDTree

SP = os.environ.get("A1_WORK", "validation/out/digitize")  # intermediate .npy arrays
PDF = Path("docs/ballin-tm100991.pdf")
NAT = Path("validation/out/digitize/a1redo")
PAGE_NO = 56
PAGE = "validation/out/digitize/a1r_600-056.png"
img = np.asarray(Image.open(PAGE).convert("L"), dtype=float) / 255.0
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


def native_bitmap() -> np.ndarray:
    """The page's own 300 dpi 1-bit image, as an ink mask.  No render, no resample."""
    NAT.mkdir(parents=True, exist_ok=True)
    if not (NAT / "native-000.png").exists():
        subprocess.run(
            [
                "pdfimages",
                "-f",
                str(PAGE_NO),
                "-l",
                str(PAGE_NO),
                "-png",
                str(PDF),
                str(NAT / "native"),
            ],
            check=True,
        )
    return ~np.array(Image.open(NAT / "native-000.png"))


NINK = native_bitmap()


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    e = np.flatnonzero(np.diff(np.concatenate(([0], mask.view(np.int8), [0]))))
    return list(zip(e[0::2], e[1::2], strict=True))


def _trace(ink, along, lo, hi, tlo, thi, maxw=5):
    """Centroid of the single narrow ink run in [lo,hi) at each position along the line.

    `along=0` traces a vertical line (one run per row), `along=1` a horizontal one.  A bare
    frame line gives exactly one narrow run; a tick merges with it and widens the run past
    `maxw`, and data or text give several -- both are skipped rather than modelled.
    """
    ts, vs = [], []
    for t in range(tlo, thi):
        seg = ink[t, lo:hi] if along == 0 else ink[lo:hi, t]
        rr = [r for r in _runs(seg) if (r[1] - r[0]) <= maxw]
        if len(rr) != 1:
            continue
        a, b = rr[0]
        ts.append(float(t))
        vs.append(lo + (a + b - 1) / 2.0)
    return np.array(ts), np.array(vs)


def _robust(t, v, tol=2.0, it=8):
    keep = np.ones(t.shape, bool)
    for _ in range(it):
        s, i = np.polyfit(t[keep], v[keep], 1)
        new = np.abs(v - (i + s * t)) < tol
        if new.sum() < 20 or (new == keep).all():
            break
        keep = new
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
        a, b, r = _robust(*_trace(ink, along, *args))
        out[nm] = (a, b)
        print(
            f"  {nm:6s} {a:10.3f} + {b:+.6f}*t   rms {r.std():.2f} px  "
            f"({np.degrees(np.arctan(b)):+.3f} deg)"
        )
    return out


def corner(hor, ver):
    a, b = hor
    c, d = ver
    x = (c + d * a) / (1.0 - d * b)
    return x, a + b * x


def homography(src, dst):
    A, B = [], []
    for (sx, sy), (dx, dy) in zip(src, dst, strict=True):
        A.append([sx, sy, 1, 0, 0, 0, -dx * sx, -dx * sy])
        B.append(dx)
        A.append([0, 0, 0, sx, sy, 1, -dy * sx, -dy * sy])
        B.append(dy)
    h = np.linalg.solve(np.asarray(A, float), np.asarray(B, float))
    return np.append(h, 1.0).reshape(3, 3)


def _apply(H, px, py):
    px = np.atleast_1d(np.asarray(px, float))
    py = np.atleast_1d(np.asarray(py, float))
    d = H[2, 0] * px + H[2, 1] * py + H[2, 2]
    return (
        (H[0, 0] * px + H[0, 1] * py + H[0, 2]) / d,
        (H[1, 0] * px + H[1, 1] * py + H[1, 2]) / d,
    )


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


def _protrusion(ink, edge, along, sign, lo, hi, maxlen=40):
    """Contiguous ink growing inward from a frame edge, per position along it."""
    a, b = edge
    out = np.zeros(hi - lo, int)
    for j, t in enumerate(range(lo, hi)):
        base = a + b * t
        k = 0
        while k < maxlen:
            q = int(round(base + sign * (3 + k)))
            on = ink[t, q] if along == 0 else ink[q, t]
            if not on:
                break
            k += 1
        out[j] = k
    return out


def _groups(p, lo, thresh=2, gap=3, maxlen=20):
    """Tick centres: ink-weighted centroid of each run of positions that protrude."""
    idx = np.flatnonzero(p >= thresh)
    out, cur = [], [int(idx[0])]
    for v in idx[1:]:
        if v - cur[-1] <= gap:
            cur.append(int(v))
        else:
            out.append(cur)
            cur = [int(v)]
    out.append(cur)
    res = []
    for g in out:
        g = np.asarray(g)
        if p[g].max() >= maxlen:  # a frame corner or a marker sitting on the line
            continue
        w = p[g].astype(float)
        res.append(float(lo + (g * w).sum() / w.sum()))
    return res


def _anchored(u, y, deg, e0, e1):
    """Fit on the interior ticks WITH the two printed frame-edge values as end anchors.

    The anchors are added only after the held-out test below has chosen `deg`; the choice
    itself is made by a fit that has never seen them."""
    uu = np.concatenate([u, [0.0, 0.0, 1.0, 1.0]])
    yy = np.concatenate([y, [e0, e0, e1, e1]])
    return np.polyfit(uu, yy, deg)


def calibrate(ink, fr):
    """Every printed tick on all four edges, and the held-out frame-edge model test.

    x runs 0..20 with a minor tick every 1.0 and one long major at 10; y runs 2..12 with a
    minor every 0.2 and majors at 4/6/8/10.  Both axes print their end values AT the frame
    lines, so fitting the interior ticks alone leaves four numbers to predict that the fit
    was never given.  That is the whole model-selection argument, and it is free."""

    def T(nm, along, sign, lo, hi):
        return np.array(_groups(_protrusion(ink, fr[nm], along, sign, lo, hi), lo))

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
    return cx, cy, alt


CX, CY, CALT = calibrate(NINK, FR)


def to_data(px, py):
    """A RECTIFIED 600 dpi pixel -> data.  The rectified frame is a rectangle by
    construction, so (px,py) -> (u,v) is exact; (u,v) is the same square the native
    homography maps onto, which is what lets the markers be measured on the render and
    calibrated on the scan.  `bridge_check` below measures the error of that step."""
    u = (np.asarray(px, float) - RFRAME[0]) / (RFRAME[2] - RFRAME[0])
    v = (RFRAME[3] - np.asarray(py, float)) / (RFRAME[3] - RFRAME[1])
    return float(np.polyval(CX, u)), float(np.polyval(CY, v))


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
np.save(SP + "/crossdots_r.npy", CROSSDOTS)
np.save(SP + "/frame.npy", np.array(RFRAME))


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


def vstrokes(lbl, L=33, region=None):
    if region is None:
        cid, cx, cy = CURVE[lbl]
        m = np.zeros((Y1 - Y0, X1 - X0), bool)
        m[cy - Y0, cx - X0] = True
    else:
        xa, ya, xb, yb = region
        m = np.zeros((Y1 - Y0, X1 - X0), bool)
        m[ya - Y0 : yb - Y0, xa - X0 : xb - X0] = INK[ya:yb, xa:xb]
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
v11 = [p for p in vstrokes(11, region=(3600, 1560, 4050, 1820)) if p[0] > 3600]
pairs, cure = [], []
for x, y, _h in v11:
    if cure and x - cure[-1][0] < 40:
        cure.append((x, y))
    else:
        if cure:
            pairs.append(cure)
        cure = [(x, y)]
pairs.append(cure)
SEEDS[11] = [(np.mean([p[0] for p in g]), np.mean([p[1] for p in g])) for g in pairs]
print(f"\ncurve 10 '0'-bowl offset from the string centre: {OFF10:.1f} px")
print("curve 11 strokes:", [(round(x), round(y)) for x, y, h in v11])
for lbl in LABELS:
    print(
        f"seeds curve {lbl:2d}: n={len(SEEDS[lbl])} "
        + " ".join(f"({x:.0f},{y:.0f})" for x, y in SEEDS[lbl])
    )


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


# ---------------------------------------------------------------- 5. refine + assemble
def bbox_centre(x0, y0, gap=4, half=36, ext=18, rows=30):
    xs0 = int(round(x0)) - half
    ys0 = int(round(y0)) - rows
    w = INK[ys0 : ys0 + 2 * rows + 1, xs0 : xs0 + 2 * half + 1]
    tall = [
        j
        for j in range(w.shape[1])
        if np.flatnonzero(w[:, j]).size
        and (np.flatnonzero(w[:, j]).max() - np.flatnonzero(w[:, j]).min() + 1) >= ext
    ]
    if not tall:
        return None
    tall = np.array(tall)
    k = int(np.argmin(np.abs(tall - half)))
    lo = hi = k
    while lo > 0 and tall[lo] - tall[lo - 1] <= gap:
        lo -= 1
    while hi < len(tall) - 1 and tall[hi + 1] - tall[hi] <= gap:
        hi += 1
    cols = tall[lo : hi + 1]
    rmin = min(np.flatnonzero(w[:, j]).min() for j in cols)
    rmax = max(np.flatnonzero(w[:, j]).max() for j in cols)
    return xs0 + (cols.min() + cols.max()) / 2.0, ys0 + (rmin + rmax) / 2.0


# The bounding-box refinement was validated on the eleven left endpoints (isolated glyphs)
# but does NOT generalise to the sloped runs: where two markers are within a glyph width it
# merges or clips their column runs and moves them by up to 10 px, wrecking the spacing.
# Structural seeds are used unrefined.
REFINE = {lbl: (False, 0) for lbl in range(1, 12)}
PX = {}
for lbl in LABELS:
    ok, gap = REFINE[lbl]
    pts = []
    for sx, sy in SEEDS[lbl]:
        pts.append(bbox_centre(sx, sy, gap=gap) if ok else (sx, sy))
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
np.save(
    SP + "/px_r.npy", np.array([[NGC[lbl], k, *PX[lbl][k]] for lbl in LABELS for k in range(7)])
)
np.save(
    SP + "/data_r.npy", np.array([[NGC[lbl], k, *DATA[lbl][k]] for lbl in LABELS for k in range(7)])
)


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
UL = float(np.polyval(CX, float(to_uv_rect(XV, 0.0)[0])))
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
    f"#   The printed vertical construction line now reads PS3/P2 = {UL:.4f}, against 0.974",
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

# per-curve two-column CSVs so `digitize.py verify` can be run on each speed line
for lbl in LABELS:
    p = f"validation/out/digitize/a1_curve{lbl:02d}.csv"
    with open(p, "w", newline="") as fh:
        fh.write(f"# Figure A1 speed line SYMBOL {lbl} = {NGC[lbl]}% NGC -- for verify only\n")
        w = csv.writer(fh)
        w.writerow(["x", "y"])
        for a, b in DATA[lbl]:
            w.writerow([f"{a:.3f}", f"{b:.4f}"])
