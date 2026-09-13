"""Extract Figure A7 (`f7`, gas-generator turbine energy) from TM-100991 pdf p.62.

Why this figure gets the most care of any in the appendix.  `f7` carries six points, five
of them packed into a knee 0.015 wide whose steepest segment runs at -526 BTU/LBM per unit
P45/P41, and a sensitivity sweep of the trimmed engine gives +2.9 % to +5.5 % of shaft
power per +1 % on `f7`, with all three published trim conditions sitting inside that knee.
One pixel of read error there is 0.047 BTU/LBM.

This is the 2026-09-11 *recalibration*, replacing the same day's first redo.  What changed:

1.  **Both** y frame values are used.  The top frame prints 52.5; the bottom frame is
    unlabelled, and the tick lattice fixes it -- the frame spans 13.03 major steps of 2.5
    and 130.3 minor steps of 0.25, so it is a lattice point and the bottom frame is
    **20.0 exactly**.  With two known frame values the held-out test can separate an even
    axis correction from an odd one, which one value cannot.  x is the same story: the
    right frame prints 0.35, the left is unlabelled and the lattice spans 15.00 majors of
    0.01, so it is **0.20 exactly**.
2.  **The minor ticks are used** -- 0.25 in y (131 of them on the left edge, every slot
    filled exactly once) and 0.002 in x (76).  Ten times the data the long ticks give, and
    enough to fit and *test* an order-5 axis map instead of guessing between order 2 and
    order 4 on twelve ticks.
3.  **The frame edges are fitted as arcs** (`digitize_native.make_sv`), not as a
    homography between four corners.  p.62's sagittae are -0.44 / -0.75 / -0.38 / -0.14 px.
4.  **The frame values are anchored exactly**, not carried as two extra data points among
    76 ticks (`fit_anchored`).  With the soft anchor the order-3 x map misses the left
    frame by 2.0 px, and marker 1 *sits on that frame*.
5.  **Seeds come from a global search before the local fit.**  Marker seeds are the cores
    of a 6x6 morphological opening -- six of them, an independent marker count -- but the
    core of a marker on a frame line is the union of glyph, frame and tick and can be 10 px
    out; handed that, the local fit slid marker 1 18 px down the joining line into a
    converged-looking wrong minimum, 0.27 BTU/LBM low.  `refine_seeds` grid-searches the
    hard-edged model over +-16 px first, and the fit then moves <= 0.66 px.
6.  Nelder-Mead is given an explicit 2 px simplex on the glyph *offset*.  Its default
    simplex is 5 % of each coordinate, which at marker 6's column of x = 2118 is a 106 px
    first step into flat, converged-looking nonsense.

What the axis order came out as, and how it was chosen: see `choose_order` -- minimise the
held-out frame error among orders whose leave-one-tick-out error is within 1.2x of the
best, ties within 1.25x to the lower order.  On the fine lattice that gives **order 5 in y**
(the quadratic in use before misses the bottom frame by 0.081-0.091 BTU/LBM) and **order 3
in x**.  The bow is even *and* odd: fitted as Legendre coefficients of the y value
function, P2 = -3.8/-4.8 px (35 sigma), P3 = -3.5/-4.1 px (27 sigma), P4 = +2.6/+2.7 px
(17 sigma), P5 = +1.1 px (6 sigma), P6 nothing.  That is the answer the single printed
frame value could not give.

Run from the repo root, with `tools/` on PYTHONPATH.  Writes overlays under
validation/out/digitize/a7redo/ and prints the rows of data/maps/f7_gg_turbine_energy.csv.
Needs `pdfimages` (poppler).

This module lives in tools/ and may use SciPy and PIL.  Nothing here is imported by
src/t700/.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

from digitize import check_against_csv
from digitize_native import (
    blend_axis,
    corner,
    curve_component,
    fit_glyphs,
    frame_curves,
    label_by_map,
    label_ticks,
    make_sv,
    make_window,
    native_bitmap,
    poly_se,
    segment_lines,
    ticks,
)

OUT = Path("validation/out/digitize/a7redo")
CSV = Path("data/maps/f7_gg_turbine_energy.csv")
PAGE = 62
N = 6

X0, X1 = 0.20, 0.35  # left frame DERIVED from the 15-step lattice span; right printed
Y0, Y1 = 20.0, 52.5  # bottom frame DERIVED from the 13-step lattice span; top printed
FRAME_ORDER, X_ORDER, Y_ORDER = 2, 3, 5  # chosen by the tests main() prints
FINE = True  # calibrate on every tick (0.002 in x, 0.25 in y), not the long ticks alone

SPEC = [
    ("left", 0, (440, 500, 340, 2600)),
    ("right", 0, (2095, 2160, 340, 2600)),
    ("top", 1, (280, 345, 500, 2100)),
    ("bottom", 1, (2600, 2700, 500, 2100)),
]

# (along, sign, lo, hi, long-tick threshold, all-tick threshold, coarse step, fine step)
EDGES = {
    "left": (0, +1, 300, 2660, 7, 3, 2.5, 0.25),
    "right": (0, -1, 300, 2670, 7, 3, 2.5, 0.25),
    "bottom": (1, -1, 450, 2140, 7, 3, 0.01, 0.002),
    "top": (1, +1, 450, 2140, 7, 3, 0.01, 0.002),
}

# frame arcs 1 or 2, x order 1/3/5, y order 3/4/5, edges blended or pooled
VARIANTS = [
    (fo, xo, yo, bl, hd)
    for fo in (1, 2)
    for xo in (1, 3, 5)
    for yo in (3, 4, 5)
    for bl in (True, False)
    for hd in (True, False)
]

MASK = {0: ("left", 7.0), 5: ("right", 7.0)}  # markers that sit ON a frame line
TMPL0 = [np.radians(-44.0), 14.9, np.radians(46.0), 15.2, 3.46, 3.52]


# --------------------------------------------------------------------------- ticks


def tickset(fc, to_sv, verbose=False, fine=None):
    """Every printed tick on the four edges, labelled, frame-edge ticks removed.

    Labelling is `label_ticks` on the long ticks (a discrete slot choice from the two frame
    values, so it cannot leak into the continuous fit that is later tested against them),
    then `label_by_map` for the fine lattice through a cubic fitted to the long ticks --
    never through the straight line between the corners, which is several pixels out
    mid-page and would silently discard a block of good minor ticks.
    """
    ink = tickset.ink
    fine = FINE if fine is None else fine
    TL = corner(fc, "top", "left", (480, 303))
    TR = corner(fc, "top", "right", (2140, 325))
    BL = corner(fc, "bottom", "left", (458, 2639))
    BR = corner(fc, "bottom", "right", (2114, 2669))
    ends = {
        "left": (BL[1], TL[1], Y0, Y1),
        "right": (BR[1], TR[1], Y0, Y1),
        "bottom": (BL[0], BR[0], X0, X1),
        "top": (TL[0], TR[0], X0, X1),
    }

    def at(nm, cs):
        cs = np.asarray(cs, float)
        if nm in ("left", "right"):
            return to_sv(np.polyval(fc[nm], cs), cs)[1]
        return to_sv(cs, np.polyval(fc[nm], cs))[0]

    P = {"_xedges": (X0, X1), "_yedges": (Y0, Y1), "_span": {}}
    for nm, (along, sign, lo, hi, th_long, th_all, step_c, step_f) in EDGES.items():
        w0, w1, v0, v1 = ends[nm]
        coarse = ticks(ink, fc[nm], along, sign, lo, hi, th_long)
        wc, vc, keep_c = label_ticks(coarse, w0, w1, v0, v1, step_c, tol=0.15)
        tc = at(nm, wc)
        pc = np.polyfit(tc, vc, 3)
        if fine:
            allt = ticks(ink, fc[nm], along, sign, lo, hi, th_all)
            t, val, keep_f = label_by_map(at(nm, allt), pc, step_f)
            n_raw, n_kept = allt.size, int(keep_f.sum())
        else:
            t, val = tc, vc
            n_raw, n_kept = coarse.size, int(keep_c.sum())
        step = step_f if fine else step_c
        # how many steps does the frame span, by a STRAIGHT line through the lattice?
        # This is the free check that fixes the unlabelled frame: it only has to pick an
        # integer, so the bow's leak into it (~0.03 steps here) cannot reach the answer.
        cl = np.polyfit(t, val, 1)
        P["_span"][nm] = (
            float((np.polyval(cl, 1.0) - np.polyval(cl, 0.0)) / step_c),
            float((np.polyval(cl, 1.0) - np.polyval(cl, 0.0)) / step_f),
        )
        interior = (val > v0 + 0.4 * step) & (val < v1 - 0.4 * step)
        P[nm] = (t[interior], val[interior])
        if verbose:
            k = np.round(val / step).astype(int)
            slots = set(range(k.min(), k.max() + 1))
            print(
                f"  {nm:6s}: {coarse.size} long ticks -> {int(keep_c.sum())} labelled "
                f"(dropped {np.round(coarse[~keep_c], 1)});  {n_raw} ticks at {step} -> "
                f"{n_kept} labelled, {interior.sum()} interior; "
                f"slots {k.min()}..{k.max()}, missing {sorted(slots - set(k.tolist()))}, "
                f"duplicates {len(k) - len(set(k.tolist()))}"
            )
    return P


# --------------------------------------------------------------------------- axis model


def _loo(t, val, order):
    """Leave-one-tick-out error: the guard against an order that buys the frame with noise."""
    e = [
        val[k] - np.polyval(np.polyfit(np.delete(t, k), np.delete(val, k), order), t[k])
        for k in range(t.size)
    ]
    return float(np.std(e))


def held_out_table(t, val, lo, hi, unit=1.0, uname="", orders=(1, 2, 3, 4, 5, 6), label=""):
    """Fit interior ticks alone; predict BOTH printed frame values; report LOO and tail CV."""
    print(f"   {label}({t.size} interior ticks)")
    print(
        "     ord  interior rms       LOO     predict low frame (err)   "
        "predict high frame (err)    tail-CV rms/max"
    )
    rows = {}
    for o in orders:
        c = np.polyfit(t, val, o)
        r = val - np.polyval(c, t)
        a, b = float(np.polyval(c, 0.0)), float(np.polyval(c, 1.0))
        loo = _loo(t, val, o)
        tc, tm = _tail(t, val, o)
        rows[o] = (float(r.std()), loo, a, b)
        print(
            f"     {o}  {r.std() / unit:10.4f}{uname} {loo / unit:9.4f}{uname}  "
            f"{a:11.5f} ({(a - lo) / unit:+8.4f}{uname})  "
            f"{b:11.5f} ({(b - hi) / unit:+8.4f}{uname})"
            f"   {tc / unit:7.3f}/{tm / unit:7.3f}"
        )
    return rows


def _tail(t, val, order, frac=0.18):
    lo, hi = t.min(), t.max()
    span = hi - lo
    inner = (t > lo + frac * span) & (t < hi - frac * span)
    c = np.polyfit(t[inner], val[inner], order)
    r = val[~inner] - np.polyval(c, t[~inner])
    return float(np.sqrt((r**2).mean())), float(np.abs(r).max())


def choose_order(rows, lo, hi, unit=1.0, uname=""):
    """A6/A8/A10's rule: minimise the held-out frame error, but only among orders whose
    leave-one-tick-out error is within 1.2x of the best -- otherwise a high order can win
    the extrapolation by fitting noise in between.  Ties within 1.25x go to the lower
    order.  Both criteria are held out; neither is a fit statistic on the fitted data."""
    errs = {o: max(abs(a - lo), abs(b - hi)) for o, (_, _, a, b) in rows.items()}
    loos = {o: v[1] for o, v in rows.items()}
    lbest = min(loos.values())
    ok = [o for o in rows if loos[o] <= 1.20 * lbest]
    order = min(ok, key=lambda k: errs[k])
    for o in sorted(ok):
        if errs[o] <= 1.25 * errs[order]:
            order = o
            break
    print(
        f"     LOO within 1.2x of best: orders {ok};  worst-frame errors "
        + ", ".join(f"{o}:{errs[o] / unit:.4f}{uname}" for o in ok)
        + f"  ->  ORDER {order}"
    )
    return order


def parity(t, val, scale, deg=6):
    """Legendre coefficients of the value function, with standard errors, in page pixels.

    Even degrees are an even distortion about mid-axis, odd degrees an odd one.  This is
    the test the A7 redo could not run: one printed frame value cannot separate them, a
    lattice of 129 minor ticks can.
    """
    u = 2 * t - 1
    A = np.polynomial.legendre.legvander(u, deg)
    coef, *_ = np.linalg.lstsq(A, val, rcond=None)
    r = val - A @ coef
    s2 = (r**2).sum() / max(t.size - deg - 1, 1)
    se = np.sqrt(np.diag(np.linalg.inv(A.T @ A) * s2))
    out = []
    for k in range(2, deg + 1):
        out.append((k, coef[k] * scale, se[k] * scale, abs(coef[k] / se[k])))
    return out


# --------------------------------------------------------------------------- markers


def cores(ink, fc, k=6, margin=8.0, minsize=30):
    """Marker seeds by morphological opening -- a SHAPE detector, not a density one.

    The joining line is ~3.5 px wide and a k x k opening erases it and every tick; where
    two strokes cross, enough ink survives in every direction.  Nothing to tune per point,
    and the number of cores is an independent marker count.  The frame band is left IN
    (markers 1 and 6 sit on it), so a core's seed is centroided over its non-frame ink.
    """
    Y, X = np.mgrid[0 : ink.shape[0], 0 : ink.shape[1]]
    dL = X - np.polyval(fc["left"], Y)
    dR = X - np.polyval(fc["right"], Y)
    dT = Y - np.polyval(fc["top"], X)
    dB = Y - np.polyval(fc["bottom"], X)
    inside = (dL > -margin) & (dR < margin) & (dT > -margin) & (dB < margin)
    off_frame = (np.abs(dL) > 4.5) & (np.abs(dR) > 4.5) & (np.abs(dT) > 4.5) & (np.abs(dB) > 4.5)
    op = ndimage.binary_opening(ink & inside, np.ones((k, k)))
    lab, n = ndimage.label(op, structure=np.ones((3, 3)))
    sz = ndimage.sum(op, lab, range(1, n + 1))
    out = []
    for i in range(1, n + 1):
        if sz[i - 1] < minsize:
            continue
        m = (lab == i) & off_frame
        if m.sum() < 8:
            m = lab == i
        out.append(
            (float((X * m).sum() / m.sum()), float((Y * m).sum() / m.sum()), float(sz[i - 1]))
        )
    out.sort(key=lambda r: r[0])
    return out


def _hard_cost(ink, X, Y, use, c, T, dirs, rad):
    """Glyph cost with no supersampling -- cheap enough to scan a grid of offsets."""
    m = np.zeros(X.shape, bool)
    for ang, half, wid, halfline in (
        (T[0], T[1], T[4], False),
        (T[2], T[3], T[4], False),
        *[(np.arctan2(d[1], d[0]), 42.0, T[5], True) for d in dirs],
    ):
        u = np.array([np.cos(ang), np.sin(ang)])
        dx, dy = X - c[0], Y - c[1]
        t = dx * u[0] + dy * u[1]
        pp = -dx * u[1] + dy * u[0]
        k = np.abs(pp) <= wid / 2
        k &= (t >= 0) & (t <= half) if halfline else (np.abs(t) <= half)
        m |= k
    sel = use & (((X - c[0]) ** 2 + (Y - c[1]) ** 2) < rad**2)
    return float(((m.astype(float) - ink) ** 2)[sel].sum())


def refine_seeds(ink, fc, seeds, lines, tmpl, span=16, rad=20.0):
    """Re-seed each glyph by a 1 px grid search of the hard-edged model over +-16 px.

    The opening core of a marker that sits ON a frame line is the union of glyph, frame
    and tick, so its centroid can be ten pixels out -- and a local optimiser handed a seed
    that far off slides down the joining line into a wrong minimum that still looks
    converged.  A grid search cannot: it evaluates every position and takes the best.
    """
    out = []
    for i in range(len(seeds)):
        x0, y0 = seeds[i]
        h = 24 + span
        xs = np.arange(int(x0) - h, int(x0) + h + 1)
        ys = np.arange(int(y0) - h, int(y0) + h + 1)
        X, Y = np.meshgrid(xs, ys)
        patch = ink[ys[0] : ys[-1] + 1, xs[0] : xs[-1] + 1].astype(float)
        use = np.ones(X.shape, bool)
        if i in MASK:
            nm, w = MASK[i]
            use &= np.abs(X - np.polyval(fc[nm], Y)) > w
        dirs = []
        if i > 0:
            dirs.append(-lines[i - 1, 2:4])
        if i < len(seeds) - 1:
            dirs.append(lines[i, 2:4])
        best, bc = None, (x0, y0)
        for dx in range(-span, span + 1):
            for dy in range(-span, span + 1):
                c = (x0 + dx, y0 + dy)
                v = _hard_cost(patch, X, Y, use, c, tmpl, dirs, rad)
                if best is None or v < best:
                    best, bc = v, c
        out.append(bc)
    return np.array(out, float)


def glyph_pass(ink, fc, seeds, lines, mask_line=False, rad=20.0, framew=None, tmpl0=None):
    def mask_for(nm, w):
        def f(X, Y):
            return np.abs(X - np.polyval(fc[nm], Y)) > w

        return f

    W = []
    for i in range(N):
        mf = None
        if i in MASK:
            nm, w = MASK[i]
            mf = mask_for(nm, w if framew is None else framew)
        W.append(
            make_window(
                ink,
                seeds[i, 0],
                seeds[i, 1],
                rad=rad,
                half=24,
                mask_fn=mf,
                seg_in=(None if (mask_line or i == 0) else lines[i - 1, 2:4]),
                seg_out=(None if (mask_line or i == N - 1) else lines[i, 2:4]),
            )
        )
    if mask_line:  # the independent check: hide the polyline instead of modelling it
        for i, w in enumerate(W):
            for j, sign in ((i - 1, -1.0), (i, +1.0)):
                if 0 <= j < N - 1:
                    d = lines[j, 2:4] * sign
                    dx = w["X"] - seeds[i, 0]
                    dy = w["Y"] - seeds[i, 1]
                    t = dx * d[0] + dy * d[1]
                    p = -dx * d[1] + dy * d[0]
                    w["use"] &= ~((t > 6) & (np.abs(p) < 6))
    t0 = tmpl0 or TMPL0
    C, T = fit_glyphs(W, seeds, t0, seglen=42.0, simplex=2.0, verbose=False)
    moved = np.hypot(*(C - seeds).T)
    assert moved.max() < 4.0, f"a glyph ran away from its seed: {np.round(moved, 2)}"
    return C, T, moved


def count_check(curve, C):
    """A marker count that shares nothing with the detector: walk the ink along each fitted
    segment in 10 px bins and look for an excess over the joining line's own baseline."""
    ys, xs = np.nonzero(curve)
    P = np.stack([xs, ys], 1).astype(float)
    out = []
    for i in range(N - 1):
        a, b = C[i], C[i + 1]
        d = b - a
        L = float(np.hypot(*d))
        d = d / L
        t = (P - a) @ d
        perp = (P - a) @ np.array([-d[1], d[0]])
        tt = t[np.abs(perp) < 12]
        nb = max(int(L // 10), 4)
        h, e = np.histogram(tt, bins=nb, range=(0.0, L))
        base = max(np.median(h[2:-2]) if nb > 6 else np.median(h), 1)
        out.append(
            [
                round((e[k] + e[k + 1]) / 2, 1)
                for k in range(nb)
                if h[k] > 1.8 * base and 12 < (e[k] + e[k + 1]) / 2 < L - 12
            ]
        )
    return out


# --------------------------------------------------------------------------- output


def overlay(ink, to_px, C, px, py):
    im = Image.fromarray((~ink * 255).astype(np.uint8)).convert("RGB")
    d = ImageDraw.Draw(im)
    tt = np.linspace(0.0, 1.0, 40001)
    for k in range(16):
        u = float(tt[np.argmin(np.abs(np.polyval(px[0], tt) - (X0 + 0.01 * k)))])
        a, b = to_px(u, 0.0), to_px(u, 1.0)
        d.line([(a[0][0], a[1][0]), (b[0][0], b[1][0])], fill=(210, 230, 255), width=1)
    for k in range(14):
        v = float(tt[np.argmin(np.abs(np.polyval(py[0], tt) - (Y0 + 2.5 * k)))])
        a, b = to_px(0.0, v), to_px(1.0, v)
        d.line([(a[0][0], a[1][0]), (b[0][0], b[1][0])], fill=(210, 230, 255), width=1)
    for x, y in C:
        d.ellipse([x - 26, y - 26, x + 26, y + 26], outline=(220, 0, 0), width=3)
        d.line([(x - 40, y), (x + 40, y)], fill=(0, 170, 0), width=1)
        d.line([(x, y - 40), (x, y + 40)], fill=(0, 170, 0), width=1)
    im.crop((300, 250, 2250, 2750)).resize((975, 1250), Image.LANCZOS).save(OUT / "a7_verify.png")
    im.crop((430, 560, 880, 1450)).resize((900, 1780), Image.LANCZOS).save(
        OUT / "a7_verify_knee.png"
    )
    im.crop((2040, 2250, 2200, 2420)).resize((800, 850), Image.LANCZOS).save(
        OUT / "a7_verify_m6.png"
    )


def fit_anchored(t, val, order, lo, hi):
    """Least squares through the ticks that passes EXACTLY through both frame values.

    Why not `blend_axis`'s two extra data points: the frame value is *known* (printed, or
    fixed by the integer lattice span) and its pixel position comes from a ~1900-row fit
    of the frame line, so it is the best-determined datum on the axis -- while two copies
    among 76 ticks give it 2.6 % of the weight.  With the soft anchor the order-3 x map
    misses the left frame by 2.0 px, and marker 1 *sits on that frame*, so the miss lands
    straight in the CSV.  Fitted as `lo + (hi-lo)t + t(1-t)q(t)`, which is exact at both
    ends by construction and leaves order-1 free parameters for the bow.
    """
    base = lo + (hi - lo) * t
    if order < 2:
        return np.array([hi - lo, lo])
    A = np.stack([t * (1 - t) * t**j for j in range(order - 1)], 1)
    q, *_ = np.linalg.lstsq(A, val - base, rcond=None)
    asc = np.polynomial.polynomial.polyadd(
        [lo, hi - lo], np.polynomial.polynomial.polymul([0.0, 1.0, -1.0], q)
    )
    return asc[::-1]


def values(P, C, to_sv, xo=None, yo=None, blend=True, hard=True):
    xo = X_ORDER if xo is None else xo
    yo = Y_ORDER if yo is None else yo
    s, v = to_sv(C[:, 0], C[:, 1])
    fit = (lambda e, o, a, b: fit_anchored(P[e][0], P[e][1], o, a, b)) if hard else None
    if blend:
        if hard:
            px = [fit(e, xo, X0, X1) for e in ("bottom", "top")]
            py = [fit(e, yo, Y0, Y1) for e in ("left", "right")]
        else:
            px = blend_axis(P["bottom"][0], P["bottom"][1], P["top"][0], P["top"][1], xo, X0, X1)
            py = blend_axis(P["left"][0], P["left"][1], P["right"][0], P["right"][1], yo, Y0, Y1)
        X = (1 - v) * np.polyval(px[0], s) + v * np.polyval(px[1], s)
        Y = (1 - s) * np.polyval(py[0], v) + s * np.polyval(py[1], v)
    else:
        tx = np.concatenate([P["bottom"][0], P["top"][0]])
        vx = np.concatenate([P["bottom"][1], P["top"][1]])
        ty = np.concatenate([P["left"][0], P["right"][0]])
        vy = np.concatenate([P["left"][1], P["right"][1]])
        if hard:
            px = (fit_anchored(tx, vx, xo, X0, X1),) * 2
            py = (fit_anchored(ty, vy, yo, Y0, Y1),) * 2
        else:
            px = (np.polyfit(np.r_[tx, 0, 0, 1, 1], np.r_[vx, X0, X0, X1, X1], xo),) * 2
            py = (np.polyfit(np.r_[ty, 0, 0, 1, 1], np.r_[vy, Y0, Y0, Y1, Y1], yo),) * 2
        X = np.polyval(px[0], s)
        Y = np.polyval(py[0], v)
    return X, Y, px, py, s, v


def main() -> int:
    ink = native_bitmap(PAGE, OUT)
    tickset.ink = ink
    print("frame arcs (printed straight; the sagitta is the page's bow):")
    fc, _ = frame_curves(ink, SPEC, order=FRAME_ORDER)
    to_sv, to_px = make_sv(fc)

    print("\ncalibration ticks:")
    P = tickset(fc, to_sv, verbose=True)
    print("\n  FREE CHECK -- how many steps does the frame span, by a straight line through")
    print("  the ticks alone?  An integer says the unlabelled frame is a lattice point:")
    for nm, (a, b) in P["_span"].items():
        stepc, stepf = EDGES[nm][6], EDGES[nm][7]
        print(
            f"    {nm:6s}: {a:8.4f} majors of {stepc:<5g} (nearest {round(a)}, margin "
            f"{abs(a - round(a)):.3f})   {b:8.3f} minors of {stepf:<5g} (nearest {round(b)})"
        )
    print(f"  => the unlabelled bottom frame is {Y1} - 13*2.5 = {Y0}; the left frame {X0}.")

    print("\nCALIBRATION MODEL SELECTION -- both frame values held out, per edge")
    orders = {}
    for axis, lo, hi, unit, uname, edges in (
        ("y = f7 (BTU/LBM)", Y0, Y1, 1.0, "", ("left", "right")),
        ("x = P45/P41", X0, X1, 1e-4, "e-4", ("bottom", "top")),
    ):
        print(f"\n  --- {axis}: printed {hi} at the high frame, {lo} derived at the low frame ---")
        rows = {}
        for nm in edges:
            t, val = P[nm]
            r = held_out_table(t, val, lo, hi, unit, uname, label=f"{nm} edge, fine lattice ")
            rows[nm] = choose_order(r, lo, hi, unit, uname)
        tt = np.concatenate([P[edges[0]][0], P[edges[1]][0]])
        vv = np.concatenate([P[edges[0]][1], P[edges[1]][1]])
        pooled = held_out_table(tt, vv, lo, hi, unit, uname, label="both edges pooled, fine ")
        # The model actually used fits each edge ALONE and blends, so the per-edge test is
        # the one that matters; pooling is reported as a cross-check and, where the two
        # edges carry a real offset from each other (x here: 2.5 px at the right frame),
        # it cannot see curvature at all -- the same defect the A10 x note describes.
        pick = choose_order(pooled, lo, hi, unit, uname)
        orders[axis[0]] = (rows[edges[0]], rows[edges[1]], pick)
        Pc = tickset(fc, to_sv, fine=False)
        tt = np.concatenate([Pc[edges[0]][0], Pc[edges[1]][0]])
        vv = np.concatenate([Pc[edges[0]][1], Pc[edges[1]][1]])
        c = held_out_table(
            tt, vv, lo, hi, unit, uname, orders=(1, 2, 3, 4), label="both edges pooled, LONG ticks "
        )
        choose_order(c, lo, hi, unit, uname)
        print("   parity of the bow (Legendre coefficients of the value function, page px):")
        for nm in edges:
            t, val = P[nm]
            sc = 2340.0 / (hi - lo) if nm in ("left", "right") else 1658.0 / (hi - lo)
            row = " ".join(
                f"P{k}{'e' if k % 2 == 0 else 'o'} {c_:+6.2f}+-{s_:.2f}({z:4.1f}s)"
                for k, c_, s_, z in parity(t, val, sc)
            )
            print(f"     {nm:6s}: {row}")
    print(f"\n  per-edge / pooled choices: y {orders['y']}, x {orders['x']}")
    print(f"  chosen: frame arcs order {FRAME_ORDER}, x order {X_ORDER}, y order {Y_ORDER}")
    assert orders["y"][:2] == (Y_ORDER, Y_ORDER), orders["y"]

    print("\nmarkers -- morphological opening (a shape detector, independent of density):")
    cs = cores(ink, fc)
    for x, y, sz in cs:
        print(f"    core size {sz:5.0f} px at ({x:8.2f}, {y:8.2f})")
    assert len(cs) == N, f"opening found {len(cs)} cores, expected {N}"
    seeds = np.array([[x, y] for x, y, _ in cs])

    cur = curve_component(ink, fc)
    ys, xs = np.nonzero(cur)
    pts = np.stack([xs, ys], 1).astype(float)

    def clear_of_frame(i, p):
        if i == 0:
            return p[:, 0] - np.polyval(fc["left"], p[:, 1]) > 13.0
        if i == N - 2:
            return p[:, 0] - np.polyval(fc["right"], p[:, 1]) < -13.0
        return np.ones(p.shape[0], bool)

    lines = segment_lines(pts, seeds, clear=26.0, band=4.5, extra=clear_of_frame)
    print("\njoining segments (total least squares, 26 px clear of every marker):")
    for i, ln in enumerate(lines):
        print(
            f"  seg{i + 1}: len {ln[4]:7.1f} px  n {int(ln[5]):5d}  rms {ln[6]:.3f} px  "
            f"dir {np.degrees(np.arctan2(ln[3], ln[2])):+8.3f} deg"
        )

    print("\nglyph fit:")
    rs = refine_seeds(ink, fc, seeds, lines, TMPL0)
    print(f"  grid-refined seeds moved {np.round(np.hypot(*(rs - seeds).T), 1)} px from the cores")
    seeds = rs
    C, T, moved = glyph_pass(ink, fc, seeds, lines)
    print(
        f"  template: {np.degrees(T[0]):+.2f} deg / {T[1]:.2f} px and {np.degrees(T[2]):+.2f} deg"
        f" / {T[3]:.2f} px, glyph width {T[4]:.2f}, line width {T[5]:.2f}"
    )
    print(f"  centres moved {np.round(moved, 2)} px from their opening-core seeds")
    C2, _, _ = glyph_pass(ink, fc, seeds, lines, mask_line=True)
    dd = np.hypot(*(C - C2).T)
    print(
        f"  polyline masked out entirely instead of modelled: max {dd.max():.2f} px, "
        f"{np.round(dd, 2)}"
    )
    var = [dd]
    for kw in (dict(rad=17.0), dict(rad=23.0), dict(framew=5.5), dict(framew=8.5)):
        Cv, _, _ = glyph_pass(ink, fc, seeds, lines, **kw)
        var.append(np.hypot(*(C - Cv).T))
        print(f"  {kw}: max {var[-1].max():.2f} px")
    glyph_px = np.maximum(np.max(var, axis=0), 0.35)

    # distance from each centre to the independently fitted joining segments either side
    dperp = []
    for i in range(N):
        for j in (i - 1, i):
            if 0 <= j < N - 1:
                c0 = lines[j, :2]
                d0 = lines[j, 2:4]
                dperp.append(float((C[i] - c0) @ np.array([-d0[1], d0[0]])))
    dperp = np.array(dperp)
    print(
        f"  centre vs adjacent fitted segments: mean {dperp.mean():+.3f} px, "
        f"sd {dperp.std():.3f}, max {np.abs(dperp).max():.2f}"
    )
    print(
        f"  count check -- interior ink excesses per segment: {count_check(cur, C)} (empty => six)"
    )

    X, Y, px, py, s, v = values(P, C, to_sv)

    print("\nFREE STRUCTURAL CHECK -- the abscissae against the 0.005 lattice (never told):")
    lat = np.array([0.200, 0.215, 0.220, 0.225, 0.230, 0.350])
    for fo in (1, 2):
        fcv, _ = frame_curves(ink, SPEC, order=fo, verbose=False)
        sv, _ = make_sv(fcv)
        Pv = tickset(fcv, sv)
        for xo in (1, 3, 5):
            Xv, _, _, _, _, _ = values(Pv, C, sv, xo=xo)
            e = Xv - lat
            print(
                f"  frame order {fo}, x order {xo}: max {np.abs(e).max():.5f} "
                f"({np.abs(e).max() * 11053:.2f} px)  sd {e.std():.5f} "
                f"({e.std() * 11053:.2f} px)  mean {e.mean():+.5f}"
            )
    print(f"  chosen model: {np.round(X - lat, 5)}")

    print("\nwhat the ORDER is worth at the six points (same centres, same geometry):")
    for nm, orders_, fmt in (("y", (2, 3, 4, 5, 6), "8.4f"), ("x", (1, 2, 3, 4, 5), "9.6f")):
        for o in orders_:
            kw = dict(yo=o) if nm == "y" else dict(xo=o)
            Xo, Yo, _, _, _, _ = values(P, C, to_sv, **kw)
            row = Yo if nm == "y" else Xo
            print(f"   {nm} order {o}: " + " ".join(f"{z:{fmt}}" for z in row))

    print("\nmodel spread -- every defensible calibration, same glyph centres:")
    cache = {}
    VX, VY = [], []
    for fo, xo, yo, bl, hd in VARIANTS:
        if fo not in cache:
            fcv, _ = frame_curves(ink, SPEC, order=fo, verbose=False)
            sv, _ = make_sv(fcv)
            cache[fo] = (tickset(fcv, sv), sv)
        Pv, sv = cache[fo]
        Xv, Yv, _, _, _, _ = values(Pv, C, sv, xo=xo, yo=yo, blend=bl, hard=hd)
        VX.append(Xv)
        VY.append(Yv)
    VX, VY = np.array(VX), np.array(VY)
    mx, my = VX.std(0), VY.std(0)
    print(f"  all {len(VARIANTS)} variants: x sd max {mx.max():.5f}   y sd max {my.max():.4f}")
    ok = [i for i, m in enumerate(VARIANTS) if m[2] in (4, 5) and m[1] in (3, 5)]
    print(
        f"  subset the held-out tests pass ({len(ok)}: y order 4/5, x order 3/5): "
        f"x sd max {VX[ok].std(0).max():.5f}   y sd max {VY[ok].std(0).max():.4f}"
    )
    Pc = tickset(fc, to_sv, fine=False)
    Xc, Yc, _, _, _, _ = values(Pc, C, to_sv)
    print(
        f"  long ticks vs all ticks: max dx {np.abs(X - Xc).max():.5f}, "
        f"max dy {np.abs(Y - Yc).max():.4f}"
    )
    mx = np.hypot(mx, np.abs(X - Xc))
    my = np.hypot(my, np.abs(Y - Yc))

    # local two-tick interpolation: bow-free, because it never fits a bow
    tw = []
    for i in range(N):
        acc = []
        for nm in ("left", "right"):
            t, val = P[nm]
            o = np.argsort(t)
            t, val = t[o], val[o]
            j = min(max(int(np.searchsorted(t, v[i])), 1), t.size - 1)
            f = (v[i] - t[j - 1]) / (t[j] - t[j - 1])
            acc.append(val[j - 1] + f * (val[j] - val[j - 1]))
        tw.append(float((1 - s[i]) * acc[0] + s[i] * acc[1]))
    tw = np.array(tw)
    print(f"  chosen model minus local two-tick interpolation: {np.round(Y - tw, 4)}")

    # ---- per-point 1 sigma
    sig = {}
    for nm, poly, order, arg in (
        ("left", py[0], Y_ORDER, v),
        ("right", py[1], Y_ORDER, v),
        ("bottom", px[0], X_ORDER, s),
        ("top", px[1], X_ORDER, s),
    ):
        t, val = P[nm]
        lo, hi = (Y0, Y1) if nm in ("left", "right") else (X0, X1)
        tt = np.concatenate([t, [0.0, 1.0]])
        vv = np.concatenate([val, [lo, hi]])
        sig[nm] = poly_se(tt, (vv - np.polyval(poly, tt)).std(), order, arg)
    cal_y = np.hypot((1 - s) * sig["left"], s * sig["right"])
    cal_x = np.hypot((1 - v) * sig["bottom"], v * sig["top"])
    hpx = np.polyval(fc["right"], C[:, 1]) - np.polyval(fc["left"], C[:, 1])
    vpx = np.polyval(fc["bottom"], C[:, 0]) - np.polyval(fc["top"], C[:, 0])
    dXdpx = np.abs(np.polyval(np.polyder(px[0]), s)) / hpx
    dYdpx = np.abs(np.polyval(np.polyder(py[0]), v)) / vpx
    sx = np.sqrt(cal_x**2 + (glyph_px * dXdpx) ** 2 + mx**2)
    sy = np.sqrt(cal_y**2 + (glyph_px * dYdpx) ** 2 + my**2)

    overlay(ink, to_px, C, px, py)

    old = np.array(
        [
            [0.199937, 48.2815],
            [0.214923, 40.3985],
            [0.219917, 39.1268],
            [0.224969, 38.3672],
            [0.229970, 37.5898],
            [0.350096, 24.6064],
        ]
    )
    print(
        "\n   k  px_x     px_y       x          y      1s_x    1s_y  |  cal_x glyph_x model_x |"
        "  cal_y glyph_y model_y |   dx      dy    (sigma)"
    )
    for i in range(N):
        print(
            f"  {i + 1:2d} {C[i, 0]:8.2f} {C[i, 1]:8.2f} {X[i]:9.6f} {Y[i]:8.4f} "
            f"{sx[i]:7.5f} {sy[i]:6.4f} | {cal_x[i]:6.5f} {glyph_px[i] * dXdpx[i]:7.5f} "
            f"{mx[i]:7.5f} | {cal_y[i]:6.4f} {glyph_px[i] * dYdpx[i]:7.4f} {my[i]:7.4f} | "
            f"{X[i] - old[i, 0]:+.6f} {Y[i] - old[i, 1]:+7.4f} "
            f"({(X[i] - old[i, 0]) / sx[i]:+.1f}, {(Y[i] - old[i, 1]) / sy[i]:+.1f})"
        )
    # The CSV is maintained by hand; these values are not. See `check_against_csv`.
    return check_against_csv(
        CSV,
        {
            "x": (X, "%.6f"),
            "y": (Y, "%.4f"),
            "sx": (sx, "%.6f"),
            "sy": (sy, "%.4f"),
        },
        label="f7 (Figure A7)",
    )


if __name__ == "__main__":
    raise SystemExit(main())
