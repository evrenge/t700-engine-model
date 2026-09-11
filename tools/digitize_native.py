"""Shared machinery for digitizing an Appendix A figure off the native 300 dpi scan.

This is `tools/digitize_a7.py`'s method, factored out so A2 and A9 can reuse it, plus the
two things A7 taught us afterwards:

* **The axis map is a polynomial and its order is chosen by a held-out test.**  Where a
  figure prints a value at a frame edge, fit the printed ticks *excluding* that value and
  predict it.  Linear vs quadratic vs cubic: whichever predicts the printed frame value
  best is the right model.  That is a held-out test, not a fit statistic.  A second,
  independent test is `tail_cv` -- fit the middle of the tick lattice and predict the two
  ends -- which needs no printed frame value at all.
* **The four frame edges are curved, not straight.**  A printed straight edge images with
  a 0.3-1.8 px sagitta on these pages, which is the same size as the errors being chased.
  So the plot coordinates are a curvilinear blend between the four *fitted arcs*
  (`make_sv`), not a homography between four corners.  On a flatbed scan there is no
  perspective to model; the defect is the paper's shape.

Nothing here resamples the page.  Rotation, keystone and bow are all carried as a map
applied to the measured points.

This module lives in tools/ and may use SciPy and PIL.  Nothing here is imported by
src/t700/.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.optimize import minimize

PDF = Path("docs/ballin-tm100991.pdf")


# --------------------------------------------------------------------------- the raster


def native_bitmap(page: int, out: Path) -> np.ndarray:
    """The page's own 300 dpi 1-bit image, as an ink mask.  No render, no resample."""
    out.mkdir(parents=True, exist_ok=True)
    if not (out / "native-000.png").exists():
        subprocess.run(
            ["pdfimages", "-f", str(page), "-l", str(page), "-png", str(PDF), str(out / "native")],
            check=True,
        )
    return ~np.array(Image.open(out / "native-000.png"))


# --------------------------------------------------------------------------- frame edges


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    e = np.flatnonzero(np.diff(np.concatenate(([0], mask.view(np.int8), [0]))))
    return list(zip(e[0::2], e[1::2], strict=True))


def _trace(ink, along, lo, hi, tlo, thi, maxw=8):
    """Centroid of the single narrow ink run in [lo,hi) at each position along the line.

    `along=0` traces a vertical edge (one run per row), `along=1` a horizontal one.
    A frame line gives exactly one run; text and data give several, and are skipped.
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


def _robust_line(t, v, tol=2.0, it=8):
    keep = np.ones(t.shape, bool)
    for _ in range(it):
        s, i = np.polyfit(t[keep], v[keep], 1)
        new = np.abs(v - (i + s * t)) < tol
        if new.sum() < 20 or (new == keep).all():
            break
        keep = new
    return keep


def frame_curves(ink, spec, order=2, verbose=True):
    """Each frame edge as a polynomial in the coordinate along it.

    Order 1 is the straight-edge model; order 2 admits the page's bow.  The sagitta
    printed here is the whole reason order 2 exists -- it is a *printed straight line*
    imaging as an arc, so it measures the page distortion directly, with no reliance on
    the tick lattice.
    """
    out, diag = {}, {}
    for nm, along, args in spec:
        t, v = _trace(ink, along, *args)
        keep = _robust_line(t, v)
        c = np.polyfit(t[keep], v[keep], order)
        r = v[keep] - np.polyval(c, t[keep])
        lo, hi = t[keep].min(), t[keep].max()
        mid = 0.5 * (lo + hi)
        sag = np.polyval(c, mid) - 0.5 * (np.polyval(c, lo) + np.polyval(c, hi))
        deg = np.degrees(np.arctan(np.polyval(np.polyder(c), mid)))
        out[nm] = c
        diag[nm] = dict(
            rms=float(r.std()), sagitta=float(sag), slope_deg=float(deg), n=int(keep.sum())
        )
        if verbose:
            print(
                f"  {nm:6s} n {keep.sum():5d}  rms {r.std():.3f} px   "
                f"sagitta {sag:+.2f} px   mid-slope {deg:+.3f} deg"
            )
    return out, diag


def make_sv(fc):
    """Curvilinear plot coordinates from the four fitted frame arcs.

    `s` is 0 on the left edge and 1 on the right, `v` is 0 on the bottom and 1 on the
    top -- exactly, by construction -- with a linear blend between.  This is the
    curved-edge generalisation of `digitize_a7.py`'s homography and it is the only piece
    of geometry ever applied to a datum.
    """

    def to_sv(px, py):
        px = np.atleast_1d(np.asarray(px, float))
        py = np.atleast_1d(np.asarray(py, float))
        xl = np.polyval(fc["left"], py)
        xr = np.polyval(fc["right"], py)
        yt = np.polyval(fc["top"], px)
        yb = np.polyval(fc["bottom"], px)
        return (px - xl) / (xr - xl), 1.0 - (py - yt) / (yb - yt)

    def to_px(s, v, iters=8):
        s = np.atleast_1d(np.asarray(s, float))
        v = np.atleast_1d(np.asarray(v, float))
        py = np.full(s.shape, float(np.polyval(fc["left"], 0.0)) * 0.0 + 1500.0)
        px = np.empty_like(py)
        for _ in range(iters):
            xl = np.polyval(fc["left"], py)
            xr = np.polyval(fc["right"], py)
            px = xl + s * (xr - xl)
            yt = np.polyval(fc["top"], px)
            yb = np.polyval(fc["bottom"], px)
            py = yt + (1.0 - v) * (yb - yt)
        return px, py

    return to_sv, to_px


def corner(fc, hor, ver, guess):
    """Intersection of two fitted frame arcs, by fixed-point iteration."""
    x, y = guess
    for _ in range(30):
        y = np.polyval(fc[hor], x)
        x = np.polyval(fc[ver], y)
    return float(x), float(y)


# --------------------------------------------------------------------------- ticks


def edge_run(ink, edge_poly, along, sign, lo, hi):
    """Length of the ink run through the frame edge, per position along it.

    A tick lengthens the run.  The frame's own width need not be known: it is removed
    afterwards as a rolling median, so a line that thickens or thins along the page
    cannot masquerade as a tick (which is what a fixed inward-probe offset lets happen).
    """
    n = hi - lo
    ln = np.zeros(n)
    H, W = ink.shape

    def on(t, q):
        if along == 0:
            return 0 <= q < W and ink[t, q]
        return 0 <= q < H and ink[q, t]

    for j, t in enumerate(range(lo, hi)):
        base = int(round(np.polyval(edge_poly, t)))
        if not on(t, base):
            for dq in (1, -1, 2, -2, 3, -3):
                if on(t, base + dq):
                    base += dq
                    break
            else:
                ln[j] = 0
                continue
        p = base
        while on(t, p + sign):
            p += sign
        m = base
        while on(t, m - sign):
            m -= sign
        ln[j] = abs(p - m) + 1
    return ln


def _rolling_median(x, half):
    return np.array(
        [np.median(x[max(0, i - half) : min(x.size, i + half + 1)]) for i in range(x.size)]
    )


def ticks(ink, edge_poly, along, sign, lo, hi, thresh, half=60, gap=2, maxw=14):
    """Tick centres along a frame edge: ink-weighted centroid of each run excess."""
    ln = edge_run(ink, edge_poly, along, sign, lo, hi)
    ex = ln - _rolling_median(ln, half)
    idx = np.flatnonzero(ex >= thresh)
    if idx.size == 0:
        return np.zeros(0)
    groups, cur = [], [int(idx[0])]
    for v in idx[1:]:
        if v - cur[-1] <= gap:
            cur.append(int(v))
        else:
            groups.append(cur)
            cur = [int(v)]
    groups.append(cur)
    out = []
    for g in groups:
        g = np.asarray(g)
        if (g[-1] - g[0] + 1) > maxw:
            continue
        w = ex[g]
        out.append(float(lo + (g * w).sum() / w.sum()))
    return np.array(out)


def label_ticks(w, w_lo, w_hi, val_lo, val_hi, step, tol=0.30):
    """Attach a printed value to each detected tick, from the two frame-edge values.

    The frame values fix only a *discrete* choice -- which lattice slot a tick occupies --
    so they cannot leak into the continuous fit that is later tested against them.  A tick
    landing further than `tol` steps from every slot (a marker's arm on the frame, a corner
    artefact) is dropped, and a slot claimed twice keeps the closer tick.

    This replaces counting.  On A2 a run of the detector that returned exactly the
    expected number of ticks contained one spurious entry and was missing one real one,
    and the resulting off-by-one labelling moved the axis by 0.03 -- invisible in the
    count, invisible in the CSV.  Assert on structure, never on count.
    """
    w = np.asarray(w, float)
    est = val_lo + (w - w_lo) * (val_hi - val_lo) / (w_hi - w_lo)
    k = np.round(est / step)
    d = np.abs(est - k * step) / step
    keep = d <= tol
    seen: dict[int, int] = {}
    for i in np.flatnonzero(keep):
        kk = int(k[i])
        if kk in seen:
            if d[i] < d[seen[kk]]:
                keep[seen[kk]] = False
                seen[kk] = i
            else:
                keep[i] = False
        else:
            seen[kk] = i
    return w[keep], k[keep] * step, keep


def label_by_map(t, poly, step, tol=0.30):
    """Label a fine tick lattice through a coarse-lattice fit, not through a straight line.

    The straight line between the two frame corners is off by several pixels in mid-page
    (that is the whole point of this exercise), which on a 16.7 px minor lattice throws
    away a contiguous block of perfectly good ticks -- 29 of A9's 141 -- and would throw
    away the wrong ones if the tolerance were widened instead.
    """
    est = np.polyval(poly, t)
    k = np.round(est / step)
    d = np.abs(est - k * step) / step
    keep = d <= tol
    seen: dict[int, int] = {}
    for i in np.flatnonzero(keep):
        kk = int(k[i])
        if kk in seen:
            if d[i] < d[seen[kk]]:
                keep[seen[kk]] = False
                seen[kk] = i
            else:
                keep[i] = False
        else:
            seen[kk] = i
    return t[keep], k[keep] * step, keep


# --------------------------------------------------------------------------- axis model


def held_out(t, val, edge_lo, edge_hi, orders=(1, 2, 3, 4), unit=1.0, uname=""):
    """Fit on interior ticks only, then predict the printed value at each frame edge."""
    rows = []
    for o in orders:
        c = np.polyfit(t, val, o)
        r = val - np.polyval(c, t)
        a, b = float(np.polyval(c, 0.0)), float(np.polyval(c, 1.0))
        rows.append((o, float(r.std()), a, b))
        ea = "     -    " if edge_lo is None else f"{a - edge_lo:+9.5f}"
        eb = "     -    " if edge_hi is None else f"{b - edge_hi:+9.5f}"
        print(
            f"    order {o}: interior rms {r.std() / unit:8.3f}{uname}   "
            f"edge0 {a:11.6f} ({ea})   edge1 {b:11.6f} ({eb})"
        )
    return rows


def tail_cv(t, val, frac=0.18, orders=(1, 2, 3, 4), unit=1.0, uname=""):
    """Fit the middle of the tick lattice, predict the two ends.  Uses no printed edge."""
    lo, hi = t.min(), t.max()
    span = hi - lo
    inner = (t > lo + frac * span) & (t < hi - frac * span)
    for o in orders:
        c = np.polyfit(t[inner], val[inner], o)
        r = val[~inner] - np.polyval(c, t[~inner])
        print(
            f"    order {o}: held-out rms {np.sqrt((r**2).mean()) / unit:8.3f}{uname}  "
            f"max {np.abs(r).max() / unit:8.3f}{uname}   (fit {inner.sum()}, test {(~inner).sum()})"
        )


def blend_axis(t_lo, v_lo, t_hi, v_hi, order, edge0, edge1):
    """Two opposite edges calibrated separately, each anchored on its printed frame values.

    Returned as a pair of polynomials; the caller blends them across the plot.  Fitting
    the two edges *together* would average away a real difference: on A2 the top and
    bottom x lattices disagree by 0.026 at the left frame.
    """
    p_lo = np.polyfit(
        np.concatenate([t_lo, [0.0, 1.0]]), np.concatenate([v_lo, [edge0, edge1]]), order
    )
    p_hi = np.polyfit(
        np.concatenate([t_hi, [0.0, 1.0]]), np.concatenate([v_hi, [edge0, edge1]]), order
    )
    return p_lo, p_hi


def axis_values(ink, fc_spec, tickspec, C, frame_order, x_order, y_order, blend=True):
    """Data values for a set of pixel centres under one complete calibration model.

    Exists so the *model* can be varied -- frame arcs straight or curved, each axis
    linear or cubic, opposite edges blended or pooled -- and the spread over defensible
    variants carried into the uncertainty as a systematic, instead of a single model's
    formal error being quoted as if the model were certain.
    """
    fc, _ = frame_curves(ink, fc_spec, order=frame_order, verbose=False)
    to_sv, _ = make_sv(fc)
    P = tickspec(fc, to_sv)
    s, v = to_sv(C[:, 0], C[:, 1])
    out = []
    for order, (nlo, nhi), (elo, ehi), arg, other in (
        (x_order, ("bottom", "top"), P["_xedges"], s, v),
        (y_order, ("left", "right"), P["_yedges"], v, s),
    ):
        if blend:
            p_lo, p_hi = blend_axis(P[nlo][0], P[nlo][1], P[nhi][0], P[nhi][1], order, elo, ehi)
            out.append((1 - other) * np.polyval(p_lo, arg) + other * np.polyval(p_hi, arg))
        else:
            t = np.concatenate([P[nlo][0], P[nhi][0], [0.0, 0.0, 1.0, 1.0]])
            w = np.concatenate([P[nlo][1], P[nhi][1], [elo, elo, ehi, ehi]])
            out.append(np.polyval(np.polyfit(t, w, order), arg))
    return out[0], out[1]


def poly_se(t_fit, resid_std, order, t_at):
    """1 sigma standard error of a polynomial fit's prediction, at `t_at`."""
    A = np.vander(np.atleast_1d(t_fit), order + 1)
    cov = np.linalg.inv(A.T @ A) * resid_std**2
    B = np.vander(np.atleast_1d(t_at), order + 1)
    return np.sqrt(np.einsum("ij,jk,ik->i", B, cov, B))


# --------------------------------------------------------------------------- glyph model

SS = 5  # supersampling per axis, for a sub-pixel coverage model


def make_window(ink, cx, cy, rad=20.0, half=24, mask_fn=None, seg_in=None, seg_out=None):
    xs = np.arange(int(cx) - half, int(cx) + half + 1)
    ys = np.arange(int(cy) - half, int(cy) + half + 1)
    X, Y = np.meshgrid(xs, ys)
    use = ((X - cx) ** 2 + (Y - cy) ** 2) < rad**2
    if mask_fn is not None:
        use &= mask_fn(X, Y)
    off = (np.arange(SS) + 0.5) / SS - 0.5
    return dict(
        X=X,
        Y=Y,
        use=use,
        data=ink[ys[0] : ys[-1] + 1, xs[0] : xs[-1] + 1].astype(float),
        OX=X[..., None, None] + off[None, None, :, None],
        OY=Y[..., None, None] + off[None, None, None, :],
        seg_in=seg_in,
        seg_out=seg_out,
    )


def _bar(w, c, ang, half, wid, halfline=False):
    u = np.array([np.cos(ang), np.sin(ang)])
    n = np.array([-u[1], u[0]])
    dx = w["OX"] - c[0]
    dy = w["OY"] - c[1]
    t = dx * u[0] + dy * u[1]
    p = dx * n[0] + dy * n[1]
    m = np.abs(p) <= wid / 2
    m &= (t >= 0) & (t <= half) if halfline else (np.abs(t) <= half)
    return m


def coverage(w, c, th1, l1, th2, l2, wg, wl, seglen=42.0):
    """`stroke1 u stroke2 u incoming half-segment u outgoing half-segment`, all radiating
    from one point, supersampled to a sub-pixel area fraction."""
    m = _bar(w, c, th1, l1, wg) | _bar(w, c, th2, l2, wg)
    if w["seg_in"] is not None:
        m |= _bar(w, c, np.arctan2(-w["seg_in"][1], -w["seg_in"][0]), seglen, wl, True)
    if w["seg_out"] is not None:
        m |= _bar(w, c, np.arctan2(w["seg_out"][1], w["seg_out"][0]), seglen, wl, True)
    return m.mean(axis=(2, 3))


def _cost_one(c, w, T, seglen):
    return float((((coverage(w, c, *T, seglen=seglen) - w["data"]) ** 2)[w["use"]]).sum())


def fit_glyphs(W, C0, tmpl0, iters=6, seglen=42.0, verbose=True, simplex=None):
    """Alternate: every centre against a fixed template, then the template against all
    centres.  The template is global -- the plotter drew the same mark at every point --
    and the score covers background as well as ink, which is what stops a near-parallel
    pair of lines from outscoring a four-armed glyph.

    `simplex`, if given, is the side of an explicit initial simplex in PIXELS and the
    centre is optimised as an OFFSET from its current value.  Nelder-Mead's default
    simplex is 5 % of each coordinate, which at an image column of 2118 is a 106 px first
    step -- outside any scoring window, where the cost is flat and converged-looking.
    Pass it (2.0 is right for these glyphs) on any figure with a marker far from the
    origin; the default `None` keeps the historical behaviour bit for bit.
    """
    C = np.array(C0, float)
    T = np.array(tmpl0, float)
    for it in range(iters):
        for i, w in enumerate(W):
            if simplex is None:
                r = minimize(
                    _cost_one,
                    C[i],
                    args=(w, T, seglen),
                    method="Nelder-Mead",
                    options=dict(xatol=1e-4, fatol=1e-8, maxiter=4000),
                )
                C[i] = r.x
                continue
            base = C[i].copy()
            r = minimize(
                lambda d, b=base, w=w, T=T: _cost_one(b + d, w, T, seglen),
                np.zeros(2),
                method="Nelder-Mead",
                options=dict(
                    xatol=1e-4,
                    fatol=1e-8,
                    maxiter=4000,
                    initial_simplex=np.array([[0.0, 0.0], [simplex, 0.0], [0.0, simplex]], float),
                ),
            )
            C[i] = base + r.x
        r = minimize(
            lambda t: float(sum(_cost_one(C[i], W[i], t, seglen) for i in range(len(W)))),
            T,
            method="Powell",
            options=dict(xtol=1e-5, ftol=1e-8, maxiter=20000, maxfev=60000),
        )
        T = r.x
        if verbose:
            print(
                f"   iter {it}: cost {r.fun:8.3f}   strokes {np.degrees(T[0]):+.2f} deg / "
                f"{T[1]:.2f} px and {np.degrees(T[2]):+.2f} deg / {T[3]:.2f} px, "
                f"glyph width {T[4]:.2f}, line width {T[5]:.2f}"
            )
    return C, T


def curve_component(ink, fc, band=4.0, margin=30.0, close=0):
    """The joining polyline and its markers, as one 8-connected blob with the frame off.

    `close` bridges a joining line the scan broke -- A9's is 1 px wide over its flat run
    and comes apart in a dozen places, which splits the curve into fragments and makes
    "largest component" pick a fragment.  Closing is used only to decide *membership*;
    the returned mask is the original ink.
    """
    Y, X = np.mgrid[0 : ink.shape[0], 0 : ink.shape[1]]
    dL = X - np.polyval(fc["left"], Y)
    dR = X - np.polyval(fc["right"], Y)
    dT = Y - np.polyval(fc["top"], X)
    dB = Y - np.polyval(fc["bottom"], X)
    inside = (dL > -margin) & (dR < margin) & (dT > -margin) & (dB < margin)
    on_frame = (np.abs(dL) < band) | (np.abs(dR) < band) | (np.abs(dT) < band) | (np.abs(dB) < band)
    m = ink & inside & ~on_frame
    g = ndimage.binary_closing(m, np.ones((close, close))) if close else m
    lab, n = ndimage.label(g, structure=np.ones((3, 3)))
    sz = ndimage.sum(g, lab, range(1, n + 1))
    return m & (lab == int(np.argmax(sz)) + 1)


def segment_lines(curve_pts, seeds, clear=24.0, band=4.5, extra=None):
    """The joining polyline, segment by segment, total least squares, clear of markers."""
    lines = []
    for i in range(len(seeds) - 1):
        a, b = np.asarray(seeds[i], float), np.asarray(seeds[i + 1], float)
        d = b - a
        L = float(np.hypot(*d))
        d = d / L
        n = np.array([-d[1], d[0]])
        t = (curve_pts - a) @ d
        perp = (curve_pts - a) @ n
        sel = (np.abs(perp) < band) & (t > clear) & (t < L - clear)
        if extra is not None:
            sel &= extra(i, curve_pts)
        pts = curve_pts[sel]
        c = pts.mean(0)
        _, _, vt = np.linalg.svd(pts - c, full_matrices=False)
        dr = vt[0]
        if dr @ d < 0:
            dr = -dr
        r = (pts - c) @ np.array([-dr[1], dr[0]])
        lines.append([c[0], c[1], dr[0], dr[1], L, int(sel.sum()), float(r.std())])
    return np.array(lines)
