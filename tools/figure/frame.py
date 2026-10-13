"""The four printed frame edges, and the map from pixels to plot coordinates.

A frame located wrongly moves every value on the page without disturbing any check that
looks only at the data. Two geometry models (recipe, trap 7):

* order 1 -- straight edges, plot coordinates by `homography` between the four corners;
* order 2 -- each edge an arc, plot coordinates by `make_sv`'s curvilinear blend.

The single-curve figures choose between them by their own held-out frame test; A1 uses
order 1, and A2, A7 and A9 order 2. The body figures 6-10 use `make_sv` on four straight
lines (`fit_edges` on 6-8, `band_line` on 9-10). Nothing here resamples the page: tilt and bow
are carried as a map applied to the measured points.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from .raster import native_bitmap

EDGE = {  # name -> (along, inward sign)
    "left": (0, +1),
    "right": (0, -1),
    "top": (1, +1),
    "bottom": (1, -1),
}


# --------------------------------------------------------------------------- tracing


def runs(mask):
    """(start, end) of every contiguous True run in a 1-D boolean array; `end` is exclusive."""
    e = np.flatnonzero(np.diff(np.concatenate(([0], np.asarray(mask).view(np.int8), [0]))))
    return list(zip(e[0::2], e[1::2], strict=True))


def trace(ink, along, lo, hi, tlo, thi, maxw=8):
    """Centroid of the single narrow ink run in [lo,hi) at each position along the line.

    `along=0` traces a vertical edge (one run per row), `along=1` a horizontal one. A bare
    frame line gives exactly one narrow run; a tick merges with it and widens the run past
    `maxw`, and data or text give several -- both are skipped rather than modelled.
    """
    ts, vs = [], []
    for t in range(tlo, thi):
        seg = ink[t, lo:hi] if along == 0 else ink[lo:hi, t]
        rr = [r for r in runs(seg) if (r[1] - r[0]) <= maxw]
        if len(rr) != 1:
            continue
        a, b = rr[0]
        ts.append(float(t))
        vs.append(lo + (a + b - 1) / 2.0)
    return np.array(ts), np.array(vs)


def robust_keep(t, v, tol=2.0, it=8):
    """Which traced points belong to the line: iterated straight fit, |resid| < `tol`."""
    keep = np.ones(t.shape, bool)
    for k in range(it):
        s, i = np.polyfit(t[keep], v[keep], 1)
        new = np.abs(v - (i + s * t)) < tol
        if new.sum() < 20 and k == 0:
            # keeping `keep` here would fit the whole contaminated trace
            raise AssertionError(
                f"only {int(new.sum())} of {t.size} traced points lie within {tol} px of a "
                f"straight fit -- the trace is not one line"
            )
        if new.sum() < 20 or (new == keep).all():
            break
        keep = new
    return keep


def fit_edges(ink, spec, order=1, verbose=True):
    """Each frame edge as a polynomial in the coordinate along it, highest power first.

    `spec` is a list of `(name, along, xlo, xhi, tlo, thi)` windows (see `auto_frame_spec`).
    Returns `(coeffs, diag)`. The order-2 sagitta in `diag` measures the page distortion
    directly, since the ink is a printed straight line.
    """
    out, diag = {}, {}
    for nm, along, xlo, xhi, tlo, thi in spec:
        t, v = trace(ink, along, xlo, xhi, tlo, thi)
        keep = robust_keep(t, v)
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


def as_lines(fc, span=None):
    """Each edge as an `(intercept, slope)` pair, for the masks and the seed search.

    Exact for an order-1 fit. For a higher order pass `span`, `(t0, t1)` per edge, to get
    the best straight line over the edge: without it the pair is the tangent at t = 0, and
    with t a pixel row running to 2600 that tangent leaves the ink entirely.
    """
    out = {}
    for nm, c in fc.items():
        c = np.atleast_1d(np.asarray(c, float))
        if c.size <= 2 or span is None:
            out[nm] = (float(c[-1]), float(c[-2]))
            continue
        t0, t1 = sorted(span[nm])
        t = np.linspace(t0, t1, 200)
        b, a = np.polyfit(t, np.polyval(c, t), 1)
        out[nm] = (float(a), float(b))
    return out


# --------------------------------------------------------------------------- corners


def corner_line(hor, ver):
    """Intersection of two straight edges, closed form. `hor`/`ver` are (intercept, slope)."""
    a, b = hor
    c, d = ver
    x = (c + d * a) / (1.0 - d * b)
    return x, a + b * x


def corner_arc(fc, hor, ver, guess):
    """Intersection of two fitted frame arcs, by fixed-point iteration."""
    x, y = guess
    for _ in range(30):
        y = np.polyval(fc[hor], x)
        x = np.polyval(fc[ver], y)
    return float(x), float(y)


def corners(fc, order=1):
    """The four corners as (TL, TR, BR, BL), by whichever formula the model admits."""
    if order == 1:
        if any(np.atleast_1d(c).size > 2 for c in fc.values()):
            # `as_lines` would return the tangent at t = 0, far off the ink
            raise AssertionError("corners(order=1) was handed a higher-order frame fit")
        fr = as_lines(fc)
        return (
            corner_line(fr["top"], fr["left"]),
            corner_line(fr["top"], fr["right"]),
            corner_line(fr["bottom"], fr["right"]),
            corner_line(fr["bottom"], fr["left"]),
        )
    # seeded from the intersection of the edges' linear parts, a pixel or two off the arcs
    lin = {nm: (float(c[-1]), float(c[-2])) for nm, c in fc.items()}
    return tuple(
        corner_arc(fc, h, v, corner_line(lin[h], lin[v]))
        for h, v in (("top", "left"), ("top", "right"), ("bottom", "right"), ("bottom", "left"))
    )


def edge_spans(corners4):
    """Each edge's extent along itself, between its own two corners `(TL, TR, BR, BL)`: rows
    for the left and right edges, columns for the top and bottom."""
    TL, TR, BR, BL = corners4
    return {
        "left": (TL[1], BL[1]),
        "right": (TR[1], BR[1]),
        "top": (TL[0], TR[0]),
        "bottom": (BL[0], BR[0]),
    }


# --------------------------------------------------------------------------- plot coords


UNIT = ((0, 1), (1, 1), (1, 0), (0, 0))
"""The unit square's corners in `(TL, TR, BR, BL)` order: u runs right, v runs up."""


def homography(src, dst):
    """The projective map taking the four `src` points to the four `dst` points."""
    A, B = [], []
    for (sx, sy), (dx, dy) in zip(src, dst, strict=True):
        A.append([sx, sy, 1, 0, 0, 0, -dx * sx, -dx * sy])
        B.append(dx)
        A.append([0, 0, 0, sx, sy, 1, -dy * sx, -dy * sy])
        B.append(dy)
    h = np.linalg.solve(np.asarray(A, float), np.asarray(B, float))
    return np.append(h, 1.0).reshape(3, 3)


def apply_h(H, px, py):
    """Apply a homography to pixel coordinates."""
    px = np.atleast_1d(np.asarray(px, float))
    py = np.atleast_1d(np.asarray(py, float))
    d = H[2, 0] * px + H[2, 1] * py + H[2, 2]
    return (
        (H[0, 0] * px + H[0, 1] * py + H[0, 2]) / d,
        (H[1, 0] * px + H[1, 1] * py + H[1, 2]) / d,
    )


def make_sv(fc):
    """Curvilinear plot coordinates from the four fitted frame arcs, as `(to_sv, to_px)`.

    `s` is 0 on the left edge and 1 on the right, `v` 0 on the bottom and 1 on the top,
    exactly by construction, with a linear blend between: the curved-edge generalisation of
    the homography.
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
        py = np.full(s.shape, 1500.0)
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


# --------------------------------------------------------------------------- masks


def dist(edge, X, Y, vertical=True):
    """Signed perpendicular distance from a straight edge `(intercept, slope)`."""
    a, b = edge
    return (
        (X - (a + b * Y)) / np.hypot(1.0, b) if vertical else (Y - (a + b * X)) / np.hypot(1.0, b)
    )


def blocked_mask(fr, allt, framew=3.6, tickpad=4.5, tickext=5.0):
    """A `mask_fn` excluding frame lines and tick marks, SYMMETRICALLY about each edge.

    Ticks protrude inward only, so masking just the inward side would strip half of a marker
    on a frame line and pull its fitted centre outward (recipe, trap 10).
    """

    def mask_fn(X, Y):
        m = np.zeros(X.shape, bool)
        for nm in ("left", "right", "top", "bottom"):
            vert = nm in ("left", "right")
            d = dist(fr[nm], X, Y, vert)
            t = Y if vert else X
            m |= np.abs(d) < framew
            for c, k, _ in allt[nm]:
                if np.abs(t - c).min() > tickpad + 3:
                    continue
                m |= (np.abs(t - c) < tickpad) & (np.abs(d) < k + tickext)
        return m

    return mask_fn


def frame_offsets(fc, X, Y):
    """Each pixel's offset from the four fitted edges, measured across each one: x minus the
    left or right edge's x at that row, y minus the top or bottom edge's y at that column.
    Inside the frame `dL` and `dT` are positive, `dR` and `dB` negative."""
    return (
        X - np.polyval(fc["left"], Y),
        X - np.polyval(fc["right"], Y),
        Y - np.polyval(fc["top"], X),
        Y - np.polyval(fc["bottom"], X),
    )


def inside_frame(offsets, margin):
    """Inside the frame, or within `margin` px outside it, from `frame_offsets`."""
    dL, dR, dT, dB = offsets
    return (dL > -margin) & (dR < margin) & (dT > -margin) & (dB < margin)


def plot_ink(ink, fc, band, margin):
    """Ink inside the frame (to `margin` px beyond it) with the frame lines themselves,
    everything within `band` px of a fitted edge, removed."""
    Y, X = np.mgrid[0 : ink.shape[0], 0 : ink.shape[1]]
    d = frame_offsets(fc, X, Y)
    dL, dR, dT, dB = d
    on_frame = (np.abs(dL) < band) | (np.abs(dR) < band) | (np.abs(dT) < band) | (np.abs(dB) < band)
    return ink & inside_frame(d, margin) & ~on_frame


# --------------------------------------------------------------------------- finding it


def _sheared_row_profile(sub, slope):
    """Row-sum of `sub` after shearing each column to undo a line of the given slope."""
    h, w = sub.shape
    xc = w / 2.0
    shift = np.rint(slope * (np.arange(w) - xc)).astype(int)
    ys = np.arange(h)[:, None] + shift[None, :]
    ok = (ys >= 0) & (ys < h)
    return np.where(ok, sub[np.clip(ys, 0, h - 1), np.arange(w)[None, :]], False).sum(axis=1)


def find_frame_box(ink, margin=150, sep_frac=0.25, slopes=None):
    """Locate the printed plot box as (left, right, top, bottom) pixel centres.

    The horizontal frames are tilted and broken: A11's bottom frame descends ten rows across
    the plot, so no single row holds half of it. Rows are therefore summed after shearing the
    page by a trial slope, keeping the slope with the strongest response. Verticals are solid
    and need no shear. See recipe traps 3 and 4.
    """
    if slopes is None:
        slopes = np.arange(-14, 15) * 0.001
    h, w = ink.shape

    def two(profile, lo, hi, minsep):
        prof = profile.astype(float).copy()
        prof[:lo] = 0.0
        prof[hi:] = 0.0
        out = []
        for _ in range(2):
            i = int(np.argmax(prof))
            if prof[i] <= 0:
                raise AssertionError("no frame candidate -- check the margin")
            out.append(i)
            prof[max(0, i - minsep) : i + minsep] = 0.0
        return sorted(out)

    lef, rig = two(ink.sum(axis=0), margin, w - margin, int(sep_frac * w))
    sub = ink[:, lef + 10 : rig - 10]

    def outermost(prof, cover):
        """Extreme rows that span the plot, not the strongest rows.

        A6's solid flat curve outscores its broken top frame on row strength. Thresholding
        on coverage of the plot WIDTH works: a footer rule is too short to qualify, a data
        line qualifies but lies inside, and a frame still spans the plot once unsheared.
        """
        q = prof.astype(float).copy()
        q[:margin] = 0.0
        q[h - margin :] = 0.0
        idx = np.flatnonzero(q >= cover * sub.shape[1])
        if idx.size == 0 or idx.max() - idx.min() < sep_frac * h:
            return None
        return int(idx.min()), int(idx.max())

    best = None
    for m in slopes:
        prof = _sheared_row_profile(sub, float(m))
        got = outermost(prof, 0.60)
        if got is None:
            continue
        a, b = got
        score = float(prof[a]) + float(prof[b])
        if best is None or score > best[0]:
            best = (score, a, b, float(m))
    if best is None:
        raise AssertionError("no horizontal frame candidates")
    _, top, bot, _m = best
    return float(lef), float(rig), float(top), float(bot)


def auto_frame_spec(page, pad=45, inset=15):
    """`spec` for `fit_edges`, derived from the page's own bitmap. See `find_frame_box`."""
    lef, rig, top, bot = find_frame_box(native_bitmap(page))
    return [
        ("left", 0, int(lef - pad), int(lef + pad), int(top + inset), int(bot - inset)),
        ("right", 0, int(rig - pad), int(rig + pad), int(top + inset), int(bot - inset)),
        ("top", 1, int(top - pad + 10), int(top + pad - 10), int(lef + inset), int(rig - inset)),
        ("bottom", 1, int(bot - pad + 10), int(bot + pad - 10), int(lef + inset), int(rig - inset)),
    ]


def curve_component(ink, fc, band=4.0, margin=30.0, close=0):
    """The joining polyline and its markers, as one 8-connected blob with the frame off.

    `close` bridges a joining line the scan broke (A9's is 1 px wide), so that "largest
    component" does not pick a fragment. Closing decides membership only; the returned mask
    is the original ink.
    """
    m = plot_ink(ink, fc, band, margin)
    g = ndimage.binary_closing(m, np.ones((close, close))) if close else m
    lab, n = ndimage.label(g, structure=np.ones((3, 3)))
    sz = ndimage.sum(g, lab, range(1, n + 1))
    return m & (lab == int(np.argmax(sz)) + 1)


def segment_lines(curve_pts, seeds, clear=24.0, band=4.5, extra=None):
    """The joining polyline, segment by segment, total least squares, clear of markers.

    Returns one row per segment: `[cx, cy, dx, dy, length, n, rms]`, where `(dx, dy)` is the
    unit direction from this marker to the next -- the glyph model needs it, because the
    printed line runs through every marker.
    """
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


def along_segments(P, C):
    """Each joining segment between consecutive centres `C`: its length, and every point of
    `P` as `(along, across)` coordinates measured from the segment's start."""
    for i in range(len(C) - 1):
        a, b = C[i], C[i + 1]
        d = b - a
        L = float(np.hypot(*d))
        d = d / L
        yield L, (P - a) @ d, (P - a) @ np.array([-d[1], d[0]])


def ink_excesses(L, t, perp, ratio, end):
    """Interior ink excesses along one segment of `along_segments` -- a marker count that
    shares nothing with the marker detector.

    The ink within 12 px of the segment is binned along it in 10 px bins; a bin holding more
    than `ratio` times the segment's own baseline, and more than `end` px from either end,
    is an excess. An empty list on every segment says no marker was missed or invented.
    """
    tt = t[np.abs(perp) < 12]
    nb = max(int(L // 10), 4)
    h, e = np.histogram(tt, bins=nb, range=(0.0, L))
    base = max(np.median(h[2:-2]) if nb > 6 else np.median(h), 1)
    return [
        round((e[k] + e[k + 1]) / 2, 1)
        for k in range(nb)
        if h[k] > ratio * base and end < (e[k] + e[k + 1]) / 2 < L - end
    ]


def print_segments(lines):
    """One line per joining segment from `segment_lines`: length, samples, rms, direction."""
    for i, ln in enumerate(lines):
        print(
            f"  seg{i + 1:2d}: len {ln[4]:6.1f} px  n {int(ln[5]):5d}  rms {ln[6]:.3f} px  "
            f"dir {np.degrees(np.arctan2(ln[3], ln[2])):+8.3f} deg"
        )


def longest_run(mask) -> int:
    """Length of the longest contiguous True run in a 1-D boolean array."""
    return max((int(b - a) for a, b in runs(mask)), default=0)


def longest_run_columns(ink, lo_frac, hi_frac, strong_frac):
    """`(left, right)`: the outermost columns whose longest vertical ink run exceeds
    `strong_frac` of the longest on the page, searched between `lo_frac` and `hi_frac` of the
    width. A vertical frame spans the plot; text and data do not. The body figures' (6-10)
    first guess at their vertical frames, refined by `band_line`."""
    w = ink.shape[1]
    m0, m1 = int(w * lo_frac), int(w * hi_frac)
    cols = np.array([longest_run(ink[:, x]) if m0 <= x < m1 else 0 for x in range(w)])
    strong = np.flatnonzero(cols > strong_frac * cols.max())
    return int(strong.min()), int(strong.max())


def band_line(ink, along, nominal, lo, hi, step, halves, max_ink, min_samples):
    """One straight frame line from the ink centroid of a narrow band about it: the body
    figures' (6-10) frame fitter.

    `along=0` fits a vertical line, column as a function of row; `along=1` a horizontal one.
    Positions `range(lo, hi, step)` along the line are sampled; pass k looks `halves[k]` px
    either side of the previous pass's line, the first about the constant `nominal`. A
    position is a sample if its band holds 1 to `max_ink` inked pixels: more is a tick or a
    crossing trace. A pass with fewer than `min_samples` samples ends the fit. Returns
    `((slope, intercept), passes completed, residual std of the last completed pass)`; the
    line is that pass's, or `(0, nominal)` if none completed. Callers choose what an
    incomplete fit means.
    """
    fitted = (0.0, float(nominal))
    res = float("nan")
    for k, half in enumerate(halves):
        ts: list[float] = []
        vs: list[float] = []
        for t in range(lo, hi, step):
            a = int(round(fitted[0] * t + fitted[1] - half))
            if along == 0:
                band = ink[t, max(a, 0) : a + 2 * half + 1]
            else:
                band = ink[max(a, 0) : a + 2 * half + 1, t]
            idx = np.flatnonzero(band)
            if idx.size and idx.size <= max_ink:
                ts.append(float(t))
                vs.append(a + float(idx.mean()))
        if len(ts) < min_samples:
            return fitted, k, res
        s, b = np.polyfit(ts, vs, 1)
        fitted = (float(s), float(b))
        res = float(np.std(np.asarray(vs) - (s * np.asarray(ts) + b)))
    return fitted, len(halves), res


def edge_uv(fr, H):
    """`uv_of(name, positions)` -- where ticks on one straight edge land in the unit square."""

    def uv_of(nm, pos):
        a, b = fr[nm]
        if nm in ("left", "right"):
            return apply_h(H, [a + b * q for q in pos], pos)
        return apply_h(H, pos, [a + b * q for q in pos])

    return uv_of


def edge_at_straight(uv_of):
    """`at(name, positions)` -- the along-edge coordinate from `edge_uv`'s map: v on the left
    and right edges, u on the top and bottom."""

    def at(nm, pos):
        q = uv_of(nm, list(np.atleast_1d(pos)))
        return q[1] if nm in ("left", "right") else q[0]

    return at


def edge_at(fc, to_sv):
    """`at(name, coords)` -- the same for the arc model, returning the along-edge coordinate."""

    def at(nm, cs):
        cs = np.asarray(cs, float)
        if nm in ("left", "right"):
            return to_sv(np.polyval(fc[nm], cs), cs)[1]
        return to_sv(cs, np.polyval(fc[nm], cs))[0]

    return at


def unit_square(fc, order):
    """`(to_uv, to_px, tick_u)` under one frame model: pixels to the unit square, back, and
    `tick_u(name, positions)` for ticks along one edge, as the along-edge coordinate.

    Order 1 maps the four straight edges by the homography between their corners; order 2
    is `make_sv`'s curvilinear blend between the four arcs.
    """
    if order == 1:
        # the homography, not `make_sv`: on straight edges they differ by up to 1.1e-3 of
        # the unit square, and every order-1 figure is calibrated through the homography
        cs = corners(fc, order=1)
        H, Hinv = homography(cs, UNIT), homography(UNIT, cs)

        def to_uv(px, py):
            return apply_h(H, px, py)

        def to_px(u, v):
            return apply_h(Hinv, u, v)

        return to_uv, to_px, edge_at_straight(edge_uv(as_lines(fc), H))
    to_uv, to_px = make_sv(fc)
    return to_uv, to_px, edge_at(fc, to_uv)
