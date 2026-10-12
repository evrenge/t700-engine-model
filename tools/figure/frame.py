"""The four printed frame edges, and the map from pixels to plot coordinates.

Everything downstream rests on this: a frame located wrongly moves every value on the
page, and it does it without disturbing any of the checks that look only at the data. A11
once acquired a bottom frame 192 px too low, which put every digitized value 0.3 too high
while the tick lattice still appeared to close.

## The two geometry models, and why both are still here

A printed straight edge does not image straight on this scan. Fitted as a quadratic, the
frame edges of these pages carry sagittae of 0.3 to 1.8 px -- p.62's are -0.44, -0.75,
-0.38 and -0.14 px -- which is the same size as the errors the whole exercise is chasing.
The paper's shape, not perspective: on a flatbed there is nothing to keystone.

So there are two models, and until 2026-09-16 they lived in two different files with two
different sets of figures behind them:

* **order 1** -- straight edges, plot coordinates by `homography` between the four
  corners. What `digitize_single_curve` used for thirteen figures, and `digitize_a1` for
  Figure A1.
* **order 2** -- each edge an arc, plot coordinates the curvilinear blend between the four
  arcs (`make_sv`). What `digitize_a7` established and A2 and A9 followed. Recipe trap 7
  once called it the right one; measured, it is right for two of the thirteen single-curve
  figures.

Both stay reachable because choosing between them per figure is a measurement, not a
refactor: it moves published numbers, so it is made deliberately, with the movement
recorded, rather than as a side effect of deleting a duplicate.

**Nothing here resamples the page.** Rotation, tilt and bow are carried as a map applied
to the measured points. Warping pixels needs interpolation, and interpolation invents ink.

This module lives in tools/ and may use SciPy and PIL. Nothing here is imported by
`src/t700/`.
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
    """(start, end) of every contiguous True run in a 1-D boolean array."""
    e = np.flatnonzero(np.diff(np.concatenate(([0], mask.view(np.int8), [0]))))
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

    `spec` is a list of `(name, along, xlo, xhi, tlo, thi)` windows -- see
    `auto_frame_spec`. Returns `(coeffs, diag)`.

    The sagitta reported for order 2 is the whole reason that order exists: it is a
    *printed straight line* imaging as an arc, so it measures the page distortion
    directly, with no reliance on the tick lattice.
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
    """Each edge as an `(intercept, slope)` pair -- what the masks and the seed search want.

    Those consumers do not care about a 2 px bow; they need a line lying along the ink. For
    an order-1 fit the pair is exact. For a higher order it must be the best straight line
    **over the edge's own span**, which is what `span` supplies as `(t0, t1)` per edge.

    Without the span the only available answer is the polynomial's tangent at t = 0, and on
    these pages t is a pixel row running to 2600, so that tangent leaves the ink entirely.
    It put A10's frame far enough out that the marker search found no ink at all and raised
    "opening found nothing". Pass the span unless you know the fit is order 1.
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


def corners(fc, order=1, guess=None):
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
    # Seed the fixed-point iteration from the straight-line intersection of the same
    # edges' linear parts. The arcs depart from those lines by a pixel or two, so the seed
    # is already inside the basin and the iteration is only polishing.
    lin = {nm: (float(c[-1]), float(c[-2])) for nm, c in fc.items()}
    g = guess or {}
    return tuple(
        corner_arc(fc, h, v, g.get((h, v)) or corner_line(lin[h], lin[v]))
        for h, v in (("top", "left"), ("top", "right"), ("bottom", "right"), ("bottom", "left"))
    )


# --------------------------------------------------------------------------- plot coords


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
    """Curvilinear plot coordinates from the four fitted frame arcs.

    `s` is 0 on the left edge and 1 on the right, `v` is 0 on the bottom and 1 on the top
    -- exactly, by construction -- with a linear blend between. This is the curved-edge
    generalisation of the homography, and where it is used it is the only piece of
    geometry ever applied to a datum.
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

    Symmetry matters, and it is trap 10: the ticks all protrude inward, so masking only the
    inward side strips ink from one half of a marker sitting on a frame line and pulls its
    fitted centre outward.
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

    The horizontal frames on these scans are **tilted and broken**, and that combination
    defeats every row-wise measure. A11's bottom frame descends ten rows across the plot
    -- 0.36 degrees -- so no single row holds more than about half of it, and simple row
    sums pick the x-axis label band or a caption rule instead.

    So the rows are summed **after shearing** the page by a trial slope, and the slope that
    maximises the response is the one the line actually has. Verticals need no such
    treatment: they are solid, and their own tilt is small over the page width.
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

        Two failures shaped this. Taking the two *strongest* rows picks Figure A6's flat
        combustor-efficiency curve instead of its top frame, because that curve is solid
        while the frame is broken. Taking the extreme rows above a fraction of the *max*
        then drops A6's top frame entirely, for the same reason.

        Thresholding on coverage of the plot WIDTH fixes both: a footer rule is too short
        to qualify, a data line qualifies but lies inside, and a frame -- however broken --
        still spans the plot once the shear has been undone.
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

    `close` bridges a joining line the scan broke -- A9's is 1 px wide over its flat run
    and comes apart in a dozen places, which splits the curve into fragments and makes
    "largest component" pick a fragment. Closing decides *membership* only; the returned
    mask is the original ink.
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
    """The joining polyline, segment by segment, total least squares, clear of markers.

    Returns one row per segment: `[cx, cy, dx, dy, length, n, rms]`, where `(dx, dy)` is the
    unit direction from this marker to the next. The glyph model needs those directions --
    the printed line runs through every marker, so a fit that does not model it is fitting
    the line as well as the glyph.
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


def longest_run(mask) -> int:
    """Length of the longest contiguous True run in a 1-D boolean array."""
    e = np.flatnonzero(np.diff(np.concatenate(([0], np.asarray(mask).view(np.int8), [0]))))
    return 0 if e.size == 0 else int((e[1::2] - e[0::2]).max())


def edge_uv(fr, H):
    """`uv_of(name, positions)` -- where ticks on one straight edge land in the unit square."""

    def uv_of(nm, pos):
        a, b = fr[nm]
        if nm in ("left", "right"):
            return apply_h(H, [a + b * q for q in pos], pos)
        return apply_h(H, pos, [a + b * q for q in pos])

    return uv_of


def edge_at(fc, to_sv):
    """`at(name, coords)` -- the same for the arc model, returning the along-edge coordinate."""

    def at(nm, cs):
        cs = np.asarray(cs, float)
        if nm in ("left", "right"):
            return to_sv(np.polyval(fc[nm], cs), cs)[1]
        return to_sv(cs, np.polyval(fc[nm], cs))[0]

    return at
