"""The marker: a generative model of the local ink, fitted sub-pixel.

A printed marker on these pages is two crossing strokes -- an `x` -- with the joining
polyline running through it. Blob detection is useless on that (trap 1 of the recipe: the
curve and its markers are one connected component), and a centroid of the local ink is
pulled off-centre by whichever arm of the line is longer. So the marker is *modelled*:
two strokes plus the two joining half-segments, all radiating from one point, supersampled
to a sub-pixel area fraction, and scored against ink **and background** over an unmasked
disc. Scoring the background is what stops a near-parallel pair of lines from outscoring a
four-armed glyph.

The template is global -- the plotter drew the same mark at every point -- so the fit
alternates: every centre against a fixed template, then the template against all centres.

## What was duplicated, and what it cost

This existed twice, and the two copies had drifted apart on exactly the two points the
recipe records as traps:

* **Trap 8, the Nelder-Mead simplex.** SciPy builds its initial simplex by perturbing each
  coordinate by 5 % of its own value. At an image column of 2269 that is a 113 px first
  step, outside any scoring window, where the cost is flat and every printed diagnostic
  looks converged -- on A6's right-hand marker the fit took that step and settled 125 px
  off the data line. The cure is to optimise the **offset** from the current centre with
  an explicit 2 px simplex. `digitize_single_curve` always did; `digitize_native` made it
  optional and defaulted it off, and of its three callers only `digitize_a7` passed it.
  A2 and A9 did not -- and those produce `f2` and `f9`, two of the highest-leverage maps
  per +1 % of their value.
* **Trap 9, the seed guard.** `digitize_single_curve` raises if a fitted centre moves more
  than 10 px from its seed, which is what makes a bad seed loud instead of plausible.
  `digitize_native` had no such check.

`simplex` therefore defaults to 2.0 and `max_move` to 10.0 here -- the guarded behaviour
is the default, and a caller that wants the unguarded one has to say so at the call site
where it can be seen and dated.

## The parameters that still differ between callers

Not unified, deliberately: unifying them would move published numbers, and that is a
separate decision from removing a duplicate implementation. Each is a caller's current
value, carried over unchanged.

    seglen     45.0 in the single-curve family, 42.0 in the native family, 32.0 in A9
    centre     Nelder-Mead tolerances: xatol 2e-3 / fatol 1e-6 / maxiter 400
               against 1e-4 / 1e-8 / 4000
    template   Powell tolerances: xtol 1e-3 / ftol 1e-5 / maxfev 1200
               against 1e-5 / 1e-8 / maxiter 20000 / maxfev 60000

This module lives in tools/ and may use SciPy and PIL. Nothing here is imported by
`src/t700/`.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

SS = 5
"""Supersampling per axis, for a sub-pixel coverage model. Both families used 5."""

TPL0 = (np.radians(-44.05), 14.88, np.radians(46.09), 15.22, 3.46, 3.52)
"""Template seed: (angle1, half-length1, angle2, half-length2, glyph width, line width).

A7's fitted template, which is the one figure whose markers stand clearest of everything
else on the page. A seed, not a constraint -- `fit_glyphs` refits it.
"""

CENTRE_CURVE = dict(xatol=2e-3, fatol=1e-6, maxiter=400)
CENTRE_NATIVE = dict(xatol=1e-4, fatol=1e-8, maxiter=4000)
TEMPLATE_CURVE = dict(xtol=1e-3, ftol=1e-5, maxfev=1200)
TEMPLATE_NATIVE = dict(xtol=1e-5, ftol=1e-8, maxiter=20000, maxfev=60000)


def windows(ink, C0, rad, mask_fn=None, pad=5):
    """One scoring window per marker: the used pixels flattened, with an SSxSS sub-grid.

    Flattened rather than kept as a dense (H, W, SS, SS) array -- the dense form was
    `digitize_native`'s and it evaluates the model at every pixel of the bounding box
    including the ones the disc excludes, which is most of them.

    `mask_fn(X, Y) -> bool array` marks pixels to EXCLUDE (frame lines, tick marks). See
    `blocked_mask`.
    """
    W = []
    off = (np.arange(SS) + 0.5) / SS - 0.5
    ox, oy = np.meshgrid(off, off, indexing="ij")
    ox, oy = ox.ravel(), oy.ravel()
    half = int(rad) + pad
    for x0, y0 in C0:
        xs = np.arange(int(x0) - half, int(x0) + half + 1)
        ys = np.arange(int(y0) - half, int(y0) + half + 1)
        X, Y = np.meshgrid(xs, ys)
        use = ((X - x0) ** 2 + (Y - y0) ** 2) < rad**2
        if mask_fn is not None:
            use &= ~mask_fn(X, Y)
        data = ink[ys[0] : ys[-1] + 1, xs[0] : xs[-1] + 1].astype(float)
        px, py = X[use].astype(float), Y[use].astype(float)
        W.append(
            dict(
                X=X,
                Y=Y,
                use=use,
                data=data,
                px=px,
                py=py,
                d=data[use],
                OX=px[:, None] + ox[None, :],
                OY=py[:, None] + oy[None, :],
                keep=np.ones(px.shape, bool),
                seg=None,
            )
        )
    return W


def _bar(w, c, ang, half, wid, halfline=False):
    """Coverage of one stroke: a bar of width `wid`, at angle `ang`, through `c`."""
    u0, u1 = np.cos(ang), np.sin(ang)
    dx = w["OX"] - c[0]
    dy = w["OY"] - c[1]
    t = dx * u0 + dy * u1
    p = -dx * u1 + dy * u0
    m = np.abs(p) <= wid / 2
    m &= (t >= 0) & (t <= half) if halfline else (np.abs(t) <= half)
    return m


def dirs(C, i):
    """Unit vectors from marker `i` towards its neighbours -- the joining segments."""
    out = []
    for j in (i - 1, i + 1):
        if 0 <= j < len(C):
            d = np.asarray(C[j], float) - np.asarray(C[i], float)
            out.append(d / np.hypot(*d))
        else:
            out.append(None)
    return out


def coverage(w, c, tpl, pin, pout, seglen, withline=True):
    """`stroke1 u stroke2 u incoming half-segment u outgoing half-segment`, all radiating
    from one point, as a sub-pixel area fraction per used pixel."""
    th1, l1, th2, l2, wg, wl = tpl
    m = _bar(w, c, th1, l1, wg) | _bar(w, c, th2, l2, wg)
    if withline:
        for p in (pin, pout):
            if p is not None:
                m |= _bar(w, c, np.arctan2(p[1], p[0]), seglen, wl, True)
    return m.mean(axis=1)


def cost(w, c, tpl, pin, pout, seglen, withline=True):
    cov = coverage(w, c, tpl, pin, pout, seglen, withline)
    return float((((cov - w["d"]) ** 2) * w["keep"]).sum())


def fit_glyphs(
    W,
    C0,
    tpl0=None,
    *,
    rounds=5,
    seglen=45.0,
    withline=True,
    simplex=2.0,
    max_move=10.0,
    fit_template=True,
    centre=None,
    template=None,
    verbose=True,
):
    """Coordinate descent: each centre against its own window, then the shared template.

    `simplex` is the side of the initial Nelder-Mead simplex in PIXELS, and the centre is
    optimised as an OFFSET from its current value.

    `max_move` raises if a centre moves further than that from its seed. `None` disables
    the check, which is trap 9.

    A window may carry its own fixed `seg` pair; otherwise the joining-segment directions
    are recomputed from the current centres each round.
    """
    centre = CENTRE_CURVE if centre is None else centre
    template = TEMPLATE_CURVE if template is None else template
    C = np.asarray(C0, float).copy()
    seed = C.copy()
    tpl = np.array(TPL0 if tpl0 is None else tpl0, float)
    cur = np.nan

    def seg_for(i):
        if W[i]["seg"] is not None:
            return W[i]["seg"]
        if not withline:
            return None, None
        return dirs(C, i)

    for it in range(rounds):
        for i in range(len(C)):
            pin, pout = seg_for(i)
            if not withline:
                pin = pout = None
            c0 = C[i].copy()
            r = minimize(
                lambda dc, i=i, pin=pin, pout=pout, c0=c0, tpl=tpl: cost(
                    W[i], c0 + dc, tpl, pin, pout, seglen, withline
                ),
                np.zeros(2),
                method="Nelder-Mead",
                options=dict(
                    initial_simplex=np.array([[0.0, 0.0], [simplex, 0.0], [0.0, simplex]], float),
                    **centre,
                ),
            )
            C[i] = c0 + r.x
            # from the seed, not per round: four 10 px rounds used to pass a 40 px walk
            walked = float(np.hypot(*(C[i] - seed[i])))
            if max_move is not None and walked > max_move:
                raise AssertionError(f"marker {i} moved {walked:.1f} px from its seed -- bad seed?")
            cur = r.fun
        if fit_template:
            D = [seg_for(i) for i in range(len(C))]
            if not withline:
                D = [(None, None)] * len(C)

            def tcost(p, D=D):
                return sum(
                    cost(W[i], C[i], p, D[i][0], D[i][1], seglen, withline) for i in range(len(C))
                )

            r = minimize(tcost, tpl, method="Powell", options=dict(template))
            tpl, cur = r.x, r.fun
        if verbose:
            print(f"    round {it + 1}: cost {cur:.3f}")
    if verbose:
        print(
            f"    template: {np.degrees(tpl[0]):+.2f} deg / {tpl[1]:.2f} px and "
            f"{np.degrees(tpl[2]):+.2f} deg / {tpl[3]:.2f} px, "
            f"glyph width {tpl[4]:.2f}, line width {tpl[5]:.2f}"
        )
    return C, tpl


def glyph_only_windows(ink, C, tpl, rad, mask_fn=None):
    """The same windows, with the joining polyline's own pixels excluded from scoring.

    One of the honest variations the uncertainty budget is built from: if the fitted centre
    barely moves when the line it sits on is hidden rather than modelled, the model of the
    line is not carrying the answer.
    """
    W = windows(ink, C, rad, mask_fn=mask_fn)
    for i, w in enumerate(W):
        pin, pout = dirs(C, i)
        drop = np.zeros(w["px"].shape, bool)
        for p in (pin, pout):
            if p is None:
                continue
            wl = tpl[5] + 3.0
            dx = w["px"] - C[i][0]
            dy = w["py"] - C[i][1]
            t = dx * p[0] + dy * p[1]
            q = -dx * p[1] + dy * p[0]
            drop |= (np.abs(q) <= wl / 2) & (t >= 0)
        w["keep"] = ~drop
    return W


def hide_segments(W, C, lines, half=6.0, reach=6.0):
    """Drop the joining polyline from the score instead of modelling it.

    The native family's independent check, written identically in `digitize_a2`,
    `digitize_a7` and `digitize_a9`. `lines[j, 2:4]` is the unit direction of segment `j`;
    segment `i-1` arrives at marker `i` and segment `i` leaves it, so the arriving one is
    walked backwards.
    """
    n = len(C)
    for i, w in enumerate(W):
        for j, sign in ((i - 1, -1.0), (i, +1.0)):
            if 0 <= j < n - 1:
                d = np.asarray(lines[j, 2:4], float) * sign
                dx = w["px"] - C[i][0]
                dy = w["py"] - C[i][1]
                t = dx * d[0] + dy * d[1]
                p = -dx * d[1] + dy * d[0]
                w["keep"] &= ~((t > reach) & (np.abs(p) < half))
    return W
