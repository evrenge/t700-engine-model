"""The axis map: pixels to printed values, with the model chosen by a held-out test.

A linear axis map cannot fit these scans. The page bows along the scan direction -- A7's y
major spacing runs 181 -> 177 -> 180.5 px top to bottom, a smooth 2 % swing worth +-4 px
against a straight line -- and that bow is a fittable, checkable polynomial rather than
noise. The question is which order, and the figures answer it themselves:

* **The held-out frame test.** Where a figure prints a value at a frame edge, fit the
  interior ticks *excluding* that value and predict it. Whichever order predicts the
  printed number best is the right model. That is a held-out test, not a fit statistic.
* **Leave-one-tick-out, as the guard.** The frame test alone over-fits: on A10's x the
  quartic wins the frame test (1.11 per-mille against linear's 1.94) while its
  leave-one-tick-out error is 40 % worse. So the winner is chosen only among orders whose
  LOO error is within 1.2x of the best, ties within 1.25x going to the lower order.
* **Tail cross-validation**, which needs no printed frame value at all: fit the middle of
  the lattice and predict the two ends.

`choose` is that rule, written once. It was written twice before -- `digitize_a7.
choose_order` and the inline block in `digitize_single_curve.calibrate_axis` -- with the
same constants and the same comment.

## Anchoring: the one place the two families still disagree

Once the order is chosen, the final map is fitted to the ticks *and* the printed frame
values. How those values enter is not a detail:

* `fit_anchored` -- **hard**. Fits `lo + (hi-lo)t + t(1-t)q(t)`, which is exact at both
  ends by construction. A7's reasoning: the frame value is *known*, and its pixel position
  comes from a ~1900-row fit of the frame line, so it is the best-determined datum on the
  axis -- while two copies among 76 ticks give it 2.6 % of the weight. "With the soft
  anchor the order-3 x map misses the left frame by 2.0 px, and marker 1 *sits on that
  frame*, so the miss lands straight in the CSV."
* `soft_anchor` -- **soft**. Appends the frame value twice among the ticks and fits an
  ordinary polynomial. What `digitize_single_curve` uses for its thirteen figures.

`digitize_fig678.bow_fit` arrived at the hard anchor independently, for the same reason,
which makes two of the five families that found it and one that did not. Both are here
because switching a figure from one to the other moves its published values, and that is a
measurement to be made and recorded per figure, not a refactor.

This module lives in tools/ and may use SciPy and PIL. Nothing here is imported by
`src/t700/`.
"""

from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------- held-out tests


def polyfit_predict(t, val, order, at):
    """Fit and evaluate, returning `(prediction, coefficients)`."""
    c = np.polyfit(t, val, order)
    return np.polyval(c, at), c


def loo(t, val, order):
    """Leave-one-tick-out error: the guard against an order that buys the frame with noise."""
    e = [
        val[k] - np.polyval(np.polyfit(np.delete(t, k), np.delete(val, k), order), t[k])
        for k in range(t.size)
    ]
    return float(np.std(e))


def tail_cv(t, val, order, frac=0.18):
    """Fit the middle of the tick lattice, predict the two ends. Uses no printed edge."""
    lo, hi = t.min(), t.max()
    span = hi - lo
    inner = (t > lo + frac * span) & (t < hi - frac * span)
    c = np.polyfit(t[inner], val[inner], order)
    r = val[~inner] - np.polyval(c, t[~inner])
    return float(np.sqrt((r**2).mean())), float(np.abs(r).max())


def choose(errs, loos, loo_tol=1.20, tie=1.25):
    """Minimise the held-out frame error among orders that also interpolate the ticks well.

    `errs` is the worst held-out frame error per order, `loos` the leave-one-tick-out error
    per order. Both are held out; neither is a fit statistic on the fitted data.
    """
    lbest = min(loos.values())
    ok = [o for o in sorted(errs) if loos[o] <= loo_tol * lbest]
    order = min(ok, key=lambda k: errs[k])
    for o in sorted(ok):
        if errs[o] <= tie * errs[order]:
            return o, ok
    return order, ok


def parity(t, val, scale, deg=6):
    """Legendre coefficients of the value function, with standard errors, in page pixels.

    Even degrees are an even distortion about mid-axis, odd degrees an odd one. This is the
    test one printed frame value cannot run and a lattice of 129 minor ticks can.
    """
    u = 2 * t - 1
    A = np.polynomial.legendre.legvander(u, deg)
    coef, *_ = np.linalg.lstsq(A, val, rcond=None)
    r = val - A @ coef
    s2 = (r**2).sum() / max(t.size - deg - 1, 1)
    se = np.sqrt(np.diag(np.linalg.inv(A.T @ A) * s2))
    return [(k, coef[k] * scale, se[k] * scale, abs(coef[k] / se[k])) for k in range(2, deg + 1)]


def poly_se(t_fit, resid_std, order, t_at):
    """1 sigma standard error of a polynomial fit's prediction, at `t_at`."""
    A = np.vander(np.atleast_1d(t_fit), order + 1)
    cov = np.linalg.inv(A.T @ A) * resid_std**2
    a = np.vander(np.atleast_1d(t_at), order + 1)
    return np.sqrt(np.einsum("ij,jk,ik->i", a, cov, a))


# --------------------------------------------------------------------------- the final fit


def fit_anchored(t, val, order, lo, hi):
    """Least squares through the ticks that passes EXACTLY through both frame values.

    `lo + (hi-lo)t + t(1-t)q(t)`: exact at both ends by construction, leaving order-1 free
    parameters for the bow.
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


def soft_anchor(t, val, sides, at, printed):
    """The ticks plus two copies of each anchored frame value. See the module docstring."""
    ta, va = list(t), list(val)
    for side in sides:
        ta += [at[side]] * 2
        va += [printed[side]] * 2
    return np.array(ta), np.array(va)


def blend_axis(t_lo, v_lo, t_hi, v_hi, order, edge0, edge1):
    """Two opposite edges calibrated separately, each softly anchored on its frame values.

    Returned as a pair of polynomials; the caller blends them across the plot. Fitting the
    two edges *together* averages away a real difference: on A2 the top and bottom x
    lattices disagree by 0.026 at the left frame.
    """
    p_lo = np.polyfit(
        np.concatenate([t_lo, [0.0, 1.0]]), np.concatenate([v_lo, [edge0, edge1]]), order
    )
    p_hi = np.polyfit(
        np.concatenate([t_hi, [0.0, 1.0]]), np.concatenate([v_hi, [edge0, edge1]]), order
    )
    return p_lo, p_hi


def blend_value(edgec, t, across):
    """The axis value at normalized coordinate `t`, blending the two ticked edges.

    `across` runs 0 at the low-side edge to 1 at the high-side edge, so the map used at a
    point is the one its own end of the plot measured. The pooled fit -- what
    `digitize_single_curve` did until 2026-09-16 -- is the average everywhere, so the two
    agree mid-plot and diverge toward the edges, by up to half the edge disagreement.

    That term stays in the uncertainty budget: blending fixes the central value but says
    nothing about how the lattice varies BETWEEN the edges, which linear interpolation only
    assumes.
    """
    lo = np.polyval(edgec[0], t)
    hi = np.polyval(edgec[1], t)
    return (1.0 - across) * lo + across * hi


# --------------------------------------------------------------------------- orchestration


def calibrate(t_lo_edge, t_hi_edge, spec, label, orders=(1, 2, 3, 4), verbose=True):
    """Pick the order by the held-out frame test, then fit -- the single-curve family's path.

    `t_lo_edge` / `t_hi_edge` are (coordinate, value) for the interior ticks read off each
    of the two edges that carry ticks for this axis -- the majors alone, or the full minor
    lattice where the page resolves one. Coordinate is u (x axis) or v (y axis);
    0 at the low-side frame, 1 at the high-side frame. The printed frame values are NOT in
    the fit used for the test -- that is the point.
    """
    from .ticks import major_step

    tt = np.concatenate([t_lo_edge[0], t_hi_edge[0]])
    vv = np.concatenate([t_lo_edge[1], t_hi_edge[1]])
    printed = {"lo": spec["lo"], "hi": spec["hi"]}
    at = {"lo": 0.0, "hi": 1.0}
    rng = spec["hi"] - spec["lo"]
    if verbose:
        print(f"  {label} axis -- held-out frame test (fit = the {tt.size} interior ticks only):")
    errs, loos, pred = {}, {}, {}
    for order in orders:
        row, worst = [], 0.0
        for side in spec["printed"]:
            p, _ = polyfit_predict(tt, vv, order, at[side])
            e = float(p) - printed[side]
            worst = max(worst, abs(e))
            row.append(f"{side} frame printed {printed[side]:g} -> {float(p):.6g} ({e:+.3g})")
        pred[order] = row
        errs[order] = worst
        loos[order] = loo(tt, vv, order)
        if verbose:
            print(
                f"    order {order}: "
                + " | ".join(row)
                + f"   [LOO {loos[order] / rng * 1000:.3f} per-mille of range]"
            )
    order, _ok = choose(errs, loos)
    if verbose:
        print(
            f"    -> chosen order {order} "
            f"({['', 'linear', 'quadratic', 'cubic', 'quartic'][order]}); "
            f"worst frame error {errs[order]:.3g} = {errs[order] / rng * 1000:.2f} per-mille "
            f"of range"
        )
        for nm, d in (("low-side edge", t_lo_edge), ("high-side edge", t_hi_edge)):
            row = []
            for side in spec["printed"]:
                p, _ = polyfit_predict(d[0], d[1], order, at[side])
                row.append(f"{side} -> {float(p):.6g} ({float(p) - printed[side]:+.3g})")
            print(f"    {nm} alone, order {order}: " + " | ".join(row))

    # How many major steps does the frame span, measured by the ticks alone? On the four
    # axes where BOTH frame values are printed it comes out within 0.01 of an integer every
    # time, so on the two where one frame is unlabelled (A8 x, A10 y) the integer count
    # fixes that frame's value and it is used as an anchor -- derived, not printed.
    c1 = np.polyfit(tt, vv, 1)
    span = (np.polyval(c1, 1.0) - np.polyval(c1, 0.0)) / major_step(spec)
    anchors = spec.get("anchor", spec["printed"])
    if verbose:
        print(
            f"    tick lattice: the frame spans {span:.4f} major steps "
            f"(nearest integer {round(span)}, printed span "
            f"{(spec['hi'] - spec['lo']) / major_step(spec):.0f})"
        )
        if set(anchors) - set(spec["printed"]):
            print(
                f"    anchoring also at the DERIVED "
                f"{sorted(set(anchors) - set(spec['printed']))} frame value(s), from that "
                f"integer step count"
            )

    ta, va = soft_anchor(tt, vv, anchors, at, printed)
    allc = {o: np.polyfit(ta, va, o) for o in orders}
    c = allc[order]
    r = va - np.polyval(c, ta)
    if verbose:
        print(
            f"    final fit over {len(ta)} points: resid rms {r.std():.4g}, "
            f"max {np.abs(r).max():.4g}  ({r.std() / rng * 1000:.3f} per-mille of range)"
        )
        print(f"    low frame -> {np.polyval(c, 0.0):.6g}   high frame -> {np.polyval(c, 1.0):.6g}")
    edgec = [
        np.polyfit(*soft_anchor(d[0], d[1], anchors, at, printed), order)
        for d in (t_lo_edge, t_hi_edge)
    ]
    return dict(
        c=c,
        order=order,
        allc=allc,
        edgec=edgec,
        resid=float(r.std()),
        errs=errs,
        loos=loos,
        pred=pred,
        n=int(tt.size),
    )


def calibrate_under(ink, fspec, cfg, forder, verbose=True):
    """Everything from the frame lines to the two axis maps, under one frame model.

    `forder` 1 is four straight edges mapped by a homography; 2 is four fitted arcs and the
    curvilinear blend between them. `digitize_single_curve.run` calls this twice and keeps
    the model the figure's own held-out frame test prefers, so the returned `score` -- the
    summed held-out frame error of both axes, in per-mille of their ranges -- is what that
    choice is made on. `digitize_multi_curve` and `digitize_c30` call it once, straight.

    It was written three times until 2026-10-12: here, and inline in those two tools.
    """
    from .frame import apply_h as _apply
    from .frame import as_lines, edge_at, edge_uv, fit_edges, homography, make_sv
    from .frame import corners as frame_corners
    from .ticks import coarse_lattice, fine_lattice_all, majors_all_edges

    if verbose:
        print("frame lines:")
    fc, frdiag = fit_edges(ink, fspec, order=forder, verbose=verbose)
    TL, TR, BR, BL = frame_corners(fc, order=forder)
    span = {
        "left": (TL[1], BL[1]),
        "right": (TR[1], BR[1]),
        "top": (TL[0], TR[0]),
        "bottom": (BL[0], BR[0]),
    }
    fr = as_lines(fc, span)
    if verbose:
        print(
            f"  corners TL({TL[0]:.3f},{TL[1]:.3f}) TR({TR[0]:.3f},{TR[1]:.3f}) "
            f"BR({BR[0]:.3f},{BR[1]:.3f}) BL({BL[0]:.3f},{BL[1]:.3f})"
        )
    H = homography([TL, TR, BR, BL], [(0, 1), (1, 1), (1, 0), (0, 0)])
    if forder == 1:
        to_uv = lambda px, py: _apply(H, px, py)  # noqa: E731
        to_px = lambda u, v: _apply(  # noqa: E731
            homography([(0, 1), (1, 1), (1, 0), (0, 0)], [TL, TR, BR, BL]), u, v
        )
        edge_coord = edge_uv(fr, H)

        def tick_u(nm, pos):
            q = edge_coord(nm, list(np.atleast_1d(pos)))
            return q[1] if nm in ("left", "right") else q[0]
    else:
        sv, px_of = make_sv(fc)
        to_uv, to_px = sv, px_of
        tick_u = edge_at(fc, sv)

    if verbose:
        print("ticks:")
    maj, allt, _ = majors_all_edges(ink, fc, cfg, verbose=verbose, order=forder)
    xs, ys = cfg["x"], cfg["y"]
    coarse_of = coarse_lattice(maj, cfg, tick_u)

    # ---- the MINOR lattice, which until 2026-09-16 this path printed and ignored.
    # The axis map was fitted on the interior majors alone, four to thirty-four per axis,
    # while the page prints five, ten or twenty times as many ticks. A7 calibrates on the
    # fine lattice and says what it buys: enough data "to fit and *test* an order-5 axis
    # map instead of guessing between order 2 and order 4 on twelve ticks". Nothing about
    # these thirteen pages made them different; the machinery simply never reached them.
    if verbose:
        print("fine lattice:")
    fine = fine_lattice_all(ink, fc, cfg, tick_u, span, coarse_of, verbose)

    if verbose:
        print("calibration:")
    CX = pick_lattice(
        coarse_of["bottom"], coarse_of["top"], fine["bottom"], fine["top"], xs, "x", verbose
    )
    CY = pick_lattice(
        coarse_of["left"], coarse_of["right"], fine["left"], fine["right"], ys, "y", verbose
    )
    score = sum(
        CA["errs"][CA["order"]] / (ax["hi"] - ax["lo"]) * 1000 for CA, ax in ((CX, xs), (CY, ys))
    )
    return dict(
        fc=fc,
        frdiag=frdiag,
        fr=fr,
        corners=(TL, TR, BR, BL),
        H=H,
        to_uv=to_uv,
        to_px=to_px,
        allt=allt,
        CX=CX,
        CY=CY,
        forder=forder,
        score=float(score),
    )


def axis_values(ink, spec, tickspec, C, frame_order, x_order, y_order, blend=True):
    """Data values for a set of pixel centres under one complete calibration model.

    Exists so the *model* can be varied -- frame edges straight or curved, each axis linear
    or cubic, opposite edges blended or pooled -- and the spread over defensible variants
    carried into the uncertainty as a systematic, instead of a single model's formal error
    being quoted as if the model were certain.

    This is the mechanism the whole package is built to make general: choosing between two
    geometry models, or two anchoring rules, is answered by running both and measuring the
    difference, not by argument.
    """
    from .frame import apply_h, corners, fit_edges, homography, make_sv

    fc, _ = fit_edges(ink, spec, order=frame_order, verbose=False)
    if frame_order == 1:
        # four straight edges map by the homography, as everywhere else in the package;
        # `make_sv` on straight edges differs from it by up to 1.1e-3 of the unit square
        H = homography(corners(fc, order=1), [(0, 1), (1, 1), (1, 0), (0, 0)])

        def to_sv(x, y):
            return apply_h(H, x, y)
    else:
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


def tail_cv_report(t, val, frac=0.18, orders=(1, 2, 3, 4), unit=1.0, uname=""):
    """`tail_cv` over a range of orders, printed. Uses no printed frame value."""
    for o in orders:
        rms, mx = tail_cv(t, val, o, frac)
        lo, hi = t.min(), t.max()
        span = hi - lo
        inner = (t > lo + frac * span) & (t < hi - frac * span)
        print(
            f"    order {o}: held-out rms {rms / unit:8.3f}{uname}  max {mx / unit:8.3f}{uname}"
            f"   (fit {inner.sum()}, test {(~inner).sum()})"
        )


def fit_soft_anchored(t, val, order, lo, hi):
    """Ordinary least squares with each printed frame value entered TWICE among the ticks.

    The weight is not a detail. `blend_axis` enters each frame value ONCE and this enters
    it twice, so the same ticks and the same order give different maps; `fit_anchored`
    gives the frame value infinite weight by making the fit exact there. All three are in
    use, on different figures, and which is right is the question B5 in the consolidation
    notes asks per figure.
    """
    return np.polyfit(
        np.concatenate([t, [0.0, 0.0, 1.0, 1.0]]),
        np.concatenate([val, [lo, lo, hi, hi]]),
        order,
    )


def pick_lattice(lo_maj, hi_maj, lo_fine, hi_fine, spec, label, verbose=True):
    """Majors alone or the full minor lattice? Let the figure's own held-out test say.

    The fine lattice is five to twenty-five times the data, and on twenty-one of the
    twenty-six axes here it predicts the printed frame value better -- A6's x by a factor
    of 44, C28's y by 9.4, C24's x by 7.8. On five it predicts it WORSE, C27's y by a factor
    of five, because more ticks also means more chances to mislabel one and those edges are
    the dirtiest on the page.

    So the choice is not made once for all thirteen figures, and it is not made by hand
    either. It is made by the same held-out frame test this project already uses to choose
    the polynomial order (trap 6): fit without the printed frame value, predict it, keep
    whichever model predicts it better. Ties go to the fine lattice, because more data is
    the better default when nothing distinguishes them.

    The weakness, stated: that test rests on one or two held-out numbers per axis. It is the
    strongest thing the page offers and it is not strong. What makes it bearable is that it
    is the SAME rule everywhere, applied by the tool rather than chosen per figure by
    whoever is reading the output.

    **A second weakness, which this choice is also the answer to.** Trap 6's guard on the
    polynomial order -- admit only orders whose leave-one-tick-out error is within 1.2x of
    the best -- goes inert on a fine lattice. With two hundred ticks the LOO error is set by
    the scatter of individual ticks and barely moves with order: on A6's x axis it runs
    0.899, 0.910, 0.906, 0.925 per-mille for orders 1 to 4, so every order clears the guard
    and the frame test alone decides. That would be a real hole if the frame test were a fit
    statistic, but it is held out, so an order that bought the interior with noise predicts
    the frame worse and loses. The five axes where the fine lattice does predict the frame
    worse -- C24 y, C27 y, C29 y, A11 x, A11 y (2026-10-12; the set moves with each
    figure's frame model) -- are exactly the ones this function sends back to their majors.
    """
    maj = calibrate(lo_maj, hi_maj, spec, label, verbose=False)
    fine = calibrate(lo_fine, hi_fine, spec, label, verbose=False)
    e_maj = maj["errs"][maj["order"]]
    e_fine = fine["errs"][fine["order"]]
    rng = spec["hi"] - spec["lo"]
    use_fine = e_fine <= e_maj
    if verbose:
        print(
            f"  {label} axis -- majors alone ({maj['n']} ticks, order {maj['order']}) predict "
            f"the printed frame to {e_maj / rng * 1000:.3f} per-mille;"
        )
        print(
            f"    the fine lattice ({fine['n']} ticks, order {fine['order']}) to "
            f"{e_fine / rng * 1000:.3f} per-mille  ->  "
            f"{'FINE LATTICE' if use_fine else 'MAJORS'}"
        )
    chosen = (lo_fine, hi_fine) if use_fine else (lo_maj, hi_maj)
    out = calibrate(chosen[0], chosen[1], spec, label, verbose=verbose)
    out["lattice"] = "fine" if use_fine else "majors"
    out["alt_err"] = float(e_maj if use_fine else e_fine)
    out["edges"] = chosen  # the two (coordinate, value) pairs the final fit actually used
    return out
