"""The axis map: pixels to printed values, its order chosen by a held-out test.

The page bows along the scan direction, so the map is a polynomial (recipe, trap 5). `choose`
picks the order: the lowest held-out frame error -- fit the interior ticks without a printed
frame value, then predict it -- among orders whose leave-one-tick-out error is within 1.2x of
the best (trap 6). `tail_cv` is a third check that needs no printed frame value.

The final fit adds the printed frame values, by one of three rules; switching a figure
between them moves its published values, so all three stay:

* `fit_anchored` -- hard, exact at both frames (A7 by default; Figures 6-8, each edge's order
  by `choose_anchored`);
* `soft_anchor` / `fit_soft_anchored` -- each frame value entered twice among the ticks
  (the single-curve figures, A1);
* `blend_axis` -- each entered once, per edge (A2, A9).
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
    per order. Ties within `tie` of the best go to the lower order.
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

    Even degrees are an even distortion about mid-axis, odd degrees an odd one: a test a
    fine tick lattice can run and one printed frame value cannot.
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
    parameters for the bow. The reason: the frame's pixel position comes from a fit of the
    whole frame line, so it is the best-determined datum on the axis.
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


def choose_anchored(t, val, lo, hi, orders=(1, 2, 3, 4)):
    """`fit_anchored`'s order by leave-one-tick-out alone, and the held-out error: the rule of
    the body figures 6-8, where both frames are anchors and none is left to hold out.

    Every order, 1 (the straight map) included, is scored by refitting without each tick and
    predicting it, so fitting the ticks does not spend the check that holds them out. The
    lowest error wins and a tie goes to the lower order. A tick on or outside a frame,
    `t(1-t) <= 1e-9`, is predicted but never fitted: the anchor fixes the map there already.
    Returns `(coefficients, loo, rms of the fit over all ticks, order)`.
    """
    t, val = np.asarray(t, float), np.asarray(val, float)
    inner = t * (1.0 - t) > 1e-9
    best = (None, float("inf"), 0)
    for order in orders:
        errs = []
        for i in range(t.size):
            m = np.ones(t.size, bool)
            m[i] = False
            if order > 1 and m.sum() < order:
                errs = None
                break
            c = fit_anchored(t[m & inner], val[m & inner], order, lo, hi)
            errs.append(val[i] - float(np.polyval(c, t[i])))
        if not errs:
            continue
        e = float(np.sqrt(np.mean(np.square(errs))))
        if e < best[1]:
            best = (fit_anchored(t[inner], val[inner], order, lo, hi), e, order)
    coef, e, order = best
    fitted = float(np.sqrt(np.mean(np.square(val - np.polyval(coef, t)))))
    return coef, e, fitted, order


def soft_anchor(t, val, sides, at, printed):
    """The ticks plus two copies of each anchored frame value."""
    ta, va = list(t), list(val)
    for side in sides:
        ta += [at[side]] * 2
        va += [printed[side]] * 2
    return np.array(ta), np.array(va)


def blend_axis(t_lo, v_lo, t_hi, v_hi, order, edge0, edge1):
    """Two opposite edges calibrated separately, each with both frame values entered once.

    Returned as a pair for `blend_value`. Fitting the two edges together would average away
    a real difference: on A2 the top and bottom x lattices disagree by 0.026 at the left frame.
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

    `across` runs 0 at the low-side edge to 1 at the high-side edge, so each point uses the
    map its own end of the plot measured. The edge disagreement stays in the uncertainty
    budget: the linear blend only assumes how the lattice varies between the edges.
    """
    lo = np.polyval(edgec[0], t)
    hi = np.polyval(edgec[1], t)
    return (1.0 - across) * lo + across * hi


# --------------------------------------------------------------------------- orchestration


def calibrate(t_lo_edge, t_hi_edge, spec, label, orders=(1, 2, 3, 4), verbose=True):
    """Pick the order by the held-out frame test, then fit: the single-curve family's path.

    `t_lo_edge` / `t_hi_edge` are `(coordinate, value)` for the interior ticks of the two
    ticked edges; the coordinate is u (x axis) or v (y axis), 0 at the low-side frame and 1
    at the high-side. The printed frame values are NOT in the fit used for the test.
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

    # The frame's span in major steps, by the ticks alone, is within 0.01 of an integer where
    # both frame values are printed; where one is not (A8 x, A10 y) that integer fixes it,
    # and it is used as a derived anchor.
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


FRAME_MARGIN = 0.90
"""Arcs replace straight edges only if they predict the printed frames 10 % better."""


def choose_frame_order(ink, fspec, cfg):
    """Straight edges or arcs? The figure decides (recipe trap 7): run both, keep whichever
    predicts the printed frame values better, summed over both axes. Asked on the fine
    lattice; on the majors alone A8 gives the opposite answer. Straight is the incumbent, so
    arcs must win by `FRAME_MARGIN`: a near-tie should not move a published number.

    Returns the chosen `forder` and both trials' `calibrate_under` results.
    """
    trials = {fo: calibrate_under(ink, fspec, cfg, fo, verbose=False) for fo in (1, 2)}
    forder = 2 if trials[2]["score"] < FRAME_MARGIN * trials[1]["score"] else 1
    print(
        f"frame model: straight {trials[1]['score']:.3f} per-mille vs arcs "
        f"{trials[2]['score']:.3f}  ->  {'ARCS' if forder == 2 else 'STRAIGHT'} "
        f"(summed held-out frame error over both axes; arcs must win by "
        f"{(1 - FRAME_MARGIN) * 100:.0f} %)"
    )
    return forder, trials


def calibrate_under(ink, fspec, cfg, forder, verbose=True):
    """Everything from the frame lines to the two axis maps, under one frame model.

    `forder` 1 is four straight edges and a homography; 2 is four fitted arcs and the
    curvilinear blend. `digitize_single_curve.run` calls this for both and keeps the model
    with the lower `score`, the summed held-out frame error of both axes in per-mille of
    their ranges.
    """
    from .frame import as_lines, edge_spans, fit_edges, unit_square
    from .frame import corners as frame_corners
    from .ticks import coarse_lattice, fine_lattice_all, majors_all_edges

    if verbose:
        print("frame lines:")
    fc, frdiag = fit_edges(ink, fspec, order=forder, verbose=verbose)
    TL, TR, BR, BL = frame_corners(fc, order=forder)
    span = edge_spans((TL, TR, BR, BL))
    fr = as_lines(fc, span)
    if verbose:
        print(
            f"  corners TL({TL[0]:.3f},{TL[1]:.3f}) TR({TR[0]:.3f},{TR[1]:.3f}) "
            f"BR({BR[0]:.3f},{BR[1]:.3f}) BL({BL[0]:.3f},{BL[1]:.3f})"
        )
    to_uv, to_px, tick_u = unit_square(fc, forder)

    if verbose:
        print("ticks:")
    maj, allt, _ = majors_all_edges(ink, fc, cfg, verbose=verbose, order=forder)
    xs, ys = cfg["x"], cfg["y"]
    coarse_of = coarse_lattice(maj, cfg, tick_u)

    # the minor lattice, five to twenty-five times the majors' data (recipe, trap 11)
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

    Lets the model be varied -- frame straight or curved, each axis's order, opposite edges
    blended or pooled -- so the spread over defensible variants enters the uncertainty as a
    systematic.
    """
    from .frame import fit_edges, unit_square

    fc, _ = fit_edges(ink, spec, order=frame_order, verbose=False)
    to_sv, _, _ = unit_square(fc, frame_order)
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


def held_out_table(t, val, lo, hi, unit=1.0, uname="", orders=(1, 2, 3, 4, 5, 6), label=""):
    """Fit interior ticks alone; predict BOTH printed frame values; report LOO and tail CV.

    `lo` or `hi` is None where that frame prints no value: the fit's prediction there is
    still shown, as a measurement of the unprinted frame, with no error beside it. Returns
    `{order: (interior rms, LOO, predicted low frame, predicted high frame)}`.
    """
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
        e_loo = loo(t, val, o)
        tc, tm = tail_cv(t, val, o)
        rows[o] = (float(r.std()), e_loo, a, b)
        ea = "unprinted" if lo is None else f"{(a - lo) / unit:+8.4f}{uname}"
        eb = "unprinted" if hi is None else f"{(b - hi) / unit:+8.4f}{uname}"
        print(
            f"     {o}  {r.std() / unit:10.4f}{uname} {e_loo / unit:9.4f}{uname}  "
            f"{a:11.5f} ({ea})  {b:11.5f} ({eb})"
            f"   {tc / unit:7.3f}/{tm / unit:7.3f}"
        )
    return rows


def point_sigma(P, fc, C, s, v, px, py, x_order, y_order, glyph_px, mx, my):
    """Per-point 1 sigma in x and y for the native family (A2, A7, A9), by quadrature.

    Three terms: the calibration, as each ticked edge's own standard error with its frame
    values appended, blended across the plot; the glyph-centre uncertainty `glyph_px`,
    converted through the local scale of the arcs; and the model spread `mx`/`my`. `P` is
    the tick lattice with its `_xedges`/`_yedges`, `px`/`py` the per-edge axis maps.
    Returns a dict of the totals `sx`, `sy` and the terms the tools print.
    """
    sig = {}
    for nm, poly, order, arg in (
        ("left", py[0], y_order, v),
        ("right", py[1], y_order, v),
        ("bottom", px[0], x_order, s),
        ("top", px[1], x_order, s),
    ):
        t, val = P[nm]
        lo, hi = P["_yedges"] if nm in ("left", "right") else P["_xedges"]
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
    return dict(sx=sx, sy=sy, cal_x=cal_x, cal_y=cal_y, dXdpx=dXdpx, dYdpx=dYdpx)


def fit_soft_anchored(t, val, order, lo, hi):
    """Ordinary least squares with each printed frame value entered TWICE among the ticks.

    Not the same map as `blend_axis` (once) or `fit_anchored` (exact) on the same ticks.
    """
    return np.polyfit(
        *soft_anchor(t, val, ("lo", "hi"), {"lo": 0.0, "hi": 1.0}, {"lo": lo, "hi": hi}), order
    )


def pick_lattice(lo_maj, hi_maj, lo_fine, hi_fine, spec, label, verbose=True):
    """Majors alone or the full minor lattice? The figure's own held-out frame test decides.

    The fine lattice usually predicts the printed frame value better; on the dirtiest edges
    more ticks means more chances to mislabel one, and it predicts worse. Ties go to the
    fine lattice. The test rests on one or two held-out numbers per axis, but it is the same
    rule everywhere. On a fine lattice the LOO guard goes inert (recipe, trap 6).
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
