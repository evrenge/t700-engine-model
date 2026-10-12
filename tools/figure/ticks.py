"""The printed tick lattice, and which printed value each tick carries.

Two detectors:

* `protrusion` + `pick_majors` probes inward from the frame and matches each candidate to
  its predicted lattice slot: the interior MAJOR ticks only.
* `edge_run` + `find_ticks` + `label_*` measures the ink run through the edge, with the
  frame's own width removed as a rolling median, so a frame line that thickens along the
  page cannot pass for a tick. It finds the MINOR lattice too (recipe, trap 11).

Match ticks to values by lattice slot, never by protrusion length: a marker on a frame line
out-protrudes a major tick, and A10's lowest left-edge major protrudes only as far as a
minor. And assert on structure, not on count: a count that matches is not a check (see the
recipe).
"""

from __future__ import annotations

import numpy as np

from .frame import EDGE, corner_arc, edge_at, edge_spans
from .frame import corners as frame_corners

# --------------------------------------------------------------------------- the lattice

MAX_TICK_WIDTH_PX = 10
"""Widest protrusion run accepted as a major tick. Matched majors are 4 px wide at the
median, 5 at the 99th percentile; wider runs are ticks merged into a thickened frame (C25's
right edge) or a marker printed on the frame."""

SLOT_TOL = 0.1
"""How far from its predicted slot a major may be matched, as a fraction of the major step.
Correct matches sit within 0.05; at 0.25 a minor tick could stand in for a major hidden
under a marker (A11, C26)."""

MINOR_RATIO_TOL = 0.02
"""Largest departure of major/minor spacing from a whole number, relative. Worst 1.7 %."""


def major_step(ax):
    """The axis's major-tick interval, whether given as `step` or implied by `majors`."""
    if "step" in ax:
        step = float(ax["step"])
    else:
        v = np.asarray(ax["majors"], float)
        d = np.diff(v)
        if d.size == 0:
            raise AssertionError("an axis needs either `step` or two or more `majors`")
        if np.ptp(d) > 1e-9 * abs(d[0]):
            raise AssertionError(f"majors are not evenly spaced: {v}")
        step = float(d[0])
    if step <= 0.0:
        # a negative step disables `_dedupe`'s off-lattice rejection without a sound
        raise AssertionError(f"major step is {step}; list `majors` ascending")
    return step


def major_values(ax):
    """Interior major-tick values of an axis: `majors` if given, else `lo + step*(k+1)`.

    The default assumes a lattice anchored on the low frame. Figure A4's is not (frame at
    WA2c = 3.0, majors printed at 4, 6, 8, 10), so such an axis lists `majors` explicitly.
    """
    if "majors" in ax:
        return np.asarray(ax["majors"], float)
    return np.array([ax["lo"] + ax["step"] * (k + 1) for k in range(ax["n"])])


# ------------------------------------------------------- detector 1: the inward probe


def protrusion(ink, edge, along, sign, lo, hi, maxlen=60):
    """Contiguous ink growing inward from a frame edge, per position along it.

    `edge` is a polynomial in the along-edge coordinate, highest power first as `np.polyfit`
    returns it, so arc edges work too. A caller holding an `(intercept, slope)` pair must
    reverse it (recipe, trap 12).
    """
    out = np.zeros(hi - lo, int)
    for j, t in enumerate(range(lo, hi)):
        base = np.polyval(edge, t)
        k = 0
        while k < maxlen:
            q = int(round(base + sign * (3 + k)))
            on = ink[t, q] if along == 0 else ink[q, t]
            if not on:
                break
            k += 1
        out[j] = k
    return out


def index_groups(idx, gap):
    """Split ascending indices into runs whose consecutive members are at most `gap` apart.
    `idx` must not be empty."""
    out, cur = [], [int(idx[0])]
    for v in idx[1:]:
        if v - cur[-1] <= gap:
            cur.append(int(v))
        else:
            out.append(cur)
            cur = [int(v)]
    out.append(cur)
    return out


def groups(p, lo, thresh, gap=3):
    """Tick candidates: ink-weighted centroid, peak protrusion, width."""
    idx = np.flatnonzero(p >= thresh)
    if idx.size == 0:
        return []
    res = []
    for g in index_groups(idx, gap):
        g = np.asarray(g)
        w = p[g].astype(float)
        res.append((float(lo + (g * w).sum() / w.sum()), float(p[g].max()), len(g)))
    return res


def pick_majors(cands, t0, t1, frac, name, max_missing=2, verbose=True):
    """Match interior major ticks to their printed values by predicted lattice slot.

    `frac` is each major's fractional position between the two frames (not necessarily even:
    A4). Returns `(centres, major tick length)`, one centre per entry of `frac` in printing
    order. A major hidden under a marker printed on the frame is NaN, to be dropped from the
    fit; more than `max_missing` of them raises.
    """
    n = len(frac)
    slots = [t0 + (t1 - t0) * f for f in frac]
    step = (t1 - t0) / (n + 1)
    cands = [g for g in cands if g[2] <= MAX_TICK_WIDTH_PX]
    win = 0.25 * abs(step)
    seed = []
    for s in slots:
        c = [g for g in cands if abs(g[0] - s) < win and 5 <= g[1] <= 25]
        if c:
            seed.append(min(c, key=lambda g: abs(g[0] - s))[1])
    if not seed:
        raise AssertionError(f"{name}: no major-tick candidates at all")
    lmaj = float(np.median(seed))
    out, missing = [], 0
    for s in slots:
        c = [
            g
            for g in cands
            if abs(g[0] - s) < SLOT_TOL * abs(step) and 0.55 * lmaj <= g[1] <= 1.7 * lmaj
        ]
        if not c:
            out.append(float("nan"))
            missing += 1
            continue
        out.append(min(c, key=lambda g: abs(g[0] - s))[0])
    if missing > max_missing:
        raise AssertionError(
            f"{name}: {missing} of {n} majors not found (L~{lmaj:.1f}) -- frame wrong?"
        )
    seen = [v for v in out if not np.isnan(v)]
    if len(set(seen)) != len(seen):
        raise AssertionError(f"{name}: duplicate tick matched")
    resid = np.array(seen) - np.array(
        [s for s, v in zip(slots, out, strict=True) if not np.isnan(v)]
    )
    if verbose:
        print(
            f"    {name:6s} {n - missing}/{n} majors, L~{lmaj:.1f} px, "
            f"slot residual {resid.std():.2f} px rms, max {np.abs(resid).max():.2f}"
            + (f"  ({missing} hidden by a marker)" if missing else "")
        )
    return out, lmaj


def majors_all_edges(ink, fc, cfg, verbose=True, order=None):
    """Every printed major on all four edges, plus every tick of any length (for masks).

    Works for either frame model. Returns `(majors, candidates, corners)`.
    """
    order = len(next(iter(fc.values()))) - 1 if order is None else order
    TL, TR, BR, BL = frame_corners(fc, order=order)
    span = edge_spans((TL, TR, BR, BL))
    maj, allt = {}, {}
    for nm in ("left", "right", "top", "bottom"):
        along, sign = EDGE[nm]
        a, b = span[nm]
        lo, hi = int(min(a, b)) + 3, int(max(a, b)) - 2
        pr = protrusion(ink, fc[nm], along, sign, lo, hi)
        cands = groups(pr, lo, 3)
        allt[nm] = cands
        ax = cfg["y" if nm in ("left", "right") else "x"]
        vals = major_values(ax)
        if nm in ("left", "right"):
            # Pixels run top to bottom, values high to low, so slot 0 is the HIGHEST
            # major. Reverse the values, then measure down from the high frame.
            frac = (ax["hi"] - vals[::-1]) / (ax["hi"] - ax["lo"])
        else:
            frac = (vals - ax["lo"]) / (ax["hi"] - ax["lo"])
        maj[nm], _ = pick_majors(cands, a, b, frac, nm, cfg.get("max_missing", 2), verbose)
    return maj, allt, (TL, TR, BR, BL)


# ------------------------------------------------- detector 2: the run through the edge


OUTWARD_CAP_PX = 6
"""How far outside the fitted edge the run may extend: the frame's own ink reaches 2 px out
at the median and 6 at the 99.9th percentile."""


def edge_run(ink, edge_poly, along, sign, lo, hi):
    """Length of the ink run through the frame edge, per position along it.

    A tick lengthens the run; the frame's own width is removed afterwards as a rolling
    median. The inward side (`sign`) is followed to its end; the outward side stops at
    `OUTWARD_CAP_PX`, so ink touching the frame from outside is not a tick.
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
        while on(t, m - sign) and abs(m - base) < OUTWARD_CAP_PX:
            m -= sign
        ln[j] = abs(p - m) + 1
    return ln


def _rolling_median(x, half):
    return np.array(
        [np.median(x[max(0, i - half) : min(x.size, i + half + 1)]) for i in range(x.size)]
    )


def find_ticks(ink, edge_poly, along, sign, lo, hi, thresh, half=60, gap=2, maxw=14):
    """Tick centres along a frame edge: ink-weighted centroid of each run excess."""
    ln = edge_run(ink, edge_poly, along, sign, lo, hi)
    ex = ln - _rolling_median(ln, half)
    idx = np.flatnonzero(ex >= thresh)
    if idx.size == 0:
        return np.zeros(0)
    out = []
    for g in index_groups(idx, gap):
        g = np.asarray(g)
        if (g[-1] - g[0] + 1) > maxw:
            continue
        w = ex[g]
        out.append(float(lo + (g * w).sum() / w.sum()))
    return np.array(out)


def _dedupe(t, est, step, tol):
    """Keep one tick per lattice slot -- the closest -- and drop anything off-lattice."""
    if step <= 0.0:
        raise AssertionError(f"lattice step is {step}; off-lattice rejection needs it positive")
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


def label_ticks(w, w_lo, w_hi, val_lo, val_hi, step, tol=0.30):
    """Attach a printed value to each detected tick, from the two frame-edge values.

    The frame values fix only a *discrete* choice -- which lattice slot a tick occupies --
    so they cannot leak into the continuous fit that is later tested against them.
    """
    w = np.asarray(w, float)
    est = val_lo + (w - w_lo) * (val_hi - val_lo) / (w_hi - w_lo)
    return _dedupe(w, est, step, tol)


def label_by_map(t, poly, step, tol=0.30):
    """Label a fine tick lattice through a coarse-lattice fit, not through a straight line.

    The straight line between the frame corners is several pixels out in mid-page, which on
    a minor lattice discards a contiguous block of good ticks (recipe, trap 11); widening the
    tolerance instead discards different ones.
    """
    return _dedupe(np.asarray(t, float), np.polyval(poly, t), step, tol)


def paired(maj, name, vals, descending):
    """Tick positions and their printed values, dropping any major hidden by a marker.

    `pick_majors` returns slot order, so slot k carries `vals[k]` on the x edges and
    `vals[-1-k]` on the y edges (`descending`). Pairing by INDEX rather than by sorting keeps
    a NaN attached to its value, so dropping it drops both.
    """
    pos = np.asarray(maj[name], float)
    v = np.asarray(vals, float)[::-1] if descending else np.asarray(vals, float)
    ok = ~np.isnan(pos)
    pos, v = pos[ok], v[ok]
    o = np.argsort(pos)
    return pos[o], v[o]


def minors_per_major(pos, major_px):
    """How many minor intervals the page prints per major, measured from the ticks.

    The major spacing over the median tick gap, which must come out a whole number; how far
    it had to be rounded is the first check. Never read it off the page by eye: A11's
    half-major ticks pass for majors (recipe, trap 11). Returns `(n, ratio, minor_px)`, or
    None if the edge carries too few ticks to tell.
    """
    pos = np.sort(np.asarray(pos, float))
    if pos.size < 4:
        return None
    gaps = np.diff(pos)
    gaps = gaps[gaps > 1.0]
    if gaps.size == 0:
        return None
    small = gaps[gaps < 0.75 * major_px]
    minor_px = float(np.median(small if small.size else gaps))
    ratio = major_px / minor_px
    return int(round(ratio)), float(ratio), minor_px


def fine_lattice(ink, edge_poly, along, sign, lo, hi, to_u, coarse, major_px, step, thresh=3):
    """Every printed tick on one edge, labelled with its value through the COARSE fit.

    `to_u` maps pixel positions along the edge to the normalized plot coordinate; `coarse`
    is a polynomial in that coordinate fitted to the majors alone (see `label_by_map`).
    Returns `(u, values, n_per_major, n_detected)` or None if no lattice resolves.
    """
    raw = find_ticks(ink, edge_poly, along, sign, lo, hi, thresh)
    got = minors_per_major(raw, major_px)
    if got is None:
        return None
    n, ratio, _minor_px = got
    if n < 1 or abs(ratio - n) > MINOR_RATIO_TOL * n:
        raise AssertionError(
            f"major/minor spacing is {ratio:.3f}, not a whole number to {MINOR_RATIO_TOL:.0%}"
        )
    u = np.asarray(to_u(raw), float)
    uk, vk, _keep = label_by_map(u, coarse, step / n)
    return uk, vk, n, int(raw.size)


def coarse_lattice(maj, cfg, tick_u):
    """Each edge's interior majors as `(coordinate, value)`, any hidden by a marker dropped.

    `tick_u(name, positions)` maps positions along an edge to the normalized plot coordinate
    under the caller's frame model. `fine_lattice_all` labels through this, and falls back
    to it.
    """
    xv, yv = major_values(cfg["x"]), major_values(cfg["y"])
    out = {}
    for nm, vals, desc in (
        ("bottom", xv, False),
        ("top", xv, False),
        ("left", yv, True),
        ("right", yv, True),
    ):
        pos, val = paired(maj, nm, vals, desc)
        out[nm] = (np.asarray(tick_u(nm, pos), float), val)
    return out


def fine_lattice_all(ink, fc, cfg, tick_u, span, coarse_of, verbose=True):
    """The minor lattice on all four edges, falling back to the majors where none resolves.

    `coarse_of` is `coarse_lattice`'s output, which the fine ticks are labelled through;
    `tick_u` is as there, and `span` is `frame.edge_spans`'s.
    """
    out = {}
    for nm in ("bottom", "top", "left", "right"):
        ax = cfg["y" if nm in ("left", "right") else "x"]
        cu, cv = coarse_of[nm]
        along, sign = EDGE[nm]
        a, b = span[nm]
        lo, hi = int(min(a, b)) + 3, int(max(a, b)) - 2
        major_px = abs(b - a) / ((ax["hi"] - ax["lo"]) / major_step(ax))
        coarse = np.polyfit(cu, cv, max(1, min(3, cu.size - 2)))
        got = fine_lattice(
            ink,
            fc[nm],
            along,
            sign,
            lo,
            hi,
            lambda raw, nm=nm: np.asarray(tick_u(nm, raw), float),
            coarse,
            major_px,
            major_step(ax),
        )
        if got is None:
            if verbose:
                print(f"  {nm:6s} no lattice resolves; this edge keeps its majors")
            out[nm] = (cu, cv)
            continue
        fu, fv, npm, ndet = got
        if verbose:
            resid = fv - np.polyval(coarse, fu)
            print(
                f"  {nm:6s} {npm:2d} minors per major, {ndet:4d} ticks detected, "
                f"{fu.size:4d} labelled ({fu.size / max(cu.size, 1):.1f}x the majors), "
                f"scatter about the coarse map {resid.std() / (major_step(ax) / npm):.3f} "
                f"of one minor interval"
            )
        out[nm] = (fu, fv)
    return out


def native_lattice(
    ink, fc, to_sv, edges, guess, xedges, yedges, fine=True, tol=0.15, verbose=False
):
    """Every printed tick on the four edges, labelled, frame-edge ticks removed: the native
    family's lattice (A2, A7, A9), in the form `axis.axis_values` reads.

    `edges[name]` is `(along, sign, lo, hi, long-tick threshold, all-tick threshold, coarse
    step, fine step)`, `guess` the seeds for the four arc corners keyed `TL`, `TR`, `BL`,
    `BR`, and `xedges`/`yedges` the frame values. `fine=False` keeps the long ticks alone
    (A2). Long ticks are labelled by `label_ticks`, the fine lattice by `label_by_map`
    through a cubic on the long ticks.
    """
    at = edge_at(fc, to_sv)
    TL = corner_arc(fc, "top", "left", guess["TL"])
    TR = corner_arc(fc, "top", "right", guess["TR"])
    BL = corner_arc(fc, "bottom", "left", guess["BL"])
    BR = corner_arc(fc, "bottom", "right", guess["BR"])
    (X0, X1), (Y0, Y1) = xedges, yedges
    ends = {
        "left": (BL[1], TL[1], Y0, Y1),
        "right": (BR[1], TR[1], Y0, Y1),
        "bottom": (BL[0], BR[0], X0, X1),
        "top": (TL[0], TR[0], X0, X1),
    }

    P = {"_xedges": xedges, "_yedges": yedges, "_span": {}}
    for nm, (along, sign, lo, hi, th_long, th_all, step_c, step_f) in edges.items():
        w0, w1, v0, v1 = ends[nm]
        coarse = find_ticks(ink, fc[nm], along, sign, lo, hi, th_long)
        wc, vc, keep_c = label_ticks(coarse, w0, w1, v0, v1, step_c, tol=tol)
        tc = at(nm, wc)
        if fine:
            pc = np.polyfit(tc, vc, 3)
            allt = find_ticks(ink, fc[nm], along, sign, lo, hi, th_all)
            t, val, keep_f = label_by_map(at(nm, allt), pc, step_f)
            n_raw, n_kept = allt.size, int(keep_f.sum())
        else:
            t, val = tc, vc
            n_raw, n_kept = coarse.size, int(keep_c.sum())
        step = step_f if fine else step_c
        # the frame's span in steps, by a straight line through the lattice, fixes an
        # unlabelled frame; it only picks an integer, so the bow (~0.03 steps on A7) cannot leak in
        cl = np.polyfit(t, val, 1)
        span = np.polyval(cl, 1.0) - np.polyval(cl, 0.0)
        P["_span"][nm] = (float(span / step_c), float(span / step_f) if fine else None)
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
