"""The printed tick lattice, and which printed value each tick carries.

The ticks are the only calibration the page offers, so how many of them a tool finds sets
how good its axis map can be. Two detectors live here, and the difference between them is
not a style choice -- it is a factor of ten in the data the axis map is fitted on:

* **`protrusion` + `pick_majors`** finds the interior MAJOR ticks by probing inward from
  the frame at a fixed offset and matching each candidate to its predicted lattice slot.
  Four to seventeen ticks per edge. What `digitize_single_curve` used for thirteen figures.
* **`edge_run` + `find_ticks` + `label_*`** measures the length of the ink run *through*
  the edge and removes the frame's own width as a rolling median, so a frame line that
  thickens or thins along the page cannot masquerade as a tick -- which is exactly what a
  fixed inward probe lets happen. It finds the MINOR lattice too: 131 y-minors on A7's
  left edge, every slot filled once, and 141 on A9. What `digitize_a7` established.

A7's docstring puts the consequence plainly: the minor lattice is "ten times the data the
long ticks give, and enough to fit and *test* an order-5 axis map instead of guessing
between order 2 and order 4 on twelve ticks."

## Two rules, both learned the hard way

**Match ticks to values by lattice SLOT, never by a protrusion-length threshold.** A
marker sitting on a frame line out-protrudes a major tick (A6, A8 and A10 each print one),
A6's flat data line protrudes the full scan length, and A10's lowest left-edge major
protrudes only as far as a minor.

**Assert on structure, not on count -- and a count that matches is not a check.** A run of
the A2 y detector returned exactly 189 minor ticks, the exactly correct number, while
holding one spurious tick and missing one real one. The off-by-one labelling moved the
axis by 0.03 in T3/T2, and every count-based check passed.

This module lives in tools/ and may use SciPy and PIL. Nothing here is imported by
`src/t700/`.
"""

from __future__ import annotations

import numpy as np

from .frame import EDGE, as_lines, corner_line

# --------------------------------------------------------------------------- the lattice


def major_step(ax):
    """The axis's major-tick interval, whether given as `step` or implied by `majors`."""
    if "step" in ax:
        return float(ax["step"])
    v = np.asarray(ax["majors"], float)
    d = np.diff(v)
    if d.size == 0:
        raise AssertionError("an axis needs either `step` or two or more `majors`")
    if np.ptp(d) > 1e-9 * abs(d[0]):
        raise AssertionError(f"majors are not evenly spaced: {v}")
    return float(d[0])


def major_values(ax):
    """Interior major-tick values of an axis.

    Defaults to `lo + step*(k+1)`, which is right whenever the lattice is anchored on the
    low frame. Figure A4 is not: its frame is at WA2c = 3.0 while its majors are printed at
    4, 6, 8, 10, so they are neither at `lo + k*step` nor evenly spaced between the frames.
    Such an axis gives `majors` explicitly.
    """
    if "majors" in ax:
        return np.asarray(ax["majors"], float)
    return np.array([ax["lo"] + ax["step"] * (k + 1) for k in range(ax["n"])])


def n_majors(ax):
    """How many interior majors the axis has."""
    return len(major_values(ax))


# ------------------------------------------------------- detector 1: the inward probe


def protrusion(ink, edge, along, sign, lo, hi, maxlen=60):
    """Contiguous ink growing inward from a straight frame edge, per position along it."""
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


def groups(p, lo, thresh, gap=3):
    """Tick candidates: ink-weighted centroid, peak protrusion, width."""
    idx = np.flatnonzero(p >= thresh)
    if idx.size == 0:
        return []
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
        w = p[g].astype(float)
        res.append((float(lo + (g * w).sum() / w.sum()), float(p[g].max()), len(g)))
    return res


def pick_majors(cands, t0, t1, frac, name, max_missing=2, verbose=True):
    """Match interior major ticks to their printed values by predicted lattice slot.

    Returns one centre per entry of `frac`, in printing order, where `frac` is each major's
    fractional position between the two frames -- not necessarily evenly spaced, because
    Figure A4's majors are not.

    A marker printed ON a frame hides the tick underneath it -- A11's first marker sits at
    y = 4.00 on the left frame, and A7 records the same thing. A hidden major is recorded
    as NaN and dropped from the fit rather than failing the extraction; the axis map is
    fitted on the majors that ARE visible, and the held-out frame test still judges it.
    Losing more than `max_missing` of them is a real failure and still raises.
    """
    n = len(frac)
    slots = [t0 + (t1 - t0) * f for f in frac]
    step = (t1 - t0) / (n + 1)
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
        c = [g for g in cands if abs(g[0] - s) < win and 0.55 * lmaj <= g[1] <= 1.7 * lmaj]
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


def majors_all_edges(ink, fc, cfg, verbose=True):
    """Every printed major on all four edges, plus every tick of any length (for masks).

    Straight-edge model only: `protrusion` probes from a line. Returns
    `(majors, candidates, corners)`.
    """
    fr = as_lines(fc)
    TL = corner_line(fr["top"], fr["left"])
    TR = corner_line(fr["top"], fr["right"])
    BL = corner_line(fr["bottom"], fr["left"])
    BR = corner_line(fr["bottom"], fr["right"])
    span = {
        "left": (TL[1], BL[1]),
        "right": (TR[1], BR[1]),
        "top": (TL[0], TR[0]),
        "bottom": (BL[0], BR[0]),
    }
    maj, allt = {}, {}
    for nm in ("left", "right", "top", "bottom"):
        along, sign = EDGE[nm]
        a, b = span[nm]
        lo, hi = int(min(a, b)) + 3, int(max(a, b)) - 2
        pr = protrusion(ink, fr[nm], along, sign, lo, hi)
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


def edge_run(ink, edge_poly, along, sign, lo, hi):
    """Length of the ink run through the frame edge, per position along it.

    A tick lengthens the run. The frame's own width need not be known: it is removed
    afterwards as a rolling median, so a line that thickens or thins along the page cannot
    masquerade as a tick.
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


def find_ticks(ink, edge_poly, along, sign, lo, hi, thresh, half=60, gap=2, maxw=14):
    """Tick centres along a frame edge: ink-weighted centroid of each run excess."""
    ln = edge_run(ink, edge_poly, along, sign, lo, hi)
    ex = ln - _rolling_median(ln, half)
    idx = np.flatnonzero(ex >= thresh)
    if idx.size == 0:
        return np.zeros(0)
    grp, cur = [], [int(idx[0])]
    for v in idx[1:]:
        if v - cur[-1] <= gap:
            cur.append(int(v))
        else:
            grp.append(cur)
            cur = [int(v)]
    grp.append(cur)
    out = []
    for g in grp:
        g = np.asarray(g)
        if (g[-1] - g[0] + 1) > maxw:
            continue
        w = ex[g]
        out.append(float(lo + (g * w).sum() / w.sum()))
    return np.array(out)


def _dedupe(t, est, step, tol):
    """Keep one tick per lattice slot -- the closest -- and drop anything off-lattice."""
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

    The straight line between the two frame corners is off by several pixels in mid-page --
    that is the whole point of this exercise -- which on a 16.7 px minor lattice throws away
    a contiguous block of perfectly good ticks (29 of A9's 141), and would throw away the
    wrong ones if the tolerance were widened instead.
    """
    return _dedupe(np.asarray(t, float), np.polyval(poly, t), step, tol)


def paired(maj, name, vals, descending):
    """Tick positions and their printed values, dropping any major hidden by a marker.

    `pick_majors` returns slot order -- low edge to high edge -- so slot k carries `vals[k]`
    on the x edges and `vals[-1-k]` on the y edges, where pixel order runs opposite to value
    order. Pairing by INDEX rather than by sorting keeps a NaN attached to the value it
    belongs to, so dropping it drops both.
    """
    pos = np.asarray(maj[name], float)
    v = np.asarray(vals, float)[::-1] if descending else np.asarray(vals, float)
    ok = ~np.isnan(pos)
    pos, v = pos[ok], v[ok]
    o = np.argsort(pos)
    return pos[o], v[o]
