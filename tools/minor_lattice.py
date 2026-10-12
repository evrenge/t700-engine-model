"""Measure each figure's MINOR tick lattice, and say whether it can be trusted.

## Why

The thirteen single-curve figures calibrate their axis maps on the interior MAJOR ticks
alone -- four to thirty-four per axis. A7 and A9 calibrate on the full minor lattice, 131
and 141 ticks, and A7's docstring says what that bought: "Ten times the data the long ticks
give, and enough to fit and *test* an order-5 axis map instead of guessing between order 2
and order 4 on twelve ticks." The same lattice is printed on all thirteen pages and none of
them uses it.

## The minor step is measured, not read

Reading "there are five minors per major" off a 1988 scan by eye is exactly the kind of
transcription this project forbids, and it is unnecessary: the lattice states its own step.
Sort the detected ticks along an edge, take the median gap, and divide the major spacing by
it. That ratio must come out an integer -- a printed lattice has a whole number of minor
intervals per major -- and how close it lands to one is the first check on the whole
measurement.

## The check that matters is SLOT OCCUPANCY, not count

`docs/notes/digitizing-recipe.md`: "assert on structure, not on count -- and a count that
matches is not a check." A run of the A2 y detector once returned exactly 189 minor ticks,
the exactly correct number, while holding one spurious tick and missing one real one; the
off-by-one labelling moved the axis by 0.03 and every count-based check passed. So what is
reported here per edge is: how many lattice slots exist between the frames, how many are
filled, how many are filled TWICE, and how far the worst tick sits from its slot.

An edge whose slots are cleanly occupied can carry a fine calibration. One whose are not
says so, and keeps its majors.

Writes nothing.

    python3 tools/minor_lattice.py           # every figure
    python3 tools/minor_lattice.py a8 a10    # just these
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from digitize_single_curve import FIGS  # noqa: E402
from figure.frame import (  # noqa: E402  # noqa: E402
    EDGE,
    as_lines,
    auto_frame_spec,
    corners,
    edge_at_straight,
    edge_uv,
    fit_edges,
    homography,
)
from figure.raster import native_bitmap  # noqa: E402
from figure.ticks import (  # noqa: E402  # noqa: E402
    find_ticks,
    groups,
    label_by_map,
    major_step,
    major_values,
    minors_per_major,
    paired,
    pick_majors,
    protrusion,
)


def label_fine(u, coarse_poly, step, tol=0.30):
    """Label the fine lattice through the COARSE fit, never through a straight line.

    This is A7's `label_by_map`, and the reason it exists is the reason the first version
    of this tool looked like a failure: the straight line between the two frame corners is
    several pixels out in mid-page -- that is the page bow, the thing being measured -- and
    on a 13 px minor lattice it throws away a contiguous block of perfectly good ticks. A9
    lost 29 of its 141 that way. Widening the tolerance instead just throws away different
    ones.
    """
    return label_by_map(np.asarray(u, float), coarse_poly, step, tol)


def report(key, verbose=True):
    cfg = FIGS[key]
    ink = native_bitmap(cfg["page"])
    fspec = cfg.get("frame_spec") or auto_frame_spec(cfg["page"])
    fc, _ = fit_edges(ink, fspec, order=1, verbose=False)
    fr = as_lines(fc)
    TL, TR, BR, BL = corners(fc, order=1)
    H = homography([TL, TR, BR, BL], [(0, 1), (1, 1), (1, 0), (0, 0)])
    at = edge_at_straight(edge_uv(fr, H))
    span = {
        "left": (TL[1], BL[1]),
        "right": (TR[1], BR[1]),
        "top": (TL[0], TR[0]),
        "bottom": (BL[0], BR[0]),
    }
    if verbose:
        print(f"\n{cfg['fig']}  (pdf p.{cfg['page']})")
        print(
            f"  {'edge':>7} {'majors':>7} {'n/maj':>6} {'ratio':>7} {'minor':>9} "
            f"{'fine':>6} {'kept':>6} {'gain':>6} {'rms slot':>9}"
        )
    out = {}
    for nm in ("bottom", "top", "left", "right"):
        along, sign = EDGE[nm]
        a, b = span[nm]
        lo, hi = int(min(a, b)) + 3, int(max(a, b)) - 2
        ax = cfg["y" if nm in ("left", "right") else "x"]
        vals = major_values(ax)
        desc = nm in ("left", "right")
        frac = (
            (ax["hi"] - vals[::-1]) / (ax["hi"] - ax["lo"])
            if desc
            else (vals - ax["lo"]) / (ax["hi"] - ax["lo"])
        )
        mpos, _lmaj = pick_majors(
            groups(protrusion(ink, fc[nm], along, sign, lo, hi), lo, 3),
            a,
            b,
            frac,
            nm,
            cfg.get("max_missing", 2),
            verbose=False,
        )
        mp, mv = paired({nm: mpos}, nm, vals, desc)
        umaj = at(nm, mp)
        major_px = abs(b - a) / ((ax["hi"] - ax["lo"]) / major_step(ax))

        ticks = find_ticks(ink, fc[nm], along, sign, lo, hi, 3)
        got = minors_per_major(ticks, major_px)
        if got is None or umaj.size < 3:
            if verbose:
                print(f"  {nm:>7} {len(vals):7d}   -- too few ticks to measure a lattice")
            continue
        n, ratio, _minor_px = got
        step = major_step(ax) / n
        # a coarse map from the majors alone, then the fine lattice labelled through it
        order = min(3, umaj.size - 2)
        coarse = np.polyfit(umaj, mv, max(order, 1))
        ufine = at(nm, ticks)
        uk, vk, keep = label_fine(ufine, coarse, step)
        resid = vk - np.polyval(coarse, uk)
        out[nm] = dict(u=uk, v=vk, n=n, step=step, order=order)
        if verbose:
            print(
                f"  {nm:>7} {len(vals):7d} {n:6d} {ratio:7.3f} {step:9.4g} "
                f"{ticks.size:6d} {int(keep.sum()):6d} "
                f"{int(keep.sum()) / max(umaj.size, 1):5.1f}x "
                f"{resid.std() / step:9.3f}"
            )
    return out


def main(argv):
    keys = [a for a in argv if a in FIGS] or list(FIGS)
    print("MINOR TICK LATTICE, MEASURED FROM THE PAGE'S OWN TICKS")
    print("ratio = major spacing / median tick gap; it must land on an integer.")
    print("`fine` ticks detected, `kept` labelled through the coarse fit, `gain` against")
    print("the majors that calibration uses today. `rms slot` is the labelled ticks'")
    print("scatter about the coarse map, in units of one minor interval.")
    for key in keys:
        report(key)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
