"""Square up a scanned figure page, so a family of curves can be separated on it.

**This resamples the page, which everything in `tools/figure/` refuses to do.** Its one
caller is `digitize_a1.py`: on a 600 dpi raster Figure A1's eleven crossing speed lines are
each one 8-connected ink component, and component labelling needs the frame square. A1's
axes are calibrated on the native scan, so nothing resampled here reaches a published number
through the calibration (see `digitize_a1.py`).

The frames are keystoned (pdf pp.58-61 and 66: the bottom line drops 10-16 px across the
plot, the top a different amount), so `rectify` fits the four frame lines, intersects them,
and maps that quadrilateral onto its enclosing rectangle.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

from figure.frame import corner_line, homography, longest_run
from figure.raster import load_gray
from figure.ticks import index_groups

# --------------------------------------------------------------------------- line finding


def long_lines(ink: np.ndarray, centre: int, half: int = 100, frac: float = 0.80) -> list[float]:
    """Centres of near-continuous ink lines crossing the column band around `centre`.

    A frame line fills nearly the whole band; a row of text does not, however much ink it
    carries. That distinction is what makes this robust where a row-sum is not.
    """
    a, b = int(centre - half), int(centre + half)
    band = ink[:, a:b]
    runs = np.array([longest_run(band[y]) for y in range(band.shape[0])])
    hits = np.flatnonzero(runs > frac * (b - a))
    if hits.size == 0:
        return []
    return [float(np.mean(g)) for g in index_groups(hits, gap=3)]


def robust_fit(xs, ys, tol: float = 12.0) -> tuple[float, float]:
    """Least-squares line with outlier rejection. Returns (slope, intercept)."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    keep = np.ones(xs.shape, bool)
    for _ in range(6):
        s, i = np.polyfit(xs[keep], ys[keep], 1)
        new = np.abs(ys - (i + s * xs)) < tol
        if new.sum() < 3 or (new == keep).all():
            break
        keep = new
    return tuple(np.polyfit(xs[keep], ys[keep], 1))


def edge_fits(ink: np.ndarray, axis: int) -> tuple[tuple[float, float], tuple[float, float]]:
    """The two frame lines perpendicular to `axis`, each as (slope, intercept).

    axis=0 gives the top and bottom lines as row(x); axis=1 the left and right as col(y).
    """
    im = ink if axis == 0 else ink.T
    cols = im.sum(axis=0)
    nz = np.flatnonzero(cols > 0.02 * im.shape[0])
    lo, hi = int(nz.min()), int(nz.max())
    ts, firsts, lasts = [], [], []
    for c in range(lo + 120, hi - 119, 60):
        ll = long_lines(im, c)
        if ll:
            ts.append(c)
            firsts.append(min(ll))
            lasts.append(max(ll))
    if len(ts) < 3:
        raise ValueError("could not find two frame lines -- is this a plot page?")
    fmed, lmed = np.median(firsts), np.median(lasts)
    kf = [i for i in range(len(ts)) if abs(firsts[i] - fmed) < 80]
    kl = [i for i in range(len(ts)) if abs(lasts[i] - lmed) < 80]
    return (
        robust_fit([ts[i] for i in kf], [firsts[i] for i in kf]),
        robust_fit([ts[i] for i in kl], [lasts[i] for i in kl]),
    )


# --------------------------------------------------------------------------- rectify


def rectify(page: Path, out: Path) -> tuple[int, int, int, int]:
    """Warp `page` so its plot frame is an exact rectangle. Returns (left, top, right, bottom)."""
    ink = load_gray(page) < 0.55
    (st, it), (sb, ib) = edge_fits(ink, 0)
    (sl, il), (sr, ir) = edge_fits(ink, 1)

    tl = corner_line((it, st), (il, sl))
    tr = corner_line((it, st), (ir, sr))
    bl = corner_line((ib, sb), (il, sl))
    br = corner_line((ib, sb), (ir, sr))
    left = round((tl[0] + bl[0]) / 2)
    right = round((tr[0] + br[0]) / 2)
    top = round((tl[1] + tr[1]) / 2)
    bottom = round((bl[1] + br[1]) / 2)

    im = Image.open(page).convert("L")
    # PIL's PERSPECTIVE maps an OUTPUT point to an input point: the rectangle to the scan.
    coeffs = homography(
        [(left, top), (right, top), (right, bottom), (left, bottom)], [tl, tr, br, bl]
    ).ravel()[:8]
    im.transform(
        im.size, Image.PERSPECTIVE, tuple(coeffs), resample=Image.BICUBIC, fillcolor=255
    ).save(out)
    return left, top, right, bottom
