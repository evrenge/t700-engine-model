"""Square up a scanned figure page, and measure its printed axis ticks.

Two things `digitize.py` assumes but this 1988 scan does not provide.

1. **The plot frame is not axis-parallel.** Pages 58-61 and 66 are keystoned: the bottom
   frame line drops 10-16 px across the plot width while the top drops a different amount.
   No single row/column pair describes the frame, so any two-point calibration is wrong
   somewhere on the page. `rectify` fits the four frame lines, intersects them for the
   four corners, and maps that quadrilateral onto its enclosing rectangle. Residual drift
   afterwards is under 1 px.

2. **`digitize.py axes` mis-locates the frame when the figure has a long flat data line.**
   `find_frame` averages every strong row in each half of the page, and a flat run of data
   is a strong row. On Figure A3 it returned top=362 against a measured 313 -- 0.009 in
   B1. `ticks` measures the printed major ticks instead, which gives a many-point
   least-squares calibration and a check on whatever `axes` reports.

Neither step invents a number: both are pixel measurements of printed ink.

    python3 tools/rectify_page.py rectify --page a3-058.png --out a3r-058.png
    python3 tools/rectify_page.py ticks --page a3r-058.png --frame 468,313,2133,2672

This module lives in tools/ and may use SciPy, Matplotlib and PIL. Nothing here is
imported by src/t700/.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image

from digitize import load_gray

# --------------------------------------------------------------------------- line finding


def longest_run(mask: np.ndarray) -> int:
    """Length of the longest run of True in a 1-D bool array."""
    edges = np.flatnonzero(np.diff(np.concatenate(([0], mask.view(np.int8), [0]))))
    return 0 if edges.size == 0 else int((edges[1::2] - edges[0::2]).max())


def _group(idx: np.ndarray, gap: int) -> list[list[int]]:
    groups: list[list[int]] = []
    cur = [int(idx[0])]
    for v in idx[1:]:
        if v - cur[-1] <= gap:
            cur.append(int(v))
        else:
            groups.append(cur)
            cur = [int(v)]
    groups.append(cur)
    return groups


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
    return [float(np.mean(g)) for g in _group(hits, gap=3)]


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


def perspective_coeffs(dst, src) -> np.ndarray:
    """Coefficients for PIL PERSPECTIVE, which maps an output point to an input point."""
    a, b = [], []
    for (bx, by), (sx, sy) in zip(dst, src, strict=True):
        a.append([bx, by, 1, 0, 0, 0, -sx * bx, -sx * by])
        b.append(sx)
        a.append([0, 0, 0, bx, by, 1, -sy * bx, -sy * by])
        b.append(sy)
    return np.linalg.solve(np.asarray(a, float), np.asarray(b, float))


def rectify(page: Path, out: Path) -> tuple[int, int, int, int]:
    """Warp `page` so its plot frame is an exact rectangle. Returns (left, top, right, bottom)."""
    ink = load_gray(page) < 0.55
    (st, it), (sb, ib) = edge_fits(ink, 0)
    (sl, il), (sr, ir) = edge_fits(ink, 1)

    def corner(srow, irow, scol, icol):
        x = (icol + scol * irow) / (1 - scol * srow)
        return x, irow + srow * x

    tl = corner(st, it, sl, il)
    tr = corner(st, it, sr, ir)
    bl = corner(sb, ib, sl, il)
    br = corner(sb, ib, sr, ir)
    left = round((tl[0] + bl[0]) / 2)
    right = round((tr[0] + br[0]) / 2)
    top = round((tl[1] + tr[1]) / 2)
    bottom = round((bl[1] + br[1]) / 2)

    im = Image.open(page).convert("L")
    coeffs = perspective_coeffs(
        [(left, top), (right, top), (right, bottom), (left, bottom)], [tl, tr, br, bl]
    )
    im.transform(
        im.size, Image.PERSPECTIVE, tuple(coeffs), resample=Image.BICUBIC, fillcolor=255
    ).save(out)
    return left, top, right, bottom


# --------------------------------------------------------------------------- ticks


def protrusion(
    ink: np.ndarray, edge: str, pos: int, lo: int, hi: int, maxlen: int = 70
) -> np.ndarray:
    """Contiguous ink length growing inward from a frame edge, per pixel along that edge.

    A major tick reaches further in than a minor one, which is how the two are told apart
    without knowing either length in advance.
    """
    n = ink.shape[0] if edge in ("left", "right") else ink.shape[1]
    out = np.zeros(n, int)
    for t in range(lo, hi):
        k = 0
        while k < maxlen:
            if edge == "left":
                on = ink[t, pos + 4 + k]
            elif edge == "right":
                on = ink[t, pos - 4 - k]
            elif edge == "bottom":
                on = ink[pos - 4 - k, t]
            else:
                on = ink[pos + 4 + k, t]
            if not on:
                break
            k += 1
        out[t] = k
    return out


def ticks(page: Path, frame: tuple[int, int, int, int], minlen: int = 6) -> dict[str, list]:
    """Tick centres just inside each frame edge, as (centre_px, protrusion_px).

    A centre whose protrusion is the full `maxlen` is a data line crossing the frame, not
    a tick -- those are reported too, so the caller can see and discard them.
    """
    ink = load_gray(page) < 0.55
    left, top, right, bottom = frame
    res = {}
    for edge, pos, lo, hi in (
        ("left", left, top + 3, bottom - 2),
        ("right", right, top + 3, bottom - 2),
        ("bottom", bottom, left + 3, right - 2),
        ("top", top, left + 3, right - 2),
    ):
        p = protrusion(ink, edge, pos, lo, hi)
        idx = np.flatnonzero(p >= minlen)
        res[edge] = (
            [(float(np.mean(g)), int(p[np.asarray(g)].max())) for g in _group(idx, 3)]
            if idx.size
            else []
        )
    return res


def calib_fit(px, values) -> tuple[float, float, np.ndarray]:
    """Least-squares value(px). Returns (intercept, slope, residuals in pixels).

    Feed it every printed major tick on the axis, with the frame lines included as the
    extreme ticks. Two points is what `--calib` takes, but two points is not what you
    should measure -- the residuals here are what tell you the map is trustworthy.
    """
    px, values = np.asarray(px, float), np.asarray(values, float)
    slope, intercept = np.polyfit(px, values, 1)
    return intercept, slope, (values - (intercept + slope * px)) / slope


# --------------------------------------------------------------------------- cli


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("rectify", help="warp the page so the plot frame is a true rectangle")
    r.add_argument("--page", type=Path, required=True)
    r.add_argument("--out", type=Path, required=True)

    t = sub.add_parser("ticks", help="measure printed tick positions on a rectified page")
    t.add_argument("--page", type=Path, required=True)
    t.add_argument("--frame", required=True, help="left,top,right,bottom from `rectify`")
    t.add_argument("--minlen", type=int, default=6, help="shortest protrusion counted as a tick")

    a = ap.parse_args(argv)

    if a.cmd == "rectify":
        left, top, right, bottom = rectify(a.page, a.out)
        print(f"wrote {a.out}")
        print(f"--frame {left},{top},{right},{bottom}")
        print("# check with `digitize.py axes`; where they disagree, the ticks decide")
        return 0

    if a.cmd == "ticks":
        frame = tuple(int(v) for v in a.frame.split(","))
        for edge, found in ticks(a.page, frame, a.minlen).items():
            print(f"{edge}: {len(found)} groups (centre, protrusion)")
            print("   ", [(round(c, 1), m) for c, m in found])
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
