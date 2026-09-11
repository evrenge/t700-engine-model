"""Digitize Figures 6, 7 and 8 -- the steady-state trim sweeps.

These are the report's operating-range validation [pdf pp.40-42], and they matter more
than their size suggests. Table B.1 gives three trim points; these give about thirty each,
across the whole fuel-flow range. Three points against eight component maps is
underdetermined -- a deviation cannot be attributed. Ninety points over-determine it, so a
map error shows as a *trend in a region* rather than as one number nobody can localize.

    Figure 6, p.40   gas generator speed, %      vs  fuel flow, lb/hr
    Figure 7, p.41   shaft horsepower            vs  fuel flow, lb/hr
    Figure 8, p.42   station 3 static pressure   vs  gas generator speed, %

Each carries **three marker series**, keyed by an in-plot legend:

    *  REAL-TIME MODEL             -- Ballin's model. Our replication target.
    ^  G.E. UNBALANCED TORQUE MODEL
    x  G.E. STATUS-81 MODEL

The two GE series are what Ballin was validating *against*; they are not what we
reproduce. They are extracted anyway, to separate files, because the gap between his model
and them is itself information -- it is the error he accepted.

Glyph separation is by shape, which is what distinguishes them: a filled disc inks most of
its bounding box, while an outline triangle and a cross ink roughly a third of theirs.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from digitize_native import native_bitmap  # noqa: E402

OUT = Path("validation/out/digitize/fig678")
REF = Path("data/reference")


@dataclass(frozen=True)
class Fig:
    num: int
    page: int
    xkey: str
    xlabel: str
    x_lo: float
    x_hi: float
    ykey: str
    ylabel: str
    y_lo: float
    y_hi: float


FIGS = [
    Fig(
        6,
        40,
        "wf_pph",
        "FUEL FLOW, lb/hr",
        100.0,
        900.0,
        "ng_pct",
        "GAS GENERATOR SPEED, percent",
        65.0,
        105.0,
    ),
    Fig(
        7,
        41,
        "wf_pph",
        "FUEL FLOW, lb/hr",
        100.0,
        900.0,
        "shp",
        "SHAFT HORSEPOWER (no unit printed)",
        -200.0,
        2000.0,
    ),
    Fig(
        8,
        42,
        "ng_pct",
        "GAS GENERATOR SPEED, percent",
        65.0,
        105.0,
        "ps3_psia",
        "STATION 3 STATIC PRESSURE, psia",
        40.0,
        280.0,
    ),
]

SERIES = {
    "realtime": "REAL-TIME MODEL -- Ballin's model, the replication target",
    "ge_unbalanced": "G.E. UNBALANCED TORQUE MODEL",
    "ge_status81": "G.E. STATUS-81 MODEL",
}


def longest_run(col: np.ndarray) -> int:
    best = cur = 0
    for v in col:
        cur = cur + 1 if v else 0
        if cur > best:
            best = cur
    return best


def find_frame(ink: np.ndarray) -> tuple[int, int, int, int]:
    """The single plot frame. Same discriminator as everywhere else in this project:
    a frame line spans the plot, text and data do not."""
    h, w = ink.shape
    m0, m1 = int(w * 0.06), int(w * 0.97)
    cols = np.array([longest_run(ink[:, x]) if m0 <= x < m1 else 0 for x in range(w)])
    strong = np.flatnonzero(cols > 0.6 * cols.max())
    left, right = int(strong.min()), int(strong.max())

    frac = ink[:, left : right + 1].mean(axis=1)
    # these pages print fainter frames than the appendix: p.40's strongest row
    # reaches only 0.437 fill, so the ladder has to go lower than elsewhere.
    for thresh in (0.8, 0.7, 0.6, 0.5, 0.4, 0.32, 0.26, 0.20):
        rows = np.flatnonzero(frac > thresh)
        if rows.size == 0:
            continue
        groups = np.split(rows, np.flatnonzero(np.diff(rows) > 3) + 1)
        if len(groups) >= 2:
            return left, int(round(groups[0].mean())), right, int(round(groups[-1].mean()))
    raise RuntimeError("could not find two horizontal frame lines")


def classify(ink: np.ndarray, frame: tuple[int, int, int, int], margin: int = 22):
    """Find every glyph inside the frame and sort it by shape.

    A filled disc inks about 78 % of its bounding box; an outline triangle and a cross ink
    roughly a third. That single number separates the replication target from the two GE
    series cleanly, and it is a property of the printed glyph rather than of any threshold
    we chose.

    Triangle against cross is the harder split, and is decided by where the ink sits: a
    triangle has a solid horizontal base, so its bottom row is nearly full; a cross has its
    ink at the corners and an empty bottom row.
    """
    left, top, right, bottom = frame
    sub = ink[top + margin : bottom - margin, left + margin : right - margin]
    lab, n = ndimage.label(sub, structure=np.ones((3, 3)))

    # First pass: every blob that could be a glyph, with its shape statistics.
    cand = []
    for i, sl in enumerate(ndimage.find_objects(lab), start=1):
        blob = lab[sl] == i
        hgt, wid = blob.shape
        area = int(blob.sum())
        if area < 20 or hgt < 6 or wid < 6 or hgt > 40 or wid > 40:
            continue
        if max(hgt, wid) / min(hgt, wid) > 1.8:
            continue
        cy, cx = ndimage.center_of_mass(blob)
        cand.append(
            {
                "x": cx + sl[1].start + left + margin,
                "y": cy + sl[0].start + top + margin,
                "fill": area / (hgt * wid),
                "base": blob[-2:, :].mean(),
            }
        )

    # Reject the in-plot legend without hard-coding where it is. The legend is *text*, and
    # text is what a data marker is not: characters sit shoulder to shoulder on a common
    # baseline, so each has several neighbours within a character's width at the same
    # height. Markers on a scatter plot do not. Counting neighbours separates them by a
    # property of what they are rather than by a rectangle someone measured once.
    X = np.array([c["x"] for c in cand])
    Y = np.array([c["y"] for c in cand])
    kept = []
    for cdt in cand:
        close = np.flatnonzero((np.abs(X - cdt["x"]) < 60.0) & (np.abs(Y - cdt["y"]) < 12.0))
        if close.size - 1 >= 3:
            continue  # three or more companions on one baseline: a word, not a datum
        kept.append(cdt)

    out: dict[str, list[dict]] = {k: [] for k in SERIES}
    for cdt in kept:
        if cdt["fill"] > 0.62:
            out["realtime"].append(cdt)
        elif cdt["base"] > 0.55:
            out["ge_unbalanced"].append(cdt)
        else:
            out["ge_status81"].append(cdt)

    stats = [(c["fill"], c["base"]) for c in kept]
    final = {k: _backbone(v) for k, v in out.items()}
    return {
        k: (np.array([(c["x"], c["y"]) for c in v]) if v else np.zeros((0, 2)))
        for k, v in final.items()
    }, stats


def _backbone(cand: list[dict]) -> list[dict]:
    """Keep only the monotone backbone of one series.

    All three figures plot a monotone relationship -- more fuel gives more speed, more
    speed gives more pressure -- so every genuine marker of a series lies on a
    non-decreasing curve. The longest non-decreasing subsequence *is* that curve, and what
    falls off it is contamination.

    This is what finally worked, after a robust polynomial fit did not. What it removes is
    exactly the two classes visible on the page: the legend's own glyph samples, which sit
    above the data near the top of the plot, and scan specks in the empty lower right. The
    inventory independently warns of one such speck on Figure 6 "with no legend shape ...
    scan dirt, not data".

    It is applied **per series**, not to the pooled set: pooling keeps one chain and so
    quietly decimates whichever series is sparser.

    It is safe because it appeals to physics rather than to a threshold -- it cannot delete
    a real point unless that point breaks monotonicity, and none does.
    """
    if len(cand) < 6:
        return cand
    cand = sorted(cand, key=lambda c: c["x"])
    ys = [c["y"] for c in cand]  # pixel y increases downward; every figure rises
    n = len(ys)
    best = [1] * n
    prev = [-1] * n
    for i in range(n):
        for j in range(i):
            if ys[j] >= ys[i] and best[j] + 1 > best[i]:
                best[i] = best[j] + 1
                prev[i] = j
    idx = int(np.argmax(best))
    chain = []
    while idx != -1:
        chain.append(idx)
        idx = prev[idx]
    return [cand[i] for i in reversed(chain)]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    REF.mkdir(parents=True, exist_ok=True)
    for f in FIGS:
        ink = native_bitmap(f.page, OUT / f"p{f.page}.pbm")
        fr = find_frame(ink)
        left, top, right, bottom = fr
        series, stats = classify(ink, fr)

        def x_of(px, f=f, left=left, right=right):
            span = f.x_hi - f.x_lo
            return f.x_lo + (np.asarray(px, float) - left) / (right - left) * span

        def y_of(py, f=f, top=top, bottom=bottom):
            g = (np.asarray(py, float) - top) / (bottom - top)
            return f.y_hi + g * (f.y_lo - f.y_hi)

        print(f"\nFigure {f.num} (pdf p.{f.page})  frame {left},{top},{right},{bottom}")
        if stats:
            fills = np.array([s[0] for s in stats])
            print(f"  {len(stats)} glyphs, fill ratio {fills.min():.2f}..{fills.max():.2f}")
        for key, pts in series.items():
            print(f"    {key:14} {pts.shape[0]:3d}")
            if pts.shape[0] == 0:
                continue
            path = REF / f"fig{f.num:02d}_{key}.csv"
            with path.open("w") as fh:
                fh.write(f"# source: TM-100991 pdf p.{f.page}, Figure {f.num}\n")
                fh.write(f"# quantity: {f.ylabel} vs {f.xlabel}; series: {SERIES[key]}\n")
                fh.write(
                    "# method: digitized -- native 300 dpi 1-bit scan, "
                    "tools/digitize_fig678.py; glyphs separated by shape\n"
                )
                fh.write(f"# {f.xkey}: {f.xlabel}, axis {f.x_lo} to {f.x_hi}\n")
                fh.write(f"# {f.ykey}: {f.ylabel}, axis {f.y_lo} to {f.y_hi}\n")
                fh.write(f"# points: {pts.shape[0]}\n")
                fh.write("# NOT transcribed. Values carry read error.\n")
                fh.write(f"{f.xkey},{f.ykey}\n")
                for px, py in pts:
                    fh.write(f"{x_of(px):.4f},{y_of(py):.4f}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
