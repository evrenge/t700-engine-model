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


FRAME_MIN_RUN_FRAC = 0.08
"""How much of the plot width a row must span in one unbroken run to be a frame line.

Measured: the frames run 144-546 px of 1688 on pdf p.40, the faint and tilted one, and
essentially the full width on pp.41-42. The two things that are not frames run 20 px (the
caption's longest letter stem) and about 25 px (a data glyph). 8 % is 135 px on p.40 --
above both by a factor of five, below the worst real frame by the same."""


def find_frame(ink: np.ndarray) -> tuple[int, int, int, int]:
    """The single plot frame. Same discriminator as everywhere else in this project:
    a frame line spans the plot, text and data do not."""
    h, w = ink.shape
    m0, m1 = int(w * 0.06), int(w * 0.97)
    cols = np.array([longest_run(ink[:, x]) if m0 <= x < m1 else 0 for x in range(w)])
    strong = np.flatnonzero(cols > 0.6 * cols.max())
    left, right = int(strong.min()), int(strong.max())

    # **Longest run per row, not fill fraction.** Fill was the criterion until 2026-09-14
    # and it put Figure 6's bottom datum on the figure's own CAPTION. That page's frames
    # are faint and tilted by ~16 px across the width, so the line's ink spreads over
    # twenty rows and no single row exceeds 0.24 fill -- while the caption, a dense line of
    # text, reaches 0.42 and won the ladder. The axis was then 2958-330 px tall against a
    # true 2769-330, so every digitized NG read high: +1.2 %NG at 90 and +2.8 %NG at 70.
    #
    # Longest run separates them structurally rather than by a threshold: a frame line is
    # continuous ink for hundreds of pixels (144-546 px on the worst page, near the full
    # width on pp.41-42), while the caption's longest run is 20 px -- one letter stem --
    # and a data glyph's is about 25. The 8 % floor sits an order of magnitude above both.
    span = right - left
    runs = np.array([longest_run(ink[y, left : right + 1]) for y in range(h)])
    rows = np.flatnonzero(runs > FRAME_MIN_RUN_FRAC * span)
    if rows.size:
        groups = np.split(rows, np.flatnonzero(np.diff(rows) > 6) + 1)
        if len(groups) >= 2:
            return (
                left,
                int(round(groups[0].mean())),
                right,
                int(round(groups[-1].mean())),
            )
    raise RuntimeError("could not find two horizontal frame lines")


FRAME_FIT_HALF_BAND_PX = 30
"""Half-width of the first-pass search band around a nominal frame position, in pixels.

Wide enough to contain the whole lean (15-19 px end to end on pp.41-42) plus the line's own
thickness; a second pass narrows to 6 px about the first fit."""


def _fit_line(samples: list[tuple[float, float]]) -> tuple[float, float]:
    a, b = np.polyfit([s[0] for s in samples], [s[1] for s in samples], 1)
    return float(a), float(b)


def fit_frames(ink: np.ndarray, frame: tuple[int, int, int, int]) -> dict:
    """Fit all four frame lines, because the plot is a parallelogram and not a rectangle.

    **These pages are sheared.** The left and right frames of pdf pp.41-42 lean by 15.7 to
    18.7 px across the plot height -- 0.92 to 1.09 percent of plot width -- and pdf p.40's
    lean the other way by 7.6 to 13.2 px. The top and bottom frames are level to 0.10
    percent of height, so unlike Figures 9-10 (open question #63) the distortion here is in
    **x**, not y.

    Uncorrected it is worth up to 8 lbm/hr on the fuel-flow axis of Figures 6 and 7 -- 17 px
    of 1710 across an 800 lbm/hr span, so 1.1 percent at 750 lbm/hr and 3.5 percent at 224,
    the low-power end where the Table B.1 vs Figure 7 disagreement of open question #46
    lives -- and 0.40 %NG on Figure 8's abscissa.

    `find_frame` returns the frame's outermost extent, so with a positive lean `left` is the
    left frame at the TOP and `right` is the right frame at the BOTTOM. The width it implies
    is therefore too large by the lean, and its origin belongs to a different row than its
    end. Fitting each line along its length removes both errors at once.

    Returns slope/intercept for each: horizontals as `y = a*x + b`, verticals as
    `x = a*y + b`.
    """
    left, top, right, bottom = frame
    out = {}
    for name, nominal, vertical in (
        ("top", top, False),
        ("bottom", bottom, False),
        ("left", left, True),
        ("right", right, True),
    ):
        half = FRAME_FIT_HALF_BAND_PX
        fitted = (0.0, float(nominal))
        for _pass in (0, 1):
            samples: list[tuple[float, float]] = []
            if vertical:
                for y in range(top + 5, bottom - 4, 4):
                    centre = fitted[0] * y + fitted[1]
                    lo = int(round(centre - half))
                    band = ink[y, max(lo, 0) : lo + 2 * half + 1]
                    idx = np.flatnonzero(band)
                    if idx.size and idx.size <= 10:
                        samples.append((float(y), lo + float(idx.mean())))
            else:
                for x in range(left + 5, right - 4, 4):
                    centre = fitted[0] * x + fitted[1]
                    lo = int(round(centre - half))
                    band = ink[max(lo, 0) : lo + 2 * half + 1, x]
                    idx = np.flatnonzero(band)
                    if idx.size and idx.size <= 10:
                        samples.append((float(x), lo + float(idx.mean())))
            if len(samples) < 50:
                break
            fitted = _fit_line(samples)
            half = 6
        out[name] = fitted
    return out


TICK_DEPTH_PX = 14
TICK_MIN_FILL = 0.6
TICK_WIDTH_PX = (2, 10)
"""What a tick is, in native-scan pixels: a stroke reaching `TICK_DEPTH_PX` inside the
frame, at least `TICK_MIN_FILL` of it inked, and between 2 and 10 px across.

The width bound is the discriminator that matters. Without it, a data glyph sitting near a
frame is a "tick": Figure 6's curve runs low on the left, and its bottom frame returned 21
detections at 13-unit spacings -- glyph debris -- against the 5 real ticks its top frame
returns. Printed glyphs are 15-25 px across on these pages and ticks are 3-6."""

TICK_RESIDUAL_TOL = {6: 2.0, 7: 2.0, 8: 0.10}
"""How far a mapped major tick may sit from its round value, in data units, before the
calibration is rejected: 2 lbm/hr of an 800 lbm/hr span, 0.10 %NG of a 40 %NG span. Both
are about 0.25 percent of the axis.

**This is a held-out check**, the same shape as the one in `digitize_a1.py`: the two frames
of each axis carry printed values and pin the map, and the interior major ticks are then
predicted rather than fitted. Measured before the shear correction the residuals ran rms
4.1-4.8 lbm/hr and 0.20-0.27 %NG; after it, 0.43-0.64 and 0.03. The tolerance is set
between those, so it fails the calibration this tool had until 2026-09-14 and passes the
one it has."""

TICK_GRID = {6: 100.0, 7: 100.0, 8: 5.0}
"""Major tick spacing on each figure's abscissa, read off the printed axis labels."""

TICK_GRID_Y = {6: 2.5, 7: 200.0, 8: 20.0}
TICK_RESIDUAL_TOL_Y = {6: 0.25, 7: 8.0, 8: 1.5}
"""The same check on the ordinate, where the ticks are minor as well as major -- Figure 6
prints one every 2.5 %NG against labels every 5.

**This is the check that caught the caption.** Figure 6's sixteen left-frame ticks read
102.91, 100.58, 98.23, 95.90, 93.59 ... under the old calibration: a 2.33 spacing, on no
grid at all, rms 0.9 %NG from the nearest 2.5. Under the fitted frames they read 102.48,
99.97, 97.45, 94.94, 92.46 -- the 2.5 grid to **rms 0.04 %NG** over sixteen held-out ticks.
Nothing else available on that page could have distinguished a bad y datum from a real
disagreement with Table B.1, and the disagreement it leaves is real and larger than before.

Tolerances sit between the two measurements, so each fails the calibration this tool had
until 2026-09-14 and passes the one it has: Figure 7 measured rms 1.95-3.45 of a 2200 span
and Figure 8 0.44-0.47 of 240, both already sound because those pages print crisp frames."""


def frame_ticks_v(
    ink: np.ndarray, line: tuple[float, float], lo: int, hi: int, inward: int
) -> np.ndarray:
    """Rows of the ticks on one vertical frame. `inward` is +1 when the interior is to the
    right of the line (the left frame) and -1 when it is to the left (the right frame)."""
    fill = []
    for y in range(lo + 4, hi - 3):
        x0 = int(round(line[0] * y + line[1]))
        a = x0 + 2 if inward > 0 else x0 - TICK_DEPTH_PX - 1
        band = ink[y, a : a + TICK_DEPTH_PX]
        if band.size < TICK_DEPTH_PX:
            fill.append(0)
            continue
        adjacent = band[0] if inward > 0 else band[-1]
        fill.append(int(band.sum()) if adjacent else 0)
    on = np.flatnonzero(np.array(fill) >= TICK_DEPTH_PX * TICK_MIN_FILL)
    if on.size == 0:
        return np.zeros(0)
    out = []
    for g in np.split(on, np.flatnonzero(np.diff(on) > 3) + 1):
        if TICK_WIDTH_PX[0] <= g.size <= TICK_WIDTH_PX[1]:
            out.append(g.mean() + lo + 4)
    return np.array(out)


def frame_ticks(
    ink: np.ndarray, line: tuple[float, float], lo: int, hi: int, inward: int
) -> np.ndarray:
    """Columns of the major ticks on one horizontal frame.

    `inward` is +1 for a frame whose interior lies below it (the top frame) and -1 for one
    whose interior lies above (the bottom frame).
    """
    fill = []
    for x in range(lo + 4, hi - 3):
        y0 = int(round(line[0] * x + line[1]))
        # The band is CONTIGUOUS with the frame. A tick is attached to it; a data glyph
        # that merely sits near it is not, and Figure 6's curve runs along its own bottom
        # frame. Requiring ink in the row adjacent to the frame is what separates them:
        # without it that frame returned 21 "ticks" at irregular spacings against the 5
        # its top frame returns.
        a = y0 + 2 if inward > 0 else y0 - TICK_DEPTH_PX - 1
        band = ink[a : a + TICK_DEPTH_PX, x]
        adjacent = band[0] if inward > 0 else band[-1]
        fill.append(int(band.sum()) if adjacent else 0)
    on = np.flatnonzero(np.array(fill) >= TICK_DEPTH_PX * TICK_MIN_FILL)
    if on.size == 0:
        return np.zeros(0)
    out = []
    for g in np.split(on, np.flatnonzero(np.diff(on) > 3) + 1):
        if TICK_WIDTH_PX[0] <= g.size <= TICK_WIDTH_PX[1]:
            out.append(g.mean() + lo + 4)
    return np.array(out)


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
    failures: list[str] = []
    for f in FIGS:
        ink = native_bitmap(f.page, OUT / f"p{f.page}.pbm")
        fr = find_frame(ink)
        left, top, right, bottom = fr
        series, stats = classify(ink, fr)

        fits = fit_frames(ink, fr)

        # The parallelogram inverse. `u` is the fraction across at this row, `v` the
        # fraction down at this column; both frames of each pair are fitted, so neither
        # axis is read against an edge position that belongs to a different part of the
        # plot. Exact for a parallelogram, and the four residuals are 0.48-1.20 px.
        def x_of(px, py, f=f, fits=fits):
            lo = fits["left"][0] * np.asarray(py, float) + fits["left"][1]
            hi = fits["right"][0] * np.asarray(py, float) + fits["right"][1]
            u = (np.asarray(px, float) - lo) / (hi - lo)
            return f.x_lo + u * (f.x_hi - f.x_lo)

        def y_of(py, px, f=f, fits=fits):
            hi = fits["top"][0] * np.asarray(px, float) + fits["top"][1]
            lo = fits["bottom"][0] * np.asarray(px, float) + fits["bottom"][1]
            v = (np.asarray(py, float) - hi) / (lo - hi)
            return f.y_hi + v * (f.y_lo - f.y_hi)

        # Held-out tick check. The frames pinned the map; the interior ticks did not.
        grid = TICK_GRID[f.num]
        tick_report = []
        for side, inward in (("top", +1), ("bottom", -1)):
            cols = frame_ticks(ink, fits[side], left, right, inward)
            if cols.size < 3:
                tick_report.append(f"{side}: {cols.size} ticks, too few to check")
                continue
            rows = fits[side][0] * cols + fits[side][1]
            vals = x_of(cols, rows)
            # Which grid the ticks are on is read off the ticks, not assumed: Figure 6's
            # bottom frame carries MINOR ticks (45 of them) where its top frame carries 8
            # major ones, and a check that insisted on the major spacing called the minor
            # ticks a calibration failure. The coarsest spacing that fits is taken, and
            # nothing finer than a fifth of the major spacing is allowed to count -- at
            # that point a random set of columns would pass.
            best = None
            for div in (1, 2, 4, 5):
                sp = grid / div
                r = vals - np.round(vals / sp) * sp
                rms = float(np.sqrt(np.mean(r**2)))
                if rms < TICK_RESIDUAL_TOL[f.num]:
                    best = (sp, rms, float(np.abs(r).max()))
                    break
            if best is None:
                sp = grid
                r = vals - np.round(vals / sp) * sp
                rms = float(np.sqrt(np.mean(r**2)))
                tick_report.append(f"{side}: {cols.size} ticks, rms {rms:.3f} FAIL")
                failures.append(
                    f"Figure {f.num} {side} frame: {cols.size} ticks land on no round grid "
                    f"down to {grid / 5:g} -- best rms {rms:.3f} against a "
                    f"{TICK_RESIDUAL_TOL[f.num]} tolerance. The x calibration is wrong, "
                    f"not the ticks"
                )
            else:
                sp, rms, mx = best
                tick_report.append(
                    f"{side}: {cols.size} ticks on a {sp:g} grid, rms {rms:.3f} max {mx:.3f} ok"
                )

        for side, inward in (("left", +1), ("right", -1)):
            rws = frame_ticks_v(ink, fits[side], top, bottom, inward)
            if rws.size < 3:
                tick_report.append(f"{side}: {rws.size} ticks, too few to check")
                continue
            cls = fits[side][0] * rws + fits[side][1]
            vals = y_of(rws, cls)
            gy = TICK_GRID_Y[f.num]
            r = vals - np.round(vals / gy) * gy
            rms = float(np.sqrt(np.mean(r**2)))
            ok = rms < TICK_RESIDUAL_TOL_Y[f.num]
            tick_report.append(
                f"{side}: {rws.size} ticks on a {gy:g} grid, rms {rms:.4f} "
                f"max {np.abs(r).max():.4f} {'ok' if ok else 'FAIL'}"
            )
            if not ok:
                failures.append(
                    f"Figure {f.num} {side} frame: {rws.size} ticks miss the {gy:g} grid by "
                    f"rms {rms:.4f} against a {TICK_RESIDUAL_TOL_Y[f.num]} tolerance -- the "
                    f"y calibration is wrong, not the ticks"
                )

        print(f"\nFigure {f.num} (pdf p.{f.page})  frame {left},{top},{right},{bottom}")
        print(
            f"  shear: left {fits['left'][0] * (bottom - top):+.2f} px, "
            f"right {fits['right'][0] * (bottom - top):+.2f} px across the height"
        )
        for line in tick_report:
            print(f"  tick check {line}")
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
                fh.write(
                    f"# shear corrected: left frame leans "
                    f"{fits['left'][0] * (bottom - top):+.2f} px across the plot height, "
                    f"right {fits['right'][0] * (bottom - top):+.2f} px; top frame "
                    f"{fits['top'][0] * (right - left):+.2f} px across the width, bottom "
                    f"{fits['bottom'][0] * (right - left):+.2f} px. The page is sheared; "
                    f"see `fit_frames`.\n"
                )
                fh.write("# NOT transcribed. Values carry read error.\n")
                fh.write(f"{f.xkey},{f.ykey}\n")
                for px, py in pts:
                    fh.write(f"{x_of(px, py):.4f},{y_of(py, px):.4f}\n")
    if failures:
        print("\nTICK CHECK FAILED -- the axis calibration does not predict the printed ticks:")
        for line in failures:
            print(f"  {line}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
