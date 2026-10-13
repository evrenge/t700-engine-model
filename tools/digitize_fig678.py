"""Digitize Figures 6, 7 and 8 -- the steady-state trim sweeps [pdf pp.40-42].

About thirty points each across the whole fuel-flow range, against Table B.1's three, so a
map error shows as a trend in a region rather than as one number nobody can localize.

    Figure 6, p.40   gas generator speed, %      vs  fuel flow, lb/hr
    Figure 7, p.41   shaft horsepower            vs  fuel flow, lb/hr
    Figure 8, p.42   station 3 static pressure   vs  gas generator speed, %

Each carries three marker series, keyed by an in-plot legend:

    *  REAL-TIME MODEL             -- Ballin's model. Our replication target.
    ^  G.E. UNBALANCED TORQUE MODEL
    x  G.E. STATUS-81 MODEL

The two GE series are what Ballin validated against, not what we reproduce; they go to
separate files because the gap between his model and them is the error he accepted.
Glyphs are separated by shape (`classify`).
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure.axis import blend_value, choose_anchored, fit_anchored  # noqa: E402
from figure.frame import band_line, longest_run, longest_run_columns, make_sv  # noqa: E402
from figure.raster import native_bitmap  # noqa: E402

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


FRAME_MIN_RUN_FRAC = 0.08
"""How much of the plot width a row must span in one unbroken run to be a frame line.

Frames run 144-546 px of 1688 on the faint, tilted p.40 and nearly the full width on
pp.41-42; the caption's longest letter stem is 20 px and a data glyph about 25. 8 % (135 px
on p.40) sits a factor of five from both."""


def find_frame(ink: np.ndarray) -> tuple[int, int, int, int]:
    """The single plot frame. Same discriminator as everywhere else in this project:
    a frame line spans the plot, text and data do not."""
    h = ink.shape[0]
    left, right = longest_run_columns(ink, 0.06, 0.97, 0.6)

    # Longest run per row, not fill fraction: p.40's faint, tilted frame spreads its ink
    # over twenty rows, and by fill the caption (dense text) wins. See FRAME_MIN_RUN_FRAC.
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


def fit_frames(ink: np.ndarray, frame: tuple[int, int, int, int]) -> dict:
    """Fit all four frame lines, because the plot is a parallelogram and not a rectangle.

    These pages are sheared in x, not y (unlike Figures 9-10, open question #63): the
    vertical frames lean up to ~19 px across the plot height while top and bottom are level.
    Uncorrected that is up to 8 lbm/hr on the fuel-flow axis, worst at the low-power end
    where open question #46 lives. `find_frame` returns the outermost extent, so its width
    is too large by the lean; fitting each line along its length removes that.

    Returns slope/intercept for each: horizontals as `y = a*x + b`, verticals as
    `x = a*y + b`.
    """
    left, top, right, bottom = frame
    out = {}
    for name, nominal, along, lo, hi in (
        ("top", top, 1, left + 5, right - 4),
        ("bottom", bottom, 1, left + 5, right - 4),
        ("left", left, 0, top + 5, bottom - 4),
        ("right", right, 0, top + 5, bottom - 4),
    ):
        # a pass short of 50 samples keeps the previous pass's line, or the nominal one
        out[name], _, _ = band_line(
            ink, along, nominal, lo, hi, 4, (FRAME_FIT_HALF_BAND_PX, 6), 10, 50
        )
    return out


TICK_DEPTH_PX = 14
TICK_MIN_FILL = 0.6
TICK_WIDTH_PX = (2, 10)
"""What a tick is, in native-scan pixels: a stroke reaching `TICK_DEPTH_PX` inside the
frame, at least `TICK_MIN_FILL` of it inked, and between 2 and 10 px across. The width
bound is what rejects a data glyph near a frame: glyphs are 15-25 px across, ticks 3-6."""

TICK_RESIDUAL_TOL = {6: 2.0, 7: 2.0, 8: 0.10}
"""How far a mapped major tick may sit from its round value, in data units, before the
calibration is rejected: about 0.25 % of each axis. A held-out check: the two printed frame
values pin the map and each interior tick is predicted without itself (`choose_anchored`).
Set between the residuals without the shear correction (rms 4.1-4.8 lbm/hr, 0.20-0.27 %NG)
and with it."""

TICK_GRID = {6: 100.0, 7: 100.0, 8: 5.0}
"""Major tick spacing on each figure's abscissa, read off the printed axis labels."""

TICK_GRID_Y = {6: 2.5, 7: 200.0, 8: 20.0}
TICK_RESIDUAL_TOL_Y = {6: 0.25, 7: 8.0, 8: 1.5}
"""The same check on the ordinate, where the ticks are minor as well as major -- Figure 6
prints one every 2.5 %NG against labels every 5. On Figure 6 it is the only thing that tells
a wrong y datum (ticks off the 2.5 grid by rms 0.9 %NG) from a real disagreement with
Table B.1; the tolerances sit between that and the fitted frames' residuals."""


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
        # The band is contiguous with the frame: a tick is attached to it, a data glyph
        # near it is not (Figure 6's curve runs along its own bottom frame).
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


# ---------------------------------------------------------------- the printed scale bows

BOW_ORDERS = (1, 2, 3, 4)
"""Orders of `axis.fit_anchored` tried per edge, chosen by leave-one-out and never by eye:
order 1 is the straight map between the two printed frame values, and 2-4 add a bow
`t(1-t)P(t)` with `P` of degree 0-2. The `t(1-t)` factor vanishes at both frames, so the two
printed frame values stay exact by construction however the ticks come out; the ticks only
say how the map bends in between, which is the one thing two endpoints cannot know."""


def classify(ink: np.ndarray, frame: tuple[int, int, int, int], margin: int = 22):
    """Find every glyph inside the frame and sort it by shape.

    A filled disc inks about 78 % of its bounding box, an outline triangle or a cross about
    a third. Triangle against cross: a triangle's solid base fills its bottom row, a cross
    has its ink at the corners and an empty bottom row.
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

    # Reject the in-plot legend without hard-coding where it is: text characters sit
    # shoulder to shoulder on a baseline, so each has several neighbours at its height.
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

    All three figures plot a monotone relationship, so every genuine marker lies on the
    longest non-decreasing subsequence; what falls off it (the legend's own glyph samples,
    scan specks) is contamination. Applied per series: pooling would keep one chain and
    decimate the sparser series.
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
        ink = native_bitmap(f.page)
        fr = find_frame(ink)
        left, top, right, bottom = fr
        series, stats = classify(ink, fr)

        fits = fit_frames(ink, fr)

        # The parallelogram inverse, `frame.make_sv` on the four fitted lines: `u` is the
        # fraction across at this row, `v` the fraction DOWN at this column (make_sv's runs up).
        to_sv, _ = make_sv(fits)

        def u_of(px, py, to_sv=to_sv):
            return to_sv(px, py)[0]

        def v_of(py, px, to_sv=to_sv):
            return 1.0 - to_sv(px, py)[1]

        # Stage one: the straight map, used only to tell each tick which round value it is.
        def x_lin(px, py, f=f):
            return f.x_lo + u_of(px, py) * (f.x_hi - f.x_lo)

        def y_lin(py, px, f=f):
            return f.y_hi + v_of(py, px) * (f.y_lo - f.y_hi)

        # Stage two: what the ticks say about the bend between the frames. Each edge gets a
        # bow term exactly when leave-one-out says it predicts ticks it was not shown;
        # nothing is bent because it looked bent.
        tick_report = []
        maps: dict[str, np.ndarray] = {}
        held_out: dict[str, float] = {}
        for side, inward in (("top", +1), ("bottom", -1), ("left", +1), ("right", -1)):
            horiz = side in ("top", "bottom")
            # each edge's map runs from the value at t = 0 to the value at t = 1, and stays
            # straight unless its ticks say otherwise
            lo_v, hi_v = (f.x_lo, f.x_hi) if horiz else (f.y_hi, f.y_lo)
            maps[side] = fit_anchored(np.zeros(0), np.zeros(0), 1, lo_v, hi_v)
            if horiz:
                pos = frame_ticks(ink, fits[side], left, right, inward)
                if pos.size < 3:
                    tick_report.append(f"{side}: {pos.size} ticks, too few to check")
                    continue
                other = fits[side][0] * pos + fits[side][1]
                vals = x_lin(pos, other)
                t = u_of(pos, other)
                grid, tol = TICK_GRID[f.num], TICK_RESIDUAL_TOL[f.num]
            else:
                pos = frame_ticks_v(ink, fits[side], top, bottom, inward)
                if pos.size < 3:
                    tick_report.append(f"{side}: {pos.size} ticks, too few to check")
                    continue
                other = fits[side][0] * pos + fits[side][1]
                vals = y_lin(pos, other)
                t = v_of(pos, other)
                grid, tol = TICK_GRID_Y[f.num], TICK_RESIDUAL_TOL_Y[f.num]

            # The grid is read off the ticks (Figure 6's bottom frame carries minors, its
            # top majors): the coarsest spacing that fits, never finer than a fifth of the
            # major, below which a random set of columns would pass.
            sp = None
            for div in (1, 2, 4, 5):
                cand = grid / div
                r = vals - np.round(vals / cand) * cand
                if float(np.sqrt(np.mean(r**2))) < tol:
                    sp = cand
                    break
            if sp is None:
                r = vals - np.round(vals / grid) * grid
                rms = float(np.sqrt(np.mean(r**2)))
                tick_report.append(f"{side}: {pos.size} ticks, rms {rms:.4f} FAIL")
                failures.append(
                    f"Figure {f.num} {side} frame: {pos.size} ticks land on no round grid "
                    f"down to {grid / 5:g} -- best rms {rms:.4f} against a {tol} "
                    f"tolerance. The calibration is wrong, not the ticks"
                )
                continue

            want = np.round(vals / sp) * sp
            resid = want - vals
            maps[side], loo, fitted, order = choose_anchored(t, want, lo_v, hi_v, BOW_ORDERS)
            held_out[side] = loo
            span = abs(hi_v - lo_v)
            shape = ""
            if resid.std() > 0:
                b = t * (1.0 - t)
                if b.std() > 0:
                    shape = f", bow correlation {float(np.corrcoef(resid, b)[0, 1]):+.2f}"
            straight = float(np.sqrt(np.mean(resid**2)))
            tick_report.append(
                f"{side}: {pos.size} ticks on a {sp:g} grid, straight rms {straight:.4f} "
                f"({straight / span * 100:.3f} %span){shape} -> "
                + (
                    f"degree {order - 2} bow, rms {fitted:.4f}, held-out {loo:.4f}"
                    if order > 1
                    else f"no bow (held-out {loo:.4f})"
                )
            )
            if loo > tol:
                failures.append(
                    f"Figure {f.num} {side} frame: even after the bow, {pos.size} ticks are "
                    f"predicted only to rms {loo:.4f} when held out, against a {tol} "
                    f"tolerance. The calibration is wrong, not the ticks"
                )

        # Stage three: the map the data goes through. Each axis carries the bow of its near
        # edge, blended to that of its far edge across the plot: the two ticked edges of
        # one axis are not the same scale. The blend weight is clipped to the frame.
        def x_of(px, py, maps=maps):
            w = np.clip(v_of(py, px), 0.0, 1.0)
            return blend_value((maps["top"], maps["bottom"]), u_of(px, py), w)

        def y_of(py, px, maps=maps):
            w = np.clip(u_of(px, py), 0.0, 1.0)
            return blend_value((maps["left"], maps["right"]), v_of(py, px), w)

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
                hx = [held_out[k] for k in ("top", "bottom") if k in held_out]
                hy = [held_out[k] for k in ("left", "right") if k in held_out]
                fh.write(
                    "# read error: after calibration the axis maps predict a tick they "
                    "were not shown to rms "
                    + (f"{max(hx):.4g} in {f.xkey}" if hx else "(x: no tick check)")
                    + " and "
                    + (f"{max(hy):.4g} in {f.ykey}" if hy else "(y: no tick check)")
                    + ", worst edge of the two. That is the scale to judge a deviation "
                    "against: a disagreement smaller than this one is not resolved by the "
                    "figure.\n"
                )
                fh.write("# NOT transcribed. Values carry read error.\n")
                fh.write(f"{f.xkey},{f.ykey}\n")
                for xv, yv in zip(
                    x_of(pts[:, 0], pts[:, 1]), y_of(pts[:, 1], pts[:, 0]), strict=True
                ):
                    fh.write(f"{xv:.4f},{yv:.4f}\n")
    if failures:
        print("\nTICK CHECK FAILED -- the axis calibration does not predict the printed ticks:")
        for line in failures:
            print(f"  {line}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
