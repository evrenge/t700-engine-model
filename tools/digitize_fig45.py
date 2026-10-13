"""Digitize Figures 4 and 5 -- the linear models' response to fuel pulses [pdf pp.30, 34].

Each figure is three stacked panels on one time axis, 0 to 16 s:

    1 FUEL  fuel flow, % of trim -- the input, one solid line
    2 NG    gas turbine speed, % of trim
    3 NP    power turbine speed, % of trim

Figure 4 overlays four traces that lie almost on one another: the Gen Hel UH-60A simulation
(solid) and the 5-DOF, 2-DOF and reduced-order linear models (dashed, dash-dot, dotted).
Figure 5 overlays two: the 6-DOF (dashed) and 5-DOF (dash-dot) linear models. The traces
coincide over most of their length and are drawn in broken line styles, so this tool does not
claim to tell them apart by style. Per column it records the ink band: the uppermost and the
lowermost line centres (each edge inset by half the measured line width). Every trace in the
panel lies between them at that time; where only one line is inked, the two are equal. It also
records every ink run, from which `validation/test_fig45.py` separates each figure's two lines.

The input panel's plateaus sit on printed tick values (100, 105, 95, 100, 102.5 % of trim), so
their digitized levels measure the page's own read error, as Figures 9-10's fuel panels do.
Works on the native 300 dpi 1-bit scan, never resampled; the panels and their skew come from
`digitize_fig910`, which owns them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from digitize_fig910 import calibrate, find_panels  # noqa: E402
from figure.raster import native_bitmap  # noqa: E402

REF = Path("data/reference")
PAGES = {4: 30, 5: 34}
RANGES = [
    ("fuel", "FUEL FLOW, % OF TRIM", 92.5, 107.5),
    ("ng", "GAS TURBINE SPEED, %", 99.0, 101.0),
    ("np", "POWER TURBINE SPEED, %", 99.0, 102.0),
]
"""Each panel's printed y range, from its own tick labels [pdf pp.30, 34]; the two figures
share them."""
T_RANGE = (0.0, 16.0)
T_MAJOR = 2.0
"""The time axis: printed 0 to 16.0, labelled every 2.0 s."""
T_END = 15.05
"""The traces stop at 15 s; past it the fuel panel holds legend text, not data."""

MARGIN_PX = 16
"""Rows kept clear of each frame, for its ticks (as `digitize_fig910.trace_line`)."""
SPECK_WINDOW = 5
"""Columns in the rolling median that removes single-column ink specks from each band edge:
a scan defect is one or two columns wide, a trace's steep segment many more."""
STRAY_WINDOW = 15
STRAY_PX = 25
"""The speck filter in `band`: its rolling median spans 15 columns (0.15 s), short enough
that a trace's sharpest corner (Figure 4's NG peak, about 2.4 px per column on the way down)
stays within 25 px of it, long enough that a five-column speck cannot move it."""
STEP_DROP = 0.5
"""A change in the fuel trace of more than this, in % of trim, between adjacent columns is a
step edge."""
INPUT_LEVELS = (105.0, 95.0, 100.0, 102.5)
"""The four plateaus after trim the input panel draws, read off its tick labels [pdf pp.30,
34]: what the digitized levels are checked against, not what they are replaced by. Before
the first step the flow is the trim, 100 % by the axis's own definition; that short segment
runs into the left frame's 100.0 tick and is not read."""
MERGE_LEVEL = 1.0
"""Adjacent segments closer than this, in % of trim, are one plateau split by a speck."""


def _runs(col: np.ndarray) -> list[tuple[int, int]]:
    ys = np.flatnonzero(col)
    if ys.size == 0:
        return []
    groups = np.split(ys, np.flatnonzero(np.diff(ys) > 1) + 1)
    return [(int(g[0]), int(g[-1])) for g in groups]


def _scalar(v) -> float:
    return float(np.ravel(v)[0])


def input_trace(ink: np.ndarray, p, x_of, v_of) -> tuple[np.ndarray, np.ndarray, float]:
    """The fuel trace as `(t, value)` per column, and its line width in pixels.

    The trace is the uppermost ink in every column: each legend sits below it."""
    ts, vs, widths = [], [], []
    for x in range(int(p.left_at_top) + MARGIN_PX, int(p.right_at_top) - MARGIN_PX):
        y0 = int(np.ceil(p.top_at(x))) + MARGIN_PX
        y1 = int(np.floor(p.bot_at(x))) - MARGIN_PX
        runs = _runs(ink[y0:y1, x])
        if not runs:
            continue
        a, b = runs[0]
        t = _scalar(x_of(x, y0 + a))
        if t > T_END:
            continue
        ts.append(t)
        vs.append(_scalar(v_of(y0 + 0.5 * (a + b), x)))
        widths.append(b - a + 1)
    w = np.asarray(widths)
    return np.asarray(ts), np.asarray(vs), float(np.median(w[w <= np.percentile(w, 75)]))


def plateaus(ts: np.ndarray, vs: np.ndarray) -> list[tuple[float, float, float]]:
    """Split the input trace at its step edges: `(start time, level, level scatter)` each,
    after the trim.

    A step's vertical edge is drawn in a single column, so the edge time is the column where
    the value jumps; the level is the median between edges, away from them."""
    vs = ndimage.median_filter(vs, SPECK_WINDOW, mode="nearest")
    jumps = np.flatnonzero(np.abs(np.diff(vs)) > STEP_DROP)
    edges = [ts[0]]
    for j in jumps:
        t = 0.5 * (ts[j] + ts[j + 1])
        if t - edges[-1] > 0.05:
            edges.append(t)
    edges.append(ts[-1])
    out = []
    for a, b in zip(edges, edges[1:], strict=False):
        sel = (ts > a + 0.03) & (ts < b - 0.03)
        if sel.sum() < 5:
            continue
        level = float(np.median(vs[sel]))
        if out and abs(level - out[-1][1]) < MERGE_LEVEL:
            continue  # a speck split one plateau: keep its first start
        out.append((a, level, float(np.std(vs[sel]))))
    if out and abs(out[0][1] - 100.0) < MERGE_LEVEL and out[0][0] < 0.5:
        out = out[1:]  # the trim segment before the first step
    return out


def column_runs(ink: np.ndarray, p) -> list[tuple[int, list[tuple[int, int]]]]:
    """Per column, the ink runs between the frames as `(top row, bottom row)`.

    Ink far from its neighbours' band is a speck, not a trace: a run more than `STRAY_PX`
    outside the rolling median of the band over `STRAY_WINDOW` columns is dropped. Figure 5 has
    one such dot above its NG traces near 8 s, five columns wide, where the traces' own dashes
    leave a gap."""
    cols = []
    for x in range(int(p.left_at_top) + MARGIN_PX, int(p.right_at_top) - MARGIN_PX):
        y0 = int(np.ceil(p.top_at(x))) + MARGIN_PX
        y1 = int(np.floor(p.bot_at(x))) - MARGIN_PX
        runs = [(y0 + a, y0 + b) for a, b in _runs(ink[y0:y1, x])]
        if runs:
            cols.append((x, runs))
    tops = ndimage.median_filter(np.array([r[0][0] for _, r in cols], float), STRAY_WINDOW)
    bots = ndimage.median_filter(np.array([r[-1][1] for _, r in cols], float), STRAY_WINDOW)
    out = []
    for (x, runs), lo_px, hi_px in zip(cols, tops, bots, strict=True):
        runs = [r for r in runs if r[1] >= lo_px - STRAY_PX and r[0] <= hi_px + STRAY_PX]
        if runs:
            out.append((x, runs))
    return out


def run_table(ink: np.ndarray, p, x_of, v_of, width: float) -> np.ndarray:
    """One row per ink run: `(column, t, centre, half span)`, the half span being how far the run
    reaches either side of its centre beyond half a line width (zero for a flat line, large on
    a steep stroke). Values in % of trim."""
    rows = []
    inset = 0.5 * (width - 1.0)
    for x, runs in column_runs(ink, p):
        for a, b in runs:
            mid = 0.5 * (a + b)
            t = _scalar(x_of(x, mid))
            if t > T_END:
                continue
            c = _scalar(v_of(mid, x))
            reach = max(0.5 * (b - a) - inset, 0.0)
            half = abs(_scalar(v_of(mid - reach, x)) - c)
            rows.append((x, t, c, half))
    return np.asarray(rows)


def band(ink: np.ndarray, p, x_of, v_of, width: float) -> np.ndarray:
    """Per column: `(t, upper, lower, runs)`, the uppermost and lowermost line centres."""
    rows = []
    inset = 0.5 * (width - 1.0)
    for x, runs in column_runs(ink, p):
        top = runs[0][0] + inset
        bot = runs[-1][1] - inset
        if bot < top:  # one line thinner than the median: its centre
            top = bot = 0.5 * (runs[0][0] + runs[-1][1])
        t = _scalar(x_of(x, top))
        if t > T_END:
            continue
        rows.append((t, _scalar(v_of(top, x)), _scalar(v_of(bot, x)), len(runs)))
    a = np.asarray(rows)
    # single-column specks: a rolling median on each edge, the run count kept as read
    a[:, 1] = ndimage.median_filter(a[:, 1], SPECK_WINDOW, mode="nearest")
    a[:, 2] = ndimage.median_filter(a[:, 2], SPECK_WINDOW, mode="nearest")
    return a


def main() -> int:
    REF.mkdir(parents=True, exist_ok=True)
    failures = []
    for fig, page in PAGES.items():
        ink = native_bitmap(page)
        panels = find_panels(ink, RANGES)
        cal = {p.key: (p, *calibrate(ink, p, T_RANGE, T_MAJOR)) for p in panels}
        p, x_of, v_of, terr = cal["fuel"]
        ts, vs, width = input_trace(ink, p, x_of, v_of)
        steps = plateaus(ts, vs)
        span = p.v_hi - p.v_lo
        levels = [s[1] for s in steps]
        if len(levels) != len(INPUT_LEVELS):
            failures.append(f"Figure {fig}: {len(levels)} plateaus after trim, the panel draws 4")
            continue
        miss = [100.0 * (got - want) / span for got, want in zip(levels, INPUT_LEVELS, strict=True)]
        worst = max(miss, key=abs)
        print(f"\nFigure {fig} (pdf p.{page}): line width {width:.1f} px, x tick {terr:.4f} s")
        for (t0, lv, _sd), m in zip(steps, miss, strict=True):
            print(f"  plateau from {t0:7.4f} s at {lv:8.4f} % of trim ({m:+.3f} %h)")
        if abs(worst) > 0.6:
            failures.append(f"Figure {fig}: a plateau misses its tick value by {worst:+.3f} %h")
        path = REF / f"fig{fig:02d}_input.csv"
        with path.open("w") as fh:
            fh.write(f"# source: TM-100991 pdf p.{page}, Figure {fig}, panel {p.label}\n")
            fh.write(
                "# quantity: the commanded fuel flow, one plateau per row after the trim "
                "(100 % of trim until the first row's start)\n"
            )
            fh.write("# method: digitized -- native 300 dpi 1-bit scan, tools/digitize_fig45.py\n")
            fh.write(
                "# t_start_s: the column where the drawn step edge sits; level_pct: the median "
                "of the trace between edges, % of trim; scatter_pct: its standard deviation "
                "about that median, % of trim\n"
            )
            fh.write(
                f"# read error: the plateaus sit on printed tick values "
                f"{', '.join(f'{v:g}' for v in INPUT_LEVELS)}; the worst digitized level "
                f"misses its value by {worst:+.3f} % of panel height\n"
            )
            fh.write(f"# points: {len(steps)}\n")
            fh.write(f"# x tick check: worst interior major tick off by {terr:.4f} s\n")
            fh.write(f"# frame tilt corrected: {p.tilt_pct_height:+.2f} % of panel height\n")
            fh.write("# NOT transcribed. Values carry read error.\n")
            fh.write("t_start_s,level_pct,scatter_pct\n")
            for t0, lv, sd in steps:
                fh.write(f"{t0:.5f},{lv:.5f},{sd:.5f}\n")
        for key in ("ng", "np"):
            p, x_of, v_of, terr = cal[key]
            b = band(ink, p, x_of, v_of, width)
            path = REF / f"fig{fig:02d}_{key}.csv"
            with path.open("w") as fh:
                fh.write(f"# source: TM-100991 pdf p.{page}, Figure {fig}, panel {p.label}\n")
                what = (
                    "the band of four nearly coincident traces: Gen Hel UH-60A simulation, "
                    "5-DOF, 2-DOF and reduced-order linear models"
                    if fig == 4
                    else "the band of two traces: 6-DOF (dashed) and 5-DOF (dash-dot) linear models"
                )
                fh.write(f"# quantity: {p.label} vs time, {what}\n")
                fh.write(
                    "# method: digitized -- native 300 dpi 1-bit scan, tools/digitize_fig45.py\n"
                )
                fh.write(
                    f"# columns: t_s; upper and lower, the uppermost and lowermost line centres "
                    f"in % of trim, each ink edge inset by half the {width:.1f} px line width; "
                    f"runs, the ink runs in the column. Ink over {STRAY_PX} px outside its "
                    f"{STRAY_WINDOW}-column neighbourhood is dropped as a speck, and a "
                    f"{SPECK_WINDOW}-column rolling median smooths upper and lower\n"
                )
                fh.write(f"# value: axis range {p.v_lo} to {p.v_hi}\n")
                fh.write(f"# points: {b.shape[0]}\n")
                fh.write(f"# x tick check: worst interior major tick off by {terr:.4f} s\n")
                fh.write(f"# frame tilt corrected: {p.tilt_pct_height:+.2f} % of panel height\n")
                fh.write("# NOT transcribed. Values carry read error.\n")
                fh.write("t_s,upper,lower,runs\n")
                for t, up, lo, n in b:
                    fh.write(f"{t:.5f},{up:.5f},{lo:.5f},{int(n)}\n")
            r = run_table(ink, p, x_of, v_of, width)
            path = REF / f"fig{fig:02d}_{key}_runs.csv"
            what = (
                "Gen Hel UH-60A simulation (solid) and the 5-DOF, 2-DOF and reduced-order "
                "linear models (dashed, dash-dot, dotted)"
                if fig == 4
                else "6-DOF (dashed) and 5-DOF (dash-dot) linear models"
            )
            with path.open("w") as fh:
                fh.write(f"# source: TM-100991 pdf p.{page}, Figure {fig}, panel {p.label}\n")
                fh.write(
                    f"# quantity: {p.label} vs time, every ink run in every column, so the "
                    f"two traces can be told apart where they cross: {what}\n"
                )
                fh.write(
                    "# method: digitized -- native 300 dpi 1-bit scan, tools/digitize_fig45.py\n"
                )
                fh.write(
                    f"# columns: column, the scan's pixel column, shared by the runs one "
                    f"column holds; t_s; centre, the run's centre in % of trim; half_span, how "
                    f"far the run reaches either side of its centre beyond half the "
                    f"{width:.1f} px line width, % of trim. Specks are dropped as for "
                    f"fig{fig:02d}_{key}.csv\n"
                )
                fh.write(f"# value: axis range {p.v_lo} to {p.v_hi}\n")
                fh.write(f"# points: {r.shape[0]}\n")
                fh.write(f"# x tick check: worst interior major tick off by {terr:.4f} s\n")
                fh.write(f"# frame tilt corrected: {p.tilt_pct_height:+.2f} % of panel height\n")
                fh.write("# NOT transcribed. Values carry read error.\n")
                fh.write("column,t_s,centre,half_span\n")
                for x, t, cc, hh in r:
                    fh.write(f"{int(x)},{t:.5f},{cc:.5f},{hh:.5f}\n")
            print(f"  {key}: {b.shape[0]} columns, x tick worst {terr:.4f} s")
    if failures:
        print("\nINPUT CHECK FAILED:", *failures, sep="\n  ")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
