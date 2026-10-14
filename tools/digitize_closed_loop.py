"""Digitize the closed-loop figures -- Gen Hel with the T700 against flight test
[pdf pp.49, 52].

Figures 12 and 14 are five stacked panels each on one time axis:

    1 INPUT      the pilot's input, in.: lateral cyclic (Figure 12), pedals (Figure 14)
    2 TORQUE     total torque of both engines, ft-lb x 10^-3
    3 NG         gas generator turbine speed, %
    4 FUEL       vehicle fuel flow, both engines, lb/hr
    5 NR         main rotor speed, %

Each panel carries two lines: "SIKORSKY GEN HEL WITH T700 AND GEARBOX (AMES VERSION)", solid,
and "1982 USAAEFA UH-60A TEST DATA", dashed [pdf pp.49, 52]. On the scan the solid line is one
long connected stroke and every dash a separate piece about 15 px long, so the two are told
apart by connected components: a component at least `SOLID_MIN_PX` wide is the solid line, a
smaller one a dash. Where the two lines cross or run together, a dash touches the solid line
and joins its component (`separate`): the column's ink is then thicker than the solid line
alone, so the solid line keeps the line-width slice on the side its clean neighbours lead to,
and the rest is the dash. Where the dashed line lies on the solid one, it is hidden for longer
than its own gaps; there it is the solid line's value, flagged `seen` 0 (`fill_hidden`).
Works on the native 300 dpi 1-bit scan, never resampled; the panels and their skew come from
`digitize_fig910`, which owns them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from digitize_fig45 import MARGIN_PX, _runs, _scalar  # noqa: E402
from digitize_fig910 import calibrate, find_panels  # noqa: E402
from figure.raster import native_bitmap  # noqa: E402

REF = Path("data/reference")
FIGURES = {
    11: {
        "page": 48,
        "ranges": [
            ("collective", "COLLECTIVE, in.", 0.0, 10.0),
            ("torque", "TOT. TORQUE, ft-lb X 10^-3", 18.0, 28.0),
            ("ng", "GAS GEN. TURBINE, %", 88.0, 93.0),
            ("fuel", "VEH. FUEL, lb/hr", 500.0, 900.0),
            ("nr", "MAIN ROTOR SPEED, %", 99.0, 102.0),
        ],
        "t_range": (2.0, 8.0),
        "t_major": 1.0,
        "legend": {"collective": (None, 4.0)},
    },
    12: {
        "page": 49,
        "ranges": [
            ("lat_cyclic", "LATERAL CYCLIC, in.", -5.0, 5.0),
            ("torque", "TOT. TORQUE, ft-lb X 10^-3", 18.0, 23.0),
            ("ng", "GAS GEN. TURBINE, %", 88.0, 91.0),
            ("fuel", "VEH. FUEL, lb/hr", 500.0, 700.0),
            ("nr", "MAIN ROTOR SPEED, %", 97.0, 99.0),
        ],
        "t_range": (0.0, 6.0),
        "t_major": 1.0,
        "legend": {"lat_cyclic": (0.5, None)},
    },
    13: {
        "page": 51,
        "ranges": [
            ("collective", "COLLECTIVE, in.", 0.0, 10.0),
            ("torque", "TOT. TORQUE, ft-lb X 10^-3", 20.0, 50.0),
            ("ng", "GAS GEN. TURBINE, %", 92.0, 96.0),
            ("fuel", "VEH. FUEL, lb/hr", 600.0, 1100.0),
            ("nr", "MAIN ROTOR SPEED, %", 98.0, 101.0),
        ],
        "t_range": (1.0, 7.0),
        "t_major": 1.0,
        "legend": {"collective": (None, 4.0)},
    },
    15: {
        "page": 53,
        "ranges": [
            ("collective", "COLLECTIVE, in.", 0.0, 10.0),
            ("accel", "NORMAL ACCELERATION, g", -1.0, 2.0),
            ("mr_torque", "MAIN ROTOR TORQUE, ft-lb X 10^-3", 0.0, 30.0),
            ("nr", "MAIN ROTOR SPEED, %", 96.0, 104.0),
        ],
        "t_range": (0.0, 6.0),
        "t_major": 1.0,
        "legend": {"collective": (6.0, None)},
    },
    14: {
        "page": 52,
        "ranges": [
            ("pedals", "PEDALS, in.", -3.0, 3.0),
            ("torque", "TOT. TORQUE, ft-lb X 10^-3", 20.0, 50.0),
            ("ng", "GAS GEN. TURBINE, %", 92.0, 97.0),
            ("fuel", "VEH. FUEL, lb/hr", 600.0, 1200.0),
            ("nr", "MAIN ROTOR SPEED, %", 98.0, 102.0),
        ],
        "t_range": (2.0, 10.0),
        "t_major": 2.0,
        "legend": {"pedals": (0.5, None)},
    },
}
"""Per figure: its page, each panel's printed y range from its own tick labels, the time axis as
printed (Figure 11 2.0 to 8.0, Figure 12 0.0 to 6.0 and Figure 13 1.0 to 7.0, each labelled
every 1.0; Figure 14 2.0 to 10.0 every 2.0), and the values beyond which a panel's ink is
legend, not data: on Figures 12 and 14 the legend fills the top panel above 1 in. and the
traces stay below 0 in.; on Figures 11 and 13 the legend lies below 3.5 in. and the traces above
4.5 in. [pdf pp.48, 49, 51, 52]."""
SOLID_MIN_PX = 40
"""A component at least this many columns wide is the solid line: the dashes are 14-15 px."""
SPECK_PX = 12
"""A component of fewer pixels is a speck, not a dash: a dash holds 40-60."""
JOIN_PX = 2.5
"""A column of the solid component more than this many pixels thicker than the solid line
alone holds a joined dash. In Figure 14's torque and fuel panels, where the two lines never
touch, no column of the solid line is thicker than that; a touching dash adds a line width,
3-4 px."""
JOIN_MIN_COLUMNS = 3
"""A joined dash thickens at least this many columns in a row: a dash is about 15 px long, and
a single thick column is a speck on the line."""
CLEAN_PX = 1.0
CLEAN_REACH_PX = 30
"""A column of the solid component no more than `CLEAN_PX` thicker than the solid line alone
holds the solid line only. Where the two lines run together a dash overlaps the line by less
than a line width, and the fused band's centre steps by a pixel at every dash: there the solid
line is interpolated between the clean columns on either side, when both lie within
`CLEAN_REACH_PX` columns (a dash and its two edges, with room)."""
END_GAP_PX = 6
END_ROWS = 12
"""A short component starting within `END_GAP_PX` columns beyond an end of the solid line, and
within `END_ROWS` rows of its height there, continues it: Figure 12's NG line breaks 4 px short
of the frame and its last 20 columns stand alone."""
NEAR_PX = 60
SPECK_LINES = 4.0
"""A dash more than `SPECK_LINES` line widths from the straight line through the other dashes
within `NEAR_PX` columns is a speck (Figure 12 has one below its rotor-speed traces near 1.9 s),
and so is a joined dash that far from them (`check_dashes`)."""
SLOPE_SPAN_PX = 3
"""A column's slope, which sets the thickness a clean column may have, is the smaller of its own
step and the step across `SLOPE_SPAN_PX` columns either side: where a dash touches the line,
the fused columns' centres jump, and only across them does the line's own slope show. Taken
from those columns alone, Figure 11's fuel line at 4.02 s, fused with a dash, looked clean."""
BRIDGE_FIT_PX = 6
"""A bridge across fused columns matches the solid line's slope on either side, fitted over
this many clean columns beyond each end."""
ZIGZAG_FRACTION = 0.5
"""A panel whose solid line is steeper than a row per column in more than this fraction of its
columns is read stroke by stroke and at its turns (`zigzag`): Figure 15's main-rotor torque, at
82 %, against at most 17 % in every other panel of Figures 11-15."""
CHAIN_PX = 25
"""A dash whose end lies within this many pixels of a kept dash's end continues the dashed line
and is kept too, however the line bends there: consecutive dashes on Figures 11-15 are 6-18 px
apart end to end, and the specks 64-97 px from any dash (Figure 15's flight-test torque, which
swings in V shapes, lost about 37 dashes to the straight-line test alone)."""
END_DASH_PX = 15
"""At a panel end, the dashed line is followed back to the solid one along its first (or last)
dash alone, the columns within this many of its end: a dash is about 15 px long, and the dashes
after it already curve (Figure 15's rotor speed)."""
HIDDEN_PX = 30
"""A stretch without dash ink longer than this many columns is not a gap in the dash pattern
(about 10 px) but the dashed line lying on the solid one."""


def components(ink: np.ndarray, p) -> tuple[np.ndarray, list]:
    """The panel's ink between its frames, labelled 8-connected, and each label's slices."""
    inside = np.zeros_like(ink, bool)
    for x in range(int(p.left_at_top) + MARGIN_PX, int(p.right_at_top) - MARGIN_PX):
        y0 = int(np.ceil(p.top_at(x))) + MARGIN_PX
        y1 = int(np.floor(p.bot_at(x))) - MARGIN_PX
        inside[y0:y1, x] = True
    lab, _ = ndimage.label(ink & inside, structure=np.ones((3, 3)))
    return lab, ndimage.find_objects(lab)


def extents(sel: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per inked column of `sel`: the column, and the top and bottom rows of its ink."""
    xs = np.flatnonzero(sel.any(axis=0))
    top = np.array([_runs(sel[:, x])[0][0] for x in xs], float)
    bot = np.array([_runs(sel[:, x])[-1][1] for x in xs], float)
    return xs, top, bot


def _slope(xs: np.ndarray, mid: np.ndarray) -> np.ndarray:
    """Rows per column, column by column: the smaller of the column's own step and the step
    across `SLOPE_SPAN_PX` columns either side."""
    own = np.gradient(mid, xs)
    k = SLOPE_SPAN_PX
    across = (np.interp(xs + k, xs, mid) - np.interp(xs - k, xs, mid)) / (2.0 * k)
    return np.where(np.abs(own) < np.abs(across), own, across)


def _end_slope(cx: np.ndarray, cy: np.ndarray, k: int, side: int) -> float:
    """The solid line's slope at clean column `cx[k]`, from the clean columns within
    `BRIDGE_FIT_PX` of it on `side` (-1 left, 1 right), the column itself included."""
    near = (side * (cx - cx[k]) >= 0) & (np.abs(cx - cx[k]) <= BRIDGE_FIT_PX)
    if near.sum() < 2:
        return float("nan")
    return float(np.polyfit(cx[near], cy[near], 1)[0])


def _bridge(x: np.ndarray, cx: np.ndarray, cy: np.ndarray, left: np.ndarray, right: np.ndarray):
    """Values at the columns `x` between clean columns `cx[left]` and `cx[right]`: the cubic
    matching the solid line's value and slope at both, straight where a slope is missing."""
    out = np.empty(x.size)
    for i, (xi, a, b) in enumerate(zip(x, left, right, strict=True)):
        x0, x1, y0, y1 = cx[a], cx[b], cy[a], cy[b]
        h = float(x1 - x0)
        u = (xi - x0) / h
        m0, m1 = _end_slope(cx, cy, a, -1), _end_slope(cx, cy, b, 1)
        if np.isnan(m0) or np.isnan(m1):
            out[i] = y0 + u * (y1 - y0)
            continue
        h00, h10 = 2 * u**3 - 3 * u**2 + 1, u**3 - 2 * u**2 + u
        h01, h11 = -2 * u**3 + 3 * u**2, u**3 - u**2
        out[i] = h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1
    return out


def separate(sel: np.ndarray, width: float) -> tuple[np.ndarray, ...]:
    """The solid component split into the solid line and the dashes joined to it, in pixels.

    A line of slope `s`, in rows per column, inks `width * sqrt(1 + s**2) + |s|` rows of a
    column, the slope taken column by column (`_slope`) so that the steep flanks of a rippling
    line count as the line, not as a join. Columns more than `JOIN_PX` thicker,
    `JOIN_MIN_COLUMNS` or more in a row, hold a joined dash: the dash is a line width in from
    the end farther from the path the clean columns draw through. Where the dashed line crosses
    without touching, a column holds two runs of ink, and the solid line is the one nearer that
    path. The solid line is read from its clean columns (`CLEAN_PX`) and bridged across the
    others between them (`CLEAN_REACH_PX`) by a cubic matching its slope on either side
    (`_bridge`); failing that, it is the line-width slice at the other end, or the column's
    centre. Returns `(columns, solid centre, solid run)` per column and `(columns, dash centre,
    dash run)` for the joined dashes."""
    xs, top, bot = extents(sel)
    mid = 0.5 * (top + bot)
    slope = _slope(xs, mid)
    run = width * np.sqrt(1.0 + slope**2) + np.abs(slope)
    thick, _ = ndimage.label(bot - top + 1.0 > run + JOIN_PX)
    sizes = np.bincount(thick.ravel())
    joined = (thick > 0) & (sizes[thick] >= JOIN_MIN_COLUMNS)
    path = np.interp(xs, xs[~joined], mid[~joined])
    upper = top + 0.5 * (run - 1.0)
    lower = bot - 0.5 * (run - 1.0)
    solid_top = np.abs(upper - path) <= np.abs(lower - path)
    centre = np.where(joined, np.where(solid_top, upper, lower), mid)
    dash = np.where(solid_top, lower, upper)[joined]
    own_c, own_n = mid.copy(), bot - top + 1.0
    for k, x in enumerate(xs):
        runs = _runs(sel[:, x])
        if len(runs) > 1:
            a, b = min(runs, key=lambda r: abs(0.5 * (r[0] + r[1]) - path[k]))
            own_c[k], own_n[k] = 0.5 * (a + b), b - a + 1.0
    clean = own_n <= run + CLEAN_PX
    centre[clean] = own_c[clean]
    cx, cy = xs[clean], own_c[clean]
    left = np.clip(np.searchsorted(cx, xs, side="right") - 1, 0, cx.size - 1)
    right = np.clip(np.searchsorted(cx, xs, side="left"), 0, cx.size - 1)
    between = (cx[left] <= xs) & (cx[right] >= xs) & (cx[right] - cx[left] <= CLEAN_REACH_PX)
    bridged = ~clean & between
    centre[bridged] = _bridge(xs[bridged], cx, cy, left[bridged], right[bridged])
    return xs, centre, run, xs[joined], dash, run[joined]


def zigzag(sel: np.ndarray, width: float) -> tuple[np.ndarray, ...]:
    """A solid line steeper than a row per column, read column by column and at its turns.

    Along a steep stroke a column's ink is the stroke's passage through it, so the line is the
    centre of that ink; where two strokes cross and a column holds two runs of ink, the line is
    the run its neighbours lead to. At a turn the column holds the end of one stroke and the
    start of the next, and its centre lies below a peak or above a trough: read that way,
    Figure 15's main-rotor torque ripple lost its printed amplitude. So at each peak, and each
    trough, standing out by at least a line width, the line is a half line width in from the
    ink's top, or bottom. No dash is cut from the solid line: where the dashed line lies on it
    (Figure 15 before 1.3 s) the turns are those of both lines together. Returns what
    `separate` does, with no joined dashes."""
    from scipy.signal import find_peaks

    xs, top, bot = extents(sel)
    inset = 0.5 * (width - 1.0)
    centre = 0.5 * (top + bot)
    runs = [_runs(sel[:, x]) for x in xs]
    single = np.array([len(r) == 1 for r in runs])
    path = np.interp(xs, xs[single], centre[single])
    for k in np.flatnonzero(~single):
        a, b = min(runs[k], key=lambda r: abs(0.5 * (r[0] + r[1]) - path[k]))
        centre[k] = 0.5 * (a + b)
    peaks, _ = find_peaks(-top, prominence=width)
    troughs, _ = find_peaks(bot, prominence=width)
    centre[peaks] = top[peaks] + inset
    centre[troughs] = bot[troughs] - inset
    empty = np.zeros(0)
    return xs, centre, bot - top + 1.0, empty, empty, empty


def steep_fraction(sel: np.ndarray) -> float:
    """The fraction of the solid line's columns where it is steeper than a row per column."""
    xs, top, bot = extents(sel)
    return float(np.mean(np.abs(np.gradient(0.5 * (top + bot), xs)) > 1.0))


def values(xs, centre, run, x_of, v_of, width: float) -> np.ndarray:
    """Pixel rows to `(t, value, half span)`, the half span being how far the run reaches
    either side of its centre beyond half a line width (zero for a flat line, large on a
    steep one)."""
    inset = 0.5 * (width - 1.0)
    rows = []
    for x, c, n in zip(xs, centre, run, strict=True):
        v = _scalar(v_of(c, x))
        reach = max(0.5 * (n - 1.0) - inset, 0.0)
        rows.append((_scalar(x_of(x, c)), v, abs(_scalar(v_of(c - reach, x)) - v)))
    return np.asarray(rows)


def fill_hidden(dx, dc, drun, sx, sc, srun, width: float):
    """The dashed line through the stretches where it lies on the solid one.

    A stretch of columns without dash ink longer than `HIDDEN_PX` whose dash ink on both
    sides (or on its one side, at a panel end) sits within a line width of the solid line is
    the dashed line hidden under the solid: there it takes the solid line's centre. At a panel
    end the dashed line may leave the solid one at a slant before its first dash shows
    (Figure 15's torque and rotor speed, Figure 14's rotor speed): there its first dash,
    extended straight back to where it meets the solid line (`END_DASH_PX`), gives the column
    where it leaves, and it is hidden from the panel end to there. Returns the columns,
    centres and runs with a `seen` flag, 1 for dash ink and 0 for hidden."""
    o = np.argsort(dx)
    dx, dc, drun = dx[o].astype(int), dc[o], drun[o]
    solid_at = dict(zip(sx.tolist(), zip(sc.tolist(), srun.tolist(), strict=True), strict=True))

    def on_solid(x: int, c: float) -> bool:
        return x in solid_at and abs(c - solid_at[x][0]) <= width

    edges = [(int(sx[0]) - 1, None)] + list(zip(dx.tolist(), dc.tolist(), strict=True))
    edges += [(int(sx[-1]) + 1, None)]
    add = []
    for (x0, c0), (x1, c1) in zip(edges, edges[1:], strict=False):
        if x1 - x0 <= HIDDEN_PX:
            continue
        ends = [on_solid(x, c) for x, c in ((x0, c0), (x1, c1)) if c is not None]
        if ends and all(ends):
            add += [x for x in range(x0 + 1, x1) if x in solid_at]
        elif c0 is None or c1 is None:
            # A panel end: follow the dashed line back from its first (or last) dash.
            x_dash, step = (int(x1), -1) if c0 is None else (int(x0), 1)
            mine = np.abs(dx - x_dash) <= END_DASH_PX
            fit = None
            if np.unique(dx[mine]).size >= 2:
                fit = tuple(float(v) for v in np.polyfit(dx[mine], dc[mine], 1))
            leave = None
            if fit is not None:
                for x in range(x_dash + step, int(x0) if c0 is None else int(x1), step):
                    if x in solid_at and abs(fit[0] * x + fit[1] - solid_at[x][0]) <= width:
                        leave = x
                        break
            if leave is not None:
                lo, hi = (int(x0) + 1, leave + 1) if c0 is None else (leave, int(x1))
                add += [x for x in range(lo, hi) if x in solid_at]
    hx = np.array(add, int)
    out_x = np.concatenate([dx, hx])
    out_c = np.concatenate([dc, [solid_at[x][0] for x in add]])
    out_r = np.concatenate([drun, [solid_at[x][1] for x in add]])
    seen = np.concatenate([np.ones(dx.size), np.zeros(hx.size)])
    o = np.argsort(out_x, kind="stable")
    return out_x[o], out_c[o], out_r[o], seen[o]


def classify(lab: np.ndarray, sl: list, v_of, legend: tuple | None) -> tuple[list, list]:
    """The labels of the solid line's pieces and of the dashes, legend and specks left out.

    A component at least `SOLID_MIN_PX` wide is solid, and so is a short one that continues it
    past an end (`END_GAP_PX`, `END_ROWS`)."""
    keep = {}
    for i, s in enumerate(sl, start=1):
        if s is None:
            continue
        ys, xs = np.nonzero(lab[s] == i)
        if ys.size < SPECK_PX:
            continue
        if legend is not None:
            v = _scalar(v_of(s[0].start + ys.mean(), s[1].start + xs.mean()))
            lo, hi = legend
            if (lo is None or v > lo) and (hi is None or v < hi):
                continue
        keep[i] = s
    solid = [i for i, s in keep.items() if s[1].stop - s[1].start >= SOLID_MIN_PX]
    grown = True
    while grown:
        grown = False
        xs, top, bot = extents(np.isin(lab, solid))
        ends = ((xs[0], 0.5 * (top[0] + bot[0]), -1), (xs[-1], 0.5 * (top[-1] + bot[-1]), 1))
        for i, s in keep.items():
            if i in solid:
                continue
            col = s[1].start if s[1].start > xs[-1] else s[1].stop - 1
            for x_end, y_end, side in ends:
                rows = np.flatnonzero(lab[:, col] == i)
                if 0 < side * (col - x_end) <= END_GAP_PX and abs(rows.mean() - y_end) <= END_ROWS:
                    solid.append(i)
                    grown = True
                    break
    return solid, [i for i in keep if i not in solid]


def legend_note(legend: tuple | None) -> str:
    """The CSV header's words for a panel's legend band."""
    if legend is None:
        return ""
    lo, hi = legend
    if hi is None:
        return f"; ink above {lo:g} is the legend"
    if lo is None:
        return f"; ink below {hi:g} is the legend"
    return f"; ink from {lo:g} to {hi:g} is the legend"


def _line_near(x: float, ox: np.ndarray, oy: np.ndarray) -> tuple[float, float] | None:
    """Slope and intercept of the straight line through the points within `NEAR_PX` columns of
    `x`; `None` with fewer than two columns of them."""
    near = np.abs(ox - x) <= NEAR_PX
    if np.unique(ox[near]).size < 2:
        return None
    a, b = np.polyfit(ox[near], oy[near], 1)
    return float(a), float(b)


def _off(x: float, y: float, fit: tuple[float, float]) -> float:
    """How far the point lies from the line `fit`, at right angles to it: on a steep stretch
    the gap in rows alone overstates it (Figure 15's torque, where the dashed line falls at
    three rows a column)."""
    a, b = fit
    return abs(y - (a * x + b)) / float(np.hypot(1.0, a))


def check_dashes(lab, dashes, jx, jc, jrun, sx, sc, width: float):
    """The dashes and joined dashes the dashed line's own ink agrees on.

    A dash is kept when the straight line through the other dashes within `NEAR_PX` columns
    passes within `SPECK_LINES` line widths of it, measured at right angles to that line
    (`_off`); with too few of them, when it lies that close to the solid line or to the joined
    dashes. A dash failing that is kept still when one of its ends lies within `CHAIN_PX` of a
    kept dash's end, so a bend in the dashed line keeps its dashes. A joined dash is dropped
    when the dashes around it put the dashed line farther than that from it: a thick stretch
    of the solid line, not a dash. Returns the kept dash labels and the kept joined dashes'
    columns, centres and runs."""
    far = SPECK_LINES * width
    cents = []
    for i in dashes:
        ys, xs = np.nonzero(lab == i)
        cents.append((float(xs.mean()), float(ys.mean())))
    cx = np.array([c[0] for c in cents])
    cy = np.array([c[1] for c in cents])
    keep = []
    for k in range(len(dashes)):
        others = np.arange(cx.size) != k
        line = _line_near(cx[k], cx[others], cy[others])
        if line is None:
            line = _line_near(cx[k], jx, jc)
            near = abs(cy[k] - np.interp(cx[k], sx, sc)) <= far
            ok = near or (line is not None and _off(cx[k], cy[k], line) <= far)
        else:
            ok = _off(cx[k], cy[k], line) <= far
        if ok:
            keep.append(k)
    ends = []
    for i in dashes:
        ys, xs = np.nonzero(lab == i)
        lo, hi = xs.min(), xs.max()
        ends.append(np.array([[lo, ys[xs == lo].mean()], [hi, ys[xs == hi].mean()]]))
    grown = True
    while grown:
        grown = False
        for k in range(len(dashes)):
            if k in keep:
                continue
            gaps = [np.hypot(*(a - b)) for j in keep for a in ends[k] for b in ends[j]]
            if gaps and min(gaps) <= CHAIN_PX:
                keep.append(k)
                grown = True
    keep.sort()
    kx, ky = cx[keep], cy[keep]
    agree = np.array(
        [
            (fit := _line_near(x, kx, ky)) is None or _off(x, c, fit) <= far
            for x, c in zip(jx, jc, strict=True)
        ],
        bool,
    )
    return [dashes[k] for k in keep], jx[agree], jc[agree], jrun[agree]


def dash_strokes(lab: np.ndarray, dashes: list[int]) -> tuple[np.ndarray, ...]:
    """Each dash read as the short straight stroke it is: per column of a dash, its ink's centre
    and run, the centres replaced by the least-squares line through that dash's own columns.
    A column's centre alone steps by whole pixels and bends at the dash's rounded ends; nothing
    is averaged across dashes or gaps. Returns `(columns, centres, runs)`."""
    out_x, out_c, out_r = [], [], []
    for i in dashes:
        xs, top, bot = extents(lab == i)
        mid = 0.5 * (top + bot)
        if xs.size >= 2:
            a, b = np.polyfit(xs, mid, 1)
            mid = a * xs + b
        out_x.append(xs)
        out_c.append(mid)
        out_r.append(bot - top + 1.0)
    if not out_x:
        return np.zeros(0, int), np.zeros(0), np.zeros(0)
    return np.concatenate(out_x), np.concatenate(out_c), np.concatenate(out_r)


def joined_strokes(jx: np.ndarray, jc: np.ndarray) -> np.ndarray:
    """The joined dashes' centres, each run of consecutive columns (a dash cut from the solid
    line) replaced by its least-squares line, as `dash_strokes` reads a free one."""
    out = jc.astype(float).copy()
    if jx.size == 0:
        return out
    order = np.argsort(jx)
    breaks = np.flatnonzero(np.diff(jx[order]) > 1) + 1
    for run in np.split(order, breaks):
        if run.size >= 2:
            a, b = np.polyfit(jx[run], jc[run], 1)
            out[run] = a * jx[run] + b
    return out


def line_width(lab: np.ndarray, ids: list[int]) -> float:
    """The solid line's width in pixels: the median run over its flatter columns."""
    sel = np.isin(lab, ids)
    w = np.array([b - a + 1 for x in np.flatnonzero(sel.any(axis=0)) for a, b in _runs(sel[:, x])])
    return float(np.median(w[w <= np.percentile(w, 50)]))


def main() -> int:
    REF.mkdir(parents=True, exist_ok=True)
    for fig, cfg in FIGURES.items():
        page = cfg["page"]
        ink = native_bitmap(page)
        print(f"\nFigure {fig} (pdf p.{page})")
        panels = []
        for p in find_panels(ink, cfg["ranges"]):
            x_of, v_of, terr = calibrate(ink, p, cfg["t_range"], cfg["t_major"])
            lab, sl = components(ink, p)
            legend = cfg["legend"].get(p.key)
            solid, dashes = classify(lab, sl, v_of, legend)
            steep = steep_fraction(np.isin(lab, solid)) > ZIGZAG_FRACTION
            panels.append((p, x_of, v_of, terr, lab, legend, solid, dashes, steep))
        # A zigzag's columns are all stroke, so its own line width cannot be read: it is the
        # pen's, the median over the figure's other panels.
        pen = float(np.median([line_width(q[4], q[6]) for q in panels if not q[8]] or [np.nan]))
        for p, x_of, v_of, terr, lab, legend, solid, dashes, steep in panels:
            width = pen if steep else line_width(lab, solid)
            read = zigzag if steep else separate
            sx, sc, srun, jx, jc, jrun = read(np.isin(lab, solid), width)
            dashes, jx, jc, jrun = check_dashes(lab, dashes, jx, jc, jrun, sx, sc, width)
            fx, fc, frun = dash_strokes(lab, dashes)
            dx, dc, drun, seen = fill_hidden(
                np.concatenate([fx, jx]),
                np.concatenate([fc, joined_strokes(jx, jc)]),
                np.concatenate([frun, jrun]),
                sx,
                sc,
                srun,
                width,
            )
            lines = {
                "gen_hel": (values(sx, sc, srun, x_of, v_of, width), None),
                "test": (values(dx, dc, drun, x_of, v_of, width), seen),
            }
            for name, what in (
                ("gen_hel", "SIKORSKY GEN HEL WITH T700 AND GEARBOX (AMES VERSION), solid"),
                ("test", "1982 USAAEFA UH-60A TEST DATA, dashed"),
            ):
                r, flag = lines[name]
                path = REF / f"fig{fig:02d}_{p.key}_{name}.csv"
                with path.open("w") as fh:
                    fh.write(f"# source: TM-100991 pdf p.{page}, Figure {fig}, panel {p.label}\n")
                    fh.write(f"# quantity: {p.label} vs time, {what}\n")
                    fh.write(
                        "# method: digitized -- native 300 dpi 1-bit scan, "
                        "tools/digitize_closed_loop.py\n"
                    )
                    fh.write(
                        f"# columns: t_s; value, the centre of the line's ink in the column; "
                        f"half_span, how far that ink reaches either side of its centre beyond "
                        f"half the {width:.1f} px line width"
                        + (
                            "; seen, 1 for dash ink, 0 where the dashed line lies on the solid "
                            "one and takes its value. Each dash's centres lie on the straight "
                            "line through its own columns"
                            if flag is not None
                            else ""
                        )
                        + f". The solid line is the components at least {SOLID_MIN_PX} px "
                        f"wide and the pieces continuing them, the dashes the rest, and a dash "
                        f"joined to the solid line is cut from it where its column is over "
                        f"{JOIN_PX} px thicker than the line; components under {SPECK_PX} px, "
                        f"and dashes over {SPECK_LINES:g} line widths from the dashes around them "
                        f"and over {CHAIN_PX} px from any kept dash's end, are specks"
                        + legend_note(legend)
                        + (
                            f". The solid line here is steeper than a row per column in over "
                            f"{ZIGZAG_FRACTION:.0%} of its columns: it is the centre of each "
                            f"column's ink along a stroke, and a half line width in from the top "
                            f"of each peak and the bottom of each trough; before the dashed line "
                            f"leaves it, the peaks and troughs are those of both lines together"
                            if steep and name == "gen_hel"
                            else ""
                        )
                        + "\n"
                    )
                    fh.write(f"# value: axis range {p.v_lo} to {p.v_hi}\n")
                    fh.write(f"# points: {r.shape[0]}\n")
                    fh.write(f"# x tick check: worst interior major tick off by {terr:.4f} s\n")
                    fh.write(
                        f"# frame tilt corrected: {p.tilt_pct_height:+.2f} % of panel height\n"
                    )
                    fh.write("# NOT transcribed. Values carry read error.\n")
                    if flag is None:
                        fh.write("t_s,value,half_span\n")
                        for t, c, h in r:
                            fh.write(f"{t:.5f},{c:.5f},{h:.5f}\n")
                    else:
                        fh.write("t_s,value,half_span,seen\n")
                        for (t, c, h), f in zip(r, flag, strict=True):
                            fh.write(f"{t:.5f},{c:.5f},{h:.5f},{int(f)}\n")
            print(
                f"  {p.key}: solid {len(solid)} component(s), {len(dashes)} dashes, "
                f"{jx.size} columns with a dash joined, {int((seen == 0).sum())} hidden, line "
                f"{width:.1f} px, x tick worst {terr:.4f} s"
                + (", read as a zigzag" if steep else "")
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
