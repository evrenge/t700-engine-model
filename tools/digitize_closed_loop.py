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


def separate(sel: np.ndarray, width: float) -> tuple[np.ndarray, ...]:
    """The solid component split into the solid line and the dashes joined to it, in pixels.

    A line of slope `s`, in rows per column, inks `width * sqrt(1 + s**2) + |s|` rows of a
    column, the slope taken column by column so that the steep flanks of a rippling line count
    as the line, not as a join. Columns more than `JOIN_PX` thicker, `JOIN_MIN_COLUMNS` or more
    in a row, hold a joined dash: the dash is a line width in from the end farther from the path
    the clean columns draw through. The solid line is read from its clean columns
    (`CLEAN_PX`) and interpolated across the others between them (`CLEAN_REACH_PX`); failing
    that, it is the line-width slice at the other end, or the column's centre. Returns
    `(columns, solid centre, solid run)` per column and `(columns, dash centre, dash run)` for
    the joined dashes."""
    xs, top, bot = extents(sel)
    mid = 0.5 * (top + bot)
    slope = np.gradient(mid)
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
    clean = bot - top + 1.0 <= run + CLEAN_PX
    cx = xs[clean]
    left = np.clip(np.searchsorted(cx, xs, side="right") - 1, 0, cx.size - 1)
    right = np.clip(np.searchsorted(cx, xs, side="left"), 0, cx.size - 1)
    between = (cx[left] <= xs) & (cx[right] >= xs) & (cx[right] - cx[left] <= CLEAN_REACH_PX)
    bridged = ~clean & between
    centre[bridged] = np.interp(xs[bridged], cx, mid[clean])
    return xs, centre, run, xs[joined], dash, run[joined]


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
    the dashed line hidden under the solid: there it takes the solid line's centre. Returns
    the columns, centres and runs with a `seen` flag, 1 for dash ink and 0 for hidden."""
    o = np.argsort(dx)
    dx, dc, drun = dx[o], dc[o], drun[o]
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


def _line_at(x: float, ox: np.ndarray, oy: np.ndarray) -> float | None:
    """The straight line through the points within `NEAR_PX` columns of `x`, at `x`; `None`
    with fewer than two columns of them."""
    near = np.abs(ox - x) <= NEAR_PX
    if np.unique(ox[near]).size < 2:
        return None
    a, b = np.polyfit(ox[near], oy[near], 1)
    return float(a * x + b)


def check_dashes(lab, dashes, jx, jc, jrun, sx, sc, width: float):
    """The dashes and joined dashes the dashed line's own ink agrees on.

    A dash is kept when the straight line through the other dashes within `NEAR_PX` columns
    passes within `SPECK_LINES` line widths of it; with too few of them, when it lies that close
    to the solid line or to the joined dashes. A joined dash is dropped when the dashes around
    it put the dashed line farther than that from it: a thick stretch of the solid line, not a
    dash. Returns the kept dash labels and the kept joined dashes' columns, centres and runs."""
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
        line = _line_at(cx[k], cx[others], cy[others])
        if line is None:
            line = _line_at(cx[k], jx, jc)
            near = abs(cy[k] - np.interp(cx[k], sx, sc)) <= far
            ok = near or (line is not None and abs(cy[k] - line) <= far)
        else:
            ok = abs(cy[k] - line) <= far
        if ok:
            keep.append(k)
    kx, ky = cx[keep], cy[keep]
    agree = np.array(
        [
            (lv := _line_at(x, kx, ky)) is None or abs(c - lv) <= far
            for x, c in zip(jx, jc, strict=True)
        ],
        bool,
    )
    return [dashes[k] for k in keep], jx[agree], jc[agree], jrun[agree]


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
        for p in find_panels(ink, cfg["ranges"]):
            x_of, v_of, terr = calibrate(ink, p, cfg["t_range"], cfg["t_major"])
            lab, sl = components(ink, p)
            legend = cfg["legend"].get(p.key)
            solid, dashes = classify(lab, sl, v_of, legend)
            width = line_width(lab, solid)
            sx, sc, srun, jx, jc, jrun = separate(np.isin(lab, solid), width)
            dashes, jx, jc, jrun = check_dashes(lab, dashes, jx, jc, jrun, sx, sc, width)
            fx, ftop, fbot = extents(np.isin(lab, dashes))
            frun = fbot - ftop + 1.0
            dx, dc, drun, seen = fill_hidden(
                np.concatenate([fx, jx]),
                np.concatenate([0.5 * (ftop + fbot), jc]),
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
                            "one and takes its value"
                            if flag is not None
                            else ""
                        )
                        + f". The solid line is the components at least {SOLID_MIN_PX} px "
                        f"wide and the pieces continuing them, the dashes the rest, and a dash "
                        f"joined to the solid line is cut from it where its column is over "
                        f"{JOIN_PX} px thicker than the line; components under {SPECK_PX} px, "
                        f"and dashes over {SPECK_LINES:g} line widths from the dashes around them, "
                        f"are specks" + legend_note(legend) + "\n"
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
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
