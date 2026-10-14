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
Glyphs are separated by shape, and taken apart where they overprint (`classify`).
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure.axis import blend_value, choose_anchored, fit_anchored  # noqa: E402
from figure.frame import find_frame_box, fit_edges, make_sv  # noqa: E402
from figure.raster import native_bitmap  # noqa: E402
from figure.ticks import find_ticks  # noqa: E402

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


FRAME_FIT_HALF_BAND_PX = 30
"""Half-width of the window `frame.fit_edges` traces each frame line in, about its nominal
position from `frame.find_frame_box`, in pixels.

Wide enough to contain the whole lean (15-19 px end to end on pp.41-42) plus the line's own
thickness; `fit_edges` keeps only the single narrow run per position and rejects samples
more than 2 px off its straight fit, so a tick or a glyph inside the window is dropped."""


def fit_frames(ink: np.ndarray, frame: tuple[int, int, int, int]) -> dict:
    """Fit all four frame lines, because the plot is a parallelogram and not a rectangle.

    These pages are sheared in x, not y (unlike Figures 9-10, open question #63): the
    vertical frames lean up to ~19 px across the plot height while top and bottom are level.
    Uncorrected that is up to 8 lbm/hr on the fuel-flow axis, worst at the low-power end
    where open question #46 lives. `frame.find_frame_box` returns one column or row per
    frame line, a single point of a leaning one; fitting each line along its length removes
    that.

    Returns slope/intercept for each: horizontals as `y = a*x + b`, verticals as
    `x = a*y + b`.
    """
    left, top, right, bottom = frame
    h = FRAME_FIT_HALF_BAND_PX
    fc, _ = fit_edges(
        ink,
        [
            ("top", 1, top - h, top + h + 1, left + 5, right - 4),
            ("bottom", 1, bottom - h, bottom + h + 1, left + 5, right - 4),
            ("left", 0, left - h, left + h + 1, top + 5, bottom - 4),
            ("right", 0, right - h, right + h + 1, top + 5, bottom - 4),
        ],
        verbose=False,
    )
    return {nm: (float(c[0]), float(c[1])) for nm, c in fc.items()}


TICK_EXCESS_PX = 7
"""How far a tick must lengthen the ink run through its frame, beyond the frame's own
rolling-median width, to be read as a labelled tick by `ticks.find_ticks`, in native-scan
pixels. Measured 2026-10-13 on all twelve edges of Figures 6-8: the long ticks (majors, and
Figure 6's 2.5 %NG half-majors) lengthen it by 8-12 px, the short minors by 2-6."""

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


# ---------------------------------------------------------------- the printed scale bows

BOW_ORDERS = (1, 2, 3, 4)
"""Orders of `axis.fit_anchored` tried per edge, chosen by leave-one-out and never by eye:
order 1 is the straight map between the two printed frame values, and 2-4 add a bow
`t(1-t)P(t)` with `P` of degree 0-2. The `t(1-t)` factor vanishes at both frames, so the two
printed frame values stay exact by construction however the ticks come out; the ticks only
say how the map bends in between, which is the one thing two endpoints cannot know."""


DISC_CORE_PX = 9.5
"""How deep inside the ink a pixel must lie to belong to a filled disc, in native-scan
pixels: the distance to the nearest background pixel. On pp.40-42 every disc on a curve
reaches 10.6 or more, the legend's sample discs 9.4-10.0. No stroke of a cross, a triangle
or a letter passes 6, and the deepest other glyph, a heavily inked cross printed over a
triangle on p.40, reaches 8.5 once its gaps are filled as pinholes."""

PINHOLE_PX = 20
"""A white speck up to this area with solid ink all round it (`PINHOLE_RING_PX`) is scan noise
inside a disc and is filled before discs are looked for: p.42 has one at a disc's centre. A
triangle's open centre stays open, and so do the pieces a cross printed over a triangle cuts
it into (p.40): the ink round those is thin strokes and paper."""

PINHOLE_RING_PX = 4
"""How far round a speck the ink must be solid for it to count as a pinhole."""

DISC_ALONE_MAX_PX = 500
"""A blob holding one disc and no more ink than this is the disc alone, read by its
centroid. Lone discs cover 428-493 px on pp.40-42; a disc touching another glyph, 520 or
more."""

LETTER_HEIGHT_PX = (26, 32)
"""Height of the legend's capital letters on pp.40-42, in pixels (28-30 measured)."""

BASELINE_TOL_PX = 3
"""How far the top and bottom rows of letters on one line may differ, in pixels."""

TRIANGLE_HOLE_PX = 30
"""Paper a lone glyph must enclose to be a triangle, in pixels: on pp.40-42 every triangle
encloses 67-97, every cross none. The bottom-row test this replaces read four triangles
whose base prints thin (pp.40, 42) as crosses."""

OUTLINE_ALONE_MAX_PX = 340
"""The same for a blob without a disc: a lone cross or triangle covers 200-330 px."""

RAY_STEP_DEG = 3.0
"""Angular spacing of the rays along which a disc's edge is found (`_disc_fit`)."""

DISC_RUNAWAY_PX = 6.0
"""How far a merged disc's fitted centre may move from its core before the fit is refused
(recipe, trap 9). Measured moves on pp.40-42 are under 3 px."""

DISC_CHECK_PX = 0.5
"""How far the disc fit may land from a lone disc's centroid, where both read the same disc;
measured 0.04-0.11 px on pp.40-42."""

EDGE_BAND_PX = 2.5
"""A ray's edge is the disc's own when it lies this close to the disc's radius. A ray that
runs on into a touching glyph or an overlapping disc is dropped."""

OUTLINE_MATCH = 0.70
"""Least fraction of a cross or triangle template's visible ink that must land on ink for a
glyph to be read in a merged blob (`_outlines`). Ink under a disc counts neither way. A cross
half behind a disc, two arms showing, scores 0.71-0.95 on pp.40-41."""

OUTLINE_VISIBLE = 0.25
"""Least fraction of a template that must lie outside the discs: two arms of a cross are
26-32 % of it."""

OUTLINE_NEW = 0.15
"""Least fraction of a template that must land on ink no glyph has explained yet, so a glyph
already read is not read again."""


def _fill_pinholes(sub: np.ndarray) -> np.ndarray:
    """`sub` with its pinholes filled: specks of up to `PINHOLE_PX` inside solid ink."""
    holes = ndimage.binary_fill_holes(sub) & ~sub
    lab, n = ndimage.label(holes)
    out = sub.copy()
    for k, sl in enumerate(ndimage.find_objects(lab), start=1):
        hole = lab[sl] == k
        g = PINHOLE_RING_PX
        if hole.sum() > PINHOLE_PX or min(sl[0].start, sl[1].start) < g:
            continue
        win = (slice(sl[0].start - g, sl[0].stop + g), slice(sl[1].start - g, sl[1].stop + g))
        h = np.pad(hole, g)
        ring = ndimage.binary_dilation(h, iterations=g) & ~h
        if sub[win][ring].all():
            out[win] |= h
    return out


def _disc_cores(depth: np.ndarray) -> list[tuple[float, float]]:
    """One point per filled disc in a blob's depth map: the centroid of each plateau of
    local maxima at least `DISC_CORE_PX` deep."""
    peak = (depth == ndimage.maximum_filter(depth, size=9)) & (depth >= DISC_CORE_PX)
    lab, n = ndimage.label(peak, structure=np.ones((3, 3)))
    pts = ndimage.center_of_mass(peak, lab, range(1, n + 1)) if n else []
    keep: list[tuple[float, float]] = []
    for p in pts:  # one disc's plateau can break into pieces a few pixels apart
        if all(np.hypot(p[0] - q[0], p[1] - q[1]) > 8.0 for q in keep):
            keep.append(p)
    return keep


def _disc_fit(filled: np.ndarray, cy: float, cx: float, r0: float):
    """Centre and radius of one filled disc, from its own free edge.

    Rays every `RAY_STEP_DEG` run out from the current centre to the first background
    pixel. An edge within `EDGE_BAND_PX` of the disc's radius `r0` is the disc's own; where a
    cross, a triangle or an overlapping disc touches it the ray runs on in ink and is
    dropped, so a touching glyph does not pull the centre. A least-squares circle goes
    through the kept edges; repeated from the new centre until it stops moving.
    """
    ang = np.deg2rad(np.arange(0.0, 360.0, RAY_STEP_DEG))
    steps = np.arange(0.25, r0 + 3.0 * EDGE_BAND_PX, 0.25)
    r, n = r0, 0
    for _ in range(8):
        iy = np.rint(cy + np.outer(np.sin(ang), steps)).astype(int)
        ix = np.rint(cx + np.outer(np.cos(ang), steps)).astype(int)
        inside = (iy >= 0) & (iy < filled.shape[0]) & (ix >= 0) & (ix < filled.shape[1])
        on = np.zeros(iy.shape, dtype=bool)
        on[inside] = filled[iy[inside], ix[inside]]
        edge = steps[np.argmax(~on, axis=1)]
        keep = ~on.all(axis=1) & (np.abs(edge - r0) <= EDGE_BAND_PX)
        n = int(keep.sum())
        if n < 20:
            break
        px = cx + edge[keep] * np.cos(ang[keep])
        py = cy + edge[keep] * np.sin(ang[keep])
        a = np.column_stack([px, py, np.ones(n)])
        d, e, f = np.linalg.lstsq(a, -(px**2 + py**2), rcond=None)[0]
        ncx, ncy = -d / 2.0, -e / 2.0
        r = float(np.sqrt(ncx**2 + ncy**2 - f))
        moved = float(np.hypot(ncx - cx, ncy - cy))
        cx, cy = float(ncx), float(ncy)
        if moved < 1e-3:
            break
    return cy, cx, r, n


def _templates(sub: np.ndarray, lone: dict, half: int = 20) -> dict:
    """Each GE glyph as this page prints it: the mean of its lone instances, centred on
    their ink centroids, kept where at least half of them carry ink. Returns, per series,
    the template and its own centroid's offset from the template centre."""
    out = {}
    for key in ("ge_unbalanced", "ge_status81"):
        patches = []
        for c in lone[key]:
            y, x = int(round(c["sy"])), int(round(c["sx"]))
            p = sub[y - half : y + half + 1, x - half : x + half + 1]
            if p.shape == (2 * half + 1, 2 * half + 1):
                patches.append(p.astype(float))
        t = np.mean(patches, axis=0) >= 0.5
        ty, tx = ndimage.center_of_mass(t)
        out[key] = (t.astype(float), ty - half, tx - half)
    return out


def _parabola(a: float, b: float, c: float) -> float:
    """Sub-pixel offset of the peak through three equally spaced samples, `b` the middle."""
    den = a - 2.0 * b + c
    return 0.0 if den == 0.0 else float(np.clip(0.5 * (a - c) / den, -0.5, 0.5))


def _outlines(resid: np.ndarray, hidden: np.ndarray, tpl: dict):
    """Crosses and triangles in a merged blob's ink once its discs are taken out.

    Greedy: the template placement whose visible ink lands best on ink, among those that
    explain ink no glyph has explained yet, is read if at least `OUTLINE_MATCH` of it does;
    its footprint is then explained, and the next is sought. Ink under a disc (`hidden`)
    counts neither for nor against a placement. Returns the glyphs, as (series, row,
    column) in the window, and the ink left unexplained.
    """
    found: list[tuple[str, float, float]] = []
    left = resid.copy()
    visible = (~hidden).astype(float)
    for _ in range(8):
        if left.sum() < 25:
            break
        best = None
        for key, (t, oy, ox) in tpl.items():
            tot = t.sum()
            on = ndimage.correlate(resid.astype(float), t, mode="constant")
            vis = ndimage.correlate(visible, t, mode="constant")
            new = ndimage.correlate(left.astype(float), t, mode="constant")
            frac = np.where(vis >= OUTLINE_VISIBLE * tot, on / np.maximum(vis, 1.0), 0.0)
            frac[new < OUTLINE_NEW * tot] = 0.0
            k = np.unravel_index(int(np.argmax(frac)), frac.shape)
            if best is None or frac[k] > best[0]:
                best = (float(frac[k]), key, k, frac, oy, ox, t)
        score, key, (py, px), frac, oy, ox, t = best
        if score < OUTLINE_MATCH:
            break
        dy = _parabola(*frac[py - 1 : py + 2, px]) if 0 < py < frac.shape[0] - 1 else 0.0
        dx = _parabola(*frac[py, px - 1 : px + 2]) if 0 < px < frac.shape[1] - 1 else 0.0
        found.append((key, py + dy + oy, px + dx + ox))
        foot = np.zeros(left.shape, dtype=bool)
        h = t.shape[0] // 2
        y0, x0 = py - h, px - h
        ys, xs = np.nonzero(t)
        ys, xs = ys + y0, xs + x0
        ok = (ys >= 0) & (ys < foot.shape[0]) & (xs >= 0) & (xs < foot.shape[1])
        foot[ys[ok], xs[ok]] = True
        left &= ~ndimage.binary_dilation(foot, iterations=2)
    return found, int(left.sum())


def classify(ink: np.ndarray, frame: tuple[int, int, int, int], margin: int = 22):
    """Find every glyph inside the frame and sort it by shape.

    A glyph standing alone is one blob. A filled disc inks about 78 % of its bounding box, an
    outline triangle or a cross about a third. Triangle against cross: a triangle encloses
    paper (`TRIANGLE_HOLE_PX`), a cross encloses none.

    Where glyphs touch, one blob holds several, and neither its centroid nor its shape is any
    one glyph's. Those blobs are taken apart (`_disc_fit`, `_outlines`): the discs are found
    by their solid cores and centred on their own free edges, then the crosses and triangles
    left beside them by templates of their lone instances on the same page.
    """
    left, top, right, bottom = frame
    sub = ink[top + margin : bottom - margin, left + margin : right - margin].astype(bool)
    filled = _fill_pinholes(sub)
    depth = ndimage.distance_transform_edt(filled)
    lab, n = ndimage.label(filled, structure=np.ones((3, 3)))

    # First pass: every blob with its shape statistics and the discs it holds.
    blobs = []
    for i, sl in enumerate(ndimage.find_objects(lab), start=1):
        blob = lab[sl] == i
        hgt, wid = blob.shape
        area = int((blob & sub[sl]).sum())
        if area < 20 or hgt < 6 or wid < 6:
            continue
        cy, cx = ndimage.center_of_mass(blob & sub[sl])
        blobs.append(
            {
                "i": i,
                "sl": sl,
                "x": cx + sl[1].start + left + margin,
                "y": cy + sl[0].start + top + margin,
                "sx": cx + sl[1].start,
                "sy": cy + sl[0].start,
                "area": area,
                "fill": area / (hgt * wid),
                "hole": int((ndimage.binary_fill_holes(blob & sub[sl]) & ~sub[sl]).sum()),
                "glyph_sized": hgt <= 40 and wid <= 40 and max(hgt, wid) / min(hgt, wid) <= 1.8,
                "cores": _disc_cores(depth[sl] * blob),
            }
        )

    # Reject the in-plot legend without hard-coding where it is: the letters of a line share
    # its top and bottom rows, so each has several neighbours printed between the same two
    # rows. Joined letters count too, however wide; a plotted glyph has no such row of
    # neighbours. The legend's three sample glyphs are left to `_backbone`.
    rows = np.array([(b["sl"][0].start, b["sl"][0].stop, b["x"]) for b in blobs], dtype=float)
    text = set()
    for b, (t0, t1, x) in zip(blobs, rows, strict=True):
        if not LETTER_HEIGHT_PX[0] <= t1 - t0 <= LETTER_HEIGHT_PX[1] or b["cores"]:
            continue
        same = (
            (np.abs(rows[:, 0] - t0) <= BASELINE_TOL_PX)
            & (np.abs(rows[:, 1] - t1) <= BASELINE_TOL_PX)
            & (np.abs(rows[:, 2] - x) < 160.0)
        )
        if np.count_nonzero(same) - 1 >= 3:
            text.add(b["i"])

    def merged(b: dict) -> bool:
        k = len(b["cores"])
        return (
            k >= 2
            or (k == 1 and b["area"] > DISC_ALONE_MAX_PX)
            or (k == 0 and b["area"] > OUTLINE_ALONE_MAX_PX)
        )

    out: dict[str, list[dict]] = {k: [] for k in SERIES}
    lone_fills = []
    for b in blobs:
        if b["i"] in text or not b["glyph_sized"] or merged(b):
            continue
        lone_fills.append(b["fill"])
        if b["fill"] > 0.62:
            out["realtime"].append(b)
        elif b["hole"] >= TRIANGLE_HOLE_PX:
            out["ge_unbalanced"].append(b)
        else:
            out["ge_status81"].append(b)

    # The disc's radius as the rays see it, from the lone discs; and each lone disc fitted the
    # same way, against its centroid, as a check on the fit.
    r0 = float(np.sqrt(np.median([b["area"] for b in out["realtime"]]) / np.pi)) + 0.5
    check = []
    for b in out["realtime"]:
        sl = b["sl"]
        own = lab[sl] == b["i"]
        (y0, x0) = b["cores"][0] if b["cores"] else (b["sy"] - sl[0].start, b["sx"] - sl[1].start)
        fy, fx, _, _ = _disc_fit(own, y0, x0, r0)
        check.append(float(np.hypot(fy + sl[0].start - b["sy"], fx + sl[1].start - b["sx"])))

    # Second pass: the merged blobs, taken apart.
    tpl = _templates(sub, out)
    report = {"merged": 0, "discs": 0, "ge_unbalanced": 0, "ge_status81": 0, "left_px": 0}
    pad = 24
    for b in blobs:
        if b["i"] in text or not (
            (b["glyph_sized"] and merged(b))
            or (not b["glyph_sized"] and (b["cores"] or b["area"] > OUTLINE_ALONE_MAX_PX))
        ):
            continue
        sl = b["sl"]
        y0 = max(sl[0].start - pad, 0)
        x0 = max(sl[1].start - pad, 0)
        win = (slice(y0, sl[0].stop + pad), slice(x0, sl[1].stop + pad))
        own = lab[win] == b["i"]
        hidden = np.zeros(own.shape, dtype=bool)
        yy, xx = np.indices(own.shape)
        report["merged"] += 1
        for cy, cx in b["cores"]:
            sy0, sx0 = cy + sl[0].start - y0, cx + sl[1].start - x0
            fy, fx, fr, _ = _disc_fit(own, sy0, sx0, r0)
            if np.hypot(fy - sy0, fx - sx0) > DISC_RUNAWAY_PX:
                raise RuntimeError(
                    f"disc at ({fx + x0 + left + margin:.0f}, {fy + y0 + top + margin:.0f}) "
                    f"ran {np.hypot(fy - sy0, fx - sx0):.1f} px from its core"
                )
            hidden |= np.hypot(yy - fy, xx - fx) <= fr + 1.5
            out["realtime"].append({"x": fx + x0 + left + margin, "y": fy + y0 + top + margin})
            report["discs"] += 1
        found, unexplained = _outlines(own & sub[win] & ~hidden, hidden, tpl)
        for key, gy, gx in found:
            out[key].append({"x": gx + x0 + left + margin, "y": gy + y0 + top + margin})
            report[key] += 1
        report["left_px"] = max(report["left_px"], unexplained)

    stats = {
        "lone_fills": lone_fills,
        "disc_fit_vs_centroid_px": max(check) if check else 0.0,
        "disc_radius_px": r0,
        **report,
    }
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
        # The plot box: the two strongest ink columns, and the outermost rows that span the
        # plot once the page's tilt is sheared out (recipe, traps 3 and 4).
        L, R, T, B = find_frame_box(ink)
        fr = (int(L), int(T), int(R), int(B))
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
                pos = find_ticks(ink, fits[side], 1, inward, left + 4, right - 3, TICK_EXCESS_PX)
                if pos.size < 3:
                    tick_report.append(f"{side}: {pos.size} ticks, too few to check")
                    continue
                other = fits[side][0] * pos + fits[side][1]
                vals = x_lin(pos, other)
                t = u_of(pos, other)
                grid, tol = TICK_GRID[f.num], TICK_RESIDUAL_TOL[f.num]
            else:
                pos = find_ticks(ink, fits[side], 0, inward, top + 4, bottom - 3, TICK_EXCESS_PX)
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
        if stats["disc_fit_vs_centroid_px"] > DISC_CHECK_PX:
            failures.append(
                f"Figure {f.num}: the disc fit lands {stats['disc_fit_vs_centroid_px']:.2f} px "
                f"from a lone disc's centroid, against {DISC_CHECK_PX}"
            )
        fills = np.array(stats["lone_fills"])
        print(f"  {fills.size} lone glyphs, fill ratio {fills.min():.2f}..{fills.max():.2f}")
        print(
            f"  {stats['merged']} merged blobs taken apart: {stats['discs']} discs, "
            f"{stats['ge_unbalanced']} triangles, {stats['ge_status81']} crosses; at most "
            f"{stats['left_px']} px of a blob's ink left unexplained"
        )
        print(
            f"  disc fit on the lone discs (radius {stats['disc_radius_px']:.2f} px): within "
            f"{stats['disc_fit_vs_centroid_px']:.3f} px of their centroids"
        )
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
