"""Digitize Figures 9 and 10 -- the open-loop fuel-step transients [pdf pp.45-46].

Each figure is six stacked panels sharing one time axis:

    1 WFPH    fuel flow, the input -- solid line only
    2 PS3     compressor discharge static pressure
    3 PCNG    percent gas generator speed
    4 T41     gas generator turbine inlet temperature
    5 T45     power turbine inlet temperature
    6 TORQ45  output torque -- the report's measure of power, since NP is held

Each data panel carries two traces that mean different things, written to separate files:

* a solid line -- Ballin's real-time model, our replication target;
* `+` markers -- the GE status-81 simulation, `REFERENCE STANDARD MODEL DATA`, the
  reference Ballin was himself validating against and does not match perfectly.

Works on the native 300 dpi 1-bit scan, never resampled. The pages are skewed, so every
frame is fitted per panel (`Panel`, `frame_tilt`), and the caption's exactly stated fuel
step checks each page (`check_fuel_panel`). Markers are found by `marker_centres`, the
model trace by `trace_line`.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure.frame import band_line, longest_run_columns, make_sv  # noqa: E402
from figure.raster import native_bitmap  # noqa: E402
from figure.ticks import find_ticks  # noqa: E402

OUT = Path("validation/out/digitize/fig910")
REF = Path("data/reference")

# y-axis ranges, transcribed from each figure's own printed tick labels
# [docs/notes/body-validation.md, sections 2.4 and 2.5]. The two figures do NOT share them.
PANEL_RANGES = {
    9: [
        ("wfph", "WFPH", 250.0, 1000.0),  # input: solid line only, no markers
        ("ps3", "PS3", 100.0, 300.0),
        ("pcng", "PCNG", 80.0, 100.0),
        ("t41", "T41", 2000.0, 3000.0),
        ("t45", "T45", 1250.0, 2250.0),
        ("torq45", "TORQ45", 100.0, 500.0),
    ],
    10: [
        ("wfph", "WFPH", 0.0, 500.0),
        ("ps3", "PS3", 0.0, 200.0),
        ("pcng", "PCNG", 60.0, 100.0),
        ("t41", "T41", 1500.0, 2500.0),
        ("t45", "T45", 1000.0, 2000.0),
        ("torq45", "TORQ45", 0.0, 400.0),
    ],
}
# The caption states the fuel step exactly [pdf pp.45-46: "from 400 to 775 lb_m per hour",
# "from 400 to 125 lb_m per hour"], so the WFPH panel is the one place where the truth behind
# a plotted curve is known: the LEVEL and FLATNESS of the commanded flow either side of the
# step check the whole page (`check_fuel_panel`).
FUEL_STEP = {9: (400.0, 775.0), 10: (400.0, 125.0)}

FUEL_LEVEL_TOL_PCT_HEIGHT = 1.0
FUEL_TILT_TOL_PCT_HEIGHT = 1.0
"""Tolerances for that check, as a percentage of panel height (a read error is a pixel
error, and the panels do not share a data range). The trace's per-column centroid scatters
0.5-1.1 px (1 sigma), 0.15-0.33 % of a 330 px panel; 1.0 % is three sigma of that, and
rejects the 11-14 px skew of an uncorrected page."""
N_PANELS = 6
T_LO, T_HI = 0.0, 5.0


@dataclass
class Panel:
    key: str
    label: str
    top: float
    bottom: float
    left: int
    right: int
    v_lo: float
    v_hi: float
    top_slope: float = 0.0
    bot_slope: float = 0.0
    """Tilt of the panel's two horizontal frame lines, in rows per column, measured at
    `left`. The pages are skewed: on pdf p.46 every frame descends 11-14 px across the panel
    width (3.3-4.2 % of panel height), on p.45 2-5 px. Ignored, that puts a fake ramp of the
    same size into every trace (open question #63); `check_fuel_panel` catches it."""

    left_slope: float = 0.0
    left_at_top: float = 0.0
    right_slope: float = 0.0
    right_at_top: float = 0.0
    """The two VERTICAL frames, fitted over this panel's own rows: columns per row, and
    the column at the panel's top row.

    Each panel needs its own time axis. The page is skewed, so the vertical frames drift
    left going down (pdf p.45: column 603.9 in WFPH, 592.6 in TORQ45); against the one
    page-wide `left`/`right` the upper panels' times read up to 38 ms late, which looks like
    a physical lag in the temperatures and is not one. Width is per panel for the same
    reason."""

    def top_at(self, x: float) -> float:
        """Row of the upper frame at column `x`."""
        return self.top + self.top_slope * (x - self.left)

    def bot_at(self, x: float) -> float:
        """Row of the lower frame at column `x`."""
        return self.bottom + self.bot_slope * (x - self.left)

    def edges(self) -> dict[str, np.ndarray]:
        """The four fitted frame lines in the form `frame.make_sv` takes: each a polynomial,
        highest power first, the verticals in the row and the horizontals in the column."""
        return {
            "left": np.array([self.left_slope, self.left_at_top - self.left_slope * self.top]),
            "right": np.array([self.right_slope, self.right_at_top - self.right_slope * self.top]),
            "top": np.array([self.top_slope, self.top - self.top_slope * self.left]),
            "bottom": np.array([self.bot_slope, self.bottom - self.bot_slope * self.left]),
        }

    @property
    def tilt_pct_height(self) -> float:
        """Mean frame drop across the panel, as a percentage of panel height."""
        rise = 0.5 * (self.top_slope + self.bot_slope) * (self.right - self.left)
        return 100.0 * rise / (self.bottom - self.top)


GLYPH_BAR_MIN_PX = 14
GLYPH_BAR_MAX_PX = 34
GLYPH_BAR_MAX_WID_PX = 20
"""Bounds on a `+` glyph's vertical stroke, in native-scan pixels. On Figure 9 the strokes
run 19-25 px tall; widths reach 17 px where the bar merges with its own horizontal arm. The
upper height bound excludes the trace's own steep segments (85-212 px)."""

LATTICE_MAX_GAP = 2.5
"""How many lattice pitches an end detection may sit from its neighbour before it is
not a sample. The genuine series is contiguous: across the ten marker-bearing panels the
largest internal gap is one pitch, and the single isolated detection sits at five."""

TRACE_MAX_STEP_PX = 8
TRACE_BACKBONE_MIN = 10
TRACE_LINK_MIN = 3
TRACE_LINK_SLOPE_PX = 5.0
TRACE_LINK_OFFSET_PX = 5.0
"""Continuity constraints on the model trace.

Between adjacent columns the trace moves a few pixels; a jump over `TRACE_MAX_STEP_PX`
splits the candidates into blocks. Blocks of `TRACE_BACKBONE_MIN` columns or more are the
backbone; a shorter block (at least `TRACE_LINK_MIN`) joins it if its endpoint reaches an
accepted block's at no more than `TRACE_LINK_SLOPE_PX` per column plus
`TRACE_LINK_OFFSET_PX`. Linked by slope, not by a minimum length: blocks between markers
are only ~14 columns, and a length threshold drops Figure 9's T41 peak (a legitimate
7-column block). Glyph debris needs 8 px per column or more to link."""


FRAME_HALF_BAND_PX = 7
"""How far either side of a nominal frame row to look for the frame's own ink: wide enough
for the up-to-14 px skew it measures, narrow enough to miss the trace where it runs along
the frame. Every frame's ink lands within 7 px of the row the panel search predicts."""

FRAME_MIN_COLUMNS = 50
"""Below this many usable columns the tilt is not measured and the panel stays level."""

VFRAME_HALF_BAND_PX = 25
VFRAME_MIN_ROWS = 40
"""The same, for a vertical frame fitted over one panel's rows. The band starts wide
because the page-wide `left` can be a dozen pixels from the line inside a given panel."""

TICK_EXCESS_PX = 8
"""How far a tick must lengthen the ink run through the lower frame, beyond the frame's own
rolling-median width, to be a major in the x tick check (`ticks.find_ticks`). Measured
2026-10-13 on all twelve panels: the whole- and half-second ticks lengthen it by 9-13 px, the
0.1 s minors by 2-7 (Figure 10's TORQ45 trace runs along its frame from 1.4 to 4.5 s and
lifts both there). At 8 every panel finds all four interior whole seconds; the band count it
replaced missed one or more on seven panels, and all four on Figure 10's T41."""


def frame_tilt(ink: np.ndarray, y0: int, left: int, right: int) -> tuple[float, float, float]:
    """Fit one horizontal frame line: `(slope, row at `left`, residual)`, all in pixels.

    The pages are skewed, so the frames are not horizontal (see `Panel.top_slope`). The
    frame is the only ink across the whole panel width at its height: take the ink centroid
    of a narrow band per column, drop columns with more than 8 inked rows (the trace
    crossing, or a tick), fit a line. The fitted row also replaces the panel search's,
    which is a centroid of the smeared tilted line.
    """
    (a, b), done, res = band_line(
        ink, 1, y0, left + 5, right - 4, 4, (FRAME_HALF_BAND_PX,), 8, FRAME_MIN_COLUMNS
    )
    if not done:
        return 0.0, float(y0), float("nan")
    return a, a * left + b, res


def vframe_fit(ink: np.ndarray, x0: float, y0: int, y1: int) -> tuple[float, float, float]:
    """Fit one vertical frame over one panel's rows: `(slope, column at y0, residual)`.

    Two passes, wide then narrow, for the same reason `frame_tilt` uses two: the starting
    guess is the page-wide extreme column and the true line can be a dozen pixels away.
    """
    (a, b), done, res = band_line(
        ink, 0, x0, y0 + 4, y1 - 3, 3, (VFRAME_HALF_BAND_PX, 6), 10, VFRAME_MIN_ROWS
    )
    if done < 2:
        return 0.0, float(x0), float("nan")
    return a, a * y0 + b, res


def find_panels(ink: np.ndarray, ranges: list) -> list[Panel]:
    """Locate the six panel frames.

    The vertical frame lines are the longest columnar runs, searched away from the page
    edges. A full-width row is not necessarily a frame (the WFPH step's two flat levels are
    full-width lines too, and pairing neighbours shifts every label by one panel), so the
    six tops are found as an arithmetic progression, which the trace lines cannot join.
    """
    h = ink.shape[0]
    left, right = longest_run_columns(ink, 0.08, 0.97, 0.7)

    frac = ink[:, left : right + 1].mean(axis=1)
    lines: np.ndarray = np.zeros(0)
    for thresh in (0.60, 0.50, 0.40, 0.32, 0.25, 0.20):
        rows = np.flatnonzero(frac > thresh)
        if rows.size == 0:
            continue
        groups = np.split(rows, np.flatnonzero(np.diff(rows) > 3) + 1)
        lines = np.array([g.mean() for g in groups if g.size])
        if lines.size >= 2 * N_PANELS:
            break
    if lines.size < 2 * N_PANELS:
        raise RuntimeError(f"only {lines.size} horizontal lines found; need 12")

    # the six tops are the best arithmetic progression among the candidates
    best = None
    for i, a in enumerate(lines):
        for b in lines[i + 1 :]:
            pitch = b - a
            if pitch < 300:
                continue
            want = a + pitch * np.arange(N_PANELS)
            if want[-1] > h:
                continue
            err = 0.0
            got = []
            for t in want:
                j = int(np.argmin(np.abs(lines - t)))
                err += abs(lines[j] - t)
                got.append(float(lines[j]))
            if len(set(got)) < N_PANELS:
                continue
            if best is None or err < best[0]:
                best = (err, got, pitch)
    if best is None:
        raise RuntimeError(f"no six-panel progression among {lines}")
    err, tops, pitch = best
    if err > 6.0 * N_PANELS:
        raise RuntimeError(f"progression fits poorly (total {err:.1f} px): {tops}")

    # Panel height is the same for all six, so take it once (the median top-to-bottom
    # distance): the nearest line below a top can be the WFPH step's level, not a frame.
    cands = []
    for t in tops:
        below = lines[(lines > t + 0.5 * pitch) & (lines < t + 0.95 * pitch)]
        for b in below:
            cands.append(float(b) - t)
    if not cands:
        raise RuntimeError("no candidate panel bottoms found")
    height = float(np.median([c for c in cands if c > 0.6 * pitch] or cands))

    out = []
    for t, (key, label, lo, hi) in zip(tops, ranges, strict=True):
        t0, b0 = int(round(t)), int(round(t + height))
        ts, ty, _ = frame_tilt(ink, t0, left, right)
        bs, by, _ = frame_tilt(ink, b0, left, right)
        # A frame the fit could not reach falls back to the other one rather than to zero:
        # the two frames of a panel are parallel to well under a pixel across the width.
        if ts == 0.0:
            ts = bs
        if bs == 0.0:
            bs = ts
        ls, lx, _ = vframe_fit(ink, left, t0, b0)
        rs, rx, _ = vframe_fit(ink, right, t0, b0)
        # The vertical frames are referenced to the panel's own top row, `ty`, not to the
        # row the fit was anchored at.
        lx += ls * (ty - t0)
        rx += rs * (ty - t0)
        out.append(Panel(key, label, ty, by, left, right, lo, hi, ts, bs, ls, lx, rs, rx))
    return out


def marker_centres(ink: np.ndarray, p: Panel, bar: int = 7, margin: int = 16) -> np.ndarray:
    """Find `+` glyph centres inside a panel.

    A `+` is the only thing with both a horizontal and a vertical stroke at one place: the
    flat trace has only the first, the steep trace only the second. `margin` (16 px) clears
    the frame ticks, which reach 13 px and would otherwise pass that test.

    The centre comes from the VERTICAL stroke alone, the horizontal one only testing
    presence: where a marker sits on the model line the horizontal stroke IS the line, and
    a centre of mass of both pulls the GE sample onto Ballin's curve.
    """
    # The rectangle spans the whole data region and the tilted frames are blanked out of it:
    # shrinking it to clear the tilt loses markers near the bottom of an axis.
    lo = int(np.floor(min(p.top_at(p.left), p.top_at(p.right)))) + margin
    hi = int(np.ceil(max(p.bot_at(p.left), p.bot_at(p.right)))) - margin
    x0, x1 = p.left + margin, p.right - margin
    sub = ink[lo:hi, x0:x1].copy()
    _rows = np.arange(lo, hi)[:, None]
    _cols = np.arange(x0, x1)[None, :]
    sub[_rows < (p.top + p.top_slope * (_cols - p.left)) + margin] = False
    sub[_rows > (p.bottom + p.bot_slope * (_cols - p.left)) - margin] = False
    horiz = ndimage.binary_opening(sub, structure=np.ones((1, bar)))
    vert = ndimage.binary_opening(sub, structure=np.ones((bar, 1)))
    horiz_near = ndimage.binary_dilation(horiz, structure=np.ones((5, 5)))
    lab, n = ndimage.label(vert, structure=np.ones((3, 3)))
    if n == 0:
        return np.zeros((0, 2))
    pts = []
    for i, sl in enumerate(ndimage.find_objects(lab), start=1):
        hgt = sl[0].stop - sl[0].start
        wid = sl[1].stop - sl[1].start
        if not (GLYPH_BAR_MIN_PX <= hgt <= GLYPH_BAR_MAX_PX) or wid > GLYPH_BAR_MAX_WID_PX:
            continue  # a steep segment of the trace, or a frame artifact -- not a bar
        blob = lab[sl] == i
        if not (blob & horiz_near[sl]).any():
            continue  # a bare vertical stroke with no arm is not a glyph
        cy, cx = ndimage.center_of_mass(blob)
        pts.append((cx + sl[1].start + x0, cy + sl[0].start + lo))
    return np.array(sorted(pts)) if pts else np.zeros((0, 2))


def snap_to_lattice(pts: np.ndarray, pitch_hint: float = 33.0) -> tuple[np.ndarray, dict]:
    """Keep one detection per sample slot, using the markers' own regular spacing.

    The reference series is sampled at a fixed interval, so a split glyph or a spurious
    blob nearby falls between lattice slots; the detection nearest each slot centre is kept.
    """
    if pts.shape[0] < 5:
        return pts, {"slots": 0, "pitch": float("nan"), "dropped": 0}
    x = np.sort(pts[:, 0])
    d = np.diff(x)
    pitch = float(np.median(d[(d > 0.6 * pitch_hint) & (d < 1.6 * pitch_hint)]))
    if not np.isfinite(pitch):
        pitch = pitch_hint

    # phase that best aligns every detection to the lattice
    phases = ((pts[:, 0] % pitch) / pitch) * 2.0 * np.pi
    phase = float(np.angle(np.exp(1j * phases).mean()) / (2.0 * np.pi) * pitch)

    slot = np.round((pts[:, 0] - phase) / pitch).astype(int)
    keep = []
    for sl in np.unique(slot):
        cand = pts[slot == sl]
        if cand.shape[0] == 1:
            keep.append(cand[0])
            continue
        # the detection closest to its slot centre is the glyph; the rest are splinters
        want = phase + sl * pitch
        keep.append(cand[np.argmin(np.abs(cand[:, 0] - want))])
    out = np.array(sorted(keep, key=lambda q: q[0]))

    # The series is sampled contiguously, so an end detection more than LATTICE_MAX_GAP
    # pitches from its neighbour is not a sample. Trimmed from both ends.
    n_trim = 0
    while out.shape[0] > 5 and out[-1, 0] - out[-2, 0] > LATTICE_MAX_GAP * pitch:
        out = out[:-1]
        n_trim += 1
    while out.shape[0] > 5 and out[1, 0] - out[0, 0] > LATTICE_MAX_GAP * pitch:
        out = out[1:]
        n_trim += 1
    return out, {
        "slots": len(out),
        "pitch": pitch,
        "dropped": len(pts) - len(out),
        "trimmed": n_trim,
    }


def trace_line(
    ink: np.ndarray, p: Panel, marks: np.ndarray, clear: float = 9.0, margin: int = 16
) -> np.ndarray:
    """Follow the solid model trace, one sample per column, markers masked out.

    A column is a sample only if, between the tilted frames less `margin`, it holds exactly
    one ink run and that run is thin: legend text, a marker, a crossing or a near-vertical
    segment is skipped, leaving a gap rather than a wrong sample. Count the runs BEFORE
    filtering by width: filtering first lets a glyph's detached bar-top win over the merged
    bar-and-line run (the spikes of open question #51).

    `margin` (16 px) clears the frame ticks, which reach 13 px at the 99th percentile; a
    wider one clips the trace where it settles near a frame (Figure 9's PCNG).

    A thin single run can still be a bar-top or a scrap of legend rule, so the candidates
    are split into connected blocks and linked by slope (see the `TRACE_*` constants).
    """
    cand = []
    for x in range(p.left + margin, p.right - margin):
        # The window follows the tilted frames; a fixed one lets the frame's own ink into
        # the column at the far end of the panel.
        y0 = int(np.ceil(p.top_at(x))) + margin
        y1 = int(np.floor(p.bot_at(x))) - margin
        col = ink[y0:y1, x]
        ys = np.flatnonzero(col)
        if ys.size == 0:
            continue
        runs = np.split(ys, np.flatnonzero(np.diff(ys) > 1) + 1)
        if len(runs) != 1 or runs[0].size > 6:
            continue  # legend text, a marker, a crossing, or a near-vertical segment
        y = float(runs[0].mean()) + y0
        if marks.size:
            near = marks[np.abs(marks[:, 0] - x) < clear]
            if near.size and np.min(np.abs(near[:, 1] - y)) < clear:
                continue
        cand.append((x, y))
    if not cand:
        return np.zeros((0, 2))

    # Continuity: split wherever adjacency breaks (a skipped column, or a jump larger than
    # the trace can make in one column), then link the blocks by slope.
    blocks: list[list[tuple[int, float]]] = [[cand[0]]]
    for (x0, y0), (x1, y1) in zip(cand, cand[1:], strict=False):
        if x1 - x0 == 1 and abs(y1 - y0) <= TRACE_MAX_STEP_PX:
            blocks[-1].append((x1, y1))
        else:
            blocks.append([(x1, y1)])

    accepted = [len(b) >= TRACE_BACKBONE_MIN for b in blocks]
    if not any(accepted):  # no backbone: take the longest block and nothing else
        accepted[max(range(len(blocks)), key=lambda i: len(blocks[i]))] = True

    def links(a: list, b: list) -> bool:
        """Do two blocks join with a slope the trace could actually have?"""
        (xa, ya), (xb, yb) = (a[-1], b[0]) if a[0][0] < b[0][0] else (a[0], b[-1])
        dx = max(abs(xb - xa), 1)
        return abs(yb - ya) <= TRACE_LINK_SLOPE_PX * dx + TRACE_LINK_OFFSET_PX

    changed = True
    while changed:  # grow outward from the backbone until nothing more attaches
        changed = False
        for i, b in enumerate(blocks):
            if accepted[i]:
                continue
            for j in (i - 1, i + 1):
                if (
                    0 <= j < len(blocks)
                    and accepted[j]
                    and len(b) >= TRACE_LINK_MIN
                    and links(blocks[j], b)
                ):
                    accepted[i] = True
                    changed = True
                    break

    kept = [q for b, ok in zip(blocks, accepted, strict=True) if ok for q in b]
    kept.sort()
    return np.array(kept) if kept else np.zeros((0, 2))


def calibrate(ink: np.ndarray, p: Panel) -> tuple:
    """Map pixels to data for one panel, and check the map against the printed ticks.

    Both frames of every axis carry a printed value, so they pin the (straight) map and the
    interior major ticks are the held-out check. `v_of` needs the column as well as the row,
    because the frames are tilted (see `Panel.top_slope`). The map is `frame.make_sv` on the
    four straight lines, not `frame.unit_square`'s homography: across these twelve panels the
    two differ by up to 2.1 px in time (Figure 9 WFPH), measured 2026-10-13.
    """
    to_sv, _ = make_sv(p.edges())

    def x_of(px, py):
        return T_LO + to_sv(px, py)[0] * (T_HI - T_LO)

    def v_of(py, px):
        return p.v_lo + to_sv(px, py)[1] * (p.v_hi - p.v_lo)

    # The ticks on the lower frame point inward (up), and nothing lies below it.
    tp = find_ticks(ink, p.edges()["bottom"], 1, -1, p.left, p.right + 1, TICK_EXCESS_PX)
    got = x_of(tp, p.bot_at(tp))
    err = float("nan")
    if got.size:
        want = np.arange(T_LO, T_HI + 1e-9, 1.0)
        near = [min(abs(got - w)) for w in want if np.any(np.abs(got - w) < 0.3)]
        if near:
            err = float(max(near))
    return x_of, v_of, err


def check_fuel_panel(fig: int, ts: np.ndarray, vs: np.ndarray, span: float) -> list[str]:
    """The WFPH panel against the caption. Returns a list of failures, empty if it passes.

    The LEVEL of each plateau against the stated flow catches a wrong axis range; the TILT
    catches page skew, since a commanded constant must be flat.
    """
    lo, hi = FUEL_STEP[fig]
    mid = 0.5 * (lo + hi)
    # the FIRST crossing of the half level, not the nearest sample to it: the post-step
    # plateau of Figure 10 sits only 137 lbm/hr below the half level and wins an argmin.
    cross = np.flatnonzero(np.diff(np.sign(vs - mid)) != 0)
    edge = float(ts[int(cross[0])]) if cross.size else float("nan")
    out: list[str] = []
    for name, mask, want in (
        ("pre-step", (ts > 0.15) & (ts < edge - 0.05), lo),
        ("post-step", ts > edge + 0.60, hi),
    ):
        w, u = ts[mask], vs[mask]
        if w.size < 20:
            out.append(f"Figure {fig} WFPH {name}: only {w.size} samples, cannot check")
            continue
        level = 100.0 * (float(np.median(u)) - want) / span
        slope = float(np.polyfit(w, u, 1)[0]) * (T_HI - T_LO) / span * 100.0
        tag = "ok"
        if abs(level) > FUEL_LEVEL_TOL_PCT_HEIGHT:
            out.append(
                f"Figure {fig} WFPH {name} level {level:+.2f} % of panel height "
                f"(median {np.median(u):.2f} against the caption's {want:.0f})"
            )
            tag = "FAIL"
        if abs(slope) > FUEL_TILT_TOL_PCT_HEIGHT:
            out.append(
                f"Figure {fig} WFPH {name} tilt {slope:+.2f} % of panel height over "
                f"the record -- a commanded constant is not flat, so the panel is skewed"
            )
            tag = "FAIL"
        print(f"    caption check {name:9} level {level:+6.2f} %h  tilt {slope:+6.2f} %h  [{tag}]")
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    REF.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []
    for page, fig in ((45, 9), (46, 10)):
        ink = native_bitmap(page)
        panels = find_panels(ink, PANEL_RANGES[fig])
        print(f"\nFigure {fig} (pdf p.{page})")
        for p in panels:
            x_of, v_of, terr = calibrate(ink, p)
            has_markers = p.key != "wfph"
            marks = np.zeros((0, 2))
            info = {"slots": 0, "pitch": float("nan")}
            if has_markers:
                marks, info = snap_to_lattice(marker_centres(ink, p))
            line = trace_line(ink, p, marks)

            for kind, pts in (("model", line), ("reference", marks)):
                if pts.shape[0] == 0:
                    continue
                path = REF / f"fig{fig:02d}_{p.key}_{kind}.csv"
                ts_all = x_of(pts[:, 0], pts[:, 1])
                vs_all = v_of(pts[:, 1], pts[:, 0])
                dropped: list[int] = []
                if ts_all.size > 8:
                    _g = np.diff(ts_all)
                    _m = float(np.median(_g))
                    if _m > 0:
                        _big = np.flatnonzero(_g > 50.0 * _m)
                        if _big.size:
                            _cut = int(_big[-1]) + 1
                            if len(ts_all) - _cut <= 8:
                                dropped = list(range(_cut, len(ts_all)))
                with path.open("w") as fh:
                    fh.write(f"# source: TM-100991 pdf p.{page}, Figure {fig}, panel {p.label}\n")
                    what = (
                        "Ballin real-time model (solid line) -- the replication target"
                        if kind == "model"
                        else "GE reference standard model data (+ markers) -- "
                        "the hardware reference Ballin was himself validating against"
                    )
                    fh.write(f"# quantity: {p.label} vs time, {what}\n")
                    fh.write(
                        "# method: digitized -- native 300 dpi 1-bit scan, "
                        "tools/digitize_fig910.py\n"
                    )
                    fh.write(
                        f"# t_s: time, seconds (0 to 5.0). The vertical frames are fitted "
                        f"over THIS panel's rows: left column {p.left_at_top:.2f} "
                        f"(slope {p.left_slope * 1e3:+.2f} m px/row), right "
                        f"{p.right_at_top:.2f} (slope {p.right_slope * 1e3:+.2f}), width "
                        f"{p.right_at_top - p.left_at_top:.2f} px against the page-wide "
                        f"{p.right - p.left}. The page is skewed; see `Panel.left_slope`.\n"
                    )
                    fh.write(
                        f"# value: {p.label}, axis range {p.v_lo} to {p.v_hi}, "
                        f"no unit printed on the figure\n"
                    )
                    fh.write(f"# points: {pts.shape[0] - len(dropped)}\n")
                    fh.write(f"# x tick check: worst interior major tick off by {terr:.4f} s\n")
                    fh.write(
                        f"# frame tilt corrected: top "
                        f"{p.top_slope * (p.right - p.left):+.2f} px, bottom "
                        f"{p.bot_slope * (p.right - p.left):+.2f} px across the panel "
                        f"({p.tilt_pct_height:+.2f} % of panel height). "
                        f"The page is skewed; see `frame_tilt`.\n"
                    )
                    if kind == "reference":
                        fh.write(
                            f"# sample lattice: pitch {info['pitch']:.2f} px "
                            f"= {x_of(p.left_at_top + info['pitch'], p.top)[0] - T_LO:.4f} s\n"
                        )
                    # Past the end of the plotted curve the walker re-acquires ink (a frame
                    # line, a legend rule, a stray glyph). Every such artifact sits behind an
                    # x-gap of more than 130x the median sample interval, every legitimate
                    # continuation behind 20-25x, so `dropped` (above) cuts at 50x. Dropped
                    # rather than kept, and recorded in the header (open question #48).
                    ts = ts_all
                    if dropped:
                        gaps = np.diff(ts)
                        med = float(np.median(gaps))
                        cut = dropped[0]
                        fh.write(
                            f"# dropped: {len(dropped)} trailing sample(s) at "
                            f"t = {', '.join(f'{ts[i]:.5f}' for i in dropped)} s, behind an "
                            f"x-gap of {gaps[cut - 1] / med:.0f}x the median sample interval "
                            f"({med:.5f} s). The curve ends at t={ts[cut - 1]:.5f} s with "
                            f"value {vs_all[cut - 1]:.5f}; "
                            f"the dropped samples read "
                            f"{', '.join(f'{vs_all[i]:.2f}' for i in dropped)}. "
                            f"Re-acquired ink past the end of the plotted curve, not data. "
                            f"See open question #48.\n"
                        )
                    fh.write("# NOT transcribed. Values carry read error.\n")
                    fh.write("t_s,value\n")
                    keep = set(range(len(pts))) - set(dropped)
                    for i in range(len(pts)):
                        if i in keep:
                            fh.write(f"{ts_all[i]:.5f},{vs_all[i]:.5f}\n")
            print(
                f"  {p.label:7} model {line.shape[0]:5d} pts | "
                f"reference {marks.shape[0]:3d} | x tick worst {terr:.4f} s | "
                f"tilt {p.tilt_pct_height:+.2f} %h"
            )
            if p.key == "wfph" and line.shape[0]:
                failures += check_fuel_panel(
                    fig,
                    x_of(line[:, 0], line[:, 1]),
                    v_of(line[:, 1], line[:, 0]),
                    abs(p.v_hi - p.v_lo),
                )
    if failures:
        print("\nCAPTION CHECK FAILED -- the digitization does not reproduce the printed step:")
        for f in failures:
            print(f"  {f}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
