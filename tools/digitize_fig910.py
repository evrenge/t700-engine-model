"""Digitize Figures 9 and 10 -- the open-loop fuel-step transients.

These are the report's transient validation evidence and the target for our real-time
model. Each figure is **six stacked panels** sharing one time axis [pdf pp.45-46]:

    1 WFPH    fuel flow, the input -- solid line only
    2 PS3     compressor discharge static pressure
    3 PCNG    percent gas generator speed
    4 T41     gas generator turbine inlet temperature
    5 T45     power turbine inlet temperature
    6 TORQ45  output torque -- the report's measure of power, since NP is held

Each data panel carries **two traces that mean different things**:

* a **solid line** -- Ballin's real-time model. This is our replication target: we are
  reproducing his model, so his trace is what our output must match.
* **`+` markers** -- the GE performance-standard status-81 simulation, labelled
  `REFERENCE STANDARD MODEL DATA`. This is the hardware-side reference Ballin was himself
  validating against, and he does not match it perfectly.

They are written to separate files. Merging them would destroy the distinction between
"what we must reproduce" and "what he was trying to reproduce".

Method is the one accumulated across Appendix A (`docs/notes/digitizing-recipe.md`):
native 300 dpi 1-bit via `pdfimages`, never resampled; per-axis polynomial order chosen by
the held-out frame-edge test; the scan bows differently on every page.

Separating the two traces: a `+` glyph has a horizontal stroke several pixels wide, while
the model trace is a thin near-vertical-gradient line. A horizontal opening isolates the
marker bars; the residual ink is the line.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from digitize_native import native_bitmap  # noqa: E402

OUT = Path("validation/out/digitize/fig910")
REF = Path("data/reference")

# y-axis ranges, transcribed from each figure's own printed tick labels
# [docs/notes/body-validation.md, sections 2.4 and 2.5].
#
# **The two figures do NOT share axis ranges.** A first pass assumed they did and applied
# Figure 9's to both; Figure 10's WFPH runs 0-500, not 250-1000, so every Figure 10 value
# came out wrong -- the fuel step read +99 % against a caption that states it exactly. The
# caption's stated step levels are what caught it, which is why that check is run.
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
# "from 400 to 125 lb_m per hour"], so the WFPH panel is the one place on either page where
# the truth behind a plotted curve is known. That makes it the calibration check for the
# whole page: the commanded flow is constant either side of the step, so both its LEVEL and
# its FLATNESS are checkable.
#
# It was declared here and used nowhere until 2026-09-14, with a comment saying the check
# "is run". It was not, and running it failed: on Figure 10 the pre-step level read +1.70 %
# and the post-step level fell 4.0 % of panel height across the record, which is the page
# skew now corrected in `frame_tilt`. A gate that cannot fail is not a gate.
FUEL_STEP = {9: (400.0, 775.0), 10: (400.0, 125.0)}

FUEL_LEVEL_TOL_PCT_HEIGHT = 1.0
FUEL_TILT_TOL_PCT_HEIGHT = 1.0
"""Tolerances for that check, both as a percentage of panel height -- the unit the error
is actually made in, since a read error is a pixel error and the panels do not share a
data range.

Set from the trace's own centroid noise, not from what happens to pass: the model trace is
4-6 px thick and its per-column centroid scatters 0.5-1.1 px (1 sigma) about a fitted
frame, which on a 330 px panel is 0.15-0.33 percent of height. 1.0 percent is three sigma
of that, and 3.3 px -- it admits ordinary line-centroid error and rejects the 11-14 px
skew that was there before."""
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
    `left`. **The pages are skewed and the frames are not horizontal**: on pdf p.46 every
    frame descends 11-14 px across the panel width, 3.3-4.2 percent of panel height, and
    on pdf p.45 2-5 px, 0.7-1.6 percent. Ignoring that put a ramp of exactly that size
    into every digitized trace on both pages -- a monotone fake decay on top of the real
    one, worst on Figure 10, which is the chop transient.

    It was found by checking the one panel whose truth is printed: the caption states the
    fuel step exactly, and the WFPH trace -- a commanded constant either side of the step
    -- read 4.0 percent of panel height lower at t=5 than at t=0 on Figure 10, matching
    the frames' own 3.35 percent to within the trace's centroid noise. `check_fuel_panel`
    below is that check, run every time rather than remembered."""

    left_slope: float = 0.0
    left_at_top: float = 0.0
    right_slope: float = 0.0
    right_at_top: float = 0.0
    """The two VERTICAL frames, fitted over this panel's own rows: columns per row, and
    the column at the panel's top row.

    **Each panel needs its own time axis.** `find_panels` returns one `left`/`right` pair
    for the whole page -- the extreme columns of the longest columnar runs -- and until
    2026-09-14 every panel was calibrated against that one pair. The pages are skewed, so
    the shared vertical frames drift left going down: on pdf p.45 the left frame sits at
    column 603.9 in the WFPH panel and 592.6 in TORQ45, and on pdf p.46 at 494.7 and
    477.6. Against a single `left` taken at the bottom, the upper panels' times read up to
    **+38 ms too large**, and the error varies from panel to panel down the page.

    That is not a rounding nuisance on a transient sampled every 7 ms. It showed up as our
    model appearing to respond *earlier* than Ballin's by an amount that grew down the
    page -- 2 ms on PCNG, 37 on T41, 49 on T45, 54 on TORQ45 -- which reads as a physical
    lag in the temperatures and is not one.

    The panel's width is per-panel too, for the same reason: 1658.9 px on Figure 9's WFPH
    against the shared 1671, a 0.75 percent stretch worth 37 ms over the 5 s record."""

    def left_at(self, y: float) -> float:
        return self.left_at_top + self.left_slope * (y - self.top)

    def right_at(self, y: float) -> float:
        return self.right_at_top + self.right_slope * (y - self.top)

    def top_at(self, x: float) -> float:
        """Row of the upper frame at column `x`."""
        return self.top + self.top_slope * (x - self.left)

    def bot_at(self, x: float) -> float:
        """Row of the lower frame at column `x`."""
        return self.bottom + self.bot_slope * (x - self.left)

    @property
    def tilt_pct_height(self) -> float:
        """Mean frame drop across the panel, as a percentage of panel height."""
        rise = 0.5 * (self.top_slope + self.bot_slope) * (self.right - self.left)
        return 100.0 * rise / (self.bottom - self.top)


GLYPH_BAR_MIN_PX = 14
GLYPH_BAR_MAX_PX = 34
GLYPH_BAR_MAX_WID_PX = 20
"""Bounds on a `+` glyph's vertical stroke, in native-scan pixels.

Measured across the five marker-bearing panels of Figure 9: the strokes run 19-25 px
tall (median 21). Widths run 4-6 px where the bar stands alone but reach 17 px where it
merges with its own horizontal arm, so the width bound is 20: at 10 it rejected eight real
markers on Figure 9's Ps3 panel and two on Figure 10's T41. The outliers at 85-212 px
are the trace's own steep segments, which the upper bound is there to exclude."""

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

The trace is one connected thin line, so between ADJACENT columns it moves at most a few
pixels; a jump of more than `TRACE_MAX_STEP_PX` between neighbouring columns is not the
trace, and splits the candidates into blocks.

A minimum block LENGTH was tried first and is wrong. Blocks are only about 14 columns
long -- the markers sit every 33 px and the mask around each removes some 18 columns -- so
a threshold anywhere near that length throws away real trace. It threw away Figure 9's T41
peak, a legitimate 7-column block holding the 2790 degR maximum, which is the single most
quoted number in the whole comparison.

So blocks are linked by SLOPE instead. Blocks of `TRACE_BACKBONE_MIN` columns or more are
the backbone; a shorter block joins it if its endpoint reaches an accepted block's
endpoint at no more than `TRACE_LINK_SLOPE_PX` per column plus a fixed
`TRACE_LINK_OFFSET_PX`. The peak block connects at 1.8 px per column. Glyph debris does
not: it is 4-6 columns wide and sits tens of pixels off the curve, needing 8 px per column
or more. The step riser needs no special case -- the levels either side of it are both
backbone, so both are kept without ever having to link across it."""


def longest_run(col: np.ndarray) -> int:
    best = cur = 0
    for v in col:
        cur = cur + 1 if v else 0
        if cur > best:
            best = cur
    return best


FRAME_HALF_BAND_PX = 7
"""How far either side of a nominal frame row to look for the frame's own ink.

The skew moves a frame by up to 14 px end to end, so the search band has to be wider than
the tilt it is measuring -- but not so wide that it catches the trace where the trace runs
along the frame. 7 px is measured: every frame's ink lands within 7 px of the row the
arithmetic-progression fit puts it at, and the residual about the fitted tilt is 0.45-1.45
px on all twenty-four frame lines."""

FRAME_MIN_COLUMNS = 50
"""Below this many usable columns the tilt is not measured and the panel stays level."""

VFRAME_HALF_BAND_PX = 25
VFRAME_MIN_ROWS = 40
"""The same, for a vertical frame fitted over one panel's rows. The band starts wide
because the page-wide `left` can be a dozen pixels from the line inside a given panel."""


def frame_tilt(ink: np.ndarray, y0: int, left: int, right: int) -> tuple[float, float, float]:
    """Fit one horizontal frame line: `(slope, row at `left`, residual)`, all in pixels.

    The frame is the only thing inked across the whole panel width at that height, so
    taking the ink centroid of a narrow vertical band per column and fitting a line
    recovers both its tilt and its true position. Columns whose band holds more than 8
    inked rows are dropped: that is the trace crossing the frame, or a tick, not the frame.

    The fitted row replaces the one the panel search produced, which was a centroid of the
    *smeared* tilted line plus a modal panel height applied uniformly. That approximation
    left the Figure 10 fuel panel reading 1.2 percent of panel height high even after the
    tilt was taken out; fitting each frame on its own removes it.
    """
    xs: list[int] = []
    ys: list[float] = []
    for x in range(left + 5, right - 4, 4):
        band = ink[y0 - FRAME_HALF_BAND_PX : y0 + FRAME_HALF_BAND_PX + 1, x]
        idx = np.flatnonzero(band)
        if idx.size and idx.size <= 8:
            xs.append(x)
            ys.append(y0 - FRAME_HALF_BAND_PX + float(idx.mean()))
    if len(xs) < FRAME_MIN_COLUMNS:
        return 0.0, float(y0), float("nan")
    a, b = np.polyfit(xs, ys, 1)
    res = float(np.std(np.asarray(ys) - (a * np.asarray(xs) + b)))
    return float(a), float(a * left + b), res


def vframe_fit(ink: np.ndarray, x0: float, y0: int, y1: int) -> tuple[float, float, float]:
    """Fit one vertical frame over one panel's rows: `(slope, column at y0, residual)`.

    Two passes, wide then narrow, for the same reason `frame_tilt` uses two: the starting
    guess is the page-wide extreme column and the true line can be a dozen pixels away.
    """
    fitted = (0.0, float(x0))
    half = VFRAME_HALF_BAND_PX
    res = float("nan")
    for _pass in (0, 1):
        ys: list[float] = []
        xs: list[float] = []
        for y in range(y0 + 4, y1 - 3, 3):
            centre = fitted[0] * y + fitted[1]
            lo = int(round(centre - half))
            band = ink[y, max(lo, 0) : lo + 2 * half + 1]
            idx = np.flatnonzero(band)
            if idx.size and idx.size <= 10:
                ys.append(float(y))
                xs.append(lo + float(idx.mean()))
        if len(ys) < VFRAME_MIN_ROWS:
            return 0.0, float(x0), float("nan")
        a, b = np.polyfit(ys, xs, 1)
        fitted = (float(a), float(b))
        res = float(np.std(np.asarray(xs) - (a * np.asarray(ys) + b)))
        half = 6
    return fitted[0], fitted[0] * y0 + fitted[1], res


def find_panels(ink: np.ndarray, ranges: list) -> list[Panel]:
    """Locate the six panel frames.

    The vertical frame lines are the longest columnar runs, searched away from the page
    edges so the scan's left-margin junk cannot pose as a frame.

    Finding the horizontal frames needs care. A row inked across the plot width is not
    necessarily a frame: panel 1's WFPH trace is a *step*, so its two flat levels are
    full-width lines too. A first attempt that paired adjacent lines by height accepted
    that trace as panel 2's top and shifted every label by one panel -- the kind of error
    that silently mislabels an entire dataset.

    The fix is structural. The six panels have equal height *and* equal pitch, so the
    tops form an arithmetic progression. Search for the progression rather than pairing
    neighbours, and the trace lines cannot join it.
    """
    h, w = ink.shape
    m0, m1 = int(w * 0.08), int(w * 0.97)
    cols = np.array([longest_run(ink[:, x]) if m0 <= x < m1 else 0 for x in range(w)])
    strong = np.flatnonzero(cols > 0.7 * cols.max())
    left, right = int(strong.min()), int(strong.max())

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

    # Panel height, like panel pitch, is the same for all six -- so take it once, as the
    # modal top-to-bottom distance, and apply it uniformly.
    #
    # Choosing each bottom independently as "the nearest line below the top" fails on the
    # WFPH panel, whose trace is a step: on Figure 10 the post-step level at 125 lbm/hr
    # sits where a frame would plausibly be, and taking it gave that panel a height of 246
    # px against ~330 everywhere else. Every Figure 10 fuel value was then wrong by 6 %,
    # which the caption's stated levels caught.
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

    A `+` is the only thing on the page with **both** a horizontal stroke and a vertical
    one at the same place. That is the discriminator:

    * the model trace where it runs flat has horizontal extent but no vertical extent;
    * the model trace where it runs steep has vertical extent but no horizontal extent;
    * a `+` has both.

    A first attempt used the horizontal opening alone and returned 146 "markers" in the
    WFPH panel -- which carries no markers at all, only a step whose two flat levels are
    long horizontal runs. Requiring both strokes removes that entire class of false
    positive.

    `margin` removes the other class. Every frame carries **tick marks**, and a tick
    crossing its own frame line presents a vertical and a horizontal stroke at the same
    place -- indistinguishable from a glyph by the test above. Ticks reach 13 px at the
    99th percentile across all twelve panels, so 16 px inside every frame clears them.

    ## The centre comes from the VERTICAL stroke, and that is the fix

    Until 2026-09-12 the centre was the centre of mass of `horiz & vert`. Where a marker
    sits **on** the model line -- which is most of them, since Ballin plotted the GE data
    against his own curve -- the horizontal stroke in that intersection *is the model
    line*, so the centre of mass was pulled onto the line and the GE sample inherited the
    model's value. That is the "GE data accounted as Ballin" half of the contamination:
    on Figure 10, three PCNG markers, and one each on Ps3, T41, T45 and TORQ45, read the
    model curve to within a pixel.

    The vertical stroke belongs only to the glyph: a `+`'s bar is symmetric about its
    centre, and where the bar merges with the line the merged run is still centred on the
    glyph. So the centre is taken from the vertical blob alone, with the horizontal stroke
    kept only as a presence test.

    The height bounds are measured, not guessed. Across the five marker-bearing panels of
    Figure 9 the vertical strokes run 19-25 px (median 21) and 4-6 px wide; the outliers
    at 85-212 px are the trace's own steep segments, which is what the upper bound excludes.
    """
    # The rectangle spans the whole data region -- shrinking it to the part no tilted
    # frame reaches costs 28 px of a 330 px panel and lost five of Figure 10's sixteen
    # TORQ45 markers, which sit near the bottom of that axis. Instead the rectangle stays
    # and the two tilted frames are blanked out of it, which is what `margin` was doing
    # before the tilt was known.
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

    The reference series is sampled at a fixed interval, so its markers sit on a lattice.
    Where the model's trace touches a glyph the detector can split it in two or add a
    spurious blob nearby; both show up as extra points *between* lattice slots. Fitting
    the lattice and taking the best detection per slot removes them without inventing
    anything -- the structure is printed on the page.

    This is the same principle the appendix work settled on: assert on structure, not on
    count. A stable count can be stably wrong.
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

    # The series is sampled contiguously, so a detection separated from the rest by
    # several empty slots is not a sample. Exactly one exists across the ten
    # marker-bearing panels: Figure 10's T45 carries a blob five pitches past the end of
    # the series, reading 1795.6 where the series ends at 1290.9, and it drew a lone `+`
    # floating above that panel in every overlay. Trimmed from both ends, because the
    # reason is symmetric.
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

    Two things had to be excluded, and the caption's stated fuel levels caught both.

    `margin` removes the frame's **tick marks**: averaging a column over the full panel
    height mixes the trace with the ticks and drags every sample toward the panel centre,
    which read Figure 9's 400 lbm/hr level as 344.

    Its value is **measured, not guessed**. Scanning the inward protrusion of ink attached
    to the top and bottom frame of all twelve panels, ticks reach at most 13 px at the 99th
    percentile (the two larger readings, 26 and 21 px, are places where the trace itself
    touches a frame). It was 20 px, and that was too much: Figure 9's PCNG trace settles at
    98.79 %, which lands on row 1177.0 of a panel whose top is 1157 -- exactly the old
    window edge. The trace was clipped out of those columns and the only ink left in them
    was a GE `+` marker at 96.95, so the tracer read the marker as the curve and drew a
    notch into the last 0.4 s of every PCNG overlay. 16 px clears the ticks and keeps the
    trace.

    The single-thin-run test removes the **legend**. It sits inside a panel in the right
    half of the page, so from about t = 2 s onward the column mean was averaging the trace
    with legend text -- the level read 775.6 correctly at t = 1.6 s and 453.6 at t = 4.6 s
    on a line that never moves. A column of the bare trace is one run of a few pixels;
    anything else is not the trace.

    ## The spikes, and why they were ours

    Until 2026-09-12 this function filtered the runs *before* counting them:
    `runs = [r for r in runs if r.size <= 6]` and then `if len(runs) != 1`. That discards
    a long run instead of rejecting the column, so a column carrying BOTH a long run and
    a short one passed, and the short one won.

    That is exactly the geometry of these panels. The GE `+` markers are drawn **on** the
    model line, and each one's vertical bar has a slightly detached top. A column through
    such a bar holds a long run -- the bar merged with the line -- plus a short, separate
    run at the bar's top. The long run was thrown away and the sample was taken 40 px
    above the trace, which on the T45 panel is 61 deg R. Those were the vertical spikes in
    every overlay plot, and they were reported as a data defect of the scan (open question
    #51) when they were a defect of this loop.

    The fix is one line: count the runs **before** filtering. A column of the bare trace
    has exactly one run and it is thin; anything else -- a marker, a crossing, legend
    text, or a near-vertical segment of the trace itself -- is skipped, leaving a gap
    rather than a wrong sample.

    ## Why thinness is not enough, and continuity is

    Counting the runs fixed the spikes that sat 40 px off the trace, but not all of them:
    a glyph's detached bar-top is itself a thin single run in its column, and so is a
    scrap of legend rule. Those survived, and on Figure 9's TORQ45 and Figure 10's T45
    they put samples on the GE curve at the very end of the record, where a local
    smoothness test cannot see them either because there is no "local" left.

    The trace is a **connected** line, which is a far stronger property than thinness. So
    the candidates are split wherever adjacency breaks -- a skipped column, or a jump
    larger than the trace can make in one column -- and the blocks are then linked by
    slope, from a backbone of the long ones outward. See the constants above for why
    linking by slope and not by block length: a length threshold discarded Figure 9's
    T41 peak.
    """
    cand = []
    for x in range(p.left + margin, p.right - margin):
        # The window follows the tilted frames. A fixed window is wrong by the tilt at the
        # far end of the panel: on pdf p.46 the frames drop 11-14 px across the width, so a
        # window pinned to the left-hand frame rows leaves only 2 px of clearance under the
        # top frame at x = right, and the frame's own ink then enters the column as "trace".
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

    # Continuity. A thin single run at a plausible y is still not proof: a glyph's
    # detached bar-top presents exactly that, and so does a scrap of legend rule. The
    # trace is a CONNECTED line, so split the candidates wherever adjacency breaks --
    # either a skipped column or a jump larger than the trace can make in one column --
    # and keep only blocks long enough to be trace rather than glyph.
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

    Both frame lines of every axis here carry a printed value, so the frames pin the map
    and the *interior major ticks* become the held-out check -- the reverse of the
    appendix test, and just as strong: fit on two numbers, predict the ones you were not
    given.

    **`v_of` needs the column as well as the row.** The frames are tilted (see
    `Panel.top_slope`), so the row that means `v_hi` depends on where along the panel you
    are. Reading the value off the left-hand frame position alone put a ramp of up to 4
    percent of panel height into every trace on pdf p.46.
    """

    def x_of(px, py):
        lo = p.left_at_top + p.left_slope * (np.asarray(py, float) - p.top)
        hi = p.right_at_top + p.right_slope * (np.asarray(py, float) - p.top)
        return T_LO + (np.asarray(px, float) - lo) / (hi - lo) * (T_HI - T_LO)

    def v_of(py, px):
        dx = np.asarray(px, float) - p.left
        top = p.top + p.top_slope * dx
        bottom = p.bottom + p.bot_slope * dx
        f = (np.asarray(py, float) - top) / (bottom - top)
        return p.v_hi + f * (p.v_lo - p.v_hi)

    # The tick band follows the lower frame too, for the same reason the trace window does.
    counts = np.zeros(p.right - p.left + 1, dtype=int)
    for i, x in enumerate(range(p.left, p.right + 1)):
        yb = p.bot_at(x)
        counts[i] = int(ink[int(round(yb)) - 16 : int(round(yb)) - 4, x].sum())
    hits = np.flatnonzero(counts >= 8)
    ticks_px = []
    if hits.size:
        for g in np.split(hits, np.flatnonzero(np.diff(hits) > 3) + 1):
            ticks_px.append(float(g.mean()) + p.left)
    got = np.array([x_of(t, p.bot_at(t)) for t in ticks_px])
    err = float("nan")
    if got.size:
        want = np.arange(T_LO, T_HI + 1e-9, 1.0)
        near = [min(abs(got - w)) for w in want if np.any(np.abs(got - w) < 0.3)]
        if near:
            err = float(max(near))
    return x_of, v_of, err


def check_fuel_panel(fig: int, ts: np.ndarray, vs: np.ndarray, span: float) -> list[str]:
    """The WFPH panel against the caption. Returns a list of failures, empty if it passes.

    Two things are checked, and the second is the one that matters. The **level** of each
    plateau against the caption's stated flow catches a wrong axis range (it caught Figure
    10 being digitized on Figure 9's 250-1000 axis). The **tilt** of each plateau catches
    a wrong axis *orientation* -- a commanded constant that is not flat is the page skew,
    and nothing else on either page can reveal it.
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
        ink = native_bitmap(page, OUT / f"p{page}.pbm")
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
                ts_all = np.array([x_of(px, py) for px, py in pts])
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
                            f"= {x_of(p.left_at_top + info['pitch'], p.top) - T_LO:.4f} s\n"
                        )
                    # The plotted curve ends before the panel does, and past its end the
                    # walker re-acquires ink -- a frame line, a legend rule, a stray glyph
                    # -- and emits it as more "samples". They are not data, and they were
                    # what put a dip in every Figure 9 PCNG overlay.
                    #
                    # The separation is unambiguous, which is why this can be a rule and
                    # not a judgment: measured over all twelve panels, every artifact sits
                    # behind an x-gap of **more than 130x** the median sample interval,
                    # while every legitimate continuation of a curve sits behind a gap of
                    # 20-25x. The threshold is set at 50x, in the empty middle.
                    #
                    # Dropped, not annotated-and-kept as before. The earlier note said
                    # "consumers must window it out", and three of them did not: the
                    # overlay plots, the whole-curve metric and the shape metrics each had
                    # to grow their own filter, and one of those filters was itself wrong
                    # for a while. Provenance is preserved by recording here exactly what
                    # was removed and why -- which is what the rule below writes.
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
                            f"value {v_of(pts[cut - 1][1], pts[cut - 1][0]):.5f}; "
                            f"the dropped samples read "
                            f"{', '.join(f'{v_of(pts[i][1], pts[i][0]):.2f}' for i in dropped)}. "
                            f"Re-acquired ink past the end of the plotted curve, not data. "
                            f"See open question #48.\n"
                        )
                    fh.write("# NOT transcribed. Values carry read error.\n")
                    fh.write("t_s,value\n")
                    keep = set(range(len(pts))) - set(dropped)
                    for i, (px, py) in enumerate(pts):
                        if i in keep:
                            fh.write(f"{x_of(px, py):.5f},{v_of(py, px):.5f}\n")
            print(
                f"  {p.label:7} model {line.shape[0]:5d} pts | "
                f"reference {marks.shape[0]:3d} | x tick worst {terr:.4f} s | "
                f"tilt {p.tilt_pct_height:+.2f} %h"
            )
            if p.key == "wfph" and line.shape[0]:
                failures += check_fuel_panel(
                    fig,
                    np.array([x_of(px, py) for px, py in line]),
                    np.array([v_of(py, px) for px, py in line]),
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
