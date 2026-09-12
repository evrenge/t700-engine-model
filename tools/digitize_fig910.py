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
# the caption states the fuel step exactly, so WFPH validates its own calibration
FUEL_STEP = {9: (400.0, 775.0), 10: (400.0, 125.0)}
N_PANELS = 6
T_LO, T_HI = 0.0, 5.0


@dataclass
class Panel:
    key: str
    label: str
    top: int
    bottom: int
    left: int
    right: int
    v_lo: float
    v_hi: float


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
        out.append(Panel(key, label, int(round(t)), int(round(t + height)), left, right, lo, hi))
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
    lo, hi = p.top + margin, p.bottom - margin
    x0, x1 = p.left + margin, p.right - margin
    sub = ink[lo:hi, x0:x1]
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
        col = ink[p.top + margin : p.bottom - margin, x]
        ys = np.flatnonzero(col)
        if ys.size == 0:
            continue
        runs = np.split(ys, np.flatnonzero(np.diff(ys) > 1) + 1)
        if len(runs) != 1 or runs[0].size > 6:
            continue  # legend text, a marker, a crossing, or a near-vertical segment
        y = float(runs[0].mean()) + p.top + margin
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
    """

    def x_of(px):
        return T_LO + (np.asarray(px, float) - p.left) / (p.right - p.left) * (T_HI - T_LO)

    def v_of(py):
        f = (np.asarray(py, float) - p.top) / (p.bottom - p.top)
        return p.v_hi + f * (p.v_lo - p.v_hi)

    band = ink[p.bottom - 16 : p.bottom - 4, p.left : p.right + 1]
    hits = np.flatnonzero(band.sum(axis=0) >= 8)
    ticks_px = []
    if hits.size:
        for g in np.split(hits, np.flatnonzero(np.diff(hits) > 3) + 1):
            ticks_px.append(float(g.mean()) + p.left)
    got = np.array([x_of(t) for t in ticks_px])
    err = float("nan")
    if got.size:
        want = np.arange(T_LO, T_HI + 1e-9, 1.0)
        near = [min(abs(got - w)) for w in want if np.any(np.abs(got - w) < 0.3)]
        if near:
            err = float(max(near))
    return x_of, v_of, err


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    REF.mkdir(parents=True, exist_ok=True)
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
                ts_all = np.array([x_of(px) for px, _ in pts])
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
                    fh.write("# t_s: time, seconds (shared axis, 0 to 5.0)\n")
                    fh.write(
                        f"# value: {p.label}, axis range {p.v_lo} to {p.v_hi}, "
                        f"no unit printed on the figure\n"
                    )
                    fh.write(f"# points: {pts.shape[0] - len(dropped)}\n")
                    fh.write(f"# x tick check: worst interior major tick off by {terr:.4f} s\n")
                    if kind == "reference":
                        fh.write(
                            f"# sample lattice: pitch {info['pitch']:.2f} px "
                            f"= {x_of(p.left + info['pitch']) - T_LO:.4f} s\n"
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
                            f"value {v_of(pts[cut - 1][1]):.5f}; the dropped samples read "
                            f"{', '.join(f'{v_of(pts[i][1]):.2f}' for i in dropped)}. "
                            f"Re-acquired ink past the end of the plotted curve, not data. "
                            f"See open question #48.\n"
                        )
                    fh.write("# NOT transcribed. Values carry read error.\n")
                    fh.write("t_s,value\n")
                    keep = set(range(len(pts))) - set(dropped)
                    for i, (px, py) in enumerate(pts):
                        if i in keep:
                            fh.write(f"{x_of(px):.5f},{v_of(py):.5f}\n")
            print(
                f"  {p.label:7} model {line.shape[0]:5d} pts | "
                f"reference {marks.shape[0]:3d} | x tick worst {terr:.4f} s"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
