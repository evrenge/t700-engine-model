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


def marker_centres(ink: np.ndarray, p: Panel, bar: int = 7, margin: int = 20) -> np.ndarray:
    """Find `+` glyph centres inside a panel.

    A `+` is the only thing on the page with **both** a horizontal stroke and a vertical
    one at the same place. That is the whole discriminator, and it is what makes this
    robust:

    * the model trace where it runs flat has horizontal extent but no vertical extent;
    * the model trace where it runs steep has vertical extent but no horizontal extent;
    * a `+` has both.

    A first attempt used the horizontal opening alone and returned 146 "markers" in the
    WFPH panel -- which carries no markers at all, only a step whose two flat levels are
    long horizontal runs. Requiring both strokes removes that entire class of false
    positive.

    `margin` then removes the other class. Every frame carries **tick marks**, and a tick
    crossing its own frame line presents a vertical stroke and a horizontal stroke at the
    same place -- indistinguishable from a glyph by the test above. The ticks protrude
    about 10 px at 300 dpi, so the interior is taken 20 px inside every frame.
    """
    sub = ink[p.top + margin : p.bottom - margin, p.left + margin : p.right - margin]
    horiz = ndimage.binary_opening(sub, structure=np.ones((1, bar)))
    vert = ndimage.binary_opening(sub, structure=np.ones((bar, 1)))
    both = ndimage.binary_dilation(horiz & vert, structure=np.ones((3, 3)))
    lab, n = ndimage.label(both, structure=np.ones((3, 3)))
    if n == 0:
        return np.zeros((0, 2))
    pts = []
    for i, sl in enumerate(ndimage.find_objects(lab), start=1):
        hgt = sl[0].stop - sl[0].start
        wid = sl[1].stop - sl[1].start
        if wid > 3 * bar or hgt > 3 * bar:
            continue  # a crossing of two long strokes, not a glyph
        blob = lab[sl] == i
        cy, cx = ndimage.center_of_mass(blob)
        pts.append((cx + sl[1].start + p.left + margin, cy + sl[0].start + p.top + margin))
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
    return out, {"slots": len(out), "pitch": pitch, "dropped": len(pts) - len(out)}


def trace_line(
    ink: np.ndarray, p: Panel, marks: np.ndarray, clear: float = 9.0, margin: int = 20
) -> np.ndarray:
    """Follow the solid model trace, one sample per column, markers masked out.

    Two things had to be excluded, and the caption's stated fuel levels caught both.

    `margin` removes the frame's **tick marks**: averaging a column over the full panel
    height mixes the trace with the ticks and drags every sample toward the panel centre,
    which read Figure 9's 400 lbm/hr level as 344.

    The single-thin-run test removes the **legend**. It sits inside a panel in the right
    half of the page, so from about t = 2 s onward the column mean was averaging the trace
    with legend text -- the level read 775.6 correctly at t = 1.6 s and 453.6 at t = 4.6 s
    on a line that never moves. A column of the bare trace is one run of a few pixels;
    anything else is not the trace.
    """
    pts = []
    for x in range(p.left + margin, p.right - margin):
        col = ink[p.top + margin : p.bottom - margin, x]
        ys = np.flatnonzero(col)
        if ys.size == 0:
            continue
        runs = np.split(ys, np.flatnonzero(np.diff(ys) > 1) + 1)
        runs = [r for r in runs if r.size <= 6]
        if len(runs) != 1:
            continue  # legend text, a marker, or a crossing -- not the bare trace
        y = float(runs[0].mean()) + p.top + margin
        if marks.size:
            near = marks[np.abs(marks[:, 0] - x) < clear]
            if near.size and np.min(np.abs(near[:, 1] - y)) < clear:
                continue
        pts.append((x, y))
    return np.array(pts) if pts else np.zeros((0, 2))


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
                    fh.write(f"# points: {pts.shape[0]}\n")
                    fh.write(f"# x tick check: worst interior major tick off by {terr:.4f} s\n")
                    if kind == "reference":
                        fh.write(
                            f"# sample lattice: pitch {info['pitch']:.2f} px "
                            f"= {x_of(p.left + info['pitch']) - T_LO:.4f} s\n"
                        )
                    fh.write("# NOT transcribed. Values carry read error.\n")
                    fh.write("t_s,value\n")
                    for px, py in pts:
                        fh.write(f"{x_of(px):.5f},{v_of(py):.5f}\n")
            print(
                f"  {p.label:7} model {line.shape[0]:5d} pts | "
                f"reference {marks.shape[0]:3d} | x tick worst {terr:.4f} s"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
