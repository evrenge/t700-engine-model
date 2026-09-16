"""Is each printed curve straight between its markers? Measured in DATA space.

## Why this exists

`t700.maps` interpolates every function table linearly, so the question that decides
whether that loses anything is not how far apart the printed markers are -- it is whether
the ink between them bends. A figure whose curve is drawn as a polyline is reproduced by
`np.interp` exactly; one whose curve genuinely bows is not.

## Why DATA space, and not pixels

A first version of this check measured the bow in pixels and was wrong for two reasons
that matter on every page in this report:

* **The frames are tilted quadrilaterals.** Measured on the scans, opposite edges differ
  in angle by up to 0.36 deg (Figure A8), so the distortion is keystone as well as
  rotation. A line that is straight in data coordinates is not straight in pixels.
* **The axis maps are polynomials.** The y map departs from linear by 0.20-0.25 % of the
  axis span on these pages. Over one short inter-marker segment that is invisible -- the
  effect scales as (segment/span)^2 -- but over a segment spanning most of the plot it is
  not. Figure C29's data line spans the whole abscissa and its quadratic x map alone bows
  a data-straight line by 6 px, which a pixel-space test reads as curvature in the ink.

Pushing the traced centreline through the figure's own calibration removes tilt, keystone,
page bow and axis nonlinearity in one step, because the calibration is exactly the thing
that models them. What is left is the ink.

## What it reports

Per inter-marker segment: the rms residual from a straight line fitted in data space (the
trace's own noise floor) and the fitted sagitta (systematic bow), both in the figure's y
units and as a multiple of that file's own quoted 1 sigma.

**Read this as drafting quality, not as information loss.** These figures are hand drawn,
and a printed symbol marks a breakpoint of the underlying table: between two consecutive
symbols the draughtsman drew a straight segment, so the chord IS the data and no denser
sampling of that segment can add anything. A sagitta above 1 sigma therefore says the drawn
line wandered, or the reproduction bowed it -- never that our table is missing a point.
That distinction was got wrong on 2026-09-16 for Figure C29, whose entire content is two
symbols and the line joining them: its 4.18 sigma was written up as evidence that Ballin's
table held unmarked breakpoints, and a table cannot hold a smooth parabola. Two points
define a line; a parabola through two points is not something a plotter or a table can
produce, and the measured departure has no kinks anywhere along it, which is the signature
of a drawn-and-reproduced straight line rather than of a polyline.

The one figure where the ink genuinely carries more than the symbols is `f1`, whose speed
lines are drawn as smooth curves rather than as segments between symbols -- which is why
tracing them was right there and is not right here.

**The statistic is rms and sagitta, never the maximum residual.** On a 1-bit scan the
maximum is dominated by ink speckle: the right frame line of pdf p.62 is straight on the
printed page and scores a 10.2 px maximum against a true sagitta of 0.61 px.

## Coverage

Runs on every figure whose CSV header publishes its frame model and axis maps, which is 15
of the 19. The other four -- f1's raw table, f2, f9, F_EC1 and F_HM7 -- carry their own
stronger checks and say so in their headers (f2's 20 abscissae land on the integer lattice
to 0.0155 without the fit being told to expect integers; F_HM7 requires every segment
between consecutive knots to lie on ink, mean coverage 0.9881).

This module lives in tools/ and may use SciPy and PIL. Nothing here is imported by
src/t700/.
"""

from __future__ import annotations

import csv
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

PDF = Path("docs/ballin-tm100991.pdf")
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
WORK = ROOT / "validation" / "out" / "digitize" / "straightness"

# Figures are DISCOVERED, not listed. Two reasons. A hardcoded list is what let
# `digitize_multi_curve` sit outside the reproducibility gate, and
# `tests/test_reproducibility_gate.py` classifies a tool as owning a data file by looking
# for a literal "data/..." path in its source -- naming one here would make this
# read-only analysis tool look like a digitizer that forgot to write its output.
SKIP_STEMS = {
    "f1_compressor_mass_flow": "raw seven-point table; the beta map is what the model "
    "loads, and tools/digitize_a1.py reports its own fit residual per speed line",
    "f1_compressor_mass_flow_beta": "derived from the raw table, not traced from a figure",
    "f2_compressor_temperature": "header does not publish the frame arcs; validated by the "
    "free integer-lattice check, 0.0155 max over 20 abscissae",
    "f9_pt_mass_flow": "header does not publish the frame arcs; validated by the tool's "
    "own verification overlay",
    "fec1_thermocouple_tau": "traced multi-curve; validated by the gap-threshold sweep, "
    "<=0.05 in y",
    "fhm7_accel_limit": "vertices, not a trace; validated by segment ink coverage, mean 0.9881",
}


def figures():
    """Every committed single-curve CSV that publishes enough header to be checked."""
    out = []
    for sub in ("maps", "schedules"):
        for path in sorted((DATA / sub).glob("*.csv")):
            if path.stem in SKIP_STEMS:
                continue
            out.append(path)
    return out


def y_unit(path) -> str:
    m = re.search(r"^# y:[^,\n]*?,\s*([^\n]+)$", header(str(path)), re.M)
    return (m.group(1).strip()[:14]) if m else ""


def header(path: str) -> str:
    return "".join(line for line in open(path) if line.startswith("#"))


def _poly(expr: str, var: str) -> np.ndarray:
    """Parse a printed polynomial into numpy coefficients, highest power first.

    Has to cope with both conventions the headers use: a frame line writes its constant
    FIRST (`x = 589.107 +0.006047 y`) and an axis map writes it LAST
    (`y = -0.0067 v^3 +0.0103 v^2 +0.1364 v -0.0200`). So terms are classified by whether
    they carry the variable, never by position.
    """
    expr = expr.split("(")[0].strip()
    c: dict[int, float] = {}
    pattern = rf"([-+]?\s*[\d.]+(?:[eE][-+]?\d+)?)\s*(?:\*\s*)?({var}(?:\^(\d+))?)?"
    for coef, hasvar, power in re.findall(pattern, expr):
        coef = coef.replace(" ", "")
        if coef in ("", "+", "-"):
            continue
        try:
            value = float(coef)
        except ValueError:
            continue
        if not hasvar:
            k = 0
        elif power:
            k = int(power)
        else:
            k = 1
        c[k] = c.get(k, 0.0) + value
    if not c:
        return np.array([0.0])
    return np.array([c.get(k, 0.0) for k in range(max(c), -1, -1)])


def parse(path: str):
    """Frame model and axis maps, out of the CSV's own header. Nothing is re-derived."""
    txt = header(path)
    frame = {}
    for name in ("left", "right", "top", "bottom"):
        m = re.search(rf"^#\s+{name}\s+([xy]) = ([^\n(]+)", txt, re.M)
        if not m:
            return None
        var = "y" if m.group(1) == "x" else "x"
        frame[name] = _poly(m.group(2), var)
    axes = {}
    for axis, var in (("x", "u"), ("y", "v")):
        # A per-edge AXIS map is written in the normalized coordinate (u, v or s); a FRAME
        # line on the same edge is written in the other pixel coordinate. Classify on the
        # expression with its parenthetical comment stripped -- the comments carry words
        # like "samples" whose trailing s matched an earlier, sloppier test for the
        # variable and made every frame line look like an axis map.
        cand = []
        for edge, expr in re.findall(
            rf"^#\s+(left|right|top|bottom)\s+{axis} = ([^\n]+)", txt, re.M
        ):
            body = expr.split("(")[0]
            for sym in (var, "s"):
                if re.search(rf"(?<![A-Za-z]){sym}(?![A-Za-z])", body):
                    cand.append((edge, body, sym))
                    break
        if len(cand) == 2:
            axes[axis] = ("blend", [_poly(b, sym) for _, b, sym in cand])
        else:
            m = re.search(rf"^#\s+{axis} = ([^\n]+)", txt, re.M)
            if not m:
                return None
            axes[axis] = ("pooled", [_poly(m.group(1), var)])
    page = int(re.search(r"pdf p\.(\d+)", txt).group(1))
    return frame, axes, page


def raster(page: int) -> np.ndarray:
    WORK.mkdir(parents=True, exist_ok=True)
    stem = WORK / f"p{page}"
    if not (stem.parent / f"p{page}-000.png").exists():
        subprocess.run(
            ["pdfimages", "-f", str(page), "-l", str(page), "-png", str(PDF), str(stem)],
            check=True,
            capture_output=True,
        )
    return ~np.array(Image.open(f"{stem}-000.png").convert("1"))


def make_sv(frame):
    """Curvilinear blend between the four fitted edges, exactly as digitize_native does.

    With straight edges this reduces to the bilinear map; with curved ones it carries the
    page bow too. Either way rotation and keystone come out, which is the whole point.
    """

    def to_sv(px, py):
        px = np.asarray(px, float)
        py = np.asarray(py, float)
        xl = np.polyval(frame["left"], py)
        xr = np.polyval(frame["right"], py)
        yt = np.polyval(frame["top"], px)
        yb = np.polyval(frame["bottom"], px)
        return (px - xl) / (xr - xl), 1.0 - (py - yt) / (yb - yt)

    def to_px(s, v, iters=10):
        s = np.asarray(s, float)
        v = np.asarray(v, float)
        py = np.full(np.shape(s), 1500.0, float)
        px = np.empty_like(py)
        for _ in range(iters):
            xl = np.polyval(frame["left"], py)
            xr = np.polyval(frame["right"], py)
            px = xl + s * (xr - xl)
            yt = np.polyval(frame["top"], px)
            yb = np.polyval(frame["bottom"], px)
            py = yt + (1.0 - v) * (yb - yt)
        return px, py

    return to_sv, to_px


def axis_value(spec, t, other):
    """Evaluate one axis map at normalized coordinate `t`, blending the two edges if the
    header publishes them separately. `other` is the coordinate across the plot."""
    kind, polys = spec
    if kind == "pooled":
        return np.polyval(polys[0], t)
    lo = np.polyval(polys[0], t)
    hi = np.polyval(polys[1], t)
    return (1.0 - other) * lo + other * hi


def col_trace(mask, x0, x1, y0, tol=90.0):
    """Per-column centreline, with a predictor that keeps advancing through gaps.

    `tol` is generous and the predictor window short, because these curves steepen fast:
    Figure A3 falls 621 rows over 96 columns between two of its markers, and a tight
    tolerance makes the trace fall behind, coast on a stale slope and never re-acquire.
    """
    y = float(y0)
    out = []
    step = 1 if x1 >= x0 else -1
    slope = 0.0
    for x in range(int(x0), int(x1) + step, step):
        col = np.nonzero(mask[:, x])[0]
        if col.size == 0:
            y += slope * step
            continue
        runs = np.split(col, np.nonzero(np.diff(col) != 1)[0] + 1)
        cen = np.array([0.5 * (r[0] + r[-1]) for r in runs])
        j = int(np.argmin(np.abs(cen - y)))
        if abs(cen[j] - y) > tol:
            y += slope * step
            continue
        out.append((x, cen[j]))
        if len(out) > 2:
            rec = np.array(out[-4:])
            slope = np.polyfit(rec[:, 0], rec[:, 1], 1)[0]
        y = cen[j] + slope * step
    return np.array(out)


def curve_mask(ink, frame, band=5.0, margin=25.0):
    """Ink inside the frame with the frame lines themselves removed."""
    Y, X = np.mgrid[0 : ink.shape[0], 0 : ink.shape[1]]
    dL = X - np.polyval(frame["left"], Y)
    dR = X - np.polyval(frame["right"], Y)
    dT = Y - np.polyval(frame["top"], X)
    dB = Y - np.polyval(frame["bottom"], X)
    inside = (dL > -margin) & (dR < margin) & (dT > -margin) & (dB < margin)
    on_frame = (np.abs(dL) < band) | (np.abs(dR) < band) | (np.abs(dT) < band) | (np.abs(dB) < band)
    return ink & inside & ~on_frame


def run_one(tag, path, yunit, clear=17.0):
    parsed = parse(path)
    if parsed is None:
        print(f"{tag:>7}  header does not publish a full calibration -- skipped")
        return None
    frame, axes, page = parsed
    ink = raster(page)
    to_sv, to_px = make_sv(frame)

    rows = [r for r in csv.DictReader(ln for ln in open(path) if not ln.startswith("#"))]
    kx = "x" if "x" in rows[0] else list(rows[0])[0]
    ky = "y" if "y" in rows[0] else list(rows[0])[1]
    ksy = [k for k in rows[0] if k.startswith("sigma_y") or k == "sy"]
    X = np.array([float(r[kx]) for r in rows])
    Y = np.array([float(r[ky]) for r in rows])
    sig = float(np.mean([float(r[ksy[0]]) for r in rows])) if ksy else float("nan")

    # marker pixel columns, from the published map
    g = np.linspace(0.0, 1.0, 400001)
    mid = 0.5
    gx = axis_value(axes["x"], g, mid)
    gy = axis_value(axes["y"], g, mid)

    def s_of(x):
        return np.interp(x, gx, g) if gx[0] < gx[-1] else np.interp(x, gx[::-1], g[::-1])

    def v_of(y):
        return np.interp(y, gy, g) if gy[0] < gy[-1] else np.interp(y, gy[::-1], g[::-1])

    PX = np.array([to_px(s_of(a), v_of(b)) for a, b in zip(X, Y, strict=True)], float)

    m = curve_mask(ink, frame)
    lab, _ = ndimage.label(m, structure=np.ones((3, 3)))
    sizes = ndimage.sum(m, lab, range(1, int(lab.max()) + 1))
    # Look near EVERY marker and keep the biggest component any of them lands on. Seeding
    # from one marker is not safe: a figure's first and last markers often sit on a frame
    # line, and Figure A6's search from its own last marker found a 30 px speck outside
    # the frame instead of the 6629 px curve.
    found = {}
    rad = 26
    for mx, my in PX:
        ys, xs = int(round(my)), int(round(mx))
        win = lab[max(ys - rad, 0) : ys + rad + 1, max(xs - rad, 0) : xs + rad + 1]
        for cand in np.unique(win[win > 0]):
            found[int(cand)] = float(sizes[int(cand) - 1])
    if not found:
        print(f"{tag:>7}  could not seed the curve -- skipped")
        return None
    cid = max(found, key=found.get)
    curve = lab == cid
    cols = np.nonzero(curve.any(axis=0))[0]
    # Seed in the MIDDLE of the curve's own column range, not at a marker: the first and
    # last markers of several figures sit on a frame line, and seeding there gives the
    # trace nowhere to go in one direction and a contaminated start row in the other.
    seed_col = int(np.clip(np.median(cols), cols.min() + 3, cols.max() - 3))
    want = float(np.interp(seed_col, PX[:, 0], PX[:, 1]))
    runs_here = np.split(
        np.nonzero(curve[:, seed_col])[0],
        np.nonzero(np.diff(np.nonzero(curve[:, seed_col])[0]) != 1)[0] + 1,
    )
    seed_row = float(min(runs_here, key=lambda r: abs(0.5 * (r[0] + r[-1]) - want)).mean())
    a = col_trace(curve, seed_col, cols.max() - 2, seed_row)
    b = col_trace(curve, seed_col, cols.min() + 2, seed_row)
    tr = np.vstack([b, a]) if len(b) and len(a) else (a if len(a) else b)
    tr = tr[np.argsort(tr[:, 0])]
    _, uniq = np.unique(tr[:, 0], return_index=True)
    tr = tr[uniq]

    # THE STEP THAT MATTERS: pixels -> data, through the figure's own calibration
    s, v = to_sv(tr[:, 0], tr[:, 1])
    dx = axis_value(axes["x"], s, v)
    dy = axis_value(axes["y"], v, s)

    print(f"\n{tag}  [pdf p.{page}]   {len(X)} markers, 1 sigma_y = {sig:.3e} {yunit}")
    print(f"{'segment':>26} {'n':>6} {'rms':>11} {'sagitta':>11} {'sigma_y':>9}")
    worst = 0.0
    for i in range(len(X) - 1):
        lo, hi = sorted((PX[i][0], PX[i + 1][0]))
        keep = (tr[:, 0] >= lo + clear) & (tr[:, 0] <= hi - clear)
        if int(keep.sum()) < 6:
            print(f"{X[i]:12.4g}-{X[i + 1]:<12.4g} {int(keep.sum()):6d}   (no clear run)")
            continue
        sx, sy = dx[keep], dy[keep]
        p = np.polyfit(sx, sy, 1)
        res = sy - np.polyval(p, sx)
        q = np.polyfit(sx, sy, 2)
        sag = abs(q[0]) * (sx.max() - sx.min()) ** 2 / 8.0
        worst = max(worst, sag / sig if sig == sig else 0.0)
        print(
            f"{X[i]:12.4g}-{X[i + 1]:<12.4g} {int(keep.sum()):6d} "
            f"{np.sqrt((res**2).mean()):11.3e} {sag:11.3e} {sag / sig:9.2f}"
        )
    print(f"{'worst bow':>26} {'':6} {'':11} {'':11} {worst:9.2f}")
    return worst


def main():
    print("STRAIGHTNESS OF EVERY PRINTED CURVE, MEASURED IN DATA SPACE")
    print("sagitta is the systematic bow; judge it against the file's own 1 sigma on y")
    worst = {}
    for path in figures():
        tag = path.stem.split("_")[0]
        w = run_one(tag, str(path), y_unit(path))
        if w is not None:
            worst[tag] = w
    print("\n" + "=" * 70)
    print("SUMMARY -- worst systematic bow per figure, in units of its own 1 sigma on y")
    for tag, w in sorted(worst.items(), key=lambda kv: -kv[1]):
        flag = "   <-- above the read error" if w > 1.0 else ""
        print(f"  {tag:>7} {w:8.2f}{flag}")
    print("\nNOT COVERED (each carries a stronger check of its own):")
    for stem, why in sorted(SKIP_STEMS.items()):
        print(f"  {stem:>30}  {why}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
