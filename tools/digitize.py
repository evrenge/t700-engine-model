"""Digitize the printed function plots of NASA TM-100991 into CSV.

Every functional relationship in Ballin's engine model (f1..f10, f_hs) and every fuel
control schedule (F_EC1, F_HM1..F_HM7) exists ONLY as a printed plot -- see
docs/notes/inventory-appendix-a.md. This tool turns those plots into data.

The figures are kind to us: they are discrete marked data points joined by straight
segments, not freehand curves. So the job is finding marker centroids, not tracing.

Workflow per figure:

    1. grid    -- overlay a labelled pixel grid so axis ticks can be read off
    2. marks   -- detect candidate marker centroids inside a region of interest
    3. extract -- convert to data coordinates and write CSV with a provenance header
    4. verify  -- re-plot the extracted points over the source image

Step 4 is not optional. A digitization nobody checked is a guess, and CLAUDE.md is
explicit that guesses do not enter this codebase.

This module lives in tools/ and may use SciPy and Matplotlib. Nothing here is imported
by src/t700/.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

REPORT = "TM-100991"


# --------------------------------------------------------------------------- calibration


@dataclass(frozen=True)
class Axis:
    """Maps pixel coordinates to data coordinates along one axis.

    Defined by two reference points -- normally two axis ticks whose printed values
    are unambiguous. Pixel coordinates increase right/down; data values may run either
    way, which the two-point form handles without a sign convention.
    """

    px_a: float
    val_a: float
    px_b: float
    val_b: float
    log: bool = False

    def to_data(self, px: np.ndarray | float) -> np.ndarray | float:
        va, vb = (
            (np.log10(self.val_a), np.log10(self.val_b))
            if self.log
            else (
                self.val_a,
                self.val_b,
            )
        )
        if self.px_b == self.px_a:
            raise ValueError("calibration reference pixels are identical")
        t = (np.asarray(px, dtype=float) - self.px_a) / (self.px_b - self.px_a)
        v = va + t * (vb - va)
        return 10.0**v if self.log else v


@dataclass(frozen=True)
class Calibration:
    x: Axis
    y: Axis

    @staticmethod
    def parse(spec: str) -> Calibration:
        """Parse 'xpx1:xval1,xpx2:xval2,ypx1:yval1,ypx2:yval2[,logx][,logy]'.

        Pixel values come from the `grid` overlay; data values are read off the printed
        tick labels.
        """
        parts = [p.strip() for p in spec.split(",") if p.strip()]
        flags = {p for p in parts if p in ("logx", "logy")}
        pairs = [p for p in parts if p not in flags]
        if len(pairs) != 4:
            raise ValueError(
                f"calibration needs 4 point:value pairs, got {len(pairs)} -- see --help"
            )
        nums = []
        for p in pairs:
            px, _, val = p.partition(":")
            nums.append((float(px), float(val)))
        return Calibration(
            x=Axis(nums[0][0], nums[0][1], nums[1][0], nums[1][1], "logx" in flags),
            y=Axis(nums[2][0], nums[2][1], nums[3][0], nums[3][1], "logy" in flags),
        )


# --------------------------------------------------------------------------- image ops


def load_gray(path: Path) -> np.ndarray:
    """Load a rendered page as a float array in [0, 1]; 0 is ink, 1 is paper."""
    return np.asarray(Image.open(path).convert("L"), dtype=float) / 255.0


def parse_roi(spec: str | None, shape: tuple[int, int]) -> tuple[int, int, int, int]:
    """Parse 'x0,y0,x1,y1' into a clipped pixel box. None means the whole page."""
    h, w = shape
    if not spec:
        return 0, 0, w, h
    x0, y0, x1, y1 = (int(round(float(v))) for v in spec.split(","))
    return max(0, x0), max(0, y0), min(w, x1), min(h, y1)


def find_marks(
    img: np.ndarray,
    roi: tuple[int, int, int, int],
    threshold: float = 0.55,
    min_area: int = 6,
    max_area: int = 600,
    max_aspect: float = 4.0,
) -> list[tuple[float, float, int]]:
    """Find candidate data-marker centroids inside `roi`.

    Returns (x_px, y_px, area) in full-page pixel coordinates, ordered left to right.

    Markers in these figures are printed glyphs -- crosses, and in the 2-D maps, digits
    identifying the parameter line. Connected-component analysis finds them; the area and
    aspect filters reject axis lines, tick marks and stray specks. The defaults suit
    300 dpi renders and want widening at 600.
    """
    x0, y0, x1, y1 = roi
    patch = img[y0:y1, x0:x1]
    ink = patch < threshold
    labels, n = ndimage.label(ink)
    if n == 0:
        return []

    out: list[tuple[float, float, int]] = []
    objects = ndimage.find_objects(labels)
    for i, sl in enumerate(objects, start=1):
        if sl is None:
            continue
        hgt = sl[0].stop - sl[0].start
        wid = sl[1].stop - sl[1].start
        area = int((labels[sl] == i).sum())
        if not (min_area <= area <= max_area):
            continue
        if max(hgt, wid) / max(1, min(hgt, wid)) > max_aspect:
            continue  # an axis line or a tick, not a marker
        cy, cx = ndimage.center_of_mass(labels == i)
        out.append((float(cx) + x0, float(cy) + y0, area))

    out.sort(key=lambda p: p[0])
    return out


def find_marks_density(
    img: np.ndarray,
    roi: tuple[int, int, int, int],
    threshold: float = 0.55,
    window: int = 21,
    min_fill: float = 0.16,
    min_sep: int = 14,
) -> list[tuple[float, float, int]]:
    """Find markers that sit ON a connecting line, by local ink density.

    This is the mode most of these figures need. Ballin's plots join their data points
    with straight segments, so the markers are not isolated blobs -- a curve and all its
    crosses form one connected component, and `find_marks` returns a single useless
    centroid for the lot.

    A marker is locally *denser* than the line running through it. Box-filter the ink
    mask, then keep local maxima above `min_fill`. `window` should be about one marker
    across, and `min_sep` a little under the closest spacing between two real points.
    """
    x0, y0, x1, y1 = roi
    ink = (img[y0:y1, x0:x1] < threshold).astype(float)
    dens = ndimage.uniform_filter(ink, size=window, mode="constant")

    peaks = (dens == ndimage.maximum_filter(dens, size=min_sep, mode="constant")) & (
        dens >= min_fill
    )
    labels, n = ndimage.label(peaks)
    if n == 0:
        return []

    out: list[tuple[float, float, int]] = []
    for cy, cx in ndimage.center_of_mass(peaks, labels, range(1, n + 1)):
        score = int(round(dens[int(round(cy)), int(round(cx))] * window * window))
        out.append((float(cx) + x0, float(cy) + y0, score))

    out.sort(key=lambda p: p[0])
    return out


def _longest_runs(mask, axis: int):
    """Longest continuous run of True along each row (axis=1) or column (axis=0)."""
    m = mask if axis == 1 else mask.T
    out = np.zeros(m.shape[0], dtype=int)
    for i, row in enumerate(m):
        if not row.any():
            continue
        edges = np.flatnonzero(np.diff(np.concatenate(([0], row.astype(np.int8), [0]))))
        out[i] = (edges[1::2] - edges[0::2]).max()
    return out


def find_frame(
    img: np.ndarray, threshold: float = 0.55, fill: float = 0.7
) -> tuple[int, int, int, int]:
    """Locate the plot frame: (left, top, right, bottom) in page pixels.

    **The page must be rectified first** -- see `tools/rectify_page.py`. On a skewed page
    this raises rather than guessing, which is the whole point of the rewrite.

    History, because it is instructive. The first version scored a row by TOTAL ink and
    averaged every strong row in each half of the page. It was wrong nearly everywhere:

    * on pdf p.62 it returned the x-axis *title* as the bottom frame -- a block of text
      out-inks a thin line;
    * on Figure A1 it averaged the eleven flat speed lines into the "top frame", off by
      420 px;
    * on Figure A3 the flat plateaux pulled it 49 px, which is 0.009 in B1.

    Every one of those errors was invisible in the resulting CSV.

    What actually discriminates a frame line is neither total ink nor longest run (ticks
    break the lines: on A3 the top frame's longest run is 994 px and the bottom's only
    351). It is the **fraction of the plot width that is inked on that row**. The frame
    spans it; text and data do not. The vertical frame lines are solid enough that a
    longest-run scan finds them directly, and they give the width to measure against.

    Raises:
        ValueError: if fewer than two full-width lines are found. On these scans that
            almost always means the page is rotated -- a tilted line puts no single row
            above the threshold. Rectify, then call this again.
    """
    h, w = img.shape
    ink = img < threshold

    col_runs = _longest_runs(ink, axis=0)
    cols = np.flatnonzero(col_runs > 0.5 * col_runs.max())
    if cols.size == 0:
        raise ValueError("no vertical frame lines found -- is this a plot page?")
    left, right = int(cols.min()), int(cols.max())

    frac = ink[:, left : right + 1].mean(axis=1)
    idx = np.flatnonzero(frac > fill)
    if idx.size == 0:
        raise ValueError(
            "no full-width horizontal line found. The page is probably rotated -- "
            "rectify it first with tools/rectify_page.py, then call find_frame again."
        )
    groups = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    if len(groups) < 2:
        raise ValueError(
            f"only one horizontal frame line found (at y={int(groups[0].mean())}). "
            f"The page is probably rotated or the opposite frame is broken -- rectify "
            f"with tools/rectify_page.py and check the render."
        )
    top = int(round(groups[0].mean()))
    bottom = int(round(groups[-1].mean()))
    return left, top, right, bottom


# --------------------------------------------------------------------------- output


def write_csv(
    path: Path,
    points: list[tuple[float, float]],
    *,
    figure: str,
    pdf_page: int,
    quantity: str,
    x_label: str,
    y_label: str,
    note: str = "",
    calibrated: str,
) -> None:
    """Write data points with the provenance header CLAUDE.md requires.

    Every data file in this project must say where its numbers came from and how they
    were obtained. `digitized` is a weaker provenance than `transcribed`, and saying so
    in the file is the point.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    stamp = calibrated
    with path.open("w", newline="") as fh:
        for line in (
            f"# source: {REPORT} pdf p.{pdf_page}, Figure {figure}",
            f"# quantity: {quantity}",
            "# method: digitized -- marker centroids from a rendered page image",
            f"# x: {x_label}",
            f"# y: {y_label}",
            f"# points: {len(points)}",
            f"# digitized: {stamp} by tools/digitize.py",
            *([f"# note: {note}"] if note else []),
            "# NOT transcribed. Values carry read error; see the verify overlay.",
        ):
            fh.write(line + "\n")
        w = csv.writer(fh)
        w.writerow(["x", "y"])
        for x, y in points:
            w.writerow([f"{x:.6g}", f"{y:.6g}"])


def verify_overlay(
    page: Path, csv_path: Path, calib: Calibration, out: Path, roi_spec: str | None
) -> int:
    """Re-plot digitized points over the source figure. Returns the point count.

    This is the check that makes a digitization trustworthy: if a marker was missed, or
    a speck was picked up as data, or the calibration is off, it is visible here and
    nowhere else.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    img = load_gray(page)
    x0, y0, x1, y1 = parse_roi(roi_spec, img.shape)

    xs, ys = [], []
    with csv_path.open() as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            xs.append(float(row["x"]))
            ys.append(float(row["y"]))

    # invert the calibration to put data points back on the page
    def inv(axis: Axis, vals: list[float]) -> np.ndarray:
        v = np.asarray(vals, dtype=float)
        va, vb = (
            (np.log10(axis.val_a), np.log10(axis.val_b))
            if axis.log
            else (
                axis.val_a,
                axis.val_b,
            )
        )
        if axis.log:
            v = np.log10(v)
        return axis.px_a + (v - va) / (vb - va) * (axis.px_b - axis.px_a)

    fig, ax = plt.subplots(figsize=(11, 8.5), dpi=130)
    ax.imshow(img, cmap="gray", vmin=0, vmax=1)
    ax.plot(
        inv(calib.x, xs),
        inv(calib.y, ys),
        "o",
        mfc="none",
        mec="#E2231A",
        mew=1.4,
        ms=11,
        label=f"{len(xs)} extracted",
    )
    ax.add_patch(
        plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="#1F7A8C", lw=1.0, ls="--")
    )
    ax.set_title(f"{csv_path.name} over {page.name}", fontsize=10)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_xlim(0, img.shape[1])
    ax.set_ylim(img.shape[0], 0)
    ax.axis("off")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return len(xs)


def grid_overlay(page: Path, out: Path, step: int = 100) -> None:
    """Write the page with a labelled pixel grid, for reading off tick coordinates."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    img = load_gray(page)
    h, w = img.shape
    fig, ax = plt.subplots(figsize=(w / 130, h / 130), dpi=170)
    ax.imshow(img, cmap="gray", vmin=0, vmax=1)
    for x in range(0, w, step):
        ax.axvline(x, color="#1F7A8C", lw=0.4, alpha=0.65)
        ax.text(x + 2, 14, str(x), color="#E2231A", fontsize=5.5, rotation=90)
    for y in range(0, h, step):
        ax.axhline(y, color="#1F7A8C", lw=0.4, alpha=0.65)
        ax.text(2, y - 3, str(y), color="#E2231A", fontsize=5.5)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def _detect(img: np.ndarray, roi: tuple[int, int, int, int], a) -> list[tuple[float, float, int]]:
    """Dispatch to the detector the figure needs."""
    if a.mode == "density":
        return find_marks_density(img, roi, a.threshold, a.window, a.min_fill, a.min_sep)
    return find_marks(img, roi, a.threshold, a.min_area, a.max_area)


# --------------------------------------------------------------------------- cli


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("grid", help="overlay a labelled pixel grid to read tick coordinates")
    g.add_argument("--page", type=Path, required=True)
    g.add_argument("--out", type=Path, required=True)
    g.add_argument("--step", type=int, default=100)

    fr = sub.add_parser("axes", help="detect the plot frame and print a ready calibration")
    fr.add_argument("--page", type=Path, required=True)
    fr.add_argument("--x0", type=float, required=True, help="data value at the left frame")
    fr.add_argument("--x1", type=float, required=True, help="data value at the right frame")
    fr.add_argument("--y0", type=float, required=True, help="data value at the BOTTOM frame")
    fr.add_argument("--y1", type=float, required=True, help="data value at the TOP frame")
    fr.add_argument("--inset", type=int, default=12, help="shrink the suggested roi by this")

    m = sub.add_parser("marks", help="list detected marker centroids in pixel coordinates")
    m.add_argument("--page", type=Path, required=True)
    m.add_argument("--roi", help="x0,y0,x1,y1 -- restrict to the plot interior")
    m.add_argument("--threshold", type=float, default=0.55)
    m.add_argument("--min-area", type=int, default=6)
    m.add_argument("--max-area", type=int, default=600)
    m.add_argument(
        "--mode",
        choices=("blob", "density"),
        default="density",
        help="density: markers lying on a joining line (most figures). blob: isolated markers only",
    )
    m.add_argument("--window", type=int, default=21, help="density: marker width in px")
    m.add_argument("--min-fill", type=float, default=0.16, help="density: peak threshold")
    m.add_argument("--min-sep", type=int, default=14, help="density: min point spacing px")

    e = sub.add_parser("extract", help="detect, convert to data coordinates, write CSV")
    e.add_argument("--page", type=Path, required=True)
    e.add_argument("--out", type=Path, required=True)
    e.add_argument(
        "--calib", required=True, help="xpx1:xval1,xpx2:xval2,ypx1:yval1,ypx2:yval2[,logx][,logy]"
    )
    e.add_argument("--roi", help="x0,y0,x1,y1")
    e.add_argument("--figure", required=True, help="e.g. A2")
    e.add_argument("--pdf-page", type=int, required=True)
    e.add_argument("--quantity", required=True, help="e.g. f2 -- compressor temperature")
    e.add_argument("--x-label", required=True)
    e.add_argument("--y-label", required=True)
    e.add_argument("--note", default="")
    e.add_argument("--threshold", type=float, default=0.55)
    e.add_argument("--min-area", type=int, default=6)
    e.add_argument("--max-area", type=int, default=600)
    e.add_argument(
        "--mode",
        choices=("blob", "density"),
        default="density",
        help="density: markers lying on a joining line (most figures). blob: isolated markers only",
    )
    e.add_argument("--window", type=int, default=21, help="density: marker width in px")
    e.add_argument("--min-fill", type=float, default=0.16, help="density: peak threshold")
    e.add_argument("--min-sep", type=int, default=14, help="density: min point spacing px")

    v = sub.add_parser("verify", help="re-plot a CSV over its source figure")
    v.add_argument("--page", type=Path, required=True)
    v.add_argument("--csv", dest="csv_path", type=Path, required=True)
    v.add_argument("--calib", required=True)
    v.add_argument("--roi")
    v.add_argument("--out", type=Path, required=True)

    a = ap.parse_args(argv)

    if a.cmd == "axes":
        img = load_gray(a.page)
        left, top, right, bottom = find_frame(img)
        i = a.inset
        print(f"# frame of {a.page.name}: left={left} top={top} right={right} bottom={bottom}")
        print(f'--calib "{left}:{a.x0},{right}:{a.x1},{bottom}:{a.y0},{top}:{a.y1}"')
        print(f"--roi {left + i},{top + i},{right - i},{bottom - i}")
        return 0

    if a.cmd == "grid":
        grid_overlay(a.page, a.out, a.step)
        print(f"wrote {a.out}")
        return 0

    if a.cmd == "marks":
        img = load_gray(a.page)
        roi = parse_roi(a.roi, img.shape)
        marks = _detect(img, roi, a)
        print(f"# {len(marks)} candidates ({a.mode}) in roi {roi} of {a.page.name}")
        print("#      x_px      y_px  score")
        for x, y, area in marks:
            print(f"{x:10.2f}{y:10.2f}{area:7d}")
        return 0

    if a.cmd == "extract":
        img = load_gray(a.page)
        roi = parse_roi(a.roi, img.shape)
        calib = Calibration.parse(a.calib)
        marks = _detect(img, roi, a)
        pts = [(float(calib.x.to_data(x)), float(calib.y.to_data(y))) for x, y, _ in marks]
        write_csv(
            a.out,
            pts,
            figure=a.figure,
            pdf_page=a.pdf_page,
            quantity=a.quantity,
            x_label=a.x_label,
            y_label=a.y_label,
            note=a.note,
        )
        print(f"wrote {a.out} -- {len(pts)} points")
        print("NOT verified. Run `digitize.py verify` and look at the overlay.")
        return 0

    if a.cmd == "verify":
        n = verify_overlay(a.page, a.csv_path, Calibration.parse(a.calib), a.out, a.roi)
        print(f"wrote {a.out} -- {n} points plotted over {a.page.name}")
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
