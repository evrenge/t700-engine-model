"""The page raster: one cache, and no resampling.

`pdfimages -list` shows each figure page of TM-100991 is a single 300 dpi 1-bit CCITT image.
That bitmap is the document; a `pdftoppm` render is a resample of it and invents greys. So
nothing resamples the page: rotation, keystone and bow are carried as a map applied to the
measured points (recipe, trap 7). `render()` is the one exception, for Figure A1.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

PDF = Path("docs/ballin-tm100991.pdf")

CACHE = Path("validation/out/digitize/pages")
"""One page cache for every tool. It is under the gitignored `validation/out/` and filled
here, as `tests/test_reproducibility_gate.py` requires of any input read from there."""


def native_bitmap(page: int) -> np.ndarray:
    """The page's own 300 dpi 1-bit image as an ink mask, True is ink. No resample.

    `convert("1")` is a no-op on the 1-bit PNGs `pdfimages` writes; it is there so that a
    page that came back greyscale is thresholded rather than bitwise-inverted into nonsense.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    stem = CACHE / f"p{page}"
    png = CACHE / f"p{page}-000.png"
    if not png.exists():
        subprocess.run(
            ["pdfimages", "-f", str(page), "-l", str(page), "-png", str(PDF), str(stem)],
            check=True,
            capture_output=True,
        )
    return ~np.array(Image.open(png).convert("1"))


def render(page: int, dpi: int, stem: Path) -> Path:
    """A `pdftoppm` raster of the page at `dpi`, cached by page and dpi. **This resamples.**

    Only Figure A1 uses it, and not for geometry: at 600 dpi each of its eleven speed lines
    is one 8-connected ink component. A1's axes are calibrated on `native_bitmap`.
    """
    stem = stem.parent / f"{stem.name}-{dpi}dpi"
    stem.parent.mkdir(parents=True, exist_ok=True)
    png = stem.parent / f"{stem.name}-{page:03d}.png"
    if not png.exists():
        subprocess.run(
            [
                "pdftoppm",
                "-f",
                str(page),
                "-l",
                str(page),
                "-r",
                str(dpi),
                "-png",
                str(PDF),
                str(stem),
            ],
            check=True,
        )
    return png


def load_gray(path: Path) -> np.ndarray:
    """A rendered page as a float array in [0, 1]; 0 is ink, 1 is paper."""
    return np.asarray(Image.open(path).convert("L"), dtype=float) / 255.0


def overlay(ink, to_px, C, grid, crops, n=20001):
    """The calibration gridlines drawn back onto the scan through `to_px`, markers ringed.

    `grid` is `((x_map, x_values), (y_map, y_values))`: each axis map as a polynomial in the
    unit-square coordinate, and the printed values to draw a line at. `crops` is a list of
    `(box, size, path)`. A picture to check by eye; nothing reads it back, so it may resize.
    """
    im = Image.fromarray((~ink * 255).astype(np.uint8)).convert("RGB")
    d = ImageDraw.Draw(im)
    tt = np.linspace(0.0, 1.0, n)
    for (poly, values), horiz in zip(grid, (False, True), strict=True):
        for val in values:
            s = float(tt[np.argmin(np.abs(np.polyval(poly, tt) - val))])
            p0, p1 = ((0.0, s), (1.0, s)) if horiz else ((s, 0.0), (s, 1.0))
            q0, q1 = to_px(*p0), to_px(*p1)
            d.line(
                [(float(q0[0][0]), float(q0[1][0])), (float(q1[0][0]), float(q1[1][0]))],
                fill=(210, 230, 255),
                width=1,
            )
    for x, y in C:
        d.ellipse([x - 26, y - 26, x + 26, y + 26], outline=(220, 0, 0), width=3)
        d.line([(x - 40, y), (x + 40, y)], fill=(0, 170, 0), width=1)
        d.line([(x, y - 40), (x, y + 40)], fill=(0, 170, 0), width=1)
    for box, size, path in crops:
        im.crop(box).resize(size, Image.LANCZOS).save(path)
