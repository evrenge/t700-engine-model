"""The page raster. One cache, one answer, no resampling.

Every digitizer starts by turning a page of the scan into an ink mask, and until
2026-09-16 four of them carried their own copy of the six lines that do it -- each with
its own cache directory, so the same page was rastered up to four times and the four
copies could in principle disagree about what ink is.

## Why `pdfimages` and not `pdftoppm`

`pdfimages -list` shows the figure pages of TM-100991 are each a single 300 dpi 1-bit
CCITT image. That bitmap *is* the document: a `pdftoppm` render at any dpi is a resample
of it, and a resample at 600 dpi is an upsample that invents intermediate greys. Sub-pixel
centroids are fine on an upsample; geometry is not. So the rule the whole package follows
is that nothing resamples the page -- rotation, keystone and page bow are carried as a
*map applied to the measured points*, never as a warp applied to the pixels.

`render()` exists for the one tool that needs a higher-dpi raster for a different reason
(Figure A1, whose eleven speed lines have to be separated as connected components), and it
says so at the call site.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

PDF = Path("docs/ballin-tm100991.pdf")

CACHE = Path("validation/out/digitize/pages")
"""One cache for every tool, keyed by page.

Previously `digitize_native` wrote `<caller's dir>/native-000.png`, `digitize_single_curve`
wrote `a6810/p{page}-000.png`, `digitize_a1` wrote its own, and `check_straightness` a
fourth. Sharing one directory means a page is rastered once per clone rather than once per
tool. It is under the gitignored `validation/out/`, and it is produced here rather than by
hand, which is what `tests/test_reproducibility_gate.py` requires of any input a tool reads
from there.
"""


def native_bitmap(page: int) -> np.ndarray:
    """The page's own 300 dpi 1-bit image, as an ink mask. No render, no resample.

    True is ink.

    The `convert("1")` is a no-op on the 1-bit PNGs `pdfimages` writes for these pages --
    PIL returns a copy when the mode already matches -- and it is here so that a page that
    ever came back as greyscale would be thresholded rather than bitwise-inverted into
    nonsense.
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

    Only Figure A1 uses it, and not for geometry: at 600 dpi each of its eleven speed
    lines is a single 8-connected ink component, which is what separates a family of
    curves that cross and crowd. A1's axes are calibrated on `native_bitmap` like every
    other figure's.
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
    `(box, size, path)`, each a region of the page saved at `size`. A picture to check by
    eye, not data: nothing reads it back, which is why it may resize.
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
