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
from PIL import Image

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


def native_bitmap(page: int, out: Path | None = None) -> np.ndarray:
    """The page's own 300 dpi 1-bit image, as an ink mask. No render, no resample.

    True is ink. `out` overrides the shared cache for a caller that wants the file
    somewhere specific; nothing about the returned array depends on it.

    The `convert("1")` is a no-op on the 1-bit PNGs `pdfimages` writes for these pages --
    PIL returns a copy when the mode already matches -- and it is here so that a page that
    ever came back as greyscale would be thresholded rather than bitwise-inverted into
    nonsense.
    """
    out = CACHE if out is None else out
    out.mkdir(parents=True, exist_ok=True)
    stem = out / f"p{page}"
    png = out / f"p{page}-000.png"
    if not png.exists():
        subprocess.run(
            ["pdfimages", "-f", str(page), "-l", str(page), "-png", str(PDF), str(stem)],
            check=True,
            capture_output=True,
        )
    return ~np.array(Image.open(png).convert("1"))


def render(page: int, dpi: int, stem: Path) -> Path:
    """A `pdftoppm` raster of the page at `dpi`, cached. **This resamples the page.**

    Only Figure A1 uses it, and not for geometry: at 600 dpi each of its eleven speed
    lines is a single 8-connected ink component, which is what separates a family of
    curves that cross and crowd. A1's axes are calibrated on `native_bitmap` like every
    other figure's.
    """
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
