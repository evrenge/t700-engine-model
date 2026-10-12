"""Shared machinery for digitizing a printed figure of TM-100991 off the source scan.

One module per stage of the pipeline; each figure script in `tools/` supplies only that
figure's own facts. `docs/notes/digitizing-recipe.md` holds the traps the docstrings here
cite by number.

    raster.py       the page bitmap: the scan's own 300 dpi 1-bit image, never resampled
    frame.py        the four printed frame edges, and the pixel -> plot-coordinate map
    ticks.py        the printed tick lattice, and each tick's printed value
    axis.py         the axis map: order chosen by held-out test, anchored, blended
    glyph.py        the marker: a generative model of the local ink, fitted
    provenance.py   the CSV and its header -- one writer, and the parser that reads it back

Where the callers still differ (frame order, anchoring rule, glyph tolerances) the difference
is a parameter here, not a second implementation: unifying it would move published numbers,
which `tools/reproduce_all.sh` checks byte for byte. Nothing here is imported by `src/t700/`.
"""
