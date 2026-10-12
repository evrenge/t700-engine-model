"""Shared machinery for digitizing a printed figure of TM-100991 off the source scan.

Until 2026-09-16 this machinery existed five times over. `digitize_native.py` held the
method A7 established and A2 and A9 reused; `digitize_single_curve.py` held a second copy
of the same ideas for its thirteen figures, which `digitize_multi_curve` and
`digitize_c30` imported; `digitize_a1.py` held a third; `digitize_fig678.py` and
`digitize_fig910.py` a fourth and fifth for the body figures; and `check_straightness.py`
a sixth of the geometry alone. They had drifted, and the drift was not visible from any
one file:

* the frame edges are fitted as **arcs** in one family and as **straight lines** in the
  other, and `docs/notes/digitizing-recipe.md` trap 7 then said the straight model was
  wrong (measured since, eleven of thirteen single-curve figures prefer it);
* the axis lattice is the full **minor** tick set in one family and the **interior
  majors** alone in the other, a factor of ten in the data the axis map is fitted on;
* the Nelder-Mead guards that trap 8 and trap 9 were written for reached one family and
  not the other, and f2 and f9 -- the two highest-leverage maps -- are in the family that
  did not get them.

So the point of this package is not tidiness. It is that an improvement made for one
figure should reach every figure that needs it, and while the code was duplicated there
was no mechanism by which it could.

**Both models stay reachable.** Consolidation must not move a number: `tools/reproduce_
all.sh` is a byte-identity gate, and a refactor that changes data silently is exactly what
it exists to catch. So where the two families genuinely differ, the difference is a
documented parameter here rather than a second implementation elsewhere, every caller
keeps the value it had, and choosing between them is a separate, measured decision with
the movement recorded.

Layout, one module per stage of the pipeline:

    raster.py       the page bitmap: the scan's own 300 dpi 1-bit image, never resampled
    frame.py        the four printed frame edges, and the pixel -> plot-coordinate map
    ticks.py        the printed tick lattice, and each tick's printed value
    axis.py         the axis map: order chosen by held-out test, anchored, blended
    glyph.py        the marker: a generative model of the local ink, fitted
    provenance.py   the CSV and its header -- one writer, and the parser that reads it back

This package lives in tools/ and may use SciPy and PIL. Nothing here is imported by
`src/t700/`, and `tests/test_architecture.py` enforces that.
"""
