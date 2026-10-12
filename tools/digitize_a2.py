"""Re-extract Figure A2 (`f2`, compressor temperature ratio) from TM-100991 pdf p.57.

Why this figure was redone.  `f2` is one of the two highest-leverage maps in the engine:
a sensitivity sweep of the trimmed model gives dSHP/df2 of -2.0 % to -7.7 % per +1 %, so
a 0.2 % calibration error is worth up to 1.5 % of shaft power.  The 2026-09-10 extraction
used a per-axis *linear* map on a perspective-rectified render, and A7 later showed that a
linear axis map cannot fit these scans.

Method, following `tools/digitize_a7.py`:

1.  The SCAN, not a render.  `pdfimages -list` says p.57 is one 300 dpi 1-bit CCITT image,
    so any higher-dpi render is an upsample.  Nothing is ever resampled.
2.  Geometry carried as a map applied to the measured points.  The four frame edges are
    fitted as *arcs*: they are printed straight and image with a -1.2 px (left) and
    -1.4 px (right) sagitta, so a straight-edge model -- a homography included -- is wrong
    by that much in mid-page.  Plot coordinates are the curvilinear blend between the arcs.
3.  Axis order chosen by a held-out test.  A2 prints a value at all four frame edges
    (0 / 20 in x, 1.1 / 3.0 in y).  Fit the interior ticks alone and predict the printed
    edge; a second test that uses no printed edge at all -- fit the middle 64 % of the
    lattice, predict the two ends -- has to agree.
4.  Markers by a generative glyph model: two strokes plus the two joining half-segments,
    all radiating from one point, supersampled, scored on ink AND background.

Run from the repo root, with `tools/` on PYTHONPATH.  Writes overlays under
validation/out/digitize/a2redo/ and prints the rows of
data/maps/f2_compressor_temperature.csv.  Needs `pdfimages` (poppler).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from figure.axis import axis_values, blend_axis, held_out_table, point_sigma
from figure.frame import (
    curve_component,
    fit_edges,
    make_sv,
    print_segments,
    segment_lines,
)
from figure.glyph import fit_on_polyline, seeds_from_lattice
from figure.provenance import calibration_lines, check_against_csv, check_calibration
from figure.raster import native_bitmap, overlay
from figure.ticks import native_lattice

OUT = Path("validation/out/digitize/a2redo")
CSV = Path("data/maps/f2_compressor_temperature.csv")
PAGE = 57
N = 20

SPEC = [
    ("left", 0, 560, 625, 320, 2610),
    ("right", 0, 2225, 2290, 320, 2610),
    ("top", 1, 265, 325, 620, 2230),
    ("bottom", 1, 2605, 2670, 620, 2230),
]

FRAME_ORDER, X_ORDER, Y_ORDER = 2, 3, 3  # chosen by the tests main() prints

VARIANTS = [
    (fo, xo, yo, bl) for fo in (1, 2) for xo in (1, 3) for yo in (1, 3, 4) for bl in (True, False)
]


# (along, sign, lo, hi, long-tick threshold, -, coarse step, -): A2 uses its long ticks alone
EDGES = {
    "left": (0, +1, 300, 2635, 7, None, 0.1, None),
    "right": (0, -1, 300, 2645, 7, None, 0.1, None),
    "bottom": (1, -1, 590, 2260, 3, None, 1.0, None),
    "top": (1, +1, 590, 2260, 3, None, 1.0, None),
}


def tickset(fc, to_sv, verbose=False):
    """Every printed tick on the four edges, labelled, frame-edge ticks removed."""
    return native_lattice(
        tickset.ink,
        fc,
        to_sv,
        EDGES,
        {"TL": (589, 291), "TR": (2255, 296), "BL": (593, 2631), "BR": (2257, 2641)},
        (0.0, 20.0),
        (1.1, 3.0),
        fine=False,
        tol=0.30,
        verbose=verbose,
    )


def report_models(P):
    for axis, lo, hi, unit, uname, edges in (
        ("y = T3/T2", 1.1, 3.0, 1e-3, " mU", ("left", "right")),
        ("x = PS3/P2", 0.0, 20.0, 1.0, "", ("bottom", "top")),
    ):
        for nm in edges:
            t, val = P[nm]
            print(f"\n  --- {axis}, {nm} edge ---")
            print(f"   frame-edge held out (printed {lo} at 0 and {hi} at 1, neither in the fit);")
            print("   tail CV: fit the middle 64 % of the lattice, predict the two ends")
            held_out_table(t, val, lo, hi, unit, uname, orders=(1, 2, 3, 4))


def glyph_pass(ink, fc, seeds, lines, mask_line=False):
    """Fit the 20 glyphs.  `mask_line` is the independent check: hide the joining polyline
    from the score entirely instead of modelling it, and see whether the centres move."""

    def frame_mask(X, Y):
        off_frame = np.abs(X - np.polyval(fc["right"], Y)) > 5.0
        off_tick = ~((Y >= 536) & (Y <= 552) & (X <= 2256))  # the 2.8 major under marker 20
        return off_frame & off_tick

    return fit_on_polyline(
        ink,
        seeds,
        lines,
        [np.radians(-45.0), 15.0, np.radians(45.0), 15.0, 3.3, 3.3],
        42.0,
        # `frame_mask` marks what to KEEP; `windows` takes what to EXCLUDE.
        {N - 1: lambda X, Y: ~frame_mask(X, Y)},
        mask_line=mask_line,
        # Trap 8's guard, adopted 2026-09-16 after measuring it: an explicit 2 px simplex
        # on the OFFSET, and a refusal of any centre more than 10 px from its seed. A2's
        # own free check -- its markers sit at integer PS3/P2 and the tool is never told
        # so -- cannot separate the two fits: worst residual 0.0224 -> 0.0218 at order 1
        # and 0.0155 -> 0.0151 at the chosen order 3, sd 0.0060 -> 0.0062. A wash, which
        # is the honest verdict. The reason to adopt it anyway is that the unguarded fit
        # depends on SciPy's default simplex not walking off the page, and on this figure
        # the markers run out to column 2256, where that first step is 113 px.
        max_move=10.0,
    )


def main() -> int:
    ink = native_bitmap(PAGE)
    tickset.ink = ink
    print("frame arcs (printed straight; the sagitta is the page's bow):")
    fc, _ = fit_edges(ink, SPEC, order=FRAME_ORDER)
    to_sv, to_px = make_sv(fc)

    print("\ncalibration ticks:")
    P = tickset(fc, to_sv, verbose=True)
    print("\nCALIBRATION MODEL SELECTION")
    report_models(P)
    print(f"\n  chosen: frame arcs order {FRAME_ORDER}, x order {X_ORDER}, y order {Y_ORDER};")
    print("  each edge fitted alone with its two printed frame values, then blended.")

    cur = curve_component(ink, fc)
    # marker 20 sits ON the right frame, so it is not in `cur`
    seeds = seeds_from_lattice(
        cur, to_sv, [k / 20.0 for k in range(1, N + 1)], {N - 1: (2257.0, 533.0)}
    )
    ys, xs = np.nonzero(cur)
    lines = segment_lines(np.stack([xs, ys], 1).astype(float), seeds)
    print("\njoining segments:")
    print_segments(lines)

    print("\nglyph fit:")
    C, _ = glyph_pass(ink, fc, seeds, lines)
    print("\nglyph fit with the joining polyline masked out entirely (independent check):")
    C2, _ = glyph_pass(ink, fc, seeds, lines, mask_line=True)
    dd = np.hypot(*(C - C2).T)
    print(f"  agreement max {dd.max():.2f} px, rms {dd.std():.2f} px:  {np.round(dd, 2)}")

    py = blend_axis(P["left"][0], P["left"][1], P["right"][0], P["right"][1], Y_ORDER, 1.1, 3.0)
    px = blend_axis(P["bottom"][0], P["bottom"][1], P["top"][0], P["top"][1], X_ORDER, 0.0, 20.0)
    s, v = to_sv(C[:, 0], C[:, 1])
    X = (1 - v) * np.polyval(px[0], s) + v * np.polyval(px[1], s)
    Y = (1 - s) * np.polyval(py[0], v) + s * np.polyval(py[1], v)

    print("\nfree structural check -- the abscissae against the integer lattice (never told):")
    for xo in (1, 2, 3):
        Xv, _ = axis_values(ink, SPEC, tickset, C, FRAME_ORDER, xo, Y_ORDER)
        e = Xv - np.arange(1, N + 1)
        print(
            f"  x order {xo}: max {np.abs(e).max():.4f} ({np.abs(e).max() * 83.2:.2f} px)  "
            f"sd {e.std():.4f} ({e.std() * 83.2:.2f} px)  mean {e.mean():+.4f}"
        )
    print(f"  chosen model: {np.round(X - np.arange(1, N + 1), 4)}")

    print("\nmodel spread -- every defensible calibration, same glyph centres:")
    var = [axis_values(ink, SPEC, tickset, C, *m) for m in VARIANTS]
    VX = np.array([a for a, _ in var])
    VY = np.array([b for _, b in var])
    print(f"  all {len(var)}: x sd max {VX.std(0).max():.4f}   y sd max {VY.std(0).max():.5f}")
    ok = [i for i, m in enumerate(VARIANTS) if m[2] in (3, 4)]
    mx, my = VX[ok].std(0), VY[ok].std(0)
    print(f"  y-order 3/4 subset ({len(ok)}): x sd max {mx.max():.4f}   y sd max {my.max():.5f}")

    # ---- per-point 1 sigma, by quadrature
    glyph_px = np.maximum(dd, 0.35)  # polyline-masked disagreement, floored at the fit noise
    unc = point_sigma(P, fc, C, s, v, px, py, X_ORDER, Y_ORDER, glyph_px, mx, my)
    sx, sy, cal_x, cal_y = unc["sx"], unc["sy"], unc["cal_x"], unc["cal_y"]
    dXdpx, dYdpx = unc["dXdpx"], unc["dYdpx"]

    overlay(
        ink,
        to_px,
        C,
        ((px[0], range(21)), (py[0], [1.1 + 0.1 * k for k in range(20)])),
        [
            ((520, 240, 2330, 2700), (905, 1230), OUT / "a2_verify.png"),
            ((2120, 420, 2320, 700), (800, 1120), OUT / "a2_verify_m20.png"),
            ((600, 2280, 900, 2460), (1200, 720), OUT / "a2_verify_m1.png"),
        ],
    )

    print(
        "\n   k        x           y      1s_x    1s_y   | cal_x  glyph_x model_x |"
        "  cal_y  glyph_y model_y"
    )
    for i in range(N):
        print(
            f"  {i + 1:2d} {X[i]:10.5f} {Y[i]:9.5f} {sx[i]:7.5f} {sy[i]:7.5f} | "
            f"{cal_x[i]:6.4f} {glyph_px[i] * dXdpx[i]:7.4f} {mx[i]:7.4f} | "
            f"{cal_y[i]:6.5f} {glyph_px[i] * dYdpx[i]:7.5f} {my[i]:7.5f}"
        )
    # The CSV is maintained by hand; these values are not. See `check_against_csv`.
    rows = check_against_csv(
        CSV,
        {
            "x": (X, "%.5f"),
            "y": (Y, "%.5f"),
            "sx": (sx, "%.5f"),
            "sy": (sy, "%.5f"),
        },
        label="f2 (Figure A2)",
    )
    # Nor is the calibration its header publishes, which `check_straightness` reads back.
    cal = check_calibration(CSV, calibration_lines(fc, px, py), label="f2 (Figure A2)")
    return max(rows, cal)


if __name__ == "__main__":
    raise SystemExit(main())
