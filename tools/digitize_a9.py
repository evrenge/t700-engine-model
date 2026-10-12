"""Re-extract Figure A9 (`f9`, power turbine corrected mass flow) from TM-100991 pdf p.64.

The CSV is hand-maintained: this tool recomputes its rows and its calibration header and
checks them, returning non-zero on any difference. Method as in `tools/digitize_a7.py`.
What is specific to this page:

*   `x` prints a value only at the RIGHT frame (0.85). The left frame's 0.30 is inferred
    from the tick count; the cubic tick fits, told nothing about it, predict it to within
    two pixels.
*   Both `y` frame edges are printed (0.25, 0.39), so `y` gets the full two-ended
    held-out test.
*   The joining polyline is 1 px wide over the flat run and the scan has broken it in a
    dozen places, so the curve is found by closing before labelling.
*   The first marker sits ON the left frame and the last ON the right frame, each with
    minor ticks underneath; both are masked, never modelled.

Run from the repo root, with `tools/` on PYTHONPATH.  Writes overlays under
validation/out/digitize/a9redo/ and prints the rows of data/maps/f9_pt_mass_flow.csv.
Needs `pdfimages` (poppler).
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

OUT = Path("validation/out/digitize/a9redo")
CSV = Path("data/maps/f9_pt_mass_flow.csv")
PAGE = 64
N = 23

SPEC = [
    ("left", 0, 435, 505, 340, 2630),
    ("right", 0, 2095, 2165, 340, 2630),
    ("top", 1, 285, 350, 500, 2100),
    ("bottom", 1, 2620, 2690, 500, 2100),
]

X0, X1, Y0, Y1 = 0.30, 0.85, 0.25, 0.39
FRAME_ORDER, X_ORDER, Y_ORDER = 2, 3, 3
FINE = True  # calibrate on every tick (0.01 in x, 0.001 in y), not on the long ticks alone

VARIANTS = [
    (fo, xo, yo, bl) for fo in (1, 2) for xo in (1, 3) for yo in (1, 3, 4) for bl in (True, False)
]

# (edge, along, sign, lo, hi, long-tick threshold, all-tick threshold, coarse step, fine step)
EDGES = {
    "left": (0, +1, 300, 2660, 6, 2, 0.01, 0.001),
    "right": (0, -1, 300, 2670, 6, 2, 0.01, 0.001),
    "bottom": (1, -1, 465, 2135, 6, 2, 0.05, 0.01),
    "top": (1, +1, 465, 2135, 6, 2, 0.05, 0.01),
}


def tickset(fc, to_sv, verbose=False, fine=None):
    """This page's tick lattice: `figure.ticks.native_lattice` with its edges and corners."""
    return native_lattice(
        tickset.ink,
        fc,
        to_sv,
        EDGES,
        {"TL": (469, 311), "TR": (2130, 321), "BL": (468, 2645), "BR": (2125, 2662)},
        (X0, X1),
        (Y0, Y1),
        fine=FINE if fine is None else fine,
        verbose=verbose,
    )


def report_models(P):
    for axis, lo, hi, edges in (
        ("y = W45C", Y0, Y1, ("left", "right")),
        ("x = Ps9/P45", None, X1, ("bottom", "top")),
    ):
        for nm in edges:
            t, val = P[nm]
            print(f"\n  --- {axis}, {nm} edge ---")
            lab = f"printed {hi} at 1" + ("" if lo is None else f", printed {lo} at 0")
            print(
                f"   frame-edge held out ({lab}; not in the fit).  edge0 unprinted = a prediction"
            )
            print("   tail CV: fit the middle 64 % of the lattice, predict the two ends")
            held_out_table(t, val, lo, hi, 1e-5, "e-5", orders=(1, 2, 3, 4))


def glyph_pass(ink, fc, seeds, lines, mask_line=False):
    def left_mask(X, Y):
        return np.abs(X - np.polyval(fc["left"], Y)) > 5.0

    def right_mask(X, Y):
        return np.abs(X - np.polyval(fc["right"], Y)) > 5.0

    return fit_on_polyline(
        ink,
        seeds,
        lines,
        [np.radians(-45.0), 15.0, np.radians(45.0), 15.0, 3.6, 2.4],
        32.0,
        # The masks mark what to KEEP; `windows` takes what to EXCLUDE.
        {0: lambda X, Y: ~left_mask(X, Y), N - 1: lambda X, Y: ~right_mask(X, Y)},
        mask_line=mask_line,
        # Recipe traps 8-9: an explicit 2 px simplex on the offset, and no centre more than
        # 10 px from its seed.
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
    print("\nCALIBRATION MODEL SELECTION (all ticks: 0.01 in x, 0.001 in y)")
    report_models(P)
    print("\n  the same test on the LONG ticks alone (0.05 in x, 0.01 in y):")
    report_models(tickset(fc, to_sv, fine=False))
    print(f"\n  chosen: frame arcs order {FRAME_ORDER}, x order {X_ORDER}, y order {Y_ORDER};")
    print("  each edge fitted alone with its printed frame values, then blended.")

    cur = curve_component(ink, fc, close=3)
    # the first and last markers sit ON the left and right frames, so are not in `cur`
    seeds = seeds_from_lattice(
        cur,
        to_sv,
        [k / (N - 1.0) for k in range(N)],
        {0: (469.4, 612.6), N - 1: (2126.0, 2468.0)},
    )
    ys, xs = np.nonzero(cur)
    lines = segment_lines(np.stack([xs, ys], 1).astype(float), seeds, clear=20.0, band=4.0)
    print("\njoining segments:")
    print_segments(lines)

    print("\nglyph fit:")
    C, _ = glyph_pass(ink, fc, seeds, lines)
    print("\nglyph fit with the joining polyline masked out entirely (independent check):")
    C2, _ = glyph_pass(ink, fc, seeds, lines, mask_line=True)
    dd = np.hypot(*(C - C2).T)
    print(f"  agreement max {dd.max():.2f} px, rms {dd.std():.2f} px:  {np.round(dd, 2)}")

    py = blend_axis(P["left"][0], P["left"][1], P["right"][0], P["right"][1], Y_ORDER, Y0, Y1)
    px = blend_axis(P["bottom"][0], P["bottom"][1], P["top"][0], P["top"][1], X_ORDER, X0, X1)
    s, v = to_sv(C[:, 0], C[:, 1])
    X = (1 - v) * np.polyval(px[0], s) + v * np.polyval(px[1], s)
    Y = (1 - s) * np.polyval(py[0], v) + s * np.polyval(py[1], v)

    print("\nfree structural check -- the abscissae against the 0.025 lattice (never told):")
    lat = X0 + 0.025 * np.arange(N)
    for xo in (1, 2, 3):
        Xv, _ = axis_values(ink, SPEC, tickset, C, FRAME_ORDER, xo, Y_ORDER)
        e = Xv - lat
        print(
            f"  x order {xo}: max {np.abs(e).max():.5f} ({np.abs(e).max() * 3018:.2f} px)  "
            f"sd {e.std():.5f} ({e.std() * 3018:.2f} px)  mean {e.mean():+.5f}"
        )
    print(f"  chosen model: {np.round(X - lat, 5)}")

    print("\nmodel spread -- every defensible calibration, same glyph centres:")
    var = [axis_values(ink, SPEC, tickset, C, *m) for m in VARIANTS]
    VX = np.array([a for a, _ in var])
    VY = np.array([b for _, b in var])
    print(f"  all {len(var)}: x sd max {VX.std(0).max():.5f}   y sd max {VY.std(0).max():.6f}")
    ok = [i for i, m in enumerate(VARIANTS) if m[2] in (3, 4) and m[1] == 3]
    mx, my = VX[ok].std(0), VY[ok].std(0)
    print(f"  order-3/4 subset ({len(ok)}): x sd max {mx.max():.5f}   y sd max {my.max():.6f}")
    Pc = tickset(fc, to_sv, fine=False)
    pyc = blend_axis(Pc["left"][0], Pc["left"][1], Pc["right"][0], Pc["right"][1], Y_ORDER, Y0, Y1)
    pxc = blend_axis(Pc["bottom"][0], Pc["bottom"][1], Pc["top"][0], Pc["top"][1], X_ORDER, X0, X1)
    Xc = (1 - v) * np.polyval(pxc[0], s) + v * np.polyval(pxc[1], s)
    Yc = (1 - s) * np.polyval(pyc[0], v) + s * np.polyval(pyc[1], v)
    print(
        f"  long ticks vs all ticks: max dx {np.abs(X - Xc).max():.5f}, "
        f"max dy {np.abs(Y - Yc).max():.6f}"
    )
    mx = np.hypot(mx, np.abs(X - Xc))
    my = np.hypot(my, np.abs(Y - Yc))

    glyph_px = np.maximum(dd, 0.35)
    unc = point_sigma(P, fc, C, s, v, px, py, X_ORDER, Y_ORDER, glyph_px, mx, my)
    sx, sy, cal_x, cal_y = unc["sx"], unc["sy"], unc["cal_x"], unc["cal_y"]
    dXdpx, dYdpx = unc["dXdpx"], unc["dYdpx"]

    overlay(
        ink,
        to_px,
        C,
        ((px[0], [X0 + 0.05 * k for k in range(12)]), (py[0], [Y0 + 0.01 * k for k in range(15)])),
        [
            ((400, 260, 2200, 2710), (900, 1225), OUT / "a9_verify.png"),
            ((430, 560, 780, 700), (1400, 560), OUT / "a9_verify_flat.png"),
            ((2050, 2380, 2200, 2560), (750, 900), OUT / "a9_verify_m23.png"),
        ],
    )

    print(
        "\n   k       x          y       1s_x    1s_y   | cal_x  glyph_x model_x |"
        "  cal_y   glyph_y model_y"
    )
    for i in range(N):
        print(
            f"  {i + 1:2d} {X[i]:9.5f} {Y[i]:9.6f} {sx[i]:7.5f} {sy[i]:7.6f} | "
            f"{cal_x[i]:6.5f} {glyph_px[i] * dXdpx[i]:7.5f} {mx[i]:7.5f} | "
            f"{cal_y[i]:7.6f} {glyph_px[i] * dYdpx[i]:7.6f} {my[i]:7.6f}"
        )
    # The CSV is maintained by hand; these values are not. See `check_against_csv`.
    rows = check_against_csv(
        CSV,
        {
            "x": (X, "%.5f"),
            "y": (Y, "%.6f"),
            "sx": (sx, "%.5f"),
            "sy": (sy, "%.6f"),
        },
        label="f9 (Figure A9)",
    )
    # Nor is the calibration its header publishes, which `check_straightness` reads back.
    cal = check_calibration(CSV, calibration_lines(fc, px, py), label="f9 (Figure A9)")
    return max(rows, cal)


if __name__ == "__main__":
    raise SystemExit(main())
