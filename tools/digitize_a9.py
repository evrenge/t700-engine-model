"""Re-extract Figure A9 (`f9`, power turbine corrected mass flow) from TM-100991 pdf p.64.

Why this figure was redone.  `f9` is the highest-leverage map in the engine: a sensitivity
sweep of the trimmed model gives +2.30 % NG and +7.46 % SHP per +1 % on `f9` at hover, so
a 0.2 % calibration error is worth 1.5 % of shaft power.  The 2026-09-10 extraction used a
per-axis *linear* map on a perspective-rectified render, and A7 later showed that a linear
axis map cannot fit these scans.

Method, following `tools/digitize_a7.py`; see `tools/figure/` for the machinery.
What is specific to this page:

*   `x` prints a value only at the RIGHT frame (0.85).  The left frame is unlabelled, so
    it is a *prediction*: the cubic tick fits, told nothing about it, put it at 0.30076 and
    0.30068 (all ticks, bottom and top edges) and 0.30050 / 0.30028 (long ticks alone).
    That is the first independent measurement of the 0.30 the 2026-09-10 header inferred
    from a tick-count argument, and it lands inside two pixels of it.
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

from figure.axis import axis_values, blend_axis, held_out, point_sigma, tail_cv_report
from figure.frame import corner_arc, curve_component, edge_at, fit_edges, make_sv, segment_lines
from figure.glyph import CENTRE_NATIVE, TEMPLATE_NATIVE, fit_glyphs, hide_segments, windows
from figure.provenance import calibration_lines, check_against_csv, check_calibration
from figure.raster import native_bitmap, overlay
from figure.ticks import find_ticks, label_by_map, label_ticks

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
    ink = tickset.ink
    at = edge_at(fc, to_sv)
    fine = FINE if fine is None else fine
    TL = corner_arc(fc, "top", "left", (469, 311))
    TR = corner_arc(fc, "top", "right", (2130, 321))
    BL = corner_arc(fc, "bottom", "left", (468, 2645))
    BR = corner_arc(fc, "bottom", "right", (2125, 2662))
    ends = {
        "left": (BL[1], TL[1], Y0, Y1),
        "right": (BR[1], TR[1], Y0, Y1),
        "bottom": (BL[0], BR[0], X0, X1),
        "top": (TL[0], TR[0], X0, X1),
    }

    P = {"_xedges": (X0, X1), "_yedges": (Y0, Y1)}
    for nm, (along, sign, lo, hi, th_long, th_all, step_c, step_f) in EDGES.items():
        w0, w1, v0, v1 = ends[nm]
        coarse = find_ticks(ink, fc[nm], along, sign, lo, hi, th_long)
        wc, vc, keep_c = label_ticks(coarse, w0, w1, v0, v1, step_c, tol=0.15)
        tc = at(nm, wc)
        pc = np.polyfit(tc, vc, 3)
        if fine:
            allt = find_ticks(ink, fc[nm], along, sign, lo, hi, th_all)
            wf, vf, keep_f = label_by_map(at(nm, allt), pc, step_f)
            t, val = wf, vf
            n_raw, n_kept = allt.size, int(keep_f.sum())
        else:
            t, val = tc, vc
            n_raw, n_kept = coarse.size, int(keep_c.sum())
        step = step_f if fine else step_c
        interior = (val > v0 + 0.4 * step) & (val < v1 - 0.4 * step)
        P[nm] = (t[interior], val[interior])
        if verbose:
            print(
                f"  {nm:6s}: {coarse.size} long ticks -> {keep_c.sum()} labelled;  "
                f"{n_raw} ticks at {step} -> {n_kept} labelled, {interior.sum()} interior"
            )
    return P


def report_models(P):
    for axis, lo, hi, edges in (
        ("y = W45C", Y0, Y1, ("left", "right")),
        ("x = Ps9/P45", None, X1, ("bottom", "top")),
    ):
        for nm in edges:
            t, val = P[nm]
            print(f"\n  --- {axis}, {nm} edge ({t.size} interior ticks) ---")
            lab = f"printed {hi} at 1" + ("" if lo is None else f", printed {lo} at 0")
            print(
                f"   frame-edge held out ({lab}; not in the fit).  edge0 unprinted = a prediction"
            )
            held_out(t, val, lo, hi, unit=1e-5, uname="e-5")
            print("   tail CV (fit the middle 64 % of the lattice, predict the two ends):")
            tail_cv_report(t, val, unit=1e-5, uname="e-5")


def seeds_from_lattice(ink, to_sv, cur):
    Y = np.mgrid[0 : ink.shape[0], 0 : ink.shape[1]][0]
    cnt = cur.sum(0)
    cols = np.flatnonzero(cnt > 0)
    yc = (cur * Y).sum(0)[cols] / cnt[cols]
    s = to_sv(cols.astype(float), yc)[0]
    seeds = []
    for k in range(N):
        j = int(np.argmin(np.abs(s - k / (N - 1.0))))
        seeds.append((float(cols[j]), float(yc[j])))
    seeds[0] = (469.4, 612.6)  # on the left frame, so not in `cur`
    seeds[-1] = (2126.0, 2468.0)  # on the right frame
    return np.array(seeds)


def glyph_pass(ink, fc, seeds, lines, mask_line=False):
    def left_mask(X, Y):
        return np.abs(X - np.polyval(fc["left"], Y)) > 5.0

    def right_mask(X, Y):
        return np.abs(X - np.polyval(fc["right"], Y)) > 5.0

    W = []
    for i in range(N):
        # The masks mark what to KEEP; `windows` takes what to EXCLUDE.
        keep = left_mask if i == 0 else (right_mask if i == N - 1 else None)
        mf = (lambda X, Y, k=keep: ~k(X, Y)) if keep is not None else None
        w = windows(ink, seeds[i : i + 1], 20.0, mask_fn=mf)[0]
        # The half-segment towards the PREVIOUS marker runs backwards along segment i-1.
        w["seg"] = (
            None if (mask_line or i == 0) else -lines[i - 1, 2:4],
            None if (mask_line or i == N - 1) else lines[i, 2:4],
        )
        W.append(w)
    if mask_line:
        hide_segments(W, seeds, lines)
    return fit_glyphs(
        W,
        seeds,
        [np.radians(-45.0), 15.0, np.radians(45.0), 15.0, 3.6, 2.4],
        rounds=6,
        seglen=32.0,
        # Trap 8's guard, adopted 2026-09-16 after measuring it. A9's own free check --
        # its abscissae fall on a 0.025 lattice the tool is never told about -- comes out
        # consistently better guarded, on the worst point at every order: 2.18 -> 2.10 px,
        # 1.99 -> 1.91, 1.51 -> 1.44 at the chosen order 3, with sd unchanged. Small, and
        # in the right direction on the map SCOPE.md puts at +7.46 % SHP per +1 %.
        simplex=2.0,
        max_move=10.0,
        centre=CENTRE_NATIVE,
        template=TEMPLATE_NATIVE,
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
    seeds = seeds_from_lattice(ink, to_sv, cur)
    ys, xs = np.nonzero(cur)
    lines = segment_lines(np.stack([xs, ys], 1).astype(float), seeds, clear=20.0, band=4.0)
    print("\njoining segments:")
    for i, ln in enumerate(lines):
        print(
            f"  seg{i + 1:2d}: len {ln[4]:6.1f} px  n {int(ln[5]):5d}  rms {ln[6]:.3f} px  "
            f"dir {np.degrees(np.arctan2(ln[3], ln[2])):+8.3f} deg"
        )

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
