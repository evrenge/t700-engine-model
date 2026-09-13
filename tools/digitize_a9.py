"""Re-extract Figure A9 (`f9`, power turbine corrected mass flow) from TM-100991 pdf p.64.

Why this figure was redone.  `f9` is the highest-leverage map in the engine: a sensitivity
sweep of the trimmed model gives +2.30 % NG and +7.46 % SHP per +1 % on `f9` at hover, so
a 0.2 % calibration error is worth 1.5 % of shaft power.  The 2026-09-10 extraction used a
per-axis *linear* map on a perspective-rectified render, and A7 later showed that a linear
axis map cannot fit these scans.

Method, following `tools/digitize_a7.py`; see `tools/digitize_native.py` for the machinery.
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
from PIL import Image, ImageDraw

from digitize import check_against_csv
from digitize_native import (
    axis_values,
    blend_axis,
    corner,
    curve_component,
    fit_glyphs,
    frame_curves,
    held_out,
    label_by_map,
    label_ticks,
    make_sv,
    make_window,
    native_bitmap,
    poly_se,
    segment_lines,
    tail_cv,
    ticks,
)

OUT = Path("validation/out/digitize/a9redo")
CSV = Path("data/maps/f9_pt_mass_flow.csv")
PAGE = 64
N = 23

SPEC = [
    ("left", 0, (435, 505, 340, 2630)),
    ("right", 0, (2095, 2165, 340, 2630)),
    ("top", 1, (285, 350, 500, 2100)),
    ("bottom", 1, (2620, 2690, 500, 2100)),
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
    fine = FINE if fine is None else fine
    TL = corner(fc, "top", "left", (469, 311))
    TR = corner(fc, "top", "right", (2130, 321))
    BL = corner(fc, "bottom", "left", (468, 2645))
    BR = corner(fc, "bottom", "right", (2125, 2662))
    ends = {
        "left": (BL[1], TL[1], Y0, Y1),
        "right": (BR[1], TR[1], Y0, Y1),
        "bottom": (BL[0], BR[0], X0, X1),
        "top": (TL[0], TR[0], X0, X1),
    }

    def at(nm, cs):
        cs = np.asarray(cs, float)
        if nm in ("left", "right"):
            return to_sv(np.polyval(fc[nm], cs), cs)[1]
        return to_sv(cs, np.polyval(fc[nm], cs))[0]

    P = {"_xedges": (X0, X1), "_yedges": (Y0, Y1)}
    for nm, (along, sign, lo, hi, th_long, th_all, step_c, step_f) in EDGES.items():
        w0, w1, v0, v1 = ends[nm]
        coarse = ticks(ink, fc[nm], along, sign, lo, hi, th_long)
        wc, vc, keep_c = label_ticks(coarse, w0, w1, v0, v1, step_c, tol=0.15)
        tc = at(nm, wc)
        pc = np.polyfit(tc, vc, 3)
        if fine:
            allt = ticks(ink, fc[nm], along, sign, lo, hi, th_all)
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
            tail_cv(t, val, unit=1e-5, uname="e-5")


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
        mf = left_mask if i == 0 else (right_mask if i == N - 1 else None)
        W.append(
            make_window(
                ink,
                seeds[i, 0],
                seeds[i, 1],
                rad=20.0,
                half=24,
                mask_fn=mf,
                seg_in=(None if (mask_line or i == 0) else lines[i - 1, 2:4]),
                seg_out=(None if (mask_line or i == N - 1) else lines[i, 2:4]),
            )
        )
    if mask_line:
        for i, w in enumerate(W):
            for j, sign in ((i - 1, -1.0), (i, +1.0)):
                if 0 <= j < N - 1:
                    d = lines[j, 2:4] * sign
                    dx = w["X"] - seeds[i, 0]
                    dy = w["Y"] - seeds[i, 1]
                    t = dx * d[0] + dy * d[1]
                    p = -dx * d[1] + dy * d[0]
                    w["use"] &= ~((t > 6) & (np.abs(p) < 6))
    return fit_glyphs(
        W, seeds, [np.radians(-45.0), 15.0, np.radians(45.0), 15.0, 3.6, 2.4], seglen=32.0
    )


def overlay(ink, to_px, C, px, py):
    im = Image.fromarray((~ink * 255).astype(np.uint8)).convert("RGB")
    d = ImageDraw.Draw(im)
    tt = np.linspace(0.0, 1.0, 20001)
    for k in range(12):
        u = float(tt[np.argmin(np.abs(np.polyval(px[0], tt) - (X0 + 0.05 * k)))])
        a, b = to_px(u, 0.0), to_px(u, 1.0)
        d.line([(a[0][0], a[1][0]), (b[0][0], b[1][0])], fill=(210, 230, 255), width=1)
    for k in range(15):
        v = float(tt[np.argmin(np.abs(np.polyval(py[0], tt) - (Y0 + 0.01 * k)))])
        a, b = to_px(0.0, v), to_px(1.0, v)
        d.line([(a[0][0], a[1][0]), (b[0][0], b[1][0])], fill=(210, 230, 255), width=1)
    for x, y in C:
        d.ellipse([x - 26, y - 26, x + 26, y + 26], outline=(220, 0, 0), width=3)
        d.line([(x - 40, y), (x + 40, y)], fill=(0, 170, 0), width=1)
        d.line([(x, y - 40), (x, y + 40)], fill=(0, 170, 0), width=1)
    im.crop((400, 260, 2200, 2710)).resize((900, 1225), Image.LANCZOS).save(OUT / "a9_verify.png")
    im.crop((430, 560, 780, 700)).resize((1400, 560), Image.LANCZOS).save(
        OUT / "a9_verify_flat.png"
    )
    im.crop((2050, 2380, 2200, 2560)).resize((750, 900), Image.LANCZOS).save(
        OUT / "a9_verify_m23.png"
    )


def main() -> int:
    ink = native_bitmap(PAGE, OUT)
    tickset.ink = ink
    print("frame arcs (printed straight; the sagitta is the page's bow):")
    fc, _ = frame_curves(ink, SPEC, order=FRAME_ORDER)
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

    sig = {}
    for nm, poly, order, arg in (
        ("left", py[0], Y_ORDER, v),
        ("right", py[1], Y_ORDER, v),
        ("bottom", px[0], X_ORDER, s),
        ("top", px[1], X_ORDER, s),
    ):
        t, val = P[nm]
        lo, hi = (Y0, Y1) if nm in ("left", "right") else (X0, X1)
        tt = np.concatenate([t, [0.0, 1.0]])
        vv = np.concatenate([val, [lo, hi]])
        sig[nm] = poly_se(tt, (vv - np.polyval(poly, tt)).std(), order, arg)
    cal_y = np.hypot((1 - s) * sig["left"], s * sig["right"])
    cal_x = np.hypot((1 - v) * sig["bottom"], v * sig["top"])
    hpx = np.polyval(fc["right"], C[:, 1]) - np.polyval(fc["left"], C[:, 1])
    vpx = np.polyval(fc["bottom"], C[:, 0]) - np.polyval(fc["top"], C[:, 0])
    dXdpx = np.abs(np.polyval(np.polyder(px[0]), s)) / hpx
    dYdpx = np.abs(np.polyval(np.polyder(py[0]), v)) / vpx
    glyph_px = np.maximum(dd, 0.35)
    sx = np.sqrt(cal_x**2 + (glyph_px * dXdpx) ** 2 + mx**2)
    sy = np.sqrt(cal_y**2 + (glyph_px * dYdpx) ** 2 + my**2)

    py1 = blend_axis(P["left"][0], P["left"][1], P["right"][0], P["right"][1], 1, Y0, Y1)
    Ylin = (1 - s) * np.polyval(py1[0], v) + s * np.polyval(py1[1], v)

    overlay(ink, to_px, C, px, py)

    print(
        "\n   k       x          y       1s_x    1s_y   | cal_x  glyph_x model_x |"
        "  cal_y   glyph_y model_y | y(linear)"
    )
    for i in range(N):
        print(
            f"  {i + 1:2d} {X[i]:9.5f} {Y[i]:9.6f} {sx[i]:7.5f} {sy[i]:7.6f} | "
            f"{cal_x[i]:6.5f} {glyph_px[i] * dXdpx[i]:7.5f} {mx[i]:7.5f} | "
            f"{cal_y[i]:7.6f} {glyph_px[i] * dYdpx[i]:7.6f} {my[i]:7.6f} | {Ylin[i]:9.6f}"
        )
    # The CSV is maintained by hand; these values are not. See `check_against_csv`.
    return check_against_csv(
        CSV,
        {
            "x": (X, "%.5f"),
            "y": (Y, "%.6f"),
            "sx": (sx, "%.5f"),
            "sy": (sy, "%.6f"),
        },
        label="f9 (Figure A9)",
    )


if __name__ == "__main__":
    raise SystemExit(main())
