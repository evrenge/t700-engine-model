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
from PIL import Image, ImageDraw

from digitize_native import (
    axis_values,
    blend_axis,
    corner,
    curve_component,
    fit_glyphs,
    frame_curves,
    held_out,
    label_ticks,
    make_sv,
    make_window,
    native_bitmap,
    poly_se,
    segment_lines,
    tail_cv,
    ticks,
)

OUT = Path("validation/out/digitize/a2redo")
CSV = Path("data/maps/f2_compressor_temperature.csv")
PAGE = 57
N = 20

SPEC = [
    ("left", 0, (560, 625, 320, 2610)),
    ("right", 0, (2225, 2290, 320, 2610)),
    ("top", 1, (265, 325, 620, 2230)),
    ("bottom", 1, (2605, 2670, 620, 2230)),
]

FRAME_ORDER, X_ORDER, Y_ORDER = 2, 3, 3  # chosen by the tests main() prints

VARIANTS = [
    (fo, xo, yo, bl) for fo in (1, 2) for xo in (1, 3) for yo in (1, 3, 4) for bl in (True, False)
]


def tickset(fc, to_sv, verbose=False):
    """Every printed tick on the four edges, labelled, frame-edge ticks removed."""
    ink = tickset.ink
    TL = corner(fc, "top", "left", (589, 291))
    TR = corner(fc, "top", "right", (2255, 296))
    BL = corner(fc, "bottom", "left", (593, 2631))
    BR = corner(fc, "bottom", "right", (2257, 2641))

    def at(nm, cs):
        cs = np.asarray(cs, float)
        if nm in ("left", "right"):
            return to_sv(np.polyval(fc[nm], cs), cs)[1]
        return to_sv(cs, np.polyval(fc[nm], cs))[0]

    raw = {
        "left": ticks(ink, fc["left"], 0, +1, 300, 2635, 7),
        "right": ticks(ink, fc["right"], 0, -1, 300, 2645, 7),
        "bottom": ticks(ink, fc["bottom"], 1, -1, 590, 2260, 3),
        "top": ticks(ink, fc["top"], 1, +1, 590, 2260, 3),
    }
    ends = {
        "left": (BL[1], TL[1], 1.1, 3.0, 0.1),
        "right": (BR[1], TR[1], 1.1, 3.0, 0.1),
        "bottom": (BL[0], BR[0], 0.0, 20.0, 1.0),
        "top": (TL[0], TR[0], 0.0, 20.0, 1.0),
    }
    P = {"_xedges": (0.0, 20.0), "_yedges": (1.1, 3.0)}
    for nm, r in raw.items():
        w, val, keep = label_ticks(r, *ends[nm])
        lo, hi, step = ends[nm][2], ends[nm][3], ends[nm][4]
        interior = (val > lo + 0.4 * step) & (val < hi - 0.4 * step)
        P[nm] = (at(nm, w[interior]), val[interior])
        if verbose:
            print(
                f"  {nm:6s}: {r.size} raw, {keep.sum()} labelled, {interior.sum()} interior;"
                f" dropped {np.round(r[~keep], 1)}"
            )
    return P


def report_models(P):
    for axis, lo, hi, unit, uname, edges in (
        ("y = T3/T2", 1.1, 3.0, 1e-3, " mU", ("left", "right")),
        ("x = PS3/P2", 0.0, 20.0, 1.0, "", ("bottom", "top")),
    ):
        for nm in edges:
            t, val = P[nm]
            print(f"\n  --- {axis}, {nm} edge ({t.size} interior ticks) ---")
            print(f"   frame-edge held out (printed {lo} at 0 and {hi} at 1, neither in the fit):")
            held_out(t, val, lo, hi, unit=unit, uname=uname)
            print("   tail CV (fit the middle 64 % of the lattice, predict the two ends):")
            tail_cv(t, val, unit=unit, uname=uname)


def seeds_from_lattice(ink, to_sv, cur):
    Y = np.mgrid[0 : ink.shape[0], 0 : ink.shape[1]][0]
    cnt = cur.sum(0)
    cols = np.flatnonzero(cnt > 0)
    yc = (cur * Y).sum(0)[cols] / cnt[cols]
    s = to_sv(cols.astype(float), yc)[0]
    seeds = []
    for k in range(1, N + 1):
        j = int(np.argmin(np.abs(s - k / 20.0)))
        seeds.append((float(cols[j]), float(yc[j])))
    seeds[-1] = (2257.0, 533.0)  # marker 20 sits ON the right frame, so it is not in `cur`
    return np.array(seeds)


def glyph_pass(ink, fc, seeds, lines, mask_line=False):
    """Fit the 20 glyphs.  `mask_line` is the independent check: hide the joining polyline
    from the score entirely instead of modelling it, and see whether the centres move."""

    def frame_mask(X, Y):
        off_frame = np.abs(X - np.polyval(fc["right"], Y)) > 5.0
        off_tick = ~((Y >= 536) & (Y <= 552) & (X <= 2256))  # the 2.8 major under marker 20
        return off_frame & off_tick

    W = []
    for i in range(N):
        W.append(
            make_window(
                ink,
                seeds[i, 0],
                seeds[i, 1],
                rad=20.0,
                half=24,
                mask_fn=frame_mask if i == N - 1 else None,
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
    return fit_glyphs(W, seeds, [np.radians(-45.0), 15.0, np.radians(45.0), 15.0, 3.3, 3.3])


def overlay(ink, to_px, C, px, py):
    im = Image.fromarray((~ink * 255).astype(np.uint8)).convert("RGB")
    d = ImageDraw.Draw(im)
    tt = np.linspace(0.0, 1.0, 20001)
    for val in range(21):
        u = float(tt[np.argmin(np.abs(np.polyval(px[0], tt) - val))])
        a, b = to_px(u, 0.0), to_px(u, 1.0)
        d.line([(a[0][0], a[1][0]), (b[0][0], b[1][0])], fill=(210, 230, 255), width=1)
    for k in range(20):
        v = float(tt[np.argmin(np.abs(np.polyval(py[0], tt) - (1.1 + 0.1 * k)))])
        a, b = to_px(0.0, v), to_px(1.0, v)
        d.line([(a[0][0], a[1][0]), (b[0][0], b[1][0])], fill=(210, 230, 255), width=1)
    for x, y in C:
        d.ellipse([x - 26, y - 26, x + 26, y + 26], outline=(220, 0, 0), width=3)
        d.line([(x - 40, y), (x + 40, y)], fill=(0, 170, 0), width=1)
        d.line([(x, y - 40), (x, y + 40)], fill=(0, 170, 0), width=1)
    im.crop((520, 240, 2330, 2700)).resize((905, 1230), Image.LANCZOS).save(OUT / "a2_verify.png")
    im.crop((2120, 420, 2320, 700)).resize((800, 1120), Image.LANCZOS).save(
        OUT / "a2_verify_m20.png"
    )
    im.crop((600, 2280, 900, 2460)).resize((1200, 720), Image.LANCZOS).save(
        OUT / "a2_verify_m1.png"
    )


def main() -> int:
    ink = native_bitmap(PAGE, OUT)
    tickset.ink = ink
    print("frame arcs (printed straight; the sagitta is the page's bow):")
    fc, _ = frame_curves(ink, SPEC, order=FRAME_ORDER)
    to_sv, to_px = make_sv(fc)

    print("\ncalibration ticks:")
    P = tickset(fc, to_sv, verbose=True)
    print("\nCALIBRATION MODEL SELECTION")
    report_models(P)
    print(f"\n  chosen: frame arcs order {FRAME_ORDER}, x order {X_ORDER}, y order {Y_ORDER};")
    print("  each edge fitted alone with its two printed frame values, then blended.")

    cur = curve_component(ink, fc)
    seeds = seeds_from_lattice(ink, to_sv, cur)
    ys, xs = np.nonzero(cur)
    lines = segment_lines(np.stack([xs, ys], 1).astype(float), seeds)
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
    sig = {}
    for nm, poly, order, arg in (
        ("left", py[0], Y_ORDER, v),
        ("right", py[1], Y_ORDER, v),
        ("bottom", px[0], X_ORDER, s),
        ("top", px[1], X_ORDER, s),
    ):
        t, val = P[nm]
        lo, hi = (1.1, 3.0) if nm in ("left", "right") else (0.0, 20.0)
        tt = np.concatenate([t, [0.0, 1.0]])
        vv = np.concatenate([val, [lo, hi]])
        sig[nm] = poly_se(tt, (vv - np.polyval(poly, tt)).std(), order, arg)
    cal_y = np.hypot((1 - s) * sig["left"], s * sig["right"])
    cal_x = np.hypot((1 - v) * sig["bottom"], v * sig["top"])
    hpx = np.polyval(fc["right"], C[:, 1]) - np.polyval(fc["left"], C[:, 1])
    vpx = np.polyval(fc["bottom"], C[:, 0]) - np.polyval(fc["top"], C[:, 0])
    dXdpx = np.abs(np.polyval(np.polyder(px[0]), s)) / hpx
    dYdpx = np.abs(np.polyval(np.polyder(py[0]), v)) / vpx
    glyph_px = np.maximum(dd, 0.35)  # polyline-masked disagreement, floored at the fit noise
    sx = np.sqrt(cal_x**2 + (glyph_px * dXdpx) ** 2 + mx**2)
    sy = np.sqrt(cal_y**2 + (glyph_px * dYdpx) ** 2 + my**2)

    py1 = blend_axis(P["left"][0], P["left"][1], P["right"][0], P["right"][1], 1, 1.1, 3.0)
    Ylin = (1 - s) * np.polyval(py1[0], v) + s * np.polyval(py1[1], v)

    overlay(ink, to_px, C, px, py)

    print(
        "\n   k        x           y      1s_x    1s_y   | cal_x  glyph_x model_x |"
        "  cal_y  glyph_y model_y | y(linear)"
    )
    for i in range(N):
        print(
            f"  {i + 1:2d} {X[i]:10.5f} {Y[i]:9.5f} {sx[i]:7.5f} {sy[i]:7.5f} | "
            f"{cal_x[i]:6.4f} {glyph_px[i] * dXdpx[i]:7.4f} {mx[i]:7.4f} | "
            f"{cal_y[i]:6.5f} {glyph_px[i] * dYdpx[i]:7.5f} {my[i]:7.5f} | {Ylin[i]:8.5f}"
        )
    print("\ncsv rows:")
    for i in range(N):
        print(f"{X[i]:.5f},{Y[i]:.5f}")
    print(f"\nheader and CSV are maintained by hand in {CSV}; compare the rows above.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
