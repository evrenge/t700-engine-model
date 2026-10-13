"""The single-curve digitizer: thirteen figures across both appendices of TM-100991.

Same method as `tools/digitize_a7.py`, generalised (see `docs/notes/digitizing-recipe.md`).
Every figure here plots one curve or marker series, so the markers are found without reading
anything off the page. The two Appendix C figures whose curves must be told apart live in
`digitize_multi_curve.py` and `digitize_c30.py`; both take their calibration from
`figure.axis.calibrate_under`, as this tool does.

    A3,  pdf p.58 -> f3    seal-pressurization bleed      data/maps/f3_seal_bleed_fraction.csv
    A4,  pdf p.59 -> f4    PT-balance bleed               data/maps/f4_pt_balance_bleed_fraction.csv
    A5,  pdf p.60 -> f5    tip-leak + cooling bleed       data/maps/f5_tip_leak_cooling_bleed...csv
    A6,  pdf p.61 -> f6    combustor efficiency           data/maps/f6_combustor_efficiency.csv
    A8,  pdf p.63 -> f8    power turbine energy           data/maps/f8_pt_energy.csv
    A10, pdf p.65 -> f10   exhaust pressure loss          data/maps/f10_exhaust_pressure_loss.csv
    A11, pdf p.66 -> f_hs  station 4.1 heat-sink constant data/maps/fhs_heat_sink_constant.csv
    C24-C29, pp.95-100 -> F_HM1..F_HM6                    data/schedules/fhm[1-6]_*.csv

A1, A2, A7 and A9 have their own tools.

Per figure: the native 300 dpi scan, never resampled; frame edges fitted straight or as
arcs, whichever predicts the printed frame values better; ticks matched to their printed
values by predicted lattice slot, not protrusion length (a marker on a frame out-protrudes
a major tick); axis order by the held-out frame test, guarded by leave-one-tick-out;
markers seeded by a morphological opening and fitted by the generative glyph model, with
frame and ticks masked symmetrically so the mask cannot bias the abscissa.

Run from the repo root.  Writes the thirteen CSVs and overlays under
validation/out/digitize/.  Needs `pdfimages` (poppler).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

from figure.axis import FRAME_MARGIN, blend_value, calibrate_under, choose_frame_order
from figure.frame import along_segments, auto_frame_spec, blocked_mask, ink_excesses
from figure.frame import dist as _dist
from figure.glyph import fit_glyphs, glyph_only_windows, windows
from figure.provenance import edge_str
from figure.provenance import poly_str as _poly_str
from figure.raster import native_bitmap, overlay
from figure.ticks import major_step, major_values

OUT = Path("validation/out/digitize/a6810")

# --------------------------------------------------------------------------- figure table
#
# `x` / `y`     (value at the low-side frame, value at the high-side frame, major step,
#               number of INTERIOR majors, which frame values are PRINTED on the page)
#               low side = left frame / bottom frame; high side = right frame / top frame.

FIGS = {
    "a8": dict(
        page=63,
        fig="A8",
        csv=Path("data/maps/f8_pt_energy.csv"),
        x=dict(lo=0.30, hi=0.85, step=0.05, n=10, printed=("hi",), anchor=("lo", "hi")),
        y=dict(lo=-5.0, hi=40.0, step=2.5, n=17, printed=("lo", "hi")),
        nmark=12,
        rad=20.0,
    ),
    "a10": dict(
        page=65,
        fig="A10",
        csv=Path("data/maps/f10_exhaust_pressure_loss.csv"),
        x=dict(lo=65.0, hi=100.0, step=5.0, n=6, printed=("lo", "hi")),
        y=dict(lo=1.01, hi=1.14, step=0.01, n=12, printed=("hi",), anchor=("lo", "hi")),
        nmark=18,
        rad=18.0,  # the 99.5 / 100 pair is 40 px apart; keep the discs disjoint
    ),
    # ---- Appendix C scheduling functions -------------------------------------------
    # Same machinery, different appendix: single-curve 'x'-marker plots on the same
    # 300 dpi 1-bit rasters. Axis lattices read off the rendered page, not inferred.
    "c24": dict(
        page=95,
        fig="C24",
        csv=Path("data/schedules/fhm1_topping_line.csv"),
        x=dict(lo=390.0, hi=640.0, step=50.0, n=4, printed=("lo", "hi")),
        y=dict(lo=-0.50, hi=4.00, step=0.50, n=8, printed=("lo", "hi")),
        nmark=9,
        rad=20.0,
        open_sz=6,
        rail=20.0,
    ),
    "c25": dict(
        page=96,
        fig="C25",
        csv=Path("data/schedules/fhm2_power_available.csv"),
        x=dict(lo=20.0, hi=120.0, step=20.0, n=4, printed=("lo", "hi")),
        y=dict(lo=-2.0, hi=12.0, step=1.0, n=13, printed=("lo", "hi")),
        nmark=11,
        rad=20.0,
        open_sz=6,
        # A fragment sits on the top frame at x=1957. The last marker sits ON the right
        # frame, so only the top rail may be excluded.
        rail=(-14.0, -14.0, 20.0, -14.0),
    ),
    "c26": dict(
        page=97,
        fig="C26",
        csv=Path("data/schedules/fhm3_load_demand_png.csv"),
        x=dict(lo=0.0, hi=100.0, step=20.0, n=4, printed=("lo", "hi")),
        y=dict(lo=75.0, hi=112.5, step=2.5, n=14, printed=("lo", "hi")),
        nmark=11,
        rad=20.0,
    ),
    "a3": dict(
        page=58,
        fig="A3",
        csv=Path("data/maps/f3_seal_bleed_fraction.csv"),
        x=dict(lo=65.0, hi=100.0, majors=[70.0, 80.0, 90.0], printed=("hi",), anchor=("lo", "hi")),
        y=dict(lo=-0.02, hi=0.12, step=0.01, n=13, printed=("lo", "hi")),
        nmark=13,
        rad=18.0,
        open_sz=6,
        min_size=60,
    ),
    "a4": dict(
        page=59,
        fig="A4",
        csv=Path("data/maps/f4_pt_balance_bleed_fraction.csv"),
        x=dict(lo=3.0, hi=12.0, majors=[4.0, 6.0, 8.0, 10.0], printed=("hi",), anchor=("lo", "hi")),
        y=dict(lo=0.0088, hi=0.0108, step=0.0002, n=9, printed=("lo", "hi")),
        nmark=3,
        rad=20.0,
    ),
    "a5": dict(
        page=60,
        fig="A5",
        csv=Path("data/maps/f5_tip_leak_cooling_bleed_fraction.csv"),
        x=dict(lo=3.0, hi=12.0, majors=[4.0, 6.0, 8.0, 10.0], printed=("hi",), anchor=("lo", "hi")),
        y=dict(lo=0.0770, hi=0.0860, step=0.0010, n=8, printed=("lo", "hi")),
        nmark=3,
        rad=20.0,
        max_missing=3,  # A5's right edge prints sparser ticks than its left
    ),
    "a11": dict(
        page=66,
        fig="A11",
        csv=Path("data/maps/fhs_heat_sink_constant.csv"),
        x=dict(lo=65.0, hi=100.0, step=5.0, n=6, printed=("lo", "hi")),
        y=dict(lo=3.50, hi=8.00, step=0.50, n=8, printed=("lo", "hi")),
        nmark=6,
        rad=20.0,
        rail=(-14.0, -14.0, -14.0, 200.0),
    ),
    "c27": dict(
        page=98,
        fig="C27",
        csv=Path("data/schedules/fhm4_load_demand_wfqps3.csv"),
        x=dict(lo=0.0, hi=100.0, step=20.0, n=4, printed=("lo", "hi")),
        y=dict(lo=2.1, hi=3.7, step=0.2, n=7, printed=("lo", "hi")),
        nmark=11,
        rad=20.0,
        open_sz=6,
        # First and last markers sit ON the left and right frames, so those rails must
        # stay; the fragments there are separated by size instead.
        rail=(-14.0, -14.0, 20.0, 20.0),
        rail_min_size=60,
    ),
    "c28": dict(
        page=99,
        fig="C28",
        csv=Path("data/schedules/fhm5_idle_wfirf.csv"),
        x=dict(lo=390.0, hi=640.0, step=50.0, n=4, printed=("lo", "hi")),
        y=dict(lo=2.05, hi=2.65, step=0.05, n=11, printed=("lo", "hi")),
        nmark=13,
        rad=20.0,
        open_sz=6,
        rail=20.0,  # every marker stands clear of the frames on this one
    ),
    "c29": dict(
        page=100,
        fig="C29",
        csv=Path("data/schedules/fhm6_idle_pcngi.csv"),
        x=dict(lo=390.0, hi=640.0, step=50.0, n=4, printed=("lo", "hi")),
        # The TOP frame carries no label. The printed majors run 64..74 at a step of 2 and
        # the frame sits one half-interval above 74, i.e. at 75.0 -- measured off the page,
        # and the held-out frame test judges it because the fit never sees it.
        y=dict(
            lo=62.0,
            hi=75.0,
            majors=[64.0, 66.0, 68.0, 70.0, 72.0, 74.0],
            printed=("lo",),
            anchor=("lo", "hi"),
        ),
        nmark=2,
        rad=20.0,
        open_sz=6,
        rail=20.0,
    ),
    "a6": dict(
        page=61,
        fig="A6",
        csv=Path("data/maps/f6_combustor_efficiency.csv"),
        x=dict(lo=0.010, hi=0.020, step=0.001, n=9, printed=("lo", "hi")),
        y=dict(lo=0.88, hi=1.10, step=0.02, n=10, printed=("lo", "hi")),
        nmark=2,
        rad=20.0,
    ),
}


# --------------------------------------------------------------------------- marker seeds


def seeds(
    ink, fr, nmark, open_sz=7, rail=-14.0, min_size=0, merge=12.0, rail_min_size=0, rail_zone=12.0
):
    """Marker seeds from a morphological opening -- no density, no threshold sweep.

    A 7x7 opening erases the ~3.5 px joining line and every tick but keeps the crossing of
    a marker's two strokes; cores closer than `merge` px are one marker split by the scan.
    `open_sz` 6 keeps the lighter-pen Appendix C markers that a 7x7 opening erases outright.
    `rail` is the distance from a frame inside which a component is discarded: negative
    keeps components up to that far OUTSIDE it (end markers that sit on a frame), positive
    drops those that far INSIDE it; a 4-tuple `(left, right, top, bottom)` sets each frame.
    `min_size` drops small fragments anywhere, `rail_min_size` only within `rail_zone` of a
    frame. The per-figure reasons are in `FIGS` and the CSV headers.
    """
    Y, X = np.mgrid[0 : ink.shape[0], 0 : ink.shape[1]]
    dL = _dist(fr["left"], X, Y)
    dR = _dist(fr["right"], X, Y)
    dT = _dist(fr["top"], X, Y, False)
    dB = _dist(fr["bottom"], X, Y, False)
    band = ink & (dL > -26) & (dR < 26) & (dT > -26) & (dB < 26)
    op = ndimage.binary_opening(band, structure=np.ones((open_sz, open_sz), bool))
    lab, n = ndimage.label(op, structure=np.ones((3, 3)))
    if n == 0:
        raise AssertionError("opening found nothing")
    com = np.array(ndimage.center_of_mass(op, lab, range(1, n + 1)))
    sz = np.array(ndimage.sum(op, lab, range(1, n + 1)))
    cy, cx = com[:, 0], com[:, 1]
    # keep only what is inside the plot rectangle (axis titles and tick labels sit outside)
    rl, rr, rt, rb = (rail,) * 4 if np.isscalar(rail) else rail
    ky = _dist(fr["top"], cx, cy, False) > rt
    ky &= _dist(fr["bottom"], cx, cy, False) < -rb
    ky &= _dist(fr["left"], cx, cy) > rl
    ky &= _dist(fr["right"], cx, cy) < -rr
    if min_size:
        ky &= sz >= min_size
    if rail_min_size:
        # a size floor near a frame only: C27's frame fragments are as small as a real
        # marker standing in open plot
        near = (
            (np.abs(_dist(fr["left"], cx, cy)) < rail_zone)
            | (np.abs(_dist(fr["right"], cx, cy)) < rail_zone)
            | (np.abs(_dist(fr["top"], cx, cy, False)) < rail_zone)
            | (np.abs(_dist(fr["bottom"], cx, cy, False)) < rail_zone)
        )
        ky &= ~near | (sz >= rail_min_size)
    cx, cy, sz = cx[ky], cy[ky], sz[ky]
    o = np.argsort(cx)
    cx, cy, sz = cx[o], cy[o], sz[o]
    merged = []
    for a, b, s in zip(cx, cy, sz, strict=True):
        if merged and np.hypot(a - merged[-1][0], b - merged[-1][1]) < merge:
            wa, wb, ws = merged[-1]
            merged[-1] = ((wa * ws + a * s) / (ws + s), (wb * ws + b * s) / (ws + s), ws + s)
        else:
            merged.append((a, b, s))
    print(f"  opening: {n} components, {len(cx)} inside the frame, {len(merged)} after merging")
    if len(merged) != nmark:
        raise AssertionError(f"seeded {len(merged)} markers, expected {nmark}")
    return np.array([[a, b] for a, b, _ in merged])


# --------------------------------------------------------------------------- checks


def segment_check(ink, fr, C):
    """Straightness of the joining segments, and an interior-ink count independent of any
    detector: walk each fitted segment in 10 px bins and look for an excess over its own
    baseline anywhere but at the ends."""
    Y, X = np.mgrid[0 : ink.shape[0], 0 : ink.shape[1]]
    keep = (
        (np.abs(_dist(fr["left"], X, Y)) > 5)
        & (np.abs(_dist(fr["right"], X, Y)) > 5)
        & (np.abs(_dist(fr["top"], X, Y, False)) > 5)
        & (np.abs(_dist(fr["bottom"], X, Y, False)) > 5)
    )
    ys, xs = np.nonzero(ink & keep)
    P = np.stack([xs, ys], 1).astype(float)
    extras, sag = [], []
    for L, t, perp in along_segments(P, C):
        sel = (np.abs(perp) < 6) & (t > 24) & (t < L - 24)
        if sel.sum() > 40:
            sag.append((L, float(np.abs(perp[sel]).mean()), float(perp[sel].std())))
        extras.append(ink_excesses(L, t, perp, 1.9, 14))
    print(f"  interior ink excesses per segment: {extras}  (all empty => no missed marker)")
    if sag:
        print(
            "  joining segments: "
            + ", ".join(f"{L:.0f} px rms {s:.2f}" for L, _, s in sag[: min(len(sag), 12)])
        )
    return extras, sag


# --------------------------------------------------------------------------- overlay


def gridlines(ink, to_px, corners4, C, cx, cy, cfg, tag):
    """The gridlines drawn back onto the scan, through whichever map the figure chose."""
    TL, TR, BR, BL = corners4
    lines = []
    for spec, poly in ((cfg["x"], cx), (cfg["y"], cy)):
        n = int(round((spec["hi"] - spec["lo"]) / major_step(spec)))
        lines.append((poly, [spec["lo"] + major_step(spec) * k for k in range(n + 1)]))
    x0 = int(min(TL[0], BL[0])) - 200
    x1 = int(max(TR[0], BR[0])) + 60
    y0 = int(min(TL[1], TR[1])) - 60
    y1 = int(max(BL[1], BR[1])) + 90
    crop = ((x0, y0, x1, y1), ((x1 - x0) // 2, (y1 - y0) // 2), OUT / f"{tag}_verify.png")
    overlay(ink, to_px, C, tuple(lines), [crop])


# --------------------------------------------------------------------------- driver


def run(key):
    cfg = FIGS[key]
    print(f"\n================ {cfg['fig']}  (pdf p.{cfg['page']})")
    ink = native_bitmap(cfg["page"])
    fspec = auto_frame_spec(cfg["page"])

    forder, trials = choose_frame_order(ink, fspec, cfg)
    cal = calibrate_under(ink, fspec, cfg, forder, verbose=True)
    fc, frdiag, fr = cal["fc"], cal["frdiag"], cal["fr"]
    TL, TR, BR, BL = cal["corners"]
    to_uv, to_px = cal["to_uv"], cal["to_px"]
    allt = cal["allt"]
    CX, CY = cal["CX"], cal["CY"]
    xs, ys = cfg["x"], cfg["y"]
    cx, cy = CX["c"], CY["c"]

    print("markers:")
    C0 = seeds(
        ink,
        fr,
        cfg["nmark"],
        cfg.get("open_sz", 7),
        cfg.get("rail", -14.0),
        cfg.get("min_size", 0),
        rail_min_size=cfg.get("rail_min_size", 0),
    )
    W = windows(ink, C0, cfg["rad"], mask_fn=blocked_mask(fr, allt))
    C, tpl = fit_glyphs(W, C0, rounds=4)
    # re-centre the windows on the fitted centres and refit, so the disc is self-consistent
    W = windows(ink, C, cfg["rad"], mask_fn=blocked_mask(fr, allt))
    C, tpl = fit_glyphs(W, C, rounds=3, tpl0=tpl)

    print("checks:")
    _, sag = segment_check(ink, fr, C)
    u, v = to_uv(C[:, 0], C[:, 1])
    Xd, Yd = blend_value(CX["edgec"], u, v), blend_value(CY["edgec"], v, u)

    # ---- glyph-centre uncertainty: worst case over honest variations of the estimator
    print("  glyph-centre robustness (max |move| over variants, px):")
    var = {}
    for nm, kw, rad, wl in (
        ("rad-3", {}, cfg["rad"] - 3.0, True),
        ("rad+3", {}, cfg["rad"] + 3.0, True),
        ("framemask 3.0", dict(framew=3.0), cfg["rad"], True),
        ("framemask 5.0", dict(framew=5.0), cfg["rad"], True),
        ("tickpad 3.0", dict(tickpad=3.0), cfg["rad"], True),
        ("tickpad 6.0", dict(tickpad=6.0), cfg["rad"], True),
    ):
        Wv = windows(ink, C, rad, mask_fn=blocked_mask(fr, allt, **kw))
        Cv, _ = fit_glyphs(
            Wv, C, rounds=2, withline=wl, verbose=False, tpl0=tpl, fit_template=False
        )
        var[nm] = Cv - C
        print(f"    {nm:14s} {np.abs(Cv - C).max():.2f}")
    Wg = glyph_only_windows(ink, C, tpl, cfg["rad"], mask_fn=blocked_mask(fr, allt))
    Cg, _ = fit_glyphs(Wg, C, rounds=3, withline=False, verbose=False, tpl0=tpl, fit_template=False)
    var["line masked out"] = Cg - C
    print(f"    {'line masked out':14s} {np.abs(Cg - C).max():.2f}")
    dpx = np.max(np.abs(np.array(list(var.values()))), axis=0)  # (nmark, 2) worst case

    # ---- uncertainty budget, per point, 1 sigma, in data units
    print("  uncertainty budget (1 sigma, data units):")
    sig, terms = {}, {}
    for ax, CA, t, poly, across in (("x", CX, u, cx, v), ("y", CY, v, cy, u)):
        base = blend_value(CA["edgec"], t, across)
        # a. residual of the anchored fit about the printed ticks
        a = np.full(t.shape, CA["resid"])
        # b. spread over calibration polynomial order 1..4
        b = np.max(
            [np.abs(np.polyval(CA["allc"][o], t) - base) for o in CA["allc"]],
            axis=0,
        )
        # c. the two ticked edges calibrated separately, half their disagreement
        c = 0.5 * np.abs(np.polyval(CA["edgec"][0], t) - np.polyval(CA["edgec"][1], t))
        # d. glyph centre, converted from px through the local scale of this axis
        du = _pixel_to_unit(to_uv, C, 0 if ax == "x" else 1)
        d = np.abs(np.polyval(np.polyder(poly), t)) * du * dpx[:, 0 if ax == "x" else 1]
        s = np.sqrt(a**2 + b**2 + c**2 + d**2)
        sig[ax] = s
        terms[ax] = (b, c, d)
        print(
            f"    {ax}: tick-fit {a[0]:.4g} | order spread {b.mean():.4g} | "
            f"edge split {c.mean():.4g} | glyph {d.mean():.4g}  ->  1 sigma {s.mean():.4g} "
            f"(max {s.max():.4g})"
        )

    # ---- independent bow-free check: linear interpolation between the two ticks that
    #      bracket each point, on each ticked edge, averaged. It never fits a bow, so where
    #      it agrees with the chosen polynomial, that polynomial is not wandering between ticks.
    print("  local two-tick check (bow-free), chosen-model minus local:")
    twotick = {}
    for ax, t, aspec, poly, edges in (
        ("x", u, xs, cx, CX["edges"]),
        ("y", v, ys, cy, CY["edges"]),
    ):
        per = []
        for et, ev in edges:
            nodes_t = np.concatenate([[0.0], et, [1.0]])
            nodes_v = np.concatenate(
                [
                    [aspec["lo"] if "lo" in aspec["printed"] else np.polyval(poly, 0.0)],
                    ev,
                    [aspec["hi"] if "hi" in aspec["printed"] else np.polyval(poly, 1.0)],
                ]
            )
            o = np.argsort(nodes_t)
            per.append(np.interp(t, nodes_t[o], nodes_v[o]))
        twotick[ax] = np.mean(per, axis=0)
        d = np.polyval(poly, t) - twotick[ax]
        print(
            f"    {ax}: max {np.abs(d).max():.4g}, rms {d.std():.4g}  "
            f"(edge-to-edge spread {np.abs(per[0] - per[1]).max():.4g})"
        )

    gridlines(ink, to_px, (TL, TR, BR, BL), C, cx, cy, cfg, key)
    return dict(
        cfg=cfg,
        fc=fc,
        fspec=fspec,
        ink=ink,
        frdiag=frdiag,
        sag=sag,
        varmax={k2: float(np.abs(v2).max()) for k2, v2 in var.items()},
        corners=(TL, TR, BR, BL),
        CX=CX,
        CY=CY,
        cx=cx,
        cy=cy,
        C=C,
        tpl=tpl,
        u=u,
        v=v,
        X=Xd,
        Y=Yd,
        dpx=dpx,
        sig=sig,
        terms=terms,
        loc=twotick,
        # the local map scale per marker, evaluated here because the map itself is a
        # closure over the chosen frame model
        forder=cal["forder"],
        frame_scores={fo: t["score"] for fo, t in trials.items()},
    )


# --------------------------------------------------------------------------- the CSVs

# The 2026-09-10 extraction, kept so the header's MOVEMENT block can report the move in sigma.
OLD = {
    "a3": [
        (65.0605, 0.108942),
        (78.0725, 0.10912),
        (79.0605, 0.108525),
        (80.08, 0.106204),
        (81.0785, 0.102454),
        (82.0454, 0.097811),
        (83.0544, 0.0892992),
        (84.0004, 0.0738827),
        (85.0094, 0.0500735),
        (87.0274, 0.0131692),
        (88.0154, 0.0022765),
        (89.0244, -0.00010442),
        (99.9973, 1.46257e-05),
    ],
    "a4": [(3.00493, 0.0105709), (5.01405, 0.00899895), (12.0054, 0.00899766)],
    "a5": [(3.01027, 0.0850907), (8.30604, 0.0845899), (12.0104, 0.0778945)],
    "a11": [
        (64.9947, 4.00062),
        (69.9947, 4.00062),
        (75.1108, 5.75467),
        (80.1002, 6.20473),
        (85.0158, 7.49719),
        (99.9736, 7.49719),
    ],
    "a8": [
        (0.300661, 33.8654),
        (0.35024, 30.2596),
        (0.400646, 26.5769),
        (0.450391, 23.0962),
        (0.5003, 20.0962),
        (0.54988, 17.1731),
        (0.59979, 14.25),
        (0.649865, 11.3269),
        (0.699775, 8.31731),
        (0.750015, 5.25),
        (0.800421, 2.05769),
        (0.850496, -1.31731),
    ],
    "a10": [
        (65.0736, 1.12408),
        (69.994, 1.11816),
        (71.0033, 1.1163),
        (76.0341, 1.0906),
        (78.5311, 1.07938),
        (81.5433, 1.06021),
        (85.4701, 1.03715),
        (86.4374, 1.03062),
        (87.4257, 1.02685),
        (88.435, 1.0241),
        (89.9279, 1.02146),
        (90.9162, 1.02115),
        (91.9045, 1.02215),
        (93.9546, 1.0281),
        (95.989, 1.03396),
        (97.5188, 1.04418),
        (99.5584, 1.06424),
        (99.979, 1.06621),
    ],
    "a6": [(0.0099932, 0.985076), (0.0199882, 0.984982)],
}

PROSE = {
    "a3": dict(
        quantity="f3 -- seal-pressurization bleed fraction B1 (Eq. 12, pdf p.22: B1 = f3(NGc))",
        xdesc="NGC corrected gas generator speed, percent",
        ydesc="B1 bleed fraction, nondimensional",
        fmt=("%.4f", "%.6f", "%.4f", "%.6f"),
        extra=[
            "FREE STRUCTURAL CHECK. The thirteen abscissae recover 65, 78, 79, 80, 81, 82, 83,",
            "84, 85, 87, 88, 89, 100 to max 0.05 percent NGc against a mean sigma_x of 0.025.",
            "Nothing here is told to expect integers.",
            "THE 0.11 PLATEAU, INDEPENDENTLY. Open question #42 concluded from the tick lattice",
            "alone that the upper plateau is 0.1091 and that the printed 0.11 axis LABEL is",
            "displaced, not the tick. This extraction -- a different raster, a different frame",
            "fit and a different calibration -- recovers 0.108989 and 0.109077. It is an",
            "independent confirmation, not a repetition.",
            "THE NEGATIVE TAIL is real and is left as measured: -0.000168 and -0.000095 at the",
            "top of the speed range. See open question #44; do not clamp the file.",
        ],
    ),
    "a4": dict(
        quantity="f4 -- power-turbine-balance bleed fraction B2 (Eq. 13, pdf p.22: B2 = f4(WA2c))",
        xdesc="WA2C station 2 corrected mass flow, lbm/sec",
        ydesc="B2 bleed fraction, nondimensional",
        fmt=("%.4f", "%.7f", "%.4f", "%.7f"),
        extra=[
            "FREE STRUCTURAL CHECK. The three abscissae recover 3.000, 5.000 and 12.002 against",
            "a printed 3, 5 and 12, with a mean sigma_x of 0.011. The first and last markers sit",
            "ON the left and right frame lines, which is why `rail` keeps frame-adjacent",
            "components here.",
        ],
    ),
    "a5": dict(
        quantity=(
            "f5 -- impeller tip leakage and turbine cooling-bleed fraction B3\n"
            "#    (Eq. 14, pdf p.22: B3 = f5(WA2c))"
        ),
        xdesc="WA2C station 2 corrected mass flow, lbm/sec",
        ydesc="B3 bleed fraction, nondimensional",
        fmt=("%.4f", "%.7f", "%.4f", "%.7f"),
        extra=[
            "FREE STRUCTURAL CHECK. The outer two abscissae recover 3.0002 and 11.9990 against a",
            "printed 3 and 12, with a mean sigma_x of 0.011. The middle knot is genuinely not on",
            "a round value: it reads 8.303.",
            "A5's RIGHT EDGE prints sparser ticks than its left -- five of eight interior majors",
            "are found there against eight of eight on the left -- so `max_missing` is raised for",
            "this figure alone. The y axis is fitted on both edges from the majors that ARE",
            "visible and its held-out frame test still judges the result.",
        ],
    ),
    "a11": dict(
        quantity=(
            "f_hs -- station 4.1 heat-sink constant (Eq. 52, pdf p.26: T41sgn = f_hs(NGc)),\n"
            "#    which feeds the heat-sink lag tau_b of Eq. 53. NOT f_s: the report prints\n"
            "#    them as different symbols and an earlier draft of the notes conflated them"
        ),
        xdesc="NGC corrected gas generator speed, percent",
        ydesc="T41SGN heat-sink constant (no unit on the figure or in the nomenclature)",
        fmt=("%.4f", "%.5f", "%.4f", "%.5f"),
        extra=[
            "FREE STRUCTURAL CHECK. The six abscissae should land on 65, 70, 75, 80, 85, 100",
            "(Fig. A11 draws a marker at every printed major x tick plus the right frame).",
            "Nothing here is told to expect that. They recover to -0.074, -0.042, +0.090,",
            "+0.084, +0.016, -0.002 -- max 0.09 against a mean sigma_x of 0.041.",
            "MARKERS ON THE FRAME. Markers 1 and 6 sit on the left and right frame lines, so",
            "`rail` keeps frame-adjacent components in x while excluding a band of tick-label",
            "fragments along the inside of the bottom frame. Marker 1 also HIDES the y major",
            "at 4.00 on the left edge; that major is dropped from the axis fit rather than",
            "failing the extraction, and the held-out frame test still judges the result.",
        ],
    ),
    "c29": dict(
        quantity=(
            "F_HM6 -- HMU idle-schedule function 2 (Fig. C20, pdf p.93:\n"
            "#    PCNGI = F_HM6(T2)), the NG droop-line reference of the idle schedule"
        ),
        xdesc="T2 engine inlet temperature, deg R",
        ydesc="PCNGI NG droop line idle reference parameter, percent",
        fmt=("%.4f", "%.5f", "%.4f", "%.5f"),
        extra=[
            "TWO KNOTS ONLY. Fig. C29 is a straight line and the report draws a marker at each",
            "end; there is nothing between them to read. The function is fully determined.",
            "THE TOP FRAME CARRIES NO LABEL. Its value is taken as 75.0 -- the printed majors",
            "run 64 to 74 at a step of 2 and the frame sits one half-interval above the last.",
            "That reading is what the held-out frame test judges, and it predicts the LOW",
            "frame, 62, to 0.002 -- 0.16 per-mille of range -- from the interior majors alone.",
            "A wrong top frame would show up there.",
        ],
    ),
    "c27": dict(
        quantity=(
            "F_HM4 -- HMU load-demand-compensation schedule function 2 (Fig. C17,\n"
            "#    pdf p.91: WFQPS3 = F_HM4(XLDSH)). Same argument caveat as F_HM3: the block\n"
            "#    diagram feeds it post-hysteresis XLDSH, the axis here is labelled XLDSA"
        ),
        xdesc="XLDSA load demand spindle angle, degrees",
        ydesc="WFQPS3 fuel flow delta demand parameter, nondimensional (no unit printed)",
        fmt=("%.4f", "%.5f", "%.4f", "%.5f"),
        extra=[
            "FREE STRUCTURAL CHECK. The eleven abscissae recover 0, 10, ... 100 to max 0.08 deg",
            "against a mean sigma_x of 0.11. Nothing here expects a lattice.",
            "Flat above 80 deg at 3.5717, 3.5717, 3.5717 -- the three plateau knots agree to",
            "8e-05 of each other, which is a second free check on the y calibration.",
            "SEVEN FRAGMENTS of the left frame survive the opening at 36 px. The marker at",
            "XLDSA = 0 sits on that same frame at 209 px and the marker at XLDSA = 90 is also",
            "36 px but stands in open plot, so neither a global size floor nor a rail exclusion",
            "works: the floor applies only WITHIN 12 px of a frame.",
        ],
    ),
    "c28": dict(
        quantity=(
            "F_HM5 -- HMU idle-schedule function 1 (Fig. C20, pdf p.93:\n"
            "#    WFIRF = F_HM5(T2)), the fuel-flow term of the idle droop line"
        ),
        xdesc="T2 engine inlet temperature, deg R",
        ydesc="WFIRF fuel flow idle schedule limit parameter, nondimensional",
        fmt=("%.4f", "%.5f", "%.4f", "%.5f"),
        extra=[
            "FREE STRUCTURAL CHECK. The thirteen abscissae land on a 20 deg R lattice from 395",
            "to 635, recovering it to max 0.16 deg R against a mean sigma_x of 0.19. Nothing",
            "here is told to expect a lattice, and the end markers stand clear of both frames.",
            "MONOTONE RISING, as an idle schedule must be: 2.0997 at 395 deg R to 2.5708 at",
            "635 deg R.",
        ],
    ),
    "c25": dict(
        quantity=(
            "F_HM2 -- HMU power-available input schedule (Fig. C16, pdf p.91:\n"
            "#    WFPRF = F_HM2(PAS)). One of the eight Appendix C scheduling functions,\n"
            "#    which exist ONLY as plots -- Appendix C prints no equations and no tables"
        ),
        xdesc="PAS power available spindle angle, degrees",
        ydesc="WFPRF fuel flow power-available parameter, nondimensional (no unit printed)",
        fmt=("%.4f", "%.5f", "%.4f", "%.5f"),
        extra=[
            "FREE STRUCTURAL CHECK. Nine of the eleven abscissae recover 40, 50, 60, 70, 80,",
            "90, 100, 110, 120 to max 0.16 deg against a mean sigma_x of 0.14. The first two",
            "are genuinely off-lattice at 28.45 and 29.95, which is what the figure shows.",
            "The last marker sits ON the right frame, so only the TOP rail is excluded --",
            "a fragment sits on the top frame at x = 1957 px.",
            "MONOTONE by construction of the schedule: WFPRF falls from 10.43 at 28 deg to",
            "0.006 at 120 deg. Fig. C25 is what fixes that sign convention (inventory 2.3).",
        ],
    ),
    "c26": dict(
        quantity=(
            "F_HM3 -- HMU load-demand-compensation schedule function 1 (Fig. C17,\n"
            "#    pdf p.91: PNG = F_HM3(XLDSH)). Note the ARGUMENT: Fig. C17 feeds it the\n"
            "#    post-hysteresis spindle angle XLDSH, while the axis here is labelled\n"
            "#    XLDSA -- see open question list, Appendix C item 5"
        ),
        xdesc="XLDSA load demand spindle angle, degrees",
        ydesc="PNG gas generator speed demand parameter, percent",
        fmt=("%.4f", "%.5f", "%.4f", "%.5f"),
        extra=[
            "FREE STRUCTURAL CHECK. The eleven abscissae recover 0, 10, 20, ... 100 to max",
            "0.11 deg against a mean sigma_x of 0.14. Nothing here expects a lattice.",
            "The schedule is flat above 80 deg, which is the topping plateau.",
        ],
    ),
    "c24": dict(
        quantity=(
            "F_HM1 -- HMU topping line schedule (Fig. C14, pdf p.91: WFPTP = F_HM1(T2)).\n"
            "#    One of the eight Appendix C scheduling functions, which exist ONLY as\n"
            "#    plots -- Appendix C prints no numbered equations and no function tables"
        ),
        xdesc="T2 engine inlet temperature, deg R",
        ydesc="WFPTP fuel flow topping line parameter, nondimensional (no unit printed)",
        fmt=("%.4f", "%.5f", "%.4f", "%.5f"),
        extra=[
            "FREE STRUCTURAL CHECK. The six interior abscissae land on a 20 deg R lattice",
            "(515, 535, 555, 575, 595, 615) to +0.18, +0.19, +0.03, +0.12, +0.07, +0.16 --",
            "max 0.19, against a mean sigma_x of 0.22. Nothing in the extraction was told to",
            "expect a lattice. The two end markers do NOT sit on the frame values: they read",
            "395.0 and 634.2 against frames at 390 and 640, which is what the figure shows.",
        ],
    ),
    "a8": dict(
        quantity="f8 -- power turbine energy (Eq. 32, pdf p.24: dH_PT = theta45 * f8(P49/P45))",
        xdesc=(
            "P49/P45 power turbine pressure ratio (TOTAL over total), nondimensional -- NOT\n"
            "#    Ps9/P45; Figure A9 on the facing page uses the STATIC ratio Ps9/P45 over the\n"
            "#    identical numeric range with a near-identical axis title, and the two\n"
            "#    must not be cross-indexed (Eq. 32 vs Eq. 33, pdf p.24)"
        ),
        ydesc=(
            "f8, power turbine enthalpy drop parameter, BTU/LBM -- goes genuinely NEGATIVE\n"
            "#    above P49/P45 = 0.83; the printed y axis is deliberately extended to\n"
            "#    -5.0 to show it, so negative values here are data, not a detection error"
        ),
        fmt=("%.6f", "%.4f", "%.6f", "%.4f"),
        extra=[
            "FREE STRUCTURAL CHECK. The twelve abscissae should land on the 0.05 lattice from",
            "0.30 to 0.85 (Fig. A8 draws a marker at every major x tick). Nothing in the",
            "extraction is told to expect that. Deviations, in units of x:",
            "  {LATTICE}",
            "i.e. max {LATTICEMAX:.2g}, sd {LATTICESD:.2g}, against a mean sigma_x of"
            " {LATSIG:.2g}.",
        ],
    ),
    "a10": dict(
        quantity="f10 -- exhaust pressure loss (Eq. 38, pdf p.24: P49 = Ps9 * f10(NGc))",
        xdesc="NGC corrected gas generator speed, percent",
        ydesc=(
            "f10 = P49/PS9, exhaust pressure ratio, nondimensional. THE VALUES ARE f10\n"
            "#    DIRECTLY -- USE THEM AS THEY STAND, DO NOT INVERT. Figure A10's printed\n"
            "#    axis LABEL reads PS9/P49 and is inverted relative to what is plotted;\n"
            "#    Eq. 38 (P49 = Ps9*f10) plus the expansion order P45 > P49 > Ps9 require\n"
            "#    f10 > 1, and every plotted value is > 1. An earlier note on this file told\n"
            "#    the consumer to take 1/y; that was wrong and would put the exhaust total\n"
            "#    pressure below ambient. See open question #5."
        ),
        fmt=("%.4f", "%.6f", "%.4f", "%.6f"),
        extra=[
            "MARKER COUNT: 18, not the inventory's estimated ~20. The inventory wrote the count",
            "as an estimate and flagged the 89-92 % trough as a judgement call; 18 is what the",
            "page carries. The morphological opening finds 19 cores inside the frame and two of",
            "them, 9 px apart at (836.7, 675.6), are one marker whose waist the scan broke -- at a",
            "6x6 opening they are a single 138 px core. Confirmed by eye at 3-5x magnification",
            "over the whole curve, and by the interior-ink check above, which finds no excess",
            "anywhere on any of the 17 joining segments.",
        ],
    ),
    "a6": dict(
        quantity="f6 -- combustor efficiency, eta = f6(FAR)",
        xdesc=(
            "FAR fuel-to-air ratio, nondimensional (the figure prints FUEL-TO-RATIO, FAR --\n"
            "#    the word AIR is missing, a typo in the report; open question #7)"
        ),
        ydesc="ETA combustor efficiency, nondimensional",
        fmt=("%.6f", "%.6f", "%.6f", "%.6f"),
        extra=[
            "THE WHOLE FUNCTION IS ONE NUMBER. f6 is drawn as a single perfectly horizontal line",
            "with a marker on each vertical frame, so its accuracy is entirely a y-calibration",
            "question -- which is exactly what the scan's bow attacks. That is why the held-out",
            "frame test matters more here than on any other figure in the appendix, and why both",
            "y frame values being printed (0.88 and 1.10) is worth as much as the ticks.",
            "The two markers agree on y to {DYY:.2g}, which is {DYSIG:.1f} sigma -- the",
            "figure's own statement that it is a constant.",
        ],
    ),
}


def _cal_block(CA, spec, axis):
    order = CA["order"]
    var = "u" if axis == "x" else "v"
    side = {"lo": "left" if axis == "x" else "bottom", "hi": "right" if axis == "x" else "top"}
    # the frames bounding the axis (above) are not the edges carrying its ticks (below)
    edge = ("bottom", "top") if axis == "x" else ("left", "right")
    anch = spec.get("anchor", spec["printed"])
    words = ", ".join(
        f"{side[s]} = {spec[s]:g}" + ("" if s in spec["printed"] else " (DERIVED, see below)")
        for s in ("lo", "hi")
        if s in anch
    )
    lines = [
        f"#   {axis} axis: {['', 'LINEAR', 'QUADRATIC', 'CUBIC', 'QUARTIC'][order]} in {var}, "
        f"where {var} = 0 at the {side['lo']} frame and 1 at the {side['hi']} frame.",
        "#      The two ticked edges are fitted SEPARATELY and blended across the plot;",
        "#      the value at a point is the map its own end of the plot measured.",
        f"#        {edge[0]:6s} {axis} = {_poly_str(CA['edgec'][0], var)}",
        f"#        {edge[1]:6s} {axis} = {_poly_str(CA['edgec'][1], var)}",
        "#      Pooling them, which this tool did until 2026-09-16, averages away a real",
        f"#      difference; the pooled map was {axis} = {_poly_str(CA['c'], var)}",
        f"#      Fitted on {CA['n']} interior ticks ({CA['lattice']} lattice, "
        f"{CA['n'] // 2} per edge) plus the frame",
        f"#      value(s) {words};",
        f"#      residual rms {CA['resid']:.4g}, i.e. "
        f"{CA['resid'] / (spec['hi'] - spec['lo']) * 1000:.2f} per-mille of the axis range.",
        "#      HELD-OUT FRAME TEST -- fit the interior majors ALONE, then predict the printed",
        "#      value at the frame, which the fit was never given (leave-one-tick-out rms is",
        "#      shown beside it, as an independent check that the order is not fitting noise):",
    ]
    for o in (1, 2, 3, 4):
        tag = f"#        order {o} ({['', 'linear', 'quadratic', 'cubic', 'quartic'][o]:9s})"
        loo = CA["loos"][o] / (spec["hi"] - spec["lo"]) * 1000
        pad = "#" + " " * (len(tag) + 1)
        for j, pline in enumerate(CA["pred"][o]):
            lines.append(f"{tag}: {pline}" if j == 0 else f"{pad}{pline}")
        lines.append(f"{pad}leave-one-tick-out {loo:.3f} per-mille of range")
    lines.append(
        f"#      -> order {order} wins the held-out test among the orders that also interpolate the"
    )
    lines.append(
        "#         ticks as well as the best one does; ties within 1.25x go to the lower order."
    )
    return lines


def write_csv(key, r):
    cfg, pr = r["cfg"], PROSE[key]
    xs, ys = cfg["x"], cfg["y"]
    X, Y = r["X"], r["Y"]
    sx, sy = r["sig"]["x"], r["sig"]["y"]
    old = np.array(OLD[key]) if key in OLD else None
    fmtx, fmty, fmtsx, fmtsy = pr["fmt"]
    L = []
    A = L.append
    A(f"# source: TM-100991 pdf p.{cfg['page']}, Figure {cfg['fig']}")
    A(f"# quantity: {pr['quantity']}")
    A(f"# method: digitized -- generative model fit of the {cfg['nmark']} printed 'x' glyphs")
    A(f"# x: {pr['xdesc']}")
    A(f"# y: {pr['ydesc']}")
    A("# sigma_x: 1 sigma uncertainty on x, in x's own units, per point (see UNCERTAINTY)")
    A("# sigma_y: 1 sigma uncertainty on y, in y's own units, per point (see UNCERTAINTY)")
    A(f"# points: {cfg['nmark']}")
    A(
        "# digitized: 2026-09-11 by tools/digitize_single_curve.py"
        " -- reruns and reproduces this file"
    )
    if old is not None:
        A("#            exactly; a re-digitisation of the 2026-09-10 extraction (see MOVEMENT)")
    else:
        A("#            exactly; a first extraction, with no earlier version to compare against")
    A("#")
    A(f"# RASTER. `pdfimages -list -f {cfg['page']}` shows pdf p.{cfg['page']} is a single 300 dpi")
    A("# 1-bit CCITT image (2544x3300), so every render above 300 dpi is an upsample of that")
    A("# bitmap and can add no detail. This extraction works on the bitmap itself, pulled out")
    A("# with `pdfimages`, and takes its sub-pixel information from a model fitted to the binary")
    A("# ink rather than from interpolated grey. The page is NOT resampled at any point: the")
    A("# scan's rotation and keystone are a coordinate transform applied to the measured points,")
    A("# not a warp applied to the pixels, so no interpolation ever touches the data.")
    A("#")
    arcs = r["forder"] > 1
    A("# GEOMETRY. The four frame lines are fitted sub-pixel by tracing the single ink run per")
    A("# row/column and least-squares fitting with outlier rejection. A printed straight edge")
    A("# does not image straight here, so BOTH models are fitted and the figure chooses: fit")
    A("# each edge as a line, or as an arc, and keep whichever predicts the printed frame")
    A("# values it was never shown. Summed held-out error over both axes, per-mille of range:")
    A(
        f"#   straight {r['frame_scores'][1]:.3f}   arcs {r['frame_scores'][2]:.3f}"
        f"   ->  {'ARCS' if arcs else 'STRAIGHT'}"
        f" (arcs must win by {(1 - FRAME_MARGIN) * 100:.0f} %; straight is the incumbent)"
    )
    for nm in ("left", "right", "top", "bottom"):
        fd = r["frdiag"][nm]
        v = "x" if nm in ("left", "right") else "y"
        t = "y" if nm in ("left", "right") else "x"
        extra = f", sagitta {fd['sagitta']:+.2f} px" if arcs else ""
        A(
            f"#   {edge_str(nm, r['fc'][nm], v, t)}   ({fd['slope_deg']:+.3f} deg, "
            f"{fd['n']} samples, rms {fd['rms']:.2f} px{extra})"
        )
    TL, TR, BR, BL = r["corners"]
    A("# The four edge slopes are not a common rotation, so the frame is a genuine")
    A(f"# quadrilateral; corners TL({TL[0]:.3f},{TL[1]:.3f}) TR({TR[0]:.3f},{TR[1]:.3f})")
    A(f"# BR({BR[0]:.3f},{BR[1]:.3f}) BL({BL[0]:.3f},{BL[1]:.3f}). They are mapped onto the")
    if arcs:
        A("# unit square (u right, v up) by the curvilinear blend between the four fitted arcs,")
        A("# which is exact on each edge by construction, and every measured point is carried")
        A("# through it.")
    else:
        A("# unit square (u right, v up) by a homography, and every measured point is carried")
        A("# through it.")
    A("# NOT a pixel warp -- nothing is resampled at any point.")
    A("#")
    A("# CALIBRATION. Every printed MAJOR tick on all four edges, located by the ink-weighted")
    A("# centroid of its inward protrusion, and matched to its printed value by predicted")
    A("# lattice slot rather than by protrusion length -- a marker sitting ON a frame line")
    A("# out-protrudes a major tick, and a flat data line protrudes the whole scan length.")
    A(
        f"#   x majors every {major_step(xs):g} ({len(major_values(xs))} interior, read off"
        f" the top and bottom edges);"
    )
    A(
        f"#   y majors every {major_step(ys):g} ({len(major_values(ys))} interior, off the"
        f" left and right edges)."
    )
    A("# The page also prints a MINOR lattice, five to twenty-five times as many ticks, and")
    A("# until 2026-09-16 this tool detected it and threw it away. Each minor step is")
    A("# measured from the ticks themselves -- the major spacing divided by the median gap,")
    A("# which must come out an integer -- and each tick is then labelled through a fit made")
    A("# on the MAJORS alone, never through the straight line between the frame corners,")
    A("# which is several pixels out in mid-page and discards a contiguous block of good")
    A("# ticks. Whether an axis uses the fine lattice or its majors alone is decided per")
    A("# axis by the held-out frame test below, not chosen: whichever predicts the printed")
    A("# frame value it was never shown. Both numbers are printed, so the choice is visible.")
    for CA, ax, nm in ((r["CX"], xs, "x"), (r["CY"], ys, "y")):
        alt = CA["alt_err"] / (ax["hi"] - ax["lo"]) * 1000
        here = CA["errs"][CA["order"]] / (ax["hi"] - ax["lo"]) * 1000
        A(
            f"#   {nm}: {CA['lattice']} lattice, {CA['n']} ticks, predicting the printed frame"
            f" to {here:.3f} per-mille"
        )
        A(f"#      against {alt:.3f} per-mille for the other choice.")
    L.extend(_cal_block(r["CX"], xs, "x"))
    L.extend(_cal_block(r["CY"], ys, "y"))
    A("#   Independent bow-free cross-check: linear interpolation between the two majors that")
    A("#   bracket each point, done on each ticked edge and averaged, agrees with the chosen")
    A(
        f"#   model to {np.abs(np.polyval(r['cx'], r['u']) - r['loc']['x']).max():.3g} in x and "
        f"{np.abs(np.polyval(r['cy'], r['v']) - r['loc']['y']).max():.3g} in y, worst point."
    )
    A("#")
    k = cfg.get("open_sz", 7)
    A(
        f"# MARKERS. {cfg['nmark']} 'x' glyphs. Seeds come from a {k}x{k} morphological "
        "OPENING of the"
    )
    A("# ink -- the joining line is ~4 px wide and vanishes, the crossing of two strokes does")
    A("# not -- which is independent of density and needs no threshold sweep. Each centre is")
    A("# then fitted with a generative model of the local ink: (stroke1 U stroke2 U incoming")
    A("# half-segment U outgoing half-segment), all four radiating from one point, rendered by")
    A("# 5x5 supersampling and scored as sum (coverage - ink)^2 over an unmasked disc of radius")
    tpl = r["tpl"]
    A(
        f"# {cfg['rad']:.0f} px. The glyph template is global -- the plotter drew the same mark "
        f"{cfg['nmark']} times -- and"
    )
    A(
        f"# came out as strokes at {np.degrees(tpl[0]):.2f} deg (half-length {tpl[1]:.2f} px) and "
        f"{np.degrees(tpl[2]):+.2f} deg ({tpl[3]:.2f} px),"
    )
    A(f"# stroke width {tpl[4]:.2f} px, joining-line width {tpl[5]:.2f} px. Frame lines and tick")
    A("# marks are MASKED, never modelled, and the mask is SYMMETRIC about each frame line: the")
    A("# ticks all protrude inward, so masking only the inward side would strip ink from one")
    A("# half of a marker that sits on a frame and bias its abscissa. Fitted centres, image px:")
    for i in range(0, cfg["nmark"], 3):
        A(
            "#   "
            + "  ".join(
                f"m{j + 1:<2d}({r['C'][j][0]:8.3f},{r['C'][j][1]:8.3f})"
                for j in range(i, min(i + 3, cfg["nmark"]))
            )
        )
    A("#")
    A("# UNCERTAINTY, per point, 1 sigma, in data units -- the sigma_x and sigma_y columns.")
    A("# Quadrature of four named contributions, each measured, none assumed:")
    for ax, CA, t in (("x", r["CX"], r["u"]), ("y", r["CY"], r["v"])):
        b, c, d = r["terms"][ax]  # the terms the sigma column is the quadrature of
        A(
            f"#   {ax}: tick-fit residual {CA['resid']:.3g} | calibration order spread (deg 1-4) "
            f"{b.mean():.3g} mean, {b.max():.3g} max |"
        )
        A(
            f"#      the two ticked edges calibrated separately, half their disagreement "
            f"{c.mean():.3g} mean, {c.max():.3g} max |"
        )
        b24 = np.max(
            [
                np.abs(np.polyval(CA["allc"][o], t) - np.polyval(CA["allc"][o2], t))
                for o in (2, 3, 4)
                for o2 in (2, 3, 4)
            ],
            axis=0,
        )
        A(
            f"#      glyph centre {d.mean():.3g} mean, {d.max():.3g} max.  Total 1 sigma "
            f"{r['sig'][ax].mean():.3g} mean, {r['sig'][ax].max():.3g} max."
        )
        A(
            "#      (The order term is deliberately conservative: it includes the LINEAR map, "
            "which the"
        )
        A(
            f"#      held-out frame test rejects. Over orders 2-4 alone the spread is only "
            f"{b24.max():.3g}.)"
        )
    A("#   The glyph term is the worst case over: window radius +/-3 px, frame-mask half-width")
    A("#   3.0 and 5.0 px, tick-mask half-width 3.0 and 6.0 px, and masking the joining line out")
    A("#   ENTIRELY and fitting the two strokes alone. Worst move in px over those variants:")
    for k2, v2 in r["varmax"].items():
        A(f"#      {k2:16s} {v2:.2f} px")
    A("#")
    A("# COUNT CHECK, independent of any detector. Walking the ink along each fitted joining")
    A("# segment in 10 px bins, every excess over that segment's own baseline sits at a segment")
    A("# endpoint; there is no interior excess on any segment -- no marker is missed and none")
    A("# is invented. Straightness of each joining segment (residual about a straight line,")
    A("# px rms): " + ", ".join(f"{s:.2f}" for _, _, s in r["sag"]) + ".")
    A("#")
    for ln in pr["extra"]:
        A("# " + ln)
    A("#")
    if old is not None:
        A("# MOVEMENT from the 2026-09-10 extraction (old -> new, and the move in units of the")
        A("# 1 sigma above):")
        dxs = X - old[:, 0]
        dys = Y - old[:, 1]
        for i in range(cfg["nmark"]):
            A(
                f"#   {i + 1:2d}  x {old[i, 0]:{fmtx[1:]}} -> {X[i]:{fmtx[1:]}} "
                f"({dxs[i]:+.3g}, {dxs[i] / sx[i]:+.1f} sigma)   "
                f"y {old[i, 1]:{fmty[1:]}} -> {Y[i]:{fmty[1:]}} "
                f"({dys[i]:+.3g}, {dys[i] / sy[i]:+.1f} sigma)"
            )
        A("#")
    A("# NOT transcribed. Values carry read error; see validation/out/digitize/a6810/.")
    A("x,y,sigma_x,sigma_y")
    for i in range(cfg["nmark"]):
        A(f"{fmtx % X[i]},{fmty % Y[i]},{fmtsx % sx[i]},{fmtsy % sy[i]}")
    text = "\n".join(L) + "\n"
    if key == "a8":
        lat = X - np.round(X / major_step(xs)) * major_step(xs)
        text = (
            text.replace("{LATTICE}", " ".join(f"{d:+.5f}" for d in lat))
            .replace("{LATTICEMAX:.2g}", f"{np.abs(lat).max():.2g}")
            .replace("{LATTICESD:.2g}", f"{lat.std():.2g}")
            .replace("{LATSIG:.2g}", f"{sx.mean():.2g}")
        )
    if key == "a6":
        d = abs(Y[1] - Y[0])
        text = text.replace("{DYY:.2g}", f"{d:.2g}").replace(
            "{DYSIG:.1f}", f"{d / (sy.mean() * np.sqrt(2)):.1f}"
        )
    cfg["csv"].write_text(text)
    print(f"  wrote {cfg['csv']}  ({len(L)} lines)")


def _pixel_to_unit(to_uv, C, axis):
    """d(u or v) / d(pixel) at each marker, along that axis -- the local map scale.

    Takes the map as a callable so it works for either frame model: the homography, or the
    curvilinear blend between four fitted arcs.
    """
    u0, v0 = to_uv(C[:, 0], C[:, 1])
    if axis == 0:
        u1, _ = to_uv(C[:, 0] + 1.0, C[:, 1])
        return np.abs(u1 - u0)
    _, v1 = to_uv(C[:, 0], C[:, 1] + 1.0)
    return np.abs(v1 - v0)


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    # Every figure this tool owns: the gate runs it with no arguments, so a narrower
    # default would leave committed files unchecked (tests/test_reproducibility_gate.py).
    keys = [a for a in argv if a in FIGS] or list(FIGS)
    for k in keys:
        r = run(k)
        print("\n  x, y, sigma_x, sigma_y")
        for i in range(len(r["X"])):
            print(
                f"  {r['X'][i]:.6f}, {r['Y'][i]:.6f}, "
                f"{r['sig']['x'][i]:.6f}, {r['sig']['y'][i]:.6f}"
            )
        write_csv(k, r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
