"""The single-curve digitizer: thirteen figures across both appendices of TM-100991.

Same method as `tools/digitize_a7.py`, generalised. Every figure here plots one curve or
one marker series per panel, which is what lets the markers be found without reading
anything off the page; the two Appendix C figures whose curves must be told apart live in
`digitize_multi_curve.py` and `digitize_c30.py`, and both import their frame location,
tick lattice and axis calibration from here unchanged.

    A3,  pdf p.58 -> f3    seal-pressurization bleed      data/maps/f3_seal_bleed_fraction.csv
    A4,  pdf p.59 -> f4    PT-balance bleed               data/maps/f4_pt_balance_bleed_fraction.csv
    A5,  pdf p.60 -> f5    tip-leak + cooling bleed       data/maps/f5_tip_leak_cooling_bleed...csv
    A6,  pdf p.61 -> f6    combustor efficiency           data/maps/f6_combustor_efficiency.csv
    A8,  pdf p.63 -> f8    power turbine energy           data/maps/f8_pt_energy.csv
    A10, pdf p.65 -> f10   exhaust pressure loss          data/maps/f10_exhaust_pressure_loss.csv
    A11, pdf p.66 -> f_hs  station 4.1 heat-sink constant data/maps/fhs_heat_sink_constant.csv
    C24-C29, pp.95-100 -> F_HM1..F_HM6                    data/schedules/fhm[1-6]_*.csv

A1, A2, A7 and A9 have their own tools: each needed something this path does not do.

Why A6, A8 and A10 were redone first, which is where the method here came from:
re-digitizing A7 showed that a LINEAR axis map cannot fit these
scans -- the page bows along the scan direction by several pixels -- and that every figure
which prints a value at a frame edge can test its own calibration model against a number
the fit was never given.  All three figures here print a value at at least three of their
four frame edges, so the test is available on every axis.

What this does:

1.  Works on the SCAN.  `pdfimages -list` shows pp.61/63/65 are each one 300 dpi 1-bit
    CCITT image, so any higher-dpi render is an upsample.  The bitmap is pulled out with
    `pdfimages` and NEVER resampled: rotation and keystone are carried as a homography
    applied to the measured points, so no interpolation touches the data.
2.  Fits the four frame lines sub-pixel and calibrates on every printed MAJOR tick on all
    four edges.  Ticks are matched to their printed values by predicted lattice slot, not
    by a protrusion-length threshold -- a marker sitting on a frame line out-protrudes a
    major tick, and on A6 the flat data line protrudes all the way across.
3.  Chooses the polynomial order per axis by the HELD-OUT FRAME TEST: fit the interior
    majors only, then predict the printed value at the frame.  Reports orders 1-4 against
    printed, guards the choice with leave-one-tick-out on the ticks themselves (without it
    a high order wins A10's x axis by fitting noise), and refits from each edge separately
    to check they agree.  The order that wins is NOT the same on any two of the six axes.
4.  Locates the markers by fitting a generative model of the local ink -- two glyph
    strokes plus the joining half-segments, all radiating from one point -- scored against
    ink AND background over an unmasked disc.  Frame lines and tick marks are masked
    (symmetrically about the frame line, so the mask cannot bias the abscissa), never
    modelled.  Seeds come from a morphological opening, which is independent of density.

Run from the repo root.  Writes the thirteen CSVs and overlays under
validation/out/digitize/.  Needs `pdfimages` (poppler).

This module lives in tools/ and may use SciPy and PIL.  Nothing here is imported by
src/t700/.
"""

from __future__ import annotations

import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
from scipy.optimize import minimize

PDF = Path("docs/ballin-tm100991.pdf")
OUT = Path("validation/out/digitize/a6810")

# --------------------------------------------------------------------------- frame finding


def _longest_run(mask, gap=8):
    """Longest run of True in a 1-D boolean array, tolerating gaps of up to `gap`.

    Gap tolerance is not a nicety. The frame lines on several of these pages are broken
    by the scan: Figure A3's horizontal frames have a longest *strictly* contiguous run of
    777 px against a 1666 px frame, so a zero-gap measure does not see them at all.
    """
    m = np.ascontiguousarray(mask).astype(np.int8)
    if not m.any():
        return 0
    if gap > 0:
        m = ndimage.binary_closing(m.astype(bool), structure=np.ones(2 * gap + 1, bool))
        m = m.astype(np.int8)
    d = np.diff(np.concatenate(([0], m, [0])))
    return int((np.flatnonzero(d == -1) - np.flatnonzero(d == 1)).max())


def _run_extent(mask, gap=8):
    """(start, end, length) of the longest gap-tolerant True run in a 1-D array."""
    m = np.ascontiguousarray(mask).astype(bool)
    if not m.any():
        return (0, 0, 0)
    if gap > 0:
        m = ndimage.binary_closing(m, structure=np.ones(2 * gap + 1, bool))
        if not m.any():  # closing erodes a run that only touches the border
            return (0, 0, 0)
    d = np.diff(np.concatenate(([0], m.astype(np.int8), [0])))
    a = np.flatnonzero(d == 1)
    b = np.flatnonzero(d == -1)
    k = int(np.argmax(b - a))
    return int(a[k]), int(b[k]), int(b[k] - a[k])


def _sheared_row_profile(sub, slope):
    """Row-sum of `sub` after shearing each column to undo a line of the given slope."""
    h, w = sub.shape
    xc = w / 2.0
    shift = np.rint(slope * (np.arange(w) - xc)).astype(int)
    ys = np.arange(h)[:, None] + shift[None, :]
    ok = (ys >= 0) & (ys < h)
    return np.where(ok, sub[np.clip(ys, 0, h - 1), np.arange(w)[None, :]], False).sum(axis=1)


def find_frame_box(ink, margin=150, sep_frac=0.25, slopes=None):
    """Locate the printed plot box as (left, right, top, bottom) pixel centres.

    The horizontal frames on these scans are **tilted and broken**, and that combination
    defeats every row-wise measure. A11's bottom frame descends ten rows across the plot
    -- 0.36 degrees -- so no single row holds more than about half of it, and simple row
    sums pick the x-axis label band or a caption rule instead. That is how A11 acquired a
    bottom frame 192 px too low, which put every digitized value 0.3 too high while the
    tick lattice still appeared to close.

    So the rows are summed **after shearing** the page by a trial slope, and the slope that
    maximises the response is the one the line actually has. Verticals need no such
    treatment: they are solid, and their own tilt is small over the page width.
    """
    if slopes is None:
        slopes = np.arange(-14, 15) * 0.001
    h, w = ink.shape

    def two(profile, lo, hi, minsep):
        prof = profile.astype(float).copy()
        prof[:lo] = 0.0
        prof[hi:] = 0.0
        out = []
        for _ in range(2):
            i = int(np.argmax(prof))
            if prof[i] <= 0:
                raise AssertionError("no frame candidate -- check the margin")
            out.append(i)
            prof[max(0, i - minsep) : i + minsep] = 0.0
        return sorted(out)

    lef, rig = two(ink.sum(axis=0), margin, w - margin, int(sep_frac * w))
    sub = ink[:, lef + 10 : rig - 10]

    def outermost(prof, cover):
        """Extreme rows that span the plot, not the strongest rows.

        Two failures shaped this. Taking the two *strongest* rows picks Figure A6's flat
        combustor-efficiency curve instead of its top frame, because that curve is solid
        while the frame is broken. Taking the extreme rows above a fraction of the *max*
        then drops A6's top frame entirely, for the same reason -- the solid data line
        sets the max.

        Thresholding on coverage of the plot WIDTH fixes both: a footer rule is too short
        to qualify, a data line qualifies but lies inside, and a frame -- however broken --
        still spans the plot once the shear has been undone.
        """
        q = prof.astype(float).copy()
        q[:margin] = 0.0
        q[h - margin :] = 0.0
        idx = np.flatnonzero(q >= cover * sub.shape[1])
        if idx.size == 0 or idx.max() - idx.min() < sep_frac * h:
            return None
        return int(idx.min()), int(idx.max())

    best = None
    for m in slopes:
        prof = _sheared_row_profile(sub, float(m))
        got = outermost(prof, 0.60)
        if got is None:
            continue
        a, b = got
        score = float(prof[a]) + float(prof[b])
        if best is None or score > best[0]:
            best = (score, a, b, float(m))
    if best is None:
        raise AssertionError("no horizontal frame candidates")
    _, top, bot, _m = best
    return float(lef), float(rig), float(top), float(bot)


def auto_frame_spec(page, pad=45, inset=15):
    """`frame_spec` for a page, derived from its own bitmap. See `find_frame_box`."""
    lef, rig, top, bot = find_frame_box(native_bitmap(page))
    return [
        ("left", 0, int(lef - pad), int(lef + pad), int(top + inset), int(bot - inset)),
        ("right", 0, int(rig - pad), int(rig + pad), int(top + inset), int(bot - inset)),
        ("top", 1, int(top - pad + 10), int(top + pad - 10), int(lef + inset), int(rig - inset)),
        ("bottom", 1, int(bot - pad + 10), int(bot + pad - 10), int(lef + inset), int(rig - inset)),
    ]


# --------------------------------------------------------------------------- figure table
#
# `frame_spec`  windows in which to trace each frame line: (name, along, xlo, xhi, tlo, thi)
#               `along=0` traces a vertical line (one ink run per row), `along=1` horizontal.
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


# --------------------------------------------------------------------------- the raster


def native_bitmap(page: int) -> np.ndarray:
    """The page's own 300 dpi 1-bit image, as an ink mask.  No render, no resample."""
    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT / f"p{page}"
    if not (OUT / f"p{page}-000.png").exists():
        subprocess.run(
            ["pdfimages", "-f", str(page), "-l", str(page), "-png", str(PDF), str(stem)],
            check=True,
        )
    return ~np.array(Image.open(OUT / f"p{page}-000.png"))


# --------------------------------------------------------------------------- frame lines


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    e = np.flatnonzero(np.diff(np.concatenate(([0], mask.view(np.int8), [0]))))
    return list(zip(e[0::2], e[1::2], strict=True))


def _trace(ink, along, lo, hi, tlo, thi, maxw=8):
    """Centroid of the single narrow ink run in [lo,hi) at each position along the line."""
    ts, vs = [], []
    for t in range(tlo, thi):
        seg = ink[t, lo:hi] if along == 0 else ink[lo:hi, t]
        rr = [r for r in _runs(seg) if (r[1] - r[0]) <= maxw]
        if len(rr) != 1:
            continue
        a, b = rr[0]
        ts.append(float(t))
        vs.append(lo + (a + b - 1) / 2.0)
    return np.array(ts), np.array(vs)


def _robust(t, v, tol=2.0, it=8):
    keep = np.ones(t.shape, bool)
    for _ in range(it):
        s, i = np.polyfit(t[keep], v[keep], 1)
        new = np.abs(v - (i + s * t)) < tol
        if new.sum() < 20 or (new == keep).all():
            break
        keep = new
    s, i = np.polyfit(t[keep], v[keep], 1)
    return float(i), float(s), (v - (i + s * t))[keep], int(keep.sum())


FRAME_DIAG: dict[str, tuple[float, float, float, int]] = {}


def frame(ink, spec) -> dict[str, tuple[float, float]]:
    """The four frame lines as (intercept, slope): verticals x(y), horizontals y(x)."""
    out = {}
    FRAME_DIAG.clear()
    for nm, along, xlo, xhi, tlo, thi in spec:
        a, b, r, n = _robust(*_trace(ink, along, xlo, xhi, tlo, thi))
        out[nm] = (a, b)
        FRAME_DIAG[nm] = (a, b, float(r.std()), n)
        print(
            f"  {nm:6s} {a:10.3f} + {b:+.6f}*t   n {n:5d}  rms {r.std():.3f} px  "
            f"({np.degrees(np.arctan(b)):+.3f} deg)"
        )
    return out


def corner(hor, ver):
    a, b = hor
    c, d = ver
    x = (c + d * a) / (1.0 - d * b)
    return x, a + b * x


def homography(src, dst):
    A, B = [], []
    for (sx, sy), (dx, dy) in zip(src, dst, strict=True):
        A.append([sx, sy, 1, 0, 0, 0, -dx * sx, -dx * sy])
        B.append(dx)
        A.append([0, 0, 0, sx, sy, 1, -dy * sx, -dy * sy])
        B.append(dy)
    h = np.linalg.solve(np.asarray(A, float), np.asarray(B, float))
    return np.append(h, 1.0).reshape(3, 3)


def _apply(H, px, py):
    px = np.atleast_1d(np.asarray(px, float))
    py = np.atleast_1d(np.asarray(py, float))
    d = H[2, 0] * px + H[2, 1] * py + H[2, 2]
    return (
        (H[0, 0] * px + H[0, 1] * py + H[0, 2]) / d,
        (H[1, 0] * px + H[1, 1] * py + H[1, 2]) / d,
    )


# --------------------------------------------------------------------------- ticks

EDGE = {  # name -> (along, inward sign)
    "left": (0, +1),
    "right": (0, -1),
    "top": (1, +1),
    "bottom": (1, -1),
}


def _protrusion(ink, edge, along, sign, lo, hi, maxlen=60):
    """Contiguous ink growing inward from a frame edge, per position along it."""
    a, b = edge
    out = np.zeros(hi - lo, int)
    for j, t in enumerate(range(lo, hi)):
        base = a + b * t
        k = 0
        while k < maxlen:
            q = int(round(base + sign * (3 + k)))
            on = ink[t, q] if along == 0 else ink[q, t]
            if not on:
                break
            k += 1
        out[j] = k
    return out


def _groups(p, lo, thresh, gap=3):
    """Tick candidates: ink-weighted centroid, peak protrusion, width."""
    idx = np.flatnonzero(p >= thresh)
    if idx.size == 0:
        return []
    out, cur = [], [int(idx[0])]
    for v in idx[1:]:
        if v - cur[-1] <= gap:
            cur.append(int(v))
        else:
            out.append(cur)
            cur = [int(v)]
    out.append(cur)
    res = []
    for g in out:
        g = np.asarray(g)
        w = p[g].astype(float)
        res.append((float(lo + (g * w).sum() / w.sum()), float(p[g].max()), len(g)))
    return res


def _pick_majors(cands, t0, t1, frac, name, max_missing=2):
    """Match interior major ticks to their printed values by predicted lattice slot.

    Thresholding on protrusion length does NOT work here.  A marker sitting on a frame
    line out-protrudes a major tick (A6, A8, A10 all print one), A6's flat data line
    protrudes the full scan length, and A10's lowest left-edge major protrudes only as far
    as a minor.  What is reliable is the SLOT: there are exactly `n` interior majors,
    evenly spaced between the two frame lines to within the scan's few-pixel bow.  So
    predict each slot, keep candidates whose protrusion is in family with the majors, and
    take the closest.  Returns one centre per entry of `frac`, in printing order, where
    `frac` is each major's fractional position between the two frames -- not necessarily
    evenly spaced, because Figure A4's majors are not.
    """
    n = len(frac)
    slots = [t0 + (t1 - t0) * f for f in frac]
    step = (t1 - t0) / (n + 1)
    win = 0.25 * abs(step)
    # protrusion length of a major: median of the nearest in-window candidate per slot,
    # over candidates that are neither a minor tick (<5) nor a runaway (>25)
    seed = []
    for s in slots:
        c = [g for g in cands if abs(g[0] - s) < win and 5 <= g[1] <= 25]
        if c:
            seed.append(min(c, key=lambda g: abs(g[0] - s))[1])
    if not seed:
        raise AssertionError(f"{name}: no major-tick candidates at all")
    lmaj = float(np.median(seed))
    # A marker printed ON a frame hides the tick underneath it -- A11's first marker sits
    # at y = 4.00 on the left frame and A7 records the same thing. A hidden major is
    # recorded as NaN and dropped from the fit rather than failing the extraction; the
    # axis map is fitted on the majors that ARE visible, and the held-out frame test still
    # judges it. Losing more than two of them is a real failure and still raises.
    out, missing = [], 0
    for s in slots:
        c = [g for g in cands if abs(g[0] - s) < win and 0.55 * lmaj <= g[1] <= 1.7 * lmaj]
        if not c:
            out.append(float("nan"))
            missing += 1
            continue
        out.append(min(c, key=lambda g: abs(g[0] - s))[0])
    if missing > max_missing:
        raise AssertionError(
            f"{name}: {missing} of {n} majors not found (L~{lmaj:.1f}) -- frame wrong?"
        )
    seen = [v for v in out if not np.isnan(v)]
    if len(set(seen)) != len(seen):
        raise AssertionError(f"{name}: duplicate tick matched")
    resid = np.array(seen) - np.array(
        [s for s, v in zip(slots, out, strict=True) if not np.isnan(v)]
    )
    print(
        f"    {name:6s} {n - missing}/{n} majors, L~{lmaj:.1f} px, "
        f"slot residual {resid.std():.2f} px rms, max {np.abs(resid).max():.2f}"
        + (f"  ({missing} hidden by a marker)" if missing else "")
    )
    return out, lmaj


def major_step(ax):
    """The axis's major-tick interval, whether given as `step` or implied by `majors`."""
    if "step" in ax:
        return float(ax["step"])
    v = np.asarray(ax["majors"], float)
    d = np.diff(v)
    if d.size == 0:
        raise AssertionError("an axis needs either `step` or two or more `majors`")
    if np.ptp(d) > 1e-9 * abs(d[0]):
        raise AssertionError(f"majors are not evenly spaced: {v}")
    return float(d[0])


def major_values(ax):
    """Interior major-tick values of an axis.

    Defaults to `lo + step*(k+1)`, which is right whenever the lattice is anchored on the
    low frame. Figure A4 is not: its frame is at WA2c = 3.0 while its majors are printed
    at 4, 6, 8, 10, so they are neither at `lo + k*step` nor evenly spaced between the
    frames. Such an axis gives `majors` explicitly.
    """
    if "majors" in ax:
        return np.asarray(ax["majors"], float)
    return np.array([ax["lo"] + ax["step"] * (k + 1) for k in range(ax["n"])])


def n_majors(ax):
    """How many interior majors the axis has."""
    return len(major_values(ax))


def ticks(ink, fr, cfg):
    """Every printed major tick on all four edges, plus every tick of any length (masks)."""
    TL = corner(fr["top"], fr["left"])
    TR = corner(fr["top"], fr["right"])
    BL = corner(fr["bottom"], fr["left"])
    BR = corner(fr["bottom"], fr["right"])
    span = {
        "left": (TL[1], BL[1]),
        "right": (TR[1], BR[1]),
        "top": (TL[0], TR[0]),
        "bottom": (BL[0], BR[0]),
    }
    maj, allt = {}, {}
    for nm in ("left", "right", "top", "bottom"):
        along, sign = EDGE[nm]
        a, b = span[nm]
        lo, hi = int(min(a, b)) + 3, int(max(a, b)) - 2
        pr = _protrusion(ink, fr[nm], along, sign, lo, hi)
        cands = _groups(pr, lo, 3)
        allt[nm] = cands
        ax = cfg["y" if nm in ("left", "right") else "x"]
        vals = major_values(ax)
        if nm in ("left", "right"):
            # Pixels run top to bottom, values high to low, so slot 0 is the HIGHEST
            # major. Reverse the values, then measure down from the high frame.
            frac = (ax["hi"] - vals[::-1]) / (ax["hi"] - ax["lo"])
        else:
            frac = (vals - ax["lo"]) / (ax["hi"] - ax["lo"])
        maj[nm], _ = _pick_majors(cands, a, b, frac, nm, cfg.get("max_missing", 2))
    return maj, allt, (TL, TR, BR, BL)


# --------------------------------------------------------------------------- calibration


def _polyfit_predict(t, val, order, at):
    c = np.polyfit(t, val, order)
    return np.polyval(c, at), c


def calibrate_axis(t_lo_edge, t_hi_edge, spec, label):
    """Pick the polynomial order by the held-out frame-edge test, then fit.

    `t_lo_edge` / `t_hi_edge` are (coordinate, value) for the interior majors read off each
    of the two edges that carry ticks for this axis.  Coordinate is u (x axis) or v (y
    axis); 0 at the low-side frame, 1 at the high-side frame.  The printed frame values are
    NOT in the fit used for the test -- that is the point.
    """
    tt = np.concatenate([t_lo_edge[0], t_hi_edge[0]])
    vv = np.concatenate([t_lo_edge[1], t_hi_edge[1]])
    printed = {"lo": spec["lo"], "hi": spec["hi"]}
    at = {"lo": 0.0, "hi": 1.0}
    rng = spec["hi"] - spec["lo"]
    orders = (1, 2, 3, 4)
    print(f"  {label} axis -- held-out frame test (fit = the {tt.size} interior majors only):")
    errs, loos, pred = {}, {}, {}
    for order in orders:
        row, worst = [], 0.0
        for side in spec["printed"]:
            p, _ = _polyfit_predict(tt, vv, order, at[side])
            e = float(p) - printed[side]
            worst = max(worst, abs(e))
            row.append(f"{side} frame printed {printed[side]:g} -> {float(p):.6g} ({e:+.3g})")
        pred[order] = row
        errs[order] = worst
        e = np.array(
            [
                vv[k] - np.polyval(np.polyfit(np.delete(tt, k), np.delete(vv, k), order), tt[k])
                for k in range(tt.size)
            ]
        )
        loos[order] = float(e.std())
        print(
            f"    order {order}: " + " | ".join(row) + f"   [LOO {loos[order] / rng * 1000:.3f} "
            f"per-mille of range]"
        )
    # Selection: the held-out frame value decides, but only among orders that also
    # INTERPOLATE the ticks as well as the best one does -- otherwise a high order can win
    # the extrapolation by accident while fitting noise in between (A10's x axis does
    # exactly that).  Ties within 1.25x go to the lower order; nothing above quartic is
    # considered, because the scan's distortion is a smooth page-scale bow and the frame
    # test itself degrades beyond quartic on all six axes here.
    lbest = min(loos.values())
    ok = [o for o in orders if loos[o] <= 1.20 * lbest]
    order = min(ok, key=lambda k: errs[k])
    for o in sorted(ok):
        if errs[o] <= 1.25 * errs[order]:
            order = o
            break
    print(
        f"    -> chosen order {order} ({['', 'linear', 'quadratic', 'cubic', 'quartic'][order]}); "
        f"worst frame error {errs[order]:.3g} = {errs[order] / rng * 1000:.2f} per-mille of range"
    )
    # each edge on its own, as an independent check of the same prediction
    for nm, d in (("low-side edge", t_lo_edge), ("high-side edge", t_hi_edge)):
        row = []
        for side in spec["printed"]:
            p, _ = _polyfit_predict(d[0], d[1], order, at[side])
            row.append(f"{side} -> {float(p):.6g} ({float(p) - printed[side]:+.3g})")
        print(f"    {nm} alone, order {order}: " + " | ".join(row))

    # final map: majors plus the printed frame values, one copy per ticked edge
    # How many major steps does the frame span, measured by the ticks alone?  On the four
    # axes here where BOTH frame values are printed it comes out within 0.01 of an integer
    # every time, so on the two where one frame is unlabelled (A8 x, A10 y) the integer
    # count fixes that frame's value and it is used as an anchor -- derived, not printed,
    # and flagged as such.
    c1 = np.polyfit(tt, vv, 1)
    span = (np.polyval(c1, 1.0) - np.polyval(c1, 0.0)) / major_step(spec)
    print(
        f"    tick lattice: the frame spans {span:.4f} major steps "
        f"(nearest integer {round(span)}, printed span "
        f"{(spec['hi'] - spec['lo']) / major_step(spec):.0f})"
    )
    anchors = spec.get("anchor", spec["printed"])
    if set(anchors) - set(spec["printed"]):
        print(
            f"    anchoring also at the DERIVED {sorted(set(anchors) - set(spec['printed']))} "
            f"frame value(s), from that integer step count"
        )

    def anchor(t, v):
        ta, va = list(t), list(v)
        for side in anchors:
            ta += [at[side]] * 2
            va += [printed[side]] * 2
        return np.array(ta), np.array(va)

    ta, va = anchor(tt, vv)
    allc = {o: np.polyfit(ta, va, o) for o in orders}
    c = allc[order]
    r = va - np.polyval(c, ta)
    print(
        f"    final fit over {len(ta)} points: resid rms {r.std():.4g}, max {np.abs(r).max():.4g}"
        f"  ({r.std() / rng * 1000:.3f} per-mille of range)"
    )
    print(f"    low frame -> {np.polyval(c, 0.0):.6g}   high frame -> {np.polyval(c, 1.0):.6g}")
    edgec = [np.polyfit(*anchor(*d), order) for d in (t_lo_edge, t_hi_edge)]
    return dict(
        c=c,
        order=order,
        allc=allc,
        edgec=edgec,
        resid=float(r.std()),
        errs=errs,
        loos=loos,
        pred=pred,
        tt=tt,
        vv=vv,
        n=int(tt.size),
    )


# --------------------------------------------------------------------------- marker seeds


def _dist(edge, X, Y, vertical=True):
    a, b = edge
    return (
        (X - (a + b * Y)) / np.hypot(1.0, b) if vertical else (Y - (a + b * X)) / np.hypot(1.0, b)
    )


def seeds(
    ink, fr, nmark, open_sz=7, rail=-14.0, min_size=0, merge=12.0, rail_min_size=0, rail_zone=12.0
):
    """Marker seeds from a morphological opening -- no density, no threshold sweep.

    The joining polyline is ~3.5 px wide; where two 3.5 px strokes cross at right angles
    the ink is ~5 px wide in every direction, so a 7x7 opening erases the line and every
    tick and leaves one core per marker.  Cores closer than 12 px are one marker split by
    the scan (A10's 70 % marker does this).

    `open_sz` and `rail` exist because the Appendix C plots are drawn with a lighter pen.
    On C24 a 7x7 opening erases the ninth marker outright -- it is not filtered out, it is
    never found -- while 6x6 keeps it and also keeps a rail of fragments along the frame
    lines. `rail` is the distance from a frame inside which a component is discarded: the
    default -14 keeps anything up to 14 px OUTSIDE the frame, which A6/A8/A10 need because
    their end markers sit on it, and a positive value discards anything within that many
    pixels INSIDE it, which C24 needs because its markers all stand at least 33 px clear.

    `min_size` drops fragments by area. A3 needs it: at the 6x6 opening its two steepest
    markers finally separate from the curve, but four fragments survive as well -- three on
    the left frame rail and one on the top frame -- and they are 36 to 48 px against 72 to
    151 px for every real marker. Size separates them where position cannot, because A3's
    first and last markers sit on those same frames.

    `rail` may also be a 4-tuple `(left, right, top, bottom)` when a page needs both at once.
    A11 does: its end markers sit ON the left and right frames, so those two must stay
    negative, while a band of tick-label fragments runs along the inside of the bottom
    frame and has to be excluded -- its real markers all stand 450 px clear of it.
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
        # A size floor that applies ONLY within `rail_zone` of a frame line. C27 needs
        # exactly this: seven fragments of its left frame survive the opening at 36 px,
        # its marker at XLDSA = 0 sits on that same frame at 209 px, and its marker at
        # XLDSA = 90 is also 36 px but stands in open plot. A global floor would drop the
        # real one; a rail exclusion would drop the other. Size near a frame separates them.
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


# --------------------------------------------------------------------------- the glyphs

SS = 5  # supersampling per axis, for a sub-pixel coverage model
HALF = 45.0  # half-segment length drawn from a marker along the joining line


def _blocked(X, Y, fr, allt, framew=3.6, tickpad=4.5, tickext=5.0):
    """Frame lines and tick marks, masked SYMMETRICALLY about each frame line.

    Symmetry matters: the ticks all protrude inward, so masking only the inward side would
    remove ink from one half of a marker that sits on a frame line and bias its abscissa.
    """
    m = np.zeros(X.shape, bool)
    for nm in ("left", "right", "top", "bottom"):
        vert = nm in ("left", "right")
        d = _dist(fr[nm], X, Y, vert)
        t = Y if vert else X
        m |= np.abs(d) < framew
        for c, k, _ in allt[nm]:
            if np.abs(t - c).min() > tickpad + 3:
                continue
            m |= (np.abs(t - c) < tickpad) & (np.abs(d) < k + tickext)
    return m


def windows(ink, fr, allt, C0, rad, **maskkw):
    """One scoring window per marker: the used pixels flattened, with a 5x5 sub-grid each."""
    W = []
    off = (np.arange(SS) + 0.5) / SS - 0.5
    ox, oy = np.meshgrid(off, off, indexing="ij")
    ox, oy = ox.ravel(), oy.ravel()
    half = int(rad) + 5
    for x0, y0 in C0:
        xs = np.arange(int(x0) - half, int(x0) + half + 1)
        ys = np.arange(int(y0) - half, int(y0) + half + 1)
        X, Y = np.meshgrid(xs, ys)
        use = ((X - x0) ** 2 + (Y - y0) ** 2) < rad**2
        use &= ~_blocked(X, Y, fr, allt, **maskkw)
        data = ink[ys[0] : ys[-1] + 1, xs[0] : xs[-1] + 1].astype(float)
        px, py = X[use].astype(float), Y[use].astype(float)
        W.append(
            dict(
                X=X,
                Y=Y,
                use=use,
                data=data,
                px=px,
                py=py,
                d=data[use],
                OX=px[:, None] + ox[None, :],
                OY=py[:, None] + oy[None, :],
                keep=np.ones(px.shape, bool),
            )
        )
    return W


def _bar(w, c, ang, half, wid, halfline=False):
    u0, u1 = np.cos(ang), np.sin(ang)
    dx = w["OX"] - c[0]
    dy = w["OY"] - c[1]
    t = dx * u0 + dy * u1
    p = -dx * u1 + dy * u0
    m = np.abs(p) <= wid / 2
    m &= (t >= 0) & (t <= half) if halfline else (np.abs(t) <= half)
    return m


def _dirs(C, i):
    """Unit vectors from marker i towards its neighbours -- the joining segments."""
    out = []
    for j in (i - 1, i + 1):
        if 0 <= j < len(C):
            d = C[j] - C[i]
            out.append(d / np.hypot(*d))
        else:
            out.append(None)
    return out


def coverage(w, c, tpl, pin, pout, withline=True):
    th1, l1, th2, l2, wg, wl = tpl
    m = _bar(w, c, th1, l1, wg) | _bar(w, c, th2, l2, wg)
    if withline:
        for p in (pin, pout):
            if p is not None:
                m |= _bar(w, c, np.arctan2(p[1], p[0]), HALF, wl, True)
    return m.mean(axis=1)


def _cost_one(w, c, tpl, pin, pout, withline=True):
    k = w["keep"]
    cov = coverage(w, c, tpl, pin, pout, withline)
    return float((((cov - w["d"]) ** 2) * k).sum())


TPL0 = (np.radians(-44.05), 14.88, np.radians(46.09), 15.22, 3.46, 3.52)  # A7's, as the seed


def fit_glyphs(W, C0, rounds=5, withline=True, verbose=True, tpl0=None, fit_template=True):
    """Coordinate descent: each centre against its own window, then the shared template."""
    C = C0.astype(float).copy()
    tpl = np.array(TPL0 if tpl0 is None else tpl0, float)
    cost = np.nan
    for it in range(rounds):
        for i in range(len(C)):
            pin, pout = _dirs(C, i)
            if not withline:
                pin = pout = None
            # Optimize the OFFSET from the current centre, with an explicit 2 px simplex.
            # Nelder-Mead's default simplex is 5 % of each coordinate, which at x ~ 2269
            # is 113 px -- far outside the scoring window, where the cost is flat, so the
            # fit walks away and stays away.  (It did exactly that on A6's right marker.)
            c0 = C[i].copy()
            r = minimize(
                lambda dc, i=i, pin=pin, pout=pout, c0=c0, tpl=tpl: _cost_one(
                    W[i], c0 + dc, tpl, pin, pout, withline
                ),
                np.zeros(2),
                method="Nelder-Mead",
                options=dict(
                    xatol=2e-3,
                    fatol=1e-6,
                    maxiter=400,
                    initial_simplex=np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]]),
                ),
            )
            if np.hypot(*r.x) > 10.0:
                raise AssertionError(f"marker {i} moved {np.hypot(*r.x):.1f} px -- bad seed?")
            C[i] = c0 + r.x
            cost = r.fun
        if fit_template:
            D = [_dirs(C, i) if withline else (None, None) for i in range(len(C))]

            def tcost(p, D=D):
                return sum(
                    _cost_one(W[i], C[i], p, D[i][0], D[i][1], withline) for i in range(len(C))
                )

            r = minimize(
                tcost, tpl, method="Powell", options=dict(xtol=1e-3, ftol=1e-5, maxfev=1200)
            )
            tpl, cost = r.x, r.fun
        if verbose:
            print(f"    round {it + 1}: cost {cost:.3f}")
    if verbose:
        print(
            f"    template: {np.degrees(tpl[0]):+.2f} deg / {tpl[1]:.2f} px and "
            f"{np.degrees(tpl[2]):+.2f} deg / {tpl[3]:.2f} px, "
            f"glyph width {tpl[4]:.2f}, line width {tpl[5]:.2f}"
        )
    return C, tpl


def glyph_only_windows(ink, fr, allt, C, tpl, rad, **maskkw):
    """Same windows, but with the joining polyline's own pixels excluded from scoring."""
    W = windows(ink, fr, allt, C, rad, **maskkw)
    for i, w in enumerate(W):
        pin, pout = _dirs(C, i)
        drop = np.zeros(w["px"].shape, bool)
        for p in (pin, pout):
            if p is None:
                continue
            wl = tpl[5] + 3.0
            dx = w["px"] - C[i][0]
            dy = w["py"] - C[i][1]
            t = dx * p[0] + dy * p[1]
            q = -dx * p[1] + dy * p[0]
            drop |= (np.abs(q) <= wl / 2) & (t >= 0)
        w["keep"] = ~drop
    return W


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
    for i in range(len(C) - 1):
        a, b = C[i], C[i + 1]
        d = b - a
        L = float(np.hypot(*d))
        d = d / L
        t = (P - a) @ d
        perp = (P - a) @ np.array([-d[1], d[0]])
        sel = (np.abs(perp) < 6) & (t > 24) & (t < L - 24)
        if sel.sum() > 40:
            sag.append((L, float(np.abs(perp[sel]).mean()), float(perp[sel].std())))
        tt = t[np.abs(perp) < 12]
        nb = max(int(L // 10), 4)
        h, e = np.histogram(tt, bins=nb, range=(0.0, L))
        base = max(np.median(h[2:-2]) if nb > 6 else np.median(h), 1)
        extras.append(
            [
                round((e[k] + e[k + 1]) / 2, 1)
                for k in range(nb)
                if h[k] > 1.9 * base and 14 < (e[k] + e[k + 1]) / 2 < L - 14
            ]
        )
    print(f"  interior ink excesses per segment: {extras}  (all empty => no missed marker)")
    if sag:
        print(
            "  joining segments: "
            + ", ".join(f"{L:.0f} px rms {s:.2f}" for L, _, s in sag[: min(len(sag), 12)])
        )
    return extras, sag


# --------------------------------------------------------------------------- overlay


def overlay(ink, fr, C, cx, cy, cfg, tag):
    im = Image.fromarray((~ink * 255).astype(np.uint8)).convert("RGB")
    d = ImageDraw.Draw(im)
    TL = corner(fr["top"], fr["left"])
    TR = corner(fr["top"], fr["right"])
    BL = corner(fr["bottom"], fr["left"])
    BR = corner(fr["bottom"], fr["right"])
    Hi = homography([(0, 1), (1, 1), (1, 0), (0, 0)], [TL, TR, BR, BL])
    tt = np.linspace(0.0, 1.0, 20001)
    for spec, poly, horiz in ((cfg["x"], cx, False), (cfg["y"], cy, True)):
        n = int(round((spec["hi"] - spec["lo"]) / major_step(spec)))
        for k in range(n + 1):
            val = spec["lo"] + major_step(spec) * k
            s = float(tt[np.argmin(np.abs(np.polyval(poly, tt) - val))])
            p0, p1 = ((0.0, s), (1.0, s)) if horiz else ((s, 0.0), (s, 1.0))
            q0, q1 = _apply(Hi, *p0), _apply(Hi, *p1)
            d.line(
                [(q0[0].item(), q0[1].item()), (q1[0].item(), q1[1].item())],
                fill=(210, 230, 255),
                width=1,
            )
    for x, y in C:
        d.ellipse([x - 26, y - 26, x + 26, y + 26], outline=(220, 0, 0), width=3)
        d.line([(x - 40, y), (x + 40, y)], fill=(0, 170, 0), width=1)
        d.line([(x, y - 40), (x, y + 40)], fill=(0, 170, 0), width=1)
    x0 = int(min(TL[0], BL[0])) - 200
    x1 = int(max(TR[0], BR[0])) + 60
    y0 = int(min(TL[1], TR[1])) - 60
    y1 = int(max(BL[1], BR[1])) + 90
    im.crop((x0, y0, x1, y1)).resize(((x1 - x0) // 2, (y1 - y0) // 2), Image.LANCZOS).save(
        OUT / f"{tag}_verify.png"
    )


# --------------------------------------------------------------------------- driver


def run(key):
    cfg = FIGS[key]
    print(f"\n================ {cfg['fig']}  (pdf p.{cfg['page']})")
    ink = native_bitmap(cfg["page"])
    print("frame lines:")
    spec = cfg.get("frame_spec") or auto_frame_spec(cfg["page"])
    fr = frame(ink, spec)
    TL = corner(fr["top"], fr["left"])
    TR = corner(fr["top"], fr["right"])
    BL = corner(fr["bottom"], fr["left"])
    BR = corner(fr["bottom"], fr["right"])
    print(
        f"  corners TL({TL[0]:.3f},{TL[1]:.3f}) TR({TR[0]:.3f},{TR[1]:.3f}) "
        f"BR({BR[0]:.3f},{BR[1]:.3f}) BL({BL[0]:.3f},{BL[1]:.3f})"
    )
    H = homography([TL, TR, BR, BL], [(0, 1), (1, 1), (1, 0), (0, 0)])
    print("ticks:")
    maj, allt, _ = ticks(ink, fr, cfg)

    # tick coordinates in the unit square
    def uv_of(nm, pos):
        if nm in ("left", "right"):
            a, b = fr[nm]
            return _apply(H, [a + b * p for p in pos], pos)
        a, b = fr[nm]
        return _apply(H, pos, [a + b * p for p in pos])

    xs = cfg["x"]
    ys = cfg["y"]
    xv = major_values(xs)
    yv = major_values(ys)

    def paired(name, vals, descending):
        """Positions and values, dropping any major hidden by a marker.

        `_pick_majors` returns slot order -- low edge to high edge -- so slot k carries
        `vals[k]` on the x edges and `vals[-1-k]` on the y edges, where pixel order runs
        opposite to value order. Pairing by INDEX rather than by sorting keeps a NaN
        attached to the value it belongs to, so dropping it drops both.
        """
        pos = np.asarray(maj[name], float)
        v = np.asarray(vals, float)[::-1] if descending else np.asarray(vals, float)
        ok = ~np.isnan(pos)
        pos, v = pos[ok], v[ok]
        o = np.argsort(pos)
        return pos[o], v[o]

    pb, xb = paired("bottom", xv, False)
    pt, xt = paired("top", xv, False)
    pl, yl = paired("left", yv, True)
    pr, yr = paired("right", yv, True)
    ub, _ = uv_of("bottom", list(pb))
    ut, _ = uv_of("top", list(pt))
    _, vl = uv_of("left", list(pl))
    _, vr = uv_of("right", list(pr))
    print("calibration:")
    CX = calibrate_axis((ub, xb), (ut, xt), xs, "x")
    CY = calibrate_axis((vl, yl), (vr, yr), ys, "y")
    cx, cy = CX["c"], CY["c"]

    print("markers:")
    C0 = seeds(
        ink,
        fr,
        cfg["nmark"],
        cfg.get("open_sz", 7),
        cfg.get("rail", -14.0),
        cfg.get("min_size", 0),
        cfg.get("merge", 12.0),
        cfg.get("rail_min_size", 0),
        cfg.get("rail_zone", 12.0),
    )
    W = windows(ink, fr, allt, C0, cfg["rad"])
    C, tpl = fit_glyphs(W, C0, rounds=4)
    # re-centre the windows on the fitted centres and refit, so the disc is self-consistent
    W = windows(ink, fr, allt, C, cfg["rad"])
    C, tpl = fit_glyphs(W, C, rounds=3, tpl0=tpl)

    print("checks:")
    extras, sag = segment_check(ink, fr, C)
    u, v = _apply(H, C[:, 0], C[:, 1])
    Xd, Yd = np.polyval(cx, u), np.polyval(cy, v)

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
        Wv = windows(ink, fr, allt, C, rad, **kw)
        Cv, _ = fit_glyphs(
            Wv, C, rounds=2, withline=wl, verbose=False, tpl0=tpl, fit_template=False
        )
        var[nm] = Cv - C
        print(f"    {nm:14s} {np.abs(Cv - C).max():.2f}")
    Wg = glyph_only_windows(ink, fr, allt, C, tpl, cfg["rad"])
    Cg, _ = fit_glyphs(Wg, C, rounds=3, withline=False, verbose=False, tpl0=tpl, fit_template=False)
    var["line masked out"] = Cg - C
    print(f"    {'line masked out':14s} {np.abs(Cg - C).max():.2f}")
    dpx = np.max(np.abs(np.array(list(var.values()))), axis=0)  # (nmark, 2) worst case

    # ---- uncertainty budget, per point, 1 sigma, in data units
    print("  uncertainty budget (1 sigma, data units):")
    sig = {}
    for ax, CA, t, poly in (("x", CX, u, cx), ("y", CY, v, cy)):
        base = np.polyval(poly, t)
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
        du = _pixel_to_unit(H, C, 0 if ax == "x" else 1)
        d = np.abs(np.polyval(np.polyder(poly), t)) * du * dpx[:, 0 if ax == "x" else 1]
        s = np.sqrt(a**2 + b**2 + c**2 + d**2)
        sig[ax] = s
        print(
            f"    {ax}: tick-fit {a[0]:.4g} | order spread {b.mean():.4g} | "
            f"edge split {c.mean():.4g} | glyph {d.mean():.4g}  ->  1 sigma {s.mean():.4g} "
            f"(max {s.max():.4g})"
        )

    # ---- independent bow-free check: local linear interpolation between the two majors
    #      that bracket each point, done on each ticked edge and averaged
    print("  local two-tick check (bow-free), chosen-model minus local:")
    twotick = {}
    for ax, t, spec, poly, edges in (
        ("x", u, xs, cx, (ub, ut)),
        ("y", v, ys, cy, (vl, vr)),
    ):
        vals = xv if ax == "x" else yv
        per = []
        for e in edges:
            nodes_t = np.concatenate([[0.0], e, [1.0]])
            nodes_v = np.concatenate(
                [
                    [spec["lo"] if "lo" in spec["printed"] else np.polyval(poly, 0.0)],
                    vals,
                    [spec["hi"] if "hi" in spec["printed"] else np.polyval(poly, 1.0)],
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

    overlay(ink, fr, C, cx, cy, cfg, key)
    return dict(
        cfg=cfg,
        key=key,
        fr=fr,
        frdiag=dict(FRAME_DIAG),
        extras=extras,
        sag=sag,
        varmax={k2: float(np.abs(v2).max()) for k2, v2 in var.items()},
        corners=(TL, TR, BR, BL),
        maj=maj,
        CX=CX,
        CY=CY,
        cx=cx,
        cy=cy,
        C=C,
        C0=C0,
        tpl=tpl,
        u=u,
        v=v,
        X=Xd,
        Y=Yd,
        dpx=dpx,
        var=var,
        sig=sig,
        loc=twotick,
        H=H,
    )


# --------------------------------------------------------------------------- the CSVs

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
    ],  # the 2026-09-10 extraction, kept so the move can be reported in sigma
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


def _poly_str(c, var="t"):
    n = len(c) - 1
    out = []
    for i, a in enumerate(c):
        p = n - i
        t = "" if p == 0 else f" {var}" if p == 1 else f" {var}^{p}"
        out.append(f"{a:+.8g}{t}")
    return " ".join(out)


def _cal_block(CA, spec, axis):
    order = CA["order"]
    var = "u" if axis == "x" else "v"
    side = {"lo": "left" if axis == "x" else "bottom", "hi": "right" if axis == "x" else "top"}
    anch = spec.get("anchor", spec["printed"])
    words = ", ".join(
        f"{side[s]} = {spec[s]:g}" + ("" if s in spec["printed"] else " (DERIVED, see below)")
        for s in ("lo", "hi")
        if s in anch
    )
    lines = [
        f"#   {axis} axis: {['', 'LINEAR', 'QUADRATIC', 'CUBIC', 'QUARTIC'][order]} in {var}, "
        f"where {var} = 0 at the {side['lo']} frame and 1 at the {side['hi']} frame.",
        f"#      {axis} = {_poly_str(CA['c'], var)}",
        f"#      Fitted on {CA['n']} interior major ticks (both ticked edges) plus the frame",
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
    A("# GEOMETRY. The four frame lines are fitted sub-pixel by tracing the single ink run per")
    A("# row/column and least-squares fitting with outlier rejection:")
    for nm in ("left", "right", "top", "bottom"):
        a, b, rms, n = r["frdiag"][nm]
        v = "x" if nm in ("left", "right") else "y"
        t = "y" if nm in ("left", "right") else "x"
        A(
            f"#   {nm:6s} {v} = {a:.3f} {b:+.6f} {t}   ({np.degrees(np.arctan(b)):+.3f} deg, "
            f"{n} samples, rms {rms:.2f} px)"
        )
    TL, TR, BR, BL = r["corners"]
    A("# The four edge slopes are not a common rotation, so the frame is a genuine")
    A(f"# quadrilateral; corners TL({TL[0]:.3f},{TL[1]:.3f}) TR({TR[0]:.3f},{TR[1]:.3f})")
    A(f"# BR({BR[0]:.3f},{BR[1]:.3f}) BL({BL[0]:.3f},{BL[1]:.3f}) are mapped onto the unit")
    A("# square (u right, v up) by a homography, and every measured point is carried through")
    A("# it. NOT `digitize.py axes`, and NOT a pixel warp -- nothing is resampled.")
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
    L.extend(_cal_block(r["CX"], xs, "x"))
    L.extend(_cal_block(r["CY"], ys, "y"))
    A("#   Independent bow-free cross-check: linear interpolation between the two majors that")
    A("#   bracket each point, done on each ticked edge and averaged, agrees with the chosen")
    A(
        f"#   model to {np.abs(np.polyval(r['cx'], r['u']) - r['loc']['x']).max():.3g} in x and "
        f"{np.abs(np.polyval(r['cy'], r['v']) - r['loc']['y']).max():.3g} in y, worst point."
    )
    A("#")
    A(f"# MARKERS. {cfg['nmark']} 'x' glyphs. Seeds come from a 7x7 morphological OPENING of the")
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
        base = np.polyval(CA["c"], t)
        b = np.max([np.abs(np.polyval(CA["allc"][o], t) - base) for o in CA["allc"]], axis=0)
        c = 0.5 * np.abs(np.polyval(CA["edgec"][0], t) - np.polyval(CA["edgec"][1], t))
        d = (
            np.abs(np.polyval(np.polyder(CA["c"]), t))
            * _pixel_to_unit(r["H"], r["C"], 0 if ax == "x" else 1)
            * r["dpx"][:, 0 if ax == "x" else 1]
        )
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


def _pixel_to_unit(H, C, axis):
    """d(u or v) / d(pixel) at each marker, along that axis -- the local map scale."""
    u0, v0 = _apply(H, C[:, 0], C[:, 1])
    if axis == 0:
        u1, _ = _apply(H, C[:, 0] + 1.0, C[:, 1])
        return np.abs(u1 - u0)
    _, v1 = _apply(H, C[:, 0], C[:, 1] + 1.0)
    return np.abs(v1 - v0)


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    # Every figure this tool owns, not a subset. It defaulted to ["a8", "a10", "a6"] --
    # the three the old filename `digitize_a6810.py` named -- while FIGS grew to thirteen,
    # and `tools/reproduce_all.sh` invokes this with no arguments. So ten committed files
    # (f3, f4, f5, f_hs and F_HM1-F_HM6) carried a header saying "reruns and reproduces
    # this file" that nothing had ever checked. Found 2026-09-14 by the rename: a full gate
    # run rewrote the provenance line in three files and left those ten untouched.
    keys = [a for a in argv if a in FIGS] or list(FIGS)
    for k in keys:
        if "--csv-only" in argv:
            with (OUT / f"{k}_result.pkl").open("rb") as fh:
                r = pickle.load(fh)
        else:
            r = run(k)
            with (OUT / f"{k}_result.pkl").open("wb") as fh:
                pickle.dump(r, fh)
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
