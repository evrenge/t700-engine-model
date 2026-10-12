"""What does the calibration model choice cost, figure by figure?

Three decisions were made differently in different corners of this repository, each for a
reason its own author wrote down, and none of them ever reached the other corner. This
tool answers all three the only way they can honestly be answered: run both, hold
everything else fixed, and print the difference in the figure's own units and in multiples
of the read error that figure already publishes.

    FRAME     order 1, four straight edges mapped by a homography  -- what the thirteen
              single-curve figures and A1 use
              order 2, four fitted arcs with a curvilinear blend   -- what A2, A7 and A9
              use. Recipe trap 7 used to call it right; measured, eleven of the
              thirteen single-curve figures prefer straight

    ANCHOR    soft, the printed frame value entered twice among the ticks -- the
              single-curve path
              once, `blend_axis` -- the native path
              hard, `fit_anchored`, exact at both frames by construction -- A7, and
              `digitize_fig678` independently: "the t(1-t) factor is the whole point,
              it vanishes at both frames"

**The glyph centres are held fixed** at the values the committed extraction found, so what
is measured here is the calibration alone. That is deliberate: the frame model also moves
the marker windows slightly, and mixing the two effects would make neither attributable.
The glyph question -- trap 8's simplex -- is measured separately, by flipping the guard on
in the tool and rerunning it.

Writes nothing. Run it, read the table, then decide per figure.

    python3 tools/compare_models.py            # every single-curve figure
    python3 tools/compare_models.py a6 a8      # just these
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from digitize_single_curve import FIGS, run  # noqa: E402
from figure.axis import blend_value, choose, fit_anchored, loo, soft_anchor  # noqa: E402
from figure.frame import (  # noqa: E402  # noqa: E402  # noqa: E402
    EDGE,
    apply_h,
    as_lines,
    auto_frame_spec,
    corners,
    edge_at,
    edge_at_straight,
    edge_uv,
    fit_edges,
    homography,
    make_sv,
)
from figure.raster import native_bitmap  # noqa: E402
from figure.ticks import (  # noqa: E402  # noqa: E402
    fine_lattice,
    major_step,
    major_values,
    majors_all_edges,
    paired,
)

UNIT = (0, 1), (1, 1), (1, 0), (0, 0)


def plot_coords(ink, spec, C, order):
    """Marker centres in the unit square, under a frame model of the given order."""
    fc, _ = fit_edges(ink, spec, order=order, verbose=False)
    if order == 1:
        H = homography(corners(fc, order=1), UNIT)
        u, v = apply_h(H, C[:, 0], C[:, 1])
        return fc, u, v
    to_sv, _ = make_sv(fc)
    s, v = to_sv(C[:, 0], C[:, 1])
    return fc, s, v


def edge_positions(ink, fc, cfg, order, fine=True):
    """Where each printed tick sits along its edge, in the unit square.

    `fine=True` uses the minor lattice on every edge where one resolves. That is not quite
    what ships: `digitize_single_curve` sends five axes back to their majors (`pick_lattice`),
    so these numbers describe a slightly different calibration. On 2026-09-19 it changed no
    frame-model choice.
    """
    maj, _allt, _cor = majors_all_edges(ink, fc, cfg, verbose=False, order=order)
    if order == 1:
        uv_of = edge_uv(as_lines(fc), homography(corners(fc, order=1), UNIT))
        at = edge_at_straight(uv_of)
    else:
        to_sv, _ = make_sv(fc)
        at = edge_at(fc, to_sv)

    TL, TR, BR, BL = corners(fc, order=order)
    span = {
        "left": (TL[1], BL[1]),
        "right": (TR[1], BR[1]),
        "top": (TL[0], TR[0]),
        "bottom": (BL[0], BR[0]),
    }
    out = {}
    for nm, ax, desc in (
        ("bottom", "x", False),
        ("top", "x", False),
        ("left", "y", True),
        ("right", "y", True),
    ):
        vals = major_values(cfg[ax])
        pos, val = paired(maj, nm, vals, desc)
        u = np.asarray(at(nm, pos), float)
        out[nm] = (u, val)
        if not fine:
            continue
        along, sign = EDGE[nm]
        a, b = span[nm]
        lo, hi = int(min(a, b)) + 3, int(max(a, b)) - 2
        major_px = abs(b - a) / ((cfg[ax]["hi"] - cfg[ax]["lo"]) / major_step(cfg[ax]))
        coarse = np.polyfit(u, val, max(1, min(3, u.size - 2)))
        got = fine_lattice(
            ink,
            fc[nm],
            along,
            sign,
            lo,
            hi,
            lambda raw, nm=nm: np.asarray(at(nm, raw), float),
            coarse,
            major_px,
            major_step(cfg[ax]),
        )
        if got is not None:
            fu, fv, _n, _d = got
            out[nm] = (fu, fv)
    return out


def axis_polys(P, spec, order, anchor):
    """The two per-edge maps for one axis, under the chosen anchoring rule."""
    lo_nm, hi_nm = ("bottom", "top") if spec["_axis"] == "x" else ("left", "right")
    printed = {"lo": spec["lo"], "hi": spec["hi"]}
    at = {"lo": 0.0, "hi": 1.0}
    sides = spec.get("anchor", spec["printed"])
    out = []
    for nm in (lo_nm, hi_nm):
        t, val = P[nm]
        if anchor == "hard":
            out.append(fit_anchored(t, val, order, spec["lo"], spec["hi"]))
        elif anchor == "once":
            out.append(
                np.polyfit(
                    np.concatenate([t, [0.0, 1.0]]),
                    np.concatenate([val, [spec["lo"], spec["hi"]]]),
                    order,
                )
            )
        else:
            out.append(np.polyfit(*soft_anchor(t, val, sides, at, printed), order))
    return out


def values(ink, spec_frame, cfg, C, frame_order, anchor, orders):
    fc, u, v = plot_coords(ink, spec_frame, C, frame_order)
    P = edge_positions(ink, fc, cfg, frame_order)
    xs = dict(cfg["x"], _axis="x")
    ys = dict(cfg["y"], _axis="y")
    px = axis_polys(P, xs, orders[0], anchor)
    py = axis_polys(P, ys, orders[1], anchor)
    return blend_value(px, u, v), blend_value(py, v, u)


def sigma_of(csv: Path):
    rows = [ln for ln in csv.read_text().splitlines() if ln and not ln.startswith("#")]
    head = [h.strip() for h in rows[0].split(",")]
    cols = np.array([[float(x) for x in r.split(",")] for r in rows[1:]])
    get = lambda n: cols[:, head.index(n)] if n in head else None  # noqa: E731
    return get("sigma_x"), get("sigma_y")


def compare(key):
    r = run(key)
    cfg, C, ink, spec = r["cfg"], r["C"], r["ink"], r["fspec"]
    orders = (r["CX"]["order"], r["CY"]["order"])
    sx, sy = sigma_of(cfg["csv"])
    base = values(ink, spec, cfg, C, 1, "soft", orders)
    rows = []
    for label, fo, anc in (
        ("arcs (frame order 2)", 2, "soft"),
        ("hard anchor", 1, "hard"),
        ("one anchor copy", 1, "once"),
        ("arcs + hard anchor", 2, "hard"),
    ):
        try:
            X, Y = values(ink, spec, cfg, C, fo, anc, orders)
        except Exception as exc:  # a variant a figure cannot support says so
            rows.append((label, None, None, None, None, str(exc)[:40]))
            continue
        dx, dy = np.abs(X - base[0]), np.abs(Y - base[1])
        rows.append(
            (
                label,
                dx.max(),
                dy.max(),
                (dx / sx).max() if sx is not None else np.nan,
                (dy / sy).max() if sy is not None else np.nan,
                "",
            )
        )
    return cfg, orders, rows


def main(argv):
    keys = [a for a in argv if a in FIGS] or list(FIGS)
    if "--verdict" in argv:
        return verdict(keys)
    print("\n" + "=" * 92)
    print("CALIBRATION MODEL COMPARISON -- glyph centres held fixed, calibration varied")
    print("=" * 92)
    table = {}
    for key in keys:
        cfg, orders, rows = compare(key)
        table[key] = (cfg, orders, rows)
    print(
        f"\n{'figure':>7} {'variant':>22} {'max dx':>11} {'max dy':>11} "
        f"{'dx/sigma':>9} {'dy/sigma':>9}"
    )
    for _key, (cfg, orders, rows) in table.items():
        print(f"{cfg['fig']:>7}  (x order {orders[0]}, y order {orders[1]})")
        for label, dx, dy, rx, ry, err in rows:
            if err:
                print(f"{'':>7} {label:>22}   {err}")
                continue
            flag = "  <-- above the read error" if max(rx, ry) > 1.0 else ""
            print(f"{'':>7} {label:>22} {dx:11.3e} {dy:11.3e} {rx:9.2f} {ry:9.2f}{flag}")
    return 0


def held_out_under(ink, spec_frame, cfg, frame_order):
    """The figure's own verdict on a frame model: how well it predicts its printed frames.

    This is the test that makes the choice a measurement rather than an argument. Each
    figure prints a value at one or both frame edges. Fit the axis map on the INTERIOR
    majors alone -- the fit never sees the frame value -- and predict it. A frame model
    that describes the page better puts the interior ticks in more nearly the right
    places, so it predicts the number it was not given more accurately.

    It is the same test `figure.axis.calibrate` already uses to choose the polynomial
    order, asked one level up.
    """
    fc, _ = fit_edges(ink, spec_frame, order=frame_order, verbose=False)
    P = edge_positions(ink, fc, cfg, frame_order)
    edges = ("left", "right", "top", "bottom")
    sag = {nm: _sagitta(ink, spec_frame, nm, frame_order) for nm in edges}
    out = {}
    for ax, (lo_nm, hi_nm) in (("x", ("bottom", "top")), ("y", ("left", "right"))):
        spec = cfg[ax]
        t = np.concatenate([P[lo_nm][0], P[hi_nm][0]])
        val = np.concatenate([P[lo_nm][1], P[hi_nm][1]])
        rng = spec["hi"] - spec["lo"]
        errs, loos = {}, {}
        for order in (1, 2, 3, 4):
            worst = 0.0
            for side in spec["printed"]:
                at = 0.0 if side == "lo" else 1.0
                pred = float(np.polyval(np.polyfit(t, val, order), at))
                worst = max(worst, abs(pred - spec[side]))
            errs[order] = worst
            loos[order] = loo(t, val, order)
        # Each frame model is judged at ITS OWN best order, chosen by the same rule
        # `figure.axis.calibrate` uses -- otherwise the comparison is rigged by whichever
        # model's order is imposed on the other.
        order, _ok = choose(errs, loos)
        out[ax] = (order, errs[order], rng)
    return out, sag


def _sagitta(ink, spec_frame, name, order):
    _fc, diag = fit_edges(ink, spec_frame, order=max(order, 2), verbose=False)
    return diag[name]["sagitta"]


def verdict(keys):
    print("\n" + "=" * 92)
    print("WHICH FRAME MODEL DOES THE FIGURE ITSELF PREFER?")
    print("Held-out frame-edge error: fit the interior majors alone, predict the printed")
    print("frame value. Lower is better. `sagitta` is the measured bow of each printed")
    print("straight edge, in pixels -- the physical quantity trap 7 started from.")
    print("=" * 92)
    print(
        f"\n{'figure':>7} {'axis':>4} {'ord s/a':>9} {'straight':>12} {'arcs':>12} "
        f"{'better':>9}   sagittae L/R/T/B (px)"
    )
    for key in keys:
        cfg = FIGS[key]
        ink = native_bitmap(cfg["page"])
        fspec = cfg.get("frame_spec") or auto_frame_spec(cfg["page"])
        h1, sag = held_out_under(ink, fspec, cfg, 1)
        h2, _ = held_out_under(ink, fspec, cfg, 2)
        sg = "  ".join(f"{sag[n]:+.2f}" for n in ("left", "right", "top", "bottom"))
        for ax in ("x", "y"):
            (o1, e1, rng), (o2, e2, _) = h1[ax], h2[ax]
            who = "arcs" if e2 < e1 else "straight"
            ratio = max(e1, e2) / min(e1, e2) if min(e1, e2) else float("inf")
            print(
                f"{cfg['fig']:>7} {ax:>4} {o1:>4}/{o2:<4} {e1 / rng * 1000:9.3f} pm "
                f"{e2 / rng * 1000:9.3f} pm {who:>9} ({ratio:.2f}x)"
                + (f"   {sg}" if ax == "x" else "")
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
