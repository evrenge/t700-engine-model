"""What does the calibration model choice cost, figure by figure? Writes nothing.

Each single-curve figure is recalibrated under each variant, everything else held fixed, and
the difference printed in the figure's own units and in multiples of its published read error.
The baseline is order 1 with the soft anchor; the variants (see `figure.axis`):

    frame     order 2: four fitted arcs and the curvilinear blend
    anchor    once (`blend_axis`), or hard, exact at both frames (`fit_anchored`)
    both      order 2 with the hard anchor

The glyph centres are held at the values the extraction finds, so only the calibration is
measured: the frame model also moves the marker windows, and mixing the two effects would make
neither attributable.

    python3 tools/compare_models.py            # every single-curve figure
    python3 tools/compare_models.py a6 a8      # just these
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from digitize_single_curve import FIGS, run  # noqa: E402
from figure.axis import blend_axis, blend_value, fit_anchored, soft_anchor  # noqa: E402
from figure.frame import edge_spans, fit_edges, unit_square  # noqa: E402
from figure.ticks import coarse_lattice, fine_lattice_all, majors_all_edges  # noqa: E402


def plot_coords(ink, spec, C, order):
    """Marker centres in the unit square, under a frame model of the given order."""
    fc, _ = fit_edges(ink, spec, order=order, verbose=False)
    to_uv, _, _ = unit_square(fc, order)
    u, v = to_uv(C[:, 0], C[:, 1])
    return fc, u, v


def edge_positions(ink, fc, cfg, order):
    """Where each printed tick sits along its edge, in the unit square.

    Uses the minor lattice wherever one resolves, whereas `pick_lattice` sends some axes back
    to their majors, so this is a slightly different calibration from the one that ships.
    """
    maj, _allt, cor = majors_all_edges(ink, fc, cfg, verbose=False, order=order)
    _, _, at = unit_square(fc, order)
    return fine_lattice_all(
        ink, fc, cfg, at, edge_spans(cor), coarse_lattice(maj, cfg, at), verbose=False
    )


def axis_polys(P, spec, order, anchor):
    """The two per-edge maps for one axis, under the chosen anchoring rule."""
    lo_nm, hi_nm = ("bottom", "top") if spec["_axis"] == "x" else ("left", "right")
    printed = {"lo": spec["lo"], "hi": spec["hi"]}
    at = {"lo": 0.0, "hi": 1.0}
    sides = spec.get("anchor", spec["printed"])
    if anchor == "once":
        return list(blend_axis(*P[lo_nm], *P[hi_nm], order, spec["lo"], spec["hi"]))
    out = []
    for nm in (lo_nm, hi_nm):
        t, val = P[nm]
        if anchor == "hard":
            out.append(fit_anchored(t, val, order, spec["lo"], spec["hi"]))
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


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
