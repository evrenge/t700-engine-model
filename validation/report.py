"""Every comparison this project can make against the report, as tables and overlay plots.

The question it answers: **what is the complete list of things the report prints, and
where do we stand against each one?** It recomputes every headline number rather than
quoting one, writes `site/assets/*.png` for the overlays in both themes, and dumps
`report.json` so the numbers can be checked against `README.md` and `SCOPE.md`
mechanically instead of by reading.

This is the only figure pipeline. `plot_validation.py` was a second one -- seven sheets,
five of them strict subsets of these, with its own hardcoded copy of `plotstyle`'s palette
and a character-identical copy of `test_whole_curve._split_strays` -- and it was deleted on
2026-09-14. Its two sheets that framed something new are absorbed: Table B.1's trims are
drawn onto the Figures 6-7 overlays, and `internal_spread` computes what its `residuals`
sheet drew.

Model runs are imported from the test modules that own them, so the report measures the
same runs the suite does. The statistics are computed here, independently of the tests'
own assertions -- an audit that shares its arithmetic with the thing it audits is not one.

Four comparison surfaces exist in the report and all four are covered:

    Table B.1   pdf p.67   21 printed numbers, the engine state at three trims
    Table 1     pdf p.31   27 printed eigenvalues, three model variants
    Appendix B  pp.68-76   297 printed matrix elements over twelve figures
    Figures 6-10           five reproducible result figures, ~8,100 digitized points

Matplotlib lives here and nowhere near `src/t700/`, per CLAUDE.md.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
for _p in (HERE, HERE.parent / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import plotstyle  # noqa: E402
from t700 import appendix_b as ab  # noqa: E402
from t700 import constants as c  # noqa: E402
from t700 import trim  # noqa: E402
from t700.engine import Ambient, frame  # noqa: E402
from t700.linear import DOF, extract  # noqa: E402
from t700.units import shp_from_torque, wf_pps_from_pph  # noqa: E402

OUT = HERE.parent / "site" / "assets"
REF = HERE.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)

THEME = plotstyle.LIGHT
"""The theme the next figure is drawn in. Rebound by `use_theme`.

Every sheet draws through the colour names below, and `main` renders the whole set once per
theme so the page can hand a reader the variant their system asks for. The cost is that each
sheet's model runs happen twice; this is a build step and not on any hot path, and splitting
a dozen sheets into compute and draw halves for the sake of it would double the file."""

OURS = THEME.ours
BALLIN = THEME.ballin
GE = THEME.ge
GE2 = THEME.ge2
MUTED = THEME.muted
GRID = THEME.grid


def use_theme(theme: plotstyle.Theme) -> None:
    """Point the module's colour names at one theme's tokens."""
    global THEME, OURS, BALLIN, GE, GE2, MUTED, GRID
    THEME = theme
    OURS, BALLIN, GE, GE2 = theme.ours, theme.ballin, theme.ge, theme.ge2
    MUTED, GRID = theme.muted, theme.grid


def fs(w: float, h: float) -> tuple[float, float]:
    """A figure size in inches, scaled so its content reads at the page's column width."""
    return (w * plotstyle.FIG_SCALE, h * plotstyle.FIG_SCALE)


def save(fig, stem: str) -> None:
    """Write `<stem>-<theme>.png` into the site's asset directory."""
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}-{THEME.name}.png", dpi=plotstyle.DPI)
    plt.close(fig)


TRIMS = ("hover", "level 80 kt", "descent 80 kt")
WF_PPH = {1: 476.3, 2: 349.3, 3: 267.7}
DQ_REQ_DNP = {1: 0.019471, 2: 0.015869, 3: 0.012950}
"""Recovered from Appendix B's own NP diagonal -- open question #6. Gen Hel's, not ours."""

TABLE_1 = {
    # [TM-100991 pdf p.31] read cell by cell from the page image; see
    # docs/notes/body-realtime.md 5.5. Blank cells are blank in the original.
    1: {
        "5dof": {"NG": -2.66, "NP": -0.565, "P3": -51.6, "P41": -4900.0, "P45": -3060.0},
        "2dof": {"NG": -2.69, "NP": -0.565},
        "red5": {"NG": -2.81, "NP": -0.565},
    },
    2: {
        "5dof": {"NG": -2.08, "NP": -0.446, "P3": -52.2, "P41": -4640.0, "P45": -4040.0},
        "2dof": {"NG": -2.23, "NP": -0.446},
        "red5": {"NG": -2.16, "NP": -0.446},
    },
    3: {
        "5dof": {"NG": -1.75, "NP": -0.357, "P3": -52.6, "P41": -4430.0, "P45": -4530.0},
        "2dof": {"NG": -1.82, "NP": -0.357},
        "red5": {"NG": -1.83, "NP": -0.357},
    },
}

STATE_B1 = {
    # [TM-100991 pdf p.67, Table B.1]
    "hover": dict(
        wf=476.3, ng=41638.0, ps3=176.34, p41=174.28, t41=2292.0, p45=37.42, t45=1632.0, shp=911.1
    ),
    "level 80 kt": dict(
        wf=349.3, ng=39768.0, ps3=142.13, p41=140.26, t41=2102.0, p45=30.66, t45=1501.0, shp=552.6
    ),
    "descent 80 kt": dict(
        wf=267.7, ng=38072.0, ps3=114.27, p41=112.77, t41=1982.0, p45=25.54, t45=1424.0, shp=302.6
    ),
}


def _style(ax, xlabel="", ylabel="", title=""):
    """Labels and title. Grid, spines, tick colours and faces all come from rcParams."""
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontfamily=[plotstyle.SANS_CONDENSED, plotstyle.SANS])
    plotstyle.mono_ticks(ax, THEME)


def _csv(path: Path, xc: str, yc: str):
    rows = list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))
    x = np.array([float(r[xc]) for r in rows])
    y = np.array([float(r[yc]) for r in rows])
    o = np.argsort(x)
    return x[o], y[o]


def _stats(dev: np.ndarray) -> dict:
    dev = np.asarray(dev, dtype=float)
    return {
        "n": int(dev.size),
        "mean": float(dev.mean()),
        "rms": float(np.sqrt((dev**2).mean())),
        "worst": float(dev[np.argmax(np.abs(dev))]),
    }


# ============================================================ Figures 6-8, the steady sweeps

FIG678 = {
    6: dict(
        file="fig06",
        xc="wf_pph",
        yc="ng_pct",
        key="ng_pct",
        xlabel="fuel flow Wf, lbm/hr",
        ylabel="gas generator speed, %NG",
        title="Figure 6 [pdf p.40] -- NG against fuel flow",
    ),
    7: dict(
        file="fig07",
        xc="wf_pph",
        yc="shp",
        key="shp",
        xlabel="fuel flow Wf, lbm/hr",
        ylabel="shaft horsepower",
        title="Figure 7 [pdf p.41] -- shaft power against fuel flow",
    ),
    8: dict(
        file="fig08",
        xc="ng_pct",
        yc="ps3_psia",
        key="ps3",
        xlabel="gas generator speed, %NG",
        ylabel="compressor discharge Ps3, psia",
        title="Figure 8 [pdf p.42] -- Ps3 against NG",
    ),
}


def _b1_markers(no: int) -> tuple[np.ndarray, np.ndarray]:
    """Table B.1's three trims plotted on Figure `no`'s own axes."""
    xs, ys = [], []
    for p in STATE_B1.values():
        ng_pct = 100.0 * p["ng"] / c.NG_DES
        xs.append(ng_pct if no == 8 else p["wf"])
        ys.append({6: ng_pct, 7: p["shp"], 8: p["ps3"]}[no])
    o = np.argsort(xs)
    return np.array(xs)[o], np.array(ys)[o]


def steady_sweeps() -> dict:
    """Overlay each of Figures 6-8 with all three printed series, and score ours."""
    import test_steady_sweeps as tss

    out = {}
    fig, axes = plt.subplots(1, 3, figsize=fs(16.5, 4.8))

    for ax, (no, spec) in zip(axes, sorted(FIG678.items()), strict=True):
        rt = _csv(REF / f"{spec['file']}_realtime.csv", spec["xc"], spec["yc"])
        s81 = _csv(REF / f"{spec['file']}_ge_status81.csv", spec["xc"], spec["yc"])
        unb = _csv(REF / f"{spec['file']}_ge_unbalanced.csv", spec["xc"], spec["yc"])

        # A dense ladder for the drawn curve, and the *reference's own abscissae* for the
        # score -- trimmed at each one and dropped where the solve will not go, which is
        # what `test_steady_sweeps` does. Interpolating a dense curve instead silently
        # extrapolates past the end of the compressor map and scores a point that is not
        # reachable: Figure 6's last marker is 100.31 %NG against `f1`'s 100 % top line.
        dense = np.arange(130.0, 815.0, 5.0)
        grid = tss._ladder(dense)
        gx = np.array(sorted(grid))
        if no == 8:
            ours_x = np.array([grid[k]["ng_pct"] for k in gx])
            ours_y = np.array([grid[k]["ps3"] for k in gx])
            inside = (rt[0] >= ours_x.min()) & (rt[0] <= ours_x.max())
            at_x, ref_x, ref_y = (
                np.interp(rt[0][inside], ours_x, ours_y),
                rt[0][inside],
                rt[1][inside],
            )
            dropped = int((~inside).sum())
        else:
            ours_x, ours_y = gx, np.array([grid[k][spec["key"]] for k in gx])
            at = tss._ladder(rt[0])
            keep = [i for i, k in enumerate(rt[0]) if k in at]
            ref_x, ref_y = rt[0][keep], rt[1][keep]
            at_x = np.array([at[k][spec["key"]] for k in ref_x])
            dropped = int(rt[0].size - len(keep))

        ax.plot(unb[0], unb[1], "s", color=GE2, ms=4, alpha=0.75, zorder=2, label="GE unbalanced")
        ax.plot(s81[0], s81[1], "^", color=GE, ms=4, alpha=0.75, zorder=2, label="GE status 81")
        ax.plot(
            rt[0],
            rt[1],
            "o",
            color=BALLIN,
            ms=5,
            mec=THEME.ground,
            mew=0.8,
            zorder=4,
            label="Ballin real-time",
        )
        ax.plot(ours_x, ours_y, "-", color=OURS, lw=2.2, zorder=3, label="our model")
        # Table B.1's three trims, on the same axes. They are printed *numbers* against a
        # digitized plot, so where the diamond misses the circle the report disagrees with
        # itself and no model can sit on both -- see `internal_spread` and SCOPE.md, "Which
        # source wins". Figure 8 carries no fuel flow, so its abscissa is NG instead.
        b1x, b1y = _b1_markers(no)
        ax.plot(
            b1x,
            b1y,
            "D",
            color=THEME.ink,
            ms=7,
            mec=THEME.ground,
            mew=1.2,
            zorder=5,
            label="Table B.1 (printed numbers)",
        )
        _style(ax, spec["xlabel"], spec["ylabel"], spec["title"])
        ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="best")

        if no == 6:
            dev = at_x - ref_y  # %NG points, an absolute difference
            unit = "%NG"
        else:
            dev = 100.0 * (at_x - ref_y) / np.abs(ref_y)
            unit = "%"
        out[f"figure_{no}"] = {
            **_stats(dev),
            "unit": unit,
            "against": "Ballin real-time",
            "printed_points": int(rt[0].size),
            "dropped_unreachable": dropped,
            "worst_at": float(ref_x[int(np.argmax(np.abs(dev)))]),
        }

    fig.tight_layout()
    save(fig, "figures-6-8-overlay")
    return out


# ========================================================= Figures 9-10, the fuel transients

PANEL_LABEL = {
    "pcng": ("gas generator speed", "%NG"),
    "ps3": ("compressor discharge Ps3", "psia"),
    "t41": ("turbine inlet T4.1", "deg R"),
    "t45": ("power turbine inlet T4.5", "deg R"),
    "torq45": ("power turbine torque", "ft*lbf"),
    "wfph": ("fuel flow (the INPUT)", "lbm/hr"),
}


def transients() -> dict:
    """Overlay all six panels of each transient figure, and score the five comparable ones.

    The sixth, `WFPH`, is the input: its two levels are printed in the caption, so drawing
    it is what shows the step our model was given is the step Ballin's was.
    """
    import test_fuel_step as tfs

    out = {}
    for no, wf_hi in ((9, 775.0), (10, 125.0)):
        tr = tfs._run(no, wf_hi, tfs.STEP_TIME[no])
        t = tr["t"]
        fig, axes = plt.subplots(2, 3, figsize=fs(16.5, 8.2))
        panels = {}

        for ax, key in zip(axes.ravel(), list(PANEL_LABEL), strict=True):
            mp = REF / f"fig{no:02d}_{key}_model.csv"
            gp = REF / f"fig{no:02d}_{key}_reference.csv"
            label, unit = PANEL_LABEL[key]
            if mp.exists():
                tb, vb = _csv(mp, "t_s", "value")
                ax.plot(tb, vb, "-", color=BALLIN, lw=1.6, zorder=3, label="Ballin real-time")
            if gp.exists():
                tg, vg = _csv(gp, "t_s", "value")
                ax.plot(tg, vg, "+", color=GE, ms=5, mew=1.2, zorder=2, label="GE reference")
            if key == "wfph":
                ax.plot(
                    t,
                    np.where(t < tfs.STEP_TIME[no], 400.0, wf_hi),
                    "-",
                    color=OURS,
                    lw=2,
                    zorder=4,
                    label="our input",
                )
            else:
                ours = tfs.PANELS[key](tr)
                ax.plot(t, ours, "-", color=OURS, lw=2, zorder=4, label="our model")
                if mp.exists():
                    late, t_last = tfs._settled_window(tb, vb)
                    if late.size:
                        mine = ours[(t > 4.0) & (t <= t_last)]
                        pre_b = vb[tb < tfs.STEP_TIME[no] - 0.05]
                        pre_u = ours[t < tfs.STEP_TIME[no] - 0.05]
                        panels[key] = {
                            "settled_dev_pct": float(
                                100.0 * (mine.mean() - late.mean()) / abs(late.mean())
                            ),
                            "initial_dev_pct": float(
                                100.0 * (pre_u.mean() - pre_b.mean()) / abs(pre_b.mean())
                            )
                            if pre_b.size
                            else None,
                            "ours_settled": float(mine.mean()),
                            "ballin_settled": float(late.mean()),
                            "ref_ends_s": float(t_last),
                        }
            _style(ax, "time, s", unit, f"{label}")
            ax.legend(frameon=False, fontsize=7.5, labelcolor=MUTED, loc="best")

        plotstyle.title(
            fig,
            f"Figure {no} [pdf p.{ {9: 45, 10: 46}[no] }] -- "
            f"fuel step 400 -> {wf_hi:.0f} lbm/hr, heat sink on",
        )
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        save(fig, f"figure-{no}-overlay")
        out[f"figure_{no}"] = panels
    return out


def whole_curve() -> dict:
    """Shape agreement over the whole record, normalised by each panel's own excursion."""
    import test_whole_curve as twc

    out = {}
    for no in (9, 10):
        t, panels = twc._run(no, heat_sink=True)
        for key in twc.PANELS:
            tb, vb = twc._ballin(no, key)
            if tb is None:
                continue
            out[f"fig{no}_{key}"] = float(twc._rms_pct(t, panels[key], tb, vb, no))
    return out


# ================================================================== Table B.1, pdf p.67


NP_TRIM_RPM = 20895.0
"""Table B.1's power turbine speed at all three trims [pdf p.67].

**Not `NP_DES`, which is 20900.** The five rpm between them is 0.024 %, which is nothing on
a speed and not quite nothing on a shaft power computed from it, and `validation/
test_trim_points.py` has always used the printed value. This file used `NP_DES` for a few
hours on 2026-09-14 and the Table B.1 rms read 0.2527 % instead of 0.2521."""


def table_b1() -> dict:
    """The complete printed engine state at three trims: 21 numbers, open loop."""
    rows, devs = [], []
    for i, name in enumerate(TRIMS, start=1):
        p = STATE_B1[name]
        wf = wf_pps_from_pph(p["wf"])
        r = trim.solve(wf, NP_TRIM_RPM, AMB)
        f = frame(r.state, wf, AMB)
        got = {
            "ng": r.state.ng_rpm,
            # Ps3 is the *static* discharge pressure Table B.1 prints; P3 is the total.
            "ps3": c.K_PS3 * r.state.p3_psia,
            "p41": r.state.p41_psia,
            "t41": f.t41_degR,
            "p45": r.state.p45_psia,
            "t45": f.t45_degR,
            "shp": shp_from_torque(f.q_pt_ftlbf, NP_TRIM_RPM),
        }
        for k, v in got.items():
            d = 100.0 * (v - p[k]) / abs(p[k])
            devs.append(d)
            rows.append(
                {
                    "trim": name,
                    "quantity": k.upper(),
                    "printed": p[k],
                    "ours": float(v),
                    "dev_pct": float(d),
                }
            )
        _ = i
    return {"rows": rows, **_stats(np.array(devs))}


def internal_spread(b1: dict) -> dict:
    """How far the report disagrees with *itself*, against how far we sit from it.

    Table B.1 prints numbers; Figures 6 and 7 print the same three operating points as
    plotted markers. They do not agree, and the gap widens as power falls -- SCOPE.md's
    "Which source wins" rule and open question #46 both turn on it, and the practical
    consequence is a discipline: **do not tune the model below the report's own internal
    spread.** That made this a governing number stated in prose and computed nowhere, so
    it is computed here, from the same digitized figures every other comparison uses.

    Ours is the deviation from Table B.1, which carries no digitizing error of ours at
    all; theirs is Table B.1 against the figure, which carries only theirs and ours of
    the figure. A row where `ours_pct` is the smaller of the two is a row where chasing
    the residual further is chasing noise in the source.
    """
    fig6 = _csv(REF / "fig06_realtime.csv", "wf_pph", "ng_pct")
    fig7 = _csv(REF / "fig07_realtime.csv", "wf_pph", "shp")
    ours = {(r["trim"], r["quantity"]): r["dev_pct"] for r in b1["rows"]}

    rows = []
    for name, p in STATE_B1.items():
        b1_ng_pct = 100.0 * p["ng"] / c.NG_DES
        for quantity, printed, curve, unit in (
            ("NG", b1_ng_pct, fig6, "%NG"),
            ("SHP", p["shp"], fig7, "%"),
        ):
            # The figures are sampled at their own abscissae -- 15 markers across Figure 7's
            # whole range -- so a trim's fuel flow almost never lands on one and the value
            # is interpolated. That interpolation is itself a source of error and it is
            # reported rather than buried: `nearest_marker_pph` is how far the trim sits
            # from the closest printed marker, and on Figure 7's steep low-power end the
            # curve climbs ~3 shp per lbm/hr, so 5 lbm/hr of it is 1.6 % of hover power.
            # Read a row whose `nearest_marker_pph` is large as the weaker of the two.
            fig = float(np.interp(p["wf"], curve[0], curve[1]))
            near = float(np.min(np.abs(curve[0] - p["wf"])))
            theirs = printed - fig if unit == "%NG" else 100.0 * (printed - fig) / abs(printed)
            rows.append(
                {
                    "trim": name,
                    "quantity": quantity,
                    "wf_pph": p["wf"],
                    "table_b1": printed,
                    "figure": fig,
                    "unit": unit,
                    "theirs_pct": float(theirs),
                    "ours_pct": float(ours[(name, quantity)]),
                    "nearest_marker_pph": near,
                }
            )
    inside = [r for r in rows if abs(r["ours_pct"]) < abs(r["theirs_pct"])]
    return {
        "rows": rows,
        "n": len(rows),
        "widest_theirs_pct": max(rows, key=lambda r: abs(r["theirs_pct"]))["theirs_pct"],
        "ours_inside_their_spread": len(inside),
        "figure_6": "pdf p.40",
        "figure_7": "pdf p.41",
        "table_b1": "pdf p.67",
    }


# =================================================== Table 1 and every other eigenvalue

DOF_FIGS = {2: DOF.TWO, 3: DOF.THREE, 5: DOF.FIVE, 6: DOF.SIX}
MODE_ORDER = ("NG", "NP", "P3", "P41", "P45")


def _our_eigs(dof: DOF, t: int):
    wf = wf_pps_from_pph(WF_PPH[t])
    r = trim.solve(wf, c.NP_DES, AMB)
    m = extract(r, wf, dof, AMB, j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[t])
    return np.linalg.eigvals(m.A), m


def _sorted_real(ev):
    """Real parts, least negative first -- the order Table 1 prints its modes in."""
    return np.sort(np.asarray(ev).real)[::-1]


def eigenvalues() -> dict:
    """Every model variant against Table 1, and against Ballin's own printed matrices.

    Two references exist and they are not the same thing.

    **Table 1 [pdf p.31]** prints 27 eigenvalues, for the 5-DOF, the 2-DOF and the
    order-reduced 5-DOF. It prints none for the 3-DOF or the 6-DOF.

    **The Appendix B matrices themselves** imply a spectrum for all twelve figures, so the
    3-DOF and 6-DOF *can* be compared by eigenvalue even though the report tabulates
    none -- with one caution that is not optional. `SCOPE.md` and `t700.appendix_b` both
    record that the 6-DOF matrices are ill-conditioned as printed: row 6 differences terms
    of order 2e5 to give order 1e3, so four printed digits carry about +/-50 there, and
    B8 comes out **unstable at +3.07 /sec**. That is lost precision in the printing, not a
    claim by the report. It is reported here because hiding it would be worse, and it is
    labelled at every point it appears.
    """
    out = {"table_1": [], "against_printed_matrices": [], "notes": {}}

    for dofno, col in ((5, "5dof"), (2, "2dof"), (None, "red5")):
        dof = DOF_FIGS.get(dofno, DOF.REDUCED_FIVE)
        for t in (1, 2, 3):
            ev = _sorted_real(_our_eigs(dof, t)[0])
            printed = TABLE_1[t][col]
            pv = _sorted_real(np.array(list(printed.values())))
            names = [k for k, _ in sorted(printed.items(), key=lambda kv: -kv[1])]
            for i, nm in enumerate(names):
                ours_v = float(ev[i])
                out["table_1"].append(
                    {
                        "model": col,
                        "trim": TRIMS[t - 1],
                        "mode": nm,
                        "printed": float(pv[i]),
                        "ours": ours_v,
                        "dev_pct": float(100.0 * (ours_v - pv[i]) / abs(pv[i])),
                    }
                )

    for dofno, dof in DOF_FIGS.items():
        for t in (1, 2, 3):
            ev, _ = _our_eigs(dof, t)
            ref = ab.find(dofno, t)
            bal = np.linalg.eigvals(ref.A)
            ours = _sorted_real(ev)
            theirs = _sorted_real(bal)
            out["against_printed_matrices"].append(
                {
                    "model": f"{dofno}-DOF",
                    "figure": ref.figure,
                    "trim": TRIMS[t - 1],
                    "ours": [float(v) for v in ours],
                    "ballin": [float(v) for v in theirs],
                    "dev_pct": [
                        float(100.0 * (a - b) / abs(b)) for a, b in zip(ours, theirs, strict=True)
                    ],
                    "ballin_complex": bool(np.max(np.abs(bal.imag)) > 1e-9),
                    "ballin_unstable": bool(np.max(bal.real) > 0.0),
                    "ours_complex": bool(np.max(np.abs(np.asarray(ev).imag)) > 1e-9),
                    "ours_unstable": bool(np.max(np.asarray(ev).real) > 0.0),
                    "ill_conditioned": dofno == 6,
                }
            )

    t1 = np.array([r["dev_pct"] for r in out["table_1"]])
    out["notes"]["table_1_stats"] = _stats(t1)
    out["notes"]["table_1_worst_row"] = max(out["table_1"], key=lambda r: abs(r["dev_pct"]))
    return out


def eigenvalue_plot(data: dict) -> None:
    """Table 1's 27 printed modes against ours, and the two spectra the report omits."""
    fig, axes = plt.subplots(1, 2, figsize=fs(15.5, 5.4), width_ratios=[1.35, 1.0])

    ax = axes[0]
    rows = data["table_1"]
    labels = [f"{r['model'][:4]} {r['trim'].split()[0][:4]} {r['mode']}" for r in rows]
    idx = np.arange(len(rows))
    ax.barh(
        idx + 0.2,
        [-r["printed"] for r in rows],
        height=0.38,
        color=BALLIN,
        zorder=3,
        label="Table 1, printed",
    )
    ax.barh(
        idx - 0.2,
        [-r["ours"] for r in rows],
        height=0.38,
        color=OURS,
        zorder=3,
        label="our Jacobian",
    )
    ax.set_xscale("log")
    ax.set_yticks(idx)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.invert_yaxis()
    for i, r in enumerate(rows):
        if abs(r["dev_pct"]) > 3.0:
            ax.text(
                max(-r["printed"], -r["ours"]) * 1.2,
                i,
                f"{r['dev_pct']:+.1f}%",
                va="center",
                fontsize=6,
                color=MUTED,
            )
    _style(
        ax,
        "|eigenvalue|, 1/sec  (log)",
        "",
        "All 27 eigenvalues Table 1 prints [pdf p.31], against ours",
    )
    ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="lower right")

    ax = axes[1]
    got = [r for r in data["against_printed_matrices"] if r["model"] in ("3-DOF", "6-DOF")]
    y, lab = [], []
    for k, r in enumerate(got):
        for a, b in zip(r["ours"], r["ballin"], strict=True):
            ax.plot([-b], [k], "o", color=BALLIN, ms=6, zorder=3)
            ax.plot([-a], [k], "o", color=OURS, ms=6, mfc="none", mew=1.6, zorder=4)
        y.append(k)
        lab.append(
            f"{r['figure']} {r['model']} {r['trim'].split()[0][:4]}"
            + ("  ill-cond." if r["ill_conditioned"] else "")
        )
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=7)
    ax.invert_yaxis()
    ax.plot([], [], "o", color=BALLIN, ms=6, label="Ballin's printed matrix")
    ax.plot([], [], "o", color=OURS, ms=6, mfc="none", mew=1.6, label="ours")
    ax.axvline(0.0, color="#b4453a", lw=1.1, ls=":", zorder=2)
    for k, r in enumerate(got):
        if r["ballin_unstable"]:
            worst = min(r["ballin"])
            ax.annotate(
                "B8 as printed is UNSTABLE at +3.07 /sec:\nfour digits carry ~+/-50 in row 6",
                xy=(-worst, k),
                xytext=(-worst * 30, k + 0.55),
                fontsize=7,
                color="#b4453a",
                arrowprops=dict(arrowstyle="->", color="#b4453a", lw=0.9),
            )
    _style(
        ax,
        "-eigenvalue, 1/sec  (symlog).  Left of the dotted line is unstable",
        "",
        "The two spectra Table 1 does NOT print, from the matrices themselves",
    )
    ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="upper left")

    fig.tight_layout()
    save(fig, "eigenvalues-all-models")


# ============================================================ Appendix B, 297 elements


def appendix_b() -> dict:
    """Every printed element against ours, as a scatter and as per-block statistics."""
    figs = {
        (2, 1): 1,
        (2, 2): 3,
        (2, 3): 5,
        (3, 1): 7,
        (3, 2): 9,
        (3, 3): 11,
        (5, 1): 2,
        (5, 2): 4,
        (5, 3): 6,
        (6, 1): 8,
        (6, 2): 10,
        (6, 3): 12,
    }
    pairs, blocks = [], {}
    a_dev, b_dev = [], []
    for (dofno, t), _ in sorted(figs.items()):
        dof = DOF_FIGS[dofno]
        wf = wf_pps_from_pph(WF_PPH[t])
        r = trim.solve(wf, c.NP_DES, AMB)
        m = extract(r, wf, dof, AMB, j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[t])
        ref = ab.find(dofno, t)
        n = len(ref.states)
        for i in range(n):
            for j in range(n):
                if ref.A[i, j] == 0.0:
                    continue
                d = 100.0 * (m.A[i, j] - ref.A[i, j]) / abs(ref.A[i, j])
                pairs.append((ref.A[i, j], m.A[i, j], f"{dofno}-DOF A"))
                a_dev.append(d)
                blocks.setdefault(f"{dofno}-DOF A", []).append(d)
        for i in range(n):
            if ref.b[i] == 0.0:
                continue
            d = 100.0 * (m.b[i] - ref.b[i]) / abs(ref.b[i])
            pairs.append((ref.b[i], m.b[i], f"{dofno}-DOF b"))
            b_dev.append(d)
            blocks.setdefault(f"{dofno}-DOF b", []).append(d)

    fig, axes = plt.subplots(1, 2, figsize=fs(14.5, 5.6))
    ax = axes[0]
    px = np.array([abs(p[0]) for p in pairs])
    py = np.array([abs(p[1]) for p in pairs])
    lim = [min(px.min(), py.min()) * 0.6, max(px.max(), py.max()) * 1.6]
    ax.plot(lim, lim, "-", color=MUTED, lw=1, zorder=2)
    for tag, col in (("A", OURS), ("b", GE)):
        sel = [k for k, p in enumerate(pairs) if p[2].endswith(tag)]
        ax.plot(
            px[sel],
            py[sel],
            "o",
            ms=4,
            color=col,
            alpha=0.65,
            zorder=3,
            label=f"{tag} elements ({len(sel)})",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    _style(
        ax, "|printed|, Appendix B", "|ours|", "All 206 non-zero printed elements, twelve figures"
    )
    ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="upper left")

    ax = axes[1]
    keys = sorted(blocks)
    ax.axvline(0.0, color=MUTED, lw=1, zorder=2)
    for k, key in enumerate(keys):
        v = np.array(blocks[key])
        ax.plot(
            v,
            np.full(v.size, k) + np.random.default_rng(0).normal(0, 0.06, v.size),
            "o",
            ms=4,
            color=OURS if key.endswith("A") else GE,
            alpha=0.6,
            zorder=3,
        )
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels(keys, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(-60, 60)
    _style(
        ax,
        "deviation from printed, %",
        "",
        "By block. The 6-DOF b is the one with a shape, and it is the heat sink's",
    )
    fig.tight_layout()
    save(fig, "appendix-b-elements")

    return {
        "A": _stats(np.array(a_dev)),
        "b": _stats(np.array(b_dev)),
        "by_block": {k: _stats(np.array(v)) for k, v in sorted(blocks.items())},
    }


# ======================================================== the closed loop, Phase 5's gate


def closed_loop() -> dict:
    """Table B.1 again, with fuel flow as an *output* -- the hardest test in the project."""
    import test_closed_loop as tcl

    rows, devs = [], []
    for name, wf, ng, ps3, q, ratio, shp in tcl.TRIMS:
        mean, ptp = tcl.settle_cycle(wf, q, ratio, name=name)
        printed = {"NG": ng, "NP": tcl.NP_TRIM_RPM, "Wf": wf, "Ps3": ps3, "shp": shp}
        for k, p in printed.items():
            d = 100.0 * (mean[k] - p) / abs(p)
            devs.append(d)
            rows.append(
                {
                    "trim": name,
                    "quantity": k,
                    "printed": p,
                    "ours": float(mean[k]),
                    "dev_pct": float(d),
                    "cycle_ptp": float(ptp[k]),
                }
            )
    return {"rows": rows, **_stats(np.array(devs))}


def phase_plane() -> dict:
    """Ps3 against NG -- the comparison that discards the time axis, and with it #37."""
    import test_fuel_step as tfs

    fig, axes = plt.subplots(1, 2, figsize=fs(13.5, 5.4))
    out = {}
    for ax, (no, wf_hi) in zip(axes, ((9, 775.0), (10, 125.0)), strict=True):
        tr = tfs._run(no, wf_hi, tfs.STEP_TIME[no])
        ng = 100.0 * tr["ng"] / c.NG_DES
        ps3 = c.K_PS3 * tr["p3"]
        tb, vb = _csv(REF / f"fig{no:02d}_pcng_model.csv", "t_s", "value")
        tp, vp = _csv(REF / f"fig{no:02d}_ps3_model.csv", "t_s", "value")
        bal_ps3 = np.interp(tb, tp, vp)
        eq = _csv(REF / "fig08_realtime.csv", "ng_pct", "ps3_psia")
        ax.plot(
            eq[0], eq[1], "--", color=MUTED, lw=1.2, zorder=2, label="Figure 8 equilibrium locus"
        )
        ax.plot(vb, bal_ps3, "-", color=BALLIN, lw=1.8, zorder=3, label="Ballin")
        ax.plot(ng, ps3, "-", color=OURS, lw=2, zorder=4, label="ours")
        _style(
            ax,
            "gas generator speed, %NG",
            "Ps3, psia",
            f"Figure {no} in the phase plane -- no time axis",
        )
        ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="best")
        lo, hi = max(vb.min(), ng.min()), min(vb.max(), ng.max())
        grid = np.linspace(lo + 0.5, hi - 0.5, 40)
        ours_i = np.interp(grid, ng[np.argsort(ng)], ps3[np.argsort(ng)])
        o = np.argsort(vb)
        bal_i = np.interp(grid, vb[o], bal_ps3[o])
        out[f"figure_{no}"] = _stats(100.0 * (ours_i - bal_i) / np.abs(bal_i))
    fig.tight_layout()
    save(fig, "phase-plane")
    return out


def table_b1_plot(b1: dict, cl: dict) -> None:
    """Open loop and closed loop against the same 21 printed numbers."""
    fig, ax = plt.subplots(figsize=fs(13.0, 5.2))
    rows = b1["rows"]
    labels = [f"{r['trim'].split()[0][:4]} {r['quantity']}" for r in rows]
    idx = np.arange(len(rows))
    ax.bar(
        idx,
        [r["dev_pct"] for r in rows],
        width=0.55,
        color=OURS,
        zorder=3,
        label="open loop, Wf prescribed",
    )
    # the closed loop reports NG / NP / Wf / Ps3 / shp; the open-loop rows are uppercased
    # station names. Only the three that exist on both sides can be overlaid -- Wf is the
    # open loop's *input*, and NP is held there.
    cl_map = {(r["trim"], r["quantity"].upper()): r["dev_pct"] for r in cl["rows"]}
    xs = [i for i, r in enumerate(rows) if (r["trim"], r["quantity"]) in cl_map]
    ax.plot(
        xs,
        [cl_map[(rows[i]["trim"], rows[i]["quantity"])] for i in xs],
        "D",
        color=BALLIN,
        ms=6,
        zorder=4,
        label="closed loop, Wf an output",
    )
    ax.axhline(0.0, color=MUTED, lw=1, zorder=2)
    for s in (0.5, -0.5):
        ax.axhline(s, color=GRID, lw=1.2, ls="--", zorder=2)
    ax.set_xticks(idx)
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    _style(
        ax,
        "",
        "deviation from printed, %",
        "Table B.1 [pdf p.67], 21 printed numbers. Dashed lines are our +/-0.5 % bound",
    )
    ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="best")
    fig.tight_layout()
    save(fig, "table-b1")


def _render(data: dict) -> None:
    print("figures 6-8 ...", flush=True)
    data["steady"] = steady_sweeps()
    print("figures 9-10 ...", flush=True)
    data["transient"] = transients()
    print("whole curve ...", flush=True)
    data["whole_curve"] = whole_curve()
    print("phase plane ...", flush=True)
    data["phase_plane"] = phase_plane()
    print("table B.1 ...", flush=True)
    data["table_b1"] = table_b1()
    print("closed loop ...", flush=True)
    data["closed_loop"] = closed_loop()
    table_b1_plot(data["table_b1"], data["closed_loop"])
    data["internal_spread"] = internal_spread(data["table_b1"])
    print("eigenvalues ...", flush=True)
    data["eigenvalues"] = eigenvalues()
    eigenvalue_plot(data["eigenvalues"])
    print("appendix B ...", flush=True)
    data["appendix_b"] = appendix_b()
    print("engine operating line ...", flush=True)
    data["engine_steady"] = sheet_engine_steady()
    print("open loop family ...", flush=True)
    data["open_loop"] = sheet_open_loop()
    print("closed loop response ...", flush=True)
    data["closed_loop_response"] = sheet_closed_loop()
    print("control schedules ...", flush=True)
    data["schedules"] = sheet_schedules()
    print("engine maps ...", flush=True)
    data["engine_maps"] = sheet_engine_maps()


# ============================================ the model's own behaviour, with no reference
#
# Nothing below is a comparison. These are the engine and its control as this replication
# runs them, across the range the report declares the model valid over -- which is what a
# reader needs in order to judge whether the agreements above are agreements about
# something physical. `SCOPE.md` records the ceiling: NG is claimed to 100 % and no further.


def sheet_engine_steady() -> dict:
    """The equilibrium operating line, station by station, over the whole trim range."""
    import test_steady_sweeps as tss

    pph = np.arange(130.0, 812.0, 4.0)
    grid = tss._ladder(pph)
    ks = np.array(sorted(grid))
    rows = []
    for k in ks:
        wf = wf_pps_from_pph(float(k))
        r = trim.solve(wf, c.NP_DES, AMB, guess=None)
        if not r.trustworthy:
            continue
        f = frame(r.state, wf, AMB)
        rows.append(
            dict(
                wf=float(k),
                ng=100.0 * r.state.ng_rpm / c.NG_DES,
                ps3=c.K_PS3 * r.state.p3_psia,
                p41=r.state.p41_psia,
                p45=r.state.p45_psia,
                t3=f.t3_degR,
                t41=f.t41_degR,
                t45=f.t45_degR,
                t49=f.t49_degR,
                wa2=f.wa2_pps,
                wa31=f.wa31_pps,
                far=f.far,
                shp=shp_from_torque(f.q_pt_ftlbf, c.NP_DES),
                pr=c.K_PS3 * r.state.p3_psia / AMB.p_amb_psia,
                bleed=100.0 * (f.wa24_bl_pps + f.wa3_bl_pps) / f.wa2_pps,
            )
        )
    g = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    g["sfc"] = g["wf"] / np.maximum(g["shp"], 1e-9)

    fig, axes = plt.subplots(2, 3, figsize=fs(16.5, 8.0))
    panels = [
        (
            "station pressures",
            "psia",
            [("Ps3", "ps3", OURS), ("P41", "p41", BALLIN), ("P45", "p45", GE)],
        ),
        (
            "station temperatures",
            "deg R",
            [
                ("T3", "t3", GE),
                ("T4.1", "t41", BALLIN),
                ("T4.5", "t45", OURS),
                ("T4.9", "t49", GE2),
            ],
        ),
        ("air path", "lbm/sec", [("Wa2 inlet", "wa2", OURS), ("Wa31 to burner", "wa31", BALLIN)]),
        ("shaft power", "hp", [("SHP", "shp", OURS)]),
        ("specific fuel consumption, above 50 hp", "lbm/hr per hp", [("SFC", "sfc", BALLIN)]),
        ("customer + cooling bleed", "% of inlet flow", [("bleed", "bleed", GE)]),
    ]
    # SFC is Wf/SHP and shaft power crosses zero near the bottom of the range -- it reaches
    # -17 hp at 68 %NG, which is the free turbine absorbing rather than delivering. Plotting
    # the quotient through that is meaningless, so the panel starts where there is power.
    usable = g["shp"] > 50.0
    for ax, (title, unit, series) in zip(axes.ravel(), panels, strict=True):
        for label, key, col in series:
            m = usable if key == "sfc" else np.ones_like(usable, dtype=bool)
            ax.plot(g["ng"][m], g[key][m], "-", color=col, lw=2, zorder=3, label=label)
        _style(ax, "gas generator speed, %NG", unit, title)
        if len(series) > 1:
            ax.legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="best")
    plotstyle.title(
        fig,
        "The equilibrium operating line, 130 to 810 lbm/hr. No reference data -- this is "
        "the model's own behaviour",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    save(fig, "engine-operating-line")

    return {
        "ng_range_pct": [float(g["ng"].min()), float(g["ng"].max())],
        "wf_range_pph": [float(g["wf"].min()), float(g["wf"].max())],
        "shp_range": [float(g["shp"].min()), float(g["shp"].max())],
        "pressure_ratio_range": [float(g["pr"].min()), float(g["pr"].max())],
        "t41_range_degR": [float(g["t41"].min()), float(g["t41"].max())],
        "bleed_pct_range": [float(g["bleed"].min()), float(g["bleed"].max())],
        "sfc_best": float(g["sfc"][g["shp"] > 50.0].min()),
        "sfc_best_at_ng": float(
            g["ng"][g["shp"] > 50.0][int(np.argmin(g["sfc"][g["shp"] > 50.0]))]
        ),
        "shp_zero_crossing_ng": float(np.interp(0.0, g["shp"], g["ng"])),
        "points": len(rows),
    }


def sheet_open_loop() -> dict:
    """A family of fuel steps from the same trim, both heat-sink configurations."""
    from t700 import realtime

    wf0 = wf_pps_from_pph(400.0)
    r0 = trim.solve(wf0, c.NP_DES, AMB)
    f0 = frame(r0.state, wf0, AMB)
    steps = [(125.0, "#9c5bd0"), (250.0, GE), (550.0, OURS), (775.0, BALLIN)]

    fig, axes = plt.subplots(1, 4, figsize=fs(17.0, 4.3))
    out = {}
    for pph, col in steps:
        st = realtime.from_trim(r0, f0.wa31_pps, f0)
        hi = wf_pps_from_pph(pph)
        tr = realtime.run(
            st,
            lambda t, v=hi: wf0 if t < 0.5 else v,
            AMB,
            duration_s=5.0,
            dt=0.007,
            q_req_ftlbf=f0.q_pt_ftlbf,
            integrate_np=False,
            heat_sink=True,
        )
        t = tr["t"]
        lab = f"400 -> {pph:.0f}"
        for ax, key, fn in (
            (axes[0], "ng", lambda v: 100.0 * v / c.NG_DES),
            (axes[1], "t41", lambda v: v),
            (axes[2], "q_pt", lambda v: v),
            (axes[3], "p3", lambda v: c.K_PS3 * v),
        ):
            ax.plot(t, fn(np.asarray(tr[key])), "-", color=col, lw=1.8, zorder=3, label=lab)
        ng = 100.0 * np.asarray(tr["ng"]) / c.NG_DES
        out[f"step_{pph:.0f}"] = {
            "ng_start": float(ng[0]),
            "ng_end": float(ng[-1]),
            "t41_peak": float(np.max(tr["t41"])),
            "t41_end": float(tr["t41"][-1]),
            "t41_overshoot": float(np.max(tr["t41"]) - tr["t41"][-1]),
        }
    for ax, ylab, title in (
        (axes[0], "%NG", "gas generator speed"),
        (axes[1], "deg R", "turbine inlet T4.1"),
        (axes[2], "ft*lbf", "power turbine torque"),
        (axes[3], "psia", "compressor discharge Ps3"),
    ):
        _style(ax, "time, s", ylab, title)
    axes[0].legend(frameon=False, fontsize=7.5, labelcolor=MUTED, loc="best", title="lbm/hr")
    plotstyle.title(
        fig,
        "Open loop: four fuel steps from the same 400 lbm/hr trim, heat sink on, "
        "power turbine held at design speed",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "open-loop-family")
    return out


def sheet_closed_loop() -> dict:
    """The governor doing its job: a load step, a collective slam, a reference change."""
    import test_closed_loop as tcl
    from t700.control import loop

    name, wf, _ng, _ps3, q, ratio, _shp = tcl.TRIMS[0]
    wf_pps = wf_pps_from_pph(wf)
    r = trim.solve(wf_pps, c.NP_DES, AMB)
    f = frame(r.state, wf_pps, AMB)
    slope = abs(tcl._dqpt_dnp(r, wf_pps)) * ratio
    base = dict(pas_deg=100.0, pcprf_pct=tcl.NP_TRIM_RPM * 100.0 / c.NP_DES)

    def run(n, at, xcpc_of, load_of, ref_of=None):
        s, rec = None, {k: [] for k in ("t", "np", "ng", "wf", "t41", "spdg", "lim")}
        for i in range(n):
            pilot = loop.Pilot(
                xcpc_pct=xcpc_of(i),
                pas_deg=100.0,
                pcprf_pct=(ref_of(i) if ref_of else base["pcprf_pct"]),
            )
            if s is None:
                s = loop.seed(r, f.wa31_pps, f, pilot, AMB)
            s, e, h, fr = loop.step(
                s, pilot, AMB, dt=tcl.DT, load=load_of(i), j_load=c.J_LOAD_UH60A, heat_sink=True
            )
            rec["t"].append(i * tcl.DT)
            rec["np"].append(s.engine.np_rpm)
            rec["ng"].append(100.0 * s.engine.ng_rpm / c.NG_DES)
            rec["wf"].append(h.wf_pph)
            rec["t41"].append(fr.t41_degR)
            rec["spdg"].append(e.spdg)
            rec["lim"].append(h.limit)
        _ = at
        return {k: (np.asarray(v) if k != "lim" else v) for k, v in rec.items()}

    xc = tcl.XCPC_PCT[name]
    cases = [
        (
            "15 % load step",
            run(
                2600,
                1400,
                lambda i: xc,
                lambda i: (
                    lambda v, s=1.15 if i >= 1400 else 1.0: q * s + slope * (v - tcl.NP_TRIM_RPM)
                ),
            ),
        ),
        (
            "collective 52.75 -> 70 %",
            run(
                1800,
                600,
                lambda i: 70.0 if i >= 600 else xc,
                lambda i: lambda v: q + slope * (v - tcl.NP_TRIM_RPM),
            ),
        ),
        (
            "speed reference 100 -> 102 %",
            run(
                2600,
                1200,
                lambda i: xc,
                lambda i: lambda v: q + slope * (v - tcl.NP_TRIM_RPM),
                lambda i: base["pcprf_pct"] * (1.02 if i >= 1200 else 1.0),
            ),
        ),
    ]

    fig, axes = plt.subplots(3, 4, figsize=fs(17.0, 10.2))
    out = {}
    for row, (label, d) in zip(axes, cases, strict=True):
        for ax, key, ylab, title, col in (
            (row[0], "np", "rpm", "power turbine speed NP", OURS),
            (row[1], "wf", "lbm/hr", "fuel flow, an OUTPUT here", BALLIN),
            (row[2], "ng", "%NG", "gas generator speed", GE),
            (row[3], "spdg", "volts", "ECU torque-motor demand SPDG", GE2),
        ):
            ax.plot(d["t"], d[key], "-", color=col, lw=1.8, zorder=3)
            _style(ax, "time, s", ylab, f"{title}")
        row[0].set_ylabel(f"{label}\n\nrpm", color=MUTED, fontsize=9)
        out[label] = {
            "np_min": float(d["np"].min()),
            "np_max": float(d["np"].max()),
            "np_settled": float(d["np"][-200:].mean()),
            "wf_min": float(d["wf"].min()),
            "wf_max": float(d["wf"].max()),
            "wf_settled": float(d["wf"][-200:].mean()),
            "limits_used": sorted(set(d["lim"])),
        }
    plotstyle.title(
        fig,
        "Closed loop: the Appendix C control system on the engine at the hover trim. "
        "Fuel flow is an output",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.968))
    save(fig, "closed-loop-response")
    return out


def sheet_schedules() -> dict:
    """All eight Appendix C scheduling functions, as digitized. They exist only as plots."""
    from t700.control import schedules as sch

    specs = [
        (
            "F_HM1",
            "Fig. C24, pdf p.95",
            "T2, deg R",
            "WFPTP",
            sch.f_hm1,
            "topping line -- the fuel ceiling",
        ),
        ("F_HM2", "Fig. C25, pdf p.96", "PAS, deg", "WFPRF", sch.f_hm2, "power available spindle"),
        (
            "F_HM3",
            "Fig. C26, pdf p.97",
            "XLDSH, deg",
            "PNG",
            sch.f_hm3,
            "load demand -> NG reference",
        ),
        (
            "F_HM4",
            "Fig. C27, pdf p.98",
            "XLDSH, deg",
            "WFQPS3",
            sch.f_hm4,
            "load demand -> fuel feedforward",
        ),
        (
            "F_HM5",
            "Fig. C28, pdf p.99",
            "T2, deg R",
            "WFIRF",
            sch.f_hm5,
            "idle schedule, fuel term",
        ),
        (
            "F_HM6",
            "Fig. C29, pdf p.100",
            "T2, deg R",
            "PCNGI",
            sch.f_hm6,
            "idle schedule, speed reference",
        ),
    ]
    fig, axes = plt.subplots(2, 4, figsize=fs(17.0, 7.6))
    out = {}
    for ax, (name, src, xl, yl, fn, what) in zip(axes.ravel(), specs, strict=False):
        cur = fn()
        xs = np.linspace(cur.x.min(), cur.x.max(), 400)
        ax.plot(xs, [float(cur(v)) for v in xs], "-", color=OURS, lw=2, zorder=3)
        ax.plot(cur.x, cur.y, "x", color=BALLIN, ms=7, mew=1.6, zorder=4)
        _style(ax, xl, yl, f"{name} -- {what}")
        ax.text(
            0.02,
            0.04,
            f"{cur.x.size} printed markers\n{src}",
            transform=ax.transAxes,
            fontsize=6.8,
            color=MUTED,
            va="bottom",
            zorder=5,
            bbox=dict(fc=THEME.ground, ec="none", alpha=0.82, pad=1.8),
        )
        out[name] = {
            "knots": int(cur.x.size),
            "source": src,
            "x_range": [float(cur.x.min()), float(cur.x.max())],
            "y_range": [float(cur.y.min()), float(cur.y.max())],
        }

    ax = axes.ravel()[6]
    m7 = sch.f_hm7()
    for k, param in enumerate(m7.params):
        line = m7.lines[k]
        ax.plot(
            line.x,
            line.y,
            "-",
            lw=1.6,
            zorder=3,
            color=plt.cm.viridis(k / max(len(m7.params) - 1, 1)),
            label=f"{param:g}",
        )
    _style(ax, "PCNGHL, %NG", "WFPAC", "F_HM7 -- acceleration limit, 7 T2 lines")
    ax.legend(
        frameon=False, fontsize=6, labelcolor=MUTED, ncol=2, title="T2, deg R", title_fontsize=6
    )
    out["F_HM7"] = {
        "lines": int(len(m7.params)),
        "source": "Fig. C30, pdf p.101",
        "params": [float(v) for v in m7.params],
    }

    ax = axes.ravel()[7]
    ec = sch.f_ec1()
    for k, param in enumerate(ec.params):
        line = ec.lines[k]
        ax.plot(
            line.x,
            line.y,
            "-",
            lw=1.8,
            zorder=3,
            color=plt.cm.plasma(0.15 + 0.7 * k / max(len(ec.params) - 1, 1)),
            label=f"{param:g}",
        )
    _style(ax, "W45R", "TAU45, sec", "F_EC1 -- thermocouple time constant")
    ax.legend(frameon=False, fontsize=6, labelcolor=MUTED, title="T45L, deg R", title_fontsize=6)
    out["F_EC1"] = {
        "lines": int(len(ec.params)),
        "source": "Fig. C23, pdf p.94",
        "params": [float(v) for v in ec.params],
    }

    plotstyle.title(
        fig,
        "Appendix C's eight scheduling functions. None exists as a table -- the fuel control "
        "is specified entirely as pictures",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    save(fig, "control-schedules")
    return out


MAP_FILE = {
    "f2": "f2_compressor_temperature.csv",
    "f3": "f3_seal_bleed_fraction.csv",
    "f4": "f4_pt_balance_bleed_fraction.csv",
    "f5": "f5_tip_leak_cooling_bleed_fraction.csv",
    "f6": "f6_combustor_efficiency.csv",
    "f7": "f7_gg_turbine_energy.csv",
    "f8": "f8_pt_energy.csv",
    "f9": "f9_pt_mass_flow.csv",
    "f10": "f10_exhaust_pressure_loss.csv",
    "f_hs": "fhs_heat_sink_constant.csv",
    "f1": "f1_compressor_mass_flow_beta.csv",
}


def _short_axis(text: str) -> str:
    """A printed axis title reduced to its symbol and unit, for a small-multiple panel.

    The provenance headers carry the full printed title plus whatever caveat the digitizer
    needed to record -- `f9`'s runs to 240 characters warning that the facing figure uses a
    different pressure ratio. Set whole on a 3x4 sheet they overrun their panels and collide
    with the neighbours. The symbol is the first token; the unit is the last comma-separated
    clause when it is short and is not the word "nondimensional", which says nothing.
    """
    head = text.split(" -- ")[0].split(" (")[0].strip()
    symbol = head.split(" ")[0].rstrip(",")
    if symbol.endswith("=") or symbol in {"f9"}:
        symbol = head.split(",")[0].split(" ")[0].rstrip(",")
    tail = head.rsplit(",", 1)[-1].strip() if "," in head else ""
    if tail and len(tail) <= 16 and tail.lower() != "nondimensional":
        return f"{symbol}, {tail}"
    return symbol


def _axis_labels(fname: str, directory: str = "maps") -> tuple[str, str]:
    """The printed axis titles, read out of the data file's own provenance header.

    Not hardcoded here. An earlier draft of this sheet typed them from memory and put
    "NGc, %" on `f4` and `f5`, whose abscissa is corrected mass flow -- a caption wrong in a
    way no test can see. The digitizers already record what the page prints; use that, and
    reduce it to what fits with `_short_axis`.
    """
    path = HERE.parent / "data" / directory / fname
    x = y = ""
    for ln in path.read_text().splitlines():
        if ln.startswith("# x:") and not x:
            x = ln[4:].strip()
        elif ln.startswith("# y:") and not y:
            y = ln[4:].strip()
        elif not ln.startswith("#"):
            break
    return (_short_axis(x), _short_axis(y))


def sheet_engine_maps() -> dict:
    """The eleven engine function tables of Appendix A, as digitized and conditioned."""
    from t700 import maps

    single = [
        ("f2", "Fig. A2, p.57", "compressor temperature, Eq. 10"),
        ("f3", "Fig. A3, p.58", "seal-pressurization bleed B1, Eq. 12"),
        ("f4", "Fig. A4, p.59", "power-turbine-balance bleed B2, Eq. 13"),
        ("f5", "Fig. A5, p.60", "tip-leak + cooling bleed B3"),
        ("f6", "Fig. A6, p.61", "combustor efficiency -- a CONSTANT"),
        ("f7", "Fig. A7, p.62", "gas generator turbine energy, Eq. 26"),
        ("f8", "Fig. A8, p.63", "power turbine energy, Eq. 32"),
        ("f9", "Fig. A9, p.64", "power turbine mass flow, Eq. 33"),
        ("f10", "Fig. A10, p.65", "exhaust pressure loss, Eq. 38"),
        ("f_hs", "Fig. A11, p.66", "station 4.1 heat-sink constant, Eq. 52"),
    ]
    fig, axes = plt.subplots(3, 4, figsize=fs(17.0, 10.0))
    out = {}
    for ax, (name, src, what) in zip(axes.ravel(), single, strict=False):
        cur = getattr(maps, name)()
        xl, yl = _axis_labels(MAP_FILE[name])
        xs = np.linspace(cur.x.min(), cur.x.max(), 400)
        ax.plot(xs, [float(cur(v)) for v in xs], "-", color=OURS, lw=2, zorder=3)
        ax.plot(cur.x, cur.y, "x", color=BALLIN, ms=6, mew=1.4, zorder=4)
        _style(ax, xl, yl, f"{name} -- {what}")
        ax.text(
            0.02,
            0.05,
            f"{cur.x.size} knots\n{src}",
            transform=ax.transAxes,
            fontsize=6.6,
            color=MUTED,
            va="bottom",
            zorder=5,
            bbox=dict(fc=THEME.ground, ec="none", alpha=0.82, pad=1.8),
        )
        out[name] = {
            "knots": int(cur.x.size),
            "source": src,
            "x_axis": xl,
            "y_axis": yl,
            "conditioning_moved": float(maps.CONDITIONING.get(name, 0.0)),
        }

    ax = axes.ravel()[10]
    m1 = maps.f1()
    for k, param in enumerate(m1.params):
        line = m1.lines[k]
        ax.plot(
            line.x,
            line.y,
            "-",
            lw=1.3,
            zorder=3,
            color=plt.cm.viridis(k / (len(m1.params) - 1)),
            label=f"{param:g}",
        )
    # f1's derived file is a beta grid and carries no `# x:` header of its own; the axes it
    # is *evaluated* on are the raw extraction's, which is what the panel plots.
    fx, fy = "Ps3/P2", "WA2c, lbm/sec"
    _style(ax, fx, fy, "f1 -- compressor map, 11 speed lines")
    ax.legend(frameon=False, fontsize=5.5, labelcolor=MUTED, ncol=2, title="%NGc", title_fontsize=6)
    out["f1"] = {
        "lines": int(len(m1.params)),
        "knots_per_line": int(m1.lines[0].x.size),
        "source": "Fig. A1, p.56",
    }

    ax = axes.ravel()[11]
    printed = maps.f1_as_printed()
    ax.plot(
        printed.lines[5].x,
        printed.lines[5].y,
        "o--",
        color=BALLIN,
        ms=5,
        lw=1.4,
        zorder=3,
        label="as printed, 7 markers",
    )
    ax.plot(
        m1.lines[5].x, m1.lines[5].y, "-", color=OURS, lw=2, zorder=4, label="on the beta grid, 56"
    )
    _style(ax, fx, fy, "f1: 89 % line, both ways")
    ax.legend(frameon=False, fontsize=7, labelcolor=MUTED, loc="best")

    plotstyle.title(
        fig,
        "Appendix A's eleven engine function tables, as digitized. Crosses are the printed "
        "markers; the line is what the model evaluates",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.958))
    save(fig, "engine-maps")
    return out


def main() -> int:
    """Every sheet, once per theme. The numbers are identical, so only the first set keeps.

    `plotstyle.available()` is reported rather than enforced: without IBM Plex installed
    Matplotlib falls back to DejaVu and the figures are still correct, just not typeset like
    the page they sit in.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    if not plotstyle.available():
        print("NOTE: IBM Plex is not installed; figures will fall back to DejaVu Sans.")
    data: dict = {}
    for theme in plotstyle.THEMES:
        print(f"\n=== {theme.name} ===", flush=True)
        use_theme(theme)
        with plotstyle.context(theme):
            _render(data if theme is plotstyle.THEMES[0] else {})
    (OUT / "report.json").write_text(json.dumps(data, indent=1))
    pngs = sorted(OUT.glob("*.png"))
    print(
        f"\nwrote {OUT}/report.json and {len(pngs)} plots "
        f"({len(pngs) // len(plotstyle.THEMES)} figures x {len(plotstyle.THEMES)} themes)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
