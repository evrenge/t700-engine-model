"""Render the validation figure set: is the model sound, and does the data look right?

Six sheets into `validation/out/plots/`. Each answers a question that a table of
percentages cannot:

    1  function-tables.png   do the eleven digitized maps look like physical functions?
    2  trim-sweep.png        does the engine track Ballin across the operating range?
    3  ps3-vs-ng.png         is the compressor side right, independently of fuel?
    4  transient-fig9.png    does the step response follow Ballin's own trace?
    5  eigenvalues.png       do the dynamics match Table 1's printed modes?
    6  residuals.png         where does the remaining disagreement sit?

Lives in `validation/`, so Matplotlib is allowed here and nowhere near `src/t700/`.

Colour follows the categorical slots in fixed order -- blue for our model, orange for
Ballin, aqua for the GE reference -- and the assignment follows the *entity*, so a series
keeps its colour on every sheet.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from t700 import constants as c  # noqa: E402
from t700 import maps, realtime, trim  # noqa: E402
from t700.engine import Ambient, frame  # noqa: E402
from t700.units import wf_pps_from_pph  # noqa: E402

OUT = Path(__file__).resolve().parent / "out" / "plots"
REF = Path(__file__).resolve().parent.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)

OURS = "#2a78d6"  # categorical slot 1 -- our model
BALLIN = "#eb6834"  # slot 2 -- the report's own model
GE = "#1baf7a"  # slot 3 -- the GE reference models
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e3e2df"

TRIMS = [  # Table B.1, pdf p.67
    ("hover", 476.3, 41638.0, 911.1),
    ("level 80 kt", 349.3, 39768.0, 552.6),
    ("descent 80 kt", 267.7, 38072.0, 302.6),
]


def _style(ax, xlabel: str, ylabel: str, title: str = "") -> None:
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_xlabel(xlabel, color=MUTED, fontsize=9)
    ax.set_ylabel(ylabel, color=MUTED, fontsize=9)
    if title:
        ax.set_title(title, color=INK, fontsize=10, loc="left", pad=8)


def _load(name: str, xk: str, yk: str):
    path = REF / name
    if not path.exists():
        return None, None
    rows = list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))
    x = np.array([float(r[xk]) for r in rows])
    y = np.array([float(r[yk]) for r in rows])
    o = np.argsort(x)
    return x[o], y[o]


def _split_strays(x, y):
    """Separate a digitized reference trace into curve and off-curve samples.

    Returns `(xc, yc, xs, ys)` -- the curve, then the strays. **Nothing is discarded
    silently and nothing is smoothed**: the strays are returned so they can be drawn,
    because they are a measured property of our own digitization (open question #51:
    61 off-curve samples in 6,632) and hiding them would misrepresent how well the
    reference is known.

    Two rules, both conservative:

    * a sample flagged by the digitizer's own `# defect:` header is not data -- the tool
      says so in the file. Five traces carry one: a trailing sample re-acquired hundreds
      of median intervals after its predecessor.
    * a sample more than six pixel quanta from the median of its four neighbours, *and*
      in a locally flat neighbourhood. The flatness guard is what keeps the genuine
      near-vertical step edge -- where a large jump is the data -- out of the count.
    """
    if x is None or len(x) < 6:
        return x, y, np.zeros(0), np.zeros(0)
    dy = np.abs(np.diff(y))
    nz = dy[dy > 0]
    quantum = float(np.median(nz)) if nz.size else 1.0
    # A 9-sample window (four neighbours each side) rather than five: the strays come in
    # runs, and the longest observed is four contiguous samples -- the plunge at the end
    # of Figure 9's PCNG trace, where four samples sit at 90.3 while the curve is at 98.8.
    # A five-sample window cannot outvote a run that long.
    half = 4
    bad = np.zeros(len(x), dtype=bool)
    for i in range(half, len(y) - half):
        nb = np.delete(y[i - half : i + half + 1], half)
        if np.ptp(nb) > 3 * quantum:
            continue
        if abs(y[i] - np.median(nb)) > 6 * quantum:
            bad[i] = True
    dt = np.diff(x)
    if dt.size and dt[-1] > 10 * np.median(dt):
        bad[-1] = True

    # A trailing BLOCK, not just a trailing sample. The digitizer's `# defect:` rule
    # catches a single re-acquired point; Figure 9's PCNG trace ends with four of them,
    # at 90.33 where the settled plateau is 98.73. The rule that covers both: once a
    # trace has settled, a contiguous run reaching the end of the record and sitting
    # more than six quanta off the plateau is re-acquisition, because a settled trace
    # cannot leave its plateau at the very end and simply stop there.
    # The guard matters: the rule presumes the trace has SETTLED. Figure 10's PCNG is
    # still decaying at the end of the record, and without this check the rule found no
    # plateau and condemned 92 perfectly good samples.
    late = x > 0.7 * x.max()
    keep = late & ~bad
    # Interquartile range, not peak-to-peak: the strays are exactly what we are trying
    # to detect, so letting them into the flatness measure blocks their own detection.
    spread = np.subtract(*np.percentile(y[keep], [75, 25])) if keep.sum() else 0.0
    if keep.sum() > 10 and abs(spread) < 10 * quantum:
        plateau = float(np.median(y[keep]))
        k = len(y) - 1
        while k >= 0 and abs(y[k] - plateau) > 6 * quantum:
            bad[k] = True
            k -= 1

    return x[~bad], y[~bad], x[bad], y[bad]


def sheet_function_tables() -> None:
    """Every digitized map, with its knots. A physical function should look like one."""
    names = ["f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "f10", "f_hs"]
    fig, axes = plt.subplots(3, 4, figsize=(15, 9))
    fig.patch.set_facecolor("white")
    for ax, nm in zip(axes.flat, names, strict=False):
        cur = getattr(maps, nm)()
        xs = np.linspace(cur.x[0], cur.x[-1], 400)
        ax.plot(xs, cur(xs), color=OURS, lw=2, zorder=2)
        ax.plot(cur.x, cur.y, "o", color=OURS, ms=5, mec="white", mew=1, zorder=3)
        _style(ax, "", "", f"{nm}   ({cur.x.size} knots)")

    m = maps.f1()
    ax = axes.flat[10]
    for line, p in zip(m.lines, m.params, strict=True):
        ax.plot(line.x, line.y, "-o", color=OURS, lw=1.2, ms=3, alpha=0.75)
        ax.annotate(f"{p:.0f}", (line.x[-1], line.y[-1]), fontsize=6, color=MUTED)
    _style(ax, "Ps3/P2", "WA2c, lbm/s", "f1  (11 speed lines x 7)")
    axes.flat[11].axis("off")
    axes.flat[11].text(
        0.0,
        0.5,
        "Every map read from the printed\nfigure, never from the text layer.\n\n"
        "Flat-then-falling on f1 is choked\nflow; f8 going negative above\n0.83 is in the report.",
        fontsize=9,
        color=MUTED,
        va="center",
    )
    fig.suptitle(
        "The eleven function tables, as digitized", color=INK, fontsize=13, x=0.09, ha="left"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT / "function-tables.png", dpi=130)
    plt.close(fig)


def sheet_trim_sweep() -> None:
    """Our engine across the whole range, against both of Ballin's datasets."""
    wfs = np.arange(150.0, 815.0, 10.0)
    res = trim.sweep(wfs)
    ok = [(w, r) for w, r in zip(wfs, res, strict=True) if r.trustworthy]
    mw = np.array([a for a, _ in ok])
    mng = np.array([100 * r.state.ng_rpm / c.NG_DES for _, r in ok])
    mshp = np.array([r.shp for _, r in ok])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.patch.set_facecolor("white")

    ax = axes[0]
    ax.plot(mw, mng, color=OURS, lw=2, label="our model", zorder=3)
    x, y = _load("fig06_realtime.csv", "wf_pph", "ng_pct")
    if x is not None:
        ax.plot(x, y, "s", color=BALLIN, ms=5, mec="white", mew=1, label="Figure 6", zorder=4)
    tb = np.array([(wf, 100 * ng / c.NG_DES) for _, wf, ng, _ in TRIMS])
    ax.plot(
        tb[:, 0], tb[:, 1], "D", color=INK, ms=9, mec="white", mew=1.5, label="Table B.1", zorder=5
    )
    _style(
        ax, "fuel flow, lbm/hr", "gas generator speed, %", "Speed: we track Table B.1, not Figure 6"
    )
    ax.legend(frameon=False, fontsize=9, labelcolor=MUTED)

    ax = axes[1]
    ax.plot(mw, mshp, color=OURS, lw=2, label="our model", zorder=3)
    x, y = _load("fig07_realtime.csv", "wf_pph", "shp")
    if x is not None:
        ax.plot(x, y, "s", color=BALLIN, ms=5, mec="white", mew=1, label="Figure 7", zorder=4)
    tb = np.array([(wf, shp) for _, wf, _, shp in TRIMS])
    ax.plot(
        tb[:, 0], tb[:, 1], "D", color=INK, ms=9, mec="white", mew=1.5, label="Table B.1", zorder=5
    )
    _style(
        ax,
        "fuel flow, lbm/hr",
        "shaft horsepower",
        "Power: Ballin's two sources differ most at low power",
    )
    ax.legend(frameon=False, fontsize=9, labelcolor=MUTED)

    fig.tight_layout()
    fig.savefig(OUT / "trim-sweep.png", dpi=130)
    plt.close(fig)


def sheet_ps3() -> None:
    wfs = np.arange(150.0, 815.0, 10.0)
    res = trim.sweep(wfs)
    ok = [r for r in res if r.trustworthy]
    mng = np.array([100 * r.state.ng_rpm / c.NG_DES for r in ok])
    mps3 = np.array([c.K_PS3 * r.state.p3_psia for r in ok])

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    fig.patch.set_facecolor("white")
    ax.plot(mng, mps3, color=OURS, lw=2, label="our model", zorder=3)
    x, y = _load("fig08_realtime.csv", "ng_pct", "ps3_psia")
    if x is not None:
        ax.plot(x, y, "s", color=BALLIN, ms=6, mec="white", mew=1, label="Figure 8", zorder=4)
    _style(
        ax,
        "gas generator speed, %",
        "station 3 static pressure, psia",
        "Compressor side: agrees to -0.43 %, sd 0.70",
    )
    ax.legend(frameon=False, fontsize=9, labelcolor=MUTED)
    fig.tight_layout()
    fig.savefig(OUT / "ps3-vs-ng.png", dpi=130)
    plt.close(fig)


def sheet_transient(fig_no: int = 9, wf_hi: float = 775.0, t_step: float = 0.539) -> None:
    """One fuel-step figure: our real-time model, Ballin's own trace, and the GE reference.

    Three series, and only one of them is the target. **Ballin's trace is what we must
    match**; the GE markers are the hardware reference he was himself validating against,
    plotted because they show how large a disagreement is ordinary between two models of
    this engine -- never as a thing to move toward. See `test_ge_reference.py`.

    Run with the **heat sink on**, which these figures were [docs/notes/heat-sink-
    configuration.md]. Until 2026-09-12 this sheet omitted the flag and so plotted the
    5-DOF model against a 6-DOF figure -- the same category error open question #47
    records, left in the plotting script after being fixed in the tests.
    """
    wf0 = wf_pps_from_pph(400.0)
    r0 = trim.solve(wf0, c.NP_DES, AMB)
    f0 = frame(r0.state, wf0, AMB)
    hi = wf_pps_from_pph(wf_hi)
    st = realtime.from_trim(r0, f0.wa31_pps, f0)
    tr = realtime.run(
        st,
        lambda t: wf0 if t < t_step else hi,
        AMB,
        duration_s=5.0,
        dt=0.007,
        q_req_ftlbf=f0.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=True,
    )

    panels = [
        ("pcng", 100 * tr["ng"] / c.NG_DES, "gas generator speed, %"),
        ("ps3", c.K_PS3 * tr["p3"], "station 3 static pressure, psia"),
        ("t41", tr["t41"], "T41, deg R"),
        ("t45", tr["t45"], "T45, deg R"),
        ("torq45", tr["q_pt"], "output torque, ft.lbf"),
    ]
    fig, axes = plt.subplots(1, 5, figsize=(19, 4.2))
    fig.patch.set_facecolor("white")
    for ax, (key, ours, lab) in zip(axes, panels, strict=True):
        x, y = _load(f"fig{fig_no:02d}_{key}_model.csv", "t_s", "value")
        xc, yc, xs, ys = _split_strays(x, y)
        if xc is not None:
            ax.plot(xc, yc, color=BALLIN, lw=2.5, alpha=0.85, label="Ballin", zorder=2)
        if len(xs):
            ax.plot(
                xs,
                ys,
                "o",
                mfc="none",
                mec=BALLIN,
                ms=5,
                mew=1.0,
                alpha=0.5,
                label=f"off-curve ({len(xs)})",
                zorder=5,
            )
        xg, yg = _load(f"fig{fig_no:02d}_{key}_reference.csv", "t_s", "value")
        if xg is not None:
            ax.plot(xg, yg, "+", color=GE, ms=7, mew=1.4, label="GE reference", zorder=3)
        ax.plot(tr["t"], ours, color=OURS, lw=2, label="our model", zorder=4)
        _style(ax, "time, s", lab, key.upper())
        ax.set_xlim(0, 5)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=MUTED, loc="best")
    direction = "increase" if wf_hi > 400.0 else "decrease"
    fig.suptitle(
        f"Figure {fig_no} - fuel step {direction} 400 to {wf_hi:.0f} lbm/hr, "
        f"heat sink on, NP integration suppressed",
        color=INK,
        fontsize=13,
        x=0.055,
        ha="left",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / f"transient-fig{fig_no}.png", dpi=130)
    plt.close(fig)


def sheet_eigenvalues() -> None:
    """Our Jacobian's modes against Table 1's printed values [pdf p.31]."""
    from t700.engine import State

    wf = wf_pps_from_pph(476.3)
    r = trim.solve(wf, 20895.0, AMB)
    x0 = r.state.as_array()
    qreq = r.frame.q_pt_ftlbf

    def deriv(x):
        f = frame(State.from_array(x), wf, AMB, q_req_ftlbf=qreq)
        return np.array([f.dng_dt, f.dnp_dt, f.dp3_dt, f.dp41_dt, f.dp45_dt])

    f0 = deriv(x0)
    A = np.zeros((5, 5))
    for i in range(5):
        h = 1e-6 * max(abs(x0[i]), 1.0)
        xp = x0.copy()
        xp[i] += h
        A[:, i] = (deriv(xp) - f0) / h
    ours = np.sort(np.linalg.eigvals(A).real)[::-1]
    printed = np.array([-0.565, -2.66, -51.6, -3060.0, -4900.0])

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor("white")
    idx = np.arange(5)
    ax.barh(idx + 0.19, -printed, height=0.34, color=BALLIN, label="Table 1", zorder=3)
    ax.barh(idx - 0.19, -ours, height=0.34, color=OURS, label="our model", zorder=3)
    ax.set_xscale("log")
    ax.set_yticks(idx)
    ax.set_yticklabels(["mode 5 (load)", "mode 4", "mode 3", "mode 2", "mode 1"])
    for i, (a, b) in enumerate(zip(ours, printed, strict=True)):
        if abs(b) > 1:
            ax.text(
                max(-a, -b) * 1.15,
                i,
                f"{abs(a - b) / abs(b) * 100:+.1f} %",
                va="center",
                fontsize=8,
                color=MUTED,
            )
    _style(
        ax,
        "|eigenvalue|, 1/sec  (log)",
        "",
        "Dynamics against Table 1.  Mode 5 depends on dQreq/dNP, from Gen Hel",
    )
    ax.legend(frameon=False, fontsize=9, labelcolor=MUTED, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "eigenvalues.png", dpi=130)
    plt.close(fig)


def sheet_residuals() -> None:
    """Where the remaining disagreement sits, against the report's own spread."""
    fig, ax = plt.subplots(figsize=(8.5, 5))
    fig.patch.set_facecolor("white")
    x6, y6 = _load("fig06_realtime.csv", "wf_pph", "ng_pct")
    x7, y7 = _load("fig07_realtime.csv", "wf_pph", "shp")

    labels, ours_dev, spread = [], [], []
    for name, wf, _ng, shp in TRIMS:
        r = trim.solve(wf_pps_from_pph(wf), 20895.0, AMB)
        ours_dev.append((r.shp - shp) / shp * 100.0)
        f7 = float(np.interp(wf, x7, y7)) if x7 is not None else shp
        spread.append(abs(shp - f7) / shp * 100.0)
        labels.append(f"{name}\n{wf:.0f} lbm/hr")

    idx = np.arange(len(labels))
    ax.bar(
        idx,
        spread,
        width=0.62,
        color=BALLIN,
        alpha=0.35,
        zorder=2,
        label="Ballin's own spread (Table B.1 vs Figure 7)",
    )
    ax.bar(
        idx, np.abs(ours_dev), width=0.3, color=OURS, zorder=3, label="our deviation from Table B.1"
    )
    for i, (a, b) in enumerate(zip(ours_dev, spread, strict=True)):
        ax.text(i, abs(a) + 0.12, f"{a:+.2f} %", ha="center", fontsize=8, color=OURS)
        ax.text(i, b + 0.12, f"{b:.2f} %", ha="center", fontsize=8, color=MUTED)
    ax.set_xticks(idx)
    ax.set_xticklabels(labels, fontsize=9, color=MUTED)
    _style(
        ax,
        "",
        "shaft power, % of Table B.1",
        "Our error sits inside the report's own internal disagreement",
    )
    ax.legend(frameon=False, fontsize=9, labelcolor=MUTED)
    fig.tight_layout()
    fig.savefig(OUT / "residuals.png", dpi=130)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    # Both transient figures. Figure 10 is the step DOWN to below-idle power and was
    # long the worst disagreement in the project; it is plotted so that stays visible.
    for fn, args in (
        (sheet_function_tables, ()),
        (sheet_trim_sweep, ()),
        (sheet_ps3, ()),
        (sheet_transient, (9, 775.0, 0.539)),
        (sheet_transient, (10, 125.0, 0.545)),
        (sheet_eigenvalues, ()),
        (sheet_residuals, ()),
    ):
        fn(*args)
        print(f"  {fn.__name__}{args if args else ''}")
    print(f"\n  written to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
