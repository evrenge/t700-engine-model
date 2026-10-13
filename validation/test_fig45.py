"""Figures 4 and 5 [pdf pp.30, 34]: the linear models' response to fuel pulses, Ballin's
curves against ours, one for one.

Figure 4 overlays the Gen Hel UH-60A simulation and the 5-DOF, 2-DOF and reduced-order linear
models, "nearly identical" [pdf p.29]: on the page the four lie on one another, so they are
read as one line, the middle of their ink, and each of ours is compared with it. Figure 5
overlays two curves that part, the 5-DOF and the 6-DOF, "with and without station 4.1
heat-sink model" [pdf p.34]: they are separated (`his_curves`), and our 5-DOF is compared with
his 5-DOF, our 6-DOF with his 6-DOF.

Both captions say only "a high-power flight condition"; it is the hover trim of Table B.1,
because Ballin's own printed hover matrices (B2, B1) reproduce his curves and the level and
descent matrices miss Figure 4's NP panel (#80). The input is the fuel trace's four plateaus at
their printed tick values, from the digitized step instants (`tools/digitize_fig45.py`). NG
and NP are % of the trim value, which is what the axes print. Every difference allows the
figure's time read error.
"""

from __future__ import annotations

import functools
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

from t700 import appendix_b as ab
from t700 import constants as c
from t700 import realtime, trim
from t700.linear import DOF
from t700.units import wf_pps_from_pph
from test_all_dof_models import AMB, DQ_REQ_DNP, NP_RPM, WF_PPH, _model
from test_fuel_step import read_reference

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
from digitize_fig45 import INPUT_LEVELS  # noqa: E402

FIGS = (4, 5)
PANELS = ("ng", "np")
TIME_TOL_S = {4: 0.0314, 5: 0.0444}
"""Each figure's time read error: its worst interior major tick, as the digitizer reports it
in every CSV header."""
LOCAL_S = 0.15
"""How far, in time, a column with only one curve inked looks for columns where both are, to
tell which curve it belongs to: Figure 5's dashes leave one curve alone for up to ten columns,
0.1 s."""
PX_PCT = {"ng": 2.0 / 707.0, "np": 3.0 / 707.0}
"""% of trim per pixel: each panel's printed range over its height, 707 px on both pages."""
GUIDE_PX = 6.0
"""Figure 5's two curves count as one where they come within this many pixels, one and a half
line widths."""
CLUSTER_PX = 4.0
"""Runs in one column closer than a line width are one stroke: a steep dash-dot line breaks
into stacked fragments."""
VOTE_S = 0.3
"""Where tracking starts, which stroke is the 5-DOF is a vote over the columns this close."""
TRACK_POINTS = 8
"""Each curve is followed by a line through its last eight points, about 0.08 s."""
RESET_S = 0.3
"""A run of columns without both curves longer than this restarts the tracking: the curves
may have crossed unseen."""


@functools.cache
def band(fig: int, key: str) -> dict[str, np.ndarray]:
    t, up = read_reference(f"fig{fig:02d}_{key}.csv", "t_s", "upper")
    _, lo = read_reference(f"fig{fig:02d}_{key}.csv", "t_s", "lower")
    _, runs = read_reference(f"fig{fig:02d}_{key}.csv", "t_s", "runs")
    return {"t": t, "upper": up, "lower": lo, "runs": runs}


@functools.cache
def steps(fig: int) -> tuple[np.ndarray, np.ndarray]:
    """`(step instants, digitized plateau levels)` of the fuel trace, after the trim."""
    return read_reference(f"fig{fig:02d}_input.csv", "t_start_s", "level_pct")


def fuel(fig: int):
    """The input, % of trim / 100: trim until the first step, then `INPUT_LEVELS`."""
    starts = steps(fig)[0]

    def u(t: float) -> float:
        k = int(np.searchsorted(starts, t, side="right"))
        return 1.0 if k == 0 else INPUT_LEVELS[k - 1] / 100.0

    return u


@functools.cache
def _hover():
    wf = wf_pps_from_pph(WF_PPH[1])
    return wf, trim.solve(wf, NP_RPM, AMB)


def linear_response(A: np.ndarray, b: np.ndarray, fig: int, dt: float = 0.001):
    """`(t, NG %, NP %)` of `x-dot = A x + b dWf` from rest at hover, zero-order hold."""
    wf0, r = _hover()
    n = A.shape[0]
    ad = expm(A * dt)
    bd = np.linalg.solve(A, (ad - np.eye(n)) @ b)
    u = fuel(fig)
    t = np.arange(0.0, 15.0 + dt / 2, dt)
    x, out = np.zeros(n), np.zeros((t.size, 2))
    for i, ti in enumerate(t):
        out[i] = x[:2]
        x = ad @ x + bd * (u(ti) - 1.0) * wf0
    return t, 100.0 * (1.0 + out[:, 0] / r.state.ng_rpm), 100.0 * (1.0 + out[:, 1] / NP_RPM)


def reduced(A: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The 5-DOF with its pressure derivatives set to zero [pdf p.28]."""
    a11, a12, a21, a22 = A[:2, :2], A[:2, 2:], A[2:, :2], A[2:, 2:]
    return a11 - a12 @ np.linalg.solve(a22, a21), b[:2] - a12 @ np.linalg.solve(a22, b[2:])


@functools.cache
def response(name: str, fig: int):
    """One named response to `fig`'s input. Ballin's are his printed hover matrices."""
    if name == "ballin 5-DOF":
        m = ab.find(5, 1)
        return linear_response(m.A, m.b, fig)
    if name == "ballin 2-DOF":
        m = ab.find(2, 1)
        return linear_response(m.A, m.b, fig)
    if name == "ballin reduced":
        return linear_response(*reduced(ab.find(5, 1).A, ab.find(5, 1).b), fig)
    if name == "nonlinear":
        return nonlinear(fig)
    dof = {"5-DOF": DOF.FIVE, "2-DOF": DOF.TWO, "reduced": DOF.REDUCED_FIVE, "6-DOF": DOF.SIX}
    m = _model(dof[name], 1)
    return linear_response(m.A, m.b, fig)


def nonlinear(fig: int, dt: float = realtime.FRAME_ENGINE_S):
    """Our real-time engine at hover with Appendix B's rotor load, heat sink off: Figure 4's
    Gen Hel trace sits on the 5-DOF curve, which has none [pdf p.27]."""
    wf0, r = _hover()
    st = realtime.from_trim(r)
    q0, dq = r.frame.q_pt_ftlbf, DQ_REQ_DNP[1]
    u = fuel(fig)
    n = int(round(15.0 / dt)) + 1
    t = np.arange(n) * dt
    ng, np_ = np.zeros(n), np.zeros(n)
    for i in range(n):
        ng[i], np_[i] = st.ng_rpm, st.np_rpm
        st, _ = realtime.step(
            st,
            u(t[i]) * wf0,
            AMB,
            dt,
            q_req_ftlbf=q0 + dq * (st.np_rpm - NP_RPM),
            j_load=c.J_LOAD_UH60A,
            integrate_np=(i % 2 == 0),  # NP every 14 ms, as in Gen Hel [pdf p.47]
            dt_np=2 * dt,
        )
    return t, 100.0 * ng / r.state.ng_rpm, 100.0 * np_ / NP_RPM


def _best(t, m, tb, target, half, tol) -> np.ndarray:
    """Per ink point, how far the response passes from the inked stroke `target +- half`,
    the smallest within the time read error: zero where it crosses the stroke."""
    out = np.full(tb.size, np.inf)
    for tau in np.linspace(-tol, tol, 9):
        v = np.interp(tb + tau, t, m)
        out = np.minimum(out, np.maximum(0.0, np.abs(v - target) - half))
    return out


@functools.cache
def runs(key: str) -> dict[str, np.ndarray]:
    """Every ink run of Figure 5's `key` panel: its scan `column`, `t`, `centre` and `half`
    span, % of trim."""
    name = f"fig05_{key}_runs.csv"
    col, t = read_reference(name, "column", "t_s")
    _, centre = read_reference(name, "column", "centre")
    _, half = read_reference(name, "column", "half_span")
    return {"column": col, "t": t, "centre": centre, "half": half}


def _predict(hist: list[tuple[float, float]], t: float) -> float:
    """A curve's value at `t`, extrapolated by a line through its last `TRACK_POINTS` points."""
    pts = np.array(hist[-TRACK_POINTS:])
    if len(pts) < 3 or np.ptp(pts[:, 0]) == 0.0:
        return float(pts[-1, 1])
    a, b = np.polyfit(pts[:, 0], pts[:, 1], 1)
    return float(a * t + b)


def _clusters(c: np.ndarray, h: np.ndarray, px: float) -> list[tuple[float, float]]:
    """One column's runs merged where they touch, as `(centre, half span)`, bottom first."""
    o = np.argsort(c)
    out: list[list[float]] = []
    for k in o:
        a, b = c[k] - h[k], c[k] + h[k]
        if out and a - out[-1][1] <= CLUSTER_PX * px:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(0.5 * (a + b), 0.5 * (b - a)) for a, b in out]


@functools.cache
def his_curves(fig: int, key: str) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Ballin's curves in one panel, each as `(t, value, half span)`.

    Figure 4: "all four", the middle of the ink in each column. Figure 5: a column's runs are
    merged into strokes (`_clusters`); where it holds two strokes apart, the one nearer
    Ballin's own printed 5-DOF (B2) response is his 5-DOF, by a vote over the columns within
    `VOTE_S` weighted by the square of their separation, so the columns where the curves are
    far apart, and B2's few-hundredths-of-a-percent offset from his ink cannot matter, decide
    the close ones until the curves cross. A column with one stroke joins the nearer of the two
    curves around it, or both where those are within `GUIDE_PX`, which is where the curves
    coincide. The values are the ink's; B2 only says which line is which, and it is his
    matrix, not ours."""
    if fig == 4:
        b = band(4, key)
        return {"all four": (b["t"], 0.5 * (b["upper"] + b["lower"]), np.zeros(b["t"].size))}
    r = runs(key)
    col, t, c, h = r["column"], r["t"], r["centre"], r["half"]
    px = PX_PCT[key]
    tg, *gs = response("ballin 5-DOF", 5)
    g = gs[0] if key == "ng" else gs[1]
    pairs, singles = [], []  # (t, lower, upper) strokes; (t, stroke)
    for idx in np.split(np.arange(t.size), np.flatnonzero(np.diff(col)) + 1):
        tt = float(np.mean(t[idx]))
        cl = _clusters(c[idx], h[idx], px)
        if len(cl) >= 2 and (cl[-1][0] - cl[-1][1]) - (cl[0][0] + cl[0][1]) > GUIDE_PX * px:
            pairs.append((tt, cl[0], cl[-1]))
        else:
            singles.extend((tt, q) for q in cl)
    pt = np.array([q[0] for q in pairs])
    lo_c = np.array([q[1][0] for q in pairs])
    up_c = np.array([q[2][0] for q in pairs])
    sep = np.array([(q[2][0] - q[2][1]) - (q[1][0] + q[1][1]) for q in pairs])
    b2 = np.interp(pt, tg, g)
    vote = np.where(np.abs(up_c - b2) < np.abs(lo_c - b2), 1.0, -1.0)
    upper_is_five = np.zeros(pt.size, bool)
    for seg in np.split(np.arange(pt.size), np.flatnonzero(np.diff(pt) > RESET_S) + 1):
        k0 = int(seg[np.argmax(sep[seg])])
        w = seg[np.abs(pt[seg] - pt[k0]) <= VOTE_S]
        upper_is_five[k0] = float(np.sum(vote[w] * sep[w] ** 2)) >= 0.0
        for order in (seg[seg > k0], seg[seg < k0][::-1]):
            h5 = [(pt[k0], up_c[k0] if upper_is_five[k0] else lo_c[k0])]
            h6 = [(pt[k0], lo_c[k0] if upper_is_five[k0] else up_c[k0])]
            for k in order:
                p5, p6 = _predict(h5, pt[k]), _predict(h6, pt[k])
                upper_is_five[k] = abs(up_c[k] - p5) + abs(lo_c[k] - p6) <= abs(lo_c[k] - p5) + abs(
                    up_c[k] - p6
                )
                h5.append((pt[k], up_c[k] if upper_is_five[k] else lo_c[k]))
                h6.append((pt[k], lo_c[k] if upper_is_five[k] else up_c[k]))
    rows5, rows6 = [], []
    for k, (tt, lo, up) in enumerate(pairs):
        rows5.append((tt, *(up if upper_is_five[k] else lo)))
        rows6.append((tt, *(lo if upper_is_five[k] else up)))
    f5, f6 = np.array(rows5), np.array(rows6)
    for tt, (cc, hh) in singles:
        near = np.abs(f5[:, 0] - tt) <= LOCAL_S
        if not near.any():
            rows5.append((tt, cc, hh))
            rows6.append((tt, cc, hh))
            continue
        a5 = float(np.interp(tt, f5[near, 0], f5[near, 1]))
        a6 = float(np.interp(tt, f6[near, 0], f6[near, 1]))
        if abs(a5 - a6) <= GUIDE_PX * px or abs(abs(cc - a5) - abs(cc - a6)) <= hh:
            rows5.append((tt, cc, hh))
            rows6.append((tt, cc, hh))
        elif abs(cc - a5) < abs(cc - a6):
            rows5.append((tt, cc, hh))
        else:
            rows6.append((tt, cc, hh))
    out = {}
    for name, rows in (("5-DOF", rows5), ("6-DOF", rows6)):
        a = np.array(sorted(rows))
        out[name] = (a[:, 0], a[:, 1], a[:, 2])
    return out


def gap(name: str, fig: int, key: str, curve: str) -> np.ndarray:
    """Per column of Ballin's `curve`, how far response `name` is from it, % of trim."""
    tb, v, half = his_curves(fig, key)[curve]
    t, ng, np_ = response(name, fig)
    return _best(t, ng if key == "ng" else np_, tb, v, half, TIME_TOL_S[fig])


def rms(v: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(v))))


# ------------------------------------------------------------------------------- the input


INPUT_LEVEL_TOL_PCT = 0.6
"""A digitized plateau against its printed tick value, % of the 15 % panel height: the
ceiling the digitizer itself enforces, Figures 9-10's read error."""


@pytest.mark.parametrize("fig", FIGS)
def test_the_input_sits_on_printed_tick_values(fig: int):
    levels = steps(fig)[1]
    miss = 100.0 * (levels - np.array(INPUT_LEVELS)) / 15.0
    assert np.max(np.abs(miss)) < INPUT_LEVEL_TOL_PCT, f"Figure {fig}: {miss} % of panel height"


def test_the_two_figures_share_one_input():
    """The same pulses, drawn twice: their step instants agree within the read error."""
    diff = np.abs(steps(4)[0] - steps(5)[0])
    assert diff.max() < max(TIME_TOL_S.values()), f"step instants differ by {diff} s"


# --------------------------------------------------------------------- the pairs, one by one


PAIRS = (
    (4, "all four", "5-DOF"),
    (4, "all four", "2-DOF"),
    (4, "all four", "reduced"),
    (4, "all four", "nonlinear"),
    (5, "5-DOF", "5-DOF"),
    (5, "6-DOF", "6-DOF"),
)
"""Ballin's curve against ours: Figure 4's one line against each of our four, Figure 5's two
curves against our two."""

CURVE_RMS_PCT = {
    "4 nonlinear": {"ng": 0.017, "np": 0.040},
    "4 linear": {"ng": 0.029, "np": 0.055},
    "5 5-DOF": {"ng": 0.007, "np": 0.035},
    "5 6-DOF": {"ng": 0.038, "np": 0.042},
}
"""Ours against Ballin's curve, rms difference, % of trim: characterizations, 1.25 times
today's worst in each group (Figure 4's three linear models share a bound)."""


def _bound(fig: int, ours: str, key: str) -> float:
    group = "nonlinear" if ours == "nonlinear" else ("linear" if fig == 4 else ours)
    return CURVE_RMS_PCT[f"{fig} {group}"][key]


@pytest.mark.parametrize("fig,his,ours", PAIRS)
@pytest.mark.parametrize("key", PANELS)
def test_ours_follows_ballins_curve(fig: int, his: str, ours: str, key: str):
    got = rms(gap(ours, fig, key, his))
    assert got < _bound(fig, ours, key), f"Figure {fig} {his} vs our {ours} {key}: {got:.4f} %"


@pytest.mark.parametrize(
    "fig,his,name",
    (
        (4, "all four", "ballin 5-DOF"),
        (4, "all four", "ballin 2-DOF"),
        (4, "all four", "ballin reduced"),
        (5, "5-DOF", "ballin 5-DOF"),
    ),
)
@pytest.mark.parametrize("key", PANELS)
def test_ballins_printed_hover_matrices_draw_his_curves(fig: int, his: str, name: str, key: str):
    """The identification (#80): his own Appendix B, at hover, lands on his curves."""
    assert rms(gap(name, fig, key, his)) < _bound(fig, "5-DOF", key), name


DECOY_MIN_RATIO = 1.5
"""The comparison must be able to fail: the wrong model sits at least this many times further
from a curve than the right one. The tightest case is Figure 5's NG 6-DOF curve, which lies
between our two models: our 6-DOF is about 1.9 times nearer to it than our 5-DOF is."""


@pytest.mark.parametrize("key", PANELS)
def test_the_comparison_tells_the_models_apart(key: str):
    """Our 6-DOF against Figure 4's line, and each of Figure 5's curves against the other
    model, sit at least `DECOY_MIN_RATIO` times further off than the right pairing."""
    assert rms(gap("6-DOF", 4, key, "all four")) > DECOY_MIN_RATIO * rms(
        gap("5-DOF", 4, key, "all four")
    )
    for ours, right, wrong in (("5-DOF", "5-DOF", "6-DOF"), ("6-DOF", "6-DOF", "5-DOF")):
        assert rms(gap(ours, 5, key, wrong)) > DECOY_MIN_RATIO * rms(gap(ours, 5, key, right))
