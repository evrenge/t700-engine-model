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
APART_PX = 8.0
"""Two curves count as apart when their centres are two line widths apart or more."""


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


def _best(t, m, tb, target, tol) -> np.ndarray:
    """Per column, the smallest |response - target| within the time read error."""
    out = np.full(tb.size, np.inf)
    for tau in np.linspace(-tol, tol, 9):
        out = np.minimum(out, np.abs(np.interp(tb + tau, t, m) - target))
    return out


@functools.cache
def his_curves(fig: int, key: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Ballin's curves in one panel, each as `(t, value)`.

    Figure 4: "all four", the middle of the ink in each column. Figure 5: where both curves are
    inked and apart, the outer is the 5-DOF and the inner the 6-DOF (the heat sink damps it);
    a column with one line inked joins the curve it continues, judged against the columns
    within `LOCAL_S` where both are seen, or both curves where none are, since they coincide
    there."""
    b = band(fig, key)
    t, up, lo = b["t"], b["upper"], b["lower"]
    mid = 0.5 * (up + lo)
    if fig == 4:
        return {"all four": (t, mid)}
    apart = (b["runs"] >= 2) & ((up - lo) > APART_PX * PX_PCT[key])
    above = mid >= 100.0
    outer = np.where(above, up, lo)
    inner = np.where(above, lo, up)
    five, six = np.full(t.size, np.nan), np.full(t.size, np.nan)
    five[apart], six[apart] = outer[apart], inner[apart]
    ta, oa, ia = t[apart], outer[apart], inner[apart]
    for j in np.flatnonzero(~apart):
        near = np.abs(ta - t[j]) <= LOCAL_S
        if not near.any():
            five[j] = six[j] = mid[j]
            continue
        o = np.interp(t[j], ta[near], oa[near])
        i = np.interp(t[j], ta[near], ia[near])
        if abs(o - i) < APART_PX * PX_PCT[key]:
            five[j] = six[j] = mid[j]
        elif abs(mid[j] - o) < abs(mid[j] - i):
            five[j] = mid[j]
        else:
            six[j] = mid[j]
    return {
        "5-DOF": (t[~np.isnan(five)], five[~np.isnan(five)]),
        "6-DOF": (t[~np.isnan(six)], six[~np.isnan(six)]),
    }


def gap(name: str, fig: int, key: str, curve: str) -> np.ndarray:
    """Per column of Ballin's `curve`, how far response `name` is from it, % of trim."""
    tb, v = his_curves(fig, key)[curve]
    t, ng, np_ = response(name, fig)
    return _best(t, ng if key == "ng" else np_, tb, v, TIME_TOL_S[fig])


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
    "5 5-DOF": {"ng": 0.012, "np": 0.036},
    "5 6-DOF": {"ng": 0.041, "np": 0.045},
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


DECOY_MIN_RATIO = 2.0
"""The comparison must be able to fail: the wrong model sits at least this many times further
from a curve than the right one."""


@pytest.mark.parametrize("key", PANELS)
def test_the_comparison_tells_the_models_apart(key: str):
    """Our 6-DOF against Figure 4's line, and each of Figure 5's curves against the other
    model, sit at least `DECOY_MIN_RATIO` times further off than the right pairing."""
    assert rms(gap("6-DOF", 4, key, "all four")) > DECOY_MIN_RATIO * rms(
        gap("5-DOF", 4, key, "all four")
    )
    for ours, right, wrong in (("5-DOF", "5-DOF", "6-DOF"), ("6-DOF", "6-DOF", "5-DOF")):
        assert rms(gap(ours, 5, key, wrong)) > DECOY_MIN_RATIO * rms(gap(ours, 5, key, right))
