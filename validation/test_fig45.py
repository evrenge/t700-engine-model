"""Figures 4 and 5 [pdf pp.30, 34]: the linear models' response to fuel pulses.

Figure 4 overlays the Gen Hel UH-60A simulation and the 5-DOF, 2-DOF and reduced-order linear
models, "nearly identical" [pdf p.29]. Figure 5 overlays the 6-DOF and 5-DOF, "with and
without station 4.1 heat-sink model" [pdf p.34]. Both captions say only "a high-power flight
condition"; it is the hover trim of Table B.1, because Ballin's own printed hover matrices
(B2, B1) land inside Figure 4's plotted band and on Figure 5's outer curve, and the level and
descent matrices miss its NP panel (#80).

`tools/digitize_fig45.py` reads each speed panel as a band: the uppermost and lowermost line
centres per column. A response is scored by how far it falls outside that band, with dash
gaps filled over `GAP_COLUMNS` and the figure's time read error allowed. In Figure 5, where
both curves are inked and apart, the outer is the 5-DOF and the inner the 6-DOF.

The input is the fuel trace's four plateaus at their printed tick values, from the digitized
step instants. NG and NP are % of the trim value, which is what the axes print.
"""

from __future__ import annotations

import functools
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage
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
GAP_COLUMNS = 21
"""The window that fills a broken line's gaps: on Figure 5 one curve runs alone for up to ten
columns where the other's dash is missing."""
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


def _best(t, m, tb, target, tol, distance) -> np.ndarray:
    """Per column, the smallest `distance` within the time read error."""
    out = np.full(tb.size, np.inf)
    for tau in np.linspace(-tol, tol, 9):
        out = np.minimum(out, distance(np.interp(tb + tau, t, m), target))
    return out


def filled(fig: int, key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """`(t, lower, upper)`: the panel's band with its dash gaps filled over `GAP_COLUMNS`."""
    b = band(fig, key)
    lo = ndimage.minimum_filter1d(b["lower"], GAP_COLUMNS)
    return b["t"], lo, ndimage.maximum_filter1d(b["upper"], GAP_COLUMNS)


def curves_5(key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """`(t, outer, inner)` of Figure 5's two curves, over the columns where both are inked and
    apart: the outer is the 5-DOF, the inner the 6-DOF."""
    b = band(5, key)
    apart = (b["runs"] >= 2) & ((b["upper"] - b["lower"]) > APART_PX * PX_PCT[key])
    above = 0.5 * (b["upper"] + b["lower"]) >= 100.0
    outer = np.where(above, b["upper"], b["lower"])
    inner = np.where(above, b["lower"], b["upper"])
    return b["t"][apart], outer[apart], inner[apart]


def _panel(name: str, fig: int, key: str) -> tuple[np.ndarray, np.ndarray]:
    t, ng, np_ = response(name, fig)
    return t, (ng if key == "ng" else np_)


def outside(name: str, fig: int, key: str) -> np.ndarray:
    """How far, in % of trim, a response falls outside the panel's gap-filled band."""
    tb, lo, up = filled(fig, key)
    t, m = _panel(name, fig, key)
    return _best(
        t,
        m,
        tb,
        (lo, up),
        TIME_TOL_S[fig],
        lambda v, g: np.maximum(0, np.maximum(g[0] - v, v - g[1])),
    )


def to_curve(name: str, key: str, which: str) -> np.ndarray:
    """Distance, % of trim, from a response to Figure 5's outer or inner curve."""
    tb, outer, inner = curves_5(key)
    t, m = _panel(name, 5, key)
    target = outer if which == "outer" else inner
    return _best(t, m, tb, target, TIME_TOL_S[5], lambda v, g: np.abs(v - g))


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
    gap = np.abs(steps(4)[0] - steps(5)[0])
    assert gap.max() < max(TIME_TOL_S.values()), f"step instants differ by {gap} s"


# ------------------------------------------------------------------------- Figure 4's band


FIG4_BAND_RMS_PCT = {"ng": 0.011, "np": 0.027}
"""Ours inside Figure 4's band, rms of the distance outside, % of trim: a characterization,
1.25 times today's worst over our 5-DOF, 2-DOF, reduced-order and nonlinear responses."""
OURS_IN_FIG4 = ("5-DOF", "2-DOF", "reduced", "nonlinear")


@pytest.mark.parametrize("name", ("ballin 5-DOF", "ballin 2-DOF", "ballin reduced"))
@pytest.mark.parametrize("key", PANELS)
def test_ballins_printed_hover_matrices_draw_figure_4(name: str, key: str):
    """The identification (#80): his own Appendix B, at hover, lands in his figure."""
    assert rms(outside(name, 4, key)) < FIG4_BAND_RMS_PCT[key], name


@pytest.mark.parametrize("name", OURS_IN_FIG4)
@pytest.mark.parametrize("key", PANELS)
def test_ours_lands_in_figure_4s_band(name: str, key: str):
    got = rms(outside(name, 4, key))
    assert got < FIG4_BAND_RMS_PCT[key], f"our {name} {key}: {got:.4f} % of trim outside"


FIG4_DECOY_MIN_PCT = 0.05
"""The band test must be able to fail: our 6-DOF, which Figure 4 does not plot, lies outside
it by several times the bound."""


@pytest.mark.parametrize("key", PANELS)
def test_figure_4s_band_rejects_the_heat_sink_model(key: str):
    assert rms(outside("6-DOF", 4, key)) > FIG4_DECOY_MIN_PCT


# ------------------------------------------------------------------------ Figure 5's curves


FIG5_CURVE_RMS_PCT = {"5-DOF": {"ng": 0.008, "np": 0.039}, "6-DOF": {"ng": 0.041, "np": 0.047}}
"""Our 5-DOF against Figure 5's outer curve and our 6-DOF against its inner one, rms over the
columns where both are inked and apart, % of trim: a characterization, 1.25 times today's.
The 6-DOF sits further off: ours is slightly more damped, and Appendix B's own heat-sink time
constants (#31) barely change that."""


@pytest.mark.parametrize("key", PANELS)
def test_ballins_printed_5_dof_is_figure_5s_outer_curve(key: str):
    assert rms(to_curve("ballin 5-DOF", key, "outer")) < FIG5_CURVE_RMS_PCT["5-DOF"][key]


@pytest.mark.parametrize("name,which", (("5-DOF", "outer"), ("6-DOF", "inner")))
@pytest.mark.parametrize("key", PANELS)
def test_ours_lands_on_figure_5s_curves(name: str, which: str, key: str):
    got = rms(to_curve(name, key, which))
    assert got < FIG5_CURVE_RMS_PCT[name][key], f"our {name} {key}: {got:.4f} % of trim"


@pytest.mark.parametrize("key", PANELS)
def test_figure_5s_curves_tell_the_two_models_apart(key: str):
    """Swapped, each model sits at least twice as far from the other's curve."""
    for name, right, wrong in (("5-DOF", "outer", "inner"), ("6-DOF", "inner", "outer")):
        assert rms(to_curve(name, key, wrong)) > 2.0 * rms(to_curve(name, key, right)), name
