"""The real-time transient model: Eqs. 69-80.

Two things are checked here that no steady-state comparison can reach: that the
quasi-steady frame reproduces the trim as a fixed point, and that the opened compressor
iteration is a *consistent discretization* -- that refining the time step converges.
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import constants as c
from t700 import maps, realtime, trim
from t700.engine import Ambient
from t700.units import wf_pps_from_pph
from test_trim_points import CASES

AMB = Ambient(14.696, 518.67)


def _seed(wf_pph: float, lag_whole_flow: bool = True):
    wf = wf_pps_from_pph(wf_pph)
    r = trim.solve(wf, 20895.0, AMB)
    st = realtime.from_trim(r, lag_whole_flow=lag_whole_flow)
    return st, wf, r.frame.q_pt_ftlbf


@pytest.mark.parametrize("wf_pph", [case.wf_pph for case in CASES] + [400.0])
def test_every_published_trim_is_a_fixed_point(wf_pph: float):
    """The trim of the differential model (Eqs. 42-45) is a fixed point of the real-time frame
    (Eqs. 74-80); a wrong algebraic reduction would drift. 400 lbm/hr is Figures 9 and 10's
    starting trim, which must be still for the half second before their step.
    """
    st, wf, qreq = _seed(wf_pph)
    tr = realtime.run(st, lambda t: wf, AMB, duration_s=2.0, q_req_ftlbf=qreq)
    for key in ("ng", "p3", "p41", "p45"):
        drift = abs(tr[key][-1] - tr[key][0]) / abs(tr[key][0])
        assert drift < 1e-9, f"{key} drifted {drift:.2e} over 2 s at a converged trim"
    # and NG to a millionth of an rpm over every frame of the first second
    first_s = tr["t"] < 1.0 + realtime.FRAME_ENGINE_S
    assert np.abs(tr["ng"][first_s] - tr["ng"][0]).max() < 1e-6
    # and still, peak to peak, over the half second a fuel step waits for
    pre = tr["ng"][tr["t"] < 0.5]
    assert pre.max() - pre.min() < 1e-6, "the engine must be still before the step"


def test_both_readings_of_eq_74_converge_to_the_same_limit():
    """Both readings of Eq. 74 [pdf p.36] are first-order discretizations with one limit (#22).

    As printed it lags the bleed, `WA31(n) = WA3(n) - WA3_bl(n-1)`; the prose lags the whole
    compressor flow. The shipped model uses the prose reading because the prose describes
    opening the loop, not because the printed equation fails. Trap: carrying `WA31` instead of
    `WA3_bl` in the printed branch makes an involution that oscillates with period two.

    `tol` is at machine precision on purpose: a refinement study at the report's 0.1 percent
    [pdf p.37] measures the solver's tolerance floor, not the discretization.
    """
    limits = {}
    for lag in (True, False):
        vals = []
        for dt_ms in (8.0, 4.0, 2.0, 1.0, 0.5, 0.25):
            st, wf0, qreq = _seed(400.0, lag_whole_flow=lag)
            tr = realtime.run(
                st,
                lambda t, w=wf0: w if t < 0.5 else wf_pps_from_pph(775.0),
                AMB,
                duration_s=1.0,
                dt=dt_ms / 1000.0,
                q_req_ftlbf=qreq,
                integrate_np=False,
                lag_whole_flow=lag,
                tol=1e-12,
            )
            vals.append(float(np.interp(1.0, tr["t"], tr["ng"])))
        diffs = np.abs(np.diff(vals))
        assert np.all(np.diff(diffs) < 0), (
            f"lag_whole_flow={lag}: successive differences {diffs} are not shrinking, so "
            f"this reading is not converging. If it is the printed one, #22's retraction "
            f"is wrong and the original verdict stands."
        )
        ratios = diffs[:-1] / diffs[1:]
        assert 1.7 < np.median(ratios) < 2.6, (
            f"lag_whole_flow={lag}: error ratios {ratios} are not first order"
        )
        assert np.all(diffs[1:] < 0.75 * diffs[:-1]), (
            f"lag_whole_flow={lag}: halving dt should roughly halve the error for a "
            f"first-order scheme; successive differences {diffs}"
        )
        limits[lag] = vals[-1]

    gap = abs(limits[True] - limits[False]) / limits[True] * 100.0
    assert gap < 0.01, (
        f"the two readings converge to different limits ({limits[True]:.1f} against "
        f"{limits[False]:.1f} rpm, {gap:.4f} %). They are two lags of one loop and should "
        f"agree in the limit; a real gap would mean one of them is not Eq. 74."
    )


def test_a_trimmed_engine_sits_still_in_both_readings_of_eq_74():
    """Neither reading of Eq. 74 changes the equilibrium, so a trimmed engine sits still in
    both. A wrong carried quantity shows here as a period-two oscillation.
    """
    for lag in (True, False):
        st, wf0, qreq = _seed(400.0, lag_whole_flow=lag)
        tr = realtime.run(
            st,
            lambda t, w=wf0: w,
            AMB,
            duration_s=2.0,
            dt=0.007,
            q_req_ftlbf=qreq,
            integrate_np=False,
            lag_whole_flow=lag,
        )
        drift = abs(tr["ng"][-1] / tr["ng"][0] - 1.0) * 100.0
        spread = float(tr["t41"].max() - tr["t41"].min())
        assert drift < 1e-6, f"lag_whole_flow={lag}: NG drifts {drift:.2e} % over 2 s"
        assert spread < 1.0, (
            f"lag_whole_flow={lag}: T41ns spans {spread:.1f} degR at a fixed fuel flow, "
            f"so the carried quantity does not match what the branch expects"
        )


def test_transient_excursions_outside_the_maps_are_counted():
    """An open-loop step leaves the digitized data (FAR past f6's plotted 0.010-0.020), and
    the clamp counter must record it. Open question #45.
    """
    maps.reset_clamps()
    st, wf0, qreq = _seed(400.0)
    tr = realtime.run(
        st,
        lambda t: wf0 if t < 0.5 else wf_pps_from_pph(775.0),
        AMB,
        duration_s=3.0,
        q_req_ftlbf=qreq,
        integrate_np=False,
    )
    report = maps.clamp_report()
    assert report, "this step is known to leave the maps; the counter must record it"
    assert tr["far"].max() > 0.020, "FAR should exceed f6's plotted range on this step"
    maps.reset_clamps()


def test_the_two_to_one_multirate_costs_little_but_is_the_reports():
    """The report's 2:1 multirate [pdf p.47]: engine every 7 ms, NP on the 14 ms rotor frame.

    On an open-loop step with NP free against a constant load (the worst case) it moves NP
    measurably and NG not at all, since Eq. 45 differences only Q_GT and Q_C. With
    `integrate_np=False`, as Figures 9 and 10 run [pdf p.39], it is a bit-identical no-op.
    """
    wf0 = wf_pps_from_pph(400.0)
    r0 = trim.solve(wf0, c.NP_DES, AMB)

    def go(multirate: bool, integrate_np: bool):
        st = realtime.from_trim(r0)
        return realtime.run(
            st,
            lambda t: wf0 if t < 0.5 else wf_pps_from_pph(775.0),
            AMB,
            duration_s=2.0,
            q_req_ftlbf=r0.frame.q_pt_ftlbf,
            heat_sink=True,
            integrate_np=integrate_np,
            multirate=multirate,
        )

    single, multi = go(False, True), go(True, True)
    d_np = float(np.abs(multi["np"] - single["np"]).max())
    assert 50.0 < d_np < 120.0, (
        f"the multirate moves NP by {d_np:.1f} rpm; near zero means the 2:1 is not applied"
    )
    assert np.array_equal(multi["ng"], single["ng"]), (
        "NP must not reach NG: Eq. 45 differences Q_GT and Q_C, neither of which sees NP"
    )

    a, b = go(False, False), go(True, False)
    assert all(np.array_equal(a[k], b[k]) for k in a), (
        "with the NP integration suppressed the multirate must be a bit-identical no-op, "
        "which is what makes every Figure 9 and 10 comparison independent of it"
    )
