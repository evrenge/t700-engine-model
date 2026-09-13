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
from t700.engine import Ambient, frame
from t700.units import wf_pps_from_pph

AMB = Ambient(14.696, 518.67)


def _seed(wf_pph: float, lag_whole_flow: bool = True):
    wf = wf_pps_from_pph(wf_pph)
    r = trim.solve(wf, 20895.0, AMB)
    f = frame(r.state, wf, AMB)
    st = realtime.from_trim(r, f.wa31_pps, f, lag_whole_flow=lag_whole_flow)
    return st, wf, f.q_pt_ftlbf


def test_a_trimmed_engine_sits_still():
    """The strongest available check that the two formulations agree.

    `trim.solve` finds the equilibrium of the differential model (Eqs. 42-45); the
    real-time frame solves the same physics with the pressures made algebraic
    (Eqs. 74-80). If the algebraic reduction were wrong, a trimmed engine would drift.
    """
    st, wf, qreq = _seed(476.3)
    tr = realtime.run(st, lambda t: wf, AMB, duration_s=2.0, q_req_ftlbf=qreq)
    for key in ("ng", "p3", "p41", "p45"):
        drift = abs(tr[key][-1] - tr[key][0]) / abs(tr[key][0])
        assert drift < 1e-9, f"{key} drifted {drift:.2e} over 2 s at a converged trim"


def test_the_opened_iteration_converges_as_the_step_shrinks():
    """Eq. 74 read as the prose states it must be first-order convergent.

    This is the test that decided open question #22. The equation as *printed* lags only
    the bleed and does not converge at all; the prose reading -- the whole compressor exit
    flow taken from the previous interval -- halves its error when the step is halved.

    **`tol` is set to machine precision here on purpose, and that is not a widened
    tolerance -- it is the opposite.** This is a numerical-convergence experiment: it
    asks whether the scheme is a discretization of a continuous system, which is only
    answerable when the time step is the *only* error source. Since 2026-09-12 the
    shipped default is the report's printed 0.1 percent [pdf p.37], and at small dt that
    tolerance floor dominates the dt error, so a refinement study run at the default
    measures the solver rather than the discretization. Replication runs use the default;
    this experiment does not.
    """
    st0, wf0, qreq = _seed(400.0)
    stepfn = lambda t: wf0 if t < 0.5 else wf_pps_from_pph(775.0)  # noqa: E731

    vals = []
    for dt_ms in (4.0, 2.0, 1.0, 0.5):
        st, _, _ = _seed(400.0)
        tr = realtime.run(
            st,
            stepfn,
            AMB,
            duration_s=1.0,
            dt=dt_ms / 1000.0,
            q_req_ftlbf=qreq,
            integrate_np=False,
            tol=1e-12,
        )
        vals.append(tr["ng"][-1])

    diffs = np.abs(np.diff(vals))
    assert np.all(np.diff(diffs) < 0), f"errors are not shrinking with dt: {diffs}"
    for coarse, fine in zip(diffs, diffs[1:], strict=False):
        assert fine < 0.75 * coarse, (
            f"halving dt should roughly halve the error for a first-order scheme; "
            f"got {coarse:.1f} then {fine:.1f}"
        )


def test_both_readings_of_eq_74_converge_to_the_same_limit():
    """**Both** readings of Eq. 74 are consistent discretizations. Open question #22.

    This test used to assert the opposite -- that the equation *as printed* "converges to
    nothing" -- and #22 was closed on it, with the report recorded as contradicting
    itself. That was our bug, not the report's.

    Eq. 74 as printed [pdf p.36] lags the **bleed**: `WA31(n) = WA3(n) - WA3_bl(n-1)`. The
    implementation carried `WA31` in that branch instead of `WA3_bl`, making the recurrence
    `WA31(n) = WA3(n) - WA31(n-1)` -- an involution whose only fixed point is `WA3/2`.
    Everything else is a period-2 oscillation, so a *trimmed* engine did not even sit still
    there (T41ns alternated ~10600 / ~2200 degR), and no convergence conclusion drawn from
    it could mean anything. The carry is now the bleed, and the two readings are simply two
    different one-frame lags of the same opened loop.

    Measured at t = 1.0 s after a 400 -> 775 lbm/hr step, refining dt from 8 ms to 0.25 ms:

        prose    44024.0  43924.7  43880.1  43860.7  43851.5  43847.1
        printed  44045.1  43936.6  43885.5  43863.3  43852.8  43847.7

    Both halve their successive difference as dt halves -- first order, as an explicit
    Euler frame with a one-step lag must be -- and they agree to **0.6 rpm (0.0014 %)** at
    the finest step. At the report's own 7 ms they differ by about 0.05 %.

    The shipped model uses the prose reading, and that has not changed: it is the reading
    whose own words describe opening the loop ("the mass flow entering the combustor for a
    given interval is approximately equal to the mass flow leaving the compressor in the
    previous interval"). But it is chosen because the prose says so, **not** because the
    printed equation fails.
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
                tol=1e-10,
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
        limits[lag] = vals[-1]

    gap = abs(limits[True] - limits[False]) / limits[True] * 100.0
    assert gap < 0.01, (
        f"the two readings converge to different limits ({limits[True]:.1f} against "
        f"{limits[False]:.1f} rpm, {gap:.4f} %). They are two lags of one loop and should "
        f"agree in the limit; a real gap would mean one of them is not Eq. 74."
    )


def test_a_trimmed_engine_sits_still_in_both_readings_of_eq_74():
    """The check that would have caught the carry bug immediately.

    Neither reading of Eq. 74 changes the *equilibrium* -- at steady state the lagged
    quantity equals its own current value, whichever quantity it is. So a trimmed engine
    must sit still in both branches. Before the carry was fixed, the printed branch
    oscillated with a period of two frames and never sat anywhere.
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


def test_fuel_step_moves_the_engine_the_right_way():
    st, wf0, qreq = _seed(400.0)
    tr = realtime.run(
        st,
        lambda t: wf0 if t < 0.5 else wf_pps_from_pph(775.0),
        AMB,
        duration_s=3.0,
        q_req_ftlbf=qreq,
        integrate_np=False,
    )
    assert tr["ng"][-1] > tr["ng"][0], "more fuel must raise gas generator speed"
    assert tr["t41"].max() > tr["t41"][0], "and raise turbine inlet temperature"
    pre = tr["ng"][tr["t"] < 0.5]
    assert pre.max() - pre.min() < 1e-6, "the engine must be still before the step"


def test_transient_excursions_outside_the_maps_are_counted():
    """An open-loop step leaves the digitized data, and the model must say so.

    FAR reaches ~0.033 against f6's plotted 0.010-0.020. That is probably faithful -- a
    1988 table processor would clamp at the table ends too -- but it must never be
    silent. Open question #45.
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


@pytest.mark.parametrize("wf_pph", [476.3, 349.3, 267.7])
def test_every_published_trim_is_a_fixed_point(wf_pph: float):
    st, wf, qreq = _seed(wf_pph)
    tr = realtime.run(st, lambda t: wf, AMB, duration_s=1.0, q_req_ftlbf=qreq)
    assert abs(tr["ng"][-1] - tr["ng"][0]) < 1e-6


def test_the_two_to_one_multirate_costs_little_but_is_the_reports():
    """`FRAME_NP_S = 0.014` was declared and referenced nowhere in `src/` until 2026-09-13.

    [pdf p.47]: "the engine model is updated twice for each rotor routine cycle, or once
    every 7 msec", with the NP degree of freedom on the 14 ms rotor frame. Every NP
    integration here ran at the 7 ms engine frame instead -- a single-rate model wearing a
    multirate constant.

    `realtime.run(multirate=True)` is now the default and advances NP on alternate frames
    over the 14 ms it spans. What that is worth, measured on an open-loop 400 -> 775
    lbm/hr step with NP free against a constant load torque -- deliberately the worst case,
    since a free turbine with no aerodynamic damping runs far off design:

        NP    up to 81.6 rpm apart, 0.39 % of NP_DES
        NG    0.000 rpm apart

    NG is untouched to the bit, which is structural rather than lucky: NP enters the gas
    generator through nothing at all. Eq. 41's power-turbine torque depends on NP, but
    Eq. 45 differences only Q_GT and Q_C.

    And with `integrate_np=False` -- the configuration Figures 9 and 10 run in [pdf p.39]
    -- the switch is a **bit-identical no-op**, which is why no transient comparison in
    this repository moved when it was implemented.
    """
    wf0 = wf_pps_from_pph(400.0)
    r0 = trim.solve(wf0, c.NP_DES, AMB)
    f0 = frame(r0.state, wf0, AMB)

    def go(multirate: bool, integrate_np: bool):
        st = realtime.from_trim(r0, f0.wa31_pps, f0)
        return realtime.run(
            st,
            lambda t: wf0 if t < 0.5 else wf_pps_from_pph(775.0),
            AMB,
            duration_s=2.0,
            q_req_ftlbf=f0.q_pt_ftlbf,
            heat_sink=True,
            integrate_np=integrate_np,
            multirate=multirate,
        )

    single, multi = go(False, True), go(True, True)
    d_np = float(np.abs(multi["np"] - single["np"]).max())
    assert 50.0 < d_np < 120.0, (
        f"the multirate moves NP by {d_np:.1f} rpm, against 81.6 on record; if this has "
        f"collapsed to zero the 2:1 is no longer being applied"
    )
    assert np.array_equal(multi["ng"], single["ng"]), (
        "NP must not reach NG: Eq. 45 differences Q_GT and Q_C, neither of which sees NP"
    )

    a, b = go(False, False), go(True, False)
    assert all(np.array_equal(a[k], b[k]) for k in a), (
        "with the NP integration suppressed the multirate must be a bit-identical no-op, "
        "which is what makes every Figure 9 and 10 comparison independent of it"
    )
