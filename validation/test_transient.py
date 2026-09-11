"""The real-time transient model: Eqs. 69-80.

Two things are checked here that no steady-state comparison can reach: that the
quasi-steady frame reproduces the trim as a fixed point, and that the opened compressor
iteration is a *consistent discretization* -- that refining the time step converges.
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import maps, realtime, trim
from t700.engine import Ambient, frame
from t700.units import wf_pps_from_pph

AMB = Ambient(14.696, 518.67)


def _seed(wf_pph: float):
    wf = wf_pps_from_pph(wf_pph)
    r = trim.solve(wf, 20895.0, AMB)
    f = frame(r.state, wf, AMB)
    return realtime.from_trim(r, f.wa31_pps), wf, f.q_pt_ftlbf


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
        )
        vals.append(tr["ng"][-1])

    diffs = np.abs(np.diff(vals))
    assert np.all(np.diff(diffs) < 0), f"errors are not shrinking with dt: {diffs}"
    for coarse, fine in zip(diffs, diffs[1:], strict=False):
        assert fine < 0.75 * coarse, (
            f"halving dt should roughly halve the error for a first-order scheme; "
            f"got {coarse:.1f} then {fine:.1f}"
        )


def test_the_printed_reading_of_eq_74_does_not_converge():
    """Records *why* the prose reading was chosen, as an executable fact.

    If this ever starts passing convergence, open question #22 must be reopened.
    """
    _, wf0, qreq = _seed(400.0)
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
            lag_whole_flow=False,
        )
        vals.append(tr["ng"][-1])

    diffs = np.abs(np.diff(vals))
    assert not np.all(np.diff(diffs) < 0), (
        "the printed reading of Eq. 74 now appears to converge -- reopen question #22"
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
