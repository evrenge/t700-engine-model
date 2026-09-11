"""The station 4.1 heat-sink model, Eqs. 48-53 [pdf pp.25-26].

Ballin ran the engine in two configurations and published results for both, so this
model is a **switch**, not a feature. What the report says about which is which:

* `T41 = T41_ns` when no heat-sink representation is used [pdf p.23, below Eq. 22].
* Table 1's eigenvalues and Appendix B figures B1-B6 (2- and 5-DOF) are heat-sink
  **off**: "Because the station 4.1 heat-sink approximation was added to both models as
  a linear lead-lag representation, **it was not included in the analysis**" [pdf p.27].
  Table 1 has no T41 row, which corroborates it.
* Appendix B figures B7-B12 (3- and 6-DOF) are heat-sink **on**: "The heat-sink model is
  contained in the three- and six-degree-of-freedom models" [pdf p.67].
* Figures 9 and 10 are heat-sink **on**. The report never writes that sentence; it is
  inferred from p.20 ("a model of the nonadiabatic energy transfer at station 4.1 was
  required"), p.38 ("Because of the use of the heat-sink model..."), and p.47's
  discussion of Figure 10, which only parses if the model is running. Recorded as an
  inference in `docs/notes/heat-sink-configuration.md`, not as a transcription.

The sixth state is **gas** temperature at station 4.1 [pdf p.29], not metal temperature.
`T_m` exists only in Eqs. 48-49 and is eliminated when they collapse into Eq. 50.
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import constants as c
from t700 import maps, realtime, trim
from t700.engine import Ambient, frame
from t700.realtime import _heat_sink
from t700.units import wf_pps_from_pph

AMB = Ambient(14.696, 518.67)
TRIMS_PPH = (300.0, 400.0, 500.0, 600.0)


def _seed(wf_pph: float):
    wf = wf_pps_from_pph(wf_pph)
    r = trim.solve(wf, c.NP_DES, AMB)
    assert r.trustworthy, f"trim at {wf_pph} lbm/hr is not on data; fix the fixture"
    f = frame(r.state, wf, AMB)
    return wf, r, f


# --------------------------------------------------------------- the transfer function


def test_dc_gain_is_exactly_one():
    """Eq. 50 evaluated at s = 0 is 1, for any coefficients.

    This is the property the whole design rests on: it is why `engine.frame` and
    `trim.solve` need no heat-sink flag, and why every Table B.1 trim comparison is
    valid for both configurations at once.
    """
    for t41_ns in (1800.0, 2200.0, 2700.0):
        t41 = t41_ns
        # Drive a constant input to convergence; the output must return the input.
        x = 0.0
        for _ in range(20000):
            t41, x = _heat_sink(t41_ns, x, t41, 8.0, 91.0, 0.007)
        assert t41 == pytest.approx(t41_ns, rel=1e-9), (
            f"steady-state gain is {t41 / t41_ns:.12f}, not 1"
        )


def test_high_frequency_gain_is_one_minus_tau_ratio():
    """Instantaneous response to a step is `k = 1 - tau_b/tau_a` [Eqs. 50-53].

    Below 1 whenever tau_b > 0, which is what "response is significantly slowed"
    [pdf p.33] means numerically: fast excursions are attenuated, slow ones pass.
    """
    t41_prev, w41, ngc = 2200.0, 6.7, 91.0
    tau_a = c.TC_T41 * t41_prev**0.5 / w41**0.8  # (51)
    tau_b = float(maps.f_hs()(ngc)) / w41  # (52), (53)
    k = 1.0 - tau_b / tau_a

    x0 = (tau_b / tau_a) * 2200.0  # memory holding T41 = T41_ns at the old level
    stepped, _ = _heat_sink(2200.0 * 1.5, x0, t41_prev, w41, ngc, 0.007)
    jump = (stepped - 2200.0) / (2200.0 * 1.5 - 2200.0)
    assert jump == pytest.approx(k, rel=1e-9)
    assert 0.0 < k < 1.0, f"k={k} is outside (0,1); the lead-lag is not attenuating"


# ------------------------------------------------------------------- the switch itself


@pytest.mark.parametrize("wf_pph", TRIMS_PPH)
def test_heat_sink_is_inert_at_steady_state(wf_pph: float):
    """A trimmed engine held at constant fuel flow lands in the same place either way.

    Measured, not asserted from theory: this is what proves no trim comparison in this
    repository can be improved -- or broken -- by flipping the switch.
    """
    wf, r, f = _seed(wf_pph)
    out = {}
    for hs in (False, True):
        st = realtime.from_trim(r, f.wa31_pps, f)
        tr = realtime.run(
            st,
            lambda _t: wf,
            AMB,
            duration_s=8.0,
            dt=0.007,
            q_req_ftlbf=f.q_pt_ftlbf,
            integrate_np=False,
            heat_sink=hs,
        )
        out[hs] = tr
    for key in ("t41", "ng", "p3", "p41", "p45", "t45"):
        a, b = out[False][key][-1], out[True][key][-1]
        assert b == pytest.approx(a, rel=1e-9), f"{key} differs at steady state: {a} vs {b}"


@pytest.mark.parametrize("wf_pph", TRIMS_PPH)
def test_switch_off_reproduces_t41_ns_exactly(wf_pph: float):
    """With the switch off, T41 is T41_ns identically [Eq. 23, pdf p.23]."""
    wf, r, f = _seed(wf_pph)
    st = realtime.from_trim(r, f.wa31_pps, f)
    tr = realtime.run(
        st,
        lambda t: wf if t < 0.5 else wf * 1.3,
        AMB,
        duration_s=2.0,
        dt=0.007,
        q_req_ftlbf=f.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=False,
    )
    assert np.array_equal(tr["t41"], tr["t41_ns"])


def test_switch_on_separates_t41_from_t41_ns_only_in_transient():
    """On, T41 departs from T41_ns during the step and rejoins it after."""
    wf, r, f = _seed(400.0)
    st = realtime.from_trim(r, f.wa31_pps, f)
    tr = realtime.run(
        st,
        lambda t: wf if t < 0.5 else wf * 1.4,
        AMB,
        duration_s=12.0,
        dt=0.007,
        q_req_ftlbf=f.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=True,
    )
    gap = tr["t41"] - tr["t41_ns"]
    during = np.abs(gap[(tr["t"] > 0.5) & (tr["t"] < 1.5)]).max()
    after = abs(gap[-1])
    assert during > 50.0, f"heat sink barely engaged: max gap {during:.1f} deg R"
    assert after < 1.0, f"heat sink left a standing offset of {after:.3f} deg R"


def test_heat_sink_attenuates_the_t41_excursion():
    """On a step up, the heat-sink run must overshoot T41 less than the bare run.

    Direction only -- the magnitude against Figure 9 lives in `validation/`.
    """
    wf, r, f = _seed(400.0)
    peaks = {}
    for hs in (False, True):
        st = realtime.from_trim(r, f.wa31_pps, f)
        tr = realtime.run(
            st,
            lambda t: wf if t < 0.539 else wf_pps_from_pph(775.0),
            AMB,
            duration_s=5.0,
            dt=0.007,
            q_req_ftlbf=f.q_pt_ftlbf,
            integrate_np=False,
            heat_sink=hs,
        )
        settled = tr["t41"][(tr["t"] > 4.0) & (tr["t"] < 4.6)].mean()
        peaks[hs] = tr["t41"][(tr["t"] > 0.539) & (tr["t"] < 4.6)].max() - settled
    assert peaks[True] < peaks[False], (
        f"overshoot did not shrink: {peaks[False]:.1f} -> {peaks[True]:.1f} deg R"
    )


def test_seeding_without_a_frame_still_builds_a_usable_state():
    """`from_trim` without a frame is the 5-DOF path and must keep working."""
    wf, r, f = _seed(400.0)
    st = realtime.from_trim(r, f.wa31_pps)
    nxt, out = realtime.step(st, wf, AMB, 0.007, q_req_ftlbf=f.q_pt_ftlbf, integrate_np=False)
    assert out.t41_degR == pytest.approx(out.t41_ns_degR)
