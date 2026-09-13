"""Structural checks on the HMU -- Appendix C, Figures C9-C22.

These test the *wiring*, not agreement with any published result: that the droop line has
the signs the figures draw, that the limit cascade is the one Figure C18 draws and in that
order, that every block with memory actually has it, and that the schedules are wired to
the signals the figures wire them to.

Appendix C has no equations, so a transcription error here is a wrong picture, not a wrong
number, and a wrong picture does not announce itself. That is what these are for.
"""

from __future__ import annotations

import pytest

from t700.control import constants as c
from t700.control import hmu, schedules

DT = 0.007
"""The report's engine frame [pdf p.47]."""

BASE = dict(ng_rpm=41638.0, ps3_psia=176.34, t2_degR=518.67, pas_deg=100.0, xcpc_pct=50.0)


def settle(n: int = 4000, dt: float = DT, **over):
    u = hmu.HMUInputs(**{**BASE, **over})
    s = hmu.seed(u)
    o = None
    for _ in range(n):
        s, o = hmu.step(s, u, dt)
    return s, o


def test_a_held_input_reaches_a_steady_state():
    """Nothing in the HMU should drift once its inputs stop moving."""
    s1, o1 = settle(4000)
    s2, o2 = settle(8000)
    assert abs(o2.wf_pph - o1.wf_pph) < 1e-6, (
        f"Wf is still moving at 28 s: {o1.wf_pph:.6f} -> {o2.wf_pph:.6f}"
    )
    assert abs(s2.tm_integrator - s1.tm_integrator) < 1e-9


def test_more_collective_asks_for_more_fuel():
    """[Fig. C17] The load demand path is a *subtraction* from the topping line, so the
    schedules must be wired so that more demand means less subtraction. Getting the sign
    of this backwards inverts the entire control and still produces plausible numbers."""
    wf = [settle(xcpc_pct=x)[1].wf_pph for x in (30.0, 40.0, 50.0, 60.0)]
    assert all(b > a for a, b in zip(wf, wf[1:], strict=False)), wf


def test_more_power_available_spindle_asks_for_more_fuel():
    """[Fig. C25] `F_HM2` falls from 10.4 at PAS 28 deg to -0.1 at 120, and it is
    subtracted -- so the schedule's own slope carries the sign."""
    wf = [settle(pas_deg=p)[1].wf_pph for p in (90.0, 100.0, 110.0)]
    assert all(b > a for a, b in zip(wf, wf[1:], strict=False)), wf


def test_speed_droop_raises_fuel_flow():
    """The whole point of a droop line: KNDRP*(NGREF - PCNGHL) with NGREF above every
    operating speed, so a falling gas generator asks for more fuel."""
    hi = settle(ng_rpm=42000.0)[1]
    lo = settle(ng_rpm=40000.0)[1]
    assert lo.wfpdm > hi.wfpdm, f"droop demand did not rise as NG fell: {lo.wfpdm} {hi.wfpdm}"


def test_the_limit_cascade_is_the_one_figure_c18_draws():
    """`HMUSEL = MAX( MIN( MAX(WFPDM, WFIDM), WFPAC ), WFPDC )`.

    The order is drawn, not inferred, and it is not symmetric: the deceleration floor
    outranks the acceleration ceiling. A min/max soup in any other order gives the same
    answer almost everywhere and a different one exactly where it matters.
    """
    for xcpc in (10.0, 30.0, 50.0, 70.0, 90.0):
        _, o = settle(xcpc_pct=xcpc)
        expect = max(min(max(o.wfpdm, o.wfidm), o.wfpac), o.wfpdc)
        assert o.hmusel == pytest.approx(expect, abs=1e-12), (xcpc, o.hmusel, expect)


def test_the_acceleration_limit_actually_binds_at_high_demand():
    """`F_HM7` exists to cap a slam acceleration. If it never binds, either it is wired to
    the wrong signal or Figure C30's extraction is wrong."""
    _, o = settle(xcpc_pct=95.0)
    assert o.limit == "accel", o.limit
    assert o.hmusel == pytest.approx(o.wfpac, abs=1e-12)


def test_the_deceleration_floor_binds_at_low_demand():
    _, o = settle(xcpc_pct=10.0)
    assert o.limit == "decel", o.limit
    assert o.hmusel == pytest.approx(o.wfpdc, abs=1e-12)


def test_the_idle_floor_is_inactive_above_flight_idle():
    """[pdf p.38] The report eliminated fuel control below flight-idle power, so `WFIDM`
    should never win in any condition the engine model is valid over. Open question #10."""
    for xcpc in (10.0, 30.0, 50.0, 70.0, 90.0):
        _, o = settle(xcpc_pct=xcpc)
        assert o.limit != "idle", (xcpc, o.wfidm, o.wfpdm)
        assert o.wfidm < o.wfpdm, (xcpc, o.wfidm, o.wfpdm)


def test_the_metering_valve_clamps_and_delays():
    """[Fig. C19] WFMIN/WFMAX bracket the valve, then a lag, then a transport delay."""
    _, o = settle(xcpc_pct=99.0, pas_deg=120.0)
    assert c.WFMIN <= o.wf_pph <= c.WFMAX
    # the delay must actually delay: settle, then step the collective and check the
    # output does not move on the first frame
    s, before = settle(xcpc_pct=40.0)
    u = hmu.HMUInputs(**{**BASE, "xcpc_pct": 60.0})
    s1, o1 = hmu.step(s, u, DT)
    s2, o2 = hmu.step(s1, u, DT)
    s3, o3 = hmu.step(s2, u, DT)
    assert o1.wfmv_pph > before.wfmv_pph, "the valve command should respond at once"
    assert o1.wf_pph == pytest.approx(before.wf_pph, rel=1e-9), (
        f"the 0.015 s transport delay should hold the delivered flow for a frame: "
        f"{before.wf_pph} -> {o1.wf_pph}"
    )
    assert o3.wf_pph > o1.wf_pph, "and then let it through"


def test_the_torque_motor_integrator_respects_its_printed_limits():
    """[Fig. C10, Table C.1] XLOLIM/XHILIM are in in/sec [nomenclature, pdf p.81]."""
    for spdg in (-5.0, 0.0, 5.0):
        s, _ = settle(spdg=spdg, n=3000)
        assert c.XLOLIM - 1e-12 <= s.tm_integrator <= c.XHILIM + 1e-12, (spdg, s.tm_integrator)


def test_the_torque_motor_is_still_at_its_null():
    """[Fig. C10] The forward path sees `SPDG - 0.44`, so with the ECU absent and SPDG at
    the null the trim must not wander. If it does, the HMU is not runnable open-loop."""
    s, o = settle(spdg=hmu.SPDG_NULL, n=6000)
    assert abs(o.tmru) < 1e-6, f"torque motor drifted to TMRU={o.tmru} with no ECU"


def test_the_deadband_is_what_stops_the_torque_motor_drifting():
    """The null only holds because of the Fig. C10 deadband: without it the integrator
    would ramp on the -31.0 current bias. This pins the mechanism, not just the result."""
    assert c.TMDB > 0.0
    _, o = settle(spdg=hmu.SPDG_NULL + 0.2, n=6000)
    assert abs(o.tmru) > 1e-6, "a displaced SPDG should move the trim"


def test_hysteresis_blocks_actually_hold():
    """[Figs. C11, C13, C15] Three blocks carry backlash. A dither smaller than the band
    must not reach the output."""
    u = hmu.HMUInputs(**BASE)
    s = hmu.seed(u)
    for _ in range(2000):
        s, _ = hmu.step(s, u, DT)
    held = s.xldsh
    small = hmu.HMUInputs(
        **{**BASE, "xcpc_pct": BASE["xcpc_pct"] + 0.4 * c.XLDHYS / c.COLLECTIVE_RIGGING_GAIN}
    )
    s2, _ = hmu.step(s, small, DT)
    assert s2.xldsh == pytest.approx(held, abs=1e-12), "backlash let a sub-band move through"


def test_every_schedule_the_hmu_needs_is_loadable():
    """Seven of the eight; `F_EC1` belongs to the ECU."""
    for fn in (
        schedules.f_hm1,
        schedules.f_hm2,
        schedules.f_hm3,
        schedules.f_hm4,
        schedules.f_hm5,
        schedules.f_hm6,
    ):
        cur = fn()
        assert cur.x.size >= 2, cur.name
    m = schedules.f_hm7()
    assert m.params.size == 7, m.params
