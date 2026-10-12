"""The HMU against Table B.1: a consistency check, since no printed trace exists.

The report's closed-loop Figures 11-15 all need Gen Hel, which it does not contain. Instead,
hold the engine at each Table B.1 trim [pdf p.67] and bisect for the collective at which the
HMU, built only from Appendix C, commands the printed fuel flow. The two halves of the report
are digitized independently, so a sign, wiring or units error (`XCPC` is percent of maximum,
not inches) shows as a collective that is off the spindle's range, out of order, or on a cam
limit rather than the droop line.
"""

from __future__ import annotations

import pytest

from t700.control import hmu, schedules
from test_trim_points import CASES, STATE_B1

DT = 0.007
"""The report's engine frame [pdf p.47]."""

STANDARD_DAY_T2 = 518.67
"""Table B.1's trims are sea level standard [pdf p.67]."""

PAS_DEG = 100.0
"""The power available spindle is not printed for these trims: a fixed choice of ours, the
same at all three. See `test_the_result_does_not_depend_on_the_unprinted_power_lever`."""

# Table B.1 [pdf p.67], from `test_trim_points`: name, NG rpm, Ps3 psia, Wf lbm/hr
TRIMS = [(case.name, case.ng_rpm, STATE_B1[case.name]["ps3"], case.wf_pph) for case in CASES]

COLLECTIVE_RANGE_PCT = (5.0, 95.0)
"""The band a commanded collective must fall inside to be believable. Ours: the report prints
no collective for these trims. See SCOPE.md."""


def settle(ng_rpm, ps3_psia, xcpc_pct, n: int = 250, t2=STANDARD_DAY_T2, pas=PAS_DEG):
    u = hmu.HMUInputs(ng_rpm=ng_rpm, ps3_psia=ps3_psia, t2_degR=t2, pas_deg=pas, xcpc_pct=xcpc_pct)
    s = hmu.seed(u)
    out = None
    for _ in range(n):
        s, out = hmu.step(s, u, DT)
    return out


def collective_for(ng_rpm, ps3_psia, wf_pph, **kw) -> tuple[float, object]:
    """Bisect for the collective at which the HMU commands `wf_pph`."""
    lo, hi = 0.0, 100.0
    mid, out = 50.0, None
    # The settle is 250 frames because the slowest HMU lag is CLLDS = 0.2 s.
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        out = settle(ng_rpm, ps3_psia, mid, **kw)
        if out.wf_pph < wf_pph:
            lo = mid
        else:
            hi = mid
    return mid, out


@pytest.mark.parametrize("name,ng,ps3,wf", TRIMS, ids=lambda v: str(v))
def test_the_hmu_commands_each_printed_trim_from_a_believable_collective(name, ng, ps3, wf):
    xcpc, out = collective_for(ng, ps3, wf)
    assert out.wf_pph == pytest.approx(wf, rel=1e-6), out.wf_pph
    lo, hi = COLLECTIVE_RANGE_PCT
    assert lo < xcpc < hi, (
        f"{name}: the HMU needs {xcpc:.1f} % collective to command the printed "
        f"{wf} lbm/hr -- outside the believable band. Appendix C and Table B.1 are "
        f"digitized independently, so this is where a sign or units error surfaces."
    )


def test_the_result_does_not_depend_on_the_unprinted_power_lever():
    """Every trim stays on the droop line, collective falling hover > level > descent, at
    any power lever: PAS shifts all three together, so the order is the model's, not the
    choice's. A trim on a cam means the cam schedule or the droop line is wrong.
    """
    for pas in (90.0, 100.0, 110.0):
        xs = []
        for _, ng, ps3, wf in TRIMS:
            x, out = collective_for(ng, ps3, wf, pas=pas)
            if out.limit != "droop":
                pytest.fail(
                    f"PAS {pas}: {ng} rpm lands on the {out.limit} limit "
                    f"(WFPDM {out.wfpdm:.3f}, WFIDM {out.wfidm:.3f}, WFPAC {out.wfpac:.3f}, "
                    f"WFPDC {out.wfpdc:.3f})"
                )
            xs.append(x)
        assert xs[0] > xs[1] > xs[2], (pas, xs)
        if pas == PAS_DEG:
            assert xs[0] - xs[2] > 10.0, (
                f"only {xs[0] - xs[2]:.1f} points of collective separate hover from descent; "
                f"the load demand path is barely doing anything"
            )


def test_inlet_temperature_moves_the_collective_the_way_the_topping_line_does():
    """[Fig. C24] `F_HM1`, the droop line's ceiling, rises with T2, so a hot day needs less
    collective for the same fuel flow. Pins the direction: a reversed x axis on C24 would be
    invisible everywhere else in the model.
    """
    name, ng, ps3, wf = TRIMS[0]
    cold_t2, hot_t2 = 470.0, 560.0
    ptp_cold = float(schedules.f_hm1()(cold_t2))
    ptp_hot = float(schedules.f_hm1()(hot_t2))
    assert ptp_hot > ptp_cold, (
        f"F_HM1 no longer rises with T2 ({ptp_cold:.3f} at {cold_t2} vs {ptp_hot:.3f} at "
        f"{hot_t2}); Figure C24 may have been re-digitized with its axis reversed"
    )
    cold, _ = collective_for(ng, ps3, wf, t2=cold_t2)
    hot, _ = collective_for(ng, ps3, wf, t2=hot_t2)
    assert hot < cold, (
        f"the topping line is {ptp_hot - ptp_cold:.3f} higher at {hot_t2} deg R, so a hot "
        f"day should need *less* collective for the same fuel flow, got {hot:.1f} % "
        f"against {cold:.1f} %"
    )
