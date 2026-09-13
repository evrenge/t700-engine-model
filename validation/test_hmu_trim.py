"""The HMU against the report -- and there is no figure to check it against.

The report's only closed-loop results are Figures 11-15, and every one of them needs the
Gen Hel UH-60A blade-element simulation, which this report consumes and does not contain.
So Phase 5 has no published trace to overlay. `SCOPE.md` records that as a gate warning
rather than a surprise.

What *is* available is a consistency check across two independently digitized halves of
the report, and it is a real one. Table B.1 [pdf p.67] prints the engine's state at three
trims -- gas generator speed, station 3 static pressure, fuel flow. The HMU is built
entirely from Appendix C: 22 block diagrams and eight scheduling functions digitized from
Figures C23-C30. Nothing in Appendix C knows about Table B.1, and nothing in Table B.1
knows about Appendix C.

So: hold the engine at each printed trim state, and ask what collective position the HMU
needs in order to command the printed fuel flow. If the two halves disagree, the answer
comes out absurd -- off the spindle's range, or in the wrong order, or pinned against a
cam limit. It does not:

    hover          Wf 476.3 lbm/hr  ->  XCPC 52.75 %   on the droop line
    level 80 kt    Wf 349.3         ->  XCPC 35.87 %   on the droop line
    descent 80 kt  Wf 267.7         ->  XCPC 23.01 %   on the droop line

All three land on the **droop line** rather than on the idle, acceleration or deceleration
cams -- which is where a governor is supposed to sit in steady flight -- and the collective
falls monotonically with power, hover to level to descent. Nothing was tuned to make that
happen; the only free variable is the collective, and it is the one the pilot holds.

This is weaker than an overlay and stronger than nothing. It would catch an inverted sign
on the droop line, a schedule wired to the wrong signal, a units error in the collective
rigging (it did catch one -- `XCPC` is percent of maximum, not inches), or a limit cascade
in the wrong order.
"""

from __future__ import annotations

import pytest

from t700.control import hmu, schedules

DT = 0.007
"""The report's engine frame [pdf p.47]."""

STANDARD_DAY_T2 = 518.67
"""Table B.1's trims are sea level standard [pdf p.67]."""

PAS_DEG = 100.0
"""The power available spindle is **not printed** for these trims -- no part of the report
states the cockpit power lever position at the Table B.1 conditions. 100 deg is a fixed
choice, held the same across all three so the comparison is between trims rather than
against an absolute. Moving it slides all three collectives together; see
`test_the_result_does_not_depend_on_the_unprinted_power_lever`."""

# Table B.1 [pdf p.67]: name, NG rpm, Ps3 psia, Wf lbm/hr
TRIMS = [
    ("hover", 41638.0, 176.34, 476.3),
    ("level 80 kt", 39768.0, 142.13, 349.3),
    ("descent 80 kt", 38072.0, 114.27, 267.7),
]

COLLECTIVE_RANGE_PCT = (5.0, 95.0)
"""The band a commanded collective must fall inside to be believable. Ours, not the
report's -- it states no collective for these trims. Declared in `SCOPE.md`."""


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
    # 40 halvings of a 100-point range is 1e-10, far past what the bisection needs; the
    # settle is 250 frames because the slowest lag in the HMU is CLLDS = 0.2 s and the
    # output stops moving at the sixth decimal by frame 200.
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


@pytest.mark.parametrize("name,ng,ps3,wf", TRIMS, ids=lambda v: str(v))
def test_every_printed_trim_sits_on_the_droop_line(name, ng, ps3, wf):
    """Not on the idle floor, the acceleration ceiling or the deceleration floor.

    A governor in steady flight sits on its droop line; the cams are for transients. If a
    printed trim landed on a cam, either the cam schedule or the droop line is wrong.
    """
    _, out = collective_for(ng, ps3, wf)
    assert out.limit == "droop", (
        f"{name}: the HMU is on the {out.limit} limit at a printed steady trim "
        f"(WFPDM {out.wfpdm:.3f}, WFIDM {out.wfidm:.3f}, WFPAC {out.wfpac:.3f}, "
        f"WFPDC {out.wfpdc:.3f})"
    )


def test_collective_falls_monotonically_from_hover_to_descent():
    """The ordering is the check that survives the unprinted power lever: PAS shifts all
    three together, so their *order* is a property of the model rather than of the choice.
    """
    xs = [collective_for(ng, ps3, wf)[0] for _, ng, ps3, wf in TRIMS]
    assert xs[0] > xs[1] > xs[2], (
        f"hover / level / descent need {xs[0]:.1f} / {xs[1]:.1f} / {xs[2]:.1f} % "
        f"collective -- not monotone in power"
    )
    assert xs[0] - xs[2] > 10.0, (
        f"only {xs[0] - xs[2]:.1f} points of collective separate hover from descent; the "
        f"load demand path is barely doing anything"
    )


def test_the_result_does_not_depend_on_the_unprinted_power_lever():
    """`PAS_DEG` is a choice, so the conclusion must not rest on it.

    Sweeping the power lever over a wide band moves every commanded collective, but the
    ordering and the droop-line finding have to survive -- otherwise this file is reporting
    the choice rather than the model.
    """
    for pas in (90.0, 100.0, 110.0):
        xs = []
        for _, ng, ps3, wf in TRIMS:
            x, out = collective_for(ng, ps3, wf, pas=pas)
            if out.limit != "droop":
                pytest.fail(f"PAS {pas}: {ng} rpm lands on the {out.limit} limit")
            xs.append(x)
        assert xs[0] > xs[1] > xs[2], (pas, xs)


def test_inlet_temperature_moves_the_collective_the_way_the_topping_line_does():
    """[Fig. C24] `F_HM1` is the only schedule in the HMU driven by inlet temperature
    alone, and it sets the ceiling of the droop line. Its slope is the thing worth pinning,
    because a reversed x axis on C24 would be invisible everywhere else in the model.

    It **rises** with T2 -- 0.24 at 395 deg R to 3.42 at 500, then roughly flat to 3.35 at
    634 -- which is the right sign for a Wf/Ps3 limit: at a hotter, less dense inlet the
    engine needs a larger fuel-to-pressure ratio for the same power, so the ceiling has to
    rise with it. A higher ceiling means the load demand has less to make up, so the
    collective needed for a *fixed* fuel flow moves opposite to the topping line.

    My first reading of this had the sign backwards and the model caught it, which is the
    argument for testing the direction rather than asserting the number.
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
