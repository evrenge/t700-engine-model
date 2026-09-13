"""Engine, ECU and HMU closed on one another, against Table B.1.

**This is the only quantitative check Phase 5 gets, and it is a real one.**

The report's closed-loop results are Figures 11-15 and every one needs the Gen Hel UH-60A
blade-element simulation, which this report consumes and does not contain. So there is no
published closed-loop trace to overlay, and that was known before Phase 5 started --
`SCOPE.md` records it as a gate warning rather than a surprise.

What there is, is Table B.1 [pdf p.67]: the complete engine state at three trims, printed
as numbers, including the power turbine speed and the shaft torque. Hold the load at the
printed torque and set the speed reference, and the loop must settle at the printed state
**with fuel flow as an output** rather than an input. That is a materially harder test than
anything in Phases 0-4, where Wf was prescribed and the engine only had to agree about what
it produced.

It is also a genuine cross-check rather than a restatement: the engine comes from the
report's body and Appendix A, the control from Appendix C's 22 block diagrams and eight
digitized schedules, and neither half knows the other exists.

## The load matters, and it is not zero

A free power turbine against a speed-*independent* load has no aerodynamic damping at all:
`dQreq/dNP = 0` leaves the governor as the only thing holding NP, which is not the system
the report models. Appendix B says so -- its NP diagonal carries a ~52 % residual against
ours precisely because `dQreq/dNP` is the other Gen Hel quantity. Open question #6 recovers
it in ratio form, **1.78 / 1.45 / 1.18 times `|dQ_PT/dNP|`** at hover / level / descent,
positive because rotor torque rises with speed. That ratio is what this file uses, together
with the load inertia `J_LOAD_UH60A` recovered the same way. Both are properties of the
airframe Appendix B was trimmed on, and both are cited, not invented.

## What the governor is supposed to absorb

The collective is a *pilot* input, and in a governed system it must not set the trim: it
sets where the ECU's integrator sits while the engine goes where the load tells it. That is
tested here explicitly, and it is the strongest statement in the file -- over a 35 to 70 %
collective sweep the settled state moves by less than 0.2 % while `SPDG` travels from -0.40
to +1.62.
"""

from __future__ import annotations

import pytest

from t700 import constants as engine_c
from t700 import maps, trim
from t700.control import loop
from t700.engine import Ambient, State, frame
from t700.units import shp_from_torque, wf_pps_from_pph

AMB = Ambient(14.696, 518.67)
DT = 0.007
"""The report's engine frame [pdf p.47]."""

SETTLE_FRAMES = 4000
"""28 s. The slowest thing in the loop is the ECU's P+I integrator against the rotor
inertia; everything is inside 0.05 % of its final value by 20 s."""

NP_TRIM_RPM = 20895.0
"""Table B.1's power turbine speed, printed at all three trims [pdf p.67]."""

CLOSED_LOOP_TOL_PCT = 1.0
"""Ours, declared in `SCOPE.md`. Generous against the measured worst of 0.74 %, which is
itself inherited: the descent shaft power is 0.72 % out open-loop too, so the control adds
essentially nothing to it."""

# name, Wf lbm/hr, NG rpm, Ps3 psia, shaft torque ft*lbf, dQreq/dNP ratio (#6), shp
TRIMS = [
    ("hover", 476.3, 41638.0, 176.34, 229.0, 1.78, 911.1),
    ("level 80 kt", 349.3, 39768.0, 142.13, 138.9, 1.45, 552.6),
    ("descent 80 kt", 267.7, 38072.0, 114.27, 76.06, 1.18, 302.6),
]

XCPC_PCT = {"hover": 52.75, "level 80 kt": 35.87, "descent 80 kt": 23.01}
"""Collective at each trim, from `test_hmu_trim.py`'s own droop-line solve.

**All three ran at the hover collective, 52.75 %, until 2026-09-13.** The settled state
barely notices -- the governor is what absorbs a wrong collective, and that is the thing
`test_the_governor_absorbs_the_collective` is about -- but the *path* there does not:

    level,   52.75 %   Wf peak +57.4 %,  NP peak +2.04 %   f6 x46
    level,   35.87 %   Wf peak  +0.6 %,  NP peak +0.03 %   none
    descent, 52.75 %   Wf peak +94.1 %,  NP peak +3.53 %   f6 x215, f1@85 x10, f1@82 x228
    descent, 23.01 %   Wf peak  +0.7 %,  NP peak +0.03 %   none

(Measured after the `ecu.seed` fix of the same day. The control-systems audit that found
this measured +78.7 / +107.9 % on the peaks against the unfixed seed, which was adding its
own +31.8 % slam on top; the f1 clamping is the part that matters and it is unchanged in
kind.)

So every level and descent run in what `SCOPE.md` calls the Phase 5 gate was leaving the
digitized compressor map -- 273 evaluations off `f1`'s 82 % speed line at descent -- and
nothing reported it, because `maps.reset_clamps()` was called only inside the acceleration
test. Found by the 2026-09-13 control-systems audit. The settled deviations against
Table B.1 are unchanged by the fix; what changes is that the model now stays inside its
own data on the way there."""


def _dqpt_dnp(result, wf_pps: float) -> float:
    """dQ_PT/dNP by central difference, the report's own method for a derivative."""
    h = 0.002 * result.state.np_rpm
    q = []
    for d in (-h, h):
        st = State(
            result.state.ng_rpm,
            result.state.np_rpm + d,
            result.state.p3_psia,
            result.state.p41_psia,
            result.state.p45_psia,
        )
        q.append(frame(st, wf_pps, AMB).q_pt_ftlbf)
    return (q[1] - q[0]) / (2.0 * h)


def settle(wf_pph, q_shaft, ratio, xcpc_pct=None, n=SETTLE_FRAMES, pcprf=None, name=None):
    """Close the loop from the open-loop trim and run it to steady state.

    `xcpc_pct` defaults to the collective that trim actually needs -- see `XCPC_PCT`.
    """
    if xcpc_pct is None:
        xcpc_pct = XCPC_PCT[name] if name in XCPC_PCT else 52.75
    wf = wf_pps_from_pph(wf_pph)
    r = trim.solve(wf, engine_c.NP_DES, AMB)
    f = frame(r.state, wf, AMB)
    slope = abs(_dqpt_dnp(r, wf)) * ratio

    def load(np_rpm: float) -> float:
        return q_shaft + slope * (np_rpm - NP_TRIM_RPM)

    pilot = loop.Pilot(
        xcpc_pct=xcpc_pct,
        pas_deg=100.0,
        pcprf_pct=pcprf if pcprf is not None else NP_TRIM_RPM * 100.0 / engine_c.NP_DES,
    )
    s = loop.seed(r, f.wa31_pps, f, pilot, AMB)
    ecu_o = hmu_o = fr = None
    for _ in range(n):
        s, ecu_o, hmu_o, fr = loop.step(
            s, pilot, AMB, dt=DT, load=load, j_load=engine_c.J_LOAD_UH60A, heat_sink=True
        )
    return s, ecu_o, hmu_o, fr


@pytest.mark.parametrize("name,wf,ng,ps3,q,ratio,shp", TRIMS, ids=lambda v: str(v))
def test_the_closed_loop_settles_on_the_printed_trim(name, wf, ng, ps3, q, ratio, shp):
    """Fuel flow is an **output** here. Nothing tells the loop what it should be."""
    s, _e, h, fr = settle(wf, q, ratio, name=name)
    got = {
        "NG": s.engine.ng_rpm,
        "NP": s.engine.np_rpm,
        "Wf": h.wf_pph,
        "Ps3": engine_c.K_PS3 * s.engine.p3_psia,
        "shp": shp_from_torque(fr.q_pt_ftlbf, s.engine.np_rpm),
    }
    want = {"NG": ng, "NP": NP_TRIM_RPM, "Wf": wf, "Ps3": ps3, "shp": shp}
    for key, w in want.items():
        dev = 100.0 * (got[key] / w - 1.0)
        assert abs(dev) < CLOSED_LOOP_TOL_PCT, (
            f"{name} {key}: closed loop settles at {got[key]:.2f} against Table B.1's "
            f"{w}, {dev:+.3f} %"
        )


@pytest.mark.parametrize("name,wf,ng,ps3,q,ratio,shp", TRIMS, ids=lambda v: str(v))
def test_every_printed_trim_governs_on_the_droop_line(name, wf, ng, ps3, q, ratio, shp):
    """Closed loop, as open loop: a governor in steady flight sits on its droop line, not
    on a cam. And the ECU should be governing speed, not limiting temperature."""
    _s, e, h, _fr = settle(wf, q, ratio, name=name)
    assert h.limit == "droop", (name, h.limit)
    assert e.limiting == "speed", (name, e.limiting, e.tsig, e.spdsf)


def test_the_governor_absorbs_the_collective():
    """The strongest claim in this file, and the one a broken loop fails first.

    A collective is a pilot input. In a governed system it must not set the trim -- the
    load does that, and the collective only moves where the ECU's integrator sits. Sweep
    it nearly twofold and the settled engine should barely notice, while `SPDG` travels.
    """
    _name, wf, ng, _ps3, q, ratio, _shp = TRIMS[0]
    seen = []
    for xcpc in (35.0, 45.0, 52.75, 60.0, 70.0):
        s, e, h, _fr = settle(wf, q, ratio, xcpc_pct=xcpc)
        seen.append((xcpc, s.engine.ng_rpm, s.engine.np_rpm, h.wf_pph, e.spdg))
    ngs = [r[1] for r in seen]
    nps = [r[2] for r in seen]
    spdgs = [r[4] for r in seen]
    assert max(ngs) - min(ngs) < 0.005 * ng, (
        f"NG moved {max(ngs) - min(ngs):.1f} rpm across the collective sweep; the "
        f"governor is not absorbing it"
    )
    assert max(nps) - min(nps) < 0.005 * NP_TRIM_RPM
    assert max(spdgs) - min(spdgs) > 1.0, (
        f"SPDG only moved {max(spdgs) - min(spdgs):.3f} across the sweep -- if the "
        f"collective is not moving the trim signal, this test is proving nothing"
    )


def test_the_speed_reference_is_what_sets_power_turbine_speed():
    """[Fig. C1] The governor drives sensed NP to `PCPRF`, so the trimmed speed is the
    cockpit reference and not `NP_DES`. Table B.1's 20895 against a design 20900 is
    99.976 %, which is a *prediction* of the loop rather than an input to it.
    """
    _name, wf, _ng, _ps3, q, ratio, _shp = TRIMS[0]
    for pcprf, expect in ((99.0, 0.99 * engine_c.NP_DES), (100.0, engine_c.NP_DES)):
        s, _e, _h, _fr = settle(wf, q, ratio, pcprf=pcprf)
        dev = 100.0 * (s.engine.np_rpm / expect - 1.0)
        assert abs(dev) < 0.1, (
            f"with PCPRF = {pcprf} % the loop should hold NP at {expect:.0f} rpm, "
            f"got {s.engine.np_rpm:.0f} ({dev:+.3f} %)"
        )


def test_a_load_step_is_rejected():
    """What a governor is for. Raise the demanded torque and NP must come back.

    This is the only dynamic statement available without Gen Hel, and it is a weak one --
    a direction and a settling time, not a trace. Recorded as such.
    """
    _name, wf, _ng, _ps3, q, ratio, _shp = TRIMS[0]
    s, _e, _h, _fr = settle(wf, q, ratio)
    np_before = s.engine.np_rpm

    r = trim.solve(wf_pps_from_pph(wf), engine_c.NP_DES, AMB)
    slope = abs(_dqpt_dnp(r, wf_pps_from_pph(wf))) * ratio
    step_q = q * 1.15

    def load(np_rpm: float) -> float:
        return step_q + slope * (np_rpm - NP_TRIM_RPM)

    pilot = loop.Pilot(
        xcpc_pct=52.75, pas_deg=100.0, pcprf_pct=NP_TRIM_RPM * 100.0 / engine_c.NP_DES
    )
    lo = np_before
    wf_after = None
    for i in range(3000):
        s, _e, h, _fr = loop.step(
            s, pilot, AMB, dt=DT, load=load, j_load=engine_c.J_LOAD_UH60A, heat_sink=True
        )
        lo = min(lo, s.engine.np_rpm)
        if i > 2000:
            wf_after = h.wf_pph
    assert lo < np_before, "a load step should droop the power turbine before recovery"
    dev = 100.0 * (s.engine.np_rpm / np_before - 1.0)
    assert abs(dev) < 0.5, f"NP did not recover after a 15 % load step: {dev:+.3f} %"
    assert wf_after > 1.05 * wf, (
        f"a 15 % load step should raise fuel flow well above {wf}, got {wf_after:.1f}"
    )


def test_the_acceleration_limit_keeps_the_engine_inside_its_digitized_envelope():
    """This is open question #45, and closing the loop answers it.

    An open-loop 400 -> 775 lbm/hr step drives the model well outside the data it was
    built from: the compressor map, the gas generator turbine energy function and the
    power turbine flow function all clamp. That has been recorded as a limitation since
    Phase 4.

    But it is an artifact of prescribing a fuel step **the control would never command**.
    `F_HM7` exists precisely to cap fuel during an acceleration, and with the loop closed a
    *harder* demand than Figure 9's -- collective slammed 52.75 -> 95 % and the power lever
    to 120 deg -- leaves `f1`, `f7` and `f9` untouched. The acceleration limit binds for
    about a hundred frames and does its job.

    So the envelope excursion is a property of the open-loop test, not of the model.
    """
    _name, wf, _ng, _ps3, q, ratio, _shp = TRIMS[0]
    wf_pps = wf_pps_from_pph(wf)
    r = trim.solve(wf_pps, engine_c.NP_DES, AMB)
    f = frame(r.state, wf_pps, AMB)
    slope = abs(_dqpt_dnp(r, wf_pps)) * ratio

    def load(np_rpm: float) -> float:
        return q + slope * (np_rpm - NP_TRIM_RPM)

    cruise = loop.Pilot(
        xcpc_pct=52.75, pas_deg=100.0, pcprf_pct=NP_TRIM_RPM * 100.0 / engine_c.NP_DES
    )
    s = loop.seed(r, f.wa31_pps, f, cruise, AMB)
    for _ in range(300):
        s, _e, _h, _fr = loop.step(
            s, cruise, AMB, dt=DT, load=load, j_load=engine_c.J_LOAD_UH60A, heat_sink=True
        )

    slam = loop.Pilot(xcpc_pct=95.0, pas_deg=120.0, pcprf_pct=NP_TRIM_RPM * 100.0 / engine_c.NP_DES)
    maps.reset_clamps()
    limits = {}
    for _ in range(1200):
        s, _e, h, _fr = loop.step(
            s, slam, AMB, dt=DT, load=load, j_load=engine_c.J_LOAD_UH60A, heat_sink=True
        )
        limits[h.limit] = limits.get(h.limit, 0) + 1

    report = {k: v for k, v in maps.clamp_report().items() if v}
    escaped = {k: v for k, v in report.items() if k.startswith(("f1", "f7", "f9"))}
    assert not escaped, (
        f"a closed-loop slam left the digitized envelope: {escaped}. The acceleration "
        f"limit should have prevented that -- open question #45."
    )
    assert limits.get("accel", 0) > 20, (
        f"the acceleration limit never bound during a slam ({limits}), so this test is "
        f"not exercising what it claims to"
    )


@pytest.mark.parametrize("name,wf,ng,ps3,q,ratio,shp", TRIMS, ids=lambda v: str(v))
def test_the_closed_loop_stays_inside_its_digitized_data(name, wf, ng, ps3, q, ratio, shp):
    """The Phase 5 gate must not reach its answer through extrapolated maps.

    Nothing checked this until 2026-09-13: `maps.reset_clamps()` was called only inside
    the acceleration test, so every level and descent run here left the compressor map
    unreported. The cause was the collective -- see `XCPC_PCT` -- and at descent it was
    273 evaluations off `f1`'s 82 % speed line, in the comparison `SCOPE.md` calls the
    gate for Phase 5.

    `f6` is the one clamp that is allowed, and only above about 590 lbm/hr: it is a
    two-point, nearly constant table (0.98504 to 0.98496) whose clamp costs 1e-4, and the
    top of the power range cannot avoid it. See `trim.TrimResult.extrapolated_tables`.
    """
    maps.reset_clamps()
    settle(wf, q, ratio, name=name)
    report = maps.clamp_report()
    maps.reset_clamps()

    offenders = {k: v for k, v in report.items() if not k.startswith("f6")}
    assert not offenders, (
        f"{name}: the closed loop reads {offenders} outside their digitized range on the "
        f"way to its settled state. Every number this file reports would then be resting "
        f"on extrapolated maps."
    )
