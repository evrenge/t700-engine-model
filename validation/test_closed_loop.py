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

import numpy as np
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

CLOSED_LOOP_TOL_PCT = 0.5
"""Ours, declared in `SCOPE.md`. This is the Phase 5 gate.

Measured worst **0.337 %**, descent fuel flow, as a cycle mean over the fifteen
comparisons. **Tightened 1.0 -> 0.5 on 2026-09-14**, never widened, and the model did not
change: what changed is the statistic. The loop settles to a limit cycle rather than to a
point (open question #61), and this file compared a *terminal sample* of it -- one phase,
depending on where the run happened to stop. That number is 0.876 %, and versions of it
have been quoted as the headline closed-loop result at 0.74, 0.696, 0.62 and 0.865 % as
the run length and the seeding moved around. The cycle mean is 0.337 % and does not move.

The worst is still descent fuel flow and it is still partly inherited: descent shaft power
is -0.900 % out open-loop. The governor converts that error rather than adding to it --
closed loop, descent shp is -0.019 % and the fuel flow carries it instead, because fuel is
an output here and the loop trims it until the torque matches."""

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


CYCLE_FRAMES = 1500
"""How many trailing frames to average a closed-loop result over: 10.5 s at 7 ms.

**The loop does not settle to a point. It settles to a limit cycle**, and a terminal
sample is a phase of that cycle rather than a steady state -- see open question #61, where
it is traced to `CH`, the printed 0.05 %NG hysteresis on the gas-generator speed sensor
[Fig. C15, pdf p.91]. Measured peak-to-peak at the three trims: NP 7.5 / 12.0 / 16.7 rpm,
NG 10.7 / 24.8 / 39.8 rpm, period 5.8 / 4.9 / 4.8 s. It survives step refinement to
1.75 ms, so it is the model's and not the integrator's.

1500 frames covers at least two full periods at every trim, so the mean is a cycle mean
and not a phase. It matters: the worst deviation from Table B.1 is **0.308 %** as a cycle
mean and **0.876 %** as a terminal sample, and this project quoted the terminal number --
0.62 %, then 0.865 % -- for as long as it has had a closed loop.
"""


def mean_over_cycle(values) -> float:
    """Mean of a trailing window that spans a whole number of limit-cycle periods."""
    return float(np.mean(np.asarray(values)[-CYCLE_FRAMES:]))


def settle(wf_pph, q_shaft, ratio, xcpc_pct=None, n=SETTLE_FRAMES, pcprf=None, name=None):
    """Close the loop from the open-loop trim and run it to steady state.

    `xcpc_pct` defaults to the collective that trim actually needs -- see `XCPC_PCT`.

    Returns the FINAL state. Use `settle_cycle` where a number is being compared against
    the report: the final state is one phase of a limit cycle (see `CYCLE_FRAMES`).
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


def settle_cycle(wf_pph, q_shaft, ratio, name=None, n=SETTLE_FRAMES):
    """As `settle`, but returning the **mean over the last whole cycles** of each output.

    Also returns the peak-to-peak of each, so a caller can say how much of a deviation is
    the cycle and how much is the model.
    """
    xcpc = XCPC_PCT.get(name, 52.75)
    wf = wf_pps_from_pph(wf_pph)
    r = trim.solve(wf, engine_c.NP_DES, AMB)
    f = frame(r.state, wf, AMB)
    slope = abs(_dqpt_dnp(r, wf)) * ratio
    pilot = loop.Pilot(
        xcpc_pct=xcpc,
        pas_deg=100.0,
        pcprf_pct=NP_TRIM_RPM * 100.0 / engine_c.NP_DES,
    )
    st = loop.seed(r, f.wa31_pps, f, pilot, AMB)
    hist: dict[str, list[float]] = {k: [] for k in ("NG", "NP", "Wf", "Ps3", "shp")}
    for _ in range(n):
        st, _e, h, fr = loop.step(
            st,
            pilot,
            AMB,
            dt=DT,
            load=lambda v: q_shaft + slope * (v - NP_TRIM_RPM),
            j_load=engine_c.J_LOAD_UH60A,
            heat_sink=True,
        )
        hist["NG"].append(st.engine.ng_rpm)
        hist["NP"].append(st.engine.np_rpm)
        hist["Wf"].append(h.wf_pph)
        hist["Ps3"].append(engine_c.K_PS3 * st.engine.p3_psia)
        hist["shp"].append(shp_from_torque(fr.q_pt_ftlbf, st.engine.np_rpm))
    mean = {k: mean_over_cycle(v) for k, v in hist.items()}
    ptp = {k: float(np.ptp(np.asarray(v)[-CYCLE_FRAMES:])) for k, v in hist.items()}
    return mean, ptp


@pytest.mark.parametrize("name,wf,ng,ps3,q,ratio,shp", TRIMS, ids=lambda v: str(v))
def test_the_closed_loop_governs_to_the_printed_trim(name, wf, ng, ps3, q, ratio, shp):
    """Fuel flow is an **output** here. Nothing tells the loop what it should be.

    Compared as a **cycle mean**, not as a terminal sample: the loop settles to a limit
    cycle rather than to a point (open question #61, `CYCLE_FRAMES`). This test was named
    `..._settles_on_the_printed_trim` and took the last frame until 2026-09-14, which is
    both the wrong statistic and the wrong word.
    """
    got, ptp = settle_cycle(wf, q, ratio, name=name)
    want = {"NG": ng, "NP": NP_TRIM_RPM, "Wf": wf, "Ps3": ps3, "shp": shp}
    for key, w in want.items():
        dev = 100.0 * (got[key] / w - 1.0)
        assert abs(dev) < CLOSED_LOOP_TOL_PCT, (
            f"{name} {key}: closed loop governs to {got[key]:.2f} against Table B.1's "
            f"{w}, {dev:+.3f} % (cycle peak-to-peak {ptp[key]:.3f})"
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


def test_the_limit_cycle_is_what_it_is_measured_to_be():
    """Pin the limit cycle, so it cannot grow or vanish unnoticed.

    It is a property of the report's control system as specified, not a defect in this
    implementation, and open question #61 sets out the three measurements behind that: it
    survives step refinement to 1.75 ms, it collapses when `CH` alone is removed, and its
    amplitude scales with the printed band width.

    `CH` is the gas-generator speed sensor's hysteresis, 0.05 %NG [Fig. C15, pdf p.91] --
    22 rpm on NG_DES. The governor cannot resolve speed more finely than that, and the loop
    cycles at 0.5 to 1.8 times the band. Every closed-loop number in this project is a mean
    over it; `CYCLE_FRAMES` is why.
    """
    want = {"hover": (7.5, 10.9), "level 80 kt": (12.0, 24.8), "descent 80 kt": (16.7, 39.8)}
    for name, wf, _ng, _ps3, q, ratio, _shp in TRIMS:
        _mean, ptp = settle_cycle(wf, q, ratio, name=name)
        np_pp, ng_pp = want[name]
        assert ptp["NP"] == pytest.approx(np_pp, abs=1.0), f"{name} NP p-p {ptp['NP']:.3f}"
        assert ptp["NG"] == pytest.approx(ng_pp, abs=2.0), f"{name} NG p-p {ptp['NG']:.3f}"
        # and it must remain small enough that the mean is the honest summary
        assert 100.0 * ptp["NG"] / _ng < 0.15, f"{name} NG cycle is {ptp['NG']:.1f} rpm"


# ---------------------------------------------- Appendix C's forward path, which nothing saw

LOAD_STEP_PCT = 15.0
NP_DIP_RPM = (193.4, 195.4)
WF_PEAK_PPH = (550.0, 558.0)
SLAM_RISE_MS = (140.0, 155.0)
"""Characterizations of the loop's *transient*, not accuracy claims -- see `SCOPE.md`,
"A third category exists". The report prints no closed-loop time history, so there is
nothing to be accurate against; what these bounds do is make Appendix C's forward-path
gains observable at all.

**They were not.** The 2026-09-14 validation-quality audit mutated 30 constants and found
nine that no test in the project detected, all in the fuel control. The mechanism is
structural: `test_hmu_trim` bisects the *collective* until the HMU commands the printed
fuel flow, so a forward-path gain change is absorbed into the collective; and the
closed-loop tests let the governor trim it out in steady state. Every steady test either
solves for a free input or integrates the error away, so no forward-path gain in Appendix C
was observable anywhere.

A transient is where a loop gain shows. Measured:

| perturbation | NP dip | Wf peak | slam rise |
|---|---|---|---|
| baseline | 194.4 rpm | 554.1 pph | 147 ms |
| `TMLG` +10 % | **189.3** | **556.9** | 147 |
| `TMGN` x10 | 193.1 | 552.8 | 147 |
| `TMGN` +10 % | 194.3 | 554.1 | 147 |
| `CLLDS` x2 | 194.4 | 554.1 | **175** |
| `TLGE` +50 % | 194.4 | 554.1 | 147 (bit-identical) |

So the load step pins `TMLG` at +10 % and `TMGN` at the order-of-magnitude level, and the
collective slam pins `CLLDS` at x2. The NP dip is the discriminator and its band is 1 % wide
about the measured 194.392 rpm, which makes it a tripwire rather than a tolerance -- any
real change to the loop will trip it and should be looked at. `TMGN` at +10 % moves the dip
by 0.07 rpm and stays out of reach even so. `TLGE` is bit-identical in every
configuration the suite runs, and `test_the_t45_harness_lag_is_structurally_unreachable`
records why."""


def test_a_load_step_response_pins_the_torque_motor_linkage():
    """A 15 % load step at hover, characterized -- so a forward-path gain cannot move freely.

    The quantities are the ones a governor is judged on and the ones a loop gain moves: how
    far NP is dragged down before the fuel catches up, and how much fuel it takes.
    """
    name, wf, _ng, _ps3, q, ratio, _shp = TRIMS[0]
    xcpc = XCPC_PCT[name]
    wf_pps = wf_pps_from_pph(wf)
    r = trim.solve(wf_pps, engine_c.NP_DES, AMB)
    f = frame(r.state, wf_pps, AMB)
    slope = abs(_dqpt_dnp(r, wf_pps)) * ratio
    pilot = loop.Pilot(
        xcpc_pct=xcpc, pas_deg=100.0, pcprf_pct=NP_TRIM_RPM * 100.0 / engine_c.NP_DES
    )
    s = loop.seed(r, f.wa31_pps, f, pilot, AMB)
    step_at, n = 1400, 2600
    np_min, wf_max = 1e9, 0.0
    for i in range(n):
        load_q = q * (1.0 + LOAD_STEP_PCT / 100.0) if i >= step_at else q
        s, _e, h, _fr = loop.step(
            s,
            pilot,
            AMB,
            dt=DT,
            load=lambda v, lq=load_q: lq + slope * (v - NP_TRIM_RPM),
            j_load=engine_c.J_LOAD_UH60A,
            heat_sink=True,
        )
        if i >= step_at:
            np_min = min(np_min, s.engine.np_rpm)
            wf_max = max(wf_max, h.wf_pph)
    dip = NP_TRIM_RPM - np_min
    assert NP_DIP_RPM[0] < dip < NP_DIP_RPM[1], f"NP dips {dip:.1f} rpm on a 15 % load step"
    assert WF_PEAK_PPH[0] < wf_max < WF_PEAK_PPH[1], f"Wf peaks at {wf_max:.1f} pph"


def test_a_collective_slam_pins_the_collective_lag():
    """`CLLDS` is the load-demand-spindle lag, and only a *moving* collective can see it.

    Every other test in the project holds the collective fixed, which is why doubling this
    constant was undetectable. Slam it from the hover setting to 85 % and the fuel command's
    10-90 % rise time is the observable: 147 ms as shipped, 175 ms at 2x.
    """
    name, wf, _ng, _ps3, q, ratio, _shp = TRIMS[0]
    wf_pps = wf_pps_from_pph(wf)
    r = trim.solve(wf_pps, engine_c.NP_DES, AMB)
    f = frame(r.state, wf_pps, AMB)
    slope = abs(_dqpt_dnp(r, wf_pps)) * ratio
    s, out, at, n = None, [], 600, 1800
    for i in range(n):
        pilot = loop.Pilot(
            xcpc_pct=85.0 if i >= at else XCPC_PCT[name],
            pas_deg=100.0,
            pcprf_pct=NP_TRIM_RPM * 100.0 / engine_c.NP_DES,
        )
        if s is None:
            s = loop.seed(r, f.wa31_pps, f, pilot, AMB)
        s, _e, h, _fr = loop.step(
            s,
            pilot,
            AMB,
            dt=DT,
            load=lambda v: q + slope * (v - NP_TRIM_RPM),
            j_load=engine_c.J_LOAD_UH60A,
            heat_sink=True,
        )
        out.append(h.wf_pph)
    post = np.asarray(out[at:])
    lo, hi = post[0], post.max()
    i10 = int(np.argmax(post > lo + 0.1 * (hi - lo)))
    i90 = int(np.argmax(post > lo + 0.9 * (hi - lo)))
    rise = (i90 - i10) * DT * 1e3
    assert SLAM_RISE_MS[0] < rise < SLAM_RISE_MS[1], (
        f"the fuel command's 10-90 % rise on a collective slam is {rise:.0f} ms"
    )


def test_the_t45_harness_lag_is_structurally_unreachable():
    """`TLGE` cannot be constrained by anything this suite runs, and the reason is the point.

    It is the lag on the T45 thermocouple harness, which feeds the ECU's **temperature
    limiter**. At every condition the project exercises the ECU is governing *speed* --
    `test_every_printed_trim_governs_on_the_droop_line` asserts exactly that -- so the T45
    path carries a signal nothing downstream acts on, and changing its time constant is
    bit-identical on every output.

    Reaching it needs an operating point where the temperature limiter takes over, which is
    a condition the report never publishes and this project has never had a reason to run.
    Recorded rather than papered over: `TLGE` is carried on the provenance of its Table C.1
    citation alone.
    """
    from t700.control import constants as control_c

    name, wf, _ng, _ps3, q, ratio, _shp = TRIMS[0]

    def settled():
        st, _e, h, fr = settle(wf, q, ratio, name=name)
        return (st.engine.ng_rpm, st.engine.np_rpm, h.wf_pph, fr.t45_degR)

    before = settled()
    original = control_c.TLGE
    try:
        control_c.TLGE = original * 1.5
        after = settled()
    finally:
        control_c.TLGE = original
    assert before == after, (
        f"TLGE now moves the settled state, {before} -> {after}. If the ECU has started "
        f"limiting on temperature somewhere in this suite, this test should become a real "
        f"comparison instead of a record of why it cannot be one."
    )
