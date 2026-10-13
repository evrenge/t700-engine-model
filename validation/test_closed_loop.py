"""Engine, ECU and HMU closed on one another, against Table B.1 [pdf p.67].

The report's closed-loop Figures 11-15 all need Gen Hel, so the Phase 5 check is this: hold
the load at Table B.1's printed torque, set the speed reference, and the loop must settle at
the printed state with fuel flow as an output. The engine (body, Appendix A) and the control
(Appendix C) are independent halves, so this is a cross-check, not a restatement.

The load is not speed-independent: `dQreq/dNP = 0` would leave the governor alone holding NP.
The slope `DQ_REQ_DNP` and the load inertia `J_LOAD_UH60A` are the ones the linear models use,
both recovered from Appendix B (open question #6): the airframe's, cited, not invented.
"""

from __future__ import annotations

import functools

import numpy as np
import pytest

from t700 import constants as engine_c
from t700 import maps, trim
from t700.control import loop
from t700.engine import Ambient
from t700.units import shp_from_torque, wf_pps_from_pph
from test_all_dof_models import DQ_REQ_DNP
from test_trim_points import CASES, STATE_B1

AMB = Ambient(14.696, 518.67)
DT = 0.007
"""The report's engine frame [pdf p.47]."""

SETTLE_FRAMES = 4000
"""28 s. The slowest thing in the loop is the ECU's P+I integrator against the rotor
inertia, settled by about 20 s."""

NP_TRIM_RPM = 20895.0
"""Table B.1's power turbine speed, printed at all three trims [pdf p.67]."""

PCPRF_PCT = NP_TRIM_RPM * 100.0 / engine_c.NP_DES
"""The cockpit speed reference that holds NP at `NP_TRIM_RPM`, percent of NP_des."""

CLOSED_LOOP_TOL_PCT = 0.5
"""Ours, declared in `SCOPE.md`: the Phase 5 gate. Applied to the cycle mean, never to a
terminal sample, which is one phase of the limit cycle (open question #61). The worst,
descent fuel flow, is partly inherited: the loop trims fuel until the torque matches, so the
open-loop descent shaft-power error shows up in fuel flow instead."""

# name, Wf lbm/hr, NG rpm, Ps3 psia, shaft torque ft*lbf, dQreq/dNP ft*lbf/rpm, shp --
# Table B.1 from `test_trim_points`, the load slope from Appendix B (open question #6)
TRIMS = [
    (b.name, b.wf_pph, b.ng_rpm, STATE_B1[b.name]["ps3"], b.q_shaft_ftlbf, DQ_REQ_DNP[n], b.shp)
    for n, b in enumerate(CASES, start=1)
]

XCPC_PCT = {"hover": 52.75, "level 80 kt": 35.87, "descent 80 kt": 23.01}
"""Collective at each trim, from `test_hmu_trim.py`'s own droop-line solve.

The settled state barely depends on it (the governor absorbs it), but the path does: at the
hover collective, the level and descent runs spike fuel flow and leave `f1`'s data on the way
there."""


CYCLE_FRAMES = 1500
"""How many trailing frames to average a closed-loop result over: 10.5 s at 7 ms.

The loop settles to a limit cycle, not a point, so a terminal sample is one phase of it
(open question #61, traced to `CH`, the 0.05 %NG hysteresis on the gas-generator speed
sensor [Fig. C15, pdf p.91]). 1500 frames covers at least two full periods at every trim.
"""


def mean_over_cycle(values) -> float:
    """Mean of a trailing window that spans a whole number of limit-cycle periods."""
    return float(np.mean(np.asarray(values)[-CYCLE_FRAMES:]))


def closed_loop(wf_pph, q_shaft, dq_req_dnp, n, pilot_of, load_scale_of=None):
    """Close the loop from the open-loop trim at `wf_pph` and step it `n` frames.

    The one runner for every closed-loop scenario here and in `report.py`. The load is
    `q_shaft * load_scale_of(i) + dq_req_dnp * (NP - NP_TRIM_RPM)`, the slope the linear
    models use (open question #6); `pilot_of(i)` is the pilot on frame `i`,
    and the loop is seeded with frame 0's. Yields `(i, state, ecu_out, hmu_out, frame_out)`
    after each frame, so each caller records what it needs.
    """
    wf = wf_pps_from_pph(wf_pph)
    r = trim.solve(wf, engine_c.NP_DES, AMB)
    s = None
    for i in range(n):
        pilot = pilot_of(i)
        if s is None:
            s = loop.seed(r, pilot, AMB)
        k = 1.0 if load_scale_of is None else load_scale_of(i)
        s, ecu_o, hmu_o, fr = loop.step(
            s,
            pilot,
            AMB,
            dt=DT,
            load=lambda v, k=k: q_shaft * k + dq_req_dnp * (v - NP_TRIM_RPM),
            j_load=engine_c.J_LOAD_UH60A,
            heat_sink=True,
        )
        yield i, s, ecu_o, hmu_o, fr


def _settle_run(wf_pph, q_shaft, dq_req_dnp, xcpc_pct, pcprf_pct, n):
    """One held-pilot run: the final `(state, ecu_out, hmu_out, frame_out)` and the history
    of the five quantities Table B.1 prints."""
    pilot = loop.Pilot(xcpc_pct=xcpc_pct, pas_deg=100.0, pcprf_pct=pcprf_pct)
    hist: dict[str, list[float]] = {k: [] for k in ("NG", "NP", "Wf", "Ps3", "shp")}
    for out in closed_loop(wf_pph, q_shaft, dq_req_dnp, n, lambda _i: pilot):
        _i, st, _e, h, fr = out
        hist["NG"].append(st.engine.ng_rpm)
        hist["NP"].append(st.engine.np_rpm)
        hist["Wf"].append(h.wf_pph)
        hist["Ps3"].append(engine_c.K_PS3 * st.engine.p3_psia)
        hist["shp"].append(shp_from_torque(fr.q_pt_ftlbf, st.engine.np_rpm))
    return out[1:], hist


_settled = functools.cache(_settle_run)
"""`_settle_run`, run once per pilot and shared. Nothing that reads it mutates it."""


def settle(
    wf_pph, q_shaft, dq_req_dnp, xcpc_pct=None, n=SETTLE_FRAMES, pcprf=None, name=None, fresh=False
):
    """Close the loop from the open-loop trim and run it to steady state.

    `xcpc_pct` defaults to the collective that trim actually needs -- see `XCPC_PCT`.

    Returns the FINAL state. Use `settle_cycle` where a number is being compared against
    the report: the final state is one phase of a limit cycle (see `CYCLE_FRAMES`).

    Cached per pilot; `fresh=True` runs it anew, for a caller that counts the run's clamps
    or patches a constant around it.
    """
    if xcpc_pct is None:
        xcpc_pct = XCPC_PCT[name] if name in XCPC_PCT else 52.75
    run = _settle_run if fresh else _settled
    return run(wf_pph, q_shaft, dq_req_dnp, xcpc_pct, PCPRF_PCT if pcprf is None else pcprf, n)[0]


def settle_cycle(wf_pph, q_shaft, dq_req_dnp, name=None, n=SETTLE_FRAMES):
    """As `settle`, but returning the mean over the last whole cycles of each output.


    Also returns the peak-to-peak of each, so a caller can say how much of a deviation is
    the cycle and how much is the model.
    """
    _final, hist = _settled(wf_pph, q_shaft, dq_req_dnp, XCPC_PCT.get(name, 52.75), PCPRF_PCT, n)
    mean = {k: mean_over_cycle(v) for k, v in hist.items()}
    ptp = {k: float(np.ptp(np.asarray(v)[-CYCLE_FRAMES:])) for k, v in hist.items()}
    return mean, ptp


@pytest.mark.parametrize("name,wf,ng,ps3,q,dq_req_dnp,shp", TRIMS, ids=lambda v: str(v))
def test_the_closed_loop_governs_to_the_printed_trim(name, wf, ng, ps3, q, dq_req_dnp, shp):
    """Fuel flow is an output here: nothing tells the loop what it should be. Compared as a
    cycle mean (open question #61, `CYCLE_FRAMES`).
    """
    got, ptp = settle_cycle(wf, q, dq_req_dnp, name=name)
    want = {"NG": ng, "NP": NP_TRIM_RPM, "Wf": wf, "Ps3": ps3, "shp": shp}
    for key, w in want.items():
        dev = 100.0 * (got[key] / w - 1.0)
        assert abs(dev) < CLOSED_LOOP_TOL_PCT, (
            f"{name} {key}: closed loop governs to {got[key]:.2f} against Table B.1's "
            f"{w}, {dev:+.3f} % (cycle peak-to-peak {ptp[key]:.3f})"
        )


@pytest.mark.parametrize("name,wf,ng,ps3,q,dq_req_dnp,shp", TRIMS, ids=lambda v: str(v))
def test_every_printed_trim_governs_on_the_droop_line(name, wf, ng, ps3, q, dq_req_dnp, shp):
    """Closed loop, as open loop: a governor in steady flight sits on its droop line, not
    on a cam. And the ECU should be governing speed, not limiting temperature."""
    _s, e, h, _fr = settle(wf, q, dq_req_dnp, name=name)
    assert h.limit == "droop", (name, h.limit)
    assert e.limiting == "speed", (name, e.limiting, e.tsig, e.spdsf)


def test_the_governor_absorbs_the_collective():
    """A collective is a pilot input and must not set the trim: the load does that, and the
    collective only moves where the ECU's integrator sits. Sweep it nearly twofold and the
    settled engine barely notices, while `SPDG` travels.
    """
    _name, wf, ng, _ps3, q, dq_req_dnp, _shp = TRIMS[0]
    seen = []
    for xcpc in (35.0, 45.0, 52.75, 60.0, 70.0):
        s, e, h, _fr = settle(wf, q, dq_req_dnp, xcpc_pct=xcpc)
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
    cockpit reference and not `NP_DES`.
    """
    _name, wf, _ng, _ps3, q, dq_req_dnp, _shp = TRIMS[0]
    for pcprf, expect in ((99.0, 0.99 * engine_c.NP_DES), (100.0, engine_c.NP_DES)):
        s, _e, _h, _fr = settle(wf, q, dq_req_dnp, pcprf=pcprf)
        dev = 100.0 * (s.engine.np_rpm / expect - 1.0)
        assert abs(dev) < 0.1, (
            f"with PCPRF = {pcprf} % the loop should hold NP at {expect:.0f} rpm, "
            f"got {s.engine.np_rpm:.0f} ({dev:+.3f} %)"
        )


def test_the_acceleration_limit_keeps_the_engine_inside_its_digitized_envelope():
    """Open question #45: an open-loop fuel step leaves the maps, but the control would never
    command one. With the loop closed, a harder demand than Figure 9's (collective 52.75 ->
    95 %, power lever 120 deg) leaves `f1`, `f7` and `f9` untouched: `F_HM7`, the
    acceleration limit, binds and caps the fuel.
    """
    _name, wf, _ng, _ps3, q, dq_req_dnp, _shp = TRIMS[0]
    cruise = loop.Pilot(xcpc_pct=52.75, pas_deg=100.0, pcprf_pct=PCPRF_PCT)
    slam = loop.Pilot(xcpc_pct=95.0, pas_deg=120.0, pcprf_pct=PCPRF_PCT)
    limits = {}
    # 300 frames at cruise, then the slam; clamps are counted from the slam on
    for i, _s, _e, h, _fr in closed_loop(
        wf, q, dq_req_dnp, 300 + 1200, lambda i: cruise if i < 300 else slam
    ):
        if i == 299:
            maps.reset_clamps()
        elif i >= 300:
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


@pytest.mark.parametrize("name,wf,ng,ps3,q,dq_req_dnp,shp", TRIMS, ids=lambda v: str(v))
def test_the_closed_loop_stays_inside_its_digitized_data(name, wf, ng, ps3, q, dq_req_dnp, shp):
    """The Phase 5 gate must not reach its answer through extrapolated maps; each trim runs
    at its own collective (see `XCPC_PCT`).
    """
    maps.reset_clamps()
    settle(wf, q, dq_req_dnp, name=name, fresh=True)  # the clamps counted are this run's
    report = maps.clamp_report()
    maps.reset_clamps()

    assert not report, (
        f"{name}: the closed loop reads {report} outside their digitized range on the "
        f"way to its settled state. Every number this file reports would then be resting "
        f"on extrapolated maps."
    )


def test_the_limit_cycle_is_what_it_is_measured_to_be():
    """Pin the limit cycle, so it cannot grow or vanish unnoticed.

    It is the report's control system as specified (open question #61): `CH`, the speed
    sensor's 0.05 %NG hysteresis [Fig. C15, pdf p.91], is finer than the governor can
    resolve, and the cycle survives step refinement and collapses when `CH` alone is removed.
    """
    want = {"hover": (13.0, 25.8), "level 80 kt": (12.7, 26.5), "descent 80 kt": (16.3, 39.3)}
    for name, wf, _ng, _ps3, q, dq_req_dnp, _shp in TRIMS:
        _mean, ptp = settle_cycle(wf, q, dq_req_dnp, name=name)
        np_pp, ng_pp = want[name]
        assert ptp["NP"] == pytest.approx(np_pp, abs=1.0), f"{name} NP p-p {ptp['NP']:.3f}"
        assert ptp["NG"] == pytest.approx(ng_pp, abs=2.0), f"{name} NG p-p {ptp['NG']:.3f}"
        # and it must remain small enough that the mean is the honest summary
        assert 100.0 * ptp["NG"] / _ng < 0.15, f"{name} NG cycle is {ptp['NG']:.1f} rpm"


# ---------------------------------------------- Appendix C's forward path, which nothing saw

LOAD_STEP_PCT = 15.0
NP_DIP_RPM = (209.4, 211.4)
WF_PEAK_PPH = (559.5, 567.5)
SLAM_RISE_MS = (140.0, 155.0)
"""Characterizations of the loop's transient, not accuracy claims (SCOPE.md's third
category): the report prints no closed-loop time history (open question #65).

Every steady test either solves for a free input or integrates the error away, so no
forward-path gain in Appendix C shows there; a transient is where a loop gain shows. The
load step pins `TMLG` (and `TMGN` only coarsely), the collective slam pins `CLLDS`. The bands
are about 1 % wide around the measured values: tripwires, not tolerances. `TLGE` is pinned in
`tests/test_ecu.py` instead (see `test_the_t45_harness_lag_is_invisible_to_the_closed_loop`)."""


def test_a_load_step_response_pins_the_torque_motor_linkage():
    """A 15 % load step at hover, characterized: how far NP is dragged down before the fuel
    catches up, how much fuel it takes, and that NP comes back.
    """
    name, wf, _ng, _ps3, q, dq_req_dnp, _shp = TRIMS[0]
    pilot = loop.Pilot(xcpc_pct=XCPC_PCT[name], pas_deg=100.0, pcprf_pct=PCPRF_PCT)
    # The dip and the peak are read over the 1200 frames after the step that NP_DIP_RPM and
    # WF_PEAK_PPH were measured on; the recovery, 3000 frames after it.
    step_at, window_end, n = 1400, 2600, 1400 + 3000
    np_min, wf_max = 1e9, 0.0
    np_before = None
    for i, s, _e, h, _fr in closed_loop(
        wf,
        q,
        dq_req_dnp,
        n,
        lambda i: pilot,
        lambda i: (1.0 + LOAD_STEP_PCT / 100.0) if i >= step_at else 1.0,
    ):
        if i == step_at - 1:
            np_before = s.engine.np_rpm
        if step_at <= i < window_end:
            np_min = min(np_min, s.engine.np_rpm)
            wf_max = max(wf_max, h.wf_pph)
    dip = NP_TRIM_RPM - np_min
    assert NP_DIP_RPM[0] < dip < NP_DIP_RPM[1], f"NP dips {dip:.1f} rpm on a 15 % load step"
    assert WF_PEAK_PPH[0] < wf_max < WF_PEAK_PPH[1], f"Wf peaks at {wf_max:.1f} pph"

    assert np_min < np_before, "a load step should droop the power turbine before recovery"
    dev = 100.0 * (s.engine.np_rpm / np_before - 1.0)
    assert abs(dev) < 0.5, f"NP did not recover after a 15 % load step: {dev:+.3f} %"
    assert h.wf_pph > 1.05 * wf, (
        f"a 15 % load step should raise fuel flow well above {wf}, got {h.wf_pph:.1f}"
    )


def test_a_collective_slam_pins_the_collective_lag():
    """`CLLDS` is the load-demand-spindle lag, and only a moving collective can see it: slam
    the collective from the hover setting to 85 % and time the fuel command's 10-90 % rise.
    """
    name, wf, _ng, _ps3, q, dq_req_dnp, _shp = TRIMS[0]
    out, at, n = [], 600, 1800
    for _i, _s, _e, h, _fr in closed_loop(
        wf,
        q,
        dq_req_dnp,
        n,
        lambda i: loop.Pilot(
            xcpc_pct=85.0 if i >= at else XCPC_PCT[name], pas_deg=100.0, pcprf_pct=PCPRF_PCT
        ),
    ):
        out.append(h.wf_pph)
    post = np.asarray(out[at:])
    lo, hi = post[0], post.max()
    i10 = int(np.argmax(post > lo + 0.1 * (hi - lo)))
    i90 = int(np.argmax(post > lo + 0.9 * (hi - lo)))
    rise = (i90 - i10) * DT * 1e3
    assert SLAM_RISE_MS[0] < rise < SLAM_RISE_MS[1], (
        f"the fuel command's 10-90 % rise on a collective slam is {rise:.0f} ms"
    )


def test_the_t45_harness_lag_is_invisible_to_the_closed_loop():
    """`TLGE`, the T45 harness lag, cannot be constrained by anything run against the engine.

    It feeds the ECU's temperature limiter, which needs an NG past the report's 100 % ceiling
    to reach, so at every trim the ECU governs speed and `TLGE` is bit-identical on every
    output. The ECU alone can reach it:
    `tests/test_ecu.py::test_the_thermocouple_lag_shapes_the_overtemperature_handover`.
    """
    from t700.control import constants as control_c

    name, wf, _ng, _ps3, q, dq_req_dnp, _shp = TRIMS[0]

    def settled():
        st, _e, h, fr = settle(wf, q, dq_req_dnp, name=name, fresh=True)  # under a patched TLGE
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
        f"limiting on temperature at a Table B.1 trim, that is a change in the engine, not "
        f"in the control, and it should be explained before this test is relaxed."
    )


GENTLE_SLAM_TO_PCT = 70.0
SLAM_NP_RISE_RPM = (536.5, 548.5)
SLAM_PEAK_PPH = (678.5, 685.5)
"""A collective slam small enough that no limit binds -- which is what makes `F_HM4` visible.

`F_HM4` is the load-demand fuel-flow feedforward, `WFQPS3 = F_HM4(XLDSH)` [Fig. C17], on the
droop line's forward path (#65). The 85 % slam above puts the Figure C18 cascade on its
acceleration ceiling, where `MIN(., WFPAC)` discards the feedforward; stopping at 70 % keeps
it on `droop`. Characterizations, as above; the NP rise is the better discriminator because
it is monotone in a change to `F_HM4`.

The load is a function of NP only, so this exercises the feedforward path alone, not a
manoeuvre: NP rises rather than dips. Ballin's load-demand compensation is an admitted limit,
"weaker than the real aircraft" [pdf p.54]; this measures the path, not its adequacy."""


def test_a_gentle_collective_slam_pins_the_load_demand_feedforward():
    """`F_HM4` moves the transient only while the Figure C18 cascade stays on the droop
    branch, so the test asserts that it does before asserting anything else."""
    name, wf, _ng, _ps3, q, dq_req_dnp, _shp = TRIMS[0]
    wf_hist, np_hist, limits, at, n = [], [], set(), 600, 1800
    for _i, s, _e, h, _fr in closed_loop(
        wf,
        q,
        dq_req_dnp,
        n,
        lambda i: loop.Pilot(
            xcpc_pct=GENTLE_SLAM_TO_PCT if i >= at else XCPC_PCT[name],
            pas_deg=100.0,
            pcprf_pct=PCPRF_PCT,
        ),
    ):
        wf_hist.append(h.wf_pph)
        np_hist.append(s.engine.np_rpm)
        limits.add(h.limit)

    assert limits == {"droop"}, (
        f"a limit bound during the gentle slam ({sorted(limits)}), which is what hides "
        f"the load-demand feedforward -- pick a smaller slam"
    )
    peak = max(wf_hist[at:])
    rise = max(np_hist[at:]) - NP_TRIM_RPM
    assert SLAM_PEAK_PPH[0] < peak < SLAM_PEAK_PPH[1], f"Wf peaks at {peak:.2f} pph"
    assert SLAM_NP_RISE_RPM[0] < rise < SLAM_NP_RISE_RPM[1], f"NP rises {rise:.1f} rpm"
