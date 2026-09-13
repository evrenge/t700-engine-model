"""Engine, ECU and HMU closed on one another -- the whole propulsion system.

Up to here the engine has been driven by a prescribed fuel flow: every result in Phases
0-4 answers "given this Wf, what does the engine do". Closing the loop makes Wf an
**output**, computed by the control from what the engine is doing, and that is a
qualitatively harder thing to get right -- a sign error that merely shifted a number
open-loop will now run away.

## What reads what, and when

    engine state (previous frame)
        |-> ECU  [Figs. C1-C8]  reads NP, TORQ45, T45, W45, P45   -> SPDG
        |-> HMU  [Figs. C9-C22] reads NG, PS3, T2, and SPDG       -> WF
        |-> engine [Eqs. 69-80] reads WF                          -> next state

Both controls read the engine's **previous** frame, which is what a real-time
implementation does and is the same one-frame opening that Eq. 74 applies inside the
engine. It is not an approximation introduced here: the report's own structure is a frame
of sensors, a frame of control, and a frame of engine.

## Why this is the only quantitative check Phase 5 gets

The report's closed-loop figures -- 11 through 15 -- all need the Gen Hel UH-60A
blade-element simulation, which this report consumes and does not contain. So there is no
published closed-loop trace to overlay.

What there is, is Table B.1 [pdf p.67]: the complete engine state at three trims, printed
as numbers, including the power turbine speed and the shaft torque. Hold the load at the
printed torque and set the speed reference, and the loop has to settle at the printed
state **with fuel flow as an output**. Nothing in Appendix C knows about Table B.1 and
nothing in Table B.1 knows about Appendix C, so that is a real cross-check and not a
restatement. `validation/test_closed_loop.py` is where it is made.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace

from t700 import constants as engine_c
from t700 import realtime
from t700.control import ecu as ecu_mod
from t700.control import hmu as hmu_mod
from t700.engine import STANDARD_DAY, Ambient


@dataclass(frozen=True)
class Pilot:
    """The three things a pilot and airframe hand the propulsion system."""

    xcpc_pct: float
    """Collective pitch, percent of maximum [nomenclature pdf p.81]."""
    pas_deg: float
    """Power available spindle -- the power lever."""
    pcprf_pct: float = 100.0
    """Reference power turbine speed, percent of design, as set by the cockpit control.

    The governor drives sensed NP to this, so it -- not `NP_DES` -- is what sets the
    trimmed power turbine speed. Table B.1's 20895 rpm against a design 20900 corresponds
    to 99.976 %, and that is a prediction of the closed loop rather than an input to it."""


@dataclass(frozen=True)
class LoopState:
    """Engine, ECU and HMU states, plus the three engine outputs the ECU senses.

    Those three are carried rather than recomputed because the controls read the previous
    frame; keeping them in the state makes the one-frame delay explicit instead of an
    accident of evaluation order.
    """

    engine: realtime.RTState
    ecu: ecu_mod.ECUState
    hmu: hmu_mod.HMUState
    t45_degR: float
    w45_pps: float
    q_pt_ftlbf: float


def seed(
    trim_result,
    wa31_pps: float,
    frame,
    pilot: Pilot,
    ambient: Ambient = STANDARD_DAY,
) -> LoopState:
    """Start the loop from a converged open-loop trim, with both controls consistent.

    The engine side is `realtime.from_trim`. The control side is `ecu.seed` and
    `hmu.seed`, which are our own defaults -- the report states no initial condition for
    any of the 22 block diagrams (open question #55).
    """
    st = realtime.from_trim(trim_result, wa31_pps, frame)
    e_in = ecu_mod.ECUInputs(
        np_rpm=trim_result.state.np_rpm,
        torq45_ftlbf=frame.q_pt_ftlbf,
        t45_degR=frame.t45_degR,
        w45_pps=frame.w45_pps,
        p45_psia=trim_result.state.p45_psia,
        pcprf_pct=pilot.pcprf_pct,
    )
    h_in = hmu_mod.HMUInputs(
        ng_rpm=trim_result.state.ng_rpm,
        ps3_psia=engine_c.K_PS3 * trim_result.state.p3_psia,
        t2_degR=ambient.t_amb_degR,
        pas_deg=pilot.pas_deg,
        xcpc_pct=pilot.xcpc_pct,
        spdg=hmu_mod.SPDG_NULL,
    )
    return LoopState(
        engine=st,
        # Both halves seeded to hold the SAME torque-motor demand. They disagreed until
        # 2026-09-13 -- the HMU held SPDG_NULL while the ECU delivered 0 on frame 1 -- and
        # the half-volt step that produced cost +31.8 % on fuel flow. See `ecu.seed`.
        ecu=ecu_mod.seed(e_in, spdg=hmu_mod.SPDG_NULL),
        hmu=hmu_mod.seed(h_in),
        t45_degR=frame.t45_degR,
        w45_pps=frame.w45_pps,
        q_pt_ftlbf=frame.q_pt_ftlbf,
    )


def step(
    state: LoopState,
    pilot: Pilot,
    ambient: Ambient = STANDARD_DAY,
    dt: float = realtime.FRAME_ENGINE_S,
    load: float | Callable[[float], float] = 0.0,
    j_load: float = 0.0,
    heat_sink: bool = True,
    integrate_np: bool = True,
    dt_np: float | None = None,
) -> tuple[LoopState, ecu_mod.ECUOutputs, hmu_mod.HMUOutputs, realtime.FrameOut]:
    """One frame of the whole propulsion system: ECU, then HMU, then engine.

    `load` is the shaft torque the airframe demands, in ft*lbf: either a constant or a
    **callable of NP**. Make it a callable whenever the result matters. A constant means
    `dQreq/dNP = 0`, and a free turbine against a speed-independent load has no
    aerodynamic damping at all -- the governor is then the only thing holding NP, which is
    not the system the report models. Appendix B says so: its NP diagonal carries a ~52 %
    residual against ours precisely because `dQreq/dNP` is a Gen Hel quantity we do not
    have, and open question #6 recovers it as roughly 1.78 / 1.45 / 1.18 times
    `|dQ_PT/dNP|` at the three trims -- positive, because rotor torque rises with speed.

    `dt_np` and `integrate_np` pass the report's 2:1 multirate through to the engine
    [pdf p.47]: a caller driving this at 7 ms who wants the shipped configuration calls it
    with `dt_np=realtime.FRAME_NP_S` on even frames and `integrate_np=False` on odd ones.
    It cannot move a settled result -- at equilibrium `dNP/dt` is zero whatever step it is
    multiplied by -- so every closed-loop comparison against Table B.1 here is
    rate-independent by construction, and only the path there responds.

    **Both arguments were documented here and did not exist** between two commits on
    2026-09-13: the patch that was meant to add them raised before writing the file, the
    follow-up fixed only this docstring, and the audit ledger recorded a fix that had not
    been made. `grep -rn dt_np` over `control/` returning nothing is what caught it. They
    exist now.
    """
    e_in = ecu_mod.ECUInputs(
        np_rpm=state.engine.np_rpm,
        torq45_ftlbf=state.q_pt_ftlbf,
        t45_degR=state.t45_degR,
        w45_pps=state.w45_pps,
        p45_psia=state.engine.p45_psia,
        pcprf_pct=pilot.pcprf_pct,
    )
    ecu_state, ecu_out = ecu_mod.step(state.ecu, e_in, dt)

    h_in = hmu_mod.HMUInputs(
        ng_rpm=state.engine.ng_rpm,
        ps3_psia=engine_c.K_PS3 * state.engine.p3_psia,
        t2_degR=ambient.t_amb_degR,
        pas_deg=pilot.pas_deg,
        xcpc_pct=pilot.xcpc_pct,
        spdg=ecu_out.spdg,
    )
    hmu_state, hmu_out = hmu_mod.step(state.hmu, h_in, dt)

    q_req = load(state.engine.np_rpm) if callable(load) else float(load)
    eng_state, frame_out = realtime.step(
        state.engine,
        hmu_out.wf_pps,
        ambient,
        dt=dt,
        q_req_ftlbf=q_req,
        j_load=j_load,
        integrate_np=integrate_np,
        heat_sink=heat_sink,
        dt_np=dt_np,
    )
    new = replace(
        state,
        engine=eng_state,
        ecu=ecu_state,
        hmu=hmu_state,
        t45_degR=frame_out.t45_degR,
        w45_pps=frame_out.w45_pps,
        q_pt_ftlbf=frame_out.q_pt_ftlbf,
    )
    return new, ecu_out, hmu_out, frame_out
