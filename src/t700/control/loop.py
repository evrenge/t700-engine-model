"""Engine, ECU and HMU closed on one another -- the whole propulsion system.

Closing the loop makes Wf an output, computed by the control from what the engine does.

    engine state (previous frame)
        |-> ECU  [Figs. C1-C8]  reads NP, TORQ45, T45, W45, P45   -> SPDG
        |-> HMU  [Figs. C9-C22] reads NG, PS3, T2, and SPDG       -> WF
        |-> engine [Eqs. 69-80] reads WF                          -> next state

Both controls read the engine's previous frame, as a real-time implementation does; it is
the same one-frame opening Eq. 74 applies inside the engine.

The report's closed-loop figures (11-15) need the Gen Hel UH-60A simulation, which it does
not contain. The check is Table B.1 [pdf p.67] instead: hold the load at the printed torque
and the loop must settle at the printed state with fuel flow as an output
(`validation/test_closed_loop.py`).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

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

    The governor drives sensed NP to this, so it, not `NP_DES`, sets the trimmed NP."""


@dataclass(frozen=True)
class LoopState:
    """Engine, ECU and HMU states, plus the three engine outputs the ECU senses.

    Carried rather than recomputed because the controls read the previous frame; this
    makes the one-frame delay explicit rather than an accident of evaluation order.
    """

    engine: realtime.RTState
    ecu: ecu_mod.ECUState
    hmu: hmu_mod.HMUState
    t45_degR: float
    w45_pps: float
    q_pt_ftlbf: float


def _ecu_inputs(
    engine: realtime.RTState, t45_degR: float, w45_pps: float, q_pt_ftlbf: float, pilot: Pilot
) -> ecu_mod.ECUInputs:
    """What the ECU senses: NP, TORQ45, T45, W45 and P45 [Figs. C1-C8]."""
    return ecu_mod.ECUInputs(
        np_rpm=engine.np_rpm,
        torq45_ftlbf=q_pt_ftlbf,
        t45_degR=t45_degR,
        w45_pps=w45_pps,
        p45_psia=engine.p45_psia,
        pcprf_pct=pilot.pcprf_pct,
    )


def _hmu_inputs(
    engine: realtime.RTState, pilot: Pilot, ambient: Ambient, spdg: float
) -> hmu_mod.HMUInputs:
    """What the HMU senses: NG, PS3 and T2, plus the pilot's levers and the ECU trim
    [Figs. C9-C22]."""
    return hmu_mod.HMUInputs(
        ng_rpm=engine.ng_rpm,
        ps3_psia=engine_c.K_PS3 * engine.p3_psia,
        t2_degR=ambient.t_amb_degR,
        pas_deg=pilot.pas_deg,
        xcpc_pct=pilot.xcpc_pct,
        spdg=spdg,
    )


def seed(
    trim_result,
    pilot: Pilot,
    ambient: Ambient = STANDARD_DAY,
) -> LoopState:
    """Start the loop from a converged open-loop trim, with both controls consistent.

    `ecu.seed` and `hmu.seed` are our own defaults: the report states no initial condition
    for the block diagrams (open question #55).
    """
    frame = trim_result.frame
    st = realtime.from_trim(trim_result)
    e_in = _ecu_inputs(st, frame.t45_degR, frame.w45_pps, frame.q_pt_ftlbf, pilot)
    h_in = _hmu_inputs(st, pilot, ambient, hmu_mod.SPDG_NULL)
    return LoopState(
        engine=st,
        # Both halves seeded to hold the same torque-motor demand; see `ecu.seed`.
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
) -> tuple[LoopState, ecu_mod.ECUOutputs, hmu_mod.HMUOutputs, realtime.FrameOut]:
    """One frame of the whole propulsion system: ECU, then HMU, then engine.

    `load` is the shaft torque the airframe demands, ft*lbf: a constant or a callable of
    NP. Make it a callable whenever the result matters. A constant means `dQreq/dNP = 0`:
    no aerodynamic damping, so the governor alone holds NP, which is not the system the
    report models (Appendix B's NP diagonal; open question #6).
    """
    e_in = _ecu_inputs(state.engine, state.t45_degR, state.w45_pps, state.q_pt_ftlbf, pilot)
    ecu_state, ecu_out = ecu_mod.step(state.ecu, e_in, dt)

    h_in = _hmu_inputs(state.engine, pilot, ambient, ecu_out.spdg)
    hmu_state, hmu_out = hmu_mod.step(state.hmu, h_in, dt)

    q_req = load(state.engine.np_rpm) if callable(load) else float(load)
    eng_state, frame_out = realtime.step(
        state.engine,
        hmu_out.wf_pps,
        ambient,
        dt=dt,
        q_req_ftlbf=q_req,
        j_load=j_load,
        heat_sink=heat_sink,
    )
    # All six fields are rewritten -- see the note in `ecu.step`.
    new = LoopState(
        engine=eng_state,
        ecu=ecu_state,
        hmu=hmu_state,
        t45_degR=frame_out.t45_degR,
        w45_pps=frame_out.w45_pps,
        q_pt_ftlbf=frame_out.q_pt_ftlbf,
    )
    return new, ecu_out, hmu_out, frame_out
