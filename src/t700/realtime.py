"""The real-time frame: Eqs. 69-80, the model Ballin actually ran.

`engine.frame` evaluates the physics with the three volume pressures as **states**
(Eqs. 42-44). That is the 5-DOF model Appendix B linearizes, and it is what the trim
solver uses. **It cannot be integrated at the report's frame time.**

Our own Jacobian at the hover trim puts the fastest mode at 4881/sec -- a time constant of
0.205 ms -- so explicit integration of five states is stable only below about 0.41 ms. The
report runs the engine at **7 ms**, seventeen times over that limit. Not a preference: it
would diverge on the first frame.

So Ballin makes the three pressures **algebraic** and solves them each frame:

* the **opened compressor mass-flow iteration** (Eq. 74) breaks the expensive outer loop
  through the compressor map by carrying one mass flow across the frame boundary;
* the **inner P3/P41 loop** (Eqs. 76, 78) is a Gauss-Seidel sweep -- "eleven arithmetic
  operations for each pass ... no relaxation algorithm is required for convergence"
  [pdf p.37];
* **P45 is solved independently** (Eq. 80), "because it is a nonlinear function of its own
  value, iterative techniques must be used" [pdf p.37].

What remains integrated is NG (tau = 18 ms) and NP (tau ~ 340 ms), both comfortably stable
at the 7 and 14 ms the report uses.

This also decodes the Conclusions [pdf p.54]: omitting the high-speed inter-volume
mass-flow dynamics "was found to be unnecessary". Those dynamics *are* the two fast modes.
They are far faster than anything a pilot or rotor can feel, and carrying them would have
cost a twentyfold smaller time step.

## Decisions recorded here, not transcribed

**The integrator is never named in the report** (open question #26). Explicit Euler is used
here: it is what a 1988 real-time model with "no iteration between time steps" [pdf p.35]
would use, and the two surviving modes are slow enough for it. Logged as a decision.

**Eq. 74 contradicts its own prose** (open question #22). As printed it lags only the
bleed, `WA31_(n) = WA3_(n) - WA3_bl_(n-1)` -- but `WA3_(n)` still depends on `P3_(n)`, so
that does not open the loop at all. The prose says the combustor inlet flow is "the mass
flow leaving the compressor in the previous interval". Only the prose reading actually
opens it. Both are implemented, selectable, so the report's own 0.1 %-per-frame criterion
can decide between them -- see `tools/` and the transient validation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from t700 import constants as c
from t700 import corrections, maps, thermo
from t700.engine import STANDARD_DAY, Ambient, State
from t700.units import BTU_TO_FTLBF, RAD_PER_SEC_TO_RPM, RPM_TO_RAD_PER_SEC

_TORQUE_SCALE = BTU_TO_FTLBF * RAD_PER_SEC_TO_RPM

FRAME_ENGINE_S = 0.007
"""Engine update interval [pdf p.47]: "updated twice for each rotor routine cycle, or once
every 7 msec"."""

FRAME_NP_S = 0.014
"""NP update interval [pdf p.47]. NP is the drive-train degree of freedom and moves at the
host rotor rate, so the model is multirate 2:1."""

MAX_STEP_S = 0.010
"""Maximum time step [pdf p.38], set by a 0.1 % maximum allowable error between steps.

Note this is *smaller* than FRAME_NP_S. The report states both; open question #33."""


@dataclass(frozen=True)
class RTState:
    """Five model states plus the one scalar the opened iteration carries across frames."""

    ng_rpm: float
    np_rpm: float
    p3_psia: float
    p41_psia: float
    p45_psia: float
    wa31_carry_pps: float
    """The opened compressor mass-flow iteration's memory (Eq. 74). Its initialization is
    not stated anywhere in the report -- open question #23."""

    def to_state(self) -> State:
        return State(self.ng_rpm, self.np_rpm, self.p3_psia, self.p41_psia, self.p45_psia)


@dataclass(frozen=True)
class FrameOut:
    """What one real-time frame produced, beyond the next state."""

    inner_iters: int
    p45_iters: int
    t41_degR: float
    t45_degR: float
    q_pt_ftlbf: float
    q_gt_ftlbf: float
    q_c_ftlbf: float
    wa2_pps: float
    w41_pps: float
    w45_pps: float
    far: float


def _inner_pressure_loop(
    p3: float,
    p41: float,
    t3: float,
    theta41: float,
    wa31: float,
    wf: float,
    tol: float = 1e-10,
    max_iter: int = 40,
) -> tuple[float, float, int]:
    """The P3/P41 fixed-point sweep, Eqs. 76 and 78.

    Gauss-Seidel exactly as the report describes it [pdf p.37]: P3 from the current P41,
    then P41 from the P3 just computed -- not a simultaneous solve. T3 and theta_41 are
    held fixed across the sweep, which is what makes it eleven arithmetic operations a
    pass; they are recomputed once per frame outside this function.
    """
    it = 0
    while it < max_iter:
        it += 1
        p3_new = 0.5 * (p41 + np.sqrt(p41 * p41 + 4.0 * c.K_DPB * t3 * wa31 * wa31))  # (76)
        w41 = c.K_WGT * p41 / np.sqrt(theta41)
        p41_new = p3_new - t3 * c.K_DPB * (w41 - wf) ** 2 / p3_new  # (78) with (77)
        if abs(p3_new - p3) < tol * p3 and abs(p41_new - p41) < tol * p41:
            p3, p41 = p3_new, p41_new
            break
        p3, p41 = p3_new, p41_new
    return float(p3), float(p41), it


def _p45_loop(
    p45: float,
    w41: float,
    b3: float,
    wa2: float,
    theta45: float,
    ps9: float,
    tol: float = 1e-10,
    max_iter: int = 40,
) -> tuple[float, int]:
    """P45 by its own iteration, Eq. 80.

    The numerator is constant over the iteration, as the report notes -- only the f9
    lookup changes, which is why a pass costs "one function-table lookup and four
    arithmetic operations" [pdf p.37].

    Eq. 79 prints the returned bleed as `B3 B4 WA2` where Eq. 44 writes `B3 K_bl WA2` for
    the identical term, and B4 is in no nomenclature (open question #30). K_bl is used.
    """
    numerator = (w41 + b3 * c.K_BL * wa2) * np.sqrt(theta45)
    it = 0
    while it < max_iter:
        it += 1
        p45_new = numerator / float(maps.f9()(ps9 / p45))
        if abs(p45_new - p45) < tol * p45:
            p45 = p45_new
            break
        p45 = p45_new
    return float(p45), it


def step(
    st: RTState,
    wf_pps: float,
    ambient: Ambient = STANDARD_DAY,
    dt: float = FRAME_ENGINE_S,
    q_req_ftlbf: float = 0.0,
    j_load: float = 0.0,
    integrate_np: bool = True,
    lag_whole_flow: bool = True,
) -> tuple[RTState, FrameOut]:
    """Advance one engine frame.

    Args:
        lag_whole_flow: how to read Eq. 74. True (default) follows the report's **prose** --
            the combustor inlet flow is the compressor exit flow of the previous interval,
            which is what actually opens the loop. False follows the equation **as
            printed**, lagging only the bleed. Open question #22; see the module docstring.
        integrate_np: False suppresses the NP integration, which is what the report does
            for the open-loop fuel steps of Figs. 9 and 10 [pdf p.39].
    """
    p2 = ambient.p_amb_psia
    t2 = ambient.t_amb_degR
    h2 = thermo.h2_from_t2(t2)  # (3)
    theta2 = corrections.theta2(t2)  # (5)
    ngc = corrections.corrected_speed(st.ng_rpm, theta2)  # (6)
    ngc_pct = 100.0 * ngc / c.NG_DES
    delta2 = corrections.delta2(p2)  # (8)

    # --- the compressor, evaluated ONCE on the frame's entering pressure ---------------
    ps3 = c.K_PS3 * st.p3_psia  # (4)
    pr = ps3 / p2
    wa2c = float(maps.f1()(pr, ngc_pct))  # (7)
    wa2 = corrections.physical_flow(wa2c, theta2, delta2)  # (9)
    t3 = t2 * float(maps.f2()(pr))  # (10)
    h3 = thermo.h3_from_t3(t3)  # (11)

    b1 = float(maps.f3()(ngc_pct))  # (12)
    b2 = float(maps.f4()(wa2c))  # (13)
    b3 = float(maps.f5()(wa2c))  # (14)
    wa3_bl = wa2 * (b3 + c.K_B3)  # (16)
    wa3 = wa2 - wa2 * (b1 + b2)  # (15), (17)

    # --- Eq. 74, the opened iteration -------------------------------------------------
    wa31 = st.wa31_carry_pps if lag_whole_flow else (wa3 - st.wa31_carry_pps)
    if wa31 <= 0.0:
        wa31 = max(wa3 - wa3_bl, 1e-6)

    # --- combustor, once per frame ----------------------------------------------------
    far = wf_pps / wa31  # (19)
    eta_b = float(maps.f6()(far))  # (20)
    h41_ns = (h3 + eta_b * far * c.HVF) / (1.0 + far)  # (21)
    t41 = thermo.t41_from_h41(h41_ns)  # (22), (23) with no heat sink
    h41 = thermo.h41_from_t41(t41)  # (24)
    theta41 = thermo.theta41_from_t41(t41)  # (25)

    # --- the two pressure solves ------------------------------------------------------
    p3, p41, inner_iters = _inner_pressure_loop(st.p3_psia, st.p41_psia, t3, theta41, wa31, wf_pps)
    w41 = c.K_WGT * p41 / np.sqrt(theta41)  # (28)

    ps9 = p2  # (37)
    p49 = ps9 * float(maps.f10()(ngc_pct))  # (38)

    dh_gt = theta41 * float(maps.f7()(st.p45_psia / p41))  # (26)
    h44 = h41 - dh_gt  # (27)
    h45 = thermo.h45_from_h44(h44)  # (29)
    t45 = thermo.t45_from_h45(h45)  # (30)
    theta45 = thermo.theta45_from_t45(t45)  # (31)

    p45, p45_iters = _p45_loop(st.p45_psia, w41, b3, wa2, theta45, ps9)

    dh_pt = theta45 * float(maps.f8()(p49 / p45))  # (32)
    w45 = float(maps.f9()(ps9 / p45)) * p45 / np.sqrt(theta45)  # (33), (34)

    # --- torques and the two surviving integrations -----------------------------------
    q_c = _TORQUE_SCALE / st.ng_rpm * (wa2 * (c.K_QC_1 * h3 - h2) + wa3 * c.K_QC_2 * h3)  # (39)
    q_gt = _TORQUE_SCALE / st.ng_rpm * w41 * dh_gt  # (40)
    q_pt = _TORQUE_SCALE / st.np_rpm * w45 * dh_pt - c.K_DAMP * RPM_TO_RAD_PER_SEC * (
        st.np_rpm - c.NP_DES
    )  # (41)

    ng = st.ng_rpm + dt * RAD_PER_SEC_TO_RPM * (q_gt - q_c) / c.J_GT  # (45)
    np_ = st.np_rpm
    if integrate_np:
        np_ = st.np_rpm + dt * RAD_PER_SEC_TO_RPM * (q_pt - q_req_ftlbf) / (
            c.J_PT + j_load
        )  # (46), (47)

    carry = wa31 if not lag_whole_flow else (wa3 - wa3_bl)
    nxt = RTState(ng, np_, p3, p41, p45, float(carry))
    out = FrameOut(
        inner_iters=inner_iters,
        p45_iters=p45_iters,
        t41_degR=t41,
        t45_degR=t45,
        q_pt_ftlbf=q_pt,
        q_gt_ftlbf=q_gt,
        q_c_ftlbf=q_c,
        wa2_pps=wa2,
        w41_pps=w41,
        w45_pps=w45,
        far=far,
    )
    return nxt, out


def from_trim(result, wa31_pps: float) -> RTState:
    """Seed a real-time state from a converged trim.

    The carried mass flow is initialized to its equilibrium value. The report never says
    how it initializes (open question #23); starting at equilibrium is the choice that
    makes a trimmed engine sit still, which is the only defensible default.
    """
    s = result.state
    return RTState(s.ng_rpm, s.np_rpm, s.p3_psia, s.p41_psia, s.p45_psia, wa31_pps)


def run(
    st: RTState,
    wf_of_t,
    ambient: Ambient = STANDARD_DAY,
    duration_s: float = 2.0,
    dt: float = FRAME_ENGINE_S,
    **kw,
) -> dict[str, np.ndarray]:
    """Integrate for `duration_s`, returning traces keyed by name.

    `wf_of_t` is a callable of time in seconds returning fuel flow in lbm/sec, so a step
    input is just a lambda.
    """
    n = int(round(duration_s / dt)) + 1
    keys = ("t", "ng", "np", "p3", "p41", "p45", "t41", "t45", "q_pt", "wf", "wa2", "far")
    tr = {k: np.zeros(n) for k in keys}
    for i in range(n):
        t = i * dt
        wf = float(wf_of_t(t))
        tr["t"][i] = t
        tr["ng"][i] = st.ng_rpm
        tr["np"][i] = st.np_rpm
        tr["p3"][i] = st.p3_psia
        tr["p41"][i] = st.p41_psia
        tr["p45"][i] = st.p45_psia
        tr["wf"][i] = wf
        st, out = step(st, wf, ambient, dt, **kw)
        tr["t41"][i] = out.t41_degR
        tr["t45"][i] = out.t45_degR
        tr["q_pt"][i] = out.q_pt_ftlbf
        tr["wa2"][i] = out.wa2_pps
        tr["far"][i] = out.far
    return tr


__all__ = ["RTState", "FrameOut", "step", "run", "from_trim", "FRAME_ENGINE_S", "MAX_STEP_S"]
