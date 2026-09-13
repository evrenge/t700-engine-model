"""One frame of the T700 gas path: Eqs. 1-49, in the report's own order.

This is the heart of the replication. Given the five states and a fuel flow, `frame()`
walks the gas path once -- inlet, compressor, bleeds, combustor, gas generator turbine,
station 4.5 mixing, power turbine, exhaust -- and returns every station value, the three
torques, and the five state derivatives.

## Why one forward pass is enough

Ballin's three pressures are **pure integrators** (Eqs. 42-44): each is a volume constant
times the integral of its own temperature times a net mass-flow imbalance. Nothing solves
for a pressure balance, so nothing iterates. That is the whole trick, and it is what makes
the model real time.

The real-time *implementation* adds an iteration on top of this for the quasi-steady case
(Eqs. 74-80, pdf pp.36-38) -- the opened compressor mass-flow loop. That belongs in the
integrator, not here. This module is the physics.

## The five states

    NG, NP    spool speeds, rpm
    P3        compressor discharge total pressure, psia
    P41       gas generator turbine inlet total pressure, psia
    P45       power turbine inlet total pressure, psia

A sixth, T41, appears in the heat-sink variant (Eqs. 48-53). Not modelled here: with no
heat-sink representation the report sets `T41 = T41_ns` [pdf p.23], which is the 5-DOF
model that Appendix B's matrices and Table B.1's trims describe.

## Units

The report's, throughout: lbm, lbm/sec, psia, deg R, rpm, ft*lbf, Btu/lbm. See
`t700.units`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from t700 import constants as c
from t700 import corrections, maps, thermo
from t700.units import BTU_TO_FTLBF, RAD_PER_SEC_TO_RPM, RPM_TO_RAD_PER_SEC


@dataclass(frozen=True)
class Ambient:
    """Flight condition. Station 1 and 2 are set equal to it [Eqs. 1-2].

    The report justifies this: *"Pressure and temperature stagnation effects have been
    found to roughly offset the aircraft inlet losses over most of the aircraft flight
    envelope"* [pdf p.22]. So there is no inlet model -- P2 and T2 *are* ambient.
    """

    p_amb_psia: float = 14.696
    t_amb_degR: float = 518.67

    @property
    def is_standard_day(self) -> bool:
        return corrections.ambient_is_standard(self.t_amb_degR, self.p_amb_psia)


STANDARD_DAY = Ambient()
"""Standard sea level, the report's trim condition [Table B.1, pdf p.67]. A frozen
singleton so it can serve as a default argument without constructing one per call."""


@dataclass(frozen=True)
class State:
    """The five-state vector, in the report's own order [pdf p.27, below Eq. 54]."""

    ng_rpm: float
    np_rpm: float
    p3_psia: float
    p41_psia: float
    p45_psia: float

    def as_array(self) -> np.ndarray:
        return np.array([self.ng_rpm, self.np_rpm, self.p3_psia, self.p41_psia, self.p45_psia])

    @staticmethod
    def from_array(v) -> State:
        return State(float(v[0]), float(v[1]), float(v[2]), float(v[3]), float(v[4]))


@dataclass(frozen=True)
class Frame:
    """Everything one pass through the gas path produces."""

    # corrected parameters
    theta2: float
    delta2: float
    ngc_pct: float
    # compressor
    ps3_psia: float
    wa2c_pps: float
    wa2_pps: float
    t3_degR: float
    h3: float
    h2: float
    # bleeds
    b1: float
    b2: float
    b3: float
    wa24_bl_pps: float
    wa3_bl_pps: float
    wa3_pps: float
    # combustor
    wa31_pps: float
    far: float
    eta_b: float
    h41_ns: float
    t41_ns_degR: float
    t41_degR: float
    h41: float
    theta41: float
    # gas generator turbine
    dh_gt: float
    h44: float
    w41_pps: float
    # station 4.5
    h45: float
    t45_degR: float
    theta45: float
    # power turbine
    ps9_psia: float
    p49_psia: float
    dh_pt: float
    w45c_pps: float
    w45_pps: float
    h49: float
    t49_degR: float
    # torques, ft*lbf
    q_c_ftlbf: float
    q_gt_ftlbf: float
    q_pt_ftlbf: float
    # derivatives of the five states
    dp3_dt: float
    dp41_dt: float
    dp45_dt: float
    dng_dt: float
    dnp_dt: float
    # diagnostics
    combustor_flow_substituted: bool
    """True when Eq. 18's pressure drop was not positive and the combustor flow and
    fuel-air ratio were substituted rather than computed. See `frame`.

    It was `combustor_dp_negative` until 2026-09-13 and missed `dp == 0` exactly, which
    takes the same `far = 0` branch with a non-negative drop. It was also computed and
    read by nothing, which is how a silent substitution stays silent."""


# The report writes 778.12 * (60 / 2*pi) * (1/N) in Eqs. 39-41. Torque is power over
# shaft speed: power in ft*lbf/sec is 778.12 * (Btu/sec), and omega is N * 2*pi/60.
_TORQUE_SCALE = BTU_TO_FTLBF * RAD_PER_SEC_TO_RPM


def frame(
    state: State,
    wf_pps: float,
    ambient: Ambient = STANDARD_DAY,
    q_req_ftlbf: float = 0.0,
    j_load: float = 0.0,
    t41_degR: float | None = None,
) -> Frame:
    """Evaluate Eqs. 1-49 once.

    Args:
        state: the five states.
        wf_pps: metered fuel flow, lbm/sec. Table B.1 quotes lbm/hr -- convert.
        ambient: flight condition.
        q_req_ftlbf: load torque on the power turbine, Qreq [Eq. 47]. An **input**: it
            comes from the Gen Hel UH-60A simulation, which this report consumes and does
            not contain.
        j_load: load inertia added to J_PT [Eq. 46]. External to the report; recovered
            from Appendix B as `constants.J_LOAD_UH60A` (open question #6, **partly**
            closed -- `J_load` is recovered, `dQreq/dNP` is still external).
        t41_degR: drive T41 instead of computing it from T41_ns. `None` -- the default --
            is Eq. 23's `T41 = T41_ns`, the 5-DOF configuration. Supplying it makes T41 an
            independent input, which is how the heat-sink (3-/6-DOF) linear models are
            extracted. It changes nothing at a trim, because Eq. 50 has unit DC gain and
            the trim value of T41 *is* T41_ns.

    Returns:
        Every station value, the three torques, and the five derivatives.
    """
    ng, np_, p3, p41, p45 = (
        state.ng_rpm,
        state.np_rpm,
        state.p3_psia,
        state.p41_psia,
        state.p45_psia,
    )

    # --- inlet and station 2, Eqs. 1-6 ------------------------------------------------
    p2 = ambient.p_amb_psia  # (1)
    t2 = ambient.t_amb_degR  # (2)
    h2 = thermo.h2_from_t2(t2)  # (3)  no intercept, unlike H3 and H41
    ps3 = c.K_PS3 * p3  # (4)
    theta2 = corrections.theta2(t2)  # (5)
    ngc = corrections.corrected_speed(ng, theta2)  # (6)

    # --- compressor, Eqs. 7-11 --------------------------------------------------------
    pr = ps3 / p2
    ngc_pct = 100.0 * ngc / c.NG_DES
    wa2c = float(maps.f1()(pr, ngc_pct))  # (7)
    delta2 = corrections.delta2(p2)  # (8)
    wa2 = corrections.physical_flow(wa2c, theta2, delta2)  # (9)
    t3 = t2 * float(maps.f2()(pr))  # (10) a ratio multiplier, not an additive rise
    h3 = thermo.h3_from_t3(t3)  # (11)

    # --- bleeds, Eqs. 12-17 -----------------------------------------------------------
    b1 = float(maps.f3()(ngc_pct))  # (12)
    b2 = float(maps.f4()(wa2c))  # (13)
    b3 = float(maps.f5()(wa2c))  # (14)
    wa24_bl = wa2 * (b1 + b2)  # (15)
    wa3_bl = wa2 * (b3 + c.K_B3)  # (16)  the K_b3 part never returns to the gas path
    wa3 = wa2 - wa24_bl  # (17)

    # --- combustor and station 4.1, Eqs. 18-25 ----------------------------------------
    # Eq. 18 inverts the combustor pressure-drop law, so it needs P3 > P41. Two guards
    # sit on it and both are SILENT SUBSTITUTIONS rather than clamps: `max(dp, 0)`
    # replaces a reversed pressure drop with zero flow, and the `far` fallback then
    # replaces an infinite fuel-air ratio with zero -- which injects fuel into the mass
    # balance (Eq. 43 carries `wa31 + wf`) and releases none of its heat, since Eq. 21's
    # `eta_b * far * HVF` term vanishes. The state is unphysical either way; what matters
    # is that the model says so instead of returning a number.
    #
    # **Latent, and measured**: zero occurrences at every trim from 100 to 800 lbm/hr and
    # across both published transients. `combustor_flow_substituted` is what makes that a
    # checked claim rather than an assumption -- see
    # `tests/test_mass_conservation.py::test_the_combustor_flow_guards_never_fire`.
    dp = p3 * (p3 - p41)
    wa31 = float(np.sqrt(max(dp, 0.0) / (c.K_DPB * t3)))  # (18)
    substituted = dp <= 0.0 or wa31 <= 0.0
    far = wf_pps / wa31 if wa31 > 0.0 else 0.0  # (19)
    eta_b = float(maps.f6()(far))  # (20)
    h41_ns = (h3 + eta_b * far * c.HVF) / (1.0 + far)  # (21)
    t41_ns = thermo.t41_from_h41(h41_ns)  # (22)
    # (23): T41 = T41_ns when no heat-sink representation is used [pdf p.23]. Passing
    # `t41_degR` drives it instead, which is what makes T41 a sixth **state** rather than
    # a computed quantity -- the 6-DOF configuration, and the form the report's own
    # derivative extraction uses: "During derivative extraction, T41 was held fixed at the
    # trim value and the change in T41 was determined for each state perturbation"
    # [pdf p.31].
    t41 = t41_ns if t41_degR is None else float(t41_degR)
    h41 = thermo.h41_from_t41(t41)  # (24)
    theta41 = thermo.theta41_from_t41(t41)  # (25)

    # --- gas generator turbine, Eqs. 26-28 --------------------------------------------
    dh_gt = theta41 * float(maps.f7()(p45 / p41))  # (26)
    h44 = h41 - dh_gt  # (27)
    w41 = c.K_WGT * p41 / np.sqrt(theta41)  # (28)

    # --- station 4.5 mixing, Eqs. 29-31 -----------------------------------------------
    h45 = thermo.h45_from_h44(h44)  # (29)
    t45 = thermo.t45_from_h45(h45)  # (30)
    theta45 = thermo.theta45_from_t45(t45)  # (31)

    # --- exhaust first, Eqs. 37-38: P49 sets the power turbine's back pressure ---------
    ps9 = ambient.p_amb_psia  # (37)
    p49 = ps9 * float(maps.f10()(ngc_pct))  # (38)  f10 > 1; see maps.f10 and question #5

    # --- power turbine, Eqs. 32-36 ----------------------------------------------------
    dh_pt = theta45 * float(maps.f8()(p49 / p45))  # (32)
    w45c = float(maps.f9()(ps9 / p45))  # (33)
    w45 = w45c * p45 / np.sqrt(theta45)  # (34)
    h49 = h45 - dh_pt  # (35)
    t49 = thermo.t49_from_h49(h49)  # (36)

    # --- torques, Eqs. 39-41 ----------------------------------------------------------
    q_c = _TORQUE_SCALE / ng * (wa2 * (c.K_QC_1 * h3 - h2) + wa3 * c.K_QC_2 * h3)  # (39)
    q_gt = _TORQUE_SCALE / ng * w41 * dh_gt  # (40)
    q_pt = _TORQUE_SCALE / np_ * w45 * dh_pt - c.K_DAMP * RPM_TO_RAD_PER_SEC * (
        np_ - c.NP_DES
    )  # (41)  damping is on speed *displacement*; the equation governs, not the prose

    # --- state derivatives, Eqs. 42-45, 47 --------------------------------------------
    # Eq. 43 as printed reads (WA31 - Wf - W41). Eq. 73 [pdf p.36] prints +Wf for the same
    # balance and Eq. 77 follows from that one, not from Eq. 43. Fuel adds mass to the
    # combustor, so +Wf is also the only physical reading. Open questions #14 and #29.
    dp3 = c.K_V3 * t3 * (wa3 - wa3_bl - wa31)  # (42)
    dp41 = c.K_V41 * t41 * (wa31 + wf_pps - w41)  # (43), sign per Eq. 73
    dp45 = c.K_V45 * t45 * (w41 - w45 + b3 * c.K_BL * wa2)  # (44)
    dng = RAD_PER_SEC_TO_RPM * (q_gt - q_c) / c.J_GT  # (45)
    j = c.J_PT + j_load  # (46)
    dnp = RAD_PER_SEC_TO_RPM * (q_pt - q_req_ftlbf) / j  # (47)

    return Frame(
        theta2=theta2,
        delta2=delta2,
        ngc_pct=ngc_pct,
        ps3_psia=ps3,
        wa2c_pps=wa2c,
        wa2_pps=wa2,
        t3_degR=t3,
        h3=h3,
        h2=h2,
        b1=b1,
        b2=b2,
        b3=b3,
        wa24_bl_pps=wa24_bl,
        wa3_bl_pps=wa3_bl,
        wa3_pps=wa3,
        wa31_pps=wa31,
        far=far,
        eta_b=eta_b,
        h41_ns=h41_ns,
        t41_ns_degR=t41_ns,
        t41_degR=t41,
        h41=h41,
        theta41=theta41,
        dh_gt=dh_gt,
        h44=h44,
        w41_pps=w41,
        h45=h45,
        t45_degR=t45,
        theta45=theta45,
        ps9_psia=ps9,
        p49_psia=p49,
        dh_pt=dh_pt,
        w45c_pps=w45c,
        w45_pps=w45,
        h49=h49,
        t49_degR=t49,
        q_c_ftlbf=q_c,
        q_gt_ftlbf=q_gt,
        q_pt_ftlbf=q_pt,
        dp3_dt=dp3,
        dp41_dt=dp41,
        dp45_dt=dp45,
        dng_dt=dng,
        dnp_dt=dnp,
        combustor_flow_substituted=bool(substituted),
    )


def residuals(state: State, wf_pps: float, ambient: Ambient = STANDARD_DAY) -> np.ndarray:
    """The four gas-generator derivatives that must vanish at trim.

    NP is excluded: Eq. 47 balances Q_PT against Qreq, and Qreq is supplied by the Gen Hel
    helicopter simulation, not by this report. At a trim we hold NP at its printed value
    and read Qreq off as whatever Q_PT comes out to be.
    """
    f = frame(state, wf_pps, ambient)
    return np.array([f.dng_dt, f.dp3_dt, f.dp41_dt, f.dp45_dt])
