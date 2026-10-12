"""One frame of the T700 gas path: Eqs. 1-49, in the report's own order.

Given the five states and a fuel flow, `frame()` walks the gas path once -- inlet,
compressor, bleeds, combustor, gas generator turbine, station 4.5 mixing, power turbine,
exhaust -- and returns every station value, the three torques and the five derivatives.

One forward pass is enough because the three pressures are pure integrators (Eqs. 42-44):
nothing solves for a pressure balance. The real-time implementation's opened iteration
(Eqs. 74-80, pdf pp.36-38) belongs in `t700.realtime`.

    NG, NP    spool speeds, rpm
    P3        compressor discharge total pressure, psia
    P41       gas generator turbine inlet total pressure, psia
    P45       power turbine inlet total pressure, psia

T41 is a sixth state only with the heat sink (Eqs. 48-53); without it `T41 = T41_ns`
[pdf p.23], the 5-DOF model of Appendix B and Table B.1. Units are the report's
(`t700.units`).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import NamedTuple

import numpy as np

from t700 import constants as c
from t700 import corrections, maps, thermo
from t700.units import BTU_TO_FTLBF, RAD_PER_SEC_TO_RPM, RPM_TO_RAD_PER_SEC


@dataclass(frozen=True)
class Ambient:
    """Flight condition. Stations 1 and 2 are set equal to it [Eqs. 1-2].

    No inlet model: stagnation effects "roughly offset the aircraft inlet losses" [pdf p.22].
    """

    p_amb_psia: float = corrections.P_STD_PSIA
    t_amb_degR: float = corrections.T_STD_DEGR


STANDARD_DAY = Ambient()
"""Standard sea level, the report's trim condition [Table B.1, pdf p.67]."""


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
    fuel-air ratio were substituted rather than computed. See `frame`."""


# The report writes 778.12 * (60 / 2*pi) * (1/N) in Eqs. 39-41. Torque is power over
# shaft speed: power in ft*lbf/sec is 778.12 * (Btu/sec), and omega is N * 2*pi/60.
_TORQUE_SCALE = BTU_TO_FTLBF * RAD_PER_SEC_TO_RPM


# The equations `frame` and `realtime.step` share, written once in the report's order;
# `realtime.step` calls them between its own pressure solves.


class Compressor(NamedTuple):
    """Station 2, the compressor and the bleeds: everything Eqs. 1-17 produce."""

    h2: float
    ps3: float
    theta2: float
    ngc_pct: float
    wa2c: float
    delta2: float
    wa2: float
    t3: float
    h3: float
    b1: float
    b2: float
    b3: float
    wa24_bl: float
    wa3_bl: float
    wa3: float


def compressor(ng_rpm: float, p3_psia: float, ambient: Ambient) -> Compressor:
    """Eqs. 1-17 from NG, P3 and the flight condition alone."""
    # --- inlet and station 2, Eqs. 1-6 ------------------------------------------------
    p2 = ambient.p_amb_psia  # (1)
    t2 = ambient.t_amb_degR  # (2)
    h2 = thermo.h2_from_t2(t2)  # (3)  no intercept, unlike H3 and H41
    ps3 = c.K_PS3 * p3_psia  # (4)
    theta2 = corrections.theta2(t2)  # (5)
    ngc = corrections.corrected_speed(ng_rpm, theta2)  # (6)

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
    return Compressor(
        h2,
        ps3,
        theta2,
        ngc_pct,
        wa2c,
        delta2,
        wa2,
        t3,
        h3,
        b1,
        b2,
        b3,
        wa24_bl,
        wa3_bl,
        wa3,
    )


def combustor(far: float, h3: float) -> tuple[float, float, float]:
    """Eqs. 20-22: `(eta_b, H41_ns, T41_ns)` from the fuel-air ratio and H3."""
    eta_b = float(maps.f6()(far))  # (20)
    h41_ns = (h3 + eta_b * far * c.HVF) / (1.0 + far)  # (21)
    t41_ns = thermo.t41_from_h41(h41_ns)  # (22)
    return eta_b, h41_ns, t41_ns


def station41(t41_degR: float) -> tuple[float, float]:
    """Eqs. 24-25: `(H41, theta41)` from T41."""
    h41 = thermo.h41_from_t41(t41_degR)  # (24)
    theta41 = thermo.theta41_from_t41(t41_degR)  # (25)
    return h41, theta41


def gg_turbine_and_station45(
    h41: float, theta41: float, p41: float, p45: float
) -> tuple[float, float, float, float, float, float]:
    """Eqs. 26-31: `(dH_GT, H44, W41, H45, T45, theta45)`."""
    # --- gas generator turbine, Eqs. 26-28 --------------------------------------------
    dh_gt = theta41 * float(maps.f7()(p45 / p41))  # (26)
    h44 = h41 - dh_gt  # (27)
    w41 = c.K_WGT * p41 / sqrt(theta41)  # (28)

    # --- station 4.5 mixing, Eqs. 29-31 -----------------------------------------------
    h45 = thermo.h45_from_h44(h44)  # (29)
    t45 = thermo.t45_from_h45(h45)  # (30)
    theta45 = thermo.theta45_from_t45(t45)  # (31)
    return dh_gt, h44, w41, h45, t45, theta45


def power_turbine(theta45: float, p45: float, p49: float, ps9: float) -> tuple[float, float, float]:
    """Eqs. 32-34: `(dH_PT, W45c, W45)`."""
    dh_pt = theta45 * float(maps.f8()(p49 / p45))  # (32)
    w45c = float(maps.f9()(ps9 / p45))  # (33)
    w45 = w45c * p45 / sqrt(theta45)  # (34)
    return dh_pt, w45c, w45


def torques(
    ng_rpm: float,
    np_rpm: float,
    comp: Compressor,
    w41: float,
    dh_gt: float,
    w45: float,
    dh_pt: float,
) -> tuple[float, float, float]:
    """Eqs. 39-41: `(Q_C, Q_GT, Q_PT)`, ft*lbf."""
    wa2, wa3, h2, h3 = comp.wa2, comp.wa3, comp.h2, comp.h3
    q_c = _TORQUE_SCALE / ng_rpm * (wa2 * (c.K_QC_1 * h3 - h2) + wa3 * c.K_QC_2 * h3)  # (39)
    q_gt = _TORQUE_SCALE / ng_rpm * w41 * dh_gt  # (40)
    q_pt = _TORQUE_SCALE / np_rpm * w45 * dh_pt - c.K_DAMP * RPM_TO_RAD_PER_SEC * (
        np_rpm - c.NP_DES
    )  # (41)  damping is on speed *displacement*; the equation governs, not the prose
    return q_c, q_gt, q_pt


def shaft_change(q_drive: float, q_load: float, inertia: float, dt: float = 1.0) -> float:
    """Eqs. 45 and 47: a spool's speed change over `dt` seconds, rpm.

    The default `dt = 1` is the derivative itself, bit for bit: `1.0 * x` is `x` exactly,
    so `frame` and `realtime.step` share one expression and one rounding order.
    """
    return dt * RAD_PER_SEC_TO_RPM * (q_drive - q_load) / inertia


def heat_sink_taus(t41_degR: float, w41_pps: float, ngc_pct: float) -> tuple[float, float]:
    """Eqs. 51-53: `(tau_a, tau_b)` of the station 4.1 heat sink [pdf p.26].

    `tau_a = M c_pm/(h A_m)` and `tau_b = M c_pm/(W_g c_pg)`, so Eq. 50 is
    `((tau_a - tau_b) s + 1) / (tau_a s + 1)`. See `realtime._heat_sink`.
    """
    tau_a = c.TC_T41 * t41_degR**0.5 / w41_pps**0.8  # (51)
    tau_b = float(maps.f_hs()(ngc_pct)) / w41_pps  # (52), (53)
    return tau_a, tau_b


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
        q_req_ftlbf: load torque on the power turbine, Qreq [Eq. 47]. An input, from the
            Gen Hel UH-60A simulation the report does not contain.
        j_load: load inertia added to J_PT [Eq. 46], recovered from Appendix B as
            `constants.J_LOAD_UH60A` (open question #6).
        t41_degR: drive T41 instead of Eq. 23's `T41 = T41_ns` (`None`, the 5-DOF
            default). How the heat-sink (3-/6-DOF) linear models are extracted; it changes
            nothing at a trim, because Eq. 50 has unit DC gain.

    Returns:
        Every station value, the three torques, and the five derivatives.
    """
    ng, np_, p3, p41, p45 = (  # np_ is NP: np is NumPy
        state.ng_rpm,
        state.np_rpm,
        state.p3_psia,
        state.p41_psia,
        state.p45_psia,
    )

    # --- inlet, compressor and bleeds, Eqs. 1-17 ----------------------------------------
    cp = compressor(ng, p3, ambient)

    # --- combustor and station 4.1, Eqs. 18-25 ----------------------------------------
    # Eq. 18 needs P3 > P41. Its two guards are silent substitutions, not clamps:
    # `max(dp, 0)` turns a reversed drop into zero flow, and the `far` fallback then puts
    # fuel into the mass balance (Eq. 43) while releasing none of its heat (Eq. 21).
    # `combustor_flow_substituted` reports it; see
    # `tests/test_mass_conservation.py::test_the_combustor_flow_guards_never_fire`.
    dp = p3 * (p3 - p41)
    wa31 = sqrt(max(dp, 0.0) / (c.K_DPB * cp.t3))  # (18)
    substituted = dp <= 0.0 or wa31 <= 0.0
    far = wf_pps / wa31 if wa31 > 0.0 else 0.0  # (19)
    eta_b, h41_ns, t41_ns = combustor(far, cp.h3)  # (20)-(22)
    # (23): T41 = T41_ns without a heat-sink representation [pdf p.23]. Passing `t41_degR`
    # makes T41 a sixth state, held fixed as the report's derivative extraction does [pdf p.31].
    t41 = t41_ns if t41_degR is None else float(t41_degR)
    h41, theta41 = station41(t41)  # (24), (25)

    # --- gas generator turbine and station 4.5 mixing, Eqs. 26-31 ---------------------
    dh_gt, h44, w41, h45, t45, theta45 = gg_turbine_and_station45(h41, theta41, p41, p45)

    # --- exhaust first, Eqs. 37-38: P49 sets the power turbine's back pressure ---------
    ps9 = ambient.p_amb_psia  # (37)
    p49 = ps9 * float(maps.f10()(cp.ngc_pct))  # (38)  f10 > 1; see maps.f10 and question #5

    # --- power turbine, Eqs. 32-36 ----------------------------------------------------
    dh_pt, w45c, w45 = power_turbine(theta45, p45, p49, ps9)  # (32)-(34)
    h49 = h45 - dh_pt  # (35)
    t49 = thermo.t49_from_h49(h49)  # (36)

    # --- torques, Eqs. 39-41 ----------------------------------------------------------
    q_c, q_gt, q_pt = torques(ng, np_, cp, w41, dh_gt, w45, dh_pt)

    # --- state derivatives, Eqs. 42-45, 47 --------------------------------------------
    # Eq. 43 as printed reads (WA31 - Wf - W41); Eq. 73 [pdf p.36] prints +Wf for the same
    # balance, Eq. 77 follows from it, and fuel adds mass. Open questions #14 and #29.
    dp3 = c.K_V3 * cp.t3 * (cp.wa3 - cp.wa3_bl - wa31)  # (42)
    dp41 = c.K_V41 * t41 * (wa31 + wf_pps - w41)  # (43), sign per Eq. 73
    dp45 = c.K_V45 * t45 * (w41 - w45 + cp.b3 * c.K_BL * cp.wa2)  # (44)
    dng = shaft_change(q_gt, q_c, c.J_GT)  # (45)
    j_pt_total = c.J_PT + j_load  # (46)
    dnp = shaft_change(q_pt, q_req_ftlbf, j_pt_total)  # (47)
    return Frame(
        theta2=cp.theta2,
        delta2=cp.delta2,
        ngc_pct=cp.ngc_pct,
        ps3_psia=cp.ps3,
        wa2c_pps=cp.wa2c,
        wa2_pps=cp.wa2,
        t3_degR=cp.t3,
        h3=cp.h3,
        h2=cp.h2,
        b1=cp.b1,
        b2=cp.b2,
        b3=cp.b3,
        wa24_bl_pps=cp.wa24_bl,
        wa3_bl_pps=cp.wa3_bl,
        wa3_pps=cp.wa3,
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
