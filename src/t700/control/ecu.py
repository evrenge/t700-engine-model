"""The T700 electrical control unit -- Appendix C, Figures C1 through C8.

The ECU governs power turbine speed and limits station 4.5 temperature. It does not meter
fuel: its one output, `SPDG`, trims the HMU through the torque motor (`hmu.SPDG_NULL`
leaves the HMU alone).

[Fig. C1, pdf p.85]

    SPDER = PCNP/(CTPL*s+1) - XQLO - PCPRF          speed error, %NP
    ET45  = T45EL - T45REF                          temperature error, deg R

Each error is compensated (`SPDSF`, `TSIG`) and a maximum error selector passes the larger
to the P+I: T4.5 limiting takes over from NP governing when its error is the bigger. Sign
trap: `SPDER` is positive when sensed NP exceeds the reference and drives `SPDG` up, and
more `SPDG` means less fuel (the HMU subtracts `TMRU`). At every Table B.1 trim the sensed
`T45EL` sits about 360-570 deg R below `T45REF`, so the temperature path is inert in steady
flight.

Single-engine, as the report draws it: Fig. C2's load-share switch is open, so `XQLO` is
zero (`engine2_torq45` lets that be tested), and Fig. C8's `ZLOLIM` is on its -1.0 branch,
the Table C.1 value.

[Fig. C3] `B4` switches the rate gain from `ZK7` to `ZK7 + ZK1` on
`1000*B6*g(SPDER) + Y > CORR`, `g` a relay on large speed errors and `Y` the integral of
`TRQL - CR`: the gain rises on a large speed error or a sustained load (open question #11).
"""

from __future__ import annotations

from dataclasses import dataclass

from t700 import constants as engine_c
from t700.control import constants as c
from t700.control import schedules
from t700.control._blocks import clamp, deadband, integrator, lag, lead_lag


@dataclass(frozen=True)
class ECUState:
    """Every continuous state in Figures C1-C8, named for the block that owns it."""

    pcnp_lag: float = 0.0
    """Sensed power turbine speed, 1/(CTPL*s + 1) [Fig. C1]."""
    trql_lag1: float = 0.0
    trql_lag2: float = 0.0
    """Sensed power turbine torque, two lags TL1 and TL2 [Fig. C2]."""
    xqlo_lag1: float = 0.0
    xqlo_lag2: float = 0.0
    """Load-share trim, 1.1/(CT7*s+1) then ZK8/(T17*s+1) [Fig. C2]. Inert single-engine."""
    rate_lag13: float = 0.0
    rate_lag14: float = 0.0
    """Governor rate compensation, CT13 then CT14 [Fig. C3]. Their difference is the
    pseudo-derivative that gives the loop its lead."""
    torque_int: float = 0.0
    """Y = integral of (TRQL - CR), limited to [YLOLIM, YHILIM] [Fig. C3]."""
    gov_leadlag: float = 0.0
    gov_lag: float = 0.0
    """Governor dynamics, (T11*s+1)/(CT2*s+1) then 1/(CT12*s+1) [Fig. C4]."""
    harness_lag: float = 0.0
    """Thermocouple harness, 1/(TLGE*s + 1) [Fig. C5]."""
    t45el_lag: float = 0.0
    """Sensed T4.5 -- a lag whose time constant is itself a lookup, F_EC1 [Fig. C5]."""
    t45comp_lag: float = 0.0
    t45comp_leadlag: float = 0.0
    """T4.5 compensation, ZK3/(CT9*s+1) then (T8*s+1)/(T10*s+1) [Fig. C6]."""
    pi_int: float = 0.0
    """P+I integrator, limited to [ZLOLIM, ZHILIM] [Fig. C7]. ZLOLIM is itself a function
    of the other engine's torque [Fig. C8]."""
    spdg_lag: float = 0.0
    """P+I output lag, 1/(CT16*s + 1) [Fig. C7]. Its state is SPDG."""


@dataclass(frozen=True)
class ECUInputs:
    """What the ECU reads. Everything but `pcprf` comes from the engine."""

    np_rpm: float
    torq45_ftlbf: float
    t45_degR: float
    w45_pps: float
    """Power turbine inlet mass flow, lbm/sec.

    The nomenclature [pdf p.80] says lb/hr, but `W45R = W45*sqrt(T45L)/P45` feeds `F_EC1`,
    whose printed x axis runs 1.0 to 15.0, and only lbm/sec lands on it. A labelling slip."""
    p45_psia: float
    pcprf_pct: float = 100.0
    """Reference power turbine speed as set by the cockpit control, percent of design."""
    engine2_torq45_ftlbf: float | None = None
    """Sensed torque of the *other* engine. `None` is the report's own single-engine
    implementation: Fig. C2's load-share switch drawn open, and Fig. C8's ZLOLIM on its
    low-torque branch."""


@dataclass(frozen=True)
class ECUOutputs:
    """Named for the report's own signals."""

    spdg: float
    spder: float
    et45: float
    spdsf: float
    tsig: float
    spdss: float
    trql: float
    xqlo: float
    b4: float
    tau45_s: float
    w45r: float
    t45el: float
    limiting: str
    """Which path won the maximum error selector: `speed` or `t45`."""


def _zlolim(engine2_torq45: float | None) -> float:
    """[Fig. C8, pdf p.88] The P+I integrator's lower limit is not a constant."""
    if engine2_torq45 is None:
        return c.ZLOLIM
    return c.ZLOLIM_HIGH_TORQUE if engine2_torq45 >= c.ENGINE2_TORQUE_THRESHOLD else c.ZLOLIM


def seed(u: ECUInputs, spdg: float = 0.0) -> ECUState:
    """A state consistent with holding `u`, as `hmu.seed` is for the HMU. Open question #55.

    Every lag starts at its input. `spdg` seeds the P+I integrator and the CT16 output lag
    to the torque-motor demand the caller wants held; it is not derivable from the ECU's own
    inputs. `loop.seed` passes `hmu.SPDG_NULL`, the demand the HMU's torque motor is seeded
    to hold. Seeding the two halves differently injects an uncommanded half-volt SPDG step
    on frame 1.
    """
    pcnp = u.np_rpm * 100.0 / engine_c.NP_DES
    trql = u.torq45_ftlbf
    t45l = u.t45_degR + c.T45COR
    return ECUState(
        pcnp_lag=pcnp,
        trql_lag1=trql,
        trql_lag2=trql,
        harness_lag=t45l,
        t45el_lag=t45l,
        pi_int=spdg,
        spdg_lag=spdg,
    )


def step(state: ECUState, u: ECUInputs, dt: float) -> tuple[ECUState, ECUOutputs]:
    """Advance the ECU one frame. Returns the new state and the signals of Figs. C1-C8."""
    # --- C-1/C-2, Fig. C1: sensed speed and the speed error ------------------------
    pcnp = u.np_rpm * 100.0 / engine_c.NP_DES
    pcnp_lag = lag(state.pcnp_lag, pcnp, c.CTPL, dt)

    # --- C-5, Fig. C2: sensed power turbine torque ---------------------------------
    trql_lag1 = lag(state.trql_lag1, u.torq45_ftlbf, c.TL1, dt)
    trql_lag2 = lag(state.trql_lag2, trql_lag1, c.TL2, dt)
    trql = trql_lag2

    # --- C-6/C-7, Fig. C2: load share. Inert with the switch drawn open ------------
    if u.engine2_torq45_ftlbf is None:
        trqer = 0.0
    else:
        trqer = u.engine2_torq45_ftlbf - trql
    xqlo_lag1 = lag(state.xqlo_lag1, c.LOAD_SHARE_LAG_GAIN * trqer, c.CT7, dt)
    xqlo_lag2 = lag(state.xqlo_lag2, c.ZK8 * (clamp(xqlo_lag1, c.DBIAS, c.CE) - c.DBIAS), c.T17, dt)
    xqlo = xqlo_lag2

    spder = pcnp_lag - xqlo - u.pcprf_pct

    # --- C-8..C-11, Fig. C3: governor rate compensation ----------------------------
    db = deadband(spder, c.CB)
    rate_lag13 = lag(state.rate_lag13, db, c.CT13, dt)
    rate_lag14 = lag(state.rate_lag14, rate_lag13, c.CT14, dt)
    pseudo_deriv = (c.ZK5 / c.CT14) * (rate_lag13 - rate_lag14)
    u_comp = c.ZK9 * spder + c.ZK10 * db + pseudo_deriv

    # --- C-12..C-15, Fig. C3: the nonlinear NP loop gain ---------------------------
    relay = c.RELAY_OUTPUT if (spder < c.NP_LOOP_RELAY_LOW or spder > c.NP_LOOP_RELAY_HIGH) else 0.0
    torque_int = integrator(state.torque_int, trql - c.CR, 1.0, c.YLOLIM, c.YHILIM, dt)
    b4 = 1.0 if (c.NP_LOOP_RELAY_GAIN * c.B6 * relay + torque_int) > c.CORR else 0.0
    spds1 = c.ZK7 * u_comp + b4 * (c.ZK1 * u_comp)

    # --- C-16, Fig. C4: governor dynamics ------------------------------------------
    gov_leadlag, gov_out = lead_lag(state.gov_leadlag, spds1, c.T11, c.CT2, dt)
    gov_lag = lag(state.gov_lag, gov_out, c.CT12, dt)
    spdsf = gov_lag

    # --- C-17..C-21, Fig. C5: thermocouple harness and sensed T4.5 -----------------
    t45e = u.t45_degR + c.T45COR
    harness_lag = lag(state.harness_lag, t45e, c.TLGE, dt)
    t45l = harness_lag
    w45r = u.w45_pps * (t45l**0.5) / u.p45_psia if u.p45_psia > 0.0 else 0.0
    tau45 = float(schedules.f_ec1()(w45r, t45l))
    t45el_lag = lag(state.t45el_lag, t45l, tau45, dt)
    t45el = t45el_lag

    # --- C-3/C-22, Figs. C1/C6: temperature error and its compensation -------------
    et45 = t45el - c.T45REF
    t45comp_lag = lag(state.t45comp_lag, c.ZK3 * et45, c.CT9, dt)
    t45comp_leadlag, tsig = lead_lag(state.t45comp_leadlag, t45comp_lag, c.T8, c.T10, dt)

    # --- C-4, Fig. C1: the maximum error selector ----------------------------------
    spdss = max(spdsf, tsig)
    limiting = "t45" if tsig > spdsf else "speed"

    # --- C-23..C-25, Figs. C7/C8: proportional plus integral -----------------------
    lo = _zlolim(u.engine2_torq45_ftlbf)
    pi_int = integrator(state.pi_int, spdss, c.XKINTG, lo, c.ZHILIM, dt)
    spdsp = c.XKPROP * spdss + pi_int
    spdg_lag = lag(state.spdg_lag, spdsp, c.CT16, dt)

    # Every field is rewritten, so construct directly: `dataclasses.replace` re-introspects
    # the field list on every call.
    new = ECUState(
        pcnp_lag=pcnp_lag,
        trql_lag1=trql_lag1,
        trql_lag2=trql_lag2,
        xqlo_lag1=xqlo_lag1,
        xqlo_lag2=xqlo_lag2,
        rate_lag13=rate_lag13,
        rate_lag14=rate_lag14,
        torque_int=torque_int,
        gov_leadlag=gov_leadlag,
        gov_lag=gov_lag,
        harness_lag=harness_lag,
        t45el_lag=t45el_lag,
        t45comp_lag=t45comp_lag,
        t45comp_leadlag=t45comp_leadlag,
        pi_int=pi_int,
        spdg_lag=spdg_lag,
    )
    out = ECUOutputs(
        spdg=spdg_lag,
        spder=spder,
        et45=et45,
        spdsf=spdsf,
        tsig=tsig,
        spdss=spdss,
        trql=trql,
        xqlo=xqlo,
        b4=b4,
        tau45_s=tau45,
        w45r=w45r,
        t45el=t45el,
        limiting=limiting,
    )
    return new, out
