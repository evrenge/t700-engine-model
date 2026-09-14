"""The T700 electrical control unit -- Appendix C, Figures C1 through C8.

The ECU governs **power turbine speed** and limits **station 4.5 temperature**. It does
not meter fuel: its entire authority is one signal, `SPDG`, which trims the HMU through
the torque motor. With `SPDG` held at `hmu.SPDG_NULL` the HMU runs alone, which is how the
engine was flown before this module existed.

## Two paths, and the larger error wins

[Fig. C1, pdf p.85]

    SPDER = PCNP/(CTPL*s+1) - XQLO - PCPRF          speed error, %NP
    ET45  = T45EL - T45REF                          temperature error, deg R

The speed error goes through rate compensation and governor dynamics to `SPDSF`; the
temperature error through its own compensation to `TSIG`. A **MAXIMUM ERROR SELECTOR**
passes the larger of the two to the proportional-plus-integral compensation. That is the
ECU's limit-arbitration point: T4.5 limiting takes over from NP governing exactly when its
error signal is the bigger number.

Note the sign. `SPDER` is *positive when sensed NP exceeds the reference*, not the other
way round, and it drives `SPDG` up; the HMU subtracts `TMRU` from its droop line, so more
`SPDG` means less fuel. A speed-topping trim, not a classical error = reference - measured.

## Why `T45REF` is 2004 deg R and nothing ever reaches it

Table B.1's three trims sit at T45 = 1632 / 1501 / 1424 deg R, so `ET45` is 350 to 570 deg
below its reference at every printed condition. The temperature path exists for
overtemperature transients and is inert in steady flight -- which is the right behaviour
for a limiter, and is why the speed path is the one to check first.

## What is single-engine here

Two things collapse in a one-engine model, and both are drawn that way in the report:

* **Load share** [Fig. C2]. The switch labelled "ONE ENGINE IMPLEMENTATION" is drawn
  *open*, so `TRQER` is zero, and `clamp(0, DBIAS, CE) - DBIAS = 0` makes `XQLO` zero too.
  The path is inert by construction rather than by our choice; `engine2_torq45` is exposed
  so that can be tested rather than assumed.
* **The P+I integrator's lower limit** [Fig. C8]. `ZLOLIM` is a two-level function of the
  *other* engine's sensed torque: -1.0 below 180 ft*lbf and -0.3 above. With no second
  engine it is the -1.0 branch, which is the value Table C.1 prints.

## The nonlinear NP loop gain

[Fig. C3] `B4` switches the rate-compensation gain between `ZK7` and `ZK7 + ZK1` -- 0.3 and
2.0, a factor of 6.7 -- on `1000*B6*g(SPDER) + Y > CORR`, where `g` is a relay on large
speed errors and `Y` integrates `TRQL - CR`. So the loop gain rises either when the speed
error is large or when the engine has been loaded for long enough, and unwinds when it is
not. It is **not** permanently 1: open question #11 measured the crossing at 0.10-0.36 s at
flight power and the unwind at >= 2 s.
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
    """Governor rate compensation, CT13 then CT14 [Fig. C3]. Their *difference* is the
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
    """Power turbine inlet mass flow, **lbm/sec**.

    The nomenclature [pdf p.80] says lb/hr. It cannot be: `W45R = W45*sqrt(T45L)/P45`
    feeds `F_EC1`, whose printed x axis runs 1.0 to 15.0, and at the three Table B.1 trims
    lbm/sec gives 8.18 / 8.42 / 8.56 while lbm/hr gives about 30000 -- 2000x off the axis.
    A labelling slip, like Table C.1 calling `CR` nondimensional (open question #11)."""
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
    to the torque-motor demand the caller wants held.

    **Why it takes an argument at all.** The default of zero was the honest choice for the
    ECU *alone*: the P+I integrator only stops when `SPDSS` is zero, so its trim value is
    not derivable from the ECU's own inputs and inventing one would be inventing a number.
    But `loop.seed` was seeding the HMU's torque motor to hold `SPDG_NULL = 0.49496` while
    leaving this at zero, so frame 1 delivered **SPDG = 0 into a torque motor seeded for
    0.495** -- a half-volt step the control never commanded. That is not an unknown
    initial condition, it is two halves of one seed contradicting each other, and the
    2026-09-13 control-systems audit measured what it costs:

        as shipped                Wf peak 627.9 pph (+31.8 %), NP +1.42 %, T41 2500.8 degR,
                                  50 `f6` clamps -- the model leaves its digitized envelope
        ECU seeded at SPDG_NULL   Wf peak 477.4 pph (+0.23 %), NP +0.03 %, T41 2293.6, none

    It decays in about 3 s, so the 4000-frame closed-loop trim tests washed it out entirely
    -- settled fuel flow moves 476.47 -> 475.95 pph -- but every transient started from
    `loop.seed` began with that slam.

    `SPDG_NULL` is not fitted: the CT16 lag has unit DC gain and `SPDSS` is zero at trim,
    and the loop settles at SPDG = 0.5038 on its own from either seed. It is the same
    principle `realtime.from_trim` already follows -- start at equilibrium, because that is
    the only choice that makes a trimmed engine sit still.
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

    # Every field of `ECUState` is rewritten here, so this is a construction and not a
    # modification; `dataclasses.replace` introspects the field list on each call to
    # rediscover that, which cost ~4 us a frame for nothing. Naming them keeps the
    # positional order from mattering.
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
