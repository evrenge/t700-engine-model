"""The T700 hydromechanical control unit -- Appendix C, Figures C9 through C22.

The HMU is the primary fuel control loop. It governs gas generator speed on a **droop
line**, applies the idle, acceleration and deceleration cam limits, and meters fuel. The
ECU only *trims* it, through the torque motor signal `SPDG`; with that trim held at its
null the HMU runs alone, which is how this module is exercised before the ECU exists.

## The droop line, and why almost everything is a subtraction

[Fig. C9, pdf p.89]

    WFPDM = WFPTP + KNDRP*(NGREF - PCNGHL) - WFPRF - max(DWFPL + TMRU, 0)

`WFPTP` is the topping line, a function of T2 alone: the ceiling. Everything else takes
fuel *away* from it. The power-available spindle schedule `F_HM2` falls from 10.4 at
PAS = 28 deg to -0.1 at 120 deg, so more spindle means less subtraction and more fuel;
the load-demand path works the same way. Reading a summing junction backwards here
inverts the whole control, which is why the signs in `docs/notes/inventory-appendix-c.md`
were read at 600 dpi and are quoted against each relation below.

Note the units: everything on the droop line is a **fuel flow parameter**, Wf/Ps3, not a
fuel flow. It becomes one only at the metering valve, `WFMV = HMUSEL * PS3L`.

## Limit arbitration is a fixed cascade, not a min-max soup

[Fig. C18, pdf p.92]

    HMUSEL = MAX( MIN( MAX(WFPDM, WFIDM), WFPAC ), WFPDC )

Idle floor first, then the acceleration ceiling, then the deceleration floor. The order is
drawn, not inferred, and it matters: the deceleration floor outranks the acceleration
ceiling.

## What this module does not do

The report's real-time implementation eliminated **fuel control below flight-idle power**
[pdf p.38], so the idle floor `WFIDM` is inactive in every condition the engine model is
valid over -- it only outbids the governor once sensed NG has drooped roughly 8 to 11
percentage points below the idle reference. It is implemented as printed anyway; see open
question #10.

There is no ECU here. `SPDG` is an input, and `SPDG_NULL` is the value at which the torque
motor sits still.

## Integration

Fixed-step explicit Euler on the lags, matching the report's real-time formulation and
CLAUDE.md's architecture constraint. Every state is named for the block that owns it.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from t700 import constants as engine_c
from t700.control import constants as c
from t700.control import schedules
from t700.control._blocks import backlash, clamp, delay, lag, lead_lag
from t700.units import wf_pps_from_pph

SPDG_NULL: float = -c.TM_INPUT_BIAS - c.TM_CURRENT_BIAS / c.TM_FORWARD_GAIN
"""The ECU trim signal at which the torque motor sits still. **Not 0.44.**

[Fig. C10, pdf p.89] The first summing junction forms `SPDG - 0.44`, and it is tempting to
call that the null. It is not: the forward path is `e * (0.04s+1)/(0.2s+1) * 564.0 - 31.0`,
so at `e = 0` it delivers **-31.0**, far outside the +/-2.0 deadband, and the integrator
ramps to its lower limit. The signal actually sits still where the forward path is zero:

    SPDG = 0.44 + 31.0/564.0 = 0.49496

which is derived from three printed constants and nothing else. That it lands within
0.01 of 0.5 is a small check on all three.

The deadband is what makes a null exist at all rather than a single exact value -- any
SPDG in [0.44 + 29/564, 0.44 + 33/564] = [0.4914, 0.4985] leaves a motor at rest at rest.
This constant is the centre of that band."""


@dataclass(frozen=True)
class HMUState:
    """Every continuous state and memory element in Figures C9-C22.

    Named for the block that owns it, so each can be traced back to a figure.
    """

    tm_leadlag: float = 0.0
    """Torque motor forward-path lead-lag [Fig. C10]."""
    tm_integrator: float = 0.0
    """Torque motor integrator, limited to [XLOLIM, XHILIM] [Fig. C10]."""
    ps3_lag: float = 0.0
    """Ps3 sensor lag, 1/(CTPS3*s + 1) [Fig. C11]."""
    ps3_hyst: float = 0.0
    """Ps3 sensor hysteresis memory, width PS3HYS [Fig. C11]."""
    load_demand_lag: float = 0.0
    """Load demand dynamics, 1/(CLLDS*s + 1) [Fig. C12]."""
    xldsh: float = 0.0
    """Collective rigging hysteresis memory, width XLDHYS [Fig. C13]."""
    pcng_hyst: float = 0.0
    """NG sensor hysteresis memory, width CH [Fig. C15]."""
    pcnghl: float = 0.0
    """NG sensor lag, 1/(CNTL*s + 1) [Fig. C15]. The sensed speed the whole HMU uses."""
    wfmv_lag: float = 0.0
    """Metering valve lag, 1/(CLMV*s + 1) [Fig. C19], in lbm/hr."""
    wf_history: tuple[float, ...] = ()
    """Metering valve transport delay buffer, e^(-0.015 s) [Fig. C19], newest last."""


@dataclass(frozen=True)
class HMUInputs:
    """What the HMU reads. `spdg` is the ECU trim; leave it at `SPDG_NULL` to run open."""

    ng_rpm: float
    ps3_psia: float
    t2_degR: float
    pas_deg: float
    """Power available spindle angle -- the cockpit power lever."""
    xcpc_pct: float
    """Helicopter collective pitch position, **percent of maximum** [nomenclature,
    pdf p.81, read off the raster]. Not inches -- at 0.914 deg per unit the rigging only
    reaches the load demand spindle's printed 0-100 deg range if this is a percentage,
    and `F_HM3`/`F_HM4` are drawn over 0-100 deg."""
    spdg: float = SPDG_NULL


@dataclass(frozen=True)
class HMUOutputs:
    """Named for the report's own signals, so the trace reads against the figures."""

    wf_pps: float
    wf_pph: float
    pcnghl: float
    ps3l_psia: float
    wfptp: float
    wfprf: float
    dwfp: float
    dwfpl: float
    tmru: float
    wfpdm: float
    wfidm: float
    wfpac: float
    wfpdc: float
    hmusel: float
    wfmv_pph: float
    limit: str
    """Which branch of the Fig. C18 cascade set `hmusel`: droop, idle, accel or decel."""


def seed(u: HMUInputs) -> HMUState:
    """A state consistent with holding `u` -- the HMU's analogue of `realtime.from_trim`.

    **The report does not say how the control initializes** (the same gap as open question
    #23 for the engine's opened iteration), and starting from zeros is not neutral here.
    Two blocks make that bite:

    * the torque motor's forward path is a **lead**-lag, so a state starting at zero sees
      a step and kicks; and
    * the integrator behind it sits in a **deadband**, which has a continuum of equilibria
      rather than one. A startup kick therefore leaves a permanent offset in `TMRU` --
      small, but it never washes out, because nothing pulls the integrator back once the
      current re-enters the band.

    So the sensible default is the one that makes a held input sit still: every lag at its
    input, the lead-lag at the value that makes its own derivative zero, the hysteresis
    memories at their inputs, and the integrator wherever the forward path is already
    inside the deadband -- which for `SPDG_NULL` is zero.
    """
    pcng = u.ng_rpm * 100.0 / engine_c.NG_DES
    xldsa = c.COLLECTIVE_RIGGING_GAIN * u.xcpc_pct + c.COLLECTIVE_RIGGING_BIAS
    png = float(schedules.f_hm3()(xldsa))
    wfqps3 = float(schedules.f_hm4()(xldsa))
    dwfp = c.KNDRP * (c.NGREF - png) + (c.LOAD_DEMAND_BIAS - wfqps3)
    # the metering valve has to be seeded too, or the loop starts with a fuel-flow step
    # from zero -- which is a slam acceleration the control never commanded
    wfptp = float(schedules.f_hm1()(u.t2_degR))
    wfprf = float(schedules.f_hm2()(u.pas_deg))
    tmru = 0.0  # the integrator is seeded at zero, and TMRU = TMLG * it
    wfpdm = wfptp + c.KNDRP * (c.NGREF - pcng) - wfprf - max(dwfp + tmru, 0.0)
    wfidm = c.KNDRP * (float(schedules.f_hm6()(u.t2_degR)) - pcng) - float(
        schedules.f_hm5()(u.t2_degR)
    )
    wfpac = float(schedules.f_hm7()(pcng, u.t2_degR))
    wfpdc = clamp(c.AWFP * pcng + c.BWFP, c.WFPDCL, c.WFPDCH)
    hmusel = max(min(max(wfpdm, wfidm), wfpac), wfpdc)
    wfmv = clamp(hmusel * u.ps3_psia, c.WFMIN, c.WFMAX)
    return HMUState(
        tm_leadlag=(u.spdg + c.TM_INPUT_BIAS),
        tm_integrator=0.0,
        ps3_lag=u.ps3_psia,
        ps3_hyst=u.ps3_psia,
        load_demand_lag=dwfp,
        xldsh=xldsa,
        pcng_hyst=pcng,
        pcnghl=pcng,
        wfmv_lag=wfmv,
        wf_history=(wfmv, wfmv, wfmv, wfmv),
    )


def step(state: HMUState, u: HMUInputs, dt: float) -> tuple[HMUState, HMUOutputs]:
    """Advance the HMU one frame. Returns the new state and the signals of Figs. C9-C22."""
    # --- C-36, Fig. C15: NG spool sensor -------------------------------------------
    pcng = u.ng_rpm * 100.0 / engine_c.NG_DES
    pcng_hyst = backlash(state.pcng_hyst, pcng, c.CH)
    pcnghl = lag(state.pcnghl, pcng_hyst, c.CNTL, dt)

    # --- C-32, Fig. C11: Ps3 sensor ------------------------------------------------
    ps3_lag = lag(state.ps3_lag, u.ps3_psia, c.CTPS3, dt)
    ps3l = backlash(state.ps3_hyst, ps3_lag, c.PS3HYS)

    # --- C-34, Fig. C13: collective to load demand spindle, UH-60A rigging ---------
    xldsa = c.COLLECTIVE_RIGGING_GAIN * u.xcpc_pct + c.COLLECTIVE_RIGGING_BIAS
    xldsh = backlash(state.xldsh, xldsa, c.XLDHYS)

    # --- C-38, Fig. C17: load demand spindle input schedule ------------------------
    # Fed XLDSH here, though Figs. C26/C27 label their x axes XLDSA -- see the
    # inventory's ambiguity 5. The figure wins: the wire carries XLDSH.
    png = float(schedules.f_hm3()(xldsh))
    wfqps3 = float(schedules.f_hm4()(xldsh))
    dwfp = c.KNDRP * (c.NGREF - png) + (c.LOAD_DEMAND_BIAS - wfqps3)

    # --- C-33, Fig. C12: load demand dynamics --------------------------------------
    load_demand_lag = lag(state.load_demand_lag, dwfp, c.CLLDS, dt)

    # --- C-28..C-31, Fig. C10: torque motor ----------------------------------------
    tm_err = (u.spdg + c.TM_INPUT_BIAS) - c.TMLVG * state.tm_integrator
    tm_leadlag, tm_out = lead_lag(state.tm_leadlag, tm_err, c.TM_LEAD, c.TM_LAG, dt)
    tm_current = tm_out * c.TM_FORWARD_GAIN + c.TM_CURRENT_BIAS
    if tm_current > c.TMDB:
        tm_db = tm_current - c.TMDB
    elif tm_current < -c.TMDB:
        tm_db = tm_current + c.TMDB
    else:
        tm_db = 0.0
    tm_integrator = clamp(state.tm_integrator + dt * c.TMGN * tm_db, c.XLOLIM, c.XHILIM)
    tmru = c.TMLG * tm_integrator

    # --- C-35/C-37, Figs. C14/C16: topping and power-available schedules ------------
    wfptp = float(schedules.f_hm1()(u.t2_degR))
    wfprf = float(schedules.f_hm2()(u.pas_deg))

    # --- C-26, Fig. C9: the droop line ---------------------------------------------
    demand_cut = max(load_demand_lag + tmru, 0.0)
    wfpdm = wfptp + c.KNDRP * (c.NGREF - pcnghl) - wfprf - demand_cut

    # --- C-41, Fig. C20: idle schedule ---------------------------------------------
    wfirf = float(schedules.f_hm5()(u.t2_degR))
    pcngi = float(schedules.f_hm6()(u.t2_degR))
    wfidm = c.KNDRP * (pcngi - pcnghl) - wfirf

    # --- C-42, Fig. C21: acceleration limit ----------------------------------------
    wfpac = float(schedules.f_hm7()(pcnghl, u.t2_degR))

    # --- C-43, Fig. C22: deceleration limit ----------------------------------------
    wfpdc = clamp(c.AWFP * pcnghl + c.BWFP, c.WFPDCL, c.WFPDCH)

    # --- C-39, Fig. C18: the limit cascade -----------------------------------------
    after_idle = max(wfpdm, wfidm)
    after_accel = min(after_idle, wfpac)
    hmusel = max(after_accel, wfpdc)
    if hmusel == wfpdc and wfpdc > after_accel:
        limit = "decel"
    elif after_accel == wfpac and wfpac < after_idle:
        limit = "accel"
    elif after_idle == wfidm and wfidm > wfpdm:
        limit = "idle"
    else:
        limit = "droop"

    # --- C-27/C-40, Figs. C9/C19: metering valve -----------------------------------
    wfmv = hmusel * ps3l
    wfmv_lag = lag(state.wfmv_lag, clamp(wfmv, c.WFMIN, c.WFMAX), c.CLMV, dt)
    history, wf_pph = delay(state.wf_history, wfmv_lag, c.FUEL_TRANSPORT_DELAY, dt)

    new = replace(
        state,
        tm_leadlag=tm_leadlag,
        tm_integrator=tm_integrator,
        ps3_lag=ps3_lag,
        ps3_hyst=ps3l,
        load_demand_lag=load_demand_lag,
        xldsh=xldsh,
        pcng_hyst=pcng_hyst,
        pcnghl=pcnghl,
        wfmv_lag=wfmv_lag,
        wf_history=history,
    )
    out = HMUOutputs(
        wf_pps=wf_pps_from_pph(wf_pph),
        wf_pph=wf_pph,
        pcnghl=pcnghl,
        ps3l_psia=ps3l,
        wfptp=wfptp,
        wfprf=wfprf,
        dwfp=dwfp,
        dwfpl=load_demand_lag,
        tmru=tmru,
        wfpdm=wfpdm,
        wfidm=wfidm,
        wfpac=wfpac,
        wfpdc=wfpdc,
        hmusel=hmusel,
        wfmv_pph=wfmv,
        limit=limit,
    )
    return new, out
