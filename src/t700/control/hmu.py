"""The T700 hydromechanical control unit -- Appendix C, Figures C9 through C22.

The HMU is the primary fuel control loop: it governs gas generator speed on a droop line,
applies the idle, acceleration and deceleration limits, and meters fuel. The ECU only trims
it, through the torque motor signal `SPDG`; held at `SPDG_NULL` the HMU runs alone.

[Fig. C9, pdf p.89]

    WFPDM = WFPTP + KNDRP*(NGREF - PCNGHL) - WFPRF - max(DWFPL + TMRU, 0)

`WFPTP`, the topping line, is the ceiling; everything else subtracts from it. Reading a
summing junction backwards inverts the whole control, so the signs are read at 600 dpi and
recorded in `docs/notes/inventory-appendix-c.md`. Everything on the droop line is a fuel
flow parameter, Wf/Ps3, which becomes a fuel flow only at the metering valve,
`WFMV = HMUSEL * PS3L`.

[Fig. C18, pdf p.92]

    HMUSEL = MAX( MIN( MAX(WFPDM, WFIDM), WFPAC ), WFPDC )

Idle floor, then acceleration ceiling, then deceleration floor. The order is drawn, and the
deceleration floor outranks the acceleration ceiling.

The real-time model eliminated fuel control below flight idle [pdf p.38], so the idle floor
`WFIDM` is implemented as printed but inactive wherever the engine model is valid (open
question #10). Drooping NG never brings it into play: `WFPDM` and `WFIDM` carry the same
gain `KNDRP` on sensed NG, so their difference does not depend on NG; T2, PAS and the
collective decide it
(`tests/test_hmu.py::test_the_droop_line_and_the_idle_floor_are_parallel_in_ng`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

from t700 import constants as engine_c
from t700.control import constants as c
from t700.control import schedules
from t700.control._blocks import backlash, clamp, deadband, delay, integrator, lag, lead_lag
from t700.units import wf_pps_from_pph

SPDG_NULL: float = -c.TM_INPUT_BIAS - c.TM_CURRENT_BIAS / c.TM_FORWARD_GAIN
"""The ECU trim signal at which the torque motor sits still. Not 0.44.

[Fig. C10, pdf p.89] The first summing junction forms `SPDG - 0.44`, but the forward path
`e * (0.04s+1)/(0.2s+1) * 564.0 - 31.0` gives -31.0 at `e = 0`, far outside the +/-2.0
deadband. It is zero at

    SPDG = 0.44 + 31.0/564.0 = 0.49496

the centre of the deadband's band of rest, [0.44 + 29/564, 0.44 + 33/564]."""


@dataclass(frozen=True)
class HMUState:
    """Every continuous state and memory element in Figures C9-C22, named for its block."""

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
    """Helicopter collective pitch position, percent of maximum [nomenclature, pdf p.81,
    read off the raster]. Not inches: only a percentage, at 0.914 deg per unit, spans the
    0-100 deg over which `F_HM3`/`F_HM4` are drawn."""
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


def _load_demand(xlds: float) -> float:
    """DWFP from the load demand spindle, Fig. C17 [C-38]."""
    png = float(schedules.f_hm3()(xlds))
    wfqps3 = float(schedules.f_hm4()(xlds))
    return c.KNDRP * (c.NGREF - png) + (c.LOAD_DEMAND_BIAS - wfqps3)


class _Demands(NamedTuple):
    wfptp: float
    wfprf: float
    wfpdm: float
    wfidm: float
    wfpac: float
    wfpdc: float
    hmusel: float
    limit: str


def _demands(pcnghl: float, t2_degR: float, pas_deg: float, dwfpl: float, tmru: float) -> _Demands:
    """The droop line, the three limits and the Fig. C18 cascade between them."""
    # --- C-35/C-37, Figs. C14/C16: topping and power-available schedules ------------
    wfptp = float(schedules.f_hm1()(t2_degR))
    wfprf = float(schedules.f_hm2()(pas_deg))

    # --- C-26, Fig. C9: the droop line ---------------------------------------------
    demand_cut = max(dwfpl + tmru, 0.0)
    wfpdm = wfptp + c.KNDRP * (c.NGREF - pcnghl) - wfprf - demand_cut

    # --- C-41, Fig. C20: idle schedule ---------------------------------------------
    wfirf = float(schedules.f_hm5()(t2_degR))
    pcngi = float(schedules.f_hm6()(t2_degR))
    wfidm = c.KNDRP * (pcngi - pcnghl) - wfirf

    # --- C-42, Fig. C21: acceleration limit ----------------------------------------
    wfpac = float(schedules.f_hm7()(pcnghl, t2_degR))

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
    return _Demands(wfptp, wfprf, wfpdm, wfidm, wfpac, wfpdc, hmusel, limit)


def seed(u: HMUInputs) -> HMUState:
    """A state consistent with holding `u` -- the HMU's analogue of `realtime.from_trim`.

    The report does not say how the control initializes (open question #55), and zeros are
    not neutral: the torque motor's lead-lag kicks on a step, and the integrator behind it
    sits in a deadband with a continuum of equilibria, so a startup kick leaves a permanent
    `TMRU` offset. Hence every lag and hysteresis memory at its input, the lead-lag at zero
    derivative, and the integrator at zero, inside the deadband for `SPDG_NULL`.
    """
    pcng = u.ng_rpm * 100.0 / engine_c.NG_DES
    xldsa = c.COLLECTIVE_RIGGING_GAIN * u.xcpc_pct + c.COLLECTIVE_RIGGING_BIAS
    dwfp = _load_demand(xldsa)
    # seed the metering valve too, or the loop starts with an uncommanded fuel-flow step
    tmru = 0.0  # the integrator is seeded at zero, and TMRU = TMLG * it
    hmusel = _demands(pcng, u.t2_degR, u.pas_deg, dwfp, tmru).hmusel
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
    # Fed XLDSH, as the Fig. C17 wire carries, though Figs. C26/C27 label their x axes
    # XLDSA (open question #12; the inventory's ambiguity 5).
    dwfp = _load_demand(xldsh)

    # --- C-33, Fig. C12: load demand dynamics --------------------------------------
    load_demand_lag = lag(state.load_demand_lag, dwfp, c.CLLDS, dt)

    # --- C-28..C-31, Fig. C10: torque motor ----------------------------------------
    tm_err = (u.spdg + c.TM_INPUT_BIAS) - c.TMLVG * state.tm_integrator
    tm_leadlag, tm_out = lead_lag(state.tm_leadlag, tm_err, c.TM_LEAD, c.TM_LAG, dt)
    tm_current = tm_out * c.TM_FORWARD_GAIN + c.TM_CURRENT_BIAS
    tm_db = deadband(tm_current, c.TMDB)
    tm_integrator = integrator(state.tm_integrator, tm_db, c.TMGN, c.XLOLIM, c.XHILIM, dt)
    tmru = c.TMLG * tm_integrator

    # --- C-26 and C-39..C-43, Figs. C9 and C18-C22: droop line, limits, cascade ------
    d = _demands(pcnghl, u.t2_degR, u.pas_deg, load_demand_lag, tmru)
    hmusel = d.hmusel

    # --- C-27/C-40, Figs. C9/C19: metering valve -----------------------------------
    wfmv = hmusel * ps3l
    wfmv_lag = lag(state.wfmv_lag, clamp(wfmv, c.WFMIN, c.WFMAX), c.CLMV, dt)
    history, wf_pph = delay(state.wf_history, wfmv_lag, c.FUEL_TRANSPORT_DELAY, dt)

    # All ten fields are rewritten -- see the note in `ecu.step`.
    new = HMUState(
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
        wfptp=d.wfptp,
        wfprf=d.wfprf,
        dwfp=dwfp,
        dwfpl=load_demand_lag,
        tmru=tmru,
        wfpdm=d.wfpdm,
        wfidm=d.wfidm,
        wfpac=d.wfpac,
        wfpdc=d.wfpdc,
        hmusel=hmusel,
        wfmv_pph=wfmv,
        limit=d.limit,
    )
    return new, out
