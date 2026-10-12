"""Each Appendix C block driven where it acts, against the values the report prints.

For control paths no scenario reaches (limits never hit, thresholds never crossed, lags no
assertion looks at), found by `tools/constant_sweep.py`. A check against `c.WFMAX` moves with
a wrong constant, so these compare with `PRINTED`, read from the page images independently
of `t700.control.constants`. Each test drives one block where it acts: a limit into its
limit, a threshold to either side, a lag through a step against the recursion `_blocks`
documents. They test the block diagram, not an engine: no printed response isolates these
paths (Figures 11-15 need Gen Hel). `UNCONSTRAINED` records what cannot be pinned this way.
"""

from __future__ import annotations

import dataclasses

import pytest

from t700 import constants as engine_c
from t700.control import constants as c
from t700.control import ecu, hmu

DT = 0.007
"""The report's engine frame [pdf p.47]."""

PRINTED = {
    # Table C.1 [pdf p.83-84]
    "AWFP": 0.05909, "BWFP": -3.927, "CB": 0.75, "CE": 13.21, "CH": 0.05, "CLLDS": 0.2,
    "CLMV": 0.03, "CNTL": 0.025, "CORR": 20.0, "CR": 20.0, "CT2": 0.088, "CT7": 1.0,
    "CT13": 1.0, "CT14": 0.10, "DBIAS": 3.25, "PS3HYS": 0.375, "T8": 0.40, "T10": 0.060,
    "T45COR": 11.0, "TLGE": 0.077, "TMDB": 2.0, "TMGN": 0.0159, "TMLG": 84.0, "TMLVG": 55.0,
    "WFMAX": 785.0, "WFMIN": 65.0, "WFPDCH": 2.10, "WFPDCL": 1.45, "XHILIM": 0.0787,
    "XLDHYS": 2.5, "XLOLIM": -0.0427, "YHILIM": 40.0, "YLOLIM": 0.0, "ZHILIM": 3.5,
    "ZK1": 1.7, "ZK5": 0.40, "ZK7": 0.30, "ZK8": 0.231, "ZK9": 0.625, "ZK10": 0.375,
    "ZLOLIM": -1.0,
    # block diagrams: Fig. C2 [p.85], C3 [p.86], C8 [p.88], C10 [p.89], C19 [p.92]
    "LOAD_SHARE_LAG_GAIN": 1.1, "NP_LOOP_RELAY_LOW": -1.0, "NP_LOOP_RELAY_HIGH": 4.0,
    "ENGINE2_TORQUE_THRESHOLD": 180.0, "ZLOLIM_HIGH_TORQUE": -0.3, "TM_INPUT_BIAS": -0.44,
    "TM_LEAD": 0.04, "TM_LAG": 0.2, "TM_FORWARD_GAIN": 564.0, "TM_CURRENT_BIAS": -31.0,
    "FUEL_TRANSPORT_DELAY": 0.015,
    # the relay chain of UNCONSTRAINED below: B6 from Table C.1 [p.83]; the 1000 gain and
    # both relays' 1 from Fig. C3 [p.86]
    "B6": 1.0, "NP_LOOP_RELAY_GAIN": 1000.0, "RELAY_OUTPUT": 1.0,
}  # fmt: skip
P = PRINTED


def test_every_printed_value_is_the_constant():
    """`PRINTED` was read off the page images independently of the constants module, so
    agreement is a check on both readings. This is what pins each value; the block tests
    below pin where it acts."""
    wrong = {n: (v, getattr(c, n, None)) for n, v in PRINTED.items() if getattr(c, n, None) != v}
    assert not wrong, f"constant differs from the printed value, as (printed, module): {wrong}"


UNCONSTRAINED = {
    "CT9": "10 ms sensor lag; at dt/tau = 0.7 its whole effect is under 0.1 % (_blocks.lag)",
    "CTPL": "10 ms sensor lag, as CT9",
    "CTPS3": "10 ms sensor lag, as CT9",
    "TL1": "10 ms sensor lag, as CT9",
    "TL2": "10 ms sensor lag, as CT9",
    "T17": "10 ms lag in the load-share path, as CT9",
    "NP_LOOP_RELAY_GAIN": "Fig. C3: 1000 x B6 x relay only has to exceed CORR = 20; any product "
    "above 20 gives the same B4 (test_the_relay_gain_only_has_to_beat_corr)",
    "B6": "as NP_LOOP_RELAY_GAIN",
    "RELAY_OUTPUT": "as NP_LOOP_RELAY_GAIN",
    "K_T49_1": "T49 is computed [Eq. 36] and read by nothing, and the report prints no T49",
    "K_T49_2": "as K_T49_1",
}
"""Pinned by value only, on purpose: the 10 ms lags by `test_control_constants.py`, the relay
chain by `PRINTED` above, `K_T49_*` by `test_constants.py`."""


def _lag(x: float, u: float, tau: float, dt: float = DT) -> float:
    """`_blocks.lag`'s forward-Euler step, restated so the printed tau goes in."""
    return x + dt * (u - x) / tau


def _hmu(n: int, u0: hmu.HMUInputs, u1: hmu.HMUInputs | None = None, dt: float = DT):
    """Seed on `u0`, then run `n` frames on `u1`; every (state, output) pair."""
    s, trace = hmu.seed(u0), []
    for _ in range(n):
        s, o = hmu.step(s, u1 or u0, dt)
        trace.append((s, o))
    return trace


def _ecu(n: int, u0: ecu.ECUInputs, u1: ecu.ECUInputs | None = None, state=None):
    s, trace = (state or ecu.seed(u0)), []
    for _ in range(n):
        s, o = ecu.step(s, u1 or u0, DT)
        trace.append((s, o))
    return trace


HMU_BASE = hmu.HMUInputs(ng_rpm=41638.0, ps3_psia=176.34, t2_degR=518.67, pas_deg=100.0,
                         xcpc_pct=50.0)  # fmt: skip
ECU_BASE = ecu.ECUInputs(np_rpm=engine_c.NP_DES, torq45_ftlbf=229.0, t45_degR=1632.0,
                         w45_pps=7.90, p45_psia=37.42, pcprf_pct=100.0)  # fmt: skip


def _ng(pcng: float) -> float:
    return pcng * engine_c.NG_DES / 100.0


# ------------------------------------------------------------------- HMU, Figs. C9-C22


def test_the_metering_valve_saturates_at_the_printed_limits():
    """[Fig. C19] WFMV = HMUSEL x PS3L, clamped to [WFMIN, WFMAX] before the CLMV lag."""
    hi = _hmu(600, dataclasses.replace(HMU_BASE, ps3_psia=600.0))[-1][1]
    lo = _hmu(600, dataclasses.replace(HMU_BASE, ps3_psia=5.0))[-1][1]
    assert hi.wfmv_pph > P["WFMAX"] and lo.wfmv_pph < P["WFMIN"], "not driven into the clamp"
    assert hi.wf_pph == pytest.approx(P["WFMAX"], rel=1e-12)
    assert lo.wf_pph == pytest.approx(P["WFMIN"], rel=1e-12)


@pytest.mark.parametrize("pcng, where", [(95.0, "line"), (110.0, "WFPDCH"), (50.0, "WFPDCL")])
def test_the_deceleration_line_and_its_printed_clamps(pcng: float, where: str):
    """[Fig. C22] WFPDC = clamp(AWFP x PCNGHL + BWFP, WFPDCL, WFPDCH)."""
    o = _hmu(1, dataclasses.replace(HMU_BASE, ng_rpm=_ng(pcng)))[-1][1]
    line = P["AWFP"] * pcng + P["BWFP"]
    want = {"line": line, "WFPDCH": P["WFPDCH"], "WFPDCL": P["WFPDCL"]}[where]
    assert P["WFPDCL"] < line < P["WFPDCH"] or where != "line"
    assert o.wfpdc == pytest.approx(want, rel=1e-12), where


def _tm_null_and_band() -> tuple[float, float]:
    """[Fig. C10] The SPDG at which the motor current is zero, and the half-width of the
    TMDB deadband expressed in SPDG."""
    null = -P["TM_INPUT_BIAS"] - P["TM_CURRENT_BIAS"] / P["TM_FORWARD_GAIN"]
    return null, P["TMDB"] / P["TM_FORWARD_GAIN"]


@pytest.mark.parametrize("k", [-1.03, -0.97, 0.97, 1.03])
def test_the_torque_motor_deadband_sits_where_the_printed_biases_put_it(k: float):
    """Inside the band the integrator never moves; three per cent outside it, it does, and
    in the direction of the error. Moving the null or the band by 5 % flips one of these."""
    null, half = _tm_null_and_band()
    u = dataclasses.replace(HMU_BASE, spdg=null + k * half)
    s = _hmu(2000, u)[-1][0]
    if abs(k) < 1.0:
        assert s.tm_integrator == 0.0, (k, s.tm_integrator)
    else:
        assert s.tm_integrator * k > 1e-9, (k, s.tm_integrator)


@pytest.mark.parametrize("spdg, limit", [(5.0, "XHILIM"), (-5.0, "XLOLIM")])
def test_the_torque_motor_saturates_at_the_printed_limits(spdg: float, limit: str):
    """[Fig. C10] TMGN/s, then a limiter at XHILIM / XLOLIM, then TMLG to TMRU."""
    s, o = _hmu(4000, dataclasses.replace(HMU_BASE, spdg=spdg))[-1]
    assert s.tm_integrator == pytest.approx(P[limit], rel=1e-12)
    assert o.tmru == pytest.approx(P["TMLG"] * P[limit], rel=1e-12)


def test_the_torque_motor_first_frames_are_the_printed_loop():
    """Two frames after a step in SPDG from the null, run through the printed Fig. C10 loop:
    bias, lead-lag (0.04s+1)/(0.2s+1), gain 564, bias -31, deadband TMDB, TMGN/s, and the
    TMLVG feedback. (Settled, the integrator rests anywhere inside the deadband: its loop
    gain is about 3.4 per 7 ms frame, so it overshoots and stops where it lands.)"""
    null, half = _tm_null_and_band()
    step = null + 20.0 * half
    u0, u1 = dataclasses.replace(HMU_BASE, spdg=null), dataclasses.replace(HMU_BASE, spdg=step)
    trace = _hmu(2, u0, u1)
    x, i_ = null + P["TM_INPUT_BIAS"], 0.0
    for s, _o in trace:
        dx = ((step + P["TM_INPUT_BIAS"] - P["TMLVG"] * i_) - x) / P["TM_LAG"]
        current = (x + P["TM_LEAD"] * dx) * P["TM_FORWARD_GAIN"] + P["TM_CURRENT_BIAS"]
        drive = current - P["TMDB"] if current > P["TMDB"] else min(current + P["TMDB"], 0.0)
        i_ = min(max(i_ + DT * P["TMGN"] * drive, P["XLOLIM"]), P["XHILIM"])
        x = x + DT * dx
        assert s.tm_integrator == pytest.approx(i_, rel=1e-12)
    assert trace[-1][0].tm_integrator != trace[0][0].tm_integrator


def test_the_ps3_sensor_hysteresis_is_the_printed_band():
    """[Fig. C11] PS3L follows a rise only after it has moved by half of PS3HYS."""
    up = _hmu(400, HMU_BASE, dataclasses.replace(HMU_BASE, ps3_psia=HMU_BASE.ps3_psia + 2.0))
    want = HMU_BASE.ps3_psia + 2.0 - P["PS3HYS"] / 2
    assert up[-1][1].ps3l_psia == pytest.approx(want, rel=1e-12)


def test_the_collective_rigging_hysteresis_is_the_printed_band():
    """[Fig. C13] XLDSH lags XLDSA by half of XLDHYS after a rise in collective."""
    s = _hmu(10, HMU_BASE, dataclasses.replace(HMU_BASE, xcpc_pct=60.0))[-1][0]
    xldsa = c.COLLECTIVE_RIGGING_GAIN * 60.0 + c.COLLECTIVE_RIGGING_BIAS
    assert s.xldsh == pytest.approx(xldsa - P["XLDHYS"] / 2, rel=1e-12)


def test_the_ng_sensor_lag_is_cntl():
    """[Fig. C15] PCNGHL = 1/(CNTL s + 1) of the hysteresis output, frame by frame."""
    u1 = dataclasses.replace(HMU_BASE, ng_rpm=_ng(95.0))
    trace = _hmu(12, HMU_BASE, u1)
    x = HMU_BASE.ng_rpm * 100.0 / engine_c.NG_DES
    for s, o in trace:
        x = _lag(x, s.pcng_hyst, P["CNTL"])
        assert o.pcnghl == pytest.approx(x, rel=1e-12)


def test_the_load_demand_lag_is_cllds():
    """[Fig. C12] DWFPL = 1/(CLLDS s + 1) of DWFP."""
    trace = _hmu(60, HMU_BASE, dataclasses.replace(HMU_BASE, xcpc_pct=60.0))
    x = hmu.seed(HMU_BASE).load_demand_lag
    for _s, o in trace:
        x = _lag(x, o.dwfp, P["CLLDS"])
        assert o.dwfpl == pytest.approx(x, rel=1e-12)


def test_the_metering_valve_lag_and_transport_delay_are_the_printed_ones():
    """[Fig. C19] 1/(CLMV s + 1), then e^(-0.015 s). At dt = 5 ms the delay is exactly three
    frames, so the fuel flow out is the lagged demand from three frames earlier."""
    dt = 0.005
    u1 = dataclasses.replace(HMU_BASE, ps3_psia=HMU_BASE.ps3_psia + 3.0)
    trace = _hmu(40, HMU_BASE, u1, dt=dt)
    x = hmu.seed(HMU_BASE).wfmv_lag
    lagged = [x] * 4
    for s, o in trace:
        x = _lag(x, min(max(o.wfmv_pph, P["WFMIN"]), P["WFMAX"]), P["CLMV"], dt)
        assert s.wfmv_lag == pytest.approx(x, rel=1e-12)
        lagged.append(x)
        frames = round(P["FUEL_TRANSPORT_DELAY"] / dt)
        assert o.wf_pph == pytest.approx(lagged[-1 - frames], rel=1e-9)


# ------------------------------------------------------------------- ECU, Figs. C1-C8


@pytest.mark.parametrize("trqer", [2.0, 5.0, 20.0])
def test_the_load_share_trim_is_the_printed_gain_and_clamp(trqer: float):
    """[Fig. C2] XQLO = ZK8 x (clamp(1.1 x TRQER, DBIAS, CE) - DBIAS) at steady state --
    below DBIAS, between, and saturated at CE."""
    u = dataclasses.replace(ECU_BASE, engine2_torq45_ftlbf=ECU_BASE.torq45_ftlbf + trqer)
    o = _ecu(3000, u)[-1][1]
    inner = min(max(P["LOAD_SHARE_LAG_GAIN"] * trqer, P["DBIAS"]), P["CE"])
    assert o.xqlo == pytest.approx(P["ZK8"] * (inner - P["DBIAS"]), abs=1e-9)


def test_the_load_share_first_lag_is_ct7():
    """[Fig. C2] 1.1/(CT7 s + 1), frame by frame from a step in the other engine's torque."""
    u1 = dataclasses.replace(ECU_BASE, engine2_torq45_ftlbf=ECU_BASE.torq45_ftlbf + 8.0)
    x = 0.0
    for s, o in _ecu(200, ECU_BASE, u1):
        x = _lag(x, P["LOAD_SHARE_LAG_GAIN"] * (u1.engine2_torq45_ftlbf - o.trql), P["CT7"])
        assert s.xqlo_lag1 == pytest.approx(x, rel=1e-12)


def _at_spder(spder: float, **over) -> ecu.ECUInputs:
    return dataclasses.replace(ECU_BASE, pcprf_pct=100.0 - spder, **over)


@pytest.mark.parametrize("spder, on", [(-1.03, 1.0), (-0.97, 0.0), (3.9, 0.0), (4.1, 1.0)])
def test_the_relay_switches_at_the_printed_thresholds(spder: float, on: float):
    """[Fig. C3] The relay is 1 below -1 and above 4 %NP, 0 between; with the torque
    integrator held at zero (an unloaded engine), B4 is the relay."""
    o = _ecu(400, _at_spder(spder, torq45_ftlbf=5.0))[-1][1]
    assert o.spder == pytest.approx(spder, abs=1e-9)
    assert o.b4 == on, (spder, o.b4)


def test_the_torque_integrator_crosses_corr_when_the_printed_numbers_say():
    """[Fig. C3] Y = integral of (TRQL - CR) from zero; B4 turns on when Y passes CORR. With
    TRQL - CR = 10 that is the first frame after 20/10 = 2 s."""
    u = dataclasses.replace(ECU_BASE, torq45_ftlbf=P["CR"] + 10.0)
    first_on = next(i for i, (_s, o) in enumerate(_ecu(400, u), start=1) if o.b4 == 1.0)
    want = int(P["CORR"] / (10.0 * DT)) + 1
    assert first_on == want, (first_on, want)


@pytest.mark.parametrize("torque, limit", [(400.0, "YHILIM"), (5.0, "YLOLIM")])
def test_the_torque_integrator_limits(torque: float, limit: str):
    s = _ecu(2000, dataclasses.replace(ECU_BASE, torq45_ftlbf=torque))[-1][0]
    assert s.torque_int == pytest.approx(P[limit], abs=1e-12)


@pytest.mark.parametrize(
    "spder, engine2, limit",
    [(8.0, None, "ZHILIM"), (-8.0, None, "ZLOLIM"), (-8.0, 179.9, "ZLOLIM"),
     (-8.0, 180.1, "ZLOLIM_HIGH_TORQUE")],
)  # fmt: skip
def test_the_pi_integrator_limits_and_their_switch(spder, engine2, limit):
    """[Figs. C7, C8] The P+I integrator stops at ZHILIM, and at ZLOLIM -- which is -0.3
    instead of -1.0 once the other engine's torque reaches 180."""
    u = _at_spder(spder, engine2_torq45_ftlbf=engine2)
    s = _ecu(3000, u)[-1][0]
    assert s.pi_int == pytest.approx(P[limit], abs=1e-12)


@pytest.mark.parametrize("spder, torque", [(2.0, 400.0), (2.0, 5.0), (-0.5, 400.0)])
def test_the_governor_steady_gain_is_the_printed_one(spder: float, torque: float):
    """[Fig. C3] Settled, the rate lags agree and the pseudo-derivative is zero, so
    SPDSF = (ZK7 + B4 ZK1)(ZK9 SPDER + ZK10 deadband(SPDER, CB))."""
    o = _ecu(4000, _at_spder(spder, torq45_ftlbf=torque))[-1][1]
    db = spder - P["CB"] if spder > P["CB"] else (spder + P["CB"] if spder < -P["CB"] else 0.0)
    want = (P["ZK7"] + o.b4 * P["ZK1"]) * (P["ZK9"] * spder + P["ZK10"] * db)
    assert o.spdsf == pytest.approx(want, rel=1e-9)


def test_the_governor_rate_lags_and_pseudo_derivative():
    """[Fig. C3] CT13 then CT14 on the deadbanded error, and ZK5/CT14 times their difference
    added to the proportional terms. SPDS1 is read back through the CT2 lead-lag state."""
    u0, u1 = _at_spder(0.0, torq45_ftlbf=400.0), _at_spder(2.0, torq45_ftlbf=400.0)
    s0 = _ecu(3000, u0)[-1][0]
    x13, x14, prev = s0.rate_lag13, s0.rate_lag14, s0
    db = 2.0 - P["CB"]
    for s, o in _ecu(80, u0, u1, state=s0):
        x13 = _lag(x13, db, P["CT13"])
        x14 = _lag(x14, x13, P["CT14"])
        assert s.rate_lag13 == pytest.approx(x13, rel=1e-12)
        assert s.rate_lag14 == pytest.approx(x14, rel=1e-12)
        u_comp = P["ZK9"] * 2.0 + P["ZK10"] * db + P["ZK5"] / P["CT14"] * (x13 - x14)
        spds1 = prev.gov_leadlag + (s.gov_leadlag - prev.gov_leadlag) * P["CT2"] / DT
        assert spds1 == pytest.approx((P["ZK7"] + o.b4 * P["ZK1"]) * u_comp, rel=1e-9)
        prev = s


def test_the_thermocouple_harness_lag_is_tlge():
    """[Fig. C5] T45L = 1/(TLGE s + 1) of T4.5 + T45COR."""
    u1 = dataclasses.replace(ECU_BASE, t45_degR=1700.0)
    x = ECU_BASE.t45_degR + P["T45COR"]
    for s, _o in _ecu(60, ECU_BASE, u1):
        x = _lag(x, 1700.0 + P["T45COR"], P["TLGE"])
        assert s.harness_lag == pytest.approx(x, rel=1e-12)


def test_the_temperature_compensation_lead_lag_is_t8_over_t10():
    """[Fig. C6] (T8 s + 1)/(T10 s + 1) after the CT9 lag, frame by frame."""
    u1 = dataclasses.replace(ECU_BASE, t45_degR=1900.0)
    s0 = _ecu(500, ECU_BASE)[-1][0]
    x = s0.t45comp_leadlag
    for s, o in _ecu(60, ECU_BASE, u1, state=s0):
        dx = (s.t45comp_lag - x) / P["T10"]
        assert o.tsig == pytest.approx(x + P["T8"] * dx, rel=1e-12)
        x = x + DT * dx
        assert s.t45comp_leadlag == pytest.approx(x, rel=1e-12)


# ------------------------------------------------------------------ the unconstrained


def test_the_unconstrained_record_names_real_constants():
    from t700 import constants as engine

    for name in UNCONSTRAINED:
        assert hasattr(c, name) or hasattr(engine, name), name


def test_the_relay_gain_only_has_to_beat_corr(monkeypatch):
    """Why NP_LOOP_RELAY_GAIN, B6 and RELAY_OUTPUT are on the record rather than pinned: B4
    is a comparison against CORR, so the relay chain's product is invisible above 20."""
    u = _at_spder(5.0, torq45_ftlbf=5.0)
    seen = set()
    for gain in (25.0, 1000.0, 4000.0):
        monkeypatch.setattr(c, "NP_LOOP_RELAY_GAIN", gain)
        seen.add(_ecu(50, u)[-1][1].b4)
    assert seen == {1.0}
