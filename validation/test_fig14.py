"""Figure 14 [pdf p.52]: Gen Hel with our engine and fuel control in it, by replaying its rotor
speed.

Figure 14 is the closed-loop response of "engine, fuel control system, and blade-element
helicopter simulation to a pedal input at hover" [pdf p.52]: Ballin's Gen Hel line (solid)
against flight test (dashed). Gen Hel itself is not in the report (#39), but its printed rotor
speed is: everything the airframe does reaches the propulsion system through that speed, and a
pedal input leaves the collective, and with it the load demand, where it was ("without the
influences of the load demand compensation", pdf p.50). So the comparison replays Ballin's
rotor speed into our closed loop as NP and predicts the three other things he prints, NG,
vehicle fuel flow and total torque, against his Gen Hel line (#81).

Readings, where the report is silent (#81): main rotor speed in % is NP in % (a rigid gearbox);
the day is standard; the trim is at his printed NG, settled in the closed loop with Appendix
B's rotor load, then the replay starts. On a standard day our engine needs about a fifth more
fuel at his NG than his trim prints, so the flight-test day was not standard, and fuel and
torque are compared in % of their own trim; total torque at the rotor would also need the
gearbox ratio, which the report does not print. Every difference allows the figure's time
read error.

His rotor speed carries a steady ripple of about 17 Hz, a few hundredths of a percent, in the
ink itself: Gen Hel's rotor, whose frequency is "the number of blades times the rotor
rotational frequency" [pdf p.50]. It is replayed as printed; smoothed out, it moves no score
here by more than a few percent of itself.
"""

from __future__ import annotations

import dataclasses
import functools

import numpy as np
import pytest

from t700 import constants as engine_c
from t700 import realtime, trim
from t700.control import loop
from t700.units import wf_pps_from_pph
from test_all_dof_models import AMB, DQ_REQ_DNP
from test_closed_loop import SETTLE_FRAMES, XCPC_PCT
from test_fig45 import _best, rms
from test_fuel_step import read_reference

DT = realtime.FRAME_ENGINE_S
PANELS = ("ng", "fuel", "torque")
TRIM_BEFORE_S = 2.6
"""Every panel is flat until the pedal moves at about 2.65 s: the trim is the median before."""
TIME_TOL_S = 0.0353
"""The figure's time read error: its worst interior major tick over the five panels, as the
digitizer writes it in each CSV header."""


@functools.cache
def his(key: str, line: str = "gen_hel") -> tuple[np.ndarray, np.ndarray]:
    """`(t, value)` of one of Ballin's lines, `gen_hel` (solid) or `test` (dashed)."""
    return read_reference(f"fig14_{key}_{line}.csv", "t_s", "value")


def his_trim(key: str, line: str = "gen_hel") -> float:
    t, v = his(key, line)
    return float(np.median(v[t < TRIM_BEFORE_S]))


@functools.cache
def our_trim() -> trim.TrimResult:
    """Our open-loop trim at his printed NG and rotor speed, on a standard day."""
    np_rpm = his_trim("nr") / 100.0 * engine_c.NP_DES
    target = his_trim("ng") / 100.0 * engine_c.NG_DES

    def miss(wf_pph: float) -> float:
        return trim.solve(wf_pps_from_pph(wf_pph), np_rpm, AMB).state.ng_rpm - target

    a, b = 450.0, 550.0
    fa, fb = miss(a), miss(b)
    while abs(fb) > 1e-3:
        a, fa, b = b, fb, b - fb * (b - a) / (fb - fa)
        fb = miss(b)
    return trim.solve(wf_pps_from_pph(b), np_rpm, AMB)


@functools.cache
def replay(np_held: bool = False) -> dict[str, np.ndarray]:
    """Our closed loop driven by his rotor speed: `t`, NG %, fuel and torque.

    Settled first with the load Appendix B's linear models carry about the trim (#6), then
    NP is set to his rotor speed every frame. `np_held` keeps NP at the trim instead: what the
    loop does without the rotor-speed input, the decoy."""
    r = our_trim()
    np_trim = r.state.np_rpm
    q_trim = r.frame.q_pt_ftlbf
    pilot = loop.Pilot(xcpc_pct=XCPC_PCT["hover"], pas_deg=100.0, pcprf_pct=his_trim("nr"))
    s = loop.seed(r, pilot, AMB)
    for _ in range(SETTLE_FRAMES):
        s, *_ = loop.step(
            s,
            pilot,
            AMB,
            dt=DT,
            load=lambda v: q_trim + DQ_REQ_DNP[1] * (v - np_trim),
            j_load=engine_c.J_LOAD_UH60A,
        )
    tn, nr = his("nr")
    t = np.arange(tn[0], tn[-1], DT)
    out = np.zeros((t.size, 3))
    for i, ti in enumerate(t):
        np_rpm = np_trim if np_held else float(np.interp(ti, tn, nr)) / 100.0 * engine_c.NP_DES
        s = dataclasses.replace(s, engine=dataclasses.replace(s.engine, np_rpm=np_rpm))
        s, _e, h, fr = loop.step(s, pilot, AMB, dt=DT, load=0.0, j_load=engine_c.J_LOAD_UH60A)
        out[i] = (100.0 * s.engine.ng_rpm / engine_c.NG_DES, h.wf_pph, fr.q_pt_ftlbf)
    return {"t": t, "ng": out[:, 0], "fuel": out[:, 1], "torque": out[:, 2]}


def ours(key: str, np_held: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """`(t, value)`: NG in %, fuel and torque in % of our trim."""
    r = replay(np_held)
    v = r[key]
    return r["t"], v if key == "ng" else 100.0 * v / v[0]


def gap(key: str, line: str = "gen_hel", np_held: bool = False) -> np.ndarray:
    """Per digitized column of his line, how far ours is from it: NG in % of design speed,
    fuel and torque in % of trim."""
    tb, vb = his(key, line)
    if key != "ng":
        vb = 100.0 * vb / his_trim(key, line)
    t, v = ours(key, np_held)
    keep = (tb >= t[0]) & (tb <= t[-1])
    return _best(t, v, tb[keep], vb[keep], np.zeros(keep.sum()), TIME_TOL_S)


# ----------------------------------------------------------------------------- the checks


FIG14_RMS = {"ng": 0.14, "fuel": 1.9, "torque": 2.1}
"""Ours against his Gen Hel line, rms: NG in % of design speed, fuel and torque in % of trim.
A characterization, 1.25 times today's value."""


@pytest.mark.parametrize("key", PANELS)
def test_ours_follows_his_gen_hel_line(key: str):
    got = rms(gap(key))
    assert got < FIG14_RMS[key], f"Figure 14 {key}: rms {got:.3f}"


DECOY_MIN_RATIO = 1.5
"""The replay must be what does it: with NP held at the trim, ours sits at least this many
times further from his line."""


@pytest.mark.parametrize("key", PANELS)
def test_the_rotor_speed_drives_the_response(key: str):
    assert rms(gap(key, np_held=True)) > DECOY_MIN_RATIO * rms(gap(key))


@pytest.mark.parametrize("key", PANELS)
def test_ours_is_nearer_gen_hel_than_the_flight_test(key: str):
    """We replicate the simulation, not the aircraft: ours sits nearer Ballin's Gen Hel line
    than his flight-test line, which "the test aircraft contains greater load demand
    compensation than the simulation" [pdf p.47] and the hover collective mismatch [pdf p.50]
    move away from it."""
    assert rms(gap(key)) < rms(gap(key, "test"))
