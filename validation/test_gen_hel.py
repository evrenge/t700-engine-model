"""Figures 12 and 14 [pdf pp.49, 52]: Gen Hel with our engine and fuel control in it, by
replaying its rotor speed.

Each figure is the closed-loop response of "engine, fuel control system, and blade-element
helicopter simulation" to one pilot input: Ballin's Gen Hel line (solid) against flight test
(dashed). Gen Hel itself is not in the report (#39), but its printed rotor speed is, and Gen Hel
drives the T700 through "a rigid-shaft drive train model" whose degree of freedom "corresponds
to that of the drive train and rotor hub" [pdf p.47]: rotor speed in % is NP in %. Everything
the airframe does reaches the propulsion system through that speed, and the flight-test inputs
are "single-axis control inputs" [pdf p.47], so a lateral cyclic or pedal input leaves the
collective, and with it the load demand, where it was ("without the influences of the load
demand compensation", pdf p.50, of Figure 14). So the comparison replays Ballin's
rotor speed into our closed loop as NP and predicts the three other things he prints, NG,
vehicle fuel flow and total torque, against his Gen Hel line (#81).

Readings, where the report is silent (#81): the day is standard; the trim is at his printed NG
and rotor speed, with the collective at which our HMU commands that trim's fuel flow, settled
in the closed loop with Appendix B's rotor load before the replay starts. On a standard day our
engine needs more fuel at his NG than his trims print, so the flight-test days were not
standard, and fuel and torque are compared in % of their own trim; total torque at the rotor
would also need the gearbox ratio, which the report does not print. Every difference allows the
figure's time read error.

His rotor speed carries a steady ripple of about 17 Hz, a few hundredths of a percent, in the
ink itself: Gen Hel's rotor, whose frequency is "the number of blades times the rotor
rotational frequency" [pdf p.50]. It is replayed as printed; smoothed out, it moves no Figure 14
score by more than a few percent of itself.
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
from test_closed_loop import SETTLE_FRAMES
from test_fig45 import _best, rms
from test_fuel_step import read_reference
from test_hmu_trim import collective_for

DT = realtime.FRAME_ENGINE_S
FIGS = (12, 14)
PANELS = ("ng", "fuel", "torque")
INPUT = {12: ("lat_cyclic", "lateral cyclic, in."), 14: ("pedals", "pedals, in.")}
"""Each figure's pilot input panel [pdf pp.49, 52]: drawn, never an input to our model."""
TRIM_BEFORE_S = {12: 2.2, 14: 2.6}
"""Every panel is flat until the input moves (about 2.3 s on Figure 12, 2.65 s on Figure 14):
the trim is the median before."""
TIME_TOL_S = {12: 0.0203, 14: 0.0353}
"""Each figure's time read error: its worst interior major tick over the five panels, as the
digitizer writes it in each CSV header."""


@functools.cache
def his(fig: int, key: str, line: str = "gen_hel") -> tuple[np.ndarray, np.ndarray]:
    """`(t, value)` of one of Ballin's lines, `gen_hel` (solid) or `test` (dashed)."""
    return read_reference(f"fig{fig:02d}_{key}_{line}.csv", "t_s", "value")


def his_trim(fig: int, key: str, line: str = "gen_hel") -> float:
    t, v = his(fig, key, line)
    return float(np.median(v[t < TRIM_BEFORE_S[fig]]))


def _np_trim_rpm(fig: int) -> float:
    return his_trim(fig, "nr") / 100.0 * engine_c.NP_DES


@functools.cache
def trim_fuel_pph(fig: int) -> float:
    """Our fuel flow, lbm/hr per engine, at his printed NG and rotor speed on a standard day."""
    target = his_trim(fig, "ng") / 100.0 * engine_c.NG_DES

    def miss(wf_pph: float) -> float:
        return trim.solve(wf_pps_from_pph(wf_pph), _np_trim_rpm(fig), AMB).state.ng_rpm - target

    a, b = 350.0, 550.0
    fa, fb = miss(a), miss(b)
    while abs(fb) > 1e-3:
        a, fa, b = b, fb, b - fb * (b - a) / (fb - fa)
        fb = miss(b)
    return b


@functools.cache
def our_trim(fig: int) -> trim.TrimResult:
    """Our open-loop trim at his printed NG and rotor speed, on a standard day."""
    return trim.solve(wf_pps_from_pph(trim_fuel_pph(fig)), _np_trim_rpm(fig), AMB)


@functools.cache
def collective(fig: int) -> float:
    """The collective, % of maximum, at which our HMU commands our trim's fuel flow."""
    r = our_trim(fig)
    ps3 = engine_c.K_PS3 * r.state.p3_psia
    return collective_for(r.state.ng_rpm, ps3, trim_fuel_pph(fig))[0]


@functools.cache
def replay(fig: int, np_held: bool = False) -> dict[str, np.ndarray]:
    """Our closed loop driven by his rotor speed: `t`, NG %, fuel and torque.

    Settled first with the load Appendix B's linear models carry about the trim (#6), then
    NP is set to his rotor speed every frame. `np_held` keeps NP at the trim instead: what the
    loop does without the rotor-speed input, the decoy."""
    r = our_trim(fig)
    np_trim = r.state.np_rpm
    q_trim = r.frame.q_pt_ftlbf
    pilot = loop.Pilot(xcpc_pct=collective(fig), pas_deg=100.0, pcprf_pct=his_trim(fig, "nr"))
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
    tn, nr = his(fig, "nr")
    t = np.arange(tn[0], tn[-1], DT)
    out = np.zeros((t.size, 3))
    for i, ti in enumerate(t):
        np_rpm = np_trim if np_held else float(np.interp(ti, tn, nr)) / 100.0 * engine_c.NP_DES
        s = dataclasses.replace(s, engine=dataclasses.replace(s.engine, np_rpm=np_rpm))
        s, _e, h, fr = loop.step(s, pilot, AMB, dt=DT, load=0.0, j_load=engine_c.J_LOAD_UH60A)
        out[i] = (100.0 * s.engine.ng_rpm / engine_c.NG_DES, h.wf_pph, fr.q_pt_ftlbf)
    return {"t": t, "ng": out[:, 0], "fuel": out[:, 1], "torque": out[:, 2]}


def ours(fig: int, key: str, np_held: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """`(t, value)`: NG in %, fuel and torque in % of our trim."""
    r = replay(fig, np_held)
    v = r[key]
    return r["t"], v if key == "ng" else 100.0 * v / v[0]


def gap(fig: int, key: str, line: str = "gen_hel", np_held: bool = False) -> np.ndarray:
    """Per digitized column of his line, how far ours is from it: NG in % of design speed,
    fuel and torque in % of trim."""
    tb, vb = his(fig, key, line)
    if key != "ng":
        vb = 100.0 * vb / his_trim(fig, key, line)
    t, v = ours(fig, key, np_held)
    keep = (tb >= t[0]) & (tb <= t[-1])
    return _best(t, v, tb[keep], vb[keep], np.zeros(keep.sum()), TIME_TOL_S[fig])


# ----------------------------------------------------------------------------- the checks


GEN_HEL_RMS = {
    12: {"ng": 0.04, "fuel": 0.5, "torque": 0.47},
    14: {"ng": 0.14, "fuel": 1.8, "torque": 2.1},
}
"""Ours against his Gen Hel line, rms: NG in % of design speed, fuel and torque in % of trim.
A characterization, 1.25 times today's value."""


@pytest.mark.parametrize("fig", FIGS)
@pytest.mark.parametrize("key", PANELS)
def test_ours_follows_his_gen_hel_line(fig: int, key: str):
    got = rms(gap(fig, key))
    assert got < GEN_HEL_RMS[fig][key], f"Figure {fig} {key}: rms {got:.3f}"


DECOY_MIN_RATIO = 1.5
"""The replay must be what does it: with NP held at the trim, ours sits at least this many
times further from his line."""


@pytest.mark.parametrize("fig", FIGS)
@pytest.mark.parametrize("key", PANELS)
def test_the_rotor_speed_drives_the_response(fig: int, key: str):
    assert rms(gap(fig, key, np_held=True)) > DECOY_MIN_RATIO * rms(gap(fig, key))


@pytest.mark.parametrize("fig", FIGS)
@pytest.mark.parametrize("key", PANELS)
def test_ours_is_nearer_gen_hel_than_the_flight_test(fig: int, key: str):
    """We replicate the simulation, not the aircraft: ours sits nearer Ballin's Gen Hel line
    than his flight-test line, which the test aircraft's greater load demand compensation
    [pdf p.47] and the collective trim mismatches [pdf pp.47, 50] move away from it."""
    assert rms(gap(fig, key)) < rms(gap(fig, key, "test"))
