"""Figures 11-14 [pdf pp.48-52]: Gen Hel with our engine and fuel control in it, by replaying
its rotor speed and collective.

Each figure is the closed-loop response of "engine, fuel control system, and blade-element
helicopter simulation" to one pilot input: Ballin's Gen Hel line (solid) against flight test
(dashed). Gen Hel itself is not in the report (#39), but its printed rotor speed is, and Gen Hel
drives the T700 through "a rigid-shaft drive train model" whose degree of freedom "corresponds
to that of the drive train and rotor hub" [pdf p.47]: rotor speed in % is NP in %. Everything
the airframe does reaches the propulsion system through that speed, except the collective,
which the fuel control reads itself (Fig. C13, pdf p.90). The flight-test inputs are
"single-axis control inputs" [pdf p.47]: Figures 12 and 14's lateral cyclic and pedals leave the
collective where it was ("without the influences of the load demand compensation", pdf p.50, of
Figure 14), and Figures 11 and 13 move it. So the comparison replays Ballin's rotor speed into
our closed loop as NP, and his collective where it moves, and predicts the three other things he
prints, NG, vehicle fuel flow and total torque, against his Gen Hel line (#81).

His collective is printed in inches, ours is "percent of maximum" [pdf p.81], and the report
prints the travel in neither. Read proportionally, 0 in. as 0 %, against the collective our HMU
needs at each figure's trim, it is 8.8 %/in. on Figure 11 and 9.6 %/in. on Figure 13; the
straight line through both trims instead gives 13.5 %/in. Each figure is compared at its own
proportional reading; `collective_scales` spans all three, and the report shows that range.

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
from t700.units import RPM_TO_RAD_PER_SEC, wf_pps_from_pph
from test_all_dof_models import AMB, DQ_REQ_DNP
from test_closed_loop import SETTLE_FRAMES
from test_fig45 import _best, rms
from test_fuel_step import read_reference
from test_hmu_trim import collective_for

DT = realtime.FRAME_ENGINE_S
FIGS = (11, 12, 13, 14)
COLLECTIVE_FIGS = (11, 13)
"""The two collective inputs: there the collective is an input to our fuel control too."""
PANELS = ("ng", "fuel", "torque")
INPUT = {
    11: ("collective", "collective, in.", (0.0, 10.0)),
    12: ("lat_cyclic", "lateral cyclic, in.", (-5.0, 5.0)),
    13: ("collective", "collective, in.", (0.0, 10.0)),
    14: ("pedals", "pedals, in.", (-3.0, 3.0)),
}
"""Each figure's pilot input panel and its printed range [pdf pp.48, 49, 51, 52]. The lateral
cyclic and the pedals are drawn only; the collective also drives our fuel control."""
TRIM_BEFORE_S = {11: 3.2, 12: 2.2, 13: 1.2, 14: 2.6}
"""Every panel is flat until the input moves (about 3.3, 2.3, 1.25 and 2.65 s on Figures 11-14):
the trim is the median before."""
TIME_TOL_S = {11: 0.0298, 12: 0.0203, 13: 0.0204, 14: 0.0353}
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
def collective_scales() -> tuple[float, float]:
    """The collective's travel, % of maximum per inch: the lowest and highest of three
    readings, since the report prints the travel in neither unit (#81). Each of Figures 11
    and 13 read as proportional, his trim inches against the collective our HMU needs there,
    and the straight line through both trims."""
    x = {f: collective(f) for f in COLLECTIVE_FIGS}
    c = {f: his_trim(f, "collective") for f in COLLECTIVE_FIGS}
    a, b = COLLECTIVE_FIGS
    readings = [x[f] / c[f] for f in COLLECTIVE_FIGS] + [(x[b] - x[a]) / (c[b] - c[a])]
    return min(readings), max(readings)


def scale_of(fig: int) -> float | None:
    """The collective scaling a figure is compared at, % of maximum per inch: proportional, his
    trim inches against the collective our HMU needs there. `XCPC` is "percent of maximum"
    [pdf p.81] and his collective axis starts at 0 in., so 0 in. reads as 0 % (#81). `None`
    where the collective does not move."""
    return collective(fig) / his_trim(fig, "collective") if fig in COLLECTIVE_FIGS else None


@functools.cache
def replay(fig: int, np_held: bool = False, scale: float | None = None) -> dict[str, np.ndarray]:
    """Our closed loop driven by his rotor speed: `t`, NG %, fuel and torque.

    Settled first with the load Appendix B's linear models carry about the trim (#6), then
    NP is set to his rotor speed every frame, and on Figures 11 and 13 the collective follows
    his, as a change from his trim times `scale`, % per inch. `np_held` keeps NP at the trim
    instead: what the loop does without the rotor-speed input, the decoy."""
    r = our_trim(fig)
    np_trim = r.state.np_rpm
    q_trim = r.frame.q_pt_ftlbf
    x0 = collective(fig)
    pilot = loop.Pilot(xcpc_pct=x0, pas_deg=100.0, pcprf_pct=his_trim(fig, "nr"))
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
    if scale is not None:
        tc, cc = his(fig, "collective")
        c0 = his_trim(fig, "collective")
    t = np.arange(tn[0], tn[-1], DT)
    out = np.zeros((t.size, 3))
    for i, ti in enumerate(t):
        np_rpm = np_trim if np_held else float(np.interp(ti, tn, nr)) / 100.0 * engine_c.NP_DES
        s = dataclasses.replace(s, engine=dataclasses.replace(s.engine, np_rpm=np_rpm))
        if scale is not None:
            xcpc = x0 + scale * (float(np.interp(ti, tc, cc)) - c0)
            pilot = dataclasses.replace(pilot, xcpc_pct=xcpc)
        s, _e, h, fr = loop.step(s, pilot, AMB, dt=DT, load=0.0, j_load=engine_c.J_LOAD_UH60A)
        out[i] = (100.0 * s.engine.ng_rpm / engine_c.NG_DES, h.wf_pph, fr.q_pt_ftlbf)
    return {"t": t, "ng": out[:, 0], "fuel": out[:, 1], "torque": out[:, 2]}


OWN = "own"
"""`scale` default: the figure's own `scale_of`."""


def ours(fig: int, key: str, np_held: bool = False, scale=OWN) -> tuple[np.ndarray, np.ndarray]:
    """`(t, value)`: NG in %, fuel and torque in % of our trim."""
    r = replay(fig, np_held, scale_of(fig) if scale == OWN else scale)
    v = r[key]
    return r["t"], v if key == "ng" else 100.0 * v / v[0]


def gap(fig: int, key: str, line: str = "gen_hel", np_held: bool = False, scale=OWN):
    """Per digitized column of his line, how far ours is from it: NG in % of design speed,
    fuel and torque in % of trim."""
    tb, vb = his(fig, key, line)
    if key != "ng":
        vb = 100.0 * vb / his_trim(fig, key, line)
    t, v = ours(fig, key, np_held, scale)
    keep = (tb >= t[0]) & (tb <= t[-1])
    return _best(t, v, tb[keep], vb[keep], np.zeros(keep.sum()), TIME_TOL_S[fig])


def held(fig: int, key: str) -> float:
    """The rms of `gap` with his inputs held at trim, NP and, where it moves, the collective:
    the decoy."""
    return rms(gap(fig, key, np_held=True, scale=0.0 if fig in COLLECTIVE_FIGS else None))


# ----------------------------------------------------------------------------- the checks


GEN_HEL_RMS = {
    11: {"ng": 0.23, "fuel": 2.4, "torque": 2.8},
    12: {"ng": 0.04, "fuel": 0.5, "torque": 0.47},
    13: {"ng": 0.49, "fuel": 2.5, "torque": 4.2},
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
"""The replay must be what does it: with his inputs held at the trim (`held`), ours sits at
least this many times further from his line."""


@pytest.mark.parametrize("fig", FIGS)
@pytest.mark.parametrize("key", PANELS)
def test_his_inputs_drive_the_response(fig: int, key: str):
    assert held(fig, key) > DECOY_MIN_RATIO * rms(gap(fig, key))


@pytest.mark.parametrize("fig", FIGS)
@pytest.mark.parametrize("key", PANELS)
def test_ours_is_nearer_gen_hel_than_the_flight_test(fig: int, key: str):
    """We replicate the simulation, not the aircraft: ours sits nearer Ballin's Gen Hel line
    than his flight-test line, which the test aircraft's greater load demand compensation
    [pdf p.47] and the collective trim mismatches [pdf pp.47, 50] move away from it."""
    assert rms(gap(fig, key)) < rms(gap(fig, key, "test"))


# ------------------------------------------- the second test: his rotor's demand in, speed out
#
# Replaying his rotor speed makes it an input, so the governor and the rotor's inertia are never
# tested. The second test runs the other way: the drive train is rigid [pdf p.47], so the torque
# the rotor demands is the engines' torque less what accelerates the inertia, and both are on
# his plots. That demand is the load, and our model predicts the rotor speed itself (#81).

RIPPLE_S = 0.12
"""His rotor speed is differentiated, and his torque smoothed, over this window: two periods of
the 17 Hz ripple in his rotor speed, so the rebuilt load carries his rotor's slow demand and not
the ripple."""
SLOPE_TRIM = {11: 2, 12: 2, 13: 1, 14: 1}
"""The Appendix B trim whose load slope (`DQ_REQ_DNP`) each figure takes: hover for the hover
cases, level flight at 80 kt for 55 and 90 kt [pdf p.67]."""
LOAD_PANELS = ("nr", "ng", "fuel")
"""What the second test predicts: his torque is its input."""
J_TOTAL = engine_c.J_PT + engine_c.J_LOAD_UH60A
"""The power turbine and its share of the rotor, as Appendix B's linear models carry it (#6)."""


def _smooth(t_src: np.ndarray, v_src: np.ndarray, t: np.ndarray, deriv: int = 0) -> np.ndarray:
    """`v_src` on the grid `t`, averaged over `RIPPLE_S` (held flat past the ends, where the
    figure's trim continues), or the time derivative of that average."""
    from scipy.ndimage import uniform_filter1d

    n = int(round(RIPPLE_S / DT)) | 1
    v = uniform_filter1d(np.interp(t, t_src, v_src), n, mode="nearest")
    return np.gradient(v, DT) if deriv else v


@functools.cache
def demand(fig: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """His rotor's torque demand per engine, ft*lbf on our engine's scale, and his NP, rpm, on
    our frame grid: `(t, q_req, np_his)`.

    `q_req = Q_eng - J_TOTAL * d(omega)/dt`, his total torque scaled to ours at the trim, since
    the gearbox ratio is not printed."""
    tn, nr = his(fig, "nr")
    tq, q = his(fig, "torque")
    t = np.arange(max(tn[0], tq[0]), min(tn[-1], tq[-1]), DT)
    q_eng = our_trim(fig).frame.q_pt_ftlbf * _smooth(tq, q, t) / his_trim(fig, "torque")
    np_his = _smooth(tn, nr, t) / 100.0 * engine_c.NP_DES
    domega = _smooth(tn, nr, t, deriv=1) / 100.0 * engine_c.NP_DES * RPM_TO_RAD_PER_SEC
    return t, q_eng - J_TOTAL * domega, np_his


@functools.cache
def load_replay(fig: int, held: bool = False) -> dict[str, np.ndarray]:
    """Our closed loop driven by his rotor's demand: `t`, NR %, NG %, fuel and torque.

    Settled at the trim as `replay` is, then the load is his demand at his speed plus Appendix
    B's slope times the difference between our speed and his, and on Figures 11 and 13 the
    collective follows his at `scale_of`. NP is our own. `held` keeps the load at the trim:
    what the loop does without his demand, the decoy."""
    r = our_trim(fig)
    np_trim = r.state.np_rpm
    q_trim = r.frame.q_pt_ftlbf
    slope = DQ_REQ_DNP[SLOPE_TRIM[fig]]
    x0, scale = collective(fig), scale_of(fig)
    pilot = loop.Pilot(xcpc_pct=x0, pas_deg=100.0, pcprf_pct=his_trim(fig, "nr"))
    s = loop.seed(r, pilot, AMB)
    at_trim = lambda v: q_trim + slope * (v - np_trim)  # noqa: E731
    for _ in range(SETTLE_FRAMES):
        s, *_ = loop.step(s, pilot, AMB, dt=DT, load=at_trim, j_load=engine_c.J_LOAD_UH60A)
    t, q_req, np_his = demand(fig)
    if scale is not None:
        tc, cc = his(fig, "collective")
        c0 = his_trim(fig, "collective")
    out = np.zeros((t.size, 4))
    for i, ti in enumerate(t):
        if held:
            load = at_trim
        else:
            load = lambda v, q=q_req[i], n=np_his[i]: q + slope * (v - n)  # noqa: E731
        if scale is not None:
            xcpc = x0 + scale * (float(np.interp(ti, tc, cc)) - c0)
            pilot = dataclasses.replace(pilot, xcpc_pct=xcpc)
        s, _e, h, fr = loop.step(s, pilot, AMB, dt=DT, load=load, j_load=engine_c.J_LOAD_UH60A)
        out[i] = (
            100.0 * s.engine.np_rpm / engine_c.NP_DES,
            100.0 * s.engine.ng_rpm / engine_c.NG_DES,
            h.wf_pph,
            fr.q_pt_ftlbf,
        )
    return {"t": t, "nr": out[:, 0], "ng": out[:, 1], "fuel": out[:, 2], "torque": out[:, 3]}


def ours_load(fig: int, key: str, held: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """`(t, value)` of the second test: NR and NG in %, fuel and torque in % of our trim."""
    r = load_replay(fig, held)
    v = r[key]
    return r["t"], v if key in ("nr", "ng") else 100.0 * v / v[0]


def gap_load(fig: int, key: str, line: str = "gen_hel", held: bool = False) -> np.ndarray:
    """Per digitized column of his line, how far the second test is from it: NR and NG in %,
    fuel and torque in % of trim."""
    tb, vb = his(fig, key, line)
    if key in ("fuel", "torque"):
        vb = 100.0 * vb / his_trim(fig, key, line)
    t, v = ours_load(fig, key, held)
    keep = (tb >= t[0]) & (tb <= t[-1])
    return _best(t, v, tb[keep], vb[keep], np.zeros(keep.sum()), TIME_TOL_S[fig])


LOAD_REPLAY_RMS = {
    11: {"nr": 0.2, "ng": 0.24, "fuel": 3.6},
    12: {"nr": 0.037, "ng": 0.046, "fuel": 0.93},
    13: {"nr": 0.2, "ng": 0.22, "fuel": 2.5},
    14: {"nr": 0.092, "ng": 0.084, "fuel": 2.8},
}
"""The second test against his Gen Hel line, rms: NR and NG in %, fuel in % of trim. A
characterization, 1.25 times today's value."""


@pytest.mark.parametrize("fig", FIGS)
@pytest.mark.parametrize("key", LOAD_PANELS)
def test_his_load_gives_his_rotor_speed(fig: int, key: str):
    got = rms(gap_load(fig, key))
    assert got < LOAD_REPLAY_RMS[fig][key], f"Figure {fig} {key}: rms {got:.3f}"


@pytest.mark.parametrize("fig", FIGS)
@pytest.mark.parametrize("key", LOAD_PANELS)
def test_his_load_drives_the_second_test(fig: int, key: str):
    """With the load held at the trim, the second test sits at least `DECOY_MIN_RATIO` times
    further from his line."""
    assert rms(gap_load(fig, key, held=True)) > DECOY_MIN_RATIO * rms(gap_load(fig, key))
