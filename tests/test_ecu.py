"""Structural checks on the ECU -- Appendix C, Figures C1-C8.

As with the HMU, these test the wiring rather than agreement with a published result:
Appendix C has no equations, so an error here is a misread picture, and a misread picture
produces plausible numbers.

The ECU is fiddlier than the HMU in exactly two places, and both get their own test: the
maximum error selector (which path is governing), and the nonlinear NP loop gain `B4`
(which is not the constant it looks like).
"""

from __future__ import annotations

import pytest

from t700 import constants as engine_c
from t700.control import constants as c
from t700.control import ecu

DT = 0.007
BASE = dict(
    np_rpm=20895.0,
    torq45_ftlbf=229.0,
    t45_degR=1632.0,
    w45_pps=7.90,
    p45_psia=37.42,
    pcprf_pct=100.0,
)


def settle(n: int = 1000, **over):
    u = ecu.ECUInputs(**{**BASE, **over})
    s = ecu.seed(u)
    o = None
    for _ in range(n):
        s, o = ecu.step(s, u, DT)
    return s, o


def test_sensed_speed_is_percent_of_design():
    """[Fig. C1] PCNP is NP as a percent of NP_DES, and the whole speed path is in those
    units -- so a reference of 100 % means NP_DES rpm, not Table B.1's 20895."""
    _, o = settle(np_rpm=engine_c.NP_DES, pcprf_pct=100.0)
    assert o.spder == pytest.approx(0.0, abs=1e-9)


def test_speed_error_is_sensed_minus_reference_not_the_other_way_round():
    """The sign that inverts the whole governor if read backwards. [Fig. C1, read at
    600 dpi] SPDER is *positive when the turbine is running fast*, and it drives SPDG up,
    and the HMU subtracts the resulting trim -- so overspeed takes fuel away."""
    _, fast = settle(np_rpm=engine_c.NP_DES * 1.02)
    _, slow = settle(np_rpm=engine_c.NP_DES * 0.98)
    assert fast.spder > 0.0 > slow.spder
    assert fast.spdg > slow.spdg, "an overspeed must command more trim, not less"


def test_the_load_share_path_is_inert_with_one_engine():
    """[Fig. C2] The switch is drawn open, and clamp(0, DBIAS, CE) - DBIAS = 0 makes the
    rest of the path zero too. Inert by construction, not by our choice."""
    _, o = settle(engine2_torq45_ftlbf=None)
    assert o.xqlo == pytest.approx(0.0, abs=1e-12)


def test_a_second_engine_moves_the_load_share_trim():
    """And with the switch closed it must not be inert, or the path is mis-wired."""
    _, same = settle(engine2_torq45_ftlbf=229.0, n=3000)
    _, more = settle(engine2_torq45_ftlbf=600.0, n=3000)
    assert same.xqlo == pytest.approx(0.0, abs=1e-9)
    assert more.xqlo != pytest.approx(0.0, abs=1e-6)


def test_the_temperature_path_is_inert_at_every_printed_trim():
    """[Table B.1] T45 is 1632 / 1501 / 1424 deg R against T45REF = 2004, so ET45 is 350
    to 570 deg below the limit. A limiter that engaged in steady flight would be wrong."""
    for t45 in (1632.0, 1501.0, 1424.0):
        _, o = settle(t45_degR=t45)
        assert o.limiting == "speed", (t45, o.tsig, o.spdsf)
        assert o.et45 < -300.0, (t45, o.et45)


def test_an_overtemperature_takes_the_selector():
    """[Fig. C1] The maximum error selector is what makes T4.5 limiting override speed
    governing, and it only does so when its error signal is the larger."""
    _, o = settle(t45_degR=2400.0, n=3000)
    assert o.et45 > 0.0
    assert o.limiting == "t45", (o.tsig, o.spdsf)
    assert o.spdss == pytest.approx(o.tsig, abs=1e-12)


def test_the_thermocouple_time_constant_is_a_lookup_not_a_constant():
    """[Fig. C5] TAU45 = F_EC1(W45R, T45L), and W45R needs W45 in lbm/sec -- the
    nomenclature's lb/hr puts it 2000x off the printed axis."""
    _, o = settle()
    lo, hi = 1.0, 15.0
    assert lo < o.w45r < hi, f"W45R {o.w45r} is outside F_EC1's printed domain"
    _, hot = settle(t45_degR=2000.0)
    assert hot.tau45_s != pytest.approx(o.tau45_s, rel=1e-6), "TAU45 is not varying"


def test_the_nonlinear_loop_gain_switches_on_load():
    """[Fig. C3] B4 is not a constant. Y integrates TRQL - CR and the gain switches when
    1000*B6*g(SPDER) + Y passes CORR -- so a loaded engine gets the high gain and an
    unloaded one does not. Open question #11 measured the crossing at 0.10-0.36 s."""
    _, loaded = settle(torq45_ftlbf=229.0, n=2000)
    assert loaded.b4 == 1.0, loaded.torq45_ftlbf if hasattr(loaded, "torq45_ftlbf") else loaded.b4
    _, idle = settle(torq45_ftlbf=5.0, n=2000)
    assert idle.b4 == 0.0, "an unloaded engine should sit on the low loop gain"


def test_a_large_speed_error_also_raises_the_loop_gain():
    """The relay path: g(SPDER) = 1 outside [-1, +4] %NP, times 1000, swamps CORR. So a
    slam transient gets the high gain immediately, without waiting for the torque
    integrator to wind up."""
    _, o = settle(torq45_ftlbf=5.0, np_rpm=engine_c.NP_DES * 0.9, n=50)
    assert o.spder < c.NP_LOOP_RELAY_LOW
    assert o.b4 == 1.0


def test_the_integrator_lower_limit_follows_the_other_engine():
    """[Fig. C8] ZLOLIM is -1.0 below 180 ft*lbf on engine 2 and -0.3 above -- so Table
    C.1's printed -1.0 is the single-engine value, not a universal constant."""
    from t700.control.ecu import _zlolim

    assert _zlolim(None) == c.ZLOLIM
    assert _zlolim(100.0) == c.ZLOLIM
    assert _zlolim(400.0) == c.ZLOLIM_HIGH_TORQUE


def test_the_integrator_respects_whichever_limit_applies():
    _, _o = settle(np_rpm=engine_c.NP_DES * 0.90, n=6000)
    s, _ = settle(np_rpm=engine_c.NP_DES * 0.90, n=6000)
    assert c.ZLOLIM - 1e-12 <= s.pi_int <= c.ZHILIM + 1e-12, s.pi_int


def test_a_held_input_does_not_make_the_ecu_diverge():
    """The P+I integrator ramps while a speed error persists -- that is its job -- but it
    must stop at its limit rather than run away."""
    s1, _ = settle(n=4000, np_rpm=engine_c.NP_DES * 0.95)
    s2, _ = settle(n=20000, np_rpm=engine_c.NP_DES * 0.95)
    assert abs(s2.pi_int) <= max(abs(c.ZLOLIM), abs(c.ZHILIM)) + 1e-12
    assert s2.pi_int == pytest.approx(c.ZHILIM, abs=1e-9) or s2.pi_int == pytest.approx(
        s1.pi_int, rel=1e-3
    )


# --------------------------------------- the temperature limiter, which nothing measured

OVERTEMP_DEGR = 2050.0
"""A T4.5 the engine model cannot reach, fed to the ECU directly.

`T45REF` is 2004 degR and `T45COR` adds 11 to the raw measurement, so the Figure C1
selector hands over above **1993 degR** of true T4.5. Table B.1's hottest trim is 1632
[pdf p.67], and reaching 1993 needs an NG well past the 100 % ceiling Ballin declares the
model valid to [pdf p.39]. So the limiter is out of reach of every *closed-loop* run in
this project, and `validation/test_closed_loop.py` records that as the reason `TLGE` moved
freely there.

It is not out of reach of the ECU, which is a block diagram with inputs. 2050 degR is
inside `F_EC1`'s own printed parameter range, whose curves run to 2060 degR [Fig. C23,
pdf p.94] -- the figure is drawn over the limiter's operating range because that is the
only place the signal matters."""

OVERTEMP_SPDG = (1.9438, 1.9478)
"""`SPDG` 8.4 s after a step to `OVERTEMP_DEGR`, from the hover trim's sensed state.

A characterization -- a tripwire on Appendix C's temperature path, not an accuracy claim;
the report prints no overtemperature transient. It closes the last two holes open question
#65 left in the ECU. Measured 1.945780, and the band is +/-0.1 %:

| | SPDG | vs band |
|---|---|---|
| as shipped | 1.945780 | -- |
| `F_EC1` x1.1 | 1.720032 | -11.6 %, far outside |
| `TLGE` x1.1 | 1.942242 | -0.18 %, 1.8x the band |
| `TLGE` x0.9 | 1.949306 | +0.18 % |

Nothing in this harness touches the engine, the thermodynamics or a digitized engine map:
every input is a literal. So the only things that can move this number are Appendix C's
own constants and `F_EC1`, which is what lets the band be this tight."""


def _overtemp(t45_hot=OVERTEMP_DEGR, n=1500, at=300):
    """Hold the hover trim's sensed state, then step T4.5 into the limiter."""
    base = dict(
        np_rpm=20895.0,
        torq45_ftlbf=229.0,
        w45_pps=8.18,
        p45_psia=37.42,
        pcprf_pct=20895.0 * 100.0 / engine_c.NP_DES,
    )
    s = ecu.seed(ecu.ECUInputs(t45_degR=1632.0, **base), spdg=0.49496)
    hot_frames, o = 0, None
    for i in range(n):
        s, o = ecu.step(s, ecu.ECUInputs(t45_degR=(t45_hot if i >= at else 1632.0), **base), DT)
        hot_frames += o.limiting == "t45"
    return o, hot_frames


def test_the_selector_hands_over_at_the_printed_reference_minus_the_correction():
    """[Figs. C1, C5] `ET45 = T45EL - T45REF` and `T45E = T45 + T45COR`, so the handover
    is at `T45REF - T45COR` of true temperature and at nothing else.

    This is a wiring test and **it does not pin either constant**: it computes the
    threshold from the same two rows it brackets, so moving both together moves the test
    with them. What pins them is `OVERTEMP_SPDG` below, where the step lands at a fixed
    2050 degR and the distance above the reference sets the whole response."""
    threshold = c.T45REF - c.T45COR
    assert _overtemp(threshold - 15.0)[1] == 0, "the limiter took over below the reference"
    assert _overtemp(threshold + 20.0)[1] > 0, "the limiter did not take over above it"


def test_the_thermocouple_lag_shapes_the_overtemperature_handover():
    """`F_EC1` and `TLGE` were observable by nothing until this test.

    Both sit on the T4.5 path: the harness lag `TLGE`, then a second lag whose own time
    constant is the `F_EC1` lookup. Downstream of them is the Figure C6 compensation and
    the Figure C1 maximum selector, so with the limiter governing they set how fast the
    ECU winds the torque motor up -- and with the limiter *not* governing they set nothing
    at all, which is why every other test in this project is blind to them.
    """
    o, hot_frames = _overtemp()
    assert o.limiting == "t45", (o.tsig, o.spdsf)
    assert hot_frames > 0
    assert OVERTEMP_SPDG[0] < o.spdg < OVERTEMP_SPDG[1], (
        f"SPDG {o.spdg:.6f} after an overtemperature step is outside "
        f"{OVERTEMP_SPDG}; the T4.5 path has moved"
    )


def test_the_sensor_time_constant_is_looked_up_inside_its_printed_domain():
    """[Fig. C23, pdf p.94] `F_EC1` is drawn over W45R 1-15 and T45L 1260-2060 degR. Both
    arguments must land inside that box at the limiter's own operating point, or the model
    is extrapolating exactly where the schedule was published to be used."""
    o, _ = _overtemp()
    assert 1.0 < o.w45r < 15.0, f"W45R {o.w45r} is outside F_EC1's printed domain"
    assert 1260.0 < o.t45el < 2060.0, f"T45L {o.t45el} is outside F_EC1's printed domain"
    assert o.tau45_s > 0.0
