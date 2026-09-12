"""Mass conservation across the engine, at steady state and through a transient.

## Why this file exists

Nothing asserted it. The three volume equations *are* the mass balances, so it is easy to
assume they close -- but the streams in them are each determined **independently**, and
that is what makes closure a real check rather than an identity:

| stream | set by | equations |
|---|---|---|
| WA2 | the compressor map `f1`, uncorrected | 7, 9 |
| WA31 | the combustor pressure-drop orifice | 18 |
| W41 | the choked turbine nozzle | 28 |
| W45 | the power turbine map `f9` | 33, 34 |

Four independent determinations, three balances. If any map or constant were wrong in a
way that broke continuity, these residuals would show it, and no other test in the
repository would.

## The overboard accounting, which is not obvious from the equations

    in   = WA2 + Wf
    out  = W45                        through the power turbine
    lost = WA2*(B1 + B2)              station 2.4 bleed, Eq. 15 -- seal air and the
                                      power turbine balance piston; neither returns
         + WA2*K_B3                   the part of Eq. 16 that never returns
         + WA2*B3*(1 - K_bl)          the diffuser bleed not reintroduced at 4.5

`K_bl = 0.7826` is the fraction of the diffuser bleed that *is* reintroduced, and it
appears in Eq. 44 as `+ B3 K_bl WA2`. So the diffuser bleed is split and only part of it
leaves; the station 2.4 bleed leaves entirely.
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import constants as c
from t700 import realtime, trim
from t700.engine import Ambient, State, frame
from t700.units import wf_pps_from_pph

AMB = Ambient(14.696, 518.67)

TRIMS_PPH = [476.3, 349.3, 267.7, 200.0, 175.0]
"""The three Table B.1 trims plus two below them. The low ones matter most here: the
seal bleed B1 is zero above NGc 89 and 0.1089 below 78, so only the last two exercise it
at all -- see `test_the_overboard_bleed_is_only_exercised_below_the_printed_trims`."""

CLOSURE_TOL = 1e-9
"""Residuals as a fraction of WA2. These are algebraic identities at a converged trim, so
the only thing between them and zero is floating point; the measured worst is ~7e-13."""


def _streams(wf_pph: float, real_time: bool = False) -> dict[str, float]:
    wf = wf_pps_from_pph(wf_pph)
    sweep = trim.sweep(list(range(130, int(wf_pph) + 1, 5)) + [wf_pph])
    trims = [r for r in sweep if r.trustworthy]
    assert trims, f"no trustworthy trim at {wf_pph} lbm/hr"
    r = trims[-1]
    f = frame(r.state, wf, AMB, q_req_ftlbf=r.q_req_ftlbf)

    if real_time:
        st = realtime.from_trim(r, f.wa31_pps, f)
        tr = realtime.run(
            st,
            lambda t: wf,
            AMB,
            duration_s=3.0,
            dt=0.007,
            q_req_ftlbf=f.q_pt_ftlbf,
            heat_sink=True,
        )
        s = State(tr["ng"][-1], tr["np"][-1], tr["p3"][-1], tr["p41"][-1], tr["p45"][-1])
        f = frame(s, wf, AMB, q_req_ftlbf=f.q_pt_ftlbf, t41_degR=tr["t41"][-1])

    lost = f.wa2_pps * (f.b1 + f.b2 + c.K_B3 + f.b3 * (1.0 - c.K_BL))
    return {
        "wa2": f.wa2_pps,
        "eq42": f.wa3_pps - f.wa2_pps * (f.b3 + c.K_B3) - f.wa31_pps,
        "eq43": f.wa31_pps + wf - f.w41_pps,
        "eq44": f.w41_pps - f.w45_pps + f.b3 * c.K_BL * f.wa2_pps,
        "global": (f.wa2_pps + wf) - f.w45_pps - lost,
        "lost_frac": lost / f.wa2_pps,
        "ngc_pct": f.ngc_pct,
        "b1": f.b1,
    }


@pytest.mark.parametrize("wf_pph", TRIMS_PPH)
@pytest.mark.parametrize("real_time", [False, True], ids=["differential", "real-time"])
def test_every_volume_balance_closes_at_steady_state(wf_pph: float, real_time: bool):
    s = _streams(wf_pph, real_time)
    for key, label in (
        ("eq42", "station 3, Eq. 42"),
        ("eq43", "station 4.1, Eq. 43"),
        ("eq44", "station 4.5, Eq. 44"),
    ):
        rel = abs(s[key]) / s["wa2"]
        assert rel < CLOSURE_TOL, (
            f"{label} does not close at {wf_pph} lbm/hr: residual {s[key]:.3e} lbm/s, "
            f"{rel:.3e} of WA2"
        )


@pytest.mark.parametrize("wf_pph", TRIMS_PPH)
@pytest.mark.parametrize("real_time", [False, True], ids=["differential", "real-time"])
def test_the_engine_conserves_mass_overall(wf_pph: float, real_time: bool):
    """In minus out minus overboard, which no individual volume equation states."""
    s = _streams(wf_pph, real_time)
    rel = abs(s["global"]) / s["wa2"]
    assert rel < CLOSURE_TOL, (
        f"global balance fails at {wf_pph} lbm/hr: (WA2 + Wf) - W45 - lost = "
        f"{s['global']:.3e} lbm/s, {rel:.3e} of WA2"
    )


def test_the_overboard_bleed_is_only_exercised_below_the_printed_trims():
    """A gap in the evidence, recorded because it bears on open question #47.

    The overboard fraction is driven by `B1`, the station 2.4 seal bleed, which Figure A3
    makes zero above NGc 89 and 0.1089 below 78. Measured here, the total overboard loss
    runs about 3 % at the two high trims and reaches 14 % at 175 lbm/hr.

    The consequence is that **`B1`'s plateau is unconstrained by any printed trim**: of
    the three Table B.1 conditions, two sit where B1 is exactly zero and the third at
    NGc 85.2 sees only 0.047, less than half the plateau. Everything below NGc 78, where
    the plateau is fully in effect, has no printed anchor at all -- and that is the region
    where Figure 10's decay disagrees with us.

    This test does not claim B1 is wrong. Figure A3 is unambiguous and its thirteen
    glyphs are reproduced to 1e-4. It records that the quantity with the most leverage on
    the one remaining transient disagreement is also the one the printed data constrains
    least, so that the next person does not have to rediscover it.
    """
    seen = {}
    for wf_pph in (476.3, 349.3, 267.7, 175.0):
        s = _streams(wf_pph)
        seen[wf_pph] = (s["ngc_pct"], s["b1"], s["lost_frac"])

    assert seen[476.3][1] == 0.0, "B1 should be zero at the hover trim"
    assert seen[349.3][1] == 0.0, "B1 should be zero at the level trim"
    assert 0.0 < seen[267.7][1] < 0.5 * 0.1089, (
        f"the descent trim should see part of the B1 ramp, not the plateau; got "
        f"{seen[267.7][1]:.5f} at NGc {seen[267.7][0]:.2f} %"
    )
    assert seen[175.0][2] > 0.12, (
        f"at 175 lbm/hr the overboard loss should exceed 12 % of WA2; got "
        f"{100 * seen[175.0][2]:.2f} %"
    )
    assert seen[476.3][2] < 0.04, (
        f"at hover it should be under 4 %; got {100 * seen[476.3][2]:.2f} %"
    )


def test_the_balances_open_during_a_transient_and_close_again():
    """The counterpart to the steady-state tests, and it is not symmetric with them.

    At equilibrium the balances close to machine precision. During a transient they must
    NOT -- if they did, nothing would be storing mass and the volume dynamics would be
    doing no work. Measured through a 400 -> 600 lbm/hr step, as a fraction of WA2:

    | t after step | differential | real-time |
    |---|---|---|
    | 7 ms | 7.5e-2 | 6.3e-2 |
    | 0.1 s | 9.8e-3 | 1.3e-2 |
    | 1.5 s | 8.0e-5 | 1.7e-3 |
    | 2.5 s | 3.3e-6 | 6.1e-4 |

    Both open by about 7 % and both decay, but for different reasons, and that is the
    quasi-steady approximation made visible. In the differential model the residual IS
    the storage: `dP3/dt` and the rest are literally these terms. In the real-time frame
    the pressures are algebraic, so nothing stores -- yet the balances still open,
    because Eqs. 76, 78 and 80 enforce the orifice, the combustor pressure drop and the
    station 4.5 continuity rather than Eqs. 42 and 43 themselves, and WA31 is a frame
    behind by Eq. 74. The real-time residual therefore decays more slowly -- 6e-4 against
    3e-6 at 2.5 s -- which is the lag working itself out rather than a volume emptying.
    """
    wf0 = wf_pps_from_pph(400.0)
    r0 = trim.solve(wf0, c.NP_DES, AMB)
    f0 = frame(r0.state, wf0, AMB)
    st = realtime.from_trim(r0, f0.wa31_pps, f0)
    tr = realtime.run(
        st,
        lambda t: wf0 if t < 0.5 else wf_pps_from_pph(600.0),
        AMB,
        duration_s=3.0,
        dt=0.007,
        q_req_ftlbf=f0.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=True,
    )
    t = np.asarray(tr["t"])

    def residuals(i: int) -> tuple[float, float, float]:
        wf = wf0 if t[i] < 0.5 else wf_pps_from_pph(600.0)
        s = State(tr["ng"][i], tr["np"][i], tr["p3"][i], tr["p41"][i], tr["p45"][i])
        f = frame(s, wf, AMB, q_req_ftlbf=f0.q_pt_ftlbf, t41_degR=tr["t41"][i])
        return (
            abs(f.wa3_pps - f.wa2_pps * (f.b3 + c.K_B3) - f.wa31_pps) / f.wa2_pps,
            abs(f.wa31_pps + wf - f.w41_pps) / f.wa2_pps,
            abs(f.w41_pps - f.w45_pps + f.b3 * c.K_BL * f.wa2_pps) / f.wa2_pps,
        )

    before = residuals(int(np.argmin(np.abs(t - 0.40))))
    assert max(before) < CLOSURE_TOL, f"the balances should be closed before the step; got {before}"

    just_after = residuals(int(np.argmax(t >= 0.5)))
    assert max(just_after) > 1e-2, (
        f"the balances should open by more than 1 % of WA2 at the step -- if they do "
        f"not, nothing is storing mass and the volume terms are inert; got {just_after}"
    )

    settled = residuals(int(np.argmin(np.abs(t - 2.9))))
    assert max(settled) < 1e-3, f"and they should close again 2.4 s after the step; got {settled}"
    assert max(settled) < 0.05 * max(just_after), (
        "the residual should fall by more than an order of magnitude between the step "
        "and 2.4 s later"
    )
