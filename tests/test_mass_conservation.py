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


# ------------------------------------------------------------------- energy, not just mass


@pytest.mark.parametrize("wf_pph", TRIMS_PPH)
def test_eq_39_is_an_exact_two_stream_work_balance(wf_pph: float):
    """Eq. 39's 0.71/0.29 split is not a fudge factor -- it places the bleed port.

    The report calls it "an empirically determined function" for the compressor
    interstage bleed [pdf p.25] and prints `K_QC1 = 0.71`, `K_QC2 = 0.29`. They sum to
    exactly 1, and that is the tell. Writing W24 for the station 2.4 bleed and h24 for
    its enthalpy, the physical two-stream work is

        WA3*(h3 - h2) + W24*(h24 - h2) = WA3*h3 + W24*h24 - WA2*h2

    since WA3 + W24 = WA2, while Eq. 39 with K_QC1 + K_QC2 = 1 is

        WA2*K_QC1*h3 + WA3*K_QC2*h3 - WA2*h2 = WA3*h3 + W24*K_QC1*h3 - WA2*h2

    The two are equal **if and only if h24 = K_QC1*h3 = 0.71 h3**. So the equation is an
    exact energy accounting that puts the station 2.4 port where the enthalpy is 71 % of
    compressor discharge -- which works out at 36-49 % of the enthalpy RISE above inlet,
    increasing with power. Verified here to 5e-13.
    """
    wf = wf_pps_from_pph(wf_pph)
    sweep = trim.sweep(list(range(130, int(wf_pph) + 1, 5)) + [wf_pph])
    r = [q for q in sweep if q.trustworthy][-1]
    f = frame(r.state, wf, AMB, q_req_ftlbf=r.q_req_ftlbf)

    w24 = f.wa2_pps * (f.b1 + f.b2)
    h24 = c.K_QC_1 * f.h3
    eq39 = f.wa2_pps * (c.K_QC_1 * f.h3 - f.h2) + f.wa3_pps * c.K_QC_2 * f.h3
    two_stream = f.wa3_pps * (f.h3 - f.h2) + w24 * (h24 - f.h2)

    assert c.K_QC_1 + c.K_QC_2 == pytest.approx(1.0, abs=1e-12), (
        "the decomposition only holds if the two coefficients sum to one"
    )
    assert abs(eq39 - two_stream) / eq39 < 1e-11, (
        f"Eq. 39 gives {eq39:.6f} Btu/s and the two-stream balance {two_stream:.6f}; "
        f"they must agree, and they do only because h24 = 0.71 h3"
    )


@pytest.mark.parametrize("wf_pph", TRIMS_PPH)
def test_the_combustor_energy_balance_closes(wf_pph: float):
    """Eq. 21 is an energy balance and must close: WA31*h3 + eta_b*Wf*HVF = W41*h41ns.

    Eq. 21 is written as `h41ns = (h3 + eta_b*FAR*HVF)/(1 + FAR)`. Multiplying through by
    `WA31*(1 + FAR) = WA31 + Wf = W41` recovers the balance, so this checks that our FAR
    and the station 4.1 mass balance are mutually consistent, not just the algebra.
    """
    wf = wf_pps_from_pph(wf_pph)
    sweep = trim.sweep(list(range(130, int(wf_pph) + 1, 5)) + [wf_pph])
    r = [q for q in sweep if q.trustworthy][-1]
    f = frame(r.state, wf, AMB, q_req_ftlbf=r.q_req_ftlbf)

    into = f.wa31_pps * f.h3 + f.eta_b * wf * c.HVF
    out = f.w41_pps * f.h41_ns
    assert abs(into - out) / into < 1e-11, f"combustor energy: in {into:.6f}, out {out:.6f} Btu/s"


@pytest.mark.parametrize("wf_pph", TRIMS_PPH)
def test_station_4_5_mixing_loses_energy_as_the_report_specifies(wf_pph: float):
    """Eq. 29 does NOT conserve energy, by construction, and this records how much.

    The report is explicit: *"At station 4.5, gases from the station 4.4 and cooling-bleed
    flow from the compressor are mixed before passing through the power turbine. The
    enthalpy of the mixed gases is proportional to enthalpy [at] station 4.4"*
    [pdf p.24] -- so Eq. 29 is `h45 = K_H45 * h44`, a single multiplicative fraction and
    not a flow-weighted mix of the two streams. A flow-weighted mix would give

        h45 = (W41*h44 + B3*K_bl*WA2*h3) / W45

    and the coefficient that would produce it is **0.974 to 0.979** across the five trims,
    strikingly constant, against the printed **0.9623**. So the printed value sits 1.1 to
    1.7 % below an energy-conserving mix, and that deficit is the entire reason the
    engine's overall energy balance does not close.

    Reproduced as printed, per CLAUDE.md, and asserted rather than corrected: the band
    below fails if the deficit ever changes, whether by our error or by a re-reading of
    `K_H45`. It affects only the power turbine -- h45 feeds T45, theta45, dH_PT and W45 --
    so it does not touch the gas generator dynamics of open question #47.
    """
    wf = wf_pps_from_pph(wf_pph)
    sweep = trim.sweep(list(range(130, int(wf_pph) + 1, 5)) + [wf_pph])
    r = [q for q in sweep if q.trustworthy][-1]
    f = frame(r.state, wf, AMB, q_req_ftlbf=r.q_req_ftlbf)

    returned = f.b3 * c.K_BL * f.wa2_pps
    h45_mixed = (f.w41_pps * f.h44 + returned * f.h3) / f.w45_pps
    deficit = (f.h45 - h45_mixed) / h45_mixed
    equivalent = h45_mixed / f.h44

    assert -0.020 < deficit < -0.010, (
        f"station 4.5 energy deficit is {100 * deficit:+.2f} % at {wf_pph} lbm/hr; "
        f"on record is -1.1 to -1.7 %. The energy-conserving coefficient here is "
        f"{equivalent:.4f} against the printed K_H45 = {c.K_H45}."
    )
    assert 0.970 < equivalent < 0.982, (
        f"the flow-weighted equivalent of K_H45 is {equivalent:.4f}; on record is "
        f"0.974-0.979 across the trims"
    )
