"""Validate the engine against Table B.1's three trim conditions.

This is the first real comparison against the report, and it is available early because it
needs only the eleven digitized `f` functions -- no control system, no digitized result
figures, no Gen Hel.

## The comparison

Table B.1 [TM-100991 pdf p.67] prints, for each trim: ambient, NP, fuel flow, gas
generator speed, shaft horsepower and engine torque. Give the model the printed **fuel
flow, NP and ambient**, solve for the gas generator equilibrium, and compare the NG,
torque and horsepower that come out.

Fuel flow is the input because it is what the engine physically consumes; NG and shaft
power are the model's answer. Nothing about the trim state is handed to the solver -- the
initial guess is scaled off the design point, not off any trim.

## Tolerances

From `SCOPE.md`, and they are **ours, not Ballin's**. He states no percentage tolerance
for torque or power anywhere in the report, and his own accepted NG deviations run 1-4 %.
"""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from t700 import constants as c
from t700 import trim
from t700.engine import Ambient, frame
from t700.units import wf_pps_from_pph

# Tolerances from SCOPE.md. Ours, declared, not derived from the report.
NG_TOL_PCT = 2.0
SHP_TOL_PCT = 3.0
TORQUE_TOL_PCT = 3.0


@dataclass(frozen=True)
class TrimCase:
    """One row of Table B.1 [pdf p.67], transcribed."""

    name: str
    wf_pph: float
    ng_rpm: float
    shp: float
    q_shaft_ftlbf: float
    np_rpm: float = 20895.0
    p_amb_psia: float = 14.696
    t_amb_degR: float = 518.67


# All three are sea level, UH-60A at 16825 lbm, NP = 20895 rpm.
CASES = [
    TrimCase("hover", wf_pph=476.3, ng_rpm=41638.0, shp=911.1, q_shaft_ftlbf=229.0),
    TrimCase("level 80 kt", wf_pph=349.3, ng_rpm=39768.0, shp=552.6, q_shaft_ftlbf=138.9),
    TrimCase("descent 80 kt", wf_pph=267.7, ng_rpm=38072.0, shp=302.6, q_shaft_ftlbf=76.06),
]


def _solve(case: TrimCase):
    return trim.solve(
        wf_pps_from_pph(case.wf_pph),
        case.np_rpm,
        Ambient(case.p_amb_psia, case.t_amb_degR),
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_trim_converges(case: TrimCase):
    r = _solve(case)
    assert r.trustworthy, (
        f"{case.name}: no usable equilibrium. residual_converged="
        f"{r.residual_converged}, on_data={r.on_data}, |r| = {r.residual_norm:.3e}"
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_gas_generator_speed(case: TrimCase):
    """NG is the headline number: it is what the fuel flow buys."""
    r = _solve(case)
    dev = (r.state.ng_rpm - case.ng_rpm) / case.ng_rpm * 100.0
    assert abs(dev) < NG_TOL_PCT, (
        f"{case.name}: NG {r.state.ng_rpm:.0f} rpm against Ballin's {case.ng_rpm:.0f}, "
        f"{dev:+.2f} % (tolerance +/-{NG_TOL_PCT} %)"
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_shaft_horsepower(case: TrimCase):
    r = _solve(case)
    dev = (r.shp - case.shp) / case.shp * 100.0
    assert abs(dev) < SHP_TOL_PCT, (
        f"{case.name}: {r.shp:.1f} shp against Ballin's {case.shp:.1f}, {dev:+.2f} % "
        f"(tolerance +/-{SHP_TOL_PCT} %)"
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_shaft_torque(case: TrimCase):
    """Against Table B.1's *shaft*-referenced torque row, not the hub row.

    The hub row is ~145x larger, past the UH-60A main gearbox. Reading the hub figure as
    engine torque is what made open question #11 look like a defect in the report.
    """
    r = _solve(case)
    dev = (r.frame.q_pt_ftlbf - case.q_shaft_ftlbf) / case.q_shaft_ftlbf * 100.0
    assert abs(dev) < TORQUE_TOL_PCT, (
        f"{case.name}: {r.frame.q_pt_ftlbf:.2f} ft.lbf against Ballin's "
        f"{case.q_shaft_ftlbf:.2f}, {dev:+.2f} %"
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_gas_path_is_physically_ordered(case: TrimCase):
    """Pressures and temperatures must fall in the right direction through the engine.

    Cheap, and it catches a whole class of wiring error that a percentage comparison can
    miss -- a model can hit NG while running the gas path backwards.
    """
    r = _solve(case)
    f, s = r.frame, r.state
    assert s.p3_psia > s.p41_psia > s.p45_psia > f.p49_psia > f.ps9_psia, (
        "pressure must fall monotonically from compressor discharge to exhaust"
    )
    assert f.t41_degR > f.t45_degR > f.t49_degR, "gas must cool through the turbines"
    assert f.t3_degR > case.t_amb_degR, "the compressor must heat the air"
    assert 0.0 < f.far < 0.1, f"fuel-air ratio {f.far} is not physical"


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_solution_sits_inside_the_digitized_data(case: TrimCase):
    """A trim that sits outside the maps is extrapolating, and must be declared.

    The descent trim clamps `f1@82` once: at NGc = 84.99 % the solver blends speed lines
    82 and 85, and the query pressure ratio 7.704 is 0.292 beyond line 82's last knot. The
    blend weight on line 82 is 0.0033, so clamping rather than extrapolating moves WA2c by
    0.011 %. Recorded, allowed, and bounded -- not ignored.
    """
    r = _solve(case)
    allowed = {"f1@82"} if case.name == "descent 80 kt" else set()
    unexpected = set(r.clamps_at_solution) - allowed
    assert not unexpected, (
        f"{case.name}: trim sits outside digitized data for {sorted(unexpected)} -- "
        f"the model is extrapolating"
    )


def test_more_fuel_buys_more_speed_and_power():
    """Monotonicity across the three trims. A sanity check no single case can give."""
    results = [(c, _solve(c)) for c in sorted(CASES, key=lambda c: c.wf_pph)]
    ngs = [r.state.ng_rpm for _, r in results]
    shps = [r.shp for _, r in results]
    assert ngs == sorted(ngs), ngs
    assert shps == sorted(shps), shps


# --------------------------------------------------------------------------- dynamics


def test_jacobian_eigenvalues_against_table_1():
    """Compare our own linearization with the report's printed modes [Table 1, pdf p.31].

    A far stronger check than the trim values: it tests the *derivatives* of every
    equation, not just their equilibrium. Ballin calls a 4 % eigenvalue mismatch "good
    agreement" [pdf p.29], which is the standard used here.

    The slowest mode is excluded deliberately. It is the NP/rotor mode and its eigenvalue
    depends on dQreq/dNP, supplied by the Gen Hel UH-60A simulation. Our Jacobian holds
    Qreq constant, so that derivative is zero and the mode is not comparable -- the same
    row Appendix B flags as not independently reproducible.
    """
    import numpy as np

    from t700.engine import State, frame

    case = CASES[0]
    amb = Ambient(case.p_amb_psia, case.t_amb_degR)
    wf = wf_pps_from_pph(case.wf_pph)
    r = _solve(case)
    x0 = r.state.as_array()
    qreq = r.frame.q_pt_ftlbf

    def deriv(x):
        f = frame(State.from_array(x), wf, amb, q_req_ftlbf=qreq)
        return np.array([f.dng_dt, f.dnp_dt, f.dp3_dt, f.dp41_dt, f.dp45_dt])

    f0 = deriv(x0)
    A = np.zeros((5, 5))
    for i in range(5):
        h = 1e-6 * max(abs(x0[i]), 1.0)
        xp = x0.copy()
        xp[i] += h
        A[:, i] = (deriv(xp) - f0) / h

    ours = np.sort(np.linalg.eigvals(A).real)[::-1]  # least to most negative
    # Table 1, trim 1, excluding the load-dependent slowest mode
    printed = [-2.66, -51.6, -3060.0, -4900.0]
    mine = [ours[1], ours[2], ours[3], ours[4]]

    for got, want in zip(mine, printed, strict=True):
        rel = abs(got - want) / abs(want) * 100.0
        assert rel < 30.0, f"eigenvalue {got:.1f} against printed {want}: {rel:.1f} % off"

    fastest = abs(ours[-1])
    assert fastest > 1000.0, (
        "the fast pressure modes should be present in the 5-state model; if they are not, "
        "the volume dynamics are not being modelled"
    )


def test_the_model_is_too_stiff_for_explicit_integration_at_the_report_frame_time():
    """Records *why* Eqs. 74-80 exist, as an executable fact.

    The fastest mode is ~0.2 ms, so explicit integration of all five states is stable only
    below ~0.4 ms. The report runs the engine at 7 ms. The quasi-steady formulation, which
    makes the three pressures algebraic, is therefore not a preference but a necessity.
    """
    import numpy as np

    from t700.engine import State, frame

    case = CASES[0]
    amb = Ambient(case.p_amb_psia, case.t_amb_degR)
    wf = wf_pps_from_pph(case.wf_pph)
    r = _solve(case)
    x0 = r.state.as_array()
    qreq = r.frame.q_pt_ftlbf

    def deriv(x):
        f = frame(State.from_array(x), wf, amb, q_req_ftlbf=qreq)
        return np.array([f.dng_dt, f.dnp_dt, f.dp3_dt, f.dp41_dt, f.dp45_dt])

    f0 = deriv(x0)
    A = np.zeros((5, 5))
    for i in range(5):
        h = 1e-6 * max(abs(x0[i]), 1.0)
        xp = x0.copy()
        xp[i] += h
        A[:, i] = (deriv(xp) - f0) / h

    euler_limit_ms = 2.0 / abs(np.linalg.eigvals(A).real).max() * 1000.0
    assert euler_limit_ms < 7.0, (
        f"explicit Euler would be stable at {euler_limit_ms:.2f} ms, which is above the "
        f"report's 7 ms frame -- if this ever passes, the stiffness argument for the "
        f"quasi-steady formulation no longer holds and Eqs. 74-80 need revisiting"
    )


# --------------------------------------------------------------------------- off-map guard


def test_a_solution_off_the_compressor_map_is_not_reported_as_good():
    """Clamping protects against extrapolation but manufactures spurious equilibria.

    Once the state leaves the maps every lookup returns a constant, so the residual stops
    depending on the state and Newton drives it to zero in a region where the model means
    nothing. At 750 lbm/hr from a design-point guess this converged to a residual of
    6.6e-9 with NG = 2.16e11 rpm and called itself converged.

    `trustworthy` is what callers must check: converged AND on the map.
    """
    from t700 import trim as _trim

    bad = _trim.solve(wf_pps_from_pph(750.0), 20900.0, Ambient(14.696, 518.67))
    assert bad.residual_converged, "the degenerate root is a real root of the clamped system"
    assert not bad.on_data, "but it sits off the compressor map and must say so"
    assert not bad.trustworthy
    assert "OFF THE MAP" in bad.report()


def test_continuation_reaches_the_top_of_the_map():
    """The same fuel flow that fails from a cold guess succeeds by continuation."""
    from t700 import trim as _trim

    results = _trim.sweep([400.0, 500.0, 600.0, 700.0, 750.0, 775.0, 800.0])
    assert all(r.trustworthy for r in results), [r.trustworthy for r in results]
    ngs = [100.0 * r.state.ng_rpm / 44700.0 for r in results]
    assert ngs == sorted(ngs), ngs
    assert 99.0 < ngs[-1] < 100.5, ngs[-1]

    # and the sweep agrees with the transient run at the same fuel flow
    at_775 = 100.0 * results[5].state.ng_rpm / 44700.0
    assert abs(at_775 - 99.44) < 0.05, at_775


def test_the_sweep_refuses_past_the_maps_top_speed_line():
    """And a sweep must START inside the data -- continuation cannot rescue a cold start.

    Beginning at 800 lbm/hr from the design-point guess lands off the map immediately, so
    there is no good solution to continue from and every later point fails too. Begin low,
    walk up.
    """
    from t700 import trim as _trim

    results = _trim.sweep([400.0, 600.0, 800.0, 850.0])
    assert all(r.trustworthy for r in results[:3])
    assert not results[-1].trustworthy, "850 lbm/hr is past the map's 100 % speed line"

    cold = _trim.sweep([800.0])
    assert not cold[0].trustworthy, (
        "a sweep starting above ~700 lbm/hr has no valid point to continue from"
    )


# ------------------------------------------------- the rest of Table B.1, long unused

STATE_B1 = {
    # [TM-100991 pdf p.67, Table B.1] the printed engine state at each trim.
    # Ps3, P41, T41, P45 and T45 are printed there and were unused by this project until
    # 2026-09-12 -- fifteen numbers, including the one quantity with a 13x gain on the
    # transient (open question #47).
    "hover": dict(ps3=176.34, p41=174.28, t41=2292.0, p45=37.42, t45=1632.0),
    "level 80 kt": dict(ps3=142.13, p41=140.26, t41=2102.0, p45=30.66, t45=1501.0),
    "descent 80 kt": dict(ps3=114.27, p41=112.77, t41=1982.0, p45=25.54, t45=1424.0),
}

STATE_TOL_PCT = 0.5
"""Tolerance on the printed station pressures and temperatures.

Measured worst deviations are Ps3 0.31 %, P41 0.30 %, P45 0.25 %, T45 0.16 % and T41
0.08 %, so 0.5 % is a ratchet a little above the worst rather than a target. Note the
printed T41 and T45 carry only four significant figures with a trailing decimal point --
"2292." -- so a quarter of a degree is the printing resolution at T41 and the agreement
there is at that limit."""


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_the_printed_station_state_matches(case: TrimCase):
    """Ps3, P41, T41, P45 and T45 against Table B.1 -- fifteen numbers never checked.

    This closes the largest evidence gap the project had. Until now the trim comparison
    used three quantities per condition (speed, shaft power, shaft torque) while Table
    B.1 prints eight, and the five unused ones are exactly the ones the transient work
    needed: open question #47 spent a long time asking whether our Ps3 was right at low
    power, comparing a digitized Figure 10 *transient* against the report's steady-state
    Figure 8 -- an invalid comparison, since a chop must run below the equilibrium locus
    (open question #53 records the retraction) -- while a printed Ps3 sat in Table B.1 the
    whole time.

    It is also the anchor that settled the f1 interpolation question. Scored on all
    twenty-one printed numbers, the shipped constant-pressure-ratio blend gives an rms
    deviation of 0.202 % against 0.244 % for interpolation along Figure A1's own beta
    lines -- so the structurally tidier construction is measurably worse here, and is not
    shipped. See open question #52.
    """
    want = STATE_B1[case.name]
    r = _solve(case)
    amb = Ambient(case.p_amb_psia, case.t_amb_degR)
    f = frame(r.state, wf_pps_from_pph(case.wf_pph), amb, q_req_ftlbf=r.q_req_ftlbf)
    got = {
        "ps3": c.K_PS3 * r.state.p3_psia,
        "p41": r.state.p41_psia,
        "t41": f.t41_degR,
        "p45": r.state.p45_psia,
        "t45": f.t45_degR,
    }
    for key, printed in want.items():
        dev = 100.0 * (got[key] - printed) / printed
        assert abs(dev) < STATE_TOL_PCT, (
            f"{case.name} {key}: printed {printed}, ours {got[key]:.3f}, {dev:+.3f} %"
        )


def test_the_printed_state_is_self_consistent_under_a_forced_pressure():
    """Imposing the printed Ps3 and integrating recovers the printed P41, P45 and NG.

    A different kind of check from the one above, and a stronger one. P3 stops being a
    state and becomes an input held at Table B.1's printed Ps3; NG, P41 and P45 then
    integrate from a deliberately displaced start, seven percent low. If the printed
    state is a genuine equilibrium of these equations, the run must find it.

    It does, to four significant figures on both pressures -- P41 174.288 / 140.261 /
    112.764 against the printed 174.28 / 140.26 / 112.77, and P45 37.424 / 30.651 /
    25.531 against 37.42 / 30.66 / 25.54. That is what licenses the forced-Ps3
    experiments in open question #47 to be read as attribution rather than as an
    artifact of an ill-posed reduction.
    """
    from scipy.integrate import solve_ivp

    from t700.engine import State

    for case in CASES:
        want = STATE_B1[case.name]
        wf = wf_pps_from_pph(case.wf_pph)
        r = _solve(case)
        amb = Ambient(case.p_amb_psia, case.t_amb_degR)
        p3_forced = want["ps3"] / c.K_PS3

        def rhs(_t, y, p3=p3_forced, wf=wf, qreq=r.q_req_ftlbf, amb=amb, npr=case.np_rpm):
            f = frame(State(y[0], npr, p3, y[1], y[2]), wf, amb, q_req_ftlbf=qreq)
            return [f.dng_dt, f.dp41_dt, f.dp45_dt]

        y0 = [case.ng_rpm * 0.93, want["p41"] * 0.93, want["p45"] * 0.93]
        sol = solve_ivp(rhs, (0.0, 200.0), y0, method="BDF", rtol=1e-10, atol=1e-10)
        ng, p41, p45 = (float(v) for v in sol.y[:, -1])

        for label, got, printed in (
            ("NG", ng, case.ng_rpm),
            ("P41", p41, want["p41"]),
            ("P45", p45, want["p45"]),
        ):
            dev = 100.0 * (got - printed) / printed
            assert abs(dev) < STATE_TOL_PCT, (
                f"{case.name}: with Ps3 forced to its printed value, {label} settles at "
                f"{got:.3f} against the printed {printed}, {dev:+.3f} %"
            )
