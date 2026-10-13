"""Validate the engine against Table B.1's three trim conditions [TM-100991 pdf p.67].

Give the model the printed fuel flow, NP and ambient, solve for the equilibrium, and compare
what comes out: NG, shaft power, torque and the printed station state. Nothing about the trim
state is handed to the solver; the initial guess is scaled off the design point.

Tolerances are ours, declared in `SCOPE.md`: the report states none for torque or power, and
Ballin's own accepted NG deviations run 1-4 %.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from t700 import constants as c
from t700 import trim
from t700.engine import Ambient, frame
from t700.linear import DOF, extract
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
    """Against Table B.1's *shaft*-referenced torque row, not the hub row, which is past the
    UH-60A main gearbox and ~145x larger (open question #11).
    """
    r = _solve(case)
    dev = (r.frame.q_pt_ftlbf - case.q_shaft_ftlbf) / case.q_shaft_ftlbf * 100.0
    assert abs(dev) < TORQUE_TOL_PCT, (
        f"{case.name}: {r.frame.q_pt_ftlbf:.2f} ft.lbf against Ballin's "
        f"{case.q_shaft_ftlbf:.2f}, {dev:+.2f} %"
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_gas_path_is_physically_ordered(case: TrimCase):
    """Pressures and temperatures must fall in the right direction through the engine: a
    model can hit NG while running the gas path backwards.
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

    The descent trim sits just above the 85 % speed line and records no clamp. It is still
    allowed `f1@82` -- line 82 read past its last knot, which a dip below 85 % would cause --
    so a small move does not fail the test; any other clamp is extrapolation.
    """
    r = _solve(case)
    allowed = {"f1@82"} if case.name == "descent 80 kt" else set()
    unexpected = set(r.clamps_at_solution) - allowed
    assert not unexpected, (
        f"{case.name}: trim sits outside digitized data for {sorted(unexpected)} -- "
        f"the model is extrapolating"
    )


# --------------------------------------------------------------------------- dynamics


EIGENVALUE_TOL_PCT = 6.0
"""Ours, declared in `SCOPE.md`. A characterization, not an accuracy claim: a regression
detector set a quarter above today's measurement at the report's +-2 % extraction step
[pdf p.28]. Ballin calls his own Table 1 columns, about 4 % apart, "in good agreement"
[pdf p.29]; ours is not that."""


def test_jacobian_eigenvalues_against_table_1():
    """Our linearization against the report's printed 5-DOF modes [Table 1, pdf p.31], at
    hover: a check on the derivatives of every equation, not just their equilibrium.

    Paired by physics, never by position in a sorted spectrum: `comparable_eigenvalues`
    drops NP structurally (`A[1,1]`), since sorting mispairs the slow modes (see
    `tests/test_linear.py::test_sorting_would_mispair_the_slow_modes`).
    """
    case = CASES[0]
    amb = Ambient(case.p_amb_psia, case.t_amb_degR)
    wf = wf_pps_from_pph(case.wf_pph)
    r = _solve(case)

    model = extract(r, wf, DOF.FIVE, amb)
    ours = np.sort(np.real(model.comparable_eigenvalues))[::-1]  # least to most negative

    # Table 1, trim 1, NP row excluded: it is built from dQreq/dNP, a Gen Hel quantity.
    # Imported here: `test_all_dof_models` imports this file, so a top-level import is circular.
    from test_all_dof_models import TABLE_1_5DOF

    printed = [TABLE_1_5DOF[1][k] for k in ("NG", "P3", "P45", "P41")]
    # at the report's +-2 % step [pdf p.28]
    on_record = [1.2, 0.2, 3.5, 2.2]

    assert ours.size == len(printed), f"expected 4 comparable modes, got {ours}"
    for got, want, expected in zip(ours, printed, on_record, strict=True):
        # magnitude convention, as `test_all_dof_models` uses: negative means our mode is
        # less damped than the printed one
        dev = (abs(got) - abs(want)) / abs(want) * 100.0
        assert abs(dev) < EIGENVALUE_TOL_PCT, (
            f"eigenvalue {got:.1f} against printed {want}: {dev:+.1f} %"
        )
        assert abs(dev - expected) < 1.0, (
            f"eigenvalue {got:.1f} against printed {want} measures {dev:+.1f} %, "
            f"{expected:+.1f} % on record. Paired by physics, not by position."
        )


def test_the_model_is_too_stiff_for_explicit_integration_at_the_report_frame_time():
    """Why Eqs. 74-80 exist: explicit integration of all five states is stable only far below
    the report's 7 ms frame, so making the three pressures algebraic is a necessity.
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
    """The degenerate root of the clamped system, and why a cold solve cannot reach it.

    Off the function tables every lookup is constant, so both powers go flat in NG while both
    torques fall as 1/NG: `dNG/dt` tends to zero as NG tends to infinity, a root of the
    residual that is not an equilibrium. `solve`'s line search refuses a step that leaves
    `f1`'s speed range widened by a tenth, so a cold solve trims on the map or fails honestly.
    A power-balance residual would remove the root too but under-weights the NG equation at
    low speed. Open question #71.
    """
    import numpy as np

    from t700 import trim as _trim
    from t700.engine import State, frame

    amb = Ambient(14.696, 518.67)
    wf = wf_pps_from_pph(160.0)

    # 1. The mechanism. Hold the pressures, sweep NG far past the map, and watch the two
    #    forms of the same imbalance diverge in behaviour.
    base = _trim.solve(wf, 20900.0, amb)
    st = base.state
    powers, torques = [], []
    for ng in (1.0e5, 1.0e7, 1.0e9):
        f = frame(State(ng, st.np_rpm, st.p3_psia, st.p41_psia, st.p45_psia), wf, amb)
        dq = f.q_gt_ftlbf - f.q_c_ftlbf
        torques.append(abs(dq))
        powers.append(abs(dq * ng * 2.0 * np.pi / 60.0 / 550.0))
    assert max(powers) - min(powers) < 1e-6 * max(powers), (
        f"the power imbalance should not depend on NG once the tables clamp: {powers}"
    )
    assert torques[-1] < 1e-3 * torques[0], (
        f"the torque imbalance must vanish as 1/NG -- that is the spurious root: {torques}"
    )
    assert max(powers) > 100.0, (
        f"and the engine is nowhere near equilibrium there: {max(powers):.1f} hp out"
    )

    # 2. And the search can no longer get there. Every cold solve across the range either
    #    trims on the map or refuses; none converges off it.
    for pph in (110.0, 115.0, 125.0, 160.0, 400.0, 750.0, 800.0, 850.0):
        r = _trim.solve(wf_pps_from_pph(pph), 20900.0, amb)
        assert r.trustworthy or not r.residual_converged, (
            f"{pph:.0f} lbm/hr converged off the map at NGc {r.frame.ngc_pct:.4g} -- the "
            f"degenerate root is back"
        )

    # 3. The report path still exists for a caller that hands in a guess already off the
    #    map, which is the one way to get there and is the escape hatch the bound leaves.
    far = State(4.0e6, 20900.0, st.p3_psia, st.p41_psia, st.p45_psia)
    off = _trim.solve(wf, 20900.0, amb, guess=far)
    assert not off.on_data, "a guess started off the map should stay off it"
    assert not off.trustworthy
    assert "OFF THE MAP" in off.report()


def test_continuation_reaches_the_top_of_the_map():
    """The same fuel flow that fails from a cold guess succeeds by continuation."""
    from t700 import trim as _trim

    results = _trim.sweep([400.0, 500.0, 600.0, 700.0, 750.0, 775.0, 800.0])
    assert all(r.trustworthy for r in results), [r.trustworthy for r in results]
    ngs = [100.0 * r.state.ng_rpm / c.NG_DES for r in results]
    assert ngs == sorted(ngs), ngs
    assert 99.0 < ngs[-1] < 100.5, ngs[-1]

    # and the sweep agrees with the transient run at the same fuel flow
    at_775 = 100.0 * results[5].state.ng_rpm / c.NG_DES
    assert abs(at_775 - 99.44) < 0.05, at_775


def test_the_sweep_refuses_past_the_maps_top_speed_line():
    """850 lbm/hr is past `f1`'s 100 % speed line and no amount of continuation reaches it:
    the solve must pin at the bound and report `residual_converged` False.

    Trap: at the top speed line `np.clip` pins the parameter and a forward-difference
    Jacobian reads a zero column; `solve` differences inward at that bound and projects trial
    steps onto it, so a cold 800 lbm/hr trims with no continuation.
    """
    from t700 import trim as _trim

    results = _trim.sweep([400.0, 600.0, 700.0, 800.0, 850.0])
    assert all(r.trustworthy for r in results[:4])
    assert not results[-1].trustworthy, "850 lbm/hr is past the map's 100 % speed line"
    assert not results[-1].residual_converged, (
        "850 must fail honestly, not converge somewhere meaningless"
    )

    # the top of the map no longer needs continuation to reach
    cold = _trim.solve(wf_pps_from_pph(800.0), 20900.0, Ambient(14.696, 518.67))
    assert cold.trustworthy, (
        f"a cold 800 lbm/hr should now trim on the map, got NGc {cold.frame.ngc_pct:.4g}"
    )
    assert abs(cold.frame.ngc_pct - 99.8904) < 0.01, cold.frame.ngc_pct

    # and past the map a cold start fails the same way a continued one does
    for pph in (850.0, 900.0):
        far = _trim.solve(wf_pps_from_pph(pph), 20900.0, Ambient(14.696, 518.67))
        assert not far.trustworthy and not far.residual_converged, (
            f"{pph:.0f} lbm/hr is past the map and must refuse, got "
            f"NGc {far.frame.ngc_pct:.4g}, converged={far.residual_converged}"
        )


# ------------------------------------------------- the rest of Table B.1, the station state

STATE_B1 = {
    # [TM-100991 pdf p.67, Table B.1] the printed engine state at each trim.
    "hover": dict(ps3=176.34, p41=174.28, t41=2292.0, p45=37.42, t45=1632.0),
    "level 80 kt": dict(ps3=142.13, p41=140.26, t41=2102.0, p45=30.66, t45=1501.0),
    "descent 80 kt": dict(ps3=114.27, p41=112.77, t41=1982.0, p45=25.54, t45=1424.0),
}

STATE_TOL_PCT = 0.5
"""Tolerance on the printed station pressures and temperatures: a ratchet a little above the
worst, declared in `SCOPE.md`. The printed T41 and T45 carry four significant figures
("2292."), so the T41 agreement is at the printing resolution."""


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_the_printed_station_state_matches(case: TrimCase):
    """Ps3, P41, T41, P45 and T45 against Table B.1: printed numbers, with no read error of
    ours in them.

    Interpolating `f1` on the beta grid is slightly worse on this table than constant
    pressure ratio and better on everything derivative; the trade is deliberate (open
    questions #43 and #60).
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

    P3 becomes an input held at Table B.1's Ps3; NG, P41 and P45 integrate from a start seven
    percent low. If the printed state is an equilibrium of these equations, the run finds it,
    which is what lets the forced-Ps3 experiments of open question #47 be read as attribution.
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
