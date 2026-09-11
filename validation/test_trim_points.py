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

from t700 import trim
from t700.engine import Ambient
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
