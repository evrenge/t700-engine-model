"""Linearizing the real-time map itself -- open question #49.

## The oddity this resolves

Table 1 [pdf p.31] prints two eigenvalue columns that ought to describe the same system:
column 4, a 2-DOF model, and column 5, an order-reduced 5-DOF. Ballin's differ -- -2.69
against -2.81 at hover, -2.23 against -2.16 at level, -1.82 against -1.83 at descent.
**Ours are identical**, to 2e-8, which is the floor a central-difference Jacobian reaches.

That identity is not a bug, and it is not luck: it is a theorem. If `g(s) = f1(s, p(s))`
with `f2(s, p) = 0` defining `p`, then implicit differentiation gives
`dg/ds = A11 - A12 A22^-1 A21` -- exactly the Schur complement that order reduction forms.
Linearizing an *exactly solved* quasi-steady system and residualizing the full Jacobian are
the same operation. Any model that solves its pressures exactly must produce two identical
columns here.

So the question is not why ours agree. It is **why Ballin's do not**, and the answer is on
pdf p.27: his 2-DOF comes from a *separately coded* nonlinear program, not the algebraic
limit. It carries the real-time numerics -- Eq. 74's opened compressor iteration, a finite
inner-loop count, explicit integration at a finite frame. His 0.12 /sec hover spread is a
measurement of what those cost.

## Putting a number on it

This file linearizes `realtime.step` as a **discrete map** -- six states, including Eq. 74's
carried mass flow, which the continuous model does not have -- and converts it with
`logm(A_d)/dt`. At hover, with the pressures converged:

    continuous 2-DOF / reduced-5        -2.444     -9.1 % against Ballin's -2.69
    discrete map at the report's 7 ms   -2.728     +1.4 %

So most of the hover gap is the frame, not the physics, and **Table 1 column 4 becomes an
independent check** rather than something we cannot address. That is the report's central
engineering claim -- that a quasi-steady real-time model is faithful enough -- quantified
from its own printed numbers.

Two further findings, both measured below:

* **The carried mass flow's mode is numerical, not physical.** It scales as 1/dt --
  -28 /sec at 14 ms, -75 at 7 ms, -3739 at 0.2 ms -- so it is an artifact of the one-frame
  opening and vanishes in the continuous limit, rather than a dynamic the engine has.
* **The descent trim keeps a ~22 % residual under every reading**, continuous or discrete,
  tight tolerance or printed. That one is not the real-time approximation; it is the
  interpolation ambiguity of open questions #31 and #43, which is worst at descent.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import logm

from t700 import constants as c
from t700 import realtime, trim
from t700.engine import Ambient, frame
from t700.linear import DOF, extract
from t700.units import wf_pps_from_pph

AMB = Ambient(14.696, 518.67)
NP_RPM = 20895.0
WF_PPH = {1: 476.3, 2: 349.3, 3: 267.7}
DQ_REQ_DNP = {1: 0.019471, 2: 0.015869, 3: 0.012950}
"""Recovered from Appendix B; see `test_all_dof_models`. Needed only for the NP mode."""

TABLE_1_2DOF = {1: -2.69, 2: -2.23, 3: -1.82}
"""Table 1 column 4 [pdf p.31] -- Ballin's separately coded quasi-steady program."""

FRAME_S = 0.007
"""The report's engine frame [pdf p.47]."""

DISCRETE_TOL_PCT = {1: 5.0, 2: 5.0, 3: 25.0}
"""Ours, declared in `SCOPE.md`. Descent is loose because it carries the #31/#43
interpolation residual, which is present in the continuous model too and is not what this
file is measuring."""


def discrete_spectrum(
    case: int, dt: float = FRAME_S, tol: float | None = None, lag_whole_flow: bool = True
) -> np.ndarray:
    """Eigenvalues of `logm(A_d)/dt`, ascending, where `A_d` linearizes one real-time frame.

    Six states: NG, NP, the three pressures, and Eq. 74's carried compressor flow. The
    pressures are states of the *map* even though they are algebraic in the continuous
    model -- the frame starts its iteration from the previous frame's values.
    """
    wf = wf_pps_from_pph(WF_PPH[case])
    r = trim.solve(wf, NP_RPM, AMB)
    f = frame(r.state, wf, AMB)
    st0 = realtime.from_trim(r, f.wa31_pps, f, lag_whole_flow=lag_whole_flow)
    q0, np0 = f.q_pt_ftlbf, r.state.np_rpm
    kw = dict(
        ambient=AMB,
        dt=dt,
        j_load=c.J_LOAD_UH60A,
        integrate_np=True,
        heat_sink=False,
        lag_whole_flow=lag_whole_flow,
    )
    if tol is not None:
        kw["tol"] = tol

    def advance(v):
        s = realtime.RTState(v[0], v[1], v[2], v[3], v[4], v[5], st0.hs_tm_degR)
        s2, _ = realtime.step(s, wf, q_req_ftlbf=q0 + DQ_REQ_DNP[case] * (v[1] - np0), **kw)
        return np.array(
            [s2.ng_rpm, s2.np_rpm, s2.p3_psia, s2.p41_psia, s2.p45_psia, s2.wa31_carry_pps]
        )

    x0 = np.array(
        [
            st0.ng_rpm,
            st0.np_rpm,
            st0.p3_psia,
            st0.p41_psia,
            st0.p45_psia,
            st0.wa31_carry_pps,
        ]
    )
    A = np.zeros((6, 6))
    for j in range(6):
        h = max(abs(x0[j]) * 2e-6, 1e-9)
        xp, xm = x0.copy(), x0.copy()
        xp[j] += h
        xm[j] -= h
        A[:, j] = (advance(xp) - advance(xm)) / (2.0 * h)
    return np.sort(np.real(np.linalg.eigvals(np.real(logm(A)) / dt)))


def continuous_ng_mode(case: int) -> float:
    wf = wf_pps_from_pph(WF_PPH[case])
    r = trim.solve(wf, NP_RPM, AMB)
    m = extract(r, wf, dof=DOF.TWO, j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[case])
    return float(np.sort(np.real(m.comparable_eigenvalues))[-1])


def test_the_two_continuous_models_are_identical_because_they_must_be():
    """The Schur complement of an exactly solved quasi-steady system *is* the reduction.

    This is the premise of the whole file: our columns 4 and 5 cannot differ, so any
    difference in Ballin's is about his numerics and not about the engine.
    """
    for case in (1, 2, 3):
        wf = wf_pps_from_pph(WF_PPH[case])
        r = trim.solve(wf, NP_RPM, AMB)
        kw = dict(j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[case])
        two = np.sort(np.real(extract(r, wf, dof=DOF.TWO, **kw).comparable_eigenvalues))
        red = np.sort(np.real(extract(r, wf, dof=DOF.REDUCED_FIVE, **kw).comparable_eigenvalues))
        assert np.allclose(two, red, atol=1e-6), (case, two, red)


@pytest.mark.parametrize("case", (1, 2, 3))
def test_the_discrete_map_moves_toward_ballins_2dof_column(case: int):
    """Table 1 column 4 is Ballin's real-time program, so compare a real-time map to it."""
    disc = discrete_spectrum(case, tol=1e-12)[-2]
    dev = 100.0 * (disc / TABLE_1_2DOF[case] - 1.0)
    assert abs(dev) < DISCRETE_TOL_PCT[case], (
        f"trim {case}: the discrete map gives {disc:.3f} /sec against Table 1 column 4's "
        f"{TABLE_1_2DOF[case]}, {dev:+.1f} %"
    )


def test_the_frame_is_what_separates_the_discrete_map_from_the_continuous_one():
    """Refine dt and the discrete map must converge on the continuous eigenvalue.

    If it did not, the difference would be a structural error rather than the cost of the
    real-time approximation -- which is the distinction this whole file rests on.
    """
    cont = continuous_ng_mode(1)
    modes = [discrete_spectrum(1, dt=dt, tol=1e-12)[-2] for dt in (0.014, 0.007, 0.0035, 0.00175)]
    errs = [abs(m - cont) for m in modes]
    assert all(b < a for a, b in zip(errs, errs[1:], strict=False)), (cont, modes)
    ratios = [a / b for a, b in zip(errs, errs[1:], strict=False)]
    assert all(1.5 < q < 3.0 for q in ratios), (
        f"halving dt should roughly halve the error for an explicit first-order scheme, "
        f"got ratios {ratios}"
    )
    assert errs[-1] < 0.05 * abs(cont), (cont, modes)


def test_the_carried_mass_flows_mode_is_numerical_and_not_physical():
    """Eq. 74's one-frame memory adds a mode the engine does not have, and it is an
    artifact of the discretization rather than a dynamic.

    The signature is in the *discrete* eigenvalue, not the continuous one: `z` stays a
    fixed fraction per **frame** as the frame shrinks -- 0.68, 0.59, 0.54, 0.49 at
    dt = 14, 7, 3.5, 1.75 ms -- rather than a fixed fraction per second. A mode that
    decays in a fixed number of frames has no continuous limit: `lambda = ln(z)/dt` runs
    to -infinity. A physical mode would hold `lambda` fixed and let `z` approach 1.
    """
    dts = (0.014, 0.007, 0.0035, 0.00175)
    modes = [discrete_spectrum(1, dt=dt, tol=1e-12)[-3] for dt in dts]
    z = [np.exp(m * dt) for m, dt in zip(modes, dts, strict=True)]
    assert all(0.35 < q < 0.80 for q in z), (
        f"the carried-flow mode's per-frame decay is no longer frame-fixed: z = {z}"
    )
    assert abs(modes[-1]) > 6.0 * abs(modes[0]), (
        f"|lambda| should grow at least as fast as 1/dt: {modes}"
    )
    slow = discrete_spectrum(1, dt=dts[0], tol=1e-12)[-2]
    z_slow = np.exp(slow * dts[0])
    assert z_slow > 0.95, (
        f"for contrast, the physical NG mode decays by only {1 - z_slow:.4f} per frame, "
        f"so it is the one with a continuous limit"
    )


def test_the_printed_iteration_tolerance_costs_a_measurable_amount():
    """The report stops its pressure loops at 0.1 percent in ten and eight passes
    [pdf p.37]. That is not free, and this says what it costs on the slow mode."""
    for case in (1, 2, 3):
        tight = discrete_spectrum(case, tol=1e-12)[-2]
        printed = discrete_spectrum(case)[-2]
        shift = 100.0 * (printed / tight - 1.0)
        assert -30.0 < shift < 0.0, (
            f"trim {case}: the printed tolerance moves the NG mode {shift:+.1f} % "
            f"({tight:.3f} -> {printed:.3f} /sec)"
        )


def test_the_descent_residual_survives_every_reading_so_it_is_not_the_frame():
    """Trim 3 is ~22 % off Table 1 continuous *and* discrete, tight tolerance and printed.
    Whatever it is, it is not the real-time approximation -- it is the interpolation
    ambiguity of #31 and #43, which is worst at descent."""
    cont = continuous_ng_mode(3)
    disc = discrete_spectrum(3, tol=1e-12)[-2]
    for name, v in (("continuous", cont), ("discrete", disc)):
        dev = 100.0 * (v / TABLE_1_2DOF[3] - 1.0)
        assert -30.0 < dev < -15.0, f"descent {name}: {v:.3f} /sec, {dev:+.1f} %"
    assert abs(disc - cont) < 0.1 * abs(cont), (
        "the frame barely moves the descent mode, which is the point: its residual is structural"
    )
