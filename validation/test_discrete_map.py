"""Linearizing the real-time map itself -- open question #49.

Table 1 [pdf p.31] prints a 2-DOF column and an order-reduced 5-DOF column that ought to
describe the same system; Ballin's differ and ours are identical. That is a theorem: if
`g(s) = f1(s, p(s))` with `f2(s, p) = 0` defining `p`, implicit differentiation gives
`dg/ds = A11 - A12 A22^-1 A21`, the Schur complement order reduction forms. His differ
because his 2-DOF is a separately coded program [pdf p.27] carrying the real-time numerics.

This file linearizes `realtime.step` as a six-state discrete map (including Eq. 74's carried
mass flow) and takes `ln(z)/dt` per mode, so Table 1 column 4 becomes an independent check of
what the frame costs. Two of the six modes are frame artifacts with no continuous limit
(P45's one-frame lag and Eq. 74's carry). Modes are identified by participation factor, never
by position in a sorted spectrum. See SCOPE.md, "What the real-time approximation costs".
"""

from __future__ import annotations

import functools

import numpy as np
import pytest

from t700 import constants as c
from t700 import realtime, trim
from t700.engine import Ambient
from t700.linear import DOF, extract
from t700.units import wf_pps_from_pph

# DQ_REQ_DNP is needed only for the NP mode. TABLE_1_2DOF is Table 1 column 4 [pdf p.31] --
# Ballin's separately coded quasi-steady program.
from test_all_dof_models import DQ_REQ_DNP, TABLE_1_2DOF, WF_PPH

AMB = Ambient(14.696, 518.67)
NP_RPM = 20895.0

FRAME_S = 0.007
"""The report's engine frame [pdf p.47]."""

DISCRETE_TOL_PCT = {1: 8.0, 2: 15.0, 3: 7.0}
"""Ours, declared in `SCOPE.md`. The loose one (level) carries the #31/#43 interpolation
residual, which the continuous model has too and which this file is not measuring."""


TANGENT = 1e-5
"""The continuous models here are tangents: the limit the discrete map converges to, and
the Jacobian the Schur-complement identity is about. The shipped extraction uses the
report's +-2 % (pdf p.28), a secant, which is neither."""

STATE_NAMES = ("NG", "NP", "P3", "P41", "P45", "WA31")
"""Order of the six states of the frame map, for participation factors."""


def participation(A: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """`(continuous eigenvalues, P)` where `P[k, i]` is state k's participation in mode i.

    Decomposed from `A_d` itself with a per-mode `ln(z)/dt`, not via `logm(A_d)`, which fails
    outright on a negative real eigenvalue. Such modes are returned as `nan`: they have no
    real logarithm and no continuous counterpart.

    `P = |V * inv(V).T|`, column-normalised: the standard participation factor, invariant to
    the units of the states (rpm, psia and lbm/s here).
    """
    z, V = np.linalg.eig(A)
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(
            (np.abs(z.imag) < 1e-12) & (z.real <= 0.0),
            np.nan,
            np.log(z.astype(complex)) / dt,
        )
    P = np.abs(V * np.linalg.inv(V).T)
    return lam, P / P.sum(axis=0, keepdims=True)


def mode_of(A: np.ndarray, state: str, dt: float = FRAME_S) -> complex:
    """The eigenvalue that `state` participates in most. Identification by physics."""
    lam, P = participation(A, dt)
    return complex(lam[int(np.argmax(P[STATE_NAMES.index(state)]))])


@functools.cache
def _frame_jacobian(
    case: int, dt: float = FRAME_S, tol: float | None = None, lag_whole_flow: bool = True
) -> np.ndarray:
    """`A_d`, the 6x6 linearization of one real-time frame.

    Six states: NG, NP, the three pressures, and Eq. 74's carried compressor flow. The
    pressures are states of the *map* even though they are algebraic in the continuous
    model -- the frame starts its iteration from the previous frame's values.
    """
    wf = wf_pps_from_pph(WF_PPH[case])
    r = trim.solve(wf, NP_RPM, AMB)
    st0 = realtime.from_trim(r, lag_whole_flow=lag_whole_flow)
    q0, np0 = r.frame.q_pt_ftlbf, r.state.np_rpm
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
    return A


def discrete_ng_mode(case: int, dt: float = FRAME_S, tol: float | None = None) -> float:
    """The frame map's NG mode, identified by participation factor rather than by position
    (see `tests/test_linear.py::test_sorting_would_mispair_the_slow_modes`)."""
    return float(np.real(mode_of(_frame_jacobian(case, dt, tol), "NG", dt)))


def continuous_ng_mode(case: int) -> float:
    wf = wf_pps_from_pph(WF_PPH[case])
    r = trim.solve(wf, NP_RPM, AMB)
    m = extract(
        r, wf, dof=DOF.TWO, j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[case], rel_step=TANGENT
    )
    return float(np.sort(np.real(m.comparable_eigenvalues))[-1])


def test_the_two_continuous_models_are_identical_because_they_must_be():
    """The Schur complement of an exactly solved quasi-steady system *is* the reduction, so
    any difference between Ballin's columns 4 and 5 is his numerics, not the engine.
    """
    for case in (1, 2, 3):
        wf = wf_pps_from_pph(WF_PPH[case])
        r = trim.solve(wf, NP_RPM, AMB)
        kw = dict(j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[case], rel_step=TANGENT)
        two = np.sort(np.real(extract(r, wf, dof=DOF.TWO, **kw).comparable_eigenvalues))
        red = np.sort(np.real(extract(r, wf, dof=DOF.REDUCED_FIVE, **kw).comparable_eigenvalues))
        assert np.allclose(two, red, atol=1e-6), (case, two, red)


@pytest.mark.parametrize("case", (1, 2, 3))
def test_the_discrete_map_moves_toward_ballins_2dof_column(case: int):
    """Table 1 column 4 is Ballin's real-time program, so compare a real-time map to it."""
    disc = discrete_ng_mode(case, tol=1e-12)
    dev = 100.0 * (disc / TABLE_1_2DOF[case] - 1.0)
    assert abs(dev) < DISCRETE_TOL_PCT[case], (
        f"trim {case}: the discrete map gives {disc:.3f} /sec against Table 1 column 4's "
        f"{TABLE_1_2DOF[case]}, {dev:+.1f} %"
    )


def test_the_frame_is_what_separates_the_discrete_map_from_the_continuous_one():
    """Refine dt and the discrete map must converge on the continuous eigenvalue, at first
    order; otherwise the difference is a structural error, not the cost of the frame.
    """
    cont = continuous_ng_mode(1)
    modes = [discrete_ng_mode(1, dt=dt, tol=1e-12) for dt in (0.014, 0.007, 0.0035, 0.00175)]
    errs = [abs(m - cont) for m in modes]
    assert all(b < a for a, b in zip(errs, errs[1:], strict=False)), (cont, modes)
    ratios = [a / b for a, b in zip(errs, errs[1:], strict=False)]
    assert all(1.5 < q < 3.0 for q in ratios), (
        f"halving dt should roughly halve the error for an explicit first-order scheme, "
        f"got ratios {ratios}"
    )
    assert errs[-1] < 0.05 * abs(cont), (cont, modes)


def test_the_opened_iterations_add_modes_the_engine_does_not_have():
    """Two of the frame map's six modes are discretization artifacts: P45's lag and Eq. 74's
    carried flow.

    A numerical mode's `z` stays a fixed fraction per frame as the frame shrinks, so
    `lambda = ln(z)/dt` grows like 1/dt and has no continuous limit; a physical mode holds
    `lambda` and lets `z` approach 1. P45's lag exists because Eq. 26 takes the entering P45,
    evaluated before Eq. 80 recomputes it: the ordering within the frame, not convergence.

    Trap: to raise an iteration cap, wrap `_p45_loop`; setting `realtime.MAX_ITER_P45` does
    nothing, because the caps are bound as default arguments at definition time.
    """
    dts = (0.014, 0.007, 0.0035)
    on_record_p45 = (0.696, 0.609, 0.554)

    for dt, z_p45 in zip(dts, on_record_p45, strict=True):
        A = _frame_jacobian(1, dt=dt, tol=1e-12)
        lam_p45 = mode_of(A, "P45", dt)
        got_p45 = abs(np.exp(lam_p45 * dt))
        assert abs(got_p45 - z_p45) < 0.02, (
            f"dt={dt * 1000:g} ms: the P45-lag mode decays {got_p45:.3f} a frame, "
            f"{z_p45:.3f} on record"
        )
        # Eq. 74's carried-flow mode is pinned by its modulus only: it sits on a collision, and
        # tiny map changes flip it between complex, negative real and positive real.
        z_carry = np.linalg.eigvals(A)[
            int(np.argmax(participation(A, dt)[1][STATE_NAMES.index("WA31")]))
        ]
        assert abs(z_carry) < 0.05, (
            f"dt={dt * 1000:g} ms: Eq. 74's mode is z={z_carry}, |z|={abs(z_carry):.4f}; a "
            f"one-frame memory must decay almost completely within a frame"
        )
        # the mode really is the state it is named for
        lam, P = participation(A, dt)
        i45 = int(np.argmax(P[STATE_NAMES.index("P45")]))
        assert P[STATE_NAMES.index("P45"), i45] > 3.0 * P[STATE_NAMES.index("WA31"), i45], (
            f"dt={dt * 1000:g} ms: the mode called the P45 lag is no longer dominated by P45"
        )

    # Eq. 74's lambda*dt must be constant, asserted only over frames the model may run at:
    # [pdf p.38] caps the step at 10 ms, so the 14 ms frame is left out here.
    carries = []
    for dt in (0.010, 0.007, 0.0035):
        A = _frame_jacobian(1, dt=dt, tol=1e-12)
        z = np.linalg.eigvals(A)[int(np.argmax(participation(A, dt)[1][STATE_NAMES.index("WA31")]))]
        carries.append(complex(np.log(complex(z)) / dt).real * dt)
    assert max(carries) - min(carries) < 0.03 * abs(np.mean(carries)), (
        f"lambda*dt for Eq. 74's mode is {carries}; a per-frame memory must hold it fixed "
        f"while lambda itself grows like 1/dt"
    )

    coarse, fine = (
        _frame_jacobian(1, dt=dts[0], tol=1e-12),
        _frame_jacobian(1, dt=dts[-1], tol=1e-12),
    )
    lo = abs(mode_of(coarse, "P45", dts[0]))
    hi = abs(mode_of(fine, "P45", dts[-1]))
    assert hi > 3.0 * lo, (
        f"P45's mode moved {lo:.1f} -> {hi:.1f} /sec over a 4x frame refinement; a frame "
        f"artifact must grow like 1/dt and a physical mode must not move"
    )

    # for contrast, the physical NG mode
    lam_ng = mode_of(coarse, "NG", dts[0])
    assert abs(np.exp(lam_ng * dts[0])) > 0.95, (
        f"the NG mode decays {1 - abs(np.exp(lam_ng * dts[0])):.4f} a frame, so it is the "
        f"one with a continuous limit"
    )


def test_the_frame_map_has_exactly_one_mode_with_no_real_logarithm():
    """At most one negative real eigenvalue, and if there is one it belongs to an opened
    iteration or a volume pressure, never to a shaft speed. Whether one exists flips with
    small map changes, so that is not asserted.
    """
    for case in (1, 2, 3):
        for tol in (None, 1e-12):
            A = _frame_jacobian(case, tol=tol)
            z = np.linalg.eigvals(A)
            negative_real = [q for q in z if abs(q.imag) < 1e-12 and q.real < 0]
            assert len(negative_real) <= 1, (
                f"trim {case}, tol={tol}: expected at most one negative real eigenvalue, "
                f"got {negative_real}"
            )
            if not negative_real:
                continue  # trim 2 at the printed tolerance has none
            _, P = participation(A, FRAME_S)
            i = int(np.argmin(np.abs(z - negative_real[0])))
            owner = STATE_NAMES[int(np.argmax(P[:, i]))]
            assert owner in ("WA31", "P3", "P41"), (
                f"trim {case}, tol={tol}: the negative real mode is dominated by {owner}. "
                f"It should be one of the opened iterations or a volume pressure -- never "
                f"a shaft speed, which would mean a physical mode had lost its logarithm."
            )
            assert abs(negative_real[0]) < 0.05, (
                f"trim {case}, tol={tol}: the negative real eigenvalue is "
                f"{negative_real[0]:.4f}; on record it is a few thousandths, i.e. gone "
                f"inside one frame. A large one would be a real oscillation."
            )


def test_the_printed_iteration_tolerance_costs_a_measurable_amount():
    """The report stops its pressure loops at 0.1 percent in ten and eight passes
    [pdf p.37]. That is not free, and this says what it costs on the slow mode."""
    for case in (1, 2, 3):
        tight = discrete_ng_mode(case, tol=1e-12)
        printed = discrete_ng_mode(case)
        shift = 100.0 * (printed / tight - 1.0)
        assert -30.0 < shift < 0.0, (
            f"trim {case}: the printed tolerance moves the NG mode {shift:+.1f} % "
            f"({tight:.3f} -> {printed:.3f} /sec)"
        )


def test_the_descent_residual_is_gone_and_the_level_one_replaced_it():
    """Pins where the `f1` interpolation residual (#43) sits against Table 1 column 4: descent
    tight, level the loose one. Changing the interpolation scheme relocates it (see
    `docs/notes/derivative-ambiguity.md`), so neither trim may drift back unnoticed.
    """
    for case, lo, hi in ((3, -7.0, 7.0), (2, 5.0, 15.0)):
        cont = continuous_ng_mode(case)
        disc = discrete_ng_mode(case, tol=1e-12)
        for name, v in (("continuous", cont), ("discrete", disc)):
            dev = 100.0 * (abs(v) / abs(TABLE_1_2DOF[case]) - 1.0)
            assert lo < dev < hi, f"trim {case} {name}: {v:.3f} /sec, {dev:+.1f} %"
