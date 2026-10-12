"""Linearizing the real-time map itself -- open question #49.

## The oddity this resolves

Table 1 [pdf p.31] prints two eigenvalue columns that ought to describe the same system:
column 4, a 2-DOF model, and column 5, an order-reduced 5-DOF. Ballin's differ -- -2.69
against -2.81 at hover, -2.23 against -2.16 at level, -1.82 against -1.83 at descent.
**Ours are identical**, to 5.0e-8 at `rel_step = 1e-5`. That was described
here as "the floor a central-difference Jacobian reaches", and it is not a floor -- it is
the agreement at one step. Sweeping the step over twelve values from 1e-3 to 1e-9, the
minimum is **8.0e-9 at rel = 3e-7**, six times smaller, with truncation error dominating
above and rounding below. The distinction matters because the number was being quoted as a
limit on how well two formulations *can* agree, which would make anything near it
unimprovable. (The exact minimum depends on the grid; what is robust is that it is
several times below the value at 1e-5.)

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

* **Two of the six modes are numerical, not physical**, and both are frame artifacts that
  vanish in the continuous limit: P45's one-frame lag (Eq. 26 takes the *entering* P45) at
  z = 0.675 / 0.590 / 0.538 per frame at dt = 14 / 7 / 3.5 ms, and Eq. 74's carried
  compressor flow as a complex pair at z = 0.154 / 0.171 / 0.183. Neither approaches 1.
  Until 2026-09-13 this file claimed one such mode and identified it positionally as
  `spectrum[-3]`; by participation factor that mode is P45's, not the carried flow's.
* **The descent trim keeps a ~22 % residual under every reading**, continuous or discrete,
  tight tolerance or printed. That one is not the real-time approximation; it is the
  interpolation ambiguity of open questions #31 and #43, which is worst at descent.
"""

from __future__ import annotations

import numpy as np
import pytest

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

DISCRETE_TOL_PCT = {1: 8.0, 2: 15.0, 3: 7.0}
"""Ours, declared in `SCOPE.md`. The loose one carries the #31/#43 interpolation residual,
which is present in the continuous model too and is not what this file is measuring.

**It was descent, and since 2026-09-13 it is level.** Putting `f1` on a beta grid
redistributed that residual: the three trims read about +6 / +14 / +6 % where they read
+1.4 / +2.5 / -21.8 %. Descent, the loose one at 25 %, is now among the tight ones; level
is the loose one. The worst case improved from 21.8 % to about 14 %, and the ambiguity moved
rather than closing. Open question #43."""


TANGENT = 1e-5
"""The continuous models here are tangents: the limit the discrete map converges to, and
the Jacobian the Schur-complement identity is about. The shipped extraction uses the
report's +-2 % (pdf p.28), a secant, which is neither."""

STATE_NAMES = ("NG", "NP", "P3", "P41", "P45", "WA31")
"""Order of the six states of the frame map, for participation factors."""


def participation(A: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """`(continuous eigenvalues, P)` where `P[k, i]` is state k's participation in mode i.

    **Decomposed from `A_d` itself, not from `logm(A_d)`.** `logm` is a function of A, so
    the two have the same eigenvectors and the same participation factors -- but `logm`
    fails outright when A has a negative real eigenvalue, and a per-mode `ln(z)/dt` does
    not. The mode that matters here is NG, whose `z` is real and positive at every trim, so
    routing it through a whole-matrix logarithm was buying a failure mode for nothing.

    It bought one on 2026-09-13: interpolating `f1` along Figure A1's construction lines
    moved a very fast mode to **z = -0.0033**, and `_real_log` refused the matrix. That
    mode decays by a factor of 300 in one frame and alternates sign -- it has no continuous
    limit and never did -- but it took the whole file down with it.

    Modes with a negative real `z` are returned as `nan`, which is the honest value: a real
    negative eigenvalue has no real logarithm and no continuous counterpart.

    `P = |V * inv(V).T|`, column-normalised -- the standard participation factor, which is
    invariant to the units of the states. That matters here: the six states are rpm, psia
    and lbm/s, so an eigenvector component says nothing on its own.

    **This file indexed the sorted spectrum positionally until 2026-09-13** -- `[-2]` for
    "the NG mode", `[-3]` for "the carried mass flow's mode" -- which is precisely what
    `tests/test_linear.py::test_sorting_would_mispair_the_slow_modes` forbids, and one of
    the two was wrong. See `test_the_opened_iterations_add_modes_the_engine_does_not_have`.
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
    return A


def discrete_ng_mode(case: int, dt: float = FRAME_S, tol: float | None = None) -> float:
    """The frame map's NG mode, identified by participation factor rather than by position.

    Every caller below used `discrete_spectrum(...)[-2]`. That happens to be right at all
    three trims and both tolerances -- checked -- but it is the positional pairing
    `tests/test_linear.py::test_sorting_would_mispair_the_slow_modes` exists to forbid,
    and its sibling `[-3]` was wrong. See `participation`.
    """
    return float(np.real(mode_of(_frame_jacobian(case, dt, tol), "NG", dt)))


def continuous_ng_mode(case: int) -> float:
    wf = wf_pps_from_pph(WF_PPH[case])
    r = trim.solve(wf, NP_RPM, AMB)
    m = extract(
        r, wf, dof=DOF.TWO, j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[case], rel_step=TANGENT
    )
    return float(np.sort(np.real(m.comparable_eigenvalues))[-1])


def test_the_two_continuous_models_are_identical_because_they_must_be():
    """The Schur complement of an exactly solved quasi-steady system *is* the reduction.

    This is the premise of the whole file: our columns 4 and 5 cannot differ, so any
    difference in Ballin's is about his numerics and not about the engine.
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
    """Refine dt and the discrete map must converge on the continuous eigenvalue.

    If it did not, the difference would be a structural error rather than the cost of the
    real-time approximation -- which is the distinction this whole file rests on.
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
    """Two of the frame map's six modes are discretization artifacts, and the file named
    the wrong state for one of them.

    ## The signature

    A numerical mode's `z` stays a fixed fraction per **frame** as the frame shrinks,
    rather than a fixed fraction per second. `lambda = ln(z)/dt` then runs to -infinity and
    the mode has no continuous limit. A physical mode holds `lambda` fixed and lets `z`
    approach 1.

    ## What was claimed, and what is there

    This test asserted one such mode, took it as `discrete_spectrum(...)[-3]`, and called
    it "Eq. 74's one-frame memory". Identified by participation factor instead of by
    position, that mode is **P45-dominated** -- P[P45] = 0.68 at the 7 ms frame against
    P[WA31] = 0.07. Eq. 74's carried flow lives in a different mode: the complex pair at
    z ~ 0.17, where P[WA31] = 0.36, shared with P3 and P41.

    The substance of the old claim survives: both are frame artifacts.

        dt (ms)    P45 lag        Eq. 74 carry
          14       z = 0.693      z = 0.0343
           7           0.605          0.0341
           3.5         0.551          0.0341

    Neither approaches 1, so neither has a continuous limit. Eq. 74's is the cleaner of
    the two: its z barely moves at all over a 4x refinement, so the memory is exactly one
    frame deep and lambda runs -240.8, -482.5, -965.0 /sec -- 1/dt to three figures.
    (Carry column restated 2026-09-15 with the traced-line `f1`. It has read a complex
    pair, a small negative real and a small positive real across three `f1` revisions; only
    the modulus has held, and only the modulus is asserted.)

    ## Why P45's lag is a mode at all, converged or not

    Because Eq. 26 takes the **entering** P45: `dh_gt = theta41 * f7(P45(n)/P41)` is
    evaluated before Eq. 80 recomputes P45. So P45 feeds forward a frame whatever the
    iteration does with it, exactly as Eq. 74 does for the compressor flow. Checked
    constructively: raising the P45 pass cap from 8 to 400 leaves this mode at -75.330
    **to the digit**, while raising the *inner* cap from 10 to 400 moves it to -81.649 --
    so it is the ordering within the frame, not the iteration's convergence.

    (Raising a cap has to be done by wrapping `_p45_loop`, not by setting
    `realtime.MAX_ITER_P45`: the caps are bound as default arguments at definition time,
    so assigning the module attribute does nothing. The audit that produced this file's
    corrections was itself caught by that.)
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
        # Eq. 74's carried-flow mode is pinned by its MODULUS, not by its sign, and the
        # difference is a correction of 2026-09-15. The sign is not a stable property of
        # this mode: the file recorded a complex pair at z ~ 0.17 until 2026-09-13, a
        # negative real z = -0.0033 after that day's `f1` interpolation change, and the
        # traced-line `f1` of 2026-09-15 puts it at a positive real +0.0343. It has also
        # been seen as a conjugate pair at z = 0.0425 +- 0.0075j under a map differing from
        # the committed one by 1.3e-04 lbm/sec. A mode whose character flips on a 1e-04
        # map change is sitting on a collision, and the sign was never the claim.
        #
        # What IS stable, and what the claim always was: |z| ~ 0.034 at every frame size,
        # so the memory is one FRAME deep and not one second deep. That makes lambda grow
        # like 1/dt and the mode has no continuous limit -- the same signature this test
        # already uses for the P45 lag, and now applied to Eq. 74 too instead of a sign.
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
            f"dt={dt * 1000:g} ms: the mode called the P45 lag is no longer dominated by "
            f"P45; this is the mis-attribution the test was rewritten to fix"
        )

    # Neither artifact has a continuous limit, and both are shown the same way: |lambda|
    # grows like 1/dt.
    #
    # Eq. 74's is the cleaner of the two, but the constancy of lambda*dt is asserted only
    # over frames the model may actually RUN at. [pdf p.38] sets a 10 ms maximum step, and
    # the 14 ms entry in `dts` is there to show frame dependence, not because anything
    # integrates there. Measured with the 2026-09-16 blended-edge `f1`:
    #
    #     dt (ms)      z        lambda /s    lambda*dt
    #       14.0    0.04396      -223.18      -3.1245     <- outside the report's limit
    #       10.0    0.03709      -329.44      -3.2944
    #        7.0    0.03596      -475.05      -3.3254
    #        3.5    0.03540      -954.60      -3.3411
    #
    # Over 10, 7 and 3.5 ms that is constant to 1.4 %. Including 14 ms it is 6.6 %, and
    # asserting over a frame the model is forbidden to use would be asserting about an
    # extrapolation. The |z| bound above still covers all three listed frames.
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
    """`np.real(logm(A))` is invalid whenever A has a negative real eigenvalue. A may have one.

    The 2026-09-13 numerical-mathematics audit found such an eigenvalue at trim 3 under the
    printed tolerance -- z = -0.0828 -- and noted that nothing guarded the operation. The
    P3/P41 stopping-rule fix removed it the same day; the `f1` interpolation change later
    that day brought one back at z = -0.0033, Eq. 74's carried mass flow; and the
    traced-line `f1` of 2026-09-15 has removed it again -- **there is now no negative real
    eigenvalue at any of the three trims under either tolerance**, and Eq. 74's mode sits
    at a positive real +0.0343.

    Whether one exists is therefore not a property worth asserting, and this test never did
    assert it: it bounds the count at one, and checks the character of any that appears.
    That is the durable statement, and it is why the file survived the mode changing sign
    three times. What WAS a defect was routing every mode through a whole-matrix logarithm
    to get at NG, whose own `z` is real and positive at every trim. `participation`
    decomposes `A_d` directly and takes `ln(z)/dt` per mode, so one mode without a real
    logarithm no longer takes the file down with it.

    This pins both halves: never more than one negative real eigenvalue, and if there is
    one it belongs to an opened iteration or a volume pressure, never to a shaft speed.
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
    """Trim 3 was ~22 % off Table 1 under every reading. It is +1.1 % now.

    The old claim, and it was right at the time: "whatever it is, it is not the real-time
    approximation -- it is the interpolation ambiguity of #31 and #43, which is worst at
    descent." Interpolating `f1` along Figure A1's own construction lines (see
    `maps.SpeedMap`) settles it: descent's continuous NG mode goes -1.409 -> **-1.841**
    against Table 1's -1.82, i.e. -22.6 % -> **+1.1 %**.

    **The ambiguity moved rather than closing.** Level went the other way, -0.3 % ->
    +11.8 %, so the worst case across the three trims improved from 22.6 % to 11.8 % and
    the residual is now at a different trim. That is exactly what
    `docs/notes/derivative-ambiguity.md` predicted: it "identifies the mechanism without
    identifying the scheme", and changing the scheme relocates the error. Open question #43
    is still open, and the report prints nothing that would close it.

    What this test now pins is the relocation, so that neither trim can quietly drift back.
    """
    for case, lo, hi in ((3, -7.0, 7.0), (2, 5.0, 15.0)):
        cont = continuous_ng_mode(case)
        disc = discrete_ng_mode(case, tol=1e-12)
        for name, v in (("continuous", cont), ("discrete", disc)):
            dev = 100.0 * (abs(v) / abs(TABLE_1_2DOF[case]) - 1.0)
            assert lo < dev < hi, f"trim {case} {name}: {v:.3f} /sec, {dev:+.1f} %"
