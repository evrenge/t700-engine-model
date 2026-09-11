"""Small-perturbation linear models: the five variants Ballin published.

Ballin extracted linear models from **two different nonlinear simulations**, each run
with and without the station 4.1 heat sink, and separately derived one by order
reduction. That is five distinct models, and they are not interchangeable -- comparing
ours against the wrong one is a category error, not a tolerance problem.

```
                          heat sink off        heat sink on
  full volume dynamics    5-DOF  B2,B4,B6      6-DOF  B8,B10,B12
  volume dyn. approx.     2-DOF  B1,B3,B5      3-DOF  B7,B9,B11
```

plus **reduced-order 5-DOF**, which appears in Table 1's third column only and is nowhere
in Appendix B.

## The one that trips people up

**The 2-DOF is not a reduction of the 5-DOF.** It is *extracted from a different
nonlinear simulation* -- the one that "approximate[s] the dynamics between the control
volumes to be instantaneous" [pdf p.67]. Only Table 1's third column is an order
reduction: "Order reduction was performed by setting the state derivatives P3-dot,
P41-dot and P45-dot to zero and solving the equations for the remaining two states"
[pdf p.28]. The same holds for 3-DOF versus 6-DOF: p.33 calls both "extracted", and the
report prints no order reduction of the 6-DOF anywhere.

So `TWO` and `REDUCED_FIVE` are genuinely different computations in Ballin's hands, and
Table 1 prints both -- at hover, -2.69 against -2.81.

**Ours coincide exactly, and that is a theorem, not a defect.** If `g(s) = f1(s, p(s))`
where `p(s)` solves `f2(s, p) = 0`, then implicit differentiation gives
`dg/ds = A11 - A12 A22^-1 A21` -- precisely the Schur complement that `REDUCED_FIVE`
forms. Linearizing an exactly-solved quasi-steady system and residualizing the full
Jacobian are the same operation. Measured, they agree to ~1e-12.

Ballin's two differ because his reduced-order model is **a separately coded nonlinear
program**, not the exact algebraic limit of his full one: it carries the opened
compressor iteration of Eq. 74, a finite inner-loop count, and the independent P45 solve,
none of which are exact. So the 0.12 /sec spread at hover measures *his* quasi-steady
numerics, and reproducing it would mean linearizing our `realtime.step` discrete map
rather than a Newton solve. That is a real and worthwhile test -- it would put a number
on what the real-time approximation costs -- but it needs a discrete-to-continuous
conversion and is deliberately **not** done here. Logged as open question #49.

The consequence for validation: **`TWO` cannot be checked independently against Table 1
column 4.** It carries the same information as `REDUCED_FIVE`, and column 5 is its proper
target.

## Which of our two model layers feeds which

Our code already mirrors Ballin's two nonlinear simulations, because it followed his
structure:

| Ballin | ours |
|---|---|
| complete nonlinear digital simulation | `engine.frame` -- pressures are states, Eqs. 42-44 |
| reduced-order nonlinear simulation | the quasi-steady system -- pressures solved algebraically |

`realtime.step` is a *discrete* implementation of the second. For linearization we want
the continuous quasi-steady system, which `_quasi_steady_deriv` builds directly: solve
the three pressures to equilibrium at fixed (NG, NP), then return the two speed
derivatives. Linearizing the discrete map instead would give a discrete-time Jacobian,
and Appendix B's matrices are continuous.

## Conventions

States and units follow the report [pdf pp.27, 32; nomenclature pp.11-13]: NG, NP in rpm,
P3/P41/P45 in psia, T41 in deg R, and the single input Wf in **lbm/sec** (Table B.1
quotes trim fuel flow in lbm/hr -- convert). So `A[i,j] = d(xi-dot)/d(xj)`.

## What is not reproducible here, and why

The **slowest mode is the NP/rotor mode** and it is not comparable at any DOF. Its
eigenvalue depends on dQreq/dNP and on the load inertia added to J_PT [Eq. 46], both
supplied by the Gen Hel UH-60A simulation, which this report consumes and does not
contain. We hold Qreq constant and default `j_load` to zero, so the derivative is zero and
the inertia is the bare turbine's. Appendix B flags the same row as not independently
reproducible. See `docs/notes/open-questions.md` #6.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from t700 import constants as c
from t700 import maps
from t700.engine import STANDARD_DAY, Ambient, Frame, State, frame


class DOF(StrEnum):
    """The five models. String-valued so `DOF("5dof")` and `DOF.FIVE` both work."""

    TWO = "2dof"
    """The quasi-steady nonlinear model, linearized. Figs. B1/B3/B5, Table 1 col. 4.

    Equals `REDUCED_FIVE` identically for us -- see the module docstring. Ballin's does
    not, because his is a separate program."""

    THREE = "3dof"
    """`TWO` plus the heat-sink state. Figs. B7/B9/B11. No printed eigenvalues exist."""

    FIVE = "5dof"
    """Extracted from the full-dynamics model. Figs. B2/B4/B6, Table 1 col. 3."""

    SIX = "6dof"
    """`FIVE` plus the heat-sink state. Figs. B8/B10/B12.

    **Never compare eigenvalues of this one.** The printed matrices are ill-conditioned
    by the report's own arithmetic -- `A(6,3)` and `A(6,4)` carry about +-50 absolute --
    and B8 is unstable as printed. Element-wise only; `appendix_b.Reference.ill_conditioned`
    flags it."""

    REDUCED_FIVE = "reduced5"
    """Order reduction of `FIVE` [pdf p.28]. Table 1 col. 5 only; not in Appendix B."""


STATES: dict[DOF, tuple[str, ...]] = {
    DOF.TWO: ("NG", "NP"),
    DOF.THREE: ("NG", "NP", "T41"),
    DOF.FIVE: ("NG", "NP", "P3", "P41", "P45"),
    DOF.SIX: ("NG", "NP", "P3", "P41", "P45", "T41"),
    DOF.REDUCED_FIVE: ("NG", "NP"),
}
"""State vectors. `FIVE` is printed [pdf p.27], `SIX` is printed [pdf p.32]. `TWO` is
stated in prose ("the two turbine-speed degrees-of-freedom", p.27) but not displayed, and
**`THREE`'s ordering is inferred, not printed** -- see `inventory-appendix-b.md` 2.4."""

HEAT_SINK = {DOF.THREE, DOF.SIX}
QUASI_STEADY = {DOF.TWO, DOF.THREE}


@dataclass(frozen=True)
class LinearModel:
    """`x-dot = A x + b Wf`, with the trim it was extracted at."""

    dof: DOF
    states: tuple[str, ...]
    A: np.ndarray
    b: np.ndarray
    x0: np.ndarray
    """The trim state the linearization is about, in the same order as `states`."""
    wf_pps: float
    np_rpm: float
    d: np.ndarray | None = None
    """Feedthrough from Wf [Eq. 67], `x = C z + d Wf` with `C = I`. Only the heat-sink
    models have one; `None` for 2-, 5- and reduced-5-DOF, which have no Wf-dot term."""

    @property
    def eigenvalues(self) -> np.ndarray:
        """Real parts, most negative first. Every mode in Table 1 is real and stable."""
        return np.sort(np.linalg.eigvals(self.A).real)

    @property
    def comparable_eigenvalues(self) -> np.ndarray:
        """Eigenvalues with the NP/rotor mode dropped -- see the module docstring.

        The NP mode is identified structurally rather than by position: column NP is
        otherwise all zeros in every printed variant, so `A[1,1]` *is* that eigenvalue.
        """
        ev = self.eigenvalues
        npq = self.A[1, 1]
        k = int(np.argmin(np.abs(ev - npq)))
        return np.delete(ev, k)


# ------------------------------------------------------------------ the two vector fields


def _full_frame(x: np.ndarray, wf: float, amb: Ambient, qreq: float, j_load: float, t41=None):
    return frame(State.from_array(x), wf, amb, q_req_ftlbf=qreq, j_load=j_load, t41_degR=t41)


def _full_deriv(x: np.ndarray, wf: float, amb: Ambient, qreq: float, j_load: float, t41=None):
    """Five derivatives from `engine.frame`: the complete nonlinear simulation.

    `t41` drives station 4.1 temperature instead of computing it, which is what makes it
    a sixth state rather than an algebraic result [pdf p.23, Eq. 23].
    """
    f = _full_frame(x, wf, amb, qreq, j_load, t41)
    return np.array([f.dng_dt, f.dnp_dt, f.dp3_dt, f.dp41_dt, f.dp45_dt])


def _solve_pressures(
    ng: float,
    np_rpm: float,
    p_guess: np.ndarray,
    wf: float,
    amb: Ambient,
    qreq: float,
    j_load: float,
    t41=None,
    tol: float = 1e-10,
    max_iter: int = 60,
) -> np.ndarray:
    """Solve P3, P41, P45 to equilibrium at fixed (NG, NP) -- the volume dynamics
    approximation, "instantaneous" between control volumes [pdf p.67].

    Newton on the three pressure residuals. This is the *continuous* statement of what
    Eqs. 74-80 compute iteratively each frame.
    """
    p = p_guess.astype(float).copy()

    def res(pv):
        return _full_deriv(np.array([ng, np_rpm, *pv]), wf, amb, qreq, j_load, t41)[2:]

    for _ in range(max_iter):
        r = res(p)
        if np.max(np.abs(r)) < tol * max(1.0, float(np.max(np.abs(p)))):
            return p
        J = np.zeros((3, 3))
        for i in range(3):
            h = 1e-6 * max(abs(p[i]), 1.0)
            pp = p.copy()
            pp[i] += h
            J[:, i] = (res(pp) - r) / h
        p = p - np.linalg.solve(J, r)
    raise RuntimeError(
        f"quasi-steady pressure solve did not converge at NG={ng:.0f}, NP={np_rpm:.0f}; "
        f"residual {np.max(np.abs(res(p))):.3e}"
    )


def _quasi_steady_deriv(
    s: np.ndarray,
    p_guess: np.ndarray,
    wf: float,
    amb: Ambient,
    qreq: float,
    j_load: float,
    t41=None,
):
    """Two speed derivatives with the pressures slaved -- the reduced-order simulation."""
    p = _solve_pressures(s[0], s[1], p_guess, wf, amb, qreq, j_load, t41)
    return _full_deriv(np.array([s[0], s[1], *p]), wf, amb, qreq, j_load, t41)[:2], p


# --------------------------------------------------------------------------- extraction


def _central_jacobian(f, x0: np.ndarray, n_out: int, rel: float) -> np.ndarray:
    """Central differences. `rel=1e-5` is inside the plateau where the 5-DOF eigenvalues
    are converged to 5 significant figures over rel = 1e-4 .. 1e-8, measured, not assumed.
    """
    A = np.zeros((n_out, len(x0)))
    for i in range(len(x0)):
        h = rel * max(abs(x0[i]), 1.0)
        xp, xm = x0.copy(), x0.copy()
        xp[i] += h
        xm[i] -= h
        A[:, i] = (f(xp) - f(xm)) / (2.0 * h)
    return A


def _heat_sink_taus(f: Frame) -> tuple[float, float]:
    """`(tau1, tau2)` at a trim -- the heat-sink lead and lag [nomenclature pdf p.13].

    Eq. 63 writes the linearized heat sink as `T41/T41_ns = (tau1 s + 1)/(tau2 s + 1)`;
    Eq. 50 gives the same transfer function built from the nonlinear model, so
    `tau2 = tau_a` (Eq. 51) and `tau1 = tau_a - tau_b` (Eqs. 52-53).

    **The report prints no numeric value for either** -- open question #31. They are
    trim-dependent by definition, "those values corresponding to the trim operating
    condition" [pdf p.31], so they can only come from our own Eqs. 50-53 at trim.
    """
    tau_a = c.TC_T41 * f.t41_ns_degR**0.5 / f.w41_pps**0.8  # (51)
    tau_b = float(maps.f_hs()(f.ngc_pct)) / f.w41_pps  # (52), (53)
    return tau_a - tau_b, tau_a


def _chen(F1, F2, G1, G2):
    """Eqs. 65-67: descriptor form to standard form.

    `F1 xdot = F2 x + G1 Wf + G2 Wfdot` is not state-space, because differentiating the
    lead term of Eq. 63 introduces `Wfdot`. With `B0 = F1^-1 G1` and `B1 = F1^-1 G2` it
    reads `xdot = A x + B0 Wf + B1 Wfdot`, and `z = x - B1 Wf` gives

        zdot = A z + (B0 + A B1) Wf        (66)
        x    = z + B1 Wf                   (67)  ->  C = I, d = B1

    which is why every printed `C` in Figures B7-B12 is the literal `[ I ]`.

    This also answers open question #32. The report writes `b = F1^-1 G1 + F F1^-1 G2`
    with `F` defined nowhere; the transformation gives `b = B0 + A B1`, so
    `F == A == F1^-1 F2`. **Derived here, not transcribed.**
    """
    A = np.linalg.solve(F1, F2)
    B0 = np.linalg.solve(F1, G1)
    B1 = np.linalg.solve(F1, G2)
    return A, B0 + A @ B1, B1


def _extract_heat_sink(dof, f_deriv, f_t41ns, z0, wf_pps, rel_step, taus, states, np_rpm):
    """Build one heat-sink linear model from its two vector fields.

    `z0` is the trim state with T41 last. `f_deriv(z, wf)` gives the `m` engine-state
    derivatives with T41 **driven**; `f_t41ns(z, wf)` gives T41_ns at the same point.

    The T41 row is Eq. 63 written out. With `cx = dT41_ns/dz` and `cw = dT41_ns/dWf`:

        tau2 T41dot + T41 = tau1 (cx.zdot + cw Wfdot) + (cx.z + cw Wf)

    Rearranged into `F1 zdot = F2 z + G1 Wf + G2 Wfdot`, `G2` comes out zero in every row
    but T41's -- which is what Appendix B prints [pdf p.33], and therefore a check on the
    derivation rather than an input to it.
    """
    tau1, tau2 = taus
    m = len(z0) - 1
    n = m + 1

    Jx = _central_jacobian(lambda z: f_deriv(z, wf_pps), z0, m, rel_step)
    cx = _central_jacobian(lambda z: np.array([f_t41ns(z, wf_pps)]), z0, 1, rel_step)[0]
    hw = rel_step * max(abs(wf_pps), 1e-6)
    Jw = (f_deriv(z0, wf_pps + hw) - f_deriv(z0, wf_pps - hw)) / (2.0 * hw)
    cw = (f_t41ns(z0, wf_pps + hw) - f_t41ns(z0, wf_pps - hw)) / (2.0 * hw)

    F1 = np.eye(n)
    F1[m, :m] = -tau1 * cx[:m]
    F1[m, m] = tau2 - tau1 * cx[m]

    F2 = np.zeros((n, n))
    F2[:m, :] = Jx
    F2[m, :m] = cx[:m]
    F2[m, m] = cx[m] - 1.0

    G1 = np.zeros(n)
    G1[:m] = Jw
    G1[m] = cw

    G2 = np.zeros(n)
    G2[m] = tau1 * cw

    A, b, d = _chen(F1, F2, G1, G2)
    return LinearModel(dof, states, A, b, z0, wf_pps, np_rpm, d=d)


def extract(
    trim_result,
    wf_pps: float,
    dof: DOF | str = DOF.FIVE,
    ambient: Ambient = STANDARD_DAY,
    j_load: float = 0.0,
    dq_req_dnp: float = 0.0,
    rel_step: float = 1e-5,
) -> LinearModel:
    """Extract one of the five linear models about a converged trim.

    Args:
        trim_result: a `trim.TrimResult`. Must be `trustworthy`; a clamped solution has a
            meaningless Jacobian, and the whole point of a linear model is its derivatives.
        wf_pps: the trim fuel flow, lbm/sec.
        dof: which of the five.
        j_load: load inertia added to J_PT [Eq. 46]. Zero means the bare power turbine.
        dq_req_dnp: slope of the load torque with power turbine speed, ft*lbf per rpm.
            `Q_req` is an **input** to this model [Eq. 47], supplied by Gen Hel; holding
            it constant means asserting `dQreq/dNP = 0`, which is not a neutral default
            but a statement that the rotor does not resist a speed change. It does. This
            completes the load specification rather than altering the engine: nothing in
            `engine.frame` changes, only what `Q_req` is at a perturbed NP. Recovered
            from Appendix B as +0.019471 / +0.015869 / +0.012950 at the three trims
            (open question #6); trim-dependent, so there is no single constant to store.

    The input column `b` is `d(x-dot)/d(Wf)`, differenced about the same trim.
    """
    dof = DOF(dof)
    if not trim_result.trustworthy:
        raise ValueError(
            "refusing to linearize an untrustworthy trim: "
            f"residual_converged={trim_result.residual_converged}, "
            f"on_data={trim_result.on_data}. A clamped solve has no usable derivatives."
        )
    amb = ambient
    full0 = trim_result.state.as_array()
    qreq0 = frame(trim_result.state, wf_pps, amb, j_load=j_load).q_pt_ftlbf
    np_rpm = float(full0[1])

    def qreq_at(np_now: float) -> float:
        """Load torque at a perturbed NP. Flat unless `dq_req_dnp` is supplied."""
        return qreq0 + dq_req_dnp * (np_now - np_rpm)

    qreq = qreq0  # the trim value, for the paths that do not perturb NP

    if dof in HEAT_SINK:
        f0 = frame(trim_result.state, wf_pps, amb, j_load=j_load)
        taus = _heat_sink_taus(f0)
        t41_0 = f0.t41_ns_degR  # at a trim T41 == T41_ns: Eq. 50 has unit DC gain

        if dof is DOF.SIX:
            z0 = np.array([*full0, t41_0])

            def f_deriv(z, wf):
                return _full_deriv(z[:5], wf, amb, qreq_at(z[1]), j_load, z[5])

            def f_t41ns(z, wf):
                return _full_frame(z[:5], wf, amb, qreq_at(z[1]), j_load, z[5]).t41_ns_degR
        else:
            z0 = np.array([full0[0], full0[1], t41_0])
            p_guess = full0[2:]

            def f_deriv(z, wf):
                return _quasi_steady_deriv(z[:2], p_guess, wf, amb, qreq_at(z[1]), j_load, z[2])[0]

            def f_t41ns(z, wf):
                q = qreq_at(z[1])
                p = _solve_pressures(z[0], z[1], p_guess, wf, amb, q, j_load, z[2])
                return _full_frame(np.array([z[0], z[1], *p]), wf, amb, q, j_load, z[2]).t41_ns_degR

        return _extract_heat_sink(
            dof, f_deriv, f_t41ns, z0, wf_pps, rel_step, taus, STATES[dof], np_rpm
        )

    if dof is DOF.FIVE:
        A = _central_jacobian(
            lambda x: _full_deriv(x, wf_pps, amb, qreq_at(x[1]), j_load), full0, 5, rel_step
        )
        hw = rel_step * max(abs(wf_pps), 1e-6)
        b = (
            _full_deriv(full0, wf_pps + hw, amb, qreq, j_load)
            - _full_deriv(full0, wf_pps - hw, amb, qreq, j_load)
        ) / (2.0 * hw)
        return LinearModel(dof, STATES[dof], A, b, full0, wf_pps, np_rpm)

    if dof is DOF.REDUCED_FIVE:
        # "setting the state derivatives P3-dot, P41-dot and P45-dot to zero and solving
        # the equations for the remaining two states" [pdf p.28]. The algebra itself is
        # not printed; this is the standard residualization and is checkable against
        # Table 1 column 5.
        m5 = extract(trim_result, wf_pps, DOF.FIVE, amb, j_load, dq_req_dnp, rel_step)
        A11, A12 = m5.A[:2, :2], m5.A[:2, 2:]
        A21, A22 = m5.A[2:, :2], m5.A[2:, 2:]
        b1, b2 = m5.b[:2], m5.b[2:]
        S = np.linalg.solve(A22, A21)
        t = np.linalg.solve(A22, b2)
        return LinearModel(dof, STATES[dof], A11 - A12 @ S, b1 - A12 @ t, full0[:2], wf_pps, np_rpm)

    # DOF.TWO -- extracted from the quasi-steady nonlinear model, NOT reduced from 5-DOF.
    s0, p0 = full0[:2], full0[2:]

    def g(s):
        d, _ = _quasi_steady_deriv(s, p0, wf_pps, amb, qreq_at(s[1]), j_load)
        return d

    A = _central_jacobian(g, s0, 2, rel_step)
    hw = rel_step * max(abs(wf_pps), 1e-6)
    bp, _ = _quasi_steady_deriv(s0, p0, wf_pps + hw, amb, qreq, j_load)
    bm, _ = _quasi_steady_deriv(s0, p0, wf_pps - hw, amb, qreq, j_load)
    return LinearModel(dof, STATES[dof], A, (bp - bm) / (2.0 * hw), s0, wf_pps, np_rpm)
