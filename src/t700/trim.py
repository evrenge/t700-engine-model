"""Find the equilibrium of the gas generator at a given fuel flow.

Ballin trims by suppressing the integrations [pdf p.28]; algebraically, solve

    dNG/dt = dP3/dt = dP41/dt = dP45/dt = 0

for (NG, P3, P41, P45) with fuel flow, NP and ambient given. NP is not solved for: Eq. 47's
Qreq comes from the Gen Hel UH-60A simulation the report does not contain, so NP is held
and Qreq read off as Q_PT, which is Table B.1's "engine torque (ref. shaft)" [pdf p.67] and so a
validation target.

Damped Newton with a differenced Jacobian, which keeps SciPy out of the core.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from t700 import constants as c
from t700 import corrections, maps
from t700.engine import STANDARD_DAY, Ambient, Frame, State, frame
from t700.units import shp_from_torque


@dataclass
class TrimResult:
    """The outcome of a trim solve, converged or not."""

    state: State
    frame: Frame
    residual_converged: bool
    """The residual fell below tolerance. Not the same as a good answer (see `on_data`);
    named so that it does not read as "it worked"."""
    iterations: int
    residual_norm: float
    clamps: dict[str, int]
    """Map clamps accumulated over the WHOLE solve, Jacobian evaluations included. Mostly
    from intermediate iterates, so not evidence about the answer."""
    clamps_at_solution: dict[str, int]
    """Map clamps from one evaluation at the converged state. The one that matters:
    nonzero means the trim itself sits outside the digitized data."""

    @property
    def on_data(self) -> bool:
        """False when the converged state lies outside the compressor map's speed range.

        Clamping manufactures equilibria: off the data every lookup is constant, the
        residual stops depending on the state, and Newton can drive it to zero where the
        model means nothing. `f1:parameter` (corrected speed itself off the map) is the
        decisive clamp but not the only one: other tables may still be extrapolated, see
        `extrapolated_tables` and `fully_on_data`.
        """
        return "f1:parameter" not in self.clamps_at_solution

    @property
    def extrapolated_tables(self) -> tuple[str, ...]:
        """Every function table being read outside its own data at this trim, sorted.

        Which tables at which fuel flows is asserted in
        `validation/test_steady_sweeps.py::test_which_tables_are_extrapolated_at_which_trims`.
        """
        return tuple(sorted(self.clamps_at_solution))

    @property
    def fully_on_data(self) -> bool:
        """No function table is being extrapolated at the converged state.

        Stricter than `trustworthy`, deliberately: the other tables degrade accuracy
        without inventing a root.
        """
        return not self.clamps_at_solution

    @property
    def trustworthy(self) -> bool:
        """Converged and the compressor map's speed range is not clamped. Check this.

        Other tables may still be extrapolated (`extrapolated_tables`): they degrade
        accuracy without inventing a root, as an `f1:parameter` clamp does.
        """
        return self.residual_converged and self.on_data

    @property
    def shp(self) -> float:
        """Shaft horsepower delivered by the power turbine at this trim."""
        return shp_from_torque(self.frame.q_pt_ftlbf, self.state.np_rpm)

    @property
    def q_req_ftlbf(self) -> float:
        """The load torque this trim implies -- Q_PT, since dNP/dt = 0 at equilibrium."""
        return self.frame.q_pt_ftlbf

    def report(self) -> str:
        s, f = self.state, self.frame
        status = "converged" if self.residual_converged else "DID NOT CONVERGE"
        if self.residual_converged and not self.on_data:
            status = "CONVERGED OFF THE MAP -- not a physical solution"
        lines = [
            f"trim {status} in {self.iterations} iterations, |r| = {self.residual_norm:.3e}",
            f"  NG   {s.ng_rpm:10.1f} rpm   ({f.ngc_pct:6.2f} % corrected)",
            f"  NP   {s.np_rpm:10.1f} rpm",
            f"  P3   {s.p3_psia:10.3f} psia  P41 {s.p41_psia:8.3f}  P45 {s.p45_psia:8.3f}",
            f"  T3   {f.t3_degR:10.1f} R     T41 {f.t41_degR:8.1f}  T45 {f.t45_degR:8.1f}",
            f"  WA2  {f.wa2_pps:10.3f} lbm/s  FAR {f.far:8.5f}",
            f"  Qpt  {f.q_pt_ftlbf:10.2f} ft.lbf -> {self.shp:8.1f} shp",
        ]
        if self.clamps_at_solution:
            lines.append(
                f"  *** AT THE SOLUTION, outside digitized data: {self.clamps_at_solution}"
            )
            if self.on_data:
                lines.append(
                    "      (trustworthy is still True: f1's speed range is not clamped, "
                    "which is the clamp that manufactures equilibria. The rest degrade "
                    "accuracy. See TrimResult.extrapolated_tables.)"
                )
        elif self.clamps:
            lines.append(
                f"  (map clamps during the search only, none at the solution: {self.clamps})"
            )
        return "\n".join(lines)


def _residuals(x: np.ndarray, np_rpm: float, wf_pps: float, amb: Ambient) -> np.ndarray:
    st = State(float(x[0]), np_rpm, float(x[1]), float(x[2]), float(x[3]))
    f = frame(st, wf_pps, amb)
    return np.array([f.dng_dt, f.dp3_dt, f.dp41_dt, f.dp45_dt])


def _jacobian(
    x: np.ndarray,
    np_rpm: float,
    wf_pps: float,
    amb: Ambient,
    r0: np.ndarray,
    ng_hi: float | None = None,
) -> np.ndarray:
    """One-sided difference Jacobian, differencing inward at an active bound.

    Steps are relative: NG is order 4e4 rpm, P45 order 40 psia. At the top of `f1`'s
    speed range a forward probe lands off the map, where the NG column reads zero and the
    Jacobian is singular; stepping inward reads the real derivative. A forward step from
    the lower bound already points into the data.
    """
    n = x.size
    jac = np.zeros((n, n))
    for i in range(n):
        h = 1e-6 * max(abs(x[i]), 1.0)
        if i == 0 and ng_hi is not None and x[i] + h > ng_hi:
            h = -h
        xp = x.copy()
        xp[i] += h
        jac[:, i] = (_residuals(xp, np_rpm, wf_pps, amb) - r0) / h
    return jac


def initial_guess(np_rpm: float, amb: Ambient) -> State:
    """A starting point good enough for Newton to find the root.

    Scaled off the design point rather than off any particular trim, so the solver is not
    quietly handed the answer it is meant to find.
    """
    p2 = amb.p_amb_psia
    return State(
        ng_rpm=0.90 * c.NG_DES,
        np_rpm=np_rpm,
        p3_psia=10.0 * p2,
        p41_psia=9.5 * p2,
        p45_psia=3.0 * p2,
    )


def solve(
    wf_pps: float,
    np_rpm: float = c.NP_DES,
    ambient: Ambient = STANDARD_DAY,
    guess: State | None = None,
    tol: float = 1e-8,
    max_iter: int = 200,
) -> TrimResult:
    """Solve for the gas generator equilibrium at a fixed fuel flow.

    The residuals have different scales (dNG/dt rpm/sec, dP/dt psia/sec), so convergence
    is judged on a relative norm, each residual over a characteristic magnitude.

    NG is confined to `f1`'s own speed range by projecting each trial step onto it, so the
    degenerate root off the map is unreachable; past the 100 % speed line the solve fails
    rather than returning a number (open question #71).
    """
    amb = ambient
    st = guess or initial_guess(np_rpm, amb)
    x = np.array([st.ng_rpm, st.p3_psia, st.p41_psia, st.p45_psia], dtype=float)

    scale = np.array([c.NG_DES, amb.p_amb_psia, amb.p_amb_psia, amb.p_amb_psia])

    # Clamps are counted in `maps.clamp_scope`s, so a caller's enclosing count survives.
    # The feasible set for NG is `f1`'s own speed range (its first and last printed speed
    # line), carried into physical rpm through Eqs. 5-6.
    ngc_lo, ngc_hi = maps.f1().param_range
    root_theta2 = float(np.sqrt(corrections.theta2(amb.t_amb_degR)))
    ng_lo = ngc_lo / 100.0 * c.NG_DES * root_theta2
    ng_hi = ngc_hi / 100.0 * c.NG_DES * root_theta2

    with maps.clamp_scope() as during:
        r = _residuals(x, np_rpm, wf_pps, amb)
        started_off_map = not (ng_lo <= x[0] <= ng_hi)
        converged = False
        it = 0
        while it < max_iter:
            it += 1
            if float(np.linalg.norm(r / scale)) < tol:
                converged = True
                break
            jac = _jacobian(x, np_rpm, wf_pps, amb, r, None if started_off_map else ng_hi)
            try:
                step = np.linalg.solve(jac, -r)
            except np.linalg.LinAlgError:
                break

            # damped line search: accept the longest step that reduces the residual
            lam = 1.0
            r_new = r
            for _ in range(40):
                trial = x + lam * step
                # keep pressures positive and ordered enough for Eq. 18 to stay real
                if trial[1] <= trial[2] or np.any(trial <= 0.0):
                    lam *= 0.5
                    continue
                # and keep NG inside the map by projection, not refusal. Off the tables
                # every lookup is constant, so the torques fall as 1/NG and dNG/dt tends to
                # zero as NG grows: a spurious root that no residual form removes, because
                # the information is absent. Projection needs no constant; rejecting the
                # step would need a margin, which is a knob.
                if not started_off_map:
                    trial[0] = float(np.clip(trial[0], ng_lo, ng_hi))
                r_new = _residuals(trial, np_rpm, wf_pps, amb)
                if np.linalg.norm(r_new / scale) < np.linalg.norm(r / scale):
                    x = trial
                    break
                lam *= 0.5
            else:
                break
            r = r_new

    st = State(float(x[0]), np_rpm, float(x[1]), float(x[2]), float(x[3]))

    # Re-evaluate cleanly at the solution: clamps during a Newton search say nothing
    # about whether the answer sits inside the data.
    with maps.clamp_scope() as at_solution:
        f = frame(st, wf_pps, amb)

    return TrimResult(
        state=st,
        frame=f,
        residual_converged=converged,
        iterations=it,
        residual_norm=float(np.linalg.norm(r / scale)),
        clamps=during,
        clamps_at_solution=at_solution,
    )


def sweep(
    wf_pph_values,
    np_rpm: float = c.NP_DES,
    ambient: Ambient = STANDARD_DAY,
) -> list[TrimResult]:
    """Trim across a range of fuel flows, continuing from each solution to the next.

    Continuation is cheaper and better conditioned than a cold start; a convenience, not
    a correctness requirement. Check `trustworthy` on each result before using it.
    """
    from t700.units import wf_pps_from_pph

    out: list[TrimResult] = []
    guess: State | None = None
    for wf_pph in wf_pph_values:
        r = solve(wf_pps_from_pph(float(wf_pph)), np_rpm, ambient, guess=guess)
        if r.trustworthy:
            guess = r.state
        out.append(r)
    return out
