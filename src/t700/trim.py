"""Find the equilibrium of the gas generator at a given fuel flow.

Ballin trims by suppressing the integrations [pdf p.28]. The same thing stated
algebraically: at equilibrium the four gas-generator derivatives vanish, and we solve

    dNG/dt = dP3/dt = dP41/dt = dP45/dt = 0

for (NG, P3, P41, P45) with fuel flow, NP and ambient given.

**NP is not solved for.** Eq. 47 balances the power turbine torque against Qreq, and Qreq
comes from the Gen Hel UH-60A simulation, which this report consumes and does not contain.
So at a trim we hold NP at its printed value and read Qreq off as whatever Q_PT comes out
to be -- which is exactly the quantity Table B.1 prints as "engine torque (ref. shaft)",
and therefore a validation target rather than an input.

Damped Newton with a numerically differenced Jacobian. Deliberately plain: this is a
numerical utility, not physics, and it keeps SciPy out of the model core.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from t700 import constants as c
from t700 import maps
from t700.engine import STANDARD_DAY, Ambient, Frame, State, frame
from t700.units import shp_from_torque

# The four unknowns, in the order the residual vector uses.
UNKNOWNS = ("NG", "P3", "P41", "P45")


@dataclass
class TrimResult:
    """The outcome of a trim solve, converged or not."""

    state: State
    frame: Frame
    residual_converged: bool
    """The residual fell below tolerance. **This is not the same as a good answer** --
    see `on_data`. Named awkwardly on purpose: a plain `converged` reads as "it worked",
    and that is exactly the mistake this class exists to prevent."""
    iterations: int
    residual_norm: float
    residuals: np.ndarray
    clamps: dict[str, int]
    """Map clamps accumulated over the WHOLE solve, Jacobian evaluations included. Mostly
    from intermediate iterates, so not evidence about the answer."""
    clamps_at_solution: dict[str, int]
    """Map clamps from one evaluation at the converged state. THIS is the one that
    matters: nonzero means the trim itself sits outside the digitized data."""

    @property
    def on_data(self) -> bool:
        """False when the converged state lies outside the compressor map's speed range.

        Clamping protects against extrapolating a map, but it also **manufactures
        equilibria**: once the state leaves the data every lookup returns a constant, the
        residual stops depending on the state, and Newton happily drives it to zero in a
        region where the model means nothing. At 750 lbm/hr this reported
        `converged=True` with a residual of 6.6e-9 and NG = 2.16e11 rpm.

        `f1:parameter` is the decisive clamp because it means corrected speed itself is
        off the map -- distinct from a query running slightly past one speed line's last
        knot, which is ordinary and happens at the descent trim.

        **It is the decisive clamp, not the only one, and this property's name oversells
        it.** `on_data` is True at trims where three other tables are being extrapolated;
        use `extrapolated_tables` to see which, and `fully_on_data` to require none. The
        distinction was invisible until the 2026-09-13 numerical-mathematics audit pointed
        out that `f9` clamped -- the very ingredient open question #45 turns on -- passes
        this test without comment.
        """
        return "f1:parameter" not in self.clamps_at_solution

    @property
    def extrapolated_tables(self) -> tuple[str, ...]:
        """Every function table being read outside its own data at this trim, sorted.

        Measured across the fuel-flow range, so the shape of it is on record rather than
        discovered per-caller:

            115 lbm/hr   f8, f9             (NGc 65.0 %, the bottom of f1's parameter range)
            125          f8
            150          f8
            200-700      none
            775          f9                 (Ps9/P45 passes f9's tabulated 0.85012)

        **Two entries this table used to carry have ceased to exist**, and both are worth
        knowing about rather than quietly dropping. `f1@65` appeared at 150 lbm/hr while the
        map was interpolated at constant abscissa; on the beta grid every speed line ends at
        its own surge limit and the blend carries the limit with it, so the query no longer
        leaves the map. `f6` appeared at every trim above about 590 lbm/hr -- the top third
        of the power range, Figure 9's endpoint included -- and it is now loaded as the
        constant Figure A6 draws, which has no domain to leave. See
        `validation/test_steady_sweeps.py::test_which_tables_are_extrapolated_at_which_trims`,
        which is where the census is asserted rather than remembered.
        """
        return tuple(sorted(self.clamps_at_solution))

    @property
    def fully_on_data(self) -> bool:
        """No function table is being extrapolated at the converged state.

        Stricter than `trustworthy`, and deliberately not what `trustworthy` means: `f6`
        clamps at every trim above ~590 lbm/hr, so requiring this would reject the top
        third of the power range including Figure 9's own endpoint.
        """
        return not self.clamps_at_solution

    @property
    def trustworthy(self) -> bool:
        """Converged *and* the compressor map's speed range is not clamped.

        This is what callers should check, and what it does **not** check is every other
        table -- see `extrapolated_tables`. The line is drawn at `f1:parameter` because
        that clamp manufactures equilibria (above), while the others degrade accuracy
        without inventing a root.
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
    x: np.ndarray, np_rpm: float, wf_pps: float, amb: Ambient, r0: np.ndarray
) -> np.ndarray:
    """Forward-difference Jacobian.

    Steps are relative, because NG is order 4e4 rpm while P45 is order 40 psia -- one
    absolute step cannot serve both.
    """
    n = x.size
    jac = np.zeros((n, n))
    for i in range(n):
        h = 1e-6 * max(abs(x[i]), 1.0)
        xp = x.copy()
        xp[i] += h
        jac[:, i] = (_residuals(xp, np_rpm, wf_pps, amb) - r0) / h
    return jac


def initial_guess(wf_pps: float, np_rpm: float, amb: Ambient) -> State:
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

    The residuals have wildly different scales -- dNG/dt is rpm/sec, dP/dt is psia/sec --
    so convergence is judged on a *relative* norm, each residual divided by a
    characteristic magnitude of its own equation.
    """
    amb = ambient
    st = guess or initial_guess(wf_pps, np_rpm, amb)
    x = np.array([st.ng_rpm, st.p3_psia, st.p41_psia, st.p45_psia], dtype=float)

    scale = np.array([c.NG_DES, amb.p_amb_psia, amb.p_amb_psia, amb.p_amb_psia])

    # Scoped rather than `reset_clamps()`: a caller measuring clamps over a run that
    # contains trims -- which `validation/` does -- used to have its count silently zeroed
    # by every trim inside it. See `maps.clamp_scope`.
    with maps.clamp_scope() as search:
        r = _residuals(x, np_rpm, wf_pps, amb)
        converged = False
        it = 0
        while it < max_iter:
            it += 1
            if float(np.linalg.norm(r / scale)) < tol:
                converged = True
                break
            jac = _jacobian(x, np_rpm, wf_pps, amb, r)
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
                r_new = _residuals(trial, np_rpm, wf_pps, amb)
                if np.linalg.norm(r_new / scale) < np.linalg.norm(r / scale):
                    x = trial
                    break
                lam *= 0.5
            else:
                break
            r = r_new
    during = search.counts

    st = State(float(x[0]), np_rpm, float(x[1]), float(x[2]), float(x[3]))

    # Re-evaluate cleanly at the solution: clamps during a Newton search say nothing
    # about whether the answer sits inside the data.
    with maps.clamp_scope() as log:
        f = frame(st, wf_pps, amb)
    at_solution = log.counts

    return TrimResult(
        state=st,
        frame=f,
        residual_converged=converged,
        iterations=it,
        residual_norm=float(np.linalg.norm(r / scale)),
        residuals=r,
        clamps=during,
        clamps_at_solution=at_solution,
    )


def sweep(
    wf_pph_values,
    np_rpm: float = c.NP_DES,
    ambient: Ambient = STANDARD_DAY,
) -> list[TrimResult]:
    """Trim across a range of fuel flows, continuing from each solution to the next.

    Continuation matters more than it looks. From a fixed design-point guess the solver
    wanders off the compressor map above about 700 lbm/hr and lands on a spurious
    equilibrium -- clamped maps make the residual insensitive to the state, so Newton can
    drive it to zero somewhere meaningless. Walking the guess along the sweep keeps every
    step inside the data, and the sweep then reaches 800 lbm/hr (NG 99.87 %) and refuses
    at 850, which is genuinely past the map's 100 % speed line.

    Check `trustworthy` on each result before using it.
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
