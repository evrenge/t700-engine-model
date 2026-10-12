"""The real-time frame: Eqs. 69-80, the model Ballin actually ran.

`engine.frame` carries the three volume pressures as states (Eqs. 42-44): the 5-DOF model
Appendix B linearizes and the trim solver uses. It cannot be integrated at the report's
7 ms frame: its fastest mode makes explicit integration stable only below about 0.4 ms. So
Ballin makes the three pressures algebraic and solves them each frame:

* the opened compressor mass-flow iteration (Eq. 74) carries one mass flow across the frame
  boundary instead of looping through the compressor map;
* the inner P3/P41 loop (Eqs. 76, 78) is a Gauss-Seidel sweep -- "eleven arithmetic
  operations for each pass ... no relaxation algorithm is required for convergence"
  [pdf p.37];
* P45 is solved independently (Eq. 80), "because it is a nonlinear function of its own
  value, iterative techniques must be used" [pdf p.37].

Only NG and NP remain integrated; their modes are hundreds of milliseconds, stable at the
report's 7 and 14 ms. The dropped fast modes are the inter-volume dynamics the Conclusions
call "unnecessary" [pdf p.54].

Two decisions, not transcriptions:

* The integrator is never named (open question #26). Explicit Euler, what a 1988 real-time
  model with "no iteration between time steps" [pdf p.35] would use.
* Eq. 74 reads two ways (open question #22). As printed it lags only the bleed,
  `WA31_(n) = WA3_(n) - WA3_bl_(n-1)`; the prose takes "the mass flow leaving the
  compressor in the previous interval", `WA31_(n) = WA31_(n-1)`. Both are implemented
  (`lag_whole_flow`); the default is the prose, the sentence that describes opening the
  loop. The two differ by less than the report's own 0.1 % criterion.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from math import sqrt
from typing import Final

import numpy as np

from t700 import constants as c
from t700 import maps
from t700.engine import (
    STANDARD_DAY,
    Ambient,
    combustor,
    compressor,
    gg_turbine_and_station45,
    heat_sink_taus,
    power_turbine,
    shaft_change,
    station41,
    torques,
)

FRAME_ENGINE_S = 0.007
"""Engine update interval [pdf p.47]: "updated twice for each rotor routine cycle, or once
every 7 msec"."""

FRAME_NP_S = 0.014
"""NP update interval [pdf p.47]. NP is the drive-train degree of freedom and moves at the
host rotor rate, so the model is multirate 2:1 -- "the engine model is updated twice for
each rotor routine cycle, or once every 7 msec". `run` defaults to it."""

MAX_STEP_S = 0.010
"""Maximum time step [pdf p.38], set by a 0.1 % maximum allowable error between steps.

Smaller than FRAME_NP_S; the report states both (open question #33). Recorded, not
enforced: the report's own configuration runs NP at 14 ms, so a check would reject it."""


@dataclass(frozen=True)
class RTState:
    """Five model states plus the one scalar the opened iteration carries across frames."""

    ng_rpm: float
    np_rpm: float
    p3_psia: float
    p41_psia: float
    p45_psia: float
    wa31_carry_pps: float
    """The opened compressor mass-flow iteration's memory (Eq. 74). Its initialization is
    not stated in the report (open question #23)."""

    hs_tm_degR: float = 0.0
    """Station 4.1 metal temperature, deg R: Eq. 48's state (see `_heat_sink` for why).
    Zero means not seeded, and `_heat_sink` starts it at equilibrium. Unused when
    `heat_sink=False`, as are the two carries below."""

    t41_carry_degR: float = 0.0
    """Previous frame's T41, deg R: gas, not metal, temperature, the report's sixth state
    [pdf p.29, p.32]. Eq. 51 needs it before this frame has computed it."""

    w41_carry_pps: float = 0.0
    """Previous frame's W41, lbm/sec, for Eqs. 51 and 53. W41 depends on T41, which the
    heat sink is computing, so it is lagged one frame as Eq. 74 opens the mass-flow loop."""


@dataclass(frozen=True)
class FrameOut:
    """What one real-time frame produced, beyond the next state."""

    inner_iters: int
    p45_iters: int
    inner_exit: Exit
    """How the P3/P41 sweep finished."""
    p45_exit: Exit
    """How the P45 iteration finished. `EXACT` below roughly 76 %NG is expected, not a
    fault -- see `_p45_loop`."""
    wa31_substituted: bool
    """True when Eq. 74's carried mass flow came back non-positive and was replaced. See
    `step`."""
    t41_degR: float
    t41_ns_degR: float
    """T41 before the heat sink. Equal to `t41_degR` when `heat_sink=False` [Eq. 23]."""
    t45_degR: float
    q_pt_ftlbf: float
    wa2_pps: float
    w45_pps: float
    far: float


# --------------------------------------------------------------------------- solver stopping rules

TOL_PRESSURE: Final = 1.0e-3
"""Relative error bound for both pressure iterations. [pdf p.37]

Printed twice on that page: the P3/P41 sweep converges "with less than 0.1 percent
error", and for P45 "eight iterations resulted in an error equal to less than 0.1 percent
of the steady-state value". Both state an error, not a convergence test, so the loops stop
on the error estimated from the iterates (`_contraction_error`), not on the step: with
contraction rho near 0.9 the two differ by rho/(1-rho), about eight times.
"""


class Exit(StrEnum):
    """How an iteration finished. Recorded because a model that stops early is a model
    running on unconverged pressures, and nothing downstream can tell from the value."""

    CONVERGED = "converged"
    """The estimated error fell below `tol` before the printed pass cap."""

    CAPPED = "capped"
    """The pass cap was reached first -- the report's "under the most extreme conditions,
    up to ten iterations may be required" [pdf p.37]. Expected during a transient."""

    EXACT = "exact"
    """Eq. 80's printed iteration did not reach the printed criterion, so the same equation
    was solved exactly by `_p45_exact`. See `_p45_loop` and open question #57."""

    DIVERGING = "diverging"
    """The iteration did not contract: successive steps grew, or held level in a cycle.
    Only the P3/P41 sweep reports it (`_p45_loop` falls back to an exact solve).
    `_stalled` decides it and does not separate divergence from a limit cycle."""


def _contraction_error(step: float, prev_step: float, value: float) -> float:
    """A-posteriori error estimate for a linearly convergent fixed-point iteration.

    For `x_(k+1) = g(x_k)` with `|g'| = rho < 1` near the fixed point,
    `|x_k - x*| <= rho/(1-rho) * |x_k - x_(k-1)|`, with rho estimated from consecutive
    steps: the estimate the report's "the iteration converges linearly" [pdf p.37] licenses.

    Infinite with no previous step (so both loops take at least two passes) or when
    `rho >= 1` (not contracting: no estimate). Zero when the step is below
    `_CONVERGED_STEP * value`: the iterate has stopped moving, and a ratio of two
    rounding-level steps is noise that would report failure at a fixed point. A held trim
    takes that branch, in one pass.
    """
    if step <= _CONVERGED_STEP * abs(value):
        return 0.0
    if prev_step <= 0.0:
        return float("inf")
    rho = step / prev_step
    if rho >= 1.0:
        return float("inf")
    return rho / (1.0 - rho) * step


def _stalled(first_step: float, last_step: float) -> bool:
    """Did the iteration fail to contract -- moving as far at the end as at the start?

    Compares the first and last steps, not any single ratio: near a fixed point the steps
    reach the rounding floor and their ratios are noise. `>=` against a relative tolerance,
    not `>`: a limit cycle (f9's expansive region against its clamped plateau) has
    `last == first` exactly, and counts as a failure to contract.
    """
    if first_step <= 0.0:
        return False
    return last_step >= first_step * (1.0 - 1e-9)


MAX_ITER_P3_P41: Final = 10
"""Pass cap for the P3/P41 sweep. [pdf p.37]  "up to ten iterations may be required",
cross-checked by the printed operation count: 11 operations a pass x 10 = 110."""

_CONVERGED_STEP: Final = 1.0e-9
"""Relative step below which an iterate is at its fixed point, whatever the step ratio says.

A noise floor for the estimator, not a tolerance: six orders below `TOL_PRESSURE`, so it
cannot mask a real failure to converge. See `_contraction_error`.
"""

MAX_ITER_P45: Final = 8
"""Pass cap for the P45 iteration. [pdf p.37]  "Eight iterations resulted in an error
equal to less than 0.1 percent of the steady-state value"."""


def _inner_pressure_loop(
    p3: float,
    p41: float,
    t3: float,
    theta41: float,
    wa31: float,
    wf: float,
    tol: float = TOL_PRESSURE,
    max_iter: int = MAX_ITER_P3_P41,
) -> tuple[float, float, int, Exit]:
    """The P3/P41 fixed-point sweep, Eqs. 76 and 78.

    Gauss-Seidel as the report describes it [pdf p.37]: P3 from the current P41, then P41
    from the P3 just computed, not a simultaneous solve. T3 and theta_41 are held across
    the sweep (recomputed once per frame outside), which makes it eleven arithmetic
    operations a pass.

    The stopping rule is the report's: "up to ten iterations may be required for
    convergence with less than 0.1 percent error, resulting in a total of 110 arithmetic
    operations" [pdf p.37]; eleven a pass times ten pins the cap. Converging harder than
    Ballin did sharpens every transient, though no trim moves (open question #25). The
    test is on the error, via `_contraction_error`; see `TOL_PRESSURE`.
    """
    # T3, theta_41 and WA31 are fixed across the sweep, so products of them alone are
    # formed once, grouped in the printed expressions' left-to-right order so every pass
    # sees the same bits.
    sqrt_theta41 = sqrt(theta41)
    dpb_t3 = t3 * c.K_DPB
    four_dpb_t3_wa31_sq = 4.0 * c.K_DPB * t3 * wa31 * wa31
    prev3 = prev41 = 0.0
    first3 = first41 = 0.0
    exit_ = Exit.CAPPED
    it = 0
    while it < max_iter:
        it += 1
        p3_new = 0.5 * (p41 + sqrt(p41 * p41 + four_dpb_t3_wa31_sq))  # (76)
        w41 = c.K_WGT * p41 / sqrt_theta41
        p41_new = p3_new - dpb_t3 * (w41 - wf) ** 2 / p3_new  # (78) with (77)
        step3, step41 = abs(p3_new - p3), abs(p41_new - p41)
        p3, p41 = p3_new, p41_new
        if it == 1:
            first3, first41 = step3, step41
        err3 = _contraction_error(step3, prev3, p3)
        err41 = _contraction_error(step41, prev41, p41)
        if err3 < tol * p3 and err41 < tol * p41:
            exit_ = Exit.CONVERGED
            break
        prev3, prev41 = step3, step41
    else:
        if _stalled(first3, prev3) or _stalled(first41, prev41):
            exit_ = Exit.DIVERGING
    return float(p3), float(p41), it, exit_


def _p45_loop(
    p45: float,
    w41: float,
    b3: float,
    wa2: float,
    theta45: float,
    ps9: float,
    tol: float = TOL_PRESSURE,
    max_iter: int = MAX_ITER_P45,
) -> tuple[float, int, Exit]:
    """P45 by its own iteration, Eq. 80.

    "Eight iterations resulted in an error equal to less than 0.1 percent of the
    steady-state value" [pdf p.37], measured under "an instantaneous step in fuel flow from
    flight idle to full power", Figure 9's own test case. The numerator is constant over the
    iteration; only the f9 lookup changes ("one function-table lookup and four arithmetic
    operations" a pass [pdf p.37]).

    Eq. 79 prints the returned bleed as `B3 B4 WA2` where Eq. 44 writes `B3 K_bl WA2` for
    the identical term, and B4 is in no nomenclature (open question #30). K_bl is used.

    Not always contractive: `g'(P45*) = dln(f9)/dln(Ps9/P45)` is f9's elasticity, which
    crosses -1 near `Ps9/P45 ~ 0.77`, and above that the fixed point repels (open question
    #45). So the printed iteration runs first, exactly as printed, and returns if it meets
    the printed criterion. If not, the same equation is solved exactly by `_p45_exact`,
    marked `Exit.EXACT` and counted as one extra pass. A deliberate departure: it
    happens below flight idle, a regime the real-time model dropped [pdf p.38] and Ballin
    never measured, where the printed method returns a wrong root (open question #57).
    """
    numerator = (w41 + b3 * c.K_BL * wa2) * sqrt(theta45)
    f9 = maps.f9()

    # --- Eq. 80 exactly as printed, first -------------------------------------------
    prev = 0.0
    it = 0
    while it < max_iter:
        it += 1
        p45_new = numerator / float(f9(ps9 / p45))
        step = abs(p45_new - p45)
        p45 = p45_new
        err = _contraction_error(step, prev, p45)
        if err < tol * p45:
            return float(p45), it, Exit.CONVERGED
        prev = step

    # --- it did not reach the printed criterion, so solve the printed equation -------
    root = _p45_exact(numerator, ps9, f9)
    f9(ps9 / root)  # one lookup, so a root off the end of f9's table still counts a clamp
    return float(root), it + 1, Exit.EXACT


def _p45_exact(numerator: float, ps9: float, f9) -> float:
    """Solve `P45 = N / f9(Ps9/P45)` exactly, on f9's own piecewise-linear table.

    With `u = Ps9/P45`, Eq. 80 becomes

        Ps9/u = N/f9(u)     <=>     f9(u) = (N/Ps9) * u

    a piecewise-linear function against a line through the origin, so the root is where the
    line crosses one of f9's segments: one division once the segment is known.
    `phi(u) = f9(u) - m*u` is strictly decreasing, because f9 is conditioned non-increasing
    (`maps.f9`) and `m = N/Ps9 > 0`, so the root is unique and a binary search over the
    knots finds its segment. Off either end f9 is its clamped constant and the root is
    `y_end/m`, still exact.
    """
    m = numerator / ps9
    if not f9.non_increasing:  # maps.f9 conditions it monotone="dec"; uniqueness needs it
        raise ValueError("f9 is not non-increasing, so Eq. 80's root is no longer unique")
    xl, yl, n = f9._xl, f9._yl, f9._n
    if yl[0] - m * xl[0] <= 0.0:  # the root sits on the left clamped extension
        return ps9 / (yl[0] / m)
    if yl[n - 1] - m * xl[n - 1] >= 0.0:  # ... or on the right one
        return ps9 / (yl[n - 1] / m)
    lo, hi = 0, n - 1
    while hi - lo > 1:
        mid = (lo + hi) >> 1
        if yl[mid] - m * xl[mid] >= 0.0:
            lo = mid
        else:
            hi = mid
    x0, y0 = xl[lo], yl[lo]
    s = (yl[lo + 1] - y0) / (xl[lo + 1] - x0)
    return ps9 / ((s * x0 - y0) / (s - m))


def _heat_sink(
    t41_ns_degR: float,
    tm_degR: float,
    t41_prev_degR: float,
    w41_prev_pps: float,
    ngc_pct: float,
    dt: float,
) -> tuple[float, float]:
    """Advance the station 4.1 heat sink one frame. Returns `(T41, next metal temp)`.

    ```
    T41       (M c_pm/(h A_m) - M c_pm/(W_g c_pg)) s + 1
    ------ = --------------------------------------------          (50)  [pdf p.26]
    T41_ns              M c_pm/(h A_m) s + 1

    M c_pm/(h A_m)     = TC_T41 * sqrt(T41) / W41^(4/5)             (51)
    T41_sgn            = f_hs(NG_c)                                 (52)
    M c_pm/(W_g c_pg)  = T41_sgn / W41                              (53)
    ```

    Write `tau_a` for Eq. 51 and `tau_b` for Eq. 53, so Eq. 50 is
    `((tau_a - tau_b) s + 1) / (tau_a s + 1)`.

    The state is the metal temperature, and that is load-bearing. Eq. 50 collapses the two
    printed heat-transfer equations [pdf p.25]:

    ```
    c_pm M dT_m/dt      = h A_m (T_gi - T_m)                        (48)
    c_pg W_g (T_gi-T_go) = c_pm M dT_m/dt                            (49)
    ```

    With `T_gi = T41_ns` and `T_go = T41`, those are exactly

    ```
    dT_m/dt = (T41_ns - T_m) / tau_a
    T41     = T41_ns - (tau_b/tau_a) * (T41_ns - T_m)
    ```

    which is what this integrates. The collapse into Eq. 50 holds only for constant
    coefficients, and Eqs. 51 and 53 vary with T41 and W41: a lead-lag memory has the
    high-frequency gain `k` baked in and goes stale when `k` moves on a fuel step, while
    `T_m` is a physical temperature whatever the coefficients do. The two agree for
    constant coefficients, so Appendix B and Table 1 are untouched; the linear models still
    carry gas temperature as the sixth state [pdf p.29].

    DC gain is exactly 1 (T41 = T41_ns at equilibrium), so the heat sink cannot move a trim
    and `engine.frame` and `trim.solve` need no heat-sink flag. The high-frequency gain is
    `k = 1 - tau_b/tau_a`: fast T41 excursions are attenuated, the "response is
    significantly slowed" of [pdf p.33].
    """
    tau_a, tau_b = heat_sink_taus(t41_prev_degR, w41_prev_pps, ngc_pct)  # (51)-(53)
    tm = tm_degR if tm_degR > 0.0 else t41_ns_degR  # unseeded: start at equilibrium
    t41 = t41_ns_degR - (tau_b / tau_a) * (t41_ns_degR - tm)  # (48) with (49)
    tm_next = tm + dt * (t41_ns_degR - tm) / tau_a  # (48)
    return t41, tm_next


def step(
    st: RTState,
    wf_pps: float,
    ambient: Ambient = STANDARD_DAY,
    dt: float = FRAME_ENGINE_S,
    q_req_ftlbf: float = 0.0,
    j_load: float = 0.0,
    integrate_np: bool = True,
    lag_whole_flow: bool = True,
    heat_sink: bool = False,
    tol: float = TOL_PRESSURE,
    dt_np: float | None = None,
) -> tuple[RTState, FrameOut]:
    """Advance one engine frame.

    Args:
        lag_whole_flow: how to read Eq. 74. True (default) follows the report's prose --
            the combustor inlet flow is the compressor exit flow of the previous interval.
            False follows the equation as printed, lagging only the bleed. Open question
            #22; see the module docstring.
        integrate_np: False suppresses the NP integration, which is what the report does
            for the open-loop fuel steps of Figs. 9 and 10 [pdf p.39].
        dt_np: the NP integration step, when it differs from the engine frame. `None`
            means `dt`, the single-rate model. `run` uses it for the report's 2:1 multirate;
            see `FRAME_NP_S`.
    """
    # --- the compressor, evaluated ONCE on the frame's entering pressure ---------------
    cp = compressor(st.ng_rpm, st.p3_psia, ambient)  # (1)-(17)

    # --- Eq. 74, the opened iteration -------------------------------------------------
    # The guard is a silent substitution, like Eq. 18's in `engine.frame`: a non-positive
    # carried flow is replaced by this frame's equilibrium estimate, floored at 1e-6 lbm/s
    # so Eq. 19 cannot divide by zero. The floor is ours, not the report's; it guards a
    # state that should not occur, and `wa31_substituted` keeps that checkable.
    wa31 = st.wa31_carry_pps if lag_whole_flow else (cp.wa3 - st.wa31_carry_pps)
    wa31_substituted = wa31 <= 0.0
    if wa31_substituted:
        wa31 = max(cp.wa3 - cp.wa3_bl, 1e-6)

    # --- combustor, once per frame ----------------------------------------------------
    far = wf_pps / wa31  # (19)
    _, _, t41_ns = combustor(far, cp.h3)  # (20)-(22)
    if heat_sink:
        t41, hs_tm = _heat_sink(
            t41_ns, st.hs_tm_degR, st.t41_carry_degR, st.w41_carry_pps, cp.ngc_pct, dt
        )  # (48)-(53)
    else:
        t41, hs_tm = t41_ns, st.hs_tm_degR  # (23), no heat-sink representation
    h41, theta41 = station41(t41)  # (24), (25)

    # --- the two pressure solves ------------------------------------------------------
    p3, p41, inner_iters, inner_exit = _inner_pressure_loop(
        st.p3_psia, st.p41_psia, cp.t3, theta41, wa31, wf_pps, tol=tol
    )

    ps9 = ambient.p_amb_psia  # (37)
    p49 = ps9 * float(maps.f10()(cp.ngc_pct))  # (38)

    dh_gt, _, w41, _, t45, theta45 = gg_turbine_and_station45(
        h41, theta41, p41, st.p45_psia
    )  # (26)-(31)

    p45, p45_iters, p45_exit = _p45_loop(st.p45_psia, w41, cp.b3, cp.wa2, theta45, ps9, tol=tol)

    dh_pt, _, w45 = power_turbine(theta45, p45, p49, ps9)  # (32)-(34)

    # --- torques and the two surviving integrations -----------------------------------
    q_c, q_gt, q_pt = torques(st.ng_rpm, st.np_rpm, cp, w41, dh_gt, w45, dh_pt)  # (39)-(41)

    ng = st.ng_rpm + shaft_change(q_gt, q_c, c.J_GT, dt)  # (45)
    np_ = st.np_rpm  # np_ is NP: np is NumPy
    if integrate_np:
        np_ = st.np_rpm + shaft_change(
            q_pt, q_req_ftlbf, c.J_PT + j_load, dt if dt_np is None else dt_np
        )  # (46), (47)

    # The carried scalar differs between the readings. Eq. 74 as printed lags the bleed,
    # `WA31(n) = WA3(n) - WA3_bl(n-1)`; carrying WA31 there instead gives an involution
    # that oscillates about WA3/2 (open question #22).
    carry = (cp.wa3 - cp.wa3_bl) if lag_whole_flow else cp.wa3_bl
    nxt = RTState(ng, np_, p3, p41, p45, float(carry), hs_tm, t41, float(w41))
    out = FrameOut(
        inner_iters=inner_iters,
        p45_iters=p45_iters,
        inner_exit=inner_exit,
        p45_exit=p45_exit,
        wa31_substituted=wa31_substituted,
        t41_degR=t41,
        t41_ns_degR=t41_ns,
        t45_degR=t45,
        q_pt_ftlbf=q_pt,
        wa2_pps=cp.wa2,
        w45_pps=w45,
        far=far,
    )
    return nxt, out


def from_trim(result, lag_whole_flow: bool = True) -> RTState:
    """Seed a real-time state from a converged trim.

    The carried mass flow starts at its equilibrium value: the report never says how it
    initializes (open question #23), and equilibrium is what makes a trimmed engine sit
    still. `lag_whole_flow` must match `step`'s, because the two readings of Eq. 74 carry
    different quantities (WA31 against WA3_bl); seeding the wrong one puts a trimmed engine
    off equilibrium on frame one.

    The heat-sink carries are seeded from the trim's own frame: at equilibrium the metal
    sits at the gas temperature (Eq. 48) and T41 = T41_ns (Eq. 49).
    """
    s, f = result.state, result.frame
    carry = f.wa31_pps if lag_whole_flow else f.wa3_bl_pps
    t41 = f.t41_ns_degR
    return RTState(
        s.ng_rpm,
        s.np_rpm,
        s.p3_psia,
        s.p41_psia,
        s.p45_psia,
        carry,
        t41,  # metal temperature: at equilibrium it sits at the gas temperature (48)
        t41,
        f.w41_pps,
    )


def run(
    st: RTState,
    wf_of_t,
    ambient: Ambient = STANDARD_DAY,
    duration_s: float = 2.0,
    dt: float = FRAME_ENGINE_S,
    multirate: bool = True,
    **kw,
) -> dict[str, np.ndarray]:
    """Integrate for `duration_s`, returning traces keyed by name.

    `wf_of_t` is a callable of time in seconds returning fuel flow in lbm/sec.

    `multirate` runs the report's 2:1 configuration [pdf p.47]: the engine every frame, NP
    on alternate frames with `dt_np = 2*dt`. A no-op when `integrate_np=False`, as in
    Figures 9 and 10. Its cost is measured by
    `validation/test_transient.py::test_the_two_to_one_multirate_costs_little_but_is_the_reports`.

    `inner_iters`, `p45_iters` and `inner_capped` are diagnostics, not model outputs: they
    show a run that sat at the pass cap. `iteration_report` summarises them.
    """
    n = int(round(duration_s / dt)) + 1
    keys = (
        "t",
        "ng",
        "np",
        "p3",
        "p41",
        "p45",
        "t41",
        "t41_ns",
        "t45",
        "q_pt",
        "wf",
        "wa2",
        "far",
        "inner_iters",
        "p45_iters",
        "inner_capped",
        "wa31_substituted",
    )
    tr = {k: np.zeros(n) for k in keys}
    want_np = kw.get("integrate_np", True)
    if multirate:
        if kw.get("dt_np") is not None:
            raise TypeError(
                "run(multirate=True) sets dt_np to 2*dt itself; pass multirate=False to "
                "choose your own NP step."
            )
        kw["dt_np"] = 2.0 * dt
    for i in range(n):
        t = i * dt
        wf = float(wf_of_t(t))
        if multirate:
            # NP advances on even frames only, over the two engine frames it spans.
            # `want_np` is read once outside the loop: assigning into `kw` from
            # `kw.get(...)` would latch False after the first odd frame.
            kw["integrate_np"] = want_np and (i % 2 == 0)
        tr["t"][i] = t
        # NG and NP are integrator states, already at t_i. The pressures are algebraic,
        # solved inside the frame from the state at t_i and returned in the next state, so
        # they are recorded after the call; before it they would lag the speeds by a frame.
        tr["ng"][i] = st.ng_rpm
        tr["np"][i] = st.np_rpm
        tr["wf"][i] = wf
        st, out = step(st, wf, ambient, dt, **kw)
        tr["p3"][i] = st.p3_psia
        tr["p41"][i] = st.p41_psia
        tr["p45"][i] = st.p45_psia
        tr["t41"][i] = out.t41_degR
        tr["t41_ns"][i] = out.t41_ns_degR
        tr["t45"][i] = out.t45_degR
        tr["q_pt"][i] = out.q_pt_ftlbf
        tr["wa2"][i] = out.wa2_pps
        tr["far"][i] = out.far
        tr["inner_iters"][i] = out.inner_iters
        tr["p45_iters"][i] = out.p45_iters
        tr["inner_capped"][i] = out.inner_exit is Exit.CAPPED
        tr["wa31_substituted"][i] = out.wa31_substituted
    return tr


def iteration_report(tr: dict[str, np.ndarray]) -> dict[str, float]:
    """How hard the two pressure solves worked over a run, and how often they gave up.

    A trace of P3 cannot say whether P3 is converged; this is `maps.clamp_report()`'s idea
    applied to the iterations. Fractions are of all frames.
    """
    n = len(tr["t"])
    return {
        "frames": float(n),
        "inner_iters_mean": float(tr["inner_iters"].mean()),
        "inner_iters_max": float(tr["inner_iters"].max()),
        "inner_capped_frac": float(tr["inner_capped"].mean()),
        "p45_iters_mean": float(tr["p45_iters"].mean()),
        "p45_iters_max": float(tr["p45_iters"].max()),
        "wa31_substituted_frac": float(tr["wa31_substituted"].mean()),
    }
