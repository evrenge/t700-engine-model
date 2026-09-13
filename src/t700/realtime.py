"""The real-time frame: Eqs. 69-80, the model Ballin actually ran.

`engine.frame` evaluates the physics with the three volume pressures as **states**
(Eqs. 42-44). That is the 5-DOF model Appendix B linearizes, and it is what the trim
solver uses. **It cannot be integrated at the report's frame time.**

Our own Jacobian at the hover trim puts the fastest mode at 4881/sec -- a time constant of
0.205 ms -- so explicit integration of five states is stable only below about 0.41 ms. The
report runs the engine at **7 ms**, seventeen times over that limit. Not a preference: it
would diverge on the first frame.

So Ballin makes the three pressures **algebraic** and solves them each frame:

* the **opened compressor mass-flow iteration** (Eq. 74) breaks the expensive outer loop
  through the compressor map by carrying one mass flow across the frame boundary;
* the **inner P3/P41 loop** (Eqs. 76, 78) is a Gauss-Seidel sweep -- "eleven arithmetic
  operations for each pass ... no relaxation algorithm is required for convergence"
  [pdf p.37];
* **P45 is solved independently** (Eq. 80), "because it is a nonlinear function of its own
  value, iterative techniques must be used" [pdf p.37].

What remains integrated is NG and NP, whose modes run **437 / 470 / 738 ms** and
**358 / 470 / 602 ms** at the three Table B.1 trims -- comfortably stable at the 7 and
14 ms the report uses. (This paragraph said "NG (tau = 18 ms)" until 2026-09-12. The 18 ms
mode is a *pressure* mode, not NG: the 5-DOF eigenvalues at hover are 0.2, 0.4, 18.0, 358
and 437 ms, and the quasi-steady approximation removes the first three. Attributing a
pressure mode to the shaft understated NG's time constant by a factor of 24.)

This also decodes the Conclusions [pdf p.54]: omitting the high-speed inter-volume
mass-flow dynamics "was found to be unnecessary". Those dynamics *are* the two fast modes.
They are far faster than anything a pilot or rotor can feel, and carrying them would have
cost a twentyfold smaller time step.

## Decisions recorded here, not transcribed

**The integrator is never named in the report** (open question #26). Explicit Euler is used
here: it is what a 1988 real-time model with "no iteration between time steps" [pdf p.35]
would use, and the two surviving modes are slow enough for it. Logged as a decision.

**Eq. 74 reads two ways, and the report is not self-contradictory** (open question #22).
As printed it lags only the bleed, `WA31_(n) = WA3_(n) - WA3_bl_(n-1)`; the prose says the
combustor inlet flow is "the mass flow leaving the compressor in the previous interval",
i.e. `WA31_(n) = WA31_(n-1)`. Both are implemented and selectable, and the shipped default
is the prose reading -- because the prose is the sentence that describes *opening* the loop,
and because Eq. 74's own left-hand side is the quantity the combustor needs.

**It is chosen on those grounds and not because the printed equation fails.** This
paragraph said the printed reading "converges to nothing", and that was our bug: the
alternative branch carried `WA31` where Eq. 74 carries `WA3_bl`, giving
`WA31_(n) = WA3_(n) - WA31_(n-1)` -- an involution whose only fixed point is `WA3/2`, so it
oscillated with a two-frame period and a trimmed engine never sat still in it. With the
bleed carried, both readings converge first-order to the same limit (43847.1 against
43847.7 rpm at dt = 0.25 ms) and differ by **at most 0.11 % on NG over a whole 400 -> 775
step at the report's own 7 ms frame**, settling to 0.0002 %. The report's 0.1 %-per-frame
criterion therefore cannot separate them either, which is the real answer to #22:
the distinction is below the tolerance the report itself works to.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final

import numpy as np

from t700 import constants as c
from t700 import corrections, maps, thermo
from t700.engine import STANDARD_DAY, Ambient, Frame, State
from t700.units import BTU_TO_FTLBF, RAD_PER_SEC_TO_RPM, RPM_TO_RAD_PER_SEC

_TORQUE_SCALE = BTU_TO_FTLBF * RAD_PER_SEC_TO_RPM

FRAME_ENGINE_S = 0.007
"""Engine update interval [pdf p.47]: "updated twice for each rotor routine cycle, or once
every 7 msec"."""

FRAME_NP_S = 0.014
"""NP update interval [pdf p.47]. NP is the drive-train degree of freedom and moves at the
host rotor rate, so the model is multirate 2:1."""

MAX_STEP_S = 0.010
"""Maximum time step [pdf p.38], set by a 0.1 % maximum allowable error between steps.

Note this is *smaller* than FRAME_NP_S. The report states both; open question #33."""


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
    not stated anywhere in the report -- open question #23."""

    hs_tm_degR: float = 0.0
    """Station 4.1 **metal** temperature, deg R -- Eq. 48's state. Unused when
    `heat_sink=False`; zero means "not seeded" and `_heat_sink` then starts it at
    equilibrium. See `_heat_sink` for why the state is this and not Eq. 50's lead-lag
    memory."""

    t41_carry_degR: float = 0.0
    """Previous frame's T41, deg R. This is the report's **sixth state** [pdf p.29,
    p.32] -- gas temperature at station 4.1, not metal temperature. Eq. 51 needs it to
    form the time constant, and this frame has not computed it yet. Unused when
    `heat_sink=False`."""

    w41_carry_pps: float = 0.0
    """Previous frame's W41, lbm/sec. Eqs. 51 and 53 both need W41, which depends on
    theta41, which depends on T41 -- the quantity the heat sink is computing. Lagged one
    frame, exactly as Eq. 74 opens the mass-flow loop. Unused when `heat_sink=False`."""

    def to_state(self) -> State:
        return State(self.ng_rpm, self.np_rpm, self.p3_psia, self.p41_psia, self.p45_psia)


@dataclass(frozen=True)
class FrameOut:
    """What one real-time frame produced, beyond the next state."""

    inner_iters: int
    p45_iters: int
    inner_exit: Exit
    """How the P3/P41 sweep finished. Discarded by `run()` until 2026-09-13, so a run that
    spent every frame at the pass cap emitted no signal at all."""
    p45_exit: Exit
    """How the P45 iteration finished. `DIVERGING` below roughly 76 %NG is expected, not a
    fault -- see `_p45_loop`."""
    t41_degR: float
    t41_ns_degR: float
    """T41 before the heat sink. Equal to `t41_degR` when `heat_sink=False` [Eq. 23]."""
    t45_degR: float
    q_pt_ftlbf: float
    q_gt_ftlbf: float
    q_c_ftlbf: float
    wa2_pps: float
    w41_pps: float
    w45_pps: float
    far: float


# --------------------------------------------------------------------------- solver stopping rules

TOL_PRESSURE: Final = 1.0e-3
"""Relative **error** bound for both pressure iterations. [pdf p.37]

Printed twice on that page: the P3/P41 sweep converges "with less than 0.1 percent
error", and for P45 "eight iterations resulted in an error equal to less than 0.1 percent
of the steady-state value". 0.1 percent is 1e-3.

**Both sentences state an error, and neither states a tolerance.** The report prints no
convergence test at all -- what it prints is a cost (eleven arithmetic operations a pass,
four plus a table lookup for P45), a pass count, and the error that pass count achieved.
Until 2026-09-13 this constant was used as a tolerance on the iterate *step*, which is a
different quantity: for a linearly convergent iteration with contraction rho the remaining
error is rho/(1-rho) times the last step, and rho measures **0.887-0.902** here, so the
step test delivered roughly eight times the number it was named for. Measured through the
Figure 9 step, the true P3 error at exit ran to a median of 0.25 % and a worst frame of
4.21 %, under a constant documented as 0.1 %.

The loops below now stop on the error itself, estimated from the iterates as
`rho/(1-rho) * step` with rho taken from consecutive steps. That introduces no constant:
rho is measured, not assumed. See `_contraction_error`.
"""


class Exit(StrEnum):
    """How an iteration finished. Recorded because a model that stops early is a model
    running on unconverged pressures, and nothing downstream can tell from the value."""

    CONVERGED = "converged"
    """The estimated error fell below `tol` before the printed pass cap."""

    CAPPED = "capped"
    """The pass cap was reached first -- the report's "under the most extreme conditions,
    up to ten iterations may be required" [pdf p.37]. Expected during a transient and
    measured: on the Figure 9 accel 142 of 287 frames exit this way."""

    DIVERGING = "diverging"
    """Successive steps grew. For P45 this is a real operating regime and not a numerical
    accident -- see `_p45_loop` and open question #45."""


def _contraction_error(step: float, prev_step: float, value: float) -> tuple[float, float]:
    """A-posteriori error estimate for a linearly convergent fixed-point iteration.

    For `x_(k+1) = g(x_k)` with `|g'| = rho < 1` near the fixed point, the standard bound
    is `|x_k - x*| <= rho/(1-rho) * |x_k - x_(k-1)|`, with rho estimated by the ratio of
    consecutive steps. Returns `(estimated_error, rho)`.

    The report describes exactly this situation -- "the iteration converges linearly"
    [pdf p.37] -- so the estimate is the one the report's own sentence licenses. It needs
    two passes before it exists, which is why both loops below take a minimum of two.

    `rho >= 1` means the steps are growing and the estimate does not apply, so the bound
    returned is infinite: an iteration that is not contracting has no a-posteriori error
    estimate at all.

    With no previous step there is likewise no estimate, and the honest return is an
    infinite bound -- which is what forces the minimum of two passes. The exception is a
    step of exactly zero: the iterate has stopped moving, so the fixed point is reached
    and no contraction factor is needed to say so. (This is the branch a trimmed engine
    takes, and it is why holding a trim costs one pass rather than two.)

    "Stopped moving" means one ulp of `value`, not exactly zero. A step at the rounding
    floor carries no information, and a *ratio* of two such steps carries less: judged on
    exact zero, a held 400 lbm/hr trim ran 2.4 % of its frames to the ten-pass cap on
    rounding noise alone. One ulp is far below anything the converging iteration produces
    -- the P45 steps at that trim measure 69 ulp, then 5, then 0.
    """
    if step <= np.spacing(value):
        return 0.0, 0.0
    if prev_step <= 0.0:
        return float("inf"), 0.0
    rho = step / prev_step
    if rho >= 1.0:
        return float("inf"), rho
    return rho / (1.0 - rho) * step, rho


def _diverged(first_step: float, last_step: float) -> bool:
    """Did the iteration move further per pass at the end than at the start?

    Deliberately not "was any single ratio >= 1". Near a fixed point the steps reach the
    rounding floor -- measured at the 400 and 775 lbm/hr trims, the P45 steps run 69 ulp,
    5 ulp, then exactly 0 -- and ratios of quantities at that level are noise. Judging one
    pass on it reported divergence on 3.5 % of the frames of a *held trim*.

    Comparing the whole run of steps is immune to that and still unambiguous where it
    matters: at the 125 lbm/hr trim the P45 steps grow 57, 88, 137, ... 7912 ulp, a ratio
    of 1.57 a pass, which is f9's elasticity at that operating point to three figures.
    """
    return last_step > first_step


MAX_ITER_P3_P41: Final = 10
"""Pass cap for the P3/P41 sweep. [pdf p.37]  "up to ten iterations may be required",
cross-checked by the printed operation count: 11 operations a pass x 10 = 110."""

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
) -> tuple[float, float, int]:
    """The P3/P41 fixed-point sweep, Eqs. 76 and 78.

    Gauss-Seidel exactly as the report describes it [pdf p.37]: P3 from the current P41,
    then P41 from the P3 just computed -- not a simultaneous solve. T3 and theta_41 are
    held fixed across the sweep, which is what makes it eleven arithmetic operations a
    pass; they are recomputed once per frame outside this function.

    **The stopping rule is the report's, not ours.** Until 2026-09-12 this ran to
    `tol=1e-10` with a cap of 40 -- seven orders of magnitude tighter than the report,
    and invented. It is printed: *"up to ten iterations may be required for convergence
    with less than 0.1 percent error, resulting in a total of 110 arithmetic
    operations"* [pdf p.37]. Eleven operations a pass times ten passes is the 110, which
    pins the cap at ten independently of the sentence.

    It is not a detail. Converging harder than Ballin did makes the pressures reach
    equilibrium within one frame where his relaxed toward it across several, which
    sharpens every transient: on the Figure 9 step, tightening from the printed rule to
    1e-10 moves the T41 overshoot from +184.9 to +195.7 degR against Ballin's +116.3.
    Steady state is untouched -- a fixed point is a fixed point -- so no trim moves.
    See open question #25 for what is still not printed: how many passes were actually
    taken per frame, which is worth another 40 degR.

    **What the stopping test measures changed on 2026-09-13**, and the difference is a
    factor of eight. It tested the iterate *step*, `|P3_(k) - P3_(k-1)| < tol*P3`, under a
    constant documented as the report's 0.1 percent *error*. For a linearly convergent
    iteration the two differ by `rho/(1-rho)`, and rho measures 0.887-0.902 here, so the
    step test delivered a median true error of 0.25 % and a worst frame of 4.21 % on the
    Figure 9 accel. It now stops on `_contraction_error`, which is the report's own
    quantity.

    The change vindicates the report's arithmetic rather than departing from it. Under the
    error test this loop runs to a median of **nine passes** on Figure 9's step and takes
    the ten-pass cap on 142 frames of 287 -- which is "under the most extreme conditions,
    up to ten iterations may be required ... resulting in a total of 110 arithmetic
    operations" [pdf p.37], the printed budget, rather than the median of *one* pass the
    step test was exiting on. And the error it then delivers on the report's own test case
    -- the flight-idle-to-full-power step -- has a median of **0.0978 %**, against a
    printed "less than 0.1 percent". Nothing was fitted to that; it falls out of the
    measured contraction.
    """
    prev3 = prev41 = 0.0
    first3 = first41 = 0.0
    exit_ = Exit.CAPPED
    it = 0
    while it < max_iter:
        it += 1
        p3_new = 0.5 * (p41 + np.sqrt(p41 * p41 + 4.0 * c.K_DPB * t3 * wa31 * wa31))  # (76)
        w41 = c.K_WGT * p41 / np.sqrt(theta41)
        p41_new = p3_new - t3 * c.K_DPB * (w41 - wf) ** 2 / p3_new  # (78) with (77)
        step3, step41 = abs(p3_new - p3), abs(p41_new - p41)
        p3, p41 = p3_new, p41_new
        if it == 1:
            first3, first41 = step3, step41
        err3, _ = _contraction_error(step3, prev3, p3)
        err41, _ = _contraction_error(step41, prev41, p41)
        if err3 < tol * p3 and err41 < tol * p41:
            exit_ = Exit.CONVERGED
            break
        prev3, prev41 = step3, step41
    else:
        if _diverged(first3, prev3) or _diverged(first41, prev41):
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
) -> tuple[float, int]:
    """P45 by its own iteration, Eq. 80.

    Same correction as the P3/P41 sweep, and the report is equally explicit for this one:
    *"Eight iterations resulted in an error equal to less than 0.1 percent of the
    steady-state value"* [pdf p.37], measured under "an instantaneous step in fuel flow
    from flight idle to full power" -- which is Figure 9's own test case.

    The numerator is constant over the iteration, as the report notes -- only the f9
    lookup changes, which is why a pass costs "one function-table lookup and four
    arithmetic operations" [pdf p.37].

    Eq. 79 prints the returned bleed as `B3 B4 WA2` where Eq. 44 writes `B3 K_bl WA2` for
    the identical term, and B4 is in no nomenclature (open question #30). K_bl is used.

    **This iteration is not always contractive, and that is a property of Eq. 80 rather
    than of the arithmetic.** Differentiating `g(P45) = N / f9(Ps9/P45)` gives
    `g'(P45*) = dln(f9)/dln(Ps9/P45)` -- the fixed point's stability is exactly f9's
    elasticity at the operating point, and nothing else. f9's elasticity crosses -1 at
    `Ps9/P45 ~ 0.77`, so **above that pressure ratio the fixed point is repelling** and no
    tolerance reaches it. That is the mechanism behind open question #45; see
    `docs/notes/open-questions.md` for the elasticity table and for the explanation it
    replaced, which had it backwards.

    So the loop reports `Exit.DIVERGING` rather than returning whatever it happened to
    hold at the cap. It still returns that value -- there is nothing better to return, and
    the report's own model had no divergence test either -- but the condition is now
    visible in `FrameOut` and in `run()`'s traces instead of being discarded.
    """
    numerator = (w41 + b3 * c.K_BL * wa2) * np.sqrt(theta45)
    prev = first = 0.0
    exit_ = Exit.CAPPED
    it = 0
    while it < max_iter:
        it += 1
        p45_new = numerator / float(maps.f9()(ps9 / p45))
        step = abs(p45_new - p45)
        p45 = p45_new
        if it == 1:
            first = step
        err, _ = _contraction_error(step, prev, p45)
        if err < tol * p45:
            exit_ = Exit.CONVERGED
            break
        prev = step
    else:
        if _diverged(first, prev):
            exit_ = Exit.DIVERGING
    return float(p45), it, exit_


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

    ## The state is the METAL temperature, and that choice is load-bearing

    Eq. 50 is a *collapse* of the two printed heat-transfer equations [pdf p.25]:

    ```
    c_pm M dT_m/dt      = h A_m (T_gi - T_m)                        (48)
    c_pg W_g (T_gi-T_go) = c_pm M dT_m/dt                            (49)
    ```

    With `T_gi = T41_ns` and `T_go = T41`, those are exactly

    ```
    dT_m/dt = (T41_ns - T_m) / tau_a
    T41     = T41_ns - (tau_b/tau_a) * (T41_ns - T_m)
    ```

    which is what this function integrates. **The collapse into Eq. 50 is only valid for
    CONSTANT coefficients**, and Eqs. 51 and 53 make `tau_a` and `tau_b` functions of T41
    and W41. Until 2026-09-12 this function carried Eq. 50's lead-lag memory `x` instead,
    with `T41 = k*T41_ns + x`. Substituting the above shows `x = (tau_b/tau_a) * T_m`, so
    **`x` has `k` baked into it**: every time the coefficients move, the stored `x` refers
    to the old `k` and is stale. On a fuel step, where `k` swings from 0.62 to 0.67 inside
    a few frames, that stale memory amplified the T41 excursion beyond the `k` that Eq. 50
    prescribes -- our Figure 9 peak ran +3.64 % and our Figure 10 minimum -2.99 %, and the
    T41 error at matched NGc reached 160 deg R.

    Integrating `T_m` has no such artifact: the state is a physical temperature and means
    the same thing whatever the coefficients do. It is identical to the lead-lag whenever
    they are constant, so nothing that depends on the constant-coefficient case moves --
    Appendix B and Table 1 are untouched, and the DC gain stays exactly 1.

    The report's *linear* models carry "gas temperature at station 4.1" as the sixth state
    [pdf p.29] and that is still what Appendix B reproduces; the choice here is only about
    which realization the nonlinear frame integrates, and Eqs. 48-49 are the printed one.

    ## Two properties worth stating, because they are load-bearing

    * **DC gain is exactly 1.** At equilibrium T41 = T41_ns, so the heat sink cannot move
      a trim. That is why `engine.frame` and `trim.solve` need no heat-sink flag and why
      Table B.1 validates both configurations.
    * **High-frequency gain is `k = 1 - tau_b/tau_a`**, which is below 1 whenever
      `tau_b > 0`. The heat sink attenuates fast excursions in T41 and passes slow ones,
      which is the "response is significantly slowed" of [pdf p.33].
    """
    tau_a = c.TC_T41 * t41_prev_degR**0.5 / w41_prev_pps**0.8  # (51)
    tau_b = float(maps.f_hs()(ngc_pct)) / w41_prev_pps  # (52), (53)
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
) -> tuple[RTState, FrameOut]:
    """Advance one engine frame.

    Args:
        lag_whole_flow: how to read Eq. 74. True (default) follows the report's **prose** --
            the combustor inlet flow is the compressor exit flow of the previous interval,
            which is what actually opens the loop. False follows the equation **as
            printed**, lagging only the bleed. Open question #22; see the module docstring.
        integrate_np: False suppresses the NP integration, which is what the report does
            for the open-loop fuel steps of Figs. 9 and 10 [pdf p.39].
    """
    p2 = ambient.p_amb_psia
    t2 = ambient.t_amb_degR
    h2 = thermo.h2_from_t2(t2)  # (3)
    theta2 = corrections.theta2(t2)  # (5)
    ngc = corrections.corrected_speed(st.ng_rpm, theta2)  # (6)
    ngc_pct = 100.0 * ngc / c.NG_DES
    delta2 = corrections.delta2(p2)  # (8)

    # --- the compressor, evaluated ONCE on the frame's entering pressure ---------------
    ps3 = c.K_PS3 * st.p3_psia  # (4)
    pr = ps3 / p2
    wa2c = float(maps.f1()(pr, ngc_pct))  # (7)
    wa2 = corrections.physical_flow(wa2c, theta2, delta2)  # (9)
    t3 = t2 * float(maps.f2()(pr))  # (10)
    h3 = thermo.h3_from_t3(t3)  # (11)

    b1 = float(maps.f3()(ngc_pct))  # (12)
    b2 = float(maps.f4()(wa2c))  # (13)
    b3 = float(maps.f5()(wa2c))  # (14)
    wa3_bl = wa2 * (b3 + c.K_B3)  # (16)
    wa3 = wa2 - wa2 * (b1 + b2)  # (15), (17)

    # --- Eq. 74, the opened iteration -------------------------------------------------
    wa31 = st.wa31_carry_pps if lag_whole_flow else (wa3 - st.wa31_carry_pps)
    if wa31 <= 0.0:
        wa31 = max(wa3 - wa3_bl, 1e-6)

    # --- combustor, once per frame ----------------------------------------------------
    far = wf_pps / wa31  # (19)
    eta_b = float(maps.f6()(far))  # (20)
    h41_ns = (h3 + eta_b * far * c.HVF) / (1.0 + far)  # (21)
    t41_ns = thermo.t41_from_h41(h41_ns)  # (22)
    if heat_sink:
        t41, hs_tm = _heat_sink(
            t41_ns, st.hs_tm_degR, st.t41_carry_degR, st.w41_carry_pps, ngc_pct, dt
        )  # (48)-(53)
    else:
        t41, hs_tm = t41_ns, st.hs_tm_degR  # (23), no heat-sink representation
    h41 = thermo.h41_from_t41(t41)  # (24)
    theta41 = thermo.theta41_from_t41(t41)  # (25)

    # --- the two pressure solves ------------------------------------------------------
    p3, p41, inner_iters, inner_exit = _inner_pressure_loop(
        st.p3_psia, st.p41_psia, t3, theta41, wa31, wf_pps, tol=tol
    )
    w41 = c.K_WGT * p41 / np.sqrt(theta41)  # (28)

    ps9 = p2  # (37)
    p49 = ps9 * float(maps.f10()(ngc_pct))  # (38)

    dh_gt = theta41 * float(maps.f7()(st.p45_psia / p41))  # (26)
    h44 = h41 - dh_gt  # (27)
    h45 = thermo.h45_from_h44(h44)  # (29)
    t45 = thermo.t45_from_h45(h45)  # (30)
    theta45 = thermo.theta45_from_t45(t45)  # (31)

    p45, p45_iters, p45_exit = _p45_loop(st.p45_psia, w41, b3, wa2, theta45, ps9, tol=tol)

    dh_pt = theta45 * float(maps.f8()(p49 / p45))  # (32)
    w45 = float(maps.f9()(ps9 / p45)) * p45 / np.sqrt(theta45)  # (33), (34)

    # --- torques and the two surviving integrations -----------------------------------
    q_c = _TORQUE_SCALE / st.ng_rpm * (wa2 * (c.K_QC_1 * h3 - h2) + wa3 * c.K_QC_2 * h3)  # (39)
    q_gt = _TORQUE_SCALE / st.ng_rpm * w41 * dh_gt  # (40)
    q_pt = _TORQUE_SCALE / st.np_rpm * w45 * dh_pt - c.K_DAMP * RPM_TO_RAD_PER_SEC * (
        st.np_rpm - c.NP_DES
    )  # (41)

    ng = st.ng_rpm + dt * RAD_PER_SEC_TO_RPM * (q_gt - q_c) / c.J_GT  # (45)
    np_ = st.np_rpm
    if integrate_np:
        np_ = st.np_rpm + dt * RAD_PER_SEC_TO_RPM * (q_pt - q_req_ftlbf) / (
            c.J_PT + j_load
        )  # (46), (47)

    # The carried scalar is a DIFFERENT quantity in the two readings, and getting that
    # wrong is what invalidated open question #22's first verdict: this line read
    # `carry = wa31 if not lag_whole_flow else (wa3 - wa3_bl)`, so the printed branch
    # computed `WA31(n) = WA3(n) - WA31(n-1)` -- an involution whose only fixed point is
    # WA3/2, which oscillates instead of converging and made the printed reading look
    # divergent. Eq. 74 lags the **bleed**: `WA31(n) = WA3(n) - WA3_bl(n-1)`.
    carry = (wa3 - wa3_bl) if lag_whole_flow else wa3_bl
    nxt = RTState(ng, np_, p3, p41, p45, float(carry), hs_tm, t41, float(w41))
    out = FrameOut(
        inner_iters=inner_iters,
        p45_iters=p45_iters,
        inner_exit=inner_exit,
        p45_exit=p45_exit,
        t41_degR=t41,
        t41_ns_degR=t41_ns,
        t45_degR=t45,
        q_pt_ftlbf=q_pt,
        q_gt_ftlbf=q_gt,
        q_c_ftlbf=q_c,
        wa2_pps=wa2,
        w41_pps=w41,
        w45_pps=w45,
        far=far,
    )
    return nxt, out


def from_trim(
    result, wa31_pps: float, frame: Frame | None = None, lag_whole_flow: bool = True
) -> RTState:
    """Seed a real-time state from a converged trim.

    The carried mass flow is initialized to its equilibrium value. The report never says
    how it initializes (open question #23); starting at equilibrium is the choice that
    makes a trimmed engine sit still, which is the only defensible default.

    `lag_whole_flow` must match the value passed to `step`, because the two readings of
    Eq. 74 carry **different quantities**: the prose reading carries WA31, the printed
    reading carries WA3_bl. Seeding the wrong one puts a trimmed engine off equilibrium on
    frame one.

    Pass `frame` -- the `engine.frame` evaluated at the same trim -- to seed the
    heat-sink carries as well. Required for `heat_sink=True` runs; harmless otherwise.
    The seed follows the same principle: at equilibrium Eq. 48's right-hand side is zero,
    so the metal sits at the gas temperature and T41 = T41_ns by Eq. 49. One line, and it
    needs no coefficients at all -- which is itself a small argument for this realization
    over the lead-lag memory it replaced, whose seed had to be computed from `tau_a` and
    `tau_b`.
    """
    s = result.state
    carry = wa31_pps
    if not lag_whole_flow:
        if frame is None:
            raise ValueError(
                "the printed reading of Eq. 74 carries WA3_bl, so `frame` is required to "
                "seed it; pass engine.frame() evaluated at the same trim"
            )
        carry = frame.wa3_bl_pps
    if frame is None:
        return RTState(s.ng_rpm, s.np_rpm, s.p3_psia, s.p41_psia, s.p45_psia, carry)

    t41 = frame.t41_ns_degR
    return RTState(
        s.ng_rpm,
        s.np_rpm,
        s.p3_psia,
        s.p41_psia,
        s.p45_psia,
        carry,
        t41,  # metal temperature: at equilibrium it sits at the gas temperature (48)
        t41,
        frame.w41_pps,
    )


def run(
    st: RTState,
    wf_of_t,
    ambient: Ambient = STANDARD_DAY,
    duration_s: float = 2.0,
    dt: float = FRAME_ENGINE_S,
    **kw,
) -> dict[str, np.ndarray]:
    """Integrate for `duration_s`, returning traces keyed by name.

    `wf_of_t` is a callable of time in seconds returning fuel flow in lbm/sec, so a step
    input is just a lambda.

    Two of the traces are diagnostics rather than model outputs. `inner_iters` and
    `p45_iters` are the pass counts, and `inner_capped` / `p45_capped` / `p45_diverging`
    are the exit conditions as 0/1 flags -- `FrameOut` carried all of this and `run()`
    dropped it until 2026-09-13, so a run that spent every frame at the pass cap, or every
    frame on a repelling P45 fixed point, produced traces indistinguishable from a
    converged one. `iteration_report` summarises them.
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
        "p45_capped",
        "p45_diverging",
    )
    tr = {k: np.zeros(n) for k in keys}
    for i in range(n):
        t = i * dt
        wf = float(wf_of_t(t))
        tr["t"][i] = t
        # NG and NP are integrator states and are already at t_i. The three pressures are
        # NOT: they are algebraic, solved *inside* the frame from the state at t_i, and
        # returned in the next state. Recording them before the call logged the previous
        # frame's solution against this frame's speeds -- a 7 ms skew, invisible at steady
        # state, worth 12.1 psia on P3 at a fuel step, and it flattered the Figure 9
        # phase-plane agreement from a true 1.42 % to an apparent 0.59 %. No test caught
        # it; found by the 2026-09-13 physics audit.
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
        tr["p45_capped"][i] = out.p45_exit is Exit.CAPPED
        tr["p45_diverging"][i] = out.p45_exit is Exit.DIVERGING
    return tr


def iteration_report(tr: dict[str, np.ndarray]) -> dict[str, float]:
    """How hard the two pressure solves worked over a run, and how often they gave up.

    The one thing a trace of P3 cannot tell you is whether P3 is converged. `maps
    .clamp_report()` already did this job for the other failure mode -- a lookup outside
    its data -- and this is the same idea for the iterations. Fractions are of all frames.
    """
    n = len(tr["t"])
    return {
        "frames": float(n),
        "inner_iters_mean": float(tr["inner_iters"].mean()),
        "inner_iters_max": float(tr["inner_iters"].max()),
        "inner_capped_frac": float(tr["inner_capped"].mean()),
        "p45_iters_mean": float(tr["p45_iters"].mean()),
        "p45_iters_max": float(tr["p45_iters"].max()),
        "p45_capped_frac": float(tr["p45_capped"].mean()),
        "p45_diverging_frac": float(tr["p45_diverging"].mean()),
    }


__all__ = [
    "FRAME_ENGINE_S",
    "MAX_STEP_S",
    "Exit",
    "FrameOut",
    "RTState",
    "from_trim",
    "iteration_report",
    "run",
    "step",
]
