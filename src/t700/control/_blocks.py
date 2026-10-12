"""The primitive blocks Appendix C draws, and nothing else.

The same few icons recur across all 22 block diagrams, so they are written once here and
the ECU and HMU modules read as the figures do. Each is a fixed-step explicit step. None
holds state: the caller owns it, so a state's name says which figure owns it.
"""

from __future__ import annotations


def clamp(v: float, lo: float, hi: float) -> float:
    """A saturation block. Appendix C draws these with their limits labelled."""
    return lo if v < lo else (hi if v > hi else v)


def lag(x: float, u: float, tau: float, dt: float) -> float:
    """One step of `1/(tau*s + 1)`. Returns the new state, which is also the output.

    `tau <= 0` passes the input straight through, so a variable time constant (Fig. C5's
    TAU45) cannot divide by zero.

    Explicit Euler puts the pole at `z = 1 - dt/tau`. Six lags (`CT9`, `CTPL`, `CTPS3`,
    `T17`, `TL1`, `TL2`) have tau = 0.010 s [Table C.1, pdf pp.83-84]: dt/tau = 0.7 at the
    7 ms frame, deadbeat at 10 ms (`realtime.MAX_STEP_S`), a negative pole at 14 ms
    (`FRAME_NP_S`); see open question #33. Kept as the model's one scheme;
    `tests/test_control_constants.py` pins it so a frame change cannot slip past.
    """
    if tau <= 0.0:
        return u
    return x + dt * (u - x) / tau


def lead_lag(x: float, u: float, lead: float, lag_t: float, dt: float) -> tuple[float, float]:
    """One step of `(lead*s + 1)/(lag*s + 1)`. Returns (new state, output).

    The state is the lag's, and the output reads the lead off the same derivative rather
    than differentiating the input -- the realization the station 4.1 heat sink uses, for
    the same reason: differentiating a sampled input amplifies its steps.
    """
    if lag_t <= 0.0:
        return u, u
    dx = (u - x) / lag_t
    return x + dt * dx, x + lead * dx


def integrator(x: float, u: float, gain: float, lo: float, hi: float, dt: float) -> float:
    """One step of a **limited** integrator, `clamp(∫gain*u dt, lo, hi)`.

    Clamping the state rather than the output makes it anti-windup. Appendix C draws
    something else: an unlimited `1/s` followed by a separate saturation [Figs. C3, C7,
    C10, pdf pp.86, 88, 89], which winds up. The two differ only while a limit binds, and
    nothing the report prints tells them apart; anti-windup is kept by decision (open
    question #59).
    """
    return clamp(x + dt * gain * u, lo, hi)


def deadband(u: float, width: float) -> float:
    """Zero inside `+/-width`, and continuous at the breakpoints outside it.

    [Figs. C3, C10] Both of Appendix C's deadbands are drawn with the outer segments
    offset to meet the axis at the breakpoints, not parallel to the input -- so the output
    is `u -/+ width`, not `u`.
    """
    if u > width:
        return u - width
    if u < -width:
        return u + width
    return 0.0


def backlash(state: float, u: float, width: float) -> float:
    """Hysteresis of total band `width`: the output follows the input, but only after the
    input has reversed by `width`.

    The report does not say whether the named constant is the whole band or the half band
    (open question #54; Figure C15's icon does not settle it); this takes it as the total
    band, the ordinary meaning. Read as a half-width, every band doubles.
    """
    half = 0.5 * width
    if u > state + half:
        return u - half
    if u < state - half:
        return u + half
    return state


def delay(
    history: tuple[float, ...], u: float, tau: float, dt: float
) -> tuple[tuple[float, ...], float]:
    """Pure transport delay `e^(-tau*s)` on a fixed step, by linear interpolation.

    [Fig. C19, pdf p.92] The fuel line's 0.015 s is about two 7 ms frames, so it is
    interpolated between samples rather than rounded to one.
    """
    if dt <= 0.0:
        return history, u
    n = max(int(tau / dt) + 2, 2)
    buf = (history + (u,))[-n:]
    if len(buf) < 2:
        return buf, u
    i = len(buf) - 1 - tau / dt
    if i <= 0:
        return buf, buf[0]
    lo_i = int(i)
    frac = i - lo_i
    hi_i = min(lo_i + 1, len(buf) - 1)
    return buf, buf[lo_i] + frac * (buf[hi_i] - buf[lo_i])
