"""The primitive blocks Appendix C draws, and nothing else.

Appendix C specifies the fuel control entirely as pictures. The same handful of icons
recurs across all 22 block diagrams -- a lag, a lead-lag, a limited integrator, a deadband,
a backlash, a saturation, a transport delay -- so they are written once here and the ECU
and HMU modules read as the figures do.

Every one is a **fixed-step explicit** step, matching the report's real-time formulation
and CLAUDE.md's architecture constraint. None of them holds state: the caller owns it, so
a state's name says which figure owns it.
"""

from __future__ import annotations


def clamp(v: float, lo: float, hi: float) -> float:
    """A saturation block. Appendix C draws these with their limits labelled."""
    return lo if v < lo else (hi if v > hi else v)


def lag(x: float, u: float, tau: float, dt: float) -> float:
    """One step of `1/(tau*s + 1)`. Returns the new state, which is also the output.

    `tau <= 0` passes the input straight through, which is what a lag of zero means and
    keeps a variable time constant (Fig. C5's TAU45) from dividing by zero.
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

    Clamping the state rather than the output is what makes it anti-windup: an integrator
    that ran free behind a clamped output would take an unbounded time to come back.
    Appendix C draws the limits inside the integrator icon, which is the same thing.
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

    **The report does not say whether the named constant is the whole band or the half
    band** -- open question #54. Figure C15's icon braces `CH` from the vertical axis to
    one branch, which argues for a half-width, but the brace measures 43 px against a
    branch separation of 58 px, so the drawing does not settle it. This takes the constant
    as the **total** band, the ordinary meaning of a hysteresis block's width parameter.
    Reading it as a half-width doubles every band; that matters most for `XLDHYS`, 2.5 deg
    on a spindle graduated in 10 deg steps.
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

    [Fig. C19, pdf p.92] The fuel line's 0.015 s against the report's 7 ms engine frame is
    about two frames -- not a rounding detail, and not resolved by the frame either, so it
    has to be interpolated between samples rather than rounded to one.
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
