"""Standard-day corrections: theta, delta, and corrected speed and flow.

Turbomachinery is described in *corrected* quantities so that one map serves every
ambient condition. Ballin follows the convention: station 2 total conditions are divided
by standard sea-level values to give theta_2 and delta_2 [pdf p.22, Eqs. 5 and 8], and
speeds and flows are corrected by them before entering the compressor map.

## Where T_std and P_std come from

**They are not printed in the report.** `T_std` and `P_std` appear in the nomenclature
[pdf pp.11, 13] and in Eqs. 5 and 8, but no value is given for either, anywhere -- that
was open question #15, and it blocked every corrected quantity.

The values below are a **decision, not a transcription**, and the distinction matters
enough to spell out:

* Table B.1 [pdf p.67] gives the sea-level trim ambient as P2 = 14.696 PSIA and
  T2 = 518.67 deg R.
* theta_2 and delta_2 must both be exactly 1.0 at standard sea level, by definition of a
  corrected quantity.
* Therefore T_std = 518.67 deg R and P_std = 14.696 psia.

These are also the standard-atmosphere sea-level values in the US customary system the
whole report uses, which is the independent confirmation. Approved as a project decision
on 2026-09-10; recorded against open question #15.

If a later reading turns up printed values that differ, these are wrong and everything
corrected by them moves. That is why they live here, named, rather than inline.
"""

from __future__ import annotations

from typing import Final

# --------------------------------------------------------------------------- decided, not cited

T_STD_DEGR: Final = 518.67
"""Standard sea-level temperature, deg R. **Decision, not a transcription** -- see the
module docstring and open question #15. Inferred from Table B.1 [pdf p.67]."""

P_STD_PSIA: Final = 14.696
"""Standard sea-level pressure, psia. **Decision, not a transcription** -- see the module
docstring and open question #15. Inferred from Table B.1 [pdf p.67]."""


# --------------------------------------------------------------------------- the ratios


def theta2(t2_degR):
    """Station 2 temperature ratio, T2/T_std, nondimensional. [pdf p.22, Eq. 5]

    Exactly 1.0 at standard sea level, which is the property that pinned T_std.
    """
    return t2_degR / T_STD_DEGR


def delta2(p2_psia):
    """Station 2 pressure ratio, P2/P_std, nondimensional. [pdf p.22, Eq. 8]

    Exactly 1.0 at standard sea level.
    """
    return p2_psia / P_STD_PSIA


# --------------------------------------------------------------------------- corrected quantities


def corrected_speed(n_rpm, theta):
    """Corrected shaft speed, N/sqrt(theta), rpm.

    Used to enter the compressor map: `f1` takes corrected gas generator speed
    [pdf p.22]. Note the report writes `NG_c`, and `docs/notes/symbols.md` keeps the
    subscript convention.
    """
    return n_rpm / theta**0.5


def corrected_flow(w_pps, theta, delta):
    """Corrected mass flow, W*sqrt(theta)/delta, lbm/sec.

    The compressor map is drawn in corrected flow, so the physical flow out of `f1` has
    to be uncorrected again before it enters any mass balance. Getting that backwards
    produces an engine that is subtly wrong only away from sea level, which is the worst
    place for an error to hide.
    """
    return w_pps * theta**0.5 / delta


def physical_flow(wc_pps, theta, delta):
    """Physical mass flow from corrected, lbm/sec. Inverse of `corrected_flow`."""
    return wc_pps * delta / theta**0.5


def ambient_is_standard(t2_degR, p2_psia, tol: float = 1e-9) -> bool:
    """True when both corrections are unity, i.e. standard sea level.

    Exists so a test can assert the Table B.1 sea-level trim really does give
    theta_2 = delta_2 = 1, which is the evidence T_std and P_std rest on.
    """
    return abs(theta2(t2_degR) - 1.0) < tol and abs(delta2(p2_psia) - 1.0) < tol
