"""Standard-day corrections: theta, delta, and corrected speed and flow.

Station 2 total conditions over standard sea-level values give theta_2 and delta_2
[pdf p.22, Eqs. 5 and 8], which correct speeds and flows for the compressor map.

`T_std` and `P_std` are in the nomenclature [pdf pp.11, 13] but no value is printed
(open question #15). The values below are a decision, not a transcription: Table B.1
[pdf p.67] gives the sea-level trim ambient as P2 = 14.696 PSIA and T2 = 518.67 deg R,
and theta_2 and delta_2 must be exactly 1.0 there. They are also the US standard
atmosphere's sea-level values.
"""

from __future__ import annotations

from typing import Final

# --------------------------------------------------------------------------- decided, not cited

T_STD_DEGR: Final = 518.67
"""Standard sea-level temperature, deg R. A decision inferred from Table B.1 [pdf p.67],
not a transcription (open question #15)."""

P_STD_PSIA: Final = 14.696
"""Standard sea-level pressure, psia. A decision inferred from Table B.1 [pdf p.67], not a
transcription (open question #15)."""


# --------------------------------------------------------------------------- the ratios


def theta2(t2_degR):
    """Station 2 temperature ratio, T2/T_std, nondimensional. [pdf p.22, Eq. 5]"""
    return t2_degR / T_STD_DEGR


def delta2(p2_psia):
    """Station 2 pressure ratio, P2/P_std, nondimensional. [pdf p.22, Eq. 8]"""
    return p2_psia / P_STD_PSIA


# --------------------------------------------------------------------------- corrected quantities


def corrected_speed(n_rpm, theta):
    """Corrected shaft speed, N/sqrt(theta), rpm; `f1` takes `NG_c` [pdf p.22]."""
    return n_rpm / theta**0.5


def physical_flow(wc_pps, theta, delta):
    """Physical mass flow from corrected, Wc*delta/sqrt(theta), lbm/sec.

    The compressor map is in corrected flow, W*sqrt(theta)/delta, so `f1`'s output must be
    uncorrected before any mass balance. Getting it backwards is wrong only off sea level.
    """
    return wc_pps * delta / theta**0.5
