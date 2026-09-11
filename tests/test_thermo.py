"""Check the thermodynamic fits against the report's own internal consistency.

The strongest available check on a transcription is the report disagreeing with itself
if you got it wrong. Ballin prints the station 4.1 enthalpy fit and its inverse with
*separately printed* coefficients. If all four numbers were read correctly, the two fits
must compose to the identity -- and four independently transcribed numbers agreeing to
four significant figures is not something a misread digit survives.
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import constants as c
from t700 import thermo


def test_station_41_fits_are_mutual_inverses():
    """K_T41_1 == 1/K_H41_1 and K_T41_2 == -K_H41_2/K_H41_1, to printed precision.

    [Table A.1 rows 11, 12, 17, 18 -- pdf p.55]
    """
    assert c.K_T41_1 == pytest.approx(1.0 / c.K_H41_1, rel=2e-4)
    assert c.K_T41_2 == pytest.approx(-c.K_H41_2 / c.K_H41_1, rel=2e-4)


def test_station_41_round_trip():
    """T41 -> H41 -> T41 returns the original within the fits' printed precision.

    0.2 deg R on a ~2100 deg R turbine inlet is 0.01%, which is the rounding in the
    printed coefficients rather than a modelling error.
    """
    t41 = np.array([1500.0, 1800.0, 2100.0, 2400.0])
    back = thermo.t41_from_h41(thermo.h41_from_t41(t41))
    assert np.allclose(back, t41, atol=0.2)


def test_enthalpy_increases_with_temperature():
    """Every enthalpy fit has positive slope. A sign error here inverts the cycle."""
    t = np.array([500.0, 1000.0, 2000.0, 3000.0])
    for fn in (thermo.h2_from_t2, thermo.h3_from_t3, thermo.h41_from_t41):
        h = fn(t)
        assert np.all(np.diff(h) > 0), f"{fn.__name__} is not increasing in T"


def test_velocity_ratio_fits_agree_between_stations():
    """theta_41 and theta_45 use numerically equal coefficients [Table A.1 rows 23-26].

    They are separate constants in the report and separate constants here. This test
    records that they currently coincide, so that if one is ever corrected the other is
    considered deliberately rather than by accident.
    """
    t = np.array([1200.0, 1600.0, 2000.0])
    assert np.array_equal(thermo.theta41_from_t41(t), thermo.theta45_from_t45(t))


def test_station_45_is_a_fraction_of_station_44():
    """H45 = K_H45 * H44 with K_H45 < 1: mixing cooling bleed lowers enthalpy."""
    assert 0.0 < c.K_H45 < 1.0
    h44 = 400.0
    assert thermo.h45_from_h44(h44) < h44


def test_backend_is_swappable_and_restores():
    """The interface exists so a compiled backend can replace it. Prove it works."""

    class Sentinel:
        name = "sentinel"

        def h2_from_t2(self, t2):
            return -1.0

    original = thermo.backend()
    try:
        thermo.set_backend(Sentinel())
        assert thermo.backend().name == "sentinel"
        assert thermo.h2_from_t2(500.0) == -1.0
    finally:
        thermo.set_backend(original)
    assert thermo.backend().name == "linear-fit"
