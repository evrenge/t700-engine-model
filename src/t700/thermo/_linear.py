"""Linear-fit thermodynamics backend. Pure NumPy.

Ballin's model has no gas property tables and no polynomials: enthalpy, temperature and
the squared critical velocity ratios relate through two-coefficient straight lines per
station, coefficients from Table A.1 (`t700.constants`). Reproducing the report means
reproducing the fits, not substituting better thermodynamics.
"""

from __future__ import annotations

from t700 import constants as c

Num = float  # or any NumPy array; every operation below is elementwise


class LinearFitThermo:
    """The report's thermodynamics, exactly as printed. Stateless."""

    name = "linear-fit"

    # ----------------------------------------------------------------- enthalpy from T

    def h2_from_t2(self, t2_degR: Num) -> Num:
        """Station 2 enthalpy, Btu/lbm. [pdf p.22, Eq. 3]  H2 = K_H2 * T2, no intercept."""
        return c.K_H2 * t2_degR

    def h3_from_t3(self, t3_degR: Num) -> Num:
        """Station 3 enthalpy, Btu/lbm. [pdf p.23]  H3 = K_H3_1*T3 + K_H3_2"""
        return c.K_H3_1 * t3_degR + c.K_H3_2

    def h41_from_t41(self, t41_degR: Num) -> Num:
        """Station 4.1 enthalpy, Btu/lbm. [pdf p.24]  H41 = K_H41_1*T41 + K_H41_2"""
        return c.K_H41_1 * t41_degR + c.K_H41_2

    def h45_from_h44(self, h44: Num) -> Num:
        """Station 4.5 enthalpy from station 4.4, Btu/lbm. [pdf p.25]  H45 = K_H45 * H44

        The report models the cooling-bleed mixing as one multiplicative fraction.
        """
        return c.K_H45 * h44

    # ----------------------------------------------------------------- T from enthalpy

    def t41_from_h41(self, h41: Num) -> Num:
        """Station 4.1 temperature, deg R. [pdf p.24]  T41 = K_T41_1*H41 + K_T41_2

        The report prints this inverse with its own coefficients; `tests/test_thermo.py`
        checks they invert the forward fit to the printed precision.
        """
        return c.K_T41_1 * h41 + c.K_T41_2

    def t45_from_h45(self, h45: Num) -> Num:
        """Station 4.5 temperature, deg R. [pdf p.25]  T45 = K_T45_1*H45 + K_T45_2

        T45 is the power turbine inlet temperature; the report has no station 4 or 5.
        """
        return c.K_T45_1 * h45 + c.K_T45_2

    def t49_from_h49(self, h49: Num) -> Num:
        """Station 4.9 temperature, deg R. [pdf p.25]  T49 = K_T49_1*H49 + K_T49_2"""
        return c.K_T49_1 * h49 + c.K_T49_2

    # ----------------------------------------------------------------- velocity ratios

    def theta41_from_t41(self, t41_degR: Num) -> Num:
        """Station 4.1 squared critical velocity ratio. [pdf p.23, Eq. 25]

        theta_41 = K_TH41_1*T41 + K_TH41_2. Feeds the turbine flow relation
        W41 = K_WGT * P41 / sqrt(theta_41) [Eq. 28], which is why Eq. 23 closes an
        algebraic loop in T41 (open question #17).
        """
        return c.K_TH41_1 * t41_degR + c.K_TH41_2

    def theta45_from_t45(self, t45_degR: Num) -> Num:
        """Station 4.5 squared critical velocity ratio. [pdf p.25]

        theta_45 = K_TH45_1*T45 + K_TH45_2; equal to the station 4.1 pair but separate
        Table A.1 rows.
        """
        return c.K_TH45_1 * t45_degR + c.K_TH45_2
