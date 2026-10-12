"""T700-GE-700 model constants, Table A.1 [TM-100991 pdf p.55].

Read from page images, never from the OCR text layer. Printed precision is kept exactly
(`0.7826`, `0.0018326`); do not tidy it.

Sub-subscript trap: `K_H41_1` and `K_H41_2` are the slope (Btu/lbm/degR) and intercept
(Btu/lbm) of one linear fit, with different units. There is no plain `K_H41`; flattening a
pair to one name silently corrupts the fit.
"""

from __future__ import annotations

from typing import Final

# --------------------------------------------------------------------------- fuel

HVF: Final = 18300.0
"""Fuel lower heating value, Btu/lbm. [Table A.1 row 1]"""

# --------------------------------------------------------------------------- inertias

J_GT: Final = 0.0445
"""Gas generator rotor polar moment of inertia, ft*lbf*sec^2. [row 2]"""

J_PT: Final = 0.062
"""Power turbine rotor polar moment of inertia, ft*lbf*sec^2. [row 3]"""

# --------------------------------------------------------------------------- bleeds

K_BL: Final = 0.7826
"""Fraction of diffuser bleed gas returned to the gas path, nondimensional. [row 4]

`bl` is b-ell, for bleed, not b-one (the OCR layer reads `gb1`). [nomenclature, pdf p.10]"""

K_B3: Final = 0.0025
"""Station 3 bleed fraction, nondimensional. [row 5]

Never returns to the gas path [body-engine-model.md, Eqs. 15-17, 44]."""

# --------------------------------------------------------------------------- gearbox

K_DAMP: Final = 0.06854
"""Helicopter gearbox damping coefficient, ft*lbf*sec/rad. [row 6]

Enters the power turbine torque balance [pdf p.25, Eq. 41] multiplied by
(NP - NP_DES), i.e. speed *displacement*, though the surrounding prose describes it as
based on change of power turbine speed. The equation governs."""

K_DPB: Final = 0.03045
"""Diffuser/burner pressure drop coefficient, lbf^2*sec^2/(lbm^2*in^4*degR). [row 7]"""

# --------------------------------------------------------------------------- enthalpy fits
# Linear enthalpy/temperature fits per station: the model's whole thermodynamics.
# Nothing outside t700.thermo may use these.

K_H2: Final = 0.239
"""Station 2 enthalpy coefficient, Btu/(lbm*degR). [row 8]  H2 = K_H2 * T2 [Eq. 3]"""

K_H3_1: Final = 0.2496
"""Station 3 enthalpy slope, Btu/(lbm*degR). [row 9]  H3 = K_H3_1*T3 + K_H3_2"""

K_H3_2: Final = -8.4
"""Station 3 enthalpy intercept, Btu/lbm. NEGATIVE. [row 10]"""

K_H41_1: Final = 0.3010
"""Station 4.1 enthalpy slope, Btu/(lbm*degR). [row 11]  H41 = K_H41_1*T41 + K_H41_2"""

K_H41_2: Final = -86.905
"""Station 4.1 enthalpy intercept, Btu/lbm. NEGATIVE. [row 12]"""

K_H45: Final = 0.9623
"""Fraction of station 4.4 enthalpy representing station 4.5, nondimensional. [row 13]

H45 = K_H45 * H44."""

# --------------------------------------------------------------------------- temperature fits
# The inverse direction, with the report's own separately printed coefficients.

K_PS3: Final = 0.956
"""Ratio of station 3 static to total pressure, nondimensional. [row 14]"""

K_QC_1: Final = 0.71
"""Compressor torque coefficient 1, nondimensional. [row 15]"""

K_QC_2: Final = 0.29
"""Compressor torque coefficient 2, nondimensional. [row 16]"""

K_T41_1: Final = 3.322
"""Station 4.1 temperature slope, lbm*degR/Btu. [row 17]  T41 = K_T41_1*H41 + K_T41_2"""

K_T41_2: Final = 288.7
"""Station 4.1 temperature intercept, degR. [row 18]"""

K_T45_1: Final = 3.519
"""Station 4.5 temperature slope, lbm*degR/Btu. [row 19]"""

K_T45_2: Final = 179.1
"""Station 4.5 temperature intercept, degR. [row 20]"""

K_T49_1: Final = 3.516
"""Station 4.9 temperature slope, lbm*degR/Btu. [row 21]"""

K_T49_2: Final = 172.3
"""Station 4.9 temperature intercept, degR. [row 22]"""

# --------------------------------------------------------------------------- velocity ratios

K_TH41_1: Final = 0.0018326
"""Station 4.1 squared-critical-velocity-ratio slope, 1/degR. [row 23]

theta_41 = K_TH41_1*T41 + K_TH41_2 [pdf p.23, Eq. 25]"""

K_TH41_2: Final = 0.0856
"""Station 4.1 velocity ratio intercept, nondimensional. [row 24]"""

K_TH45_1: Final = 0.0018326
"""Station 4.5 squared-critical-velocity-ratio slope, 1/degR. [row 25]

Equal to K_TH41_1 but its own row in the report, so kept separate."""

K_TH45_2: Final = 0.0856
"""Station 4.5 velocity ratio intercept, nondimensional. [row 26]"""

# --------------------------------------------------------------------------- volumes

K_V3: Final = 0.97
"""Station 3 control volume constant, lbf/(in^2*lbm*degR). [row 27]

Printed as `lbf/in^2 * lbm * degR`; Eq. 42 forces lbf/(in^2*lbm*degR)."""

K_V41: Final = 6.17
"""Station 4.1 control volume constant, lbf/(in^2*lbm*degR). [row 28]"""

K_V45: Final = 13.63
"""Station 4.5 control volume constant, lbf/(in^2*lbm*degR). [row 29]"""

K_WGT: Final = 0.0876
"""Gas generator turbine flow coefficient, lbm*in^2/(lbf*sec). [row 30]

W41 = K_WGT * P41 / sqrt(theta_41) [pdf p.23, Eq. 28]"""

# --------------------------------------------------------------------------- design speeds

NG_DES: Final = 44700.0
"""Gas generator design speed, rpm. [row 31]

Also the `NGDES` of the NG sensor gain 100/NGDES [Fig. C15], undefined there (open
question #9)."""

NP_DES: Final = 20900.0
"""Power turbine design speed, rpm. [row 32]

Consistent with the Appendix B trim condition NP = 20895 rpm [pdf p.67]."""

# --------------------------------------------------------------------------- derived, not printed

J_LOAD_UH60A: Final = 0.57352
"""Load inertia added to J_PT [Eq. 46], for the UH-60A drivetrain of the Appendix B trims.

DERIVED, not printed: Table A.1 gives no load inertia (open question #6). With
`j_load = 0` our Appendix B `b(NP)` is Ballin's times (J_PT + j_load) / J_PT, and the
hover trim gives that ratio as 10.2503; the other two trims agree, as an inertia must.

    J_PT + j_load = 10.2503 * J_PT = 0.63552
    j_load        =  9.2503 * J_PT = 0.57352

A property of the airframe, not the engine: opt-in for `t700.linear.extract` and
`engine.frame`, whose default stays 0.
"""

# --------------------------------------------------------------------------- heat sink

TC_T41: Final = 0.29
"""Station 4.1 heat-sink time constant coefficient. [row 33]

Units printed lbm^(4/5)*sec^(9/5)/degR^(1/2) in Table A.1 and the nomenclature [pdf p.12],
but Eq. 51 [pdf p.26] requires sec^(1/5): a typesetting error (open question #4). The value
is unaffected; do not rescale 0.29 to fit the printed units."""
