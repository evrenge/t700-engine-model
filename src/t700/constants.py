"""T700-GE-700 model constants.

Every value here is transcribed from Table A.1 [TM-100991 pdf p.55], read from a page
image at 300-600 dpi. Not one of them came from the scan's OCR text layer, which renders
`K` as `g` and drops minus signs, nor from any handbook.

Printed precision is preserved exactly. `0.7826` is not rounded to `0.783`, and
`0.0018326` keeps all five significant figures. If a value looks over-specified, that is
Ballin's choice, not ours to tidy.

**The sub-subscript trap.** The report writes some constants with a subscript AND a
sub-subscript: `K_H3` with a further `1` or `2`. These are sixteen distinct constants,
not eight. Each pair is the slope and intercept of one linear fit and **the two members
have different units** -- `K_H41_1` is Btu/lbm/degR while `K_H41_2` is Btu/lbm. Flattening
a pair to a single name silently corrupts the fit. They are named `_1` and `_2` here.
"""

from __future__ import annotations

from typing import Final

# Source for every constant in this module.
SOURCE: Final = "TM-100991 pdf p.55, Table A.1"

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

Read `bl` as b-ell, for *bleed* -- NOT b-one. The OCR text layer renders this `gb1`,
which is how it entered an early draft of the notes as `Kb1`. Confirmed at 600 dpi
against the nomenclature [pdf p.10]."""

K_B3: Final = 0.0025
"""Station 3 bleed fraction, nondimensional. [row 5]

This component of the compressor discharge bleed never returns to the gas path
[body-engine-model.md, Eqs. 15-17, 44]."""

# --------------------------------------------------------------------------- gearbox

K_DAMP: Final = 0.06854
"""Helicopter gearbox damping coefficient, ft*lbf*sec/rad. [row 6]

Enters the power turbine torque balance [pdf p.25, Eq. 41] multiplied by
(NP - NP_DES), i.e. speed *displacement*, though the surrounding prose describes it as
based on change of power turbine speed. The equation governs."""

K_DPB: Final = 0.03045
"""Diffuser/burner pressure drop coefficient, lbf^2*sec^2/(lbm^2*in^4*degR). [row 7]"""

# --------------------------------------------------------------------------- enthalpy fits
# Linear enthalpy/temperature relations, per station. This is the whole of the model's
# "thermodynamics": there are no property tables and no polynomials, which is what the
# real-time budget bought. See t700.thermo -- nothing outside that package may use these.

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

H45 = K_H45 * H44 -- stations 4.4 and 4.5 differ only by cooling-bleed mixing."""

# --------------------------------------------------------------------------- temperature fits
# The inverse direction, with the report's own separately-printed coefficients.
# K_T41_* is the exact inverse of K_H41_* to the printed precision; see tests.

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

Numerically equal to K_TH41_1, and printed as its own row. Kept separate because they
are separate constants in the report; a future revision could differ."""

K_TH45_2: Final = 0.0856
"""Station 4.5 velocity ratio intercept, nondimensional. [row 26]"""

# --------------------------------------------------------------------------- volumes

K_V3: Final = 0.97
"""Station 3 control volume constant, lbf/(in^2*lbm*degR). [row 27]

The printed units read `lbf/in^2 * lbm * degR`; Eq. 42 forces them to be
lbf/(in^2*lbm*degR). One of the three constants that make this model real-time:
P3 = K_V3 * integral of T3*(net mass flow) dt. Pressure is a pure integrator, so
nothing iterates on it."""

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

Also the `NGDES` used by the fuel control's NG sensor gain 100/NGDES [Fig. C15], where
it is referenced but never defined -- that was open question #9."""

NP_DES: Final = 20900.0
"""Power turbine design speed, rpm. [row 32]

Consistent with the Appendix B trim condition NP = 20895 rpm [pdf p.67]."""

# --------------------------------------------------------------------------- heat sink

TC_T41: Final = 0.29
"""Station 4.1 heat-sink time constant coefficient. [row 33]

**The report contradicts itself on this one's units** (open question #4). Printed as
lbm^(4/5)*sec^(9/5)/degR^(1/2) in both Table A.1 and the nomenclature [pdf p.12], but
Eq. 51 [pdf p.26] requires sec^(1/5). The glyph is unambiguously a 9 at 600 dpi, so it
is Ballin's typesetting error, printed twice.

The VALUE is unaffected and is what the model uses. Do not rescale 0.29 to make the
printed units work."""
