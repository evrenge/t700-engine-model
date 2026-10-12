"""Units and conversion factors.

Internal units are the report's, US customary (lbm, lbm/sec, psia, deg R, hp, ft*lbf, rpm,
in^2), so every printed table compares with model output with no conversion in the loop.
Nothing else in the codebase computes a conversion inline.

Each constant is either cited to the report or marked `definitional` (exact by definition
of the unit, nothing to cite).
"""

from __future__ import annotations

import math

# --------------------------------------------------------------------------- from the report

# [TM-100991 pdf p.25, Eqs. 39, 40, 41] mechanical equivalent of heat, ft*lbf per Btu.
# A bare number in the torque equations, not a Table A.1 row (open question #16). Do not
# substitute the modern 778.169: reproducing Ballin means his constant.
BTU_TO_FTLBF = 778.12

# --------------------------------------------------------------------------- definitional

SEC_PER_MIN = 60.0  # definitional
RAD_PER_REV = 2.0 * math.pi  # definitional

# rpm -> rad/sec: the explicit factor 60/(2*pi) of the torque equations [pdf p.25, Eqs. 39-41].
RPM_TO_RAD_PER_SEC = RAD_PER_REV / SEC_PER_MIN
RAD_PER_SEC_TO_RPM = SEC_PER_MIN / RAD_PER_REV

LBM_PER_HR_TO_LBM_PER_SEC = 1.0 / 3600.0  # definitional


def wf_pps_from_pph(wf_pph: float) -> float:
    """Fuel flow lbm/hr to lbm/sec.

    Easy to confuse: Table B.1 trims are lbm/hr, while `Wf` in the equations and as the
    linear-model input is lbm/sec [nomenclature pdf p.13].
    """
    return wf_pph * LBM_PER_HR_TO_LBM_PER_SEC


def shp_from_torque(q_ftlbf: float, n_rpm: float) -> float:
    """Shaft horsepower from torque (ft*lbf) and speed (rpm).

    Uses the definitional 550 ft*lbf/sec per horsepower, not the report's 778.12 Btu
    factor. Reproduces Table B.1: 911.1 hp at 20895 rpm is 229.0 ft*lbf [pdf p.67].
    """
    return q_ftlbf * n_rpm * RPM_TO_RAD_PER_SEC / 550.0
