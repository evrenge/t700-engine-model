"""Units and conversion factors.

Internal units are **the report's units** -- US customary: lbm, lbm/sec, psia, deg R,
hp, ft*lbf, rpm, in^2. Not SI. The reason (CLAUDE.md) is that every printed table and
figure can then be compared to model output literally, with no conversion sitting in the
loop waiting to be wrong.

So this module is small on purpose. It holds the one conversion the report itself uses
inside its equations, plus the boundary conversions for user-facing input and output.
Nothing else in the codebase computes a conversion inline.

Two kinds of constant live here and they are not the same kind of thing:

* **From the report.** Cited to a page. `BTU_TO_FTLBF` is Ballin's own 778.12, printed
  as a bare number three times. Use his value, not a modern handbook's 778.169, because
  the point is to reproduce his arithmetic.
* **Definitional.** Exact by definition of the unit (60 seconds in a minute). These carry
  no citation because there is nothing to cite; they are marked `definitional`.
"""

from __future__ import annotations

import math

# --------------------------------------------------------------------------- from the report

# [TM-100991 pdf p.25, Eqs. 39, 40, 41] mechanical equivalent of heat, ft*lbf per Btu.
# Printed as a bare number in all three torque equations; it is NOT a Table A.1 row and
# NOT a nomenclature entry (open question #16). The modern value is 778.169 -- do not
# substitute it. Reproducing Ballin means reproducing Ballin's constant.
BTU_TO_FTLBF = 778.12

# --------------------------------------------------------------------------- definitional

SEC_PER_MIN = 60.0  # definitional
RAD_PER_REV = 2.0 * math.pi  # definitional

# rpm -> rad/sec. The report's torque equations carry this as the explicit factor
# 60/(2*pi) applied to a shaft speed in rpm [pdf p.25, Eqs. 39-41].
RPM_TO_RAD_PER_SEC = RAD_PER_REV / SEC_PER_MIN
RAD_PER_SEC_TO_RPM = SEC_PER_MIN / RAD_PER_REV

LBM_PER_HR_TO_LBM_PER_SEC = 1.0 / 3600.0  # definitional


def wf_pps_from_pph(wf_pph: float) -> float:
    """Fuel flow lbm/hr to lbm/sec.

    The report quotes fuel flow both ways and they are easy to confuse: Table B.1 trims
    are lbm/hr (476.3 at hover), while `Wf` in the equations and as the linear-model
    control input is lbm/sec [nomenclature pdf p.13].
    """
    return wf_pph * LBM_PER_HR_TO_LBM_PER_SEC


def shp_from_torque(q_ftlbf: float, n_rpm: float) -> float:
    """Shaft horsepower from torque (ft*lbf) and speed (rpm).

    Uses the definitional 33000 ft*lbf/min per horsepower, NOT the report's 778.12 --
    these are different constants for different jobs and conflating them is a classic
    unit error in this kind of model.

    This reproduces Table B.1's engine torque row exactly: 911.1 hp at 20895 rpm gives
    229.0 ft*lbf, matching the printed 229.0 [pdf p.67].
    """
    return q_ftlbf * n_rpm * RPM_TO_RAD_PER_SEC / 550.0
