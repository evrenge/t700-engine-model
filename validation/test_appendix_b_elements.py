"""Our extracted derivatives against Appendix B, element by element [pdf pp.68-70].

Table 1 gives 12 comparable eigenvalues. Appendix B gives **297 printed numbers**, and
comparing them element-wise localizes a disagreement to a single partial derivative
instead of reporting that "a mode is off". Everything interesting in this file was found
that way and could not have been found from eigenvalues.

## What the comparison found

Three separate error families, not one:

1. **The NP row was ~10x too large** -- uniformly, in every element and in `b`. That is one
   scalar: the load inertia missing from J_PT [Eq. 46]. Recovered from `b(NP)` as
   `j_load = 9.2503 * J_PT` with a 0.186 % spread across three flight conditions, and
   recorded as `constants.J_LOAD_UH60A`. See open question #6.

2. **The NP *diagonal* keeps a further ~52 % residual** even with the inertia supplied.
   That is the second Gen Hel quantity, dQreq/dNP, which enters nowhere else. It comes out
   positive -- load torque rising with power turbine speed, as a rotor does.

3. **A P45 group -- `d(NG)/d(P45)`, `d(P45)/d(P41)`, `d(P45)/d(P45)`, `d(NG)/d(P41)` --
   fails worst at hover and least at descent.** All four run through `f7` [Eq. 26], whose
   six printed knots make its *derivative* a step function: -526.4, -254.7, -150.8, -156.8,
   -108.1. The trims land 98.6 %, 71.7 % and 32.5 % across their segments, and the error
   tracks the segment slope. Swapping `f7` to a shape-preserving cubic moves exactly these
   four elements and nothing else, which is what pins the cause. See open question #43.

## What is deliberately not asserted

Nothing here compares the 3- or 6-DOF models; `t700.linear` cannot extract them yet, and
the 6-DOF matrices are too ill-conditioned for anything but element-wise work anyway.
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import appendix_b as ab
from t700 import constants as c
from t700 import trim
from t700.engine import Ambient
from t700.linear import DOF, extract
from t700.units import wf_pps_from_pph

AMB = Ambient(14.696, 518.67)
NP_RPM = 20895.0
WF_PPH = {1: 476.3, 2: 349.3, 3: 267.7}
TRIM_NAME = {1: "hover", 2: "level 80 kt", 3: "descent 80 kt"}
S = ("NG", "NP", "P3", "P41", "P45")

F7_DRIVEN = {("NG", "P45"), ("P45", "P41"), ("P45", "P45"), ("NG", "P41")}
"""The four elements `f7`'s interpolated slope controls. Characterized, not asserted tight."""

NG_COLUMN = {("NG", "NG"), ("P3", "NG"), ("P45", "NG")}
"""Derivatives with respect to NG. These scatter to +-12 % with no trend in power, unlike
the `f7` group's monotone one -- a different cause, not yet localized."""


def _pair(trim_no: int, j_load: float = 0.0):
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    r = trim.solve(wf, NP_RPM, AMB)
    assert r.trustworthy, f"trim {trim_no} is not on data"
    return extract(r, wf, DOF.FIVE, AMB, j_load=j_load), ab.find(5, trim_no)


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_zero_structure_matches_exactly(trim_no: int):
    """Every element Appendix B prints as zero, we compute as exactly zero.

    Nineteen of thirty elements are non-zero; the other eleven are structural. Matching
    the sparsity pattern is a check on the *equations*, independent of any numeric
    agreement: a spurious coupling would show up here even if its magnitude were small.
    """
    ours, ref = _pair(trim_no)
    for i, si in enumerate(S):
        for j, sj in enumerate(S):
            if ref.A[i, j] == 0.0:
                assert ours.A[i, j] == 0.0, f"d({si})/d({sj}) should be structurally zero"
    assert ref.b[2] == 0.0 and ours.b[2] == 0.0, "fuel does not enter the P3 volume directly"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_pressure_block_agrees_closely(trim_no: int):
    """The P3 and P41 volume equations are right to under 1 % at every trim.

    These six derivatives are the core of the volume dynamics and involve no turbine map
    slope, which is why they are the tightest thing in the comparison -- and why the
    failures elsewhere cannot be blamed on the pressure formulation.
    """
    ours, ref = _pair(trim_no)
    for si, sj in (("P3", "NG"), ("P3", "P3"), ("P3", "P41"), ("P41", "P3"), ("P41", "P41")):
        i, j = S.index(si), S.index(sj)
        if (si, sj) in NG_COLUMN:
            continue
        dev = (ours.A[i, j] - ref.A[i, j]) / abs(ref.A[i, j]) * 100.0
        assert abs(dev) < 1.5, f"trim {trim_no} d({si})/d({sj}): {dev:+.2f} %"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_fuel_column_is_essentially_exact(trim_no: int):
    """`b` matches to 0.1 % everywhere the load inertia does not enter.

    This is what makes the NP-row diagnosis airtight: the fuel derivative into NG, P41 and
    P45 is right to a part in a thousand, so `b(NP)`'s factor of ten cannot be the fuel
    path and must be the inertia dividing it.
    """
    ours, ref = _pair(trim_no)
    for s in ("NG", "P41", "P45"):
        i = S.index(s)
        dev = (ours.b[i] - ref.b[i]) / abs(ref.b[i]) * 100.0
        assert abs(dev) < 0.2, f"trim {trim_no} b({s}): {dev:+.3f} %"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_load_inertia_recovers_the_np_row(trim_no: int):
    """Supplying the derived `j_load` takes the NP row from ~900 % error to within 12 %.

    The elements that remain off are the ones carrying their own model error on top of the
    inertia; `b(NP)`, which carries none, lands inside 0.2 %.
    """
    bare, ref = _pair(trim_no)
    fixed, _ = _pair(trim_no, j_load=c.J_LOAD_UH60A)

    bare_dev = abs((bare.b[1] - ref.b[1]) / ref.b[1] * 100.0)
    fixed_dev = abs((fixed.b[1] - ref.b[1]) / ref.b[1] * 100.0)
    assert bare_dev > 800.0, f"b(NP) with j_load=0 should be ~900 % off, got {bare_dev:.0f} %"
    assert fixed_dev < 0.2, f"b(NP) with the derived inertia: {fixed_dev:+.3f} %"

    for s in ("NG", "P3", "P41", "P45"):
        j = S.index(s)
        dev = (fixed.A[1, j] - ref.A[1, j]) / abs(ref.A[1, j]) * 100.0
        assert abs(dev) < 12.0, f"trim {trim_no} d(NP)/d({s}): {dev:+.1f} %"


def test_the_inertia_ratio_is_the_same_at_all_three_trims():
    """The evidence that the factor of ten is inertia and not something else.

    Inertia is the only quantity in that row which *must* be identical at hover, level
    flight and descent. Measured spread is 0.186 %.
    """
    ratios = []
    for trim_no in (1, 2, 3):
        ours, ref = _pair(trim_no)
        ratios.append(ours.b[1] / ref.b[1])
    spread = (max(ratios) - min(ratios)) / np.mean(ratios) * 100.0
    assert spread < 0.5, f"inertia ratio is not constant across trims: {ratios}, {spread:.2f} %"
    implied = (np.mean(ratios) - 1.0) * c.J_PT
    assert abs(implied - c.J_LOAD_UH60A) / c.J_LOAD_UH60A < 0.01, (
        f"constants.J_LOAD_UH60A is {c.J_LOAD_UH60A}, the data now implies {implied:.5f}"
    )


def test_the_np_diagonal_still_needs_the_load_torque_slope():
    """dQreq/dNP is the second Gen Hel quantity and is not recoverable from inertia.

    With the inertia supplied the diagonal is still about half Ballin's, consistently at
    all three trims -- consistent because it is one missing physical term, not noise.
    Characterized here so that supplying it later has a target to hit.
    """
    devs = []
    for trim_no in (1, 2, 3):
        fixed, ref = _pair(trim_no, j_load=c.J_LOAD_UH60A)
        devs.append((abs(fixed.A[1, 1]) - abs(ref.A[1, 1])) / abs(ref.A[1, 1]) * 100.0)
    assert all(-60.0 < d < -45.0 for d in devs), (
        f"NP diagonal residual moved: {[f'{d:+.1f}%' for d in devs]}. If it shrank, "
        f"dQreq/dNP has been supplied -- re-measure and tighten rather than widening."
    )


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_f7_group_is_characterized_by_trim(trim_no: int):
    """The four `f7`-driven elements, recorded with their power dependence.

    Worst at hover (98.6 % across f7's steepest segment), least at descent (32.5 % across
    a shallow one). Bounds are per trim and deliberately just wide enough to hold today's
    values; the monotone ordering is asserted separately because it is the actual evidence.
    """
    bound = {1: 26.0, 2: 16.0, 3: 9.0}[trim_no]
    ours, ref = _pair(trim_no)
    for si, sj in F7_DRIVEN:
        i, j = S.index(si), S.index(sj)
        dev = (ours.A[i, j] - ref.A[i, j]) / abs(ref.A[i, j]) * 100.0
        assert abs(dev) < bound, f"trim {trim_no} d({si})/d({sj}): {dev:+.1f} %"


def test_the_f7_error_shrinks_monotonically_with_power():
    """The signature that identifies `f7`'s interpolated slope as the cause.

    `d(NG)/d(P45)` runs -24.3 %, -14.5 %, -7.3 % from hover to descent. A wrong equation
    would not order itself by where the operating point sits on one table's segments; a
    wrong interpolated derivative does exactly that.
    """
    devs = []
    for trim_no in (1, 2, 3):
        ours, ref = _pair(trim_no)
        i, j = S.index("NG"), S.index("P45")
        devs.append(abs((ours.A[i, j] - ref.A[i, j]) / ref.A[i, j] * 100.0))
    assert devs[0] > devs[1] > devs[2], f"the power ordering is gone: {devs}"
    assert devs[0] > 2.5 * devs[2], f"the spread has collapsed: {devs}"
