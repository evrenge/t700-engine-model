"""Our extracted derivatives against Appendix B, element by element [pdf pp.68-70].

Appendix B prints 297 numbers; comparing them element-wise localizes a disagreement to a
single partial derivative, which a mode cannot. Three error families, each found that way:

1. The NP row is ~10x too large without the load inertia missing from J_PT [Eq. 46],
   recovered from `b(NP)` as `constants.J_LOAD_UH60A` (open question #6).
2. The NP diagonal keeps a ~52 % residual even with it: dQreq/dNP, the second Gen Hel
   quantity, which enters nowhere else.
3. The four elements through `f7` [Eq. 26] track where each trim sits on `f7`'s piecewise
   linear segments, whose derivative is a step function (open question #43).

Most bounds here are characterizations (SCOPE.md's third category); see SCOPE.md,
"Appendix B, element by element".
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import appendix_b as ab
from t700 import constants as c
from t700 import maps, trim
from t700.engine import Ambient
from t700.linear import DOF, extract
from t700.units import wf_pps_from_pph
from test_all_dof_models import DQ_REQ_DNP, WF_PPH
from test_all_dof_models import TRIM as TRIM_NAME

AMB = Ambient(14.696, 518.67)
NP_RPM = 20895.0
S = ("NG", "NP", "P3", "P41", "P45")

TANGENT_STEP = 1e-5
"""For the mechanism tests: knot proximity in `f1` and `f7` (#43) is a statement about local
slopes. The shipped extraction uses the report's +-2 % (pdf p.28, #56), which averages
across those knots and removes most of what these tests describe."""

F7_DRIVEN = {("NG", "P45"), ("P45", "P41"), ("P45", "P45"), ("NG", "P41")}
"""The four elements `f7`'s interpolated slope controls. Characterized, not asserted tight."""

NG_COLUMN = {("NG", "NG"), ("P3", "NG"), ("P45", "NG")}
"""Derivatives with respect to NG: they carry the `f1` speed-line error (#43), not `f7`'s, and
show no monotone trend in power."""


DOF_NO = {DOF.TWO: 2, DOF.THREE: 3, DOF.FIVE: 5, DOF.SIX: 6}


def _pair(
    dof: DOF,
    trim_no: int,
    j_load: float = 0.0,
    rel_step: float | None = None,
    dq_req_dnp: float = 0.0,
):
    """Our `dof` model at Table B.1 trim `trim_no`, and the Appendix B figure it answers to.

    `j_load` defaults to the bare power turbine; the heat-sink models are compared with
    `c.J_LOAD_UH60A` supplied, because Appendix B's matrices carry the UH-60A drivetrain.
    """
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    r = trim.solve(wf, NP_RPM, AMB)
    assert r.trustworthy, f"trim {trim_no} is not on data"
    kw = {} if rel_step is None else {"rel_step": rel_step}
    ours = extract(r, wf, dof, AMB, j_load=j_load, dq_req_dnp=dq_req_dnp, **kw)
    return ours, ab.find(DOF_NO[dof], trim_no)


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_zero_structure_matches_exactly(trim_no: int):
    """Every element Appendix B prints as zero, we compute as exactly zero: a check on the
    equations, independent of numeric agreement, that catches a small spurious coupling.
    """
    ours, ref = _pair(DOF.FIVE, trim_no)
    for i, si in enumerate(S):
        for j, sj in enumerate(S):
            if ref.A[i, j] == 0.0:
                assert ours.A[i, j] == 0.0, f"d({si})/d({sj}) should be structurally zero"
    assert ref.b[2] == 0.0 and ours.b[2] == 0.0, "fuel does not enter the P3 volume directly"


PRESSURE_BLOCK_TOL_PCT = 3.0
"""The 5-DOF P3/P41 block, per element, at the report's +-2 % step [pdf p.28]. Eq. 18's
square root makes the secant visibly non-tangent here, and B2/B4/B6 sit with the tangent
(#56)."""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_pressure_block_agrees_closely(trim_no: int):
    """The P3 and P41 volume equations involve no turbine map slope, so they are among the
    tightest elements, and failures elsewhere cannot be blamed on the pressure formulation.
    """
    ours, ref = _pair(DOF.FIVE, trim_no)
    for si, sj in (("P3", "NG"), ("P3", "P3"), ("P3", "P41"), ("P41", "P3"), ("P41", "P41")):
        i, j = S.index(si), S.index(sj)
        if (si, sj) in NG_COLUMN:
            continue
        dev = (ours.A[i, j] - ref.A[i, j]) / abs(ref.A[i, j]) * 100.0
        assert abs(dev) < PRESSURE_BLOCK_TOL_PCT, f"trim {trim_no} d({si})/d({sj}): {dev:+.2f} %"


FUEL_COLUMN_TOL_PCT = 0.2
"""The 5-DOF `b` into NG, P41 and P45."""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_fuel_column_is_essentially_exact(trim_no: int):
    """`b` is essentially exact wherever the load inertia does not enter, so `b(NP)`'s factor
    of ten cannot be the fuel path and must be the inertia dividing it.
    """
    ours, ref = _pair(DOF.FIVE, trim_no)
    for s in ("NG", "P41", "P45"):
        i = S.index(s)
        dev = (ours.b[i] - ref.b[i]) / abs(ref.b[i]) * 100.0
        assert abs(dev) < FUEL_COLUMN_TOL_PCT, f"trim {trim_no} b({s}): {dev:+.3f} %"


NP_ROW_TOL_PCT = 17.0
"""The 5-DOF NP row with the derived inertia, at the report's +-2 % step."""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_load_inertia_recovers_the_np_row(trim_no: int):
    """Supplying the derived `j_load` takes the NP row from ~900 % error to within the band;
    `b(NP)`, which carries no model error of its own, lands within a fraction of a percent.
    """
    bare, ref = _pair(DOF.FIVE, trim_no)
    fixed, _ = _pair(DOF.FIVE, trim_no, c.J_LOAD_UH60A)

    bare_dev = abs((bare.b[1] - ref.b[1]) / ref.b[1] * 100.0)
    fixed_dev = abs((fixed.b[1] - ref.b[1]) / ref.b[1] * 100.0)
    assert bare_dev > 800.0, f"b(NP) with j_load=0 should be ~900 % off, got {bare_dev:.0f} %"
    assert fixed_dev < 0.25, f"b(NP) with the derived inertia: {fixed_dev:+.3f} %"

    for s in ("NG", "P3", "P41", "P45"):
        j = S.index(s)
        dev = (fixed.A[1, j] - ref.A[1, j]) / abs(ref.A[1, j]) * 100.0
        assert abs(dev) < NP_ROW_TOL_PCT, f"trim {trim_no} d(NP)/d({s}): {dev:+.1f} %"


def test_the_inertia_ratio_is_the_same_at_all_three_trims():
    """The factor of ten is inertia: the only quantity in that row which must be identical
    at hover, level flight and descent.
    """
    ratios = []
    for trim_no in (1, 2, 3):
        ours, ref = _pair(DOF.FIVE, trim_no)
        ratios.append(ours.b[1] / ref.b[1])
    spread = (max(ratios) - min(ratios)) / np.mean(ratios) * 100.0
    assert spread < 0.5, f"inertia ratio is not constant across trims: {ratios}, {spread:.2f} %"
    implied = (np.mean(ratios) - 1.0) * c.J_PT
    assert abs(implied - c.J_LOAD_UH60A) / c.J_LOAD_UH60A < 0.01, (
        f"constants.J_LOAD_UH60A is {c.J_LOAD_UH60A}, the data now implies {implied:.5f}"
    )


F7_GROUP_TOL_PCT = {1: 8.0, 2: 12.0, 3: 5.0}
"""Per trim, the four `f7`-driven elements of the 5-DOF `A`, at the report's +-2 % step."""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_f7_group_is_characterized_by_trim(trim_no: int):
    """The four `f7`-driven elements, per trim, with bounds just wide enough to hold today's
    values. The power ordering, asserted separately on the tangent, is the evidence.
    """
    bound = F7_GROUP_TOL_PCT[trim_no]
    ours, ref = _pair(DOF.FIVE, trim_no)
    for si, sj in F7_DRIVEN:
        i, j = S.index(si), S.index(sj)
        dev = (ours.A[i, j] - ref.A[i, j]) / abs(ref.A[i, j]) * 100.0
        assert abs(dev) < bound, f"trim {trim_no} d({si})/d({sj}): {dev:+.1f} %"


def test_the_f7_error_shrinks_monotonically_with_power():
    """The signature that identifies `f7`'s interpolated slope as the cause.

    On the tangent, `d(NG)/d(P45)`'s error falls from hover to descent: a wrong equation
    would not order itself by where the trim sits on one table's segments. At the report's
    +-2 % step the ordering goes, because the secant averages across the knot.
    """
    devs = []
    for trim_no in (1, 2, 3):
        ours, ref = _pair(DOF.FIVE, trim_no, rel_step=TANGENT_STEP)
        i, j = S.index("NG"), S.index("P45")
        devs.append(abs((ours.A[i, j] - ref.A[i, j]) / ref.A[i, j] * 100.0))
    assert devs[0] > devs[1] > devs[2], f"the power ordering is gone: {devs}"
    assert devs[0] > 2.5 * devs[2], f"the spread has collapsed: {devs}"


# --------------------------------------------------------------- knot proximity


def _knot_distance_f7(trim_no: int) -> float:
    """How far this trim's P45/P41 sits from the nearest `f7` knot, as a fraction of its
    segment. 0 means on a knot, where a linear interpolant's derivative is undefined."""
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    r = trim.solve(wf, NP_RPM, AMB)
    arg = r.state.p45_psia / r.state.p41_psia
    x = maps.f7().x
    k = int(np.searchsorted(x, arg))
    frac = (arg - x[k - 1]) / (x[k] - x[k - 1])
    return min(frac, 1.0 - frac)


def _knot_distance_f1(trim_no: int) -> float:
    """Same, for NGc against `f1`'s speed lines."""
    from t700.engine import frame

    wf = wf_pps_from_pph(WF_PPH[trim_no])
    r = trim.solve(wf, NP_RPM, AMB)
    ngc = frame(r.state, wf, AMB).ngc_pct
    p = maps.f1().params
    k = int(np.searchsorted(p, ngc))
    frac = (ngc - p[k - 1]) / (p[k] - p[k - 1])
    return min(frac, 1.0 - frac)


def _err(trim_no: int, si: str, sj: str, rel_step: float | None = None) -> float:
    ours, ref = _pair(DOF.FIVE, trim_no, rel_step=rel_step)
    i, j = S.index(si), S.index(sj)
    return abs((ours.A[i, j] - ref.A[i, j]) / ref.A[i, j] * 100.0)


def test_hover_is_the_double_case_that_pins_the_mechanism():
    """One operating point, two tables: hover sits on an `f7` knot (the worst `f7`-driven
    errors, `test_the_f7_error_shrinks_monotonically_with_power`) and mid-segment on `f1`,
    where its NG column stays small. See `docs/notes/derivative-ambiguity.md`, #43.
    """
    assert _knot_distance_f7(1) < 0.05, "hover no longer sits on an f7 knot"
    assert _knot_distance_f1(1) > 0.30, "hover no longer sits mid-segment on f1"

    for si, sj in NG_COLUMN:
        assert _err(1, si, sj, TANGENT_STEP) < 5.0, (
            f"hover d({si})/d({sj}) is {_err(1, si, sj, TANGENT_STEP):.1f} %"
        )


def test_the_ng_column_error_is_redistributed_not_removed():
    """`f1` is piecewise linear between speed lines, so d/dNGc jumps at each one.

    Which trim carries the NG-column error depends on the interpolation scheme, which the
    report does not state (open question #43, `docs/notes/derivative-ambiguity.md`). This
    pins the current ordering and that the spread stays bounded.
    """
    d = {t: _knot_distance_f1(t) for t in (1, 2, 3)}
    assert d[1] > d[3] > d[2], f"knot distances changed: {d}"
    for si, sj in NG_COLUMN:
        e = {t: _err(t, si, sj, TANGENT_STEP) for t in (1, 2, 3)}
        # Level and descent trade places between elements; only "hover is best" holds on
        # every element, so only that is asserted.
        assert e[1] < e[2] and e[1] < e[3], (
            f"d({si})/d({sj}): hover {e[1]:.1f} % should be the best, against level "
            f"{e[2]:.1f} % and descent {e[3]:.1f} %; look at the f1 interpolation."
        )
        assert max(e.values()) < 25.0, f"d({si})/d({sj}) spread widened: {e}"


# ------------------------------------------------- the heat-sink models, B7-B12


SIX_DOF_T41_COLUMN_TOL_PCT = 1.0
"""The 6-DOF T41 column, per element."""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_six_dof_t41_column_is_close(trim_no: int):
    """How the five engine equations respond to T41, the 6-DOF's new physics: a check on
    Eq. 23's coupling of T41 into h41 and theta41.
    """
    ours, ref = _pair(DOF.SIX, trim_no, c.J_LOAD_UH60A)
    t41 = ours.states.index("T41")
    for i, s in enumerate(ours.states):
        if s == "T41" or ref.A[i, t41] == 0.0:
            continue
        dev = (ours.A[i, t41] - ref.A[i, t41]) / abs(ref.A[i, t41]) * 100.0
        assert abs(dev) < SIX_DOF_T41_COLUMN_TOL_PCT, f"trim {trim_no} d({s})/d(T41): {dev:+.2f} %"


def test_the_six_dof_loses_the_p3_to_np_coupling_the_five_dof_has():
    """The 5-DOF prints a non-zero `A(NP,P3)`; the 6-DOF prints zero at all three trims. With
    T41 a state, P3 reaches power-turbine torque only through T41, so ours must be zero too.
    """
    assert ab.find(5, 1).A[1, 2] != 0.0
    for trim_no in (1, 2, 3):
        ours, ref = _pair(DOF.SIX, trim_no, c.J_LOAD_UH60A)
        assert ref.A[1, 2] == 0.0, "B8/B10/B12 should print zero for d(NP)/d(P3)"
        assert ours.A[1, 2] == 0.0, f"trim {trim_no}: we still couple P3 to NP directly"


T41_ROW = (0.90, 0.96)
T41_ROW_SPREAD_PCT = 1.0
"""The 6-DOF T41 row less its NG element: ratio band, and max-min spread over the mean."""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_t41_row_is_uniformly_low_by_the_lead_lag_ratio(trim_no: int):
    """Characterization, on the tangent: the T41 row sits uniformly low, for one reason.

    The row scales with `tau1/tau2`, and uniformity across three trims says the cause is the
    time constants, not the structure (open question #31, `test_heat_sink_taus.py`). The NG
    element is excluded: it is built through the NG column and carries the `f1` error too.
    Not at the report's +-2 % step, where Eq. 18's secant offsets the deficit. If this
    moves, `f_hs` or `TC_T41` has: re-measure and tighten.
    """

    ours, ref = _pair(DOF.SIX, trim_no, c.J_LOAD_UH60A, TANGENT_STEP)
    t41 = ours.states.index("T41")
    ng = ours.states.index("NG")
    ratios = [
        ours.A[t41, j] / ref.A[t41, j]
        for j in range(len(ours.states))
        if ref.A[t41, j] != 0.0 and j != ng
    ]
    assert all(T41_ROW[0] < r < T41_ROW[1] for r in ratios), f"T41 row ratios moved: {ratios}"
    spread = (max(ratios) - min(ratios)) / np.mean(ratios) * 100.0
    assert spread < T41_ROW_SPREAD_PCT, (
        f"the row is no longer uniformly scaled: spread {spread:.2f} %"
    )


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_t41_row_ng_element_carries_the_f1_error_on_top(trim_no: int):
    """`A(T41, NG)` is the one element of the T41 row built through the NG column, so it
    carries the `f1` error on top of the lead-lag ratio: bounded loosely, on its own.
    """
    ours, ref = _pair(DOF.SIX, trim_no, c.J_LOAD_UH60A)
    t41, ng = ours.states.index("T41"), ours.states.index("NG")
    r_ng = ours.A[t41, ng] / ref.A[t41, ng]
    assert 0.75 < r_ng < 1.15, f"trim {trim_no} A(T41,NG) ratio {r_ng:.3f}"


# --------------------------------------------------- the load-torque slope, and its limits
# `DQ_REQ_DNP` is `test_all_dof_models.DQ_REQ_DNP`, recovered from the elements it then
# corrects; the test below says so rather than claiming a win.


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_supplying_the_load_slope_closes_the_np_diagonal_circularly(trim_no: int):
    """With `dQreq/dNP` supplied the NP diagonal matches: circular, since the slope came
    from this element, so only a plumbing guard. The non-circular half: the slope moves the
    NP diagonal and nothing else, because Q_req enters only dNP/dt and depends only on NP.
    """
    bare, ref = _pair(DOF.FIVE, trim_no, c.J_LOAD_UH60A)
    fixed, _ = _pair(DOF.FIVE, trim_no, c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[trim_no])

    # Without the slope it sits ~52 % low, which NP_DIAGONAL_RESIDUAL_PCT holds in
    # `test_the_heat_sink_np_diagonals_carry_the_same_missing_load_slope`.
    after = (abs(fixed.A[1, 1]) - abs(ref.A[1, 1])) / abs(ref.A[1, 1]) * 100.0
    assert abs(after) < 1.0, f"NP diagonal is {after:+.1f} % with it"

    # The structural half: nothing else may move.
    for i in range(5):
        for j in range(5):
            if (i, j) == (1, 1):
                continue
            assert fixed.A[i, j] == pytest.approx(bare.A[i, j], rel=1e-9), (
                f"d({S[i]})/d({S[j]}) moved; Q_req enters only dNP/dt and depends only on NP"
            )
    assert np.allclose(fixed.b, bare.b, rtol=1e-9)


# ======================================== the 3- and 6-DOF blocks, and the rest of the 5-DOF
#
# Characterizations in SCOPE.md's third sense: bounds fitted round the current measurement,
# so a regression fails and each block's shape is on record rather than hidden in an aggregate.

THREE_DOF_B_GAS = (0.89, 0.99)
THREE_DOF_B_T41 = (0.74, 0.88)
"""The 3-DOF fuel column against B7/B9/B11, ratio ours/printed, at the report's +-2 % step.

Low for one diagnosed reason: Appendix B was generated with a TC_T41 of about 0.32 against
the printed 0.29, and `b(T41)` is a near-cancellation that amplifies it
(`test_heat_sink_taus.py`, open questions #31 and #73). Hence a band of its own for T41."""

SIX_DOF_B_GAS = (0.93, 0.96)
SIX_DOF_B_T41 = (0.93, 0.97)
"""The 6-DOF fuel column against B8/B10/B12, ratio ours/printed; `b(P3)` is printed zero.

The 5-DOF fuel column is exact, and the heat sink is the only difference between the two
models, so the whole deficit here is the heat sink's: the time constants of open question
#31. With the ones Appendix B was generated with, it matches (`test_heat_sink_taus.py`)."""

SIX_DOF_PRESSURE = (1.00, 1.02)
"""Rows P3 and P41 against columns P3 and P41 in B8/B10/B12, at the report's +-2 % step. As
tight as the 5-DOF pressure block: promoting T41 to a state does not disturb the volumes."""

SIX_DOF_REST = (0.80, 1.26)
"""Everything left in the 6-DOF `A`: the NG, NP and P45 rows outside the T41 column, less
`A(P45,P3)` (below) and the Gen Hel `A(NP,NP)` diagonal. The same `f1`/`f7` derivative
spread the 5-DOF carries (open questions #31 and #43). A tripwire for a sign flip or a
dropped term, not evidence of agreement."""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_three_dof_fuel_column_against_b7_b9_b11(trim_no: int):
    """Nine printed numbers, the 3-DOF fuel column."""
    ours, ref = _pair(DOF.THREE, trim_no, c.J_LOAD_UH60A)
    for i, s in enumerate(ref.states):
        assert ref.b[i] != 0.0, f"B{ref.figure}: b({s}) is printed zero, which it was not"
        ratio = ours.b[i] / ref.b[i]
        lo, hi = THREE_DOF_B_T41 if s == "T41" else THREE_DOF_B_GAS
        assert lo < ratio < hi, f"{TRIM_NAME[trim_no]} b({s}) is {ratio:.4f} of printed"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_six_dof_fuel_column_is_low_where_the_five_dof_one_is_exact(trim_no: int):
    """Fifteen printed numbers, asserted in two groups because they are two: the four
    gas-path states move together and T41 sits apart.
    """
    ours, ref = _pair(DOF.SIX, trim_no, c.J_LOAD_UH60A)
    p3 = ref.states.index("P3")
    assert ref.b[p3] == 0.0, "B8/B10/B12 print b(P3) as zero"
    # ours is a `np.linalg.solve` residual at ~1e-11 rather than a literal zero, as the
    # off-row `d` entries are -- compare it against an entry that should be non-zero.
    assert abs(ours.b[p3]) < 1e-9 * np.abs(ours.b).max(), f"b(P3) is {ours.b[p3]:.3g}"
    for i, s in enumerate(ref.states):
        if s == "P3":
            continue
        ratio = ours.b[i] / ref.b[i]
        lo, hi = SIX_DOF_B_T41 if s == "T41" else SIX_DOF_B_GAS
        assert lo < ratio < hi, f"{TRIM_NAME[trim_no]} b({s}) is {ratio:.4f} of printed"


def test_the_six_dof_fuel_column_falls_monotonically_with_power():
    """Every entry is lowest at descent and highest at hover: five independent entries
    agreeing on the ordering is one mechanism, not five coincidences.
    """
    by_trim = {}
    for trim_no in (1, 2, 3):
        ours, ref = _pair(DOF.SIX, trim_no, c.J_LOAD_UH60A)
        by_trim[trim_no] = {
            s: ours.b[i] / ref.b[i] for i, s in enumerate(ref.states) if ref.b[i] != 0.0
        }
    for s in by_trim[1]:
        assert by_trim[1][s] > by_trim[2][s] > by_trim[3][s], (
            f"b({s}) is not monotone in power: "
            f"{by_trim[1][s]:.4f} / {by_trim[2][s]:.4f} / {by_trim[3][s]:.4f}"
        )


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_six_dof_pressure_block_survives_promoting_t41(trim_no: int):
    """The volume dynamics must not care that T41 became a state. Twelve numbers."""
    ours, ref = _pair(DOF.SIX, trim_no, c.J_LOAD_UH60A)
    idx = [ref.states.index(s) for s in ("P3", "P41")]
    for i in idx:
        for j in idx:
            assert ref.A[i, j] != 0.0
            ratio = ours.A[i, j] / ref.A[i, j]
            assert SIX_DOF_PRESSURE[0] < ratio < SIX_DOF_PRESSURE[1], (
                f"{TRIM_NAME[trim_no]} A({ref.states[i]},{ref.states[j]}) {ratio:.4f}"
            )


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_rest_of_the_six_dof_a_matrix_is_bounded_at_all(trim_no: int):
    """The rest of the 6-DOF `A`. A tripwire, not a tolerance."""
    ours, ref = _pair(DOF.SIX, trim_no, c.J_LOAD_UH60A)
    t41 = ref.states.index("T41")
    p3, p45, np_ = ref.states.index("P3"), ref.states.index("P45"), ref.states.index("NP")
    for i, si in enumerate(ref.states):
        if i == t41:
            continue
        for j, sj in enumerate(ref.states):
            # A(NP,NP) is the Gen Hel element: it carries dQreq/dNP, which this file
            # extracts without by convention -- see the circular-closure test above.
            if j == t41 or ref.A[i, j] == 0.0 or (i, j) in {(p45, p3), (np_, np_)}:
                continue
            ratio = ours.A[i, j] / ref.A[i, j]
            assert SIX_DOF_REST[0] < ratio < SIX_DOF_REST[1], (
                f"{TRIM_NAME[trim_no]} A({si},{sj}) is {ratio:.3f} of printed"
            )


def test_the_six_dof_p45_p3_element_is_a_cancellation_and_is_excluded_for_that_reason():
    """`A(P45,P3)`, excluded above: the 6-DOF prints it forty to sixty times smaller than the
    5-DOF does, a residual of cancelling terms printed to four figures (the trap
    `t700.appendix_b` names for row 6). Ours collapses too, which is what is asserted;
    bounding the ratio would be bounding Ballin's rounding.
    """
    for trim_no, ref5_val in ((1, -269.5), (2, -255.1), (3, -261.9)):
        ours, ref = _pair(DOF.SIX, trim_no, c.J_LOAD_UH60A)
        p3, p45 = ref.states.index("P3"), ref.states.index("P45")
        assert abs(ref.A[p45, p3]) < 0.05 * abs(ref5_val), (
            f"B{ref.figure} A(P45,P3) is no longer the small residual this records"
        )
        assert abs(ours.A[p45, p3]) < 0.10 * abs(ref5_val), (
            f"trim {trim_no}: ours is {ours.A[p45, p3]:.4g}, which is not a collapse "
            f"against the 5-DOF's {ref5_val}"
        )


P3_COLUMN_TOL_PCT = 8.0
"""`A(NG,P3)` and `A(P45,P3)` in the 5-DOF, the P3 column outside the pressure block, at the
report's +-2 % step: Eq. 18's secant widens them as it does the pressure block."""

THREE_DOF_NG_COLUMN = (0.76, 1.12)
"""The 3-DOF NG column against B7/B9/B11, ratio ours/printed.

The `f1` derivative ambiguity of open question #43 and, in the T41 row, the report's TC_T41
inconsistency (`test_heat_sink_taus.py`). Descent `A(T41,NG)` is excluded: it has its own
band, `test_all_dof_models.DESCENT_T41_NG` (#73), and one element must not carry two bands
that disagree."""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_p3_column_outside_the_pressure_block(trim_no: int):
    """How NG and P45 respond to compressor discharge pressure. Six printed numbers."""
    ours, ref = _pair(DOF.FIVE, trim_no)
    j = S.index("P3")
    for si in ("NG", "P45"):
        i = S.index(si)
        assert ref.A[i, j] != 0.0
        dev = (ours.A[i, j] - ref.A[i, j]) / abs(ref.A[i, j]) * 100.0
        assert abs(dev) < P3_COLUMN_TOL_PCT, f"trim {trim_no} d({si})/d(P3): {dev:+.3f} %"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_three_dof_ng_column_against_b7_b9_b11(trim_no: int):
    """The 3-DOF NG column. Descent's `A(T41,NG)` is skipped here and bounded by its own
    test; see `THREE_DOF_NG_COLUMN`.
    """
    ours, ref = _pair(DOF.THREE, trim_no, c.J_LOAD_UH60A)
    j = ref.states.index("NG")
    for i, si in enumerate(ref.states):
        if ref.A[i, j] == 0.0:
            continue
        if trim_no == 3 and si == "T41":
            continue  # bounded by its own test; see THREE_DOF_NG_COLUMN
        ratio = ours.A[i, j] / ref.A[i, j]
        assert THREE_DOF_NG_COLUMN[0] < ratio < THREE_DOF_NG_COLUMN[1], (
            f"{TRIM_NAME[trim_no]} A({si},NG) is {ratio:.4f} of printed"
        )


NP_DIAGONAL_RESIDUAL_PCT = (-56.0, -50.0)
"""`A(NP,NP)` without `dQreq/dNP`, in the 5-, 3- and 6-DOF variants. The residual is the
same in every variant, as a term missing from the NP equation alone must be: the power
turbine state is decoupled [pdf p.29]."""


@pytest.mark.parametrize("dof_no,dof", [(5, DOF.FIVE), (3, DOF.THREE), (6, DOF.SIX)])
def test_the_heat_sink_np_diagonals_carry_the_same_missing_load_slope(dof_no, dof):
    """The NP diagonal without dQreq/dNP, the second Gen Hel quantity: about half Ballin's
    at every trim, one missing physical term. If it shrank, dQreq/dNP has been supplied:
    re-measure and tighten rather than widening.
    """

    for trim_no in (1, 2, 3):
        ours, ref = _pair(dof, trim_no, c.J_LOAD_UH60A)
        assert ref.dof == dof_no
        i = ref.states.index("NP")
        dev = (abs(ours.A[i, i]) - abs(ref.A[i, i])) / abs(ref.A[i, i]) * 100.0
        assert NP_DIAGONAL_RESIDUAL_PCT[0] < dev < NP_DIAGONAL_RESIDUAL_PCT[1], (
            f"B{ref.figure} A(NP,NP) residual is {dev:+.2f} % without dQreq/dNP"
        )
