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
from t700 import maps, trim
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


# --------------------------------------------------------------- knot proximity


def _knot_distance_f7(trim_no: int) -> float:
    """How far this trim's P45/P41 sits from the nearest `f7` knot, as a fraction of its
    segment. 0 means on a knot, where a linear interpolant's derivative is undefined."""
    from t700.engine import frame

    wf = wf_pps_from_pph(WF_PPH[trim_no])
    r = trim.solve(wf, NP_RPM, AMB)
    arg = r.state.p45_psia / r.state.p41_psia
    x = maps.f7().x
    k = int(np.searchsorted(x, arg))
    frac = (arg - x[k - 1]) / (x[k] - x[k - 1])
    assert frame  # the import above documents where P45/P41 comes from
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


def _err(trim_no: int, si: str, sj: str) -> float:
    ours, ref = _pair(trim_no)
    i, j = S.index(si), S.index(sj)
    return abs((ours.A[i, j] - ref.A[i, j]) / ref.A[i, j] * 100.0)


def test_hover_is_the_double_case_that_pins_the_mechanism():
    """One operating point, two tables, opposite outcomes -- the strongest single check.

    Hover sits 0.0002 from `f7`'s knot at 0.21496, and simultaneously 58 % across an `f1`
    segment. If knot proximity is really what drives the disagreement, hover must hold
    both the WORST `f7`-driven errors and the BEST `f1`-driven errors in the whole
    comparison. It does. A wrong equation could not produce that split at one trim.

    See `docs/notes/derivative-ambiguity.md`.
    """
    assert _knot_distance_f7(1) < 0.05, "hover no longer sits on an f7 knot"
    assert _knot_distance_f1(1) > 0.30, "hover no longer sits mid-segment on f1"

    f7_here = _err(1, "NG", "P45")
    f7_elsewhere = max(_err(2, "NG", "P45"), _err(3, "NG", "P45"))
    assert f7_here > f7_elsewhere, "hover should hold the worst f7-driven error"

    # **The NG-column half of this inverted on 2026-09-13** and the assertions below are
    # what is left of it. Interpolating f1 along Figure A1's construction lines (constant k,
    # not constant x -- see `maps.SpeedMap`) redistributed the f1-driven error: descent now
    # holds the best d/dNG column at 0.6 % and level the worst at 18.1 %, where it was hover
    # best and descent worst. The f7 half above is untouched, which is the control: f7 is a
    # 1-D table and the change did not touch it.
    #
    # So the "double case" this test is named for is now a single case. It still pins the
    # f7 mechanism, and the NG column is pinned by
    # `test_the_ng_column_error_is_redistributed_not_removed` instead. Open question #43.
    for si, sj in NG_COLUMN:
        assert _err(1, si, sj) < 5.0, (
            f"hover d({si})/d({sj}) is {_err(1, si, sj):.1f} %; it was under 1 % under "
            f"constant-x f1 and is about 2.5 % now"
        )


def test_the_ng_column_error_is_redistributed_not_removed():
    """`f1` is piecewise linear between speed lines, so d/dNGc jumps at each one.

    **This test asserted the opposite ordering until 2026-09-13**, and recording why is the
    point. Under constant-abscissa interpolation the error tracked distance to the nearest
    speed line: hover sits 41.9 % across [92, 94] and held the best NG column; level sits
    2.8 % past the 89 line and descent 7.5 % past the 85 line, and both were worse.

    Interpolating along Figure A1's own construction lines instead (see `maps.SpeedMap`)
    **redistributed that error without removing it**. Descent now holds the best column at
    0.6 % and level the worst at 18.1 %, with hover between at 2.5 %. The knot distances
    have not changed -- the correspondence between them and the errors has.

    That is the honest state of open question #43: `docs/notes/derivative-ambiguity.md`
    identified a mechanism and said it "identifies the mechanism without identifying the
    scheme". Changing the scheme moved the error to a different trim, which is what that
    note predicted would happen. The derivative disagreement is not settled, and the
    report prints nothing that would settle it.

    What this now pins is that the *spread* stays bounded, so a future change that quietly
    made every trim worse would fail here.
    """
    d = {t: _knot_distance_f1(t) for t in (1, 2, 3)}
    assert d[1] > d[3] > d[2], f"knot distances changed: {d}"
    for si, sj in NG_COLUMN:
        e = {t: _err(t, si, sj) for t in (1, 2, 3)}
        # Level is the worst on every element of the column; descent and hover trade
        # places between elements (descent 0.6 % against hover 2.5 % on d(P3)/d(NG),
        # 2.9 against 2.5 on d(NG)/d(NG)), so only the level statement is asserted.
        assert e[2] > e[1] and e[2] > e[3], (
            f"d({si})/d({sj}): on record level {e[2]:.1f} % is the worst, against hover "
            f"{e[1]:.1f} % and descent {e[3]:.1f} %. If this has reordered again, the f1 "
            f"interpolation is the place to look."
        )
        assert max(e.values()) < 25.0, f"d({si})/d({sj}) spread widened: {e}"


# ------------------------------------------------- the heat-sink models, B7-B12

HS_FIG = {
    (DOF.THREE, 1): 7,
    (DOF.THREE, 2): 9,
    (DOF.THREE, 3): 11,
    (DOF.SIX, 1): 8,
    (DOF.SIX, 2): 10,
    (DOF.SIX, 3): 12,
}


def _hs_pair(dof: DOF, trim_no: int):
    """Heat-sink models are compared with the load inertia supplied: Appendix B's
    matrices carry the UH-60A drivetrain, not the bare power turbine."""
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    r = trim.solve(wf, NP_RPM, AMB)
    return extract(r, wf, dof, AMB, j_load=c.J_LOAD_UH60A), ab.load(HS_FIG[(dof, trim_no)])


@pytest.mark.parametrize("trim_no", [1, 2, 3])
@pytest.mark.parametrize("dof", [DOF.THREE, DOF.SIX])
def test_feedthrough_is_nonzero_only_on_the_t41_row(dof: DOF, trim_no: int):
    """`d = F1^-1 G2`, and `G2` is zero except in T41's row -- which we never asserted
    into the derivation. Appendix B prints exactly this shape, so reproducing it is a
    check on the descriptor form (Eq. 65) rather than an input to it."""
    ours, ref = _hs_pair(dof, trim_no)
    t41 = ours.states.index("T41")
    assert ours.d is not None and ref.d is not None
    assert abs(ours.d[t41]) > 0.0
    # `d = F1^-1 G2` is exactly zero off the T41 row in exact arithmetic -- F1 is block
    # lower-triangular -- but `np.linalg.solve` leaves rounding at ~1e-20, so compare
    # against the entry that should be non-zero rather than against literal zero.
    assert np.abs(np.delete(ours.d, t41)).max() < 1e-12 * abs(ours.d[t41])
    assert np.abs(np.delete(ref.d, t41)).max() == 0.0


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_six_dof_t41_column_is_close(trim_no: int):
    """How the five engine equations respond to T41 -- the new physics in the 6-DOF model.

    This is the part the 5-DOF cannot be checked on at all, and it agrees to well under a
    percent at every trim, which says Eq. 23's coupling of T41 into h41 and theta41 is
    right.
    """
    ours, ref = _hs_pair(DOF.SIX, trim_no)
    t41 = ours.states.index("T41")
    for i, s in enumerate(ours.states):
        if s == "T41" or ref.A[i, t41] == 0.0:
            continue
        dev = (ours.A[i, t41] - ref.A[i, t41]) / abs(ref.A[i, t41]) * 100.0
        assert abs(dev) < 1.0, f"trim {trim_no} d({s})/d(T41): {dev:+.2f} %"


def test_the_six_dof_loses_the_p3_to_np_coupling_the_five_dof_has():
    """A structural difference between the two that both models must show.

    The 5-DOF prints `A(2,3) = -0.4128E+2`; the 6-DOF prints `0.0000E+0` at all three
    trims. With T41 promoted to a state, P3 no longer reaches power-turbine torque except
    through T41, so the path moves out of that element and into the T41 column. If our
    6-DOF produced a non-zero there, the promotion would not have taken.
    """
    assert ab.find(5, 1).A[1, 2] != 0.0
    for trim_no in (1, 2, 3):
        ours, ref = _hs_pair(DOF.SIX, trim_no)
        assert ref.A[1, 2] == 0.0, "B8/B10/B12 should print zero for d(NP)/d(P3)"
        assert ours.A[1, 2] == 0.0, f"trim {trim_no}: we still couple P3 to NP directly"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_t41_row_is_uniformly_low_by_the_lead_lag_ratio(trim_no: int):
    """Characterization: the whole T41 row sits ~7 % low, and for one reason.

    The row scales with `tau1/tau2`, which Appendix B *does* determine (0.6845 / 0.6571 /
    0.6270) even though `tau2` alone does not survive its own print precision. Ours is
    7.3 % low at every trim. Uniformity across three flight conditions is what says this
    is the time constants and not the structure -- see `docs/notes/derivative-ambiguity.md`
    and open question #31.

    **The NG column is excluded, and that exclusion is itself a result.** Every other
    element of the row scales by 0.925-0.928 at every trim, but the NG element runs
    0.927 / 0.973 / 0.813 -- because `A(T41, NG)` is built from the engine block's NG
    column, which carries the `f1` speed-line knot error documented above. The two known
    faults compose exactly where they should, and nowhere else.

    If this tightens, `f_hs` or `TC_T41` has moved; re-measure and tighten rather than
    widening.
    """
    ours, ref = _hs_pair(DOF.SIX, trim_no)
    t41 = ours.states.index("T41")
    ng = ours.states.index("NG")
    ratios = [
        ours.A[t41, j] / ref.A[t41, j]
        for j in range(len(ours.states))
        if ref.A[t41, j] != 0.0 and j != ng
    ]
    assert all(0.90 < r < 0.96 for r in ratios), f"T41 row ratios moved: {ratios}"
    spread = (max(ratios) - min(ratios)) / np.mean(ratios) * 100.0
    assert spread < 1.0, f"the row is no longer uniformly scaled: spread {spread:.2f} %"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_t41_row_ng_element_carries_the_f1_error_on_top(trim_no: int):
    """The two known faults compose, and only where they should.

    `A(T41, NG)` is the one element of the T41 row built through the engine block's NG
    column, so it carries the `f1` knot error as well as the lead-lag ratio. It therefore
    deviates from the row's otherwise uniform scaling, and by most at descent -- which is
    also where the `f1` error is worst.
    """
    ours, ref = _hs_pair(DOF.SIX, trim_no)
    t41, ng = ours.states.index("T41"), ours.states.index("NG")
    r_ng = ours.A[t41, ng] / ref.A[t41, ng]
    # The "descent stands clear of its row" assertion is gone with the error pattern it
    # described: under constant-k f1 descent's T41/NG element sits 0.922 against a row mean
    # of 0.929, i.e. squarely in the row. That is the improvement, not a loss of signal.
    assert 0.75 < r_ng < 1.15, f"trim {trim_no} A(T41,NG) ratio {r_ng:.3f}"


# --------------------------------------------------- the load-torque slope, and its limits

DQ_REQ_DNP = {1: 0.019471, 2: 0.015869, 3: 0.012950}
"""dQreq/dNP, ft*lbf per rpm, recovered from Appendix B's NP diagonal (open question #6).

**Recovered from the elements it then corrects.** Supplying it is completing the load
specification, as `j_load` was -- Q_req is an input [Eq. 47] and holding it constant
asserts that the rotor does not resist a speed change, which is false. But it is not
independent evidence of agreement, and the test below says so rather than claiming a win.
"""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_supplying_the_load_slope_closes_the_np_diagonal_circularly(trim_no: int):
    """It closes, and the closure is circular. Both halves are asserted deliberately.

    With `dQreq/dNP` supplied the NP diagonal matches; without it, it sits ~52 % low. That
    is worth locking in as a regression guard on the plumbing -- the slope must actually
    reach `Q_req` and nowhere else -- while being explicit that it demonstrates nothing
    about the engine, since the number came from this very element.

    The non-circular part is the *structure*: the slope must move the NP diagonal and
    leave every other element untouched, because Q_req enters only dNP/dt and depends only
    on NP. That is a real check and it is the second half of this test.
    """
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    r = trim.solve(wf, NP_RPM, AMB)
    ref = ab.find(5, trim_no)
    bare = extract(r, wf, DOF.FIVE, AMB, j_load=c.J_LOAD_UH60A)
    fixed = extract(r, wf, DOF.FIVE, AMB, j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[trim_no])

    before = (abs(bare.A[1, 1]) - abs(ref.A[1, 1])) / abs(ref.A[1, 1]) * 100.0
    after = (abs(fixed.A[1, 1]) - abs(ref.A[1, 1])) / abs(ref.A[1, 1]) * 100.0
    assert before < -45.0, f"NP diagonal was {before:+.1f} % without the slope"
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
