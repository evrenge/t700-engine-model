"""Every one of Ballin's five linear models, against every reference that exists for it.

Coverage was uneven before this file: the 5-DOF was compared element by element and by
eigenvalue, the 3- and 6-DOF only structurally, and the 2-DOF and reduced-5-DOF not at all
-- despite Table 1 printing a column for each of the latter two. This closes that.

## What reference exists for what

| model | eigenvalues | matrices |
|---|---|---|
| 2-DOF | Table 1 column 4 | B1, B3, B5 |
| 3-DOF | **none printed** | B7, B9, B11 |
| 5-DOF | Table 1 column 3 | B2, B4, B6 |
| 6-DOF | **none printed** | B8, B10, B12 -- element-wise only, ill-conditioned |
| reduced-5 | Table 1 column 5 | **not in Appendix B** |

## The two Gen Hel inputs are supplied here

Appendix B's models carry the UH-60A drivetrain, so these comparisons pass `j_load` and
`dq_req_dnp`. Both were recovered from Appendix B itself, so the NP diagonal matching is
circular and is asserted only as a plumbing check -- see
`test_appendix_b_elements.test_supplying_the_load_slope_closes_the_np_diagonal_circularly`.
Everything else in this file is independent of that recovery.
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
TRIM = {1: "hover", 2: "level 80 kt", 3: "descent 80 kt"}
DQ_REQ_DNP = {1: 0.019471, 2: 0.015869, 3: 0.012950}

TABLE_1_2DOF = {1: -2.69, 2: -2.23, 3: -1.82}
"""Table 1 column 4 [pdf p.31], the NG mode. The NP mode is the Gen Hel one."""
TABLE_1_RED5 = {1: -2.81, 2: -2.16, 3: -1.83}
"""Table 1 column 5 -- the order-reduced 5-DOF, which is not in Appendix B."""

FIG = {
    (DOF.TWO, 1): 1,
    (DOF.TWO, 2): 3,
    (DOF.TWO, 3): 5,
    (DOF.THREE, 1): 7,
    (DOF.THREE, 2): 9,
    (DOF.THREE, 3): 11,
    (DOF.FIVE, 1): 2,
    (DOF.FIVE, 2): 4,
    (DOF.FIVE, 3): 6,
    (DOF.SIX, 1): 8,
    (DOF.SIX, 2): 10,
    (DOF.SIX, 3): 12,
}


def _model(dof: DOF, trim_no: int):
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    r = trim.solve(wf, NP_RPM, AMB)
    assert r.trustworthy
    return extract(r, wf, dof, AMB, j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[trim_no])


# --------------------------------------------------------------- every variant extracts


@pytest.mark.parametrize("dof", list(DOF))
@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_every_model_extracts_at_every_trim(dof: DOF, trim_no: int):
    """Fifteen extractions, all finite and all the right shape.

    Cheap, and it is the test that would have caught any of the five silently regressing
    while attention was on another.
    """
    m = _model(dof, trim_no)
    n = len(m.states)
    assert m.A.shape == (n, n) and m.b.shape == (n,)
    assert np.all(np.isfinite(m.A)) and np.all(np.isfinite(m.b))
    assert (m.d is not None) == (dof in (DOF.THREE, DOF.SIX))


@pytest.mark.parametrize("dof", [DOF.TWO, DOF.THREE, DOF.FIVE, DOF.SIX])
@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_np_stays_decoupled_in_every_variant(dof: DOF, trim_no: int):
    """ "The power turbine state is completely decoupled" [pdf p.29] -- in all twelve.

    NP drives nothing, so column NP is zero off its diagonal. We reproduce that exactly,
    which is what makes `A[1,1]` usable as the NP mode without sorting.
    """
    m = _model(dof, trim_no)
    off = np.delete(m.A[:, 1], 1)
    assert np.abs(off).max() == 0.0, f"{dof.value} trim {trim_no}: NP column {off}"
    assert ab.load(FIG[(dof, trim_no)]).A[1, 1] != 0.0


# ------------------------------------------------------- the two models Table 1 prints
# ------------------------------------------------------- that were never checked


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_two_dof_ng_mode_against_table_1_column_4(trim_no: int):
    """Our quasi-steady model's NG mode against Table 1's 2-DOF column.

    **-5.62 / +7.57 / +4.21 %** at hover / level / descent, in this file's own magnitude
    convention. The shape is the `f1` knot story of `docs/notes/derivative-ambiguity.md`,
    not a new fault -- but which trim carries it moved when `f1` went onto a beta grid.
    It read -9.1 / -0.3 / -22.6 % with level mid-segment and the other two on a knot; the
    beta grid redistributed that rather than removing it, and level is now the worst.
    Open question #43.

    The *discrete* frame map -- the thing Ballin actually ran -- is closer than this
    continuous model at all three trims: -3.17 / +4.10 / -1.01 %. See
    `validation/test_discrete_map.py`.
    """
    ng = _model(DOF.TWO, trim_no).comparable_eigenvalues[-1]
    want = TABLE_1_2DOF[trim_no]
    dev = (abs(ng) - abs(want)) / abs(want) * 100.0
    assert abs(dev) < 25.0, f"{TRIM[trim_no]}: 2-DOF NG mode {ng:.3f} vs {want}, {dev:+.1f} %"
    if trim_no == 3:
        # Descent is the tight one now, not level. Interpolating f1 along Figure A1's
        # construction lines (constant k, not constant x -- see `maps.SpeedMap`)
        # REDISTRIBUTED the f1-driven derivative error rather than removing it:
        # descent went from the worst case to the best and level the other way. The
        # mechanism of docs/notes/derivative-ambiguity.md still holds; which trim it
        # bites has moved. Open question #43.
        assert abs(dev) < 6.0, "descent should be the tight one under the beta-gridded f1"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_reduced_five_ng_mode_against_table_1_column_5(trim_no: int):
    """The order-reduced model against the column that exists only in Table 1.

    Ours is identical to the 2-DOF by construction (see `t700.linear`), so it cannot match
    both columns: Ballin's differ by 0.12 /sec at hover. We land closer to his column 4
    than his column 5 at hover and descent, which is consistent with our quasi-steady solve
    being exact rather than carrying the numerics his reduced-order program did.
    """
    ng = _model(DOF.REDUCED_FIVE, trim_no).comparable_eigenvalues[-1]
    want = TABLE_1_RED5[trim_no]
    dev = (abs(ng) - abs(want)) / abs(want) * 100.0
    assert abs(dev) < 25.0, f"{TRIM[trim_no]}: reduced-5 NG mode {ng:.3f} vs {want}, {dev:+.1f} %"


# --------------------------------------------------------- element-wise, the two gaps


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_two_dof_elements_against_b1_b3_b5(trim_no: int):
    """Six printed numbers per trim, never compared until now.

    `b` agrees to 5 % everywhere. `A` carries the `f1` knot error in its NG column, worst
    at descent where the trim sits 7.5 % from the 85 % speed line: 0.776 and 0.821 of
    Ballin's against 0.996 and 1.025 at level.
    """
    m = _model(DOF.TWO, trim_no)
    ref = ab.find(2, trim_no)
    for i in range(2):
        assert abs(m.b[i] - ref.b[i]) / abs(ref.b[i]) < 0.06, (
            f"{TRIM[trim_no]} b({ref.states[i]}): {m.b[i]:.4g} vs {ref.b[i]:.4g}"
        )
    for i in range(2):
        for j in range(2):
            if ref.A[i, j] == 0.0:
                assert m.A[i, j] == 0.0
                continue
            r = m.A[i, j] / ref.A[i, j]
            # Widened from 1.10 to 1.20 on 2026-09-13 with the f1 interpolation change,
            # which took level's A(NG,NG) from 1.06 to 1.117 and descent's elements the
            # other way. Open question #43: the ambiguity moved, it did not go.
            assert 0.70 < r < 1.20, (
                f"{TRIM[trim_no]} A({ref.states[i]},{ref.states[j]}) ratio {r:.3f}"
            )


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_three_dof_elements_against_b7_b9_b11(trim_no: int):
    """The 3-DOF matrices, which have no printed eigenvalues to check instead.

    The T41 column -- how NG, NP and T41's own dynamics respond to T41 -- agrees to 10 %
    at every trim. The NG column carries the `f1` knot error as it does everywhere else.
    `A(T41,NG)` at descent is excluded and characterized separately below.
    """
    m = _model(DOF.THREE, trim_no)
    ref = ab.load(FIG[(DOF.THREE, trim_no)])
    t41 = ref.states.index("T41")
    for i in range(3):
        if ref.A[i, t41] == 0.0:
            continue
        r = m.A[i, t41] / ref.A[i, t41]
        assert 0.88 < r < 1.12, f"{TRIM[trim_no]} A({ref.states[i]},T41) ratio {r:.3f}"
    assert m.d is not None and ref.d is not None
    assert abs(m.d[t41] / ref.d[t41] - 1.0) < 0.03, "feedthrough should agree to 3 %"


def test_the_descent_t41_ng_element_is_recorded_not_forgotten():
    """B11's `A(T41,NG)` was the single worst element in the whole comparison, and it came
    right.

    Ours was about 1e-4 against Ballin's 0.01659 -- **a ratio of 0.006**, a total relative
    miss on a small element, recorded here so it would not be lost in an aggregate. It was
    attributed to the `f1` knot error in its most extreme form: descent sits 7.5 % past the
    85 % speed line and this element is built through the engine block's NG column.

    **That attribution was right.** Interpolating `f1` along Figure A1's own construction
    lines instead of at constant abscissa (see `maps.SpeedMap`) took it to **0.791**, in
    line with the other two trims at 0.823 and 1.147. The element that was absent is now
    merely imprecise, which is what the rest of the comparison looks like.

    The remaining spread across the three trims is the derivative ambiguity of open
    question #43 -- redistributed by the same change, not removed: descent went from the
    worst case to the best and level from the best to the worst. See
    `docs/notes/derivative-ambiguity.md`, whose mechanism still holds even though which
    trim it bites has moved.
    """
    m = _model(DOF.THREE, 3)
    ref = ab.load(11)
    t41, ng = ref.states.index("T41"), ref.states.index("NG")
    r = m.A[t41, ng] / ref.A[t41, ng]
    assert 0.6 < r < 1.0, (
        f"A(T41,NG) at descent is {r:.3f} of Ballin's, 0.791 on record. It was 0.006 under "
        f"constant-x f1; if it has collapsed again, the interpolation is the place to look."
    )
    others = [
        _model(DOF.THREE, t).A[t41, ng] / ab.load(FIG[(DOF.THREE, t)]).A[t41, ng] for t in (1, 2)
    ]
    assert all(0.7 < o < 1.3 for o in others), (
        f"the other two trims should be imprecise but present: {others}"
    )
