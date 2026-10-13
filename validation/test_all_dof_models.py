"""Every one of Ballin's five linear models, against every reference that exists for it.

| model | eigenvalues | matrices |
|---|---|---|
| 2-DOF | Table 1 column 4 | B1, B3, B5 |
| 3-DOF | none printed | B7, B9, B11 |
| 5-DOF | Table 1 column 3 | B2, B4, B6 |
| 6-DOF | none printed | B8, B10, B12 -- element-wise only, ill-conditioned |
| reduced-5 | Table 1 column 5 | not in Appendix B |

Appendix B's models carry the UH-60A drivetrain, so these comparisons pass `j_load` and
`dq_req_dnp`, both recovered from Appendix B itself: the NP diagonal is circular and checked
only as plumbing (see `test_appendix_b_elements`). Everything else here is independent of it.
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
from test_trim_points import CASES, EIGENVALUE_TOL_PCT

AMB = Ambient(14.696, 518.67)
NP_RPM = 20895.0
WF_PPH = {n: case.wf_pph for n, case in enumerate(CASES, start=1)}
"""Table B.1's fuel flows by trim number, from `test_trim_points.CASES`."""
TRIM = {n: case.name for n, case in enumerate(CASES, start=1)}
DQ_REQ_DNP = {1: 0.019471, 2: 0.015869, 3: 0.012950}
"""dQreq/dNP, ft*lbf per rpm, recovered from Appendix B's NP diagonal (open question #6).
Gen Hel's, not ours. It completes the load specification [Eq. 47] but, recovered from the
elements it then corrects, is not evidence of agreement; see
`test_appendix_b_elements.test_supplying_the_load_slope_closes_the_np_diagonal_circularly`.
"""

TABLE_1_5DOF = {
    # [TM-100991 pdf p.31] read cell by cell from the page image; see
    # docs/notes/body-realtime.md 5.5. Blank cells are blank in the original.
    1: {"NG": -2.66, "NP": -0.565, "P3": -51.6, "P41": -4900.0, "P45": -3060.0},
    2: {"NG": -2.08, "NP": -0.446, "P3": -52.2, "P41": -4640.0, "P45": -4040.0},
    3: {"NG": -1.75, "NP": -0.357, "P3": -52.6, "P41": -4430.0, "P45": -4530.0},
}
"""Table 1 column 3 [pdf p.31], the 5-DOF, every mode. Columns 4 and 5 print the same NP
mode, so only their NG modes are kept below."""
TABLE_1_2DOF = {1: -2.69, 2: -2.23, 3: -1.82}
"""Table 1 column 4 [pdf p.31], the NG mode. The NP mode is the Gen Hel one."""
TABLE_1_RED5 = {1: -2.81, 2: -2.16, 3: -1.83}
"""Table 1 column 5 -- the order-reduced 5-DOF, which is not in Appendix B."""

TABLE_1 = {
    t: {
        DOF.FIVE: TABLE_1_5DOF[t],
        DOF.TWO: {"NG": TABLE_1_2DOF[t], "NP": TABLE_1_5DOF[t]["NP"]},
        DOF.REDUCED_FIVE: {"NG": TABLE_1_RED5[t], "NP": TABLE_1_5DOF[t]["NP"]},
    }
    for t in (1, 2, 3)
}
"""All 27 printed modes [pdf p.31], by trim and model. Each column prints the same NP mode."""

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


def participation(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """`(eigenvalues, P)`: `P[k, i]` is state k's participation in mode i, `|V * inv(V).T|`
    column-normalised, so it does not depend on the states' units."""
    z, V = np.linalg.eig(A)
    P = np.abs(V * np.linalg.inv(V).T)
    return z, P / P.sum(axis=0, keepdims=True)


def modes_by_state(m) -> dict[str, float]:
    """Each mode of `m`, named by the state that dominates it: by physics, never by position
    in a sorted spectrum (see `test_linear.py::test_sorting_would_mispair_the_slow_modes`).
    NP is `A[1,1]` exactly: its column is zero off the diagonal."""
    z, P = participation(m.A)
    k_np = m.states.index("NP")
    picks = {s: int(np.argmax(P[k])) for k, s in enumerate(m.states) if k != k_np}
    assert len(set(picks.values())) == len(picks), f"two states claim one mode: {picks}"
    out = {s: float(np.real(z[i])) for s, i in picks.items()}
    out["NP"] = float(m.A[k_np, k_np])
    return out


def _model(dof: DOF, trim_no: int):
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    r = trim.solve(wf, NP_RPM, AMB)
    assert r.trustworthy
    return extract(r, wf, dof, AMB, j_load=c.J_LOAD_UH60A, dq_req_dnp=DQ_REQ_DNP[trim_no])


# --------------------------------------------------------------- every variant extracts


@pytest.mark.parametrize("dof", [DOF.TWO, DOF.THREE, DOF.FIVE, DOF.SIX])
@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_np_stays_decoupled_in_every_variant(dof: DOF, trim_no: int):
    """ "The power turbine state is completely decoupled" [pdf p.29]: column NP is exactly
    zero off its diagonal, which makes `A[1,1]` usable as the NP mode without sorting.
    """
    m = _model(dof, trim_no)
    off = np.delete(m.A[:, 1], 1)
    assert np.abs(off).max() == 0.0, f"{dof.value} trim {trim_no}: NP column {off}"


# ------------------------------------------------ Table 1, each column's NG mode


NG_MODE_TOL_PCT = {DOF.FIVE: 2.0, DOF.TWO: 6.0, DOF.REDUCED_FIVE: 3.0}
"""Ours, declared in `SCOPE.md`: a characterization set a quarter above today's measurement,
not an accuracy claim."""


@pytest.mark.parametrize("dof", [DOF.FIVE, DOF.TWO, DOF.REDUCED_FIVE], ids=str)
@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_each_model_matches_its_own_table_1_column(dof: DOF, trim_no: int):
    """Table 1's three columns are made three ways [pdf pp.27-28]: the 5-DOF and 2-DOF by
    +-2 % central differences of two different nonlinear simulations, the reduced-order
    5-DOF by setting the 5-DOF's pressure derivatives to zero. `extract` makes each of ours
    the same way, so each NG mode is held to its own column (open question #49).
    """
    ng = modes_by_state(_model(dof, trim_no))["NG"]
    want = TABLE_1[trim_no][dof]["NG"]
    dev = (abs(ng) - abs(want)) / abs(want) * 100.0
    assert abs(dev) < NG_MODE_TOL_PCT[dof], (
        f"{TRIM[trim_no]}: {dof} NG mode {ng:.3f} vs {want}, {dev:+.2f} %"
    )


# --------------------------------------------------------- element-wise, the two gaps


TWO_DOF_A = (0.70, 1.20)
TWO_DOF_B_TOL_PCT = 6.0
"""2-DOF against B1/B3/B5: `A` as a ratio to the printed value, `b` as a deviation."""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_two_dof_elements_against_b1_b3_b5(trim_no: int):
    """Six printed numbers per trim. `A`'s NG column carries the `f1` derivative error (#43),
    hence its wider band.
    """
    m = _model(DOF.TWO, trim_no)
    ref = ab.find(2, trim_no)
    for i in range(2):
        assert abs(m.b[i] - ref.b[i]) / abs(ref.b[i]) * 100.0 < TWO_DOF_B_TOL_PCT, (
            f"{TRIM[trim_no]} b({ref.states[i]}): {m.b[i]:.4g} vs {ref.b[i]:.4g}"
        )
    for i in range(2):
        for j in range(2):
            if ref.A[i, j] == 0.0:
                assert m.A[i, j] == 0.0
                continue
            r = m.A[i, j] / ref.A[i, j]
            assert TWO_DOF_A[0] < r < TWO_DOF_A[1], (
                f"{TRIM[trim_no]} A({ref.states[i]},{ref.states[j]}) ratio {r:.3f}"
            )


THREE_DOF_T41_COLUMN = (0.88, 1.12)
THREE_DOF_FEEDTHROUGH_TOL_PCT = 7.0
"""3-DOF against B7/B9/B11: the T41 column of `A` as a ratio, the T41 feedthrough `d` as a
deviation. With T41 held fixed during extraction [pdf p.31], `d` carries the report's TC_T41
inconsistency in full (`test_heat_sink_taus.py`)."""


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_three_dof_elements_against_b7_b9_b11(trim_no: int):
    """The 3-DOF T41 column and feedthrough; the report prints no 3-DOF eigenvalues."""
    m = _model(DOF.THREE, trim_no)
    ref = ab.load(FIG[(DOF.THREE, trim_no)])
    t41 = ref.states.index("T41")
    for i in range(3):
        if ref.A[i, t41] == 0.0:
            continue
        r = m.A[i, t41] / ref.A[i, t41]
        assert THREE_DOF_T41_COLUMN[0] < r < THREE_DOF_T41_COLUMN[1], (
            f"{TRIM[trim_no]} A({ref.states[i]},T41) ratio {r:.3f}"
        )
    assert m.d is not None and ref.d is not None
    dev = abs(m.d[t41] / ref.d[t41] - 1.0) * 100.0
    assert dev < THREE_DOF_FEEDTHROUGH_TOL_PCT, f"{TRIM[trim_no]} feedthrough {dev:.2f} %"


DESCENT_T41_NG = (0.6, 1.0)
"""`A(T41,NG)` at the descent trim, 3-DOF, as a ratio to B11's printed value.

A band of its own: this element moves with every `f1` revision while the trim under it barely
moves, a pure derivative effect (open question #73). `THREE_DOF_NG_COLUMN` in
`test_appendix_b_elements` skips it rather than assert it twice."""


def test_the_descent_t41_ng_element_is_recorded_not_forgotten():
    """B11's `A(T41,NG)`, built through the NG column where descent sits far from the 85 %
    speed line. Constant-abscissa `f1` interpolation collapses it to near zero; along Figure
    A1's construction lines (`maps.SpeedMap`) it is merely imprecise (#43, #73).
    """
    m = _model(DOF.THREE, 3)
    ref = ab.load(11)
    t41, ng = ref.states.index("T41"), ref.states.index("NG")
    r = m.A[t41, ng] / ref.A[t41, ng]
    assert DESCENT_T41_NG[0] < r < DESCENT_T41_NG[1], (
        f"A(T41,NG) at descent is {r:.3f} of Ballin's; near zero means the f1 interpolation "
        f"has gone back to constant x."
    )
    # The other two trims, imprecise but present, are bounded tighter by
    # `test_appendix_b_elements.py::test_the_three_dof_ng_column_against_b7_b9_b11`.


@pytest.mark.parametrize("t", [1, 2, 3])
@pytest.mark.parametrize("dof", [DOF.FIVE, DOF.TWO, DOF.REDUCED_FIVE], ids=str)
def test_every_printed_table_1_mode(dof, t):
    """All 27 modes of Table 1 [pdf p.31] within `EIGENVALUE_TOL_PCT` (SCOPE.md), paired by
    dominant state. The 9 NP modes are circular: `DQ_REQ_DNP` comes from Appendix B."""
    ours = modes_by_state(_model(dof, t))
    dev = {k: 100.0 * (abs(ours[k]) - abs(v)) / abs(v) for k, v in TABLE_1[t][dof].items()}
    bad = {k: round(d, 2) for k, d in dev.items() if abs(d) >= EIGENVALUE_TOL_PCT}
    assert not bad, f"{dof} trim {t}: {bad} % against Table 1"
