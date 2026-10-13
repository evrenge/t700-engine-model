"""The five linear-model variants, and the structure that makes them comparable.

The eigenvalue comparison is re-derived here rather than re-read: it is easy to record with
numbers no code produces, or with the NG and NP modes paired against each other's rows.
"""

from __future__ import annotations

import functools

import numpy as np
import pytest

from t700 import trim
from t700.engine import Ambient
from t700.linear import DOF, STATES, extract
from t700.units import wf_pps_from_pph
from test_appendix_b import TABLE_1_5DOF

AMB = Ambient(14.696, 518.67)
NP_RPM = 20895.0

# Table B.1 [pdf p.67] fuel flows; Table 1 [pdf p.31] eigenvalues by mode, transcribed once
# in `test_appendix_b.py`.
TRIMS = {
    "hover": (476.3, TABLE_1_5DOF[1]),
    "level": (349.3, TABLE_1_5DOF[2]),
    "descent": (267.7, TABLE_1_5DOF[3]),
}


@functools.cache
def _model(name: str, dof: DOF = DOF.FIVE, rel_step: float = 1e-5):
    wf = wf_pps_from_pph(TRIMS[name][0])
    r = trim.solve(wf, NP_RPM, AMB)
    assert r.trustworthy
    return extract(r, wf, dof, AMB, rel_step=rel_step)


@pytest.mark.parametrize("name", TRIMS)
def test_sorting_would_mispair_the_slow_modes(name: str):
    """Guard the method, not just the numbers.

    At hover our NP mode is more damped than our NG mode, while Ballin's is far less, so
    pairing sorted eigenvalues against Table 1's row order swaps the two slow modes. This
    fails if positional pairing comes back.
    """
    m = _model(name)
    ev = np.sort(np.linalg.eigvals(m.A).real)
    np_mode = m.A[1, 1]
    assert m.comparable_eigenvalues.size == 4
    assert np_mode not in m.comparable_eigenvalues
    # Hover only: at descent NP is the slowest mode, so positional pairing works there by
    # coincidence -- the coincidence this test warns about.
    if name == "hover":
        assert ev[-1] != np_mode, (
            "NP is no longer the slowest mode at hover; positional pairing now works by "
            "coincidence everywhere, re-read this test"
        )


@pytest.mark.parametrize("name", TRIMS)
def test_eigenvalues_are_step_independent(name: str):
    """A Jacobian that moves with the differencing step is measuring the step.

    The system spans four decades, 4900 to 0.5 /sec, where one-sided differences go soft.
    Central differences are converged over rel = 1e-4 .. 1e-8.

    The report's step, +/-2 % of equilibrium [pdf p.28], is three decades above that plateau:
    an amplitude matched to real rotor-speed excursions, so a secant rather than the
    derivative. The last assertion pins the separation so the two cannot be conflated (open
    question #56).
    """
    wf = wf_pps_from_pph(TRIMS[name][0])
    r = trim.solve(wf, NP_RPM, AMB)
    ref = np.sort(np.linalg.eigvals(extract(r, wf, DOF.FIVE, AMB, rel_step=1e-5).A).real)
    for rel in (1e-4, 1e-6, 1e-7, 1e-8):
        ev = np.sort(np.linalg.eigvals(extract(r, wf, DOF.FIVE, AMB, rel_step=rel).A).real)
        assert np.allclose(ev, ref, rtol=1e-3), f"rel={rel:g} moved the spectrum: {ev}"

    report = np.sort(np.linalg.eigvals(extract(r, wf, DOF.FIVE, AMB, rel_step=0.02).A).real)
    shifts = [100.0 * (b / a - 1.0) for a, b in zip(ref, report, strict=True)]
    worst = max(shifts, key=abs)
    # Recorded values. The gap is a second-derivative effect, so an `f1` revision moves only
    # a trim whose perturbation lands on the changed part: hover sits short of its knee on
    # the flat choked run, which every revision draws the same; descent sits past its knee.
    expected = {"hover": 21.6, "level": -7.4, "descent": -1.5}[name]
    assert abs(worst - expected) < 1.0, (
        f"{name}: the report's +/-2 % perturbation [pdf p.28] moves the 5-DOF spectrum by "
        f"{worst:+.1f} % against rel=1e-5, {expected:+.1f} % on record. If this has "
        f"collapsed, open question #56 can be closed."
    )


@pytest.mark.parametrize("name", TRIMS)
def test_two_dof_equals_reduced_five_at_a_tangent_step(name: str):
    """With exact derivatives, linearizing the quasi-steady system *is* residualizing the
    full Jacobian: dg/ds = A11 - A12 A22^-1 A21 by implicit differentiation.

    Checked at rel 1e-6, where both finite differences are converged, then at `_model`'s
    default 1e-5 as a growth alarm. Neither is the shipped extraction step, `linear.py`'s
    `rel_step` (the report's two percent): a central difference that wide carries
    truncation, and the quasi-steady `g(s)`, which slaves the pressures, curves differently
    in NG from the 5-DOF `f(x)`, so the two routes part there, NG's diagonal included, as
    Ballin's Table 1 columns 4 and 5 do (open question #49).
    """
    a = _model(name, DOF.TWO, rel_step=1e-6)
    b = _model(name, DOF.REDUCED_FIVE, rel_step=1e-6)
    assert np.allclose(a.A, b.A, rtol=1e-6, atol=1e-9), f"\n{a.A}\n{b.A}"
    assert np.allclose(a.b, b.b, rtol=1e-6, atol=1e-9)

    # `_model`'s default step, not the shipped one: a growth alarm, so a failure points at
    # the step rather than the data.
    a5 = _model(name, DOF.TWO)
    b5 = _model(name, DOF.REDUCED_FIVE)
    gap = float(np.max(np.abs(a5.A - b5.A) / np.maximum(np.abs(a5.A), 1e-30)))
    assert gap < 2e-3, (
        f"{name}: at rel_step=1e-5 the two routes to the 2-DOF A differ by "
        f"{gap:.2e}; the step is no longer adequate."
    )


@pytest.mark.parametrize("name", TRIMS)
def test_five_dof_modes_against_table_1(name: str):
    """The two pressure modes, then the NG mode, against Table 1.

    The pressure modes are the well-conditioned part. The 16 % bound holds hover's slower
    one, the -14.9 % outlier; if it tightens, tighten this. The NG mode, the one a pilot
    feels and our worst comparison, is a characterization: its bound is just wide enough to
    hold today's values, so it cannot drift unnoticed.
    """
    m = _model(name)
    want = sorted([TRIMS[name][1]["P41"], TRIMS[name][1]["P45"]])
    got = sorted(m.comparable_eigenvalues[:2])
    for o, w in zip(got, want, strict=True):
        dev = (abs(o) - abs(w)) / abs(w) * 100.0
        assert abs(dev) < 16.0, f"{name}: {o:.1f} against {w:.1f}, {dev:+.1f} %"

    ng = m.comparable_eigenvalues[-1]
    want = TRIMS[name][1]["NG"]
    dev = (abs(ng) - abs(want)) / abs(want) * 100.0
    assert abs(dev) < 25.0, f"{name}: NG mode {ng:.3f} against {want}, {dev:+.1f} %"


@pytest.mark.parametrize("dof", [DOF.THREE, DOF.SIX])
def test_heat_sink_variants_extract_with_the_right_shape(dof: DOF):
    """These two carry a `Wf-dot` feedthrough. Shape and structure only: the comparison
    against Figures B7-B12 is in `validation/test_appendix_b_elements.py`."""
    wf = wf_pps_from_pph(476.3)
    r = trim.solve(wf, NP_RPM, AMB)
    m = extract(r, wf, dof, AMB)
    n = len(STATES[dof])
    assert m.A.shape == (n, n) and m.b.shape == (n,)
    assert m.states[-1] == "T41", "T41 is the appended state [pdf pp.29, 32]"
    assert m.d is not None, "the heat-sink models have a Wf-dot feedthrough [Eq. 67]"
    t41 = m.states.index("T41")
    assert abs(m.d[t41]) > 0.0
    assert np.abs(np.delete(m.d, t41)).max() < 1e-12 * abs(m.d[t41])


@pytest.mark.parametrize("dof", [DOF.TWO, DOF.FIVE, DOF.REDUCED_FIVE])
def test_non_heat_sink_variants_have_no_feedthrough(dof: DOF):
    """`d` exists only where a `Wf-dot` term had to be absorbed -- Eq. 65's G2 comes from
    differentiating the heat-sink lead, so without the heat sink there is nothing to
    absorb. Appendix B prints `C` and `d` only on B7-B12, and this mirrors that."""
    wf = wf_pps_from_pph(476.3)
    r = trim.solve(wf, NP_RPM, AMB)
    assert extract(r, wf, dof, AMB).d is None


def test_extract_refuses_an_untrustworthy_trim():
    """A clamped solve has no usable derivatives; linearizing it is the whole trap."""
    wf = wf_pps_from_pph(900.0)
    r = trim.solve(wf, NP_RPM, AMB)
    if r.trustworthy:
        pytest.skip("900 lbm/hr now trims on data; pick a hotter fixture")
    with pytest.raises(ValueError, match="untrustworthy"):
        extract(r, wf, DOF.FIVE, AMB)


def test_state_names_match_the_report():
    """5-DOF [pdf p.27] and 6-DOF [pdf p.32] are printed; 3-DOF ordering is inferred."""
    assert STATES[DOF.FIVE] == ("NG", "NP", "P3", "P41", "P45")
    assert STATES[DOF.SIX] == ("NG", "NP", "P3", "P41", "P45", "T41")
    assert STATES[DOF.TWO] == STATES[DOF.REDUCED_FIVE] == ("NG", "NP")


def test_dof_is_addressable_by_string():
    """`extract(..., dof="2dof")` must work -- this is the switch the project asked for."""
    assert DOF("5dof") is DOF.FIVE
    assert DOF("reduced5") is DOF.REDUCED_FIVE
    assert {d.value for d in DOF} == {"2dof", "3dof", "5dof", "6dof", "reduced5"}
