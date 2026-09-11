"""The five linear-model variants, and the structure that makes them comparable.

These tests exist because the eigenvalue comparison was twice recorded wrong in this
repository -- once with numbers no committed code produces, once with the NG and NP modes
paired against each other's rows. Both were caught by re-deriving rather than re-reading.
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import trim
from t700.engine import Ambient
from t700.linear import DOF, STATES, extract
from t700.units import wf_pps_from_pph

AMB = Ambient(14.696, 518.67)
NP_RPM = 20895.0

# Table B.1 [pdf p.67] fuel flows; Table 1 [pdf p.31] eigenvalues by mode.
TRIMS = {
    "hover": (476.3, dict(NG=-2.66, NP=-0.565, P3=-51.6, P41=-4900.0, P45=-3060.0)),
    "level": (349.3, dict(NG=-2.08, NP=-0.446, P3=-52.2, P41=-4640.0, P45=-4040.0)),
    "descent": (267.7, dict(NG=-1.75, NP=-0.357, P3=-52.6, P41=-4430.0, P45=-4530.0)),
}


def _model(name: str, dof: DOF = DOF.FIVE):
    wf = wf_pps_from_pph(TRIMS[name][0])
    r = trim.solve(wf, NP_RPM, AMB)
    assert r.trustworthy
    return extract(r, wf, dof, AMB)


@pytest.mark.parametrize("name", TRIMS)
def test_np_column_is_zero_off_the_diagonal(name: str):
    """The structure Appendix B prints, and the reason NP can be identified at all.

    "The power turbine state is completely decoupled from the other states" [pdf p.29].
    NP drives nothing: no other state's derivative depends on it. So `A[1,1]` is exactly
    an eigenvalue, and that -- not sorting -- is how the NP mode must be found.
    """
    A = _model(name).A
    off = np.delete(A[:, 1], 1)
    assert np.abs(off).max() == 0.0, f"NP column is not decoupled: {off}"
    ev = np.sort(np.linalg.eigvals(A).real)
    assert np.abs(ev - A[1, 1]).min() < 1e-6 * abs(A[1, 1])


@pytest.mark.parametrize("name", TRIMS)
def test_sorting_would_mispair_the_slow_modes(name: str):
    """Guard the *method*, not just the numbers.

    At hover and descent our NP mode is more damped than our NG mode, while Ballin's is
    far less. Pairing sorted eigenvalues against Table 1's row order therefore swaps them
    and reports the two slow modes as roughly +-5 % when they are -14 % and -23 %. This
    test fails if anyone reintroduces positional pairing.
    """
    m = _model(name)
    ev = np.sort(np.linalg.eigvals(m.A).real)
    np_mode = m.A[1, 1]
    assert m.comparable_eigenvalues.size == 4
    assert np_mode not in m.comparable_eigenvalues
    if name in ("hover", "descent"):
        assert ev[-1] != np_mode, (
            "NP is no longer the *slowest* mode, so positional pairing would now "
            "coincidentally work here; re-read this test before trusting it"
        )


@pytest.mark.parametrize("name", TRIMS)
def test_eigenvalues_are_step_independent(name: str):
    """A Jacobian that moves with the differencing step is measuring the step.

    The system spans four decades, 4900 to 0.5 /sec, which is exactly where one-sided
    differences go soft. Central differences are converged over rel = 1e-4 .. 1e-8.
    """
    wf = wf_pps_from_pph(TRIMS[name][0])
    r = trim.solve(wf, NP_RPM, AMB)
    ref = np.sort(np.linalg.eigvals(extract(r, wf, DOF.FIVE, AMB, rel_step=1e-5).A).real)
    for rel in (1e-4, 1e-6, 1e-7, 1e-8):
        ev = np.sort(np.linalg.eigvals(extract(r, wf, DOF.FIVE, AMB, rel_step=rel).A).real)
        assert np.allclose(ev, ref, rtol=1e-3), f"rel={rel:g} moved the spectrum: {ev}"


@pytest.mark.parametrize("name", TRIMS)
def test_two_dof_equals_reduced_five_identically(name: str):
    """Linearizing the exact quasi-steady system *is* residualizing the full Jacobian.

    dg/ds = A11 - A12 A22^-1 A21 by implicit differentiation. Ballin's two differ (-2.69
    vs -2.81 at hover) because his reduced-order model is a separate program with the
    opened Eq. 74 iteration and finite loop counts; ours is a Newton solve to 1e-10, so
    it lands exactly on the Schur complement. See open question #49.
    """
    a = _model(name, DOF.TWO)
    b = _model(name, DOF.REDUCED_FIVE)
    assert np.allclose(a.A, b.A, rtol=1e-6, atol=1e-9), f"\n{a.A}\n{b.A}"
    assert np.allclose(a.b, b.b, rtol=1e-6, atol=1e-9)


@pytest.mark.parametrize("name", TRIMS)
def test_five_dof_fast_modes_match_table_1(name: str):
    """The two pressure modes are the well-conditioned part of the comparison.

    Within 3 % everywhere except hover's slower one, which is the -14.9 % outlier. Bound
    set at 16 % to hold that one; **if it tightens, tighten this** rather than leaving
    slack that hides a regression.
    """
    m = _model(name)
    want = sorted([TRIMS[name][1]["P41"], TRIMS[name][1]["P45"]])
    got = sorted(m.comparable_eigenvalues[:2])
    for o, w in zip(got, want, strict=True):
        dev = (abs(o) - abs(w)) / abs(w) * 100.0
        assert abs(dev) < 16.0, f"{name}: {o:.1f} against {w:.1f}, {dev:+.1f} %"


@pytest.mark.parametrize("name", TRIMS)
def test_ng_mode_is_recorded_not_forgotten(name: str):
    """Characterization. The NG mode is the one a pilot feels and our worst comparison.

    -14.0 % hover, +2.3 % level, -22.6 % descent -- no monotone pattern, good in the
    middle and poor either side. Recorded so it cannot drift unnoticed; the bound is
    deliberately just wide enough to hold today's values.
    """
    m = _model(name)
    ng = m.comparable_eigenvalues[-1]
    want = TRIMS[name][1]["NG"]
    dev = (abs(ng) - abs(want)) / abs(want) * 100.0
    assert abs(dev) < 25.0, f"{name}: NG mode {ng:.3f} against {want}, {dev:+.1f} %"


@pytest.mark.parametrize("dof", [DOF.THREE, DOF.SIX])
def test_heat_sink_variants_refuse_clearly(dof: DOF):
    """Not implemented, and it must say why rather than return something plausible."""
    wf = wf_pps_from_pph(476.3)
    r = trim.solve(wf, NP_RPM, AMB)
    with pytest.raises(NotImplementedError, match="Chen|heat-sink"):
        extract(r, wf, dof, AMB)


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
