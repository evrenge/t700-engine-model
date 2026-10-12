"""Verify Appendix B's transcription against Table 1 -- the report checking itself.

Appendix B prints 297 matrix elements [pdf pp.68-76]. Table 1 prints 27 eigenvalues
[pdf p.31] that were computed from six of those twelve matrices. Neither derives from the
other in this repository, so agreement is an end-to-end check on the transcription:
90 numbers read off six pages reproducing 15 numbers read off a seventh.

None of this involves our model: it tests only that we read the report correctly.
Comparing our own derivatives against these matrices is `validation/`'s job.
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import appendix_b as ab

# Table 1 [pdf p.31], by trim and mode. Blank cells are blank in the original.
TABLE_1_5DOF = {
    1: dict(NG=-2.66, NP=-0.565, P3=-51.6, P41=-4900.0, P45=-3060.0),
    2: dict(NG=-2.08, NP=-0.446, P3=-52.2, P41=-4640.0, P45=-4040.0),
    3: dict(NG=-1.75, NP=-0.357, P3=-52.6, P41=-4430.0, P45=-4530.0),
}
TABLE_1_2DOF = {1: [-2.69, -0.565], 2: [-2.23, -0.446], 3: [-1.82, -0.357]}


def test_all_twelve_load():
    refs = ab.all_references()
    assert len(refs) == 12
    assert sum(r.A.size + r.b.size + (0 if r.d is None else r.d.size) for r in refs) == 297


@pytest.mark.parametrize("trim", [1, 2, 3])
def test_five_dof_eigenvalues_reproduce_table_1(trim: int):
    """The strongest available check on 90 transcribed numbers.

    Every residual is print precision (worst 3.33e-3 relative): Appendix B gives 4
    significant figures, Table 1 gives 3. Table 1 truncates rather than rounds in two places
    (-2.66843 printed "-2.66", -2.08692 "-2.08"), which is why the bound is 4e-3, not 2e-3.
    """
    ev = np.sort(np.linalg.eigvals(ab.find(5, trim).A).real)
    want = np.array(sorted(TABLE_1_5DOF[trim].values()))
    for got, exp in zip(ev, want, strict=True):
        assert abs(got - exp) / abs(exp) < 4e-3, (
            f"trim {trim}: {got:.5f} against Table 1's {exp:.4g}"
        )


@pytest.mark.parametrize("trim", [1, 2, 3])
def test_two_dof_is_lower_triangular_so_its_diagonal_is_its_spectrum(trim: int):
    """A(1,2) = 0 exactly, so Table 1's 2-DOF column is the diagonal of B1/B3/B5, to Table
    1's printing precision (-1.816 printed as "-1.82")."""
    r = ab.find(2, trim)
    diag = sorted(np.diag(r.A))
    for got, exp in zip(diag, sorted(TABLE_1_2DOF[trim]), strict=True):
        assert abs(got - exp) / abs(exp) < 3e-3, f"trim {trim}: {got:.5f} vs {exp}"


@pytest.mark.parametrize("trim", [1, 2, 3])
def test_np_is_decoupled_in_every_variant(trim: int):
    """ "The power turbine state is completely decoupled from the other states" [pdf p.29].

    Column NP is zero off the diagonal in all twelve figures. This is the structural fact
    that lets the NP mode be identified without sorting -- see `t700.linear`.
    """
    for dof in (2, 3, 5, 6):
        r = ab.find(dof, trim)
        off = np.delete(r.A[:, 1], 1)
        assert np.abs(off).max() == 0.0, f"{r.figure}: NP column not decoupled: {off}"


@pytest.mark.parametrize("trim", [1, 2, 3])
def test_the_np_eigenvalue_is_identical_across_all_four_variants(trim: int):
    """A(2,2) is the same number in the 2-, 3-, 5- and 6-DOF models at each trim: NP is
    decoupled, so its mode cannot depend on the rest of the model. A transcription error in
    any of the four breaks this."""
    vals = [ab.find(dof, trim).A[1, 1] for dof in (2, 3, 5, 6)]
    assert len(set(vals)) == 1, f"trim {trim}: A(2,2) differs across variants: {vals}"
    # Appendix B prints 4 figures (-0.3567), Table 1 prints 3 (-0.357).
    np_mode = TABLE_1_5DOF[trim]["NP"]
    assert abs(vals[0] - np_mode) / abs(np_mode) < 3e-3


@pytest.mark.parametrize("fig", [7, 8, 9, 10, 11, 12])
def test_feedthrough_is_nonzero_only_on_the_t41_row(fig: int):
    """`d = F1^-1 G2` and `G2` [pdf p.33] is zero in every row but T41's.

    This is also one of the three pieces of evidence pinning the 3-DOF state ordering:
    rows 1 and 2 of the printed 3-DOF `d` are zero and only row 3 is not.
    """
    r = ab.load(fig)
    assert r.d is not None
    t41 = r.states.index("T41")
    assert r.d[t41] != 0.0
    assert np.abs(np.delete(r.d, t41)).max() == 0.0


def test_the_six_dof_instability_is_recorded_not_repaired():
    """B8 is unstable as printed, at about +3 /sec. Do not "fix" it.

    Row 6 of `F1^-1 F2` is a difference of terms of order 2e5 giving order 1e3, so the
    printed four digits carry about +-50 absolute in A(6,3) and A(6,4); perturbing those
    within their printed rounding moves this mode from -0.99 to +7.3. The instability is
    lost precision, not a claim by the report.
    """
    r = ab.load(8)
    assert r.ill_conditioned
    assert np.linalg.eigvals(r.A).real.max() > 1.0
    for trim in (2, 3):
        assert np.linalg.eigvals(ab.find(6, trim).A).real.max() < 0.0, (
            "trims 2 and 3 are stable as printed; only trim 1 is not"
        )


def test_provenance_metadata_survives_the_round_trip():
    """Every figure knows its page, its model, and whether its states are printed."""
    for r in ab.all_references():
        assert 67 <= r.page <= 76
        assert r.dof in (2, 3, 5, 6)
        assert r.heat_sink == (r.dof in (3, 6))
        assert len(r.states) == r.dof
        assert r.states[:2] == ("NG", "NP")
        # Only the 3-DOF ordering is inferred [inventory-appendix-b.md 2.4].
        assert r.states_are_printed == (r.dof != 3)
        assert (r.d is None) == (not r.heat_sink)


def test_matrices_are_read_only():
    """Loaded arrays are cached and shared; a caller must not be able to mutate them."""
    with pytest.raises(ValueError):
        ab.load(2).A[0, 0] = 0.0
