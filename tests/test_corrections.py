"""Check the standard-day corrections, including the evidence T_std and P_std rest on.

T_STD_DEGR and P_STD_PSIA are the only two numbers in the model core that are a decision
rather than a transcription (open question #15); these tests pin the reasoning.
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import corrections as corr


def test_table_b1_sea_level_trim_gives_unity_corrections():
    """The evidence for T_std and P_std.

    Table B.1 [pdf p.67] prints the sea-level trim ambient as P2 = 14.696 PSIA,
    T2 = 518.67 deg R. A corrected quantity is unity at standard sea level by definition,
    so those two printed numbers pin the two unprinted ones.
    """
    assert corr.theta2(518.67) == pytest.approx(1.0)
    assert corr.delta2(14.696) == pytest.approx(1.0)


def test_corrections_move_the_right_way():
    """Hot day: theta > 1. High altitude: delta < 1."""
    assert corr.theta2(560.0) > 1.0
    assert corr.theta2(480.0) < 1.0
    assert corr.delta2(10.0) < 1.0


def test_physical_flow_uncorrects_by_delta_over_root_theta():
    """Uncorrecting must apply W = Wc*delta/sqrt(theta), the inverse of the correction.
    An inverted correction is invisible at sea level and wrong everywhere else."""
    wc = np.array([2.5, 7.5, 12.0])
    for t2, p2 in ((518.67, 14.696), (560.0, 10.9), (470.0, 4.4)):
        th, de = corr.theta2(t2), corr.delta2(p2)
        assert np.allclose(corr.physical_flow(wc, th, de), wc * de / np.sqrt(th))


def test_physical_flow_is_not_the_identity_off_standard_day():
    """Guard against someone 'simplifying' the correction away."""
    wc = 10.0
    th, de = corr.theta2(560.0), corr.delta2(10.9)
    assert corr.physical_flow(wc, th, de) != pytest.approx(wc, rel=1e-3)
