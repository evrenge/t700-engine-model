"""The digitizing library's guards against silent wrong answers (review 2026-09-19, #10-#19).

Each case is one the library used to accept and carry on from with a wrong number.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))

from figure import frame, provenance, ticks  # noqa: E402


def test_parse_poly_refuses_a_foreign_variable():
    assert provenance.parse_poly("466.717 +0.000853938042 y", "y").tolist() == [
        0.000853938042,
        466.717,
    ]
    with pytest.raises(ValueError, match="'u'"):
        provenance.parse_poly("+0.5503129 u +0.29944711", "y")  # was [0.84976]


def test_the_frame_is_read_from_its_own_lines_whatever_the_header_order(tmp_path):
    src = ROOT / "data" / "maps" / "f3_seal_bleed_fraction.csv"
    lines = src.read_text().splitlines(keepends=True)
    first_frame = next(i for i, ln in enumerate(lines) if ln.startswith("#   left   x = "))
    axis_maps = [i for i, ln in enumerate(lines) if ln.startswith("#        left   x = ")]
    assert axis_maps, "the fixture header has changed shape"
    rest = [ln for i, ln in enumerate(lines) if i not in axis_maps]
    moved = rest[:first_frame] + [lines[i] for i in axis_maps] + rest[first_frame:]
    reordered = tmp_path / src.name
    reordered.write_text("".join(moved))
    want, got = provenance.parse_calibration(src), provenance.parse_calibration(reordered)
    assert want is not None and got is not None
    for name in ("left", "right", "top", "bottom"):
        assert np.array_equal(want[0][name], got[0][name]), name


def test_edge_str_refuses_an_intercept_slope_pair():
    assert provenance.edge_str("left", (0.000854, 466.717), "x", "y").startswith("left")
    with pytest.raises(ValueError, match="highest power first"):
        provenance.edge_str("left", (466.717, 0.000854), "x", "y")


@pytest.mark.parametrize("ax", [{"majors": [10.0, 8.0, 6.0]}, {"step": -2.0}])
def test_a_descending_lattice_is_refused(ax):
    with pytest.raises(AssertionError, match="ascending"):
        ticks.major_step(ax)


def test_corners_refuses_an_arc_fit_as_straight():
    arc = {nm: np.array([1e-6, 1e-3, 100.0]) for nm in ("left", "right", "top", "bottom")}
    with pytest.raises(AssertionError, match="higher-order"):
        frame.corners(arc, order=1)


def test_robust_keep_refuses_a_trace_that_is_not_one_line():
    t = np.arange(200.0)
    assert frame.robust_keep(t, 3.0 + 0.01 * t).all()
    with pytest.raises(AssertionError, match="not one line"):
        frame.robust_keep(t, 40.0 * np.sin(t))


def test_a_minor_tick_cannot_stand_in_for_a_hidden_major():
    """Slots every 100 px. The 300 px major is hidden; a minor sits 20 px from its slot
    (0.2 of the step, inside the old 0.25 window). The 500 px major is inside a 129 px run
    of merged ink, as on C25's right edge."""
    t0, t1, frac = 0.0, 600.0, [1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6]
    cands = [(100.0, 9.0, 4), (200.0, 9.0, 4), (320.0, 5.0, 4), (400.0, 9.0, 4), (495.0, 15.0, 129)]
    out, _ = ticks.pick_majors(cands, t0, t1, frac, "test", max_missing=2, verbose=False)
    assert out[:2] == [100.0, 200.0] and out[3] == 400.0
    assert np.isnan(out[2]), "the minor 20 px off was taken as the hidden major"
    assert np.isnan(out[4]), "the merged run was taken as a tick"
