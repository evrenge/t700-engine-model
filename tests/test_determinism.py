"""The core must produce bit-identical output for identical input.

Not "close". Identical. A real-time engine model that drifts between runs cannot be
regression-tested against digitized reference traces at all, because every deviation
becomes ambiguous: the model, or the machine?

`test_imports.py` enforces this structurally by banning clocks and RNGs. This file
checks it behaviourally, and grows as the model does.
"""

from __future__ import annotations

import numpy as np

from t700 import thermo, units


def test_thermo_is_bit_identical_across_calls():
    t = np.linspace(500.0, 3000.0, 257)
    first = np.concatenate([thermo.h41_from_t41(t), thermo.theta41_from_t41(t)])
    second = np.concatenate([thermo.h41_from_t41(t), thermo.theta41_from_t41(t)])
    assert np.array_equal(first, second), "identical input gave different output"
    assert first.dtype == np.float64


def test_scalar_and_array_paths_agree_exactly():
    """A scalar call and an array call must give the same bits, not merely close values.

    They diverge the moment someone writes a fast path, and a silent scalar/array
    mismatch would show up much later as an unreproducible validation failure.
    """
    values = [600.0, 1234.5, 2100.0]
    scalar = np.array([thermo.h41_from_t41(v) for v in values])
    array = thermo.h41_from_t41(np.array(values))
    assert np.array_equal(scalar, array)


def test_unit_conversions_round_trip_exactly_where_they_should():
    """lbm/hr <-> lbm/sec is exact in binary for these magnitudes."""
    for pph in (125.0, 267.7, 349.3, 400.0, 476.3, 775.0):
        assert units.wf_pph_from_pps(units.wf_pps_from_pph(pph)) == pph
