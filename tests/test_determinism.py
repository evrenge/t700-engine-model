"""The core must produce bit-identical output for identical input.

Not "close". Identical. A real-time engine model that drifts between runs cannot be
regression-tested against digitized reference traces at all, because every deviation
becomes ambiguous: the model, or the machine?

`test_architecture.py` enforces this structurally by banning clocks and RNGs. This file
checks it behaviourally.

**It did not grow as the model did.** Until 2026-09-13 it exercised `thermo` (two calls)
and `units` (a round trip) and nothing else -- no `engine.frame`, no `realtime.step` or
`run`, no `trim.solve`, no `control.loop`. Those are the only places where the things that
could actually break determinism live: `maps._clamps`, the `_cached` table loaders,
`thermo.set_backend`, and the carried state of the opened iterations. CLAUDE.md lists
determinism as `live` and names this file as half the mechanism, so for the stateful half
of the model that claim was resting on the import ban alone. Found by the 2026-09-13
code-quality audit.
"""

from __future__ import annotations

import numpy as np

from t700 import constants as c
from t700 import maps, realtime, thermo, trim, units
from t700.control import loop
from t700.engine import Ambient, frame
from t700.units import wf_pps_from_pph


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


# --------------------------------------------------------------------------------------
# The stateful half. Every one of these runs the model twice with something in between
# that a non-deterministic implementation would notice, and demands bit-identical output.
# --------------------------------------------------------------------------------------

AMB = Ambient(14.696, 518.67)
WF0 = wf_pps_from_pph(400.0)


def _figure_9_trace():
    r = trim.solve(WF0, c.NP_DES, AMB)
    f = frame(r.state, WF0, AMB)
    st = realtime.from_trim(r, f.wa31_pps, f)
    return realtime.run(
        st,
        lambda t: WF0 if t < 0.5 else wf_pps_from_pph(775.0),
        AMB,
        duration_s=1.0,
        q_req_ftlbf=f.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=True,
    )


def test_a_real_time_run_is_bit_identical_to_itself():
    a, b = _figure_9_trace(), _figure_9_trace()
    for key in a:
        assert np.array_equal(a[key], b[key]), f"trace {key!r} differs between runs"


def test_a_real_time_run_is_bit_identical_across_an_intervening_clamp_count():
    """`maps._clamps` is the one mutable module global CLAUDE.md sanctions in the core.

    The guarantee that makes it acceptable is that no model output can depend on it.
    `tests/test_maps.py` checks that for a trim; this checks it for a whole transient,
    which is where the counter is actually hammered.
    """
    a = _figure_9_trace()
    maps.reset_clamps()
    for _ in range(500):
        maps.f9()(0.95)  # outside f9's table, so every call counts
    b = _figure_9_trace()
    for key in a:
        assert np.array_equal(a[key], b[key]), f"trace {key!r} moved with the clamp counter"
    maps.reset_clamps()


def test_a_trim_is_bit_identical_across_an_intervening_solve():
    """The cached loaders and the continuation seed are both process-global."""
    first = trim.solve(WF0, c.NP_DES, AMB).state
    trim.sweep([800.0, 700.0, 600.0])  # a different ladder, in between
    second = trim.solve(WF0, c.NP_DES, AMB).state
    assert first == second, f"{first}\n{second}"


def test_the_closed_loop_is_bit_identical_to_itself():
    """The control carries 22 continuous states and 4 memory elements, and three of the
    blocks are path-dependent (backlash). If anything in the core reads a clock, a run
    this long is where it shows."""

    def go():
        r = trim.solve(WF0, c.NP_DES, AMB)
        f = frame(r.state, WF0, AMB)
        pilot = loop.Pilot(xcpc_pct=52.75, pas_deg=100.0, pcprf_pct=100.0)
        s = loop.seed(r, f.wa31_pps, f, pilot, AMB)
        out = []
        for _ in range(400):
            s, _e, h, fr = loop.step(s, pilot, AMB, dt=0.007, load=f.q_pt_ftlbf, heat_sink=True)
            out.append((s.engine.ng_rpm, s.engine.np_rpm, h.wf_pph, fr.t41_degR))
        return out

    assert go() == go(), "the closed loop is not reproducible"
