"""The core must produce bit-identical output for identical input.

Not "close". Identical. A real-time engine model that drifts between runs cannot be
regression-tested against digitized reference traces at all, because every deviation
becomes ambiguous: the model, or the machine?

`test_architecture.py` enforces this structurally by banning clocks and RNGs. This file
checks it behaviourally.

**It did not grow as the model did.** Until 2026-09-13 it exercised `thermo` (two calls)
and `units` (a round trip) and nothing else -- no `engine.frame`, no `realtime.step` or
`run`, no `trim.solve`, no `control.loop`. Those are the only places where the things that
could actually break determinism live: `maps._clamps`, the cached table loaders,
`thermo.set_backend`, and the carried state of the opened iterations. CLAUDE.md lists
determinism as `live` and names this file as half the mechanism, so for the stateful half
of the model that claim was resting on the import ban alone. Found by the 2026-09-13
code-quality audit.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from t700 import constants as c
from t700 import maps, realtime, thermo, trim
from t700.control import loop
from t700.engine import Ambient
from t700.units import wf_pps_from_pph


def test_scalar_and_array_paths_agree_exactly():
    """A scalar call and an array call must give the same bits, not merely close values.

    They diverge the moment someone writes a fast path, and a silent scalar/array
    mismatch would show up much later as an unreproducible validation failure.
    """
    values = [600.0, 1234.5, 2100.0]
    scalar = np.array([thermo.h41_from_t41(v) for v in values])
    array = thermo.h41_from_t41(np.array(values))
    assert np.array_equal(scalar, array)
    assert array.dtype == np.float64
    assert thermo.theta41_from_t41(np.array(values)).dtype == np.float64


# --------------------------------------------------------------------------------------
# The stateful half. Every one of these runs the model twice with something in between
# that a non-deterministic implementation would notice, and demands bit-identical output.
# --------------------------------------------------------------------------------------

AMB = Ambient(14.696, 518.67)
WF0 = wf_pps_from_pph(400.0)


def _figure_9_trace():
    r = trim.solve(WF0, c.NP_DES, AMB)
    st = realtime.from_trim(r)
    return realtime.run(
        st,
        lambda t: WF0 if t < 0.5 else wf_pps_from_pph(775.0),
        AMB,
        duration_s=1.0,
        q_req_ftlbf=r.frame.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=True,
    )


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


def _closed_loop():
    r = trim.solve(WF0, c.NP_DES, AMB)
    pilot = loop.Pilot(xcpc_pct=52.75, pas_deg=100.0, pcprf_pct=100.0)
    s = loop.seed(r, pilot, AMB)
    out = []
    for _ in range(400):
        s, _e, h, fr = loop.step(s, pilot, AMB, dt=0.007, load=r.frame.q_pt_ftlbf, heat_sink=True)
        out.append((s.engine.ng_rpm, s.engine.np_rpm, h.wf_pph, fr.t41_degR))
    return out


def test_the_closed_loop_is_bit_identical_to_itself():
    """The control carries 22 continuous states and 4 memory elements, and three of the
    blocks are path-dependent (backlash). If anything in the core reads a clock, a run
    this long is where it shows."""
    assert _closed_loop() == _closed_loop(), "the closed loop is not reproducible"


def _digest() -> str:
    """SHA-256 over the bits of a trim, a real-time transient and a closed-loop run."""
    h = hashlib.sha256()
    h.update(repr(trim.solve(WF0, c.NP_DES, AMB).state).encode())
    trace = _figure_9_trace()
    for key in sorted(trace):
        h.update(key.encode())
        h.update(np.ascontiguousarray(trace[key], dtype=np.float64).tobytes())
    h.update(np.array(_closed_loop(), dtype=np.float64).tobytes())
    return h.hexdigest()


def test_output_is_bit_identical_across_processes():
    """Every test above runs both trials in one process, where `hash()` and set order are
    fixed. A planted `hash()`-dependent factor gave three different NG values in three
    processes with this file green (audit, 2026-09-19). Two fresh interpreters with
    different hash seeds must match this one."""
    root = Path(__file__).resolve().parent.parent
    code = "import test_determinism as t; print(t._digest())"
    digests = {"this process": _digest()}
    for seed in ("0", "4242"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        env["PYTHONPATH"] = os.pathsep.join([str(root / "src"), str(root / "tests")])
        out = subprocess.run(
            [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
        )
        digests[f"PYTHONHASHSEED={seed}"] = out.stdout.strip()
    assert len(set(digests.values())) == 1, f"output depends on the process: {digests}"
