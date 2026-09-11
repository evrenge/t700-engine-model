"""Compare the real-time model against Ballin's published fuel-step transients.

Figures 9 and 10 [pdf pp.45-46] are the report's open-loop transient evidence: a step up
from 400 to 775 lbm/hr and a step down from 400 to 125 lbm/hr, both with "the power
turbine speeds ... held constant by suppressing the NP integration" [pdf p.39].

The comparison target is **Ballin's own model trace** (the solid line), not the GE
reference markers. We are replicating his model, so his output is what ours must match;
the GE data is the hardware reference he was himself validating against, and he does not
match it perfectly either.

## What these tests establish, and what they do not

They compare the **initial trim** and the **settled final state** on each panel. They do
not yet compare the shape of the transition, which needs the traces resampled onto a
common time base and a decision about how to treat the step-time uncertainty (the step
time is not printed; ours measures 0.539 s and 0.545 s -- open question #17).

Tolerances here are **ours**. The report states no percentage tolerance for any transient
quantity; its only numeric transient claim is a 1-2 % NG overestimate on the open-loop
step [pdf pp.45-46].
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from t700 import constants as c
from t700 import maps, realtime, trim
from t700.engine import Ambient, frame
from t700.units import wf_pps_from_pph

REF = Path(__file__).resolve().parent.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)

# Figure, final fuel flow, measured step time
STEPS = [(9, 775.0, 0.539), (10, 125.0, 0.545)]

# Panels comparable to a model output, and how to form that output from a run.
PANELS = {
    "pcng": lambda tr: 100.0 * tr["ng"] / c.NG_DES,
    "ps3": lambda tr: c.K_PS3 * tr["p3"],
    "t41": lambda tr: tr["t41"],
    "t45": lambda tr: tr["t45"],
    "torq45": lambda tr: tr["q_pt"],
}

INITIAL_TOL_PCT = 1.5
"""The trim is the same physics the steady-state tests already check, so it should agree
tightly -- this is really a check that the figures' own y calibration is sound."""

FINAL_TOL_PCT = {9: 5.0, 10: 20.0}
"""Figure 10 gets a much looser bound, and deliberately: see
`test_step_down_runs_off_the_bottom_of_the_maps`, which measures why."""


UNTRUSTED = {(10, "torq45")}
"""Panels excluded from the comparison, with the evidence.

Both figures begin at the same 400 lbm/hr trim, so their t = 0 values must agree with each
other regardless of any model. Five of six panels do, to within 0.62 %. Figure 10's
`TORQ45` disagrees with Figure 9's by **5.10 %** -- and that panel also yielded only 16
markers against 45 everywhere else. The reference data is wrong there, not the model, so
using it would be measuring our digitizer.
"""


def _reference(fig: int, key: str) -> tuple[np.ndarray, np.ndarray]:
    path = REF / f"fig{fig:02d}_{key}_model.csv"
    rows = list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))
    t = np.array([float(r["t_s"]) for r in rows])
    v = np.array([float(r["value"]) for r in rows])
    o = np.argsort(t)
    return t[o], v[o]


def _run(fig: int, wf_hi: float, t_step: float):
    wf0 = wf_pps_from_pph(400.0)
    r0 = trim.solve(wf0, c.NP_DES, AMB)
    f0 = frame(r0.state, wf0, AMB)
    hi = wf_pps_from_pph(wf_hi)
    st = realtime.from_trim(r0, f0.wa31_pps)
    return realtime.run(
        st,
        lambda t: wf0 if t < t_step else hi,
        AMB,
        duration_s=5.0,
        dt=0.007,
        q_req_ftlbf=f0.q_pt_ftlbf,
        integrate_np=False,
    )


@pytest.mark.parametrize("fig,wf_hi,t_step", STEPS, ids=lambda v: str(v))
@pytest.mark.parametrize("key", list(PANELS))
def test_initial_state_matches_the_figure(fig: int, wf_hi: float, t_step: float, key: str):
    """Before the step, our engine and Ballin's should be at the same trim."""
    if (fig, key) in UNTRUSTED:
        pytest.skip(f"fig {fig} {key}: reference data untrusted, see UNTRUSTED")
    path = REF / f"fig{fig:02d}_{key}_model.csv"
    if not path.exists():
        pytest.skip(f"{path.name} not digitized")
    tb, vb = _reference(fig, key)
    tr = _run(fig, wf_hi, t_step)
    ours = PANELS[key](tr)

    pre_us = ours[tr["t"] < t_step - 0.05]
    pre_them = vb[tb < t_step - 0.05]
    if pre_them.size == 0:
        pytest.skip("no pre-step reference samples")
    dev = (pre_us.mean() - pre_them.mean()) / abs(pre_them.mean()) * 100.0
    assert abs(dev) < INITIAL_TOL_PCT, (
        f"fig {fig} {key}: initial {pre_us.mean():.1f} against Ballin's "
        f"{pre_them.mean():.1f}, {dev:+.2f} %"
    )


@pytest.mark.parametrize("fig,wf_hi,t_step", STEPS, ids=lambda v: str(v))
@pytest.mark.parametrize("key", list(PANELS))
def test_settled_state_matches_the_figure(fig: int, wf_hi: float, t_step: float, key: str):
    if (fig, key) in UNTRUSTED:
        pytest.skip(f"fig {fig} {key}: reference data untrusted, see UNTRUSTED")
    path = REF / f"fig{fig:02d}_{key}_model.csv"
    if not path.exists():
        pytest.skip(f"{path.name} not digitized")
    tb, vb = _reference(fig, key)
    tr = _run(fig, wf_hi, t_step)
    ours = PANELS[key](tr)

    late_us = ours[tr["t"] > 4.0]
    late_them = vb[(tb > 4.0) & (tb < 4.6)]
    if late_them.size == 0:
        pytest.skip("no settled reference samples")
    dev = (late_us.mean() - late_them.mean()) / abs(late_them.mean()) * 100.0
    assert abs(dev) < FINAL_TOL_PCT[fig], (
        f"fig {fig} {key}: settled {late_us.mean():.1f} against Ballin's "
        f"{late_them.mean():.1f}, {dev:+.2f} %"
    )


def test_step_up_stays_inside_the_compressor_map():
    """Figure 9 keeps the engine on data, which is why it agrees to a few percent."""
    maps.reset_clamps()
    tr = _run(9, 775.0, 0.539)
    ngc = 100.0 * tr["ng"] / c.NG_DES
    assert ngc.min() > 88.0, "the step up should not approach the map's lower edge"
    assert ngc.max() < 101.0
    rep = maps.clamp_report()
    assert not any(k.startswith("f1@") and v > 50 for k, v in rep.items()), (
        f"the compressor map should barely be extrapolated on the step up: {rep}"
    )
    maps.reset_clamps()


def test_step_down_runs_off_the_bottom_of_the_maps():
    """Figure 10 disagrees far more, and this measures why rather than asserting it.

    The step to 125 lbm/hr drives NGc down to about 68 %, against a compressor map whose
    lowest speed line is 65 %, and drives the fuel-air ratio to 0.005 against an `f6`
    table that starts at 0.010. The model spends the whole late transient extrapolating
    at the bottom of its own data, so the 9-16 % deviations there are a statement about
    the digitized envelope, not about the equations.

    Figure 10's reference data is also the weaker of the two: its `WFPH` settles 4.35 %
    from the value its caption states, against 0.23 % for Figure 9, and its `TORQ45` panel
    yielded only 16 markers against 45 elsewhere and is not used.
    """
    maps.reset_clamps()
    tr = _run(10, 125.0, 0.545)
    ngc = 100.0 * tr["ng"] / c.NG_DES
    rep = maps.clamp_report()

    assert ngc.min() < 70.0, "the step down should reach the bottom of the speed range"
    assert tr["far"].min() < 0.010, "and below f6's tabulated fuel-air ratio"
    assert rep.get("f1@65", 0) > 100, (
        "the lowest compressor speed line should be heavily extrapolated here; if it is "
        "not, the explanation for Figure 10's disagreement needs revisiting"
    )
    maps.reset_clamps()


def test_the_two_figures_step_in_opposite_directions():
    """Guards against the two runs being accidentally identical."""
    up = _run(9, 775.0, 0.539)
    down = _run(10, 125.0, 0.545)
    assert up["ng"][-1] > up["ng"][0]
    assert down["ng"][-1] < down["ng"][0]


@pytest.mark.parametrize("key", list(PANELS) + ["wfph"])
def test_the_two_figures_agree_at_their_shared_trim(key: str):
    """A model-free check on the reference data itself.

    Figures 9 and 10 both start from the same 400 lbm/hr trim, so their t = 0 values must
    agree with each other whatever our model says. This caught Figure 10's TORQ45 panel
    disagreeing with Figure 9's by 5.10 % while every other panel agreed to 0.62 % -- and
    that is how we know the fault is in that panel's digitization rather than in the
    engine.
    """
    vals = []
    for fig in (9, 10):
        path = REF / f"fig{fig:02d}_{key}_model.csv"
        if not path.exists():
            pytest.skip(f"{path.name} not digitized")
        t, v = _reference(fig, key)
        vals.append(v[t < 0.48].mean())
    dev = (vals[1] - vals[0]) / abs(vals[0]) * 100.0
    if key == "torq45":
        assert abs(dev) > 2.0, (
            "Figure 10's TORQ45 panel is on file as disagreeing with Figure 9's by 5.1 %. "
            "If it now agrees, it has been re-digitized and should leave UNTRUSTED."
        )
        return
    assert abs(dev) < 1.5, (
        f"{key}: the two figures' shared trim disagrees by {dev:+.2f} % -- one of the "
        f"two digitizations is wrong"
    )


def test_transient_shape_gap_is_characterized_not_forgotten():
    """Record the shape mismatch so it cannot regress silently, and so it flags when fixed.

    The endpoint comparisons above pass at 0.9-3.9 % while the *shape* is visibly wrong:
    our engine responds faster than Ballin's and spikes turbine temperature far harder.
    Plotting found this; percentages did not.

    The leading suspect is the heat-sink model (Eqs. 48-53), which is not implemented --
    see open question #47. This is a **characterization** test: it asserts the gap that
    exists today. When the heat sink lands it should FAIL, and that failure is the signal
    to re-measure and tighten it.
    """
    tr = _run(9, 775.0, 0.539)
    t41 = tr["t41"]
    settled = t41[tr["t"] > 4.0].mean()
    overshoot = t41.max() - settled

    assert overshoot > 250.0, (
        f"T41 overshoot is now {overshoot:.0f} degR against the {400.8:.0f} on record. "
        f"If the heat-sink model has been added, re-measure against Ballin's +120.6 and "
        f"tighten this bound -- do not simply relax it."
    )
    assert overshoot < 600.0, "overshoot has grown; something regressed"
