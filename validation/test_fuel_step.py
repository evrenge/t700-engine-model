"""Compare the real-time model against Ballin's published fuel-step transients.

Figures 9 and 10 [pdf pp.45-46] are the report's open-loop transient evidence: a step up
from 400 to 775 lbm/hr and a step down from 400 to 125 lbm/hr, both with "the power
turbine speeds ... held constant by suppressing the NP integration" [pdf p.39].

The comparison target is **Ballin's own model trace** (the solid line), not the GE
reference markers. We are replicating his model, so his output is what ours must match;
the GE data is the hardware reference he was himself validating against, and he does not
match it perfectly either.

## Which configuration these run

**Heat sink on.** Figures 9 and 10 were generated with the station 4.1 heat-sink model
active -- inferred, not printed; the evidence and its strength are laid out in
`docs/notes/heat-sink-configuration.md`. Running our 5-DOF (`T41 = T41ns`) model against
them is a category error, and was the leading cause of transients that ran 1.4-3.5x too
fast. Table 1 and Appendix B figures B1-B6 are the opposite case, heat sink **off**, and
are compared elsewhere against the 5-DOF model.

Because Eq. 50 has unit DC gain, the switch cannot move a trim -- so the initial-state
tests below are configuration-independent by construction, and only the transient shape
and settled state respond to it.

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


def _run(fig: int, wf_hi: float, t_step: float, heat_sink: bool = True):
    wf0 = wf_pps_from_pph(400.0)
    r0 = trim.solve(wf0, c.NP_DES, AMB)
    f0 = frame(r0.state, wf0, AMB)
    hi = wf_pps_from_pph(wf_hi)
    st = realtime.from_trim(r0, f0.wa31_pps, f0)
    return realtime.run(
        st,
        lambda t: wf0 if t < t_step else hi,
        AMB,
        duration_s=5.0,
        dt=0.007,
        q_req_ftlbf=f0.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=heat_sink,
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
    # The 4.6 upper bound is load-bearing: fig09_t41_model.csv and fig09_torq45_model.csv
    # each carry a spurious trailing sample at t ~ 4.9 (see their headers, and open
    # question #48). Widening this window silently corrupts the comparison.
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

    The step to 125 lbm/hr drives NGc down to about 71 %, against a compressor map whose
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

    # 70.94 % with the heat sink on, against 67.64 % without it -- the lead-lag slows the
    # decay, so the engine no longer plunges as deep, and Ballin's own trace bottoms at
    # 74.2 %. Still below the 75 % where `f1` has real speed lines either side, so the
    # extrapolation argument below survives; it is simply less severe than it was.
    assert ngc.min() < 72.0, "the step down should still reach the bottom of the speed range"
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


def _t41_overshoot(t: np.ndarray, v: np.ndarray, t_step: float) -> tuple[float, float]:
    """Peak T41 above the settled value. Windowed at t < 4.6 -- see `late_them` above."""
    settled = v[(t > 4.0) & (t < 4.6)].mean()
    return v[(t > t_step) & (t < 4.6)].max() - settled, settled


def test_heat_sink_closes_most_of_the_t41_overshoot_gap():
    """The heat sink is the larger part of the shape mismatch, but not all of it.

    History, because the numbers only mean something against it. With `T41 = T41ns` our
    overshoot was 400.8 degR against Ballin's 117.4 -- 3.41x. That was recorded as a
    characterization test with a note saying that if the heat-sink model ever landed, the
    bound was to be **re-measured and tightened, not relaxed**. It landed on 2026-09-11
    and this is that re-measurement: 196.1 degR, 1.67x.

    So Eqs. 48-53 account for roughly 65 % of the peak error and the remaining 1.67x is
    unexplained. Candidates not yet eliminated, in no order: the compressor map's
    digitized transient envelope, the `lag_whole_flow` reading of Eq. 74 (open question
    #22), and the possibility that Figures 9/10 used the NASA-Lewis test-engine function
    set rather than the specification set -- the report states that substitution for
    Tables 2/3 [pdf p.39] but says nothing about these figures.

    This stays a characterization test. If it fails low, something improved and the bound
    should be tightened again rather than widened.
    """
    tr = _run(9, 775.0, 0.539)
    ours, _ = _t41_overshoot(tr["t"], tr["t41"], 0.539)
    tb, vb = _reference(9, "t41")
    theirs, _ = _t41_overshoot(tb, vb, 0.539)
    ratio = ours / theirs

    assert 1.4 < ratio < 2.0, (
        f"T41 overshoot is {ours:.0f} degR against Ballin's {theirs:.0f}, a ratio of "
        f"{ratio:.2f}x. On record is 1.67x (196.1 vs 117.4) with the heat sink on, down "
        f"from 3.41x without it. Below this band something improved -- re-measure and "
        f"tighten. Above it, something regressed."
    )


def test_the_heat_sink_is_what_closed_it():
    """Guard the attribution itself, not just the number.

    The claim above is that Eqs. 48-53 are responsible. That is checkable directly by
    running the same step both ways, and it is worth checking: if the heat sink were
    silently disabled, `test_heat_sink_closes_most_of_the_t41_overshoot_gap` alone would
    fail with no indication of why.
    """
    bare = _run(9, 775.0, 0.539, heat_sink=False)
    sunk = _run(9, 775.0, 0.539, heat_sink=True)
    off, _ = _t41_overshoot(bare["t"], bare["t41"], 0.539)
    on, _ = _t41_overshoot(sunk["t"], sunk["t41"], 0.539)
    assert on < 0.6 * off, (
        f"the heat sink should roughly halve the T41 overshoot; it went {off:.0f} -> {on:.0f} degR"
    )
