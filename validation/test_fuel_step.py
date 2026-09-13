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
time is not printed; ours measures 0.539 s and 0.545 s -- open question #37).

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

FINAL_TOL_PCT = {9: 5.0, 10: 8.0}
"""Figure 10 still gets a looser bound, but far less loose than it needed before.

It was 20.0 while the settled T41 ran +11.5 % and T45 +17.4 %. Adopting the report's own
printed convergence criterion [pdf p.37] brought those to +3.2 % and +3.4 %, and the
worst remaining panel is Ps3 at +6.9 %. **Tightened, never widened** -- and this pair
belongs in `SCOPE.md` per CLAUDE.md rather than inline here."""


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
    """Figure 10 runs outside its own digitized data, and that is where its error is.

    ## The history, because it is the point of this test

    Three values of the NGc minimum have been on record, and the differences are entirely
    in how hard the pressure loops were converged:

    | inner P3/P41 stopping rule | NGc bottoms at |
    |---|---|
    | `tol=1e-10` on the step, cap 40 -- invented, pre-2026-09-12 | 70.94 % |
    | `tol=1e-3` on the **step** -- the printed number, the wrong quantity | 74.64 % |
    | `tol=1e-3` on the **error** -- what the report states [pdf p.37] | **69.92 %** |

    The middle row looked like a triumph: Ballin's own trace bottoms at 74.2 %. It was
    not. Testing the iterate step rather than the error under-converges the loop by a
    factor of `rho/(1-rho)` ~ 8 (see `realtime.TOL_PRESSURE`), and the resulting lag in
    the pressures damped the plunge into near-agreement. Converging the loop to the error
    the report actually prints puts the floor back at 69.92 %, close to where the
    seven-orders-too-tight version had it -- which is the giveaway, since 1e-10 and a
    correctly applied 1e-3 should agree, and they do.

    ## What this costs, and where it goes

    Figure 10's whole-curve RMS moved pcng 2.35 -> 4.21 %, t41 5.18 -> 5.78 %, t45
    7.62 -> 13.15 %. Figure 9, the accel, improved on four of five panels over the same
    change. The asymmetry is the finding: **the chop's error is the map extrapolation,
    not the numerics.**

    Measured over the run, with the loop converged:

        f1's 65 % speed line clamped   389 times   (it is the lower bracket below 80 %NGc)
        f6 clamped                     425 times   (FAR falls to 0.00527 against a table
                                                    starting at 0.00999)
        f9 clamped                     722 times
        f8 clamped                      38 times

    `f1` has no digitized speed line between 65 and 80 % -- a 15-point hole -- so
    everything below 80 %NGc interpolates across it, and below 74 % the 65 line is being
    asked for pressure ratios outside its own range as well. The chop now spends its whole
    floor in that band. **Fixing `f1`'s low-speed interpolation is the next piece of work**,
    and it is the thing that should bring Figure 10 back.

    Figure 10's reference data is also the weaker of the two: its `WFPH` settles 4.35 %
    from the value its caption states, against 0.23 % for Figure 9, and its `TORQ45`
    panel yielded only 16 markers against 45 elsewhere and is not used.
    """
    maps.reset_clamps()
    tr = _run(10, 125.0, 0.545)
    ngc = 100.0 * tr["ng"] / c.NG_DES
    rep = maps.clamp_report()

    # Ours bottoms at 69.92 %, Ballin's at 74.2 %. The window is tight on both sides on
    # purpose: a floor that rises back toward 74 % would most likely mean the pressure
    # loops have stopped converging again rather than that the model improved, and a
    # deeper one would mean the fuel cut is no longer being followed at all.
    assert 68.5 < ngc.min() < 71.5, (
        f"NGc bottoms at {ngc.min():.2f} %, against Ballin's 74.2 % and our 69.92 % on "
        f"record. If this has risen, check `realtime.TOL_PRESSURE`'s stopping rule before "
        f"believing the model got better."
    )
    assert tr["far"].min() < 0.010, "and below f6's tabulated fuel-air ratio"
    assert rep.get("f1@65", 0) > 100, (
        "f1's 65 % line is the lower bracket below 80 % NGc and is asked for pressure "
        "ratios outside its own range; if that stops happening, the extrapolation "
        "caveat on every Figure 10 number needs revisiting"
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
    """The T41 overshoot against Figure 9, and a ratchet on it.

    Two model changes moved this, both of them replacements of an invention by something
    printed, and neither aimed at the figure:

    * 2026-09-12, the pressure iterations were given the report's own 0.1 percent
      convergence criterion [pdf p.37] in place of our invented 1e-10.
    * 2026-09-12, the heat sink was rebuilt on Eqs. 48-49 -- integrating the station 4.1
      METAL temperature -- instead of Eq. 50's collapsed lead-lag memory. The collapse is
      only valid for constant coefficients, and Eqs. 51 and 53 make them functions of T41
      and W41. See `realtime._heat_sink`.

    The overshoot went 3.41x Ballin's without the heat sink, 1.67x with it and the old
    realization, and is now BELOW his. The band is a ratchet: leaving it in either
    direction means re-measuring, and the lower bound exists so an improvement is noticed
    rather than silently absorbed.
    """
    tb, vb = _reference(9, "t41")
    pre_them = float(np.median(vb[tb < 0.45]))
    late_them = float(np.median(vb[(tb > 3.5) & (tb < 4.5)]))
    theirs = float(vb.max()) - late_them

    tr = _run(9, 775.0, 0.539)
    ours_trace = PANELS["t41"](tr)
    late_us = float(np.median(ours_trace[tr["t"] > 4.0]))
    ours = float(ours_trace.max()) - late_us

    ratio = ours / theirs
    assert 0.6 < ratio < 1.2, (
        f"T41 overshoot is {ours:.0f} degR against Ballin's {theirs:.0f}, a ratio of "
        f"{ratio:.2f}x. On record is 0.87x after the Eqs. 48-49 rebuild, down from 1.67x "
        f"with the collapsed lead-lag and 3.41x with no heat sink at all. Outside this "
        f"band, re-measure -- and if below, tighten."
    )
    assert pre_them > 0.0


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
