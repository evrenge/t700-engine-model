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
time is not printed; measured at 0.5232 s and 0.5217 s -- open question #37).

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
STEP_TIME = {9: 0.5232, 10: 0.5217}
"""When the fuel steps, seconds, **measured from the WFPH panel's own 50 % crossing**.

The caption states both flow levels and not the instant [pdf pp.45-46], so this is read off
the figure -- open question #37, closed 2026-09-14 by this measurement.

It read 0.539 and 0.545 until then, eyeballed, and the two figures were assumed to differ.
They do not: measured on the deskewed panels with each panel's own time axis (#63), Figure
9 steps at **0.5232 s** and Figure 10 at **0.5217 s**, agreeing to 1.5 ms, as one test
setup should. The 16-23 ms correction is not a nicety on a transient sampled every 7 ms --
it was the whole of the apparent lag between our temperatures and Ballin's. See
`test_whole_curve`."""

STEPS = [(9, 775.0, STEP_TIME[9]), (10, 125.0, STEP_TIME[10])]

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


UNTRUSTED: set[tuple[int, str]] = set()
"""Panels excluded from the comparison, with the evidence. **Currently none.**

Both figures begin at the same 400 lbm/hr trim, so their t = 0 values must agree with each
other regardless of any model. `(10, "torq45")` sat here from 2026-09-12 to 2026-09-14 on
the evidence that it disagreed with Figure 9's by **5.10 %** while the other five panels
agreed to 0.62 %, so "the reference data is wrong there, not the model".

**That was our own page skew.** pdf p.46 is scanned at a slight rotation and every panel
frame on it descends 11-14 px across the panel width; the digitizer mapped rows to values
through the frame positions at the panel's *left* edge only, which laid a monotone ramp of
up to 4.2 % of panel height onto every trace. TORQ45 has the smallest data excursion of the
six panels relative to its axis, so it showed the artifact most. With the frames fitted
along their length the two figures agree on TORQ45 to **0.93 %**, and the panel is back in
the comparison.

What remains true of it is the marker count -- 13 slots against 45 on every other panel --
but that is the GE status-81 reference series, not the model trace this file compares
against, and it is recorded in the CSV header rather than here.
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
    tr = _run(9, 775.0, STEP_TIME[9])
    ngc = 100.0 * tr["ng"] / c.NG_DES
    assert ngc.min() > 88.0, "the step up should not approach the map's lower edge"
    assert ngc.max() < 101.0
    rep = maps.clamp_report()
    assert not any(k.startswith("f1@") and v > 50 for k, v in rep.items()), (
        f"the compressor map should barely be extrapolated on the step up: {rep}"
    )
    maps.reset_clamps()


def test_step_down_runs_off_the_bottom_of_the_maps():
    """Figure 10 still reads `f1` outside its data, and the floor is now Ballin's.

    ## The history, because it is the point of this test

    Four values of the NGc minimum have been on record, and the differences are entirely in
    how the two pressure solves were converged:

    | pressure solves | NGc bottoms at |
    |---|---|
    | inner `tol=1e-10` on the step, cap 40 -- invented, pre-2026-09-12 | 70.94 % |
    | inner `tol=1e-3` on the **step** -- the printed number, wrong quantity | 74.64 % |
    | inner on the **error**, Eq. 80 as printed | 69.92 % |
    | inner on the error, **Eq. 80 solved** when its iteration fails | **74.09 %** |

    Against Ballin's own 74.18 %. The second row looked like a triumph and was not: testing
    the iterate step rather than the error under-converges the inner loop by a factor of
    `rho/(1-rho)` ~ 8, and the lag damped the plunge into near-agreement. The third row is
    what the inner loop honestly gives with Eq. 80 left as printed. The fourth is what both
    solves converged gives, and it lands 0.09 %NG from Ballin.

    So the deep plunge was **Eq. 80's fixed-point iteration**, which does not converge where
    f9's elasticity is below -1, and not the compressor map. That matters because the third
    row was attributed to `f1` at the time, and on that attribution the whole-curve ratchet
    was raised from 8 to 14 %. It is back at 8.

    ## And the `f1` extrapolation was never real either

    With `f1` interpolated at constant abscissa, this run clamped the 65 % speed line **389
    times** -- every frame below 80 %NGc, asking it for pressure ratios up to 5.44 against
    its own last knot at 3.753 -- and that was recorded as open question #58, "f1's 15-point
    data hole".

    Interpolating along Figure A1's own printed construction lines instead (see
    `maps.SpeedMap`) takes that to **zero**. Each speed line ends at its own surge limit, so
    blending the two lines knot by knot blends their limits too, and the blended line's
    right edge runs 4.726 at 70 %NGc and 5.504 at 74 -- comfortably outside the chop's worst
    query of 5.439. The model never leaves the compressor map. The extrapolation was an
    artifact of the evaluation, not a property of the data.

    What is left is `f8` and `f9`, both near the top of their own ranges, and the fuel-air
    ratio falling to 0.00527 -- which is no longer an extrapolation of anything either,
    since `f6` is the constant Figure A6 draws and a constant has no domain to leave.

    The fuel-air ratio still falls to 0.00527. That is no longer an extrapolation of
    anything: `f6` is loaded as the constant Figure A6 draws, and a constant has no domain
    to leave. See `maps.f6`.

    Figure 10's reference data is also the weaker of the two: its `WFPH` settles 4.35 %
    from the value its caption states, against 0.23 % for Figure 9, and its `TORQ45`
    panel yielded only 16 markers against 45 elsewhere and is not used.
    """
    maps.reset_clamps()
    tr = _run(10, 125.0, STEP_TIME[10])
    ngc = 100.0 * tr["ng"] / c.NG_DES
    rep = maps.clamp_report()

    # Ours bottoms at 74.09 %, Ballin's at 74.18. The window is tight on both sides: a
    # floor that fell back toward 70 % would mean a pressure solve has stopped converging,
    # and one that rose would mean the fuel cut is no longer being followed.
    assert 73.0 < ngc.min() < 75.0, (
        f"NGc bottoms at {ngc.min():.2f} %, against Ballin's 74.18 % and our 74.09 % on "
        f"record. If this has fallen, check that both pressure solves still converge "
        f"before believing the model changed physically."
    )
    assert tr["far"].min() < 0.010, "the fuel-air ratio still falls below f6's plotted range"
    assert "f1@65" not in rep and "f1" not in rep, (
        f"the chop is reading f1 outside its data again: {rep}. Under constant-k "
        f"interpolation the blended line carries both speed lines' surge limits, so a "
        f"query inside both is inside the blend and nothing should clamp."
    )
    assert "f6" not in rep, "f6 is a constant now and cannot be clamped -- see maps.f6"
    maps.reset_clamps()


def test_the_two_figures_step_in_opposite_directions():
    """Guards against the two runs being accidentally identical."""
    up = _run(9, 775.0, STEP_TIME[9])
    down = _run(10, 125.0, STEP_TIME[10])
    assert up["ng"][-1] > up["ng"][0]
    assert down["ng"][-1] < down["ng"][0]


@pytest.mark.parametrize("key", list(PANELS) + ["wfph"])
def test_the_two_figures_agree_at_their_shared_trim(key: str):
    """A model-free check on the reference data itself.

    Figures 9 and 10 both start from the same 400 lbm/hr trim, so their t = 0 values must
    agree with each other whatever our model says. This caught Figure 10's TORQ45 panel
    disagreeing with Figure 9's by 5.10 % while every other panel agreed to 0.62 %, which
    put that panel on `UNTRUSTED` for two days. The cause turned out to be the page skew the
    digitizer corrects as of 2026-09-14, and all six panels now agree: worst 0.93 % on
    TORQ45, 0.25 % on the other five.
    """
    vals = []
    for fig in (9, 10):
        path = REF / f"fig{fig:02d}_{key}_model.csv"
        if not path.exists():
            pytest.skip(f"{path.name} not digitized")
        t, v = _reference(fig, key)
        vals.append(v[t < 0.48].mean())
    dev = (vals[1] - vals[0]) / abs(vals[0]) * 100.0
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

    The overshoot went 3.17x Ballin's without the heat sink, 1.67x with it and the old
    realization, and is now BELOW his at **0.83x**. (3.41x and 0.87x were on record until
    2026-09-13; both moved with that day's stopping-rule correction and neither was swept
    up. Measured today with this file's own `_t41_overshoot` definition: ours 100.1 degR
    against Ballin's 120.8 with the heat sink, 383.4 without.)

    The band is a ratchet: leaving it in either
    direction means re-measuring, and the lower bound exists so an improvement is noticed
    rather than silently absorbed.
    """
    tb, vb = _reference(9, "t41")
    pre_them = float(np.median(vb[tb < 0.45]))
    late_them = float(np.median(vb[(tb > 3.5) & (tb < 4.5)]))
    theirs = float(vb.max()) - late_them

    tr = _run(9, 775.0, STEP_TIME[9])
    ours_trace = PANELS["t41"](tr)
    late_us = float(np.median(ours_trace[tr["t"] > 4.0]))
    ours = float(ours_trace.max()) - late_us

    ratio = ours / theirs
    assert 0.6 < ratio < 1.2, (
        f"T41 overshoot is {ours:.0f} degR against Ballin's {theirs:.0f}, a ratio of "
        f"{ratio:.2f}x. On record is 0.83x after the Eqs. 48-49 rebuild, down from 1.67x "
        f"with the collapsed lead-lag and 3.17x with no heat sink at all. Outside this "
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
    bare = _run(9, 775.0, STEP_TIME[9], heat_sink=False)
    sunk = _run(9, 775.0, STEP_TIME[9], heat_sink=True)
    off, _ = _t41_overshoot(bare["t"], bare["t41"], STEP_TIME[9])
    on, _ = _t41_overshoot(sunk["t"], sunk["t41"], STEP_TIME[9])
    assert on < 0.6 * off, (
        f"the heat sink should roughly halve the T41 overshoot; it went {off:.0f} -> {on:.0f} degR"
    )
