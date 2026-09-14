"""Whole-curve error against Ballin's transients, and the heat-sink configuration test.

## Why this file exists

Until 2026-09-12 the transient comparison measured three points per panel -- the pre-step
plateau, each trace's own extremum, and the settled level. Those are the quantities that
need no time alignment, which is why they were chosen, but they cannot see a wrong
settling time or a wrong path between the points, and they made the agreement look several
times better than it is. Against Figure 9 they reported 1-4 %; the whole curve is 2.4-5.3 %
(it was 6-11 % when this was written, before the heat sink was rebuilt on Eqs. 48-49).

So this integrates the error over the whole record, normalised by the excursion each panel
actually makes so the numbers are comparable between panels and between figures.

## Alignment

Our step fires at 0.539 s. Ballin's is bracketed to [0.539, 0.551] by two independent
readings of his own figure: his Wf trace has no samples between 0.539 and 0.566 -- the
vertical riser, which a line tracer cannot sample -- and his PCNG trace is flat at 0.530
and already rising at 0.551. Errors are computed at the midpoint. The bracket is worth
about 0.6 % of RMS on the worst panel, which is small against the errors themselves.
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
from test_fuel_step import STEP_TIME

REF = Path(__file__).resolve().parent.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)
WF0 = wf_pps_from_pph(400.0)
CHOP_STEP = STEP_TIME[10]
"""Our step and Ballin's are the same instant, per figure, and `test_fuel_step.STEP_TIME`
is where it is measured and explained.

There used to be two numbers and a shift between them -- ours at 0.539 for both figures,
his at 0.545, "midpoint of the bracket [0.539, 0.551]". Both were eyeballed off an
uncorrected time axis, and the residue showed up as our temperatures leading his by 22-49
ms while PCNG and Ps3 led by 2. Stepping at the measured instant removes it: the best-fit
time offset per panel is now -24 to +10 ms with mixed signs, i.e. within three engine
frames and with no systematic lead left to explain."""

STEPS = {9: 775.0, 10: 125.0}
PANELS = ("pcng", "ps3", "t41", "t45", "torq45")
UNTRUSTED: set[tuple[int, str]] = set()  # see test_fuel_step.UNTRUSTED


def _split_strays(x, y):
    """Drop off-curve samples: more than six pixel quanta from the median of the eight
    neighbours, in a locally flat neighbourhood. The flatness guard keeps a genuine
    near-vertical segment, where a large jump IS the data, out of the count.

    The trailing-sample and trailing-block rules that used to live here are gone: the
    digitizer drops re-acquired ink at source as of 2026-09-12 (open question #51), and
    the trailing-block rule had begun condemning the last sixteen genuine samples of
    Figure 10's T41 trace. This is the only copy of the rule; `plot_validation.py` held a
    character-identical second one until it was deleted on 2026-09-14, which is the same
    duplication #51 was raised about.
    """
    dy = np.abs(np.diff(y))
    nz = dy[dy > 0]
    quantum = float(np.median(nz)) if nz.size else 1.0
    bad = np.zeros(len(x), dtype=bool)
    half = 4
    for i in range(half, len(y) - half):
        nb = np.delete(y[i - half : i + half + 1], half)
        if np.ptp(nb) > 3 * quantum:
            continue
        if abs(y[i] - np.median(nb)) > 6 * quantum:
            bad[i] = True
    return x[~bad], y[~bad]


def _ballin(fig: int, key: str):
    path = REF / f"fig{fig:02d}_{key}_model.csv"
    if not path.exists():
        return None, None
    rows = list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))
    t = np.array([float(r["t_s"]) for r in rows])
    v = np.array([float(r["value"]) for r in rows])
    o = np.argsort(t)
    return _split_strays(t[o], v[o])


def _run(fig: int, heat_sink: bool):
    r0 = trim.solve(WF0, c.NP_DES, AMB)
    f0 = frame(r0.state, WF0, AMB)
    hi = wf_pps_from_pph(STEPS[fig])
    st = realtime.from_trim(r0, f0.wa31_pps, f0)
    tr = realtime.run(
        st,
        lambda t: WF0 if t < STEP_TIME[fig] else hi,
        AMB,
        duration_s=5.0,
        dt=0.007,
        q_req_ftlbf=f0.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=heat_sink,
    )
    t = np.asarray(tr["t"])
    return t, {
        "pcng": 100.0 * np.asarray(tr["ng"]) / c.NG_DES,
        "ps3": c.K_PS3 * np.asarray(tr["p3"]),
        "t41": np.asarray(tr["t41"]),
        "t45": np.asarray(tr["t45"]),
        "torq45": np.asarray(tr["q_pt"]),
    }


STEP_JUMP_TOL_PCT = 12.0
"""How far our motion across the riser may differ from Ballin's, as a percent of the
panel's own excursion. Declared in `SCOPE.md`.

Generous by design: the riser's two endpoints are the *last* sample before the edge and the
*first* after it, and where they fall depends on the tracer, so a 3 ms pitch on a curve
moving 600 degR in 90 ms puts tens of degrees of quantisation into the comparison on its
own. It is a check that the edge is the right size, not a precision measurement of it."""

RISER_WINDOW_S = 0.30
"""How far either side of the step to look for the riser hole, in seconds."""


def riser(tb: np.ndarray, t_step: float) -> tuple[float, float] | None:
    """The sampling hole at the step: the times either side of it, or None.

    **A line tracer cannot follow a vertical edge.** Where Ballin's curve moves almost
    instantaneously at the fuel step, the walker in `tools/digitize_fig910.py` finds no
    single thin run per column and emits nothing, so every reference trace has a gap there
    -- 64 to 165 ms wide, against a 3 ms median sampling pitch. `np.interp` then draws a
    **straight chord** across it, and on Figure 9's T41 panel that chord spans 603.97 degR,
    99.9 percent of the panel's entire excursion.

    Comparing our near-vertical edge against that fabricated ramp measures the tracer, not
    the model. It contributed **96.5 %** of the sum of squares on Figure 9's T41 and 93.9 %
    on its T45, which is what made the two best panels in the set look like the two worst.
    Found by the 2026-09-14 validation-quality audit.

    The widest gap *near the step*, not the widest in the record: three panels have a wider
    one elsewhere, where the curve runs flat and the walker lost it against a frame line.
    Panels whose curve is slow enough to trace through the step -- Figure 10's T41 and
    PCNG -- have no hole there, and excluding their widest near-step gap costs 0.1 % and
    0.2 % respectively, so the rule needs no special case.
    """
    mid = 0.5 * (tb[:-1] + tb[1:])
    near = np.flatnonzero(np.abs(mid - t_step) < RISER_WINDOW_S)
    if near.size == 0:
        return None
    i = int(near[np.argmax(np.diff(tb)[near])])
    return float(tb[i]), float(tb[i + 1])


def _rms_pct(t, ours, tb, vb, fig: int, exclude_riser: bool = True):
    """RMS |ours - Ballin| over the common window, as a percent of Ballin's excursion.

    No time shift. There was one, `BALLIN_STEP - OUR_STEP` = 6 ms, because our step time
    and his were two different eyeballed guesses; both runs now step at the instant
    measured from his own WFPH panel, so the traces share an origin by construction.

    The riser hole is excluded -- see `riser`. What that region *does* carry is the size and
    the time of the edge, both readable even though the path between its endpoints is not,
    and `test_the_step_edge_itself_matches` asserts those separately. `exclude_riser=False`
    reproduces the number this file reported before 2026-09-14.
    """
    lo, hi = float(tb.min()), min(float(tb.max()), float(t.max()))
    grid = np.linspace(lo, hi, 800)
    d = np.interp(grid, t, ours) - np.interp(grid, tb, vb)
    if exclude_riser:
        gap = riser(tb, STEP_TIME[fig])
        if gap is not None:
            d = d[(grid < gap[0]) | (grid > gap[1])]
    return 100.0 * float(np.sqrt((d**2).mean())) / float(np.ptp(vb))


def test_the_heat_sink_configuration_is_the_better_fit():
    """Figures 9 and 10 were run with the heat sink ON -- measured, not inferred.

    `docs/notes/heat-sink-configuration.md` argues this from four printed facts, the
    strongest being pdf p.47's sentence about Figure 10: "the real-time heat-sink model
    constants were made independent of the direction of power change", which only parses
    if the heat sink is running in the figure under discussion.

    This is the measurement behind that inference, and it is not close: heat sink on wins
    on **every one of the nine comparable panels**, and the mean whole-curve RMS is
    **3.79 %** against **19.12 %**. Figure 10's PCNG alone goes from 30.55 % to 2.35 %,
    and the worst panel off is Figure 10's T45 at 34.53 %.
    """
    means = {}
    for heat_sink in (True, False):
        errs = []
        for fig in STEPS:
            t, ours = _run(fig, heat_sink)
            for key in PANELS:
                if (fig, key) in UNTRUSTED:
                    continue
                tb, vb = _ballin(fig, key)
                if tb is None:
                    continue
                errs.append(_rms_pct(t, ours[key], tb, vb, fig))
        means[heat_sink] = float(np.mean(errs))
    assert means[True] < means[False], (
        f"heat sink on gives mean whole-curve RMS {means[True]:.2f} % against "
        f"{means[False]:.2f} % off; if that ever reverses, the configuration inference "
        f"in docs/notes/heat-sink-configuration.md is wrong"
    )
    assert means[True] < 0.75 * means[False], (
        f"the margin should be large, not marginal: {means[True]:.2f} % against "
        f"{means[False]:.2f} %"
    )


WHOLE_CURVE_RMS_CEILING_PCT = 4.0
"""Ceiling on any single panel's whole-curve RMS, as a percent of its own excursion.

**This is a ratchet, not a tolerance.** Set just above the worst panel measured, so a
regression fails while an improvement is free. Lower it whenever the model improves.
Declared in `SCOPE.md`, and `tests/test_scope_tolerances.py` checks that the two agree.

Current, over ten panels, with the riser hole excluded (see `riser`):

| panel | Fig. 9 | Fig. 10 |
|---|---|---|
| PCNG | 1.474 | 1.447 |
| PS3 | 0.760 | 1.698 |
| T41 | **0.706** | **3.820** |
| T45 | **0.836** | 2.937 |
| TORQ45 | 0.673 | 2.421 |

Mean **1.677 %**, worst 3.820 on Figure 10's T41. **Every Figure 9 panel is under 1.5 %**
and the best two in the set are its TORQ45 at 0.673 and T41 at 0.706. What is left lives
entirely on Figure 10, the chop.

(This read "Figure 9's two temperature panels are now the best two in the set" until
2026-09-14. They are second and fourth: TORQ45 at 0.673 beats T41's 0.706 and T45's 0.836,
and the table directly above says so. A sentence that contradicts the table it sits under
is the easiest kind of claim to leave standing, which is why the numbers here are now
regenerated by `validation/report.py` rather than typed.)

History, including the one time it went the wrong way and the three times the *measurement*
was the thing that moved rather than the model:

* 17.0 (worst 15.8 %) -> **8.0** (worst 7.62 %) when the heat sink was rebuilt on Eqs.
  48-49. Nine-panel mean 9.78 % -> 3.79 %.
* 8.0 -> **14.0** -> **8.0**, both on 2026-09-13, recorded rather than quietly absorbed
  because the diagnosis that justified the raise turned out to be half right. The P3/P41
  stopping test was corrected to measure the error the report states rather than the
  iterate step; Figure 9 improved on four of five panels and Figure 10's chop degraded to
  t45 13.15 %, attributed to `f1`'s data hole. **It was not `f1`** -- Eq. 80's fixed-point
  iteration does not converge where f9's elasticity is below -1, and solving Eq. 80 where
  the printed iteration fails took t45 back to 7.72 % with Figure 9 untouched.
* **8.0 -> 5.0 on 2026-09-14, model unchanged.** Two defects in our own digitization of
  pdf pp.45-46: the panels are skewed, so every trace carried a ramp of up to 4.2 % of
  panel height (#63); and each panel's time axis was read against one page-wide pair of
  vertical frame columns though those frames drift left going down the page, so the six
  panels' clocks disagreed by up to 38 ms. The second is what produced the long-standing
  appearance of our temperatures leading Ballin's -- best-fit offsets of +2 ms on PCNG,
  +37 on T41, +49 on T45, +54 on TORQ45, a lag growing down the page. After: -24 to +10 ms
  with mixed signs. Figure 10's T45 went 7.238 -> 4.713.
* **5.0 -> 4.0 the same day, model unchanged again.** The metric was drawing a straight
  chord across the reference traces' riser hole and scoring our near-vertical edge against
  it -- 96.5 % of the sum of squares on Figure 9's T41. See `riser`. Excluding it takes the
  mean 2.651 -> 1.677 and moves Figure 9's T41 and T45 off the bottom of the set entirely:
  3.712 -> 0.706 and 3.353 -> 0.836, from the worst two panels to second and fourth best.

Three of those four moves were our measurement and not our model, which is worth stating
plainly: this number was for a long time a worse description of the model than it looked.
"""


@pytest.mark.parametrize("fig", sorted(STEPS))
@pytest.mark.parametrize("key", PANELS)
def test_whole_curve_error_does_not_regress(fig: int, key: str):
    if (fig, key) in UNTRUSTED:
        pytest.skip(f"fig {fig} {key}: reference data untrusted, see test_fuel_step.UNTRUSTED")
    tb, vb = _ballin(fig, key)
    if tb is None:
        pytest.skip("not digitized")
    t, ours = _run(fig, heat_sink=True)
    rms = _rms_pct(t, ours[key], tb, vb, fig)
    assert rms < WHOLE_CURVE_RMS_CEILING_PCT, (
        f"fig {fig} {key}: whole-curve RMS {rms:.2f} % of excursion, past the "
        f"{WHOLE_CURVE_RMS_CEILING_PCT} % ratchet"
    )


def test_figure_10_has_not_settled_by_the_end_of_its_record():
    """A guard against calling a mid-transient value a settled one, which we did.

    Every "settled" number reported for Figure 10 before 2026-09-12 was measured over
    t = 4.0-4.6 s. Neither trace is settled there: Ballin's PCNG is still falling at
    -2.45 %NG/s and ours at -1.46 %NG/s.
    """
    r0 = trim.solve(WF0, c.NP_DES, AMB)
    f0 = frame(r0.state, WF0, AMB)
    st = realtime.from_trim(r0, f0.wa31_pps, f0)
    tr = realtime.run(
        st,
        lambda t: WF0 if t < CHOP_STEP else wf_pps_from_pph(125.0),
        AMB,
        duration_s=6.0,
        dt=0.007,
        q_req_ftlbf=f0.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=True,
    )
    t = np.asarray(tr["t"])
    pcng = 100.0 * np.asarray(tr["ng"]) / c.NG_DES
    window = (t >= 3.5) & (t <= 4.5)
    slope = float(np.polyfit(t[window], pcng[window], 1)[0])
    assert slope < -1.0, (
        f"our Figure 10 run is still decaying at {slope:.2f} %NG/s over 3.5-4.5 s; "
        f"nothing measured in that window may be called a settled value"
    )


def test_eq_80s_fixed_point_is_repelling_below_flight_idle_and_the_frame_reaches_it_anyway():
    """Open question #45, closed 2026-09-13, and its recorded mechanism was wrong twice.

    ## What was recorded

    "At 125 lbm/hr the printed 0.1 percent P45 criterion settles on a false root" -- the
    Figure 10 chop run past 25 s settled at 75.7 %NG where the differential model trims at
    67.0 %, and the explanation on file was *"a loose tolerance on a clamped, therefore
    nearly flat, iteration function is exactly how a false fixed point appears."*

    ## Why the mechanism was backwards

    Eq. 80 iterates `g(P45) = N / f9(Ps9/P45)` with N constant over the pass, so

        g'(P45*) = dln(f9) / dln(Ps9/P45)

    exactly -- the fixed point's multiplier **is** f9's elasticity at the operating point,
    and nothing else enters. Measured on our own digitized f9:

        Ps9/P45   0.60    0.70    0.75    0.77    0.80    0.82    0.849
        dln f9/dln x   -0.253  -0.755  -0.947  -1.062  -1.258  -1.351  -1.568

    It crosses -1 at about **0.77**, so above that ratio the fixed point is **repelling**
    and no tolerance whatsoever reaches it. The settled chop sits at Ps9/P45 = 0.8492 and
    spends 96 % of its frames above 0.77. That point is *inside* f9's data -- the table
    clamps above about 0.86, where the elasticity is exactly 0.000, and that flat clamped
    region is what the old explanation was describing. It is the wrong region: the
    operating point is below it, and the map there is expansive rather than flat.

    ## Why it was also the wrong loop

    The false root is gone **from the chop path**, and correcting P45 is not what removed
    it. (This paragraph said "the false root is gone" without that qualifier for part of
    2026-09-13. Held at its own 125 lbm/hr trim the shipped model still settles at 75.743
    %NG, and which root it reaches depends on the **parity** of `MAX_ITER_P45` and nothing
    else: 75.743 at caps 2, 4, 6, 8, 20 and 67.610 at 3, 5, 7, 9, 21. See
    `test_the_sub_idle_root_is_selected_by_the_parity_of_the_pass_count`.) On 2026-09-13 the
    **P3/P41** stopping test was corrected to measure the error the report states rather
    than the iterate step (see `realtime.TOL_PRESSURE`). With the inner loop converged, the
    125 lbm/hr chop settles at **66.895 %NG against a 67.039 % differential trim, -0.21 %**
    -- and identically at tol 1e-3, 1e-6 and 1e-9, to five figures. The P45 tolerance never
    mattered. The old table in this docstring, which showed 75.745 %NG at the printed
    tolerance and 66.895 at 1e-6, was reading the inner loop's convergence through the P45
    knob it happened to be varying, because `step(tol=...)` sets both.

    ## Why a repelling fixed point does not blow the run up

    Each frame restarts the iteration from the *previous frame's* P45, which is already at
    the fixed point to within rounding. Eight passes amplify that by 1.57^8 ~ 37, which
    leaves it at rounding. Divergence is detectable only on frames whose entering P45 is
    genuinely far off -- 0.4 % of the run, just after the step -- and `_p45_loop` reports
    those as `Exit.DIVERGING` rather than discarding them.

    ## The two explanations are complementary

    Calling the recorded wording "backwards" was itself wrong, and this docstring did that
    for part of 2026-09-13. f9's data ends at Ps9/P45 = 0.85012; the true equilibrium sits
    at 0.84903, inside by 0.00109, elasticity -1.5685 and therefore repelling. The clamped
    plateau immediately beyond has elasticity **0.000** -- strongly attracting -- and that
    is what "a clamped, therefore nearly flat, iteration function invites a false fixed
    point" was describing. The false root sits on it, at 0.8692. The iterate is expelled
    from a repelling root and captured by an attracting plateau 0.13 % away. Both halves
    are real and they are one mechanism.
    """
    f9 = maps.f9()

    def elasticity(x: float, h: float = 1e-6) -> float:
        return (np.log(float(f9(x + h))) - np.log(float(f9(x - h)))) / (
            np.log(x + h) - np.log(x - h)
        )

    assert elasticity(0.70) > -1.0, "f9's elasticity at 0.70 should be inside the unit circle"
    assert elasticity(0.80) < -1.0, "f9's elasticity at 0.80 should be outside it"
    crossing = next(x / 1000 for x in range(700, 860) if elasticity(x / 1000) < -1.0)
    assert 0.75 < crossing < 0.79, (
        f"f9's elasticity crosses -1 at Ps9/P45 = {crossing:.3f}; on record it is about "
        f"0.77, and that crossing is the whole of open question #45's mechanism"
    )

    r0 = trim.solve(WF0, c.NP_DES, AMB)
    f0 = frame(r0.state, WF0, AMB)

    def settle(tol: float) -> tuple[float, float]:
        st = realtime.from_trim(r0, f0.wa31_pps, f0)
        tr = realtime.run(
            st,
            lambda t: WF0 if t < CHOP_STEP else wf_pps_from_pph(125.0),
            AMB,
            duration_s=120.0,
            dt=0.007,
            q_req_ftlbf=f0.q_pt_ftlbf,
            integrate_np=False,
            heat_sink=True,
            tol=tol,
        )
        t = np.asarray(tr["t"])
        late = t > 110.0
        ng = 100.0 * float(np.median(np.asarray(tr["ng"])[late])) / c.NG_DES
        ratio = AMB.p_amb_psia / float(np.median(np.asarray(tr["p45"])[late]))
        return ng, ratio

    trims = [r for r in trim.sweep(list(range(400, 120, -25)) + [125.0]) if r.trustworthy]
    differential = 100.0 * trims[-1].state.ng_rpm / c.NG_DES

    printed, ratio = settle(realtime.TOL_PRESSURE)
    assert ratio > 0.77, (
        f"the settled chop sits at Ps9/P45 = {ratio:.4f}, below the -1 elasticity "
        f"crossing; the repelling regime this test is about is no longer being entered"
    )
    assert abs(printed - differential) / differential < 0.01, (
        f"at the printed criterion the 125 lbm/hr run settles at {printed:.3f} %NG and "
        f"the differential model trims at {differential:.3f} %. These are the same "
        f"physics and must agree. They did not until the P3/P41 loop was made to stop on "
        f"the error rather than the step, which is what removed the 75.7 % false root."
    )

    converged, _ = settle(1e-9)
    assert abs(converged - printed) / printed < 1e-3, (
        f"the settled speed must not depend on the tolerance any more: {printed:.3f} %NG "
        f"printed against {converged:.3f} % at 1e-9"
    )


def test_the_sub_idle_root_no_longer_depends_on_the_pass_count():
    """It used to depend on the *parity* of `MAX_ITER_P45` and nothing else.

    Held at its own 125 lbm/hr trim, the model settled at 75.743 %NG for an even cap and
    67.610 for an odd one, against a 67.039 differential trim -- magnitude irrelevant out
    to 21 passes. At 150 lbm/hr the even branch sat 7.07 % low. That was Eq. 80's
    fixed-point iteration failing to reach its root wherever f9's elasticity is below -1,
    and landing wherever the parity put it.

    `_p45_loop` now falls back to `_p45_bisect` when the printed iteration does not meet
    the printed criterion, so the root is found rather than guessed at. Every cap from 2 to
    21 gives the same answer, to the digit.
    """
    original = realtime._p45_loop
    r0 = trim.solve(wf_pps_from_pph(125.0), c.NP_DES, AMB)
    f0 = frame(r0.state, wf_pps_from_pph(125.0), AMB)
    differential = 100.0 * r0.state.ng_rpm / c.NG_DES

    def settle_with_cap(cap: int) -> float:
        def capped(*a, **k):
            return original(*a[:6], k.get("tol", realtime.TOL_PRESSURE), cap)

        realtime._p45_loop = capped
        try:
            st = realtime.from_trim(r0, f0.wa31_pps, f0)
            tr = realtime.run(
                st,
                lambda t: wf_pps_from_pph(125.0),
                AMB,
                duration_s=40.0,
                dt=0.007,
                q_req_ftlbf=f0.q_pt_ftlbf,
                integrate_np=False,
                heat_sink=True,
            )
        finally:
            realtime._p45_loop = original
        t = np.asarray(tr["t"])
        pcng = 100.0 * np.asarray(tr["ng"]) / c.NG_DES
        return float(np.median(pcng[t > 35.0]))

    settled = {cap: settle_with_cap(cap) for cap in (2, 3, 4, 5, 8, 9, 20, 21)}
    spread = max(settled.values()) - min(settled.values())
    assert spread < 0.01, (
        f"the settled speed still moves with the pass cap: {settled}. Spread {spread:.4f} "
        f"%NG. Parity dependence is the signature of Eq. 80 not reaching its root."
    )
    for cap, got in settled.items():
        assert abs(got - differential) / differential < 0.001, (
            f"cap {cap} settles at {got:.3f} %NG against a differential trim of "
            f"{differential:.3f}. These are the same physics and must agree."
        )


def test_the_frame_map_agrees_with_the_differential_model_below_flight_idle():
    """It did not, by up to 7.07 %NG, and the cause was Eq. 80's iteration alone.

    On record before 2026-09-13, chopping from 400 lbm/hr and holding:

        lbm/hr   differential   frame settled   deviation
           175       78.509 %        78.509 %      +0.00 %
           150       74.009          68.781        -7.07 %
           140       71.727          68.153        -4.98 %
           130       68.897          67.520        -2.00 %
           125       67.039          66.895        -0.21 %

    Identical at tol 1e-3 and 1e-9 and dead still at the end, so not a convergence
    artifact in the outer loop -- two formulations of the same physics settling in
    different places.

    **All of it was Eq. 80.** Established constructively: replacing its fixed-point
    iteration with a bracketed root-find on the same residual -- same equation, same f9
    data, convergence guaranteed by the bracket rather than by contraction -- put every one
    of those flows on the differential trim to 0.00 %. `_p45_loop` now does exactly that
    when the printed iteration fails to meet the printed criterion, and the table above
    reads +0.003 % or better everywhere.

    The boundary is f9's elasticity crossing. At 175 lbm/hr Ps9/P45 = 0.727, below the
    -1 crossing at ~0.77, so the printed iteration converges on its own in 2.0 passes and
    the fallback never fires. Below about 160 lbm/hr it fires every frame.
    """
    r0 = trim.solve(WF0, c.NP_DES, AMB)
    f0 = frame(r0.state, WF0, AMB)

    for pph in (175.0, 160.0, 150.0, 140.0, 130.0, 125.0):
        sub = trim.solve(wf_pps_from_pph(pph), c.NP_DES, AMB)
        differential = 100.0 * sub.state.ng_rpm / c.NG_DES
        st = realtime.from_trim(r0, f0.wa31_pps, f0)
        tr = realtime.run(
            st,
            lambda t, g=pph: WF0 if t < CHOP_STEP else wf_pps_from_pph(g),
            AMB,
            duration_s=60.0,
            dt=0.007,
            q_req_ftlbf=f0.q_pt_ftlbf,
            integrate_np=False,
            heat_sink=True,
        )
        t = np.asarray(tr["t"])
        late = (100.0 * np.asarray(tr["ng"]) / c.NG_DES)[t > 55.0]
        assert late.max() - late.min() < 0.01, (
            f"{pph:.0f} lbm/hr has not settled by 55 s (spread {late.max() - late.min():.4f})"
        )
        dev = 100.0 * (float(np.median(late)) - differential) / differential
        assert abs(dev) < 0.01, (
            f"{pph:.0f} lbm/hr: the frame settles {dev:+.5f} % from the differential trim "
            f"of {differential:.3f} %NG. The two formulations are the same physics; a "
            f"deviation here means Eq. 80 is not reaching its root again."
        )


# --------------------------------------------------------- why the transient cannot be tightened


def test_ps3_leverage_on_the_speed_derivative_is_what_it_was_measured_to_be():
    """The structural reason open question #47 cannot be closed from the report.

    `dNG/dt` is the difference of two nearly equal torques, so a small error in station 3
    static pressure is hugely amplified. Ps3 enters twice and both routes pull the same
    way: `T3 = T2*f2(Ps3/P2)` (Eq. 10) into the compressor torque (39), and
    `WA31 = sqrt(P3(P3 - P41)/(K_dpb T3))` (Eq. 18) into the whole core flow -- and Eq. 18
    is the strong one, because P41 tracks P3 closely so a small move in P3 moves the
    pressure *drop* by much more.

    **What is invariant here is the sensitivity, not the gain**, and this test pinned the
    wrong one until 2026-09-13. It asserted a *ratio* -- the percent change in `dNG/dt`
    per percent change in Ps3 -- whose denominator is the local deceleration rate, which
    depends entirely on where on the trajectory 76 %NG happens to fall. Correcting the
    P3/P41 stopping rule moved that: at 76 %NG the chop now runs at -2948 rpm/s where the
    under-converged one ran at -880, so the same physics reports a gain of 4.1 instead of
    13.8. Nothing about the amplification changed.

    The absolute sensitivity did not move at all. Measured on the Figure 10 chop, before
    and after that correction:

        76 %NG   -170.0  ->  -169.9 rpm/s per psia of Ps3
        80 %NG   -163.5  ->  -163.4
        86 %NG   -129.9  ->  -133.9

    -- 0.1 % on the first two, across a change that moved the trajectory's floor by 4.3
    percentage points of NG. That is the quantity to pin, and this test pins it.

    The amplification is then that sensitivity divided by whatever `dNG/dt` is locally,
    and it is large wherever the chop is slow: 13.8 at 76 %NG on the old trajectory,
    10.5 at 80 %NG on the new one. Either way, 1 to 2 % on Ps3 -- as well as a digitized
    figure can be read -- is tens of percent on the rate.

    That is why the remaining deceleration disagreement is not localisable. Every input
    has been verified against Ballin's own printed data -- WA31 to 0.8 % with no map in
    the loop, Ps3 to 1.7 % dynamically and 1.0 % against Figure 8, `f2` digitized to
    0.025 % per point with knots every 1.0 in pressure ratio, `f3` confirmed against
    Figure A3, Eq. 39 shown to be an exact two-stream balance -- and 1 to 2 % on Ps3,
    which is as well as a digitized figure can be read, *is* 15 to 30 % on the derivative.

    Constructive proof of the attribution: pinning Ps3 to Ballin's own trace and
    integrating removes the low-speed deficit entirely, taking the rate ratio at 76 %NG
    from 1.29 to 1.00. It also breaks the high-speed end, 1.01 to 0.87 at 86 %NG, because
    there our computed Ps3 is the more accurate of the two. So Ps3 accounts for the whole
    residual in both directions, and neither reading is good enough to do better.

    This test pins the leverage itself, so the argument stays measured.
    """
    wf0 = wf_pps_from_pph(400.0)
    r0 = trim.solve(wf0, c.NP_DES, AMB)
    f0 = frame(r0.state, wf0, AMB)
    st = realtime.from_trim(r0, f0.wa31_pps, f0)
    tr = realtime.run(
        st,
        lambda t: wf0 if t < CHOP_STEP else wf_pps_from_pph(125.0),
        AMB,
        duration_s=5.0,
        dt=0.007,
        q_req_ftlbf=f0.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=True,
    )
    from t700.engine import State

    pcng = 100.0 * np.asarray(tr["ng"]) / c.NG_DES

    # rpm/s per psia of Ps3, measured on record at the three speeds the chop passes
    # 86 %NG moved -133.9 -> -156.9 with the f1 interpolation change of 2026-09-13; the
    # two lower speeds barely moved, which is consistent with the change being in how the
    # speed lines are blended rather than in the physics.
    ON_RECORD = {76.0: -171.7, 80.0: -163.5, 86.0: -156.9}

    for target, expected in ON_RECORD.items():
        i = int(np.argmin(np.abs(pcng - target)))
        assert abs(pcng[i] - target) < 0.5, (
            f"the chop no longer passes {target:.0f} %NG (floor {pcng.min():.2f} %); the "
            f"trajectory has changed shape and this measurement needs redoing"
        )

        def dng_with_ps3(ps3: float, i: int = i) -> float:
            s = State(tr["ng"][i], tr["np"][i], ps3 / c.K_PS3, tr["p41"][i], tr["p45"][i])
            fr = frame(s, tr["wf"][i], AMB, q_req_ftlbf=f0.q_pt_ftlbf, t41_degR=tr["t41"][i])
            return fr.dng_dt

        ps3 = c.K_PS3 * tr["p3"][i]
        sensitivity = (dng_with_ps3(ps3 * 1.01) - dng_with_ps3(ps3)) / (0.01 * ps3)
        assert abs(sensitivity - expected) / abs(expected) < 0.05, (
            f"d(dNG/dt)/dPs3 at {target:.0f} %NG measures {sensitivity:.1f} rpm/s/psia "
            f"against {expected:.1f} on record. This is the invariant behind open "
            f"question #47's error budget; if it has moved, the budget needs redoing."
        )


@pytest.mark.parametrize("fig", sorted(STEPS))
@pytest.mark.parametrize("key", PANELS)
def test_the_step_edge_itself_matches(fig: int, key: str):
    """What the riser region does carry: the size of the edge, and when it happened.

    `_rms_pct` excludes the hole because the *path* across it is drawn by `np.interp` and
    not by Ballin. Its two endpoints are his, though, so the jump between them is a real
    measurement of how far the quantity moved while the fuel changed, and the midpoint is a
    real measurement of when. Asserting those keeps the fastest part of the transient under
    test rather than merely excluded.

    The bound is on the *jump*, normalised by the panel's excursion, because that is the
    quantity the riser measures; the edge time is bounded at 40 ms, which is the widest
    riser (165 ms on Figure 9's TORQ45) divided by four and about six engine frames.
    """
    tb, vb = _ballin(fig, key)
    gap = riser(tb, STEP_TIME[fig])
    if gap is None:
        pytest.skip(f"fig {fig} {key}: no sampling gap near the step")
    t0, t1 = gap
    i = int(np.flatnonzero(tb == t0)[0])
    theirs = vb[i + 1] - vb[i]
    if abs(theirs) < 0.05 * np.ptp(vb):
        pytest.skip(f"fig {fig} {key}: no riser, the tracer followed the edge ({theirs:.2f})")

    t, ours = _run(fig, True)
    o = ours[key]
    mine = float(np.interp(t1, t, o) - np.interp(t0, t, o))
    dev = 100.0 * (mine - theirs) / np.ptp(vb)
    assert abs(dev) < STEP_JUMP_TOL_PCT, (
        f"fig {fig} {key}: across the riser Ballin moves {theirs:.2f} and we move "
        f"{mine:.2f}, {dev:+.2f} % of the panel's {np.ptp(vb):.2f} excursion"
    )


def test_the_control_volumes_cannot_be_constrained_by_any_transient_panel():
    """Doubling K_V3, K_V41, K_V45 or K_DAMP leaves every reference point bit-identical.

    Recorded rather than fixed, because it is a property of the report's own formulation
    and not a gap we can close from Figures 9 and 10.

    **The volumes.** The real-time frame solves the three pressures algebraically [Eqs.
    76-80] -- that is what the quasi-steady reduction on pdf p.36 *is*. The control volumes
    are the coefficients of the pressure derivatives, and the derivatives are what got
    removed, so `K_V3`, `K_V41` and `K_V45` appear nowhere in `realtime.py`. They survive
    only in `engine.frame`, i.e. in the continuous 5-DOF model, which is what Appendix B
    compares against. Those three constants are therefore pinned by the Appendix B
    element comparison and by nothing else in the project.

    **The damping.** `K_DAMP` (Eq. 41) multiplies `NP - NP_DES`, and Figures 9 and 10 are
    run with `integrate_np=False` at NP = NP_DES, so that term is identically zero on every
    transient frame. At the three Table B.1 trims NP is 20895 against NP_DES 20900, so it
    contributes 0.036 ft*lbf -- 0.016 to 0.048 % of Q_PT -- and `K_DAMP` would have to be
    wrong by about 60x before the +/-3 % torque comparison noticed.

    Found by the 2026-09-14 validation-quality audit, which mutated 30 constants and
    measured which tests noticed. This test exists so the four that cannot be noticed here
    are on record as such, rather than appearing to be covered.
    """
    base = {}
    for fig in STEPS:
        t, ours = _run(fig, True)
        base[fig] = {k: v.copy() for k, v in ours.items()}

    for name in ("K_V3", "K_V41", "K_V45", "K_DAMP"):
        original = getattr(c, name)
        try:
            setattr(c, name, original * 2.0)
            for fig in STEPS:
                _t, ours = _run(fig, True)
                for key in PANELS:
                    assert np.array_equal(base[fig][key], ours[key]), (
                        f"doubling {name} moved fig {fig} {key}. That is a real change in "
                        f"the model and this test's premise -- that the transient cannot "
                        f"see these four constants -- no longer holds."
                    )
        finally:
            setattr(c, name, original)
