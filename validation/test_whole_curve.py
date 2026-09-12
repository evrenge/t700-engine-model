"""Whole-curve error against Ballin's transients, and the heat-sink configuration test.

## Why this file exists

Until 2026-09-12 the transient comparison measured three points per panel -- the pre-step
plateau, each trace's own extremum, and the settled level. Those are the quantities that
need no time alignment, which is why they were chosen, but they cannot see a wrong
settling time or a wrong path between the points, and they made the agreement look several
times better than it is. Against Figure 9 they reported 1-4 %; the whole curve is 6-11 %.

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
from t700 import realtime, trim
from t700.engine import Ambient, frame
from t700.units import wf_pps_from_pph

REF = Path(__file__).resolve().parent.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)
WF0 = wf_pps_from_pph(400.0)
OUR_STEP = 0.539
BALLIN_STEP = 0.545
"""Midpoint of the bracket [0.539, 0.551]; see the module docstring."""

STEPS = {9: 775.0, 10: 125.0}
PANELS = ("pcng", "ps3", "t41", "t45", "torq45")
UNTRUSTED = {(10, "torq45")}


def _split_strays(x, y):
    """Drop off-curve samples. One rule, matching `plot_validation._split_strays`.

    The trailing-sample and trailing-block rules that used to live here are gone: the
    digitizer drops re-acquired ink at source as of 2026-09-12. See that function.
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
        lambda t: WF0 if t < OUR_STEP else hi,
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


def _rms_pct(t, ours, tb, vb):
    """RMS |ours - Ballin| over the common window, as a percent of Ballin's excursion."""
    shift = BALLIN_STEP - OUR_STEP
    lo, hi = float(tb.min()), min(float(tb.max()), float(t.max()) - shift)
    grid = np.linspace(lo, hi, 800)
    return (
        100.0
        * float(np.sqrt(((np.interp(grid - shift, t, ours) - np.interp(grid, tb, vb)) ** 2).mean()))
        / float(np.ptp(vb))
    )


def test_the_heat_sink_configuration_is_the_better_fit():
    """Figures 9 and 10 were run with the heat sink ON -- measured, not inferred.

    `docs/notes/heat-sink-configuration.md` argues this from four printed facts, the
    strongest being pdf p.47's sentence about Figure 10: "the real-time heat-sink model
    constants were made independent of the direction of power change", which only parses
    if the heat sink is running in the figure under discussion.

    This is the measurement behind that inference, and it is not close: heat sink on wins
    on **every one of the nine comparable panels**, and the mean whole-curve RMS is 9.8 %
    against 18.5 %. Figure 10's PCNG alone goes from 30.6 % to 5.6 %.
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
                errs.append(_rms_pct(t, ours[key], tb, vb))
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


WHOLE_CURVE_RMS_CEILING_PCT = 8.0
"""Ceiling on any single panel's whole-curve RMS, as a percent of its own excursion.

**This is a ratchet, not a tolerance.** Set just above the worst panel measured, so a
regression fails while an improvement is free. Lower it whenever the model improves, and
never raise it. Belongs in `SCOPE.md` per CLAUDE.md.

History, each step a replacement of an invention by something printed:
17.0 (worst panel 15.8 %) -> **8.0** (worst 7.0 %) when the heat sink was rebuilt on
Eqs. 48-49. Mean over the nine panels went 9.78 % -> 3.63 %.
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
    rms = _rms_pct(t, ours[key], tb, vb)
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
        lambda t: WF0 if t < OUR_STEP else wf_pps_from_pph(125.0),
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


def test_the_printed_p45_tolerance_admits_a_false_equilibrium_below_flight_idle():
    """At 125 lbm/hr the printed 0.1 percent P45 criterion settles on a false root.

    Found 2026-09-12 while rebuilding the heat sink, and it is a property of the pressure
    solve rather than of the heat sink. Run the Figure 10 chop out past ~25 s -- twenty
    times Ballin's 4.5 s of record -- and the real-time frame settles at **75.7 %NG with
    T41 1617 degR**, while the differential model trims at **67.0 %NG, T41 1762**. Two
    roots, and the frame picks the wrong one.

    It is the P45 *tolerance*, not the pass cap, and tightening only that fixes it:

    | P45 tol | cap | settles at |
    |---|---|---|
    | 1e-3 (printed) | 8 | 75.745 %NG |
    | 1e-3 (printed) | 20 | 75.745 %NG |
    | 1e-6 | 40 | 66.895 %NG |
    | 1e-10 | 200 | 66.895 %NG |

    **The printed criterion is kept.** 125 lbm/hr is below flight idle, and the report
    states that fuel control below flight-idle power was one of the features eliminated
    from the real-time model [pdf p.38] -- so this condition is outside the envelope the
    criterion was chosen for, and `f9` is being asked for pressure ratios outside its
    table 136,704 times in a 200 s run there. A loose tolerance on a clamped, therefore
    nearly flat, iteration function is exactly how a false fixed point appears. Open
    question #45.

    This test pins the behaviour so it stays known rather than being rediscovered, and
    checks the diagnosis: with P45 converged tightly, the two formulations agree.
    """
    r0 = trim.solve(WF0, c.NP_DES, AMB)
    f0 = frame(r0.state, WF0, AMB)

    def settle(tol: float) -> float:
        st = realtime.from_trim(r0, f0.wa31_pps, f0)
        tr = realtime.run(
            st,
            lambda t: WF0 if t < OUR_STEP else wf_pps_from_pph(125.0),
            AMB,
            duration_s=120.0,
            dt=0.007,
            q_req_ftlbf=f0.q_pt_ftlbf,
            integrate_np=False,
            heat_sink=True,
            tol=tol,
        )
        t = np.asarray(tr["t"])
        return 100.0 * float(np.median(np.asarray(tr["ng"])[t > 110.0])) / c.NG_DES

    trims = [r for r in trim.sweep(list(range(400, 120, -25)) + [125.0]) if r.trustworthy]
    differential = 100.0 * trims[-1].state.ng_rpm / c.NG_DES

    printed = settle(realtime.TOL_PRESSURE)
    assert abs(printed - 75.7) < 0.5, (
        f"under the printed criterion the 125 lbm/hr run settles at {printed:.2f} %NG; "
        f"the false root on record is 75.75 %. If this has moved, re-derive the table "
        f"in this docstring."
    )

    converged = settle(1e-9)
    assert abs(converged - differential) / differential < 0.01, (
        f"with P45 converged the real-time frame settles at {converged:.2f} %NG and the "
        f"differential model trims at {differential:.2f} %; these are the same physics "
        f"and must agree, which is the evidence that the false root is a tolerance "
        f"artifact and not a disagreement between the two formulations"
    )


# --------------------------------------------------------- why the transient cannot be tightened


def test_ps3_has_about_thirteenfold_leverage_on_the_speed_derivative():
    """The structural reason open question #47 cannot be closed from the report.

    `dNG/dt` is the difference of two nearly equal torques, so a small error in station 3
    static pressure is hugely amplified. Ps3 enters twice and both routes pull the same
    way: `T3 = T2*f2(Ps3/P2)` (Eq. 10) into the compressor torque (39), and
    `WA31 = sqrt(P3(P3 - P41)/(K_dpb T3))` (Eq. 18) into the whole core flow -- and Eq. 18
    is the strong one, because P41 tracks P3 closely so a small move in P3 moves the
    pressure *drop* by much more.

    Measured at 76 %NG on the Figure 10 chop: substituting Ballin's plotted Ps3, which
    differs from ours by 2.2 %, changes `dNG/dt` by 29 %. A gain near thirteen.

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
        lambda t: wf0 if t < OUR_STEP else wf_pps_from_pph(125.0),
        AMB,
        duration_s=5.0,
        dt=0.007,
        q_req_ftlbf=f0.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=True,
    )
    pcng = 100.0 * np.asarray(tr["ng"]) / c.NG_DES
    i = int(np.argmin(np.abs(pcng - 76.0)))

    from t700.engine import State

    def dng_with_ps3(ps3: float) -> float:
        s = State(tr["ng"][i], tr["np"][i], ps3 / c.K_PS3, tr["p41"][i], tr["p45"][i])
        return frame(s, tr["wf"][i], AMB, q_req_ftlbf=f0.q_pt_ftlbf, t41_degR=tr["t41"][i]).dng_dt

    ps3 = c.K_PS3 * tr["p3"][i]
    base = dng_with_ps3(ps3)
    perturbed = dng_with_ps3(ps3 * 1.01)
    gain = abs((perturbed - base) / base) / 0.01

    assert 8.0 < gain < 20.0, (
        f"dNG/dt gain on Ps3 measures {gain:.1f} at 76 %NG; on record is about 13. If this "
        f"has moved a long way, the error budget in open question #47 needs redoing."
    )
