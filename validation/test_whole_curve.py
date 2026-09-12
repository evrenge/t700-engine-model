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
    """Drop off-curve samples. Same three rules as `plot_validation._split_strays`."""
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
    dt = np.diff(x)
    if dt.size and dt[-1] > 10 * np.median(dt):
        bad[-1] = True
    late = x > 0.7 * x.max()
    keep = late & ~bad
    if keep.sum() > 10:
        spread = abs(np.subtract(*np.percentile(y[keep], [75, 25])))
        if spread < 10 * quantum:
            plateau = float(np.median(y[keep]))
            k = len(y) - 1
            while k >= 0 and abs(y[k] - plateau) > 6 * quantum:
                bad[k] = True
                k -= 1
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


WHOLE_CURVE_RMS_CEILING_PCT = 17.0
"""Ceiling on any single panel's whole-curve RMS, as a percent of its own excursion.

**This is a ratchet, not a tolerance.** It is set just above the worst panel measured
(Figure 10's T45 at 15.8 %) so that a regression fails while an improvement is free. It
does not mean 17 % is acceptable -- it is not; see open question #47. Lower it whenever
the model improves, and never raise it. Belongs in `SCOPE.md` per CLAUDE.md.
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
    -2.45 %NG/s and ours at -1.46 %NG/s. Run out to 60 s, ours reaches 66.89 %NG, which
    is our own differential trim's 67.04 % -- so the two formulations do agree, and the
    apparent 12 % gap between them was this mistake, not a defect.
    """
    r0 = trim.solve(WF0, c.NP_DES, AMB)
    f0 = frame(r0.state, WF0, AMB)
    st = realtime.from_trim(r0, f0.wa31_pps, f0)
    tr = realtime.run(
        st,
        lambda t: WF0 if t < OUR_STEP else wf_pps_from_pph(125.0),
        AMB,
        duration_s=60.0,
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

    late = t > 50.0
    settled = float(np.median(pcng[late]))
    trims = [r for r in trim.sweep(list(range(400, 120, -25)) + [125.0]) if r.trustworthy]
    differential = 100.0 * trims[-1].state.ng_rpm / c.NG_DES
    assert abs(settled - differential) / differential < 0.01, (
        f"the real-time frame settles at {settled:.2f} %NG and the differential model "
        f"trims at {differential:.2f} %; these are the same physics and must agree"
    )
