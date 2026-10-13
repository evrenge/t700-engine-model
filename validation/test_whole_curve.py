"""Whole-curve error against Ballin's Figure 9/10 transients, and the heat-sink configuration.

Point checks (plateau, extremum, settled level) cannot see a wrong settling time or path, so
this integrates the error over the whole record, normalised by each panel's own excursion so
panels and figures are comparable. Our step and Ballin's are the same measured instant
(`test_fuel_step.STEP_TIME`), so no time shift is applied; the sampling hole at the step is
excluded (see `riser`, open question #66).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import test_fuel_step as tfs
from t700 import constants as c
from t700 import maps, realtime, trim
from t700.engine import Ambient, frame
from t700.units import wf_pps_from_pph
from test_fuel_step import STEP_TIME, read_reference

REF = Path(__file__).resolve().parent.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)
WF0 = wf_pps_from_pph(400.0)
CHOP_STEP = STEP_TIME[10]
"""Figure 10's step instant, ours and Ballin's alike; measured in `test_fuel_step.STEP_TIME`."""

STEPS = {9: 775.0, 10: 125.0}
PANELS = ("pcng", "ps3", "t41", "t45", "torq45")
UNTRUSTED: set[tuple[int, str]] = set()  # see test_fuel_step.UNTRUSTED


def _split_strays(x, y):
    """Drop off-curve samples: more than six pixel quanta from the median of the eight
    neighbours, in a locally flat neighbourhood. The flatness guard keeps a genuine
    near-vertical segment, where a large jump IS the data, out of the count. The digitizer
    drops re-acquired ink at source (open question #51).
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
    return _split_strays(*read_reference(path, "t_s", "value"))


def _run(fig: int, heat_sink: bool):
    return _panels(tfs._run(STEPS[fig], STEP_TIME[fig], heat_sink))


def _panels(tr):
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
panel's own excursion. Declared in `SCOPE.md`. Generous by design: where the riser's two
endpoints fall depends on the tracer, so this checks the edge is the right size, no more."""

RISER_WINDOW_S = 0.30
"""How far either side of the step to look for the riser hole, in seconds."""


def riser(tb: np.ndarray, t_step: float) -> tuple[float, float] | None:
    """The sampling hole at the step: the times either side of it, or None.

    A line tracer cannot follow a vertical edge, so every reference trace has a gap at the
    fuel step and `np.interp` draws a straight chord across it; scoring our edge against that
    chord measures the tracer, not the model (#66). The widest gap near the step, not in the
    whole record: some panels lose the curve against a frame line elsewhere.
    """
    mid = 0.5 * (tb[:-1] + tb[1:])
    near = np.flatnonzero(np.abs(mid - t_step) < RISER_WINDOW_S)
    if near.size == 0:
        return None
    i = int(near[np.argmax(np.diff(tb)[near])])
    return float(tb[i]), float(tb[i + 1])


def _rms_pct(t, ours, tb, vb, fig: int, exclude_riser: bool = True):
    """RMS |ours - Ballin| over the common window, as a percent of Ballin's excursion.

    No time shift: both runs step at the instant measured from his own WFPH panel. The riser
    hole is excluded (see `riser`); `test_the_step_edge_itself_matches` asserts the edge
    separately. `exclude_riser=False` scores the chord too.
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
    """Figures 9 and 10 were run with the heat sink on: it wins on mean whole-curve RMS by a
    wide margin. The documentary argument, from pdf p.47 among others, is in
    `docs/notes/heat-sink-configuration.md`.
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

A ratchet, not a tolerance: set just above the worst panel, so a regression fails and an
improvement is free; lower it when the model improves. Declared in `SCOPE.md`. Current values
per panel are `whole_curve` in `results/report.json`.
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


def test_eq_80s_fixed_point_is_repelling_below_flight_idle_and_the_frame_reaches_it_anyway():
    """Eq. 80's fixed point repels below flight idle, and the chop still settles on the
    differential trim (open question #45).

    Eq. 80 iterates `g(P45) = N / f9(Ps9/P45)` with N constant over the pass, so
    `g'(P45*) = dln(f9) / dln(Ps9/P45)`: the multiplier is f9's elasticity, which crosses -1
    near Ps9/P45 = 0.77. Above that the root repels, and the clamped plateau just beyond f9's
    data (elasticity 0) attracts a false root. Each frame restarts from the previous frame's
    P45, already at the root, so the run does not blow up; where the printed iteration
    fails, `_p45_loop` solves Eq. 80 exactly (`Exit.EXACT`). The settled speed must not
    depend on `tol`, which `step(tol=...)` sets for the P3/P41 and P45 loops alike.
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

    def settle(tol: float) -> tuple[float, float]:
        tr = tfs.fuel_step(125.0, CHOP_STEP, duration_s=120.0, tol=tol)
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
        f"physics and must agree."
    )

    converged, _ = settle(1e-9)
    assert abs(converged - printed) / printed < 1e-3, (
        f"the settled speed must not depend on the tolerance any more: {printed:.3f} %NG "
        f"printed against {converged:.3f} % at 1e-9"
    )


def test_the_sub_idle_root_no_longer_depends_on_the_pass_count():
    """Held at 125 lbm/hr, the settled speed must not depend on the P45 pass cap: dependence
    on its parity is the signature of Eq. 80's iteration not reaching its root. `_p45_loop`
    falls back to the exact solve `_p45_exact` when the printed iteration misses the printed
    criterion.
    """
    original = realtime._p45_loop
    r0 = trim.solve(wf_pps_from_pph(125.0), c.NP_DES, AMB)
    differential = 100.0 * r0.state.ng_rpm / c.NG_DES

    def settle_with_cap(cap: int) -> float:
        def capped(*a, **k):
            return original(*a[:6], k.get("tol", realtime.TOL_PRESSURE), cap)

        realtime._p45_loop = capped
        try:  # held at its own trim: a step from 125 to 125 lbm/hr
            tr = tfs.fuel_step(125.0, 0.0, duration_s=40.0, wf0_pph=125.0)
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
    """Chopping from 400 lbm/hr and holding, the frame map settles on the differential trim
    at every sub-idle flow (#57). Past the f9 elasticity crossing (see above) the printed
    Eq. 80 iteration cannot reach its root, and the exact fallback is what makes them agree.
    """

    # The differential trims come by continuation and are checked `trustworthy`: a cold
    # `solve` can converge to a real root of the clamped system off the map (open question
    # #71, `test_trim_points.py::test_a_solution_off_the_compressor_map_is_not_reported_as_good`).
    ladder = (175.0, 160.0, 150.0, 140.0, 130.0, 125.0)
    subs = trim.sweep(ladder, c.NP_DES, AMB)
    for pph, sub in zip(ladder, subs, strict=True):
        assert sub.trustworthy, (
            f"{pph:.0f} lbm/hr did not trim on the map: converged="
            f"{sub.residual_converged}, on_data={sub.on_data}, "
            f"NGc={sub.frame.ngc_pct:.4g} %. Continuation is supposed to prevent this."
        )
        differential = 100.0 * sub.state.ng_rpm / c.NG_DES
        tr = tfs.fuel_step(pph, CHOP_STEP, duration_s=60.0)
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

    `dNG/dt` is the difference of two nearly equal torques, and Ps3 enters it twice, pulling
    the same way: through T3 (Eq. 10) into the compressor torque (Eq. 39), and through WA31
    (Eq. 18) into the core flow -- the strong route, since P41 tracks P3 and the pressure drop
    moves much more. So 1-2 % on Ps3, as well as a digitized figure can be read, is tens of
    percent on the rate. Pins the absolute sensitivity (rpm/s per psia), not the gain, whose
    denominator is the local deceleration rate and moves with the trajectory.
    """
    f0 = trim.solve(WF0, c.NP_DES, AMB).frame
    tr = tfs._run(125.0, CHOP_STEP, True)  # Figure 10's chop
    from t700.engine import State

    pcng = 100.0 * np.asarray(tr["ng"]) / c.NG_DES

    # rpm/s per psia of Ps3, measured at three speeds the chop passes
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
    """The riser's two endpoints are Ballin's, so the jump between them is a real measurement
    of how far the quantity moved while the fuel changed. Bounded as a percent of the panel's
    excursion, it keeps the fastest part of the transient under test rather than excluded.
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
    """Doubling K_V3, K_V41, K_V45 or K_DAMP leaves every reference point bit-identical: a
    property of the report's formulation, recorded so these four are not taken as covered.

    The real-time frame solves the pressures algebraically [Eqs. 76-80, pdf p.36], so the
    volumes appear only in the continuous model and are pinned by Appendix B alone. `K_DAMP`
    (Eq. 41) multiplies `NP - NP_DES`, identically zero with `integrate_np=False` at NP_DES.
    `f4` (Figure A4) is likewise seen by no comparison against the report (not asserted
    here); only `tools/reproduce_all.sh` ties it to the figure.
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
                # uncached: these runs must see the patched constant
                _t, ours = _panels(tfs.fuel_step(STEPS[fig], STEP_TIME[fig]))
                for key in PANELS:
                    assert np.array_equal(base[fig][key], ours[key]), (
                        f"doubling {name} moved fig {fig} {key}. That is a real change in "
                        f"the model and this test's premise -- that the transient cannot "
                        f"see these four constants -- no longer holds."
                    )
        finally:
            setattr(c, name, original)
