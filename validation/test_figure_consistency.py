"""Are the report's own result figures mutually consistent, and are we on all of them?

Trap: a transient (Figures 9-10) must not be compared with an equilibrium sweep (Figures 6-8)
at matched NG. During a chop T41 falls and the choked station 4.1 nozzle passes the same flow
at a lower P41, so P3 sits below its equilibrium value by construction (#46). Valid
comparisons are trace against trace, settled state against the sweep, and the phase plane.

Ballin ties the figures to one condition in text [pdf p.39]: "Gas generator speed is
overestimated by 1 to 2 percent; this is reflected in the trim differences between the
real-time model and the status-81 model in figure 6." See SCOPE.md, "Comparing a transient
to a steady-state figure" and "The figures' own read error".
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import test_fuel_step as tfs
from t700 import constants as c
from t700 import trim
from t700.engine import Ambient, frame
from t700.units import shp_from_torque, wf_pps_from_pph
from test_fuel_step import read_reference

REF = Path(__file__).resolve().parent.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)

PHASE_PLANE_TOL_PCT = 2.5
BOTH_SIDES_TOL_PCT = 0.5
READ_ERROR_CEILING_PCT_FS = 0.6
"""All three ours, declared in `SCOPE.md`; the report states no transient tolerance. The
read-error ceiling is a fraction of each panel's full scale, not of value -- see
`PANEL_SPAN`."""

# Full-scale span of each panel's y axis, from `tools/digitize_fig910.py`'s PANEL_RANGES.
# `calibrate()` there maps pixels to values by a straight line between the two frame rows,
# so a read error is a *pixel* offset: a fixed fraction of the span, not of the value.
PANEL_SPAN = {
    9: dict(wfph=750.0, ps3=200.0, pcng=20.0, t41=1000.0, t45=1000.0, torq45=400.0),
    10: dict(wfph=500.0, ps3=200.0, pcng=40.0, t41=1000.0, t45=1000.0, torq45=400.0),
}

# [pdf p.43, Table 2] NASA-Lewis test conditions, the only non-standard ambient printed;
# it belongs to Table 3. Here so `test_figures_9_and_10_are_sea_level_standard` can rule it out.
T2_WF_PPH = np.array([140.1, 297.2, 372.0, 458.4, 560.6, 694.4])
T2_P2_PSIA = np.array([14.37, 14.17, 14.16, 14.09, 14.02, 13.92])
T2_T2_DEGR = np.array([516.7, 515.6, 508.3, 508.0, 507.2, 507.2])

PRE_STEP_KEYS = ("pcng", "ps3", "t41", "t45", "torq45", "wfph")
"""Panels whose pre-step plateau is measured from the reference files, never transcribed
here: a literal copy goes stale silently when the digitizer changes."""


def _pre_step(fig: int) -> dict[str, float]:
    """Both figures begin at the same 400 lbm/hr trim, so these are two reads of one state."""
    t_step = tfs.STEP_TIME[fig]
    out = {}
    for key in PRE_STEP_KEYS:
        t, v = read_reference(f"fig{fig:02d}_{key}_model.csv")
        out[key] = float(v[t < t_step - 0.05].mean())
    return out


# The WFPH panel is the input, and the caption prints both of its levels.
WFPH_PRINTED = {9: (400.0, 775.0), 10: (400.0, 125.0)}


# --------------------------------------------------------------------- the reference itself


@pytest.mark.parametrize("fig", [9, 10])
def test_the_figures_own_read_error_is_measured_not_estimated(fig: int):
    """The WFPH panel's true value is printed, so digitizing it measures our read error.

    In percent of full scale, not of value: the digitizer maps pixels linearly between two
    frame rows, so the error is a pixel offset. Quoted per value, one offset reads very
    differently on a wide WFPH axis and the narrow PCNG one.
    """
    lo, hi = WFPH_PRINTED[fig]
    t_step = tfs.STEP_TIME[fig]
    span = PANEL_SPAN[fig]["wfph"]
    t, v = read_reference(f"fig{fig:02d}_wfph_model.csv")
    pre = v[t < t_step - 0.05].mean()
    post = v[t > t_step + 0.35].mean()
    for got, want, lbl in ((pre, lo, "pre-step"), (post, hi, "post-step")):
        dev_fs = 100.0 * (got - want) / span
        assert abs(dev_fs) < READ_ERROR_CEILING_PCT_FS, (
            f"fig {fig} WFPH {lbl}: printed {want}, digitized {got:.2f}, "
            f"{dev_fs:+.3f} % of the {span:.0f} full scale "
            f"({100 * (got / want - 1.0):+.2f} % of value) -- the digitizer has drifted; "
            f"this panel's answer is known exactly"
        )


def test_the_two_figures_read_one_state_to_within_a_pixel():
    """Figures 9 and 10 both start at the 400 lbm/hr trim, on axes of different span, so the
    difference between their pre-step plateaus is pure read error, on every panel (#63).
    """
    a, b = _pre_step(9), _pre_step(10)
    for key in PRE_STEP_KEYS:
        span = max(PANEL_SPAN[9][key], PANEL_SPAN[10][key])
        dev_fs = 100.0 * (b[key] - a[key]) / span
        assert abs(dev_fs) < READ_ERROR_CEILING_PCT_FS, (
            f"{key}: fig 9 reads {a[key]:.3f}, fig 10 reads {b[key]:.3f}, "
            f"{dev_fs:+.3f} % of the {span:.0f} full scale. Two reads of one state cannot "
            f"disagree by more than the read error, so either the digitizer has drifted or "
            f"the floor in SCOPE.md is too tight."
        )


def test_ballins_own_cross_reference_between_figures_9_and_6_holds():
    """[pdf p.39] "overestimated by 1 to 2 percent; this is reflected in ... figure 6."

    The trim gap between his model and status-81 is the same quantity in Figure 9 and in
    Figure 6 at the same fuel flow.
    """
    ng9_b, ng9_g = _pre_step(9)["pcng"], 89.759  # Ballin's line, GE '+' markers
    gap_transient = 100.0 * (ng9_b / ng9_g - 1.0)

    wf_b, nb = read_reference("fig06_realtime.csv")
    wf_g, ng = read_reference("fig06_ge_status81.csv")
    b = float(np.interp(400.0, wf_b, nb))
    g = float(np.interp(400.0, wf_g, ng))
    gap_steady = 100.0 * (b / g - 1.0)

    assert 1.0 <= gap_transient <= 2.0, f"Figure 9 trim gap {gap_transient:+.2f} %"
    assert abs(gap_transient - gap_steady) < 0.5, (
        f"Figure 9 gap {gap_transient:+.2f} % against Figure 6 gap {gap_steady:+.2f} % -- "
        f"Ballin ties these together in text, so they should not diverge"
    )


def test_figures_9_and_10_are_sea_level_standard():
    """Their captions do not say so (Figures 6-8's do); Figures 9-10 use the NASA-Lewis
    dynamometer [pdf p.39], so Table 2's Lewis condition [pdf p.43] is a candidate (#53).

    Our standard-day trim fits Ballin's pre-step plateau better than the Lewis condition on
    the four trusted channels (TORQ45 left out). A supporting leg only: the load-bearing
    evidence is Ballin's own cross-reference, checked above.
    """
    p2 = float(np.interp(400.0, T2_WF_PPH, T2_P2_PSIA))
    t2 = float(np.interp(400.0, T2_WF_PPH, T2_T2_DEGR))
    wf0 = wf_pps_from_pph(400.0)

    rms = {}
    for label, amb in (("standard", AMB), ("lewis", Ambient(p2, t2))):
        r = trim.solve(wf0, c.NP_DES, amb)
        f = frame(r.state, wf0, amb)
        ours = dict(
            pcng=100.0 * r.state.ng_rpm / c.NG_DES,
            ps3=f.ps3_psia,
            t41=f.t41_degR,
            t45=f.t45_degR,
        )
        ref9 = _pre_step(9)
        dev = [100.0 * (ours[k] / ref9[k] - 1.0) for k in ours]
        rms[label] = float(np.sqrt(np.mean(np.square(dev))))

    assert rms["standard"] < rms["lewis"], (
        f"sea-level standard rms {rms['standard']:.2f} % against the Table 2 Lewis "
        f"condition {rms['lewis']:.2f} % -- the Lewis condition now fits BETTER, which "
        f"would reopen the question this row closed. The documentary evidence on pdf "
        f"p.39 would then be in conflict with the measurement."
    )
    assert rms["standard"] < 0.6, f"standard-day rms {rms['standard']:.2f} %"


# --------------------------------------------------------------------- transient vs steady


def _ballin_phase_plane(fig: int, t_step: float) -> tuple[np.ndarray, np.ndarray]:
    """Ballin's Ps3 against his NG after the step, NG ascending (Figure 10's falls, so it
    is reversed for interpolation)."""
    tb_n, vb_n = read_reference(f"fig{fig:02d}_pcng_model.csv")
    tb_p, vb_p = read_reference(f"fig{fig:02d}_ps3_model.csv")
    mb = tb_n > t_step
    bn = vb_n[mb]
    bp = np.interp(tb_n[mb], tb_p, vb_p)
    if fig == 10:
        bn, bp = bn[::-1], bp[::-1]
    return bn, bp


@pytest.mark.parametrize(
    "fig,wf_hi,t_step,grid",
    [
        (9, 775.0, tfs.STEP_TIME[9], (92.0, 94.0, 96.0, 98.0)),
        (10, 125.0, tfs.STEP_TIME[10], (78.0, 80.0, 82.0, 84.0, 86.0, 88.0, 90.0)),
    ],
)
def test_the_phase_plane_trajectory_matches(fig, wf_hi, t_step, grid):
    """Ps3 against NG -- the comparison that discards the time axis, and with it the
    unprinted step time (#37).

    Only Figure 10's panel is evidence: Figure 9's Ps3(NG) is nearly a straight chord between
    the two trims (see `test_the_phase_plane_is_not_a_degenerate_comparison`). The metric is
    blind to the volume constants `K_V3` and `K_V41`, which drop out of the algebraic
    pressure solve, and on Figure 9 it is not converged in dt.
    """
    bn, bp = _ballin_phase_plane(fig, t_step)
    tr = tfs._run(wf_hi, t_step, True)
    our_n = 100.0 * tr["ng"] / c.NG_DES
    our_p = c.K_PS3 * tr["p3"]

    mo = tr["t"] > t_step
    on, op = our_n[mo], our_p[mo]
    if fig == 10:  # NG falls, so reverse for interpolation
        on, op = on[::-1], op[::-1]

    for n in grid:
        b = float(np.interp(n, bn, bp))
        o = float(np.interp(n, on, op))
        dev = 100.0 * (o / b - 1.0)
        assert abs(dev) < PHASE_PLANE_TOL_PCT, (
            f"fig {fig} at {n:.0f} %NG: Ps3 ours {o:.2f} against Ballin's {b:.2f}, "
            f"{dev:+.2f} % -- the state trajectory has moved, not just the rate"
        )


def test_the_phase_plane_is_not_a_degenerate_comparison():
    """Figure 10's trajectory is curved; Figure 9's is not. Only the first is evidence: a
    model matching the two endpoint trims already reproduces the straight chord.
    """
    for fig, t_step, grid, straight in (
        (9, tfs.STEP_TIME[9], (92.0, 94.0, 96.0, 98.0), True),
        (10, tfs.STEP_TIME[10], (78.0, 80.0, 82.0, 84.0, 86.0, 88.0, 90.0), False),
    ):
        bn, bp = _ballin_phase_plane(fig, t_step)
        lo, hi = min(grid), max(grid)
        p_lo, p_hi = float(np.interp(lo, bn, bp)), float(np.interp(hi, bn, bp))
        chord = max(
            abs(
                100.0
                * (
                    (p_lo + (n - lo) / (hi - lo) * (p_hi - p_lo)) / float(np.interp(n, bn, bp))
                    - 1.0
                )
            )
            for n in grid
        )
        if straight:
            assert chord < 1.0, (
                f"fig {fig}: the chord now misses Ballin's trace by {chord:.2f} %, so this "
                f"panel has become informative"
            )
        else:
            assert chord > 4.0, (
                f"fig {fig}: the chord misses Ballin's trace by only {chord:.2f} %, so the "
                f"phase-plane agreement here is no longer distinguishable from hitting the "
                f"two endpoint trims -- the conclusion drawn from it does not hold"
            )


@pytest.mark.parametrize(
    "fig,wf_hi,t_step", [(9, 775.0, tfs.STEP_TIME[9]), (10, 125.0, tfs.STEP_TIME[10])]
)
def test_the_transient_sits_off_the_steady_locus_in_the_right_direction(fig, wf_hi, t_step):
    """An accel runs above the equilibrium Ps3 line, a chop below it: the sign of a physical
    effect, so a statement about the model rather than the digitizer.
    """
    ng8, ps8 = read_reference("fig08_realtime.csv")
    tr = tfs._run(wf_hi, t_step, True)
    our_n = 100.0 * tr["ng"] / c.NG_DES
    our_p = c.K_PS3 * tr["p3"]

    # sample mid-transient, where the excursion is largest
    m = (tr["t"] > t_step + 0.1) & (tr["t"] < t_step + 0.4)
    n, p = our_n[m].mean(), our_p[m].mean()
    assert ng8.min() < n < ng8.max()
    excursion = 100.0 * (p / float(np.interp(n, ng8, ps8)) - 1.0)

    if fig == 9:
        assert excursion > 1.0, f"accel should run above the steady locus, got {excursion:+.2f} %"
    else:
        assert excursion < -1.0, f"chop should run below the steady locus, got {excursion:+.2f} %"


def test_we_match_both_sides_of_the_figure_6_to_9_gap():
    """At 775 lbm/hr the sweeps (Figures 6-8) and Figure 9 disagree because Figure 9 has not
    settled when its record ends. We match both at once, and our own unsettled deficit
    accounts for the gap.
    """
    t_end = 4.46  # where Ballin's Figure 9 record stops
    r775 = trim.sweep(np.arange(400.0, 776.0, 25.0), c.NP_DES, AMB)[-1]
    assert r775.trustworthy
    f775 = frame(r775.state, wf_pps_from_pph(775.0), AMB)
    # from the sweep's own first rung: a cold 400 lbm/hr trim, which is what `fuel_step` solves
    tr = tfs.fuel_step(775.0, tfs.STEP_TIME[9], duration_s=40.0)

    wf6, ng6 = read_reference("fig06_realtime.csv")
    wf7, shp7 = read_reference("fig07_realtime.csv")
    ng8, ps8 = read_reference("fig08_realtime.csv")
    ng_f6 = float(np.interp(775.0, wf6, ng6))

    cases = {
        "NG": (
            ng_f6,
            read_reference("fig09_pcng_model.csv")[1][-1],
            100.0 * r775.state.ng_rpm / c.NG_DES,
            float(np.interp(t_end, tr["t"], 100.0 * tr["ng"] / c.NG_DES)),
        ),
        "Ps3": (
            float(np.interp(ng_f6, ng8, ps8)),
            float(np.interp(t_end, *read_reference("fig09_ps3_model.csv"))),
            f775.ps3_psia,
            float(np.interp(t_end, tr["t"], c.K_PS3 * tr["p3"])),
        ),
        "shp": (
            float(np.interp(775.0, wf7, shp7)),
            shp_from_torque(
                float(np.interp(t_end, *read_reference("fig09_torq45_model.csv"))), c.NP_DES
            ),
            shp_from_torque(f775.q_pt_ftlbf, c.NP_DES),
            shp_from_torque(float(np.interp(t_end, tr["t"], tr["q_pt"])), c.NP_DES),
        ),
    }

    for name, (steady, transient, our_steady, our_transient) in cases.items():
        gap = 100.0 * (steady / transient - 1.0)
        assert gap > 0.4, f"{name}: the sweep/transient gap has vanished ({gap:+.2f} %)"
        d_s = 100.0 * (our_steady / steady - 1.0)
        d_t = 100.0 * (our_transient / transient - 1.0)
        assert abs(d_s) < BOTH_SIDES_TOL_PCT, f"{name} vs the sweep: {d_s:+.2f} %"
        assert abs(d_t) < BOTH_SIDES_TOL_PCT, f"{name} vs the transient: {d_t:+.2f} %"
        # and our own unsettled deficit must account for a real share of the gap
        mine = 100.0 * (our_steady / our_transient - 1.0)
        assert 0.85 < mine / gap < 1.10, (
            f"{name}: our unsettled deficit is {mine:+.2f} % of a {gap:+.2f} % gap "
            f"({100 * mine / gap:.0f} %); the settling explanation has stopped holding"
        )
