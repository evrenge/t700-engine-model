"""Are the report's own result figures mutually consistent, and are we on all of them?

This file exists because of a mistake. Comparing Figure 10's transient Ps3 against
Figure 8's steady-state sweep at matched gas generator speed produced 3-6 % deviations,
and those were written up as the report disagreeing with itself. **That comparison is
invalid.** Figures 6-8 are equilibrium sweeps; Figures 9-10 are transients. During a chop
fuel is cut, T41 falls, and the choked station 4.1 nozzle then passes the same flow at a
lower P41 -- so P3 *must* sit below its equilibrium value at that speed. The excursion is
the reason a transient model exists.

Ballin also settles the question in text. Of Figure 9 he writes [pdf p.39]:

    "Gas generator speed is overestimated by 1 to 2 percent; this is reflected in the
    trim differences between the real-time model and the status-81 model in figure 6."

He cross-references the transient figure's trim to the steady figure's, which he could
only do if they share a condition. The tests below check that he was right, and that we
sit on both sides of every apparent disagreement at once -- which is the operative point,
because a single model cannot match two genuinely contradictory datasets.

See `SCOPE.md` -> "Comparing a transient to a steady-state figure" for the tolerances and
for the measured read error of Figures 9-10.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from t700 import constants as c
from t700 import realtime, trim
from t700.engine import Ambient, frame
from t700.units import shp_from_torque, wf_pps_from_pph

REF = Path(__file__).resolve().parent.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)

PHASE_PLANE_TOL_PCT = 2.5
BOTH_SIDES_TOL_PCT = 3.0
SHARED_TRIM_TOL_PCT = 1.5
READ_ERROR_CEILING_PCT = 2.5
"""All four declared in `SCOPE.md`; none is derived from the report, which states no
transient tolerance at all."""

# [pdf p.43, Table 2] NASA-Lewis test conditions. The only non-standard ambient printed
# in the report -- and it belongs to Table 3, not to Figures 9-10. Kept here because
# `test_figures_9_and_10_are_sea_level_standard` has to rule it out with numbers.
T2_WF_PPH = np.array([140.1, 297.2, 372.0, 458.4, 560.6, 694.4])
T2_P2_PSIA = np.array([14.37, 14.17, 14.16, 14.09, 14.02, 13.92])
T2_T2_DEGR = np.array([516.7, 515.6, 508.3, 508.0, 507.2, 507.2])

# Ballin's pre-step plateau, digitized from both figures' own panels. Both figures begin
# at the same 400 lbm/hr trim, so these are two independent reads of one state.
PRE_STEP = {
    9: dict(pcng=90.903, ps3=156.953, t41=2188.33, t45=1563.46, torq45=175.20),
    10: dict(pcng=91.235, ps3=157.543, t41=2193.50, t45=1571.79, torq45=184.11),
}

# The WFPH panel is the input, and the caption prints both of its levels.
WFPH_PRINTED = {9: (400.0, 775.0, 0.539), 10: (400.0, 125.0, 0.545)}


def _trace(name: str) -> tuple[np.ndarray, np.ndarray]:
    path = REF / name
    rows = list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))
    x_key, y_key = list(rows[0])[:2]
    x = np.array([float(r[x_key]) for r in rows])
    y = np.array([float(r[y_key]) for r in rows])
    o = np.argsort(x)
    return x[o], y[o]


def _run(wf_hi_pph: float, t_step: float, duration_s: float = 5.0):
    """Figures 9 and 10: a step from the 400 lbm/hr trim, NP integration suppressed."""
    wf0 = wf_pps_from_pph(400.0)
    r0 = trim.solve(wf0, c.NP_DES, AMB)
    f0 = frame(r0.state, wf0, AMB)
    hi = wf_pps_from_pph(wf_hi_pph)
    st = realtime.from_trim(r0, f0.wa31_pps, f0)
    return realtime.run(
        st,
        lambda t: wf0 if t < t_step else hi,
        AMB,
        duration_s=duration_s,
        dt=0.007,
        q_req_ftlbf=f0.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=True,
    )


# --------------------------------------------------------------------- the reference itself


@pytest.mark.parametrize("fig", [9, 10])
def test_the_figures_own_read_error_is_measured_not_estimated(fig: int):
    """The WFPH panel's true value is printed, so digitizing it measures our read error.

    This is the floor under every other comparison against these two figures. It comes
    out at 1.67-1.83 % at the 400 lbm/hr level -- which is larger than several deviations
    that have been treated as findings in this project.
    """
    lo, hi, t_step = WFPH_PRINTED[fig]
    t, v = _trace(f"fig{fig:02d}_wfph_model.csv")
    pre = v[t < t_step - 0.05].mean()
    post = v[t > t_step + 0.35].mean()
    for got, want, lbl in ((pre, lo, "pre-step"), (post, hi, "post-step")):
        dev = 100.0 * (got / want - 1.0)
        assert abs(dev) < READ_ERROR_CEILING_PCT, (
            f"fig {fig} WFPH {lbl}: printed {want}, digitized {got:.2f}, {dev:+.2f} % -- "
            f"the digitizer has drifted; this panel's answer is known exactly"
        )


def test_figures_9_and_10_share_one_trim_and_agree_on_it():
    """Both figures start at 400 lbm/hr, so their pre-step states must agree.

    `torq45` is excluded: Figure 10's panel yielded 16 markers against 45 everywhere else
    and disagrees with Figure 9's by 5.1 %. That exclusion is recorded in
    `test_fuel_step.py`'s UNTRUSTED and this test is where the 5.1 % is pinned, so the
    exclusion cannot quietly become permanent without evidence.
    """
    worst = 0.0
    for key in ("pcng", "ps3", "t41", "t45"):
        a, b = PRE_STEP[9][key], PRE_STEP[10][key]
        dev = abs(100.0 * (a / b - 1.0))
        worst = max(worst, dev)
        assert dev < SHARED_TRIM_TOL_PCT, f"{key}: fig 9 {a} vs fig 10 {b}, {dev:.2f} %"
    assert worst < SHARED_TRIM_TOL_PCT
    tq = abs(100.0 * (PRE_STEP[9]["torq45"] / PRE_STEP[10]["torq45"] - 1.0))
    assert tq > 2.0, (
        f"fig 9/10 torq45 now agree to {tq:.2f} %; if the digitizer has been fixed, drop "
        f"the UNTRUSTED exclusion in test_fuel_step.py"
    )


def test_ballins_own_cross_reference_between_figures_9_and_6_holds():
    """[pdf p.39] "overestimated by 1 to 2 percent; this is reflected in ... figure 6."

    The trim gap between his model and status-81 should be the same quantity in the
    transient figure and in the steady sweep. It is: +1.27 % in Figure 9, +1.03 % in
    Figure 6 at the same fuel flow, against his stated 1-2 %.
    """
    ng9_b, ng9_g = PRE_STEP[9]["pcng"], 89.759  # Ballin's line, GE '+' markers
    gap_transient = 100.0 * (ng9_b / ng9_g - 1.0)

    wf_b, nb = _trace("fig06_realtime.csv")
    wf_g, ng = _trace("fig06_ge_status81.csv")
    b = float(np.interp(400.0, wf_b, nb))
    g = float(np.interp(400.0, wf_g, ng))
    gap_steady = 100.0 * (b / g - 1.0)

    assert 1.0 <= gap_transient <= 2.0, f"Figure 9 trim gap {gap_transient:+.2f} %"
    assert abs(gap_transient - gap_steady) < 0.5, (
        f"Figure 9 gap {gap_transient:+.2f} % against Figure 6 gap {gap_steady:+.2f} % -- "
        f"Ballin ties these together in text, so they should not diverge"
    )


def test_figures_9_and_10_are_sea_level_standard():
    """Their captions do not say so; Figures 6-8's do. Decided by experiment.

    Table 2 [pdf p.43] is the only non-standard ambient the report prints (P2 13.92-14.37
    psia, T2 507-517 deg R), and Figures 9-10 use "the dynamometer used for testing of the
    NASA-Lewis experimental engine" [pdf p.39] -- so the Lewis test-cell condition is a
    live possibility. It loses: our trim at sea-level standard fits Ballin's pre-step
    plateau roughly four times better.

    That is consistent with Table 2 belonging to Table 3 alone, which also swapped in
    Lewis-derived compressor and turbine functions "in place of the standard functions"
    and is therefore a different model, not merely a different day.
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
        dev = [100.0 * (ours[k] / PRE_STEP[9][k] - 1.0) for k in ours]
        rms[label] = float(np.sqrt(np.mean(np.square(dev))))

    assert rms["standard"] < rms["lewis"] / 2.0, (
        f"sea-level standard rms {rms['standard']:.2f} % against the Table 2 Lewis "
        f"condition {rms['lewis']:.2f} % -- the discrimination has weakened, so the "
        f"ambient behind Figures 9-10 is no longer settled by this evidence"
    )
    assert rms["standard"] < 1.0


# --------------------------------------------------------------------- transient vs steady


@pytest.mark.parametrize(
    "fig,wf_hi,t_step,grid",
    [
        (9, 775.0, 0.539, (92.0, 94.0, 96.0, 98.0)),
        (10, 125.0, 0.545, (78.0, 80.0, 82.0, 84.0, 86.0, 88.0, 90.0)),
    ],
)
def test_the_phase_plane_trajectory_matches(fig, wf_hi, t_step, grid):
    """Ps3 against NG -- the comparison that needs no time axis.

    Every time-domain comparison against these figures carries the unprinted step time
    (open question #17) and a rate comparison on top of the state error. Plotting Ps3
    against NG discards both: what is left is the thermodynamic path the engine takes
    through the state space, which is set by the turbine flow function and the station 4.1
    temperature dynamics. This is the strongest evidence in the project that the gas path
    is right, and it is the metric that localised the remaining Figure 10 error to the
    torque balance rather than to the gas path.
    """
    tb_n, vb_n = _trace(f"fig{fig:02d}_pcng_model.csv")
    tb_p, vb_p = _trace(f"fig{fig:02d}_ps3_model.csv")
    tr = _run(wf_hi, t_step)
    our_n = 100.0 * tr["ng"] / c.NG_DES
    our_p = c.K_PS3 * tr["p3"]

    mb, mo = tb_n > t_step, tr["t"] > t_step
    bn = vb_n[mb]
    bp = np.interp(tb_n[mb], tb_p, vb_p)
    on, op = our_n[mo], our_p[mo]
    if fig == 10:  # NG falls, so reverse for interpolation
        bn, bp, on, op = bn[::-1], bp[::-1], on[::-1], op[::-1]

    for n in grid:
        b = float(np.interp(n, bn, bp))
        o = float(np.interp(n, on, op))
        dev = 100.0 * (o / b - 1.0)
        assert abs(dev) < PHASE_PLANE_TOL_PCT, (
            f"fig {fig} at {n:.0f} %NG: Ps3 ours {o:.2f} against Ballin's {b:.2f}, "
            f"{dev:+.2f} % -- the state trajectory has moved, not just the rate"
        )


@pytest.mark.parametrize("fig,wf_hi,t_step", [(9, 775.0, 0.539), (10, 125.0, 0.545)])
def test_the_transient_sits_off_the_steady_locus_in_the_right_direction(fig, wf_hi, t_step):
    """An accel runs above the equilibrium Ps3 line, a chop below it. Both models do.

    This is the test that would have caught the invalid comparison in the first place. It
    asserts the sign of a physical effect, not a magnitude, so it is a statement about the
    model rather than about the digitizer.
    """
    ng8, ps8 = _trace("fig08_realtime.csv")
    tr = _run(wf_hi, t_step)
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
    """At 775 lbm/hr the sweep and the transient disagree -- because one has not settled.

    Figure 6 reads 99.69 %NG, Figure 9 reads 98.88 at the end of its record; Figure 8 and
    Figure 7 show the same pattern on Ps3 and shp. Ballin's Figure 9 trace is still
    climbing at +0.14 %NG/s when the record ends, so the two figures are not describing
    the same instant.

    The operative fact is that we match **both**, simultaneously, to 1.92 % worst: a
    single model cannot fit two mutually contradictory datasets, so there is no
    contradiction to resolve. Between 39 and 62 % of each gap is reproduced by our own
    model as an unsettled transient; the remainder is inside the 1.8 % read error measured
    above.
    """
    t_end = 4.46  # where Ballin's Figure 9 record stops
    sweep = trim.sweep(np.arange(400.0, 776.0, 25.0), c.NP_DES, AMB)
    r0, r775 = sweep[0], sweep[-1]
    assert r775.trustworthy
    wf0, hi = wf_pps_from_pph(400.0), wf_pps_from_pph(775.0)
    f0 = frame(r0.state, wf0, AMB)
    f775 = frame(r775.state, hi, AMB)
    st = realtime.from_trim(r0, f0.wa31_pps, f0)
    tr = realtime.run(
        st,
        lambda t: wf0 if t < 0.539 else hi,
        AMB,
        duration_s=40.0,
        dt=0.007,
        q_req_ftlbf=f0.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=True,
    )

    wf6, ng6 = _trace("fig06_realtime.csv")
    wf7, shp7 = _trace("fig07_realtime.csv")
    ng8, ps8 = _trace("fig08_realtime.csv")
    ng_f6 = float(np.interp(775.0, wf6, ng6))

    cases = {
        "NG": (
            ng_f6,
            _trace("fig09_pcng_model.csv")[1][-1],
            100.0 * r775.state.ng_rpm / c.NG_DES,
            float(np.interp(t_end, tr["t"], 100.0 * tr["ng"] / c.NG_DES)),
        ),
        "Ps3": (
            float(np.interp(ng_f6, ng8, ps8)),
            float(np.interp(t_end, *_trace("fig09_ps3_model.csv"))),
            f775.ps3_psia,
            float(np.interp(t_end, tr["t"], c.K_PS3 * tr["p3"])),
        ),
        "shp": (
            float(np.interp(775.0, wf7, shp7)),
            shp_from_torque(float(np.interp(t_end, *_trace("fig09_torq45_model.csv"))), c.NP_DES),
            shp_from_torque(f775.q_pt_ftlbf, c.NP_DES),
            shp_from_torque(float(np.interp(t_end, tr["t"], tr["q_pt"])), c.NP_DES),
        ),
    }

    for name, (steady, transient, our_steady, our_transient) in cases.items():
        gap = 100.0 * (steady / transient - 1.0)
        assert gap > 0.5, f"{name}: the sweep/transient gap has vanished ({gap:+.2f} %)"
        d_s = 100.0 * (our_steady / steady - 1.0)
        d_t = 100.0 * (our_transient / transient - 1.0)
        assert abs(d_s) < BOTH_SIDES_TOL_PCT, f"{name} vs the sweep: {d_s:+.2f} %"
        assert abs(d_t) < BOTH_SIDES_TOL_PCT, f"{name} vs the transient: {d_t:+.2f} %"
        # and our own unsettled deficit must account for a real share of the gap
        mine = 100.0 * (our_steady / our_transient - 1.0)
        assert 0.3 < mine / gap < 0.9, (
            f"{name}: our unsettled deficit is {mine:+.2f} % of a {gap:+.2f} % gap "
            f"({100 * mine / gap:.0f} %); the settling explanation has stopped holding"
        )
