"""Compare the real-time model against Ballin's published fuel-step transients.

Figures 9 and 10 [pdf pp.45-46] are the report's open-loop transient evidence: a step up
from 400 to 775 lbm/hr and a step down from 400 to 125 lbm/hr, both with "the power
turbine speeds ... held constant by suppressing the NP integration" [pdf p.39].

The target is Ballin's own model trace (the solid line), not the GE reference markers he was
himself validating against. This file compares the initial trim and the settled state on
each panel; the whole curve is `test_whole_curve.py`.

Heat sink on: Figures 9 and 10 were generated with the station 4.1 heat sink active,
inferred rather than printed (see `docs/notes/heat-sink-configuration.md`). Table 1 and
B1-B6 are heat sink off. Eq. 50 has unit DC gain, so the switch cannot move a trim.

Tolerances are ours (SCOPE.md); the report's only numeric transient claim is a 1-2 % NG
overestimate on the open-loop step [pdf pp.45-46].
"""

from __future__ import annotations

import csv
import functools
from pathlib import Path

import numpy as np
import pytest

from t700 import constants as c
from t700 import maps, realtime, trim
from t700.engine import Ambient
from t700.units import wf_pps_from_pph

REF = Path(__file__).resolve().parent.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)

# Figure, final fuel flow, measured step time
STEP_TIME = {9: 0.5232, 10: 0.5217}
"""When the fuel steps, seconds, measured from the WFPH panel's own 50 % crossing on the
deskewed panels (#37, #63). The caption states both flow levels but not the instant
[pdf pp.45-46]. A few ms matter here: a transient sampled every 7 ms."""

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
"""The trim is the same physics the steady tests check, so this is really a check that the
figures' own y calibration is sound. See SCOPE.md."""

FINAL_TOL_PCT = {9: 1.0, 10: 1.0}
"""Settled state against the figure, both figures. See SCOPE.md and `_settled_window`."""


UNTRUSTED: set[tuple[int, str]] = set()
"""Panels excluded from the comparison, with the evidence. Currently none: a panel that
disagrees between Figures 9 and 10 at their shared 400 lbm/hr trim is the place to look
(#63)."""


def read_reference(
    name: str | Path, xc: str | None = None, yc: str | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Two columns of a digitized CSV under `data/reference/`, sorted by the first.

    `xc` and `yc` name the columns and default to the file's first two. The one reader for
    the reference data: every file that compares against a figure reads it through this.
    """
    path = REF / name
    rows = list(csv.DictReader(ln for ln in path.open() if not ln.startswith("#")))
    if xc is None:
        xc, yc = list(rows[0])[:2]
    x = np.array([float(r[xc]) for r in rows])
    y = np.array([float(r[yc]) for r in rows])
    o = np.argsort(x)
    return x[o], y[o]


def _reference(fig: int, key: str) -> tuple[np.ndarray, np.ndarray]:
    return read_reference(f"fig{fig:02d}_{key}_model.csv", "t_s", "value")


def _settled_window(tb, vb):
    """The reference's own settled samples, and the time its record actually stops.

    The 4.6 upper bound is load-bearing: `fig09_t41_model.csv` and `fig09_torq45_model.csv`
    each carry a spurious trailing sample at t ~ 4.9 (open question #48). Callers average
    ours only up to the returned stop time (about 4.47 s, where Ballin's frame ends), since
    on Figure 10 nothing has settled by then (#68).
    """
    sel = (tb > 4.0) & (tb < 4.6)
    if not sel.any():
        return vb[sel], 0.0
    return vb[sel], float(tb[sel].max())


def fuel_step(
    wf_hi: float,
    t_step: float,
    heat_sink: bool = True,
    duration_s: float = 5.0,
    wf0_pph: float = 400.0,
    **kw,
) -> dict[str, np.ndarray]:
    """Figures 9 and 10's experiment: a fuel step from a trim, NP integration suppressed.

    From the `wf0_pph` trim at NP_des to `wf_hi` lbm/hr at `t_step` seconds, on the
    report's 7 ms frame; a hold is a step to the same flow. `kw` goes to `realtime.run`
    (`tol`, for one). Uncached, so use it where a run's side effects are measured -- clamp
    counts -- or a constant or function is patched around it; `_run` is the cached one.
    """
    wf0 = wf_pps_from_pph(wf0_pph)
    r0 = trim.solve(wf0, c.NP_DES, AMB)
    hi = wf_pps_from_pph(wf_hi)
    st = realtime.from_trim(r0)
    return realtime.run(
        st,
        lambda t: wf0 if t < t_step else hi,
        AMB,
        duration_s=duration_s,
        dt=0.007,
        q_req_ftlbf=r0.frame.q_pt_ftlbf,
        integrate_np=False,
        heat_sink=heat_sink,
        **kw,
    )


@functools.cache
def _run(wf_hi: float, t_step: float, heat_sink: bool) -> dict[str, np.ndarray]:
    """`fuel_step` from 400 lbm/hr over 5 s, run once and shared by every comparison.

    The arrays are read-only, so no caller can alter what another reads. Pass `heat_sink`
    positionally: the cache is keyed on the call's spelling.
    """
    tr = fuel_step(wf_hi, t_step, heat_sink)
    for v in tr.values():
        v.flags.writeable = False
    return tr


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
    tr = _run(wf_hi, t_step, True)
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
    tr = _run(wf_hi, t_step, True)
    ours = PANELS[key](tr)

    late_them, t_last = _settled_window(tb, vb)
    if late_them.size == 0:
        pytest.skip("no settled reference samples")  # fig 10 TORQ45 ends before t = 4.0
    late_us = ours[(tr["t"] > 4.0) & (tr["t"] <= t_last)]
    dev = (late_us.mean() - late_them.mean()) / abs(late_them.mean()) * 100.0
    assert abs(dev) < FINAL_TOL_PCT[fig], (
        f"fig {fig} {key}: settled {late_us.mean():.1f} against Ballin's "
        f"{late_them.mean():.1f}, {dev:+.2f} %"
    )


def test_step_up_stays_inside_the_compressor_map():
    """Figure 9 keeps the engine on data, which is why it agrees to a few percent."""
    maps.reset_clamps()
    tr = fuel_step(775.0, STEP_TIME[9])  # uncached: the clamps it counts are this run's
    ngc = 100.0 * tr["ng"] / c.NG_DES
    assert ngc.min() > 88.0, "the step up should not approach the map's lower edge"
    assert ngc.max() < 101.0
    rep = maps.clamp_report()
    assert not any(k.startswith("f1@") and v > 50 for k, v in rep.items()), (
        f"the compressor map should barely be extrapolated on the step up: {rep}"
    )
    maps.reset_clamps()


def test_step_down_runs_off_the_bottom_of_the_maps():
    """Figure 10's chop: our speed at Ballin's last printed sample is within the panel's read
    error of his, and the run never leaves `f1`.

    Compare at his last sample, never at our minimum: neither trace has settled when his
    record stops (#68). The run's own minimum is a tripwire on the pressure solves: Eq. 80's
    fixed-point iteration does not converge where f9's elasticity is below -1, and an
    unconverged solve plunges NG too deep. Interpolating `f1` along Figure A1's construction
    lines (`maps.SpeedMap`) blends the speed lines' surge limits too, so nothing clamps (#58).
    `f6` is a constant (`maps.f6`) and cannot clamp.
    """
    maps.reset_clamps()
    tr = fuel_step(125.0, STEP_TIME[10])  # uncached: the clamps it counts are this run's
    ngc = 100.0 * tr["ng"] / c.NG_DES
    rep = maps.clamp_report()

    tb, vb = _reference(10, "pcng")
    ours_at_end = float(np.interp(tb[-1], tr["t"], ngc))
    gap = ours_at_end - vb[-1]
    assert abs(gap) < 0.24, (
        f"at Ballin's last sample (t = {tb[-1]:.3f} s) we read {ours_at_end:.3f} %NGc "
        f"against his {vb[-1]:.3f}, a gap of {gap:+.3f} %NG, outside this panel's read error"
    )

    # The run's own minimum (its last frame): a minimum near 70 % means a pressure solve has
    # stopped converging.
    assert 73.4 < ngc.min() < 74.2, (
        f"NGc bottoms at {ngc.min():.2f} % over the full 5 s: check that both pressure "
        f"solves still converge (not a comparison against Ballin)"
    )
    assert tr["far"].min() < 0.010, "the fuel-air ratio still falls below f6's plotted range"
    assert "f1@65" not in rep and "f1" not in rep, (
        f"the chop is reading f1 outside its data again: {rep}"
    )
    assert "f6" not in rep, "f6 is a constant and cannot be clamped -- see maps.f6"
    maps.reset_clamps()

    # Neither trace has settled by the end of Ballin's record, so nothing measured there may
    # be called a settled value.
    window = (tr["t"] >= 3.5) & (tr["t"] <= 4.5)
    slope = float(np.polyfit(tr["t"][window], ngc[window], 1)[0])
    assert slope < -1.0, (
        f"our Figure 10 run is still decaying at {slope:.2f} %NG/s over 3.5-4.5 s; "
        f"nothing measured in that window may be called a settled value"
    )


def _t41_overshoot(t: np.ndarray, v: np.ndarray, t_step: float) -> float:
    """Peak T41 above the settled value. Windowed at t < 4.6 -- see `_settled_window`."""
    settled = v[(t > 4.0) & (t < 4.6)].mean()
    return v[(t > t_step) & (t < 4.6)].max() - settled


def test_heat_sink_closes_most_of_the_t41_overshoot_gap():
    """The T41 overshoot against Figure 9, as a ratio to Ballin's, and a ratchet on it.

    The heat sink integrates the station 4.1 metal temperature (Eqs. 48-49), not Eq. 50's
    collapsed lead-lag, which is valid only for constant coefficients; Eqs. 51 and 53 make
    them functions of T41 and W41 (see `realtime._heat_sink`). Leaving the band in either
    direction means re-measuring; the lower bound makes an improvement visible.
    """
    tb, vb = _reference(9, "t41")
    late_them = float(np.median(vb[(tb > 3.5) & (tb < 4.5)]))
    theirs = float(vb.max()) - late_them

    tr = _run(775.0, STEP_TIME[9], True)
    ours_trace = PANELS["t41"](tr)
    late_us = float(np.median(ours_trace[tr["t"] > 4.0]))
    ours = float(ours_trace.max()) - late_us

    ratio = ours / theirs
    assert 0.6 < ratio < 1.2, (
        f"T41 overshoot is {ours:.0f} degR against Ballin's {theirs:.0f}, a ratio of "
        f"{ratio:.2f}x. Outside this band, re-measure -- and if below, tighten."
    )


def test_the_heat_sink_is_what_closed_it():
    """The same step with and without Eqs. 48-53, so a silently disabled heat sink fails
    here with its cause named.
    """
    bare = _run(775.0, STEP_TIME[9], False)
    sunk = _run(775.0, STEP_TIME[9], True)
    off = _t41_overshoot(bare["t"], bare["t41"], STEP_TIME[9])
    on = _t41_overshoot(sunk["t"], sunk["t41"], STEP_TIME[9])
    assert on < 0.6 * off, (
        f"the heat sink should roughly halve the T41 overshoot; it went {off:.0f} -> {on:.0f} degR"
    )
