"""Figures 6, 7 and 8 -- the steady-state sweeps [pdf pp.40-42], against Ballin's model
and against the two GE models he validated against.

Every other comparison evaluates the model at Table B.1's three trims; a sweep catches a model
that is right only there. Trims come by continuation (`trim.sweep`): a cold-started trim
fails near the top of the range.

The target is Ballin, and only Ballin. The GE series (status-81 and unbalanced torque) give a
scale -- how far apart careful models of this engine sit -- not a verdict; closeness to GE
would be a symptom, not evidence. No model value derives from them: `src/t700/` never reads
`data/reference/` (`tests/test_architecture.py`).
"""

from __future__ import annotations

from functools import cache, lru_cache
from pathlib import Path

import numpy as np
import pytest

from t700 import constants as c
from t700 import trim
from t700.engine import Ambient, frame
from t700.units import shp_from_torque, wf_pps_from_pph
from test_fuel_step import read_reference
from test_trim_points import CASES

REF = Path(__file__).resolve().parent.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)
FIG6_RMS_PCT_NG = 0.35
FIG6_WORST_PCT_NG = 1.25
"""Deviation from Figure 6's real-time series, in %NG. Both declared in `SCOPE.md`. The rms
is the meaningful bound; the worst is one marker on a steep stretch of the curve, i.e.
abscissa error rather than a speed disagreement."""

FIG7_RMS_SHP = 2.5
FIG7_WORST_SHP = 6.5
"""Deviation from Figure 7's real-time series, in shaft horsepower. Both declared in
`SCOPE.md`, set above the measured values to catch a regression."""

NP_RPM = c.NP_DES
"""NP_des, Table A.1 [pdf p.55]. Figures 6-8 run at design speed (#35, read off Figure 7
itself), Ballin's curves and the GE models' alike."""


@lru_cache(maxsize=2)
def _rungs(np_rpm: float) -> tuple[tuple[float, object], ...]:
    """A dense continuation ladder, built once per speed and independent of any query.

    Chaining through the reference flows is not enough: Figure 7's ~50 lb/hr steps are too
    coarse to carry the solution through the top of the compressor map.
    """
    out, guess = [], None
    for pph in np.arange(120.0, 815.0, 10.0):
        r = trim.solve(wf_pps_from_pph(float(pph)), np_rpm, AMB, guess=guess)
        if r.trustworthy:
            guess = r.state
            out.append((float(pph), r.state))
    return tuple(out)


def _ladder(fuel_flows, np_rpm: float = NP_RPM):
    """Trim at each requested flow, seeded from the nearest lower rung."""
    rungs = _rungs(np_rpm)
    out = {}
    for pph in sorted(float(v) for v in fuel_flows):
        lower = [st for x, st in rungs if x <= pph]
        r = trim.solve(wf_pps_from_pph(pph), np_rpm, AMB, guess=lower[-1] if lower else None)
        if not r.trustworthy:
            continue
        f = frame(r.state, wf_pps_from_pph(pph), AMB)
        out[pph] = dict(
            ng_pct=100.0 * r.state.ng_rpm / c.NG_DES,
            shp=shp_from_torque(f.q_pt_ftlbf, np_rpm),
            ps3=f.ps3_psia,
        )
    return out


def _ps3_at_ng(ng_pct: np.ndarray, rungs: int, np_rpm: float):
    """Our Ps3 at the given %NG, read off an `rungs`-point ladder over 150-800 lb/hr.

    Returns `(mask, values)`: Figure 8 carries no fuel flow, so the comparison is made in NG
    and only where the ladder spans it.
    """
    pts = list(_ladder(np.linspace(150, 800, rungs), np_rpm).values())
    ng = np.array([p["ng_pct"] for p in pts])
    ps = np.array([p["ps3"] for p in pts])
    o = np.argsort(ng)
    m = (ng_pct >= ng[o].min()) & (ng_pct <= ng[o].max())
    return m, np.interp(ng_pct[m], ng[o], ps[o])


def test_figure_8_ps3_against_ng_tracks_across_the_range():
    """The strongest independent check here: Ps3 against %NG has no fuel flow in it, so it
    exercises the compressor, the station-3 volume and the bleeds without the combustor.
    """
    x, y = read_reference("fig08_realtime.csv", "ng_pct", "ps3_psia")
    inside, ours = _ps3_at_ng(x, 60, NP_RPM)
    assert inside.sum() >= 20, "the sweep no longer spans Figure 8"
    dev = (ours - y[inside]) / y[inside] * 100.0
    assert np.abs(dev).max() < 3.0, f"worst Ps3 deviation {np.abs(dev).max():.1f} %"
    assert abs(dev.mean()) < 1.0, f"mean Ps3 deviation {dev.mean():+.2f} %"


def test_figure_7_shaft_power_tracks_across_the_range():
    """Shaft power against fuel flow over every printed point, in shaft horsepower.

    In shp, not percent: the low-power points fall to 9 shp, where the axis read error alone
    is 15 % of the value, so a percentage there measures the denominator (#83).
    """
    x, y = read_reference("fig07_realtime.csv", "wf_pph", "shp")
    got = _ladder(x)
    missing = [v for v in x if v not in got]
    assert not missing, f"part of Figure 7 is now unreachable: {missing}"
    dev = np.array([got[k]["shp"] - v for k, v in zip(x, y, strict=True)])
    rms = float(np.sqrt(np.mean(dev**2)))
    assert rms < FIG7_RMS_SHP, f"rms deviation from Figure 7 {rms:.2f} shp over {dev.size} points"
    assert np.abs(dev).max() < FIG7_WORST_SHP, (
        f"worst deviation from Figure 7 {np.abs(dev).max():.2f} shp"
    )


def test_we_now_track_figure_6_across_its_whole_range():
    """Figure 6's NG against fuel flow, over every one of its points.

    Trap: pdf p.40's frames are faint and tilted, and a digitizer that finds a frame by row
    fill anchors the axis on the denser caption instead, reading every speed high; that was
    the whole of #46's low-power half (#64). Ballin's top point sits on `f1`'s 100 % top
    line, which ours reaches at a hair less fuel; `maps.F1_TOP_ALLOWANCE_PCT` covers it (#77).
    """
    x, y = read_reference("fig06_realtime.csv", "wf_pph", "ng_pct")
    got = _ladder(x)
    missing = [v for v in x if v not in got]
    assert not missing, f"part of Figure 6 is now unreachable: {missing}"
    assert _ladder([800.0])[800.0]["ng_pct"] > 99.5, "we should still reach ~100 %NG"

    pairs = [(k, got[k]["ng_pct"] - v) for k, v in zip(x, y, strict=True) if k in got]
    dev = np.array([p[1] for p in pairs])

    rms = float(np.sqrt(np.mean(dev**2)))
    assert rms < FIG6_RMS_PCT_NG, (
        f"rms deviation from Figure 6 {rms:.3f} %NG over {dev.size} points"
    )
    assert np.abs(dev).max() < FIG6_WORST_PCT_NG, (
        f"worst deviation {dev[np.argmax(np.abs(dev))]:+.3f} %NG at "
        f"{[p[0] for p in pairs][int(np.argmax(np.abs(dev)))]:.1f} lb/hr"
    )
    assert abs(dev.mean()) < 0.2, f"mean offset {dev.mean():+.3f} %NG"


@pytest.mark.parametrize("pph,ng_b1", sorted((case.wf_pph, case.ng_rpm) for case in CASES))
def test_the_sweep_still_passes_through_table_b1(pph: float, ng_b1: float):
    """The sweep and the trim tests agree where they overlap, so the Figure 6 result is
    evidence about the report rather than about a drifting sweep.
    """
    got = _ladder([pph])
    dev = (got[pph]["ng_pct"] * c.NG_DES / 100.0 - ng_b1) / ng_b1 * 100.0
    assert abs(dev) < 0.2, f"{pph} lb/hr: NG {dev:+.3f} % from Table B.1"


def test_which_tables_are_extrapolated_at_which_trims():
    """Which tables each trim reads outside their data (`TrimResult.clamps_at_solution`).

    `trustworthy` checks only `f1:parameter`, deliberately: tightening it to `fully_on_data`
    would reject the top of the power range over the f9 clamp there. `f6` is the constant
    Figure A6 [pdf p.61] draws and has no domain to leave (`maps.f6`); `f1` stays on its data
    because blended speed lines blend their surge limits too (#58).
    """
    expected = {
        115.0: ("f8", "f9"),
        125.0: ("f8",),
        150.0: ("f8",),
        200.0: (),
        300.0: (),
        400.0: (),
        550.0: (),
        600.0: (),
        700.0: (),
        775.0: ("f9",),
    }
    # Continued, not cold-started: near the top of the range a cold trim from the
    # design-point guess lands off the map.
    ladder = sorted(expected)
    results = dict(zip(ladder, trim.sweep(ladder), strict=True))

    for pph, tables in expected.items():
        r = results[pph]
        assert r.trustworthy, f"{pph:.0f} lbm/hr no longer trims"
        assert r.extrapolated_tables == tables, (
            f"{pph:.0f} lbm/hr extrapolates {r.extrapolated_tables}, on record "
            f"{tables}. A change here moves which results carry an extrapolation caveat."
        )
        assert r.fully_on_data == (not tables)

    top = results[775.0]
    assert top.trustworthy and not top.fully_on_data, (
        "the f9 clamp near the top of the range must not make a trim untrustworthy -- see "
        "the docstring for why that line is drawn where it is"
    )


# ------------------------------------------------- against the two GE models, same figures


def _curve(fig: int, xs: np.ndarray):
    """Our model sampled at the given abscissae; returns (mask, values)."""
    if fig == 8:
        return _ps3_at_ng(xs, 80, NP_RPM)
    key = "ng_pct" if fig == 6 else "shp"
    ours = _ladder(xs, NP_RPM)
    got = [ours.get(float(v), {}).get(key) for v in xs]
    m = np.array([v is not None for v in got])
    return m, np.array([v for v in got if v is not None])


SPEC = {
    6: ("fig06", "wf_pph", "ng_pct", "abs"),
    7: ("fig07", "wf_pph", "shp", "pct"),
    8: ("fig08", "ng_pct", "ps3_psia", "pct"),
}


def _spreads(fig: int, series: str):
    """(Ballin - GE, ours - GE) on GE's own abscissae, in %NG or percent."""
    stem, xc, yc, mode = SPEC[fig]
    bx, by = read_reference(f"{stem}_realtime.csv", xc, yc)
    gx, gy = read_reference(f"{stem}_{series}.csv", xc, yc)
    inb = (gx >= bx.min()) & (gx <= bx.max())
    bal = np.interp(gx[inb], bx, by)
    m, ov = _curve(fig, gx[inb])
    if mode == "abs":
        return bal - gy[inb], ov - gy[inb][m]
    return (bal - gy[inb]) / np.abs(gy[inb]) * 100, (ov - gy[inb][m]) / np.abs(gy[inb][m]) * 100


@pytest.mark.parametrize("fig", [6, 7, 8])
@pytest.mark.parametrize("series", ["ge_status81", "ge_unbalanced"])
def test_our_distance_from_ge_is_comparable_to_ballins(fig: int, series: str):
    """We sit no further from the GE models than Ballin's own model does, give or take a
    point and a half. Not a claim that we beat him: a deviation from Ballin of a percent or
    two is within the scatter of his own validation.
    """
    d_bal, d_our = _spreads(fig, series)
    assert d_our.size > 3, "too few overlapping points to judge"
    assert abs(d_our.mean()) < abs(d_bal.mean()) + 1.5, (
        f"fig {fig} vs {series}: ours {d_our.mean():+.2f} against Ballin's "
        f"{d_bal.mean():+.2f} from GE"
    )


def test_on_figure_6_we_sit_where_ballin_sits_and_he_says_where_that_is():
    """[pdf p.39] "Gas generator speed is overestimated by 1 to 2 percent; this is reflected
    in the trim differences between the real-time model and the status-81 model in figure 6."

    We and Ballin sit the same distance above GE status-81: the displacement is his, and we
    inherit it.
    """
    d_bal, d_our = _spreads(6, "ge_status81")
    assert abs(d_our.mean() - d_bal.mean()) < 0.25, (
        f"ours {d_our.mean():+.3f} %NG from GE against Ballin's {d_bal.mean():+.3f} -- we "
        f"should inherit his displacement, not have our own"
    )
    for who, d in (("Ballin", d_bal), ("ours", d_our)):
        assert 0.5 < d.mean() < 2.0, (
            f"{who} sits {d.mean():+.3f} %NG above GE status-81; pdf p.39 says the "
            f"real-time model overestimates NG by 1 to 2 percent"
        )


GE_SCATTER_FLOOR = {6: 1.5, 7: 10.0, 8: 8.0}
"""How far apart GE's own two models must stay, per figure (%NG on Figure 6, percent on 7
and 8), for the scatter argument to hold. Set below the measured worst."""


def test_the_two_ge_models_disagree_with_each_other_substantially():
    """GE's own two models differ from each other by far more than we differ from Ballin:
    whatever "agreement" meant in 1988, it did not mean a few tenths of a percent.
    """
    for fig in (6, 7, 8):
        stem, xc, yc, mode = SPEC[fig]
        sx, sy = read_reference(f"{stem}_ge_status81.csv", xc, yc)
        ux, uy = read_reference(f"{stem}_ge_unbalanced.csv", xc, yc)
        inb = (ux >= sx.min()) & (ux <= sx.max())
        if inb.sum() < 3:
            continue
        s = np.interp(ux[inb], sx, sy)
        d = (s - uy[inb]) if mode == "abs" else (s - uy[inb]) / np.abs(uy[inb]) * 100
        assert np.abs(d).max() > GE_SCATTER_FLOOR[fig], (
            f"fig {fig}: the two GE models now agree to {np.abs(d).max():.2f}; "
            f"if so, the 1988 scatter argument needs re-reading"
        )


# ------------------------------------------------ the report's own printed claim about Fig. 6

FUEL_OVERESTIMATE_BAND = (81.0, 86.0)
FUEL_OVERESTIMATE_CEILING_PCT = 8.0
"""[pdf p.39] "The steady-state operation of the real-time model is bounded by the two
analysis models over the range of operation except between 81 and 86 percent of gas
generator speed. In this area, the real-time model tends to overestimate fuel consumption
by as much as five percent."


Read from the page raster, not the text layer. Checked on Figure 6's digitized series (#67).
A harsh currency: the NG-vs-Wf curve is steep through that band, so 1 %NG of speed error is
roughly 9 % of fuel. The ceiling is ours, declared in SCOPE.md.
"""


@cache
def _fig6_series():
    """Fuel flow as a function of NG for all three Figure 6 series, plus ours. Cached: the
    three tests that read it only read it."""
    out = {}
    for key in ("realtime", "ge_status81", "ge_unbalanced"):
        wf, ng = read_reference(f"fig06_{key}.csv", "wf_pph", "ng_pct")
        o = np.argsort(ng)
        out[key] = (ng[o], wf[o])
    grid = np.arange(125.0, 851.0, 2.5)
    ok = [(x, r) for x, r in zip(grid, trim.sweep(grid, NP_RPM, AMB), strict=True) if r.trustworthy]
    out["ours"] = (
        np.array([100.0 * r.state.ng_rpm / c.NG_DES for _, r in ok]),
        np.array([float(x) for x, _ in ok]),
    )
    return out


def _wf_at(series, ng: float) -> float:
    n, w = series
    return float(np.interp(ng, n, w))


def fuel_overestimate() -> tuple[np.ndarray, np.ndarray]:
    """Fuel flow over GE status-81's, in %, ours and Ballin's real-time line, every 0.5 %NG
    across `FUEL_OVERESTIMATE_BAND`. `validation/report.py` exports both peaks from this."""
    s = _fig6_series()
    lo, hi = FUEL_OVERESTIMATE_BAND
    ngs = np.arange(lo, hi + 0.01, 0.5)
    ours = np.array(
        [100.0 * (_wf_at(s["ours"], g) / _wf_at(s["ge_status81"], g) - 1.0) for g in ngs]
    )
    theirs = np.array(
        [100.0 * (_wf_at(s["realtime"], g) / _wf_at(s["ge_status81"], g) - 1.0) for g in ngs]
    )
    return ours, theirs


def test_the_printed_fuel_overestimate_band_is_reproduced():
    """We overestimate fuel consumption in 81-86 %NG against GE status-81, the standard the
    sentence is about, with a peak within 2 points of Ballin's own line. The same
    disagreement as open question #46, in the report's own currency.
    """
    ours, theirs = fuel_overestimate()

    assert ours.max() > 3.0, f"no overestimate in the printed band: peak {ours.max():+.2f} %"
    assert ours.max() < FUEL_OVERESTIMATE_CEILING_PCT, (
        f"we overestimate fuel by {ours.max():+.2f} % in 81-86 %NG against the report's "
        f"'as much as five percent' and its own line's {theirs.max():+.2f} %"
    )
    assert abs(ours.max() - theirs.max()) < 2.0, (
        f"our peak {ours.max():+.2f} % against Ballin's {theirs.max():+.2f} % -- we should "
        f"reproduce his excursion, not have our own"
    )


def test_above_the_band_we_need_less_fuel_than_status_81_as_the_report_says():
    """[pdf p.39] "At higher gas generator speeds, the real-time model displays the
    characteristic of the GE unbalanced torque model, requiring less fuel ... than the GE
    status-81 model."
    """
    s = _fig6_series()
    for ng in (88.0, 90.0, 95.0):
        ours = 100.0 * (_wf_at(s["ours"], ng) / _wf_at(s["ge_status81"], ng) - 1.0)
        theirs = 100.0 * (_wf_at(s["realtime"], ng) / _wf_at(s["ge_status81"], ng) - 1.0)
        assert ours < -2.0, f"{ng} %NG: we need {ours:+.2f} % fuel against status-81"
        assert abs(ours - theirs) < 2.0, f"{ng} %NG: ours {ours:+.2f} vs Ballin {theirs:+.2f}"


def test_we_leave_the_ge_envelope_no_more_often_than_ballin_does():
    """[pdf p.39] "bounded by the two analysis models over the range of operation except
    between 81 and 86 percent of gas generator speed."

    Only testable where both GE series have data (the unbalanced series starts near 83 %NG).
    Ballin's own line leaves the envelope at a speed his sentence does not exempt, so the
    bound is against his count, not the printed claim.
    """

    s = _fig6_series()
    lo = max(s["ge_status81"][0].min(), s["ge_unbalanced"][0].min())
    hi = min(s["ge_status81"][0].max(), s["ge_unbalanced"][0].max(), s["ours"][0].max())
    out_b = out_o = 0
    for ng in np.arange(np.ceil(lo), np.floor(hi) + 0.01, 1.0):
        a, b = _wf_at(s["ge_status81"], ng), _wf_at(s["ge_unbalanced"], ng)
        band = (min(a, b), max(a, b))
        out_b += not (band[0] <= _wf_at(s["realtime"], ng) <= band[1])
        out_o += not (band[0] <= _wf_at(s["ours"], ng) <= band[1])
    assert out_o <= out_b + 3, (
        f"we leave the GE envelope at {out_o} sampled speeds against Ballin's {out_b}"
    )
