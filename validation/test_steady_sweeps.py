"""Figures 6, 7 and 8 -- the steady-state sweeps [pdf pp.40-42].

**These were digitized on day one and no test touched them until 2026-09-11.** Every other
comparison in this repository -- Table B.1, Table 1, all twelve Appendix B matrices --
evaluates the model at the *same three operating points*. Three points cannot tell a model
that is right from one that is right there, which is exactly the failure a sweep catches.

The sweeps are Ballin's own real-time model over the full fuel-flow range: Figure 6 is
%NG against Wf (29 points, 136 to 810 lb/hr), Figure 7 is shaft horsepower against Wf
(15 points), Figure 8 is Ps3 against %NG (27 points).

## Continuation is not optional

A cold-started trim fails above about 750 lb/hr and succeeds when seeded from the next
lower solution. The first version of this comparison omitted that and reported the model
as unable to cover a quarter of Figure 6's range; it covers all of it but the last point.
`trim.sweep` exists for this reason, and these tests use the same ladder.

## What each figure actually tests

Figure 8 is the strongest independent check in the repository: it relates Ps3 to NG with
**no fuel flow in it at all**, so it validates the compressor and the station-3 volume
without touching the combustor. It agrees to 1.6 % worst case over 24 points.

Figure 6 does not, and the reason is a conflict inside the report -- see
`test_figure_6_offset_is_the_table_b1_conflict`.
"""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from t700 import constants as c
from t700 import trim
from t700.engine import Ambient, frame
from t700.units import shp_from_torque, wf_pps_from_pph

REF = Path(__file__).resolve().parent.parent / "data" / "reference"
AMB = Ambient(14.696, 518.67)
NP_RPM = c.NP_DES
"""Open question #35: the sweeps' power turbine speed is not printed. NP_des is the
reading the report's own text points at, and NG is decoupled from NP at steady state --
sweeping 16000-22000 rpm moves the NG-vs-Wf curve by under 0.01 %NG -- so the choice
cannot explain any deviation below."""


def _reference(name: str, xc: str, yc: str):
    rows = list(csv.DictReader(ln for ln in (REF / name).open() if not ln.startswith("#")))
    x = np.array([float(r[xc]) for r in rows])
    y = np.array([float(r[yc]) for r in rows])
    o = np.argsort(x)
    return x[o], y[o]


@lru_cache(maxsize=1)
def _rungs() -> tuple[tuple[float, object], ...]:
    """A dense continuation ladder, built once and independent of any query.

    Chaining through the reference flows themselves is not enough: Figure 7's points jump
    702 -> 750 -> 801 lb/hr, and a 50 lb/hr step is too coarse to carry the solution
    through the top of the compressor map. The ladder must be fine and the queries must
    hang off it, not the other way round.
    """
    out, guess = [], None
    for pph in np.arange(120.0, 815.0, 10.0):
        r = trim.solve(wf_pps_from_pph(float(pph)), NP_RPM, AMB, guess=guess)
        if r.trustworthy:
            guess = r.state
            out.append((float(pph), r.state))
    return tuple(out)


def _ladder(fuel_flows):
    """Trim at each requested flow, seeded from the nearest lower rung."""
    rungs = _rungs()
    out = {}
    for pph in sorted(float(v) for v in fuel_flows):
        lower = [st for x, st in rungs if x <= pph]
        r = trim.solve(wf_pps_from_pph(pph), NP_RPM, AMB, guess=lower[-1] if lower else None)
        if not r.trustworthy:
            continue
        f = frame(r.state, wf_pps_from_pph(pph), AMB)
        out[pph] = dict(
            ng_pct=100.0 * r.state.ng_rpm / c.NG_DES,
            shp=shp_from_torque(f.q_pt_ftlbf, NP_RPM),
            ps3=f.ps3_psia,
        )
    return out


def test_the_sweep_covers_the_whole_printed_range_but_the_last_point():
    """We trim everywhere Ballin's Figure 6 plots except its final point.

    That point is 100.31 %NG and `f1`'s top speed line is 100 %, so trimming it would mean
    extrapolating the compressor map off its printed data. Refusing is correct behaviour,
    not a shortfall -- but it is a real boundary and it is stated here rather than left to
    be rediscovered.
    """
    x, _ = _reference("fig06_realtime.csv", "wf_pph", "ng_pct")
    got = _ladder(x)
    missing = [v for v in x if v not in got]
    assert len(missing) <= 1, f"more of Figure 6 is now unreachable: {missing}"
    if missing:
        assert missing[0] > 800.0, f"the unreachable point moved down to {missing[0]:.0f} lb/hr"
        assert _ladder([800.0])[800.0]["ng_pct"] > 99.5, "we should still reach ~100 %NG"


def test_figure_8_ps3_against_ng_tracks_across_the_range():
    """The strongest independent check here: no fuel flow appears in it.

    Ps3 against %NG exercises the compressor map, the station-3 volume and the bleed
    schedules without the combustor, so it cannot be rescued by a compensating fuel-path
    error. 24 points, worst 1.6 %.
    """
    x, y = _reference("fig08_realtime.csv", "ng_pct", "ps3_psia")
    got = _ladder(np.linspace(150, 800, 60))
    ng = np.array([v["ng_pct"] for v in got.values()])
    ps = np.array([v["ps3"] for v in got.values()])
    o = np.argsort(ng)
    inside = (x >= ng[o].min()) & (x <= ng[o].max())
    assert inside.sum() >= 20, "the sweep no longer spans Figure 8"
    dev = (np.interp(x[inside], ng[o], ps[o]) - y[inside]) / y[inside] * 100.0
    assert np.abs(dev).max() < 3.0, f"worst Ps3 deviation {np.abs(dev).max():.1f} %"
    assert abs(dev.mean()) < 1.0, f"mean Ps3 deviation {dev.mean():+.2f} %"


def test_figure_7_shaft_power_tracks_except_at_the_lowest_power():
    """Shaft power against fuel flow: within 5 % over 13 of 15 points.

    The exception is the low-power end, +11.3 % at 224 lb/hr, decaying monotonically to
    -0.7 % by 700. That is the same conditioning the descent trim has -- shaft power is a
    small difference of large enthalpies there, so upstream error is amplified -- and the
    same end where Table B.1 and Figure 7 disagree with *each other* by 5.49 % (#46).
    """
    x, y = _reference("fig07_realtime.csv", "wf_pph", "shp")
    got = _ladder(x)
    dev = np.array(
        [(got[k]["shp"] - v) / abs(v) * 100.0 for k, v in zip(x, y, strict=True) if k in got]
    )
    within = int((np.abs(dev) < 5.0).sum())
    assert within >= 12, f"only {within} of {dev.size} points within 5 %"
    assert np.abs(dev).max() < 13.0, f"worst SHP deviation {np.abs(dev).max():.1f} %"
    assert dev[0] > dev[-1], "the low-power end should be the worst; the trend has changed"


def test_figure_6_offset_is_the_table_b1_conflict():
    """Figure 6 is the one sweep we do not track, and it is a conflict in the report.

    We sit a systematic -1.38 %NG below Figure 6 across its whole range. That is not our
    error against the report -- it is the report disagreeing with itself, logged as open
    question #46: at Wf = 267.7 / 349.3 / 476.3 Table B.1 sits -1.50 / -1.18 / -0.83 %NG
    below Figure 6, and our sweep reproduces -1.55 / -1.07 / -0.83 at those same flows.
    We match Table B.1 -- printed numbers, no digitization -- to 0.014 / 0.132 / 0.026 %.

    What this test adds to #46 is that the conflict is a **whole-curve** systematic, not
    three points: it runs from -2.0 %NG at 157 lb/hr to -0.4 % at 709, narrowing
    monotonically with power. #46 recorded only the three trim conditions.

    It asserts the offset's shape rather than its absence, so that if the model ever
    stopped matching Table B.1 this would move.
    """
    x, y = _reference("fig06_realtime.csv", "wf_pph", "ng_pct")
    got = _ladder(x)
    pairs = [(k, got[k]["ng_pct"] - v) for k, v in zip(x, y, strict=True) if k in got]
    flows = np.array([p[0] for p in pairs])
    dev = np.array([p[1] for p in pairs])

    assert dev.max() < 0.0, "the offset should be negative everywhere -- we sit below Fig. 6"
    assert -2.5 < dev.mean() < -0.9, f"mean offset {dev.mean():+.2f} %NG"
    lo = dev[flows < 250.0].mean()
    hi = dev[flows > 600.0].mean()
    assert lo < hi - 0.8, (
        f"the offset should narrow with power: {lo:+.2f} %NG below 250 lb/hr against "
        f"{hi:+.2f} above 600. If it no longer does, #46's reading has changed."
    )


@pytest.mark.parametrize("pph,ng_b1", [(267.7, 38072.0), (349.3, 39768.0), (476.3, 41638.0)])
def test_the_sweep_still_passes_through_table_b1(pph: float, ng_b1: float):
    """The sweep and the trim tests must agree where they overlap.

    Trivial-looking, and it is the tie that makes the Figure 6 result interpretable: it is
    only evidence about the report's internal conflict if we are genuinely on Table B.1
    at the same flows, rather than drifting and coincidentally landing nearby.
    """
    got = _ladder([pph])
    dev = (got[pph]["ng_pct"] * c.NG_DES / 100.0 - ng_b1) / ng_b1 * 100.0
    assert abs(dev) < 0.2, f"{pph} lb/hr: NG {dev:+.3f} % from Table B.1"
