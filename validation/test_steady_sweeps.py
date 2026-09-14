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
FIG6_RMS_PCT_NG = 0.35
FIG6_WORST_PCT_NG = 1.25
"""Deviation from Figure 6's real-time series, in %NG. Both declared in `SCOPE.md`.

Measured over the 29 digitized points: mean **+0.000**, rms **0.234**, worst **1.006**.
Before the page-40 frame fix this was a whole-curve systematic running -2.0 %NG at the low
end to -0.4 at the high.

The rms is the meaningful bound; the worst is one glyph. It sits at 249.2 lb/hr, where the
figure's own local slope is 0.060 %NG per lb/hr -- so 1.006 %NG is **16.7 lb/hr of abscissa
error on one marker**, not a speed disagreement. The next worst is 0.585 at 142.2 lb/hr,
which on that steep part of the curve is 2.1 lb/hr."""

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


def test_the_sweep_covers_the_whole_printed_range():
    """We trim at every point Ballin's Figure 6 plots. All twenty-nine.

    **This test was called `..._but_the_last_point` and said we did not** -- "that point is
    100.31 %NG and `f1`'s top speed line is 100 %, so trimming it would mean extrapolating
    the compressor map off its printed data". That was true when it was written and the
    page-40 frame fix of 2026-09-14 (#64) removed it: every digitized speed on that figure
    had read high by `(105 - NG) * 0.081`, and with the frame found correctly the last
    marker is **806.0 lb/hr at 99.98 %NG**, inside the map. We reach it at 99.97.

    The assertion still tolerates one unreachable point, because the boundary is real even
    though this figure no longer crosses it -- but the name and the docstring now describe
    what happens rather than what used to.
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


def test_we_now_track_figure_6_across_its_whole_range():
    """Figure 6 was the one sweep we did not track. **It was our digitizer.**

    This test asserted the opposite until 2026-09-14: "we sit a systematic -1.38 %NG below
    Figure 6 across its whole range ... that is not our error against the report, it is the
    report disagreeing with itself, logged as open question #46". The disagreement was
    entirely ours.

    `find_frame` in `tools/digitize_fig678.py` located a horizontal frame by the fraction of
    a row that is inked. pdf p.40 prints faint frames and is tilted by ~16 px across the
    plot, so the bottom frame's ink spreads over twenty rows and no row exceeds 0.24 fill --
    while the figure's **caption**, a dense line of text 190 px lower, reaches 0.42 and won.
    Figure 6's NG axis was therefore 2958-330 px tall against a true 2768-323, and every
    digitized speed read high by (105 - NG) * 0.081: +1.2 %NG at 90, +2.8 %NG at 70.

    That is the whole of #46's low-power half. Against Table B.1's three printed trims,
    Figure 6 read +0.83 / +1.18 / +1.50 %NG high; it now reads **-0.043 / -0.024 / -0.033**,
    rms 1.201 -> **0.034 %NG**. The two are the same dataset and always were.

    The check that settles it is held out and is in the tool: the sixteen minor ticks on
    each vertical frame land on Figure 6's printed 2.5 %NG grid to **rms 0.045 %NG**, where
    the caption-anchored map put them at 102.91, 100.58, 98.23 -- a 2.33 spacing, on no
    grid at all.
    """
    x, y = _reference("fig06_realtime.csv", "wf_pph", "ng_pct")
    got = _ladder(x)
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


def test_which_tables_are_extrapolated_at_which_trims():
    """The map extrapolation at the operating point: counted since day one, never surfaced.

    `TrimResult.clamps_at_solution` has always recorded this and `trustworthy` has always
    ignored all of it but `f1:parameter`, so a trim reading three tables outside their
    data reported clean. The 2026-09-13 engine-physics and numerical-mathematics audits
    both landed on it; this is the shape, measured:

        115 lbm/hr   f8, f9          NGc 65.4 %, just above the bottom of f1's range
                                     (110 lbm/hr lands AT NGc 65.000 and is rejected --
                                      `f1:parameter` clamps, so `trustworthy` is False.
                                      It trimmed at 65.01 under the pre-2026-09-14 map;
                                      the beta grid moved it a hundredth of a percent the
                                      other side of the edge, which is the edge doing its
                                      job rather than a change in the physics.)
        125          f8
        150          f8              (it read `f1@65` too until 2026-09-13, under
                                        constant-abscissa interpolation)
        200-725      none
        750 and up   f9              Ps9/P45 passes f9's tabulated 0.85012

    **`f6` used to appear here at every trim above about 590 lbm/hr** -- the top third of
    the power range, Figure 9's 775 lbm/hr endpoint included -- and it does not any more.
    Figure A6 [pdf p.61] draws a single horizontal line across a y axis spanning 0.88 to
    1.10, with an `x` at each frame edge and nothing between; the two digitized endpoints
    differ by 7.3e-5 against a 4.06e-4 read error apiece, which is 0.13 sigma and 0.46 px.
    The slope was the page's skew. As of 2026-09-13 `f6` is loaded as the constant it is,
    and **a constant has no domain to leave**, so what was the most-reported clamp in the
    project has ceased to exist rather than been suppressed. See `maps.f6`.

    **`f1` no longer appears in this census at any trim.** It was clamped at 150 lbm/hr --
    the 65 % speed line asked for pressure ratios past its own last knot while it brackets
    the 80 % line from below -- and that was open question #58, "f1's 15-point data hole",
    believed to be where Figure 10's residual lived. On the beta grid each speed line ends
    at its own surge limit and blending two lines blends their limits too, so the query
    stays on the map: **zero of 715 frames** on either published transient. #58 closed
    2026-09-14 and the row records that the hole was real as a statement about the data and
    irrelevant as a diagnosis.

    `trustworthy` is deliberately not tightened to `fully_on_data`: that would reject the
    top third of the power range over an `f6` clamp worth 1e-4.
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
    # Continued, not cold-started: above ~725 lbm/hr a cold trim from the design-point
    # guess lands off the map, which `test_the_sweep_refuses_past_the_maps_top_speed_line`
    # is about. The ladder is the one `trim.sweep` exists for.
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
