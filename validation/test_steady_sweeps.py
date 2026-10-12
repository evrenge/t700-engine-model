"""Figures 6, 7 and 8 -- the steady-state sweeps [pdf pp.40-42], against Ballin's model
and against the two GE models he validated against.

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

Figure 6 did not, and the reason was taken for a conflict inside the report; it was our
digitizer -- see `test_we_now_track_figure_6_across_its_whole_range`.

## The two GE models on the same figures

**The target is Ballin, and only Ballin.** Closeness to GE is not evidence that we are
right, and this file must not be read as claiming otherwise. We are replicating one
model; agreeing with a different model of the same engine is at best a coincidence and
at worst a symptom -- if our numbers ever sat closer to GE than to Ballin, the first
thing to suspect would be that we had imported something from the wrong source.

Audited 2026-09-12, because that suspicion deserves a check rather than a disclaimer:
**no number in the model derives from the GE or NASA-Lewis series.** The GE traces exist
only under `data/reference/`, nothing in `src/t700/` reads that directory at all
(`tests/test_architecture.py` now enforces it), and every transient assertion in
`test_fuel_step.py` targets `fig??_*_model.csv` -- Ballin's own curve -- never
`*_reference.csv`. The GE series has never touched a model value.

### What the GE series legitimately provides: a scale, not a verdict

Every deviation in this project is quoted against Ballin without any sense of how large a
deviation is normal between two careful models of the same engine. Figures 6, 7 and 8
each plot **three** series: Ballin's real-time model and two General Electric models --
the status-81 performance standard and an unbalanced-torque variant. Measuring all three
pairs gives that scale: Ballin's own mean distance from GE status-81 runs 0.45-1.97 %
with worst cases to 8.5 %, and from the unbalanced-torque model 1.4-5.3 % with worst
cases to 17 %. Ours from Ballin is 0.43-1.51 % mean.

That is a ruler, not a grade. It says a 1 % steady-state spread is ordinary among models
of this engine; it says nothing about whether our 1 % is in the right place.

### Where it does bear on a question: #46

On Figure 6 we sit -0.08 %NG from GE where Ballin sits +0.89 %NG. The useful content is
not "we are closer" but that **Figure 6 carries an offset of its own**, which is exactly
the disagreement open question #46 records between Table B.1 and Figure 6. Three sources,
and Table B.1 and GE fall together. That removes the presumption that the gap to Figure 6
is automatically our error -- it does not turn the gap into a credential.
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


@lru_cache(maxsize=2)
def _rungs(np_rpm: float) -> tuple[tuple[float, object], ...]:
    """A dense continuation ladder, built once per speed and independent of any query.

    Chaining through the reference flows themselves is not enough: Figure 7's points jump
    702 -> 750 -> 801 lb/hr, and a 50 lb/hr step is too coarse to carry the solution
    through the top of the compressor map. The ladder must be fine and the queries must
    hang off it, not the other way round. The GE comparison below builds its own at its
    own speed.
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
    """The strongest independent check here: no fuel flow appears in it.

    Ps3 against %NG exercises the compressor map, the station-3 volume and the bleed
    schedules without the combustor, so it cannot be rescued by a compensating fuel-path
    error. 24 points, worst 1.6 %.
    """
    x, y = read_reference("fig08_realtime.csv", "ng_pct", "ps3_psia")
    inside, ours = _ps3_at_ng(x, 60, NP_RPM)
    assert inside.sum() >= 20, "the sweep no longer spans Figure 8"
    dev = (ours - y[inside]) / y[inside] * 100.0
    assert np.abs(dev).max() < 3.0, f"worst Ps3 deviation {np.abs(dev).max():.1f} %"
    assert abs(dev.mean()) < 1.0, f"mean Ps3 deviation {dev.mean():+.2f} %"


def test_figure_7_shaft_power_tracks_except_at_the_lowest_power():
    """Shaft power against fuel flow: within 5 % over 13 of 15 points.

    The exception is the low-power end, +11.3 % at 224 lb/hr, decaying monotonically to
    -0.7 % by 700. That is the same conditioning the descent trim has -- shaft power is a
    small difference of large enthalpies there, so upstream error is amplified -- and the
    same end where Table B.1 and Figure 7 disagree with *each other* most -- though after the
    2026-09-14 digitizer corrections that disagreement is 0.734 % and not the 5.49 % this
    docstring used to quote, so it no longer covers the residual (#46; SCOPE.md, "Which
    source wins", and `report.py::internal_spread`, which computes it).
    """
    x, y = read_reference("fig07_realtime.csv", "wf_pph", "shp")
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

    ## And we trim at every point it plots, all twenty-nine

    A test named `..._but_the_last_point` said we did not -- "that point is 100.31 %NG and
    `f1`'s top speed line is 100 %, so trimming it would mean extrapolating the compressor
    map off its printed data". That was true when it was written and the same frame fix
    (#64) removed it: with the frame found correctly the last marker is **806.0 lb/hr at
    99.98 %NG**, inside the map. We reach it at 99.97. One unreachable point is still
    tolerated, because the boundary is real even though this figure no longer crosses it.
    """
    x, y = read_reference("fig06_realtime.csv", "wf_pph", "ng_pct")
    got = _ladder(x)
    missing = [v for v in x if v not in got]
    assert len(missing) <= 1, f"more of Figure 6 is now unreachable: {missing}"
    if missing:
        assert missing[0] > 800.0, f"the unreachable point moved down to {missing[0]:.0f} lb/hr"
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


# ------------------------------------------------- against the two GE models, same figures

NP_GE_RPM = 20895.0
"""The speed the GE comparison has trimmed at since it was written: Table B.1's [pdf p.67].
The sweeps against Ballin above use `NP_RPM`, NP_des (#35). Which one this comparison
should use is not yet decided; the two are kept apart so that merging the files moved
nothing."""


def _curve(fig: int, xs: np.ndarray):
    """Our model sampled at the given abscissae; returns (mask, values)."""
    if fig == 8:
        return _ps3_at_ng(xs, 80, NP_GE_RPM)
    key = "ng_pct" if fig == 6 else "shp"
    ours = _ladder(xs, NP_GE_RPM)
    got = [ours.get(float(v), {}).get(key) for v in xs]
    m = np.array([v is not None for v in got])
    return m, np.array([v for v in got if v is not None])


SPEC = {
    6: ("fig06", "wf_pph", "ng_pct", "abs"),
    7: ("fig07", "wf_pph", "shp", "pct"),
    8: ("fig08", "ng_pct", "ps3_psia", "pct"),
}


def _spreads(fig: int, series: str):
    """(Ballin - GE, ours - GE) on GE's own abscissae, in %NG or per cent."""
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
    """We sit no further from the GE standard than Ballin's own model does, give or take.

    "Give or take" is deliberate: the bound is that our mean distance does not exceed his
    by more than 1.5 percentage points (or 1.5 %NG on Figure 6). It is not a claim that we
    beat him -- on Figure 7's status-81 comparison he is slightly closer. The point is that
    the two are the same size, so a deviation of ours from Ballin of a per cent or two is
    within the scatter of his own validation rather than evidence of a defect.
    """
    d_bal, d_our = _spreads(fig, series)
    assert d_our.size > 3, "too few overlapping points to judge"
    assert abs(d_our.mean()) < abs(d_bal.mean()) + 1.5, (
        f"fig {fig} vs {series}: ours {d_our.mean():+.2f} against Ballin's "
        f"{d_bal.mean():+.2f} from GE"
    )


def test_on_figure_6_we_sit_where_ballin_sits_and_he_says_where_that_is():
    """Ballin's own sentence about Figure 6, reproduced.

    This test asserted the opposite claim until 2026-09-14 -- "Figure 6 is the one sweep we
    do not track, at -1.38 %NG ... against GE status-81 ours is -0.08 %NG mean and Ballin's
    is +0.89 %NG. He is the one displaced from GE there, not us." The -1.38 was our own
    digitizer reading Figure 6's y axis off the figure's caption; see
    `test_steady_sweeps.test_we_now_track_figure_6_across_its_whole_range`.

    With the axis right, we and Ballin are the same curve: **ours -- Ballin is +0.034 %NG
    mean and 0.220 %NG worst** over the fifteen shared points. Both then sit **+0.96 and
    +1.00 %NG above GE status-81**, which is the number the report prints for itself
    [pdf p.39]: "Gas generator speed is overestimated by 1 to 2 percent; this is reflected
    in the trim differences between the real-time model and the status-81 model in figure
    6." So the displacement from the hardware standard is his, we inherit it exactly, and
    the report says so in text.

    That is a far stronger statement than the one it replaces, and it is what makes the
    residual Figure 6 comparison interpretable at all.
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
"""How far apart GE's own two models must stay, per figure, for the scatter argument to
hold. Measured worst: **1.83 %NG on Figure 6, 15.1 % on Figure 7, 11.0 % on Figure 8.**

Figure 6's floor was 2.0 until 2026-09-14 and the measurement was 1.98. Both series on that
page are digitized against the same y axis, so the page-40 frame fix scaled their
difference down with everything else -- 1.98 -> 1.83, an 8 % reduction that is exactly the
axis correction and not a change in the report."""


def test_the_two_ge_models_disagree_with_each_other_substantially():
    """The scale that makes every other number here readable.

    GE's own two models differ from one another by far more than we differ from Ballin:
    1.83 %NG, 15.1 % and 11.0 % against our 0.22 %NG worst on Figure 6. Whatever
    "agreement" meant in 1988, it did not mean a few tenths of a per cent.
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
analysis models over the range of operation **except between 81 and 86 percent of gas
generator speed. In this area, the real-time model tends to overestimate fuel consumption
by as much as five percent.**"

Read from the page raster, not the text layer.

**This is the report's sharpest printed steady-state claim about itself, and nothing
checked it until 2026-09-14** -- found by that day's validation-quality audit. It is
checkable entirely from series already digitized on Figure 6, and it is a harsher currency
than the one this project usually reports in: the NG-vs-Wf curve is steep through that
band, so 1 %NG of speed error is roughly 9 % of fuel.

The ceiling is 8 % against a measured peak of 6.17 %, and against Ballin's own 5.82 %.
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


def test_the_printed_fuel_overestimate_band_is_reproduced():
    """We overestimate fuel consumption in 81-86 %NG, by about as much as the report says.

    Measured against GE status-81, which is the performance standard the sentence is
    about:

    | %NG | Ballin | ours |
    |---|---|---|
    | 81 | +0.92 % | +1.25 % |
    | 82 | +5.99 % | +6.31 % |
    | 83 | +5.49 % | +5.93 % |
    | 84 | +2.32 % | +6.16 % |
    | 85 | +0.72 % | +0.98 % |
    | 86 | -1.87 % | -1.52 % |

    So the excursion is real, positive, confined to the band the report names, and peaks at
    6.31 % against his 5.99 % and the printed "as much as five percent". We run 0.2 to 5.1
    points above his line inside the band, worst at 82 %NG.

    (Remeasured 2026-09-19. This table read 6.17 / 5.82, peaking at 84 %NG, from before the
    minor-lattice recalibration; the peak has moved down the band and the widest gap to his
    line has gone 3.8 -> 5.1 points. Our peak is still under the 8 % ceiling and still within
    2 points of his, which is what the assertions below require.)

    **This is the same disagreement as open question #46, in the report's own currency.**
    #46 closed when Figure 6's y axis stopped being anchored on the figure's caption; before
    that fix this comparison read +7.5 to +20.0 % across the band.
    """
    s = _fig6_series()
    lo, hi = FUEL_OVERESTIMATE_BAND
    ngs = np.arange(lo, hi + 0.01, 0.5)
    ours = np.array(
        [100.0 * (_wf_at(s["ours"], g) / _wf_at(s["ge_status81"], g) - 1.0) for g in ngs]
    )
    theirs = np.array(
        [100.0 * (_wf_at(s["realtime"], g) / _wf_at(s["ge_status81"], g) - 1.0) for g in ngs]
    )

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

    Measured at 88 / 90 / 95 %NG: Ballin -4.69 / -8.06 / -12.80 %, ours -5.07 / -9.25 /
    -13.13 %. Same sign, same order, and monotone in both.
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

    Only testable where **both** GE series have data -- Figure 6's unbalanced-torque series
    starts at NG 83.03 %, so below that there is one bound and not an envelope. Over
    83-100 %NG at whole-percent steps: **Ballin's own line leaves the envelope at 2 of 16
    points and ours at 4.** His two are 84 % -- inside the exception he names -- and 98 %,
    which his sentence does not exempt. So the printed claim does not hold at the top of
    the range for his own line on the digitized data, and this test says so rather than
    asserting something the report does not deliver.
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
