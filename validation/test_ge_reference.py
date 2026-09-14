"""Ours, Ballin's model, and the two GE models he validated against [Figs. 6-8].

**The target is Ballin, and only Ballin.** Closeness to GE is not evidence that we are
right, and this file must not be read as claiming otherwise. We are replicating one
model; agreeing with a different model of the same engine is at best a coincidence and
at worst a symptom -- if our numbers ever sat closer to GE than to Ballin, the first
thing to suspect would be that we had imported something from the wrong source.

Audited 2026-09-12, because that suspicion deserves a check rather than a disclaimer:
**no number in the model derives from the GE or NASA-Lewis series.** The GE traces exist
only under `data/reference/`, nothing in `src/t700/` reads that directory at all
(`tests/test_imports.py` now enforces it), and every transient assertion in
`test_fuel_step.py` targets `fig??_*_model.csv` -- Ballin's own curve -- never
`*_reference.csv`. The GE series has never touched a model value.

## What the GE series legitimately provides: a scale, not a verdict

Every deviation in this project is quoted against Ballin without any sense of how large a
deviation is normal between two careful models of the same engine. Figures 6, 7 and 8
each plot **three** series: Ballin's real-time model and two General Electric models --
the status-81 performance standard and an unbalanced-torque variant. Measuring all three
pairs gives that scale: Ballin's own mean distance from GE status-81 runs 0.45-1.97 %
with worst cases to 8.5 %, and from the unbalanced-torque model 1.4-5.3 % with worst
cases to 17 %. Ours from Ballin is 0.43-1.51 % mean.

That is a ruler, not a grade. It says a 1 % steady-state spread is ordinary among models
of this engine; it says nothing about whether our 1 % is in the right place.

## Where it does bear on a question: #46

On Figure 6 we sit -0.08 %NG from GE where Ballin sits +0.89 %NG. The useful content is
not "we are closer" but that **Figure 6 carries an offset of its own**, which is exactly
the disagreement open question #46 records between Table B.1 and Figure 6. Three sources,
and Table B.1 and GE fall together. That removes the presumption that the gap to Figure 6
is automatically our error -- it does not turn the gap into a credential.
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
NP_RPM = 20895.0


def _read(name: str, xc: str, yc: str):
    rows = list(csv.DictReader(ln for ln in (REF / name).open() if not ln.startswith("#")))
    x = np.array([float(r[xc]) for r in rows])
    y = np.array([float(r[yc]) for r in rows])
    o = np.argsort(x)
    return x[o], y[o]


@lru_cache(maxsize=1)
def _rungs():
    out, guess = [], None
    for pph in np.arange(120.0, 815.0, 10.0):
        r = trim.solve(wf_pps_from_pph(float(pph)), NP_RPM, AMB, guess=guess)
        if r.trustworthy:
            guess = r.state
            out.append((float(pph), r.state))
    return tuple(out)


def _ours(pph: float):
    lower = [st for x, st in _rungs() if x <= pph]
    r = trim.solve(wf_pps_from_pph(pph), NP_RPM, AMB, guess=lower[-1] if lower else None)
    if not r.trustworthy:
        return None
    f = frame(r.state, wf_pps_from_pph(pph), AMB)
    return dict(
        ng=100.0 * r.state.ng_rpm / c.NG_DES,
        shp=shp_from_torque(f.q_pt_ftlbf, NP_RPM),
        ps3=f.ps3_psia,
    )


def _curve(fig: int, xs: np.ndarray):
    """Our model sampled at the given abscissae; returns (mask, values)."""
    if fig == 8:
        pts = [p for p in (_ours(float(v)) for v in np.linspace(150, 800, 80)) if p]
        ng = np.array([p["ng"] for p in pts])
        ps = np.array([p["ps3"] for p in pts])
        o = np.argsort(ng)
        m = (xs >= ng[o].min()) & (xs <= ng[o].max())
        return m, np.interp(xs[m], ng[o], ps[o])
    key = "ng" if fig == 6 else "shp"
    got = [(_ours(float(v)) or {}).get(key) for v in xs]
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
    bx, by = _read(f"{stem}_realtime.csv", xc, yc)
    gx, gy = _read(f"{stem}_{series}.csv", xc, yc)
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
        sx, sy = _read(f"{stem}_ge_status81.csv", xc, yc)
        ux, uy = _read(f"{stem}_ge_unbalanced.csv", xc, yc)
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


def _fig6_series():
    """Fuel flow as a function of NG for all three Figure 6 series, plus ours."""
    out = {}
    for key in ("realtime", "ge_status81", "ge_unbalanced"):
        wf, ng = _read(f"fig06_{key}.csv", "wf_pph", "ng_pct")
        o = np.argsort(ng)
        out[key] = (ng[o], wf[o])
    grid = np.arange(125.0, 851.0, 2.5)
    gx, gn, guess = [], [], None
    for x in grid:
        r = trim.solve(wf_pps_from_pph(float(x)), c.NP_DES, AMB, guess=guess)
        if not r.trustworthy:
            continue
        guess = r.state
        gx.append(float(x))
        gn.append(100.0 * r.state.ng_rpm / c.NG_DES)
    out["ours"] = (np.array(gn), np.array(gx))
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
    | 81 | +0.71 % | +1.08 % |
    | 82 | +5.82 % | +6.13 % |
    | 83 | +5.50 % | +5.91 % |
    | 84 | +2.40 % | +6.17 % |
    | 85 | -0.99 % | +1.07 % |
    | 86 | -1.84 % | -1.34 % |

    So the excursion is real, positive, confined to the band the report names, and peaks at
    6.17 % against his 5.82 % and the printed "as much as five percent". We run 0.4 to 3.8
    points above his line inside the band, worst at 84 %NG.

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
