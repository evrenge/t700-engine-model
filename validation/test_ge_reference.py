"""Ours, Ballin's model, and the two GE models he validated against [Figs. 6-8].

Every deviation in this project has been quoted against Ballin without a scale. Figures 6,
7 and 8 each plot **three** series: Ballin's real-time model and two General Electric
models -- the status-81 performance standard and an unbalanced-torque variant. Comparing
all three answers the question that gives the others meaning: **how far was Ballin himself
from the standard he was validating against?**

The answer is that he was further from GE than we are from him. His mean distance from GE
status-81 runs 0.45-1.97 % with worst cases to 8.5 %, and from the unbalanced-torque model
1.4-5.3 % with worst cases to 17 %. Ours from Ballin is 0.43-1.51 % mean.

## The result that matters most

On Figure 6, **we are closer to GE than Ballin is**: -0.08 %NG mean against his +0.89 %NG.
That reframes the -1.38 %NG "disagreement" recorded in open question #46 -- it exists
because Ballin's Figure 6 curve sits about 0.9 %NG above GE while we sit on it. Since we
also match Table B.1 to 0.13 %, the reading is that **Table B.1 and GE agree, and Figure 6
is the outlier of the three**, rather than that we drift.

This does not make Figure 6 wrong; it makes it a third source with its own offset, and it
removes the presumption that any gap to it is our error.
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


def test_on_figure_6_we_are_closer_to_ge_than_ballin_is():
    """The finding that reframes open question #46.

    Figure 6 is the one sweep we do not track, at -1.38 %NG. Against GE status-81 -- the
    hardware performance standard, not a model -- ours is -0.08 %NG mean and Ballin's is
    +0.89 %NG. He is the one displaced from GE there, not us, and we also sit on Table B.1
    to 0.13 %. So the gap to Figure 6 is not evidence that we drift.

    Asserted rather than noted because if it ever reverses, #46's whole reading changes.
    """
    d_bal, d_our = _spreads(6, "ge_status81")
    assert abs(d_our.mean()) < abs(d_bal.mean()), (
        f"ours {d_our.mean():+.2f} %NG from GE against Ballin's {d_bal.mean():+.2f} %NG"
    )
    assert abs(d_our.mean()) < 0.5, f"ours should sit on GE: {d_our.mean():+.2f} %NG"


def test_the_two_ge_models_disagree_with_each_other_substantially():
    """The scale that makes every other number here readable.

    GE's own two models differ from one another by far more than we differ from Ballin.
    Whatever "agreement" meant in 1988, it did not mean a few tenths of a per cent.
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
        assert np.abs(d).max() > 2.0, (
            f"fig {fig}: the two GE models now agree to {np.abs(d).max():.2f}; "
            f"if so, the 1988 scatter argument needs re-reading"
        )
