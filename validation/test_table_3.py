"""Table 3 [pdf p.44]: the NASA-Lewis test engine, Ballin's model and ours, at Table 2's six
test conditions [pdf p.43].

"For the six cases, the inlet and exhaust conditions, fuel flow, and load conditions were
specified so that internal states could be compared. [...] The model power turbine speed was
allowed to vary in order to achieve trim with the test load torque" [pdf p.39]. So each case
is a trim at Table 2's P2, T2 and Wf, with NP solved until the power turbine's torque equals
the printed load torque, and the exhaust set by the printed P49.

A characterization, not a pass/fail target (#34): Ballin's model column was made with
"compressor-mass-flow and turbine-energy functions derived from this test engine [...] in
place of the standard functions" [pdf p.39], which the report never prints. Ours runs the
standard functions of Appendix A. The difference between his column and ours is therefore
mostly those maps, and nothing here may be tuned toward either column.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import pytest
from scipy.optimize import brentq

from t700 import constants as c
from t700 import maps, trim
from t700.engine import Ambient, frame
from t700.units import wf_pps_from_pph


@dataclass(frozen=True)
class Case:
    """One row of Table 2 [pdf p.43], read from the page image."""

    wf_pph: float
    p2_psia: float
    t2_degR: float
    p49_psia: float
    load_ftlbf: float


TABLE_2 = {
    1: Case(140.1, 14.37, 516.7, 14.37, 30.1),
    2: Case(297.2, 14.17, 515.6, 14.43, 90.1),
    3: Case(372.0, 14.16, 508.3, 14.46, 148.3),
    4: Case(458.4, 14.09, 508.0, 14.60, 206.5),
    5: Case(560.6, 14.02, 507.2, 14.63, 274.3),
    6: Case(694.4, 13.92, 507.2, 14.72, 360.8),
}
"""Table 2 [pdf p.43]: Wf lbm/hr, P2 psia, T2 degR, P49 psia, load torque ft*lbf."""

QUANTITIES = ("NG", "NP", "WA2", "P3", "T3", "T45")
"""Table 3's columns [pdf p.44]: NG and NP in %, WA2 lbm/sec, P3 psia, T3 and T45 degR."""

TABLE_3 = {
    "engine": {
        1: (65.9, 52.6, 3.20, 58.0, 832.0, 1413.0),
        2: (84.7, 95.7, 5.16, 113.1, 1026.0, 1577.0),
        3: (87.7, 95.7, 6.16, 139.0, 1081.0, 1626.0),
        4: (90.4, 95.7, 6.92, 161.1, 1127.0, 1731.0),
        5: (92.6, 95.7, 7.66, 184.8, 1173.0, 1838.0),
        6: (95.9, 95.7, 8.50, 211.9, 1228.0, 1974.0),
    },
    "ballin": {
        1: (70.0, 83.8, 3.59, 63.4, 855.9, 1372.0),
        2: (86.4, 100.6, 5.52, 119.0, 1043.0, 1555.0),
        3: (88.1, 97.8, 6.27, 141.7, 1084.0, 1629.0),
        4: (90.7, 97.3, 7.00, 163.3, 1130.0, 1734.0),
        5: (92.8, 98.2, 7.73, 186.4, 1178.0, 1855.0),
        6: (96.1, 99.4, 8.53, 213.8, 1233.0, 2009.0),
    },
}
"""Table 3 [pdf p.44], read from the page image: the experimental engine's rows and the
"Real-time model" rows, in `QUANTITIES` order. Cases 2-6 print T3 and T45 with a trailing
point and no decimal digit; that is the report's own formatting."""


@dataclass(frozen=True)
class Ours:
    values: tuple[float, ...]
    """In `QUANTITIES` order and units."""
    np_rpm: float
    ps9_psia: float
    trustworthy: bool
    clamps: dict[str, int]


def _ambient(case: Case, ps9: float) -> Ambient:
    return Ambient(case.p2_psia, case.t2_degR, ps9_psia=ps9)


def _trim(case: Case, np_rpm: float, ps9: float):
    return trim.solve(wf_pps_from_pph(case.wf_pph), np_rpm, _ambient(case, ps9))


@functools.cache
def run(n: int) -> Ours:
    """Table 2 case `n`, trimmed with NP free and P49 as printed.

    Eq. 38 makes P49 = Ps9 * f10(NGc), so the printed P49 fixes Ps9 once NGc is known: a
    fixed point, a few passes. Inside it, NP is the root of Q_PT - load torque, which falls
    monotonically with NP (Eq. 41).
    """
    case = TABLE_2[n]
    ps9 = case.p49_psia / 1.05
    np_rpm = c.NP_DES
    for _ in range(50):

        def excess(x: float, ps9: float = ps9) -> float:
            return _trim(case, x, ps9).frame.q_pt_ftlbf - case.load_ftlbf

        np_rpm = brentq(excess, 0.3 * c.NP_DES, 1.3 * c.NP_DES, xtol=1e-9)
        r = _trim(case, np_rpm, ps9)
        new = case.p49_psia / float(maps.f10()(r.frame.ngc_pct))
        if abs(new - ps9) < 1e-12:
            break
        ps9 = new
    else:
        raise RuntimeError(f"case {n}: P49 fixed point did not converge")
    r = _trim(case, np_rpm, ps9)
    with maps.clamp_scope() as clamps:  # the solution's own reads, not the search's
        frame(r.state, wf_pps_from_pph(case.wf_pph), _ambient(case, ps9))
    f = r.frame
    values = (
        100.0 * r.state.ng_rpm / c.NG_DES,
        100.0 * np_rpm / c.NP_DES,
        f.wa2_pps,
        r.state.p3_psia,
        f.t3_degR,
        f.t45_degR,
    )
    return Ours(values, np_rpm, ps9, r.trustworthy, dict(clamps))


def deviation(n: int, column: str) -> dict[str, float]:
    """Ours against one printed column of case `n`, % of the printed value."""
    printed = TABLE_3[column][n]
    return {
        q: 100.0 * (o - p) / abs(p)
        for q, o, p in zip(QUANTITIES, run(n).values, printed, strict=True)
    }


@pytest.mark.parametrize("n", sorted(TABLE_2))
def test_each_case_trims_on_the_printed_conditions(n: int):
    """The printed load torque and P49 are met, and the trim stands on the maps' data."""
    case, ours = TABLE_2[n], run(n)
    r = _trim(case, ours.np_rpm, ours.ps9_psia)
    assert ours.trustworthy, f"case {n}: trim not trustworthy, clamps {ours.clamps}"
    assert abs(r.frame.q_pt_ftlbf - case.load_ftlbf) < 1e-6
    p49 = ours.ps9_psia * float(maps.f10()(r.frame.ngc_pct))
    assert abs(p49 - case.p49_psia) < 1e-9


CHARACTERIZED = (2, 3, 4, 5, 6)
"""Case 1 is reported, not bounded: "the model is not valid for such low power turbine
speeds" [pdf p.39], and its NP sits far from design."""

TABLE_3_CHAR_PCT = {"NG": 4.0, "NP": 10.0, "WA2": 17.0, "P3": 14.0, "T3": 4.0, "T45": 11.0}
"""Ours against Ballin's "Real-time model" column, cases 2-6, per quantity: a
characterization, the smallest whole percent at least 1.25 times today's worst. Not an
accuracy claim: his column ran the test engine's own maps (see the module docstring)."""


@pytest.mark.parametrize("n", CHARACTERIZED)
def test_ours_against_ballins_column_is_characterized(n: int):
    dev = deviation(n, "ballin")
    bad = {q: round(d, 2) for q, d in dev.items() if abs(d) >= TABLE_3_CHAR_PCT[q]}
    assert not bad, f"case {n}: {bad} % against Table 3's model column"


def test_the_standard_maps_make_a_stronger_engine_than_the_prototype():
    """At every case 2-6, against Ballin's column, ours passes more air at a higher P3, runs a
    cooler T45 and settles at a higher NP for the same fuel and load torque: the standard
    maps describe a better engine than the prototype's, which is the report's own reason for
    swapping them [pdf p.39]. A sign flip here would mean the difference is no longer the
    maps alone."""
    higher, lower = ("NG", "NP", "WA2", "P3", "T3"), ("T45",)
    for n in CHARACTERIZED:
        dev = deviation(n, "ballin")
        assert all(dev[q] > 0 for q in higher), f"case {n}: {dev}"
        assert all(dev[q] < 0 for q in lower), f"case {n}: {dev}"
