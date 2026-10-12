"""The heat-sink time constants Appendix B was generated with, recovered from its own numbers.

Open questions #31 and #73, diagnosed 2026-10-12. Eqs. 50-53 [pdf p.26], `TC_T41 = 0.29`
[p.55] and the f_hs plateau at 7.50 [Fig. A11, p.66] are all as printed, and they give time
constants 14 % (tau1) and 9 % (tau2) shorter than the ones Appendix B was generated with.
Expressed as inputs, the whole difference is TC_T41: Appendix B behaves as if it were 0.320,
and as if f_hs were the plotted one to 0.4 %. The replication keeps the printed 0.29; this
file pins the inconsistency.

With T41 held fixed during extraction [pdf p.31], `dT41_ns/dT41` is zero, and writing
`q = A(T41,NG) / (1 + tau1 A(NG,NG))` the Chen form gives, per trim:

    B0    = d / tau1 + tau1 Jw q                  Jw = b(NG) - A(NG,T41) d
    tau2  = -1 / (A(T41,T41) - tau1 A(NG,T41) q)

so both come from one printed 3-DOF block alone. Recovered from B7/B9/B11, tau1/tau2 is
0.6712 / 0.6446 / 0.6202; derivative-ambiguity.md's feedthrough estimator on the 6-DOF
blocks gives 0.6712 / 0.6450 / 0.6182.
"""

from __future__ import annotations

import pytest
from scipy.optimize import brentq

from t700 import appendix_b as ab
from t700 import constants as c
from t700 import linear, maps, trim
from t700.engine import frame
from t700.linear import DOF
from t700.units import wf_pps_from_pph
from test_all_dof_models import AMB, DQ_REQ_DNP, FIG, NP_RPM, TRIM, WF_PPH, _model

TAU1_RATIO = (0.84, 0.88)
TAU2_RATIO = (0.89, 0.93)
"""Ours over Appendix B's. Measured 0.862 / 0.855 / 0.852 and 0.907 / 0.906 / 0.910."""

SIX_DOF_SWAP_TOL_PCT = 0.5
"""The 6-DOF fuel column -- gas-path `b` and `d(T41)` -- with the taus recovered from the
3-DOF. Measured worst 0.37 %; with our own taus it runs 5-6 % low."""

IMPLIED_TC_T41 = (0.315, 0.325)
IMPLIED_F_HS_TOL_PCT = 1.0
"""Eqs. 51 and 53 inverted at our trims for the recovered taus. Measured TC_T41 0.3198 /
0.3201 / 0.3188 against the printed 0.29; f_hs within 0.12 / 0.16 / 0.37 % of the plot."""

THREE_DOF_SWAP_TOL_PCT = 2.0
"""The 3-DOF T41 row -- `A(T41,T41)`, `b(T41)`, `d` -- with those taus, hover and descent, on
the tangent. Measured worst 1.7 % (4.0 % at hover with the shipped +-2 % step, whose secant
through the settled pressures moves `b(T41)`). Level is not asserted: A(T41,T41) 1.07 and
b(T41) 1.15 there, open under #73."""


def _parts(m):
    ng, t41 = m.states.index("NG"), m.states.index("T41")
    a, b, d = m.A, m.b, m.d[t41]
    return d, b[t41] - a[t41, t41] * d, a[ng, ng], a[ng, t41], a[t41, ng], a[t41, t41], b[ng]


def _taus(m) -> tuple[float, float]:
    d, b0, j11, j13, a31, a33, b1 = _parts(m)
    jw = b1 - j13 * d

    def f(x):
        return d / x + x * jw * a31 / (1 + x * j11) - b0

    # two further roots near 0.6-1.4 s give tau2 under 0.1 s; the bracket excludes them
    tau1 = brentq(f, 1.5, 6.0)
    return tau1, -1.0 / (a33 - tau1 * j13 * a31 / (1 + tau1 * j11))


def _ours(trim_no):
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    return linear._heat_sink_taus(frame(trim.solve(wf, NP_RPM, AMB).state, wf, AMB))


def _recovered(trim_no):
    return _taus(ab.load(FIG[(DOF.THREE, trim_no)]))


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_recovery_returns_our_own_taus_from_our_matrices(trim_no: int):
    assert _taus(_model(DOF.THREE, trim_no)) == pytest.approx(_ours(trim_no), rel=1e-6)


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_appendix_b_was_generated_with_longer_time_constants(trim_no: int):
    (t1, t2), (b1, b2) = _ours(trim_no), _recovered(trim_no)
    r1, r2 = t1 / b1, t2 / b2
    assert TAU1_RATIO[0] < r1 < TAU1_RATIO[1], f"{TRIM[trim_no]} tau1 ours/printed {r1:.3f}"
    assert TAU2_RATIO[0] < r2 < TAU2_RATIO[1], f"{TRIM[trim_no]} tau2 ours/printed {r2:.3f}"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_difference_is_tc_t41_alone(trim_no: int):
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    f = frame(trim.solve(wf, NP_RPM, AMB).state, wf, AMB)
    tau1, tau2 = _recovered(trim_no)
    tc = tau2 * f.w41_pps**0.8 / f.t41_ns_degR**0.5  # (51) inverted
    f_hs = (tau2 - tau1) * f.w41_pps  # (52)-(53) inverted
    plotted = float(maps.f_hs()(f.ngc_pct))
    assert IMPLIED_TC_T41[0] < tc < IMPLIED_TC_T41[1], f"{TRIM[trim_no]} TC_T41 {tc:.4f}"
    dev = (f_hs / plotted - 1.0) * 100.0
    assert abs(dev) < IMPLIED_F_HS_TOL_PCT, f"{TRIM[trim_no]} f_hs {dev:+.2f} % off the plot"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_those_time_constants_reproduce_the_six_dof_fuel_column(trim_no: int, monkeypatch):
    """Recovered from the 3-DOF block, checked on the 6-DOF one: the two printed blocks
    agree on the heat sink, and the time constants are the whole of the 6-DOF column's
    deficit."""
    taus = _recovered(trim_no)
    monkeypatch.setattr(linear, "_heat_sink_taus", lambda _f: taus)
    ours, printed = _model(DOF.SIX, trim_no), ab.load(FIG[(DOF.SIX, trim_no)])
    t41 = printed.states.index("T41")
    pairs = [(f"b({s})", ours.b[i], printed.b[i]) for i, s in enumerate(printed.states)]
    pairs = [p for p in pairs if p[2] != 0.0 and p[0] != "b(T41)"]
    pairs.append(("d(T41)", ours.d[t41], printed.d[t41]))
    for name, o, p in pairs:
        dev = (o / p - 1.0) * 100.0
        assert abs(dev) < SIX_DOF_SWAP_TOL_PCT, f"{TRIM[trim_no]} {name} {dev:+.2f} %"


@pytest.mark.parametrize("trim_no", [1, 3])
def test_those_time_constants_reproduce_the_three_dof_t41_row(trim_no: int, monkeypatch):
    taus = _recovered(trim_no)
    monkeypatch.setattr(linear, "_heat_sink_taus", lambda _f: taus)
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    ours = linear.extract(
        trim.solve(wf, NP_RPM, AMB),
        wf,
        DOF.THREE,
        AMB,
        j_load=c.J_LOAD_UH60A,
        dq_req_dnp=DQ_REQ_DNP[trim_no],
        rel_step=1e-5,
    )
    printed = ab.load(FIG[(DOF.THREE, trim_no)])
    t = printed.states.index("T41")
    for name, o, p in (
        ("A(T41,T41)", ours.A[t, t], printed.A[t, t]),
        ("b(T41)", ours.b[t], printed.b[t]),
        ("d(T41)", ours.d[t], printed.d[t]),
    ):
        dev = (o / p - 1.0) * 100.0
        assert abs(dev) < THREE_DOF_SWAP_TOL_PCT, f"{TRIM[trim_no]} {name} {dev:+.2f} %"
