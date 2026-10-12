"""The heat-sink time constants Appendix B was generated with, recovered from its own numbers.

Open questions #31 and #73, diagnosed 2026-10-12. The 3-DOF fuel column is a near-cancellation,
`b(T41) = B0 + A(T41,T41) d`, and its T41 element runs 0.87 / 0.98 / 0.85 of printed because
our heat-sink time constants are short. Eqs. 50-53 [pdf p.26], `TC_T41 = 0.29` [p.55] and the
f_hs plateau at 7.50 [Fig. A11, p.66] are all as printed; the printed 3- and 6-DOF matrices
agree with each other on time constants those printed inputs do not give. Expressed as
inputs, Appendix B behaves as if TC_T41 were about 0.314 and the plateau about 7.07.

Writing `q = A(T41,NG) / (1 + tau1 A(NG,NG))`, the Chen form gives, per trim:

    B0       = d / tau1 + tau1 Jw q                  Jw = b(NG) - A(NG,T41) d
    tau2     = D + tau1 c3,   D = (c3 - 1) / (A(T41,T41) - tau1 A(NG,T41) q)

so tau1 comes from a printed 3-DOF block alone, and tau2 needs only `c3 = dT41_ns/dT41`,
taken from our model (it is 0.06; a 10 % error in it moves tau2 by about 1 %).
"""

from __future__ import annotations

import pytest
from scipy.optimize import brentq

from t700 import appendix_b as ab
from t700 import linear, trim
from t700.engine import frame
from t700.linear import DOF
from t700.units import wf_pps_from_pph
from test_all_dof_models import AMB, FIG, NP_RPM, TRIM, WF_PPH, _model

TAU1_RATIO = (0.84, 0.88)
TAU2_RATIO = (0.91, 0.95)
"""Ours over Appendix B's. Measured 0.862 / 0.855 / 0.852 and 0.924 / 0.925 / 0.935."""

SWAP_ROW_TOL_PCT = 2.0
"""6-DOF T41 row with the recovered taus. Measured within 0.4 % (hover, level), 1.8 %
(descent); with our own taus it is 7 % low."""


def _parts(m):
    ng, t41 = m.states.index("NG"), m.states.index("T41")
    a, b, d = m.A, m.b, m.d[t41]
    return d, b[t41] - a[t41, t41] * d, a[ng, ng], a[ng, t41], a[t41, ng], a[t41, t41], b[ng]


def _tau1(m) -> float:
    d, b0, j11, j13, a31, _a33, b1 = _parts(m)
    jw = b1 - j13 * d

    def f(x):
        return d / x + x * jw * a31 / (1 + x * j11) - b0

    # two further roots near 0.6-1.4 s give tau2 under 0.1 s; the bracket excludes them
    return brentq(f, 1.5, 6.0)


def _k(m, tau1):
    _d, _b0, j11, j13, a31, a33, _b1 = _parts(m)
    return a33 - tau1 * j13 * a31 / (1 + tau1 * j11)


def _tau2(m, tau1, c3):
    return (c3 - 1) / _k(m, tau1) + tau1 * c3


def _ours(trim_no):
    wf = wf_pps_from_pph(WF_PPH[trim_no])
    taus = linear._heat_sink_taus(frame(trim.solve(wf, NP_RPM, AMB).state, wf, AMB))
    return taus, _model(DOF.THREE, trim_no)


def _recovered(trim_no):
    (t1, t2), ours = _ours(trim_no)
    k = _k(ours, t1)
    c3 = (1 + t2 * k) / (1 + t1 * k)
    printed = ab.load(FIG[(DOF.THREE, trim_no)])
    tau1 = _tau1(printed)
    return tau1, _tau2(printed, tau1, c3)


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_the_recovery_returns_our_own_tau1_from_our_matrices(trim_no: int):
    (t1, _t2), ours = _ours(trim_no)
    assert _tau1(ours) == pytest.approx(t1, rel=1e-6)


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_appendix_b_was_generated_with_longer_time_constants(trim_no: int):
    (t1, t2), _ = _ours(trim_no)
    b1, b2 = _recovered(trim_no)
    r1, r2 = t1 / b1, t2 / b2
    assert TAU1_RATIO[0] < r1 < TAU1_RATIO[1], f"{TRIM[trim_no]} tau1 ours/printed {r1:.3f}"
    assert TAU2_RATIO[0] < r2 < TAU2_RATIO[1], f"{TRIM[trim_no]} tau2 ours/printed {r2:.3f}"


@pytest.mark.parametrize("trim_no", [1, 2, 3])
def test_those_time_constants_make_the_six_dof_t41_row_agree(trim_no: int, monkeypatch):
    """Recovered from the 3-DOF block, checked on the 6-DOF one: the two printed blocks
    agree on the heat sink, and the time constants are the whole of the 6-DOF T41 row's
    7 % deficit."""
    taus = _recovered(trim_no)
    monkeypatch.setattr(linear, "_heat_sink_taus", lambda _f: taus)
    ours, printed = _model(DOF.SIX, trim_no), ab.load(FIG[(DOF.SIX, trim_no)])
    t41, ng = printed.states.index("T41"), printed.states.index("NG")
    for j, s in enumerate(printed.states):
        if j == ng or printed.A[t41, j] == 0.0:
            continue
        dev = (ours.A[t41, j] / printed.A[t41, j] - 1.0) * 100.0
        assert abs(dev) < SWAP_ROW_TOL_PCT, f"{TRIM[trim_no]} A(T41,{s}) {dev:+.2f} %"
