"""The flight-test days of Figures 11-14, inferred from Ballin's trims [pdf pp.48-52] (#82).

The report prints neither the day nor the gearbox ratio, but each figure's trim prints three
engine numbers: NG, vehicle fuel flow and total torque at the rotor. Our engine, trimmed at his
NG and rotor speed on a given day, has one fuel flow and one torque; his fuel must match ours
directly, his torque through one gearbox ratio shared by all four figures. That is eight
equations; a pressure and a temperature per figure and the ratio are nine unknowns, so one
assumption is needed, and two are fitted:

- **A, standard temperature**: 59 deg F on every day, a pressure per figure and the ratio, five
  unknowns;
- **B, one temperature**: a temperature common to the four days, a pressure per figure and the
  ratio, six unknowns.

On a standard day our engine needs about a fifth more fuel and a fifth more torque at his NG
than his trims print, on every figure, so most of what separates them is pressure. The fits
are evidence about the boundary conditions he flew, never a constant of the model: the
comparisons stay on a standard day until the reading is confirmed (#81).
"""

from __future__ import annotations

import functools

import numpy as np
from scipy.optimize import least_squares

from t700.engine import Ambient
from test_gen_hel import FIGS, his_trim, our_trim, trim_fuel_pph

T_STD_DEGR = 518.67
"""59 deg F, the standard day's temperature [Table B.1, pdf p.67]."""


def _residuals(pressures, ratio, t_degr) -> np.ndarray:
    """Per figure: our trim fuel over his, and our torque times the ratio over his, less one."""
    out = []
    for f, p in zip(FIGS, pressures, strict=True):
        day = Ambient(float(p), float(t_degr))
        fuel = 2.0 * trim_fuel_pph(f, day) / his_trim(f, "fuel")
        torque = 2.0 * ratio * our_trim(f, day).frame.q_pt_ftlbf / (1000.0 * his_trim(f, "torque"))
        out += [fuel - 1.0, torque - 1.0]
    return np.array(out)


@functools.cache
def fit(common_temperature: bool) -> dict:
    """The least-squares day: `pressures` psia per figure, `ratio`, `t_degr`, and the
    `residuals`, fractions; reading B when `common_temperature`, else A."""
    n = len(FIGS)
    if common_temperature:
        start = [13.2] * n + [89.0, 535.0]

        def res(x):
            return _residuals(x[:n], x[n], x[n + 1])

    else:
        start = [12.2] * n + [82.0]

        def res(x):
            return _residuals(x[:n], x[n], T_STD_DEGR)

    s = least_squares(res, start, diff_step=1e-4)
    t = s.x[n + 1] if common_temperature else T_STD_DEGR
    return {"pressures": s.x[:n], "ratio": s.x[n], "t_degr": float(t), "residuals": s.fun}


def days(common_temperature: bool) -> dict[int, Ambient]:
    """Each figure's inferred day."""
    f = fit(common_temperature)
    return {
        fig: Ambient(float(p), f["t_degr"]) for fig, p in zip(FIGS, f["pressures"], strict=True)
    }


def rms_pct(v: np.ndarray) -> float:
    return float(100.0 * np.sqrt(np.mean(np.square(v))))


# ----------------------------------------------------------------------------- the checks

FIT_RMS_PCT = {False: 0.8, True: 0.1}
"""How well each reading reproduces his eight trim numbers, rms %: a characterization, 1.25
times today's value."""


def test_a_standard_temperature_day_fits_his_trims():
    got = rms_pct(fit(False)["residuals"])
    assert got < FIT_RMS_PCT[False], f"reading A misses his trims by {got:.2f} % rms"


def test_one_common_temperature_fits_his_trims():
    got = rms_pct(fit(True)["residuals"])
    assert got < FIT_RMS_PCT[True], f"reading B misses his trims by {got:.2f} % rms"


PRESSURE_SPREAD_PSIA = 0.6
"""The four days' pressures agree within this, under either reading: the figures share one
test series [pdf p.47]. A characterization, 1.25 times today's widest spread."""


def test_the_four_days_agree_in_pressure():
    for common in (False, True):
        p = fit(common)["pressures"]
        assert np.ptp(p) < PRESSURE_SPREAD_PSIA, f"reading {'B' if common else 'A'}: {p}"
