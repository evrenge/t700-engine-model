"""Appendix C's eight scheduling functions, loaded from `data/schedules/`.

The fuel control has no equations in the report -- it is one constants table and thirty
figures, of which eight are function plots. These are those eight. Seven are 1-D; `F_HM7`
is the only 2-D one, and it was the hardest thing in Phase 1 (open question #50).

They are the same *kind* of object as the engine's `f1`-`f10`: a printed plot, digitized,
linear between knots, clamped outside. So they reuse `maps.Curve`, `maps.SpeedMap` and the
loaders unchanged, including the clamp counting -- a schedule driven off the end of its
printed range is exactly as interesting as a compressor map driven off the end of its, and
`maps.clamp_report()` reports both.

Naming follows the report: `F_HM1` is the HMU's topping line, `F_EC1` the ECU's
thermocouple time constant. `docs/notes/inventory-appendix-c.md` traces which block
consumes each one.
"""

from __future__ import annotations

from t700 import maps

_SCHEDULES = maps.SCHEDULE_DIR


@maps._cached
def f_hm1() -> maps.Curve:
    """WFPTP = F_HM1(T2) -- HMU topping line, the ceiling of the droop line.

    [pdf p.95, Fig. C24] Consumed by Fig. C14, p.91. x is T2 in deg R; y is the topping
    fuel-flow parameter Wf/Ps3.
    """
    return maps.load_curve("fhm1_topping_line.csv", "F_HM1", directory=_SCHEDULES)


@maps._cached
def f_hm2() -> maps.Curve:
    """WFPRF = F_HM2(PAS) -- power-available spindle input schedule.

    [pdf p.96, Fig. C25] Consumed by Fig. C16, p.91. Falls monotonically from 10.4 at
    PAS = 28 deg to -0.1 at 120 deg: it is *subtracted* from the topping line, so more
    power-available spindle means less subtraction and more fuel.
    """
    return maps.load_curve("fhm2_power_available.csv", "F_HM2", directory=_SCHEDULES)


@maps._cached
def f_hm3() -> maps.Curve:
    """PNG = F_HM3(XLDSH) -- load-demand compensation, gas generator speed demand.

    [pdf p.97, Fig. C26] Consumed by Fig. C17, p.91. x is the hysteresis-filtered load
    demand spindle angle in deg; y is a speed demand in percent NG.
    """
    return maps.load_curve("fhm3_load_demand_png.csv", "F_HM3", directory=_SCHEDULES)


@maps._cached
def f_hm4() -> maps.Curve:
    """WFQPS3 = F_HM4(XLDSH) -- load-demand compensation, fuel flow delta demand.

    [pdf p.98, Fig. C27] Consumed by Fig. C17, p.91.
    """
    return maps.load_curve("fhm4_load_demand_wfqps3.csv", "F_HM4", directory=_SCHEDULES)


@maps._cached
def f_hm5() -> maps.Curve:
    """WFIRF = F_HM5(T2) -- idle schedule, fuel flow reference.

    [pdf p.99, Fig. C28] Consumed by Fig. C20, p.93.
    """
    return maps.load_curve("fhm5_idle_wfirf.csv", "F_HM5", directory=_SCHEDULES)


@maps._cached
def f_hm6() -> maps.Curve:
    """PCNGI = F_HM6(T2) -- idle schedule, gas generator speed reference.

    [pdf p.100, Fig. C29] Consumed by Fig. C20, p.93. Two knots only: the printed curve
    is a straight line.
    """
    return maps.load_curve("fhm6_idle_pcngi.csv", "F_HM6", directory=_SCHEDULES)


@maps._cached
def f_hm7() -> maps.SpeedMap:
    """WFPAC = F_HM7(PCNGHL, T2) -- HMU acceleration fuel limit. The only 2-D schedule.

    [pdf p.101, Fig. C30] Consumed by Fig. C21, p.93. Seven curves in T2, and they
    **cross**: after the bundle at PCNGHL 85-95 the order top to bottom is 4, 3, 5, 6, 7.
    That is why it needed its own digitizer (`tools/digitize_c30.py`, open question #50)
    and why rank in y is not identity here.

    Argument order follows `SpeedMap`: the swept variable first, the curve parameter
    second, so this is called `f_hm7()(pcnghl_pct, t2_degR)`.
    """
    return maps.load_speed_map(
        "fhm7_accel_limit.csv",
        "F_HM7",
        param="t2_degR",
        xcol="pcnghl_pct",
        ycol="wfpac",
        directory=_SCHEDULES,
    )


@maps._cached
def f_ec1() -> maps.SpeedMap:
    """TAU45 = F_EC1(W45R, T45L) -- ECU thermocouple sensor time constant, sec.

    [pdf p.94, Fig. C23] Consumed by Fig. C5, p.87. Belongs to the ECU, not the HMU; it
    is here because it is one of the eight and the loader is the same.
    """
    return maps.load_speed_map(
        "fec1_thermocouple_tau.csv",
        "F_EC1",
        param="t45l",
        xcol="w45r",
        ycol="tau45_s",
        directory=_SCHEDULES,
    )
