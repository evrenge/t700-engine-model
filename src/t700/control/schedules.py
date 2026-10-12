"""Appendix C's eight scheduling functions, loaded from `data/schedules/`.

Seven are 1-D; `F_HM7` is the only 2-D one (open question #50). Like the engine's
`f1`-`f10` they are digitized plots, linear between knots and clamped outside, so they
reuse the `maps` loaders, clamp counting included. `docs/notes/inventory-appendix-c.md`
traces which block consumes each one.
"""

from __future__ import annotations

from functools import cache

from t700 import maps

_SCHEDULES = maps.SCHEDULE_DIR


@cache
def f_hm1() -> maps.Curve:
    """WFPTP = F_HM1(T2) -- HMU topping line, the ceiling of the droop line.

    [pdf p.95, Fig. C24] Consumed by Fig. C14, p.91. x is T2 in deg R; y is the topping
    fuel-flow parameter Wf/Ps3.
    """
    return maps.load_curve("fhm1_topping_line.csv", "F_HM1", directory=_SCHEDULES)


@cache
def f_hm2() -> maps.Curve:
    """WFPRF = F_HM2(PAS) -- power-available spindle input schedule.

    [pdf p.96, Fig. C25] Consumed by Fig. C16, p.91. Falls monotonically and is subtracted
    from the topping line, so more power-available spindle means more fuel.
    """
    return maps.load_curve("fhm2_power_available.csv", "F_HM2", directory=_SCHEDULES)


@cache
def f_hm3() -> maps.Curve:
    """PNG = F_HM3(XLDSH) -- load-demand compensation, gas generator speed demand.

    [pdf p.97, Fig. C26] Consumed by Fig. C17, p.91. x is the hysteresis-filtered load
    demand spindle angle in deg; y is a speed demand in percent NG.
    """
    return maps.load_curve("fhm3_load_demand_png.csv", "F_HM3", directory=_SCHEDULES)


@cache
def f_hm4() -> maps.Curve:
    """WFQPS3 = F_HM4(XLDSH) -- load-demand compensation, fuel flow delta demand.

    [pdf p.98, Fig. C27] Consumed by Fig. C17, p.91.
    """
    return maps.load_curve("fhm4_load_demand_wfqps3.csv", "F_HM4", directory=_SCHEDULES)


@cache
def f_hm5() -> maps.Curve:
    """WFIRF = F_HM5(T2) -- idle schedule, fuel flow reference.

    [pdf p.99, Fig. C28] Consumed by Fig. C20, p.93.
    """
    return maps.load_curve("fhm5_idle_wfirf.csv", "F_HM5", directory=_SCHEDULES)


@cache
def f_hm6() -> maps.Curve:
    """PCNGI = F_HM6(T2) -- idle schedule, gas generator speed reference.

    [pdf p.100, Fig. C29] Consumed by Fig. C20, p.93. Two knots only: the printed curve
    is a straight line.
    """
    return maps.load_curve("fhm6_idle_pcngi.csv", "F_HM6", directory=_SCHEDULES)


@cache
def f_hm7() -> maps.SpeedMap:
    """WFPAC = F_HM7(PCNGHL, T2) -- HMU acceleration fuel limit. The only 2-D schedule.

    [pdf p.101, Fig. C30] Consumed by Fig. C21, p.93. Seven curves in T2 that cross near
    PCNGHL 85-95, so rank in y is not identity (`tools/digitize_c30.py`, open question #50).

    Called `f_hm7()(pcnghl_pct, t2_degR)`: swept variable first, as in `SpeedMap`.
    """
    return maps.load_speed_map(
        "fhm7_accel_limit.csv",
        "F_HM7",
        param="t2_degR",
        xcol="pcnghl_pct",
        ycol="wfpac",
        directory=_SCHEDULES,
    )


@cache
def f_ec1() -> maps.SpeedMap:
    """TAU45 = F_EC1(W45R, T45L) -- ECU thermocouple sensor time constant, sec.

    [pdf p.94, Fig. C23] Consumed by Fig. C5, p.87. Belongs to the ECU; here because the
    loader is the same.
    """
    return maps.load_speed_map(
        "fec1_thermocouple_tau.csv",
        "F_EC1",
        param="t45l",
        xcol="w45r",
        ycol="tau45_s",
        directory=_SCHEDULES,
    )
