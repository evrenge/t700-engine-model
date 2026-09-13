"""Function tables: loading the digitized figures and evaluating them.

Every functional relationship in Ballin's engine exists only as a printed plot -- there
are no numeric tables in the report for `f1`..`f10` or `f_hs`. `tools/digitize.py`
turned those plots into the CSVs under `data/maps/`; this module is what the model uses
to read them.

## Interpolation is a decision, not a transcription

The body calls these "function table lookups" [pdf p.37] but **never states the
interpolation scheme** (open question #43). Implemented here as **linear between
breakpoints**, for two reasons:

* the figures are drawn as discrete knots joined by *straight segments*, which is what a
  linear table looks like when plotted; and
* a 1988 real-time model on a 10 ms frame would not have afforded spline evaluation.

Both are inference. If a later reading of the report states the scheme, every value
*between* knots moves.

## Off the end of a map

Clamping, and **counted**. A model that wanders outside its data during a transient is
inventing physics, and the dangerous version of that is the one nobody notices. Every
evaluation outside the digitized range increments a counter, so a run can report how
often it left the data it was built from. See `clamp_report()`.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_DATA = Path(__file__).resolve().parent.parent.parent / "data"
DATA_DIR = _DATA / "maps"
SCHEDULE_DIR = _DATA / "schedules"
"""The engine's function tables and the fuel control's scheduling functions are the same
*kind* of object -- a printed plot, digitized, linear between knots -- so they share
`Curve`, `SpeedMap` and the loaders. They live in different directories because they
belong to different phases and different appendices."""

_clamps: Counter[str] = Counter()


def clamp_report() -> dict[str, int]:
    """How many evaluations fell outside each map's data range, since the last reset.

    A nonzero count is not automatically wrong -- a transient may legitimately touch the
    edge -- but it must be reported alongside any result, because outside the knots the
    model is extrapolating a constant rather than following the report.
    """
    return dict(_clamps)


def reset_clamps() -> None:
    """Zero the clamp counters. Call at the start of a run you intend to report."""
    _clamps.clear()


def _read_csv(path: Path) -> tuple[list[str], dict[str, np.ndarray], list[str]]:
    """Return (column names, columns as arrays, header comment lines)."""
    header = [ln.rstrip("\n") for ln in path.read_text().splitlines() if ln.startswith("#")]
    with path.open() as fh:
        reader = csv.DictReader(ln for ln in fh if not ln.startswith("#"))
        names = list(reader.fieldnames or [])
        rows = list(reader)
    cols = {n: np.array([float(r[n]) for r in rows]) for n in names}
    return names, cols, header


@dataclass(frozen=True)
class Curve:
    """A one-dimensional function table, linear between knots, clamped outside.

    `name` is used in the clamp report, so make it the report's own symbol.
    """

    name: str
    x: np.ndarray
    y: np.ndarray
    source: str = ""

    def __post_init__(self) -> None:
        if self.x.size < 2:
            raise ValueError(f"{self.name}: need at least two knots, got {self.x.size}")
        if not np.all(np.diff(self.x) > 0):
            raise ValueError(f"{self.name}: knots are not strictly increasing in x")

    @property
    def domain(self) -> tuple[float, float]:
        return float(self.x[0]), float(self.x[-1])

    def __call__(self, xq):
        xa = np.asarray(xq, dtype=float)
        lo, hi = self.domain
        outside = int(np.count_nonzero((xa < lo) | (xa > hi)))
        if outside:
            _clamps[self.name] += outside
        return np.interp(xa, self.x, self.y)  # np.interp clamps at both ends


@dataclass
class SpeedMap:
    """A two-dimensional map: z = f(x, parameter).

    `f1` is the only one in the engine -- corrected compressor mass flow against static
    pressure ratio, parameterised by corrected gas generator speed. Each speed line is
    its own `Curve` with its own x range, so evaluation interpolates *along* the two
    bracketing lines first and *between* them second. That order matters: the lines do
    not share breakpoints, so there is no rectangular grid to interpolate on.
    """

    name: str
    params: np.ndarray  # the parameter value of each line, ascending
    lines: list[Curve] = field(default_factory=list)
    source: str = ""

    @property
    def param_range(self) -> tuple[float, float]:
        return float(self.params[0]), float(self.params[-1])

    def __call__(self, xq: float, pq: float) -> float:
        lo, hi = self.param_range
        if pq < lo or pq > hi:
            _clamps[f"{self.name}:parameter"] += 1
        p = float(np.clip(pq, lo, hi))

        j = int(np.searchsorted(self.params, p))
        if j <= 0:
            return float(self.lines[0](xq))
        if j >= len(self.params):
            return float(self.lines[-1](xq))

        p0, p1 = self.params[j - 1], self.params[j]
        z0 = float(self.lines[j - 1](xq))
        z1 = float(self.lines[j](xq))
        w = 0.0 if p1 == p0 else (p - p0) / (p1 - p0)
        return z0 + w * (z1 - z0)


def load_speed_map(
    filename: str,
    name: str,
    param: str,
    xcol: str,
    ycol: str,
    phys: Physics | None = None,
    directory: Path | None = None,
) -> SpeedMap:
    """Load a 2-D map whose rows are grouped by a parameter column.

    `phys` conditions each speed line separately, which is the only sense a per-line
    constraint has: the lines do not share breakpoints. The reported move is the largest
    over all lines.
    """
    path = (directory or DATA_DIR) / filename
    _, cols, header = _read_csv(path)
    src = next((h[len("# source:") :].strip() for h in header if h.startswith("# source:")), "")
    params = np.unique(cols[param])
    lines = []
    moved = 0.0
    for p in params:
        m = cols[param] == p
        order = np.argsort(cols[xcol][m])
        line = Curve(name=f"{name}@{p:g}", x=cols[xcol][m][order], y=cols[ycol][m][order])
        if phys is not None:
            line, d = condition(line, phys)
            moved = max(moved, d)
        lines.append(line)
    if phys is not None:
        CONDITIONING[name] = moved
    return SpeedMap(name=name, params=params, lines=lines, source=src)


# --------------------------------------------------------------------------- physical conditioning


@dataclass(frozen=True)
class Physics:
    """What physics requires of a function, declared per table.

    These plots are pen-plotter drawings of an underlying numeric table, scanned in 2009
    and read back with a pixel counter. What we want is Ballin's *table*, not the ink --
    and where the ink and physics disagree at the level of read noise, physics wins.

    This is emphatically **not** tuning. Nothing here is adjusted to make a result match:
    the constraints are declared from what the quantity *is* (a bleed fraction cannot be
    negative; a choked mass flow cannot rise with back pressure), applied before any
    comparison, and every change is measured and reported. The CSVs under `data/` stay
    exactly as digitized, so the provenance chain back to the page is unbroken.
    """

    nonnegative: bool = False
    at_most_one: bool = False
    at_least_one: bool = False
    monotone: str | None = None  # "inc" or "dec"
    why: str = ""


def _isotonic(y: np.ndarray, increasing: bool) -> np.ndarray:
    """Least-squares monotone fit (pool adjacent violators).

    Chosen over a running min/max because it moves the data as little as possible: a
    single noisy knot is averaged with its neighbour rather than dragging the whole tail.
    """
    v = y.astype(float).copy() if increasing else -y.astype(float)
    n = v.size
    vals = np.zeros(n)
    wts = np.zeros(n)
    idx = -1
    for i in range(n):
        idx += 1
        vals[idx] = v[i]
        wts[idx] = 1.0
        while idx > 0 and vals[idx - 1] > vals[idx]:
            vals[idx - 1] = (wts[idx] * vals[idx] + wts[idx - 1] * vals[idx - 1]) / (
                wts[idx] + wts[idx - 1]
            )
            wts[idx - 1] += wts[idx]
            idx -= 1
    out = np.zeros(n)
    j = 0
    for k in range(idx + 1):
        out[j : j + int(wts[k])] = vals[k]
        j += int(wts[k])
    return out if increasing else -out


def condition(curve: Curve, phys: Physics) -> tuple[Curve, float]:
    """Apply declared physical constraints. Returns the curve and the largest change."""
    y = curve.y.astype(float).copy()
    if phys.monotone in ("inc", "dec"):
        y = _isotonic(y, increasing=(phys.monotone == "inc"))
    if phys.nonnegative:
        y = np.maximum(y, 0.0)
    if phys.at_most_one:
        y = np.minimum(y, 1.0)
    if phys.at_least_one:
        y = np.maximum(y, 1.0)
    moved = float(np.abs(y - curve.y).max())
    return Curve(name=curve.name, x=curve.x, y=y, source=curve.source), moved


CONDITIONING: dict[str, float] = {}
"""Largest change physical conditioning made to each table, in data units. Report it
alongside any result -- if a number here is ever large, the digitization is wrong, not
the physics."""


def load_curve(
    filename: str,
    name: str,
    xcol: str = "x",
    ycol: str = "y",
    phys: Physics | None = None,
    directory: Path | None = None,
) -> Curve:
    """Load a 1-D function table, optionally conditioned. Defaults to data/maps/."""
    path = (directory or DATA_DIR) / filename
    _, cols, header = _read_csv(path)
    src = next((h[len("# source:") :].strip() for h in header if h.startswith("# source:")), "")
    order = np.argsort(cols[xcol])
    c = Curve(name=name, x=cols[xcol][order], y=cols[ycol][order], source=src)
    if phys is not None:
        c, moved = condition(c, phys)
        CONDITIONING[name] = moved
    return c


# --------------------------------------------------------------------------- the engine's tables


def _cached(fn):
    """Load once. These are read-only tables; re-reading them per frame is waste."""
    store: dict[str, object] = {}

    def wrapper():
        if "v" not in store:
            store["v"] = fn()
        return store["v"]

    wrapper.__name__ = fn.__name__
    wrapper.__doc__ = fn.__doc__
    return wrapper


@_cached
def f1() -> SpeedMap:
    """Compressor corrected mass flow, WA2c = f1(Ps3/P2, NGc). [pdf p.22, Eq. 7; Fig. A1]

    The only 2-D map in the engine, and the one the real-time scheme is built around
    avoiding: closing the compressor loop would mean evaluating this several times per
    frame, which is why Ballin opened it.
    """
    return load_speed_map(
        "f1_compressor_mass_flow.csv",
        "f1",
        "ngc_pct",
        "ps3_p2",
        "wa2c_lbm_per_s",
        phys=Physics(
            nonnegative=True,
            monotone="dec",
            why="Along one speed line a compressor passes less corrected flow as it works "
            "against a higher pressure ratio -- that descent toward surge is what the "
            "characteristic is. Two of the 66 digitized segments rise instead, by 4e-4 "
            "and 1.7e-3 lbm/sec on lines spanning 0.17 and 0.26, which is pen width. "
            "The cross-line ordering physics also requires -- more corrected speed "
            "passes more corrected flow at a given ratio -- is NOT imposed, because the "
            "digitized map already satisfies it at every one of 40 sampled ratios; "
            "`tests/test_maps.py` asserts that rather than conditioning it away.",
        ),
    )


@_cached
def f2() -> Curve:
    """Compressor temperature ratio, T3/T2 = f2(Ps3/P2). [Fig. A2]"""
    return load_curve(
        "f2_compressor_temperature.csv",
        "f2",
        phys=Physics(
            at_least_one=True,
            monotone="inc",
            why="A compressor raises total temperature, so T3/T2 > 1, and raising the "
            "pressure ratio must raise it further.",
        ),
    )


@_cached
def f3() -> Curve:
    """Seal-pressurization bleed fraction, B1 = f3(NGc). [Fig. A3]

    Its upper plateau is 0.1091, not 0.110 -- the printed `0.11` tick label on p.58 sits
    ~19 px below its own tick mark. See open question #42.
    """
    return load_curve(
        "f3_seal_bleed_fraction.csv",
        "f3",
        phys=Physics(
            nonnegative=True,
            at_most_one=True,
            monotone="dec",
            why="A seal bleed fraction cannot be negative, and Figure A3 shows one shape: "
            "flat plateau, monotone descent, flat at exactly zero above NGc 89. The "
            "digitized tail wobbles +/-1e-4 around zero, which is read noise on a "
            "pen-plotter line, not a bleed that reopens.",
        ),
    )


@_cached
def f4() -> Curve:
    """Power-turbine-balance bleed fraction. [Fig. A4]"""
    return load_curve(
        "f4_pt_balance_bleed_fraction.csv",
        "f4",
        phys=Physics(
            nonnegative=True,
            at_most_one=True,
            why="A bleed fraction. No monotonicity is imposed: the report draws a 3-point "
            "piecewise line and nothing physical forbids a turning point.",
        ),
    )


@_cached
def f5() -> Curve:
    """Impeller tip leakage and turbine cooling bleed fraction. [Fig. A5]"""
    return load_curve(
        "f5_tip_leak_cooling_bleed_fraction.csv",
        "f5",
        phys=Physics(nonnegative=True, at_most_one=True, why="A bleed fraction."),
    )


@_cached
def f6() -> Curve:
    """Combustor efficiency against fuel-air ratio. [Fig. A6]

    Two knots -- the function is one number, about 0.985, and the read-error bound on it
    is 5e-4 (the scan's vertical bow). The printed abscissa reads "FUEL-TO-RATIO, FAR";
    the word AIR is missing on the page (open question #7, closed as a report typo).
    """
    return load_curve(
        "f6_combustor_efficiency.csv",
        "f6",
        phys=Physics(
            nonnegative=True,
            at_most_one=True,
            why="An efficiency. The report draws a horizontal line, so this is one number "
            "to within its 5e-4 read error.",
        ),
    )


@_cached
def f7() -> Curve:
    """Gas generator turbine energy parameter against P45/P41. [pdf p.24, Eq. 26; Fig. A7]"""
    return load_curve(
        "f7_gg_turbine_energy.csv",
        "f7",
        phys=Physics(
            nonnegative=True,
            monotone="dec",
            why="Less expansion across the turbine extracts less energy, so the parameter "
            "falls as P45/P41 rises.",
        ),
    )


@_cached
def f8() -> Curve:
    """Power turbine energy parameter against P49/P45. [Eq. 32; Fig. A8]

    Goes genuinely negative above P49/P45 ~= 0.83; the printed axis is deliberately
    extended to -5.0. Note the argument is total/total -- `f9` next door takes
    static/total and the two must never be cross-indexed.
    """
    return load_curve(
        "f8_pt_energy.csv",
        "f8",
        phys=Physics(
            monotone="dec",
            why="Same as f7. NOT constrained non-negative: the report deliberately extends "
            "the axis to -5.0 and the function genuinely goes negative above 0.83.",
        ),
    )


@_cached
def f9() -> Curve:
    """Power turbine corrected mass flow against Ps9/P45. [Eq. 33; Fig. A9]

    Argument is **static** station 9 over total station 45 -- different from `f8`.
    """
    return load_curve(
        "f9_pt_mass_flow.csv",
        "f9",
        phys=Physics(
            nonnegative=True,
            monotone="dec",
            why="Corrected mass flow through the power turbine cannot rise as back "
            "pressure rises. The single +1.0e-4 rise sits at the choked end where the "
            "curve is flat, between the two markers the 2026-09-11 re-digitization "
            "reads at 0.372143 and 0.372244 -- the first of which sits ON the left "
            "frame line.",
        ),
    )


@_cached
def f10() -> Curve:
    """Exhaust pressure loss, f10 = P49/Ps9, against NGc. [pdf p.24, Eq. 38; Fig. A10]

    **The figure's axis LABEL is inverted, not its values.** Figure A10's ordinate is
    printed `PS9/P49` while Eq. 38 uses f10 as a multiplier on Ps9 to give P49, which
    requires P49/Ps9. Open question #5 originally resolved this by digitizing as printed
    and taking the reciprocal at the point of use. **That was wrong**, and the physics
    settles it:

    * Eq. 37 sets Ps9 = P_amb, so P49 = P_amb * f10.
    * Expansion through the power turbine requires P45 > P49 > Ps9, so P49/Ps9 > 1.
    * The plotted values are 1.021-1.124, all greater than 1.
    * Taking them as printed gives P49 = 1.02..1.12 x P_amb -- the exhaust total exceeds
      ambient, which is what lets the engine exhaust at all.
    * Taking the reciprocal gives P49 = 0.89..0.98 x P_amb -- exhaust total *below*
      ambient, which is physically impossible for a running engine.

    Corroboration from the neighbouring maps: `f8` takes P49/P45 over 0.3007-0.8505 and
    `f9` takes Ps9/P45 over 0.3003-0.8500. Those spans coincide, so P49 and Ps9 are close
    in magnitude and their ratio is near 1 -- exactly the 1.02-1.12 plotted.

    So the CSV holds f10 directly and there is no inversion anywhere. The defect is the
    printed axis label, and it joins the report's other typesetting errors.
    """
    return load_curve(
        "f10_exhaust_pressure_loss.csv",
        "f10",
        phys=Physics(
            at_least_one=True,
            why="f10 = P49/Ps9 and the exhaust total must exceed the ambient static it "
            "discharges into, or the engine cannot exhaust. See open question #5.",
        ),
    )


@_cached
def f_hs() -> Curve:
    """Station 4.1 heat-sink constant against NGc. [pdf p.26, Eq. 52; Fig. A11]

    Note this is `f_hs`, the lookup, not `f_s` -- the heat-sink transfer function of
    Eqs. 48-53 that contains it. The report prints them as different symbols and an
    earlier draft of the notes conflated them.
    """
    return load_curve(
        "fhs_heat_sink_constant.csv",
        "f_hs",
        phys=Physics(nonnegative=True, why="A heat-sink time-constant coefficient."),
    )
