"""Function tables: loading the digitized figures and evaluating them.

Every function in Ballin's engine (`f1`..`f10`, `f_hs`) exists only as a printed plot. The
digitizers under `tools/` turned them into the CSVs in `data/maps/`
(`docs/notes/digitizing-recipe.md`); this module reads and evaluates them.

The report calls these "function table lookups" [pdf p.37] but never states the
interpolation scheme (open question #43). Linear between breakpoints is an inference: the
figures are knots joined by straight segments, and a 1988 real-time model could not afford
splines. Off the end of a table the value is clamped and counted (`clamp_report()`), so a
run can report how often it left its data.
"""

from __future__ import annotations

import csv
from bisect import bisect_left, bisect_right
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path

import numpy as np

from ._data import DATA_ROOT

_INF = float("inf")
_NINF = float("-inf")

DATA_DIR = DATA_ROOT / "maps"
SCHEDULE_DIR = DATA_ROOT / "schedules"
"""The fuel control's scheduling functions: the same kind of object as the engine's tables,
so they share `Curve`, `SpeedMap` and the loaders."""

_clamps: Counter[str] = Counter()
"""The clamp counter: mutable module state, sanctioned in CLAUDE.md, diagnostic only.

A clamp is a property of an evaluation, which a component cannot hand back to whoever reads
the result three calls up; threading a log through every signature would carry something
no equation reads. The cost is that `frame()` is impure. No model equation reads the
counter, so no model output can depend on it
(`tests/test_maps.py::test_the_clamp_counter_cannot_influence_any_model_output`).
"""


@contextmanager
def clamp_scope() -> Iterator[dict[str, int]]:
    """Measure clamps over a block without destroying an enclosing measurement.

        with maps.clamp_scope() as counts:
            ...
        counts        # only what happened inside the block, filled in on exit

    The enclosing counter is restored on exit with the block's counts added, so nesting
    composes, where `reset_clamps()` would zero the outer count (a run that counts
    contains trims that count).
    """
    global _clamps
    outer = _clamps
    inner: Counter[str] = Counter()
    _clamps = inner
    counts: dict[str, int] = {}
    try:
        yield counts
    finally:
        counts.update(inner)
        outer.update(inner)
        _clamps = outer


def clamp_report() -> dict[str, int]:
    """How many evaluations fell outside each map's data range, since the last reset.

    Not automatically wrong, but report it with any result: outside the knots the model
    extrapolates a constant rather than following the report.
    """
    return dict(_clamps)


def reset_clamps() -> None:
    """Zero the clamp counters. Call at the start of a run you intend to report.

    Prefer `clamp_scope()` where a measurement may be nested inside another.
    """
    _clamps.clear()


_PROPAGATES = "propagating it silently makes the failure surface hundreds of frames later"


def _non_finite(name: str, what: str, consequence: str) -> ValueError:
    """The error a function table raises when asked for NaN or an infinity."""
    return ValueError(
        f"{name} was asked for a non-finite {what}. A NaN reaching a function table means "
        f"something upstream has already failed; {consequence}."
    )


def _read_csv(path: Path) -> tuple[dict[str, np.ndarray], list[str], str]:
    """Return (columns as arrays, in file order; header comment lines; the `# source:` line)."""
    header = [ln.rstrip("\n") for ln in path.read_text().splitlines() if ln.startswith("#")]
    with path.open() as fh:
        reader = csv.DictReader(ln for ln in fh if not ln.startswith("#"))
        names = list(reader.fieldnames or [])
        rows = list(reader)
    cols = {n: np.array([float(r[n]) for r in rows]) for n in names}
    src = next((h[len("# source:") :].strip() for h in header if h.startswith("# source:")), "")
    return cols, header, src


@dataclass(frozen=True)
class Curve:
    """A one-dimensional function table, linear between knots, clamped outside.

    `name` is used in the clamp report, so make it the report's own symbol.
    """

    name: str
    x: np.ndarray
    y: np.ndarray
    source: str = ""
    constant: bool = False
    """Set by `condition` when `Physics.constant` declares the figure a horizontal line.

    A constant has no domain to leave, so evaluations outside `[x[0], x[-1]]` are not
    counted as clamps. `f6` is the only one."""

    def __post_init__(self) -> None:
        if self.x.size < 2:
            raise ValueError(f"{self.name}: need at least two knots, got {self.x.size}")
        if not np.all(np.diff(self.x) > 0):
            raise ValueError(f"{self.name}: knots are not strictly increasing in x")
        # `frozen=True` protects the references, not the arrays, and the cached loaders
        # share one `Curve` per table: lock the data, or `maps.f2().y *= 1.02` rewrites the
        # table for the whole process (and only then raises).
        self.x.flags.writeable = False
        self.y.flags.writeable = False
        # Scalar-lookup caches: every model call site asks for one abscissa, and the NumPy
        # round-trip costs far more than the arithmetic on Python floats. Copies of the
        # read-only `x` and `y`, so they cannot drift.
        set_ = object.__setattr__
        set_(self, "_xl", self.x.tolist())
        set_(self, "_yl", self.y.tolist())
        set_(self, "_n", int(self.x.size))
        set_(self, "_lo", float(self.x[0]))
        set_(self, "_hi", float(self.x[-1]))
        # For `realtime._p45_exact`: Eq. 80's root is unique only where the table cannot rise.
        set_(self, "non_increasing", bool(np.all(np.diff(self.y) <= 0.0)))

    @property
    def domain(self) -> tuple[float, float]:
        return self._lo, self._hi

    def __call__(self, xq):
        if isinstance(xq, float | int):
            return self._at(xq)
        xa = np.asarray(xq, dtype=float)
        if not np.all(np.isfinite(xa)):
            raise _non_finite(self.name, f"abscissa ({xa})", _PROPAGATES)
        lo, hi = self._lo, self._hi
        if not self.constant:
            outside = int(np.count_nonzero((xa < lo) | (xa > hi)))
            if outside:
                _clamps[self.name] += outside
        return np.interp(xa, self.x, self.y)  # np.interp clamps at both ends

    def _at(self, x: float) -> float:
        """One abscissa, clamped at both ends -- `np.interp`'s arithmetic, bit for bit.

        Reproduces NumPy's `arr_interp` for one scalar: the bracket is the largest `j` with
        `x[j] <= x`, the ends return `y[0]` and `y[-1]`, a query exactly on a knot returns
        that knot's ordinate, and the interior is `slope*(x - x[j]) + y[j]` with the slope
        formed per call. Keep that form; it is what makes the bits match.
        """
        if x != x or x == _INF or x == _NINF:
            raise _non_finite(self.name, f"abscissa ({x})", _PROPAGATES)
        if not self.constant and (x < self._lo or x > self._hi):
            _clamps[self.name] += 1
        xl = self._xl
        yl = self._yl
        j = bisect_right(xl, x) - 1
        if j < 0:
            return yl[0]
        if j >= self._n - 1:
            return yl[-1]
        xj = xl[j]
        if xj == x:
            return yl[j]
        yj = yl[j]
        return (yl[j + 1] - yj) / (xl[j + 1] - xj) * (x - xj) + yj


@dataclass
class SpeedMap:
    """A two-dimensional map: z = f(x, parameter), one `Curve` per parameter value.

    `f1` is the only one in the engine: corrected compressor mass flow against static
    pressure ratio, one line per corrected gas generator speed.

    Lines are blended at constant beta, the standard compressor-map coordinate (0 at the
    choked end of a speed line, 1 at surge): knot k names the same beta on every line, so
    the two bracketing lines are blended knot by knot, x and y together, and the query is
    evaluated on the blended line. Never blend at constant x: each line ends at its own
    surge point, so the shorter line would answer with its clamped end value and blend in a
    fictitious extrapolation (open question #58).

    Figure A1 prints the beta lines itself: dotted construction lines join the k-th marker
    of all eleven speed lines, beta = k/6 (`f1_as_printed`). The speed lines are drawn as
    polylines, so linear interpolation between markers reproduces Ballin's table. Blending
    at the knots is exact
    (`tests/test_maps.py::test_blending_at_the_printed_beta_values_is_exact`).
    """

    name: str
    params: np.ndarray  # the parameter value of each line, ascending
    lines: list[Curve] = field(default_factory=list)
    source: str = ""

    def __post_init__(self) -> None:
        self.params.flags.writeable = False  # see Curve.__post_init__
        # Scalar evaluation, as in `Curve.__post_init__`: cache what is fixed at load.
        self._pl: list[float] = self.params.tolist()
        self._plo: float = self._pl[0]
        self._phi: float = self._pl[-1]
        self._param_key: str = f"{self.name}:parameter"
        self._beta_grid: bool = len({line.x.size for line in self.lines}) == 1

    @property
    def param_range(self) -> tuple[float, float]:
        return self._plo, self._phi

    @property
    def beta_grid(self) -> bool:
        """True when every speed line carries the same number of knots, so knot k names the
        same beta on every line."""
        return self._beta_grid

    def __call__(self, xq: float, pq: float) -> float:
        # A NaN passes every range test and bisects past the last line, so it would
        # silently return the top speed line with an empty clamp report: raise instead.
        if xq != xq or xq == _INF or xq == _NINF or pq != pq or pq == _INF or pq == _NINF:
            raise _non_finite(
                self.name,
                f"argument (x={xq}, param={pq})",
                "silently returning the top speed line hides that",
            )
        lo = self._plo
        hi = self._phi
        if pq < lo or pq > hi:
            _clamps[self._param_key] += 1
        p = lo if pq < lo else (hi if pq > hi else float(pq))

        pl = self._pl
        j = bisect_left(pl, p)
        if j <= 0:
            return float(self.lines[0](xq))
        if j >= len(pl):
            return float(self.lines[-1](xq))

        p0 = pl[j - 1]
        p1 = pl[j]
        w = 0.0 if p1 == p0 else (p - p0) / (p1 - p0)

        if not self._beta_grid:  # unequal beta grids: no correspondence to interpolate on
            z0 = float(self.lines[j - 1](xq))
            z1 = float(self.lines[j](xq))
            return z0 + w * (z1 - z0)

        # Blend the two speed lines at constant beta, then evaluate the blended line (see
        # the class docstring). The blended line is never materialized: `xb[k] = om*a.x[k]
        # + w*b.x[k]` is formed on demand inside the same bisection `np.interp` runs, with
        # the array form's arithmetic element for element, so the bits are identical.
        a, b = self.lines[j - 1], self.lines[j]
        om = 1.0 - w
        ax, ay, bx, by = a._xl, a._yl, b._xl, b._yl
        n = a._n
        if xq < om * ax[0] + w * bx[0] or xq > om * ax[n - 1] + w * bx[n - 1]:
            _clamps[self.name] += 1
        # NumPy's own bisection: the largest k with xb[k] <= xq, or -1 / n-1 off the ends.
        klo, khi = 0, n
        while klo < khi:
            mid = (klo + khi) >> 1
            if xq >= om * ax[mid] + w * bx[mid]:
                klo = mid + 1
            else:
                khi = mid
        k = klo - 1
        if k < 0:
            return om * ay[0] + w * by[0]
        if k >= n - 1:
            return om * ay[n - 1] + w * by[n - 1]
        xk = om * ax[k] + w * bx[k]
        yk = om * ay[k] + w * by[k]
        if xk == xq:
            return yk
        xk1 = om * ax[k + 1] + w * bx[k + 1]
        yk1 = om * ay[k + 1] + w * by[k + 1]
        return (yk1 - yk) / (xk1 - xk) * (xq - xk) + yk


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
    cols, _, src = _read_csv(path)
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

    The plots are pen-plotter drawings of Ballin's numeric table; where ink and physics
    disagree at the level of read noise, physics wins. Not tuning: constraints come from
    what the quantity is (a bleed fraction cannot be negative), are applied before any
    comparison, and every change is recorded in `CONDITIONING`. The CSVs stay as digitized.
    """

    nonnegative: bool = False
    at_most_one: bool = False
    at_least_one: bool = False
    monotone: str | None = None  # "inc" or "dec"
    constant: bool = False
    """The figure draws a horizontal line and the printed spread is read noise.

    Declared only where the ink says so: `f6`, Figure A6 [pdf p.61], one horizontal line
    whose digitized endpoints differ by far less than the read error; the residual slope is
    the page's skew. A two-point table is how a function-table processor writes a constant."""
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
    if phys.constant:
        y = np.full_like(y, float(np.mean(y)))
    if phys.monotone in ("inc", "dec"):
        y = _isotonic(y, increasing=(phys.monotone == "inc"))
    if phys.nonnegative:
        y = np.maximum(y, 0.0)
    if phys.at_most_one:
        y = np.minimum(y, 1.0)
    if phys.at_least_one:
        y = np.maximum(y, 1.0)
    moved = float(np.abs(y - curve.y).max())
    return (
        Curve(name=curve.name, x=curve.x, y=y, source=curve.source, constant=phys.constant),
        moved,
    )


CONDITIONING: dict[str, float] = {}
"""Largest change physical conditioning made to each table, in data units. Report it with
any result: a large number means the digitization is wrong, not the physics.

Mutable module state sanctioned in CLAUDE.md, but write-once: the loaders are cached, so
each table writes its entry once at load."""


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
    cols, _, src = _read_csv(path)
    order = np.argsort(cols[xcol])
    c = Curve(name=name, x=cols[xcol][order], y=cols[ycol][order], source=src)
    if phys is not None:
        c, moved = condition(c, phys)
        CONDITIONING[name] = moved
    return c


# --------------------------------------------------------------------------- the engine's tables


# Each loader is `functools.cache`d: the read-only tables load once.


@cache
def f1() -> SpeedMap:
    """Compressor corrected mass flow, WA2c = f1(Ps3/P2, NGc). [pdf p.22, Eq. 7; Fig. A1]

    Loads the beta-gridded map, `f1_compressor_mass_flow_beta.csv`: the eleven speed lines
    re-expressed by `tools/digitize_a1.py` on a common grid of 56 beta values, beta being
    the fraction of each line's own pressure-ratio span. A departure from replication
    (open question #60): the printed markers carry digitizing jitter that would propagate
    into every interpolated value. `f1_as_printed()` loads Ballin's seven-point table. The
    fit lives in the digitizer because it needs SciPy.
    """
    return load_speed_map(
        "f1_compressor_mass_flow_beta.csv",
        "f1",
        "ngc_pct",
        "ps3_p2",
        "wa2c_lbm_per_s",
        phys=Physics(
            nonnegative=True,
            monotone="dec",
            why="Along one speed line a compressor passes less corrected flow as it works "
            "against a higher pressure ratio. Applied to the densified map for the same "
            "reason it is applied to the raw one; the shape-preserving fit cannot "
            "introduce a rise the anchors do not have, so this conditions nothing.",
        ),
    )


@cache
def f1_as_printed() -> SpeedMap:
    """`f1` on the seven beta values Figure A1 prints. [Fig. A1]

    The faithful transcription: Ballin's own table, interpolated linearly as a 1988
    function-table processor did. `f1()` loads the densified extension and says why.
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


@cache
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


@cache
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


@cache
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


@cache
def f5() -> Curve:
    """Impeller tip leakage and turbine cooling bleed fraction. [Fig. A5]"""
    return load_curve(
        "f5_tip_leak_cooling_bleed_fraction.csv",
        "f5",
        phys=Physics(nonnegative=True, at_most_one=True, why="A bleed fraction."),
    )


@cache
def f6() -> Curve:
    """Combustor efficiency against fuel-air ratio. [Fig. A6]

    A constant (open question #62): Figure A6 draws one horizontal line, and the two
    digitized endpoints differ by far less than the read error; the slope was the page's
    skew. Conditioned to the mean, recorded in `CONDITIONING`; the CSV keeps the ink. A
    constant has no domain to leave, so FAR beyond the table (open question #45) is not an
    extrapolation.

    The printed abscissa reads "FUEL-TO-RATIO, FAR"; AIR is missing on the page (open
    question #7).
    """
    return load_curve(
        "f6_combustor_efficiency.csv",
        "f6",
        phys=Physics(
            nonnegative=True,
            at_most_one=True,
            constant=True,
            why="An efficiency, and Figure A6 draws it as a single horizontal line. The "
            "two endpoints differ by 0.13 sigma of the read error, which is the page's "
            "own skew rather than a slope.",
        ),
    )


@cache
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


@cache
def f8() -> Curve:
    """Power turbine energy parameter against P49/P45. [Eq. 32; Fig. A8]

    Goes genuinely negative above P49/P45 ~= 0.83; the printed axis extends to -5.0. The
    argument is total/total, while `f9` takes static/total: never cross-index them.
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


@cache
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


@cache
def f10() -> Curve:
    """Exhaust pressure loss, f10 = P49/Ps9, against NGc. [pdf p.24, Eq. 38; Fig. A10]

    The axis label is inverted, not the values (open question #5). Figure A10 prints
    `PS9/P49`, but Eq. 38 multiplies Ps9 (= P_amb, Eq. 37) by f10 to give P49, and the
    plotted values, all above 1, put the exhaust total above ambient as a running engine
    requires; the reciprocal would put it below. The CSV holds f10 directly; nothing is
    inverted.
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


@cache
def f_hs() -> Curve:
    """Station 4.1 heat-sink constant against NGc. [pdf p.26, Eq. 52; Fig. A11]

    `f_hs` is the lookup, not `f_s`, the heat-sink transfer function of Eqs. 48-53 that
    contains it; the report prints different symbols.
    """
    return load_curve(
        "fhs_heat_sink_constant.csv",
        "f_hs",
        phys=Physics(nonnegative=True, why="A heat-sink time-constant coefficient."),
    )
