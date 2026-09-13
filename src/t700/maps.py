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
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ._data import DATA_ROOT

DATA_DIR = DATA_ROOT / "maps"
SCHEDULE_DIR = DATA_ROOT / "schedules"
"""The engine's function tables and the fuel control's scheduling functions are the same
*kind* of object -- a printed plot, digitized, linear between knots -- so they share
`Curve`, `SpeedMap` and the loaders. They live in different directories because they
belong to different phases and different appendices."""

_clamps: Counter[str] = Counter()
"""The one piece of mutable module state in the model core, and it is diagnostic only.

CLAUDE.md sanctions exactly one exception to "no global mutable state" -- the thermo
backend interface -- so this is a second, and it is recorded here rather than left to be
rediscovered (the 2026-09-13 engine-physics audit found it).

**Why it is kept.** A clamp is a property of an *evaluation*, not of a component: `f6` is
asked for a fuel-air ratio outside its table by the combustor, which has no way to return
that fact to whoever is interpreting the result three call levels up. Threading a log
through `frame`, `step`, `run`, `trim.solve` and every map call would put diagnostic
plumbing in the signature of every function in the model to carry something no model
equation reads.

**What it costs, bounded.** `frame()` is impure: calling it twice runs the counter up.
Nothing downstream of a map lookup reads the counter, so no model output can depend on
it -- `tests/test_maps.py::test_the_clamp_counter_cannot_influence_any_model_output`
asserts that behaviourally. And clamp accounting is not reentrant, which is what
`clamp_scope` exists to fix: nesting a measurement inside another one used to have the
inner `reset_clamps()` destroy the outer one's count.
"""


class ClampLog:
    """What fell outside a map's data during one scoped measurement."""

    __slots__ = ("counts",)

    def __init__(self) -> None:
        self.counts: dict[str, int] = {}


@contextmanager
def clamp_scope() -> Iterator[ClampLog]:
    """Measure clamps over a block without destroying an enclosing measurement.

    `trim.solve` counts clamps over the whole Newton search and then again at the solution
    alone, and `validation/` counts them over a run that contains trims. Each of those was
    calling `reset_clamps()`, so an outer count was silently zeroed by an inner one.

        with maps.clamp_scope() as log:
            ...
        log.counts        # only what happened inside the block

    The enclosing counter is restored on exit with the block's counts added, so nesting
    composes.
    """
    global _clamps
    outer = _clamps
    inner: Counter[str] = Counter()
    _clamps = inner
    log = ClampLog()
    try:
        yield log
    finally:
        log.counts = dict(inner)
        outer.update(inner)
        _clamps = outer


def clamp_report() -> dict[str, int]:
    """How many evaluations fell outside each map's data range, since the last reset.

    A nonzero count is not automatically wrong -- a transient may legitimately touch the
    edge -- but it must be reported alongside any result, because outside the knots the
    model is extrapolating a constant rather than following the report.
    """
    return dict(_clamps)


def reset_clamps() -> None:
    """Zero the clamp counters. Call at the start of a run you intend to report.

    Prefer `clamp_scope()` where a measurement may be nested inside another.
    """
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
    constant: bool = False
    """Set by `condition` when `Physics.constant` declares the figure a horizontal line.

    A constant has no domain to leave, so evaluations outside `[x[0], x[-1]]` are not
    counted as clamps: there is nothing being extrapolated. `f6` is the only one, and it is
    why "f6 clamped N times" stopped appearing in clamp reports on 2026-09-13. The report
    asked for FAR up to 0.0327 against a table ending at 0.0200, and the answer was the
    same number at both."""

    def __post_init__(self) -> None:
        if self.x.size < 2:
            raise ValueError(f"{self.name}: need at least two knots, got {self.x.size}")
        if not np.all(np.diff(self.x) > 0):
            raise ValueError(f"{self.name}: knots are not strictly increasing in x")
        # `frozen=True` protects the *references*, not what they point at, and the loaders
        # are cached -- so every caller shares one `Curve` per table and
        # `maps.f2().y[:] = ...` silently rewrote the report's data for the whole process.
        # The dangerous spelling was `maps.f2().y *= 1.02`: numpy applies the in-place
        # multiply first and the dataclass guard fires afterwards, so the "protected" form
        # corrupted the table and then raised. Found by the 2026-09-13 code-quality audit.
        self.x.flags.writeable = False
        self.y.flags.writeable = False

    @property
    def domain(self) -> tuple[float, float]:
        return float(self.x[0]), float(self.x[-1])

    def __call__(self, xq):
        xa = np.asarray(xq, dtype=float)
        if not np.all(np.isfinite(xa)):
            raise ValueError(
                f"{self.name} was asked for a non-finite abscissa ({xa}). A NaN reaching a "
                f"function table means something upstream has already failed; propagating "
                f"it silently makes the failure surface hundreds of frames later."
            )
        lo, hi = self.domain
        if not self.constant:
            outside = int(np.count_nonzero((xa < lo) | (xa > hi)))
            if outside:
                _clamps[self.name] += outside
        return np.interp(xa, self.x, self.y)  # np.interp clamps at both ends


@dataclass
class SpeedMap:
    """A two-dimensional map: z = f(x, parameter).

    `f1` is the only one in the engine -- corrected compressor mass flow against static
    pressure ratio, parameterised by corrected gas generator speed.

    ## Interpolation is at constant BETA, along the figure's own beta lines

    Beta is the standard compressor-map coordinate: **0 at the choked end of a speed line,
    1 at surge**, so that a given beta names corresponding points on every speed line no
    matter how different their pressure-ratio ranges are. Interpolating between speed lines
    at constant beta is what performance models do, and it is what this class does.

    **Ballin drew the beta lines.** Figure A1 prints six dotted construction lines, each
    joining the k-th marker of all eleven speed lines, plus a seventh vertical one joining
    their left ends -- and those are beta lines. `tools/digitize_a1.py` tests every
    extracted point against them, and the result is a perfect **11 x 7 grid**: each speed
    line carries markers k = 0..6 at the same seven beta values,

        beta = k / 6,   beta = 0 at the choked left end, beta = 1 at surge,

    with the knee (the end of the flat choked extension) at beta = 1/6. So the report gives
    the beta grid rather than leaving it to be invented, which is why this needs no
    parameterisation of ours.

    ## The seven values are unevenly spaced, and that is correct rather than a defect

    They look alarming: on an average speed line the first interval, beta 0 to 1/6, spans
    **9.15** in pressure ratio and the remaining five together span **1.8**. But the first
    interval is the *choked* part of the line, where corrected flow does not depend on back
    pressure, and it is flat -- measured over all eleven lines, y varies by at most **0.0048
    lbm/s across it, 0.15 %**. The five crowded intervals cover the unchoked run, where y
    varies by **3.0 to 6.4 %**. Ballin put the resolution where the curvature is.

    Checked rather than argued: inserting an extra beta line anywhere inside the choked
    segment, at any skew, moves `f1` by at most **4.3e-16**. Nothing there needs resolving,
    so nothing is lost by not resolving it. Equally spaced beta lines would move resolution
    away from the only part of the line that bends, and would replace a correspondence the
    figure prints with one of ours.

    ## Why seven, and not more

    A real compressor characteristic is smooth, so representing a speed line by six chords
    must in principle lose something -- which is an argument for digitizing intermediate
    beta values off the printed curve. **The printed curve does not have any.** Read at
    600 dpi, Figure A1's speed lines are drawn as *polylines*: straight runs with a visible
    corner at each marker, not splines through them. The plotter drew Ballin's table, and
    the table has seven points per line.

    That makes linear interpolation between the seven a reproduction of his model rather
    than an approximation of it -- a 1988 function-table processor interpolates linearly,
    and the polyline on the page is what that produces. Adding beta values would mean
    inventing intermediate points the report does not contain, and would make this map
    *smoother than the one being replicated*. If that is ever wanted it is a departure to
    be argued and recorded, like Eq. 80's solver, not a digitizing job.

    (The evidence is visual and bounded: the digit markers are about 50 px tall on chords
    of 50 to 83 px, so they occupy most of each segment and a tracer cannot separate the
    line's own bow from glyph strokes. What can be said is that the corners are visible at
    the markers and the segments between them read straight.)

    Evaluation blends the two bracketing speed lines at equal beta -- x and y together --
    and evaluates the query on the blended line. Blending the seven printed beta values
    directly is exact: both coordinates are piecewise linear in beta with the same
    breakpoints on every line, so resampling each line onto a finer beta grid first and
    blending there gives an identical answer. Measured, against 50 points per line:
    **2.4e-16**. A coarser resampling is *worse*, not better -- at 20 points it misses the
    knots and loses 3.2e-4 -- so the seven printed values are both the cheapest and the most
    accurate grid available. See
    `tests/test_maps.py::test_blending_at_the_printed_beta_values_is_exact`.

    ## Interpolating at constant pressure ratio manufactures an extrapolation

    That is what this class did until 2026-09-13. Each speed line ends at its own last
    marker, which is the surge limit for that speed: the 65 % line stops at Ps3/P2 = 3.753
    where the 80 % line reaches 6.671, because a compressor at 65 % corrected speed cannot
    make a pressure ratio of 6. Asking both lines for the same x therefore asks the shorter
    one for a pressure ratio it cannot reach, takes its clamped end value, and blends a real
    number with a fictitious one. Over the Figure 10 chop that happened on **389 frames**,
    every one below 80 %NGc, up to 45 % past the 65 % line's last knot -- and it was
    recorded as open question #58, "f1's 15-point data hole".

    **There is no hole.** Blending along the construction lines, the map's own right edge
    runs 3.753 at 65 %NGc, 4.726 at 70, 5.504 at 74 and 6.671 at 80, while the chop's worst
    query is 5.439 at 74.09 %NGc: **zero of 715 frames fall outside the map.** The
    extrapolation was an artifact of the evaluation, not a property of the data.

    (The 65-to-80 gap is still fifteen points where every other gap is two or three, so the
    blend there is a longer secant than elsewhere. That is a real limit of what the report
    printed, and it is not extrapolation.)

    The previous docstring also claimed the *order* of the two interpolations mattered. It
    did not, under constant-x evaluation -- both orders agreed to 2.9e-16. Under constant-k
    the question does not arise: there is one blended line and one evaluation on it.
    """

    name: str
    params: np.ndarray  # the parameter value of each line, ascending
    lines: list[Curve] = field(default_factory=list)
    source: str = ""

    def __post_init__(self) -> None:
        self.params.flags.writeable = False  # see Curve.__post_init__

    @property
    def param_range(self) -> tuple[float, float]:
        return float(self.params[0]), float(self.params[-1])

    @property
    def beta_grid(self) -> bool:
        """True when every speed line carries the same number of beta values, so knot k of
        one line names the same beta as knot k of the next. `f1` is an 11 x 7 grid, with
        beta = k/6 from the choked end to surge."""
        return len({line.x.size for line in self.lines}) == 1

    def __call__(self, xq: float, pq: float) -> float:
        # `nan < lo or nan > hi` is False, `np.clip(nan, lo, hi)` is nan, and
        # `np.searchsorted(params, nan)` returns len(params) -- so a NaN speed used to
        # select the TOP speed line and return the 100 % value, with an empty clamp
        # report. A corrupted speed yielding maximum compressor flow is the worst
        # available failure mode. Found by the 2026-09-13 numerical-mathematics audit.
        if not (np.isfinite(xq) and np.isfinite(pq)):
            raise ValueError(
                f"{self.name} was asked for a non-finite argument (x={xq}, param={pq}). "
                f"A NaN reaching a function table means something upstream has already "
                f"failed; silently returning the top speed line hides that."
            )
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
        w = 0.0 if p1 == p0 else (p - p0) / (p1 - p0)

        if not self.beta_grid:  # unequal beta grids: no correspondence to interpolate on
            z0 = float(self.lines[j - 1](xq))
            z1 = float(self.lines[j](xq))
            return z0 + w * (z1 - z0)

        # Blend the two speed lines at constant beta -- along the figure's own printed
        # beta lines -- then evaluate the blended line. See the class docstring.
        a, b = self.lines[j - 1], self.lines[j]
        xb = (1.0 - w) * a.x + w * b.x
        yb = (1.0 - w) * a.y + w * b.y
        if xq < xb[0] or xq > xb[-1]:
            _clamps[self.name] += 1
        return float(np.interp(xq, xb, yb))


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
    constant: bool = False
    """The figure draws a horizontal line and the printed spread is read noise.

    Declared only where the *ink* says so: `f6` is Figure A6 [pdf p.61], a single
    horizontal line across a y axis spanning 0.88 to 1.10, with an `x` marker at each
    frame edge and nothing between them. Its two digitized endpoints differ by 7.3e-5
    against a 1-sigma read error of 4.06e-4 each -- **0.13 sigma on the difference**, and
    0.46 px at 200 dpi. The page's four frame edges are fitted at +0.209, +0.155, -0.109
    and +0.150 degrees, four *different* angles, so the frame is a genuine quadrilateral
    and a horizontal line across it does not come back horizontal. That skew is what the
    residual slope is.

    A two-point table is how you write a constant for a function-table processor that
    wants endpoints. Reading a slope out of it is reading the paper, not the model."""
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
"""Largest change physical conditioning made to each table, in data units. Report it
alongside any result -- if a number here is ever large, the digitization is wrong, not
the physics.

The second piece of mutable module state in the core, and unlike `_clamps` it is
**write-once**: the loaders are cached, so each table writes its entry exactly once at
load and nothing touches it again. Named in CLAUDE.md alongside `_clamps` after the
2026-09-13 code-quality audit pointed out that it was a third global nobody had counted."""


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

    **This loads the beta-gridded map**, `f1_compressor_mass_flow_beta.csv`: the same eleven
    speed lines re-expressed on a common beta grid of 56 values, produced by
    `tools/digitize_a1.py` from the raw seven-point extraction. `f1_as_printed()` loads the
    seven-point file.

    Beta is the fraction of a speed line's own pressure-ratio span -- 0 at the choked left
    end, 1 at surge -- **not** the six dotted construction lines Figure A1 prints, although
    those are beta lines and are Ballin's own. The markers carry digitizing noise: as a
    fraction of each line's span the k=1 markers run 0.785, 0.817, 0.831, 0.832, 0.814,
    0.811, 0.855 ... -- non-monotone -- while the spans themselves are smooth. Using the
    markers as the correspondence propagates that jitter into every interpolated value.

    **Using this file is a departure from replication and is recorded as one** (open
    question #60). Figure A1's speed lines are *drawn* as polylines, so the raw file is
    Ballin's table and linear interpolation between those points is his model. The
    characteristic they sample is smooth; the corners are an artifact of plotting.

    Worth, measured against the printed-marker correspondence: Table B.1's rms
    0.2471 -> 0.2439 %, the worst Table 1 NG mode 8.18 -> 7.56 %, transients and the
    zero-extrapolation property unchanged, at 56 stored points instead of 97. Against the
    seven-point original: worst Table 1 NG mode 11.80 -> 7.56 %.

    The fit lives in the digitizer, not here: it needs SciPy, the core may not import it,
    and a derived data file carries its own provenance header and is checked by the
    reproducibility gate like every other.
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


@_cached
def f1_as_printed() -> SpeedMap:
    """`f1` on the seven beta values Figure A1 actually prints. [Fig. A1]

    The faithful transcription: seven markers per speed line, interpolated linearly, which
    is Ballin's own table and what a 1988 function-table processor did with it. `f1()`
    loads the densified extension instead and says why.
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

    **This is a constant, and as of 2026-09-13 it is loaded as one.** Figure A6 draws a
    single horizontal line at 0.985 across a y axis spanning 0.88 to 1.10, with an `x`
    marker at each frame edge and no other ink. The two digitized endpoints, 0.985038 and
    0.984965, differ by 7.3e-5 against a 1-sigma read error of 4.06e-4 apiece -- 0.13 sigma
    on the difference, and 0.46 px at the 200 dpi the figure was read at. The page's four
    frame edges sit at four *different* angles (+0.209, +0.155, -0.109, +0.150 degrees), so
    the frame is a quadrilateral and a printed horizontal line does not come back
    horizontal. The slope was the paper.

    Conditioned to the mean, **0.985002**, which moves each endpoint by 3.65e-5 and is
    recorded in `CONDITIONING`. The CSV is untouched -- it still carries the ink as
    measured, so the provenance chain back to the page is unbroken and the reproducibility
    gate is unaffected.

    Two consequences worth stating. The model asks for FAR up to 0.0327 against a table
    ending at 0.0200 [open question #45], and **that stops being an extrapolation at all**:
    a constant has no domain to leave. And the combustor efficiency stops being a source of
    transient error by construction, which the 2026-09-12 measurement had already shown it
    was not -- relaxing the clamp to linear extrapolation moved the T41 peak by -0.2 degR
    against a 112 degR gap.

    The printed abscissa reads "FUEL-TO-RATIO, FAR"; the word AIR is missing on the page
    (open question #7, closed as a report typo).
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
