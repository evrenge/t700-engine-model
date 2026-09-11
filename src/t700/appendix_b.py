"""Appendix B's twelve linear models, as printed [pdf pp.67-76].

Appendix B is Table B.1 plus twelve bare matrix figures -- **no eigenvalues, no time
constants, no transfer functions and no mode descriptions appear anywhere in it**. It is
the richest comparison surface in the report: 297 printed numbers against the 12 our
eigenvalue check uses.

These are **transcribed**, not digitized. Every value is a printed number read from the
page image, so they carry no read error in the sense a plot digitization does -- only the
report's own 4-significant-figure rounding, which matters in exactly one place (see
`ill_conditioned`).

## The layout

```
                          heat sink off        heat sink on
  full volume dynamics    5-DOF  B2,B4,B6      6-DOF  B8,B10,B12
  volume dyn. approx.     2-DOF  B1,B3,B5      3-DOF  B7,B9,B11
```

Trim conditions are Table B.1's three, in order: hover, level 80 kt, descent 80 kt.

## Two traps

**The 6-DOF matrices do not support an eigenvalue comparison.** Row 6 of `A` is row 6 of
`F1^-1 F2`, a difference of terms of order 2e5 producing a result of order 1e3, so the
printed four digits carry roughly +-50 absolute in `A(6,3)` and `A(6,4)`. Perturbing those
two *within their printed rounding* moves the slow mode anywhere from -0.99 to +7.3 /sec.
As transcribed, B8 is unstable at +3.07 /sec. That is lost precision, not a claim by the
report, and `ill_conditioned` flags it.

**The 3-DOF state ordering is inferred**, not printed -- `{NG, NP, T41}`, from three
pieces of evidence set out in `inventory-appendix-b.md` 2.4. `states_are_printed` says so
per figure. Anything depending on that ordering should say it is doing so.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent.parent.parent / "data" / "linear"


@dataclass(frozen=True)
class Reference:
    """One printed linear model."""

    figure: str
    """`"B2"`, as the report labels it."""

    page: int
    dof: int
    trim: int
    """1, 2 or 3 -- Table B.1's hover, level 80 kt, descent 80 kt."""

    heat_sink: bool
    states: tuple[str, ...]
    states_are_printed: bool
    """False only for the 3-DOF models, whose ordering is inferred."""

    A: np.ndarray
    b: np.ndarray
    d: np.ndarray | None
    """Feedthrough from Wf [Eq. 67]. Present only on the heat-sink models; `C = [I]` is
    printed literally in all six and is not stored."""

    @property
    def ill_conditioned(self) -> bool:
        """True for the 6-DOF models. Compare these element by element only."""
        return self.dof == 6

    def __repr__(self) -> str:
        hs = "heat sink" if self.heat_sink else "no heat sink"
        return f"<{self.figure}: {self.dof}-DOF trim {self.trim}, {hs}, pdf p.{self.page}>"


def _header(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        if ":" in line:
            k, _, v = line[1:].partition(":")
            k = k.strip()
            if k and not k.startswith(" "):
                out.setdefault(k, v.strip())
    return out


@cache
def load(figure: str | int) -> Reference:
    """Load one figure. `load("B2")`, `load("b02")` and `load(2)` are the same thing."""
    n = int(str(figure).lower().lstrip("b"))
    if not 1 <= n <= 12:
        raise ValueError(f"Appendix B has figures B1..B12; got {figure!r}")
    path = DATA / f"b{n:02d}.csv"
    head = _header(path)
    states = tuple(head["states"].split()[0].split(","))
    k = len(states)

    A = np.zeros((k, k))
    b = np.zeros(k)
    d = np.zeros(k)
    idx = {s: i for i, s in enumerate(states)}
    with path.open() as fh:
        for row in csv.DictReader(ln for ln in fh if not ln.startswith("#")):
            v = float(row["value"])
            if row["matrix"] == "A":
                A[idx[row["row"]], idx[row["col"]]] = v
            elif row["matrix"] == "b":
                b[idx[row["row"]]] = v
            else:
                d[idx[row["row"]]] = v

    model = head["model"]
    dof = int(model.split("-DOF")[0])
    A.setflags(write=False)
    b.setflags(write=False)
    d.setflags(write=False)
    return Reference(
        figure=f"B{n}",
        page=int(head["source"].split("pdf p.")[1].split(",")[0]),
        dof=dof,
        trim=int(model.split("trim condition ")[1].split(" ")[0]),
        heat_sink=head["heat_sink"].startswith("on"),
        states=states,
        states_are_printed="INFERRED" not in head["states"],
        A=A,
        b=b,
        d=d if head["heat_sink"].startswith("on") else None,
    )


def find(dof: int, trim: int) -> Reference:
    """The figure for one model at one trim, e.g. `find(5, 1)` is B2."""
    for n in range(1, 13):
        r = load(n)
        if r.dof == dof and r.trim == trim:
            return r
    raise ValueError(f"no Appendix B figure for {dof}-DOF at trim {trim}")


def all_references() -> list[Reference]:
    """All twelve, in figure order."""
    return [load(n) for n in range(1, 13)]
