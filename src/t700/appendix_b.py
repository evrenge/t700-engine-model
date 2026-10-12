"""Appendix B's twelve linear models, as printed [pdf pp.67-76].

Table B.1 plus twelve bare matrix figures, with no eigenvalues, time constants or mode
descriptions. Transcribed from the page images, so the only error is the report's own
4-significant-figure rounding, which matters in one place (`ill_conditioned`).

```
                          heat sink off        heat sink on
  full volume dynamics    5-DOF  B2,B4,B6      6-DOF  B8,B10,B12
  volume dyn. approx.     2-DOF  B1,B3,B5      3-DOF  B7,B9,B11
```

Trim conditions are Table B.1's three, in order: hover, level 80 kt, descent 80 kt.

Two traps:

* The 6-DOF matrices do not support an eigenvalue comparison. Row 6 of `A` is row 6 of
  `F1^-1 F2`, terms of order 2e5 cancelling to order 1e3, so the printed four digits carry
  about +-50 in `A(6,3)` and `A(6,4)`; within that rounding the slow mode moves from -0.99
  to +7.3 /sec (B8 as transcribed is unstable). Lost precision, not a claim by the report.
* The 3-DOF state ordering `{NG, NP, T41}` is inferred, not printed
  (`inventory-appendix-b.md` 2.4); `states_are_printed` says so per figure.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np

from ._data import DATA_ROOT

DATA = DATA_ROOT / "linear"


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
            if k:
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
