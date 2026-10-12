"""Every tolerance named in `SCOPE.md` must equal the value the code actually asserts.

CLAUDE.md lists "never silently widen a tolerance" as enforced by "tolerances live in
`SCOPE.md`, never inline in a test file". That was not a mechanism and the description was
not even accurate: every tolerance *is* a literal beside its assertion, and `SCOPE.md`
narrates it. Nothing read `SCOPE.md`, so the two could drift -- and had. On 2026-09-14 the
docs-and-ledger audit found `SCOPE.md` advertising the discrete-map bound as "+/-5 %
(hover, level), +/-25 % (descent)" against a shipped `{1: 8.0, 2: 15.0, 3: 7.0}`: hover and
level **widened** 5 -> 8 and 5 -> 15 with the document still quoting 5, and the loose trim
moved from descent to level, inverting SCOPE's stated rationale as well.

This file is that mechanism. It parses the tolerance table out of `SCOPE.md`, resolves each
`` `module.SYMBOL` `` reference, and requires the numbers the document prints to be exactly
the numbers the symbol holds. Widening is still allowed -- this project has widened a bound
deliberately and recorded why -- but it can no longer happen in only one of the two places.

It deliberately does **not** try to check the prose ("measured 0.62 % worst"), which goes
stale constantly and cannot be parsed reliably. Only the Tolerance column is enforced.

## Three holes closed on 2026-09-19

The first version resolved `_REF.search(basis)` -- the FIRST symbol a row names -- and
checked only that every value the symbol holds is printed in the document. That left three
ways for a bound to escape:

* **A row naming two symbols only had one checked.** Figure 6's row names both an rms and a
  worst bound, and `FIG6_WORST_PCT_NG` was never resolved. It was also written without its
  module prefix, so it could not have matched anyway. Every fully-qualified reference in a
  row is resolved now, and the document spells both out.
* **A number printed in the document that no symbol holds was invisible.** Widening the code
  was caught; inventing a tighter number in the document was not. Both directions are
  checked now, which is what makes the row a statement about the code rather than beside it.
* **The headline tolerance table named no symbols at all**, so the four bounds a reader most
  likely takes as "our declared accuracy" -- NG, gas path, fuel flow, torque -- were
  enforced by nothing. They carry references now and are checked like the rest.

## The fourth, closed 2026-10-12

A row naming *no* symbol was dropped silently: seven Appendix B rows named a test function,
and SCOPE.md advertised +/-99.5 % against an asserted 0.2 % with the suite green. Rows are
now read per tolerance table, and a row in one that names no symbol fails.
"""

from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCOPE = ROOT / "SCOPE.md"

_REF = re.compile(r"`(test_[a-z0-9_]+)\.([A-Z][A-Z0-9_]*)`")
_NUM = re.compile(r"(?<![\d+/.])-?\d+(?:\.\d+)?")
"""A number, with its sign when it has one.

The minus is taken only where it cannot be something else: the lookbehind rejects it after
a digit or a dot (`ratio 0.90-0.96` is a range, not a negative) and after `+` or `/`
(`+/-2.5` is a magnitude). Without that, the parser read `-56 to -50 %` as `{50, 56}` and
a symbol holding `(-56.0, -50.0)` failed against a document that printed it correctly --
which is this test wrong about the code rather than the code wrong about the document."""


_TOLERANCE_HEADERS = {"tolerance", "bound", "steady state"}
"""The second header cell of each tolerance table in SCOPE.md, lowercased."""

EXPECTED_TABLES = 3
EXPECTED_ROWS = 37


def _parse() -> tuple[list[tuple[str, str, str]], list[str], int]:
    """Rows of every tolerance table, the rows among them naming no symbol, and the table
    count. A row whose Tolerance cell starts with "see" points elsewhere and is skipped."""
    rows: list[tuple[str, str, str]] = []
    unnamed: list[str] = []
    tables = 0
    in_table = in_tolerance = False
    for line in SCOPE.read_text().splitlines():
        if not line.startswith("|"):
            in_table = False
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if not in_table:  # a header row
            in_table = True
            in_tolerance = len(cells) > 1 and cells[1].lower() in _TOLERANCE_HEADERS
            tables += in_tolerance
            continue
        if not in_tolerance or cells[0].startswith("---"):
            continue
        # the headline table is | quantity | steady | transient | basis |
        tolerance = " ".join(cells[1:-1])
        if _REF.search(cells[-1]):
            rows.append((cells[0], tolerance, cells[-1]))
        elif not tolerance.lower().startswith("see"):
            unnamed.append(cells[0])
    return rows, unnamed, tables


def _load(module: str):
    """Import the validation module by its real name.

    Not `spec_from_file_location` under an alias: these modules define dataclasses, and
    `dataclasses` resolves annotations through `sys.modules[cls.__module__]`, which an
    unregistered alias is not in.
    """
    path = ROOT / "validation" / f"{module}.py"
    if not path.exists():
        pytest.fail(f"SCOPE.md names `{module}` and validation/{module}.py does not exist")
    if str(ROOT / "validation") not in sys.path:
        sys.path.insert(0, str(ROOT / "validation"))
    return importlib.import_module(module)


def _numbers(value) -> set[float]:
    if isinstance(value, dict):
        return {float(v) for v in value.values()}
    if isinstance(value, (tuple, list, set)):
        return {float(v) for v in value}
    return {float(value)}


ROWS, UNNAMED, TABLES = _parse()


def test_the_tolerance_tables_were_found():
    assert TABLES == EXPECTED_TABLES, f"found {TABLES} SCOPE.md tolerance tables"
    assert len(ROWS) == EXPECTED_ROWS, (
        f"SCOPE.md's tolerance tables have {len(ROWS)} checked rows, expected "
        f"{EXPECTED_ROWS}. Update EXPECTED_ROWS only if a row was added or removed on purpose."
    )


def test_every_tolerance_row_names_a_symbol():
    assert not UNNAMED, (
        f"SCOPE.md tolerance rows naming no `module.SYMBOL`, so nothing checks them: "
        f"{UNNAMED}. Hoist the bound into a named constant and name it in the row."
    )


@pytest.mark.parametrize("row", ROWS, ids=[r[0][:40] for r in ROWS])
def test_every_named_tolerance_matches_the_code(row: tuple[str, str, str]):
    comparison, tolerance, basis = row
    have: set[float] = set()
    named: list[str] = []
    for module, symbol in _REF.findall(basis):
        mod = _load(module)
        assert hasattr(mod, symbol), (
            f"SCOPE.md row {comparison!r} names `{module}.{symbol}` and it does not exist"
        )
        have |= _numbers(getattr(mod, symbol))
        named.append(f"{module}.{symbol}")
    printed = {float(n) for n in _NUM.findall(tolerance)}

    missing = have - printed
    assert not missing, (
        f"{comparison}: {', '.join(named)} holds {sorted(have)} and SCOPE.md's Tolerance "
        f"column prints {sorted(printed)}. The values {sorted(missing)} appear in the code "
        f"and not in the document. Update SCOPE.md in the same commit -- a tolerance that "
        f"moves in one place only is exactly what this test exists to stop."
    )
    unbacked = printed - have
    assert not unbacked, (
        f"{comparison}: SCOPE.md's Tolerance column prints {sorted(unbacked)}, which none "
        f"of {', '.join(named)} holds. A number in the document that no code asserts is a "
        f"tolerance the reader will believe and the suite will not enforce. Either name the "
        f"symbol that holds it or take the number out of the column."
    )
