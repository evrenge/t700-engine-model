"""Every tolerance named in `SCOPE.md` must equal the value the code actually asserts.

Parses each tolerance table of `SCOPE.md`, resolves every `` `module.SYMBOL` `` a row names,
and requires the Tolerance column to print exactly the numbers those symbols hold, in both
directions. Widening is allowed; doing it in only one of the two places is not. A row that
names no symbol fails. The prose ("measured ... worst") is deliberately not checked.
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

The lookbehind rejects a minus after a digit or dot (`0.90-0.96` is a range) and after `+`
or `/` (`+/-2.5` is a magnitude), so `-56 to -50 %` reads as `{-56, -50}`."""


_SPLIT = re.compile(r"(?<!\\)\|")
"""A cell boundary: markdown splits a table row on unescaped pipes only (`\\|` is a literal
one, as in a `<!--rj:...\\|abs-->` marker, see `tests/test_live_numbers.py`)."""

_TOLERANCE_HEADERS = {"tolerance", "bound", "steady state"}
"""The second header cell of each tolerance table in SCOPE.md, lowercased."""

EXPECTED_TABLES = 3
EXPECTED_ROWS = 43


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
        cells = [c.strip() for c in _SPLIT.split(line.strip().strip("|"))]
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
    """Import the validation module by its real name, not under a `spec_from_file_location`
    alias: dataclasses resolve annotations through `sys.modules[cls.__module__]`."""
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
        f"{comparison}: {', '.join(named)} holds {sorted(have)}, SCOPE.md prints "
        f"{sorted(printed)}; {sorted(missing)} missing there. Update SCOPE.md in the same commit."
    )
    unbacked = printed - have
    assert not unbacked, (
        f"{comparison}: SCOPE.md prints {sorted(unbacked)}, which none of {', '.join(named)} "
        f"holds; name the symbol that holds it or take the number out."
    )
