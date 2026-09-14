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


def _table_rows() -> list[tuple[str, str, str]]:
    """Every `| comparison | tolerance | basis |` row that names a code symbol."""
    out = []
    for line in SCOPE.read_text().splitlines():
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 3 or cells[0].startswith("---"):
            continue
        if _REF.search(cells[2]):
            out.append((cells[0], cells[1], cells[2]))
    return out


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


ROWS = _table_rows()


def test_the_tolerance_table_was_found():
    assert len(ROWS) >= 8, f"only {len(ROWS)} SCOPE.md tolerance rows name a symbol"


@pytest.mark.parametrize("row", ROWS, ids=[r[0][:40] for r in ROWS])
def test_every_named_tolerance_matches_the_code(row: tuple[str, str, str]):
    comparison, tolerance, basis = row
    module, symbol = _REF.search(basis).groups()
    mod = _load(module)
    assert hasattr(mod, symbol), (
        f"SCOPE.md row {comparison!r} names `{module}.{symbol}` and it does not exist"
    )
    have = _numbers(getattr(mod, symbol))
    printed = {float(n) for n in _NUM.findall(tolerance)}
    missing = have - printed
    assert not missing, (
        f"{comparison}: {module}.{symbol} holds {sorted(have)} and SCOPE.md's Tolerance "
        f"column prints {sorted(printed)}. The values {sorted(missing)} appear in the code "
        f"and not in the document. Update SCOPE.md in the same commit -- a tolerance that "
        f"moves in one place only is exactly what this test exists to stop."
    )
