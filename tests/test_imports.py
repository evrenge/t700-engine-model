"""Enforce the architecture constraints in CLAUDE.md.

Until this file existed, three of the rules in CLAUDE.md's enforcement table were marked
`advisory` -- which is to say, they were wishes. These tests make three of them real:

  * the model core imports NumPy and the standard library, nothing else
  * gas property evaluation happens only inside `t700.thermo`
  * the core does not read a clock or a random number generator

The provenance rule -- every number cites a report page -- stays advisory. It cannot be
mechanised, which is precisely why it is the one to be most careful about.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent.parent / "src" / "t700"

ALLOWED_THIRD_PARTY = {"numpy"}
"""The complete list. Adding to it is an architecture decision, not a convenience."""

BANNED_IN_CORE = {"scipy", "matplotlib", "pandas", "numba", "sympy", "PIL"}
"""Analysis and plotting live in tools/ and validation/. Numba is deferred (SCOPE.md)."""

NONDETERMINISTIC = {"random", "time", "datetime", "secrets", "uuid", "os.urandom"}
"""Same input must give bit-identical output, run to run and machine to machine."""


def core_modules() -> list[Path]:
    return sorted(SRC.rglob("*.py"))


def top_level_imports(path: Path) -> set[str]:
    """Every top-level package name imported by one module."""
    tree = ast.parse(path.read_text(), filename=str(path))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                names.add(alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level:  # relative import, stays inside the package
                continue
            if node.module:
                names.add(node.module.split(".")[0])
    return names


def test_core_has_files():
    """Guard against the other tests passing vacuously on an empty tree."""
    assert core_modules(), f"no python files under {SRC}"


@pytest.mark.parametrize("path", core_modules(), ids=lambda p: p.name)
def test_core_imports_numpy_and_stdlib_only(path: Path):
    stdlib = set(sys.stdlib_module_names)
    for name in top_level_imports(path):
        if name in ("t700", "__future__") or name in stdlib:
            continue
        assert name in ALLOWED_THIRD_PARTY, (
            f"{path.relative_to(SRC.parent.parent)} imports {name!r}. The model core "
            f"takes NumPy and the stdlib only; analysis code belongs in tools/ or "
            f"validation/. See CLAUDE.md."
        )


@pytest.mark.parametrize("path", core_modules(), ids=lambda p: p.name)
def test_core_does_not_import_analysis_libraries(path: Path):
    found = top_level_imports(path) & BANNED_IN_CORE
    assert not found, (
        f"{path.name} imports {sorted(found)}. These make the core non-portable or "
        f"non-deterministic; they belong in tools/ or validation/."
    )


@pytest.mark.parametrize("path", core_modules(), ids=lambda p: p.name)
def test_core_is_deterministic_by_construction(path: Path):
    found = top_level_imports(path) & NONDETERMINISTIC
    assert not found, (
        f"{path.name} imports {sorted(found)}. The core must give bit-identical output "
        f"for identical input; no clocks and no RNG."
    )


GAS_PROPERTY_CONSTANTS = (
    "K_H2",
    "K_H3_1",
    "K_H3_2",
    "K_H41_1",
    "K_H41_2",
    "K_H45",
    "K_T41_1",
    "K_T41_2",
    "K_T45_1",
    "K_T45_2",
    "K_T49_1",
    "K_T49_2",
    "K_TH41_1",
    "K_TH41_2",
    "K_TH45_1",
    "K_TH45_2",
)


@pytest.mark.parametrize("path", core_modules(), ids=lambda p: p.name)
def test_gas_properties_only_inside_thermo(path: Path):
    """No enthalpy or temperature-fit arithmetic outside `t700.thermo`.

    This is the rule that lets the backend be swapped without touching component code.
    A component that reaches for `K_H41_1` directly has silently welded the linear fits
    into itself.

    `constants.py` is exempt: it is where the values are declared. `thermo/` is exempt:
    it is the door.
    """
    rel = path.relative_to(SRC)
    if rel.parts[0] == "thermo" or rel.name == "constants.py":
        return

    source = path.read_text()
    used = [k for k in GAS_PROPERTY_CONSTANTS if k in source]
    assert not used, (
        f"{path.name} references gas property constants {used}. Go through "
        f"t700.thermo instead -- see CLAUDE.md, 'All gas property evaluation goes "
        f"through t700.thermo'."
    )
