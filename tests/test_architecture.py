"""Enforce the architecture constraints in CLAUDE.md.

  * the model core imports NumPy and the standard library, nothing else
  * gas property evaluation happens only inside `t700.thermo`
  * the core does not read a clock or a random number generator

An AST walk over `import` statements answers "what does this module import"; the rules are
about "what can this module reach", so each route around the import walk is closed
separately below.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent.parent / "src" / "t700"

ALLOWED_THIRD_PARTY = {"numpy"}
"""The complete list. Adding to it is an architecture decision, not a convenience."""

FORBIDDEN_MODULES = (
    "numpy.random",
    "os.urandom",
    "os.getrandom",
    "secrets",
    "random",
    "time",
    "datetime",
    "uuid",
    "timeit",
    "hashlib",
)
"""Dotted module paths the core may not import by any spelling: no clocks and no RNG.

`hashlib` because a digest over floats is a plausible route to a "deterministic" value that
is not; `timeit` because `timeit.default_timer` is a clock that does not import `time`.
Matched on the full dotted path: `top_level_imports` reduces `from numpy import random` to
the allowed `"numpy"`.
"""

FORBIDDEN_ATTRIBUTES = {
    # RNG, reachable through the one third-party package the core is allowed
    ("np", "random"): "NumPy's RNG is still an RNG; the core must be reproducible",
    ("numpy", "random"): "NumPy's RNG is still an RNG; the core must be reproducible",
    # clocks and entropy that ride in on stdlib modules the core legitimately uses
    ("os", "urandom"): "entropy source",
    ("os", "times"): "a clock",
    ("os", "getrandom"): "entropy source",
    ("time", "time"): "a clock",
    ("time", "perf_counter"): "a clock",
    ("time", "monotonic"): "a clock",
    # dynamic import, which routes around the AST import walk entirely
    ("importlib", "import_module"): "a dynamic import; the import walk cannot see it",
    ("ctypes", "CDLL"): "loads an arbitrary shared library",
    ("ctypes", "cdll"): "loads an arbitrary shared library",
}
"""Dotted names the import walk structurally cannot catch.

`np.random.default_rng()` needs only the allowed `import numpy`; `os` is stdlib and cannot
be banned outright, so its entropy and clock members are listed here.
"""

FORBIDDEN_CALLS = {
    "__import__": "a dynamic import; the import walk cannot see it",
    "eval": "arbitrary code, invisible to every check in this file",
    "exec": "arbitrary code, invisible to every check in this file",
}
"""Bare builtins with the same property, checked by call target rather than by attribute."""


def core_modules() -> list[Path]:
    return sorted(SRC.rglob("*.py"))


def dotted_imports(path: Path) -> set[str]:
    """Every full dotted path an import introduces, not just its first component.

    `from numpy import random` yields `numpy` and `numpy.random`, so a forbidden module is
    caught however it is spelled and whatever it is renamed to.
    """
    tree = ast.parse(path.read_text(), filename=str(path))
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                out.add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level or not node.module:
                continue  # relative import, stays inside the package
            out.add(node.module)
            for alias in node.names:
                out.add(f"{node.module}.{alias.name}")
    return out


def import_aliases(path: Path) -> dict[str, str]:
    """Local name -> module: `import numpy.random as r` gives `{"r": "numpy.random"}`."""
    tree = ast.parse(path.read_text(), filename=str(path))
    out: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                out[alias.asname or alias.name.split(".")[0]] = alias.name
    return out


def dotted_attributes(path: Path) -> set[tuple[str, str]]:
    """Every `<name>.<attr>` in one module, as pairs, with import aliases resolved.

    `import numpy as _np; _np.random.default_rng()` comes back as `("np", "random")`: on
    spelling alone, a two-line rename defeats the check. The local spelling is emitted too,
    because `FORBIDDEN_ATTRIBUTES` is keyed by the conventional name (`np`, `os`) and a
    module may reach a forbidden attribute through a name it never imported.
    """
    tree = ast.parse(path.read_text(), filename=str(path))
    aliases = import_aliases(path)
    canonical = {"numpy": "np"}
    out: set[tuple[str, str]] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
            name = node.value.id
            out.add((name, node.attr))
            module = aliases.get(name)
            if module:
                out.add((module, node.attr))
                out.add((canonical.get(module, module.split(".")[0]), node.attr))
    return out


def bare_calls(path: Path) -> set[str]:
    """Every `name(...)` call in one module, by callee name."""
    tree = ast.parse(path.read_text(), filename=str(path))
    return {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }


def top_level_imports(path: Path) -> set[str]:
    """Every top-level package name imported by one module."""
    return {n.split(".")[0] for n in dotted_imports(path)}


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
            f"{path.relative_to(SRC.parent.parent)} imports {name!r}; the core takes NumPy "
            f"and the stdlib only, analysis code belongs in tools/ or validation/."
        )


@pytest.mark.parametrize("path", core_modules(), ids=lambda p: p.name)
def test_core_does_not_import_a_forbidden_module_under_any_spelling(path: Path):
    """`from numpy import random` is `numpy.random`, whatever the first component says."""
    found = sorted(
        name
        for name in dotted_imports(path)
        for banned in FORBIDDEN_MODULES
        if name == banned or name.startswith(banned + ".")
    )
    assert not found, f"{path.name} imports {found}; see FORBIDDEN_MODULES."


@pytest.mark.parametrize("path", core_modules(), ids=lambda p: p.name)
def test_core_does_not_reach_a_clock_or_an_rng_through_an_allowed_module(path: Path):
    """The import ban is necessary and not sufficient -- see `FORBIDDEN_ATTRIBUTES`."""
    found = dotted_attributes(path) & set(FORBIDDEN_ATTRIBUTES)
    assert not found, "; ".join(
        f"{path.name} uses {a}.{b} -- {FORBIDDEN_ATTRIBUTES[(a, b)]}" for a, b in sorted(found)
    )


@pytest.mark.parametrize("path", core_modules(), ids=lambda p: p.name)
def test_core_does_not_import_or_execute_dynamically(path: Path):
    """`__import__`, `eval` and `exec` make every other check in this file advisory."""
    found = bare_calls(path) & set(FORBIDDEN_CALLS)
    assert not found, "; ".join(
        f"{path.name} calls {name}() -- {FORBIDDEN_CALLS[name]}" for name in sorted(found)
    )


# --------------------------------------------------------------------------------------
# A static analyser that silently matches nothing passes every test above. These run each
# detector against source that should trip it.
# --------------------------------------------------------------------------------------

DECOYS = [
    ("import scipy\n", "top_level_imports", "scipy"),
    ("from scipy import ndimage\n", "top_level_imports", "scipy"),
    ("def f():\n    import matplotlib\n", "top_level_imports", "matplotlib"),
    ("import numpy as np\nrng = np.random.default_rng()\n", "dotted", ("np", "random")),
    # aliased spellings
    ("import numpy as _np\nrng = _np.random.default_rng()\n", "dotted", ("np", "random")),
    ("import numpy\nrng = numpy.random.default_rng()\n", "dotted", ("np", "random")),
    ("import os as _os\nx = _os.urandom(8)\n", "dotted", ("os", "urandom")),
    ("from numpy import random\n", "dotted_import", "numpy.random"),
    ("import numpy.random as npr\n", "dotted_import", "numpy.random"),
    ("from numpy.random import default_rng\n", "dotted_import", "numpy.random"),
    ("m = __import__('scipy')\n", "call", "__import__"),
]


@pytest.mark.parametrize("source,kind,expected", DECOYS, ids=lambda v: str(v)[:40])
def test_each_detector_fires_on_source_that_should_trip_it(tmp_path, source, kind, expected):
    """A detector that matches nothing passes silently. Prove each one still matches."""
    decoy = tmp_path / "decoy.py"
    decoy.write_text(source)
    found = {
        "top_level_imports": top_level_imports,
        "dotted": dotted_attributes,
        "dotted_import": dotted_imports,
        "call": bare_calls,
    }[kind](decoy)
    assert expected in found, f"{kind} did not see {expected!r} in:\n{source}"


def test_every_banned_entry_is_owned_by_a_detector_that_can_see_it():
    """Every entry of the three banned sets has a shape its detector can match."""
    for a, b in FORBIDDEN_ATTRIBUTES:
        assert "." not in a and "." not in b, (
            f"FORBIDDEN_ATTRIBUTES key ({a!r}, {b!r}) is dotted; `dotted_attributes` "
            f"yields single-component pairs, so it can never match."
        )
    for name in FORBIDDEN_CALLS:
        assert "." not in name, (
            f"FORBIDDEN_CALLS entry {name!r} is dotted; `bare_calls` yields `ast.Name` "
            f"callees only, so it can never match. An attribute call belongs in "
            f"FORBIDDEN_ATTRIBUTES."
        )
    for name in FORBIDDEN_MODULES:
        assert name and not name.startswith("."), f"bad FORBIDDEN_MODULES entry {name!r}"


def _declares_or_evaluates_gas_properties(path: Path) -> bool:
    """`constants.py` declares the fit values and `thermo/` is the door; both are exempt."""
    rel = path.relative_to(SRC)
    # The top-level file only: `control/constants.py` has no business naming a gas fit.
    return rel.parts[0] == "thermo" or rel == Path("constants.py")


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
    """No gas-property fit constant is named outside `t700.thermo`, so the backend can be
    swapped without touching component code."""
    if _declares_or_evaluates_gas_properties(path):
        return

    source = path.read_text()
    used = [k for k in GAS_PROPERTY_CONSTANTS if k in source]
    assert not used, f"{path.name} references gas property constants {used}; use t700.thermo."


def _gas_property_values() -> dict[float, str]:
    """The numeric value of every gas-property fit constant, keyed by value."""
    from t700 import constants as c

    out: dict[float, str] = {}
    for name in GAS_PROPERTY_CONSTANTS:
        v = float(getattr(c, name))
        out[v] = name
        out[-v] = name  # a sign-flipped intercept is the same fit
    return out


@pytest.mark.parametrize("path", core_modules(), ids=lambda p: p.name)
def test_gas_property_constants_are_not_inlined_as_literals(path: Path):
    """The same fit written as bare numbers must not pass: `(h41 + 86.905) / 0.301` is
    `K_T41` by hand, invisible to the identifier check above.

    A second layer, not a proof: a rearranged fit (`h41 * 3.3223 + ...`) still escapes.
    """
    if _declares_or_evaluates_gas_properties(path):
        return

    values = _gas_property_values()
    tree = ast.parse(path.read_text(), filename=str(path))
    hits = sorted(
        {
            (node.value, values[node.value])
            for node in ast.walk(tree)
            if isinstance(node, ast.Constant)
            and isinstance(node.value, float)
            and node.value in values
        }
    )
    assert not hits, f"{path.name} inlines gas-property fit values {hits}; use t700.thermo."


def test_core_never_reads_the_reference_data():
    """`src/t700/` must not read `data/reference/`, the traces it is validated against
    (Ballin's curves and the GE models he checked against)."""
    for path in core_modules():
        source = path.read_text()
        for needle in ("data/reference", "reference/fig", '"reference"', "'reference'"):
            assert needle not in source, (
                f"{path.name} names {needle!r}; the core must not read validation data."
            )
