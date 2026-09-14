"""Enforce the architecture constraints in CLAUDE.md.

Until this file existed, three of the rules in CLAUDE.md's enforcement table were marked
`advisory` -- which is to say, they were wishes. These tests make three of them real:

  * the model core imports NumPy and the standard library, nothing else
  * gas property evaluation happens only inside `t700.thermo`
  * the core does not read a clock or a random number generator

The provenance rule -- every number cites a report page -- stays advisory. It cannot be
mechanised, which is precisely why it is the one to be most careful about.

The 2026-09-13 dependency audit went looking for ways to satisfy these tests and break the
rules anyway, and found several. They are closed below and each is documented where it
sits, because the interesting artifact is the *shape* of the hole rather than the patch: an
AST walk over `import` statements answers "what does this module import", and three of the
four rules above are really about "what can this module reach".
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

NONDETERMINISTIC = {"random", "time", "datetime", "secrets", "uuid", "timeit", "hashlib"}
"""Same input must give bit-identical output, run to run and machine to machine.

`"os.urandom"` sat in this set until 2026-09-13 and could never have matched: the walk
below yields top-level package names only, so the entry a reader would have pointed at as
proof the RNG ban was covered was the one entry that did nothing. `os` cannot be banned
outright -- it is the standard library -- so its nondeterministic members are caught by
`FORBIDDEN_ATTRIBUTES` instead, which is where `os.urandom` actually lives now.

`hashlib` is here because `hash()`-like digests over floats are a plausible route to a
"deterministic" value that is not; `timeit` because `timeit.default_timer` is a clock that
does not import `time`.
"""

FORBIDDEN_MODULES = (
    "numpy.random",
    "os.urandom",
    "os.getrandom",
    "secrets",
    "random",
)
"""Dotted module paths the core may not import **by any spelling**.

`top_level_imports` reduces everything to its first component, which is the right
question for "is this package allowed" and the wrong one for `numpy.random`: numpy is the
one permitted third-party package, so `from numpy import random`,
`import numpy.random as npr` and `from numpy.random import default_rng` all reduced to
`"numpy"` and passed. The 2026-09-13 code-quality audit demonstrated all three against the
detectors added earlier the same day -- the RNG ban was the headline of that morning's fix
and three of its four spellings still escaped.

`FORBIDDEN_ATTRIBUTES` catches the fourth, `np.random.default_rng()`, at the use site.
This catches the rest at the import, which is strictly stronger: a module that cannot be
imported cannot be reached under an alias either.
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

The 2026-09-13 dependency audit enumerated the ways a module could satisfy
`test_core_imports_numpy_and_stdlib_only` and still break the rule it stands for. The
sharpest was `np.random.default_rng()`: NumPy is the one allowed third-party package, so
an RNG was reachable in the core with the import ban fully satisfied. No live violation
existed -- `src/` was clean on every vector the audit tried -- but a mechanism with a hole
that size is a wish wearing a mechanism's clothes.
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

    `import numpy.random as npr` yields `numpy.random`; `from numpy.random import
    default_rng` yields `numpy.random` and `numpy.random.default_rng`; `from numpy import
    random` yields `numpy` and `numpy.random`. So a forbidden module is caught however it
    is spelled and whatever it is renamed to.
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
    """Local name -> the module it actually names.

    `import numpy as _np` gives `{"_np": "numpy"}`; `import os` gives `{"os": "os"}`;
    `import numpy.random as r` gives `{"r": "numpy.random"}`. Used to resolve an attribute
    access back to the module it reaches, whatever the module was renamed to locally.
    """
    tree = ast.parse(path.read_text(), filename=str(path))
    out: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                out[alias.asname or alias.name.split(".")[0]] = alias.name
    return out


def dotted_attributes(path: Path) -> set[tuple[str, str]]:
    """Every `<name>.<attr>` in one module, as pairs -- with import aliases resolved.

    It matches on the spelling *after* resolving the local name through the module's own
    imports, so `np.random.default_rng()`, `numpy.random.default_rng()` and
    `import numpy as _np; _np.random.default_rng()` all come back as `("np", "random")`.

    **The alias resolution is the fix for a real hole.** Until 2026-09-14 this walker was
    "deliberately shallow: it matches on the *spelling*", and a two-line rename defeated it
    -- `import numpy as _np` then `_np.random.default_rng()` passed the whole import suite.
    CLAUDE.md claimed `FORBIDDEN_ATTRIBUTES` closed "the other dotted routes"; it closed the
    `np.`-spelled one. Determinism survived only because `tests/test_determinism.py` catches
    anything that reaches output -- an RNG that merely logged would have slipped both
    layers. Found by the 2026-09-14 docs-and-ledger audit, which broke the rule to check it.

    Both spellings are emitted, not just the canonical one, because `FORBIDDEN_ATTRIBUTES`
    is keyed by the conventional local name (`np`, `os`) and a module may also reach a
    forbidden attribute through a name it never imported.
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


@pytest.mark.parametrize("path", core_modules(), ids=lambda p: p.name)
def test_core_does_not_import_a_forbidden_module_under_any_spelling(path: Path):
    """`from numpy import random` is `numpy.random`, whatever the first component says."""
    found = sorted(
        name
        for name in dotted_imports(path)
        for banned in FORBIDDEN_MODULES
        if name == banned or name.startswith(banned + ".")
    )
    assert not found, (
        f"{path.name} imports {found}. The first component of those is an allowed "
        f"package, which is why the import walk passed them -- see FORBIDDEN_MODULES."
    )


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
# The checks above are static analysers, and a static analyser that silently matches
# nothing passes every test in this file. `"os.urandom"` sat in NONDETERMINISTIC for three
# days doing exactly that. These run each detector against source that *should* trip it.
# --------------------------------------------------------------------------------------

DECOYS = [
    ("import scipy\n", "top_level_imports", "scipy"),
    ("from scipy import ndimage\n", "top_level_imports", "scipy"),
    ("def f():\n    import matplotlib\n", "top_level_imports", "matplotlib"),
    ("import numpy as np\nrng = np.random.default_rng()\n", "dotted", ("np", "random")),
    # the alias routes that defeated this walker until 2026-09-14
    ("import numpy as _np\nrng = _np.random.default_rng()\n", "dotted", ("np", "random")),
    ("import numpy\nrng = numpy.random.default_rng()\n", "dotted", ("np", "random")),
    ("import os as _os\nx = _os.urandom(8)\n", "dotted", ("os", "urandom")),
    ("import time as _t\nx = _t.monotonic()\n", "dotted", ("time", "monotonic")),
    ("from numpy import random\n", "dotted_import", "numpy.random"),
    ("import numpy.random as npr\n", "dotted_import", "numpy.random"),
    ("from numpy.random import default_rng\n", "dotted_import", "numpy.random"),
    ("import numpy.random\n", "dotted_import", "numpy.random"),
    ("import os\nx = os.urandom(8)\n", "dotted", ("os", "urandom")),
    ("import os\nt = os.times()\n", "dotted", ("os", "times")),
    (
        "import importlib\nm = importlib.import_module('scipy')\n",
        "dotted",
        ("importlib", "import_module"),
    ),
    ("import ctypes\nlib = ctypes.CDLL('libm.so.6')\n", "dotted", ("ctypes", "CDLL")),
    ("m = __import__('scipy')\n", "call", "__import__"),
    ("eval('1+1')\n", "call", "eval"),
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


@pytest.mark.parametrize("name", sorted(NONDETERMINISTIC))
def test_no_nondeterministic_entry_is_spelled_so_it_can_never_match(name: str):
    """NONDETERMINISTIC is compared against top-level names, so a dotted entry is dead.

    This is the check that would have caught `"os.urandom"` sitting in that set for three
    days doing nothing.
    """
    assert "." not in name, (
        f"NONDETERMINISTIC entry {name!r} is dotted, and the import walk yields top-level "
        f"names only -- it can never match. Put it in FORBIDDEN_ATTRIBUTES."
    )


def test_every_banned_entry_is_owned_by_a_detector_that_can_see_it():
    """The sibling check for the other three sets, and it was missing.

    `test_no_banned_name_is_unreachable_by_construction` parametrized over all four sets
    and asserted nothing for three of them -- 8 of its 15 cases ran an empty body,
    including the `FORBIDDEN_ATTRIBUTES` half its own docstring claimed to check. Found by
    the 2026-09-13 code-quality audit.
    """
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
    # `rel.name == "constants.py"` was the test until 2026-09-13, which also exempted
    # `control/constants.py` -- a file in a different appendix that has no business
    # naming a gas-property fit. The exemption is for the one file that declares the
    # values, so it names that one file.
    if rel.parts[0] == "thermo" or rel == Path("constants.py"):
        return

    source = path.read_text()
    used = [k for k in GAS_PROPERTY_CONSTANTS if k in source]
    assert not used, (
        f"{path.name} references gas property constants {used}. Go through "
        f"t700.thermo instead -- see CLAUDE.md, 'All gas property evaluation goes "
        f"through t700.thermo'."
    )


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
    """The same fit, written as a bare number, must not pass.

    `test_gas_properties_only_inside_thermo` matches **identifiers**, so it enforces a
    naming convention rather than the rule the convention stands for. The 2026-09-14
    docs-and-ledger audit broke the rule without tripping it: `K_T41_1 * h41 + K_T41_2`
    rewritten as `(h41 + 86.905) / 0.301` sat in `engine.py` with the whole suite green,
    defeating the thermo rule and the provenance rule in one line. That is the "detector
    fooled by its own subject matter" pattern this project has hit before.

    This closes the literal route: every float in a core module is compared against the
    *values* of the fit constants. It cannot catch an obfuscated arithmetic rearrangement
    -- `(h41 + 86.905) / 0.301` is caught because 86.905 and 0.301 are both literally in
    `constants.py`, but `h41 * 3.3223 + ...` written to four places would not be -- so it
    is a second layer, not a proof. Both layers together are still weaker than the reading
    CLAUDE.md asks for, which is why the gas-property rule stays worth challenging on sight.
    """
    rel = path.relative_to(SRC)
    if rel.parts[0] == "thermo" or rel == Path("constants.py"):
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
    assert not hits, (
        f"{path.name} contains the numeric value of a gas-property fit constant as a bare "
        f"literal: {hits}. That is the fit welded in by hand, and it defeats both the "
        f"thermo rule and the provenance rule. Go through t700.thermo."
    )


def test_core_never_reads_the_reference_data():
    """`src/t700/` must not consume anything under `data/reference/`.

    That directory holds the traces we validate *against* -- Ballin's published curves and
    the two GE models he himself was checking. A model that read its own validation data
    would not be a model. The GE series in particular is a different engine model, and
    agreeing with it is not evidence of anything; see `validation/test_ge_reference.py`.

    This is the mechanism behind the 2026-09-12 audit recorded in that file.
    """
    for path in core_modules():
        source = path.read_text()
        for needle in ("data/reference", "reference/fig", '"reference"', "'reference'"):
            assert needle not in source, (
                f"{path.name} names {needle!r}. The model core must not read validation "
                f"data; comparisons live in validation/."
            )
