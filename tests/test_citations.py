"""Every test or code name that a document or a comment cites must exist.

A sentence naming a test goes on naming it after the test is renamed, merged or deleted, and
nothing else notices: docstrings here have cited tests that did not exist. This test reads

- the markdown (`README.md`, `SCOPE.md`, `CLAUDE.md`, `site/README.md`, `docs/notes/*.md`)
  and `site/template.html`, whole;
- the comments (by `tokenize`) and the string statements -- docstrings, including the ones
  under constants (by `ast`) -- of every `.py` file in `src/`, `tests/`, `validation/` and
  `tools/`,

and resolves three citation forms against an AST index of those four folders. Nothing is
imported: importing a validation module runs SciPy.

    <file>.py::<name>, <dir>/<file>.py::<name>, test_<module>::<name>
                          <name> is defined in that file (any file whose path ends so)
    test_<module>.<name>  the same, for a module whose name starts with `test_`
    `test_<name>`         in backticks, bare: a test function of that name exists somewhere,
                          or a test module does (the docs name `test_trim_points` the file)

"Defined" means a module-level function, class or assignment, or `<Class>.<method>` (also
written `<Class>::<method>`). A placeholder that is deliberately not a name goes in `ALLOWED`
with its reason. Missed: any other citation form (`maps.clamp_scope()`, a test named in plain
prose without backticks), and a citation that resolves to the wrong existing thing.
"""

from __future__ import annotations

import ast
import io
import re
import tokenize
from collections.abc import Iterator
from functools import cache
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
CODE = ("src", "tests", "validation", "tools")
DOCS = sorted(
    [ROOT / "README.md", ROOT / "SCOPE.md", ROOT / "CLAUDE.md", ROOT / "site" / "README.md"]
    + list((ROOT / "docs" / "notes").glob("*.md"))
    + [ROOT / "site" / "template.html"]
)

ALLOWED: dict[str, str] = {
    "test_x.SYMBOL": "CLAUDE.md's placeholder for the form test_scope_tolerances.py resolves",
}
"""Exact citation strings that are not names, each with the reason it is not."""

_FILE_REF = re.compile(
    r"(?<![\w/.-])((?:[\w.-]+/)*(?:\w+\.py|test_\w+))::(\w+(?:(?:::|\.)\w+)?)"
)
_MODULE_REF = re.compile(r"(?<![\w.])(test_\w+)\.([A-Za-z_]\w*)")
_BARE_TEST = re.compile(r"`(test_\w+)(?:\[[^`\]]*\])?`")


def _bound(node: ast.AST) -> Iterator[str]:
    """Names a module-level statement binds: definitions and assignment targets, including
    those inside `if`/`try`/`with` blocks at module level."""
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
        yield node.name
        if isinstance(node, ast.ClassDef):
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    yield f"{node.name}.{item.name}"
        return
    targets: list[ast.expr] = []
    if isinstance(node, ast.Assign):
        targets = node.targets
    elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
        targets = [node.target]
    for t in targets:
        for n in ast.walk(t):
            if isinstance(n, ast.Name):
                yield n.id
    for field in ("body", "orelse", "finalbody", "handlers"):
        for child in getattr(node, field, ()):
            yield from _bound(child)


def _python_files() -> list[Path]:
    return sorted(
        p
        for folder in CODE
        for p in (ROOT / folder).rglob("*.py")
        if "__pycache__" not in p.parts
    )


@cache
def index() -> tuple[dict[str, frozenset[str]], frozenset[str]]:
    """Names defined per file (keyed by path relative to the root), and every test function
    and test module name in the tree."""
    names: dict[str, frozenset[str]] = {}
    tests: set[str] = set()
    for path in _python_files():
        tree = ast.parse(path.read_text(), filename=str(path))
        names[path.relative_to(ROOT).as_posix()] = frozenset(
            n for stmt in tree.body for n in _bound(stmt)
        )
        if path.stem.startswith("test_"):
            tests.add(path.stem)
        tests |= {
            n.name
            for n in ast.walk(tree)
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
            and n.name.startswith("test_")
        }
    return names, frozenset(tests)


def _defines(file: str, name: str) -> bool:
    """Whether any indexed file whose path ends in `file` defines `name`."""
    if not file.endswith(".py"):
        file += ".py"
    names, _ = index()
    hits = [v for k, v in names.items() if k == file or k.endswith("/" + file)]
    return any(name.replace("::", ".") in v for v in hits)


def broken_citations(text: str) -> list[tuple[int, str]]:
    """Every citation in `text` that does not resolve, as (line within `text`, citation)."""
    _, tests = index()
    bad: list[tuple[int, str]] = []

    def report(m: re.Match, ok: bool) -> None:
        if not ok and m.group(0).strip("`") not in ALLOWED:
            bad.append((text.count("\n", 0, m.start()) + 1, m.group(0)))

    for m in _FILE_REF.finditer(text):
        report(m, _defines(m.group(1), m.group(2)))
    for m in _MODULE_REF.finditer(text):
        if m.group(2) != "py" and not text.startswith("::", m.end()):
            report(m, _defines(m.group(1), m.group(2)))
    for m in _BARE_TEST.finditer(text):
        report(m, m.group(1) in tests)
    return bad


def _python_prose(path: Path) -> Iterator[tuple[int, str]]:
    """Each comment and each string statement in a module, with its first line."""
    source = path.read_text()
    for tok in tokenize.generate_tokens(io.StringIO(source).readline):
        if tok.type == tokenize.COMMENT:
            yield tok.start[0], tok.string
    for node in ast.walk(ast.parse(source, filename=str(path))):
        if (
            isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
        ):
            yield node.lineno, node.value.value


def _sources() -> Iterator[tuple[Path, int, str]]:
    for doc in DOCS:
        if doc.exists():
            yield doc, 1, doc.read_text()
    for path in _python_files():
        for line, text in _python_prose(path):
            yield path, line, text


def test_every_cited_name_exists():
    broken = [
        f"{path.relative_to(ROOT)}:{first + line - 1}: {cite}"
        for path, first, text in _sources()
        for line, cite in broken_citations(text)
    ]
    assert not broken, "citations of names that do not exist:\n  " + "\n  ".join(broken)


def test_the_index_sees_what_the_citations_point_at():
    """A broken index would pass everything: the checker must resolve each form it reads."""
    good = (
        "tests/test_maps.py::test_the_clamp_counter_cannot_influence_any_model_output, "
        "`validation/report.py::internal_spread`, test_trim_points.EIGENVALUE_TOL_PCT, "
        "test_fuel_step._run, `test_every_cited_name_exists`, "
        "test_trim_points::test_jacobian_eigenvalues_against_table_1"
    )
    assert broken_citations(good) == []


@pytest.mark.parametrize(
    "decoy",
    [
        "see `test_a_test_nobody_ever_wrote`",
        "pinned by test_trim_points.NO_SUCH_TOL_PCT",
        "removed with test_ge_reference._read",
        "tests/test_maps.py::test_a_clamp_test_that_is_not_there",
        "report.py::no_such_function",
        "test_linear.py::NoClass::method",
    ],
)
def test_a_made_up_name_fails(decoy: str):
    assert broken_citations(decoy) == [(1, decoy.split()[-1])], decoy


def test_the_placeholders_are_still_used():
    """An allow-list entry whose text no longer appears anywhere is dead weight that would
    excuse the next real citation to read the same."""
    corpus = "\n".join(text for _, _, text in _sources())
    unused = [k for k in ALLOWED if k not in corpus]
    assert not unused, f"ALLOWED entries no source contains any more: {unused}"
