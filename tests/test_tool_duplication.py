"""No function body may appear twice in `tools/`, `validation/`, `tests/` or `src/`, whether
in two files or in one.

Copied digitizers let a fix made for one figure miss the others (`docs/notes/digitizing-recipe.md`
traps 7 and 8). So: one implementation per concept, in `tools/figure/`, and a second copy is a
test failure.

Compared is the AST of each function body with the docstring and `-> type` annotations
stripped and local names normalised, so a copy that renamed a loop variable or gained an
annotation still matches. Missed: a copy that has already been edited. This catches a copy
before it diverges, the only moment it is cheap.
"""

from __future__ import annotations

import ast
import hashlib
from collections import defaultdict
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FOLDERS = ("tools", "validation", "tests", "src")

MIN_STATEMENTS = 3
"""Below this, two functions coinciding is a coincidence rather than a copy. Counted with
the `def` itself, so 3 is a two-line body."""

ALLOWED: dict[str, str] = {}
"""Deliberate duplicates, keyed `path:function` so one entry excuses one copy, each with the
reason it cannot be shared: sharing would couple two things that should not be."""


def _locals_of(fn: ast.FunctionDef) -> set[str]:
    """Argument and assigned names -- everything whose spelling is the author's choice."""
    names = {a.arg for a in fn.args.args + fn.args.kwonlyargs + fn.args.posonlyargs}
    if fn.args.vararg:
        names.add(fn.args.vararg.arg)
    if fn.args.kwarg:
        names.add(fn.args.kwarg.arg)
    for node in ast.walk(fn):
        if isinstance(node, ast.Name) and isinstance(node.ctx, (ast.Store, ast.Del)):
            names.add(node.id)
        elif isinstance(node, ast.ExceptHandler) and node.name:
            names.add(node.name)
    return names


class _Normalise(ast.NodeTransformer):
    def __init__(self, local_names: set[str]) -> None:
        self.locals = local_names
        self.seen: dict[str, str] = {}

    def _rename(self, name: str) -> str:
        if name not in self.locals:
            return name
        if name not in self.seen:
            self.seen[name] = f"_v{len(self.seen)}"
        return self.seen[name]

    def visit_Name(self, node: ast.Name) -> ast.Name:
        node.id = self._rename(node.id)
        return node

    def visit_arg(self, node: ast.arg) -> ast.arg:
        node.arg = self._rename(node.arg)
        node.annotation = None
        return node


def _fingerprint(fn: ast.FunctionDef) -> tuple[int, str]:
    fn = ast.parse(ast.unparse(fn)).body[0]
    fn.returns = None
    body = fn.body
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        if isinstance(body[0].value.value, str):
            fn.body = body[1:] or [ast.Pass()]
    fn = _Normalise(_locals_of(fn)).visit(fn)
    fn.name = "_f"
    n = sum(1 for node in ast.walk(fn) if isinstance(node, ast.stmt))
    return n, hashlib.sha256(ast.dump(fn, include_attributes=False).encode()).hexdigest()[:16]


def _copies(files: dict[str, str], allowed: dict[str, str]) -> list[list[str]]:
    """Every group of two or more functions sharing one body, as `path:function` entries."""
    seen: dict[str, list[str]] = defaultdict(list)
    for rel, text in files.items():
        for node in ast.walk(ast.parse(text)):
            if isinstance(node, ast.FunctionDef):
                where = f"{rel}:{node.name}"
                n, key = _fingerprint(node)
                if n >= MIN_STATEMENTS and where not in allowed:
                    seen[key].append(where)
    return [places for places in seen.values() if len(places) > 1]


def test_no_function_body_appears_twice_in_the_project():
    files = {
        str(path.relative_to(ROOT)): path.read_text()
        for folder in FOLDERS
        for path in sorted((ROOT / folder).rglob("*.py"))
    }
    unread = [f for f in FOLDERS if not any(rel.startswith(f + "/") for rel in files)]
    assert not unread, f"no Python file found under {unread}: the scan no longer reaches them"
    dups = _copies(files, ALLOWED)
    assert not dups, (
        "the same function body appears more than once:\n  "
        + "\n  ".join("  ==  ".join(places) for places in dups)
        + "\nKeep one and import it; under tools/, it belongs in tools/figure/."
    )


FINGERPRINT_DECOYS = [
    # the same body with different local names must still collide
    (
        "def f(pos, vals):\n"
        "    out = []\n"
        "    for k in pos:\n"
        "        out.append(vals[k])\n"
        "    return out\n",
        "def g(positions, values):\n"
        "    acc = []\n"
        "    for j in positions:\n"
        "        acc.append(values[j])\n"
        "    return acc\n",
        True,
    ),
    # the other way a checker fails: flagging everything
    (
        "def f(x):\n    y = x + 1\n    z = y * 2\n    return z\n",
        "def g(x):\n    y = x - 1\n    z = y / 2\n    return z\n",
        False,
    ),
    # a `-> type` annotation on one of two copies must not hide them
    (
        "def f(x) -> float:\n    y = x + 1\n    return y\n",
        "def g(x):\n    y = x + 1\n    return y\n",
        True,
    ),
]


@pytest.mark.parametrize(
    "src_a, src_b, same", FINGERPRINT_DECOYS, ids=["rename", "logic", "annotation"]
)
def test_the_fingerprint_sees_through_spelling_and_not_through_logic(src_a, src_b, same):
    """Decoys, because a checker that cannot fail proves nothing."""
    na, ka = _fingerprint(ast.parse(src_a).body[0])
    _, kb = _fingerprint(ast.parse(src_b).body[0])
    assert (ka == kb) == same, (
        f"the fingerprint {'no longer sees' if same else 'now invents'} a copy in:\n"
        f"{src_a}\n{src_b}"
    )
    assert na >= MIN_STATEMENTS or not same, f"{na} statements is below the copy threshold"


def test_a_two_line_copy_inside_one_file_is_found():
    """Decoy: the shortest body that counts, copied within a single file."""
    body = "    q = uv_of(nm, pos)\n    return q[1]\n"
    src = f"def at(nm, pos):\n{body}\n\ndef tick_u(nm, pos):\n{body}"
    assert _copies({"tools/a.py": src}, {}) == [["tools/a.py:at", "tools/a.py:tick_u"]], (
        "a two-line copy, or a copy within one file, is no longer found"
    )


def test_an_exception_excuses_one_copy_and_not_every_function_of_that_name():
    """Decoy: an `ALLOWED` entry names a file, so the other copies still collide."""
    src = "def f(x):\n    y = x + 1\n    return y\n"
    files = {"tools/a.py": src, "tools/b.py": src, "validation/c.py": src}
    assert _copies(files, {"tools/a.py:f": "decoy"}) == [["tools/b.py:f", "validation/c.py:f"]]
