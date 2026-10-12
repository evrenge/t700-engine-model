"""No function body may appear twice in `tools/`, `validation/`, `tests/` or `src/`, whether
in two files or in one.

## Why this is a test and not a note

Every digitizer in this repository grew out of the one before it, by copying. That is how
five separate implementations of one pipeline came to exist, and the cost was not the
duplication itself -- it was that **an improvement made for one figure did not reach the
others**, and nothing said so. Three examples, all found on 2026-09-16 and all recorded in
`docs/notes/digitizing-recipe.md` as traps the project had already learned:

* trap 7, the frame edges are arcs and a homography is wrong: fixed for A2, A7 and A9,
  never propagated to the thirteen figures of the single-curve path;
* trap 8, Nelder-Mead's default simplex walks a marker off the page: fixed in one copy of
  `fit_glyphs`, left optional and defaulted off in the other, so `f2` and `f9` -- the two
  highest-leverage maps in the engine -- were fitted without it;
* the anchoring rule of the axis map, which `digitize_a7` argued out in its docstring and
  `digitize_fig678` rediscovered independently, while the thirteen-figure path kept the
  weaker form.

None of those was a bug anyone could see in one file. They were only visible by reading two
files side by side and noticing they disagreed, which is exactly what nobody does. So the
rule is mechanised: one implementation per concept, in `tools/figure/`, and a second copy
of anything is a test failure.

## What it actually compares

The AST of each function body, with the docstring and the `-> type` annotations stripped
and local names normalised, so a copy that drifted in its prose, gained an annotation or
renamed a loop variable still matches. `v` against `vv` in three copies of `paired` is
precisely the difference that hid the third copy from a first, name-sensitive version of
this check; a return annotation on one of two `_rungs` hid them from the second.

## What it misses

**A copy that was edited.** Two functions that do the same thing by different code are
invisible here, and that is the dangerous case, because it is what "improved in one place"
produces: the copies stop being identical the moment one of them is fixed. This test
catches the copy *before* it diverges, which is the only moment it is cheap. After that it
takes a reader.

It also says nothing about whether the one surviving implementation is right.
"""

from __future__ import annotations

import ast
import hashlib
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FOLDERS = ("tools", "validation", "tests", "src")

MIN_STATEMENTS = 3
"""Below this, two functions coinciding is a coincidence rather than a copy. Counted with
the `def` itself, so 3 is a two-line body: at 4, a two-line helper had been copied into
three tools unseen."""

ALLOWED: dict[str, str] = {}
"""Deliberate duplicates, keyed `path:function` so that one entry excuses one copy and not
every function of that name, each with the reason it cannot be shared. Empty, and the bar
for adding an entry is that sharing the code would couple two things that should not be."""


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
        + "\nKeep one and import it; under tools/, it belongs in tools/figure/. A second copy "
        "is how an improvement made in one place fails to reach the others."
    )


def test_the_detector_finds_a_copy_it_is_shown():
    """A decoy, because a checker that cannot fail proves nothing.

    The same body, spelled with different local names in two notional files, must still
    collide -- that is the case a name-sensitive fingerprint missed.
    """
    a = ast.parse(
        "def f(pos, vals):\n"
        "    out = []\n"
        "    for k in pos:\n"
        "        out.append(vals[k])\n"
        "    return out\n"
    ).body[0]
    b = ast.parse(
        "def g(positions, values):\n"
        "    acc = []\n"
        "    for j in positions:\n"
        "        acc.append(values[j])\n"
        "    return acc\n"
    ).body[0]
    na, ka = _fingerprint(a)
    nb, kb = _fingerprint(b)
    assert na >= MIN_STATEMENTS and ka == kb, (
        "the duplicate detector no longer sees through a rename, so a copied function "
        "with tidied variable names would pass unnoticed"
    )


def test_the_detector_does_not_call_two_different_functions_a_copy():
    """The other way a checker fails: flagging everything."""
    a = ast.parse("def f(x):\n    y = x + 1\n    z = y * 2\n    return z\n").body[0]
    b = ast.parse("def g(x):\n    y = x - 1\n    z = y / 2\n    return z\n").body[0]
    assert _fingerprint(a)[1] != _fingerprint(b)[1]


def test_a_return_annotation_does_not_hide_a_copy():
    """Decoy for the blind spot that hid the two `_rungs` in `validation/`."""
    a = ast.parse("def f(x) -> float:\n    y = x + 1\n    return y\n").body[0]
    b = ast.parse("def g(x):\n    y = x + 1\n    return y\n").body[0]
    assert _fingerprint(a)[1] == _fingerprint(b)[1], (
        "a `-> type` annotation on one of two copies makes them look different again"
    )


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
