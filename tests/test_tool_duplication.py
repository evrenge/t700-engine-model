"""No two files under `tools/` may contain the same function twice.

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

The AST of each function body, with the docstring stripped and local names normalised, so
a copy that drifted in its prose or renamed a loop variable still matches. `v` against `vv`
in three copies of `paired` is precisely the difference that hid the third copy from a
first, name-sensitive version of this check.

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

TOOLS = Path(__file__).resolve().parent.parent / "tools"

MIN_STATEMENTS = 4
"""Below this, two functions coinciding is a coincidence rather than a copy."""

ALLOWED: dict[str, str] = {}
"""Deliberate duplicates, each with the reason it cannot be shared. Empty, and the bar for
adding an entry is that sharing the code would couple two things that should not be."""


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
    body = fn.body
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        if isinstance(body[0].value.value, str):
            fn.body = body[1:] or [ast.Pass()]
    fn = _Normalise(_locals_of(fn)).visit(fn)
    fn.name = "_f"
    n = sum(1 for node in ast.walk(fn) if isinstance(node, ast.stmt))
    return n, hashlib.sha256(ast.dump(fn, include_attributes=False).encode()).hexdigest()[:16]


def test_no_function_body_appears_in_two_files_under_tools():
    seen: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for path in sorted(TOOLS.rglob("*.py")):
        rel = str(path.relative_to(TOOLS.parent))
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.FunctionDef):
                n, key = _fingerprint(node)
                if n >= MIN_STATEMENTS and node.name not in ALLOWED:
                    seen[key].append((rel, node.name))
    dups = {key: places for key, places in seen.items() if len({f for f, _ in places}) > 1}
    assert not dups, (
        "the same function body appears in more than one file under tools/:\n  "
        + "\n  ".join("  ==  ".join(f"{f}:{n}" for f, n in places) for places in dups.values())
        + "\nMove it into tools/figure/ and import it. A second copy is how an improvement "
        "made for one figure fails to reach the others."
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
