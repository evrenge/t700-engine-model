"""Every markdown table in the repository's own documents must actually render.

A table needs a `|---|---|` separator row directly under its header. Without one, GitHub
and every other CommonMark renderer shows the raw pipes as a paragraph -- the content is
still *there*, so nothing fails, no link breaks, and the document reads fine in an editor.
It is invisible until someone opens the published page.

That is what happened on 2026-09-14. `SCOPE.md` was rewritten by a script that lifted two
tolerance tables out of the old file with

    [ln for ln in block.splitlines() if ln.startswith("| ") and "---" not in ln]

-- which drops the separator along with the rules it was meant to exclude -- and emitted
them without putting it back. Both tables shipped unrendered, and
`tests/test_scope_tolerances.py` did not notice because it parses rows, not layout.

Cell counts are checked too: a row with more cells than its header silently loses the
overflow in some renderers and shifts every column in others.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DOCS = sorted(
    [ROOT / "README.md", ROOT / "SCOPE.md", ROOT / "CLAUDE.md", ROOT / "site" / "README.md"]
    + list((ROOT / "docs" / "notes").glob("*.md"))
)

_SPLIT = re.compile(r"(?<!\\)\|")
_SEP = re.compile(r"^\s*\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)*\|?\s*$")


def _blocks(text: str):
    """Contiguous runs of lines starting with a pipe, with their 1-based line numbers."""
    run: list[tuple[int, str]] = []
    for n, line in enumerate(text.splitlines(), start=1):
        if line.lstrip().startswith("|"):
            run.append((n, line))
        elif run:
            yield run
            run = []
    if run:
        yield run


def _cells(line: str) -> int:
    return len(_SPLIT.split(line.strip().strip("|")))


@pytest.mark.parametrize("doc", DOCS, ids=lambda p: str(p.relative_to(ROOT)))
def test_every_table_has_a_separator_row(doc: Path):
    if not doc.exists():
        pytest.skip(f"{doc.name} not present")
    bad = []
    for run in _blocks(doc.read_text()):
        if len(run) < 2:
            continue  # a lone piped line is prose, not a table
        if not _SEP.match(run[1][1]):
            bad.append((run[0][0], run[0][1][:70]))
    assert not bad, (
        f"{doc.name}: table header not followed by a `|---|---|` separator, so it renders "
        f"as raw pipes: {bad}"
    )


@pytest.mark.parametrize("doc", DOCS, ids=lambda p: str(p.relative_to(ROOT)))
def test_every_table_row_matches_its_header(doc: Path):
    if not doc.exists():
        pytest.skip(f"{doc.name} not present")
    bad = []
    for run in _blocks(doc.read_text()):
        if len(run) < 2 or not _SEP.match(run[1][1]):
            continue
        want = _cells(run[0][1])
        for n, line in run[2:]:
            got = _cells(line)
            if got != want:
                bad.append((n, got, want))
    assert not bad, (
        f"{doc.name}: rows whose cell count differs from the header (line, got, expected): {bad}"
    )
