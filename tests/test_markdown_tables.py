"""Every markdown table in the repository's own documents must actually render.

A table needs a `|---|---|` separator row directly under its header; without one, CommonMark
renderers show the raw pipes as a paragraph, and the document still reads fine in an editor.
A row with more cells than its header loses the overflow in some renderers and shifts every
column in others.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DOCS = sorted(
    [ROOT / "README.md", ROOT / "SCOPE.md", ROOT / "CLAUDE.md"]
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
def test_every_table_has_a_separator_and_rows_that_match_its_header(doc: Path):
    if not doc.exists():
        pytest.skip(f"{doc.name} not present")
    no_sep, bad = [], []
    for run in _blocks(doc.read_text()):
        if len(run) < 2:
            continue  # a lone piped line is prose, not a table
        if not _SEP.match(run[1][1]):
            no_sep.append((run[0][0], run[0][1][:70]))
            continue
        want = _cells(run[0][1])
        for n, line in run[2:]:
            got = _cells(line)
            if got != want:
                bad.append((n, got, want))
    assert not no_sep, (
        f"{doc.name}: table header not followed by a `|---|---|` separator, so it renders "
        f"as raw pipes: {no_sep}"
    )
    assert not bad, (
        f"{doc.name}: rows whose cell count differs from the header (line, got, expected): {bad}"
    )
