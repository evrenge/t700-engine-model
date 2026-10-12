"""The open-questions ledger must be structurally sound and must count itself correctly.

`docs/notes/open-questions.md` has two conventions: five columns, one row per question; and
the first word of the Status cell is the verdict. Appending a closure note without changing
that word is the recurring failure, and makes the ledger under-report its own progress. The
header's tally is derived here rather than trusted.
"""

from __future__ import annotations

import re
from pathlib import Path

LEDGER = Path(__file__).resolve().parent.parent / "docs" / "notes" / "open-questions.md"

VERDICTS = ("closed", "partly closed", "open")

# A cell may hold an escaped pipe; markdown splits on the unescaped ones.
_SPLIT = re.compile(r"(?<!\\)\|")


def rows() -> list[tuple[int, list[str]]]:
    out = []
    for line in LEDGER.read_text().splitlines():
        if not line.startswith("| "):
            continue
        cells = _SPLIT.split(line)
        if len(cells) < 3 or not cells[1].strip().isdigit():
            continue
        out.append((int(cells[1].strip()), cells))
    return out


def verdict_of(status: str) -> str:
    s = status.strip().lstrip("*").strip().lower()
    for v in ("partly closed", "closed", "open"):
        if s.startswith(v):
            return v
    return f"UNRECOGNISED: {s[:40]!r}"


def test_no_row_contradicts_its_own_verdict():
    """A cell whose verdict word is `open` must not announce a closure further in.

    A closed row may quote the word "open" while narrating its history, so only the open
    direction is checked.
    """
    offenders = []
    for n, c in rows():
        status = c[5]
        if verdict_of(status) != "open":
            continue
        body = status.strip()
        # skip the leading verdict clause itself
        tail = body[body.find(".") + 1 :] if "." in body else body
        if re.search(r"\bCLOSED\b", tail):
            offenders.append(n)
    assert not offenders, (
        f"rows {offenders} lead with `open` but announce CLOSED further into the cell. "
        f"Move the verdict to the front; the header says the first word is the verdict."
    )


def test_row_numbers_are_unique():
    ns = [n for n, _ in rows()]
    dupes = {n for n in ns if ns.count(n) > 1}
    assert not dupes, f"duplicate row numbers {sorted(dupes)} -- never renumber, but never collide"


def test_the_header_tally_matches_the_rows():
    """The header's counts are derived here, not trusted. That needs every Status cell to
    open with its verdict, so that is asserted first."""
    bad = [(n, verdict_of(c[5])) for n, c in rows() if verdict_of(c[5]).startswith("UNRECOG")]
    assert not bad, (
        f"the first word of the Status cell is the verdict, and these do not have one: {bad}"
    )
    counts = {v: 0 for v in VERDICTS}
    for _, c in rows():
        counts[verdict_of(c[5])] += 1
    header = LEDGER.read_text().split("| # |")[0]
    m = re.search(r"(\d+) logged — (\d+) closed,\s*\n?>?\s*(\d+) partly closed, (\d+) open", header)
    assert m, "the header tally line is missing or has changed shape"
    stated = dict(
        zip(
            ("logged", "closed", "partly closed", "open"),
            (int(g) for g in m.groups()),
            strict=True,
        )
    )
    parsed = {"logged": len(rows()), **counts}
    assert stated == parsed, f"header says {stated}, parsed {parsed}"


def test_every_closed_row_says_what_was_decided():
    """A verdict word is not an answer: a closed row states what was decided.

    Eight words is deliberately low. Row 28's whole answer is "moot: `R` is never
    applied" and it is complete; the bar is a finding, not a length.
    """
    thin = []
    for n, c in rows():
        if verdict_of(c[5]) != "closed":
            continue
        words = len(re.sub(r"[*`]", "", c[5]).split())
        if words < 8:
            thin.append((n, words))
    assert not thin, (
        f"closed rows stating a verdict but not a finding (row, words): {thin}; say what "
        f"was decided."
    )
