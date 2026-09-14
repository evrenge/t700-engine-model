"""The open-questions ledger must be structurally sound and must count itself correctly.

`docs/notes/open-questions.md` decides what this project believes is known. It has two
conventions and both were violated inside 48 hours of being restated in its own header:

* **Five columns, one row per question.** Row 58's 2026-09-14 closure was written into the
  *Question* column instead of the Status column, giving that row six cells. Markdown
  renders it without complaint.
* **The first word of the Status cell is the verdict.** Rows 57 and 58 carried "open" as
  their first word with "CLOSED" a thousand characters later. Appending a closure note
  without changing the leading word is the recurring failure, and it makes the ledger
  under-report its own progress -- which is how the header's tally came to be wrong by
  three closed, three open and one whole row.

The header's tally line was maintained by hand and drifted. It is checked here instead.
Found by the 2026-09-14 docs-and-ledger audit.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

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


def test_the_ledger_has_rows():
    assert len(rows()) > 50, "the ledger did not parse; the table format has changed"


def test_every_row_has_exactly_five_cells():
    """Six cells is not a formatting nit: it moves the verdict out of the Status column."""
    bad = [(n, len(c) - 2) for n, c in rows() if len(c) != 7]
    assert not bad, f"rows with the wrong cell count (number, cells): {bad}"


def test_every_status_cell_opens_with_a_verdict():
    bad = [(n, verdict_of(c[5])) for n, c in rows() if verdict_of(c[5]).startswith("UNRECOG")]
    assert not bad, (
        f"the first word of the Status cell is the verdict, and these do not have one: {bad}"
    )


def test_no_row_contradicts_its_own_verdict():
    """A cell whose verdict word is `open` must not announce a closure further in.

    This is the specific defect rows 57 and 58 carried. A closed row may of course quote
    the word "open" while narrating its own history, so only the open direction is checked.
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


@pytest.mark.parametrize("verdict", VERDICTS)
def test_the_header_tally_matches_the_rows(verdict: str):
    """The header states three counts. They are derived here, not trusted.

    Until 2026-09-14 the header said "56 logged -- 40 closed, 5 partly closed, 11 open"
    against an actual 56 / 37 / 5 / 14, and named row 57 as closed when its verdict word
    said open and row 58 as the next piece of work when it had shipped.
    """
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
    assert stated["logged"] == len(rows()), (
        f"header says {stated['logged']} logged, parsed {len(rows())}"
    )
    assert stated[verdict] == counts[verdict], (
        f"header says {stated[verdict]} {verdict}, parsed {counts[verdict]}"
    )


def test_every_closed_row_says_what_was_decided():
    """A verdict word is not an answer. `closed` alone tells a reader nothing.

    The records behind these rows were a separate 26,000-word file for half a day. It was
    folded back in on 2026-09-14: a finished project needs the decision, not the
    investigation that reached it, and the investigation is in git history. What this
    checks is that the fold left every closed row self-contained -- eleven of the
    forty-nine had a bare `**closed**` as their entire status, and eight more stated a
    verdict without a finding; all nineteen were written from the record before it went.

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
        f"closed rows stating a verdict but not a finding (row, words): {thin}. "
        f"Say what was decided; the record behind it is no longer in the tree."
    )


def test_nothing_still_points_at_the_deleted_record_file():
    assert "closed-questions.md" not in LEDGER.read_text(), (
        "the ledger points at closed-questions.md, which was folded back in and deleted"
    )
