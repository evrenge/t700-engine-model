"""A measurement that README.md or SCOPE.md quotes from report.json must be today's value.

Measurements move whenever the model or a digitization does, and nothing re-reads the
sentences that quote them: README and SCOPE have both quoted numbers that had changed. So a
live number carries an invisible marker straight after it naming its source in
`site/assets/report.json`, which `validation/report.py` writes:

    rms 1.94 %<!--rj:eigenvalues.notes.table_1_stats.rms-->
    worst 0.81 %<!--rj:phase_plane.figure_9.worst|abs-->

- The path is dotted, each segment a dict key or an integer list index -- the lookup
  `tools/build_site.py` uses for the template's `{{rj:path|fmt}}`, in `tools/report_json.py`.
- The number is the last one before the marker on its line, with at most a unit and `**`
  between them. It must equal the value rounded to as many decimals as are printed. A
  printed sign (`+`, `-`, `−`) must match; an unsigned number is positive.
- `|abs` compares magnitudes, for a value quoted without its sign. A number printed with a
  sign never takes it. Inside a markdown table, write it `\\|abs`, as a table cell needs.

Missed: a number with no marker -- only what is marked is checked -- and the sentence around
a marked number, which can still say the wrong thing about the right value.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))

from report_json import load, lookup  # noqa: E402

DOCS = (ROOT / "README.md", ROOT / "SCOPE.md")

_MARKER = re.compile(r"<!--rj:([^|>\\]+?)(\\?\|abs)?-->")
_ANY_MARKER = re.compile(r"<!--\s*rj")
_NUMBER = re.compile(r"(?<![\w.])([+\-−]?)(\d[\d,]*(?:\.\d+)?)\s*[^\s\d*<>|]*\s*\**$")


def stale(text: str, data: dict) -> list[str]:
    """Every marker in `text` whose number is not the report.json value, as a message."""
    bad = []
    for m in _MARKER.finditer(text):
        line = text.count("\n", 0, m.start()) + 1
        path, absolute = m.group(1), bool(m.group(2))
        before = text[text.rfind("\n", 0, m.start()) + 1 : m.start()]
        n = _NUMBER.search(before)
        if n is None:
            bad.append(f"line {line}: no number just before {m.group(0)}")
            continue
        sign, digits = n.group(1), n.group(2).replace(",", "")
        try:
            value = lookup(data, path)
        except (KeyError, IndexError, TypeError, ValueError):
            bad.append(f"line {line}: {path} is not in report.json")
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            bad.append(f"line {line}: {path} holds {value!r}, not a number")
            continue
        if absolute:
            if sign:
                bad.append(f"line {line}: {n.group(0).strip()} is signed, so it cannot take abs")
                continue
            value = abs(value)
        decimals = len(digits.partition(".")[2])
        expected = format(value, f".{decimals}f")
        printed = -float(digits) if sign in ("-", "−") else float(digits)
        if printed != float(expected):
            bad.append(f"line {line}: prints {sign}{digits}, report.json {path} is {expected}")
    return bad


@pytest.mark.parametrize("doc", DOCS, ids=lambda p: p.name)
def test_every_marked_number_is_current(doc: Path):
    text = doc.read_text()
    assert _MARKER.search(text), f"{doc.name} has no live-number markers; the check is idle"
    assert len(_ANY_MARKER.findall(text)) == len(_MARKER.findall(text)), (
        f"{doc.name} has an `<!--rj` marker this test cannot parse, so it checks nothing"
    )
    bad = stale(text, load())
    assert not bad, (
        f"{doc.name} quotes numbers report.json no longer holds -- rerun validation/report.py "
        f"if the model moved, then correct the text:\n  " + "\n  ".join(bad)
    )


DATA = {"a": {"rms": 1.9372, "worst": -0.806}, "rows": [{"v": 7.64}, {"v": 16.38}]}


@pytest.mark.parametrize(
    "text",
    [
        "rms 1.94 %<!--rj:a.rms-->",
        "rms **1.9 %**<!--rj:a.rms-->",
        "worst −0.81 %<!--rj:a.worst-->",
        "worst -0.806<!--rj:a.worst-->",
        "worst 0.81 %<!--rj:a.worst|abs-->",
        "| worst 0.81 %<!--rj:a.worst\\|abs--> |",
        "from 7.6<!--rj:rows.0.v-->–16.4 rpm<!--rj:rows.1.v-->",
    ],
)
def test_a_current_number_passes(text: str):
    assert stale(text, DATA) == []


@pytest.mark.parametrize(
    "text",
    [
        "rms 1.95 %<!--rj:a.rms-->",  # the wrong number
        "rms 1.93 %<!--rj:a.rms-->",  # truncated, not rounded
        "worst 0.81 %<!--rj:a.worst-->",  # the sign dropped, with no abs
        "worst +0.81 %<!--rj:a.worst-->",  # the wrong sign
        "worst −0.81 %<!--rj:a.worst|abs-->",  # a signed number taking abs
        "rms 1.94 %<!--rj:a.mean-->",  # a path report.json does not have
        "rms 1.94 %<!--rj:rows.2.v-->",  # an index past the end
        "rms %<!--rj:a.rms-->",  # no number to check
    ],
)
def test_a_stale_number_fails(text: str):
    assert len(stale(text, DATA)) == 1, text
