"""One value out of `results/report.json`, addressed by a dotted path.

Used by `tests/test_live_numbers.py`, which holds the `<!--rj:path-->` markers in `README.md`
and `SCOPE.md` to that file.
"""

from __future__ import annotations

import json
from pathlib import Path

REPORT_JSON = Path(__file__).resolve().parent.parent / "results" / "report.json"


def load() -> dict:
    return json.loads(REPORT_JSON.read_text())


def lookup(data, path: str):
    """The value at `path`, e.g. `closed_loop.rows.1.cycle_ptp`: each dot-separated segment is
    a dict key, or an integer index into a list. A key that contains a dot cannot be named."""
    value = data
    for segment in path.split("."):
        value = value[int(segment)] if isinstance(value, list) else value[segment]
    return value
