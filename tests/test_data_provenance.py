"""Every file in data/ must say where its numbers came from.

A machine cannot tell a cited value from a plausible one, but it can check that the
citation exists, is well formed and points at a real page of the report. Whether the values
are right is for the verify overlays.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parent.parent / "data"
PDF_PAGES = 104  # docs/ballin-tm100991.pdf

REQUIRED = ("# source:", "# quantity:", "# method:")
"""Fields every data file must carry. Axes are checked by column name instead, since a 2-D
map has three or more columns."""


def data_files() -> list[Path]:
    return sorted(DATA.rglob("*.csv"))


def header_of(path: Path) -> list[str]:
    return [ln.rstrip("\n") for ln in path.read_text().splitlines() if ln.startswith("#")]


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_has_required_provenance_fields(path: Path):
    head = "\n".join(header_of(path))
    missing = [f for f in REQUIRED if f not in head]
    assert not missing, f"{path.name} is missing {missing}. See CLAUDE.md, the provenance rule."


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_cites_a_real_page_of_the_report(path: Path):
    """The source line must name TM-100991 and a page inside the document."""
    src = next((h for h in header_of(path) if h.startswith("# source:")), "")
    assert "TM-100991" in src, f"{path.name}: source does not name the report"
    m = re.search(r"pdf p\.(\d+)", src)
    assert m, f"{path.name}: source does not cite a pdf page -- got {src!r}"
    page = int(m.group(1))
    assert 1 <= page <= PDF_PAGES, f"{path.name}: cites pdf p.{page}, outside the document"


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_method_is_declared_and_honest(path: Path):
    """`digitized` is weaker provenance than `transcribed`, and must say so."""
    method = next((h for h in header_of(path) if h.startswith("# method:")), "")
    assert any(w in method for w in ("digitized", "transcribed", "derived")), method
    if "digitized" in method:
        head = "\n".join(header_of(path))
        assert "NOT transcribed" in head, (
            f"{path.name}: digitized data must carry the read-error disclaimer"
        )


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_every_column_is_documented(path: Path):
    """Each column name must appear in the header comments."""
    head = "\n".join(header_of(path))
    with path.open() as fh:
        reader = csv.reader(ln for ln in fh if not ln.startswith("#"))
        columns = next(reader)
    undocumented = [c for c in columns if c not in head]
    assert not undocumented, (
        f"{path.name}: columns {undocumented} are not described in the header. "
        f"Name every column and give its units."
    )


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_declared_point_count_matches_the_rows(path: Path):
    declared = next((h for h in header_of(path) if h.startswith("# points:")), None)
    if declared is None:
        pytest.skip(f"{path.name} declares no point count")
    n = int(re.search(r"(\d+)", declared).group(1))
    with path.open() as fh:
        rows = list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))
    assert len(rows) == n, f"{path.name}: header says {n} points, file holds {len(rows)}"


def _is_number(v: str) -> bool:
    try:
        float(v)
    except (TypeError, ValueError):
        return False
    return True


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_values_are_finite_numbers(path: Path):
    """Numeric columns hold finite numbers, and no column mixes numbers with text.

    Per column: `data/linear/` stores Appendix B's matrices long, `matrix,row,col,value`,
    with three label columns. A mixed column (a stray `N/A`, a truncated write, a header row
    duplicated into the body) is the failure guarded against.
    """
    with path.open() as fh:
        rows = list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))
    if not rows:
        pytest.skip(f"{path.name} holds no rows")

    for key in rows[0]:
        vals = [r[key] for r in rows]
        numeric = [_is_number(v) for v in vals]
        if any(numeric) and not all(numeric):
            bad = [(i, v) for i, (v, ok) in enumerate(zip(vals, numeric, strict=True)) if not ok]
            pytest.fail(
                f"{path.name}: column {key!r} mixes numbers and text; first offenders {bad[:3]}"
            )
        if not all(numeric):
            continue  # a label column, e.g. the state names in data/linear/
        for i, v in enumerate(vals):
            f = float(v)
            assert f == f and abs(f) != float("inf"), f"{path.name} row {i}: {key}={v}"
