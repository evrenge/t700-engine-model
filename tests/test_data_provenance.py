"""Every file in data/ must say where its numbers came from.

CLAUDE.md's provenance rule is the one that cannot be fully mechanised -- a machine
cannot tell a cited constant from a plausible one. But it *can* check that the citation
exists, is well formed, and points at a real page of the report. That is worth having:
the failure mode this catches is a CSV that arrives with no header during a long session
and is indistinguishable from a good one six weeks later.

What it cannot check is whether the values are right. That is what the verify overlays
are for.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parent.parent / "data"
PDF_PAGES = 104  # docs/ballin-tm100991.pdf

REQUIRED = ("# source:", "# quantity:", "# method:")
"""Fields every data file must carry. Axis documentation is checked separately, by
column name, because a 2-D map documents three or more columns and cannot use a fixed
`# x:` / `# y:` pair."""


def data_files() -> list[Path]:
    return sorted(DATA.rglob("*.csv"))


def header_of(path: Path) -> list[str]:
    return [ln.rstrip("\n") for ln in path.read_text().splitlines() if ln.startswith("#")]


def test_there_are_data_files():
    """Guard against every test below passing vacuously."""
    assert data_files(), "no CSVs under data/ -- has digitizing produced anything?"


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
    """Each column name must appear in the header comments.

    This replaces a fixed `# x:` / `# y:` check, which quietly assumed every data file is
    a 2-column curve. `f1` is a 2-D map with five columns and documents each by name --
    that is better provenance, not worse, and the test now says so.
    """
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


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_values_are_finite_numbers(path: Path):
    with path.open() as fh:
        rows = list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))
    for i, row in enumerate(rows):
        for key, val in row.items():
            f = float(val)
            assert f == f and abs(f) != float("inf"), f"{path.name} row {i}: {key}={val}"


def test_f10_is_stored_as_printed_not_inverted():
    """Open question #5's *final* resolution, pinned.

    Figure A10's ordinate is labelled PS9/P49 but every plotted value exceeds 1, and
    Eq. 38 (P49 = Ps9 * f10) with the expansion order P45 > P49 > Ps9 requires f10 > 1.
    So **the axis label is the defect and the values are f10 directly** -- no inversion
    anywhere.

    An earlier resolution of #5 said to store as printed and take the reciprocal at use.
    That was wrong: it would put the exhaust total pressure below ambient, and no running
    engine can do that. If someone ever 'helpfully' inverts the file, every value drops
    below 1 and this fails.
    """
    path = DATA / "maps" / "f10_exhaust_pressure_loss.csv"
    if not path.exists():
        pytest.skip("f10 not digitized yet")
    with path.open() as fh:
        rows = list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))
    ys = [float(r["y"]) for r in rows]
    assert all(y > 1.0 for y in ys), (
        "f10 values are P49/Ps9 and must all exceed 1 -- the printed axis LABEL is "
        "inverted, the values are not. Do not invert this file. See open question #5."
    )
