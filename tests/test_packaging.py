"""The distribution must contain the data the model reads.

Three of the nine findings in the 2026-09-13 dependency audit were the same defect seen
from different sides: `data/` lived outside `src/`, nothing declared it as package data,
and `maps.py` resolved it as `<package>/../../data` -- a path that exists in a source
checkout and nowhere else. `pip install .` therefore produced a `t700` whose first
`maps.f1()` raised `FileNotFoundError`, and the sdist shipped `tests/` without a single
CSV those tests read.

Nothing noticed for one reason: every way the project was ever exercised keeps the source
tree on `sys.path`. `pip install -e .` does, and so does pytest's `pythonpath = ["src"]`.
The bug was invisible to the entire test suite by construction, which is why these checks
read the *declarations* rather than the running import -- a test that imports `t700` and
finds its data is, in this tree, guaranteed to pass whether or not the wheel works.

What is not checked here: that a built wheel actually installs and runs. That needs a
build and a clean virtualenv, which is a minute of work rather than a millisecond, and it
belongs in whatever eventually runs CI. It was done by hand on 2026-09-13 -- wheel built,
installed into an empty venv outside the source tree, `trim.solve` run from `/tmp` -- and
`docs/notes/environment.md` records the result.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src" / "t700"
DATA = ROOT / "data"


@pytest.fixture(scope="module")
def pyproject() -> dict:
    with open(ROOT / "pyproject.toml", "rb") as fh:
        return tomllib.load(fh)


def test_the_data_directory_is_declared_as_package_data(pyproject: dict):
    """Without this mapping the wheel is a model with no maps in it."""
    tool = pyproject["tool"]["setuptools"]
    assert "t700.data" in tool["packages"], (
        "`t700.data` is not in [tool.setuptools].packages, so the repository's data/ "
        "directory is not part of the distribution. A wheel built from this tree would "
        "install a t700 that cannot load a single function table."
    )
    assert tool["package-dir"]["t700.data"] == "data", (
        "the t700.data package must be rooted at the repository's data/ directory"
    )
    globs = tool["package-data"]["t700.data"]
    assert any("csv" in g for g in globs), f"package-data for t700.data is {globs}"


def test_every_committed_data_file_is_covered_by_the_package_data_globs(pyproject: dict):
    """A glob that covers 61 of 62 CSVs fails in exactly one place, months later."""
    globs = pyproject["tool"]["setuptools"]["package-data"]["t700.data"]
    covered = {p for g in globs for p in DATA.glob(g)}
    on_disk = {p for p in DATA.rglob("*") if p.is_file() and p.suffix == ".csv"}
    assert on_disk, "no CSVs under data/ -- this test would pass vacuously"
    missing = sorted(p.relative_to(DATA) for p in on_disk - covered)
    assert not missing, (
        f"these committed data files match no package-data glob {globs} and would be "
        f"absent from a built wheel: {missing}"
    )


def test_the_sdist_manifest_ships_the_data_the_tests_read():
    """setuptools' defaults ship `tests/` and no CSV, which is worse than shipping neither."""
    manifest = ROOT / "MANIFEST.in"
    assert manifest.exists(), (
        "MANIFEST.in is gone. Without it setuptools falls back to its defaults, which "
        "shipped the test suite and none of the data it reads."
    )
    text = manifest.read_text()
    assert "recursive-include data" in text, "the sdist would carry no data/"
    assert "ballin-tm100991.pdf" in text, (
        "the sdist would carry no source document, and every citation in it is `pdf p.NN` "
        "against that one file"
    )


def test_no_core_module_resolves_data_by_walking_up_from_its_own_file():
    """The original defect, as a shape rather than as a path.

    `Path(__file__).parent.parent.parent / "data"` is the repository root's data/ only in
    a source checkout. `src/t700/_data.py` probes the two layouts that actually occur and
    raises a legible error otherwise; it is the one place allowed to know where data is.
    """
    for path in sorted(SRC.rglob("*.py")):
        if path.name == "_data.py":
            continue
        source = path.read_text()
        assert "parent.parent.parent" not in source, (
            f"{path.name} walks up from __file__ to find data/. That path exists only in "
            f"a source checkout -- import DATA_ROOT from t700._data instead."
        )


def test_installing_the_dev_extra_is_enough_to_run_the_suite(pyproject: dict):
    """`testpaths` includes validation/, and those modules import SciPy at module level.

    pytest collection is all-or-nothing: a module that fails to import during collection
    takes the run down with it, `tests/` included. So `pip install -e .[dev]` followed by
    `pytest` -- two lines anyone would try -- collected nothing at all.
    """
    extras = pyproject["project"]["optional-dependencies"]
    dev = " ".join(extras["dev"])
    testpaths = pyproject["tool"]["pytest"]["ini_options"]["testpaths"]
    if "validation" not in testpaths:
        pytest.skip("validation/ is no longer collected by default")
    assert "t700[tools]" in dev, (
        f"testpaths is {testpaths} and validation/ imports SciPy, which is declared in "
        f"the `tools` extra. `dev` is {extras['dev']} and does not pull it in, so "
        f"`pip install -e .[dev] && pytest` dies during collection."
    )
