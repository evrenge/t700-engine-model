"""The distribution must contain the data the model reads.

`data/` lives outside `src/`. Every way the project is exercised (`pip install -e .`,
pytest's `pythonpath = ["src"]`) keeps the source tree on `sys.path`, so a test that imports
`t700` and finds its data passes whether or not the wheel works. These checks read the
declarations instead. Not checked: that a built wheel installs and runs in a clean venv
(see `docs/notes/environment.md`).
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
        "`t700.data` is not in [tool.setuptools].packages, so a wheel would carry no data/"
    )
    assert tool["package-dir"]["t700.data"] == "data", (
        "the t700.data package must be rooted at the repository's data/ directory"
    )


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
    assert manifest.exists(), "MANIFEST.in is gone; setuptools' defaults ship no data/"
    text = manifest.read_text()
    assert "recursive-include data" in text, "the sdist would carry no data/"
    assert "ballin-tm100991.pdf" in text, (
        "the sdist would carry no source document, and every citation in it is `pdf p.NN` "
        "against that one file"
    )


def test_no_core_module_resolves_data_by_walking_up_from_its_own_file():
    """Walking up from `__file__` finds data/ only in a source checkout. `src/t700/_data.py`
    probes the layouts that occur and is the one place allowed to know where data is."""
    for path in sorted(SRC.rglob("*.py")):
        if path.name == "_data.py":
            continue
        source = path.read_text()
        assert "parent.parent.parent" not in source, (
            f"{path.name} walks up from __file__ to find data/; import t700._data.DATA_ROOT."
        )


def test_installing_the_dev_extra_is_enough_to_run_the_suite(pyproject: dict):
    """`testpaths` includes validation/, whose modules import SciPy at module level, and one
    failed import during collection takes the whole run down."""
    extras = pyproject["project"]["optional-dependencies"]
    dev = " ".join(extras["dev"])
    testpaths = pyproject["tool"]["pytest"]["ini_options"]["testpaths"]
    if "validation" not in testpaths:
        pytest.skip("validation/ is no longer collected by default")
    assert "t700[tools]" in dev, (
        f"testpaths is {testpaths} and validation/ needs SciPy from the `tools` extra; "
        f"`dev` is {extras['dev']} and does not pull it in."
    )
