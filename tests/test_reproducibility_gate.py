"""The data-reproducibility gate, `bash tools/reproduce_all.sh`, must be able to fail.

If no digitizer runs, `data/` stays clean and a verdict read from `git status` alone says
"every data file reproduces byte for byte". These tests read the script rather than run it
(that takes minutes and needs rendered pages). They cannot prove the digitizers correct,
only that a failure to run them cannot be reported as success.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
TOOLS = ROOT / "tools"
DATA = ROOT / "data"
SCRIPT = TOOLS / "reproduce_all.sh"


def data_tools() -> set[str]:
    """Tools that own a committed file under `data/`, found rather than hardcoded so that a
    new tool cannot sit outside the gate unnoticed.

    "Owns" rather than "writes": three of them check instead (see the next test). The glob
    does not reach the library in `tools/figure/`; `rectify_page` writes nothing.
    """
    out = set()
    for path in TOOLS.glob("*.py"):
        if path.stem == "rectify_page":
            continue
        if re.search(r'Path\("data/|"data/\w+/', path.read_text()):
            out.add(path.stem)
    return out


def test_every_data_tool_either_writes_its_file_or_checks_it():
    """A data tool must write its file or check it with `check_against_csv`.

    A tool that only prints its rows runs, succeeds and proves nothing. The hand-maintained
    `f2`, `f7` and `f9` are why some tools check rather than write.
    """
    for name in sorted(data_tools()):
        writes, checks = _writes_and_checks((TOOLS / f"{name}.py").read_text())
        assert writes or checks, (
            f"tools/{name}.py names a path under data/ but neither writes nor checks it; "
            f"the gate would prove nothing about that file."
        )


_WRITES = {"write", "write_text", "writerow", "writerows", "to_csv", "write_csv"}


def _writes_and_checks(source: str) -> tuple[bool, bool]:
    """Whether the tool *calls* a writer, and whether it calls `check_against_csv`.

    Read from the syntax tree, so a comment or docstring naming either does not count.
    """
    writes = checks = False
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", "")
        if name == "check_against_csv":
            checks = True
        elif name in _WRITES or (name == "save" and getattr(f.value, "id", "") == "np"):
            writes = True
        elif name == "open":
            pos = 0 if isinstance(f, ast.Attribute) else 1  # path.open(mode) vs open(p, mode)
            mode = node.args[pos] if len(node.args) > pos else None
            mode = next((k.value for k in node.keywords if k.arg == "mode"), mode)
            writes |= isinstance(mode, ast.Constant) and "w" in str(mode.value)
    return writes, checks


@pytest.fixture(scope="module")
def source() -> str:
    assert SCRIPT.exists(), f"{SCRIPT} is gone; the reproducibility gate has no script"
    return SCRIPT.read_text()


def _code(source: str) -> str:
    """The script without its comments, so a needle mentioned in prose cannot satisfy a guard.

    A `#` after whitespace starts a comment in bash unless it is quoted; cutting a quoted one
    too can only hide a needle, never invent one, so it errs toward failing.
    """
    lines = (ln for ln in source.splitlines() if not ln.lstrip().startswith("#"))
    return "\n".join(re.sub(r"\s#.*$", "", ln) for ln in lines)


_NO_POPPLER = (
    "the gate no longer checks that poppler is on PATH; a machine without it gets a "
    "FileNotFoundError per tool and no indication of what to install"
)


@pytest.mark.parametrize(
    "needle, why",
    [
        (
            "import numpy",
            "the script no longer verifies that its python3 can import NumPy. Without that "
            "check, running it outside the t700-dev container crashes every digitizer and "
            "reports success, which is the 2026-09-12 defect.",
        ),
        ("pdftoppm", _NO_POPPLER),
        ("pdfimages", _NO_POPPLER),
    ],
)
def test_the_gate_checks_its_prerequisites_before_doing_anything(source, needle, why):
    """NumPy and poppler are checked before the digitizer loop.

    Without poppler each tool fails with a FileNotFoundError that reads as a broken
    digitizer. Searched in the script's code only: its header comment names both binaries.
    """
    code = _code(source)
    assert needle in code, why
    assert code.index(needle) < code.index("for entry in "), (
        f"the {needle} check must come before the digitizer loop"
    )


def test_a_failed_digitizer_poisons_the_verdict(source: str):
    """A clean `data/` after a crash means nothing ran, not that the files are current."""
    assert re.search(r"failed=\$\(\(failed \+ 1\)\)", source), (
        "the script no longer counts failed digitizers; a crash would again be "
        "indistinguishable from a byte-identical reproduction"
    )
    # the phrase is also in the header comment; the last occurrence is the success echo
    success = source.rindex("reproduces byte for byte")
    check = source.index('if [ "$failed" -ne 0 ]')
    assert check < success, (
        "the failure count must be tested before the success message can be reached"
    )


def test_no_gate_tool_reads_a_gitignored_input_it_does_not_produce():
    """A tool that names a path under the gitignored `validation/out/` must also contain the
    command that creates it, or the gate runs only where the intermediates already exist."""
    ignored = (ROOT / ".gitignore").read_text()
    assert "validation/out/" in ignored, "validation/out/ is tracked now; revisit this test"

    for path in sorted(TOOLS.glob("*.py")):
        text = path.read_text()
        reads = set(re.findall(r'"(validation/out/[\w./-]*)"', text))
        if not reads:
            continue
        # Follow one hop of tools-local imports (`figure.raster` rasters for all of them).
        # Both spellings resolve: a sibling module under tools/ and a `figure` submodule.
        reachable = text
        for mod in re.findall(r"(?:^|\n)(?:from|import) ([\w.]+)", text):
            for candidate in (
                TOOLS / f"{mod.replace('.', '/')}.py",
                TOOLS / f"{mod}.py",
            ):
                if candidate.exists():
                    reachable += candidate.read_text()
                    break
        assert "pdftoppm" in reachable or "pdfimages" in reachable, (
            f"{path.name} reads {sorted(reads)} under the gitignored validation/out/ and "
            f"neither it nor anything it imports rasters a page, so it cannot run from a "
            f"fresh clone"
        )


def test_digitize_a1_cannot_rewrite_a_data_file_as_a_side_effect_of_being_imported():
    """`digitize_a1` works at module level, so it refuses to be imported rather than rewrite
    `data/maps/f1_compressor_mass_flow.csv` on import."""
    source = (TOOLS / "digitize_a1.py").read_text()
    assert '__name__ != "__main__"' in source, (
        "tools/digitize_a1.py no longer refuses to be imported; importing it rewrites data."
    )
    guard = source.index('__name__ != "__main__"')
    writes = source.index('"data/maps')
    assert guard < writes, "the refusal must come before anything that writes"


def test_every_data_writing_tool_is_rerun_by_the_gate(source: str):
    """A tool outside the loop is a committed data file with no reproducibility check."""
    loop = source.split("for entry in", 1)[1].split("do", 1)[0]
    listed = set(re.findall(r"digitize_\w+", loop))
    writers = data_tools()
    missing = writers - listed
    assert not missing, (
        f"these tools write committed files under data/ but are not rerun by the gate, "
        f"so their output is unverified: {sorted(missing)}"
    )
    for name in listed:
        assert (TOOLS / f"{name}.py").exists(), (
            f"the gate reruns tools/{name}.py, which does not exist"
        )


def test_the_gate_reaches_every_committed_data_file(source: str):
    """`digitize_single_curve.main()` with no arguments, as the gate calls it, runs every
    figure in `FIGS`: listing a tool is not running the part that writes your file."""
    src = (TOOLS / "digitize_single_curve.py").read_text()
    start = src.index("FIGS = {")
    figs = src[start : src.index("\n}\n", start)]
    owned = set(re.findall(r'^    "([a-z0-9]+)":', figs, re.M))
    assert len(owned) >= 13, f"FIGS parse looks wrong, found {sorted(owned)}"

    default = re.search(r"keys = \[a for a in argv if a in FIGS\] or (.+)", src)
    assert default, "main()'s default figure list has changed shape; re-read it"
    assert "list(FIGS)" in default.group(1), (
        f"digitize_single_curve.main() defaults to {default.group(1).strip()}, not every "
        f"figure in FIGS; the gate runs it with no arguments."
    )


def test_no_committed_data_file_names_a_tool_that_does_not_exist():
    """A provenance header is a claim about a file on disk. Renames must carry it."""
    bad = []
    for csv in sorted(DATA.rglob("*.csv")):
        for m in re.finditer(r"tools/(\w+)\.py", csv.read_text()):
            if not (TOOLS / f"{m.group(1)}.py").exists():
                bad.append(f"{csv.relative_to(DATA.parent)} -> tools/{m.group(1)}.py")
    assert not bad, (
        f"committed data files cite tools that no longer exist: {sorted(set(bad))}. "
        f"Rerun the digitizers; never hand-edit the header."
    )
