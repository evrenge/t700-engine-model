"""The data-reproducibility gate must be able to fail.

`bash tools/reproduce_all.sh` is one of only two mechanisms CLAUDE.md marks `live` behind
the rule that a generated data file is never hand-edited. On 2026-09-12 it was found to be
incapable of failing:

  * it invoked `python3`, and run from the Bazzite host that interpreter has no NumPy, so
    every digitizer died on `import numpy`;
  * its verdict was computed from `git status --porcelain data/` alone -- and `data/` was
    clean *because nothing had written to it*.

So a total failure to execute printed "every data file reproduces byte for byte". Every
number downstream of `data/` had been resting on that sentence.

These tests are static: they read the script rather than run it, because running it takes
several minutes and needs the rendered page images. They check the two properties whose
absence made the gate vacuous. They cannot prove the digitizers are correct -- only that a
failure to run them can no longer be reported as success.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "tools" / "reproduce_all.sh"

TOOLS = Path(__file__).resolve().parent.parent / "tools"
DATA = Path(__file__).resolve().parent.parent / "data"


def data_tools() -> set[str]:
    """Tools that own a committed file under `data/`, found rather than hardcoded.

    A hardcoded list is what let `digitize_multi_curve` sit outside the gate: it writes the
    seven Appendix C schedules and nobody noticed it was unlisted. `digitize.py` and
    `digitize_native.py` are shared libraries -- they name no output path, so they are
    excluded automatically rather than by exception.

    "Owns" rather than "writes", because three of them do not write. See
    `test_every_data_tool_either_writes_its_file_or_checks_it`.
    """
    out = set()
    for path in TOOLS.glob("*.py"):
        if path.stem in ("digitize", "digitize_native", "rectify_page"):
            continue
        if re.search(r'Path\("data/|"data/\w+/', path.read_text()):
            out.add(path.stem)
    return out


def test_every_data_tool_either_writes_its_file_or_checks_it():
    """A tool that only *prints* its rows puts a human eye between figure and data.

    `digitize_a2`, `digitize_a7` and `digitize_a9` contained no write call at all. They
    printed their extraction and ended "header and CSV are maintained by hand in ...;
    compare the rows above" -- so `tools/reproduce_all.sh` ran them, they succeeded, they
    touched nothing, `git status data/` stayed clean, and the gate reported "every data
    file reproduces byte for byte" for `f2`, `f7` and `f9`. Demonstrated by the 2026-09-13
    accuracy audit by corrupting a value in `f9_pt_mass_flow.csv` and watching the gate
    pass. That is the 2026-09-12 defect CLAUDE.md records as closed, alive for three of
    the eleven engine maps -- `f9` among them, which is the map Eq. 80's repelling fixed
    point depends on.

    They now call `digitize.check_against_csv` and return non-zero on any difference, so
    the CSVs stay hand-maintained and their numbers stop being unverified. This test is
    what keeps that true: a data tool must either write its file or check it.

    It also replaces a detector that was fooled by its own subject matter --
    `data_tools()` matches a `data/` path anywhere in the source, and in those three files
    the only such path was inside the "maintained by hand" message.
    """
    for name in sorted(data_tools()):
        text = (TOOLS / f"{name}.py").read_text()
        writes = bool(
            re.search(r"\.write\(|\.write_text\(|writerow|to_csv|np\.save|open\([^)]*[\"']w", text)
        )
        checks = "check_against_csv" in text
        assert writes or checks, (
            f"tools/{name}.py names a path under data/ but neither writes it nor checks "
            f"it against the extraction. The reproducibility gate will run it, it will "
            f"succeed, and it will prove nothing about that file."
        )


@pytest.fixture(scope="module")
def source() -> str:
    assert SCRIPT.exists(), f"{SCRIPT} is gone; the reproducibility gate has no script"
    return SCRIPT.read_text()


def test_the_gate_checks_its_interpreter_before_doing_anything(source: str):
    """Refusing an interpreter without NumPy is what makes the rest of the run meaningful."""
    assert "import numpy" in source, (
        "the script no longer verifies that its python3 can import NumPy. Without that "
        "check, running it outside the t700-dev container crashes every digitizer and "
        "reports success, which is the 2026-09-12 defect."
    )
    guard = source.index("import numpy")
    loop = source.index("for entry in ")
    assert guard < loop, "the NumPy check must come before the digitizer loop"


def test_a_failed_digitizer_poisons_the_verdict(source: str):
    """A clean `data/` after a crash means nothing ran, not that the files are current."""
    assert re.search(r"failed=\$\(\(failed \+ 1\)\)", source), (
        "the script no longer counts failed digitizers; a crash would again be "
        "indistinguishable from a byte-identical reproduction"
    )
    # the phrase also appears in the script's own header comment describing the defect,
    # so take the last occurrence -- the echo that actually declares success
    success = source.rindex("reproduces byte for byte")
    check = source.index('if [ "$failed" -ne 0 ]')
    assert check < success, (
        "the failure count must be tested before the success message can be reached"
    )


def test_the_gate_checks_for_poppler_before_doing_anything(source: str):
    """The prerequisite a Python dependency list cannot express.

    Every digitizer shells out to `pdftoppm` or `pdfimages`. Without them the tools fail
    one at a time with a FileNotFoundError naming a binary, which reads as a broken
    digitizer rather than a missing package -- and the string "poppler" appeared nowhere
    in this repository until 2026-09-13.
    """
    assert "pdftoppm" in source and "pdfimages" in source, (
        "the gate no longer checks that poppler is on PATH; a machine without it gets a "
        "FileNotFoundError per tool and no indication of what to install"
    )
    guard = min(source.index("pdftoppm"), source.index("pdfimages"))
    loop = source.index("for entry in ")
    assert guard < loop, "the poppler check must come before the digitizer loop"


def test_no_gate_tool_reads_a_gitignored_input_it_does_not_produce():
    """A gate that only runs on the one machine that still has the intermediates.

    `digitize_a1.py` read `validation/out/digitize/a1r_600-056.png` -- gitignored, and
    produced by no tool in the gate. On a fresh clone it failed loudly rather than falsely
    passing, which is the better of the two failures, but CLAUDE.md names this gate as one
    of only two `live` mechanisms and it was unrunnable anywhere else.

    The rule this enforces: a tool that names a path under `validation/out/` must also
    contain the command that creates it.
    """
    ignored = (Path(__file__).resolve().parent.parent / ".gitignore").read_text()
    assert "validation/out/" in ignored, "validation/out/ is tracked now; revisit this test"

    for path in sorted(TOOLS.glob("*.py")):
        text = path.read_text()
        reads = set(re.findall(r'"(validation/out/[\w./-]*)"', text))
        if not reads:
            continue
        # A tool may raster the page itself or delegate to a sibling that does --
        # `digitize_native.native_bitmap` is shared by five of them -- so follow one hop
        # of tools-local imports before deciding it cannot produce its own inputs.
        reachable = text
        for sibling in re.findall(r"(?:^|\n)(?:from|import) (\w+)", text):
            candidate = TOOLS / f"{sibling}.py"
            if candidate.exists():
                reachable += candidate.read_text()
        assert "pdftoppm" in reachable or "pdfimages" in reachable, (
            f"{path.name} reads {sorted(reads)} under the gitignored validation/out/ and "
            f"neither it nor anything it imports rasters a page, so it cannot run from a "
            f"fresh clone"
        )


def test_digitize_a1_cannot_rewrite_a_data_file_as_a_side_effect_of_being_imported():
    """It is the one tool whose work happens at module level.

    `python3 -c "import digitize_a1"` used to rerun the whole Figure A1 extraction and
    overwrite `data/maps/f1_compressor_mass_flow.csv`. Every other tool in tools/ has an
    `if __name__ == "__main__"` guard; this one has nothing to put behind such a guard, so
    it refuses the import instead.
    """
    source = (TOOLS / "digitize_a1.py").read_text()
    assert '__name__ != "__main__"' in source, (
        "tools/digitize_a1.py no longer refuses to be imported. Everything it does "
        "happens at module level, including rewriting a committed data file."
    )
    guard = source.index('__name__ != "__main__"')
    writes = source.index('"data/maps')
    assert guard < writes, "the refusal must come before anything that writes"


def test_every_data_writing_tool_is_rerun_by_the_gate(source: str):
    """A tool outside the loop is a committed data file with no reproducibility check.

    This is how `data/schedules/` -- the seven Appendix C scheduling functions -- went
    unchecked: `digitize_multi_curve` writes them and was never added to the loop.
    """
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
    """Listing a tool is not the same as running the part of it that writes your file.

    `digitize_single_curve.py` owns thirteen figures. Its `main()` defaulted to
    `["a8", "a10", "a6"]` -- the three the old filename `digitize_a6810.py` named -- while
    `FIGS` grew to thirteen, and the gate invokes it with no arguments. So ten committed
    files (`f3`, `f4`, `f5`, `f_hs`, `F_HM1`-`F_HM6`) carried a header reading "reruns and
    reproduces this file" that nothing had ever executed.

    The previous check only asked whether each data-writing tool appears in the gate's
    loop. This one asks whether the *default invocation* covers everything the tool owns,
    which is the question that matters: a tool in the list that silently does a third of
    its job leaves the gate reporting success over files it never touched.

    Found 2026-09-14, when renaming the tool rewrote the provenance line in three files
    and left the other ten reading the old name after a full gate run.
    """
    src = (TOOLS / "digitize_single_curve.py").read_text()
    figs = src[src.index("FIGS = {") : src.index("EDGE = {")]
    owned = set(re.findall(r'^    "([a-z0-9]+)":', figs, re.M))
    assert len(owned) >= 13, f"FIGS parse looks wrong, found {sorted(owned)}"

    default = re.search(r"keys = \[a for a in argv if a in FIGS\] or (.+)", src)
    assert default, "main()'s default figure list has changed shape; re-read it"
    assert "list(FIGS)" in default.group(1), (
        f"digitize_single_curve.main() defaults to {default.group(1).strip()} rather than "
        f"every figure in FIGS. The gate runs it with no arguments, so anything outside "
        f"that default is a committed file whose provenance header is never executed."
    )


def test_no_committed_data_file_names_a_tool_that_does_not_exist(source: str):
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
