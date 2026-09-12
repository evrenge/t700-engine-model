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


def data_writing_tools() -> set[str]:
    """Tools that write a committed file under `data/`, found rather than hardcoded.

    A hardcoded list is what let `digitize_appc_multi` sit outside the gate: it writes the
    seven Appendix C schedules and nobody noticed it was unlisted. `digitize.py` and
    `digitize_native.py` are shared libraries -- they name no output path, so they are
    excluded automatically rather than by exception.
    """
    out = set()
    for path in TOOLS.glob("*.py"):
        if re.search(r'Path\("data/|"data/\w+/', path.read_text()):
            out.add(path.stem)
    return out


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


def test_every_data_writing_tool_is_rerun_by_the_gate(source: str):
    """A tool outside the loop is a committed data file with no reproducibility check.

    This is how `data/schedules/` -- the seven Appendix C scheduling functions -- went
    unchecked: `digitize_appc_multi` writes them and was never added to the loop.
    """
    loop = source.split("for entry in", 1)[1].split("do", 1)[0]
    listed = set(re.findall(r"digitize_\w+", loop))
    writers = data_writing_tools()
    missing = writers - listed
    assert not missing, (
        f"these tools write committed files under data/ but are not rerun by the gate, "
        f"so their output is unverified: {sorted(missing)}"
    )
    for name in listed:
        assert (TOOLS / f"{name}.py").exists(), (
            f"the gate reruns tools/{name}.py, which does not exist"
        )
