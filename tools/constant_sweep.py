"""Which model constants does the test suite actually notice? Change each one and find out.

Every numeric constant in `t700.constants` and `t700.control.constants` is moved by +5 %
(a zero by `ZERO_STEP`), the whole suite runs, and the failing tests are recorded; the ones
nothing caught are retried at -5 %, because a limit can bind in one direction only. A
constant is

  * **behavioural** if a failing test does not name it -- a model output moved;
  * **value-only** if every failing test reads it straight off its module (`k.CORR`), as a
    transcription check does, which says nothing about whether the model depends on it;
  * **unnoticed** if nothing fails.

It runs in scratch git worktrees of HEAD, one per worker, so the real tree is never written.

    python3 tools/constant_sweep.py                 # all constants, 4 workers
    python3 tools/constant_sweep.py --workers 2 CR CORR
    python3 tools/constant_sweep.py --keep          # leave the worktrees and the raw results

A long job; each worker is one single-threaded pytest. The worktrees go in the system temp
directory: set TMPDIR to a disk path where /tmp is RAM.
"""

from __future__ import annotations

import argparse
import json
import os
import queue
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FILES = ("src/t700/constants.py", "src/t700/control/constants.py")
ASSIGN = re.compile(r"^([A-Z][A-Z0-9_]*)\s*(?::[^=]*)?=\s*(-?[0-9][0-9_.]*(?:[eE][-+]?\d+)?)\b")
ZERO_STEP = {"YLOLIM": 2.0}
"""Additive steps for constants whose value is zero, sized to their companion: YLOLIM is
the lower end of a 0-40 range, so 2.0 is 5 % of it. Anything else at zero moves by 0.05."""


def constants(tree: Path, only: set[str]) -> list[tuple[str, str, float]]:
    out = []
    for f in FILES:
        for line in (tree / f).read_text().splitlines():
            m = ASSIGN.match(line)
            if m and (not only or m.group(1) in only):
                out.append((f, m.group(1), float(m.group(2).replace("_", ""))))
    return out


def run_case(tree: Path, f: str, name: str, value: float, factor: float) -> list[str]:
    path = tree / f
    original = path.read_text()
    new = value + ZERO_STEP.get(name, 0.05) if value == 0 else value * factor
    pattern = rf"^({re.escape(name)}\s*(?::[^=]*)?=\s*)(-?[0-9][0-9_.]*(?:[eE][-+]?\d+)?)"
    text, n = re.subn(pattern, lambda m: m.group(1) + repr(new), original, count=1, flags=re.M)
    if n != 1:
        raise RuntimeError(f"{name}: assignment not found in {f}")
    path.write_text(text)
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    try:
        r = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", "-o", "addopts=", "-rf"]
            + ["-p", "no:cacheprovider"],
            cwd=tree,
            env=env,
            capture_output=True,
            text=True,
            timeout=900,
        )
        return sorted({ln.split()[1] for ln in r.stdout.splitlines() if ln.startswith("FAILED ")})
    except subprocess.TimeoutExpired:
        return ["TIMEOUT"]
    finally:
        path.write_text(original)


def sweep(trees: list[Path], items, factor: float) -> dict[str, list[str]]:
    q: queue.Queue = queue.Queue()
    for it in items:
        q.put(it)
    results: dict[str, list[str]] = {}
    lock = threading.Lock()

    def worker(tree: Path) -> None:
        while True:
            try:
                f, name, value = q.get_nowait()
            except queue.Empty:
                return
            failed = run_case(tree, f, name, value, factor)
            with lock:
                results[name] = failed
                print(f"  x{factor:g}  {name:26s} {len(failed):4d} failed", flush=True)

    threads = [threading.Thread(target=worker, args=(t,)) for t in trees]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    return results


def names_it(tree: Path, nodeid: str, name: str) -> bool:
    """Whether the failing test's own function reads the constant off a module (`k.CORR`).

    Attribute access, not the bare word: a block test keyed on the printed value
    (`PRINTED["CORR"]`) is behavioural, and a docstring may name anything."""
    path, _, rest = nodeid.partition("::")
    func = rest.split("[")[0].split("::")[-1]
    src = (tree / path).read_text()
    i = src.find(f"def {func}(")
    if i < 0:
        return False
    m = re.search(r"^(def |class |@pytest)", src[i + 1 :], re.M)
    body = src[i : i + 1 + m.start()] if m else src[i:]
    return re.search(rf"\.{re.escape(name)}\b", body) is not None


def classify(tree: Path, name: str, failed: list[str]) -> str:
    if not failed:
        return "unnoticed"
    if any(not names_it(tree, f, name) for f in failed if f != "TIMEOUT"):
        return "behavioural"
    return "value-only"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*", help="only these constants")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("JOBS", "4")))
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args(argv)

    # a polite stop still restores and removes the worktrees; SIGKILL cannot, and leaves
    # them for `git worktree prune` once the scratch directory is deleted
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))
    scratch = Path(tempfile.mkdtemp(prefix="t700-constant-sweep-"))
    trees = [scratch / f"wt{i}" for i in range(max(1, args.workers))]
    for t in trees:
        subprocess.run(
            ["git", "worktree", "add", "-q", "--detach", str(t), "HEAD"], cwd=ROOT, check=True
        )
    try:
        items = constants(trees[0], set(args.names))
        print(f"{len(items)} constants, {len(trees)} workers, scratch {scratch}")
        up = sweep(trees, items, 1.05)
        quiet = [it for it in items if classify(trees[0], it[1], up[it[1]]) != "behavioural"]
        down = sweep(trees, quiet, 0.95)
        verdict = {}
        for _f, name, _v in items:
            classes = {classify(trees[0], name, up[name])}
            if name in down:
                classes.add(classify(trees[0], name, down[name]))
            verdict[name] = (
                "behavioural"
                if "behavioural" in classes
                else ("value-only" if "value-only" in classes else "unnoticed")
            )
        for cls in ("unnoticed", "value-only"):
            names = sorted(n for n, v in verdict.items() if v == cls)
            print(f"\n{cls} ({len(names)}): {' '.join(names)}")
        print(f"\nbehavioural: {sum(v == 'behavioural' for v in verdict.values())}")
        if args.keep:
            (scratch / "results.json").write_text(json.dumps({"up": up, "down": down}, indent=1))
            print(f"kept: {scratch}")
    finally:
        if not args.keep:
            for t in trees:
                subprocess.run(["git", "worktree", "remove", "--force", str(t)], cwd=ROOT)
            shutil.rmtree(scratch, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
