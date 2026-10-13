# CLAUDE.md — T700 Replication

## What this is

A Python replication of:

> Ballin, M. G., *A High Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
> NASA TM-100991, July 1988.

A component-level model of the General Electric T700-GE-700 turboshaft engine and its
fuel control system. The source document lives at `docs/ballin-tm100991.pdf` and is
called **the report** throughout this file and in commit messages.

The goal is **replication, not improvement**. Where the report is dated, simplified, or
arguably wrong, reproduce it as written and record the objection in
`docs/notes/open-questions.md`. Improvements come after the replication validates.

## Commands

Everything runs inside the `t700-dev` distrobox (Fedora): the Bazzite host is immutable and
has no NumPy. Either `distrobox enter t700-dev` first, or prefix with
`distrobox enter t700-dev --`. `$HOME` is shared, so the project path is the same inside.

| Task | Command |
|---|---|
| Fast unit tests | `pytest tests -q` |
| Report-comparison tests | `pytest validation -q` |
| Everything | `pytest` |
| Lint and format | `ruff check --fix . && ruff format .` |
| Find which page holds a term | `grep -rin "<term>" docs/extracted/` |
| Read that page (the real source) | `pdftoppm -f NN -l NN -r 200 -png docs/ballin-tm100991.pdf /tmp/p` |
| Compare model against the report | `/validate` |
| Prove every digitized file is current | `bash tools/reproduce_all.sh` |

Page images `docs/extracted/p-0NN.png` (pages 8-101) are gitignored, so a fresh clone or a
worktree has only the `.txt` text layer, which the OCR corollary forbids transcribing from.
Render what you need with the `pdftoppm` line above.

**The edit hook.** `.claude/settings.json` runs two PostToolUse hooks after every Edit or
Write of a `.py` file: `ruff check --fix` then `ruff format` on that file (silent; it never
blocks), then, for a file matching `*/t700/*.py` or `*/tests/*.py`, `pytest tests -q -x` in
the container, in the session's project directory. A failing test run exits 2, and Claude
Code reports its last 25 lines as a blocking error: the edit is already made, so fix the
failure before going on.

## The provenance rule

This is the rule that matters most here. Every number that enters the codebase — a
constant, a map point, a schedule breakpoint, a gain, a limit, a time constant — carries
a citation to the report:

```python
# [TM-100991 pdf p.55, Table A.1] fuel lower heating value
HVF_BTU_PER_LBM = 18300.0
```

**Citations use PDF page numbers**, written `pdf p.NN`. The scan has no usable printed
page numbers in its text layer, and the cover says "95 p" against 104 PDF pages, so the
two numberings differ — never mix them.

Data files carry the same, in a header block: source page, figure or table number, and
how the value was obtained (`transcribed` / `digitized` / `derived`, with the derivation).

Never supply a value from general turbomachinery knowledge, from another engine model,
from a T700 spec sheet, or from a plausible-looking guess. If the report does not give it:

1. Search the report again — it is very likely in Appendix A (constants and function
   tables), Appendix C (fuel control), or on a figure axis. `docs/notes/index.md` has the
   section map.
2. If it is genuinely absent: add an entry to `docs/notes/open-questions.md`, mark the
   site `# UNVERIFIED:` with what is missing, and **stop and ask**.

A model that runs on invented numbers is worse than one that does not run, because it
looks like it works.

### The OCR corollary

The text layer of this 1988 scan corrupts symbols and digits: `K` reads as `g` or `i_`,
subscripts flatten or vanish, `lbf` reads as `Ibf`/`lb!`, minus signs drop. **Never
transcribe a number from the text layer.** Use it to locate; render the page (command
above) and read the PNG. `docs/notes/index.md` records the observed damage patterns.

## Units

Internal units are **the report's units** — US customary: lbm, lbm/s, psia, °R, hp,
ft·lbf, rpm, in². Not SI. Every printed table and figure can then be compared to model
output literally, with no conversion in the loop waiting to be wrong. Conversions happen
only at user-facing I/O boundaries; reversing this would be expensive.

- `t700.units` holds every conversion factor and nothing else computes one inline.
- Any quantity whose unit is not obvious carries a suffix: `t45_degR`, `wf_pps`,
  `ps3_psia`, `q_req_ftlbf`.

## Naming follows the report

Variables use the report's own symbols, so the code reads against the page: `NG`, `NP`,
`WF`, `T45`, `PS3`, `QREQ` become `ng`, `np_`, `wf`, `t45`, `ps3`, `qreq`; constants keep
the report's capitals (`K_V41`, `TC_T41`). `docs/notes/symbols.md` is the authoritative
symbol table, read from page images.

Four name collisions, resolved once here so they stay resolved:

| Report symbol | Python name | Why |
|---|---|---|
| NP | `np_` | `np` is NumPy |
| h (convective coefficient) | `h_conv` if ever needed | collides with the enthalpies `h2`, `h3`, `h41`. Unused today: Eq. 51 gives only the grouping as a time constant |
| M (heat-sink metal mass) | `m_metal` if ever needed | bare `m` reads as a generic mass. Unused today, as for `h` |
| J (inertia) | `J_GT`, `J_PT`, `J_LOAD_UH60A` | never bare `j` |

Two typesetting traps the OCR layer cannot show, both load-bearing:

- **Station indices are set full size** (`T41` is T-4-1 inline); **true subscripts are set
  small** (`H41_ns`). Different things.
- **Sixteen constants carry a subscript on a subscript** — `K_H41₁` and `K_H41₂` are the
  slope and intercept of one linear fit and have **different units**. There is no plain
  `K_H41`. Flattening a pair to a single name silently corrupts the fit.

## Architecture constraints

- **`src/t700/` imports NumPy and the stdlib. Nothing else.** No SciPy, Matplotlib,
  pandas, or Numba in the model core. Analysis, digitizing, and plotting live in `tools/`
  and `validation/`, where SciPy and Matplotlib are welcome.
- **All gas property evaluation goes through `t700.thermo`.** No inline `cp`, `gamma`,
  `h`, or T-from-h arithmetic anywhere else, ever. The backend is pure NumPy; a C++
  lookup-table backend is deferred (`SCOPE.md` → Deferred), and this interface is what lets
  it drop in without touching component code.
- **Fixed-step explicit integration only in the core**, matching the report's real-time
  formulation. An adaptive solver would hide the behaviour we are trying to reproduce.
  SciPy adaptive integration may be used in `validation/` to confirm the fixed-step scheme
  is converged — never to run the model.
- **Determinism.** No wall-clock reads, no RNG, no reliance on dict ordering in the core.
  Same input → bit-identical output, run to run and machine to machine.
- Components are plain functions or small classes: explicit state in, derivatives out. No
  global mutable state. The table loaders in `maps.py`, `control/schedules.py` and
  `appendix_b.py` are `functools.cache`d: each table loads once and is read-only after.
  Beyond those, three sanctioned exceptions:
  - the thermo backend (`thermo.set_backend`);
  - `maps.CONDITIONING`, write-once: the cached loaders record how far physical
    conditioning moved each table, and nothing writes it after load;
  - `maps._clamps`, the live one: a diagnostic counter of which tables were read outside
    their data. Nothing downstream of a lookup reads it, and
    `tests/test_maps.py::test_the_clamp_counter_cannot_influence_any_model_output` checks
    that behaviourally. The cost is that `frame()` is impure. Use `maps.clamp_scope()`
    rather than `reset_clamps()` wherever a measurement may nest inside another.
- No abstraction added "for later" (`advisory`). The thermo backend interface is the one
  exception.

## Helpers

- `.claude/agents/doc-extract.md` for anything that comes out of the PDF: it returns values
  with page citations. Do not page through 100 pages by hand into context.
- `.claude/agents/physics-review.md`: run before calling a component module done.
- `/validate`: the model against the digitized reference data, as a deviation table and plots.
- **Digitizing anything: `tools/figure/` is the library, one script per figure above it.**
  Read `docs/notes/digitizing-recipe.md` first — twelve traps, each of which cost a wrong
  CSV that looked right. Never add a second copy of something already in there; a test
  fails if you do.
- Before believing a map, run `tools/check_straightness.py` (is the printed curve straight
  between its markers, in data space) and `tools/compare_models.py` (what the frame model and
  the anchoring rule cost, per figure). Both write nothing.
- `tools/constant_sweep.py` audits the tests: it moves each model constant by ±5 % in a
  scratch worktree and reports which ones no test notices. A long job (one suite per constant).
- `tools/model_digest.py` hashes every model output. Run it before and after a refactor and
  `--compare` the two files to prove that no number moved.

## Non-negotiables, and what actually enforces them

A rule that lives only in this file is a wish. Each one below names its mechanism. A rule
marked `advisory` has nothing behind it but judgment — those are the ones that slip late
in a long session, so mechanize any that proves to matter.

| Rule | Enforced by | What the mechanism still misses |
|---|---|---|
| Every number cites a report page | `advisory` | everything — challenge any bare constant on sight |
| Never tune a constant to make a validation test pass | `advisory` | a diff touching `data/` in the same breath as a failing test is the tell |
| Never silently widen a tolerance | `tests/test_scope_tolerances.py` — parses `SCOPE.md`'s three tolerance tables, resolves each `` `test_x.SYMBOL` `` a row names, fails if the two numbers differ or a row names no symbol | nothing structural. Widening is allowed; doing it in one place only is not |
| No SciPy / Matplotlib / pandas / Numba under `src/t700/` | `tests/test_architecture.py` — parses every core module's imports | |
| Gas properties only via `t700.thermo` | `tests/test_architecture.py`, in two layers — any module outside `thermo/` that *names* a `K_H*`/`K_T*`/`K_TH*` constant, or contains one of their *values* as a bare float | **a rearranged fit** whose printed values no longer appear (`h41 * 3.3223 + ...`). This rule wants reading, not only running |
| Core is deterministic | `tests/test_architecture.py` structurally (clocks, RNGs, `FORBIDDEN_ATTRIBUTES`, and `dotted_attributes` to resolve import aliases), with decoys proving each detector fires; `tests/test_determinism.py` behaviourally, within one process and across fresh interpreters with different hash seeds | **a second machine.** Nothing runs the model on other hardware or other library builds |
| No global mutable state in the core | `advisory`, with three named exceptions — the thermo backend, `maps._clamps`, and the write-once `maps.CONDITIONING` — beside the read-only cached loaders | the clamp counter is bounded behaviourally instead: no model output can depend on it, over a trim (`test_maps.py`) and over a whole transient (`test_determinism.py`) |
| No force-push | `permissions.ask` | plain `git push` is authorized; force-push still prompts |
| Never hand-edit a generated data file | `bash tools/reproduce_all.sh` — reruns every digitizer and diffs | **three of the eleven engine maps are not generated**: `f2`, `f7` and `f9` are hand-maintained CSVs whose headers carry provenance no generator writes, so their tools *check* rather than write. **Run it inside the container, and read its verdict rather than its exit code.** A hand-added annotation is wiped by the next rerun — put it in the tool |
| One implementation per concept | `tests/test_tool_duplication.py` — AST fingerprint of every function body in `tools/`, `validation/`, `tests/` and `src/`, local names and annotations normalised, copies within one file included, with decoys | **the copy that has already been edited.** It catches a copy before it diverges, which is the only moment that is cheap; after that it takes a reader |
| A test or name a sentence cites exists | `tests/test_citations.py` — finds every `<file>.py::<name>`, `test_<module>.<name>` and backticked bare `test_<name>` in the markdown and every comment and docstring under `src/`, `tests/`, `validation/`, `tools/`, and resolves each against an AST index of those folders. Placeholders are an `ALLOWED` list, each with its reason; decoys prove a made-up name fails | **any other citation form** — `maps.clamp_scope()`, a test named in prose without backticks — and a citation that resolves to the *wrong* existing name. It proves the name exists, not that it still says what the sentence claims |
| A measurement the docs quote is today's value | `tests/test_live_numbers.py` — a live number in `README.md` or `SCOPE.md` carries an invisible `<!--rj:<path>-->` marker naming its source in `results/report.json`, and must equal that value at the printed decimals, sign included unless the marker adds `abs`; decoys prove a wrong number fails. | **a number nobody marked.** Marking is opt-in, so a newly quoted measurement is unguarded until someone marks it, and one with no report.json source (frame cost, read-error floors, the discrete-map eigenvalues) cannot be marked at all: describe it in words instead. Nor the sentence around a marked number: a current value can still be described wrongly |
| Python stays formatted, linted and passing the fast tests | the PostToolUse hooks (see Commands) | ruff's verdict is discarded, so a lint it cannot fix passes silently; the tests run in the session's project directory, not in a worktree elsewhere; `validation/` runs only when you run it |

The two that remain `advisory` are the two that matter most, and they cannot be
mechanised — a machine cannot tell a cited constant from a plausible one. Be
correspondingly careful. Every mechanism here was added after its rule had been broken;
trust one no further than its third column says.

## Communication

Terse, no preamble. Report deviations with numbers, not adjectives — "NG is 2.4% high
against a 1% tolerance", not "close but slightly off". Say plainly what was not checked.

**Writing in the tree** (docs, docstrings, comments). Say what a thing is and the one
non-obvious reason, with its `pdf p.NN`, then stop. No history ("used to", "on 2026-09-14",
superseded findings): git and the ledger keep it. One home per explanation, a pointer
elsewhere ("see #63"). A measured number is marked live (`<!--rj:path-->` in README/SCOPE) or
described in words; commit messages are dated records and may carry numbers. Quote the report
exactly as the page image prints it, with the page and `[sic]` for its typos. CSV headers,
report.json text, `Physics(why=...)` and text a test searches for are data, not prose: never
reword them.
