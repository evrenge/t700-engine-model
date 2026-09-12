# CLAUDE.md — T700 Replication

## What this is

A Python replication of:

> Ballin, M. G., *A High-Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
> NASA TM-100991, July 1988.

A component-level model of the General Electric T700-GE-700 turboshaft engine and its
fuel control system. The source document lives at `docs/ballin-tm100991.pdf` and is
called **the report** throughout this file and in commit messages.

The goal is **replication, not improvement**. Where the report is dated, simplified, or
arguably wrong, reproduce it as written and record the objection in
`docs/notes/open-questions.md`. Improvements come after the replication validates.

## Commands

Everything runs inside the container. Either `distrobox enter t700-dev` first, or prefix
with `distrobox enter t700-dev --`.

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

Pre-rendered page images already exist for **every page from 8 to 101** (94 files) as
`docs/extracted/p-0NN.png`. Render anything outside that range with the `pdftoppm` line
above.

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
transcribe a number from the text layer.** Use it to locate; render the page and read it:

```bash
pdftoppm -f 55 -l 55 -r 200 -png docs/ballin-tm100991.pdf /tmp/p   # then Read the PNG
```

`docs/notes/index.md` records the observed damage patterns.

## Units

Internal units are **the report's units** — US customary: lbm, lbm/s, psia, °R, hp,
ft·lbf, rpm, in². Not SI.

The reason: every printed table and figure can then be compared to model output
literally, with no conversion sitting in the loop waiting to be wrong. Conversions happen
only at user-facing I/O boundaries.

- `t700.units` holds every conversion factor and nothing else computes one inline.
- Any quantity whose unit is not obvious carries a suffix: `t45_degR`, `wf_pps`,
  `ps3_psia`, `q_ftlbf`.

This decision is reversible but expensive to reverse. Raise it before Phase 3, not after.

## Naming follows the report

Variables use the report's own symbols, so the code reads against the page: `NG`, `NP`,
`WF`, `T45`, `PS3`, `QREQ`, `LDS` become `ng`, `np_`, `wf`, `t45`, `ps3`, `qreq`, `lds`.
`docs/notes/symbols.md` is the authoritative symbol table — 106 entries, read from page
images, complete.

Four name collisions, resolved once here so they stay resolved:

| Report symbol | Python name | Why |
|---|---|---|
| NP | `np_` | `np` is NumPy |
| h (convective coefficient) | `h_conv` | collides with the enthalpies `h2`, `h3`, `h41` |
| m (heat-sink metal mass) | `m_metal` | bare `m` reads as a generic mass |
| J (inertia) | `j_gt`, `j_pt`, ... | never bare `j` |

Two typesetting traps the OCR layer cannot show, both load-bearing:

- **Station indices are set full size** (`T41` is T-4-1 inline); **true subscripts are set
  small** (`H41_ns`). Different things.
- **Sixteen constants carry a subscript on a subscript** — `K_H41₁` and `K_H41₂` are the
  slope and intercept of one linear fit and have **different units**. There is no plain
  `K_H41`. Flattening a pair to a single name silently corrupts the fit.

## Architecture constraints

- **`src/t700/` imports NumPy and the stdlib. Nothing else.** No SciPy, Matplotlib,
  pandas, or Numba in the model core. Analysis, digitizing, and plotting live in `tools/`
  and `validation/`, where SciPy and Matplotlib are welcome. A test in `tests/` enforces
  this by walking the imports.
- **All gas property evaluation goes through `t700.thermo`.** No inline `cp`, `gamma`,
  `h`, or T-from-h arithmetic anywhere else, ever. The current backend is pure NumPy; a
  C++ lookup-table backend is planned (see `SCOPE.md` → Deferred), and this interface is
  the reason it can be dropped in later without touching component code.
- **Fixed-step explicit integration only in the core**, matching the report's real-time
  formulation. An adaptive solver would hide the behaviour we are trying to reproduce.
  SciPy adaptive integration may be used in `validation/` to confirm the fixed-step scheme
  is converged — never to run the model.
- **Determinism.** No wall-clock reads, no RNG, no reliance on dict ordering in the core.
  Same input → bit-identical output, run to run and machine to machine.
- Components are plain functions or small classes: explicit state in, derivatives out. No
  global mutable state, no hidden caches.
- No abstraction added "for later" (`advisory`). The thermo backend interface is the one
  exception and it is already specified above.

## Helpers

- `.claude/agents/doc-extract.md` — use this for anything that comes out of the PDF. It
  returns values with page citations. Do not page through 100 pages by hand into context.
- `.claude/agents/physics-review.md` — run before calling a component module done.
- `/validate` — runs the model against digitized reference data and prints a deviation
  table plus overlay plots.
- Tests: `pytest`. Fast unit tests in `tests/`; report-comparison tests in `validation/`.

## Non-negotiables, and what actually enforces them

A rule that lives only in this file is a wish. Each one below names its mechanism. A rule
marked `advisory` has nothing behind it but judgment — those are the ones that slip late
in a long session, so mechanize any that proves to matter.

| Rule | Enforced by |
|---|---|
| Every number cites a report page | `advisory` — the rule that matters most here has the weakest enforcement; challenge any bare constant on sight |
| Never tune a constant to make a validation test pass | `advisory` — a diff touching `data/` in the same breath as a failing test is the tell |
| Never silently widen a tolerance | Tolerances live in `SCOPE.md`, never inline in a test file |
| No SciPy / Matplotlib / pandas / Numba under `src/t700/` | **live** — `tests/test_imports.py` parses every core module's imports |
| Gas properties only via `t700.thermo` | **live** — `tests/test_imports.py` fails any module outside `thermo/` that names a `K_H*`/`K_T*`/`K_TH*` constant |
| Core is deterministic | **live** — `tests/test_imports.py` bans clocks and RNGs structurally; `tests/test_determinism.py` checks it behaviourally |
| No force-push | `permissions.ask`. A plain `git push` was denied outright until 2026-09-11, when the remote was created and publishing was authorized; force-push still prompts |
| Never hand-edit a generated data file | **live** — `bash tools/reproduce_all.sh` reruns every digitizer and diffs; a hand-added note is silently wiped by the next rerun, which happened once to the Figure 9 defect annotations. Put it in the tool. **Run it inside the container**, and read its verdict, not its exit alone: until 2026-09-12 it computed the verdict from `git status data/` only, so run from the host — no NumPy — every digitizer crashed on import, `data/` stayed clean *because nothing had written to it*, and a total failure to execute printed "every data file reproduces byte for byte". It now refuses to run on an interpreter without NumPy and a failed tool poisons the verdict. A gate that cannot fail is not a gate. |
| Python stays formatted and linted | PostToolUse hook: `ruff check --fix` then `ruff format` |

Those tests are written and passing as of 2026-09-10, so three rules that were wishes
are now mechanisms. The two that remain `advisory` are the two that matter most, and
they are the two that cannot be mechanised — a machine cannot tell a cited constant from
a plausible one. Be correspondingly careful.

## Communication

Terse, no preamble. Report deviations with numbers, not adjectives — "NG is 2.4% high
against a 1% tolerance", not "close but slightly off". Say plainly what was not checked.

## Environment

Work inside the `t700-dev` distrobox (Fedora) — the Bazzite host is immutable and has no
NumPy:

```
distrobox enter t700-dev
```

`$HOME` is shared with the host, so the project path is identical inside and out.
