# The environment, and the provenance of the source scan

Two things a replication has as inputs and does not usually write down: the machine that
produced the numbers, and the exact copy of the document they were read from. Both were
unrecorded here until the 2026-09-13 dependency audit
(2026-09-13, findings 7 and 8) pointed out that the only
description of the environment was the sentence "work inside the `t700-dev` distrobox"
in CLAUDE.md — and that container existed on exactly one machine.

## The source document

Every citation in this repository is written `pdf p.NN` **against this specific scan**.
The PDF page numbering is an artifact of the scan's front matter, not of the report
(printed = pdf − 14 across the appendices; see `index.md`), so a different digitization
of TM-100991 would renumber every citation in the codebase silently. The checksum below
is what makes a citation checkable.

| | |
|---|---|
| Path | `docs/ballin-tm100991.pdf` |
| SHA-256 | `55af702919aa809871fcc2035f063d22388ca8d797ae18c7879ac0e4c7cc82e7` |
| Size | 3 028 321 bytes |
| Pages | 104 |
| Page size | 610 × 792 pt |
| PDF producer | `Pdf.Capture version 5.2`, ModDate 2003-05-23 |

Its embedded metadata is empty — no title, author, subject or keywords — so nothing
inside the file identifies it. What identifies it is printed on the NTRS cover sheet,
**read off the raster at 300 dpi, not from the text layer**, which corrupts the
accession number to `NBB-2637B`:

| | |
|---|---|
| Report number | NASA TM-100991 |
| NTRS accession | **N88-26378** |
| Performing organization report | A-88151 [pdf p.103] |
| Work unit | 505-61-51 [pdf p.103] |
| Subject category / microfiche | `G3/08 0156166`, CSCL 21E, Unclas |
| Collation | "(NASA) 95 p" |

**The URL it was downloaded from is not recorded and cannot be recovered.** The checksum
stands in for it: any future copy can be checked against this one, and if it differs,
every `pdf p.NN` in the codebase must be re-checked before it is trusted.

## The environment

Recorded 2026-09-13 from the `t700-dev` container, which is what produced every number
in `README.md` and `SCOPE.md`.

| | |
|---|---|
| OS | Fedora 44 |
| CPython | 3.14.7 (`python3-3.14.7-1.fc44.x86_64`) |
| NumPy | 2.4.6 |
| SciPy | 1.16.2 |
| Matplotlib | 3.10.8 |
| Pillow | 12.3.0 |
| pytest | 8.4.2 |
| ruff | 0.16.6 |
| poppler-utils | 26.01.0 (`poppler-utils-26.01.0-3.fc44.x86_64`) |

`Containerfile` rebuilds it; `requirements-recorded.txt` is the pinned list it installs.
`pyproject.toml` keeps lower bounds instead, because those express what the code needs
rather than what it happened to run on.

### poppler is a prerequisite no Python dependency list can express

Every digitizer under `tools/` shells out to `pdftoppm`, `pdfimages` or `pdfinfo`, so
`bash tools/reproduce_all.sh` needs poppler on PATH. The string "poppler" appeared
nowhere in this repository until 2026-09-13, and the gate would have failed on a fresh
machine with a `FileNotFoundError` naming `pdftoppm` and no indication of what to
install. It is now checked up front by the gate, alongside NumPy.

`src/t700/` needs none of it. The model reads the committed CSVs and nothing else.

### Version drift, measured rather than assumed

The dependency audit found the container (numpy 2.4.6, scipy 1.16.2, pytest 8.4.2) three
minor versions behind a fresh `pip install` of the same `pyproject.toml` (numpy 2.5.3,
scipy 1.18.1, pytest 9.1.1) — every dependency being a bare lower bound. That is a real
reproducibility hole, and it is also the reason a lock file alone would not have closed
it: nothing was *comparing* the two.

Measured on 2026-09-13, against the wheel built from this tree and installed into a clean
venv that resolved numpy 2.5.3:

    trim.solve(wf_pps=400/3600).state
      numpy 2.4.6  ng_rpm=40667.95858675376  p3_psia=163.66043457653322
      numpy 2.5.3  ng_rpm=40667.95858675376  p3_psia=163.66043457653322

Bit-identical on all five states. That is one trim, not a proof; it is the check that
was run, and it is cheap enough to rerun whenever the pins move.

## What this does not cover

No CI runs any of it. The container definition and the pinned list are a record of an
environment, not a gate that fails when the environment drifts from them — and nothing
here notices if `requirements-recorded.txt` and the machine you are on disagree.
