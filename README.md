# t700-engine-model

A Python replication of the real-time turboshaft engine model in:

> Ballin, M. G., *A High Fidelity Real-Time Simulation of a Small Turboshaft Engine*,
> NASA TM-100991, July 1988.

The subject is the **General Electric T700-GE-700**, the engine on the UH-60A Black Hawk.
Ballin built a component-level thermodynamic model that had to compute each frame inside a
**10 ms** step [pdf p.38] while a pilot flew a helicopter simulation around it; the report is
largely the record of how he made an implicit engine model run in real time, and what that
cost. The report is included as `docs/ballin-tm100991.pdf` (NTRS accession N88-26378), and
everything here is derived from it.

## Results

The numbers below come from `results/report.json`, which `validation/report.py` writes by
re-running the model, together with the comparison figures in `results/`;
`tests/test_live_numbers.py` fails if a number is no longer current.

| | |
|---|---|
| **Table B.1, the printed engine state** | **rms 0.22 %<!--rj:table_b1.rms-->** over 21<!--rj:table_b1.n--> numbers; worst shaft power −0.81 %<!--rj:table_b1.worst--> at descent |
| **Table 1, printed eigenvalues** | all **27**<!--rj:eigenvalues.notes.table_1_stats.n--> modes across three model variants: **rms 1.9 %<!--rj:eigenvalues.notes.table_1_stats.rms-->**, worst 4.7 %<!--rj:eigenvalues.notes.table_1_stats.worst--> (ours more damped) on the 2-DOF NG mode at hover, modes paired by dominant state, with the report's own ±2 % extraction step [pdf p.28] |
| **Appendix B, 297 printed elements** | zero structure exact; all 206 non-zero elements compared one by one. `A` rms 6.3 %<!--rj:appendix_b.A.rms--> over 158<!--rj:appendix_b.A.n-->, `b` rms 6.4 %<!--rj:appendix_b.b.rms--> over 42<!--rj:appendix_b.b.n-->; the heat-sink rows carry the report's own TC_T41 inconsistency (open question #31) |
| **Closed loop vs Table B.1** | **worst 0.25 %<!--rj:closed_loop.worst-->** over NG, NP, Wf, Ps3 and shaft power at three trims, **with fuel flow as an output**, averaged over the governor's own limit cycle |
| Figure 6, NG vs fuel flow | mean -0.01 %NG<!--rj:steady.figure_6.mean-->, rms 0.23<!--rj:steady.figure_6.rms-->, worst −1.02<!--rj:steady.figure_6.worst--> over 29<!--rj:steady.figure_6.n--> of 29<!--rj:steady.figure_6.printed_points--> printed points |
| Figure 7, shaft power vs fuel flow | mean +0.18 %<!--rj:steady.figure_7.mean-->, rms 0.66 %<!--rj:steady.figure_7.rms-->, worst +2.17<!--rj:steady.figure_7.worst--> over 15<!--rj:steady.figure_7.n--> of 15<!--rj:steady.figure_7.printed_points--> |
| Figure 8, Ps3 vs NG | mean +0.07 %<!--rj:steady.figure_8.mean-->, rms 0.26 %<!--rj:steady.figure_8.rms-->, worst +0.64<!--rj:steady.figure_8.worst--> over 27<!--rj:steady.figure_8.n--> of 27<!--rj:steady.figure_8.printed_points--> |
| Figures 9 and 10, fuel transients | whole-curve rms **0.62<!--rj:whole_curve_summary.min-->–3.80 %<!--rj:whole_curve_summary.max-->** of each panel's excursion, mean 1.64 %<!--rj:whole_curve_summary.mean--> over ten panels |
| The report's own claim about Figure 6 | reproduced: in the 81–86 %NG band the report names, we overestimate fuel by **+6.34 %<!--rj:fuel_overestimate.ours_peak_pct-->**, against its own line's +5.98 %<!--rj:fuel_overestimate.ballin_peak_pct--> and its printed "as much as five percent" [pdf p.39] |
| Figures 4 and 5, the linear models in time | at hover, Ballin's curves against ours, rms in % of trim: his 5-DOF against our 5-DOF 0.028<!--rj:figures_4_5.pairs.figure 5: 5-DOF vs our 5-DOF.np--> in NP, his 6-DOF (heat sink) against ours 0.033<!--rj:figures_4_5.pairs.figure 5: 6-DOF vs our 6-DOF.np-->, ours slightly more damped; Figure 4's linear models against our 5-DOF 0.025<!--rj:figures_4_5.pairs.figure 4: linear vs our 5-DOF.np-->, and his Gen Hel simulation against our nonlinear engine 0.045<!--rj:figures_4_5.pairs.figure 4: Gen Hel vs our nonlinear.np-->, whose rotor-blade lag dynamics the linear rotor load leaves out [pdf p.29] |
| Figures 11–14, Gen Hel's closed loop | his printed rotor speed, and his collective where it moves, replayed into our engine and fuel control; rms against his Gen Hel line in NG (%NG), fuel flow and torque (% of trim): collective at 90 kt 0.18<!--rj:gen_hel.figure 11.rms vs his gen hel.ng-->, 1.9<!--rj:gen_hel.figure 11.rms vs his gen hel.fuel-->, 2.2<!--rj:gen_hel.figure 11.rms vs his gen hel.torque-->; lateral cyclic at 55 kt 0.030<!--rj:gen_hel.figure 12.rms vs his gen hel.ng-->, 0.40<!--rj:gen_hel.figure 12.rms vs his gen hel.fuel-->, 0.37<!--rj:gen_hel.figure 12.rms vs his gen hel.torque-->; collective at hover 0.38<!--rj:gen_hel.figure 13.rms vs his gen hel.ng-->, 2.0<!--rj:gen_hel.figure 13.rms vs his gen hel.fuel-->, 3.3<!--rj:gen_hel.figure 13.rms vs his gen hel.torque-->; pedals at hover 0.11<!--rj:gen_hel.figure 14.rms vs his gen hel.ng-->, 1.4<!--rj:gen_hel.figure 14.rms vs his gen hel.fuel-->, 1.6<!--rj:gen_hel.figure 14.rms vs his gen hel.torque--> (#81) |
| Figures 11–15, the second test | his rotor's demand rebuilt from his torque and rotor speed (a rigid drive train, pdf p.47) is the load, and our model makes the rotor speed: rms against his Gen Hel rotor speed 0.155<!--rj:gen_hel_load.figure 11.rms vs his gen hel.nr--> %, 0.029<!--rj:gen_hel_load.figure 12.rms vs his gen hel.nr--> %, 0.158<!--rj:gen_hel_load.figure 13.rms vs his gen hel.nr--> %, 0.073<!--rj:gen_hel_load.figure 14.rms vs his gen hel.nr--> % on Figures 11–14, and 0.89<!--rj:gen_hel_15.rms vs his gen hel.nr--> % on Figure 15, a large collective input, from its main-rotor torque (#81) |
| Table 3, the NASA-Lewis test engine | a characterization, not a target (#34): on the standard maps, at Table 2's conditions, ours passes +12.0 %<!--rj:table_3.stats_cases_2_6.WA2.ours_vs_ballin_pct.mean--> more air and runs a −7.5 %<!--rj:table_3.stats_cases_2_6.T45.ours_vs_ballin_pct.mean--> cooler T45 than Ballin's model (means, cases 2-6), which ran the test engine's own maps [pdf p.39] |
| Frame cost | tens of microseconds per frame, against the report's 10 ms budget |

## What is in here

```
src/t700/        the model: NumPy and the standard library, nothing else
  engine.py      the gas path, Eqs. 1-49, in the report's own order
  realtime.py    the real-time frame, Eqs. 69-80 -- the model Ballin ran
  trim.py, linear.py   equilibrium solver; the Jacobian, by the report's own method
  maps.py        the eleven function tables, loaded and interpolated
  constants.py   Table A.1, every value cited to a page
  control/       the fuel control: ECU, HMU, the closed loop, 73 constants, 8 schedules
  thermo/        gas properties behind one interface
data/            11 engine maps, 8 control schedules, the 297 Appendix B elements, and
                 the transient and sweep traces to validate against
tools/           the digitizers; tools/reproduce_all.sh reruns them and diffs every CSV
tools/figure/    the one library they share: each stage of reading a printed figure, once
tests/           fast unit tests, and the architecture rules as executable checks
validation/      comparisons against the report, and the figure set
docs/notes/      the report transcribed (equations, symbols, appendix inventories), and
                 the open-questions ledger
results/         report.json and the comparison figures, written by validation/report.py
```

## The rules it was built under

`CLAUDE.md` carries the full set and names what enforces each. Two matter most, and no
machine can check either:

> **Every number carries a citation to a page of the report.** Never a value from general
> turbomachinery knowledge, from another engine model, or from a plausible-looking guess.
> If the report does not give it, that goes in the open-questions ledger and the work stops
> until it is decided.

**Never tune a constant to make a comparison pass.** Shaft power at the hover and level trims
sits at −0.060 %<!--rj:table_b1.rows.6.dev_pct--> and −0.077 %<!--rj:table_b1.rows.13.dev_pct-->
from Table B.1 with no data value adjusted to get there: where a comparison was off, the
fault was in our own pixels-to-numbers map, not in the report's data.

## Status

Phases 0–6 of `SCOPE.md` are complete. The engine trims and runs transients in **either of
the two configurations Ballin published**, with or without the station 4.1 heat-sink model,
and the Appendix C fuel control is closed around it, so fuel flow is an **output** of the
model rather than an input to it.

`SCOPE.md` sets out what is not claimed. In short: Figures 11–15 need the Gen Hel UH-60A
blade-element simulation, which the report consumes and does not contain, so Figures 11–15
are compared by replaying what it prints: rotor speed, collective and the rotor's demand; Table 3 was
produced with functions from a different engine and is evidence of method, never pass/fail;
and nothing above 100 %NG is claimed, because Ballin claims nothing there either [pdf p.39].

As of 2026-10-14, `docs/notes/open-questions.md` has logged 78 questions: 75 closed, 3 open
(#60, #74, #82). Each row states what was asked and what was decided.

## Running it

The model needs NumPy. Analysis and digitizing also want SciPy, Matplotlib and Pillow, and
**poppler** (`pdftoppm`, `pdfimages`), which is not a Python package; the digitizers shell
out to it to raster the scan.

```bash
pip install -e '.[tools,dev]'   # [dev] pulls in [tools]: validation/ imports SciPy
pytest
PYTHONPATH=src python validation/report.py   # every comparison, as figures and report.json
```

A plain `pip install .` also works; the wheel carries `data/`. `Containerfile` builds the
environment the published numbers were produced under, and `docs/notes/environment.md`
records its versions.

## Licence and provenance

This implementation is MIT; see `LICENSE`. The report is a NASA technical memorandum, not
subject to US copyright under 17 U.S.C. §105. `LICENSE` writes out two qualifications: §105
is US-only, and the report reproduces at least one figure it does not own. Figure 2 is
captioned "(From ref. 6.)" [pdf p.18], a General Electric training guide; no number here
comes from it. The scan is pinned by checksum in `docs/notes/environment.md`: every
citation here is `pdf p.NN` against *this* digitization, and another copy may number its
pages differently.
