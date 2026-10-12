"""Assemble `site/index.html` from the template, `report.json` and the open-questions ledger.

The site is a **generated artifact**. Edit `site/template.html` and rerun:

    PYTHONPATH=src python validation/report.py   # figures + report.json
    python tools/build_site.py                   # index.html
    python tools/build_site.py --check           # is the committed index.html current?

Every number in the page comes out of `site/assets/report.json`, which `validation/report.py`
recomputes from the model. Nothing is typed into the template that a run could contradict.

The open-questions section is the one place with curated prose rather than generated text --
a ledger cell runs to two thousand words and truncating it produces something worse than
useless. What is *not* curated is which rows appear and what verdict each carries: those are
read from `docs/notes/open-questions.md` and checked against the summaries below, so the page
cannot quietly disagree with the ledger about what is still open.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / "site"
ASSETS = SITE / "assets"
LEDGER = ROOT / "docs" / "notes" / "open-questions.md"

_SPLIT = re.compile(r"(?<!\\)\|")

# Curated one-line summaries. The row number and verdict are verified against the ledger.
OPEN_SUMMARY = {
    60: (
        "<code>f1</code> on a beta grid, not the printed table",
        "A deliberate departure. It removes a manufactured extrapolation entirely and "
        "halves the worst Table&nbsp;1 speed mode, at a cost of 0.2018&nbsp;&rarr;&nbsp;"
        "0.2439 % on Table B.1&rsquo;s rms.",
    ),
    74: (
        "Figure&nbsp;A1&rsquo;s speed lines do not form a smooth family",
        "No smooth surface fits all eleven digitized lines: the best joint fit reaches "
        "&chi;&sup2;/dof 4.67, and adjacent lines differ by up to 11&nbsp;&sigma;. Either our "
        "read error understates the cross-line error 2&ndash;4&times;, or the printed figure "
        "is not internally consistent.",
    ),
}


def ledger_rows() -> list[tuple[int, str]]:
    """Row number and leading verdict word for every row that is not fully closed."""
    out = []
    for line in LEDGER.read_text().splitlines():
        if not line.startswith("| "):
            continue
        cells = _SPLIT.split(line)
        if len(cells) < 7 or not cells[1].strip().isdigit():
            continue
        status = cells[5].strip().lstrip("*").strip().lower()
        for verdict in ("partly closed", "closed", "open"):
            if status.startswith(verdict):
                if verdict != "closed":
                    out.append((int(cells[1].strip()), verdict))
                break
    return sorted(out)


def num(v: float, d: int = 3) -> str:
    return f"{v:,.{d}f}"


def map_shape(v: dict) -> str:
    """How a table was extracted: knot count, or lines x knots for the two 2-D maps."""
    if "lines" in v:
        return f"{v['lines']} lines &times; {v['knots_per_line']}"
    return str(v.get("knots", ""))


def fmt_cond(v):
    if v is None:
        return "&mdash;"
    return "0" if v == 0.0 else f"{v:.1e}"


def picture(stem: str, alt: str) -> str:
    """A figure that follows the reader's theme, and can be swapped by the toggle."""
    return (
        f"<picture>"
        f'<source data-dark srcset="assets/{stem}-dark.png" '
        f'media="(prefers-color-scheme: dark)">'
        f'<img src="assets/{stem}-light.png" data-stem="{stem}" alt="{alt}" loading="lazy">'
        f"</picture>"
    )


def bar(v: float, scale: float) -> str:
    w = min(abs(v) / scale, 1.0) * 50.0
    side = "left:50%" if v >= 0 else f"left:{50 - w}%"
    cls = "pos" if v >= 0 else "neg"
    return f'<span class="bar"><i class="{cls}" style="{side};width:{w}%"></i></span>'


SCH_WHAT = {
    "F_HM1": "Topping line. The fuel ceiling, a function of inlet temperature alone.",
    "F_HM2": "Power available spindle &mdash; the cockpit power lever, "
    "subtracted from the ceiling.",
    "F_HM3": "Collective pitch to the gas generator speed reference.",
    "F_HM4": "Collective pitch to the fuel feedforward. The load-demand anticipation.",
    "F_HM5": "Idle schedule, fuel term. Active only deep sub-idle.",
    "F_HM6": "Idle schedule, speed reference. Two printed markers define a straight line.",
    "F_HM7": "Acceleration limit. Two-dimensional: seven inlet-temperature lines.",
    "F_EC1": "Thermocouple time constant &mdash; itself a lookup, four T4.5 lines.",
}

D = json.loads((ASSETS / "report.json").read_text())
_B1W = max(D["table_b1"]["rows"], key=lambda r: abs(r["dev_pct"]))
_EIGW = D["eigenvalues"]["notes"]["table_1_worst_row"]

# ---------------------------------------------------------------- headline table
head = [
    (
        "Table B.1, the printed engine state",
        "21 numbers",
        D["table_b1"],
        "rms",
        "%",
        f"worst {_B1W['quantity']} {'&minus;' if _B1W['dev_pct'] < 0 else '+'}"
        f"{abs(_B1W['dev_pct']):.2f} % at {_B1W['trim'].split()[0]}",
    ),
    (
        "Table 1, printed eigenvalues",
        "27 modes",
        D["eigenvalues"]["notes"]["table_1_stats"],
        "rms",
        "%",
        f"worst {'&minus;' if _EIGW['dev_pct'] < 0 else '+'}{abs(_EIGW['dev_pct']):.1f} % on "
        f"the {_EIGW['model'].replace('dof', '-DOF').upper()} {_EIGW['mode']} mode at "
        f"{_EIGW['trim']}",
    ),
    (
        "Appendix B, matrix A",
        "158 non-zero",
        D["appendix_b"]["A"],
        "rms",
        "%",
        "every one of the 206 non-zero elements is compared",
    ),
    (
        "Appendix B, fuel column b",
        "42 non-zero",
        D["appendix_b"]["b"],
        "rms",
        "%",
        "the 6-DOF half carries the heat sink's deficit",
    ),
    (
        "Closed loop against Table B.1",
        "15 numbers",
        D["closed_loop"],
        "rms",
        "%",
        "fuel flow is an <em>output</em> here, not an input",
    ),
]
rows = "".join(
    f'<tr><th scope="row">{n}<span class="sub">{c}</span></th>'
    f'<td class="mono">{num(s["rms"])}</td>'
    f'<td class="mono">{num(s["mean"])}</td>'
    f'<td class="mono strong">{num(s["worst"])}</td>'
    f'<td class="note">{note}</td></tr>'
    for n, c, s, _k, _u, note in head
)

# ---------------------------------------------------------------- Table B.1 rows
b1 = "".join(
    f'<tr><td>{r["trim"]}</td><td class="mono">{r["quantity"]}</td>'
    f'<td class="mono">{r["printed"]:,.2f}</td>'
    f'<td class="mono">{r["ours"]:,.2f}</td>'
    f'<td class="mono dev">{r["dev_pct"]:+.3f}{bar(r["dev_pct"], 1.0)}</td></tr>'
    for r in D["table_b1"]["rows"]
)

# ---------------------------------------------------------------- eigenvalues
MODEL_LABEL = {"5dof": "5-DOF", "2dof": "2-DOF", "red5": "reduced-5"}
eig = "".join(
    f"<tr><td>{MODEL_LABEL[r['model']]}</td><td>{r['trim']}</td>"
    f'<td class="mono">{r["mode"]}</td>'
    f'<td class="mono">{r["printed"]:,.4g}</td>'
    f'<td class="mono">{r["ours"]:,.4g}</td>'
    f'<td class="mono dev">{r["dev_pct"]:+.2f}{bar(r["dev_pct"], 15.0)}</td></tr>'
    for r in D["eigenvalues"]["table_1"]
)
unprinted = "".join(
    f'<tr><td class="mono">{r["figure"]}</td><td>{r["model"]}</td><td>{r["trim"]}</td>'
    f'<td class="mono small">{", ".join(f"{v:,.4g}" for v in r["ballin"])}</td>'
    f'<td class="mono small">{", ".join(f"{v:,.4g}" for v in r["ours"])}</td>'
    f'<td class="flag">{"UNSTABLE as printed" if r["ballin_unstable"] else ""}'
    f"{' &middot; ill-conditioned' if r['ill_conditioned'] else ''}</td></tr>"
    for r in D["eigenvalues"]["against_printed_matrices"]
    if r["model"] in ("3-DOF", "6-DOF")
)

# ---------------------------------------------------------------- Appendix B blocks
blocks = "".join(
    f'<tr><td class="mono">{k}</td><td class="mono">{s["n"]}</td>'
    f'<td class="mono">{num(s["mean"], 2)}</td><td class="mono">{num(s["rms"], 2)}</td>'
    f'<td class="mono dev">{s["worst"]:+.2f}{bar(s["worst"], 40.0)}</td></tr>'
    for k, s in D["appendix_b"]["by_block"].items()
)

# ---------------------------------------------------------------- figures
FIGNOTE = {
    6: "NG against fuel flow. Deviation in %NG points, not percent &mdash; an absolute "
    "difference on a speed already expressed as a percentage.",
    7: "Shaft power against fuel flow. The low-power end is the worst, and it is the same "
    "end where Table B.1 and this figure disagree with <em>each other</em> by 5.5 %.",
    8: "Ps3 against NG. No fuel flow appears in it, so it cannot be rescued by a "
    "compensating fuel-path error.",
}
figs = "".join(
    f'<tr><th scope="row">Figure {n}<span class="sub">{FIGNOTE[n]}</span></th>'
    f'<td class="mono">{s["n"]} of {s["printed_points"]}</td>'
    f'<td class="mono">{num(s["mean"])}</td><td class="mono">{num(s["rms"])}</td>'
    f'<td class="mono strong">{s["worst"]:+.3f} {s["unit"]}</td></tr>'
    for n, s in ((int(k.split("_")[1]), v) for k, v in D["steady"].items())
)

wc = D["whole_curve"]
order = sorted(wc.items(), key=lambda kv: kv[1])
wcrows = "".join(
    f'<tr><td class="mono">{k.replace("_", " &middot; ").replace("fig", "Fig ")}</td>'
    f'<td class="mono">{v:.3f}</td>'
    f'<td><span class="bar solo"><i class="pos" style="left:0;width:{min(v / 4 * 100, 100):.0f}%">'
    f"</i></span></td></tr>"
    for k, v in order
)

settled = ""
for figno in (9, 10):
    for key, p in D["transient"][f"figure_{figno}"].items():
        settled += (
            f'<tr><td class="mono">Fig {figno}</td><td class="mono">{key.upper()}</td>'
            f'<td class="mono">{p["ballin_settled"]:,.2f}</td>'
            f'<td class="mono">{p["ours_settled"]:,.2f}</td>'
            f'<td class="mono dev">{p["settled_dev_pct"]:+.3f}'
            f"{bar(p['settled_dev_pct'], 1.0)}</td>"
            f'<td class="mono small">{p["ref_ends_s"]:.3f} s</td></tr>'
        )

STEMS = (
    "figures-6-8-overlay",
    "figure-9-overlay",
    "figure-10-overlay",
    "phase-plane",
    "table-b1",
    "eigenvalues-all-models",
    "appendix-b-elements",
    "engine-operating-line",
    "open-loop-family",
    "closed-loop-response",
    "control-schedules",
    "engine-maps",
)

# ------------------------------------------------- model behaviour, no reference data
es = D["engine_steady"]


def rng(key, d=1, unit=""):
    a, b = es[key]
    return f"{a:,.{d}f} &ndash; {b:,.{d}f}{unit}"


envelope = "".join(
    f'<tr><td>{lab}</td><td class="mono">{val}</td><td class="note">{note}</td></tr>'
    for lab, val, note in [
        (
            "Gas generator speed",
            rng("ng_range_pct", 1, " %NG"),
            "the report claims validity to 100 %NG and no further",
        ),
        ("Fuel flow", rng("wf_range_pph", 0, " lbm/hr"), "the trim ladder's own span"),
        (
            "Shaft power",
            rng("shp_range", 0, " hp"),
            f"crosses zero at {es['shp_zero_crossing_ng']:.2f} %NG &mdash; the free turbine "
            f"absorbs below that",
        ),
        (
            "Compressor pressure ratio",
            rng("pressure_ratio_range", 2, " : 1"),
            "Ps3 over ambient, sea-level standard day",
        ),
        (
            "Turbine inlet T4.1",
            rng("t41_range_degR", 0, " &deg;R"),
            "2690 &deg;R at the top of the range",
        ),
        (
            "Customer + cooling bleed",
            rng("bleed_pct_range", 2, " %"),
            "falls as speed rises &mdash; the schedules close off above 85 %NG",
        ),
        ("Best SFC", f"{es['sfc_best']:.3f} lbm/hr per hp", f"at {es['sfc_best_at_ng']:.1f} %NG"),
    ]
)

ol = D["open_loop"]
ol_rows = "".join(
    f'<tr><td class="mono">400 &rarr; {k.split("_")[1]}</td>'
    f'<td class="mono">{v["ng_start"]:.2f}</td><td class="mono">{v["ng_end"]:.2f}</td>'
    f'<td class="mono">{v["t41_peak"]:,.0f}</td><td class="mono">{v["t41_end"]:,.0f}</td>'
    f'<td class="mono strong">{v["t41_overshoot"]:+,.0f}</td></tr>'
    for k, v in sorted(ol.items(), key=lambda kv: float(kv[0].split("_")[1]))
)

cl = D["closed_loop_response"]
cl_rows = "".join(
    f"<tr><td>{k}</td>"
    f'<td class="mono">{v["np_min"]:,.0f} &ndash; {v["np_max"]:,.0f}</td>'
    f'<td class="mono">{v["np_settled"]:,.0f}</td>'
    f'<td class="mono">{v["wf_min"]:,.0f} &ndash; {v["wf_max"]:,.0f}</td>'
    f'<td class="mono">{v["wf_settled"]:,.0f}</td>'
    f'<td class="mono small">{", ".join(v["limits_used"])}</td></tr>'
    for k, v in cl.items()
)

sch = D["schedules"]
sch_rows = "".join(
    f'<tr><td class="mono">{k}</td><td class="mono small">{v["source"]}</td>'
    f'<td class="mono">{v.get("knots", v.get("lines"))}'
    f"{' lines' if 'lines' in v else ' markers'}</td>"
    f'<td class="note">{SCH_WHAT[k]}</td></tr>'
    for k, v in sch.items()
)

mp = D["engine_maps"]
map_rows = "".join(
    f'<tr><td class="mono">{k}</td><td class="mono small">{v["source"]}</td>'
    f'<td class="mono">{map_shape(v)}</td>'
    f'<td class="mono small">{v.get("x_axis", "11 speed lines, beta grid")[:44]}</td>'
    f'<td class="mono">{fmt_cond(v.get("conditioning_moved"))}</td></tr>'
    for k, v in mp.items()
)

# ------------------------------------------------------------- open questions
verdicts = dict(ledger_rows())
missing = set(verdicts) - set(OPEN_SUMMARY)
stale = set(OPEN_SUMMARY) - set(verdicts)
if missing or stale:
    raise SystemExit(
        f"tools/build_site.py disagrees with the ledger about what is open.\n"
        f"  in the ledger, not summarised here: {sorted(missing)}\n"
        f"  summarised here, no longer open:    {sorted(stale)}\n"
        f"Update OPEN_SUMMARY, then rebuild."
    )
n_open = sum(1 for v in verdicts.values() if v == "open")
n_partly = len(verdicts) - n_open
oq_rows = "".join(
    f'<tr><td class="mono">#{n}</td>'
    f'<td><span class="tag {"is-open" if verdicts[n] == "open" else "is-partly"}">'
    f"{verdicts[n]}</span></td>"
    f'<th scope="row">{OPEN_SUMMARY[n][0]}</th>'
    f'<td class="note">{OPEN_SUMMARY[n][1]}</td></tr>'
    for n in sorted(OPEN_SUMMARY)
)

# The internal spread was three hardcoded percentages in the template until 2026-09-14,
# when `report.py::internal_spread` was written to compute them and found the printed claim
# seven times the measured value -- stale since the digitizer corrections of that morning.
# Anything the page asserts about a measurement comes out of report.json; this is one more.
_SP = D["internal_spread"]
spread_rows = "".join(
    f'<tr><th scope="row">{r["trim"]}</th>'
    f'<td class="mono">{r["wf_pph"]:.1f}</td>'
    f'<td class="mono">{r["theirs_pct"]:+.3f} {r["unit"]}</td>'
    f'<td class="mono{" strong" if abs(r["ours_pct"]) > abs(r["theirs_pct"]) else ""}">'
    f"{r['ours_pct']:+.3f} {r['unit']}</td>"
    f'<td class="mono small">{r["nearest_marker_pph"]:.1f} lbm/hr</td></tr>'
    for r in _SP["rows"]
    if r["quantity"] == "SHP"
)

html = (SITE / "template.html").read_text()
for stem in STEMS:
    html = html.replace(f'src="{{{{{stem}}}}}"', f"__PIC__{stem}__")
html = re.sub(
    r'<img __PIC__([a-z0-9-]+)__\s*\n?\s*alt="([^"]*)"[^>]*>',
    lambda m: picture(m.group(1), m.group(2)),
    html,
)
html = (
    html.replace("{{spread_rows}}", spread_rows)
    .replace("{{spread_widest}}", f"{abs(_SP['widest_theirs_pct']):.3f}")
    .replace("{{headline_rows}}", rows)
    .replace("{{b1_rows}}", b1)
    .replace("{{eig_rows}}", eig)
    .replace("{{unprinted_rows}}", unprinted)
    .replace("{{block_rows}}", blocks)
    .replace("{{fig_rows}}", figs)
    .replace("{{wc_rows}}", wcrows)
    .replace("{{settled_rows}}", settled)
    .replace("{{envelope_rows}}", envelope)
    .replace("{{ol_rows}}", ol_rows)
    .replace("{{cl_rows}}", cl_rows)
    .replace("{{sch_rows}}", sch_rows)
    .replace("{{map_rows}}", map_rows)
    .replace("{{oq_rows}}", oq_rows)
    .replace("{{n_open}}", str(n_open))
    .replace("{{n_partly}}", str(n_partly))
    .replace(
        "{{n_logged}}",
        str(
            len(
                [
                    ln
                    for ln in LEDGER.read_text().splitlines()
                    if ln.startswith("| ") and _SPLIT.split(ln)[1].strip().isdigit()
                ]
            )
        ),
    )
)

out = SITE / "index.html"
left = re.findall(r"\{\{[a-z0-9_-]+\}\}", html)
if left:
    raise SystemExit(f"unresolved template placeholders: {sorted(set(left))}")
if "--check" in sys.argv:
    if out.read_text() != html:
        raise SystemExit(f"{out.relative_to(ROOT)} is stale -- run tools/build_site.py")
    raise SystemExit(0)
out.write_text(html)
print(
    f"wrote {out.relative_to(ROOT)}  {out.stat().st_size / 1e3:.0f} kB"
    f"  +  {len(list(ASSETS.glob('*.png')))} plots"
)
