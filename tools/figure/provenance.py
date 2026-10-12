"""The CSV header: one format, written once and read once.

Every digitized file carries its provenance in `#` comment lines -- which page, which
figure, how the numbers were obtained, and the full calibration that produced them. That
header is not decoration. `tools/check_straightness.py` reads the frame model and the axis
maps back out of it and re-derives the figure's coordinates from them, so the header is a
**machine interface**, and it was being written by four different pieces of hand-rolled
f-string formatting and parsed by one set of regexes that had to cope with all of them.

It showed. Two of the eleven engine maps -- `f2` and `f9` -- were excluded from the
straightness audit for no reason except that their hand-written headers omitted the frame
polynomials, while `f7`, written by hand in the same family on the same day, included them.
Since 2026-10-12 they publish them: `calibration_lines` writes the text and the two tools
run `check_calibration` on it.
And the parser carries a scar: an earlier, sloppier test for the axis variable matched the
trailing `s` of the word "samples" in a parenthetical comment, which made every frame line
look like an axis map.

So the formatter and the parser live here, next to each other, where a change to one that
breaks the other is visible in a single file.

## The two conventions, which are both load-bearing

A **frame line** writes its constant FIRST, because it is read as "this edge sits at x =
589.107, sloping by +0.006047 per row":

    #   left   x = 589.107 +0.006047 y   (+0.346 deg, 2295 samples, rms 0.31 px)

An **axis map** writes its constant LAST, in descending powers, because it is a polynomial
in the normalized plot coordinate:

    #   left   y = -0.0067 v^3 +0.0103 v^2 +0.1364 v -0.0200

`parse_poly` classifies terms by whether they carry the variable, never by position, so it
reads both.

This module lives in tools/ and may use SciPy and PIL. Nothing here is imported by
`src/t700/`.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- writing


def poly_str(c, var="t"):
    """A polynomial as printed text, highest power first, constant last."""
    n = len(c) - 1
    out = []
    for i, a in enumerate(c):
        p = n - i
        t = "" if p == 0 else f" {var}" if p == 1 else f" {var}^{p}"
        out.append(f"{a:+.8g}{t}")
    return " ".join(out)


def edge_str(name, coeffs, var, along):
    """One frame edge as printed text: constant first, then the terms in `along`.

    `coeffs` is highest-power-first, as `np.polyfit` returns. Order 1 prints
    `x = 589.107 +0.006047 y`; order 2 adds `-1.23456789e-06 y^2`, which is how the arc
    model publishes itself -- and publishing it is what lets `check_straightness` audit the
    figure at all.

    **Nine significant figures on the non-constant terms, not six decimal places.** A
    quadratic coefficient here is of order 1e-6 and multiplies a pixel coordinate squared,
    so `%.6f` rounds it to one digit and reconstructs the edge two pixels out -- which is
    larger than the bow the term exists to describe.
    """
    c = np.atleast_1d(np.asarray(coeffs, float))
    if c.size > 1 and abs(c[-2]) > 0.1:
        # frames are within a degree of square; a large slope is an (intercept, slope) pair
        raise ValueError(f"{name}: slope {c[-2]:.4g} -- coefficients not highest power first?")
    terms = [f"{c[-1]:.3f}"]
    for p in range(1, len(c)):
        a = c[len(c) - 1 - p]
        t = f" {along}" if p == 1 else f" {along}^{p}"
        terms.append(f"{a:+.9g}{t}")
    return f"{name:6s} {var} = " + " ".join(terms)


def calibration_lines(fc, px, py, xvar="s", yvar="v"):
    """A figure's frame edges and its two blended axis maps, as the header lines
    `parse_calibration` reads back.

    `px` and `py` are `(low-side edge, high-side edge)` polynomial pairs -- bottom then top
    for x, left then right for y -- and that order is printed, because the blend reads the
    first as the side where the cross coordinate is 0.
    """
    lines = [
        edge_str(nm, fc[nm], *(("x", "y") if nm in ("left", "right") else ("y", "x")))
        for nm in ("left", "right", "top", "bottom")
    ]
    for axis, var, sides, polys in (
        ("x", xvar, ("bottom", "top"), px),
        ("y", yvar, ("left", "right"), py),
    ):
        lines += [
            f"{nm:6s} {axis} = {poly_str(c, var)}" for nm, c in zip(sides, polys, strict=True)
        ]
    return lines


# --------------------------------------------------------------------------- reading


def header(path) -> str:
    """The `#` comment block of a CSV, as one string."""
    return "".join(line for line in open(path) if line.startswith("#"))


def parse_poly(expr: str, var: str) -> np.ndarray:
    """Parse a printed polynomial into numpy coefficients, highest power first.

    Copes with both conventions in the module docstring by classifying each term on
    whether it carries the variable, never on where it sits. The parenthetical comment is
    stripped first: it carries words like "samples" whose letters would otherwise be read
    as the variable.
    """
    expr = expr.split("(")[0].strip()
    foreign = set(re.findall(r"(?<![\d.])[A-Za-z_]+", expr)) - {var}
    if foreign:
        # parsing `0.55 u + 0.30` for `y` used to sum the two into a constant, silently
        raise ValueError(f"{expr!r} names {sorted(foreign)}, expected only {var!r}")
    c: dict[int, float] = {}
    pattern = rf"([-+]?\s*[\d.]+(?:[eE][-+]?\d+)?)\s*(?:\*\s*)?({var}(?:\^(\d+))?)?"
    for coef, hasvar, power in re.findall(pattern, expr):
        coef = coef.replace(" ", "")
        if coef in ("", "+", "-"):
            continue
        try:
            value = float(coef)
        except ValueError:
            continue
        if not hasvar:
            k = 0
        elif power:
            k = int(power)
        else:
            k = 1
        c[k] = c.get(k, 0.0) + value
    if not c:
        raise ValueError(f"no polynomial in {expr!r}")
    return np.array([c.get(k, 0.0) for k in range(max(c), -1, -1)])


def parse_calibration(path):
    """Frame model and axis maps, out of the CSV's own header. Nothing is re-derived.

    Returns `(frame, axes, page)` or `None` if the header does not publish a full
    calibration. `axes[a]` is `("blend", [poly_lo, poly_hi])` where the two ticked edges
    were fitted separately, or `("pooled", [poly])` where one map covers the plot.
    """
    return _calibration_of(header(path))


def _calibration_of(txt: str):
    frame = {}
    for name in ("left", "right", "top", "bottom"):
        # the axis maps share this prefix (`left x = ... u`); take the line in the frame's
        # own variable, not the first match
        for m in re.finditer(rf"^#\s+{name}\s+([xy]) = ([^\n(]+)", txt, re.M):
            var = "y" if m.group(1) == "x" else "x"
            if re.search(rf"(?<![A-Za-z]){var}(?![A-Za-z])", m.group(2)):
                frame[name] = parse_poly(m.group(2), var)
                break
        else:
            return None
    axes = {}
    for axis, var in (("x", "u"), ("y", "v")):
        cand = []
        for edge, expr in re.findall(
            rf"^#\s+(left|right|top|bottom)\s+{axis} = ([^\n]+)", txt, re.M
        ):
            body = expr.split("(")[0]
            for sym in (var, "s"):
                if re.search(rf"(?<![A-Za-z]){sym}(?![A-Za-z])", body):
                    cand.append((edge, body, sym))
                    break
        if len(cand) == 2:
            axes[axis] = ("blend", [parse_poly(b, sym) for _, b, sym in cand])
        else:
            m = re.search(rf"^#\s+{axis} = ([^\n]+)", txt, re.M)
            if not m:
                return None
            axes[axis] = ("pooled", [parse_poly(m.group(1), var)])
    page = int(re.search(r"pdf p\.(\d+)", txt).group(1))
    return frame, axes, page


# --------------------------------------------------------------------------- checking


def check_against_csv(csv_path, columns, label=""):
    """Compare computed columns against a **hand-maintained** CSV. 0 if identical, 1 if not.

    Three of the eleven engine maps -- `f2` (Figure A2), `f7` (Figure A7) and `f9` (Figure
    A9) -- have CSVs whose header and rows are written by hand, because the header carries
    a provenance narrative no generator produces. Their tools printed the rows and said
    "compare the rows above", which means a human eye was the only thing between the
    extraction and the committed data.

    **That made `tools/reproduce_all.sh` vacuous for those three maps**, and it was
    demonstrated rather than argued by the 2026-09-13 accuracy audit: corrupt a value in
    `f9_pt_mass_flow.csv`, commit it so `data/` is clean, run the gate, and it prints
    "VERDICT: every data file reproduces byte for byte" and exits 0. `f9` is the map
    Eq. 80's whole repelling-fixed-point story rests on.

    So the tools compare instead of printing, and return non-zero on any difference. The
    CSVs stay hand-maintained; what changes is that the numbers in them are checked against
    the extraction on every gate run.

    `columns` maps CSV column name -> (values, format), where the format is the one that
    reproduces that column's printed precision. The comparison is on the formatted text, so
    it is exactly the check a writer would perform.
    """
    path = Path(csv_path)
    lines = [ln.rstrip("\n") for ln in path.read_text().splitlines()]
    body = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
    if not body:
        print(f"FAIL {label}: {path} has no data rows")
        return 1
    names = [n.strip() for n in body[0].split(",")]
    rows = [[v.strip() for v in ln.split(",")] for ln in body[1:]]

    missing = [n for n in columns if n not in names]
    if missing:
        print(f"FAIL {label}: {path} has no column(s) {missing}; header is {names}")
        return 1

    n_expected = len(next(iter(columns.values()))[0])
    if len(rows) != n_expected:
        print(f"FAIL {label}: {path} has {len(rows)} rows, the tool computed {n_expected}")
        return 1

    bad = 0
    for i, row in enumerate(rows):
        for name, (values, fmt) in columns.items():
            got = fmt % values[i]
            want = row[names.index(name)]
            if got != want:
                bad += 1
                print(f"  row {i + 1} column {name}: csv has {want}, the tool computes {got}")
    if bad:
        print(
            f"FAIL {label}: {bad} value(s) in {path} differ from what the tool extracts "
            f"today. The CSV is maintained by hand, so either the hand edit was wrong or "
            f"the tool has changed -- decide which, and say so in the commit."
        )
        return 1
    print(f"ok  {label}: all {len(rows)} rows of {path} match the extraction")
    return 0


def check_calibration(csv_path, lines, label=""):
    """The calibration a **hand-maintained** header publishes, against the tool's own.

    `lines` is `calibration_lines(...)` for what the tool computed today. What is compared
    is what `parse_calibration` reads out of the header -- the thing `check_straightness`
    uses -- so a stale line, a missing one, or a second line for the same edge that the
    parser would take first all fail. 0 if identical, 1 if not.
    """
    got = parse_calibration(csv_path)
    want = _calibration_of("".join(f"#   {ln}\n" for ln in lines) + "# pdf p.0\n")
    same = (
        got is not None
        and all(np.array_equal(got[0][nm], want[0][nm]) for nm in want[0])
        and all(
            got[1][a][0] == want[1][a][0]
            and len(got[1][a][1]) == len(want[1][a][1])
            and all(np.array_equal(p, q) for p, q in zip(got[1][a][1], want[1][a][1], strict=True))
            for a in ("x", "y")
        )
    )
    if not same:
        print(
            f"FAIL {label}: the calibration {csv_path}'s header publishes is not the one the "
            f"tool computes today. The header is maintained by hand; the tool's lines are:"
        )
        for ln in lines:
            print(f"#   {ln}")
        return 1
    print(f"ok  {label}: the header publishes the tool's frame arcs and axis maps")
    return 0
