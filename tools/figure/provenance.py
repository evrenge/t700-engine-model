"""The CSV header: one format, written and read in one place.

Every digitized file carries its provenance and full calibration in `#` lines.
`tools/check_straightness.py` reads the frame model and axis maps back out of it, so the
header is a machine interface; the formatter and the parser live together here.

Two conventions, both load-bearing:

    #   left   x = 589.107 +0.006047 y   (...)      frame edge: constant FIRST
    #   left   y = -0.0067 v^3 ... -0.0200          axis map: constant LAST, descending

`parse_poly` classifies terms by whether they carry the variable, so it reads both.
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

    `coeffs` is highest power first, as `np.polyfit` returns. Nine significant figures on
    the non-constant terms, not `%.6f`: a quadratic coefficient is ~1e-6, and `%.6f` rounds
    it to one digit and reconstructs the edge two pixels out.
    """
    c = np.atleast_1d(np.asarray(coeffs, float))
    if c.size > 1 and abs(c[-2]) > 0.1:
        # frames are within a degree of square; a large slope is an (intercept, slope) pair
        # (recipe, trap 12)
        raise ValueError(f"{name}: slope {c[-2]:.4g} -- coefficients not highest power first?")
    terms = [f"{c[-1]:.3f}"]
    for p in range(1, len(c)):
        a = c[len(c) - 1 - p]
        t = f" {along}" if p == 1 else f" {along}^{p}"
        terms.append(f"{a:+.9g}{t}")
    return f"{name:6s} {var} = " + " ".join(terms)


def calibration_lines(fc, px, py, xvar="s", yvar="v"):
    """The frame edges and the two blended axis maps, as header lines `parse_calibration` reads.

    `px` and `py` are `(low-side edge, high-side edge)` polynomial pairs -- bottom then top
    for x, left then right for y -- printed in that order, because the blend reads the first
    as the side where the cross coordinate is 0.
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

    Terms are classified by whether they carry the variable, not by position. The
    parenthetical comment is stripped first: its words ("samples") would read as the variable.
    """
    expr = expr.split("(")[0].strip()
    foreign = set(re.findall(r"(?<![\d.])[A-Za-z_]+", expr)) - {var}
    if foreign:
        # else `0.55 u + 0.30` parsed for `y` sums to a constant, silently
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
    """Compare computed columns with a hand-maintained CSV (f2, f7, f9). 0 if equal, 1 if not.

    `columns` maps column name -> (values, printf format). The comparison is on the formatted
    text, which is exactly the check a writer would perform.
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
    """The calibration a hand-maintained header publishes, against the tool's own.

    `lines` is `calibration_lines(...)` as computed today. The header is compared as
    `parse_calibration` reads it, so a stale line, a missing one, or a duplicate the parser
    would read first all fail. 0 if identical, 1 if not.
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
