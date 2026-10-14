"""Every figure draws each data source the same way: Ballin orange, flight test grey dashed, GE
green, ours blue, as `validation/plotstyle.py` sets them.

The figures are drawn by `validation/report.py` alone, so the rule is checked there,
statically. A colour reaches a plot call only through `plotstyle.line`, `plotstyle.points`,
`plotstyle.guide`, `plotstyle.colour` or `plotstyle.shades`, or as one of the theme's neutral
inks for text, frames and guides. A hex literal, a colormap or a colour letter in a format
string is refused.

Missed: a role used for the wrong source, such as `plotstyle.line("ballin")` on a line of ours.
The checker cannot tell whose data a line carries; the label beside it must be read.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPORT = ROOT / "validation" / "report.py"

COLOUR_KEYWORDS = {
    "color",
    "c",
    "colors",
    "mfc",
    "mec",
    "markerfacecolor",
    "markeredgecolor",
    "facecolor",
    "edgecolor",
    "fc",
    "ec",
    "labelcolor",
}
NEUTRAL = {
    "MUTED",
    "THEME.ink",
    "THEME.ground",
    "THEME.muted",
    "THEME.faint",
    "THEME.grid",
    "'none'",
}
"""Inks that carry no data: text, frames, guides, and a hollow marker's face."""
ROLE_CALLS = ("plotstyle.colour(", "plotstyle.shades(")
HEX = re.compile(r"#[0-9a-fA-F]{3}(?:[0-9a-fA-F]{3})?(?:[0-9a-fA-F]{2})?")
COLOUR_LETTER = re.compile(r"[bgrcmykw]")


def violations(source: str) -> list[str]:
    """Each place in `source` where a colour reaches a plot without its role."""
    out = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if HEX.fullmatch(node.value):
                out.append(f"line {node.lineno}: hex colour {node.value}")
        elif isinstance(node, ast.Attribute):
            if ast.unparse(node).startswith(("plt.cm", "matplotlib.cm", "cm.")):
                out.append(f"line {node.lineno}: colormap {ast.unparse(node)}")
        elif isinstance(node, ast.Call):
            fn = ast.unparse(node.func)
            for kw in node.keywords:
                if kw.arg not in COLOUR_KEYWORDS:
                    continue
                v = ast.unparse(kw.value)
                # Inside a role call, a bare name is one of `plotstyle.shades`'s colours.
                shade = fn.startswith("plotstyle.") and isinstance(kw.value, ast.Name)
                if not (shade or v in NEUTRAL or v.startswith(ROLE_CALLS)):
                    out.append(f"line {node.lineno}: {fn}(..., {kw.arg}={v})")
            if fn.endswith(".plot"):
                for a in node.args:
                    if isinstance(a, ast.Constant) and isinstance(a.value, str):
                        if COLOUR_LETTER.search(a.value):
                            out.append(f"line {node.lineno}: format string {a.value!r}")
    return out


def test_every_figure_takes_its_colours_from_the_roles():
    found = violations(REPORT.read_text())
    assert not found, "colours outside plotstyle's roles:\n" + "\n".join(found)


DECOYS = (
    'ax.plot(x, y, color="#c94f1d")',
    "ax.plot(x, y, color=plt.cm.viridis(0.5))",
    "ax.plot(x, y, color=ORANGE)",
    'ax.plot(x, y, "r--")',
    'ax.bar(x, y, facecolor="tab:green")',
)


def test_the_checker_catches_each_decoy():
    for d in DECOYS:
        assert violations(d), d


def test_the_roles_are_fixed():
    """Ballin orange, flight test grey dashed, GE green, ours blue, and an unstyled series can
    only come out blue."""
    sys.path.insert(0, str(ROOT / "validation"))
    import plotstyle as ps

    th = ps.LIGHT
    assert ps.line("ballin")["color"] == th.ballin
    assert ps.line("ours")["color"] == th.ours
    assert ps.points("ge")["color"] == th.ge
    flight = ps.line("flight")
    assert (flight["color"], flight["ls"]) == (th.muted, "--")
    assert all(ps.line(s)["ls"] == "-" for s in ("ballin", "ge", "ours"))
    cycle = ps.rc(th)["axes.prop_cycle"].by_key()["color"]
    assert set(cycle) == {th.ours}
