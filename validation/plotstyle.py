"""Matplotlib styled to match the report page, so the figures belong to the document.

`validation/report.py` writes the figures that `site/index.html` embeds. Before this module
they were Matplotlib's defaults on a white ground -- DejaVu Sans, a blue-orange-green cycle
that happened to match the page's accents only because the page was built around them, and
a fixed light background that glared inside a dark-themed page.

Three things it fixes.

**One palette, named once.** The page defines its colours as CSS custom properties; the same
values live here, so `--ours` is one hex in both places and a series keeps its colour from a
table cell to a plot mark.

**The page's own typefaces.** IBM Plex Sans Condensed for titles, Plex Sans for labels, Plex
Mono for tick numbers -- the superfamily the page is set in, which Fedora packages and the
`Containerfile` installs. If they are missing Matplotlib falls back to DejaVu silently;
`available()` reports that rather than raising, because a missing font is not a reason to
stop a validation run.

**Two themes.** Every figure is rendered twice, light and dark, and the page picks with a
`<picture>` media query. The dark variant is not an inversion: the grounds come from the
page's dark tokens and the three series colours are lifted for contrast on them, exactly as
the CSS does.
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

SANS = "IBM Plex Sans"
SANS_CONDENSED = "IBM Plex Sans Condensed"
MONO = "IBM Plex Mono"


@dataclass(frozen=True)
class Theme:
    """One column of the page's token table."""

    name: str
    ground: str
    surface: str
    ink: str
    muted: str
    faint: str
    grid: str
    ours: str
    ballin: str
    ge: str
    ge2: str


LIGHT = Theme(
    name="light",
    ground="#ffffff",
    surface="#f6f7f9",
    ink="#16191d",
    muted="#5d6673",
    faint="#8b94a1",
    grid="#e2e6ec",
    ours="#2a78d6",
    ballin="#d95a22",
    ge="#12996a",
    ge2="#7c5ccc",
)
DARK = Theme(
    name="dark",
    ground="#1b1f26",
    surface="#20252d",
    ink="#e7eaef",
    muted="#9aa4b2",
    faint="#727c8a",
    grid="#2e343e",
    ours="#5d9dea",
    ballin="#f0834d",
    ge="#2ec18c",
    ge2="#a98ae4",
)
THEMES = (LIGHT, DARK)


def available() -> bool:
    """True when the page's own faces are installed; False means a silent DejaVu fallback."""
    names = {f.name for f in fm.fontManager.ttflist}
    return {SANS, MONO}.issubset(names)


def rc(theme: Theme) -> dict:
    """The full rcParams for one theme. Applied with `plt.rc_context`, never globally."""
    family = [SANS, "DejaVu Sans", "sans-serif"]
    return {
        "figure.facecolor": theme.ground,
        "figure.edgecolor": theme.ground,
        "savefig.facecolor": theme.ground,
        "savefig.edgecolor": theme.ground,
        "axes.facecolor": theme.ground,
        "axes.edgecolor": theme.grid,
        "axes.labelcolor": theme.muted,
        "axes.titlecolor": theme.ink,
        "axes.linewidth": 0.9,
        "axes.grid": True,
        "axes.axisbelow": True,
        "axes.titlesize": 10.5,
        "axes.titleweight": "semibold",
        "axes.titlelocation": "left",
        "axes.titlepad": 9.0,
        "axes.labelsize": 9.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "grid.color": theme.grid,
        "grid.linewidth": 0.8,
        "grid.alpha": 1.0,
        "xtick.color": theme.faint,
        "ytick.color": theme.faint,
        "xtick.labelcolor": theme.muted,
        "ytick.labelcolor": theme.muted,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "text.color": theme.ink,
        "font.family": family,
        "font.size": 9.5,
        "mathtext.fontset": "dejavusans",
        "legend.frameon": False,
        "legend.fontsize": 8.0,
        "legend.labelcolor": theme.muted,
        "legend.handlelength": 1.6,
        "legend.handletextpad": 0.6,
        "legend.columnspacing": 1.2,
        "figure.titlesize": 11.5,
        "figure.titleweight": "semibold",
        "lines.solid_capstyle": "round",
        "axes.prop_cycle": matplotlib.cycler(color=[theme.ours, theme.ballin, theme.ge, theme.ge2]),
    }


def title(fig, text: str, theme: Theme | None = None) -> None:
    """A figure title set like the page's section headings: condensed, left, flush.

    `theme` defaults to whatever `rc_context` is currently in, read back out of rcParams,
    so a caller inside `plotstyle.context(...)` need not pass it.
    """
    ink = theme.ink if theme is not None else matplotlib.rcParams["text.color"]
    fig.suptitle(
        text,
        color=ink,
        fontsize=11.5,
        fontweight=600,
        x=0.010,
        ha="left",
        fontfamily=[SANS_CONDENSED, SANS, "DejaVu Sans"],
    )


def mono_ticks(ax, theme: Theme) -> None:
    """Tick labels in the mono face, as every number on the page is."""
    for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        lab.set_fontfamily([MONO, "DejaVu Sans Mono", "monospace"])
        lab.set_fontsize(7.6)
        lab.set_color(theme.muted)


def context(theme: Theme):
    """`with plotstyle.context(theme):` around every figure this project draws."""
    return plt.rc_context(rc(theme))
