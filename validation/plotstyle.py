"""The Matplotlib style of the figures `validation/report.py` writes: one light theme, the IBM
Plex faces. Missing fonts fall back to DejaVu silently; `available()` reports it rather than
raising.
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

FIG_SCALE = 0.74
"""Every figure is drawn this much smaller in inches, at a correspondingly higher dpi, so the
content is larger relative to the frame, and the 0.8 pt grid lines survive a downscale."""

DPI = 190
"""Chosen with `FIG_SCALE`: a 16.5 in sheet becomes 12.2 in at 190 dpi, 2320 px wide."""

SANS = "IBM Plex Sans"
SANS_CONDENSED = "IBM Plex Sans Condensed"
MONO = "IBM Plex Mono"


@dataclass(frozen=True)
class Theme:
    """The colours of one theme."""

    name: str
    ground: str
    ink: str
    muted: str
    faint: str
    grid: str
    ours: str
    ballin: str
    ge: str


LIGHT = Theme(
    name="light",
    ground="#f7f8fa",
    ink="#14171c",
    muted="#4d5563",
    faint="#626a77",
    grid="#dadee5",
    ours="#2364b5",
    ballin="#c94f1d",
    ge="#0d8a5f",
)

# ------------------------------------------------------------------ one look per data source
#
# Every figure draws a source the same way: Ballin orange, flight test grey dashed, GE green,
# ours blue. Curves are lines, printed or computed numbers are markers in the source's colour.
# Orange, green and grey dashed never mark anything of ours.

SOURCES = ("ballin", "flight", "ge", "ours")
LINE_WIDTH = {"ballin": 1.2, "flight": 1.0, "ge": 1.2, "ours": 1.6}
FAMILY = ("-", "--", ":", "-.")
"""Several curves of one source in one panel: its colour, told apart by these line styles in
order, each labelled; `shades` beyond four."""


def colour(source: str, theme: Theme = LIGHT) -> str:
    """The colour of a source: `ballin`, `flight`, `ge` or `ours`."""
    return {
        "ballin": theme.ballin,
        "flight": theme.muted,
        "ge": theme.ge,
        "ours": theme.ours,
    }[source]


def line(source: str, theme: Theme = LIGHT, **kw) -> dict:
    """`plot` keywords for a curve of `source`; the flight test is always dashed."""
    out = {
        "color": colour(source, theme),
        "lw": LINE_WIDTH[source],
        "ls": "--" if source == "flight" else "-",
    }
    return out | kw


def points(source: str, theme: Theme = LIGHT, **kw) -> dict:
    """`plot` keywords for numbers of `source` drawn as markers: printed ones filled, ours
    hollow, GE's as plus signs unless a second GE set needs another shape."""
    c = colour(source, theme)
    out = {"ls": "none", "color": c, "ms": 5.0}
    if source == "ours":
        out |= {"marker": "o", "mfc": "none", "mew": 1.4}
    elif source == "ge":
        out |= {"marker": "+", "mew": 1.2}
    else:
        out |= {"marker": "o", "mfc": c, "mec": c}
    return out | kw


def guide(theme: Theme = LIGHT, **kw) -> dict:
    """`plot`/`axhline` keywords for a reference line that is no data: zero, 1:1, a bound."""
    return {"color": theme.faint, "lw": 0.8} | kw


def shades(source: str, n: int, theme: Theme = LIGHT) -> list[str]:
    """`n` shades of a source's colour, light to dark, for a family of more than four
    curves of one source: mixed with the ground from 60 % down to none."""
    rgb = matplotlib.colors.to_rgb(colour(source, theme))
    ground = matplotlib.colors.to_rgb(theme.ground)
    out = []
    for k in range(n):
        f = 0.6 * (1.0 - k / max(n - 1, 1))
        out.append(
            matplotlib.colors.to_hex(
                [(1.0 - f) * a + f * b for a, b in zip(rgb, ground, strict=True)]
            )
        )
    return out


def available() -> bool:
    """True when the IBM Plex faces are installed; False means a silent DejaVu fallback."""
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
        "axes.linewidth": 1.0,
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
        "grid.linewidth": 0.85,
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
        "lines.linewidth": 2.0,
        "lines.solid_capstyle": "round",
        "axes.prop_cycle": matplotlib.cycler(color=[theme.ours]),
    }


def title(fig, text: str) -> None:
    """A figure title: condensed, left, flush.

    The ink is read back out of rcParams, so a caller inside `plotstyle.context(...)` gets its
    theme's.
    """
    fig.suptitle(
        text,
        color=matplotlib.rcParams["text.color"],
        fontsize=11.5,
        fontweight=600,
        x=0.010,
        ha="left",
        fontfamily=[SANS_CONDENSED, SANS, "DejaVu Sans"],
    )


def mono_ticks(ax, theme: Theme) -> None:
    """Tick labels in the mono face."""
    for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        lab.set_fontfamily([MONO, "DejaVu Sans Mono", "monospace"])
        lab.set_fontsize(7.6)
        lab.set_color(theme.muted)


def context(theme: Theme):
    """`with plotstyle.context(theme):` around every figure this project draws."""
    return plt.rc_context(rc(theme))
