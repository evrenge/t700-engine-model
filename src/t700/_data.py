"""Where `data/` is, resolved once, for both source trees and installed packages.

Two layouts, probed in this order:

  * source checkout -- `<repo>/src/t700/_data.py` next to `<repo>/data/`
  * installed distribution -- `t700/data/` inside the package, which `pyproject.toml`
    maps from the repository's `data/`

Checkout first, so a stale installed wheel never supplies the maps while the digitizers
rewrite the tree. This is the one module allowed to know where data is.
"""

from __future__ import annotations

from pathlib import Path

_HERE = Path(__file__).resolve().parent

_CANDIDATES = (
    _HERE.parent.parent / "data",  # source checkout: <repo>/data
    _HERE / "data",  # installed: site-packages/t700/data
)


def _resolve() -> Path:
    for candidate in _CANDIDATES:
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(
        "t700 cannot find its data directory. Looked in:\n  "
        + "\n  ".join(str(c) for c in _CANDIDATES)
        + "\nA source checkout keeps it at <repo>/data; an installed distribution keeps "
        "it inside the package. Neither is present, so the function tables, the "
        "schedules and the Appendix B matrices cannot be loaded."
    )


DATA_ROOT = _resolve()
"""The directory holding `maps/`, `schedules/`, `linear/`, `reference/`."""
