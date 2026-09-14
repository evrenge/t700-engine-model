"""Where `data/` is, resolved once, for both source trees and installed packages.

`maps.py` and `appendix_b.py` each built their own path as
`Path(__file__).parent.parent.parent / "data"`. That is the repository root's `data/`
directory, which exists only in a source checkout: `pip install .` (non-editable) put
`t700` in site-packages and the first call to `maps.f1()` raised `FileNotFoundError`.
An editable install and pytest's `pythonpath = ["src"]` both leave the source tree on
the path, which is why nothing here ever noticed.

Two layouts are supported, and nothing else:

  * **source checkout** -- `<repo>/src/t700/_data.py` next to `<repo>/data/`
  * **installed distribution** -- `t700/data/` inside the package, which `pyproject.toml`
    produces by mapping the repository's `data/` directory onto the `t700.data` package

They are probed in that order, so a checkout always reads the checkout's data even if a
copy is also installed -- otherwise a stale wheel could silently supply the maps while
the digitizers rewrite the tree you are looking at.
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
