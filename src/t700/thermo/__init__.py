"""Gas property evaluation. The single door.

**No module outside this package may compute an enthalpy, a temperature from an
enthalpy, or a critical velocity ratio.** That rule (CLAUDE.md) is what lets the backend
be replaced without touching a line of component code, and it is enforced by
`tests/test_architecture.py`.

The current backend is `LinearFitThermo` -- Ballin's two-coefficient straight lines, in
pure NumPy. A C++ lookup-table backend can be installed with `set_backend()` when
profiling justifies it, though see `_linear.py`: there is nothing here to look up, so
such a backend would buy speed and nothing else.

    from t700 import thermo
    h = thermo.h41_from_t41(2100.0)
"""

from __future__ import annotations

from ._linear import LinearFitThermo

__all__ = [
    "LinearFitThermo",
    "backend",
    "h2_from_t2",
    "h3_from_t3",
    "h41_from_t41",
    "h45_from_h44",
    "set_backend",
    "t41_from_h41",
    "t45_from_h45",
    "t49_from_h49",
    "theta41_from_t41",
    "theta45_from_t45",
]

_backend = LinearFitThermo()


def backend() -> LinearFitThermo:
    """The active backend. Its `.name` identifies it in run provenance."""
    return _backend


def set_backend(new_backend) -> None:
    """Install a different property backend.

    Exists so a compiled backend can be dropped in later. Changing it changes model
    output, so a run that uses a non-default backend must say so in its results.
    """
    global _backend
    _backend = new_backend


def h2_from_t2(t2_degR):
    """Station 2 enthalpy, Btu/lbm."""
    return _backend.h2_from_t2(t2_degR)


def h3_from_t3(t3_degR):
    """Station 3 enthalpy, Btu/lbm."""
    return _backend.h3_from_t3(t3_degR)


def h41_from_t41(t41_degR):
    """Station 4.1 enthalpy, Btu/lbm."""
    return _backend.h41_from_t41(t41_degR)


def h45_from_h44(h44):
    """Station 4.5 enthalpy from station 4.4, Btu/lbm."""
    return _backend.h45_from_h44(h44)


def t41_from_h41(h41):
    """Station 4.1 temperature, deg R."""
    return _backend.t41_from_h41(h41)


def t45_from_h45(h45):
    """Station 4.5 (power turbine inlet) temperature, deg R."""
    return _backend.t45_from_h45(h45)


def t49_from_h49(h49):
    """Station 4.9 temperature, deg R."""
    return _backend.t49_from_h49(h49)


def theta41_from_t41(t41_degR):
    """Station 4.1 squared critical velocity ratio."""
    return _backend.theta41_from_t41(t41_degR)


def theta45_from_t45(t45_degR):
    """Station 4.5 squared critical velocity ratio."""
    return _backend.theta45_from_t45(t45_degR)
