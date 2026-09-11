"""Python replication of NASA TM-100991 (Ballin, 1988).

A component-level real-time model of the General Electric T700-GE-700 turboshaft engine
and its fuel control system.

This package imports NumPy and the standard library and nothing else. SciPy, Matplotlib
and pandas belong in `tools/` and `validation/`; a test enforces the boundary.
"""

from __future__ import annotations

__version__ = "0.0.1"

REPORT = "NASA TM-100991 (Ballin, 1988)"
"""The source document. Every constant in this package cites a page of it."""
