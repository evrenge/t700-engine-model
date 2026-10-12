"""Check the fuel control constants against the printed table's structure.

Table C.1 is printed across two pages -- 28 rows then 29 -- so the count is a real check
on the transcription, not a tautology. If someone drops a row while editing, or a
duplicate creeps in, the count moves.
"""

from __future__ import annotations

import math

from t700.control import constants as k

TABLE_C1_ROWS = 57
"""28 rows on pdf p.83 plus 29 on pdf p.84."""

FIGURE_ONLY_CONSTANTS = 16
"""19 numbers appear only on block diagrams. Three are not new constants here: the
ZLOLIM below-threshold value duplicates the Table C.1 row, and 100/NGDES and 1/3600 are
a design speed and a unit conversion, which live in t700.constants and t700.units."""


def _names():
    """Split the module's public constants by where they came from.

    Table C.1 entries keep the report's own symbols, which carry no underscore.
    Figure-only constants were named by us and all contain one.
    """
    public = [n for n in dir(k) if n.isupper() and not n.startswith("_")]
    table = [n for n in public if "_" not in n]
    figure = [n for n in public if "_" in n]
    return table, figure


def test_table_c1_row_count():
    table, _ = _names()
    assert len(table) == TABLE_C1_ROWS, (
        f"expected {TABLE_C1_ROWS} Table C.1 constants (28 on pdf p.83 + 29 on p.84), "
        f"found {len(table)}: {sorted(table)}"
    )


def test_figure_only_constants_are_present():
    _, figure = _names()
    assert len(figure) == FIGURE_ONLY_CONSTANTS, sorted(figure)


def test_every_constant_is_a_number():
    table, figure = _names()
    for n in table + figure:
        assert isinstance(getattr(k, n), (int, float)), n


def test_time_constants_are_positive():
    """A negative lag is an unstable pole, not a slow one.

    Only the time constants pinned nowhere else: the printed ones are pinned by value in
    `test_control_blocks.py::PRINTED`, the 10 ms ones by the test below."""
    for n in ("CT12", "CT16", "T11"):
        assert getattr(k, n) > 0.0, n


def test_the_ten_millisecond_lags_are_still_lags_at_the_report_frame():
    """Six blocks at tau = 0.010 s against a 7 ms frame: dt/tau = 0.7. See `_blocks.lag`.

    Explicit Euler puts the pole at `1 - dt/tau`: z = +0.300, an effective time constant
    41.9 % short. At `realtime.MAX_STEP_S` (10 ms) the pole is zero and the blocks respond in
    one frame; at `FRAME_NP_S` (14 ms) it is negative and they alternate. The report names
    both frames (open question #33). Recorded, not removed: explicit Euler is the model's
    scheme. This test stops a change of frame from making them deadbeat quietly.
    """
    from t700 import realtime

    ten_ms = {
        "CT9": k.CT9,
        "CTPL": k.CTPL,
        "CTPS3": k.CTPS3,
        "T17": k.T17,
        "TL1": k.TL1,
        "TL2": k.TL2,
    }
    for name, tau in ten_ms.items():
        assert tau == 0.010, f"{name} is {tau}, not the 0.010 s this test is about"

    dt = realtime.FRAME_ENGINE_S
    pole = 1.0 - dt / 0.010
    assert abs(pole - 0.300) < 1e-12, f"dt/tau is no longer 0.7: pole {pole}"
    tau_eff = -dt / math.log(pole)
    assert abs(100.0 * (tau_eff / 0.010 - 1.0) + 41.9) < 0.1, (
        f"the effective time constant is {tau_eff * 1000:.3f} ms, "
        f"{100 * (tau_eff / 0.010 - 1):+.1f} % against -41.9 % on record"
    )
    assert pole > 0.0, "the 10 ms lags are deadbeat or worse at the shipped frame"

    # the two frames the report itself names, for the record
    assert 1.0 - realtime.MAX_STEP_S / 0.010 == 0.0, "10 ms is exactly deadbeat for these"
    assert 1.0 - realtime.FRAME_NP_S / 0.010 < 0.0, "14 ms puts the pole negative"
