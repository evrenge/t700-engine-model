"""Check the fuel control constants against the printed table's structure.

Table C.1 is printed across two pages -- 28 rows then 29 -- so the count is a real check
on the transcription, not a tautology. If someone drops a row while editing, or a
duplicate creeps in, the count moves.

The ordering assertions matter more than they look. Every limit pair in this control
system is a low/high bracket, and a transposed pair is the kind of error that produces a
model which runs, saturates in the wrong direction, and looks merely "detuned".
"""

from __future__ import annotations

import math

import pytest

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
    public = [n for n in public if n != "SOURCE_TABLE"]
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


@pytest.mark.parametrize(
    "low,high",
    [
        ("XLOLIM", "XHILIM"),
        ("YLOLIM", "YHILIM"),
        ("ZLOLIM", "ZHILIM"),
        ("WFMIN", "WFMAX"),
        ("WFPDCL", "WFPDCH"),
        ("NP_LOOP_RELAY_LOW", "NP_LOOP_RELAY_HIGH"),
    ],
)
def test_limit_pairs_are_ordered(low: str, high: str):
    assert getattr(k, low) < getattr(k, high), f"{low} must be below {high}"


def test_time_constants_are_positive():
    """A negative lag is an unstable pole, not a slow one."""
    for n in (
        "CLMV",
        "CLLDS",
        "CNTL",
        "CT2",
        "CT7",
        "CT9",
        "CT12",
        "CT13",
        "CT14",
        "CT16",
        "CTPL",
        "CTPS3",
        "T8",
        "T10",
        "T11",
        "T17",
        "TL1",
        "TL2",
        "TLGE",
        "TM_LEAD",
        "TM_LAG",
        "FUEL_TRANSPORT_DELAY",
    ):
        assert getattr(k, n) > 0.0, n


def test_the_three_known_negative_constants():
    """These are negative in the report and must stay negative.

    A dropped minus sign is the single most likely OCR failure on this scan, so the ones
    that are genuinely negative are pinned explicitly.
    """
    assert k.BWFP < 0
    assert k.XLOLIM < 0
    assert k.ZLOLIM < 0
    assert k.TM_INPUT_BIAS < 0
    assert k.TM_CURRENT_BIAS < 0
    assert k.ZLOLIM_HIGH_TORQUE < 0


def test_ocr_traps_have_the_corrected_values():
    """Two values the text layer reports wrongly; both read from the page image.

    Recorded as tests so a future 'tidy-up' against the OCR text cannot silently undo
    them.
    """
    assert k.CORR == 20.0, "OCR text layer says 2010; the page says 20.0"
    assert k.T10 == 0.060, "OCR text layer says 0.050; the page says 0.060"


def test_zlolim_switch_matches_table_value():
    """Fig. C8 switches ZLOLIM on engine-2 torque; the low branch matches Table C.1.

    With the one-engine implementation switch open, engine 2 torque is 0, which is below
    the 180 ft*lbf threshold -- so the Table C.1 value is the one in force. That the two
    agree is a consistency check on both readings.
    """
    assert k.ZLOLIM == -1.0
    assert k.ZLOLIM_HIGH_TORQUE == -0.3
    assert k.ENGINE2_TORQUE_THRESHOLD == 180.0


def test_the_ten_millisecond_lags_are_still_lags_at_the_report_frame():
    """Six blocks at tau = 0.010 s against a 7 ms frame: dt/tau = 0.7. See `_blocks.lag`.

    Explicit Euler puts the pole at `1 - dt/tau`, so these six sit at z = +0.300 with an
    effective time constant of 5.814 ms -- **41.9 % short**. At `realtime.MAX_STEP_S`
    (10 ms) the pole is exactly zero and the blocks respond in one frame; at `FRAME_NP_S`
    (14 ms) it is negative and they alternate instead of lagging. Both of those frames are
    named by the report (open question #33), so neither degenerate case is hypothetical.

    Isolated, the discretization is worth at most 0.094 % on the closed-loop settled state
    (see `_blocks.lag`), because every affected block is a sensor lag far above the loop
    bandwidth. The condition is recorded, not removed: explicit Euler is the scheme the
    model uses everywhere, by the architecture rule in CLAUDE.md.

    This test exists so that a change of frame cannot make them deadbeat quietly.
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
    assert pole > 0.0, (
        "the 10 ms lags have gone deadbeat or worse at the shipped frame; they are no "
        "longer first-order lags and every Appendix C sensor path is affected"
    )

    # the two frames the report itself names, for the record
    assert 1.0 - realtime.MAX_STEP_S / 0.010 == 0.0, "10 ms is exactly deadbeat for these"
    assert 1.0 - realtime.FRAME_NP_S / 0.010 < 0.0, "14 ms puts the pole negative"
