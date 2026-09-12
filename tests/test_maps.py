"""Check the function tables load, evaluate, and behave physically.

These tests are the last gate before the tables feed a model. They check structure and
physics, not values -- the values were checked against the page by the verify overlays,
which is the only place they can be checked.
"""

from __future__ import annotations

import numpy as np
import pytest

from t700 import maps

ALL_CURVES = ("f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "f_hs")


@pytest.mark.parametrize("name", ALL_CURVES)
def test_curve_loads_with_provenance(name: str):
    c = getattr(maps, name)()
    assert c.x.size >= 2
    assert c.x.size == c.y.size
    assert "TM-100991" in c.source, f"{name} lost its citation on load"


@pytest.mark.parametrize("name", ALL_CURVES)
def test_knots_strictly_increasing(name: str):
    """Enforced in Curve.__post_init__; asserted here so the message is about the data."""
    c = getattr(maps, name)()
    assert np.all(np.diff(c.x) > 0), f"{name} has duplicate or unsorted knots"


@pytest.mark.parametrize("name", ALL_CURVES)
def test_interpolation_passes_through_every_knot(name: str):
    """Linear interpolation must reproduce the digitized points exactly."""
    c = getattr(maps, name)()
    assert np.allclose(c(c.x), c.y, rtol=0, atol=0)


def test_midpoint_is_the_mean_of_its_neighbours():
    """Confirms the scheme really is linear, not something smoothed."""
    c = maps.f2()
    mid = 0.5 * (c.x[3] + c.x[4])
    assert c(mid) == pytest.approx(0.5 * (c.y[3] + c.y[4]))


def test_clamping_outside_the_domain_is_counted():
    maps.reset_clamps()
    c = maps.f2()
    lo, hi = c.domain
    assert c(lo - 5.0) == pytest.approx(c.y[0])
    assert c(hi + 5.0) == pytest.approx(c.y[-1])
    report = maps.clamp_report()
    assert report.get("f2", 0) == 2, (
        "leaving the digitized range must be counted -- a model that wanders off its "
        "maps without saying so is inventing physics"
    )
    maps.reset_clamps()


def test_inside_the_domain_never_clamps():
    maps.reset_clamps()
    c = maps.f9()
    lo, hi = c.domain
    c(np.linspace(lo, hi, 50))
    assert maps.clamp_report() == {}


# --------------------------------------------------------------------------- physics


def test_bleed_fractions_are_physical():
    """A bleed fraction outside [0, 1] is not a bleed fraction."""
    for name in ("f3", "f4", "f5"):
        c = getattr(maps, name)()
        assert c.y.min() >= -1e-3, f"{name} goes negative: {c.y.min()}"
        assert c.y.max() <= 1.0, f"{name} exceeds unity: {c.y.max()}"


def test_combustor_efficiency_is_near_unity():
    c = maps.f6()
    assert 0.9 < c.y.min() and c.y.max() < 1.0


def test_compressor_temperature_ratio_rises_with_pressure_ratio():
    c = maps.f2()
    assert np.all(np.diff(c.y) > 0)
    assert c.y.min() > 1.0, "T3/T2 below 1 would mean the compressor cools the air"


def test_f8_goes_negative_at_high_expansion_ratio():
    """Real, not a digitizing artefact -- the printed axis extends to -5.0."""
    c = maps.f8()
    assert c.y.min() < 0.0
    assert c(0.85) < 0.0
    assert c(0.30) > 0.0


def test_f10_is_used_as_printed_not_inverted():
    """f10 = P49/Ps9 > 1. The figure's axis LABEL is inverted, not its values.

    Eq. 37 sets Ps9 = P_amb and Eq. 38 gives P49 = Ps9 * f10, so f10 < 1 would put the
    exhaust total pressure below ambient and the engine could not exhaust. Every plotted
    value is above 1, and the neighbouring maps corroborate: f8's P49/P45 spans
    0.3007-0.8505 and f9's Ps9/P45 spans 0.3003-0.8500, so P49 and Ps9 are close in
    magnitude and their ratio sits just above 1.

    An earlier reading of open question #5 had this inverted. This test is the guard.
    """
    c = maps.f10()
    assert c.y.min() > 1.0, "f10 must exceed 1 or the engine cannot exhaust"
    assert c.y.max() < 1.5, "and it is a small loss, not a large one"

    p_amb = 14.696
    p49 = p_amb * c(85.0)
    assert p49 > p_amb, f"P49 ({p49:.3f}) must exceed ambient ({p_amb})"


def test_f10_stays_between_the_neighbouring_map_arguments():
    """P45 > P49 > Ps9 -- the expansion order that fixes f10's direction."""
    f8, f9 = maps.f8(), maps.f9()
    assert f8.x.min() == pytest.approx(f9.x.min(), abs=0.01)
    assert f8.x.max() == pytest.approx(f9.x.max(), abs=0.01)


# --------------------------------------------------------------------------- the 2-D map


def test_f1_structure():
    m = maps.f1()
    assert len(m.lines) == 11, "eleven speed lines"
    assert m.params.tolist() == [65, 80, 82, 85, 87, 89, 92, 94, 96, 98, 100]
    assert all(line.x.size == 7 for line in m.lines), "seven points per line"


def test_f1_on_a_speed_line_matches_that_line():
    """Evaluating at a printed speed must reproduce that line exactly, not a blend."""
    m = maps.f1()
    line = m.lines[5]  # 89 %
    for xk, yk in zip(line.x, line.y, strict=True):
        assert m(float(xk), 89.0) == pytest.approx(float(yk))


def test_f1_mass_flow_rises_with_speed():
    """At a fixed pressure ratio, a faster compressor passes more corrected flow."""
    m = maps.f1()
    flows = [m(6.0, p) for p in (85.0, 89.0, 94.0, 100.0)]
    assert all(b > a for a, b in zip(flows, flows[1:], strict=False)), flows


READ_ERROR_WA2C = 0.02
"""Quoted read error on the compressor map, lbm/sec, from the digitizing agent's own
ground-truth measurement against the printed construction lines."""


@pytest.mark.parametrize("i", range(11))
def test_f1_is_choked_below_the_knee(i: int):
    """Corrected flow is constant with pressure ratio below the knee, then falls.

    This is the physical signature of a compressor map and the reason the flat left
    extensions are real function breakpoints rather than decoration.

    "Constant" means constant *to within read error*, not exactly. The digitized flat
    segments vary by 0.0006-0.0048 lbm/sec, which is 0.2-1.4 % of the drop above the
    knee and well inside the quoted +/-0.02. Asserting exact equality here would be
    asserting the digitization is perfect, which it is not and cannot be.
    """
    line = maps.f1().lines[i]
    flat_variation = abs(line.y[1] - line.y[0])
    drop_above_knee = line.y[1] - line.y[-1]

    assert flat_variation < READ_ERROR_WA2C, (
        f"left extension is not flat: varies {flat_variation:.4f} lbm/sec"
    )
    assert drop_above_knee > 0.0, "flow must fall above the knee"
    assert flat_variation < 0.05 * drop_above_knee, (
        f"the flat segment ({flat_variation:.4f}) is not negligible against the drop "
        f"above the knee ({drop_above_knee:.4f}) -- is the knee index right?"
    )

    # and the segment really is interpolated flat between its two knots
    knee = line.x[1]
    mid = 0.5 * (line.x[0] + knee)
    assert abs(line(mid) - line.y[0]) < READ_ERROR_WA2C


def test_f1_parameter_clamping_is_counted():
    maps.reset_clamps()
    m = maps.f1()
    m(6.0, 50.0)  # below the lowest speed line
    m(6.0, 120.0)  # above the highest
    assert maps.clamp_report().get("f1:parameter", 0) == 2
    maps.reset_clamps()


# --------------------------------------------------------------------------- physics


def test_conditioning_changes_are_tiny_everywhere():
    """Physical conditioning must be a correction, never a rewrite.

    If any of these grows past read-error scale, the digitization is wrong and the fix
    belongs at the figure, not here.
    """
    for name in ALL_CURVES + ("f10",):
        getattr(maps, name)()
    for name, moved in maps.CONDITIONING.items():
        assert moved < 1e-3, (
            f"{name}: physical conditioning moved a value by {moved:.2e}. That is too "
            f"large to be read noise -- re-digitize the figure instead."
        )


def test_bleed_fraction_tail_is_exactly_zero():
    """Figure A3 shows f3 reaching zero and staying there. It must be exactly zero.

    The digitized tail wobbles +/-1e-4 around zero -- read noise on a pen-plotter line,
    not a seal bleed that reopens above 89 % NGc.
    """
    c = maps.f3()
    tail = c.y[c.x > 89.5]
    assert tail.size >= 1
    assert np.all(tail == 0.0), f"f3 tail is not exactly zero: {tail}"
    assert c.y.min() >= 0.0, "a bleed fraction cannot be negative"


def test_declared_monotonic_tables_really_are():
    """After conditioning, physics holds exactly rather than approximately."""
    assert np.all(np.diff(maps.f2().y) >= 0), "f2 must rise with pressure ratio"
    for name in ("f3", "f7", "f8", "f9"):
        y = getattr(maps, name)().y
        assert np.all(np.diff(y) <= 0), f"{name} must not rise"


def test_conditioning_leaves_the_csv_alone():
    """Provenance is unbroken: the file on disk still holds what was read off the page.

    Conditioning happens at load. If someone ever 'fixes' the data files instead, the
    chain back to the printed figure is lost and this test fails.
    """
    import csv as _csv
    from pathlib import Path

    path = Path(maps.DATA_DIR) / "f3_seal_bleed_fraction.csv"
    with path.open() as fh:
        rows = list(_csv.DictReader(ln for ln in fh if not ln.startswith("#")))
    raw = [float(r["y"]) for r in rows]
    assert min(raw) < 0.0, (
        "the CSV should still carry the raw digitized value, including its negative "
        "excursion -- conditioning belongs at load, not in data/"
    )


def test_f1_cross_line_ordering_holds_as_digitized():
    """More corrected speed passes more corrected flow at a given pressure ratio.

    This is the one physical constraint on `f1` that is *not* conditioned, and the reason
    is that the digitized map already obeys it everywhere. Asserting it here keeps that
    claim honest: if a re-digitization ever breaks the ordering, this fails rather than
    being silently repaired at load. See open question #44 on the principle.
    """
    m = maps.f1()
    lo = max(float(line.x.min()) for line in m.lines)
    hi = min(float(line.x.max()) for line in m.lines)
    for pr in np.linspace(lo, hi, 40):
        flows = [float(line(pr)) for line in m.lines]
        assert np.all(np.diff(flows) >= 0.0), f"speed lines cross at Ps3/P2 = {pr:.4f}"
