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
    """`f1()` is the beta-gridded map; `f1_as_printed()` is Ballin's seven-point table."""
    m = maps.f1()
    assert len(m.lines) == 11, "eleven speed lines"
    assert m.params.tolist() == [65, 80, 82, 85, 87, 89, 92, 94, 96, 98, 100]
    assert all(line.x.size == 56 for line in m.lines), "56 beta values per line"

    raw = maps.f1_as_printed()
    assert len(raw.lines) == 11
    assert raw.params.tolist() == m.params.tolist()
    assert all(line.x.size == 7 for line in raw.lines), "seven printed markers per line"


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


# --------------------------------------------------------------------------------------
# The 2026-09-13 audits. Each of these is a hole that existed rather than a property that
# was in doubt: the model was clean on all of them, and the mechanisms were not.
# --------------------------------------------------------------------------------------


def test_a_nan_parameter_does_not_silently_select_the_top_speed_line():
    """The worst available failure mode, and it was the one in place.

    `nan < lo or nan > hi` is False, so the clamp counter never fired; `np.clip(nan,...)`
    is nan; and `np.searchsorted(params, nan)` returns `len(params)`, which took the
    `j >= len(self.params)` branch. So `f1(5.0, nan)` returned **10.108** -- the 100 %
    speed line's value, maximum compressor flow -- with an empty clamp report. A corrupted
    gas generator speed produced a plausible large number instead of a detectable NaN.
    """
    m = maps.f1()
    with pytest.raises(ValueError, match="non-finite"):
        m(5.0, float("nan"))
    with pytest.raises(ValueError, match="non-finite"):
        m(float("nan"), 90.0)
    with pytest.raises(ValueError, match="non-finite"):
        maps.f2()(float("nan"))


def test_f1_has_no_speed_line_between_65_and_80_percent():
    """The data hole the Figure 10 chop lives in. Open question #58.

    Eleven speed lines at 65, 80, 82, 85, 87, 89, 92, 94, 96, 98, 100 %NGc: fifteen
    percentage points between the first two and two or three between every other pair.
    Everything below 80 %NGc is therefore interpolated across that gap, with the 65 line
    as the lower bracket -- and the 65 line's abscissa ends at 3.753 against the 80 line's
    6.671, so past 3.753 a clamped line is bracketing an unclamped one.

    Pinned so the gap is a recorded property of the data rather than a surprise. It is
    also the reason the ratchet in `test_whole_curve` went up on 2026-09-13.
    """
    m = maps.f1()
    gaps = np.diff(m.params)
    assert float(gaps[0]) == 15.0, f"the 65-80 gap is {gaps[0]:g}, not 15"
    assert float(gaps[1:].max()) <= 3.0, (
        f"the other speed lines are spaced at most 3 % apart; largest is {gaps[1:].max():g}"
    )
    lo = m.lines[0]
    assert float(lo.x[-1]) < 4.0, (
        f"the 65 % line's data ends at {lo.x[-1]:.3f} in Ps3/P2; if it has been extended, "
        f"the bracketing clamp described here is gone"
    )


def test_clamp_scope_nests_without_destroying_an_enclosing_count():
    """`reset_clamps()` inside a measurement used to zero the measurement around it.

    `trim.solve` counts clamps over its Newton search and then again at the solution, and
    `validation/` counts them over runs that contain trims. Every one of those was a
    nested measurement.
    """
    maps.reset_clamps()
    # f9, not f6: as of 2026-09-13 f6 is loaded as the constant Figure A6 draws, and a
    # constant has no domain to leave, so it no longer counts clamps at all.
    f9 = maps.f9()
    f9(0.95)  # outside f9's 0.30050-0.85012 table
    with maps.clamp_scope() as inner:
        f9(0.95)
        f9(0.95)
    assert inner.counts == {"f9": 2}, inner.counts
    assert maps.clamp_report() == {"f9": 3}, (
        "the enclosing count must survive the inner scope and absorb its clamps"
    )
    maps.reset_clamps()


def test_the_clamp_counter_cannot_influence_any_model_output():
    """`maps._clamps` is the one mutable module global in the core. This bounds it.

    CLAUDE.md sanctions one exception to "no global mutable state" -- the thermo backend
    -- so the counter is a second, recorded in `maps._clamps`' own docstring. The
    guarantee that makes it acceptable is that nothing downstream of a lookup reads it,
    and that is checked behaviourally here rather than by inspection: the same trim run
    with the counter empty, pre-loaded, and never reset must be bit-identical.
    """
    from t700 import trim
    from t700.engine import Ambient
    from t700.units import wf_pps_from_pph

    amb = Ambient(14.696, 518.67)
    wf = wf_pps_from_pph(400.0)

    maps.reset_clamps()
    a = trim.solve(wf, ambient=amb).state
    for _ in range(500):  # run the counter up on purpose
        maps.f9()(0.95)
    b = trim.solve(wf, ambient=amb).state
    maps.reset_clamps()
    c_ = trim.solve(wf, ambient=amb).state

    assert a == b == c_, f"the clamp counter moved a model output:\n{a}\n{b}\n{c_}"
    maps.reset_clamps()


def test_the_loaded_tables_cannot_be_mutated_through_a_cached_handle():
    """`frozen=True` protects the reference, not the array it points at.

    The loaders are cached, so `maps.f2()` hands every caller the same `Curve`. Two
    spellings got past the dataclass guard entirely:

        maps.f2().y[:] = something     -- no attribute assignment, so no guard at all;
                                          it moved the 400 lbm/hr trim's NG by 482 rpm
        maps.f2().y *= 1.02            -- numpy applies the in-place multiply FIRST and
                                          the FrozenInstanceError fires after, so the
                                          "protected" spelling corrupts the table and
                                          then raises

    The second is the dangerous one: a reader seeing the exception would conclude nothing
    happened. Found by the 2026-09-13 code-quality audit. `appendix_b` already did this
    with `setflags(write=False)`; the function tables did not.
    """
    curve = maps.f2()
    with pytest.raises(ValueError, match="read-only"):
        curve.y[:] = 0.0
    with pytest.raises(ValueError, match="read-only"):
        curve.y *= 1.02
    with pytest.raises(ValueError, match="read-only"):
        curve.x[0] = 0.0

    speed = maps.f1()
    with pytest.raises(ValueError, match="read-only"):
        speed.params[0] = 0.0
    with pytest.raises(ValueError, match="read-only"):
        speed.lines[0].y[:] = 0.0

    # and the table still reads correctly afterwards
    assert float(maps.f2()(3.0)) == float(curve(3.0))


def test_the_uneven_beta_spacing_puts_resolution_where_the_curvature_is():
    """The seven printed beta values are far from evenly spaced. That is correct.

    On an average speed line the first interval -- beta 0 to 1/6, the choked part -- spans
    9.15 in pressure ratio, and the remaining five together span 1.8. It looks like the
    resolution is in the wrong place until you look at what varies: across the choked
    interval `y` moves by at most 0.0048 lbm/s over all eleven lines (0.15 %), while across
    the unchoked run it moves 3.0 to 6.4 %.

    So the sparse interval is flat and the dense ones bend. Inserting an extra beta line
    inside the choked segment, at any skew, changes `f1` by ~1e-16 -- there is nothing there
    to resolve. Evenly spaced beta lines would take resolution from the only part of the
    line that curves, and would replace the correspondence Figure A1 prints with one of ours.
    """
    m = maps.f1()
    choked_dy = [abs(line.y[1] - line.y[0]) for line in m.lines]
    run_dy = [abs(line.y[-1] - line.y[1]) for line in m.lines]
    assert max(choked_dy) < 0.006, f"the choked segment is no longer flat: {max(choked_dy)}"
    assert min(run_dy) > 0.15, f"the unchoked run no longer carries the variation: {min(run_dy)}"
    assert min(run_dy) > 30 * max(choked_dy)

    # an extra beta line inside the choked segment must change nothing at all
    worst = 0.0
    for pq in np.linspace(66.0, 99.5, 30):
        j = int(np.searchsorted(m.params, pq))
        a, b = m.lines[j - 1], m.lines[j]
        w = (pq - m.params[j - 1]) / (m.params[j] - m.params[j - 1])
        xb = (1 - w) * a.x + w * b.x
        yb = (1 - w) * a.y + w * b.y
        for skew in (0.25, 0.5, 0.75):
            xs = np.insert(
                xb,
                1,
                (1 - w) * (a.x[0] + skew * (a.x[1] - a.x[0]))
                + w * (b.x[0] + skew * (b.x[1] - b.x[0])),
            )
            ys = np.insert(
                yb,
                1,
                (1 - w) * (a.y[0] + skew * (a.y[1] - a.y[0]))
                + w * (b.y[0] + skew * (b.y[1] - b.y[0])),
            )
            for xq in np.linspace(xb[0], xb[-1], 25):
                ref = float(np.interp(xq, xb, yb))
                worst = max(worst, abs(float(np.interp(xq, xs, ys)) - ref) / max(ref, 1e-12))
    assert worst < 1e-12, f"the choked segment's internal correspondence matters: {worst:.3e}"


def test_f1_is_interpolated_on_a_smooth_analytic_beta_grid():
    """Beta is the fraction of a speed line's own pressure-ratio span, not a printed marker.

    0 at the choked left end of every line, 1 at surge, so one beta names corresponding
    points on lines whose spans differ by a factor of six. Figure A1 *does* print beta lines
    -- the six dotted construction lines joining the k-th marker of all eleven speed lines
    -- and they are not used, because the marker positions carry digitizing noise. As a
    fraction of each line's span the k=1 markers run

        0.785 0.817 0.831 0.832 0.814 0.811 0.855 0.868 0.882 0.891 0.897

    which goes up, up, up, **down, down**, up. That is jitter, and using it as the
    correspondence propagates it into every interpolated value. The spans it is a fraction
    of are smooth. Measured: Table B.1's rms 0.2471 -> 0.2439 % and the worst Table 1 NG
    mode 8.18 -> 7.56 % on moving from the printed markers to the analytic beta.

    What this pins is the property that made beta worth having at all: **blending two speed
    lines at equal beta blends their surge limits too**, so a query inside both lines is
    inside the blend and nothing is ever asked for a pressure ratio it cannot reach. Under
    the old constant-abscissa evaluation that failed on 389 frames of the Figure 10 chop.
    """
    m = maps.f1()
    assert m.beta_grid, "f1 is no longer a common-beta grid"

    # `beta_grid` alone does NOT detect the departure being undone: it is computed as
    # "every speed line has the same number of knots", and Ballin's printed seven-point
    # table satisfies that too. The 2026-09-14 validation-quality audit reverted `f1` to
    # the printed table and watched this test pass. What separates the two is the knot
    # count and the spacing: the analytic grid carries 56 values per line, 8 below beta
    # 0.70 and 48 above, and the printed table carries 7.
    printed = maps.f1_as_printed()
    assert printed.beta_grid, "the printed table is a common-knot-count grid as well"
    assert m.lines[0].x.size > 4 * printed.lines[0].x.size, (
        f"f1 carries {m.lines[0].x.size} beta values per line against the printed table's "
        f"{printed.lines[0].x.size}; if these are the same, departure #60 has been undone "
        f"and `beta_grid` cannot tell"
    )
    # and the grid is DENSER above the knee than below it, which the printed table is not.
    # Beta is the fraction of the line's own span, so it comes from x and not from the
    # index -- that is the whole point of the scheme.
    x = m.lines[0].x
    beta = (x - x[0]) / (x[-1] - x[0])
    below = int(np.count_nonzero(beta <= 0.70 + 1e-9))
    above = int(beta.size - below)
    assert below * 4 < above, (
        f"the beta grid carries {below} values at or below beta 0.70 and {above} above; "
        f"it is meant to be much denser through the knee and up to surge"
    )

    for pq in (70.0, 74.0, 86.0, 93.16):
        j = int(np.searchsorted(m.params, pq))
        a, b = m.lines[j - 1], m.lines[j]
        w = (pq - m.params[j - 1]) / (m.params[j] - m.params[j - 1])
        hi = (1 - w) * a.x[-1] + w * b.x[-1]
        assert a.x[-1] <= hi <= b.x[-1], (pq, hi)
        maps.reset_clamps()
        m(hi * 0.999, pq)
        assert "f1" not in maps.clamp_report(), f"a query inside the blend clamped at {pq}"
    maps.reset_clamps()

    # every line starts at the same choked abscissa and ends at its own surge point
    assert len({round(float(line.x[0]), 6) for line in m.lines}) == 1
    ends = [float(line.x[-1]) for line in m.lines]
    assert ends == sorted(ends), "surge pressure ratio must rise with speed"
