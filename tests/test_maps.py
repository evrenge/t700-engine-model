"""Check the function tables load, evaluate, and behave physically.

These tests are the last gate before the tables feed a model. They check structure and
physics, not values -- the values were checked against the page by the verify overlays,
which is the only place they can be checked.
"""

from __future__ import annotations

import re

import numpy as np
import pytest

from t700 import maps
from t700._data import DATA_ROOT

ALL_CURVES = ("f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "f_hs")


@pytest.mark.parametrize("name", ALL_CURVES)
def test_curve_loads_with_provenance(name: str):
    c = getattr(maps, name)()
    assert c.x.size == c.y.size
    assert "TM-100991" in c.source, f"{name} lost its citation on load"


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


def test_f6_is_loaded_as_the_constant_figure_a6_draws():
    """Open question #62's departure, pinned so it cannot be undone without a test failing.

    Nothing detected it before 2026-09-14. `test_combustor_efficiency_is_near_unity` above
    passes on the digitized slope too -- the two printed endpoints differ by less than a
    tenth of their own read error -- so every range check in this file is blind to whether
    `f6` is a constant or a line. That is the same hole
    `test_f1_is_interpolated_on_a_smooth_analytic_beta_grid` exists to close for #60, found the
    same way: revert the departure and watch the suite stay green.

    Four things are asserted, and the last two are what keep the provenance honest.

    **The expected value is computed from the CSV, not frozen into this file.** It used to
    be `0.5 * (0.985038 + 0.984965)`, a pair of numbers from an extraction two revisions
    old: the 2026-09-16 blend of the two ticked edges redistributed those endpoints while
    leaving their mean untouched, so the assertion went on passing against numbers the file
    no longer held, and only the fine-lattice recalibration made it fail. A frozen number
    that survives the data changing underneath it is not pinning anything. What is actually
    being asserted -- that the loader takes the MEAN of the two measured endpoints, and
    that those endpoints are the same number to within the page's own read error -- is
    asserted directly instead.
    """
    c = maps.f6()
    assert c.y.size >= 2
    assert len(set(c.y.tolist())) == 1, (
        f"f6 is not constant: {c.y}. Figure A6 draws one horizontal line, and the tilt on "
        f"the page is the page (open question #62). If this is deliberate, move #62 to "
        f"closed and say why here."
    )
    # the CSV still carries the ink as measured: the conditioning is what makes it
    # constant, so the chain back to the page is unbroken and the reproducibility gate
    # still has something to reproduce.
    cols, _, _ = maps._read_csv(DATA_ROOT / "maps" / "f6_combustor_efficiency.csv")
    ys = cols["y"].tolist()
    sy = cols["sigma_y"].tolist()
    assert len(set(ys)) == 2, f"the f6 CSV should still hold the two measured points: {ys}"
    # the value is the mean of the two measured endpoints, not either one of them
    assert c.y[0] == pytest.approx(0.5 * (ys[0] + ys[1]), abs=1e-12), c.y[0]
    # and #62's premise: those two endpoints ARE one number, as far as the page can say
    gap = abs(ys[0] - ys[1]) / min(sy)
    assert gap < 0.5, (
        f"f6's two printed endpoints now differ by {gap:.2f} of their own 1 sigma read "
        f"error. Open question #62 rests on them being the same number; if the page really "
        f"does draw a slope, #62 is answered and this departure has to go."
    )
    assert maps.CONDITIONING["f6"] == pytest.approx(abs(ys[0] - c.y[0]), rel=1e-6)


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


def test_f1_parameter_clamping_is_counted():
    maps.reset_clamps()
    m = maps.f1()
    m(6.0, 50.0)  # below the lowest speed line
    m(6.0, 120.0)  # above the highest
    assert maps.clamp_report().get("f1:parameter", 0) == 2
    maps.reset_clamps()


# --------------------------------------------------------------------------- physics


CONDITIONING_TOL_SIGMA = 3.0
"""How far physical conditioning may move a value, in that map's OWN read error.

Was an absolute 1e-3 for every map until 2026-09-16, which is not what the docstring below
says and not what it measured. These maps carry y in units that differ by four orders of
magnitude, so one absolute number means wildly different things:

    f9   1e-3 = 35.49 sigma_y     f3   1e-3 =  5.92 sigma_y
    f6   1e-3 =  2.46 sigma_y     f1   1e-3 =  0.09 sigma_y

It was caught by f1. Blending `f1`'s two ticked edges instead of pooling them -- which
makes its eleven physically-flat extensions flatter, 0.573 px mean end-to-end difference
against 0.409 px, a free check the calibration is never told about -- left the flat pair of
one speed line 2.3e-03 apart, so the monotone step moves each end by 1.15e-03. Against the
absolute bound that FAILED. Against f1's own read error it is **0.098 sigma**, a tenth of
the noise the value is known to. The bound was wrong, not the data.

Three sigma is stricter than the old bound for three of the four maps that condition
anything, and the honest reading of "read-error scale" for the fourth."""


def _read_error(name: str) -> float | None:
    """That map's mean quoted 1 sigma on y, from its own CSV. Found, not listed."""
    # `f_hs` is stored as `fhs_heat_sink_constant.csv`, so try the name with and without
    # its underscore rather than keeping a name-to-file table that can go stale.
    stems = {name, name.replace("_", "")}
    for sub in ("maps", "schedules"):
        for stem in stems:
            for path in (DATA_ROOT / sub).glob(f"{stem}_*.csv"):
                cols, header, _ = maps._read_csv(path)
                col = [c for c in cols if c in ("sigma_y", "sy")]
                if col:
                    return float(np.mean(cols[col[0]]))
                # `f1` quotes its read error as a table in the header rather than as
                # columns, because it varies by marker index and speed line rather than by
                # point. Take the largest, which is the class the conditioning acts on.
                head = "\n".join(header)
                quoted = re.findall(r"^#\s+k=\S[^\n]*?(\d+\.\d+)\s*$", head, re.M)
                if quoted:
                    return max(float(q) for q in quoted)
    return None


def test_conditioning_changes_are_tiny_everywhere():
    """Physical conditioning must be a correction, never a rewrite.

    If any of these grows past read-error scale, the digitization is wrong and the fix
    belongs at the figure, not here. "Read-error scale" means that map's own quoted
    sigma -- see `CONDITIONING_TOL_SIGMA` for why it cannot be one absolute number.
    """
    for name in ALL_CURVES + ("f10",):
        getattr(maps, name)()
    maps.f1_as_printed()
    for name, moved in maps.CONDITIONING.items():
        sigma = _read_error(name)
        assert sigma, f"{name}: no quoted read error to judge its conditioning against"
        assert moved < CONDITIONING_TOL_SIGMA * sigma, (
            f"{name}: physical conditioning moved a value by {moved:.2e}, which is "
            f"{moved / sigma:.2f} of its own 1 sigma read error ({sigma:.2e}). That is too "
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
    for name in ("f3", "f7", "f8", "f9"):
        y = getattr(maps, name)().y
        assert np.all(np.diff(y) <= 0), f"{name} must not rise"


def test_conditioning_leaves_the_csv_alone():
    """Provenance is unbroken: the file on disk still holds what was read off the page.

    Conditioning happens at load. If someone ever 'fixes' the data files instead, the
    chain back to the printed figure is lost and this test fails.
    """
    cols, _, _ = maps._read_csv(maps.DATA_DIR / "f3_seal_bleed_fraction.csv")
    raw = cols["y"].tolist()
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


F1_MARKER_TOL_PCT = 1.0
F1_MARKER_BIAS_PCT = 0.15


def test_the_printed_f1_markers_lie_on_the_traced_map():
    """`f1_compressor_mass_flow.csv` against the beta map traced off the same figure.

    The beta map is fitted mostly to the traced line, so the two are near-independent reads
    of Figure A1. Measured worst 0.67 % at one marker, mean +0.038 % over all 77. Until
    2026-10-12 nothing in pytest read the seven-point file: +5 % on its whole wa2c column
    failed no test.
    """
    raw, traced = maps.f1_as_printed(), maps.f1()
    devs = []
    for lr, lt in zip(raw.lines, traced.lines, strict=True):
        dev = (lr.y - np.interp(lr.x, lt.x, lt.y)) / lr.y * 100.0
        assert np.abs(dev).max() < F1_MARKER_TOL_PCT, f"{lr.name}: {np.round(dev, 2)} %"
        devs.append(dev)
    bias = float(np.concatenate(devs).mean())
    assert abs(bias) < F1_MARKER_BIAS_PCT, f"mean marker deviation {bias:+.3f} %"


def test_blending_at_the_printed_beta_values_is_exact():
    """Resampling every line onto a finer beta grid that keeps the knots changes nothing.

    Both coordinates are piecewise linear in beta with shared breakpoints, so the blended
    line is the same polyline either way. Cited by `maps.SpeedMap`'s docstring.
    """
    raw = maps.f1_as_printed()
    beta = np.linspace(0.0, 1.0, raw.lines[0].x.size)
    fine = np.union1d(np.linspace(0.0, 1.0, 50), beta)
    dense = maps.SpeedMap(
        name=raw.name,
        params=raw.params.copy(),
        lines=[
            maps.Curve(name=ln.name, x=np.interp(fine, beta, ln.x), y=np.interp(fine, beta, ln.y))
            for ln in raw.lines
        ],
    )
    worst = 0.0
    with maps.clamp_scope():
        for pq in np.linspace(65.0, 100.0, 71):
            for xq in np.linspace(1.0, 16.0, 61):
                worst = max(worst, abs(raw(xq, pq) - dense(xq, pq)))
    assert worst < 1e-12, f"resampled blend differs by {worst:.2e}"


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
    lo = maps.f1().lines[0]
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
    assert inner == {"f9": 2}, inner
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
    drop = [line.y[1] - line.y[-1] for line in m.lines]  # signed: flow FALLS above the knee
    assert max(choked_dy) < 0.006, f"the choked segment is no longer flat: {max(choked_dy)}"
    assert min(drop) > 0.15, f"the unchoked run no longer carries the variation: {min(drop)}"
    assert min(drop) > 30 * max(choked_dy)


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
    # The knot counts that separate this grid from the printed table (56 against 7 per
    # line) are pinned by `test_f1_structure`.
    # The grid is DENSER above the knee than below it, which the printed table is not.
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
