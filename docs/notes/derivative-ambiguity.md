# Why our Appendix B derivatives disagree, and why no interpolation fixes it

Established 2026-09-11, from the element-wise comparison in
`validation/test_appendix_b_elements.py`.

**Every remaining disagreement between our extracted linear model and Appendix B --
after the two Gen Hel terms are accounted for -- is explained by one thing: the trim
sits at or near a knot of the function table that controls that derivative.**

## The mechanism

Our function tables are interpolated linearly between knots (open question #43, a
recorded decision, not a transcription). A piecewise-linear interpolant reproduces the
*value* between knots to whatever the knots support, but its *derivative* is a step
function -- constant within a segment, discontinuous at every knot. Appendix B compares
derivatives.

So a trim that lands mid-segment has an unambiguous derivative, and one that lands on a
knot has a derivative that is genuinely undefined in our representation: it is one number
approaching from the left and another from the right.

## The evidence

Two independent families, controlled by two different tables, and both follow the rule.

| trim | distance to nearest `f7` knot | `f7`-driven errors | distance to nearest `f1` speed line | `f1`-driven errors |
|---|---|---|---|---|
| hover | **1.4 %** of the segment | **-24.3 %, +15.0 %** | 41.9 % | +0.6 %, -0.5 %, -0.1 % |
| level 80 kt | 28.3 % | -14.5 %, +2.6 % | **2.8 %** | -7.5 %, +5.2 %, +9.5 % |
| descent 80 kt | 32.5 % | -7.3 %, -3.0 % | **7.5 %** | +11.2 %, -12.2 %, -9.6 % |

Hover is the cleanest case because it is extreme in both directions at once: it sits
**0.0002** from `f7`'s knot at 0.21496 -- the worst `f7` errors in the comparison -- and
simultaneously 58 % across an `f1` segment, giving the *best* `f1` derivatives in the
comparison. Level and descent are the mirror image. A single cause producing opposite
outcomes at the same operating point is much harder to explain away than a trend.

Magnitude also needs the size of the discontinuity, not just the distance. Descent sits
further from its `f1` knot than level (7.5 % against 2.8 %) but has larger errors, because
the slope jump across the 85 line is -27.1 % against +7.6 % across the 89 line.

Signs check out too. At level we sit 2.8 % *above* the 89 line and use the segment above,
whose slope is 7.6 % steeper than the one below: our `|A(1,1)|` comes out too large,
-3.935 against Ballin's -3.659, i.e. -7.5 %. At descent we use a segment 27 % shallower
than the one below: `|A(1,1)|` too small, -3.282 against -3.696, i.e. +11.2 %.

## Why neither candidate fix works

**A smoother interpolant does not settle it.** Swapping `f7` alone to a shape-preserving
cubic moves exactly the four `f7`-driven elements and nothing else -- which is what proves
the causal link -- but it is nearly exact at level (-14.5 % -> +0.6 %), overshoots at hover
(-24.3 % -> +22.5 %) and worsens descent (-7.3 % -> -12.2 %). The implied true slopes are
about -423 at hover (between linear's -526 and the cubic's -330), -222 at level and -146 at
descent (below both).

**A larger extraction perturbation does not settle it either.** The report says the
matrices were extracted by perturbation [pdf p.31] and never states the step, and a step
large enough to straddle a knot averages the staircase automatically. Sweeping it, the rms
over five affected elements bottoms at:

| trim | best relative step | dNG | as %NGc | rms before -> after |
|---|---|---|---|---|
| hover | 3e-3 | 125 rpm | 0.279 % | 12.8 -> **1.3** |
| level 80 kt | 1e-2 | 398 rpm | 0.891 % | 8.8 -> 5.0 |
| descent 80 kt | 2e-2 | 761 rpm | 1.703 % | 9.3 -> 2.8 |

Hover's tenfold improvement is real and is exactly the straddle effect -- 125 rpm is what
it takes to cross that `f7` knot 0.0002 away. But **the best step is not the same at the
three trims**, and it cannot be, because the distance to the controlling knot is different
at each. A single step size is therefore not the explanation, and tuning one per trim
would be fitting, not modelling.

## What this means for the comparison

The residual is a property of **how finely the report printed its function tables**, not
of our equations. Supporting evidence: the trim *values* match Table B.1 to 0.13 % on NG
and 0.05 % on shaft power while the *derivatives* run 5-25 % out, which is precisely the
signature of an interpolation mismatch and not of a wrong equation. The zero structure
matches exactly, the P3 and P41 volume derivatives agree to under 1.5 %, and `b` agrees to
0.1 % wherever inertia does not enter.

**Do not change `maps.Curve` on this evidence.** It identifies the mechanism without
identifying the scheme, and changing the interpolation moves every value between knots in
every result in the repository, to fix derivatives that are ambiguous by construction.

What would settle it: the report stating its interpolation scheme or its extraction step
size. Neither is printed. Recorded against open question #43.

---

# Addendum: what the 3- and 6-DOF models added (2026-09-11)

## The structure is confirmed

The heat-sink linear models were derived rather than transcribed -- descriptor form
Eq. 65, Chen transformation Eqs. 66-67 -- and three printed features fall out of that
derivation without being put in, which is the check that it is the right one:

* **`G2` is zero in every row but T41's**, so `d = F1^-1 G2` is too. Appendix B prints
  exactly that in all six heat-sink figures.
* **`C = [ I ]`**, printed literally in all six, is what `z = x - B1 Wf` gives.
* **The 6-DOF `A(2,3)` is zero** where the 5-DOF's is `-0.4128E+2`. With T41 a separate
  state, P3 no longer reaches power-turbine torque except through T41, so that path moves
  out of the element and into the T41 column. We reproduce both.

It also answers **open question #32**: the report writes `b = F1^-1 G1 + F F1^-1 G2` with
`F` defined nowhere, and the transformation gives `b = B0 + A B1`, so `F == A == F1^-1 F2`.
Derived, not transcribed.

## tau1/tau2 is recoverable; tau2 is not

Open question #31 records that no numeric value of either time constant is printed. The
6-DOF T41 row is linear in `(1/tau2, tau1/tau2)`, so both look recoverable -- but they are
not equally determined:

* `tau2`'s independent signature in that row is worth about **18**, in elements of about
  **198 400**. The report prints four significant figures, i.e. +-50. **tau2 is below its
  own print precision** and a least-squares fit for it returns noise -- run naively it
  produces *negative* time constants at all three trims, which is how this was noticed.
* `tau1/tau2` rides the dominant term and is well determined: 0.6845 / 0.6571 / 0.6270,
  with the six row elements agreeing to 0.74 % at hover. An independent estimator, the
  feedthrough `d`, gives 0.6712 / 0.6450 / 0.6182.

**Ours is 0.6378 / 0.6084 / 0.5810 -- uniformly 7.3 % low**, ratio 0.932 / 0.926 / 0.927.
The uniformity is the evidence: a structural error would not be constant across three
flight conditions.

Equivalently our `tau_b/tau_a` is about 14 % too high, and

    tau_b/tau_a = f_hs(NGc) / (TC_T41 * sqrt(T41) * W41^0.2)

so the candidates are `f_hs` (digitized, Fig. A11) and `TC_T41 = 0.29` (transcribed,
Table A.1). **They cannot be separated on this data**: all three trims land on `f_hs`'s
flat top -- NGc 93.2 %, 89.1 % and 85.2 %, and the plateau runs 85 % to 100 % at 7.497 --
so `f_hs` is the same number at every one of them and is indistinguishable from a constant.
Separating them needs a trim below 85 % NGc, and the report gives none.

Recorded against open questions #31 and #4. Do not adjust either value on this evidence.
