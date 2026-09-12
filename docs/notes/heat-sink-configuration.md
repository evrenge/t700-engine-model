# Which results used the heat-sink model, and which did not

Ballin ran the engine **in two configurations** and published results for both. Comparing
our output against a figure requires knowing which one produced it; getting this wrong is
a category error, not a tolerance problem.

Extracted 2026-09-11. Every row cites a page; the inference rows are marked as inferences.

## The table

| Report result | pdf page | Heat sink | Basis |
|---|---|---|---|
| Table 1, eigenvalues | p.31 | **off** | **printed**, p.27 |
| Appendix B, figs B1-B6 (2- and 5-DOF) | pp.68-70 | **off** | **printed**, p.67 |
| Appendix B, figs B7-B12 (3- and 6-DOF) | pp.71-76 | **on** | **printed**, p.67 |
| Figure 5 (linear doublet, 6-DOF vs 5-DOF) | p.34 | both | **printed** caption |
| Figures 9, 10 (fuel-step transients) | pp.45, 46 | **on** | **inferred** — see below |
| Table B.1 trims | p.67 | either | DC gain is 1; see "Why trims are exempt" |

## What is printed

Heat sink is optional in the formulation [pdf p.23, below Eq. 22]:

> The temperature at station 4.1 is expressed as a transfer function with slowly varying
> coefficients. (See the following subsection.) **If no heat-sink representation is used,
> T41 = T41ns.**

The linear analysis deliberately excluded it [pdf p.27, above Eq. 54]:

> Because the station 4.1 heat-sink approximation was added to both models as a linear
> lead-lag representation, **it was not included in the analysis.** The resulting five
> degrees-of-freedom represented were pressures at stations 3, 4.1, and 4.5 and the two
> turbine speeds.

Appendix B states the partition outright [pdf p.67]:

> The two- and three-degree-of-freedom models approximate the dynamics between the
> control volumes to be instantaneous; the five- and six-degree-of-freedom models contain
> complete dynamics. **The heat-sink model is contained in the three- and
> six-degree-of-freedom models.**

## What is inferred, and how strongly

**Figures 9 and 10 are heat-sink on.** The report never writes this sentence. Three
printed facts support it, and one measurement does:

1. [pdf p.20] "During the validation effort for the real-time simulation, some expansion
   of the model was found to be necessary... **a model of the nonadiabatic energy
   transfer at station 4.1 was required.**"
2. [pdf p.38] "**Because of the use of the heat-sink model**, station 4.1 temperature must
   be an independent input to the function table."
3. [pdf p.47], discussing Figure 10: "**Under a simplifying assumption, the real-time
   heat-sink model constants were made independent of the direction of power change.**"
   That sentence only parses if the heat sink is running in the figure under discussion.
4. Figure 9's T41 panel shows a spike decaying to a lower settled value. A `T41 = T41ns`
   run has **no mechanism** to produce that decay — the shape itself is the lead-lag.

**Settled by measurement, 2026-09-12, and it is not close.** The earlier line here quoted
three-point errors, which flatter. Integrating the error over the whole record instead,
normalised by each panel's own excursion, heat sink **on** wins on **every one of the nine
comparable panels** of Figures 9 and 10:

| | mean whole-curve RMS | worst panel |
|---|---|---|
| heat sink **on** | **3.79 %** | Fig. 10 T45, 7.62 % |
| heat sink off | 19.12 % | Fig. 10 T45, 34.53 % |

Figure 10's PCNG alone goes from 30.55 % to 2.35 %. (Re-measured 2026-09-12. This table
read 9.8 / 18.5 % with the worst-off panel named as PCNG; the "on" row moved when the heat
sink was rebuilt on Eqs. 48-49, the "off" row when the Figure 9/10 traces were re-digitized,
and the worst panel off is T45, not PCNG. The conclusion -- on wins on all nine -- is
unchanged and the margin is now far larger.) `validation/test_whole_curve.py` asserts
the margin so the inference stays measured rather than argued, and fails if it reverses.

**Caution on one sentence.** [pdf p.20] also says "No modeling of compressor surge,
heat-sink losses, or exhaust pressure losses was attempted." That refers to the
**predecessor NASA-Lewis hybrid model**, not Ballin's. Do not cite it as evidence the
heat sink is off.

## The sixth state is gas temperature, not metal temperature

[pdf p.29]:

> Inclusion of station 4.1 heat-sink effects in the small perturbation model requires the
> addition of an extra state. **This state, gas temperature at station 4.1,** is related
> to temperature without heat-sink dynamics by [Eq. 63]

`T_m` appears only in Eqs. 48-49 [pdf p.25] and is **eliminated** when those two equations
collapse into the transfer function Eq. 50 [pdf p.26]. The lead-lag carries the metal's
thermal inertia implicitly, and **Appendix B's sixth state is the gas temperature**, so the
linear models must be built on `T41`, not on `T_m`.

**This paragraph used to end "A model carrying a `T_m` state will not reproduce Appendix B."
That is false, and it is the one sentence here that could undo correct work.** The two are
not alternatives. `src/t700/linear.py` builds the linear models on the gas-temperature
state, exactly as p.29 requires, while the *nonlinear* real-time frame integrates the metal
temperature per Eqs. 48-49 (`RTState.hs_tm_degR`) -- which is what the report prints and
what removed the transient artifact, since Eq. 50's collapse is only valid for constant
coefficients and Eqs. 51/53 make sure they are not. All 297 Appendix B elements and Table
1's eigenvalues reproduce with that arrangement in place.

State vectors, both printed:

* 5-state [pdf p.27, below Eq. 54]: `{NG, NP, P3, P41, P45}`
* 6-state [pdf p.32, below Eq. 65]: `{NG, NP, P3, P41, P45, T41}`

The 2- and 3-state vectors are **not printed anywhere**. The 3-state is inferred as
`{NG, NP, T41}`.

## Why trims are exempt

Eq. 50 has **unit DC gain** by construction, so at equilibrium `T41 = T41_ns` regardless
of the switch. Every Table B.1 trim comparison is therefore valid for both configurations
simultaneously, and no trim result in this repository can be improved — or damaged — by
flipping it. `tests/test_heat_sink.py::test_heat_sink_is_inert_at_steady_state` measures
this rather than trusting the algebra: the two configurations agree to better than 1e-12
relative, against the test's own 1e-9 bound. (Stated as "~1e-15" until 2026-09-12; measured
2e-14 to 3e-13 depending on how long the frame is run before the comparison.)

## Consequences for our comparisons

* **Table 1 is the right target for our 5-DOF model** and the wrong target for a
  heat-sink-enabled one. Our existing eigenvalue comparison was correctly configured.
* **Figures 9 and 10 must be run heat-sink on.** Comparing our 5-DOF output against them
  was mis-targeted, and that was about **two thirds** of the transient error (open question
  #47). The rest was two further defects, both since fixed: an invented `tol=1e-10` on the
  pressure iterations where the report prints 0.1 percent in ten and eight passes, and
  carrying Eq. 50's collapsed lead-lag instead of integrating Eqs. 48-49. The Figure 9 T41
  overshoot went 3.41x -> 1.67x -> **0.87x** across those three steps.
* **The report prints no 3-DOF or 6-DOF eigenvalues anywhere.** If we want them we must
  compute them from the B7-B12 matrices ourselves.

## Not checked

* Whether Appendix C's fuel control couples to the heat sink.
* Which compressor/turbine function set produced Figures 9/10. The report says the
  NASA-Lewis test-engine functions replaced the specification functions for the Tables 2/3
  comparison [pdf p.39], but says nothing about Figures 9/10. Do not assume.
