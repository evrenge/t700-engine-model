---
description: Run the model against the digitized reference data and report deviations
---

Validate the current model against the report's data.

1. Run the validation suite: `pytest validation/ -v` inside the `t700-dev` distrobox.
2. For every reference case in `data/reference/`, run the model and compute the
   deviation against the reference trace.
3. Print a table: case, quantity, model value, report value, deviation, tolerance from
   `SCOPE.md`, pass/fail.
4. Write overlay plots (model vs report) to `validation/out/`.
5. Summarize: what passes, what fails, and for each failure your best hypothesis for the
   modelling difference behind it.

Rules, from CLAUDE.md:

- Do **not** adjust any constant, map value, or tolerance to make a case pass. A
  deviation is a finding to report, not a number to tune. If you believe a constant is
  wrong, say so with the report page that supports you and stop.
- Report deviations smaller than the digitization error recorded for that figure as
  "within read error", not as a pass on physics.
- If reference data for a case is missing, say so rather than skipping quietly.

$ARGUMENTS
