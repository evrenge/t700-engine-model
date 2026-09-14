# `site/` — the published report

`index.html` and everything in `assets/` are **generated**. Do not hand-edit them; the next
build overwrites both. Edit `template.html` instead, then:

```bash
distrobox enter t700-dev
PYTHONPATH=src python validation/report.py   # 24 figures + assets/report.json
python tools/build_site.py                   # index.html
```

`validation/report.py` re-runs the model for every number and every figure, so a rebuild is
a fresh measurement rather than a re-render. It writes each figure twice — light and dark,
styled by `validation/plotstyle.py` from the same palette the page uses — and the page picks
with a `<picture>` media query, with a toggle that overrides it.

`tools/build_site.py` fills the template from `assets/report.json` and cross-checks the
open-questions section against `docs/notes/open-questions.md`: if a row opens or closes and
the summary here is not updated, the build fails rather than publishing a stale claim.

Published by `.github/workflows/pages.yml` on any push that touches this directory. No
repository setting is needed — the workflow's `configure-pages` step enables Pages itself on
the first run.
