# `site/` — the published report

`index.html` and everything in `assets/` are **generated**. Do not hand-edit them; the next
build overwrites both. Edit `template.html` instead, then, inside the `t700-dev` container:

```bash
PYTHONPATH=src python validation/report.py   # 12 figures, light and dark, + assets/report.json
python tools/build_site.py                   # index.html
python tools/build_site.py --check           # is the committed index.html current?
```

`validation/report.py` re-runs the model for every number and every figure, so a rebuild is
a fresh measurement rather than a re-render. Each figure is drawn twice, light and dark,
styled by `validation/plotstyle.py` from the page's palette; the page picks one with a
`<picture>` media query, and a toggle overrides it.

`tools/build_site.py` fills the template from `assets/report.json`: named placeholders, and
`{{rj:<path>|<fmt>}}` anywhere for a number quoted in prose, so no measured number is typed
into the template. It stops on any placeholder left unfilled. Which open questions the page
lists, and their verdicts, are read from `docs/notes/open-questions.md` and checked against
the build's `OPEN_SUMMARY`, so the page cannot disagree with the ledger about what is open.
`tests/test_site_html.py` runs the `--check`.

`.github/workflows/pages.yml` publishes the committed `site/` on a push to `main` that
touches it (or the workflow), or on a manual run; it does not rebuild. No repository setting
is needed: its `configure-pages` step enables Pages on the first run.
