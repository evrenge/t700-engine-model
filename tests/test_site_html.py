"""Structural checks on the generated site, for the defects that shipped on 2026-09-14.

`site/index.html` is built by `tools/build_site.py` from `site/template.html`. Neither has a
runtime here -- there is no browser in the test environment -- so these are static checks,
and they are deliberately written against two specific mistakes rather than as a general
HTML validator.

**The page was not a document.** The template began life as an Artifact, which is wrapped in
`<!doctype html><head>...</head><body>` by the host at publish time. Served raw from GitHub
Pages it had no doctype, no `<html>`, no `<head>` and no `<body>`, so every browser rendered
it in quirks mode.

**The theme toggle could not move the figures.** Each figure is a `<picture>` with a
`<source media="(prefers-color-scheme: dark)">`. A `<picture>` resolves its `<source>` before
the `<img>` and keeps that choice, so the toggle's `img.src = ...` changed nothing visible
while the source still matched: on a dark-preferring machine the page chrome went light and
all twelve figures stayed dark. The fix drives `source.media` instead, and this file checks
that the script still does.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / "site"
INDEX = SITE / "index.html"

pytestmark = pytest.mark.skipif(not INDEX.exists(), reason="site/index.html not built")


def html() -> str:
    return INDEX.read_text()


def test_the_page_is_a_complete_html_document():
    s = html()
    low = s.lower()
    assert low.lstrip().startswith("<!doctype html>"), "no doctype: renders in quirks mode"
    for tag in ("<html", "<head", "<body", "</body>", "</html>"):
        assert tag in low, f"missing {tag}"
    assert re.search(r"<html[^>]*\blang=", low), "<html> needs a lang attribute"
    assert low.index("<head") < low.index("<body"), "head must precede body"
    assert low.index("<style") < low.index("<body"), "styles belong in the head"


def test_every_figure_offers_both_themes():
    s = html()
    pictures = re.findall(r"<picture>.*?</picture>", s, re.S)
    assert pictures, "no <picture> elements -- the figures are not theme-aware"
    for pic in pictures:
        assert "data-dark" in pic, "the dark <source> needs the data-dark hook the script uses"
        assert "-dark.png" in pic and "-light.png" in pic, pic[:120]
        assert 'media="(prefers-color-scheme: dark)"' in pic, pic[:120]
        assert "<img " in pic and "alt=" in pic, "every figure needs alt text"


def test_the_toggle_drives_the_picture_source_and_not_the_img():
    """The specific bug: `img.src` cannot override a `<source>` that still matches."""
    s = html()
    script = s[s.index("<script>") : s.index("</script>")]
    assert "source[data-dark]" in script, (
        "the toggle must select the <source> elements; setting img.src leaves every figure "
        "on whatever the media query already chose"
    )
    assert re.search(r"\.media\s*=", script), "the toggle must rewrite source.media"
    assert "not all" in script, "ruling the dark source out needs a non-matching media query"
    assert not re.search(r"\bimg\.src\s*=", script), (
        "img.src assignment is the defect this test exists for"
    )


def test_the_stamp_is_removed_when_the_reader_returns_to_system():
    script = html()
    script = script[script.index("<script>") : script.index("</script>")]
    assert "removeAttribute" in script, (
        "apply(null) must clear data-theme, or the page can never follow the system again"
    )


def test_no_asset_is_referenced_that_is_not_committed():
    s = html()
    missing = [
        ref
        for ref in set(re.findall(r'(?:src|srcset)="(assets/[^"]+)"', s))
        if not (SITE / ref).exists()
    ]
    assert not missing, f"referenced but not present: {sorted(missing)}"


def test_the_build_is_current():
    """`index.html` must not be older than the template it is generated from."""
    tpl = SITE / "template.html"
    assert INDEX.stat().st_mtime >= tpl.stat().st_mtime - 1, (
        "site/template.html is newer than site/index.html -- run tools/build_site.py"
    )
