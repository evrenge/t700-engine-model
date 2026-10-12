"""Static checks on `site/index.html`, built by `tools/build_site.py` from `site/template.html`.

There is no browser here, so these target specific defects rather than validate HTML: the
page must be a complete document (no doctype means quirks mode), and the theme toggle must
drive each `<picture>`'s dark `<source media>`, because a `<picture>` keeps the `<source>`
it resolved and `img.src = ...` changes nothing visible.
"""

from __future__ import annotations

import re
import subprocess
import sys
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
    """`index.html` is exactly what `build_site.py` makes from the tree as it stands (by
    content: git does not keep mtimes)."""
    tool = SITE.parent / "tools" / "build_site.py"
    run = subprocess.run([sys.executable, str(tool), "--check"], capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr


# ------------------------------------------------------------------------------ navigation
# A link to a missing id, or a section missing from the table of contents, shows only by
# clicking.

_ANCHOR = re.compile(r'href="#([a-z0-9-]+)"')
_ID = re.compile(r'\bid="([a-z0-9-]+)"')


def test_every_in_page_link_resolves():
    s = html()
    ids = set(_ID.findall(s))
    broken = sorted({a for a in _ANCHOR.findall(s) if a not in ids})
    assert not broken, f"links to ids the page does not define: {broken}"


def test_the_table_of_contents_lists_every_section():
    """A section missing from the nav is reachable only by scrolling past everything above."""
    s = html()
    sections = set(re.findall(r'<section id="(s\d+)"', s))
    nav = s[s.index('<nav class="toc">') : s.index("</nav>")]
    listed = set(_ANCHOR.findall(nav))
    missing = sorted(sections - listed, key=lambda x: int(x[1:]))
    assert not missing, f"sections with no table-of-contents entry: {missing}"
