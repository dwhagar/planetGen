"""UX.119: a Document written as HTML (`describe.html`) reads in both themes and at phone width.

One fixture Document (facts, ratings, composition and prose, with a very long unbroken name) is rendered, put into the
site's own sector page and checked in a real browser: every kind of text keeps contrast, the tier colours read, nothing
runs past the edge of a phone screen, and the quantities are `qty` spans the unit preference can change.
"""

import json

import pytest

from planetgen.describe import html, model
from tests.map_site_support import SECTOR_ID, map_site  # noqa: F401
from tests.test_web_browser_colors import CONTRAST
from tests.test_web_browser_fixture_maps import gl_browser, page  # noqa: F401

LONG_NAME = "Pneumonoultramicroscopicsilicovolcanoconiosis" * 3


def _document():
    return model.document("system", "system:1", "Fixture", [model.section(
        "planet", "planet:7", "Fixture II", [
            model.facts("Planet Data", [
                model.fact("distance", "Distance", model.qty("length", 1.496e11)),
                model.fact("mass", "Mass", model.qty("mass", 5.97e24)),
                model.fact("temperature", "Temperature", model.qty("temperature", 288.0)),
                model.fact("name", "Long name", LONG_NAME),
                model.fact("radius", "Radius", model.qty("radius", 6.371e6)),
                model.fact("host", "Host", model.ref("star:1", "Fixture A"))]),
            model.ratings({"name": "PHI-6", "score": 0.62, "tier": "Yellow"},
                          [{"domain": "Pressure", "tier": "Blue", "score": 1.0, "text": ["Ideal."]},
                           {"domain": "Temperature", "tier": "Green", "score": 0.8, "text": ["Cool."]},
                           {"domain": "Air", "tier": "Yellow", "score": 0.5, "text": ["Thin."]},
                           {"domain": "Radiation", "tier": "Red", "score": 0.1, "text": ["Dangerous."]}],
                          ["a mask"]),
            model.composition("Composition", [{"name": "iron", "share": 0.32}, {"name": "silicates", "share": 0.68}]),
            model.prose("summary", "A world of ", model.strong("rock"), " and ", model.qty("length", 3.0e8), ".")],
        [model.section("moon", "moon:1", "Fixture II a", [model.prose("summary", "A small moon.")])])])


INJECT = """
(markup) => {
  const box = document.createElement('div');
  box.id = 'doc119';
  box.className = 'panel';
  box.style.cssText = 'position:absolute;left:0;top:0;right:0;z-index:9999;background:var(--bg);padding:8px';
  box.innerHTML = markup;
  document.body.append(box);
}
"""


def _check_contrast(page, os_scheme, theme):
    page.evaluate("""() => {
      const mark = (selector, id) => { const el = document.querySelector(selector); el.id = id; };
      mark('#doc119 .doc-facts tbody th', 'c-row-head'); mark('#doc119 .doc-facts td', 'c-cell');
      mark('#doc119 .doc-ratings .tier-blue', 'c-blue'); mark('#doc119 .doc-ratings .tier-green', 'c-green');
      mark('#doc119 .doc-ratings .tier-yellow', 'c-yellow'); mark('#doc119 .doc-ratings .tier-red', 'c-red');
      mark('#doc119 .doc-prose', 'c-prose'); mark('#doc119 .doc-facts caption', 'c-caption');
    }""")
    failures = []
    for ident in ("c-row-head", "c-cell", "c-blue", "c-green", "c-yellow", "c-red", "c-prose", "c-caption"):
        ratio = page.evaluate(CONTRAST, ident)
        if ratio < 4.5:
            failures.append(f"{ident}: {ratio:.2f}")
    assert not failures, f"{os_scheme}/{theme}: " + "; ".join(failures)


@pytest.mark.parametrize("os_scheme,theme", [("dark", None), ("light", None), ("light", "dark"), ("dark", "light")])
def test_the_document_reads_in_both_themes(page, map_site, os_scheme, theme):  # noqa: F811
    page.emulate_media(color_scheme=os_scheme)
    page.goto(f"{map_site}/sector/{SECTOR_ID}", wait_until="load")
    if theme:
        page.evaluate("t => document.documentElement.setAttribute('data-theme', t)", theme)
    page.evaluate(INJECT, html.render(_document(), link=lambda ref: f"/{ref.replace(':', '/')}"))
    _check_contrast(page, os_scheme, theme)
    spans = page.evaluate("""() => [...document.querySelectorAll('#doc119 .qty')].map(
        (el) => [el.dataset.q, el.dataset.si, el.textContent.length > 0])""")
    assert [q for q, _si, _text in spans] == ["length", "mass", "temperature", "length"]
    assert all(float(si) > 0 and text for _q, si, text in spans)
    assert page.evaluate("() => document.querySelectorAll('#doc119 a[href=\"/star/1\"]').length") == 1


@pytest.mark.viewport(390, 800)
def test_the_document_stays_inside_a_phone_screen(page, map_site):  # noqa: F811
    page.goto(f"{map_site}/sector/{SECTOR_ID}", wait_until="load")
    before = page.evaluate("() => document.documentElement.scrollWidth")  # the sector page's own, not the Document's
    page.evaluate(INJECT, html.render(_document()))
    widths = page.evaluate("""() => {
      const box = document.getElementById('doc119').getBoundingClientRect();
      const wide = [...document.querySelectorAll('#doc119 *')].filter(
        (el) => el.closest('.table-scroll') === null && el.getBoundingClientRect().right > box.right + 1);
      return {right: box.right, view: window.innerWidth, doc: document.documentElement.scrollWidth,
              wide: wide.map((el) => el.tagName + '.' + el.className)};
    }""")
    assert widths["right"] <= widths["view"] + 1, json.dumps(widths)
    assert widths["doc"] <= max(before, widths["view"]) + 1, json.dumps({**widths, "before": before})
    assert not widths["wide"], json.dumps(widths)
