# tests/test_life_text_a11y_browser.py

"""
UX.257: the Life section of an object page at phone width, in both themes, checked with axe-core on the site's own stylesheet.

`tests/test_web_a11y.py` runs every page, but its sample bodies carry no life text, so the two faults Planets 1 found (the readings
table 438 px wide on a 390 px phone, and the species link in the life text 2.49:1 light and 1.91:1 dark against the words around it)
never showed there. The section here is the markup of `templates/item.html`: a paragraph of the life text with a species link and a
help link, and the readings table in its scroll region. A change that drops the region or the underline fails here.
"""

import os
import re

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

pytestmark = pytest.mark.browser

HERE = os.path.dirname(os.path.abspath(__file__))
STATIC = os.path.join(HERE, "..", "html", "static")
AXE_JS = os.path.join(HERE, "vendor", "axe-core", "axe.min.js")
TEMPLATE = os.path.join(HERE, "..", "planetgen", "web", "templates", "item.html")
PHONE = {"width": 390, "height": 844}
RULES = ["link-in-text-block", "scrollable-region-focusable", "color-contrast"]


def _read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


SECTION = """
<main id="main"><section class="panel"><h2>About a planet</h2>
<p class="life-glance"><strong>Life.</strong> Microbes: rare, stable, secure. Conditions suit
<a class="help-link" href="#c1">Class 1a</a> (oxygenic light-eaters) best. <a class="help-link" href="#h">What these words mean</a></p>
<h3>Life</h3>
<p class="life-summary">Life here is <a class="help-link" href="#c1">Class 1a</a>, at the stage of small mobile animals.</p>
<div class="table-scroll" tabindex="0"><table class="life-readings">
<caption class="sr-only">Each kind of life on this planet</caption>
<thead><tr><th scope="col">Kind</th><th scope="col">How much</th><th scope="col">Trend</th><th scope="col">Security</th><th scope="col">Conditions last</th></tr></thead>
<tbody>
<tr><th scope="row">Flora</th><td>Occasional (3 of 9)</td><td>Declining</td><td>Precarious</td><td>about 0.4 billion years</td></tr>
<tr><th scope="row">Fauna</th><td>Rare (1 of 9)</td><td>Stable</td><td>Secure</td><td>about 0.4 billion years</td></tr>
</tbody></table></div>
<p class="sheet-help life-help-links">Help: <a class="help-link" href="#s">stages and animal grades</a>, <a class="help-link" href="#u">suitability</a>.</p>
<div class="prose"><p>The dominant species is the <a href="/object/species:1">Velarix</a>, who farm the shallow seas.</p></div>
</section></main>
"""


@pytest.fixture(scope="module")
def browser():
    try:
        playwright = sync_api.sync_playwright().start()
    except Exception as exc:  # noqa: BLE001 -- any failure here means "no browser"
        pytest.skip(f"Playwright could not start ({exc})")
    try:
        chromium = playwright.chromium.launch()
    except Exception as exc:  # noqa: BLE001
        playwright.stop()
        pytest.skip(f"No Chromium for Playwright ({str(exc).splitlines()[0]})")
    try:
        yield chromium
    finally:
        chromium.close()


def _open(browser, theme):
    page = browser.new_page(viewport=PHONE, color_scheme=theme)
    css = _read(os.path.join(STATIC, "style.css"))
    page.set_content(f'<!doctype html><html lang="en" data-theme="{theme}"><head><meta charset="utf-8"><title>Life</title>'
                     f'<meta name="viewport" content="width=device-width, initial-scale=1"><style>{css}</style></head><body>{SECTION}</body></html>')
    return page


def test_the_template_wraps_the_readings_table_in_a_focusable_scroll_region():
    template = _read(TEMPLATE)
    assert re.search(r'<div class="table-scroll" tabindex="0"><table class="life-readings">.*?</table></div>', template, re.S)


@pytest.mark.parametrize("theme", ["light", "dark"])
def test_the_life_section_passes_axe_and_fits_a_phone(browser, theme):
    page = _open(browser, theme)
    try:
        page.evaluate(_read(AXE_JS))
        result = page.evaluate("""(rules) => axe.run(document, {runOnly: {type: "rule", values: rules}, resultTypes: ["violations"]})""", RULES)
        problems = [f"{v['id']}: {[n['target'] for n in v['nodes'][:4]]}" for v in result["violations"]]
        assert not problems, problems
        assert page.evaluate("document.documentElement.scrollWidth") <= PHONE["width"]
        table_width = page.evaluate("document.querySelector('.life-readings').getBoundingClientRect().width")
        region_width = page.evaluate("document.querySelector('.table-scroll').getBoundingClientRect().width")
        assert table_width > region_width or table_width <= PHONE["width"]       # wide enough to need the region, or it fits
    finally:
        page.close()


@pytest.mark.parametrize("selector", [".prose p a", ".life-glance a.help-link", ".life-summary a.help-link", ".life-help-links a.help-link"])
def test_the_links_in_the_running_text_of_the_life_section_are_underlined(browser, selector):
    page = _open(browser, "light")
    try:
        assert "underline" in page.evaluate("s => getComputedStyle(document.querySelector(s)).textDecorationLine", selector)
    finally:
        page.close()
