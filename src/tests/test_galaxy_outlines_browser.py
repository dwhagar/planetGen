# tests/test_galaxy_outlines_browser.py

"""
UX.135: what a browser draws from the galaxy family outlines (`static/galaxy/*.svg`, the `.galaxy-outline` rule of
`static/style.css`), on the fixture site of `tests/habitability_site.py`: each of the six is drawn in the page's text colour in the
light and the dark theme, as a line drawing (most of its box is the page behind it), and the six differ in shape, so no family is
told from another by anything but its drawing. The markup is the one `galaxy_outline(id)` gives a page, put on a page of the site.
"""

import io
import re

import numpy as np
import pytest
from PIL import Image

from planetgen.web.lib import galaxy_outlines

from tests.map_site_support import map_site  # noqa: F401
from tests.test_web_browser_fixture_maps import gl_browser, page  # noqa: F401 -- fixtures the lanes' own tests share

pytestmark = pytest.mark.browser

THEMES = ("light", "dark")
IDS = tuple(outline.id for outline in galaxy_outlines.OUTLINES)
SHAPE_PX = 96
"""int: How big the outlines are drawn to compare them."""
MAX_SHAPE_OVERLAP = 0.8
"""float: The most two outlines' pixels may overlap (intersection over union) and still count as different shapes."""
MAX_INK_SHARE = 0.35
"""float: The most of its box an outline may cover with its line: a line drawing, not a filled shape."""


def _parse(color):
    return tuple(float(value) for value in re.findall(r"[\d.]+", color)[:3])


def _open(page, base, theme):
    """The site's help page, in `theme`, with every outline put on it as `galaxy_outline` would, `SHAPE_PX` across."""
    page.goto(f"{base}/help/habitability", wait_until="load")
    page.emulate_media(color_scheme=theme)
    page.evaluate("t => document.documentElement.setAttribute('data-theme', t)", theme)
    markup = {outline_id: galaxy_outlines.outline_html(outline_id, lambda path: f"/static/{path}") for outline_id in IDS}
    page.evaluate("""([markup, size]) => {
        const holder = document.createElement('div');
        holder.id = 'outline-test';
        holder.style.cssText = 'position:fixed;top:0;left:0;background:var(--bg);padding:8px;z-index:9999;display:flex;gap:8px';
        for (const [id, html] of Object.entries(markup)) {
            const cell = document.createElement('div');
            cell.innerHTML = html;
            const svg = cell.firstChild;
            svg.style.width = size + 'px';
            svg.style.height = size + 'px';
            holder.appendChild(cell);
        }
        document.body.appendChild(holder);
    }""", [markup, SHAPE_PX])
    page.wait_for_function("document.querySelectorAll('#outline-test svg[data-galaxy-outline]').length === 6")


def _outline(page, outline_id):
    return page.locator(f'#outline-test svg[data-galaxy-outline="{outline_id}"]')


def _shot(locator):
    return np.asarray(Image.open(io.BytesIO(locator.screenshot())).convert("RGB")).astype(float)


def _solid(page, outline_id, colour, tries=30):
    """How many pixels of the outline are `colour`, once the file behind its `<use>` has loaded and painted (polled, as a file arrives late)."""
    count = 0
    for _ in range(tries):
        count = int((np.abs(_shot(_outline(page, outline_id)) - np.array(colour)).max(axis=2) <= 14).sum())
        if count:
            break
        page.wait_for_timeout(100)
    return count


def _ink(page, outline_id):
    """A boolean mask of the pixels of the outline's box that are not the page behind it."""
    background = _parse(page.evaluate("getComputedStyle(document.querySelector('#outline-test')).backgroundColor"))
    for _ in range(30):                          # polled, as the file behind the `<use>` may arrive late
        mask = np.abs(_shot(_outline(page, outline_id)) - np.array(background)).max(axis=2) > 40
        if mask.any():
            break
        page.wait_for_timeout(100)
    return mask


@pytest.mark.parametrize("theme", THEMES)
def test_each_outline_is_drawn_in_the_text_colour_as_a_line_drawing_with_its_name(page, map_site, theme):
    _open(page, map_site, theme)
    text = _parse(page.evaluate("getComputedStyle(document.body).color"))
    for outline_id in IDS:
        svg = _outline(page, outline_id)
        assert svg.get_attribute("role") == "img" and svg.get_attribute("aria-label") == galaxy_outlines.BY_ID[outline_id].label
        assert _parse(svg.evaluate("s => getComputedStyle(s).color")) == text, (outline_id, theme)
        solid = _solid(page, outline_id, text)
        assert solid >= 60, (outline_id, theme, solid)               # the file's line was drawn in the text colour
        share = _ink(page, outline_id).mean()
        assert 0.02 <= share <= MAX_INK_SHARE, (outline_id, theme, round(float(share), 3))


def test_the_six_outlines_differ_so_no_family_is_told_from_another_by_anything_but_its_drawing(page, map_site):
    _open(page, map_site, "light")
    masks = {outline_id: _ink(page, outline_id) for outline_id in IDS}
    for first_index, first in enumerate(IDS):
        for second in IDS[first_index + 1:]:
            overlap = (masks[first] & masks[second]).sum() / (masks[first] | masks[second]).sum()
            assert overlap < MAX_SHAPE_OVERLAP, (first, second, round(float(overlap), 3))
