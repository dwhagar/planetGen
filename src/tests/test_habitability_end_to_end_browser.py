# tests/test_habitability_end_to_end_browser.py

"""
TEST.300: the habitability system checked together in a browser, on a real database (`tests/habitability_site.py`): a
sector with a seeded system that has Blue or Green worlds beside a sector whose system has none, every body's PHI-6
rating stored by the generator, the Galaxy Map's block totals and habitable-system rows recorded.

Each lane tests its own piece against fixture data (UX.100 the symbols, MAP.143 the colouring, MAP.220 the markers,
UX.181 the help links); these four check that the pieces agree on what the generator stored:

1. the Planets and Moons lists draw every body's stored rating as the right symbol, in the colours of both themes;
2. a sector is coloured by its habitable-world count, the stage the map reads carries the stored count and the
   legend's scale holds it;
3. the habitable-world marker is off until its kind is switched on, finds the system with the right count (and not the
   one without), and goes when the kind is switched off;
4. a symbol opens the help page at its own section, and the page's levels table is the same one the score reads.
"""

import io
import re
import time

import numpy as np
import pytest
from PIL import Image

from planetgen.db import store
from planetgen.physics import habitability
from planetgen.physics.habitability_world import TIER_NAMES
from planetgen.web.lib import habitability_ref

from tests.habitability_site import BARREN_CENTER_PC, HABITABLE_CENTER_PC, SYSTEM_OFFSET_PC, habitable_galaxy  # noqa: F401
from tests.test_web_browser_fixture_maps import (  # noqa: F401 -- fixtures and the map drivers the lanes' own tests share
    MARKERS, _markers_until, _open_galaxy, _open_menu, _screen_of, _settle, gl_browser, page)

pytestmark = pytest.mark.browser

THEMES = ("light", "dark")
TIER_WORDS = {name: index for index, name in enumerate(TIER_NAMES)}
"""dict: Tier name -> the stored tier index (0 Blue to 3 Red)."""


def _stored_tiers(config, table):
    """`{name: overall tier name}` of every planet or moon with a stored overall score."""
    conn = store.get_connection(config)
    try:
        return {row["name"]: habitability.tier(row["phi4"])
                for row in conn.execute(f"SELECT name, phi4 FROM {table} WHERE phi4 IS NOT NULL").fetchall()}
    finally:
        conn.close()


def _row_count(config, table):
    conn = store.get_connection(config)
    try:
        return conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
    finally:
        conn.close()


def _listed(page, base, kind, rows):
    """`{name: the seven symbols' words}` of every row of a list. The table draws only the rows in view (and fetches them as they
    come into view), so this scrolls it from the top to the bottom and gathers each row as it is drawn."""
    page.goto(f"{base}/lists/{kind}", wait_until="load")
    page.wait_for_selector("tbody tr td:first-child a", timeout=30000)
    gathered = {}
    deadline = time.monotonic() + 60
    while time.monotonic() < deadline:
        state = page.evaluate("""() => {
            const scroller = document.querySelector('.datatable-scroll');
            const drawn = Array.from(document.querySelectorAll('tbody tr')).filter((row) => row.querySelector('td:first-child a')).map((row) => [
                row.querySelector('td:first-child a').textContent.trim(),
                Array.from(row.querySelectorAll('svg.hab-symbol')).map((s) => s.getAttribute('aria-label'))]);
            const before = scroller.scrollTop;
            scroller.scrollTop = before + scroller.clientHeight / 2;
            return {drawn, atBottom: scroller.scrollTop === before};
        }""")
        gathered.update({name: symbols for name, symbols in state["drawn"]})
        if len(gathered) >= rows and state["atBottom"]:
            break
        page.wait_for_timeout(100)
    assert len(gathered) == rows, (len(gathered), rows)
    return gathered


def _parse(color):
    return tuple(float(v) for v in re.findall(r"[\d.]+", color)[:3])


def _luminance(rgb):
    def channel(v):
        v /= 255.0
        return v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4
    r, g, b = (channel(v) for v in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _contrast(a, b):
    high, low = sorted((_luminance(a), _luminance(b)), reverse=True)
    return (high + 0.05) / (low + 0.05)


# --- 1. The lists ----------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("kind,table", (("planet", "planets"), ("moon", "moons")))
def test_the_lists_draw_each_bodys_stored_rating_as_the_symbol_of_its_tier(page, habitable_galaxy, kind, table):
    stored = _stored_tiers(habitable_galaxy.config, table)
    listed = _listed(page, habitable_galaxy.base_url, kind, _row_count(habitable_galaxy.config, table))
    missing = {name: listed.get(name) for name in stored if not listed.get(name)}
    assert not missing, (missing, len(listed), page.locator("tbody tr").count(), page.locator("tbody").inner_html()[:1500])
    assert {name for name, symbols in listed.items() if symbols} == set(stored) and len(stored) > 3
    assert all(not symbols for name, symbols in listed.items() if name not in stored)     # a body with no score shows none
    for name in stored:
        symbols = listed[name]
        assert len(symbols) == 7, (name, symbols)
        assert symbols[0] == f"Overall habitability: {stored[name]}", (name, symbols[0], stored[name])
        assert all(re.fullmatch(r"[A-Za-z ]+: (Blue|Green|Yellow|Red)", words) for words in symbols), symbols
    if kind == "planet":
        assert {tier for tier in stored.values()} & set(habitability.POTENTIALLY_HABITABLE_TIERS), "the seeded system has no habitable planet"
        habitable_names = {name for _kind, _id, name in habitable_galaxy.habitable_bodies if _kind == "planet"}
        assert {name for name, tier in stored.items() if tier in habitability.POTENTIALLY_HABITABLE_TIERS} >= habitable_names


def test_the_tier_colours_read_in_both_themes(page, habitable_galaxy):
    """The four tiers draw four colours that are told apart and stand off the page in each theme, and the ones that
    mean 'habitable' look cooler than the ones that do not (Blue and Green against Yellow and Red)."""
    page.goto(f"{habitable_galaxy.base_url}/lists/planet", wait_until="load")
    page.wait_for_selector("svg.hab-symbol", state="attached")
    seen = {}
    for theme in THEMES:
        page.emulate_media(color_scheme=theme)
        page.evaluate("t => document.documentElement.setAttribute('data-theme', t)", theme)
        colors = page.evaluate("""() => { const o = {}; document.querySelectorAll('svg.hab-symbol').forEach((s) => {
            const tier = Array.from(s.classList).find((c) => c.startsWith('tier-'));
            o[tier] = getComputedStyle(s).color; }); return o; }""")
        background = _parse(page.evaluate("getComputedStyle(document.body).backgroundColor"))
        assert set(colors) == {"tier-blue", "tier-green", "tier-yellow", "tier-red"}, (theme, colors)
        rgb = {tier: _parse(value) for tier, value in colors.items()}
        for tier, value in rgb.items():
            assert _contrast(value, background) >= 3.0, (theme, tier, value, background)    # WCAG 1.4.11, a graphic that carries meaning
        tiers = sorted(rgb)
        for index, first in enumerate(tiers):
            for second in tiers[index + 1:]:
                assert sum(abs(a - b) for a, b in zip(rgb[first], rgb[second])) > 60, (theme, first, second)
        blue, green, red = rgb["tier-blue"], rgb["tier-green"], rgb["tier-red"]
        assert blue[2] > blue[0] and green[1] > green[0] and red[0] > red[1], (theme, rgb)
        seen[theme] = rgb
    assert seen["light"] != seen["dark"]


# --- 2. The sector colouring ------------------------------------------------------------------------------------------

def _shot(page):
    return np.asarray(Image.open(io.BytesIO(page.locator("#galaxymap3d-canvas").screenshot())).convert("RGB")).astype(float)


def _patch(image, page, center_pc, half=8):
    """The mean colour of the patch of the map around a sector's center."""
    spot = page.evaluate("(p) => document.querySelector('#galaxymap3d-canvas').galaxyScreenOf(...p)", list(center_pc))
    box = page.locator("#galaxymap3d-canvas").bounding_box()
    x, y = int(spot["x"] - box["x"]), int(spot["y"] - box["y"])
    return image[y - half:y + half + 1, x - half:x + half + 1].reshape(-1, 3).mean(axis=0)


def _legend_number(text):
    return float(text.replace(",", ""))


def test_a_sector_is_coloured_by_its_habitable_worlds_and_the_legend_holds_the_stored_count(page, habitable_galaxy):
    # The stage the map reads carries the stored count for the block that holds both sectors.
    _open_galaxy(page, habitable_galaxy.base_url, "")
    stage = page.evaluate("fetch('/galaxy/stage').then((r) => r.json())")
    blocks = [child for child in stage["children"] if child["stats"]["systems"] > 0]
    assert sum(child["stats"]["habitable_worlds"] for child in blocks) == habitable_galaxy.habitable_worlds
    assert sum(child["stats"]["systems"] for child in blocks) == 2

    # The three-sector window round the habitable sector, the info panel closed, the menu shut.
    page.evaluate("([cell, g]) => document.querySelector('#galaxymap3d-canvas').galaxyOpenWindow(cell, g)",
                  [{"m": 1, "shell": 1, "band": 0, "wedge": 0}, 3])
    _settle(page)
    close = page.locator("#galaxymap3d-info .map-info-close")
    if close.is_visible():
        close.click()
        page.wait_for_timeout(300)
    before = _shot(page)
    _open_menu(page)
    page.select_option("#galaxymap3d-color-by", "habitable")
    _open_menu_closed(page)
    _settle(page)
    after = _shot(page)

    legend = page.locator("#galaxymap3d-legend")
    assert legend.is_visible()
    assert "Habitable worlds" in legend.locator(".galaxy-legend-title").inner_text()
    low, high = (_legend_number(legend.locator(f".galaxy-legend-{end}").inner_text()) for end in ("low", "high"))
    per_sector = habitable_galaxy.habitable_worlds / 2.0       # the block holds two charted sectors, one with the worlds
    assert low <= per_sector <= high, (low, per_sector, high)
    assert "color=habitable" in page.url

    # The habitable sector turns green and the one without worlds does not.
    green_shift = _hue_shift(before, after, page, HABITABLE_CENTER_PC)
    barren_shift = _hue_shift(before, after, page, BARREN_CENTER_PC)
    assert green_shift > barren_shift + 3.0, (green_shift, barren_shift)


def _open_menu_closed(page):
    if page.locator("#galaxymap3d-menu").evaluate("m => m.open"):
        page.locator("#galaxymap3d-menu > summary").click()


def _hue_shift(before, after, page, center_pc):
    """How much greener (green less blue, 0 to 255) the patch round a sector got between two screenshots."""
    first, second = _patch(before, page, center_pc), _patch(after, page, center_pc)
    return float((second[1] - second[2]) - (first[1] - first[2]))


# --- 3. The markers ---------------------------------------------------------------------------------------------------

def test_the_habitable_world_marker_finds_the_system_with_the_stored_count_and_goes_with_its_toggle(page, habitable_galaxy):
    _open_galaxy(page, habitable_galaxy.base_url, "")
    page.evaluate("([cell, g]) => document.querySelector('#galaxymap3d-canvas').galaxyOpenWindow(cell, g)",
                  [{"m": 1, "shell": 1, "band": 0, "wedge": 0}, 3])
    _settle(page)
    page.wait_for_timeout(500)
    button = page.locator('#galaxymap3d-kinds button[data-layer="marker-habitable"]')
    _open_menu(page)
    assert button.count() == 1 and button.get_attribute("aria-pressed") == "false"
    assert not [e for e in (page.evaluate(MARKERS) or []) if e["kind"] == "marker-habitable"]    # off until switched on

    button.click()
    assert "show=marker-habitable" in page.url
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-habitable" for e in entries or []))
    named = {e["name"] for e in found if e["kind"] == "marker-habitable"}
    assert named == {habitable_galaxy.habitable_system_name}, named
    assert habitable_galaxy.barren_system_name not in named

    # A click on it names the system and the stored count.
    _open_menu_closed(page)
    page.wait_for_timeout(300)
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    system = tuple(c + o for c, o in zip(HABITABLE_CENTER_PC, SYSTEM_OFFSET_PC))
    spot = _screen_of(page, system)
    assert spot
    page.mouse.move(*spot)
    page.mouse.move(spot[0] + 1, spot[1])
    page.wait_for_timeout(200)
    count = habitable_galaxy.habitable_worlds
    words = f"{habitable_galaxy.habitable_system_name}, {count} habitable world{'' if count == 1 else 's'}"
    assert words in page.locator("#galaxymap3d-tooltip").inner_text()
    assert "Habitable worlds" in page.locator("#galaxymap3d-layer-legend").inner_text()

    # Switched off again, the marker is gone and the address forgets it.
    _open_menu(page)
    button.click()
    assert button.get_attribute("aria-pressed") == "false" and "marker-habitable" not in page.url
    page.wait_for_timeout(400)
    assert not [e for e in (page.evaluate(MARKERS) or []) if e["kind"] == "marker-habitable"]


# --- 4. The help page -------------------------------------------------------------------------------------------------

def test_a_symbol_opens_the_help_page_at_its_section_and_the_page_lists_the_tiers_the_score_reads(page, habitable_galaxy):
    page.goto(f"{habitable_galaxy.base_url}/lists/planet", wait_until="load")
    page.wait_for_selector("svg.hab-symbol", state="attached")
    for anchor in ("overall", "domain-radiation"):
        link = page.locator(f'a.help-link[href="/help/habitability#{anchor}"]').first
        link.click()
        page.wait_for_url(f"**/help/habitability#{anchor}")
        target = page.locator(f"#{anchor}")
        assert target.count() == 1
        assert page.evaluate("(id) => document.querySelector(':target') === document.getElementById(id)", anchor)
        box = target.bounding_box()
        assert box and box["y"] >= 0, "the section sits under the sticky header"
        page.go_back()
        page.wait_for_selector("svg.hab-symbol", state="attached")

    # The levels table is the one the score reads: every cell of every row comes from the limits the score uses.
    page.goto(f"{habitable_galaxy.base_url}/help/habitability", wait_until="load")
    text = " ".join(page.evaluate("document.body.textContent").split())
    rows = [(domain, row) for domain, variables in habitability_ref.domains() for row in variables]
    assert len(rows) >= 10
    for domain, (variable, *cells) in rows:
        assert " ".join(variable.split()) in text, (domain, variable)
        for cell in cells:
            assert not cell or " ".join(cell.split()) in text, (domain, variable, cell)
    assert "Ultraviolet light at the ground" in text
