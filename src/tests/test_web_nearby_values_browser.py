# tests/test_web_nearby_values_browser.py

"""
UX.228: What's nearby's physical-value controls in a real browser, on a real database (`tests/scored_site.py`): a gravity range typed into the
form lists the planets in the sphere whose stored gravity lies in it (chips, a Values column, the form filled again), a chip's cross takes the
range away, and an order by a value lists that kind by it. Skipped without Playwright and Chromium.
"""

import json
import re

import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.db import query, store  # noqa: E402
from planetgen.web.lib.pagination import PAGE_SIZE  # noqa: E402
from tests.scored_site import scored_galaxy  # noqa: E402,F401 -- fixture
from tests.test_web_search_values_browser import _fill, browser, page  # noqa: E402,F401 -- fixtures and helper

pytestmark = pytest.mark.browser

DISTANCE = 50


@pytest.fixture(scope="module")
def placed(scored_galaxy):
    """The sample galaxy with its systems given places in the galaxy frame, a parsec or two apart (the fixture leaves them unplaced)."""
    conn = store.get_connection(scored_galaxy.config)
    try:
        conn.execute("UPDATE sectors SET center_x_pc = 2000.0, center_y_pc = 100.0, center_z_pc = 5.0, galactic_radius_pc = 2002.5 WHERE name = 'Scored'")
        for index, row in enumerate(conn.execute("SELECT id FROM star_systems ORDER BY id").fetchall()):
            conn.execute("UPDATE star_systems SET position_x_mpc = 0, position_y_mpc = 0, position_z_mpc = 0, "
                         "galaxy_x_pc = ?, galaxy_y_pc = 100.0, galaxy_z_pc = 5.0 WHERE id = ?", (2000.0 + 1.5 * index, row["id"]))
        conn.commit()
    finally:
        conn.close()
    return scored_galaxy


def _planets(config):
    conn = query.open_readonly(config)
    try:
        return {row["name"]: dict(row) for row in conn.execute("SELECT name, gravity_g FROM planets WHERE gravity_g IS NOT NULL").fetchall()}
    finally:
        conn.close()


def _place(site, page):
    systems = json.loads(page.request.get(f"{site.base_url}/api/systems?limit=1").text())["items"]
    return f"system:{systems[0]['id']}"


def _names(page):
    """The names in the results table (a plain table, not a virtual one)."""
    rows = page.locator("#nearby-results tbody tr")
    return [rows.nth(i).locator("td").nth(1).inner_text().strip() for i in range(rows.count())]


def test_a_gravity_range_typed_in_the_nearby_form_lists_the_planets_in_it(placed, page):
    scored_galaxy = placed
    planets = _planets(scored_galaxy.config)
    place = _place(scored_galaxy, page)
    sphere = json.loads(page.request.get(f"{scored_galaxy.base_url}/api/near?from={place}&distance={DISTANCE}&kinds=planet&limit=200").text())
    assert "rows" in sphere, sphere
    inside = {row["name"] for row in sphere["rows"]}
    gravities = sorted(planets[name]["gravity_g"] for name in inside if name in planets)
    assert len(gravities) > 8, "the sample galaxy has too few planets in the sphere"
    least, most = round(gravities[len(gravities) // 4] - 1e-4, 4), round(gravities[3 * len(gravities) // 4] + 1e-4, 4)
    want = {name for name in inside if least <= planets[name]["gravity_g"] <= most}
    assert want and want != inside and len(want) <= PAGE_SIZE
    assert page.goto(f"{scored_galaxy.base_url}/nearby?place={place}&distance={DISTANCE}&kinds=planet", wait_until="load").status == 200
    page.locator("details.search-values > summary").click()
    _fill(page, "planet_min_gravity", least)
    _fill(page, "planet_max_gravity", most)
    page.get_by_role("button", name="Search").click()
    page.wait_for_url(re.compile(r"planet_max_gravity="))
    assert "_unit=" not in page.url and "planet_other_gas" not in page.url          # a short address
    assert set(_names(page)) == want
    assert "values" in page.locator("#nearby-results thead").inner_text().lower()
    assert "Surface gravity" in page.locator("#nearby-results tbody tr").first.inner_text()
    chips = page.locator(".active-filters .filter-chip").all_inner_texts()
    assert any(text.startswith(f"Planet surface gravity: {least:g}–{most:g} g₀") for text in chips), chips
    assert page.locator('sl-input[name="planet_min_gravity"] input').input_value() == f"{least:g}"
    page.locator(".active-filters .filter-chip", has_text="surface gravity").locator("a").click()
    page.wait_for_url(re.compile(r"^(?!.*planet_min_gravity).*/nearby"))
    names = _names(page)
    assert set(names) <= inside and len(names) == min(PAGE_SIZE, len(inside)) and len(names) > len(want)       # the page of the whole sphere again


def test_an_order_by_a_value_lists_the_planets_by_it(placed, page):
    scored_galaxy = placed
    planets = _planets(scored_galaxy.config)
    place = _place(scored_galaxy, page)
    assert page.goto(f"{scored_galaxy.base_url}/nearby?place={place}&distance={DISTANCE}&planet_value_sort=-gravity", wait_until="load").status == 200
    names = _names(page)
    gravities = [planets[name]["gravity_g"] for name in names]
    assert len(gravities) > 3 and gravities == sorted(gravities, reverse=True)
    assert page.locator('select[name="planet_value_sort"]').input_value() == "-gravity"
    assert any("ordered by surface gravity, high to low" in text for text in page.locator(".active-filters .filter-chip").all_inner_texts())
