# tests/test_web_object_list_values_browser.py

"""
UX.228: the object lists' physical-value controls in a real browser, on a real database (`tests/scored_site.py`): a gravity range typed into a
planet list's form keeps the planets whose stored gravity lies in it (chips, a Values column, the form filled again), a chip's cross takes only
its range away, a column of a physical value sorts the list by it, and an order by a value alone sorts it. Skipped without Playwright and Chromium.
"""

import re

import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.db import query  # noqa: E402
from tests.scored_site import scored_galaxy  # noqa: E402,F401 -- fixture
from tests.test_web_search_values_browser import _fill, browser, page  # noqa: E402,F401 -- fixtures and helper

pytestmark = pytest.mark.browser


def _planets(config):
    conn = query.open_readonly(config)
    try:
        return {row["name"]: dict(row) for row in conn.execute(
            "SELECT name, gravity_g, surface_temperature_k FROM planets WHERE gravity_g IS NOT NULL "
            "AND surface_temperature_k IS NOT NULL").fetchall()}
    finally:
        conn.close()


def _rows(page):
    rows = page.locator("#object-list tbody tr[data-index]")
    rows.first.wait_for(state="attached", timeout=15000)
    return rows


def _listed(page):
    """Every name in the list: the table draws only the rows in view, so it is scrolled through, to its count (read again once the table has settled)."""
    rows = _rows(page)
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(300)
    total = int(re.match(r"\d+", page.locator("#object-list .datatable-count-line").inner_text()).group())
    names = set()
    for _ in range(300):
        names |= {rows.nth(i).locator("td").first.inner_text().strip() for i in range(rows.count())} - {"…"}      # "…" is a row still loading
        if len(names) >= total:
            break
        page.locator("#object-list .datatable-scroll").evaluate("el => { el.scrollTop += el.clientHeight - 40; }")
        page.wait_for_timeout(100)
    assert len(names) == total, (len(names), total)
    return names


def _first_names(page, count=5):
    rows = _rows(page)
    return [rows.nth(i).locator("td").first.inner_text().strip() for i in range(min(count, rows.count()))]


def test_a_gravity_range_typed_in_a_planet_lists_form_keeps_the_planets_in_it(scored_galaxy, page):
    planets = _planets(scored_galaxy.config)
    assert len(planets) > 3
    gravities = sorted(p["gravity_g"] for p in planets.values())
    least, most = round(gravities[1] - 1e-3, 4), round(gravities[-2] + 1e-3, 4)
    want = {name for name, p in planets.items() if least <= p["gravity_g"] <= most}
    assert want
    assert page.goto(f"{scored_galaxy.base_url}/lists/planet", wait_until="load").status == 200
    assert "values" not in page.locator("#object-list thead").inner_text().lower()          # no value asked for: no Values column
    page.locator("details.object-list-values > summary, .object-list-values summary").first.click()
    _fill(page, "planet_min_gravity", least)
    _fill(page, "planet_max_gravity", most)
    page.get_by_role("button", name="Apply").click()
    page.wait_for_url(re.compile(r"planet_max_gravity="))
    assert "_unit=" not in page.url
    listed = _listed(page)
    assert want <= listed and {name for name in listed if name in planets} == want
    assert "values" in page.locator("#object-list thead").inner_text().lower()
    assert "Surface gravity" in _rows(page).first.inner_text()
    chips = page.locator(".active-filters .filter-chip").all_inner_texts()
    assert any(text.startswith(f"surface gravity: {least:g}–{most:g} g₀") for text in chips), chips
    assert page.locator('sl-input[name="planet_min_gravity"] input').input_value() == f"{least:g}"
    page.locator(".active-filters .filter-chip", has_text="surface gravity").locator("a").click()
    page.wait_for_url(re.compile(r"^(?!.*planet_min_gravity).*/lists/planet"))
    assert len(_listed(page)) >= len(want)


def test_a_physical_value_column_sorts_the_list_and_survives_a_range(scored_galaxy, page):
    planets = _planets(scored_galaxy.config)
    page.goto(f"{scored_galaxy.base_url}/lists/planet?planet_min_gravity=0", wait_until="load")
    page.get_by_role("columnheader", name=re.compile("Surface temperature")).get_by_role("link").click()
    page.wait_for_url(re.compile(r"planet_sort=surface_temperature"))
    assert "planet_min_gravity=0" in page.url                 # the range stays while the table is sorted
    names = [name for name in _first_names(page, 8) if name in planets]
    temperatures = [planets[name]["surface_temperature_k"] for name in names]
    assert len(temperatures) >= 3 and temperatures == sorted(temperatures)


def test_an_order_by_a_value_alone_sorts_the_list_by_it(scored_galaxy, page):
    planets = _planets(scored_galaxy.config)
    assert page.goto(f"{scored_galaxy.base_url}/lists/planet?planet_value_sort=-gravity", wait_until="load").status == 200
    names = [name for name in _first_names(page, 8) if name in planets]
    gravities = [planets[name]["gravity_g"] for name in names]
    assert len(gravities) >= 3 and gravities == sorted(gravities, reverse=True)
    assert page.locator('select[name="planet_value_sort"]').input_value() == "-gravity"
