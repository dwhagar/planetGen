# tests/test_web_systems_values_browser.py

"""
UX.228: the Systems list's "Filter by what the system holds" in a real browser, on a real database (`tests/scored_site.py`): a gravity range typed into
the form keeps the systems that hold a planet in it (chips, the form filled again), a chip's cross takes the range away, and the table keeps the range while
it is sorted. Skipped without Playwright and Chromium.
"""

import re

import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.db import query  # noqa: E402
from tests.scored_site import scored_galaxy  # noqa: E402,F401 -- fixture
from tests.test_web_search_values_browser import _fill, browser, page  # noqa: E402,F401 -- fixtures and helper

pytestmark = pytest.mark.browser


def _rows(config):
    conn = query.open_readonly(config)
    try:
        return [dict(row) for row in conn.execute(
            "SELECT ss.name AS system_name, MAX(p.gravity_g) AS peak FROM star_systems ss JOIN planets p ON p.star_system_id = ss.id "
            "WHERE p.gravity_g IS NOT NULL GROUP BY ss.id, ss.name").fetchall()]
    finally:
        conn.close()


def _listed(page):
    """Every system name in the list: the table draws only the rows in view, so it is scrolled through, to its count."""
    rows = page.locator("#all-systems tbody tr[data-index]")
    rows.first.wait_for(state="attached", timeout=15000)
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(300)
    total = int(re.match(r"\d+", page.locator("#all-systems .datatable-count-line").inner_text()).group())
    names = set()
    for _ in range(200):
        names |= {rows.nth(i).locator("td").first.inner_text().strip() for i in range(rows.count())} - {"…"}
        if len(names) >= total:
            break
        page.locator("#all-systems .datatable-scroll").evaluate("el => { el.scrollTop += el.clientHeight - 40; }")
        page.wait_for_timeout(100)
    assert len(names) == total, (len(names), total)
    return names


def test_a_gravity_range_typed_in_the_systems_form_keeps_the_systems_with_a_planet_in_it(scored_galaxy, page):
    systems = _rows(scored_galaxy.config)
    assert len(systems) >= 3
    tops = sorted(row["peak"] for row in systems)
    cut = round(tops[len(tops) // 2] - 1e-4, 4)
    want = {row["system_name"] for row in systems if row["peak"] >= cut}
    assert want and len(want) < len(systems)
    assert page.goto(f"{scored_galaxy.base_url}/systems", wait_until="load").status == 200
    page.locator("details.search-values > summary").click()
    _fill(page, "planet_min_gravity", cut)
    page.get_by_role("button", name="Apply", exact=True).click()
    page.wait_for_url(re.compile(r"planet_min_gravity="))
    assert "_unit=" not in page.url and "planet_other_gas" not in page.url
    listed = _listed(page)
    assert {name for name in listed if name in {row["system_name"] for row in systems}} == want
    chips = page.locator(".active-filters .filter-chip").all_inner_texts()
    assert any(text.startswith(f"Planet surface gravity: ≥ {cut:g} g₀") for text in chips), chips
    assert page.locator('sl-input[name="planet_min_gravity"] input').input_value() == f"{cut:g}"
    # sorting keeps the range
    page.get_by_role("columnheader", name=re.compile("Stars")).get_by_role("link").first.click()
    page.wait_for_url(re.compile(r"systems_sort=stars"))
    assert f"planet_min_gravity={cut:g}" in page.url
    assert _listed(page) == listed
    page.locator(".active-filters .filter-chip", has_text="surface gravity").locator("a").click()
    page.wait_for_url(re.compile(r"^(?!.*planet_min_gravity).*/systems"))
    assert len(_listed(page)) > len(want)
