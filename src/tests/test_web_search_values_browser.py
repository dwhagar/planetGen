# tests/test_web_search_values_browser.py

"""
UX.227: the search page's physical-value controls in a real browser, on a real database (`tests/scored_site.py`): a gravity and a dose
range typed into the form give the planets whose stored gravity and dose lie in them, as the chips and the Values column say, and the
unit menu turns the same pressure into the same planets. Skipped without Playwright and Chromium.
"""

import re

import pytest

pytest.importorskip("playwright.sync_api")
from playwright import sync_api  # noqa: E402

from planetgen.db import query  # noqa: E402
from tests.scored_site import scored_galaxy  # noqa: E402,F401 -- fixture

pytestmark = pytest.mark.browser


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
        playwright.stop()


@pytest.fixture()
def page(browser):
    context = browser.new_context(viewport={"width": 1280, "height": 900}, reduced_motion="reduce")
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.errors = errors
    try:
        yield page
    finally:
        context.close()
    assert not errors, f"uncaught exceptions: {errors}"


def _planets(config):
    conn = query.open_readonly(config)
    try:
        return [dict(row) for row in conn.execute(
            "SELECT name, gravity_g, surface_dose_msv_yr, atmospheric_pressure_pa FROM planets "
            "WHERE gravity_g IS NOT NULL AND surface_dose_msv_yr IS NOT NULL").fetchall()]
    finally:
        conn.close()


def _fill(page, name, value):
    page.locator(f'sl-input[name="{name}"] input').fill(str(value))


def _listed(page, panel="planets"):
    """The names in a result panel: the table draws only the rows in view, so the panel is scrolled through, to its count."""
    rows = page.locator(f"#search-{panel} tbody tr[data-index]")
    rows.first.wait_for(state="attached", timeout=15000)
    total = int(re.match(r"\d+", page.locator(f"#search-{panel} .datatable-count-line").inner_text()).group())
    names = set()
    for _ in range(200):
        names |= {rows.nth(i).locator("td").first.inner_text().strip() for i in range(rows.count())}
        if len(names) >= total:
            break
        page.locator(f"#search-{panel} .datatable-scroll").evaluate("el => { el.scrollTop += el.clientHeight - 40; }")
        page.wait_for_timeout(100)
    assert len(names) == total, (len(names), total)
    return names


def test_a_gravity_and_a_dose_range_typed_in_the_form_find_the_right_planets(scored_galaxy, page):
    planets = _planets(scored_galaxy.config)
    assert len(planets) > 3
    gravities = sorted(p["gravity_g"] for p in planets)
    doses = sorted(p["surface_dose_msv_yr"] for p in planets)
    least, most, limit = round(gravities[1] - 1e-3, 4), round(gravities[-2] + 1e-3, 4), round(doses[len(doses) // 2] + 1e-3, 4)
    want = {p["name"] for p in planets if least <= p["gravity_g"] <= most and p["surface_dose_msv_yr"] <= limit}
    assert want
    assert page.goto(f"{scored_galaxy.base_url}/search", wait_until="load").status == 200
    page.locator("details.search-values > summary").click()
    _fill(page, "planet_min_gravity", least)
    _fill(page, "planet_max_gravity", most)
    _fill(page, "planet_max_surface_dose", limit)
    page.get_by_role("button", name="Search").click()
    page.wait_for_url(re.compile(r"planet_max_surface_dose="))
    assert f"planet_min_gravity={least}" in page.url and f"planet_max_gravity={most}" in page.url
    assert "_unit=" not in page.url and "planet_other_gas" not in page.url          # a unit that says nothing is not in the address
    assert _listed(page) == want
    chips = page.locator(".active-filters .filter-chip").all_inner_texts()
    assert any(text.startswith(f"Planet surface gravity: {least:g}–{most:g} g₀") for text in chips), chips
    assert any(text.startswith(f"Planet surface radiation dose: ≤ {limit:g} mSv/yr") for text in chips), chips
    cell = page.locator("#search-planets tbody tr[data-index]").first.inner_text()
    assert "Surface gravity" in cell and "Surface radiation dose" in cell
    # the form comes back filled, and a chip's cross takes only its range away
    assert page.locator('sl-input[name="planet_min_gravity"] input').input_value() == f"{least:g}"
    page.locator(".active-filters .filter-chip", has_text="Planet surface gravity").locator("a").click()
    page.wait_for_url(re.compile(r"^(?!.*planet_min_gravity).*planet_max_surface_dose="))
    wider = {p["name"] for p in planets if p["surface_dose_msv_yr"] <= limit}
    assert _listed(page) == wider


def test_the_unit_menu_gives_the_same_planets_for_the_same_pressure(scored_galaxy, page):
    conn = query.open_readonly(scored_galaxy.config)
    try:
        rows = [dict(r) for r in conn.execute("SELECT name, atmospheric_pressure_pa AS pa FROM planets WHERE atmospheric_pressure_pa IS NOT NULL").fetchall()]
    finally:
        conn.close()
    pressures = sorted(r["pa"] for r in rows)
    assert len(pressures) > 3
    middle = pressures[len(pressures) // 2]
    bound = round(middle / 1000.0 * (1 - 1e-6), 6)      # kPa
    want = {r["name"] for r in rows if r["pa"] >= bound * 1000.0}
    assert want
    assert page.goto(f"{scored_galaxy.base_url}/search", wait_until="load").status == 200
    page.locator("details.search-values > summary").click()
    _fill(page, "planet_min_surface_pressure", bound)
    page.select_option('select[name="planet_surface_pressure_unit"]', "kPa")
    page.get_by_role("button", name="Search").click()
    page.wait_for_url(re.compile(r"planet_surface_pressure_unit=kPa"))
    assert _listed(page) == want
    assert any("kPa" in text for text in page.locator(".active-filters .filter-chip").all_inner_texts())


def test_an_order_by_a_value_sorts_the_planets_by_it(scored_galaxy, page):
    conn = query.open_readonly(scored_galaxy.config)
    try:
        planets = {r["name"]: r["gravity_g"] for r in conn.execute("SELECT name, gravity_g FROM planets WHERE gravity_g IS NOT NULL").fetchall()}
    finally:
        conn.close()
    assert page.goto(f"{scored_galaxy.base_url}/search?planet_min_gravity=0&planet_value_sort=-gravity", wait_until="load").status == 200
    rows = page.locator("#search-planets tbody tr[data-index]")
    rows.first.wait_for(state="attached", timeout=15000)
    shown = [planets[rows.nth(i).locator("td").first.inner_text().strip()] for i in range(rows.count())]
    assert shown == sorted(shown, reverse=True) and len(shown) > 3      # the first rows in view are the largest
    assert page.locator('select[name="planet_value_sort"]').input_value() == "-gravity"
