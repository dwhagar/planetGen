"""
The data table in a real browser (UX.41): the phenomena list scrolls through
every row while holding a few dozen in the page, sorts from its headers,
filters from its menus, keeps the address bar in step, and still works with
scripts off.

The page is the real Flask view over a fake API of 230 phenomena, in headless
Chromium, so it runs wherever Playwright and Chromium are installed (MySQL or
not). The pure logic is in `tests/js/datatablestate.test.mjs`; the server
side in `test_web_system_phen.py` and `test_api_phenomena_table.py`.
"""

import threading
import time
from urllib.parse import parse_qs, urlparse

import pytest
from werkzeug.serving import make_server

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.api.config import Config  # noqa: E402
from planetgen.api.limiter import PAGE_LIMITS_OFF  # noqa: E402
from planetgen.web.app import create_app  # noqa: E402
from planetgen.web.lib import apiclient  # noqa: E402

DB = "planetgen_table_fixture"
VIEWPORT = {"width": 1280, "height": 900}
TOTAL = 230
KINDS = [("nebula", "emission_nebula"), ("nebula", "dark_nebula"), ("rogue_planet", "terrestrial"),
         ("rogue_planet", "gas giant"), ("black_hole", "quiescent")]
ROWS = [
    {"id": i, "type": KINDS[i % len(KINDS)][0], "descriptor": KINDS[i % len(KINDS)][1],
     "name": f"Phenomenon {i:03d}", "radius_ly": 1.5 + i if KINDS[i % len(KINDS)][0] == "nebula" else None,
     "sector_id": None, "sector_name": None, "placed": bool(i % 2)}
    for i in range(TOTAL)
]


def _fake_phenomena(db, limit=None, offset=None, sort=None, descending=False, types=(), descriptors=(),
                    placed=None, facets=False):
    """What `apiclient.get_phenomena` answers from `ROWS`."""
    rows = [row for row in ROWS if (not types or row["type"] in types)
            and (not descriptors or row["descriptor"] in descriptors)]
    key = {"radius": "radius_ly", "sector": "sector_name"}.get(sort or "name", sort or "name")
    rows = sorted(rows, key=lambda row: (row[key] is None, row[key] if row[key] is not None else 0, row["name"]),
                  reverse=bool(descending))
    body = {"items": rows[offset:offset + limit], "total": len(rows), "limit": limit, "offset": offset}
    if facets:
        body["facets"] = {}
        for column, own in (("type", types), ("descriptor", descriptors)):
            pool = [row for row in ROWS if (column == "type" or not types or row["type"] in types)
                    and (column == "descriptor" or not descriptors or row["descriptor"] in descriptors)]
            body["facets"][column] = [
                {"value": value, "count": sum(1 for row in pool if row[column] == value)}
                for value in sorted({row[column] for row in pool})]
    return body


class _SiteConfig(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "table-fixture-secret"
    RATELIMIT_PAGES = PAGE_LIMITS_OFF


@pytest.fixture(scope="module")
def table_site():
    """The base URL of a local server for the Phenomena list, on fixture data."""
    patch = pytest.MonkeyPatch()
    try:
        patch.setattr(apiclient, "get_phenomena", _fake_phenomena)
        patch.setattr(apiclient, "auth_me", lambda cookie_header: None)
        patch.setattr(apiclient, "get_population_status", lambda db: dict(apiclient.POPULATION_NONE))
        app = create_app(_SiteConfig)
        app.testing = True
        server = make_server("127.0.0.1", 0, app, threaded=True)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            yield f"http://127.0.0.1:{server.server_port}"
        finally:
            server.shutdown()
            thread.join(timeout=10)
    finally:
        patch.undo()


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
    context = browser.new_context(viewport=VIEWPORT, reduced_motion="reduce", color_scheme="dark")
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.errors = errors
    try:
        yield page
    finally:
        context.close()
    assert not errors, f"uncaught exceptions: {errors}"


def _open(page, base_url, query=""):
    assert page.goto(f"{base_url}/phenomena{query}", wait_until="load").status == 200
    page.locator("[data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)


def _first_names(page, count=3):
    rows = page.locator("tbody tr[data-index]")
    rows.first.wait_for(state="attached")
    return [rows.nth(i).locator("td").first.inner_text() for i in range(min(count, rows.count()))]


def _count_line(page):
    return page.locator(".datatable-count-line").inner_text()


def _wait_count(page, text):
    page.locator(".datatable-count-line", has_text=text).wait_for(state="attached", timeout=10000)


def test_the_table_holds_a_window_of_rows_not_every_row(page, table_site):
    _open(page, table_site)
    assert _count_line(page) == f"{TOTAL} phenomena"
    assert _first_names(page, 2) == ["Phenomenon 000", "Phenomenon 001"]
    assert 5 < page.locator("tbody tr[data-index]").count() < 60
    # The pager is for scripts-off visitors.
    assert not page.locator(".datatable-pager .pagination").is_visible()
    assert page.locator("table").get_attribute("aria-rowcount") == str(TOTAL + 1)


def test_scrolling_brings_in_the_rest_of_the_rows(page, table_site):
    _open(page, table_site)
    scroller = page.locator(".datatable-scroll")
    scroller.evaluate("el => { el.scrollTop = el.scrollHeight; }")
    last = page.locator("tbody tr[data-index]", has_text=f"Phenomenon {TOTAL - 1:03d}")
    last.wait_for(state="attached", timeout=10000)
    assert "datatable-pending" not in (last.get_attribute("class") or "")
    assert page.locator("tbody tr[data-index]").count() < 60
    # Scrolled halfway: rows from the middle pages, fetched on demand.
    scroller.evaluate("el => { el.scrollTop = el.scrollHeight / 2; }")
    shown = "[...document.querySelectorAll('tbody tr[data-index]:not(.datatable-pending)')].map(tr => +tr.dataset.index)"
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        indexes = page.evaluate(shown)
        if indexes and 40 < min(indexes) and max(indexes) < TOTAL - 40:
            break
        page.wait_for_timeout(100)
    middle = min(indexes)
    assert 40 < middle < TOTAL - 40


def test_a_header_sorts_both_ways_and_the_address_follows(page, table_site):
    _open(page, table_site)
    header = page.locator("th[data-col='name']")
    assert header.get_attribute("aria-sort") == "ascending"
    page.locator("th[data-col='name'] .datatable-sort").click()
    page.locator("th[data-col='name'][aria-sort='descending']").wait_for(state="attached")
    page.locator("tbody tr[data-index='0']", has_text=f"Phenomenon {TOTAL - 1:03d}").wait_for(state="attached")
    assert parse_qs(urlparse(page.url).query) == {"sort": ["name"], "order": ["desc"]}

    page.locator("th[data-col='type'] .datatable-sort").click()
    page.locator("th[data-col='type'][aria-sort='ascending']").wait_for(state="attached")
    assert page.locator("th[data-col='name']").get_attribute("aria-sort") == "none"
    page.locator("tbody tr[data-index='0']", has_text="Black Hole").wait_for(state="attached")
    assert parse_qs(urlparse(page.url).query) == {"sort": ["type"]}

    # The address gives the same table on a reload (the server renders it).
    page.reload()
    page.locator("[data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert page.locator("th[data-col='type']").get_attribute("aria-sort") == "ascending"
    assert page.locator("tbody tr[data-index='0'] td").nth(1).inner_text() == "Black Hole"


def test_filters_narrow_the_rows_and_the_other_menu(page, table_site):
    _open(page, table_site)
    kind = page.locator(".datatable-facet[data-param='type']")
    kind.locator("summary").click()
    kind.get_by_label("Rogue Planet").check()
    rogue = sum(1 for row in ROWS if row["type"] == "rogue_planet")
    _wait_count(page, f"{rogue} phenomena match")
    assert parse_qs(urlparse(page.url).query) == {"type": ["rogue_planet"]}
    for name in _first_names(page, 5):
        assert int(name.split()[-1]) % 5 in (2, 3)
    # The Descriptor menu now counts only rogue planets and hides the rest.
    descriptor = page.locator(".datatable-facet[data-param='descriptor']")
    descriptor.locator("summary").click()
    assert descriptor.locator("li:not([hidden])").count() == 2
    assert descriptor.locator("li", has_text="Gas giant").inner_text().endswith(f"({rogue // 2})")

    descriptor.get_by_label("Gas giant").check()
    _wait_count(page, f"{rogue // 2} phenomena match")
    assert parse_qs(urlparse(page.url).query) == {"type": ["rogue_planet"], "descriptor": ["gas giant"]}

    page.locator(".datatable-clear").click()
    _wait_count(page, f"{TOTAL} phenomena")
    assert "match" not in _count_line(page)
    assert urlparse(page.url).query == ""
    assert page.locator(".datatable-clear").is_hidden()
    assert descriptor.locator("li:not([hidden])").count() == len(KINDS)


def test_a_filtered_address_loads_filtered(page, table_site):
    _open(page, table_site, "?type=black_hole&sort=name&order=desc")
    holes = sum(1 for row in ROWS if row["type"] == "black_hole")
    assert _count_line(page) == f"{holes} phenomena match"
    assert page.locator(".datatable-facet[data-param='type'] input:checked").count() == 1
    assert _first_names(page, 1)[0] == "Phenomenon 229"


def test_the_table_scrolls_from_the_keyboard(page, table_site):
    _open(page, table_site)
    scroller = page.locator(".datatable-scroll")
    assert scroller.get_attribute("tabindex") == "0"
    scroller.focus()
    for _ in range(3):
        page.keyboard.press("PageDown")
    page.wait_for_timeout(300)
    assert scroller.evaluate("el => el.scrollTop") > 0


def test_without_scripts_it_is_a_sortable_filterable_paged_table(browser, table_site):
    context = browser.new_context(viewport=VIEWPORT, java_script_enabled=False)
    page = context.new_page()
    try:
        page.goto(f"{table_site}/phenomena", wait_until="load")
        assert page.locator("tbody tr").count() == 50
        assert page.locator(".datatable-pager .pagination").is_visible()
        page.get_by_role("link", name="Next").click()
        assert page.locator("tbody tr").first.locator("td").first.inner_text() == "Phenomenon 050"
        page.get_by_role("link", name="Name").click()
        assert page.locator("th[data-col='name']").get_attribute("aria-sort") == "descending"
        page.locator(".datatable-facet[data-param='type'] summary").click()
        page.locator(".datatable-facet[data-param='type']").get_by_label("Black Hole").check()
        page.get_by_role("button", name="Apply filters").click()
        assert page.locator("tbody tr").count() == 46
        assert parse_qs(urlparse(page.url).query)["type"] == ["black_hole"]
    finally:
        context.close()
