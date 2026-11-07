# tests/test_web_uncharted_sectors_browser.py

"""
MAP.399 in a real browser: the Sectors list tells generated sectors from uncharted placeholders (the rows a drifting rogue planet makes in a cell nothing was
generated in). A placeholder is named from its cell and never by its raw hex designation, says it is uncharted, shows its region and no density, and the header counts
the generated sectors with the uncharted number beside it; its page says it is uncharted and in which region. Needs Playwright with Chromium and the MySQL test server,
like test_web_a11y.py; skipped without them.

Run with: pytest tests/test_web_uncharted_sectors_browser.py
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.db import chart_ledger, store  # noqa: E402
from planetgen.galaxy.spherical import provisional_sector_designation  # noqa: E402
from planetgen.names import naming_key  # noqa: E402
from tests.test_web_a11y import (  # noqa: E402,F401 -- the browser and axe fixtures, and the checks
    SCHEMES, VIEWPORTS, _axe_blocking, _overflow, axe_source, browser,
)
from tests.test_web_settings_browser import _page, config, server  # noqa: E402,F401

pytestmark = pytest.mark.browser

KEY = "0123abcd"
CHARTED = [(9, 2, 5), (9, 2, 6)]
UNCHARTED = [(9, 2, 7), (10, 3, 4), (11, 1, 2)]


@pytest.fixture
def sectors(mysql_config, server):
    """Two charted sector rows and three placeholders, under a naming key."""
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            for address in CHARTED + UNCHARTED:
                sector_id = store.uncharted_sector_id(conn, address, 4.0)
                if address in CHARTED:
                    chart_ledger.record(conn, address, store.sector_run_id(conn), 1)
                    conn.execute("UPDATE sectors SET name = ? WHERE id = ?", (f"Charted {address[2]}", sector_id))
    finally:
        conn.close()
    naming_key.set_resolver(lambda: KEY)
    yield mysql_config
    naming_key.set_resolver(None)


def _rows(page):
    return page.evaluate("""() => [...document.querySelectorAll('#sectors table tbody tr')].map((tr) =>
        [...tr.querySelectorAll('td')].map((td) => td.textContent.trim()))""")


@pytest.mark.parametrize("scheme", SCHEMES)
@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_the_list_names_labels_and_places_the_placeholders_and_counts_them_apart(browser, axe_source, sectors, server, viewport, scheme):
    context, page, problems, base_url = _page(browser, server, viewport, scheme)
    try:
        assert page.goto(base_url + "/sectors?coverage=all&empty=1", wait_until="load").status == 200
        page.wait_for_selector("#sectors table tbody tr")
        badges = page.locator(".page-subhead .badge, .badges .badge").all_inner_texts()
        assert "2 sectors" in badges and "3 empty" in badges, badges     # nothing is detected in these three cells
        rows = _rows(page)
        names = {naming_key.uncharted_sector_name(address, KEY): address for address in UNCHARTED}
        shown = {row[0]: row for row in rows}
        assert set(names) <= set(shown), (sorted(names), sorted(shown))
        for name, address in names.items():
            cells = shown[name]
            assert cells[1] == "Uncharted" and cells[2] in ("–", "-", "None"), cells
            assert cells[3].endswith(" Region"), cells
            assert provisional_sector_designation(*address) not in "".join("".join(row) for row in rows)
        for row in rows:
            if row[0].startswith("Charted"):
                assert row[1] != "Uncharted" and row[2] not in ("–", "-", "None"), row
        problems += _overflow(page) + _axe_blocking(page, axe_source)
    finally:
        context.close()
    assert not problems, f"at {viewport}/{scheme}:\n  " + "\n  ".join(problems)


def test_the_page_of_a_placeholder_says_it_is_uncharted_and_in_which_region(browser, sectors, server):
    context, page, problems, base_url = _page(browser, server)
    try:
        name = naming_key.uncharted_sector_name(UNCHARTED[0], KEY)
        page.goto(base_url + "/sectors?coverage=all&empty=1", wait_until="load")
        page.wait_for_selector("#sectors table tbody tr")
        page.get_by_role("link", name=name).click()
        page.wait_for_selector("h1")
        assert name in page.locator("h1").inner_text()
        badges = page.locator(".badges .badge").all_inner_texts()
        assert "Uncharted" in badges and any(text.endswith(" Region") for text in badges), badges
        assert "uncharted sector" in page.locator("main").inner_text().lower()
    finally:
        context.close()
    assert not problems, "\n".join(problems)


@pytest.mark.parametrize("scheme", SCHEMES)
@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_the_page_of_a_placeholder_has_a_what_is_placed_here_panel_that_fits_and_passes_axe(browser, axe_source, sectors, server, viewport, scheme):
    """MAP.401: the panel is there on every placeholder's page, says the sector is uncharted, and fits the screen."""
    context, page, problems, base_url = _page(browser, server, viewport, scheme)
    try:
        name = naming_key.uncharted_sector_name(UNCHARTED[0], KEY)
        page.goto(base_url + "/sectors?coverage=all&empty=1", wait_until="load")
        page.wait_for_selector("#sectors table tbody tr")
        page.get_by_role("link", name=name).click()
        page.wait_for_selector("#uncharted-contents")
        panel = page.locator("#uncharted-contents")
        assert "Detected by long-range scans" in panel.inner_text()
        assert "long-range scans" in panel.inner_text()
        problems += _overflow(page) + _axe_blocking(page, axe_source)
    finally:
        context.close()
    assert not problems, f"at {viewport}/{scheme}:\n  " + "\n  ".join(problems)


def test_the_page_of_a_placeholder_shows_the_sector_map_locked_to_its_cell(browser, sectors, server):
    """MAP.404: the Sector Map is on the page, locked to the cell, and carries the link to the cell's scene."""
    context, page, problems, base_url = _page(browser, server)
    try:
        name = naming_key.uncharted_sector_name(UNCHARTED[0], KEY)
        page.goto(base_url + "/sectors?coverage=all&empty=1", wait_until="load")
        page.wait_for_selector("#sectors table tbody tr")
        page.get_by_role("link", name=name).click()
        page.wait_for_selector("#map canvas")
        assert "Sector Map" in page.locator("#map").inner_text()
        assert "/galaxy/uncharted/" in page.content()      # the map can open what the scans detected in the cell
    finally:
        context.close()
    assert not problems, "\n".join(problems)
