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
from tests.served import make_quiet_server

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.api.config import Config  # noqa: E402
from planetgen.api.limiter import PAGE_LIMITS_OFF  # noqa: E402
from planetgen.db import query as queryDb  # noqa: E402
from planetgen.web import sector_page  # noqa: E402
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
     "sector_id": None, "sector_name": None, "placed": bool(i % 2), "scattered": False}
    for i in range(TOTAL)
]


def _fake_phenomena(db, limit=None, offset=None, sort=None, descending=False, types=(), descriptors=(),
                    placed=None, facets=False, after=None, start=None):
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


UNCHARTED = [{"id": 5000 + i, "star_type": "G2V", "population": "old", "luminosity_sol": 1.0 + i, "temperature_k": 3000.0 + 10 * i,
              "age_gy": 1.0, "designation": "7.2.3", "shell_index": 7, "band_index": 2, "wedge_index": 3, "x": 1.0,
              "y": 2.0, "z": 3.0, "local_x": 0.1, "local_y": 0.2, "local_z": 0.3, "member_count": 0,
              "system_mass_solar": None} for i in range(130)]


def _fake_uncharted(db, limit=50, offset=0, sort=None, descending=False, after=None, start=None):
    """What `apiclient.get_uncharted_systems` answers from `UNCHARTED`, for the temperature sort and its jump box
    (PERF.135): a start past the last star answers the last page, as the API does."""
    rows = sorted(UNCHARTED, key=lambda row: row["temperature_k"], reverse=bool(descending))
    body = {"limit": limit, "offset": offset, "past_end": False}
    if start:
        kept = [row for row in rows if (row["temperature_k"] <= float(start) if descending
                                        else row["temperature_k"] >= float(start))]
        if kept:
            rows = kept
        else:
            rows, offset, body["past_end"] = rows[-limit:], 0, True
    body.update(items=rows[offset:offset + limit], total=len(rows))
    return body


SECTORS = [
    {"id": i, "name": f"Sector {i:03d}", "system_count": i % 7, "edge_ly": 10.0, "placed": bool(i % 3),
     "center_x_pc": 100.0 if i % 2 else -100.0, "center_y_pc": 50.0, "galactic_radius_ly": 100.0 + i}
    for i in range(130)
]
STANDALONE = [{"id": 1000 + i, "name": f"Drifting {i:03d}", "is_binary": i % 2, "star_summary": "G2V",
               "sector_id": None, "sector_name": None, "quadrant": None} for i in range(130)]


def _fake_sectors(db, limit=None, offset=None, sort=None, descending=False, quadrants=(), facets=False):
    rows = [row for row in SECTORS if not quadrants or ("I" if row["center_x_pc"] > 0 else "II") in quadrants
            or ("unplaced" in quadrants and not row["placed"])]
    key = {"name": "name", "systems": "system_count", "distance": "galactic_radius_ly"}.get(sort or "distance", "name")
    rows = sorted(rows, key=lambda row: (row[key], row["name"]), reverse=bool(descending))
    body = {"items": rows[offset:offset + limit], "total": len(rows), "limit": limit, "offset": offset}
    if facets:
        body["facets"] = {"quadrant": [{"value": "I", "count": 65}, {"value": "II", "count": 65}]}
    return body


def _fake_systems(db, star_type=None, sector_id=None, limit=None, offset=None, sort=None, descending=False,
                  binary=None, placement=None, octants=(), facets=False, after=None, start=None):
    rows = [row for row in STANDALONE if binary is None or bool(row["is_binary"]) == binary]
    key = {"binary": "is_binary"}.get(sort or "name", "name")
    rows = sorted(rows, key=lambda row: (row[key], row["name"]), reverse=bool(descending))
    if start and (sort or "name") == "name":  # PERF.75: the jump-by-value box, as the real reader applies it
        rows = [row for row in rows if (row["name"].casefold() < queryDb._after_prefix(start.casefold())
                                        if descending else row["name"].casefold() >= start.casefold())]
    if after is not None:  # PERF.75: the page that starts after the row with this id (the real reader keys on the sort)
        ids = [row["id"] for row in rows]
        offset = ids.index(int(after)) + 1 if int(after) in ids else offset
    body = {"items": rows[offset:offset + limit], "total": len(rows), "limit": limit, "offset": offset}
    if facets:
        body["facets"] = {"placement": [{"value": "standalone", "count": len(STANDALONE)}],
                          "binary": [{"value": "yes", "count": 65}, {"value": "no", "count": 65}], "octant": []}
    return body


POLITY = {"id": 2, "name": "Velar Concord", "government": "Concord", "color": "#3366cc", "reach_ly": 80.0,
          "species_id": 3, "species_name": "Velar", "era": "ancient", "civilization_age_years": 3.4e6,
          "capital_system_id": 5, "capital_name": "Kepler", "system_count": 120}
OWNED = [{"id": 2000 + i, "name": f"Holding {i:03d}", "distance_ly": float(i)} for i in range(120)]


def _fake_polity(db, polity_id, limit=None, offset=None, sort=None, descending=False):
    if polity_id != POLITY["id"]:
        raise apiclient.NotFoundError("no such polity")
    rows = sorted(OWNED, key=lambda row: row["name" if sort == "name" else "distance_ly"], reverse=bool(descending))
    return {**POLITY, "systems": rows[offset:offset + limit], "limit": limit, "offset": offset}


def _contents_system(i):
    return {
        "id": 3000 + i, "name": f"Star {i:03d}", "is_binary": 0, "binary_type": None,
        "stars": [{"star_type": "G2V", "temperature_k": 5772.0, "radius_km": 696000.0, "luminosity_w": 3.8e26}],
        "quadrant": "+x+y+z", "location": f"Fixture Sector -- nearest: Star {(i + 1) % 130:03d} (1.0 ly)",
        "center_distance_ly": float(i), "position_x_mpc": 1.0, "position_y_mpc": 0.5, "position_z_mpc": -0.5,
    }


_ROGUE = {"type": "rogue_planet", "descriptor": "jupiter", "class": "J", "radius_ly": 0.0, "octant": "+x+y+z",
          "offset_x_ly": 0.1, "offset_y_ly": 0.0, "offset_z_ly": 0.0,
          "nearest": [{"id": 3001, "name": "Star 001", "distance_ly": 0.5}]}
SECTOR = {
    "id": 5, "name": "Fixture Sector", "edge_mpc": 3066.0, "edge_ly": 10.0, "shell_index": 5, "band_index": 1,
    "wedge_index": 20, "placed": True, "center_x_pc": 500.0, "center_y_pc": 200.0, "center_z_pc": 10.0,
    "wiki_url": None, "systems": [_contents_system(i) for i in range(130)],
    "phenomena": [
        {**_ROGUE, "id": 21, "name": "Drifter", "distance_ly": 0.1},
        {**_ROGUE, "id": 22, "name": "Wanderer", "distance_ly": 0.2},
        {"id": 3, "type": "nebula", "name": "Veil", "descriptor": "emission", "radius_ly": 2.5, "distance_ly": 4.0,
         "offset_x_ly": 1.0, "offset_y_ly": 0.0, "offset_z_ly": 0.0},
    ],
    "neighbors": [],
}


def _fake_sector(db, sector_id):
    if int(sector_id) != SECTOR["id"]:
        raise apiclient.NotFoundError("no such sector")
    return SECTOR


def _fake_search(db, texts, tags, sizes=None, limit=None, offsets=None, panels=None):
    offset = (offsets or {}).get("sectors", 0)
    rows = [{"id": i, "name": f"Found {i:03d}", "edge_mpc": 3.07} for i in range(120)]
    panel = {"rows": rows[offset:offset + limit], "total": 120, "total_capped": False, "limit": limit,
             "offset": offset, "truncated": False}
    results = {name: None for name in ("sectors", "systems", "stars", "planets", "moons", "belts", "phenomena")}
    results["sectors"] = panel
    return {"facets": {}, "autocomplete": {"sectors": [], "systems": [], "stars": [], "planets": [], "moons": []},
            "facet_labels": {}, "results": results}


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
        patch.setattr(apiclient, "get_sectors", _fake_sectors)
        patch.setattr(apiclient, "get_systems", _fake_systems)
        patch.setattr(apiclient, "get_uncharted_systems", _fake_uncharted)
        patch.setattr(apiclient, "auth_me", lambda cookie_header: None)
        patch.setattr(apiclient, "get_population_status",
                      lambda db: {"generated": True, "species": True, "polities": True, "territories": False})
        patch.setattr(apiclient, "get_polity", _fake_polity)
        patch.setattr(apiclient, "get_search", _fake_search)
        patch.setattr(apiclient, "get_sector", _fake_sector)
        patch.setattr(apiclient, "get_sector_facilities", lambda db, sector_id: [])
        patch.setattr(apiclient, "get_galaxy_shape", lambda db: None)
        patch.setattr(sector_page, "fetch_tiles", lambda db, keys, known_stamp=None: {
            "stamp": "0" * 16, "tiles": {}, "density": {}, "edge_pc": 3.066, "has_shape": False})
        app = create_app(_SiteConfig)
        app.testing = True
        server = make_quiet_server("127.0.0.1", 0, app, threaded=True)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            yield f"http://127.0.0.1:{server.server_port}"
        finally:
            server.shutdown()
            thread.join(timeout=10)
    finally:
        patch.undo()


ADMIN = {"username": "boss", "must_change_credentials": False}
KEYS = [{"id": i, "label": f"key {i:03d}", "created_at": f"2026-01-{1 + i % 28:02d} 00:00:00", "last_used_at": None,
         "revoked_at": "2026-02-01 00:00:00" if i % 5 == 0 else None} for i in range(1, 121)]
JOBS = [{"id": f"20260101-000000-{i:08d}", "kind": "galaxy", "title": f"Job {i:03d}", "status": "done",
         "state": "done", "control": None, "live": False, "started_at": None, "created_at": "2026-01-01 00:00:00",
         "finished_at": None, "seconds": 5.0, "age_seconds": 5.0, "workers": None, "holder": None,
         "web_job_id": None, "children": [], "tasks": [],
         "totals": {"tasks": 4, "queued": 0, "running": 0, "done": 4, "failed": 0, "cancelled": 0,
                    "work_seconds": 0, "eta_seconds": None}} for i in range(120)]


def _fake_work(cookie_header, limit=None, offset=None):
    return {"available": True, "status": {"load": {"text": "0 / 0 / 0"}}, "items": JOBS[offset:offset + limit],
            "total": len(JOBS)}


@pytest.fixture(scope="module")
def admin_site():
    """A local server whose visitor is a signed-in admin, over fixture API keys and jobs."""
    patch = pytest.MonkeyPatch()
    try:
        patch.setattr(apiclient, "auth_me", lambda cookie_header: ADMIN)
        patch.setattr(apiclient, "auth_list_api_keys", lambda cookie_header: KEYS)
        patch.setattr(apiclient, "admin_work", _fake_work)
        patch.setattr(apiclient, "get_population_status",
                      lambda db: {"generated": True, "species": True, "polities": True, "territories": False})
        app = create_app(_SiteConfig)
        app.testing = True
        server = make_quiet_server("127.0.0.1", 0, app, threaded=True)
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
    # The Descriptor menu now counts only rogue planets and hides the rest (the counts follow the rows, UX.95).
    descriptor = page.locator(".datatable-facet[data-param='descriptor']")
    descriptor.locator("summary").click()
    descriptor.locator("li", has_text="Gas giant").locator(".datatable-count", has_text=f"({rogue // 2})") \
        .wait_for(state="attached", timeout=10000)
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
    deadline = time.monotonic() + 10
    while descriptor.locator("li:not([hidden])").count() != len(KINDS) and time.monotonic() < deadline:
        page.wait_for_timeout(100)
    assert descriptor.locator("li:not([hidden])").count() == len(KINDS)


def test_the_rows_show_before_the_filter_counts_arrive(page, table_site):
    """UX.95: a slow count request holds the menus ("...") but not the rows."""
    held = []
    page.route("**/table/**facets=1**", lambda route: held.append(route))   # kept, not answered yet
    _open(page, table_site)
    kind = page.locator(".datatable-facet[data-param='type']")
    kind.locator("summary").click()
    kind.get_by_label("Rogue Planet").check()
    rogue = sum(1 for row in ROWS if row["type"] == "rogue_planet")
    _wait_count(page, f"{rogue} phenomena match")                       # the rows and the count line are in
    deadline = time.monotonic() + 10
    while not held and time.monotonic() < deadline:
        page.wait_for_timeout(50)
    assert held, "the counts were asked for in a request of their own"
    assert kind.locator(".datatable-count").first.inner_text() == "..."   # the counts are not
    for route in held:
        route.continue_()
    kind.locator(".datatable-count", has_text="(").first.wait_for(state="attached", timeout=10000)


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


def test_a_table_keeps_its_sort_and_filter_across_a_reload(page, table_site):
    assert page.goto(f"{table_site}/systems", wait_until="load").status == 200
    page.locator("#systems-table[data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert page.locator("#systems-table .datatable-count-line").inner_text() == "130 systems"

    page.locator("#systems-table th[data-col='binary'] .datatable-sort").click()
    page.locator("#systems-table th[data-col='binary'][aria-sort='ascending']").wait_for(state="attached")
    assert parse_qs(urlparse(page.url).query) == {"systems_sort": ["binary"]}

    facet = page.locator("#systems-table .datatable-facet", has_text="Stars")
    facet.locator("summary").click()
    facet.get_by_label("Binary").check()
    page.locator("#systems-table .datatable-count-line", has_text="65 systems match") \
        .wait_for(state="attached", timeout=10000)
    query = parse_qs(urlparse(page.url).query)
    assert query["systems_sort"] == ["binary"] and query["binary"] == ["yes"]

    # A reload renders the table the way it was left.
    page.reload()
    page.locator("#systems-table[data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert page.locator("#systems-table th[data-col='binary']").get_attribute("aria-sort") == "ascending"
    assert page.locator("#systems-table .datatable-count-line").inner_text() == "65 systems match"


def test_the_systems_jump_box_begins_the_list_at_a_typed_name(page, table_site):
    """PERF.75: "Starting at" begins the Systems list at the first name from the text on, with the count of what is
    left, in the address, and it goes away under a sort that has no box and comes back empty."""
    assert page.goto(f"{table_site}/systems", wait_until="load").status == 200
    page.locator("#systems-table[data-enhanced='true']").wait_for(state="attached", timeout=15000)
    box = page.locator("#systems-table .datatable-jump input[type='text']")
    assert box.is_visible()
    box.fill("Drifting 100")
    page.locator("#systems-table .datatable-jump button").click()
    page.locator("#systems-table .datatable-count-line", has_text="30 systems from \u201cDrifting 100\u201d") \
        .wait_for(state="attached")
    page.locator("#systems-table tbody tr[data-index='0']", has_text="Drifting 100").wait_for(state="attached")
    assert parse_qs(urlparse(page.url).query) == {"systems_from": ["Drifting 100"]}

    # A sort with no box drops the start and hides the box; the name sort again begins at the top.
    page.locator("#systems-table th[data-col='binary'] .datatable-sort").click()
    page.locator("#systems-table th[data-col='binary'][aria-sort='ascending']").wait_for(state="attached")
    page.locator("#systems-table .datatable-count-line", has_text="130 systems").wait_for(state="attached")
    assert not page.locator("#systems-table .datatable-jump").is_visible()
    assert parse_qs(urlparse(page.url).query) == {"systems_sort": ["binary"]}

    # The address gives the same list on a reload, and the link beside the box goes back to the top.
    page.goto(f"{table_site}/systems?systems_from=Drifting+100", wait_until="load")
    page.locator("#systems-table[data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert page.locator("#systems-table .datatable-jump input[type='text']").input_value() == "Drifting 100"
    page.locator("#systems-table tbody tr[data-index='0']", has_text="Drifting 100").wait_for(state="attached")
    page.locator("#systems-table .datatable-jump-clear").click()
    page.locator("#systems-table .datatable-count-line", has_text="130 systems").wait_for(state="attached")
    assert parse_qs(urlparse(page.url).query) == {}


def test_the_uncharted_jump_box_begins_the_list_at_a_typed_value_and_says_when_it_lay_past_the_end(page, table_site):
    """PERF.135: the Uncharted stars table's "Starting at" box takes a temperature under the temperature sort; a value
    past the last star shows the last page and says so, and the line goes when the start is cleared."""
    assert page.goto(f"{table_site}/systems?uncharted=1&uncharted_sort=temperature", wait_until="load").status == 200
    page.locator("#uncharted_systems-table[data-enhanced='true']").wait_for(state="attached", timeout=15000)
    box = page.locator("#uncharted_systems-table .datatable-jump input[type='text']")
    assert box.is_visible() and "Kelvin" in box.get_attribute("placeholder")
    box.fill("4000")
    page.locator("#uncharted_systems-table .datatable-jump button").click()
    page.locator("#uncharted_systems-table .datatable-count-line", has_text="30 uncharted stars from \u201c4000\u201d") \
        .wait_for(state="attached")
    page.locator("#uncharted_systems-table tbody tr[data-index='0']", has_text="star 5100").wait_for(state="attached")
    assert parse_qs(urlparse(page.url).query)["uncharted_from"] == ["4000"]

    box.fill("99999")
    page.locator("#uncharted_systems-table .datatable-jump button").click()
    line = page.locator("#uncharted_systems-table .datatable-count-line", has_text="Nothing from \u201c99999\u201d on")
    line.wait_for(state="attached")
    assert "showing the last 50" in line.inner_text() and "50 uncharted stars" in line.inner_text()
    page.locator("#uncharted_systems-table tbody tr[data-index='0']", has_text="star 5080").wait_for(state="attached")

    page.locator("#uncharted_systems-table .datatable-jump-clear").click()
    page.locator("#uncharted_systems-table .datatable-count-line", has_text="130 uncharted stars").wait_for(state="attached")
    assert "Nothing" not in page.locator("#uncharted_systems-table .datatable-count-line").inner_text()


def test_without_scripts_the_jump_box_is_a_plain_form(browser, table_site):
    context = browser.new_context(java_script_enabled=False, viewport=VIEWPORT)
    try:
        page = context.new_page()
        assert page.goto(f"{table_site}/systems?systems_sort=name&systems_order=desc", wait_until="load").status == 200
        page.locator("#systems-table .datatable-jump input[type='text']").fill("Drifting 05")
        page.locator("#systems-table .datatable-jump button").click()
        page.wait_for_load_state("load")
        query = parse_qs(urlparse(page.url).query)
        assert query["systems_from"] == ["Drifting 05"] and query["systems_order"] == ["desc"]
        assert "from \u201cDrifting 05\u201d" in page.locator("#systems-table .datatable-count-line").inner_text()
        assert page.locator("#systems-table tbody tr td a").first.inner_text() == "Drifting 059"
    finally:
        context.close()


def test_a_politys_systems_scroll_and_sort_through_the_tables_own_route(page, table_site):
    assert page.goto(f"{table_site}/polities/2", wait_until="load").status == 200
    page.locator("[data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert page.locator(".datatable-count-line").inner_text() == "120 systems"
    assert _first_names(page, 2) == ["Holding 000", "Holding 001"]
    assert page.locator("th[data-col='distance']").get_attribute("aria-sort") == "ascending"
    page.locator(".datatable-scroll").evaluate("el => { el.scrollTop = el.scrollHeight; }")
    page.locator("tbody tr[data-index]", has_text="Holding 119").wait_for(state="attached", timeout=10000)

    page.locator("th[data-col='distance'] .datatable-sort").click()
    page.locator("th[data-col='distance'][aria-sort='descending']").wait_for(state="attached")
    page.locator("tbody tr[data-index='0']", has_text="Holding 119").wait_for(state="attached")
    assert parse_qs(urlparse(page.url).query) == {"systems_sort": ["distance"], "systems_order": ["desc"]}


def test_the_sector_page_runs_its_map_and_its_table_together(page, table_site):
    """The map script finds the Contents buttons through the table's "rows drawn" event; the
    page fixture fails this on any script error from either."""
    assert page.goto(f"{table_site}/sector/5", wait_until="load").status == 200
    page.locator("#sector-contents [data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
    page.locator("#galaxymap3d-canvas, canvas").first.wait_for(state="attached", timeout=15000)


def test_a_sectors_contents_scroll_filter_and_keep_their_map_buttons(page, table_site):
    # The map is not under test here, and its software-rendered frames starve the page of time.
    page.route("**/galaxymap3d.js*", lambda route: route.abort())
    assert page.goto(f"{table_site}/sector/5", wait_until="load").status == 200
    page.locator("#sector-contents [data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
    count = page.locator("#sector-contents .datatable-count-line")
    assert count.inner_text() == "132 items"  # 130 systems, the nebula and the folded rogue planets
    assert _first_names(page, 2) == ["Star 000", "Star 001"]
    # A Location cell keeps its links.
    assert page.locator("#sector-contents tbody tr[data-index='0'] td").nth(4).locator("a").count() == 1
    page.locator("#sector-contents .datatable-scroll").evaluate("el => { el.scrollTop = el.scrollHeight; }")
    page.locator("#sector-contents tbody tr[data-index]", has_text="2 rogue planets").wait_for(
        state="attached", timeout=10000)
    # The folded row's link lists the planets one by one, each with its (hidden) map button.
    page.locator("#sector-contents tbody a", has_text="2 rogue planets").click()
    page.locator("#sector-contents [data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert count.inner_text() == "2 items match"
    assert _first_names(page, 2) == ["Drifter", "Wanderer"]
    buttons = page.locator("#sector-contents tbody [data-map-target]")
    assert buttons.count() == 2 and buttons.first.get_attribute("data-map-target") == "rogue_planet:21"
    assert parse_qs(urlparse(page.url).query) == {"contents_type": ["Rogue Planet"]}



def _sign_in(page, base):
    page.context.add_cookies([{"name": "pg_admin_session", "value": "signed-in", "url": base}])


def test_the_api_keys_scroll_filter_and_keep_their_revoke_buttons(page, admin_site):
    _sign_in(page, admin_site)
    assert page.goto(f"{admin_site}/admin?keys_sort=label", wait_until="load").status == 200
    page.locator("#api-keys [data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert page.locator("#api-keys .datatable-count-line").inner_text() == "120 keys"
    page.locator("#api-keys .datatable-scroll").evaluate("el => { el.scrollTop = el.scrollHeight; }")
    page.locator("#api-keys tbody tr[data-index]", has_text="key 119").wait_for(state="attached", timeout=10000)
    form = page.locator("#api-keys tbody tr[data-index]", has_text="key 119").locator("form")
    assert form.get_attribute("method") == "post" and form.get_attribute("action") == "/admin"
    names = form.locator("input[type=hidden]").evaluate_all("els => els.map(e => e.name)")
    assert names == ["action", "key_id", "csrf_token"]
    page.locator("#api-keys summary", has_text="Status").click()
    page.locator("#api-keys input[value='revoked']").check()
    page.locator("#api-keys .datatable-count-line", has_text="24 keys match").wait_for(state="attached", timeout=10000)


def test_a_table_with_no_sortable_columns_still_scrolls(page, admin_site):
    _sign_in(page, admin_site)
    assert page.goto(f"{admin_site}/admin/queue", wait_until="load").status == 200
    page.locator("#jobs [data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert page.locator("#jobs .datatable-count-line").inner_text() == "120 jobs"
    page.locator("#jobs .datatable-scroll").evaluate("el => { el.scrollTop = el.scrollHeight; }")
    page.locator("#jobs tbody tr[data-index]", has_text="Job 119").wait_for(state="attached", timeout=10000)
    assert page.locator("#jobs th a.datatable-sort").count() == 0


def test_search_results_scroll_through_every_match(page, table_site):
    assert page.goto(f"{table_site}/search?sector_q=Found", wait_until="load").status == 200
    page.locator("#search-sectors [data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert page.locator("#search-sectors .datatable-count-line").inner_text() == "120 sectors"
    page.locator("#search-sectors .datatable-scroll").evaluate("el => { el.scrollTop = el.scrollHeight; }")
    page.locator("#search-sectors tbody tr[data-index]", has_text="Found 119").wait_for(state="attached", timeout=10000)
    assert page.url.endswith("sector_q=Found")


# --- UX.117: the Apply filters button never goes away ------------------------------------------------

def _wait_held(page, held, count=1):
    deadline = time.monotonic() + 10
    while len(held) < count and time.monotonic() < deadline:
        page.wait_for_timeout(50)
    assert len(held) >= count, "the table never asked for its rows"


def test_apply_filters_stays_while_the_first_rows_load_and_after(page, table_site):
    """The first rows are held back: the button is there while they load and still there once the table is live."""
    held = []
    page.route("**/table/**", lambda route: held.append(route))
    assert page.goto(f"{table_site}/phenomena", wait_until="load").status == 200
    apply = page.locator(".datatable-apply")
    _wait_held(page, held)
    assert apply.is_visible()
    for route in held:
        route.continue_()
    page.locator("[data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
    assert apply.is_visible()


def test_apply_filters_works_when_the_first_load_fails(page, table_site):
    """A failed first fetch leaves the server's table; the button then submits the form to the server."""
    failing = [True]
    page.route("**/table/**", lambda route: route.abort() if failing[0] else route.continue_())
    assert page.goto(f"{table_site}/phenomena", wait_until="load").status == 200
    apply = page.locator(".datatable-apply")
    assert apply.is_visible()
    kind = page.locator(".datatable-facet[data-param='type']")
    kind.locator("summary").click()
    kind.get_by_label("Rogue Planet").check()
    failing[0] = False                      # the server answers from here on
    with page.expect_navigation(timeout=10000):
        apply.click()
    # The server was asked for the filtered page: its address carries the filter (the test site renders every
    # row whatever the address says, so the count line is not what this checks).
    assert parse_qs(urlparse(page.url).query) == {"type": ["rogue_planet"]}
    assert page.locator(".datatable-facet[data-param='type'] input:checked").count() == 1


def test_apply_filters_retries_a_reload_that_failed(page, table_site):
    """The table is live, a later fetch fails, and the same button asks again."""
    _open(page, table_site)
    page.route("**/table/**", lambda route: route.abort())
    kind = page.locator(".datatable-facet[data-param='type']")
    kind.locator("summary").click()
    kind.get_by_label("Rogue Planet").check()
    apply = page.locator(".datatable-apply")
    assert apply.is_visible()
    page.unroute("**/table/**")
    apply.click()
    rogue = sum(1 for row in ROWS if row["type"] == "rogue_planet")
    _wait_count(page, f"{rogue} phenomena match")
    assert parse_qs(urlparse(page.url).query) == {"type": ["rogue_planet"]}


def test_every_filtered_table_keeps_its_apply_button(page, table_site, admin_site):
    for base, path, root in ((table_site, "/sector/5", "#sector-contents"), (admin_site, "/admin?keys_sort=label", "#api-keys")):
        if base is admin_site:
            _sign_in(page, admin_site)
        page.route("**/galaxymap3d.js*", lambda route: route.abort())
        assert page.goto(f"{base}{path}", wait_until="load").status == 200
        page.locator(f"{root} [data-datatable][data-enhanced='true']").wait_for(state="attached", timeout=15000)
        assert page.locator(f"{root} .datatable-apply").is_visible(), root
