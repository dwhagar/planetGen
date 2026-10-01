# tests/test_galaxy_urls_combined.py

"""
TEST.50: Galaxy Map URLs combined and abused. `/galaxy` with the
drill-down's `?at=`, `?p=` and `?sector=` (read by
`static/galaxystages.js` in the browser, not by the server) together with
the server's own `?quadrant=`, `?page=`, `?pick=` and `?course=`;
`?course=` naming objects that were deleted; and `/galaxy/locate` (and
the API's `/api/galaxy/locate` behind it) fed unicode, very long input,
NaN/inf-looking text, FULLTEXT operators and ambiguous names, plus
out-of-range coordinates for `/api/galaxy/cell` and `/galaxy/stage`.

Each answer must be a page, a redirect or a clean 4xx -- never a 500 or
a traceback. The apps here run with `TESTING` set, so an unhandled
exception fails the test outright rather than hiding behind a 500 page.

The first group fakes the data layer with `test_web_galaxy.py`'s own
fixtures; the rest run in-process against a real throwaway database and
are skipped without a MySQL server.
"""

import re
from urllib.parse import quote

import pytest

from stellarObjects import _db

import web  # noqa: F401 -- puts src/html/lib on sys.path
import apiclient  # noqa: E402

# Imported fixtures (see test_bughunt_api_gaps.py for why this works).
from tests.test_web_galaxy import _place_sector, _scene, app, client, db_client, fake  # noqa: F401,E402


def _assert_clean_page(response, statuses=(200,)):
    """A rendered page (or redirect) with an allowed status and no
    traceback or raw exception text in it."""
    assert response.status_code in statuses, response.status_code
    text = response.get_data(as_text=True)
    assert "Traceback" not in text and "Internal Server Error" not in text
    return text


def _scene_without_volatile(html):
    scene = _scene(html)
    scene.pop("csrf", None)
    scene.pop("initial", None)  # tile-cache bookkeeping differs between a first and a later request
    return scene


# --- ?at= / ?p= / ?sector= with the server's own parameters (faked data) ---------------

DRILL_COMBOS = [
    "at=243.7.14.0&p=1,2,3&sector=R5.L1.S20",
    "at=243.7.14.0&sector=R5.L1.S20",
    "p=4,5&sector=R0.L0.S0",
    "at=1.0.0.0&p=&sector=",
    "at=nan.inf.0.0&p=nan,inf,-inf,1e999&sector=R-1.L99999999999999999999.S0",
    "at=243.7.14.0&at=27.1.2.3&p=1&p=2&sector=a&sector=b",
    "at=%00&p=%00&sector=%00",
    "at=" + "9" * 5000 + "&p=" + ",".join(["1"] * 2000) + "&sector=" + "R" * 5000,
    "at=%E2%98%83&p=%F0%9F%9A%80&sector=%E5%90%8D",
]


@pytest.mark.parametrize("query", DRILL_COMBOS)
def test_drill_parameters_are_left_to_the_browser(client, fake, query):
    """The server ignores the drill-down's own parameters: any mix of
    them -- even nonsense -- renders the same map the bare URL does."""
    plain = _scene_without_volatile(_assert_clean_page(client.get("/galaxy")))
    html = _assert_clean_page(client.get(f"/galaxy?{query}"))
    assert _scene_without_volatile(html) == plain


@pytest.mark.parametrize("extra", ["quadrant=II", "quadrant=II&page=2", "quadrant=nonsense&page=-4",
                                   "pick=from&to=system:1", "course=system:1"])
def test_drill_parameters_with_the_servers_own(client, fake, extra):
    html = _assert_clean_page(client.get(f"/galaxy?at=243.7.14.0&p=1,2&sector=R5.L1.S20&{extra}"))
    assert "galaxymap3d-data" in html


def test_drill_parameters_are_never_echoed_raw(client, fake):
    payload = '"><script>alert(1)</script>'
    html = _assert_clean_page(client.get("/galaxy", query_string={"at": payload, "p": payload, "sector": payload}))
    assert "<script>alert(1)</script>" not in html


def test_drill_parameters_with_a_course(client, fake, monkeypatch):
    """`?course=` is the one the server reads; the drill-down parameters
    beside it change nothing about the course it draws."""
    from web import nav_page

    course = {"scope": "sector", "navUrl": "/nav?from=system:1&to=system:2", "points": [],
              "sector": {"id": 9, "name": "Nine", "ring": 5, "layer": 1, "slot": 20}}
    asked = []
    monkeypatch.setattr(nav_page, "galaxy_course", lambda f, t: asked.append((f, t)) or course)
    html = _assert_clean_page(client.get("/galaxy?course=system:1,system:2&at=243.7.14.0&p=1&sector=R1.L1.S1"))
    assert asked == [("system:1", "system:2")]
    assert _scene(html)["course"]["sector"]["id"] == 9


# --- /galaxy/locate input (faked API) -----------------------------------------------------

@pytest.mark.parametrize("q", [
    "Ærøskøbing", "名前", "🚀 Rocket", "‮evil", "éé", "nan", "inf", "-inf", "1e999",
    "NaN,NaN,NaN", "1e308,-1e308,0", "-99999999999999999999999", "x" * 10000, "%", "_", "\\", "'", '"',
])
def test_locate_endpoint_passes_odd_input_through_trimmed(client, fake, monkeypatch, q):
    """Whatever is typed reaches the API as text, cut to 200 characters,
    and the answer is JSON."""
    calls = []
    monkeypatch.setattr(apiclient, "get_galaxy_locate", lambda db, term: calls.append(term) or [])
    response = client.get("/galaxy/locate", query_string={"q": q})
    assert response.status_code == 200
    assert response.get_json() == {"matches": []}
    assert calls == [q.strip()[:200]]


def test_locate_endpoint_not_found_from_the_api_is_json(client, fake, monkeypatch):
    def missing(db, q):
        raise apiclient.NotFoundError("no such database")

    monkeypatch.setattr(apiclient, "get_galaxy_locate", missing)
    response = client.get("/galaxy/locate?q=Vega")
    # JSON like its sibling endpoints (/galaxy/tiles, /galaxy/stage), not
    # the HTML 404 page.
    assert response.status_code == 404
    assert response.get_json() == {"error": "no such database"}


# --- Real database ------------------------------------------------------------------------

def _system_in(mysql_config, sector_id, name):
    """Adds one star system named `name` to `sector_id` (raw rows are
    enough for the locate and course lookups)."""
    from stellarObjects.config import SystemConfig
    from stellarObjects.spaceSector import SpaceSector
    from stellarObjects.systemData import StarSystem

    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    holder = SpaceSector("Holder " + name, edge_ly=10.0)
    holder.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    holder_id = _db.save_sector(holder, config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (holder_id,)).fetchone()["id"]
            conn.execute("UPDATE star_systems SET name = ?, sector_id = ? WHERE id = ?", (name, sector_id, system_id))
            conn.execute("UPDATE stars SET name = ? WHERE star_system_id = ?", (name, system_id))
            conn.execute("DELETE FROM sectors WHERE id = ?", (holder_id,))
    finally:
        conn.close()
    return system_id


def _delete(mysql_config, table, row_id):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.execute(f"DELETE FROM {table} WHERE id = ?", (row_id,))
    finally:
        conn.close()


def test_course_to_deleted_objects_is_a_clean_404(db_client, mysql_config):
    sector_id = _place_sector(mysql_config, "Course Sector", address=(5, 1, 20))
    gone = _system_in(mysql_config, sector_id, "Gone Away")
    kept = _system_in(mysql_config, sector_id, "Still Here")
    _delete(mysql_config, "star_systems", gone)

    for course in (f"system:{gone},system:{kept}", f"system:{kept},system:{gone}", f"{gone},{kept}",
                   f"nebula:999999999,system:{kept}", f"system:{kept},black_hole:999999999",
                   f"system:{10 ** 25},system:{kept}", f"system:{kept},bogus:1", f"system:x,system:{kept}"):
        _assert_clean_page(db_client.get(f"/galaxy?course={quote(course)}"), statuses=(404,))
        # The drill-down parameters beside it don't change that.
        _assert_clean_page(db_client.get(f"/galaxy?course={quote(course)}&at=243.7.14.0&p=1&sector=R5.L1.S20"),
                           statuses=(404,))


def test_course_whose_sector_was_deleted_still_renders(db_client, mysql_config):
    """Deleting the sector detaches its systems (`sector_id` NULL): the
    course's endpoints still exist, so the map renders -- with or
    without a course drawn -- rather than failing."""
    sector_id = _place_sector(mysql_config, "Vanishing Sector", address=(6, 0, 7))
    first = _system_in(mysql_config, sector_id, "Left Behind")
    second = _system_in(mysql_config, sector_id, "Also Left")
    _delete(mysql_config, "sectors", sector_id)

    html = _assert_clean_page(db_client.get(f"/galaxy?course=system:{first},system:{second}&at=243.7.14.0"))
    assert "galaxymap3d-data" in html
    # The same system at both ends.
    _assert_clean_page(db_client.get(f"/galaxy?course=system:{first},system:{first}"))


def test_real_page_with_every_parameter_at_once(db_client, mysql_config):
    sector_id = _place_sector(mysql_config, "Busy Sector", address=(5, 1, 20))
    system_id = _system_in(mysql_config, sector_id, "Busy Star")
    query = (f"/galaxy?at=243.7.14.0&p=1,2&sector=R5.L1.S20&quadrant=III&page=999"
             f"&pick=from&to=system:{system_id}&course=system:{system_id},system:{system_id}")
    _assert_clean_page(db_client.get(query))


@pytest.mark.parametrize("q", [
    "Ærøskøbing", "名前", "🚀", "‮", "é", "nan", "inf", "-inf", "1e999", "NaN,NaN,NaN",
    "1e308,-1e308,0", "-99999999999999999999999", "%", "_", "\\", "'", '"', "*", "+Vega -Vega",
    "\"Vega", "Vega*", "(Vega)", "~Vega", "<Vega>", "@3", "a b c d e f g h i j k l m n o p",
    "x" * 300, " ".join(["Vega"] * 200), "\x00", "Vega\x00Prime",
])
def test_locate_odd_input_against_the_real_database(db_client, mysql_config, q):
    """Unicode, FULLTEXT boolean operators, LIKE wildcards, numbers that
    look like coordinates and very long input: matches or none, never an
    error. Through the page (cut to 200 characters) and straight to the
    API (no cut)."""
    _place_sector(mysql_config, "Vega Prime", address=(5, 1, 20))
    response = db_client.get("/galaxy/locate", query_string={"q": q})
    assert response.status_code == 200, response.get_data(as_text=True)[:300]
    assert isinstance(response.get_json()["matches"], list)
    response = db_client.get("/api/galaxy/locate", query_string={"q": q})
    assert response.status_code == 200, response.get_data(as_text=True)[:300]
    assert isinstance(response.get_json()["matches"], list)


def test_locate_very_long_input_straight_to_the_api(db_client, mysql_config):
    _place_sector(mysql_config, "Vega Prime", address=(5, 1, 20))
    for q in ("Vega " * 2000, "x" * 20000):
        response = db_client.get("/api/galaxy/locate", query_string={"q": q})
        assert response.status_code == 200
        assert response.get_json() == {"matches": []} or all(
            "Vega" in m["name"] for m in response.get_json()["matches"])


def test_locate_ambiguous_names(db_client, mysql_config):
    """Several sectors and systems share a name's words: exact names
    first, then names starting with it, then the rest, sectors and
    systems mixed; a sector without an address and a system outside one
    aren't offered."""
    vega = _place_sector(mysql_config, "Vega", address=(5, 1, 20))
    _place_sector(mysql_config, "Vega Prime", address=(5, 1, 21))
    _place_sector(mysql_config, "Old Vega", address=(5, 1, 22))
    _system_in(mysql_config, vega, "Vega Minor")
    unplaced = _place_sector(mysql_config, "Vega Unplaced", address=(5, 1, 23))
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE sectors SET ring_index = NULL, layer_index = NULL, ring_slot_index = NULL"
                         " WHERE id = ?", (unplaced,))
    finally:
        conn.close()

    matches = db_client.get("/galaxy/locate?q=vega").get_json()["matches"]
    names = [m["name"] for m in matches]
    assert names[0] == "Vega"
    assert set(names) == {"Vega", "Vega Prime", "Old Vega", "Vega Minor"}
    assert names.index("Old Vega") > max(names.index("Vega Prime"), names.index("Vega Minor"))
    minor = next(m for m in matches if m["name"] == "Vega Minor")
    assert (minor["kind"], minor["sector_id"], minor["ring"], minor["slot"]) == ("system", vega, 5, 20)

    # Still ambiguous after the exact match is gone: no 500, the rest
    # remain. Asked of the API itself: the page's answer comes through
    # `lib/pagecache.py`, which a raw SQL delete doesn't invalidate.
    _delete(mysql_config, "sectors", vega)
    names = [m["name"] for m in db_client.get("/api/galaxy/locate?q=vega").get_json()["matches"]]
    assert set(names) == {"Vega Prime", "Old Vega"}  # Vega Minor lost its sector


# --- Out-of-range coordinates and stage keys ---------------------------------------------------

@pytest.mark.parametrize("query", [
    "x=inf&y=0&z=0", "x=0&y=-inf&z=0", "x=0&y=0&z=nan", "x=1e309&y=0&z=0", "x=1e308&y=1e308&z=1e308",
    "x=-1e300&y=1e300&z=0", "x=1,2&y=0&z=0", "x=%E2%91%A0&y=0&z=0",
    f"ring=0&layer={10 ** 30}&slot=0", f"ring=5&layer=0&slot={10 ** 30}", "ring=-5&layer=0&slot=0",
])
def test_cell_out_of_range_coordinates_are_400(db_client, mysql_config, query):
    _db.get_connection(mysql_config).close()
    response = db_client.get(f"/api/galaxy/cell?{query}")
    assert response.status_code == 400, response.get_data(as_text=True)[:300]
    assert "error" in response.get_json()


def test_cell_far_past_the_galaxy_by_ring_is_described_not_refused(db_client, mysql_config):
    """Pinned as found: a ring far past any galaxy (10**30) is still an
    integer address, so the cell is described (about 4e30 pc out, with
    `in_galaxy` null before a plan) rather than refused -- a clean answer
    either way, never an overflow."""
    _db.get_connection(mysql_config).close()
    response = db_client.get(f"/api/galaxy/cell?ring={10 ** 30}&layer=0&slot=0")
    assert response.status_code in (200, 400)
    if response.status_code == 200:
        assert response.get_json()["in_galaxy"] is None


@pytest.mark.parametrize("at", [
    "nan.0.0.0", "243.inf.0.0", "243.-1.0.0", "243.7.14.0.0", "243.7.14",
    "243.7.14.0;drop", "%E2%91%A1.0.0.0", "",
])
def test_stage_out_of_range_keys_against_the_real_database(db_client, mysql_config, at):
    _db.get_connection(mysql_config).close()
    response = db_client.get(f"/galaxy/stage?at={at}")
    assert response.status_code in ((200,) if at == "" else (400,)), response.get_data(as_text=True)[:300]
    assert isinstance(response.get_json(), dict)
    if response.status_code == 400:
        assert re.search(r"\S", response.get_json()["error"])


def test_stage_far_past_the_galaxy_is_an_empty_block(db_client, mysql_config):
    """Pinned as found: a block ring far past any galaxy is a well-formed
    key, so it answers an empty stage rather than a 400."""
    _db.get_connection(mysql_config).close()
    response = db_client.get(f"/galaxy/stage?at=243.{10 ** 30}.0.0")
    assert response.status_code in (200, 400)
    if response.status_code == 200:
        assert response.get_json()["children"] == []
