"""
UX.231: the line "in Vang Zarh Region" in the At a glance block of a sector, a system and an object page. The API answers carry
`region` (`{name, parent}`) for a position inside a charted region and `None` for one outside every region, in a cell with no generated
sector, or while no naming key is in force; a page draws the line from it and draws nothing where it is `None`.
"""

import pytest

from planetgen.api import routes
from planetgen.db import region_stats, store
from planetgen.galaxy import named_regions, spherical
from planetgen.names import naming_key
from planetgen.web import glance, item_pages
from tests.test_region_stats import _addresses
from tests.test_region_wire_names import EDGE, KEY, galaxy, keyed  # noqa: F401 -- fixtures

NAMED = {"name": "Vang Zarh Region", "parent": "Vang Reach"}


def test_a_charted_region_gives_its_line_and_no_region_gives_none():
    assert item_pages.region_row({"region": NAMED}) == {"label": "Region", "html": "in Vang Zarh Region"}
    assert item_pages.region_row({"region": None}) is None
    assert item_pages.region_row({}) is None
    assert item_pages.region_row({"region": {"name": None, "parent": None}}) is None


def test_the_name_is_escaped():
    assert item_pages.region_row({"region": {"name": "<b>X</b> Region"}})["html"] == "in &lt;b&gt;X&lt;/b&gt; Region"


def test_a_sector_and_a_system_page_list_the_line_last_and_only_when_there_is_a_region():
    sector = {"systems": [], "phenomena": [], "region": NAMED}
    assert glance.sector_rows(sector, [], {})[-1]["label"] == "Region"
    assert str(glance.sector_rows(sector, [], {})[-1]["html"]) == "in Vang Zarh Region"
    assert all(row["label"] != "Region" for row in glance.sector_rows({**sector, "region": None}, [], {}))
    system = {"stars": [], "planets": [], "region": NAMED}
    assert glance.system_rows(system, [], lambda *a, **k: "")[-1]["label"] == "Region"
    assert all(row["label"] != "Region" for row in glance.system_rows({**system, "region": None}, [], lambda *a, **k: ""))


def test_an_object_context_gets_the_line_after_its_other_rows():
    context = {"glance": [{"label": "Mass", "html": "1 kg"}]}
    assert [row["label"] for row in item_pages.with_region(context, {"region": NAMED})["glance"]] == ["Mass", "Region"]
    assert item_pages.with_region(context, {"region": None}) is context



@pytest.mark.db
def test_a_position_inside_a_charted_cell_is_named_and_one_outside_is_not(galaxy, keyed):
    conn = store.get_connection(galaxy)
    try:
        charted = _addresses()[0]
        inside = spherical.sector_position_pc(*charted, EDGE)
        found = routes.region_of(conn, inside)
        cell = region_stats.region_cell(charted)
        assert found == {"name": named_regions.region_name(cell, KEY, EDGE),
                         "parent": named_regions.family_name(named_regions.family_of(cell), KEY, EDGE)}
        assert found["name"].endswith(" Region")
        empty = next(address for address in ((200, 150, 0), (210, 150, 0), (220, 150, 0))
                     if region_stats.region_cell(address) not in {region_stats.region_cell(a) for a in _addresses()})
        assert routes.region_of(conn, spherical.sector_position_pc(*empty, EDGE)) is None
        assert routes.region_of(conn, None) is None
    finally:
        conn.close()


@pytest.mark.db
def test_without_a_naming_key_no_region_is_named(galaxy):
    conn = store.get_connection(galaxy)
    try:
        assert routes.region_of(conn, spherical.sector_position_pc(*_addresses()[0], EDGE)) is None
    finally:
        conn.close()


# --- The pages, on the small generated galaxy of test_web_a11y -------------------------------------

from tests.publicids import pid  # noqa: E402
from tests.test_web_a11y import site_app, site_db  # noqa: E402,F401 -- fixtures


@pytest.fixture(scope="module")
def charted(site_db):
    """The sample galaxy with its regions stored, and the printed IDs of a sector, a system and a planet in it."""
    conn = store.get_connection(site_db)
    try:
        with conn:
            region_stats.rebuild(conn)
        sector = conn.execute("SELECT s.id FROM sectors s JOIN star_systems ss ON ss.sector_id = s.id WHERE s.center_x_pc IS NOT NULL"
                              " ORDER BY s.id LIMIT 1").fetchone()
        system = conn.execute("SELECT ss.id FROM star_systems ss JOIN planets p ON p.star_system_id = ss.id WHERE ss.sector_id = ?"
                              " ORDER BY ss.id LIMIT 1", (sector["id"],)).fetchone()
        planet = conn.execute("SELECT id FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 1", (system["id"],)).fetchone()
        return {"sector": pid("sector", sector["id"], conn), "system": pid("system", system["id"], conn),
                "planet": pid("planet", planet["id"], conn)}
    finally:
        conn.close()


def _pages(site_app, charted):
    client = site_app.test_client()
    return {"sector": client.get(f"/sector/{charted['sector']}").get_data(as_text=True),
            "system": client.get(f"/system/{charted['system']}").get_data(as_text=True),
            "planet": client.get(f"/object/planet:{charted['planet']}").get_data(as_text=True)}


@pytest.mark.db
def test_sector_system_and_object_pages_say_which_region_they_are_in(site_app, charted, keyed):
    pages = _pages(site_app, charted)
    for kind, page in pages.items():
        assert '<th scope="row">Region</th><td>in ' in page and " Region</td>" in page, kind
    names = {page.split('<th scope="row">Region</th><td>in ', 1)[1].split("</td>", 1)[0] for page in pages.values()}
    assert len(names) == 1, "a sector, its system and the system's planet are in one region"
