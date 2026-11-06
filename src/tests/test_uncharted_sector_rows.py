# tests/test_uncharted_sector_rows.py

"""
MAP.399: the Sectors list tells generated sectors from uncharted placeholders (the rows `store.uncharted_sector_id` makes where a drifting rogue planet lands in a
cell nothing was generated in). A placeholder is named from its cell under the naming key (`naming_key.uncharted_sector_name`, never the raw hex designation), says it is
uncharted, carries its region, and the header counts the generated sectors with the uncharted number beside it.
"""

import pytest

from planetgen.web.app import create_app
from planetgen.api.config import Config
from planetgen.db import chart_ledger, query, store
from planetgen.galaxy.spherical import provisional_sector_designation
from planetgen.names import naming_key

pytestmark = pytest.mark.db

KEY = "0123abcd"
EDGE_PC = 4.0
CHARTED = [(9, 2, 5), (9, 2, 6)]
UNCHARTED = [(9, 2, 7), (10, 3, 4), (11, 1, 2)]


@pytest.fixture
def keyed():
    naming_key.set_resolver(lambda: KEY)
    yield KEY
    naming_key.set_resolver(None)


@pytest.fixture
def galaxy(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            for address in CHARTED + UNCHARTED:
                sector_id = store.uncharted_sector_id(conn, address, EDGE_PC)
                if address in CHARTED:
                    chart_ledger.record(conn, address, store.sector_run_id(conn), 1)
                    conn.execute("UPDATE sectors SET name = ? WHERE id = ?", (f"Charted {address[2]}", sector_id))
    finally:
        conn.close()
    return mysql_config


@pytest.fixture
def client(galaxy):
    class TestConfig(Config):
        MYSQL_CONFIG = galaxy
        WRITE_MYSQL_CONFIG = galaxy
        CONTROL_MYSQL_CONFIG = galaxy

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def test_the_list_marks_the_placeholders_and_counts_them_beside_the_generated_sectors(galaxy, client, keyed):
    body = client.get("/api/sectors?limit=50").get_json()
    by_name = {item["name"]: item for item in body["items"]}
    assert body["total"] == 5 and body["uncharted_total"] == 3
    assert sorted(item["uncharted"] for item in body["items"]) == [False, False, True, True, True]
    assert {name for name, item in by_name.items() if not item["uncharted"]} == {"Charted 5", "Charted 6"}


def test_a_placeholder_is_named_from_its_cell_not_from_its_hex_designation_and_carries_its_region(galaxy, client, keyed):
    items = [item for item in client.get("/api/sectors?limit=50").get_json()["items"] if item["uncharted"]]
    assert len(items) == 3
    for item in items:
        address = (item["shell_index"], item["band_index"], item["wedge_index"])
        assert item["designation"] == provisional_sector_designation(*address)
        assert item["name"] == naming_key.uncharted_sector_name(address, KEY) and item["name"] != item["designation"]
        assert " " in item["name"] and item["region"]["name"].endswith(" Region") and item["region"]["parent"].endswith(" Reach")
    assert len({item["name"] for item in items}) == 3


def test_without_a_naming_key_a_placeholder_says_uncharted_with_its_designation(galaxy, client):
    items = [item for item in client.get("/api/sectors?limit=50").get_json()["items"] if item["uncharted"]]
    assert all(item["name"] == f"Uncharted sector {item['designation']}" and item["region"] is None for item in items)


def test_the_sector_detail_of_a_placeholder_is_uncharted_and_a_generated_one_is_not(galaxy, client, keyed):
    items = client.get("/api/sectors?limit=50").get_json()["items"]
    placeholder = next(item for item in items if item["uncharted"])
    generated = next(item for item in items if not item["uncharted"])
    detail = client.get(f"/api/sectors/{placeholder['id']}").get_json()
    assert detail["uncharted"] and detail["name"] == placeholder["name"] and detail["region"] == placeholder["region"]
    assert not client.get(f"/api/sectors/{generated['id']}").get_json()["uncharted"]


def test_the_uncharted_cell_answer_names_the_cell_and_its_region(galaxy, client, keyed):
    shell, band, wedge = UNCHARTED[0]
    cell = client.get(f"/api/galaxy/uncharted?shell={shell}&band={band}&wedge={wedge}").get_json()
    assert cell["name"] == naming_key.uncharted_sector_name(UNCHARTED[0], KEY) and cell["region"]["name"].endswith(" Region")
    assert cell["sector_id"] is None, "a placeholder row is not a generated sector"


def test_neighbouring_cells_have_different_first_words_and_a_key_changes_every_name():
    first = [naming_key.uncharted_sector_name((9, 2, wedge), KEY) for wedge in range(1, 9)]
    assert len(set(first)) == 8 and len({name.split()[0] for name in first}) >= 6
    assert naming_key.uncharted_sector_name((9, 2, 5), KEY) != naming_key.uncharted_sector_name((9, 2, 5), "0123abce")


def test_the_uncharted_count_is_the_placed_rows_without_a_ledger_cell(galaxy):
    conn = store.get_connection(galaxy)
    try:
        assert query.count_uncharted_sectors(conn) == 3
    finally:
        conn.close()
