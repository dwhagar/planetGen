# tests/test_db_collation.py

"""
TEST.13: collation collisions. Names that `utf8mb4_unicode_ci` holds
equal but that differ byte for byte -- case ("Vega"/"vega"/"VEGA"),
accents ("Vega"/"Véga"), an ignorable character -- go through the name
registries (`reserve_system_names`, `confirm_system_names`,
`reserve_sector_name`) as one name: each later holder gets a Greek
letter or a diminutive, the first keeps its own spelling, and no two
rows end up with names the database (and so the web) can't tell apart.
"Æron" is not "Aeron" to the collation, so both stay bare.
`name_in_use` and the renames see the same equality.

Then a database created with a non-default collation (`utf8mb4_general_ci`,
`latin1`, `utf8mb4_bin`) still saves, loads and searches a sector -- the
tables declare their own collation, so nothing mixes collations (error
1267).

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import uuid

import pymysql
import pytest

from planetgen.db import query
from planetgen.db import store
from planetgen.db.store import MySQLConfig
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem

from tests.conftest import _test_server_kwargs


def _system(name):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    system.name = name
    return system, cfg


def _sector(name, system_names, rogue_names=()):
    sector = SpaceSector(name, edge_ly=11.5)
    for index, system_name in enumerate(system_names):
        system, cfg = _system(system_name)
        sector.add_system(system, position=(float(index), 0.0, 0.0), system_config=cfg)
    for index, rogue_name in enumerate(rogue_names):
        sector.add_phenomenon(RoguePlanet(SystemConfig(), name=rogue_name), "rogue-planet",
                              position=(0.0, float(index) + 0.5, 0.0))
    return sector


def _names(conn, table):
    return [r["name"] for r in conn.execute(f"SELECT name FROM {table} ORDER BY id").fetchall()]


def _no_collation_duplicates(conn):
    for table in store.UNIQUE_NAME_TABLES:
        duplicates = conn.execute(f"SELECT name FROM {table} GROUP BY name HAVING COUNT(*) > 1").fetchall()
        assert duplicates == [], (table, duplicates)


# ---------------------------------------------------------------------------
# The registries
# ---------------------------------------------------------------------------

def test_case_variants_saved_one_by_one_are_one_name(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        for name in ("Vega", "vega", "VEGA"):
            system, cfg = _system(name)
            with conn:
                store.insert_star_system(conn, system, cfg)
        assert _names(conn, "star_systems") == ["Alpha Vega", "Beta vega", "Gamma VEGA"]
        registry = conn.execute("SELECT base_name, occurrence_count FROM system_name_registry").fetchall()
        assert [(r["base_name"], r["occurrence_count"]) for r in registry] == [("Vega", 3)]
        _no_collation_duplicates(conn)
    finally:
        conn.close()


def test_accent_variants_in_one_sector_are_one_name(mysql_config):
    # Before the fix, the batched reservation keyed names by casefold():
    # "Vega" and "Véga" both came out "Beta ...", equal to the database.
    sector_id = store.save_sector(_sector("Lyra Reach", ["Vega", "Véga", "VÉGA"], ["vegà"]), config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        assert _names(conn, "star_systems") == ["Alpha Vega", "Beta Véga", "Gamma VÉGA"]
        assert _names(conn, "rogue_planets") == ["Delta vegà"]
        assert conn.execute("SELECT COUNT(*) AS n FROM system_name_registry").fetchone()["n"] == 1
        _no_collation_duplicates(conn)
        loaded = store.load_sector(conn, sector_id)
        assert sorted(e.star_system.name for e in loaded.entries) == ["Alpha Vega", "Beta Véga", "Gamma VÉGA"]
    finally:
        conn.close()


def test_accent_variant_of_a_stored_name_renames_the_holder_in_its_own_spelling(mysql_config):
    store.save_sector(_sector("First Reach", ["Vega"]), config=mysql_config)
    store.save_sector(_sector("Second Reach", ["Véga"]), config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        assert _names(conn, "star_systems") == ["Alpha Vega", "Beta Véga"]
        # The single star shares its system's name and follows the rename.
        assert _names(conn, "stars") == ["Alpha Vega", "Beta Véga"]
        _no_collation_duplicates(conn)
    finally:
        conn.close()


def test_a_spelling_python_folds_apart_still_shares_its_registry_row(mysql_config):
    # A zero-width space is ignored by the collation but not by
    # `_name_key`: the two names still land on one registry row.
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            reserved = store.reserve_system_names(conn, ["Vega", "Vega​"])
        assert [r[0] for r in reserved] == ["Alpha Vega", "Beta Vega​"]
        assert conn.execute("SELECT occurrence_count FROM system_name_registry").fetchone()["occurrence_count"] == 2
    finally:
        conn.close()


@pytest.mark.parametrize("sector_name", ["VEGA", "Véga"])
def test_system_after_a_differently_spelled_sector_gets_a_diminutive(mysql_config, sector_name):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.insert_sector(conn, SpaceSector(name=sector_name))
        system, cfg = _system("Vega")
        with conn:
            store.insert_star_system(conn, system, cfg)
        assert _names(conn, "sectors") == [sector_name]
        assert _names(conn, "star_systems") == ["Little Vega"]
    finally:
        conn.close()


def test_sector_after_a_differently_spelled_system_renames_the_system(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        system, cfg = _system("Vega")
        with conn:
            store.insert_star_system(conn, system, cfg)
        with conn:
            store.insert_sector(conn, SpaceSector(name="Véga"))
        assert _names(conn, "sectors") == ["Véga"]
        assert _names(conn, "star_systems") == ["Little Vega"]
    finally:
        conn.close()


def test_sector_variants_keep_the_first_holders_spelling(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        for name in ("Vega", "Véga", "vega"):
            with conn:
                store.insert_sector(conn, SpaceSector(name=name))
        assert _names(conn, "sectors") == ["Alpha Vega", "Beta Véga", "Gamma vega"]
        _no_collation_duplicates(conn)
    finally:
        conn.close()


def test_ae_ligature_is_a_different_name(mysql_config):
    store.save_sector(_sector("Ligature Reach", ["Aeron", "Æron"]), config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        assert _names(conn, "star_systems") == ["Aeron", "Æron"]
        assert conn.execute("SELECT COUNT(*) AS n FROM system_name_registry").fetchone()["n"] == 2
    finally:
        conn.close()


def test_name_in_use_and_renames_see_collation_equal_names(mysql_config):
    store.save_sector(_sector("Rename Reach", ["Vega", "Altair"], ["Drifter"]), config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        vega = conn.execute("SELECT id FROM star_systems WHERE name = 'Vega'").fetchone()["id"]
        altair = conn.execute("SELECT id FROM star_systems WHERE name = 'Altair'").fetchone()["id"]
        rogue = conn.execute("SELECT id FROM rogue_planets").fetchone()["id"]
        # The API's 409 check: another row's name in any spelling is taken.
        for spelling in ("VEGA", "véga", "Vega "):
            assert store.name_in_use(conn, spelling, exclude=[("star_systems", altair)]) == "star_systems"
        assert store.name_in_use(conn, "DRIFTER") == "rogue_planets"
        assert store.name_in_use(conn, "rename reach") == "sectors"
        assert store.name_in_use(conn, "Æltair") is None
        star = conn.execute("SELECT id FROM stars WHERE star_system_id = ?", (vega,)).fetchone()["id"]
        own = [("star_systems", vega), ("stars", star)]
        assert store.name_in_use(conn, "VÉGA", exclude=own) is None

        # A respelling of its own name is a real rename, star included.
        with conn:
            assert store.rename_star_system(conn, vega, "VÉGA")
            assert store.rename_phenomenon(conn, "rogue_planets", rogue, "drifter")
        assert conn.execute("SELECT name FROM star_systems WHERE id = ?", (vega,)).fetchone()["name"] == "VÉGA"
        assert conn.execute("SELECT name FROM stars WHERE id = ?", (star,)).fetchone()["name"] == "VÉGA"
        assert _names(conn, "rogue_planets") == ["drifter"]
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# A database whose default collation isn't the tables'
# ---------------------------------------------------------------------------

@pytest.fixture(params=[
    "CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci",
    "CHARACTER SET utf8mb4 COLLATE utf8mb4_bin",
    "CHARACTER SET latin1",
])
def odd_collation_config(request, _mysql_server_available):
    kwargs = _test_server_kwargs()
    db_name = f"planetgen_test_{uuid.uuid4().hex[:16]}"
    admin = pymysql.connect(**kwargs)
    try:
        with admin.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{db_name}` {request.param}")
        admin.commit()
    finally:
        admin.close()
    config = MySQLConfig(database=db_name, **kwargs)
    try:
        yield config
    finally:
        store.close_pool(config)
        admin = pymysql.connect(**kwargs)
        try:
            with admin.cursor() as cur:
                cur.execute(f"DROP DATABASE IF EXISTS `{db_name}`")
            admin.commit()
        finally:
            admin.close()


def test_a_database_with_another_default_collation_reads_and_writes(odd_collation_config, monkeypatch, tmp_path):
    from tests.test_fuzz_web_routes import make_app

    config = odd_collation_config
    sector_id = store.save_sector(_sector("Véga Reach", ["Vega", "Véga", "Æron"], ["Drifter"]), config=config)
    conn = store.get_connection(config)
    try:
        collations = {r["c"] for r in conn.execute(
            "SELECT TABLE_COLLATION AS c FROM information_schema.TABLES"
            " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_TYPE = 'BASE TABLE'").fetchall()}
        assert collations == {"utf8mb4_unicode_ci"}

        loaded = store.load_sector(conn, sector_id)
        assert loaded.name == "Véga Reach"
        assert sorted(e.star_system.name for e in loaded.entries) == ["Alpha Vega", "Beta Véga", "Æron"]
        assert [s["name"] for s in query.list_sectors(conn)] == ["Véga Reach"]
        assert query.sector_detail(conn, sector_id) is not None
        for row in conn.execute("SELECT id, name FROM star_systems").fetchall():
            assert query.system_detail(conn, row["id"])["name"] == row["name"]
        assert store.name_in_use(conn, "VEGA REACH") == "sectors"
        assert store.name_in_use(conn, "beta vega") == "star_systems"

        tags = {facet: set() for facet in query.SEARCH_TAG_FACETS}
        found = query.search(conn, {"system_q": "vega", "sector_q": "véga"}, tags)["results"]
        assert sorted(r["name"] for r in found["systems"]["rows"]) == ["Alpha Vega", "Beta Véga"]
        assert [r["name"] for r in found["sectors"]["rows"]] == ["Véga Reach"]
    finally:
        conn.close()

    monkeypatch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path / "tiles"))
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    client = make_app(config).test_client()
    page = client.get("/search?system_q=V%C3%A9ga&sector_q=vega")
    assert page.status_code == 200
    body = page.get_data(as_text=True)
    assert "Alpha Vega" in body and "Beta Véga" in body
    for path in ("/", "/sectors", "/systems", f"/sector/{sector_id}", "/phenomena", "/api/search?system_q=VEGA"):
        assert client.get(path).status_code == 200, path
