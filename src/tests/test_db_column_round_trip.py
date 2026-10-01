# tests/test_db_column_round_trip.py

"""
TEST.11: every column of every content table is written by the code that
saves its rows and read back by the code that loads them. A small galaxy
is built through the real `generate.py` commands (`rich_galaxy_support`),
then every column in information_schema must hold a value in some row
(so a new column the INSERT forgot, left NULL, fails), and must be read
by the loaders and the API's detail endpoints (so a column nothing ever
loads fails). The exceptions below each say why.
"""

import pytest

from stellarObjects import _db
from tests.column_read_support import columns_read, tracking_reads
from tests.rich_galaxy_support import rich_galaxy  # noqa: F401  (fixture)

pytestmark = pytest.mark.db

BOOKKEEPING_TABLES = {
    "schema_migrations", "id_blocks", "system_name_registry", "sector_name_registry", "population_state",
    "orbit_simulation_state", "galaxy_column", "galaxy_layer", "bright_star_blocks", "bright_stars",
    "nearest_systems",
}
"""Tables no object loader reads row by row: migration and id bookkeeping,
the name registries, generation progress, and the galaxy skeleton and
bright-star scatter (read in bulk by the Galaxy Map's tile queries,
covered by test_galaxymap3d.py/test_bright_star_scatter.py)."""

CHANCE = "chance"
"""Reason for a column whose value in this galaxy depends on the draw."""

CONTAINMENT = ("inside_nebula_id", "inside_remnant_id")
SECTOR_PLACEMENT = ("sector_id", "center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc", "quadrant")

NULL_IN_THIS_GALAXY = {
    # Only something sitting inside a nebula or supernova remnant.
    **{(table, column): CHANCE for table in _db.CONTAINABLE_TABLES for column in CONTAINMENT},
    # A quasar only exists placed at a galaxy's own center; the CLI's is standalone.
    **{("quasars", column): "standalone quasar" for column in SECTOR_PLACEMENT + ("jet_length_ly",)},
    # The CLI can't place a supernova remnant, and whether it has a core is chance.
    **{("supernova_remnants", column): "standalone remnant" for column in SECTOR_PLACEMENT},
    **{("supernova_remnants", column): CHANCE for column in
       ("compact_remnant_kind", "compact_remnant_black_hole_id", "compact_remnant_neutron_star_id")},
    # `plan --no-bright-stars`: the scatter is slow and tested on its own.
    ("galaxy_shape", "bright_star_min_luminosity_sol"): "no scatter", ("galaxy_shape", "bright_star_seed"): "no scatter",
    # Facilities on stars, moons and asteroid fields, and the field's galaxy position.
    **{("facilities", column): "facility hosts used" for column in
       ("star_id", "moon_id", "asteroid_field_id", "sector_id", "center_x_pc", "center_y_pc", "center_z_pc",
        "galactic_radius_pc")},
    # A life-bearing moon (rare) and its paragraphs and spectrum.
    **{("moons", column): CHANCE for column in
       ("atm_density", "atm_molar_density", "scale_height_km", "life_chemical")},
    **{("moon_evolutionary_paragraphs", column): CHANCE for column in
       ("id", "moon_id", "position", "paragraph")},
    **{("moon_reflection_spectrum", column): CHANCE for column in
       ("id", "moon_id", "spectrum_type", "position", "value")},
    # The CLI places an interstellar comet near the galaxy's sectors only sometimes.
    **{("interstellar_comets", column): CHANCE for column in SECTOR_PLACEMENT},
    ("species", "civilization_age_years"): "no civilization", ("species", "era"): "no civilization",
    ("star_systems", "system_flavor_text"): CHANCE,
    ("star_systems", "runaway_class"): CHANCE, ("star_systems", "runaway_speed_kms"): CHANCE,
    ("sectors", "wiki_url"): "never published", ("star_systems", "wikijs_url"): "never published",
    ("star_systems", "mediawiki_url"): "never published",
    ("system_configs", "age"): "only with --age",
}
"""`(table, column)` allowed to be NULL in every row of the seeded galaxy,
with why: each is set only by something this small galaxy doesn't
contain (the cases are exercised in their own tests)."""

NEVER_READ = {
    # Row order in a child table: the loaders read them ORDER BY position.
    **{(table, "position"): "sort key" for table in (
        "asteroid_belt_composition", "asteroid_field_composition", "comet_composition",
        "interstellar_comet_composition", "planet_evolutionary_paragraphs", "planet_reflection_spectrum",
        "moon_evolutionary_paragraphs", "moon_reflection_spectrum")},
    # Child tables' own ids and parent keys: matched in the WHERE clause.
    **{(table, "id"): "child row id" for table in (
        "asteroid_belt_composition", "asteroid_field_composition", "comet_composition",
        "interstellar_comet_composition", "planet_evolutionary_paragraphs", "planet_reflection_spectrum",
        "moon_evolutionary_paragraphs", "moon_reflection_spectrum", "system_config_slots", "system_configs")},
    **{(table, parent): "parent key" for table, parent in (
        ("comet_composition", "comet_id"), ("interstellar_comet_composition", "comet_id"),
        ("asteroid_field_composition", "field_id"), ("planet_evolutionary_paragraphs", "planet_id"),
        ("planet_reflection_spectrum", "planet_id"), ("moon_evolutionary_paragraphs", "moon_id"),
        ("moon_reflection_spectrum", "moon_id"), ("system_config_slots", "config_id"),
        ("system_owners", "star_system_id"))},
    ("star_systems", "schema_version"): "the version that wrote the row, for diagnosis",
    ("galaxy_shape", "id"): "singleton key", ("galaxy_shape", "bright_star_seed"): "only to repeat a scatter",
    ("facilities", "galactic_radius_pc"): "an index column; pages place a facility by its center",
    # Itemized copies of the parent's `composition_summary`, which the
    # pages show; nothing reads the items back yet.
    ("asteroid_field_composition", "component"): "itemized summary",
    ("asteroid_field_composition", "concentration"): "itemized summary",
    ("interstellar_comet_composition", "component"): "itemized summary",
    # Empty unless a moon bears life (see NULL_IN_THIS_GALAXY).
    ("moon_evolutionary_paragraphs", "paragraph"): CHANCE,
    ("moon_reflection_spectrum", "spectrum_type"): CHANCE, ("moon_reflection_spectrum", "value"): CHANCE,
    **{(table, stamp): "row timestamps (the Galaxy Map's change feed reads modified_at in bulk)"
       for table in ("star_systems", "facilities") for stamp in ("created_at", "modified_at")},
    ("polities", "created_at"): "row timestamp", ("species", "created_at"): "row timestamp",
}
"""`(table, column)` the loaders and detail endpoints don't read, with why."""


def _client(config):
    from api.app import create_app
    from api.config import Config

    class TestConfig(Config):
        MYSQL_CONFIG = config
        WRITE_MYSQL_CONFIG = config
        CONTROL_MYSQL_CONFIG = config
        RATELIMIT_ENABLED = False

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def _content_columns(conn):
    rows = conn.execute(
        "SELECT c.TABLE_NAME AS t, c.COLUMN_NAME AS c FROM information_schema.COLUMNS c"
        " JOIN information_schema.TABLES t ON t.TABLE_SCHEMA = c.TABLE_SCHEMA AND t.TABLE_NAME = c.TABLE_NAME"
        " WHERE c.TABLE_SCHEMA = DATABASE() AND t.TABLE_TYPE = 'BASE TABLE'"
        " ORDER BY c.TABLE_NAME, c.ORDINAL_POSITION"
    ).fetchall()
    return [(row["t"], row["c"]) for row in rows if row["t"] not in BOOKKEEPING_TABLES]


def test_every_column_is_written(rich_galaxy):
    conn = _db.get_connection(rich_galaxy)
    try:
        always_null = [
            (table, column) for table, column in _content_columns(conn)
            if conn.execute(f"SELECT 1 FROM `{table}` WHERE `{column}` IS NOT NULL LIMIT 1").fetchone() is None
        ]
    finally:
        conn.close()
    unexpected = sorted(set(always_null) - set(NULL_IN_THIS_GALAXY))
    assert not unexpected, f"NULL in every row (never written?): {unexpected}"
    stale = sorted(key for key, why in NULL_IN_THIS_GALAXY.items() if why != CHANCE and key not in always_null)
    assert not stale, f"written after all; drop from NULL_IN_THIS_GALAXY: {stale}"


def test_every_column_is_read_back(rich_galaxy, monkeypatch):
    client = _client(rich_galaxy)
    conn = _db.get_connection(rich_galaxy)
    try:
        def ids(table):
            return [row["id"] for row in conn.execute(f"SELECT id FROM {table} ORDER BY id").fetchall()]

        import queryDb

        with tracking_reads(monkeypatch) as reads:
            pages = [f"/api/sectors/{i}" for i in ids("sectors")]
            pages += [f"/api/systems/{i}{part}" for i in ids("star_systems")
                      for part in ("", "/sections", "/owner", "/facilities")]
            pages += [f"/api/phenomena/{label}/{i}" for table, label, *_ in queryDb._PHENOMENON_TABLES
                      for i in ids(table)]
            pages += [f"/api/facilities/{i}" for i in ids("facilities")]
            pages += [f"/api/species/{i}" for i in ids("species")]
            pages += [f"/api/polities/{i}" for i in ids("polities")]
            pages += ["/api/galaxy/shape", "/api/population", "/api/territories"]
            for page in pages:
                assert client.get(page).status_code == 200, page
            for sector_id in ids("sectors"):
                _db.load_sector(conn, sector_id)
            for system_id in ids("star_systems"):
                _db.load_star_system(conn, system_id)
            columns = _content_columns(conn)
            read = columns_read(reads, {table for table, _column in columns})
    finally:
        conn.close()
    unread = sorted((table, column) for table, column in columns if column not in read[table])
    unexpected = sorted(set(unread) - set(NEVER_READ))
    assert not unexpected, f"never read back: {unexpected}"
    stale = sorted(key for key, why in NEVER_READ.items() if why != CHANCE and key not in unread)
    assert not stale, f"read after all; drop from NEVER_READ: {stale}"
