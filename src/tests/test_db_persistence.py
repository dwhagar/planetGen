# tests/test_db_persistence.py

"""
End-to-end save-path tests for `stellarObjects._db.insert_star_system`,
specifically that a freshly generated system's planets and moons land in
their respective schema-v2 tables (`planets` vs. `moons`, see
`schema.sql`'s "v2" header note) rather than sharing one table the way
schema v1 did. There's otherwise no test coverage of `_db.py`'s save path
at all (`docs/database-schema.md`'s own "no read path yet" status,
historical, since resolved), so this is deliberately a real
generate-then-save-then-query round trip rather than a unit test against
a hand-built object graph -- it's the same shape of bug (an `INSERT`'s
column list and value tuple silently drifting out of count/order) that
unit tests calling `insert_moon` directly with hand-picked arguments
would be less likely to catch.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
they're skipped, not failed, when no MySQL test server is configured/
reachable.
"""

import math

import pytest

from stellarObjects import _db
from stellarObjects import physical_constants as pc
from stellarObjects.config import SystemConfig
from stellarObjects.doubleStar import BinaryStarProxy
from stellarObjects.planetPhysics import calculate_orbital_period_years
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem
from stellarObjects.utils import circular_orbital_speed_kms, minimum_update_interval_years, orbital_position_au


def _drop_v17_phenomenon_columns(conn):
    """
    Drops the v17 `galactic_orbital_*` columns from the six pre-existing
    exotic-phenomenon tables (`black_holes`, `neutron_stars`, `nebulae`,
    `supernova_remnants`, `rogue_planets`, `interstellar_comets`) -- shared
    by every `test_migrate_vN_to_vN+1_*` test below that simulates a
    database older than v16 (the version these tables were introduced in).
    `mysql_config` always bootstraps a fresh database at the CURRENT
    schema (today, v17), so those tables and their v17 columns already
    exist even though a real pre-v16 database would have neither -- left
    undropped, `migrate_database`'s replayed `_migrate_v16_to_v17` would
    try to `ADD COLUMN` ones that already exist and fail with a duplicate-
    column error. See `schema.sql`'s "v17" header note.
    """
    for table in ("black_holes", "neutron_stars", "nebulae", "supernova_remnants",
                  "rogue_planets", "interstellar_comets"):
        conn.execute(
            f"ALTER TABLE {table} "
            "DROP COLUMN galactic_orbital_speed_kms, DROP COLUMN galactic_orbital_period_gy, "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years"
        )


def _drop_v18_phenomenon_columns(conn):
    """
    Drops the v18 galaxy-frame placement columns from `nebulae`/
    `asteroid_fields` -- same "already exists on a freshly-bootstrapped
    test database, so migrate_database's replayed _migrate_v17_to_v18
    would hit a duplicate-column error" reasoning
    `_drop_v17_phenomenon_columns` gives for its own six tables. See
    `schema.sql`'s "v18" header note.

    The named `chk_{table}_placement` CHECK constraint must be dropped
    first (`DROP CONSTRAINT`, portable across MySQL 8.0.19+/MariaDB --
    MySQL's own `DROP CHECK` alias is not): MySQL refuses `DROP COLUMN` on
    a column a CHECK still references (error 3959), and a fresh test
    database already has this constraint (`schema.sql`'s own
    `CREATE TABLE` body), unlike a real pre-v18 database, which never had
    it either.
    """
    for table in ("nebulae", "asteroid_fields"):
        conn.execute(
            f"ALTER TABLE {table} "
            f"DROP CONSTRAINT chk_{table}_placement, "
            f"DROP INDEX idx_{table}_galactic_radius_pc, "
            "DROP COLUMN center_x_pc, DROP COLUMN center_y_pc, "
            "DROP COLUMN center_z_pc, DROP COLUMN galactic_radius_pc"
        )


def _make_system_with_moons_and_belt():
    """Retries generation (bounded) until a system with at least one moon
    and one asteroid belt comes out -- MOONS/ASTEROID_BELT=True bias
    generation heavily toward both, but neither is unconditionally
    guaranteed for every planet/system."""
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        cfg.ASTEROID_BELT = True
        # Pinned False: left at its own default, `BINARY_SYSTEM` now rolls
        # real chance (StarSystem._should_generate_binary), and this
        # helper's only caller counts expected rows from `system.planets`
        # alone -- an S-type (wide) pair's own `secondary_planets` would
        # add real DB rows that count wouldn't account for.
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        planets = [obj for obj in system.planets if obj.body_type != "a"]
        belts = [obj for obj in system.planets if obj.body_type == "a"]
        if belts and any(p.moons for p in planets):
            return system, cfg
    pytest.fail("could not generate a system with both a moon and an asteroid belt")


def test_insert_star_system_splits_planets_and_moons_into_their_own_tables(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
    finally:
        conn.close()

    expected_planets = [obj for obj in system.planets if obj.body_type != "a"]
    expected_belts = [obj for obj in system.planets if obj.body_type == "a"]
    expected_moon_names = sorted(moon.name for planet in expected_planets for moon in planet.moons)
    assert expected_moon_names, "test fixture must actually contain moons"

    conn = _db.get_connection(mysql_config)
    try:
        db_planets = conn.execute("SELECT * FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert len(db_planets) == len(expected_planets)
        assert "is_moon" not in db_planets[0].keys()
        assert "parent_planet_id" not in db_planets[0].keys()

        db_moons = conn.execute("SELECT * FROM moons WHERE star_system_id = ?", (system_id,)).fetchall()
        assert sorted(row["name"] for row in db_moons) == expected_moon_names

        planet_ids_by_name = {row["name"]: row["id"] for row in db_planets}
        for planet in expected_planets:
            moon_rows = [row for row in db_moons if row["planet_id"] == planet_ids_by_name[planet.name]]
            assert sorted(row["name"] for row in moon_rows) == sorted(m.name for m in planet.moons)

        belts = conn.execute("SELECT * FROM asteroid_belts WHERE star_system_id = ?", (system_id,)).fetchall()
        assert len(belts) == len(expected_belts)

        version_row = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()
        assert version_row["version"] == _db.SCHEMA_VERSION
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Read path (_db.load_star_system / load_sector / load_system_config).
# Covers the gaps this module's own docstring called out as still missing:
# binary systems, lifespan_gy round-tripping, and foreign-key integrity,
# plus save_sector/insert_sector and insert_system_config's SLOTS child rows.
# ---------------------------------------------------------------------------

def test_load_star_system_round_trips_single_star_system_exactly(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert str(reloaded) == str(system)
    assert reloaded.planet_count == system.planet_count
    assert reloaded.belt_count == system.belt_count
    assert reloaded.moon_count == system.moon_count
    assert reloaded.hab_count == system.hab_count
    assert reloaded.m_count == system.m_count
    # Shared back-references: the SAME star/system_config instance everywhere.
    for obj in reloaded.planets:
        assert obj.system_config is reloaded.system_config
        if obj.body_type != "a":
            assert obj.star is reloaded.star


def test_load_star_system_round_trips_binary_system_exactly(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.BINARY_SYSTEM = True
    cfg.WIDE_BINARY = False  # this test specifically asserts BinaryStarProxy-only behavior
    cfg.PLANETS = False
    system = StarSystem(system_config=cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert isinstance(reloaded.star, BinaryStarProxy)
    assert len(reloaded.stars) == 2
    assert str(reloaded) == str(system)
    assert reloaded.star.mass == system.star.mass
    assert reloaded.star.luminosity == system.star.luminosity
    # The DB schema has exactly one system_config_id per star_systems row,
    # so the secondary star's generation-time-only deep-copied config
    # (see TODO.md's resolved "secondary-star system_config asymmetry"
    # question) is naturally collapsed to the one shared config on reload.
    assert reloaded.secondary_star.system_config is reloaded.system_config
    assert reloaded.primary_star.system_config is reloaded.system_config


def test_lifespan_gy_null_round_trips_to_infinite_lifespan(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M2VII"  # white dwarf: lifespan == float('inf')
    cfg.PLANETS = False
    system = StarSystem(system_config=cfg)
    assert system.star.lifespan == float('inf')

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        row = conn.execute(
            "SELECT lifespan_gy FROM stars WHERE star_system_id = ?", (system_id,)
        ).fetchone()
        assert row["lifespan_gy"] is None  # NULL in storage, never the JSON "Infinity" token

        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert reloaded.star.lifespan == float('inf')
    assert str(reloaded) == str(system)


def test_galactic_orbit_fields_round_trip_exactly(mysql_config):
    """
    Regression test for the v10 `galactic_orbital_speed_kms`/
    `galactic_orbital_period_gy` columns, for both a single star and a
    binary pair -- the same "INSERT column list vs. value tuple drift"
    bug class `test_orbital_motion_fields_round_trip_exactly` guards
    against for the v9 columns.
    """
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert reloaded.star.galactic_orbital_speed_kms == pytest.approx(system.star.galactic_orbital_speed_kms)
    assert reloaded.star.galactic_orbital_period_gy == pytest.approx(system.star.galactic_orbital_period_gy)

    binary_cfg = SystemConfig()
    binary_cfg.STAR_TYPE = "G2V"
    binary_cfg.BINARY_SYSTEM = True
    binary_cfg.WIDE_BINARY = False  # pin to close/P-type -- these assertions are BinaryStarProxy-specific
    binary_cfg.PLANETS = False
    binary_system = StarSystem(system_config=binary_cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            binary_system_id = _db.insert_star_system(conn, binary_system, binary_cfg)
        reloaded_binary = _db.load_star_system(conn, binary_system_id)
    finally:
        conn.close()

    assert reloaded_binary.star.galactic_orbital_speed_kms == pytest.approx(binary_system.star.galactic_orbital_speed_kms)
    assert reloaded_binary.star.galactic_orbital_period_gy == pytest.approx(binary_system.star.galactic_orbital_period_gy)


def test_star_motion_fields_round_trip_exactly(mysql_config):
    """
    Regression test for the v13 star-motion columns and the v14 binary
    mutual-orbit position columns, for both a single star and a binary
    pair -- the same "INSERT column list vs. value tuple drift" bug class
    `test_orbital_motion_fields_round_trip_exactly`/
    `test_galactic_orbit_fields_round_trip_exactly` guard against for the
    v9/v10 columns.
    """
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert reloaded.star.galactic_orbital_phase_deg == pytest.approx(system.star.galactic_orbital_phase_deg)
    assert reloaded.star.galactic_min_update_interval_years == pytest.approx(
        system.star.galactic_min_update_interval_years
    )

    binary_cfg = SystemConfig()
    binary_cfg.STAR_TYPE = "G2V"
    binary_cfg.BINARY_SYSTEM = True
    binary_cfg.WIDE_BINARY = False  # pin to close/P-type -- these assertions are BinaryStarProxy-specific
    binary_cfg.PLANETS = False
    binary_system = StarSystem(system_config=binary_cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            binary_system_id = _db.insert_star_system(conn, binary_system, binary_cfg)
        reloaded_binary = _db.load_star_system(conn, binary_system_id)
    finally:
        conn.close()

    proxy = binary_system.star
    reloaded_proxy = reloaded_binary.star
    assert reloaded_proxy.galactic_orbital_phase_deg == pytest.approx(proxy.galactic_orbital_phase_deg)
    assert reloaded_proxy.galactic_min_update_interval_years == pytest.approx(
        proxy.galactic_min_update_interval_years
    )
    # Both stars of a binary pair, and the proxy, always share one galactic phase.
    assert reloaded_binary.primary_star.galactic_orbital_phase_deg == pytest.approx(proxy.galactic_orbital_phase_deg)
    assert reloaded_binary.secondary_star.galactic_orbital_phase_deg == pytest.approx(
        proxy.galactic_orbital_phase_deg
    )

    assert reloaded_proxy.binary_mutual_orbital_period_years == pytest.approx(
        proxy.binary_mutual_orbital_period_years
    )
    assert reloaded_proxy.binary_mutual_orbital_speed_kms == pytest.approx(proxy.binary_mutual_orbital_speed_kms)
    assert reloaded_proxy.binary_mutual_orbital_inclination_deg == pytest.approx(
        proxy.binary_mutual_orbital_inclination_deg
    )
    assert reloaded_proxy.binary_mutual_orbital_ascending_node_deg == pytest.approx(
        proxy.binary_mutual_orbital_ascending_node_deg
    )
    assert reloaded_proxy.binary_mutual_orbital_phase_deg == pytest.approx(proxy.binary_mutual_orbital_phase_deg)
    assert reloaded_proxy.binary_mutual_min_update_interval_years == pytest.approx(
        proxy.binary_mutual_min_update_interval_years
    )
    assert reloaded_proxy.binary_mutual_position_x == pytest.approx(proxy.binary_mutual_position_x)
    assert reloaded_proxy.binary_mutual_position_y == pytest.approx(proxy.binary_mutual_position_y)
    assert reloaded_proxy.binary_mutual_position_z == pytest.approx(proxy.binary_mutual_position_z)

    # The stored position must actually match orbital_position_au at the
    # stored separation/inclination/ascending_node/phase, not just survive
    # a round trip unchanged.
    expected = orbital_position_au(
        proxy.binary_separation_au, proxy.binary_mutual_orbital_inclination_deg,
        proxy.binary_mutual_orbital_ascending_node_deg, proxy.binary_mutual_orbital_phase_deg,
    )
    assert (proxy.binary_mutual_position_x, proxy.binary_mutual_position_y, proxy.binary_mutual_position_z) == \
        pytest.approx(expected)


def test_insert_star_system_respects_foreign_keys(mysql_config):
    # MySQL/InnoDB enforces every foreign key eagerly, at INSERT time --
    # unlike SQLite (which needs a separate `PRAGMA foreign_key_check`
    # pass after the fact to catch a constraint left unenforced mid-
    # transaction), an insert violating a foreign key here would already
    # have raised `pymysql.err.IntegrityError` above rather than
    # completing silently. This test's real assertion is simply that the
    # insert-then-reload round trip above completes without one.
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_star_system(conn, system, cfg)
    finally:
        conn.close()


def test_save_sector_and_load_sector_round_trip(mysql_config):
    sector = SpaceSector("Persisted Sector", edge_ly=20.0)

    system_a, cfg_a = _make_system_with_moons_and_belt()
    sector.add_system(system_a, position=(1.0, 2.0, 3.0), system_config=cfg_a)

    cfg_b = SystemConfig()
    cfg_b.STAR_TYPE = "M5V"
    system_b = StarSystem(system_config=cfg_b)
    sector.add_system(system_b, position=(-4.0, 0.0, 5.5), system_config=cfg_b)

    sector_id = _db.save_sector(sector, config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        reloaded_sector = _db.load_sector(conn, sector_id)
    finally:
        conn.close()

    assert reloaded_sector.name == "Persisted Sector"
    assert reloaded_sector.edge_ly == pytest.approx(20.0)
    assert len(reloaded_sector) == 2

    reloaded_by_name = {entry.star_system.star.name: entry for entry in reloaded_sector.entries}
    for original_system, original_position in ((system_a, (1.0, 2.0, 3.0)), (system_b, (-4.0, 0.0, 5.5))):
        entry = reloaded_by_name[original_system.star.name]
        assert str(entry.star_system) == str(original_system)
        # Round-tripped through ly -> milliparsecs -> ly, so exact equality
        # isn't expected -- only floating-point-close.
        assert entry.position == pytest.approx(original_position)


def test_orbital_motion_fields_round_trip_exactly(mysql_config):
    """
    Regression test for the exact bug class this module's own docstring
    warns about (an `INSERT`'s column list and value tuple silently
    drifting out of count): `insert_moon`'s `VALUES` clause was originally
    missing one placeholder relative to its column list when the v9
    orbital-motion columns were added, which pymysql surfaced as a generic
    "not all arguments converted during string formatting" `TypeError` --
    a symptom easy to mistake for something else. Confirms every
    orbital-motion field survives a save/load round trip exactly, for both
    planets and moons.
    """
    system, cfg = _make_system_with_moons_and_belt()
    planets = [obj for obj in system.planets if obj.body_type != "a"]
    assert any(p.moons for p in planets), "test fixture must actually contain moons"

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        reloaded = _db.load_star_system(conn, system_id)
    finally:
        conn.close()

    reloaded_planets_by_name = {p.name: p for p in reloaded.planets if p.body_type != "a"}
    for planet in planets:
        reloaded_planet = reloaded_planets_by_name[planet.name]
        assert reloaded_planet.orbital_inclination_deg == pytest.approx(planet.orbital_inclination_deg)
        assert reloaded_planet.orbital_ascending_node_deg == pytest.approx(planet.orbital_ascending_node_deg)
        assert reloaded_planet.orbital_phase_deg == pytest.approx(planet.orbital_phase_deg)
        assert reloaded_planet.position_x == pytest.approx(planet.position_x)
        assert reloaded_planet.position_y == pytest.approx(planet.position_y)
        assert reloaded_planet.position_z == pytest.approx(planet.position_z)
        assert reloaded_planet.orbital_speed_kms == pytest.approx(planet.orbital_speed_kms)
        assert reloaded_planet.min_update_interval_years == pytest.approx(planet.min_update_interval_years)
        assert reloaded_planet.rotation_period_hours == pytest.approx(planet.rotation_period_hours)

        reloaded_moons_by_name = {m.name: m for m in reloaded_planet.moons}
        for moon in planet.moons:
            reloaded_moon = reloaded_moons_by_name[moon.name]
            assert reloaded_moon.orbital_inclination_deg == pytest.approx(moon.orbital_inclination_deg)
            assert reloaded_moon.orbital_ascending_node_deg == pytest.approx(moon.orbital_ascending_node_deg)
            assert reloaded_moon.orbital_phase_deg == pytest.approx(moon.orbital_phase_deg)
            assert reloaded_moon.position_x == pytest.approx(moon.position_x)
            assert reloaded_moon.position_y == pytest.approx(moon.position_y)
            assert reloaded_moon.position_z == pytest.approx(moon.position_z)
            assert reloaded_moon.orbital_speed_kms == pytest.approx(moon.orbital_speed_kms)
            assert reloaded_moon.min_update_interval_years == pytest.approx(moon.min_update_interval_years)
            assert reloaded_moon.rotation_period_hours == pytest.approx(moon.rotation_period_hours)


# ---------------------------------------------------------------------------
# Orbital motion updates (stellarObjects._db.advance_orbital_phases /
# get_orbit_update_elapsed_years) -- see updateOrbits.py.
# ---------------------------------------------------------------------------

def test_get_orbit_update_elapsed_years_is_none_before_first_update(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        assert _db.get_orbit_update_elapsed_years(conn) is None
    finally:
        conn.close()


def test_advance_orbital_phases_applies_the_correct_delta_and_leaves_other_fields_untouched(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_star_system(conn, system, cfg)

        before = conn.execute(
            "SELECT id, orbital_phase_deg, orbital_inclination_deg, orbital_ascending_node_deg, "
            "rotation_period_hours, period_years, distance_km, "
            "position_x_km, position_y_km, position_z_km, orbital_speed_kms FROM planets"
        ).fetchall()

        counts = _db.advance_orbital_phases(conn, elapsed_years=0.5)
        assert counts["planets"] > 0
        assert counts["moons"] > 0
        assert counts["stars"] > 0

        after = conn.execute(
            "SELECT id, orbital_phase_deg, orbital_inclination_deg, orbital_ascending_node_deg, "
            "rotation_period_hours, period_years, distance_km, "
            "position_x_km, position_y_km, position_z_km, orbital_speed_kms FROM planets"
        ).fetchall()
    finally:
        conn.close()

    after_by_id = {row["id"]: row for row in after}
    for row in before:
        updated = after_by_id[row["id"]]
        # Orientation/rotation are fixed at generation time -- only phase
        # (and position, in lockstep with it) moves. Speed is constant
        # around a circular orbit.
        assert updated["orbital_inclination_deg"] == pytest.approx(row["orbital_inclination_deg"])
        assert updated["orbital_ascending_node_deg"] == pytest.approx(row["orbital_ascending_node_deg"])
        assert updated["rotation_period_hours"] == pytest.approx(row["rotation_period_hours"])
        assert updated["orbital_speed_kms"] == pytest.approx(row["orbital_speed_kms"])

        expected_phase = (row["orbital_phase_deg"] + (0.5 / row["period_years"]) * 360) % 360
        assert updated["orbital_phase_deg"] == pytest.approx(expected_phase, abs=1e-6)

        expected_x_au, expected_y_au, expected_z_au = orbital_position_au(
            row["distance_km"] / pc.AU_TO_KM, row["orbital_inclination_deg"],
            row["orbital_ascending_node_deg"], expected_phase,
        )
        assert updated["position_x_km"] == pytest.approx(expected_x_au * pc.AU_TO_KM, abs=1e-3)
        assert updated["position_y_km"] == pytest.approx(expected_y_au * pc.AU_TO_KM, abs=1e-3)
        assert updated["position_z_km"] == pytest.approx(expected_z_au * pc.AU_TO_KM, abs=1e-3)

    # get_orbit_update_elapsed_years should now report ~0 elapsed time
    # (the call above just set last_updated_at to NOW()), not None.
    conn = _db.get_connection(mysql_config)
    try:
        elapsed = _db.get_orbit_update_elapsed_years(conn)
    finally:
        conn.close()
    assert elapsed is not None
    assert 0 <= elapsed < 0.01


def test_advance_orbital_phases_advances_binary_mutual_orbit_position_in_lockstep(mysql_config):
    """
    Regression test for the v14 addition: `advance_orbital_phases` must
    move `binary_mutual_orbital_phase_deg` *and* recompute
    `binary_mutual_position_x/y/z_km` from that new phase in the same
    pass -- the same "position has no independent update of its own"
    treatment `test_advance_orbital_phases_applies_the_correct_delta_and_leaves_other_fields_untouched`
    already verifies for planets/moons.
    """
    binary_cfg = SystemConfig()
    binary_cfg.STAR_TYPE = "G2V"
    binary_cfg.BINARY_SYSTEM = True
    binary_cfg.WIDE_BINARY = False  # pin to close/P-type -- these assertions are BinaryStarProxy-specific
    binary_cfg.PLANETS = False
    binary_system = StarSystem(system_config=binary_cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            binary_id = _db.insert_star_system(conn, binary_system, binary_cfg)

        before = conn.execute(
            "SELECT binary_separation_km, binary_mutual_orbital_period_years, "
            "binary_mutual_orbital_inclination_deg, binary_mutual_orbital_ascending_node_deg, "
            "binary_mutual_orbital_phase_deg "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()

        # Big enough elapsed time to clearly move the (fast) mutual orbit,
        # regardless of its randomly generated period.
        elapsed_years = before["binary_mutual_orbital_period_years"] * 137.25
        _db.advance_orbital_phases(conn, elapsed_years=elapsed_years)

        after = conn.execute(
            "SELECT binary_mutual_orbital_phase_deg, binary_mutual_position_x_km, "
            "binary_mutual_position_y_km, binary_mutual_position_z_km "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()
    finally:
        conn.close()

    expected_phase = (
        before["binary_mutual_orbital_phase_deg"]
        + (elapsed_years / before["binary_mutual_orbital_period_years"]) * 360
    ) % 360
    assert after["binary_mutual_orbital_phase_deg"] == pytest.approx(expected_phase, abs=1e-6)
    assert after["binary_mutual_orbital_phase_deg"] != pytest.approx(before["binary_mutual_orbital_phase_deg"])

    expected_x_au, expected_y_au, expected_z_au = orbital_position_au(
        before["binary_separation_km"] / pc.AU_TO_KM, before["binary_mutual_orbital_inclination_deg"],
        before["binary_mutual_orbital_ascending_node_deg"], expected_phase,
    )
    assert after["binary_mutual_position_x_km"] == pytest.approx(expected_x_au * pc.AU_TO_KM, abs=1e-3)
    assert after["binary_mutual_position_y_km"] == pytest.approx(expected_y_au * pc.AU_TO_KM, abs=1e-3)
    assert after["binary_mutual_position_z_km"] == pytest.approx(expected_z_au * pc.AU_TO_KM, abs=1e-3)


def test_advance_orbital_phases_rejects_negative_elapsed_years(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError):
            _db.advance_orbital_phases(conn, elapsed_years=-1.0)
    finally:
        conn.close()


def test_migrate_v8_to_v9_adds_orbital_motion_columns(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing v8
    database by tearing the v9 through v15 additions back out (a
    real v8 database would have none of them), then confirms
    `migrate_database` puts everything back and reports the database as
    fully current again (not just v9 -- `migrate_database` applies every
    step up to `SCHEMA_VERSION` in one call, so a v8 database now lands on
    v15 directly).
    """
    conn = _db.get_connection(mysql_config)
    try:
        for table in ("planets", "moons"):
            conn.execute(
                f"ALTER TABLE {table} "
                f"DROP COLUMN orbital_inclination_deg, DROP COLUMN orbital_ascending_node_deg, "
                f"DROP COLUMN orbital_phase_deg, "
                f"DROP COLUMN position_x_km, DROP COLUMN position_y_km, DROP COLUMN position_z_km, "
                f"DROP COLUMN orbital_speed_kms, DROP COLUMN min_update_interval_years, "
                f"DROP COLUMN rotation_period_hours"
            )
        conn.execute("DROP TABLE orbit_simulation_state")
        conn.execute(
            "ALTER TABLE stars "
            "DROP COLUMN galactic_orbital_speed_kms, DROP COLUMN galactic_orbital_period_gy, "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years, "
            "DROP COLUMN wide_binary_a_crit_km"
        )
        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_galactic_orbital_speed_kms, DROP COLUMN binary_galactic_orbital_period_gy, "
            "DROP COLUMN binary_galactic_orbital_phase_deg, DROP COLUMN binary_galactic_min_update_interval_years, "
            "DROP COLUMN binary_mutual_orbital_period_years, DROP COLUMN binary_mutual_orbital_speed_kms, "
            "DROP COLUMN binary_mutual_orbital_inclination_deg, DROP COLUMN binary_mutual_orbital_ascending_node_deg, "
            "DROP COLUMN binary_mutual_orbital_phase_deg, DROP COLUMN binary_mutual_min_update_interval_years, "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (9, 10, 11, 12, 13, 14, 15, 16, 17, 18)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (8)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 8
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION == 18

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        planet_columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM planets").fetchall()}
        assert {
            "orbital_inclination_deg", "orbital_ascending_node_deg",
            "orbital_phase_deg", "rotation_period_hours",
            "position_x_km", "position_y_km", "position_z_km", "orbital_speed_kms",
            "min_update_interval_years",
        } <= planet_columns
        assert _db.get_orbit_update_elapsed_years(conn) is None  # table exists, no row yet

        star_columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM stars").fetchall()}
        assert {
            "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
            "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
        } <= star_columns
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == 18


def test_migrate_v9_to_v10_adds_galactic_orbit_columns(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing v9
    database by tearing the v10 through v15 additions back out (the new
    `stars`/`star_systems`/`planets`/`moons`/`asteroid_belts` columns and the
    `schema_migrations` v10/v11/v12/v13/v14/v15 rows) before calling
    `migrate_database`, and confirms it puts everything back and reports
    the database as current again, without disturbing the v9
    orbital-motion columns already in place.
    """
    conn = _db.get_connection(mysql_config)
    try:
        conn.execute(
            "ALTER TABLE stars "
            "DROP COLUMN galactic_orbital_speed_kms, DROP COLUMN galactic_orbital_period_gy, "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years, "
            "DROP COLUMN wide_binary_a_crit_km"
        )
        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_galactic_orbital_speed_kms, DROP COLUMN binary_galactic_orbital_period_gy, "
            "DROP COLUMN binary_galactic_orbital_phase_deg, DROP COLUMN binary_galactic_min_update_interval_years, "
            "DROP COLUMN binary_mutual_orbital_period_years, DROP COLUMN binary_mutual_orbital_speed_kms, "
            "DROP COLUMN binary_mutual_orbital_inclination_deg, DROP COLUMN binary_mutual_orbital_ascending_node_deg, "
            "DROP COLUMN binary_mutual_orbital_phase_deg, DROP COLUMN binary_mutual_min_update_interval_years, "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        for table in ("planets", "moons"):
            conn.execute(
                f"ALTER TABLE {table} "
                f"DROP COLUMN position_x_km, DROP COLUMN position_y_km, DROP COLUMN position_z_km, "
                f"DROP COLUMN orbital_speed_kms, DROP COLUMN min_update_interval_years"
            )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (10, 11, 12, 13, 14, 15, 16, 17, 18)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (9)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 9
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION == 18

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        star_columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM stars").fetchall()}
        assert {
            "galactic_orbital_speed_kms", "galactic_orbital_period_gy",
            "galactic_orbital_phase_deg", "galactic_min_update_interval_years",
        } <= star_columns

        star_system_columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM star_systems").fetchall()}
        assert {
            "binary_galactic_orbital_speed_kms", "binary_galactic_orbital_period_gy",
            "binary_galactic_orbital_phase_deg", "binary_galactic_min_update_interval_years",
            "binary_mutual_orbital_period_years", "binary_mutual_orbital_speed_kms",
            "binary_mutual_orbital_inclination_deg", "binary_mutual_orbital_ascending_node_deg",
            "binary_mutual_orbital_phase_deg", "binary_mutual_min_update_interval_years",
            "binary_mutual_position_x_km", "binary_mutual_position_y_km", "binary_mutual_position_z_km",
        } <= star_system_columns

        planet_columns = {row["Field"] for row in conn.execute("SHOW COLUMNS FROM planets").fetchall()}
        assert {
            "position_x_km", "position_y_km", "position_z_km", "orbital_speed_kms",
            "min_update_interval_years",
        } <= planet_columns
        # v9's orbital-motion columns are untouched by these migration steps.
        assert "orbital_phase_deg" in planet_columns
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == 18


def test_migrate_v10_to_v11_adds_and_backfills_position_columns(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing
    v10 database by tearing the v11 through v15 additions back out (the new
    `planets`/`moons`/`stars`/`star_systems`/`asteroid_belts` columns and the
    `schema_migrations` v11/v12/v13/v14/v15 rows) before calling
    `migrate_database`, and confirms it not only adds the columns back but
    backfills them with real derived values (not an arbitrary placeholder)
    from each row's own already-stored `distance_km`/`period_years`/
    orbital-motion columns -- see `_migrate_v10_to_v11`'s docstring.
    """
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)

        for table in ("planets", "moons"):
            conn.execute(
                f"ALTER TABLE {table} "
                f"DROP COLUMN position_x_km, DROP COLUMN position_y_km, DROP COLUMN position_z_km, "
                f"DROP COLUMN orbital_speed_kms, DROP COLUMN min_update_interval_years"
            )
        # A "v10" database also predates v13/v14/v15's stars/star_systems/asteroid_belts columns.
        conn.execute(
            "ALTER TABLE stars "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years, "
            "DROP COLUMN wide_binary_a_crit_km"
        )
        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_galactic_orbital_phase_deg, DROP COLUMN binary_galactic_min_update_interval_years, "
            "DROP COLUMN binary_mutual_orbital_period_years, DROP COLUMN binary_mutual_orbital_speed_kms, "
            "DROP COLUMN binary_mutual_orbital_inclination_deg, DROP COLUMN binary_mutual_orbital_ascending_node_deg, "
            "DROP COLUMN binary_mutual_orbital_phase_deg, DROP COLUMN binary_mutual_min_update_interval_years, "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (11, 12, 13, 14, 15, 16, 17, 18)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (10)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 10
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION == 18

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute(
            "SELECT distance_km, orbital_inclination_deg, orbital_ascending_node_deg, orbital_phase_deg, "
            "period_years, position_x_km, position_y_km, position_z_km, orbital_speed_kms "
            "FROM planets WHERE star_system_id = ?", (system_id,)
        ).fetchall()
        assert rows, "test fixture must actually contain planets"
        for row in rows:
            expected_x_au, expected_y_au, expected_z_au = orbital_position_au(
                row["distance_km"] / pc.AU_TO_KM, row["orbital_inclination_deg"],
                row["orbital_ascending_node_deg"], row["orbital_phase_deg"],
            )
            assert row["position_x_km"] == pytest.approx(expected_x_au * pc.AU_TO_KM, abs=1e-3)
            assert row["position_y_km"] == pytest.approx(expected_y_au * pc.AU_TO_KM, abs=1e-3)
            assert row["position_z_km"] == pytest.approx(expected_z_au * pc.AU_TO_KM, abs=1e-3)

            expected_speed = (2 * math.pi * row["distance_km"]) / (row["period_years"] * pc.SECONDS_PER_YEAR)
            assert row["orbital_speed_kms"] == pytest.approx(expected_speed, rel=1e-9)
    finally:
        conn.close()


def test_migrate_v11_to_v12_adds_and_backfills_min_update_interval_years(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing
    v11 database by tearing the v12 through v15 additions back out (the
    new `planets`/`moons.min_update_interval_years` column, `stars`/
    `star_systems`/`asteroid_belts`' v13/v14/v15 columns, and the
    `schema_migrations` v12/v13/v14/v15 rows) before calling
    `migrate_database`, and confirms it not only adds the column back but
    backfills it with a real derived value (not an arbitrary placeholder)
    from each row's own already-stored `period_years` -- see
    `_migrate_v11_to_v12`'s docstring. Scoped to `planets`/`moons` only for
    v12: `stars`/`star_systems` didn't gain a column there (see
    `_migrate_v11_to_v12`'s docstring for why), only later at v13.
    """
    system, cfg = _make_system_with_moons_and_belt()

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)

        for table in ("planets", "moons"):
            conn.execute(f"ALTER TABLE {table} DROP COLUMN min_update_interval_years")
        conn.execute(
            "ALTER TABLE stars "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years, "
            "DROP COLUMN wide_binary_a_crit_km"
        )
        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_galactic_orbital_phase_deg, DROP COLUMN binary_galactic_min_update_interval_years, "
            "DROP COLUMN binary_mutual_orbital_period_years, DROP COLUMN binary_mutual_orbital_speed_kms, "
            "DROP COLUMN binary_mutual_orbital_inclination_deg, DROP COLUMN binary_mutual_orbital_ascending_node_deg, "
            "DROP COLUMN binary_mutual_orbital_phase_deg, DROP COLUMN binary_mutual_min_update_interval_years, "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (12, 13, 14, 15, 16, 17, 18)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (11)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 11
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION == 18

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        planet_rows = conn.execute(
            "SELECT period_years, min_update_interval_years "
            "FROM planets WHERE star_system_id = ?", (system_id,)
        ).fetchall()
        assert planet_rows, "test fixture must actually contain planets"
        for row in planet_rows:
            expected = minimum_update_interval_years(row["period_years"])
            assert row["min_update_interval_years"] == pytest.approx(expected, rel=1e-9)

        moon_rows = conn.execute(
            "SELECT period_years, min_update_interval_years "
            "FROM moons WHERE star_system_id = ?", (system_id,)
        ).fetchall()
        assert moon_rows, "test fixture must actually contain moons"
        for row in moon_rows:
            expected = minimum_update_interval_years(row["period_years"])
            assert row["min_update_interval_years"] == pytest.approx(expected, rel=1e-9)
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == 18


def test_migrate_v12_to_v13_adds_and_backfills_star_motion_columns(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing
    v12 database by tearing the v13 through v15 additions back out (the new
    `stars`/`star_systems`/`asteroid_belts` columns and the
    `schema_migrations` v13/v14/v15 rows) before calling `migrate_database`,
    and confirms it not only adds the
    columns back but backfills every real-derived one (the two
    `*_min_update_interval_years` guards, plus the binary pair's
    `binary_mutual_orbital_period_years`/`_speed_kms`) with real values
    from each row's own already-stored `galactic_orbital_period_gy`/
    `binary_separation_km`/`binary_effective_mass_kg` -- see
    `_migrate_v12_to_v13`'s docstring. The phase/orientation columns with
    no derivable "correct" value (`galactic_orbital_phase_deg` and the
    three `binary_mutual_orbital_{inclination,ascending_node,phase}_deg`
    columns) get the arbitrary `0` placeholder instead.
    """
    system, cfg = _make_system_with_moons_and_belt()
    binary_config = SystemConfig()
    binary_config.BINARY_SYSTEM = True
    binary_config.WIDE_BINARY = False  # pin to close/P-type -- these assertions are BinaryStarProxy-specific
    binary_system = StarSystem(binary_config)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
            binary_id = _db.insert_star_system(conn, binary_system, binary_config)

        conn.execute(
            "ALTER TABLE stars "
            "DROP COLUMN galactic_orbital_phase_deg, DROP COLUMN galactic_min_update_interval_years, "
            "DROP COLUMN wide_binary_a_crit_km"
        )
        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_galactic_orbital_phase_deg, DROP COLUMN binary_galactic_min_update_interval_years, "
            "DROP COLUMN binary_mutual_orbital_period_years, DROP COLUMN binary_mutual_orbital_speed_kms, "
            "DROP COLUMN binary_mutual_orbital_inclination_deg, DROP COLUMN binary_mutual_orbital_ascending_node_deg, "
            "DROP COLUMN binary_mutual_orbital_phase_deg, DROP COLUMN binary_mutual_min_update_interval_years, "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (13, 14, 15, 16, 17, 18)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (12)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 12
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION == 18

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        star_rows = conn.execute(
            "SELECT galactic_orbital_phase_deg, galactic_orbital_period_gy, galactic_min_update_interval_years "
            "FROM stars WHERE star_system_id = ?", (system_id,)
        ).fetchall()
        assert star_rows, "test fixture must actually contain a star"
        for row in star_rows:
            assert row["galactic_orbital_phase_deg"] == 0
            expected = minimum_update_interval_years(row["galactic_orbital_period_gy"] * 1e9)
            assert row["galactic_min_update_interval_years"] == pytest.approx(expected, rel=1e-9)

        binary_row = conn.execute(
            "SELECT binary_galactic_orbital_phase_deg, binary_galactic_orbital_period_gy, "
            "binary_galactic_min_update_interval_years, binary_separation_km, binary_effective_mass_kg, "
            "binary_mutual_orbital_period_years, binary_mutual_orbital_speed_kms, "
            "binary_mutual_orbital_inclination_deg, binary_mutual_orbital_ascending_node_deg, "
            "binary_mutual_orbital_phase_deg, binary_mutual_min_update_interval_years "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()
        assert binary_row["binary_galactic_orbital_phase_deg"] == 0
        expected_galactic_guard = minimum_update_interval_years(
            binary_row["binary_galactic_orbital_period_gy"] * 1e9
        )
        assert binary_row["binary_galactic_min_update_interval_years"] == pytest.approx(
            expected_galactic_guard, rel=1e-9
        )

        expected_mutual_period = calculate_orbital_period_years(
            binary_row["binary_separation_km"] / pc.AU_TO_KM, binary_row["binary_effective_mass_kg"]
        )
        assert binary_row["binary_mutual_orbital_period_years"] == pytest.approx(expected_mutual_period, rel=1e-9)

        expected_mutual_speed = circular_orbital_speed_kms(
            binary_row["binary_separation_km"] / pc.AU_TO_KM, expected_mutual_period
        )
        assert binary_row["binary_mutual_orbital_speed_kms"] == pytest.approx(expected_mutual_speed, rel=1e-9)

        assert binary_row["binary_mutual_orbital_inclination_deg"] == 0
        assert binary_row["binary_mutual_orbital_ascending_node_deg"] == 0
        assert binary_row["binary_mutual_orbital_phase_deg"] == 0

        expected_mutual_guard = minimum_update_interval_years(expected_mutual_period)
        assert binary_row["binary_mutual_min_update_interval_years"] == pytest.approx(
            expected_mutual_guard, rel=1e-9
        )
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == 18


def test_migrate_v13_to_v14_adds_and_backfills_binary_mutual_position(mysql_config):
    """
    `mysql_config` yields a fresh database, which `get_connection` already
    creates at the current schema (v15) -- so this simulates an existing
    v13 database by tearing the v14 *and* v15 additions back out (the new
    `star_systems.binary_mutual_position_x_km`/`_y_km`/`_z_km` columns, the
    v15 wide-binary columns, and the `schema_migrations` v14/v15 rows)
    before calling `migrate_database`, and confirms it not only adds the
    columns back but backfills them with a real derived value (not an
    arbitrary placeholder) from each row's own already-stored
    `binary_separation_km` and v13
    `binary_mutual_orbital_{inclination,ascending_node,phase}_deg` columns
    -- see `_migrate_v13_to_v14`'s docstring.
    """
    binary_config = SystemConfig()
    binary_config.BINARY_SYSTEM = True
    binary_config.WIDE_BINARY = False  # pin to close/P-type -- these assertions are BinaryStarProxy-specific
    binary_system = StarSystem(binary_config)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            binary_id = _db.insert_star_system(conn, binary_system, binary_config)

        conn.execute(
            "ALTER TABLE star_systems "
            "DROP COLUMN binary_mutual_position_x_km, DROP COLUMN binary_mutual_position_y_km, "
            "DROP COLUMN binary_mutual_position_z_km, "
            "DROP COLUMN binary_configuration, DROP COLUMN binary_eccentricity, "
            "DROP COLUMN binary_periapsis_km, DROP COLUMN binary_apoapsis_km"
        )
        conn.execute("ALTER TABLE stars DROP COLUMN wide_binary_a_crit_km")
        conn.execute(
            "ALTER TABLE asteroid_belts "
            "DROP FOREIGN KEY fk_asteroid_belts_star, DROP INDEX idx_asteroid_belts_star_id, "
            "DROP COLUMN star_id"
        )
        _drop_v17_phenomenon_columns(conn)
        _drop_v18_phenomenon_columns(conn)
        conn.execute("DELETE FROM schema_migrations WHERE version IN (14, 15, 16, 17, 18)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (13)")
        conn.commit()

        version_before = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()["version"]
        assert version_before == 13
    finally:
        conn.close()

    version_after = _db.migrate_database(mysql_config)
    assert version_after == _db.SCHEMA_VERSION == 18

    conn = _db.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute(
            "SELECT binary_separation_km, binary_mutual_orbital_inclination_deg, "
            "binary_mutual_orbital_ascending_node_deg, binary_mutual_orbital_phase_deg, "
            "binary_mutual_position_x_km, binary_mutual_position_y_km, binary_mutual_position_z_km "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()
        expected_x_au, expected_y_au, expected_z_au = orbital_position_au(
            row["binary_separation_km"] / pc.AU_TO_KM, row["binary_mutual_orbital_inclination_deg"],
            row["binary_mutual_orbital_ascending_node_deg"], row["binary_mutual_orbital_phase_deg"],
        )
        assert row["binary_mutual_position_x_km"] == pytest.approx(expected_x_au * pc.AU_TO_KM, abs=1e-3)
        assert row["binary_mutual_position_y_km"] == pytest.approx(expected_y_au * pc.AU_TO_KM, abs=1e-3)
        assert row["binary_mutual_position_z_km"] == pytest.approx(expected_z_au * pc.AU_TO_KM, abs=1e-3)
    finally:
        conn.close()

    # Idempotent: running it again against an already-current database is a no-op.
    assert _db.migrate_database(mysql_config) == 18


def test_insert_system_config_round_trips_slots_child_rows(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.SLOTS = [
        None,
        {"type": "planet", "planet_class": "M", "moons": 2},
        {"type": "asteroid_belt", "planet_class": None, "moons": None},
    ]

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            config_id = _db.insert_system_config(conn, cfg)
        reloaded = _db.load_system_config(conn, config_id)
    finally:
        conn.close()

    assert reloaded.SLOTS == cfg.SLOTS
    assert reloaded.STAR_TYPE == cfg.STAR_TYPE
