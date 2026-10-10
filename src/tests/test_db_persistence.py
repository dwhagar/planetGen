# tests/test_db_persistence.py

"""
End-to-end save-path tests for `planetgen.db.store.insert_star_system`,
specifically that a freshly generated system's planets and moons land in
their respective schema-v2 tables (`planets` vs. `moons`, see
`schema.sql`'s "v2" header note) rather than sharing one table the way
schema v1 did. There's otherwise no test coverage of `store.py`'s save path
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
import time

import pytest

from planetgen.db import store
from planetgen import tuning
from planetgen.physics import constants as pc
from planetgen.generation.comet import Comet
from planetgen.generation.phenomena.compact_remnant import BlackHole
from planetgen.generation.config import SystemConfig
from planetgen.generation.binary import BinaryStarProxy
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from planetgen.physics.orbits import orbital_position_au












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


def _make_system_with_comets():
    """Retries generation (bounded) until a system with at least one
    elliptical AND at least one parabolic comet comes out -- COMETS=True
    biases generation heavily toward having some, but not toward any
    particular orbit_type mix (see tuning.COMET_PARABOLIC_CHANCE)."""
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.COMETS = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        orbit_types = {c.orbit_type for c in system.comets}
        if {"elliptical", "parabolic"} <= orbit_types:
            return system, cfg
    pytest.fail("could not generate a system with both an elliptical and a parabolic comet")


def test_insert_star_system_persists_comets_in_their_own_table(mysql_config):
    system, cfg = _make_system_with_comets()

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)

        db_comets = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert len(db_comets) == system.comet_count

        star_row = conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'single'", (system_id,)
        ).fetchone()

        db_comets_by_name = {row["name"]: row for row in db_comets}
        for comet in system.comets:
            row = db_comets_by_name[comet.name]
            assert row["orbit_type"] == comet.orbit_type
            # A single star DOES get a real stars.id (unlike a 'close'
            # binary's merged proxy, which has none) -- see insert_comet's
            # docstring and insert_star_system's own primary_star_id handling.
            assert row["star_id"] == star_row["id"]
            assert math.isclose(row["perihelion_distance_km"], comet.perihelion_distance_au * pc.AU_TO_KM, rel_tol=1e-9)
            assert math.isclose(row["eccentricity"], comet.eccentricity, rel_tol=1e-9)
            assert row["composition_summary"] == comet.get_composition_summary()

            comp_rows = conn.execute(
                "SELECT component FROM comet_composition WHERE comet_id = ? ORDER BY position", (row["id"],)
            ).fetchall()
            assert [r["component"] for r in comp_rows] == list(comet.composition)

            if comet.orbit_type == "elliptical":
                assert row["orbital_period_years"] is not None
                assert row["mean_anomaly_deg"] is not None
                assert row["parabolic_mean_anomaly"] is None
            else:
                assert row["orbital_period_years"] is None
                assert row["mean_anomaly_deg"] is None
                assert row["parabolic_mean_anomaly"] is not None
    finally:
        conn.close()


def test_load_star_system_round_trips_comets_with_composition(mysql_config):
    system, cfg = _make_system_with_comets()

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        reloaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert reloaded.comet_count == system.comet_count
    original_by_name = {c.name: c for c in system.comets}
    for comet in reloaded.comets:
        original = original_by_name[comet.name]
        assert comet.orbit_type == original.orbit_type
        assert comet.period_class == original.period_class
        assert math.isclose(comet.perihelion_distance_au, original.perihelion_distance_au, rel_tol=1e-9)
        assert math.isclose(comet.eccentricity, original.eccentricity, rel_tol=1e-9)
        assert math.isclose(comet.distance_au, original.distance_au, rel_tol=1e-9)
        assert math.isclose(comet.orbital_speed_kms, original.orbital_speed_kms, rel_tol=1e-9)
        assert comet.composition == original.composition
        assert comet.system_config is reloaded.system_config
    assert str(reloaded) == str(system)


def _make_wide_binary_system_with_comets():
    """Retries generation (bounded) until a 'wide' (S-type) binary comes
    out with at least one comet on EACH star -- needed to exercise
    insert_star_system's/load_star_system's star_id-based grouping for
    comets (see insert_comet's and load_star_system's own docstrings),
    which _make_system_with_comets' single-star fixture above can't."""
    for _ in range(50):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.COMETS = True
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = True
        system = StarSystem(system_config=cfg)
        if system.binary_type == "wide" and system.comets and system.secondary_comets:
            return system, cfg
    pytest.fail("could not generate a wide binary with comets on both stars")


def test_insert_and_load_star_system_round_trips_comets_for_a_wide_binary(mysql_config):
    """
    A 'wide' binary's two stars each get their own, independently-rolled
    comet population (StarSystem._generate_comets), disambiguated in the
    database by each row's own star_id (see insert_comet's docstring) --
    unlike the single-star case test_insert_star_system_persists_comets_in_their_own_table
    and test_load_star_system_round_trips_comets_with_composition cover
    above, this is the only place that star_id-based split is actually
    exercised through a real insert/load round trip, rather than just the
    in-memory to_dict/from_dict path StarSystem's own serialization tests
    (test_serialization.py) cover.
    """
    system, cfg = _make_wide_binary_system_with_comets()

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)

        primary_row = conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'primary'", (system_id,)
        ).fetchone()
        secondary_row = conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'secondary'", (system_id,)
        ).fetchone()

        db_comets = conn.execute("SELECT name, star_id FROM comets WHERE star_system_id = ?", (system_id,)).fetchall()
        primary_names = {row["name"] for row in db_comets if row["star_id"] == primary_row["id"]}
        secondary_names = {row["name"] for row in db_comets if row["star_id"] == secondary_row["id"]}
        assert primary_names == {c.name for c in system.comets}
        assert secondary_names == {c.name for c in system.secondary_comets}
        # Every row's star_id must land in exactly one of the two buckets
        # above -- confirms there's no third value (e.g. a stray NULL,
        # only ever correct for a 'close' binary's merged proxy) hiding.
        assert primary_names | secondary_names == {row["name"] for row in db_comets}

        reloaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert {c.name for c in reloaded.comets} == {c.name for c in system.comets}
    assert {c.name for c in reloaded.secondary_comets} == {c.name for c in system.secondary_comets}
    assert reloaded.comet_count == system.comet_count


def _make_close_binary_system_with_comets():
    """Retries generation (bounded) until a 'close' (P-type) binary comes
    out with at least one comet -- a close pair's comets are bound to the
    merged BinaryStarProxy (StarSystem._generate_comets), not either
    individually-stored star row, so insert_comet stores star_id as NULL
    for them (see that function's own docstring); untested by any other
    comet fixture above, all of which use a real star_id."""
    for _ in range(50):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.COMETS = True
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = False
        system = StarSystem(system_config=cfg)
        if system.binary_type == "close" and system.comets:
            return system, cfg
    pytest.fail("could not generate a close binary with comets")


def test_insert_and_load_star_system_round_trips_comets_for_a_close_binary(mysql_config):
    system, cfg = _make_close_binary_system_with_comets()

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)

        db_comets = conn.execute("SELECT name, star_id FROM comets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert len(db_comets) == system.comet_count
        assert all(row["star_id"] is None for row in db_comets)

        reloaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()

    assert {c.name for c in reloaded.comets} == {c.name for c in system.comets}
    assert reloaded.secondary_comets == []
    assert reloaded.comet_count == system.comet_count


def test_insert_star_system_splits_planets_and_moons_into_their_own_tables(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
    finally:
        conn.close()

    expected_planets = [obj for obj in system.planets if obj.body_type != "a"]
    expected_belts = [obj for obj in system.planets if obj.body_type == "a"]
    expected_moon_names = sorted(moon.name for planet in expected_planets for moon in planet.moons)
    assert expected_moon_names, "test fixture must actually contain moons"

    conn = store.get_connection(mysql_config)
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
        assert version_row["version"] == store.SCHEMA_VERSION
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Read path (store.load_star_system / load_sector / load_system_config).
# Covers the gaps this module's own docstring called out as still missing:
# binary systems, lifespan_gy round-tripping, and foreign-key integrity,
# plus save_sector/insert_sector and insert_system_config's SLOTS child rows.
# ---------------------------------------------------------------------------

def test_load_star_system_round_trips_single_star_system_exactly(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        reloaded = store.load_star_system(conn, system_id)
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

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        reloaded = store.load_star_system(conn, system_id)
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

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        row = conn.execute(
            "SELECT lifespan_gy FROM stars WHERE star_system_id = ?", (system_id,)
        ).fetchone()
        assert row["lifespan_gy"] is None  # NULL in storage, never the JSON "Infinity" token

        reloaded = store.load_star_system(conn, system_id)
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

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        reloaded = store.load_star_system(conn, system_id)
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

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            binary_system_id = store.insert_star_system(conn, binary_system, binary_cfg)
        reloaded_binary = store.load_star_system(conn, binary_system_id)
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

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        reloaded = store.load_star_system(conn, system_id)
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

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            binary_system_id = store.insert_star_system(conn, binary_system, binary_cfg)
        reloaded_binary = store.load_star_system(conn, binary_system_id)
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

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.insert_star_system(conn, system, cfg)
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

    sector_id = store.save_sector(sector, config=mysql_config)

    conn = store.get_connection(mysql_config)
    try:
        reloaded_sector = store.load_sector(conn, sector_id)
    finally:
        conn.close()

    assert reloaded_sector.name == "Persisted Sector"
    assert reloaded_sector.edge_ly == pytest.approx(20.0)
    assert len(reloaded_sector) == 2

    reloaded_by_name = {entry.star_system.name: entry for entry in reloaded_sector.entries}
    for original_system, original_position in ((system_a, (1.0, 2.0, 3.0)), (system_b, (-4.0, 0.0, 5.5))):
        entry = reloaded_by_name[original_system.name]
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

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        reloaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()

    # Matched positionally, not by name: `insert_moon`/`insert_planet` persist
    # each body's list position via `orbital_index`, and `load_star_system`
    # loads both planets and moons `ORDER BY orbital_index`, so `reloaded`
    # preserves the original ordering exactly. Name-keyed dicts would break
    # here -- moon (and planet) names come from `generate_phoneme_salad_name`,
    # which picks names independently per body with no check against its own
    # siblings, so two moons on the same planet can occasionally collide on
    # the same generated name and silently overwrite each other in a dict.
    reloaded_planets = [p for p in reloaded.planets if p.body_type != "a"]
    assert len(reloaded_planets) == len(planets)
    for planet, reloaded_planet in zip(planets, reloaded_planets):
        assert reloaded_planet.orbital_inclination_deg == pytest.approx(planet.orbital_inclination_deg)
        assert reloaded_planet.orbital_ascending_node_deg == pytest.approx(planet.orbital_ascending_node_deg)
        assert reloaded_planet.orbital_phase_deg == pytest.approx(planet.orbital_phase_deg)
        assert reloaded_planet.position_x == pytest.approx(planet.position_x)
        assert reloaded_planet.position_y == pytest.approx(planet.position_y)
        assert reloaded_planet.position_z == pytest.approx(planet.position_z)
        assert reloaded_planet.orbital_speed_kms == pytest.approx(planet.orbital_speed_kms)
        assert reloaded_planet.min_update_interval_years == pytest.approx(planet.min_update_interval_years)
        assert reloaded_planet.rotation_period_hours == pytest.approx(planet.rotation_period_hours)

        assert len(reloaded_planet.moons) == len(planet.moons)
        for moon, reloaded_moon in zip(planet.moons, reloaded_planet.moons):
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
# Orbital motion updates (planetgen.db.store.advance_orbital_phases /
# get_orbit_update_elapsed_years) -- see planetgen.cli.orbits.
# ---------------------------------------------------------------------------

def test_get_orbit_update_elapsed_years_is_none_before_first_update(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        assert store.get_orbit_update_elapsed_years(conn) is None
    finally:
        conn.close()


def test_advance_orbital_phases_applies_the_correct_delta_and_leaves_other_fields_untouched(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.insert_star_system(conn, system, cfg)

        before = conn.execute(
            "SELECT id, orbital_phase_deg, orbital_inclination_deg, orbital_ascending_node_deg, "
            "rotation_period_hours, period_years, distance_km, "
            "position_x_km, position_y_km, position_z_km, orbital_speed_kms FROM planets"
        ).fetchall()

        clock = store.orbit_clock(conn, 0.5)
        counts = store.advance_orbital_phases(conn, clock)
        assert counts["planets"] > 0
        assert counts["moons"] > 0
        store.finish_orbit_update(conn, clock)

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

    # The clock now stands at the end of the step (half a year from now,
    # since no update had run), not None.
    conn = store.get_connection(mysql_config)
    try:
        epoch = store.get_orbit_epoch_unix(conn)
    finally:
        conn.close()
    assert epoch == pytest.approx(clock.end_unix, abs=1.0)


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

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            binary_id = store.insert_star_system(conn, binary_system, binary_cfg)

        before = conn.execute(
            "SELECT binary_separation_km, binary_mutual_orbital_period_years, "
            "binary_mutual_orbital_inclination_deg, binary_mutual_orbital_ascending_node_deg, "
            "binary_mutual_orbital_phase_deg "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()

        # Big enough elapsed time to clearly move the (fast) mutual orbit,
        # regardless of its randomly generated period.
        elapsed_years = before["binary_mutual_orbital_period_years"] * 137.25
        store.advance_orbital_phases(conn, store.orbit_clock(conn, elapsed_years))

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


def test_advance_orbital_phases_advances_both_binary_members_barycenter_offsets(mysql_config):
    """
    Regression test for the v18 addition: both `binary_primary_position_*`
    and `binary_secondary_position_*` must move together with the freshly-
    advanced `binary_mutual_position_*` -- not just one star sitting fixed
    while the other orbits it -- and `secondary_position == primary_position
    + mutual_position` must keep holding after every advance, not just at
    generation time.
    """
    binary_cfg = SystemConfig()
    binary_cfg.STAR_TYPE = "G2V"
    binary_cfg.BINARY_SYSTEM = True
    binary_cfg.WIDE_BINARY = False
    binary_cfg.PLANETS = False
    binary_system = StarSystem(system_config=binary_cfg)

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            binary_id = store.insert_star_system(conn, binary_system, binary_cfg)

        before = conn.execute(
            "SELECT binary_mutual_orbital_period_years, binary_primary_position_x_km, "
            "binary_secondary_position_x_km, binary_secondary_mass_fraction "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()

        elapsed_years = before["binary_mutual_orbital_period_years"] * 137.25
        counts = store.advance_orbital_phases(conn, store.orbit_clock(conn, elapsed_years))
        assert counts["binary_mutual_orbits"] > 0

        after = conn.execute(
            "SELECT binary_mutual_position_x_km, binary_mutual_position_y_km, binary_mutual_position_z_km, "
            "binary_primary_position_x_km, binary_primary_position_y_km, binary_primary_position_z_km, "
            "binary_secondary_position_x_km, binary_secondary_position_y_km, binary_secondary_position_z_km "
            "FROM star_systems WHERE id = ?", (binary_id,)
        ).fetchone()
    finally:
        conn.close()

    # Both members actually moved -- neither is the trivial "stayed put" case.
    assert after["binary_primary_position_x_km"] != pytest.approx(before["binary_primary_position_x_km"])
    assert after["binary_secondary_position_x_km"] != pytest.approx(before["binary_secondary_position_x_km"])

    mu = before["binary_secondary_mass_fraction"]
    assert after["binary_primary_position_x_km"] == pytest.approx(-mu * after["binary_mutual_position_x_km"], abs=1e-3)
    assert after["binary_primary_position_y_km"] == pytest.approx(-mu * after["binary_mutual_position_y_km"], abs=1e-3)
    assert after["binary_primary_position_z_km"] == pytest.approx(-mu * after["binary_mutual_position_z_km"], abs=1e-3)
    assert after["binary_secondary_position_x_km"] == pytest.approx(
        after["binary_primary_position_x_km"] + after["binary_mutual_position_x_km"], abs=1e-3
    )
    assert after["binary_secondary_position_y_km"] == pytest.approx(
        after["binary_primary_position_y_km"] + after["binary_mutual_position_y_km"], abs=1e-3
    )
    assert after["binary_secondary_position_z_km"] == pytest.approx(
        after["binary_primary_position_z_km"] + after["binary_mutual_position_z_km"], abs=1e-3
    )


def test_advance_orbital_phases_recomputes_star_and_planet_reflex_offsets(mysql_config):
    """
    Regression test for the v18 addition: `stars.reflex_offset_*_km`/
    `planets.reflex_offset_*_km` must be recomputed from each hosted
    body's freshly-advanced position on every call (they have no
    independent guard interval of their own -- see
    `advance_orbital_phases`' own docstring), not just set once at
    generation time.
    """
    system = None
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        candidate = StarSystem(system_config=cfg)
        real_planets = [p for p in candidate.planets if p.body_type != "a"]
        if real_planets and any(p.moons for p in real_planets):
            system = candidate
            break
    if system is None:
        pytest.fail("could not generate a single-star system with planets and at least one moon")

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, system.system_config)

        star_id = conn.execute(
            "SELECT id FROM stars WHERE star_system_id = ? AND role = 'single'", (system_id,)
        ).fetchone()["id"]
        before_offset = conn.execute(
            "SELECT reflex_offset_x_km FROM stars WHERE id = ?", (star_id,)
        ).fetchone()["reflex_offset_x_km"]

        # Big enough elapsed time to clearly move every planet, regardless
        # of its randomly generated period.
        max_period = conn.execute(
            "SELECT MAX(period_years) AS p FROM planets WHERE star_id = ?", (star_id,)
        ).fetchone()["p"]
        counts = store.advance_orbital_phases(conn, store.orbit_clock(conn, max_period * 137.25))
        assert counts["star_reflex_offsets"] > 0
        assert counts["planet_reflex_offsets"] > 0

        star_row = conn.execute(
            "SELECT mass_kg, reflex_offset_x_km, reflex_offset_y_km, reflex_offset_z_km "
            "FROM stars WHERE id = ?", (star_id,)
        ).fetchone()
        planet_rows = conn.execute(
            "SELECT mass_kg, position_x_km, position_y_km, position_z_km FROM planets WHERE star_id = ?",
            (star_id,),
        ).fetchall()
    finally:
        conn.close()

    # The offset actually changed -- not accidentally still the
    # generation-time value.
    assert star_row["reflex_offset_x_km"] != pytest.approx(before_offset)

    expected_x = -sum(
        (p["mass_kg"] / (star_row["mass_kg"] + p["mass_kg"])) * p["position_x_km"] for p in planet_rows
    )
    assert star_row["reflex_offset_x_km"] == pytest.approx(expected_x, abs=1e-3)


def test_advance_orbital_phases_rejects_negative_elapsed_years(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError):
            store.advance_orbital_phases(conn, store.orbit_clock(conn, -1.0))
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Comet orbital motion updates (planetgen.db.store.advance_comet_orbits) --
# a separate Python-loop call from advance_orbital_phases above, since a
# comet's position isn't a linear function of elapsed time the way a
# circular planet/moon orbit's is -- see that function's own docstring.
# ---------------------------------------------------------------------------

def test_advance_comet_orbits_advances_elliptical_mean_anomaly_and_recomputes_position(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.COMETS = False
    system = StarSystem(system_config=cfg)
    comet = Comet(cfg, primary_mass_solar=system.star.mass / pc.SOLAR_MASS_TO_KG, orbit_type="elliptical")
    comet.mean_anomaly_deg = 10.0
    comet.update_orbital_state()
    system.comets = [comet]

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)

        before = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
        step_years = before["orbital_period_years"] * 0.05  # 18 degrees of mean anomaly

        updated = store.advance_comet_orbits(conn, store.orbit_clock(conn, step_years))
        assert updated == 1

        after = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
    finally:
        conn.close()

    assert math.isclose(after["mean_anomaly_deg"], 28.0, rel_tol=1e-9)  # 10 + 18
    assert after["parabolic_mean_anomaly"] is None
    # Moved further from perihelion (mean anomaly 10 -> 28, still well
    # short of aphelion) -- distance must have increased.
    assert after["distance_km"] > before["distance_km"]
    # Position is a pure function of distance/orientation/anomaly, so it
    # must have moved along with mean_anomaly_deg/distance_km.
    assert (after["position_x_km"], after["position_y_km"], after["position_z_km"]) != (
        before["position_x_km"], before["position_y_km"], before["position_z_km"]
    )


def test_advance_comet_orbits_advances_parabolic_mean_anomaly_without_wrapping(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.COMETS = False
    system = StarSystem(system_config=cfg)
    comet = Comet(cfg, primary_mass_solar=system.star.mass / pc.SOLAR_MASS_TO_KG, orbit_type="parabolic")
    comet.parabolic_mean_anomaly = -1.0  # still approaching perihelion
    comet.update_orbital_state()
    system.comets = [comet]

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)

        before = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
        assert before["mean_anomaly_deg"] is None
        assert before["min_update_interval_years"] is None

        # A large elapsed time -- since a parabolic anomaly doesn't wrap
        # (unlike mean_anomaly_deg's MOD 360), this should NOT be clamped
        # or wrapped, just added linearly.
        updated = store.advance_comet_orbits(conn, store.orbit_clock(conn, 50.0))
        assert updated == 1

        after = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
    finally:
        conn.close()

    assert after["mean_anomaly_deg"] is None
    assert after["parabolic_mean_anomaly"] > before["parabolic_mean_anomaly"]
    assert after["distance_km"] > before["distance_km"]


def test_advance_comet_orbits_skips_elliptical_rows_below_min_update_interval(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.COMETS = False
    system = StarSystem(system_config=cfg)
    comet = Comet(cfg, primary_mass_solar=system.star.mass / pc.SOLAR_MASS_TO_KG, orbit_type="elliptical")
    system.comets = [comet]

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)

        before = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
        # Comfortably below this comet's own min_update_interval_years floor.
        tiny_elapsed = before["min_update_interval_years"] / 2
        updated = store.advance_comet_orbits(conn, store.orbit_clock(conn, tiny_elapsed))
        assert updated == 0

        after = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchone()
    finally:
        conn.close()

    assert after["mean_anomaly_deg"] == before["mean_anomaly_deg"]
    assert after["distance_km"] == before["distance_km"]


def test_advance_comet_orbits_rejects_negative_elapsed_years(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError):
            store.advance_comet_orbits(conn, store.orbit_clock(conn, -1.0))
    finally:
        conn.close()


















def test_insert_system_config_round_trips_slots_child_rows(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.SLOTS = [
        None,
        {"type": "planet", "planet_class": "M", "moons": 2},
        {"type": "asteroid_belt", "planet_class": None, "moons": None},
    ]

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            config_id = store.insert_system_config(conn, cfg)
        reloaded = store.load_system_config(conn, config_id)
    finally:
        conn.close()

    assert reloaded.SLOTS == cfg.SLOTS
    assert reloaded.STAR_TYPE == cfg.STAR_TYPE












def test_modified_at_tracks_edits_but_not_orbit_ticks(mysql_config):
    """
    v27: `modified_at` moves when a row is edited (MySQL's `ON UPDATE`)
    or one of a system's child rows changes (`touch_star_system`), but
    NOT when the orbit update ticks the simulation clock forward
    -- see `schema.sql`'s "v27" header note.
    """
    binary_cfg = SystemConfig()
    binary_cfg.STAR_TYPE = "G2V"
    binary_cfg.BINARY_SYSTEM = True
    binary_cfg.WIDE_BINARY = False
    binary_cfg.PLANETS = False
    binary_system = StarSystem(system_config=binary_cfg)

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, binary_system, binary_cfg)
            bh_id = store.insert_black_hole(conn, BlackHole(SystemConfig()))

        def modified(table, row_id):
            return conn.execute(f"SELECT modified_at FROM {table} WHERE id = ?", (row_id,)).fetchone()["modified_at"]

        system_before = modified("star_systems", system_id)
        bh_before = modified("black_holes", bh_id)
        time.sleep(0.05)

        clock = store.orbit_clock(conn, 1e5)
        counts = store.advance_orbital_phases(conn, clock)
        assert counts["binary_mutual_orbits"] > 0
        motion = store.advance_galactic_positions(conn, clock)
        assert motion["counts"]["star_systems"] > 0
        assert motion["counts"]["black_holes"] > 0
        assert modified("star_systems", system_id) == system_before
        assert modified("black_holes", bh_id) == bh_before

        with conn:
            conn.execute("UPDATE star_systems SET name = ? WHERE id = ?", ("Renamed For Test", system_id))
        system_renamed = modified("star_systems", system_id)
        assert system_renamed > system_before

        time.sleep(0.05)
        with conn:
            store.touch_star_system(conn, system_id)
        assert modified("star_systems", system_id) > system_renamed
    finally:
        conn.close()




def _sector_with_one_system(name):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    sector = SpaceSector(name, edge_ly=11.5)
    sector.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    return sector














def test_runaway_flags_round_trip(mysql_config):
    system = StarSystem(SystemConfig())
    system.runaway_class, system.runaway_speed_kms = "runaway", 55.0
    system_id = store.save_system(system, SystemConfig(), config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        loaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()
    assert (loaded.runaway_class, loaded.runaway_speed_kms) == ("runaway", 55.0)
    assert StarSystem.from_dict(loaded.to_dict()).runaway_class == "runaway"




def test_classed_phenomena_round_trip_their_contents(mysql_config):
    from planetgen.generation.phenomena.asteroid_field import AsteroidField
    from planetgen.generation.phenomena.nebula import Nebula

    nebula = Nebula(SystemConfig(), nebula_class="Q")
    field = AsteroidField(SystemConfig())
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            nebula_id = store.insert_nebula(conn, nebula)
            field_id = store.insert_asteroid_field(conn, field)
        row = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        assert (row["nebula_class"], row["nebula_type"]) == ("Q", "dark")
        assert row["temperature_k"] == pytest.approx(nebula.temperature_k)
        field_row = conn.execute("SELECT * FROM asteroid_fields WHERE id = ?", (field_id,)).fetchone()
        assert field_row["field_class"] == field.field_class
    finally:
        conn.close()


def test_innermost_container_picks_the_smallest_cloud_holding_the_point():
    big = {"column": "inside_nebula_id", "id": 1, "center": (0.0, 0.0, 0.0), "radius_pc": 10.0}
    small = {"column": "inside_remnant_id", "id": 2, "center": (1.0, 0.0, 0.0), "radius_pc": 2.0}
    assert store.innermost_container((1.5, 0.0, 0.0), [big, small]) is small
    assert store.innermost_container((5.0, 0.0, 0.0), [big, small]) is big
    assert store.innermost_container((20.0, 0.0, 0.0), [big, small]) is None
    # A nebula only nests in a larger cloud, and never in itself.
    assert store.innermost_container((1.0, 0.0, 0.0), [big, small], own_radius_pc=3.0) is big
    assert store.innermost_container((0.0, 0.0, 0.0), [big], own_radius_pc=0.0,
                                   own=("inside_nebula_id", 1)) is None


def test_a_nebula_holds_only_the_points_inside_its_shape():
    """GEN.75: a point in the nebula's sphere but outside its shape is not inside it."""
    class Slab:
        def contains(self, point):
            return abs(point[2]) < 0.1

    nebula = {"column": "inside_nebula_id", "id": 1, "center": (0.0, 0.0, 0.0), "radius_pc": 10.0,
              "shape": lambda: Slab()}
    assert store.innermost_container((5.0, 0.0, 0.5), [nebula]) is nebula
    assert store.innermost_container((5.0, 0.0, 5.0), [nebula]) is None
    assert store.innermost_container((11.0, 0.0, 0.0), [nebula]) is None


def _point_inside(shape):
    """A point inside `shape`: a metaball centre when one is (the usual
    case), else the nearest to the middle of a coarse grid. The shape's warp
    can carry every centre out of it (about 1 draw in 40), which made this
    helper's first version raise StopIteration at random."""
    centres = ([v / shape.scale for v in c] for c in shape.centres)
    found = next((point for point in centres if shape.contains(point)), None)
    if found is not None:
        return found
    steps = [i / 4.0 for i in range(-8, 9)]
    grid = sorted(([x, y, z] for x in steps for y in steps for z in steps), key=lambda p: sum(v * v for v in p))
    return next(point for point in grid if shape.contains(point))


def _placed_nebula(conn, sector_id, around_pc, radius_ly, nebula_class="D"):
    """A nebula placed so `around_pc` lies inside its shape (GEN.75): a
    nebula holds only the points within its shape, not its whole sphere."""
    from planetgen.generation.phenomena.nebula import Nebula
    nebula = Nebula(SystemConfig(), nebula_class=nebula_class)
    nebula.radius_ly = radius_ly
    shape = nebula.get_shape()
    from planetgen.physics.units import ly_to_pc

    inside = _point_inside(shape)
    x, y, z = (around_pc[i] - inside[i] * ly_to_pc(radius_ly) for i in range(3))
    return store.insert_nebula(conn, nebula, sector_id=sector_id, placement={
        "center_x_pc": x, "center_y_pc": y, "center_z_pc": z, "galactic_radius_pc": math.hypot(x, y, z),
    })


def _system_point_pc(conn, sector_id):
    """The first system of `sector_id` in galaxy-frame parsecs."""
    sector = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc FROM sectors WHERE id = ?",
                          (sector_id,)).fetchone()
    system = conn.execute("SELECT position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems"
                          " WHERE sector_id = ?", (sector_id,)).fetchone()
    return store.local_to_galaxy_pc(
        (sector["center_x_pc"], sector["center_y_pc"], sector["center_z_pc"]),
        (system["position_x_mpc"] / 1000.0, system["position_y_mpc"] / 1000.0, system["position_z_mpc"] / 1000.0))


def test_systems_inside_a_nebula_point_at_it_and_the_innermost_wins(mysql_config):
    sector_id = store.save_sector(_sector_with_one_system("Cloudy"), config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0,
    })
    far_id = store.save_sector(_sector_with_one_system("Clear"), config=mysql_config, galaxy_position={
        "center_x_pc": 400.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 400.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            big = _placed_nebula(conn, sector_id, _system_point_pc(conn, sector_id), radius_ly=100.0, nebula_class="E")
        row = conn.execute("SELECT inside_nebula_id, inside_remnant_id FROM star_systems WHERE sector_id = ?",
                           (sector_id,)).fetchone()
        assert (row["inside_nebula_id"], row["inside_remnant_id"]) == (big, None)
        far = conn.execute("SELECT inside_nebula_id FROM star_systems WHERE sector_id = ?", (far_id,)).fetchone()
        assert far["inside_nebula_id"] is None

        with conn:
            small = _placed_nebula(conn, sector_id, _system_point_pc(conn, sector_id), radius_ly=3.0, nebula_class="C")
        row = conn.execute("SELECT inside_nebula_id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()
        assert row["inside_nebula_id"] == small
        nested = conn.execute("SELECT inside_nebula_id FROM nebulae WHERE id = ?", (small,)).fetchone()
        assert nested["inside_nebula_id"] == big

        with conn:
            conn.execute("DELETE FROM nebulae WHERE id = ?", (small,))
            store.refresh_containment(conn, [sector_id])
        row = conn.execute("SELECT inside_nebula_id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()
        assert row["inside_nebula_id"] == big
        from planetgen.db import query as queryDb
        inside = queryDb.sector_detail(conn, sector_id)["systems"][0]["inside"]
        assert (inside["type"], inside["id"], inside["class"]) == ("nebula", big, "E")
        assert inside["descriptor"] in ("diffuse", "emission", "reflection", "planetary", "dark")
    finally:
        conn.close()


def test_containment_counts_what_entered_and_left_a_nebula(mysql_config):
    """GEN.107: the orbit update's summary counts nebula entries and exits
    from `refresh_containment`."""
    sector_id = store.save_sector(_sector_with_one_system("Drifting"), config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            cloud = _placed_nebula(conn, sector_id, _system_point_pc(conn, sector_id), radius_ly=3.0, nebula_class="C")
        home = conn.execute("SELECT center_x_pc FROM nebulae WHERE id = ?", (cloud,)).fetchone()["center_x_pc"]
        with conn:
            conn.execute("UPDATE nebulae SET center_x_pc = center_x_pc + 50 WHERE id = ?", (cloud,))
            left = store.refresh_containment(conn, [sector_id])
        assert left == {"entered_nebula": 0, "left_nebula": 1, "entered_remnant": 0, "left_remnant": 0}
        with conn:
            conn.execute("UPDATE nebulae SET center_x_pc = ? WHERE id = ?", (home, cloud))
            entered = store.refresh_containment(conn, [sector_id])
        assert entered == {"entered_nebula": 1, "left_nebula": 0, "entered_remnant": 0, "left_remnant": 0}
        with conn:
            assert store.refresh_containment(conn, [sector_id]) == dict.fromkeys(store.CONTAINMENT_COUNTS, 0)
    finally:
        conn.close()




def _named_system(name):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    system.name = name
    return system, cfg


def _remnant_with_core(name):
    from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
    for _ in range(200):
        remnant = SupernovaRemnant(SystemConfig(), name=name)
        if remnant.compact_remnant is not None:
            return remnant
    pytest.fail("no supernova remnant with a detectable core")


def test_phenomena_share_the_system_name_registry(mysql_config):
    """GEN.13 (v40): a nebula whose name clashes with a system's is
    decorated like a second system would be, and the first holder is
    renamed whichever kind it is."""
    from planetgen.generation.phenomena.nebula import Nebula

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system, cfg = _named_system("Kelvaro")
            system_id = store.insert_star_system(conn, system, cfg)
            nebula = Nebula(SystemConfig(), name="Kelvaro")
            nebula_id = store.insert_nebula(conn, nebula)
        assert conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"] == "Alpha Kelvaro"
        assert conn.execute("SELECT name FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()["name"] == "Beta Kelvaro"
        assert nebula.name == "Beta Kelvaro"

        with conn:
            first = store.insert_rogue_planet(conn, __import__(
                "planetgen.generation.phenomena.rogue", fromlist=["RoguePlanet"]).RoguePlanet(SystemConfig(), name="Ossandre"))
            system, cfg = _named_system("Ossandre")
            store.insert_star_system(conn, system, cfg)
        assert conn.execute("SELECT name FROM rogue_planets WHERE id = ?", (first,)).fetchone()["name"] == "Alpha Ossandre"
        assert system.name == "Beta Ossandre"
        assert store.name_in_use(conn, "Alpha Ossandre") == "rogue_planets"
    finally:
        conn.close()


def test_a_renamed_remnant_takes_its_core_along(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            remnant_id = store.insert_supernova_remnant(conn, _remnant_with_core("Thessavel"))
            store.insert_nebula(conn, __import__(
                "planetgen.generation.phenomena.nebula", fromlist=["Nebula"]).Nebula(SystemConfig(), name="Thessavel"))
        row = conn.execute("SELECT * FROM supernova_remnants WHERE id = ?", (remnant_id,)).fetchone()
        assert row["name"] == "Alpha Thessavel"
        core_table = "black_holes" if row["compact_remnant_black_hole_id"] else "neutron_stars"
        core_id = row["compact_remnant_black_hole_id"] or row["compact_remnant_neutron_star_id"]
        core = conn.execute(f"SELECT name FROM {core_table} WHERE id = ?", (core_id,)).fetchone()
        assert core["name"] == "Alpha Thessavel Core"
        registry = conn.execute(
            "SELECT occurrence_count, first_object_table FROM system_name_registry WHERE base_name = 'Thessavel'"
        ).fetchone()
        # The core isn't registered on its own.
        assert (registry["occurrence_count"], registry["first_object_table"]) == (2, "supernova_remnants")
    finally:
        conn.close()


def test_comets_carry_designations_that_follow_their_host(mysql_config):
    from planetgen.generation.comet import PERIODIC_COMET_MAX_PERIOD_YEARS

    system, cfg = _make_system_with_comets()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        rows = conn.execute(
            "SELECT name, orbit_type, orbital_period_years FROM comets WHERE star_system_id = ? ORDER BY id",
            (system_id,),
        ).fetchall()
        for n, row in enumerate(rows, start=1):
            periodic = row["orbit_type"] == "elliptical" and row["orbital_period_years"] < PERIODIC_COMET_MAX_PERIOD_YEARS
            assert row["name"] == f"{'P' if periodic else 'C'}/{system.name}-{n}"
        with conn:
            store.rename_star_system(conn, system_id, "Neraloth")
        renamed = conn.execute("SELECT name FROM comets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()
        assert [r["name"][2:] for r in renamed] == [f"Neraloth-{n}" for n in range(1, len(rows) + 1)]
    finally:
        conn.close()


def test_interstellar_comets_and_asteroid_fields_are_designated_by_sector(mysql_config):
    from planetgen.generation.phenomena.asteroid_field import AsteroidField
    from planetgen.generation.phenomena.rogue import InterstellarComet

    sector_id = store.save_sector(_sector_with_one_system("Designated"), config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.insert_interstellar_comet(conn, InterstellarComet(SystemConfig()), sector_id=sector_id)
            store.insert_interstellar_comet(conn, InterstellarComet(SystemConfig()), sector_id=sector_id)
            field = AsteroidField(SystemConfig())
            store.insert_asteroid_field(conn, field, sector_id=sector_id)
        names = [r["name"] for r in conn.execute(
            "SELECT name FROM interstellar_comets WHERE sector_id = ? ORDER BY id", (sector_id,)).fetchall()]
        sector_name = conn.execute("SELECT name FROM sectors WHERE id = ?", (sector_id,)).fetchone()["name"]
        assert names == [f"I/{sector_name}-1", f"I/{sector_name}-2"]
        assert field.name == f"AF {field.field_class}-{sector_name}-01"
    finally:
        conn.close()




def test_system_grid_nearest_matches_brute_force():
    import random as _random
    rng = _random.Random(7)
    systems = [(i, (rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(-5, 5))) for i in range(300)]
    grid = store._SystemGrid(systems)
    for _ in range(50):
        point = (rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(-5, 5))
        brute = sorted((math.dist(point, p), i) for i, p in systems if math.dist(point, p) <= 4.0)[:3]
        assert grid.nearest(point) == brute


def test_galaxy_to_local_undoes_local_to_galaxy():
    from planetgen.galaxy.geometry import galaxy_to_local_pc, local_to_galaxy_pc
    center = (30.0, -40.0, 2.0)
    offset = (1.2, -0.7, 0.3)
    assert galaxy_to_local_pc(center, local_to_galaxy_pc(center, offset)) == pytest.approx(offset)


def _sector_with_systems(name, positions_ly):
    sector = SpaceSector(name, edge_ly=13.0)
    for position in positions_ly:
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = False
        cfg.BINARY_SYSTEM = False
        sector.add_system(StarSystem(system_config=cfg), position=position, system_config=cfg)
    return sector


def test_nearest_systems_cross_sector_boundaries(mysql_config):
    """UX.18 (v41): a system at a sector's edge lists a system just
    across the boundary once that sector is generated."""
    from planetgen.db import query as queryDb
    from planetgen.physics.units import pc_to_ly

    edge_ly = pc_to_ly(4.0)
    first = store.save_sector(_sector_with_systems("Westmark", [(edge_ly / 2 - 0.5, 0.0, 0.0), (-5.0, 0.0, 0.0)]),
                            config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        border_id, far_id = [r["id"] for r in conn.execute(
            "SELECT id FROM star_systems WHERE sector_id = ? ORDER BY position_x_mpc DESC", (first,)).fetchall()]
        before = queryDb.nearest_systems(conn, "star_systems", [border_id])[border_id]
        assert [n["id"] for n in before] == [far_id]
    finally:
        conn.close()

    second = store.save_sector(_sector_with_systems("Eastmark", [(-edge_ly / 2 + 0.5, 0.0, 0.0)]),
                             config=mysql_config, galaxy_position={
        "center_x_pc": 104.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 104.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        across = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (second,)).fetchone()["id"]
        after = queryDb.nearest_systems(conn, "star_systems", [border_id])[border_id]
        assert [n["id"] for n in after] == [across, far_id]
        assert after[0]["distance_ly"] == pytest.approx(1.0, abs=0.01)
        own = queryDb.nearest_systems(conn, "star_systems", [across])[across]
        assert own[0]["id"] == border_id
        detail = queryDb.sector_detail(conn, first)
        by_id = {system["id"]: system for system in detail["systems"]}
        assert by_id[border_id]["nearest"][0]["id"] == across
    finally:
        conn.close()


def test_sectors_linked_later_get_the_links_they_would_have_had_at_save(mysql_config):
    """PERF.45: sectors saved with `link_neighbors=False` and linked at the end
    have the nearest systems a brute-force search gives, the neighbour that was
    already linked included."""
    import math
    from planetgen.physics.units import milliparsecs_to_ly, pc_to_ly

    edge_ly = pc_to_ly(4.0)

    def place(center_x):
        return {"center_x_pc": center_x, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": center_x}

    first = store.save_sector(_sector_with_systems("Linked", [(edge_ly / 2 - 0.5, 0.0, 0.0), (-5.0, 0.0, 0.0)]),
                              config=mysql_config, galaxy_position=place(100.0))
    saved = [
        store.save_sector(_sector_with_systems("LaterB", [(-edge_ly / 2 + 0.5, 0.0, 0.0), (0.0, 3.0, 0.0)]),
                          config=mysql_config, galaxy_position=place(104.0), link_neighbors=False),
        store.save_sector(_sector_with_systems("LaterC", [(-edge_ly / 2 + 1.5, 1.0, 0.0), (2.0, 2.0, 2.0)]),
                          config=mysql_config, galaxy_position=place(108.0), link_neighbors=False),
    ]
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM nearest_systems WHERE sector_id IN (?, ?)",
                            tuple(saved)).fetchone()["n"] == 0
    finally:
        conn.close()
    steps = []
    assert store.link_sector_neighbors(mysql_config, saved, lambda *step: steps.append(step)) == 2
    # UX.83: three steps a sector (containment, nearest systems, the neighbours' lists), ending full.
    assert [label for _done, _total, label in steps[:3]] == ["containment", "nearest systems", "neighbours' lists"]
    assert steps[-1][:2] == (6, 6) and all(total == 6 for _done, total, _label in steps)
    assert [done for done, _total, _label in steps] == sorted(done for done, _total, _label in steps)

    conn = store.get_connection(mysql_config)
    try:
        systems = {}
        for row in conn.execute(
                "SELECT s.id, s.position_x_mpc AS x, s.position_y_mpc AS y, s.position_z_mpc AS z, c.center_x_pc AS cx"
                " FROM star_systems s JOIN sectors c ON c.id = s.sector_id").fetchall():
            systems[row["id"]] = (pc_to_ly(row["cx"]) + milliparsecs_to_ly(row["x"]), milliparsecs_to_ly(row["y"]),
                                  milliparsecs_to_ly(row["z"]))
        limit_ly = pc_to_ly(store.NEAREST_SYSTEMS_SEARCH_PC)
        for system_id, point in systems.items():
            expected = sorted((math.dist(point, other), other_id) for other_id, other in systems.items()
                              if other_id != system_id and math.dist(point, other) <= limit_ly)
            stored = conn.execute(
                "SELECT neighbor_system_id AS n FROM nearest_systems WHERE object_table = 'star_systems'"
                " AND object_id = ? ORDER BY neighbor_rank", (system_id,)).fetchall()
            assert [row["n"] for row in stored] == [other_id for _d, other_id in expected][:store.NEAREST_SYSTEMS_COUNT]
    finally:
        conn.close()


def test_phenomena_store_their_octant_and_nearest_systems(mysql_config):
    from planetgen.db import query as queryDb
    sector_id = store.save_sector(_sector_with_systems("Octmark", [(1.0, 1.0, 1.0)]), config=mysql_config,
                                galaxy_position={
        "center_x_pc": 0.0, "center_y_pc": 200.0, "center_z_pc": 0.0, "galactic_radius_pc": 200.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            # Local +X points away from the galactic axis (+y here), local +Y along -x.
            nebula_id = _placed_nebula(conn, sector_id, (0.5, 200.5, -0.5), radius_ly=0.5)
            store.refresh_nearest_systems(conn, [sector_id])
        quadrant = conn.execute("SELECT quadrant FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()["quadrant"]
        from planetgen.galaxy.sector import classify_octant
        assert quadrant == classify_octant((0.5, -0.5, -0.5))[0]
        listed = [p for p in queryDb.phenomena_near_sector(conn, sector_id) if p["type"] == "nebula"][0]
        assert listed["octant"] == quadrant
        assert len(listed["nearest"]) == 1
    finally:
        conn.close()




def _bright_row(ring=3, layer=0, slot=1, luminosity_sol=800.0):
    return (ring, layer, slot, 12000, 3000, 0, "young", "B2V", "V", 1.4e31, 3.0e6, 22000.0,
            luminosity_sol * 3.828e26,
            0.02, 0.03, 7.0, 0.03, 42)


def test_bright_stars_store_and_clear(mysql_config):
    """Bright-star pre-placement storage (v43): bulk insert, per-sector
    lookup brightest first, fill link, and a plan re-run clearing it."""
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc,"
                         " bulge_scale_radius_pc, bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude,"
                         " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
                         " expected_system_count_at_density_1, outer_ring_index)"
                         " VALUES (1, 1, 1, 1, 1, 2, 0.2, 0.3, 1, 0, 1, 4, 10, 5)")
            assert store.bright_star_scatter_settings(conn) is None
            written = store.insert_bright_stars(conn, [_bright_row(luminosity_sol=600.0), _bright_row(),
                                                     _bright_row(slot=2)], batch_size=2)
            store.record_bright_star_scatter(conn, 500.0, 7)
        assert written == 3
        assert store.bright_star_scatter_settings(conn) == (500.0, 7)
        found = store.bright_stars_for_sector(conn, 3, 0, 1)
        assert len(found) == 2 and found[0]["luminosity_w"] > found[1]["luminosity_w"]

        system, cfg = _named_system("Beaconholm")
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
            store.mark_bright_star_filled(conn, found[0]["id"], system_id)
        assert len(store.bright_stars_for_sector(conn, 3, 0, 1)) == 1

        store.clear_bright_stars(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 0
        assert store.bright_star_scatter_settings(conn) is None
    finally:
        conn.close()


def test_bright_star_web_queries(mysql_config):
    """The Stats page's placed/filled counts, the sector page's per-cell and
    per-box listings (unfilled only by default) and the Generate page's
    scatter status all agree with what was stored."""
    from planetgen.db import stats
    from planetgen.db import query as queryDb

    conn = store.get_connection(mysql_config)
    try:
        assert stats.bright_star_counts(conn) == {"placed": 0, "filled": 0, "unfilled": 0}
        status = queryDb.bright_star_scatter_status(conn)
        assert status["scattered"] is False and status["min_luminosity_sol"] is None and status["seed"] is None
        assert status["default_min_luminosity_sol"] == tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL == 3000.0
        assert queryDb.bright_stars_in_sector(conn, 3, 0, 1) == []
        with conn:
            conn.execute("INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc,"
                         " bulge_scale_radius_pc, bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude,"
                         " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
                         " expected_system_count_at_density_1, outer_ring_index)"
                         " VALUES (1, 1, 1, 1, 1, 2, 0.2, 0.3, 1, 0, 1, 4, 10, 5)")
            store.insert_bright_stars(conn, [_bright_row(luminosity_sol=600.0), _bright_row(),
                                           _bright_row(slot=2)], batch_size=2)
            store.record_bright_star_scatter(conn, 100.0, 9)
        listed = queryDb.bright_stars_in_sector(conn, 3, 0, 1)
        assert [s["luminosity_sol"] for s in listed] == pytest.approx([800.0, 600.0], rel=1e-2)
        assert listed[0]["x"] == pytest.approx(12.0) and listed[0]["yerkes_class"] == "V"
        assert listed[0]["system_id"] is None

        system, cfg = _named_system("Lanternfall")
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
            store.mark_bright_star_filled(conn, listed[0]["id"], system_id)
        assert stats.bright_star_counts(conn) == {"placed": 3, "filled": 1, "unfilled": 2}
        assert [s["id"] for s in queryDb.bright_stars_in_sector(conn, 3, 0, 1)] == [listed[1]["id"]]
        every = queryDb.bright_stars_in_sector(conn, 3, 0, 1, unfilled_only=False)
        assert [(s["id"], s["system_id"]) for s in every] == [(listed[0]["id"], system_id), (listed[1]["id"], None)]

        lo, hi = (11.0, 2.0, -1.0), (13.0, 4.0, 1.0)
        boxed = queryDb.galaxy_bright_stars_in_box(conn, lo, hi, 4.0)
        assert system_id in [s["system_id"] for s in boxed]
        assert queryDb.galaxy_bright_stars_in_box(conn, lo, hi, 4.0, unfilled_only=True) \
            == [s for s in boxed if s["system_id"] is None]

        status = queryDb.bright_star_scatter_status(conn)
        assert status == {"scattered": True, "min_luminosity_sol": 100.0, "seed": 9,
                          "default_min_luminosity_sol": 3000.0}

        # A re-scatter empties the table and restarts the ids, so the id
        # span still counts exactly.
        store.clear_bright_stars(conn)
        with conn:
            store.insert_bright_stars(conn, [_bright_row(slot=s) for s in range(5)])
        assert stats.bright_star_counts(conn) == {"placed": 5, "filled": 0, "unfilled": 5}
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"] == 5
    finally:
        conn.close()


def _bright_row_at(z_pc, luminosity_sol, slot=1):
    """A bright star `z_pc` above the plane, addressed as the scatter would."""
    from planetgen.db.query import layer_index_at, ring_index_at

    return (ring_index_at(12.37, 4.0), layer_index_at(z_pc, 4.0), slot, 12000, 3000, int(z_pc * 1000), "young", "B2V", "V", 1.4e31, 3.0e6, 22000.0,
            luminosity_sol * 3.828e26, 0.02, 0.03, 7.0, 0.03, 42)


def test_bright_stars_in_a_tall_box_span_the_heights_it_reaches(mysql_config):
    """GEN.117: a crowd of luminous stars on the plane no longer crowds out
    the dimmer old giants above and below it."""
    from planetgen.db import query as queryDb

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc,"
                         " bulge_scale_radius_pc, bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude,"
                         " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
                         " expected_system_count_at_density_1, outer_ring_index)"
                         " VALUES (1, 1, 1, 1, 1, 2, 0.2, 0.3, 1, 0, 1, 4, 10, 5)")
            # One star per layer (4 pc apart), so every address is its own.
            rows = [_bright_row_at(4.0 * n, 5000.0 + n) for n in range(60)]
            rows += [_bright_row_at(400.0 + 4.0 * n, 1500.0) for n in range(10)]
            rows += [_bright_row_at(-400.0 - 4.0 * n, 1400.0) for n in range(10)]
            store.insert_bright_stars(conn, rows)
        lo, hi = (0.0, 0.0, -1000.0), (100.0, 100.0, 1000.0)
        picked = queryDb.galaxy_bright_stars_in_box(conn, lo, hi, 4.0, limit=20)
        heights = sorted(star["z"] for star in picked)
        assert len(picked) == 20
        assert sum(z > 250 for z in heights) == 5 and sum(z < -250 for z in heights) == 5, heights
        assert [s["luminosity_sol"] for s in picked] == sorted((s["luminosity_sol"] for s in picked), reverse=True)
        # A box inside the plane's band picks the brightest as before.
        flat = queryDb.galaxy_bright_stars_in_box(conn, (0.0, 0.0, -100.0), (100.0, 100.0, 100.0), 4.0, limit=20)
        assert len(flat) == 20 and all(abs(s["z"]) < 100 for s in flat)
        # Few stars off the plane leave their room to the plane.
        tall = queryDb.galaxy_bright_stars_in_box(conn, (0.0, 0.0, 300.0), (100.0, 100.0, 1000.0), 4.0, limit=20)
        assert len(tall) == 10 and all(s["z"] > 250 for s in tall)
    finally:
        conn.close()


def test_every_population_is_listed_beside_the_brighter_young_stars(mysql_config):
    """A crowd of luminous young stars no longer crowds out the dimmer stars
    of the bulge and the thick disk, in the galaxy-wide sample or in a tile
    (whether the tile reads its own addresses or walks the luminosity
    index)."""
    from planetgen.db import query as queryDb

    def row(z_pc, luminosity_sol, population, slot):
        stored = list(_bright_row_at(z_pc, luminosity_sol, slot))
        stored[6] = population
        return tuple(stored)

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc,"
                         " bulge_scale_radius_pc, bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude,"
                         " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
                         " expected_system_count_at_density_1, outer_ring_index)"
                         " VALUES (1, 1, 1, 1, 1, 2, 0.2, 0.3, 1, 0, 1, 4, 10, 5)")
            rows = [row(4.0 * n, 50000.0 + n, "young", 1) for n in range(40)]
            rows += [row(4.0 * n, 1500.0 + n, "bulge", 2) for n in range(10)]
            rows += [row(-4.0 * n, 1400.0 + n, "old", 3) for n in range(10)]
            rows += [row(400.0 + 4.0 * n, 1300.0 + n, "old", 4) for n in range(10)]
            store.insert_bright_stars(conn, rows)

        def populations(stars):
            return sorted({star["population"] for star in stars})

        sample = queryDb.galaxy_brightest_stars(conn, 8)
        assert populations(sample) == ["bulge", "old", "young"], populations(sample)
        assert sum(star["population"] == "young" for star in sample) <= 4 + 4
        assert [s["luminosity_sol"] for s in sample] == sorted((s["luminosity_sol"] for s in sample), reverse=True)

        for lo, hi in [((0.0, 0.0, -200.0), (100.0, 100.0, 200.0)), ((-5000.0, -5000.0, -200.0), (5000.0, 5000.0, 200.0))]:
            picked = queryDb.galaxy_bright_stars_in_box(conn, lo, hi, 4.0, limit=12)
            assert len(picked) == 12
            assert {"young", "bulge", "old"} <= set(populations(picked)), (lo, populations(picked))
            assert sum(star["population"] == "young" for star in picked) <= 6
        # Past the row budget the galaxy-wide sample answers instead, still sharing the picks out.
        old_budget = queryDb.GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET
        queryDb.GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET = 0
        try:
            lo, hi = (0.0, 0.0, -200.0), (100.0, 100.0, 200.0)
            sampled = queryDb.galaxy_bright_stars_in_box(conn, lo, hi, 4.0, limit=12,
                                                         brightest=lambda: queryDb.galaxy_brightest_stars(conn, 60))
            assert {"young", "bulge", "old"} <= set(populations(sampled)), populations(sampled)
        finally:
            queryDb.GALAXY_TILE_BRIGHT_STAR_ROW_BUDGET = old_budget
    finally:
        conn.close()












