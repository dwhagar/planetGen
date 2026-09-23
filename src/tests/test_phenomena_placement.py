# tests/test_phenomena_placement.py

"""
Tests for schema v18's galaxy-frame placement of nebulae/asteroid fields --
`stellarObjects._db.compute_phenomenon_placement`/`insert_nebula`/
`insert_asteroid_field`/`save_phenomenon`, and the read side
(`queryDb.phenomena_near_sector`/`galaxy_placed_phenomena`) `html/
lib/starmap.py`'s Sector Map and `html/lib/galaxymap.py`'s Galaxy Map
build on. See `schema.sql`'s "v18" header note for the full design: a
nebula/asteroid field gets a real galaxy-frame sphere (`center_x/y/z_pc`,
`radius_ly`) centered near a given sector, not a sector-relative offset,
since it's frequently far larger than any one sector's own cube.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
they're skipped, not failed, when no MySQL test server is configured/
reachable.
"""

import math

import pytest

import queryDb
from stellarObjects import _db
from stellarObjects.asteroidFieldData import AsteroidField
from stellarObjects.config import SystemConfig
from stellarObjects.nebulaData import Nebula
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.utils import ly_to_pc, mpc_to_pc, pc_to_ly


def _place_sector(mysql_config, name, center_pc, edge_ly=11.5):
    """Persists an empty, galaxy-placed sector at an exact `center_pc` --
    a plain synthetic placement (no real Fibonacci-sphere/shell address),
    which is all `compute_phenomenon_placement`/`phenomena_near_sector`
    need (they only ever read `sectors.center_x/y/z_pc`/`edge_mpc`, never
    `shell_index`/`shell_slot_index`)."""
    sector = SpaceSector(name, edge_ly=edge_ly)
    cx, cy, cz = center_pc
    galaxy_position = {
        "center_x_pc": cx, "center_y_pc": cy, "center_z_pc": cz,
        "galactic_radius_pc": math.sqrt(cx * cx + cy * cy + cz * cz),
        "vertices_pc": {"inner": [], "outer": []},
    }
    return _db.save_sector(sector, config=mysql_config, galaxy_position=galaxy_position)


def test_compute_phenomenon_placement_center_falls_within_sector_cube_half_extent(mysql_config):
    sector_id = _place_sector(mysql_config, "Placement Test Sector", (100.0, -50.0, 25.0), edge_ly=11.5)
    edge_pc = ly_to_pc(11.5)

    conn = _db.get_connection(mysql_config)
    try:
        placement = _db.compute_phenomenon_placement(conn, sector_id)
    finally:
        conn.close()

    for axis_value, sector_center in zip(
        (placement["center_x_pc"], placement["center_y_pc"], placement["center_z_pc"]),
        (100.0, -50.0, 25.0),
    ):
        assert abs(axis_value - sector_center) <= edge_pc / 2 + 1e-9

    expected_radius = math.sqrt(
        placement["center_x_pc"] ** 2 + placement["center_y_pc"] ** 2 + placement["center_z_pc"] ** 2
    )
    assert placement["galactic_radius_pc"] == pytest.approx(expected_radius)


def test_compute_phenomenon_placement_rejects_an_unplaced_sector(mysql_config):
    sector = SpaceSector("Unplaced Sector", edge_ly=11.5)
    sector_id = _db.save_sector(sector, config=mysql_config)  # no galaxy_position

    conn = _db.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError, match="no galaxy placement"):
            _db.compute_phenomenon_placement(conn, sector_id)
    finally:
        conn.close()


def test_compute_phenomenon_placement_rejects_a_nonexistent_sector(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError):
            _db.compute_phenomenon_placement(conn, 999999)
    finally:
        conn.close()


def test_save_phenomenon_persists_placement_and_nearest_sector_for_nebula_and_asteroid_field(mysql_config):
    sector_id = _place_sector(mysql_config, "Nebula Home Sector", (10.0, 20.0, 30.0))

    cfg = SystemConfig()
    nebula = Nebula(cfg)
    field = AsteroidField(cfg)

    nebula_id = _db.save_phenomenon(nebula, cfg, "nebula", config=mysql_config, sector_id=sector_id)
    field_id = _db.save_phenomenon(field, cfg, "asteroid-field", config=mysql_config, sector_id=sector_id)

    conn = _db.get_connection(mysql_config)
    try:
        nebula_row = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        field_row = conn.execute("SELECT * FROM asteroid_fields WHERE id = ?", (field_id,)).fetchone()
    finally:
        conn.close()

    for row in (nebula_row, field_row):
        assert row["sector_id"] == sector_id
        assert row["center_x_pc"] is not None
        assert row["center_y_pc"] is not None
        assert row["center_z_pc"] is not None
        assert row["galactic_radius_pc"] == pytest.approx(
            math.sqrt(row["center_x_pc"] ** 2 + row["center_y_pc"] ** 2 + row["center_z_pc"] ** 2)
        )


def test_save_phenomenon_without_sector_id_leaves_it_unplaced(mysql_config):
    cfg = SystemConfig()
    nebula = Nebula(cfg)
    nebula_id = _db.save_phenomenon(nebula, cfg, "nebula", config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
    finally:
        conn.close()

    assert row["sector_id"] is None
    assert row["center_x_pc"] is None
    assert row["center_y_pc"] is None
    assert row["center_z_pc"] is None
    assert row["galactic_radius_pc"] is None


def test_phenomena_near_sector_finds_a_nebula_placed_at_that_sector(mysql_config):
    sector_id = _place_sector(mysql_config, "Nearby Sector", (200.0, 0.0, 0.0))

    cfg = SystemConfig()
    nebula = Nebula(cfg)
    nebula.radius_ly = 5.0  # small and deterministic, for a clean "definitely inside" case
    nebula_id = _db.save_phenomenon(nebula, cfg, "nebula", config=mysql_config, sector_id=sector_id)

    conn = _db.get_connection(mysql_config)
    try:
        matches = queryDb.phenomena_near_sector(conn, sector_id)
    finally:
        conn.close()

    assert len(matches) == 1
    match = matches[0]
    assert match["id"] == nebula_id
    assert match["type"] == "nebula"
    assert match["name"] == nebula.name
    assert match["descriptor"] == nebula.nebula_type
    assert match["radius_ly"] == pytest.approx(5.0)
    # offset_*_ly must actually be the phenomenon's real galaxy-frame
    # offset from the sector's own center, not a placeholder.
    expected_distance = math.sqrt(
        match["offset_x_ly"] ** 2 + match["offset_y_ly"] ** 2 + match["offset_z_ly"] ** 2
    )
    assert match["distance_ly"] == pytest.approx(expected_distance)


def test_phenomena_near_sector_excludes_a_small_far_nebula_but_includes_a_huge_one(mysql_config):
    # Two sectors far enough apart (500 pc ~= 1630 ly) that a small nebula
    # placed at sector A's own position can't plausibly reach sector B's
    # cube -- but a nebula whose own radius exceeds that separation must
    # still be found from B too (the whole point of a real sphere-vs-cube
    # overlap test rather than trusting the "nearest sector" link alone).
    sector_a = _place_sector(mysql_config, "Sector A", (0.0, 0.0, 0.0))
    sector_b = _place_sector(mysql_config, "Sector B", (500.0, 0.0, 0.0))

    cfg = SystemConfig()
    small_nebula = Nebula(cfg)
    small_nebula.radius_ly = 5.0
    small_id = _db.save_phenomenon(small_nebula, cfg, "nebula", config=mysql_config, sector_id=sector_a)

    huge_nebula = Nebula(cfg, name="Huge Nebula")
    huge_nebula.radius_ly = 3000.0  # bigger than the whole 500 pc (~1630 ly) separation
    huge_id = _db.save_phenomenon(huge_nebula, cfg, "nebula", config=mysql_config, sector_id=sector_a)

    conn = _db.get_connection(mysql_config)
    try:
        near_a = {m["id"] for m in queryDb.phenomena_near_sector(conn, sector_a)}
        near_b = {m["id"] for m in queryDb.phenomena_near_sector(conn, sector_b)}
    finally:
        conn.close()

    assert near_a == {small_id, huge_id}
    assert near_b == {huge_id}  # the small one is real but too far away


def test_phenomena_near_sector_is_empty_for_an_unplaced_sector(mysql_config):
    sector = SpaceSector("Unplaced Sector", edge_ly=11.5)
    sector_id = _db.save_sector(sector, config=mysql_config)

    conn = _db.get_connection(mysql_config)
    try:
        assert queryDb.phenomena_near_sector(conn, sector_id) == []
    finally:
        conn.close()


def test_galaxy_placed_phenomena_returns_only_placed_ones(mysql_config):
    sector_id = _place_sector(mysql_config, "Galaxy Dot Sector", (42.0, 42.0, 42.0))

    cfg = SystemConfig()
    placed_nebula = Nebula(cfg)
    placed_id = _db.save_phenomenon(placed_nebula, cfg, "nebula", config=mysql_config, sector_id=sector_id)
    unplaced_field = AsteroidField(cfg)
    _db.save_phenomenon(unplaced_field, cfg, "asteroid-field", config=mysql_config)  # no sector_id

    conn = _db.get_connection(mysql_config)
    try:
        placed = queryDb.galaxy_placed_phenomena(conn)
    finally:
        conn.close()

    ids = {(p["id"], p["type"]) for p in placed}
    assert (placed_id, "nebula") in ids
    assert len(placed) == 1  # the unplaced asteroid field must not appear

    entry = next(p for p in placed if p["id"] == placed_id and p["type"] == "nebula")
    assert entry["galactic_radius_pc"] == pytest.approx(
        math.sqrt(entry["x"] ** 2 + entry["y"] ** 2 + entry["z"] ** 2)
    )


def test_sector_detail_includes_nearby_phenomena(mysql_config):
    sector_id = _place_sector(mysql_config, "Detail Sector", (7.0, 8.0, 9.0))
    cfg = SystemConfig()
    field = AsteroidField(cfg)
    field_id = _db.save_phenomenon(field, cfg, "asteroid-field", config=mysql_config, sector_id=sector_id)

    conn = _db.get_connection(mysql_config)
    try:
        detail = queryDb.sector_detail(conn, sector_id)
    finally:
        conn.close()

    assert any(p["id"] == field_id and p["type"] == "asteroid_field" for p in detail["phenomena"])


def test_pc_ly_round_trip_used_by_placement_math():
    # Sanity-check the conversion helpers this module's own math leans on --
    # not new code, but a cheap guard against a future accidental change to
    # either direction silently breaking the placement/overlap math above.
    for value in (0.1, 1.0, 11.5, 500.0):
        assert pc_to_ly(ly_to_pc(value)) == pytest.approx(value)
        assert mpc_to_pc(value * 1000) == pytest.approx(value)


# ---------------------------------------------------------------------------
# v21: sector-level exotic phenomena persistence --
# stellarObjects._db._galaxy_placement_from_sector_offset/insert_black_hole/
# insert_neutron_star's new sector_id/placement params, insert_sector
# persisting SpaceSector.phenomena for all seven phenomenon types, and
# save_phenomenon now actually wiring sector_id through for
# supernova-remnant/rogue-planet/comet (previously always NULL). See
# schema.sql's "v21" header note.
# ---------------------------------------------------------------------------

from stellarObjects.compactRemnant import BlackHole, NeutronStar
from stellarObjects.roguePlanetData import InterstellarComet, RoguePlanet
from stellarObjects.supernovaRemnantData import SupernovaRemnant


def test_galaxy_placement_from_sector_offset_converts_ly_offset_to_absolute_pc():
    galaxy_position = {"center_x_pc": 10.0, "center_y_pc": -5.0, "center_z_pc": 2.0}
    offset_ly = (1.0, 0.0, 0.0)

    placement = _db._galaxy_placement_from_sector_offset(galaxy_position, offset_ly)

    expected_offset_pc = ly_to_pc(1.0)
    assert placement["center_x_pc"] == pytest.approx(10.0 + expected_offset_pc)
    assert placement["center_y_pc"] == pytest.approx(-5.0)
    assert placement["center_z_pc"] == pytest.approx(2.0)
    assert placement["galactic_radius_pc"] == pytest.approx(
        math.sqrt(placement["center_x_pc"] ** 2 + placement["center_y_pc"] ** 2 + placement["center_z_pc"] ** 2)
    )


def test_galaxy_placement_from_sector_offset_returns_none_for_an_unplaced_sector():
    assert _db._galaxy_placement_from_sector_offset(None, (1.0, 2.0, 3.0)) is None


def test_insert_black_hole_and_neutron_star_persist_sector_id_and_placement(mysql_config):
    sector_id = _place_sector(mysql_config, "Compact Remnant Home Sector", (10.0, 20.0, 30.0))
    placement = {
        "center_x_pc": 10.5, "center_y_pc": 20.0, "center_z_pc": 30.0,
        "galactic_radius_pc": math.sqrt(10.5 ** 2 + 20.0 ** 2 + 30.0 ** 2),
    }

    cfg = SystemConfig()
    bh = BlackHole(cfg)
    ns = NeutronStar(cfg)

    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            bh_id = _db.insert_black_hole(conn, bh, sector_id=sector_id, placement=placement)
            ns_id = _db.insert_neutron_star(conn, ns, sector_id=sector_id, placement=placement)

        bh_row = conn.execute("SELECT * FROM black_holes WHERE id = ?", (bh_id,)).fetchone()
        ns_row = conn.execute("SELECT * FROM neutron_stars WHERE id = ?", (ns_id,)).fetchone()
    finally:
        conn.close()

    for row in (bh_row, ns_row):
        assert row["sector_id"] == sector_id
        assert row["center_x_pc"] == pytest.approx(10.5)
        assert row["galactic_radius_pc"] == pytest.approx(placement["galactic_radius_pc"])


def test_insert_black_hole_without_sector_id_leaves_it_unplaced(mysql_config):
    bh = BlackHole(SystemConfig())
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            bh_id = _db.insert_black_hole(conn, bh)
        row = conn.execute("SELECT * FROM black_holes WHERE id = ?", (bh_id,)).fetchone()
    finally:
        conn.close()

    assert row["sector_id"] is None
    assert row["center_x_pc"] is None
    assert row["galactic_radius_pc"] is None


def test_save_phenomenon_now_wires_sector_id_for_supernova_remnant_rogue_planet_and_comet(mysql_config):
    # Previously a gap: save_phenomenon accepted sector_id but never passed
    # it through to insert_supernova_remnant/insert_rogue_planet/
    # insert_interstellar_comet -- these tables' own sector_id column
    # stayed permanently NULL. Now fixed; verify it actually persists.
    sector_id = _place_sector(mysql_config, "Reserved Column Sector", (1.0, 1.0, 1.0))
    cfg = SystemConfig()

    remnant = SupernovaRemnant(cfg)
    planet = RoguePlanet(cfg)
    comet = InterstellarComet(cfg)

    remnant_id = _db.save_phenomenon(remnant, cfg, "supernova-remnant", config=mysql_config, sector_id=sector_id)
    planet_id = _db.save_phenomenon(planet, cfg, "rogue-planet", config=mysql_config, sector_id=sector_id)
    comet_id = _db.save_phenomenon(comet, cfg, "comet", config=mysql_config, sector_id=sector_id)

    conn = _db.get_connection(mysql_config)
    try:
        remnant_row = conn.execute(
            "SELECT sector_id FROM supernova_remnants WHERE id = ?", (remnant_id,)
        ).fetchone()
        planet_row = conn.execute("SELECT sector_id FROM rogue_planets WHERE id = ?", (planet_id,)).fetchone()
        comet_row = conn.execute(
            "SELECT sector_id FROM interstellar_comets WHERE id = ?", (comet_id,)
        ).fetchone()
    finally:
        conn.close()

    assert remnant_row["sector_id"] == sector_id
    assert planet_row["sector_id"] == sector_id
    assert comet_row["sector_id"] == sector_id


def test_save_phenomenon_computes_real_placement_for_a_sector_linked_black_hole(mysql_config):
    # save_phenomenon's own compute_phenomenon_placement jitter path (for
    # phenomenonGen.py's standalone --sector-id use), now also reachable
    # for black-hole/neutron-star, not just nebula/asteroid-field.
    sector_id = _place_sector(mysql_config, "BH Sector-ID Sector", (50.0, 0.0, 0.0), edge_ly=11.5)
    edge_pc = ly_to_pc(11.5)

    bh = BlackHole(SystemConfig())
    bh_id = _db.save_phenomenon(bh, SystemConfig(), "black-hole", config=mysql_config, sector_id=sector_id)

    conn = _db.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT * FROM black_holes WHERE id = ?", (bh_id,)).fetchone()
    finally:
        conn.close()

    assert row["sector_id"] == sector_id
    assert abs(row["center_x_pc"] - 50.0) <= edge_pc / 2 + 1e-9


def test_insert_sector_persists_every_phenomenon_type_with_correct_placement(mysql_config):
    """
    Full pipeline test: a SpaceSector with a star system and all seven
    exotic phenomenon types, saved via save_sector -- confirms every type
    lands in its own table, linked by sector_id, with black-hole/
    neutron-star/nebula/asteroid-field additionally getting a real
    galaxy-frame position converted from their own sector-relative offset
    (not an independently re-randomized jitter).
    """
    from stellarObjects.asteroidFieldData import AsteroidField
    from stellarObjects.nebulaData import Nebula
    from stellarObjects.systemData import StarSystem

    sector = SpaceSector("Full Pipeline Sector", edge_ly=40.0)
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    system = StarSystem(system_config=cfg)
    sector.add_system(system, system_config=cfg)

    black_hole_entry = sector.add_phenomenon(BlackHole(SystemConfig()), "black-hole")
    neutron_star_entry = sector.add_phenomenon(NeutronStar(SystemConfig()), "neutron-star")
    nebula_entry = sector.add_phenomenon(Nebula(SystemConfig()), "nebula")
    field_entry = sector.add_phenomenon(AsteroidField(SystemConfig()), "asteroid-field")
    remnant_entry = sector.add_phenomenon(SupernovaRemnant(SystemConfig()), "supernova-remnant")
    planet_entry = sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    comet_entry = sector.add_phenomenon(InterstellarComet(SystemConfig()), "comet")

    galaxy_position = {
        "center_x_pc": 100.0, "center_y_pc": -40.0, "center_z_pc": 5.0,
        "galactic_radius_pc": math.sqrt(100.0 ** 2 + 40.0 ** 2 + 5.0 ** 2),
        "vertices_pc": {"inner": [], "outer": []},
    }
    sector_id = _db.save_sector(sector, config=mysql_config, galaxy_position=galaxy_position)

    conn = _db.get_connection(mysql_config)
    try:
        bh_row = conn.execute("SELECT * FROM black_holes WHERE sector_id = ?", (sector_id,)).fetchone()
        ns_row = conn.execute("SELECT * FROM neutron_stars WHERE sector_id = ?", (sector_id,)).fetchone()
        nebula_row = conn.execute("SELECT * FROM nebulae WHERE sector_id = ?", (sector_id,)).fetchone()
        field_row = conn.execute("SELECT * FROM asteroid_fields WHERE sector_id = ?", (sector_id,)).fetchone()
        remnant_row = conn.execute(
            "SELECT * FROM supernova_remnants WHERE sector_id = ?", (sector_id,)
        ).fetchone()
        planet_row = conn.execute("SELECT * FROM rogue_planets WHERE sector_id = ?", (sector_id,)).fetchone()
        comet_row = conn.execute(
            "SELECT * FROM interstellar_comets WHERE sector_id = ?", (sector_id,)
        ).fetchone()
    finally:
        conn.close()

    assert bh_row["name"] == black_hole_entry.phenomenon.name
    assert ns_row["name"] == neutron_star_entry.phenomenon.name
    assert nebula_row["name"] == nebula_entry.phenomenon.name
    assert field_row["name"] == field_entry.phenomenon.name
    assert remnant_row["name"] == remnant_entry.phenomenon.name
    assert planet_row["name"] == planet_entry.phenomenon.name
    assert comet_row["name"] == comet_entry.phenomenon.name

    # The four galaxy-placeable types get a real position converted from
    # their own sector-relative offset (galaxy_position + offset, NOT an
    # independent random jitter within the cube's half-extent, which is
    # what compute_phenomenon_placement's own jitter would give instead).
    for row, entry in (
        (bh_row, black_hole_entry), (ns_row, neutron_star_entry),
        (nebula_row, nebula_entry), (field_row, field_entry),
    ):
        expected = _db._galaxy_placement_from_sector_offset(galaxy_position, entry.position)
        assert row["center_x_pc"] == pytest.approx(expected["center_x_pc"])
        assert row["center_y_pc"] == pytest.approx(expected["center_y_pc"])
        assert row["center_z_pc"] == pytest.approx(expected["center_z_pc"])
        assert row["galactic_radius_pc"] == pytest.approx(expected["galactic_radius_pc"])


def test_insert_sector_links_phenomena_without_placement_for_an_unplaced_sector(mysql_config):
    # A sector generated via sectorGen.py's own standalone CLI (no galaxy
    # position) still links its phenomena by sector_id -- they just have
    # no galaxy-frame placement to compute.
    sector = SpaceSector("Unplaced Pipeline Sector", edge_ly=40.0)
    black_hole_entry = sector.add_phenomenon(BlackHole(SystemConfig()), "black-hole")

    sector_id = _db.save_sector(sector, config=mysql_config)  # no galaxy_position

    conn = _db.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT * FROM black_holes WHERE sector_id = ?", (sector_id,)).fetchone()
    finally:
        conn.close()

    assert row is not None
    assert row["name"] == black_hole_entry.phenomenon.name
    assert row["center_x_pc"] is None


def test_galaxy_placed_phenomena_includes_black_holes_and_neutron_stars(mysql_config):
    sector_id = _place_sector(mysql_config, "Compact Remnant Galaxy Dot Sector", (7.0, 7.0, 7.0))

    bh = BlackHole(SystemConfig())
    bh.has_accretion_disk = True
    ns = NeutronStar(SystemConfig())

    bh_id = _db.save_phenomenon(bh, SystemConfig(), "black-hole", config=mysql_config, sector_id=sector_id)
    ns_id = _db.save_phenomenon(ns, SystemConfig(), "neutron-star", config=mysql_config, sector_id=sector_id)

    conn = _db.get_connection(mysql_config)
    try:
        placed = queryDb.galaxy_placed_phenomena(conn)
    finally:
        conn.close()

    by_id_and_type = {(p["id"], p["type"]): p for p in placed}
    assert (bh_id, "black_hole") in by_id_and_type
    assert (ns_id, "neutron_star") in by_id_and_type
    assert by_id_and_type[(bh_id, "black_hole")]["descriptor"] == "accreting"
    assert by_id_and_type[(bh_id, "black_hole")]["radius_ly"] == 0
    assert by_id_and_type[(ns_id, "neutron_star")]["descriptor"] == ns.pulsar_type


def test_phenomena_near_sector_finds_a_black_hole_placed_at_that_sector(mysql_config):
    sector_id = _place_sector(mysql_config, "Compact Remnant Nearby Sector", (300.0, 0.0, 0.0))

    bh = BlackHole(SystemConfig())
    bh_id = _db.save_phenomenon(bh, SystemConfig(), "black-hole", config=mysql_config, sector_id=sector_id)

    conn = _db.get_connection(mysql_config)
    try:
        matches = queryDb.phenomena_near_sector(conn, sector_id)
    finally:
        conn.close()

    match = next(m for m in matches if m["id"] == bh_id and m["type"] == "black_hole")
    assert match["name"] == bh.name
    assert match["radius_ly"] == 0
    # A point object (radius_ly=0) still shows up when its own center
    # falls within the sector's bounding sphere -- compute_phenomenon_placement's
    # own jitter always keeps it within the cube's half-extent.
    assert match["distance_ly"] >= 0


# ---------------------------------------------------------------------------
# supernova_remnant support in list_phenomena/count_phenomena/
# phenomenon_detail (html/phenomena.py's listing, html/phenomenon.py's
# detail page) -- added alongside this diagram's own AU-scale feature,
# since supernova_remnants previously had no web page at all. It has its
# own real radius_ly (unlike black_hole/neutron_star) but NO galaxy-frame
# placement columns (unlike the other four phenomenon types), which is
# exactly what queryDb._SUPERNOVA_REMNANT_TABLE is kept separate from
# _PHENOMENON_TABLES for -- see that constant's own docstring, and
# test_navigation.py's test_nav_between_rejects_a_real_supernova_remnant_endpoint
# for the NAV side of that same distinction.
# ---------------------------------------------------------------------------

def test_list_phenomena_includes_a_supernova_remnant_always_unplaced(mysql_config):
    remnant = SupernovaRemnant(SystemConfig())
    conn = _db.get_connection(mysql_config)
    try:
        remnant_id = _db.insert_supernova_remnant(conn, remnant)
        conn.commit()
        rows = queryDb.list_phenomena(conn)
    finally:
        conn.close()

    row = next(r for r in rows if r["id"] == remnant_id and r["type"] == "supernova_remnant")
    assert row["name"] == remnant.name
    assert row["descriptor"] == remnant.morphology
    assert row["radius_ly"] == pytest.approx(remnant.radius_ly)
    assert row["placed"] is False


def test_count_phenomena_includes_supernova_remnants(mysql_config):
    # A core-collapse progenitor can leave a detectable compact remnant of
    # its own (SupernovaRemnant.__init__'s own probabilistic roll), which
    # would insert a SECOND, genuinely-standalone-shaped black_holes/
    # neutron_stars row (star_id NULL) that count_phenomena's own
    # black_holes/neutron_stars branch can't distinguish from a real
    # standalone one -- retry for a guaranteed Type Ia progenitor (never
    # leaves anything behind) so this test's "+1" is deterministic.
    remnant = SupernovaRemnant(SystemConfig())
    for _ in range(50):
        if remnant.compact_remnant is None:
            break
        remnant = SupernovaRemnant(SystemConfig())
    else:
        pytest.fail("could not generate a supernova remnant with no embedded compact remnant")

    conn = _db.get_connection(mysql_config)
    try:
        before = queryDb.count_phenomena(conn)
        _db.insert_supernova_remnant(conn, remnant)
        conn.commit()
        after = queryDb.count_phenomena(conn)
    finally:
        conn.close()

    assert after == before + 1


def test_phenomenon_detail_returns_a_supernova_remnants_own_columns(mysql_config):
    remnant = SupernovaRemnant(SystemConfig())
    conn = _db.get_connection(mysql_config)
    try:
        remnant_id = _db.insert_supernova_remnant(conn, remnant)
        conn.commit()
        detail = queryDb.phenomenon_detail(conn, "supernova_remnant", remnant_id)
    finally:
        conn.close()

    assert detail["name"] == remnant.name
    assert detail["morphology"] == remnant.morphology
    assert detail["progenitor_type"] == remnant.progenitor_type
    assert detail["radius_ly"] == pytest.approx(remnant.radius_ly)
    assert detail["type"] == "supernova_remnant"
    assert detail["sector_name"] is None
