# tests/test_galactic_motion.py

"""The correlative update's galactic motion (GEN.6): everything
moves along its galactic orbit and is refiled into the sector it drifts
into, with octant, location, containment and nearest systems following."""

import math

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.frame import local_to_galaxy_pc
from planetgen.galaxy.spherical import band_center_pc, sector_address_at, sector_position_pc, wedge_longitude_bounds
from planetgen.galaxy.galactic_orbit import calculate_galactic_orbit, circular_orbit
from planetgen.galaxy.sector import SpaceSector
from planetgen.galaxy.straight_line import scatter_position_mpc, straight_line_position_pc
from planetgen.generation.system import StarSystem
from planetgen.physics import constants as physical_constants
from planetgen.physics.position import update_interval_s
from planetgen.physics.units import pc_to_ly
from tests.toy_shape import toy_shape

EDGE_PC = 4.0
SHELL = 10


def test_rotation_follows_galactic_phase():
    x, y, z = store._rotate_about_axis((10.0, 0.0, 1.0), math.pi / 2)
    assert (x, y, z) == pytest.approx((0.0, 10.0, 1.0))
    assert store._galactic_turn(125e6, 0.25) == pytest.approx(math.pi)
    assert store._galactic_turn(250e6, 0.25) == 0.0  # a whole orbit lands where it started (GEN.108)
    assert store._galactic_turn(1.0, None) == 0.0


def _place(wedge, positions_ly, name):
    center = sector_position_pc(SHELL, 0, wedge, EDGE_PC)
    sector = SpaceSector(name, edge_ly=pc_to_ly(EDGE_PC))
    for position in positions_ly:
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = False
        cfg.BINARY_SYSTEM = False
        sector.add_system(StarSystem(system_config=cfg), position=position, system_config=cfg)
    return sector, {"center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
                    "galactic_radius_pc": math.hypot(*center), "shell_index": SHELL, "band_index": 0,
                    "wedge_index": wedge}


def test_a_system_that_drifts_across_a_boundary_is_refiled(mysql_config):
    # A system just short of wedge 0's leading edge (local +Y is the
    # direction of travel), and a second sector ahead of it in wedge 1.
    start, end = wedge_longitude_bounds(SHELL, 0, 0)
    half_ly = pc_to_ly(band_center_pc(SHELL, 0, EDGE_PC)[0] * math.sin((end - start) / 2))
    sector, position = _place(0, [(0.0, half_ly - 0.3, 0.0), (0.0, -half_ly + 0.5, 0.0)], "Trailing")
    first = store.save_sector(sector, config=mysql_config, galaxy_position=position)
    sector, position = _place(1, [(0.0, 0.0, 0.0)], "Leading")
    second = store.save_sector(sector, config=mysql_config, galaxy_position=position)

    conn = store.get_connection(mysql_config)
    try:
        mover, stay = [r["id"] for r in conn.execute(
            "SELECT id FROM star_systems WHERE sector_id = ? ORDER BY position_y_mpc DESC", (first,)).fetchall()]
        # Run far enough to carry the mover about a light-year forward at the speed it stores (the stepper moves a system by its velocity, GEN.340).
        speed_kms = store._speed_kms(conn.execute("SELECT velocity_x_kms, velocity_y_kms, velocity_z_kms FROM star_systems WHERE id = ?", (mover,)).fetchone())
        elapsed = 1.0 * physical_constants.LIGHTYEAR_M / (speed_kms * 1000.0) / physical_constants.SECONDS_PER_YEAR
        with conn:
            motion = store.advance_galactic_positions(conn, store.orbit_clock(conn, elapsed))
            rewritten = store.refresh_after_motion(conn, motion["sectors"])["locations"]
        assert motion["refiled"] >= 1 and {first, second} <= motion["sectors"]
        row = conn.execute("SELECT * FROM star_systems WHERE id = ?", (mover,)).fetchone()
        assert row["sector_id"] == second
        assert row["location"].startswith("Leading")
        assert rewritten >= 1
        # It sits where the turn put it, in its new sector's frame.
        center = sector_position_pc(SHELL, 0, 1, EDGE_PC)
        point = local_to_galaxy_pc(center, (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000,
                                            row["position_z_mpc"] / 1000))
        low, high = wedge_longitude_bounds(SHELL, 0, 1)
        assert low <= math.atan2(point[1], point[0]) <= high
        assert conn.execute("SELECT sector_id FROM star_systems WHERE id = ?", (stay,)).fetchone()["sector_id"] == first
        near = conn.execute("SELECT sector_id FROM nearest_systems WHERE object_table = 'star_systems'"
                            " AND object_id = ? LIMIT 1", (mover,)).fetchone()
        assert near["sector_id"] == second
    finally:
        conn.close()


def test_phenomena_and_facilities_move_and_the_quasar_stays(mysql_config):
    sector, position = _place(3, [(0.0, 0.0, 0.0)], "Drift")
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position=position)
    conn = store.get_connection(mysql_config)
    try:
        from planetgen.generation.phenomena.nebula import Nebula
        center = (position["center_x_pc"], position["center_y_pc"], position["center_z_pc"])
        with conn:
            nebula = Nebula(SystemConfig())
            nebula_id = store.insert_nebula(conn, nebula, sector_id=sector_id, placement={
                "center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
                "galactic_radius_pc": math.hypot(*center)})
            facility_id = store.add_facility(conn, "Drifter", "station", "standalone", "space", sector_id)
            speed_kms = store._speed_kms(conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone())
            motion = store.advance_galactic_positions(conn, store.orbit_clock(conn, 1e3))
        assert motion["moved"] >= 3
        row = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        # A thousand years is a short way: the stepper carries the nebula its speed times the time (GEN.340; it was the angle its period sweeps).
        flown = math.dist((row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]), center)
        assert flown == pytest.approx(speed_kms * physical_constants.SECONDS_PER_YEAR * 1e3 / physical_constants.PARSEC_M * 1000.0, rel=0.02)
        moved = conn.execute("SELECT center_x_pc FROM facilities WHERE id = ?", (facility_id,)).fetchone()
        assert moved["center_x_pc"] != pytest.approx(center[0], abs=1e-12)
    finally:
        conn.close()


def test_a_standalone_facility_stores_the_rotation_curves_velocity_and_turns_it_with_its_position(mysql_config):
    """GEN.125: a stand-alone facility has a galactic velocity, the rotation curve's tangent at its
    place, which `advance_galactic_positions` turns by the same angle as its position."""
    sector, position = _place(3, [(0.0, 0.0, 0.0)], "Drift")
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position=position)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            facility_id = store.add_facility(conn, "Drifter", "station", "standalone", "space", sector_id)
        before = conn.execute("SELECT center_x_pc, center_y_pc, velocity_x_kms, velocity_y_kms, velocity_z_kms"
                              " FROM facilities WHERE id = ?", (facility_id,)).fetchone()
        speed = math.hypot(before["velocity_x_kms"], before["velocity_y_kms"])
        radius = math.hypot(before["center_x_pc"], before["center_y_pc"])
        expected, _period = calculate_galactic_orbit(pc_to_ly(radius))
        assert speed == pytest.approx(expected) and speed > 0.0 and before["velocity_z_kms"] == 0.0
        # Counterclockwise from galactic north: perpendicular to the radius, a positive turn.
        assert before["center_x_pc"] * before["velocity_x_kms"] + before["center_y_pc"] * before["velocity_y_kms"] \
            == pytest.approx(0.0, abs=1e-6 * speed * radius)
        assert before["center_x_pc"] * before["velocity_y_kms"] - before["center_y_pc"] * before["velocity_x_kms"] > 0
        with conn:
            store.advance_galactic_positions(conn, store.orbit_clock(conn, 1e7))
        after = conn.execute("SELECT center_x_pc, center_y_pc, velocity_x_kms, velocity_y_kms, velocity_z_kms"
                             " FROM facilities WHERE id = ?", (facility_id,)).fetchone()
        position_turn = math.atan2(after["center_y_pc"], after["center_x_pc"]) - math.atan2(before["center_y_pc"],
                                                                                           before["center_x_pc"])
        velocity_turn = math.atan2(after["velocity_y_kms"], after["velocity_x_kms"]) - math.atan2(
            before["velocity_y_kms"], before["velocity_x_kms"])
        assert position_turn != 0.0 and velocity_turn == pytest.approx(position_turn)
        assert math.hypot(after["velocity_x_kms"], after["velocity_y_kms"]) == pytest.approx(speed)
    finally:
        conn.close()


def test_facilities_on_a_body_keep_a_zero_galactic_velocity(mysql_config):
    from tests.test_facilities import _ids, _system_with_worlds
    system_id = _system_with_worlds(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        ids = _ids(conn, system_id)
        with conn:
            facility_id = store.add_facility(conn, "Ring", "station", "orbital", "planet", ids["giant"], phase_deg=0.0)
        row = conn.execute("SELECT velocity_x_kms, velocity_y_kms, velocity_z_kms FROM facilities WHERE id = ?",
                           (facility_id,)).fetchone()
        assert (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]) == (0.0, 0.0, 0.0)
    finally:
        conn.close()


def test_orbital_facilities_advance_their_phase(mysql_config):
    from tests.test_facilities import _ids, _system_with_worlds
    system_id = _system_with_worlds(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        ids = _ids(conn, system_id)
        with conn:
            facility_id = store.add_facility(conn, "Ring", "station", "orbital", "planet", ids["giant"], phase_deg=0.0)
        row = conn.execute("SELECT orbit_period_years, orbit_distance_km FROM facilities WHERE id = ?",
                           (facility_id,)).fetchone()
        period = row["orbit_period_years"]
        # A facility is due once it has moved its scale's threshold (100,000 km round a planet), so one close round a
        # small giant is not due after a quarter turn: the step is whole turns more than that takes, plus the quarter
        # (the whole turns wrap away, GEN.108), which keeps the phase it checks (TEST.251).
        speed_ms = row["orbit_distance_km"] * 1000.0 * 2 * math.pi / (period * physical_constants.SECONDS_PER_YEAR)
        due_years = update_interval_s(speed_ms, "planetary") / physical_constants.SECONDS_PER_YEAR
        step = (math.ceil(due_years / period) + 0.25) * period
        with conn:
            assert store.advance_facility_orbits(conn, store.orbit_clock(conn, step)) == 1
        phase = conn.execute("SELECT orbit_phase_deg FROM facilities WHERE id = ?", (facility_id,)).fetchone()
        assert phase["orbit_phase_deg"] == pytest.approx(90.0)
    finally:
        conn.close()


def test_update_orbits_runs_end_to_end(mysql_config, monkeypatch, capsys):
    """The whole correlative update against a real current-schema
    database: a generated sector with planets and phenomena, run twice
    (the first run only sets the clock)."""
    from planetgen.cli import orbits as updateOrbits

    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    sector = SpaceSector("Clockwork", edge_ly=pc_to_ly(EDGE_PC))
    sector.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    _sector, position = _place(5, [], "unused")
    store.save_sector(sector, config=mysql_config, galaxy_position=position)

    argv = ["planetgen.cli.orbits", "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
            "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password or "",
            "--mysql-database", mysql_config.database]
    monkeypatch.setattr("sys.argv", argv)
    updateOrbits.main()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            # Back-date the last run and the system's own clock: due now.
            conn.execute("UPDATE orbit_simulation_state SET last_updated_at = NOW() - INTERVAL 400 DAY")
            conn.execute("UPDATE star_systems SET epoch_unix = NULL, next_update_due = NULL")
    finally:
        conn.close()
    updateOrbits.main()
    out = capsys.readouterr().out
    assert "Moved 1 object(s) along their galactic orbits; 0 changed sector." in out


def test_a_hypervelocity_star_flies_a_straight_line_while_a_bound_star_orbits(mysql_config):
    sector, position = _place(3, [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], "Runaway")
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position=position)
    center = (position["center_x_pc"], position["center_y_pc"], position["center_z_pc"])
    velocity = (30.0, -40.0, 5.0)
    elapsed = 1e5
    conn = store.get_connection(mysql_config)
    try:
        flyer, bound = [r["id"] for r in conn.execute(
            "SELECT id FROM star_systems WHERE sector_id = ? ORDER BY id", (sector_id,)).fetchall()]
        with conn:
            conn.execute("UPDATE stars SET galactic_orbital_period_gy = 1.0")
            conn.execute("UPDATE star_systems SET runaway_class = 'hypervelocity', runaway_speed_kms = 50.0,"
                         " velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ? WHERE id = ?", (*velocity, flyer))
            start = conn.execute("SELECT position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems"
                                 " WHERE id = ?", (flyer,)).fetchone()
            motion = store.advance_galactic_positions(conn, store.orbit_clock(conn, elapsed))
        assert motion["moved"] >= 2
        row = conn.execute("SELECT * FROM star_systems WHERE id = ?", (flyer,)).fetchone()
        assert (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]) == pytest.approx(velocity)
        orbiting = conn.execute("SELECT velocity_x_kms, velocity_y_kms FROM star_systems WHERE id = ?",
                                (bound,)).fetchone()
        assert (orbiting["velocity_x_kms"], orbiting["velocity_y_kms"]) != pytest.approx(velocity[:2])
        found = conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors WHERE id = ?",
                             (row["sector_id"],)).fetchone()
        now = local_to_galaxy_pc(sector_position_pc(found["shell_index"], found["band_index"],
                                                    found["wedge_index"], EDGE_PC),
                                 (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000,
                                  row["position_z_mpc"] / 1000))
    finally:
        conn.close()
    origin = local_to_galaxy_pc(center, tuple(start[k] / 1000 for k in
                                              ("position_x_mpc", "position_y_mpc", "position_z_mpc")))
    expected = straight_line_position_pc(origin, velocity, elapsed * physical_constants.SECONDS_PER_YEAR)
    assert now == pytest.approx(expected, abs=1e-3)
    assert math.dist(now, origin) > 0.5  # it did go somewhere


def test_straight_line_position_is_the_start_plus_velocity_times_time():
    km_per_pc = physical_constants.PARSEC_M / 1000.0
    assert straight_line_position_pc((1.0, 2.0, 3.0), (km_per_pc, 0.0, -km_per_pc), 10.0) == pytest.approx(
        (11.0, 2.0, -7.0))
    assert straight_line_position_pc((1.0, 2.0, 3.0), (5.0, 5.0, 5.0), 0.0) == (1.0, 2.0, 3.0)
    assert straight_line_position_pc((0.0, 0.0, 0.0), (km_per_pc, 0.0, 0.0), -3.0) == pytest.approx((-3.0, 0.0, 0.0))


def test_a_scatter_row_flies_on_from_its_epoch_and_other_kinds_stay():
    row = {"kind": "hypervelocity-star", "position_x_mpc": 1000.0, "position_y_mpc": 0.0, "position_z_mpc": 0.0,
           "velocity_x_kms": physical_constants.PARSEC_M / 1000.0, "velocity_y_kms": 0.0, "velocity_z_kms": 0.0,
           "epoch_unix": 100.0}
    assert scatter_position_mpc(row, 105.0) == pytest.approx((6000.0, 0.0, 0.0))
    assert scatter_position_mpc(dict(row, epoch_unix=None), 105.0) == (1000.0, 0.0, 0.0)
    assert scatter_position_mpc(dict(row, kind="neutron-star"), 105.0) == (1000.0, 0.0, 0.0)


SCATTER_RADIUS_PC = 8000.0


def _seed_scatter_grid(mysql_config):
    """The galaxy's stored cell size (a point is refiled by it) and one point of every moving scatter kind."""
    shape = toy_shape(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
                               bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=2000.0,
                            config=mysql_config)


def _insert_scatter_point(conn, table, point_pc, kind=None, velocity=None, built=False):
    shell, band, wedge = sector_address_at(point_pc, EDGE_PC)
    x, y, z = (c * 1000.0 for c in point_pc)
    v = velocity or (None, None, None)
    if table == "bright_stars":
        return conn.execute(
            "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc,"
            " position_z_mpc, velocity_x_kms, velocity_y_kms, velocity_z_kms, population, star_type, yerkes_class,"
            " mass_kg, radius_km, temperature_k, luminosity_w, age_gy, initial_mass_sol, seed)"
            " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 'young', 'B2V', 'V', 1e31, 1e6, 20000, 1e30, 0.1, 5.0, 1)",
            (shell, band, wedge, x, y, z, *v)).lastrowid
    return conn.execute(
        "INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, position_x_mpc, position_y_mpc,"
        " position_z_mpc, velocity_x_kms, velocity_y_kms, velocity_z_kms, seed, built_at)"
        f" VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1, {'CURRENT_TIMESTAMP' if built else 'NULL'})",
        (shell, band, wedge, kind, x, y, z, *v)).lastrowid


def test_every_kind_of_scatter_point_moves_and_is_refiled_and_built_or_nucleus_points_stay(mysql_config):
    """GEN.223: a scatter point is a point in space like any other; nothing in the scatter is left standing."""
    _seed_scatter_grid(mysql_config)
    point = (SCATTER_RADIUS_PC, 0.0, 0.0)
    elapsed = 1e6
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            ids = {
                "bright_stars": _insert_scatter_point(conn, "bright_stars", point),
                "black-hole": _insert_scatter_point(conn, "phenomenon_scatter", point, "black-hole"),
                "neutron-star": _insert_scatter_point(conn, "phenomenon_scatter", point, "neutron-star"),
            }
            built = _insert_scatter_point(conn, "phenomenon_scatter", point, "black-hole", built=True)
            nucleus = _insert_scatter_point(conn, "phenomenon_scatter", point, "quasar")
            motion = store.advance_scatter_points(conn, store.orbit_clock(conn, elapsed))
            pot = store._galactic_potential(conn)
            rows = {
                "bright_stars": conn.execute("SELECT * FROM bright_stars WHERE id = ?",
                                             (ids["bright_stars"],)).fetchone(),
                **{kind: conn.execute("SELECT * FROM phenomenon_scatter WHERE id = ?", (ids[kind],)).fetchone()
                   for kind in ("black-hole", "neutron-star")},
            }
            stayed = [conn.execute("SELECT position_x_mpc, epoch_unix FROM phenomenon_scatter WHERE id = ?",
                                   (row_id,)).fetchone() for row_id in (built, nucleus)]
    finally:
        conn.close()
    assert motion["counts"] == {"bright_stars": 1, "phenomenon_scatter": 2} and motion["moved"] == 3
    _speed, period_gy = circular_orbit(pot, pc_to_ly(SCATTER_RADIUS_PC))
    angle = 2 * math.pi * (elapsed / (period_gy * 1e9))
    for name, row in rows.items():
        moved = (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0, row["position_z_mpc"] / 1000.0)
        assert moved == pytest.approx((point[0] * math.cos(angle), point[0] * math.sin(angle), 0.0), abs=0.5), name   # the stepper's circle, not the exact turn: a part in 1e4 of the radius
        assert row["epoch_unix"] is not None and row["next_update_due"] is not None, name
        assert math.hypot(row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]) > 1.0, name
        assert sector_address_at(moved, EDGE_PC) == (row["shell_index"], row["band_index"], row["wedge_index"]), name
        assert (row["shell_index"], row["band_index"], row["wedge_index"]) != sector_address_at(point, EDGE_PC), name
    assert all(row["position_x_mpc"] == pytest.approx(point[0] * 1000.0) and row["epoch_unix"] is None
               for row in stayed)


def test_a_hypervelocity_scatter_point_flies_a_straight_line(mysql_config):
    _seed_scatter_grid(mysql_config)
    point = (SCATTER_RADIUS_PC, 0.0, 0.0)
    velocity = (400.0, -300.0, 20.0)
    elapsed = 1e5
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "phenomenon_scatter", point, "hypervelocity-star", velocity)
            store.advance_scatter_points(conn, store.orbit_clock(conn, elapsed))
            row = conn.execute("SELECT * FROM phenomenon_scatter WHERE id = ?", (row_id,)).fetchone()
    finally:
        conn.close()
    expected = straight_line_position_pc(point, velocity, elapsed * physical_constants.SECONDS_PER_YEAR)
    assert (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0,
            row["position_z_mpc"] / 1000.0) == pytest.approx(expected, abs=1e-3)
    assert (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]) == pytest.approx(velocity)


def test_a_systems_galaxy_point_and_cell_key_follow_every_save_and_move(mysql_config):
    """DB.34: `star_systems.galaxy_*_pc` is the sector's center plus the offset, and `cell_key` (MySQL's) follows it."""
    from planetgen.galaxy import cellkey
    sector, position = _place(3, [(0.0, 0.0, 0.0), (0.3, -0.2, 0.1)], "Keyed")
    store.save_sector(sector, config=mysql_config, galaxy_position=position)
    center = (position["center_x_pc"], position["center_y_pc"], position["center_z_pc"])

    def check(conn):
        rows = conn.execute("SELECT sector_id, position_x_mpc, position_y_mpc, position_z_mpc, galaxy_x_pc, galaxy_y_pc,"
                            " galaxy_z_pc, cell_key FROM star_systems").fetchall()
        assert len(rows) == 2
        for row in rows:
            owner = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc FROM sectors WHERE id = ?",
                                 (row["sector_id"],)).fetchone()
            expected = local_to_galaxy_pc((owner["center_x_pc"], owner["center_y_pc"], owner["center_z_pc"]),
                                          tuple(row[f"position_{a}_mpc"] / 1000.0 for a in "xyz"))
            point = (row["galaxy_x_pc"], row["galaxy_y_pc"], row["galaxy_z_pc"])
            assert point == pytest.approx(expected, abs=1e-6)
            assert row["cell_key"] == cellkey.cell_key(*point)

    conn = store.get_connection(mysql_config)
    try:
        check(conn)
        before = conn.execute("SELECT galaxy_x_pc FROM star_systems ORDER BY id").fetchall()
        with conn:
            conn.execute("UPDATE stars SET galactic_orbital_period_gy = 1.0")
            store.advance_galactic_positions(conn, store.orbit_clock(conn, 1e6))
        check(conn)
        after = conn.execute("SELECT galaxy_x_pc FROM star_systems ORDER BY id").fetchall()
        assert [r["galaxy_x_pc"] for r in after] != pytest.approx([r["galaxy_x_pc"] for r in before], abs=1e-9)
    finally:
        conn.close()
    assert center  # the sector was placed


def test_a_neutron_star_the_galaxy_ejected_flies_in_a_straight_line_and_the_nucleus_stays(mysql_config):
    """GEN.238: a black hole or neutron star with a velocity and no galactic orbit moves like a hypervelocity star, by
    its velocity from its epoch, and keeps the velocity; the galaxy's own black hole has neither and stays."""
    from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
    sector, position = _place(3, [(0.0, 0.0, 0.0)], "Ejected")
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position=position)
    center = (position["center_x_pc"], position["center_y_pc"], position["center_z_pc"])
    placement = {"center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
                 "galactic_radius_pc": math.hypot(*center)}
    velocity = (300.0, -120.0, 200.0)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            neutron_star = NeutronStar(SystemConfig(), galactic_center_dist_ly=pc_to_ly(math.hypot(*center)))
            neutron_star.eject(velocity)
            star_id = store.insert_neutron_star(conn, neutron_star, sector_id=sector_id, placement=placement)
            hole = BlackHole(SystemConfig(), mass_class="supermassive")
            hole_id = store.insert_black_hole(conn, hole, sector_id=sector_id, placement=placement)
            store.advance_galactic_positions(conn, store.orbit_clock(conn, 1e6))
        row = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc, velocity_x_kms, velocity_y_kms, velocity_z_kms,"
                           " galactic_orbital_period_gy FROM neutron_stars WHERE id = ?", (star_id,)).fetchone()
        assert row["galactic_orbital_period_gy"] is None
        expected = straight_line_position_pc(center, velocity, 1e6 * physical_constants.SECONDS_PER_YEAR)
        assert (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]) == pytest.approx(expected, rel=1e-9)
        assert math.dist(expected, center) > 100.0       # 300 km/s for a million years is hundreds of parsecs
        assert (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]) == pytest.approx(velocity)
        still = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc FROM black_holes WHERE id = ?",
                             (hole_id,)).fetchone()
        assert (still["center_x_pc"], still["center_y_pc"], still["center_z_pc"]) == pytest.approx(center)
    finally:
        conn.close()
