# tests/test_sky_sphere.py

"""
VIEW.22 (System sky S4): `planetgen.db.sky_sphere.stars_in_sphere`, the stars of a sphere above a luminosity floor from the charted systems, the
pre-placed bright stars, and the scattered phenomena, as arrays with the cell levels the sphere reached.

The rows are written straight into the tables (no generation), so each test knows exactly what is stored, and the answer is checked against a brute
force over the same numbers: every row, its distance worked out by `math.dist`, its luminosity against the floor.
"""

import math
import random
import time

import numpy as np
import pytest

from planetgen.db import sky_sphere, store
from planetgen.db.sky_sphere import BRIGHT, CHARTED, SCATTER, SkySphereError
from planetgen.galaxy import spherical
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.physics.units import ly_to_pc
from planetgen import tuning
from tests import bright_star_rows
from tests.toy_shape import toy_shape

SOL = 3.828e26
EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)


@pytest.fixture
def conn(mysql_config):
    connection = store.get_connection(mysql_config)
    try:
        yield connection
    finally:
        connection.close()


def add_system(conn, point, stars, velocity=(1.0, 2.0, 3.0), epoch=1.0e9, name="S"):
    """A charted system at `point` (parsecs) with one star for each `(luminosity_sol, temperature_k)` of `stars`; returns `(system id, [star ids])`."""
    config_id = conn.execute("INSERT INTO system_configs (name) VALUES (NULL)").lastrowid
    cur = conn.execute("INSERT INTO star_systems (system_config_id, name, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc, velocity_x_kms, velocity_y_kms,"
                 " velocity_z_kms, epoch_unix) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", (config_id, name, *point, *velocity, epoch))
    system_id = cur.lastrowid
    ids = []
    for index, (luminosity, temperature) in enumerate(stars):
        cur = conn.execute(
            "INSERT INTO stars (star_system_id, role, name, star_type, yerkes_class, mass_kg, radius_km, temperature_k, luminosity_w, age_gy, lifespan_gy,"
            " habitable_zone_inner_km, habitable_zone_outer_km, system_perimeter_km, heliosphere_radius_km, galactic_orbital_period_gy,"
            " galactic_orbital_phase_deg, galactic_min_update_interval_years, snow_line_km) VALUES (?, 'single', ?, 'G2V', 'V', 2e30, 7e5, ?, ?, 4.0,"
            " 10.0, 1.0, 2.0, 3.0, 4.0, 0.2, 0.0, 1.0, 1.0)", (system_id, f"{name} {index}", temperature, luminosity * SOL))
        ids.append(cur.lastrowid)
    return system_id, ids


def add_bright(conn, point, luminosity, system_id=None, velocity=None, cluster_id=None, address=(1, 0, 0), temperature=9000.0):
    """An unbuilt bright star at `point` (parsecs) of `luminosity` solar; returns its id."""
    values = dict(luminosity_w=luminosity * SOL, temperature_k=temperature)
    row = bright_star_rows.star_row(address, tuple(1000.0 * value for value in point), **values)
    store.insert_bright_stars(conn, [row])
    star_id = conn.execute("SELECT MAX(id) AS i FROM bright_stars").fetchone()["i"]
    if system_id is not None:
        conn.execute("UPDATE bright_stars SET star_system_id = ? WHERE id = ?", (system_id, star_id))
    if velocity is not None:
        conn.execute("UPDATE bright_stars SET velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ?, epoch_unix = 5.0e8 WHERE id = ?", (*velocity, star_id))
    if cluster_id is not None:
        conn.execute("UPDATE bright_stars SET cluster_id = ? WHERE id = ?", (cluster_id, star_id))
    return star_id


def add_scatter(conn, point, luminosity, kind="black-hole", built=False, velocity=(None, None, None)):
    """A scattered phenomenon at `point` (parsecs); returns its id."""
    cur = conn.execute("INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, position_x_mpc, position_y_mpc, position_z_mpc, velocity_x_kms,"
                 " velocity_y_kms, velocity_z_kms, epoch_unix, seed, map_luminosity_sol, built_at) VALUES (1, 0, 0, ?, ?, ?, ?, ?, ?, ?, 7.0e8, 1, ?, "
                 + ("CURRENT_TIMESTAMP" if built else "NULL") + ")", (kind, *(1000.0 * value for value in point), *velocity, luminosity))
    return cur.lastrowid


def add_cluster(conn):
    return conn.execute("INSERT INTO galaxy_components (kind, name, source, center_x_pc, center_y_pc, center_z_pc, velocity_x_kms, velocity_y_kms,"
                        " velocity_z_kms, epoch_unix, mass_solar, half_mass_radius_pc, profile, king_w0, tidal_radius_pc, age_gy, metallicity_fe_h, seed)"
                        " VALUES ('cluster', 'Test cluster', 'generated', 0, 0, 0, 0, 0, 0, 0, 2e5, 5, 'king', 7, 40, 11, -1.5, 1)").lastrowid


def keys(found):
    return set(zip(found.source.tolist(), found.ident.tolist()))


# --- A sphere against a brute force over the same rows ---------------------------------------------------------------------------------

def test_a_sphere_holds_exactly_the_rows_a_brute_force_finds(conn):
    rng = random.Random(7)
    spot = lambda: tuple(round(rng.uniform(-40.0, 40.0), 3) for _ in range(3))     # bright stars and phenomena are stored in millparsecs
    expected = {}
    for index in range(60):
        point, lum = spot(), round(10 ** rng.uniform(-1, 4), 3)
        system_id, ids = add_system(conn, point, [(lum, 5800.0)], name=f"S{index}")
        expected[(CHARTED, ids[0])] = (point, lum)
    for _ in range(120):
        point, lum = spot(), round(10 ** rng.uniform(-1, 4), 3)
        expected[(BRIGHT, add_bright(conn, point, lum))] = (point, lum)
    for _ in range(40):
        point, lum = spot(), round(10 ** rng.uniform(-1, 4), 3)
        expected[(SCATTER, add_scatter(conn, point, lum))] = (point, lum)
    conn.commit()
    center, radius, floor = (3.0, -2.0, 5.0), 30.0, 50.0
    found = sky_sphere.stars_in_sphere(conn, center, radius, floor)
    want = {key: value for key, value in expected.items() if math.dist(center, value[0]) <= radius and value[1] * SOL >= floor * SOL}
    assert keys(found) == set(want) and len(want) > 8
    assert len(found) == len(want) and found.count(CHARTED) + found.count(BRIGHT) + found.count(SCATTER) == len(want)
    for index in range(len(found)):
        point, lum = want[(int(found.source[index]), int(found.ident[index]))]
        assert found.distance_pc[index] == pytest.approx(math.dist(center, point), rel=1e-12)
        assert found.positions_pc()[index] == pytest.approx(point)
        assert found.luminosity_sol[index] == pytest.approx(lum, rel=1e-12)
    assert list(found.distance_pc) == sorted(found.distance_pc)


def test_the_rows_come_nearest_first_then_by_source_and_id_and_the_same_every_time(conn):
    for source in (SCATTER, BRIGHT, CHARTED):
        {CHARTED: lambda: add_system(conn, (5.0, 0.0, 0.0), [(100.0, 5000.0)]), BRIGHT: lambda: add_bright(conn, (5.0, 0.0, 0.0), 100.0),
         SCATTER: lambda: add_scatter(conn, (5.0, 0.0, 0.0), 100.0)}[source]()
    add_bright(conn, (1.0, 0.0, 0.0), 100.0)
    conn.commit()
    first = sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 10.0, 1.0)
    second = sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 10.0, 1.0)
    assert first.source.tolist() == [BRIGHT, CHARTED, BRIGHT, SCATTER] and first.distance_pc.tolist() == [1.0, 5.0, 5.0, 5.0]
    for name in sky_sphere._FIELDS:
        np.testing.assert_array_equal(getattr(first, name), getattr(second, name))


# --- The floor, the radius and the places a row can be ------------------------------------------------------------------------------------

def test_a_star_at_exactly_the_floor_is_in_and_one_a_little_below_is_out(conn):
    _system, at = add_system(conn, (1.0, 0.0, 0.0), [(100.0, 5000.0)])
    _system, below = add_system(conn, (2.0, 0.0, 0.0), [(99.99, 5000.0)])
    bright_at, bright_below = add_bright(conn, (3.0, 0.0, 0.0), 100.0), add_bright(conn, (4.0, 0.0, 0.0), 99.99)
    scatter_at, scatter_below = add_scatter(conn, (5.0, 0.0, 0.0), 100.0), add_scatter(conn, (6.0, 0.0, 0.0), 99.99)
    conn.commit()
    found = sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 10.0, 100.0)
    assert keys(found) == {(CHARTED, at[0]), (BRIGHT, bright_at), (SCATTER, scatter_at)}


def test_only_the_stars_above_the_floor_of_a_system_are_returned(conn):
    _system, ids = add_system(conn, (1.0, 1.0, 1.0), [(500.0, 9000.0), (0.5, 3000.0), (60.0, 6000.0)])
    conn.commit()
    found = sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 5.0, 50.0)
    assert sorted(found.ident.tolist()) == sorted([ids[0], ids[2]])
    assert sorted(found.luminosity_sol.tolist()) == pytest.approx([60.0, 500.0])
    assert set(found.temperature_k.tolist()) == {9000.0, 6000.0} and set(found.system_id.tolist()) == {found.system_id[0]}
    assert sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 5.0, 1000.0).source.size == 0


def test_a_row_on_the_surface_is_in_and_one_past_it_is_out_and_a_zero_radius_holds_only_the_center(conn):
    on = add_bright(conn, (3.0, 4.0, 0.0), 10.0)
    add_bright(conn, (3.0, 4.0, 0.001), 10.0)
    center = add_bright(conn, (0.0, 0.0, 0.0), 10.0)
    conn.commit()
    assert keys(sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 5.0, 1.0)) == {(BRIGHT, on), (BRIGHT, center)}
    assert keys(sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 0.0, 1.0)) == {(BRIGHT, center)}


def test_what_is_built_is_read_where_it_was_built_and_never_twice(conn):
    system_id, ids = add_system(conn, (2.0, 0.0, 0.0), [(300.0, 6000.0)])
    unbuilt = add_bright(conn, (1.0, 0.0, 0.0), 300.0)
    add_bright(conn, (2.0, 0.0, 0.0), 300.0, system_id=system_id)     # the bright star that became the charted system's star
    unlit = add_scatter(conn, (3.0, 0.0, 0.0), 300.0)
    add_scatter(conn, (3.5, 0.0, 0.0), 300.0, built=True)
    conn.execute("INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, position_x_mpc, position_y_mpc, position_z_mpc, seed)"
                 " VALUES (1, 0, 0, 'quasar', 4000, 0, 0, 1)")        # no map luminosity: nothing to see
    conn.commit()
    found = sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 10.0, 1.0)
    assert keys(found) == {(CHARTED, ids[0]), (BRIGHT, unbuilt), (SCATTER, unlit)}


def test_a_row_beyond_the_box_of_cell_keys_is_found_by_the_branch_for_rows_with_no_key(conn):
    far = (300000.0, 10.0, -5.0)
    inside = add_bright(conn, far, 40.0)
    add_bright(conn, (-300000.0, 0.0, 0.0), 40.0)
    add_scatter(conn, (300000.0, 10.0, 200.0), 40.0)
    conn.commit()
    assert conn.execute("SELECT cell_key FROM bright_stars WHERE id = ?", (inside,)).fetchone()["cell_key"] is None
    found = sky_sphere.stars_in_sphere(conn, far, 20.0, 1.0)
    assert keys(found) == {(BRIGHT, inside)}


def test_an_empty_galaxy_gives_empty_arrays_of_the_right_types_and_no_cells(conn):
    found = sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 100.0, 1.0)
    assert len(found) == 0 and found.positions_pc().shape == (0, 3) and found.velocities_kms().shape == (0, 3)
    assert found.source.dtype == np.int8 and found.ident.dtype == np.int64 and found.x_pc.dtype == np.float64 and found.kind.dtype == object
    assert found.levels.cells == () and found.levels.complete_above_sol == math.inf and found.levels.untouched == 0


# --- What each row carries ----------------------------------------------------------------------------------------------------------------

def test_the_velocity_epoch_age_cluster_and_kind_of_each_source(conn):
    cluster = add_cluster(conn)
    _system, ids = add_system(conn, (1.0, 0.0, 0.0), [(100.0, 5000.0)], velocity=(10.0, -20.0, 30.0), epoch=2.5e9)
    placed = add_bright(conn, (2.0, 0.0, 0.0), 100.0)
    moving = add_bright(conn, (3.0, 0.0, 0.0), 100.0, velocity=(4.0, 5.0, 6.0), cluster_id=cluster)
    runaway = add_scatter(conn, (4.0, 0.0, 0.0), 100.0, kind="hypervelocity-star", velocity=(900.0, 0.0, 0.0))
    still = add_scatter(conn, (5.0, 0.0, 0.0), 100.0, kind="neutron-star")
    conn.commit()
    found = sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 10.0, 1.0)
    row = {(int(s), int(i)): n for n, (s, i) in enumerate(zip(found.source, found.ident))}
    charted, a, b, c, d = (row[(CHARTED, ids[0])], row[(BRIGHT, placed)], row[(BRIGHT, moving)], row[(SCATTER, runaway)], row[(SCATTER, still)])
    assert found.velocities_kms()[charted].tolist() == [10.0, -20.0, 30.0] and found.epoch_unix[charted] == 2.5e9
    assert found.age_gy[charted] == 4.0 and found.lifespan_gy[charted] == 10.0 and found.cluster_id[charted] == 0 and found.kind[charted] == ""
    assert np.isnan(found.velocities_kms()[a]).all() and np.isnan(found.epoch_unix[a])      # stored NULL until the system is built
    assert found.velocities_kms()[b].tolist() == [4.0, 5.0, 6.0] and found.epoch_unix[b] == 5.0e8 and found.cluster_id[b] == cluster
    assert found.velocities_kms()[c].tolist() == [900.0, 0.0, 0.0] and found.kind[c] == "hypervelocity-star" and found.epoch_unix[c] == 7.0e8
    assert np.isnan(found.velocities_kms()[d]).all() and found.kind[d] == "neutron-star"
    assert np.isnan(found.temperature_k[c]) and np.isnan(found.age_gy[c]) and found.system_id[charted] > 0 and found.system_id[a] == 0
    assert found.temperature_k[a] == 9000.0 and found.age_gy[a] == bright_star_rows.DEFAULTS["age_gy"]


def test_a_star_that_never_ends_has_an_infinite_lifespan(conn):
    star = add_bright(conn, (1.0, 0.0, 0.0), 100.0)
    conn.execute("UPDATE bright_stars SET lifespan_gy = NULL WHERE id = ?", (star,))
    conn.commit()
    assert sky_sphere.stars_in_sphere(conn, (0.0, 0.0, 0.0), 5.0, 1.0).lifespan_gy.tolist() == [math.inf]


# --- Bad asks ------------------------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("center, radius, floor", [
    ((0.0, 0.0), 1.0, 1.0), ((0.0, 0.0, math.nan), 1.0, 1.0), ((0.0, 0.0, 0.0), -1.0, 1.0), ((0.0, 0.0, 0.0), math.inf, 1.0),
    ((0.0, 0.0, 0.0), 1.0, -1.0), ((0.0, 0.0, 0.0), 1.0, math.nan), ((0.0, "x", 0.0), 1.0, 1.0), ((0.0, 0.0, 0.0), None, 1.0)])
def test_a_center_a_radius_or_a_floor_that_is_not_a_finite_number_is_refused(conn, center, radius, floor):
    with pytest.raises(SkySphereError):
        sky_sphere.stars_in_sphere(conn, center, radius, floor)


# --- The cell levels --------------------------------------------------------------------------------------------------------------------------

def test_the_cell_levels_name_each_sector_the_sphere_reaches_and_how_far_each_was_drawn(mysql_config):
    shape = toy_shape(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0, bulge_amplitude=2.0, arm_count=2,
                      pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=2000.0, config=mysql_config,
                            cell_star_ceiling=1e12)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE galaxy_shape SET bright_star_min_luminosity_sol = 800.0, bright_star_seed = 11 WHERE id = 1")
        center = sector_position_pc(3, 0, 4, EDGE_PC)
        probe = sky_sphere.cell_levels(conn, center, 2.0 * EDGE_PC)
        assert len(probe.cells) > 5 and (3, 0, 4) in [cell.address for cell in probe.cells]
        assert probe.scatter_level_sol == 800.0 and probe.untouched == len(probe.cells) and probe.filled == 0
        assert probe.complete_above_sol == 800.0 and all(cell.complete_above_sol == 800.0 for cell in probe.cells)
        addresses = [cell.address for cell in probe.cells]
        drawn, deeper, filled = addresses[0], addresses[1], addresses[2]
        for address, level, mass in ((drawn, 300.0, 2.0), (deeper, 40.0, None), (filled, 0.0, 1.0)):
            conn.execute("INSERT INTO sector_stats (shell_index, band_index, wedge_index, bright_level_sol, bright_mass_sol) VALUES (?, ?, ?, ?, ?)",
                         (*address, level, mass))
        conn.commit()
        found = sky_sphere.stars_in_sphere(conn, center, 2.0 * EDGE_PC, 1.0)
        levels = {cell.address: cell for cell in found.levels.cells}
        assert sorted(levels) == sorted(addresses)
        assert (levels[drawn].level_sol, levels[drawn].mass_sol, levels[drawn].complete_above_sol) == (300.0, 2.0, 300.0)
        assert (levels[deeper].level_sol, levels[deeper].mass_sol, levels[deeper].complete_above_sol) == (40.0, None, 40.0)
        assert (levels[filled].level_sol, levels[filled].mass_sol, levels[filled].complete_above_sol) == (0.0, 1.0, 0.0)
        assert found.levels.untouched == len(addresses) - 3 and found.levels.filled == 1
        assert found.levels.complete_above_sol == 800.0           # an untouched sector is complete only to the scatter's level
        for address in addresses[3:]:
            conn.execute("INSERT INTO sector_stats (shell_index, band_index, wedge_index, bright_level_sol) VALUES (?, ?, ?, 25.0)", address)
        conn.commit()
        assert sky_sphere.cell_levels(conn, center, 2.0 * EDGE_PC).complete_above_sol == 300.0      # now the shallowest drawn sector is the limit
    finally:
        conn.close()


def test_every_sector_that_holds_a_returned_star_is_among_the_cells_reached(mysql_config):
    shape = toy_shape(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0, bulge_amplitude=2.0, arm_count=2,
                      pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=2000.0, config=mysql_config,
                            cell_star_ceiling=1e12)
    conn = store.get_connection(mysql_config)
    try:
        rng = random.Random(3)
        center = sector_position_pc(3, 0, 4, EDGE_PC)
        for _ in range(80):
            point = tuple(c + rng.uniform(-1.5 * EDGE_PC, 1.5 * EDGE_PC) for c in center)
            add_bright(conn, point, 100.0, address=spherical.sector_address_at(point, EDGE_PC))
        conn.commit()
        found = sky_sphere.stars_in_sphere(conn, center, 1.2 * EDGE_PC, 1.0)
        reached = {cell.address for cell in found.levels.cells}
        held = {tuple(row.values()) for row in conn.execute(
            "SELECT shell_index, band_index, wedge_index FROM bright_stars WHERE id IN (%s)" % ", ".join(str(i) for i in found.ident.tolist())).fetchall()}
        assert found.source.size > 10 and held <= reached
    finally:
        conn.close()


# --- Cost ----------------------------------------------------------------------------------------------------------------------------------

def test_a_sphere_of_a_crowded_galaxy_is_read_by_the_cell_index_and_the_time_is_measured(conn):
    """30,000 bright stars, 6,000 charted systems and 3,000 phenomena in a 400 pc box: a 60 pc sphere costs a fraction of a second, and the answer is a brute force's."""
    rng = random.Random(11)
    rows, points = [], []
    for _ in range(30000):
        point = tuple(rng.uniform(-200.0, 200.0) for _ in range(3))
        lum = 10 ** rng.uniform(0, 4)
        points.append((BRIGHT, point, lum))
        rows.append(bright_star_rows.star_row((1, 0, 0), tuple(1000.0 * v for v in point), luminosity_w=lum * SOL))
    store.insert_bright_stars(conn, rows)
    config_id = conn.execute("INSERT INTO system_configs (name) VALUES (NULL)").lastrowid
    system_rows = []
    for index in range(6000):
        point = tuple(rng.uniform(-200.0, 200.0) for _ in range(3))
        system_rows.append((config_id, f"S{index}", *point))
    conn.executemany("INSERT INTO star_systems (system_config_id, name, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc) VALUES (?, ?, ?, ?, ?)", system_rows)
    ids = [row["id"] for row in conn.execute("SELECT id FROM star_systems ORDER BY id").fetchall()]
    star_rows = []
    for system_id, (_c, _n, *point) in zip(ids, system_rows):
        lum = 10 ** rng.uniform(0, 4)
        points.append((CHARTED, tuple(point), lum))
        star_rows.append((system_id, "x", lum * SOL))
    conn.executemany(
        "INSERT INTO stars (star_system_id, role, name, star_type, yerkes_class, mass_kg, radius_km, temperature_k, luminosity_w, age_gy, habitable_zone_inner_km,"
        " habitable_zone_outer_km, system_perimeter_km, heliosphere_radius_km, galactic_orbital_period_gy, galactic_orbital_phase_deg,"
        " galactic_min_update_interval_years, snow_line_km) VALUES (?, 'single', ?, 'G2V', 'V', 2e30, 7e5, 5800, ?, 4.0, 1, 2, 3, 4, 0.2, 0, 1, 1)", star_rows)
    for _ in range(3000):
        point = tuple(rng.uniform(-200.0, 200.0) for _ in range(3))
        lum = 10 ** rng.uniform(0, 4)
        points.append((SCATTER, point, lum))
        add_scatter(conn, point, lum)
    conn.commit()
    center, radius, floor = (10.0, 20.0, -30.0), 60.0, 300.0
    started = time.perf_counter()
    found = sky_sphere.stars_in_sphere(conn, center, radius, floor)
    seconds = time.perf_counter() - started
    want = [point for source, point, lum in points if math.dist(center, point) <= radius and lum * SOL >= floor * SOL]
    assert len(found) == len(want) > 100
    print(f"\nsky sphere: {len(found)} rows of {len(points)} in a 60 pc sphere, {seconds * 1000:.0f} ms")
    assert seconds < 5.0
