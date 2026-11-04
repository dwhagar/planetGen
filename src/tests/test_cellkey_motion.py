# tests/test_cellkey_motion.py

"""
The cell key follows an object across a cube boundary, and the point columns and the key round-trip bit for bit (TEST.166).

`test_cellkey.py` proves MySQL's key equals the Python key on a scratch table. This file proves it on the real tables and
through the real movers: the orbit update steps an object a few hundredths of a parsec over the face of a 0.25 pc cube and
the stored key is the new cube's; a step that stays inside the cube leaves the key alone; and a point written to any
positioned table reads back as the same bits with the same key.
"""

import math
import re
import struct

import pytest

from planetgen.db import store
from planetgen.galaxy import cellkey
from planetgen.galaxy.frame import galaxy_to_local_pc
from planetgen.galaxy.galactic_orbit import circular_orbit
from planetgen.galaxy.sector import SpaceSector
from planetgen.galaxy.spherical import sector_address_at, sector_position_pc
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants as physical_constants
from planetgen.physics.units import pc_to_ly
from tests.toy_shape import toy_shape
from tests.rich_galaxy_support import rich_galaxy  # noqa: F401  (the fixture)

EDGE_PC = 4.0
SHELL = 10
STEP_PC = 0.04
"""float: How far a crossing step moves an object: more than the 0.01 pc it starts short of the face, far less than a cube."""
MARGIN_PC = 0.01
"""float: How far short of the face a crossing object starts."""
STAY_PC = 0.004
"""float: A step that must stay inside its cube (it starts at the cube's center)."""
RADIUS_PC = 8000.0


def _cube(point):
    return tuple(cellkey.cell_index(c) for c in point)


def _centered(point):
    """The center of the cube holding `point`."""
    return tuple(cellkey.cell_index(c) * cellkey.CELL_PC - cellkey.REACH_PC + cellkey.CELL_PC / 2 for c in point)


def _tangent(point):
    """The unit vector an orbit moves a point along (counterclockwise about the galaxy's axis)."""
    radius = math.hypot(point[0], point[1])
    return (-point[1] / radius, point[0] / radius, 0.0)


def _just_short_of_a_face(near, direction):
    """
    A point in the cube of `near`, `MARGIN_PC` short of the face a step along `direction` reaches first (the axis the
    direction leans on most), at the cube's center on the other axes; and that axis.
    """
    start = list(_centered(near))
    axis = max(range(3), key=lambda a: abs(direction[a]))
    cube = cellkey.cell_index(near[axis])
    if direction[axis] > 0:
        start[axis] = (cube + 1) * cellkey.CELL_PC - cellkey.REACH_PC - MARGIN_PC
    else:
        start[axis] = cube * cellkey.CELL_PC - cellkey.REACH_PC + MARGIN_PC
    return tuple(start), axis


def _years_for(distance_pc, radius_pc, period_gy):
    """The years an orbit of `period_gy` takes to carry a point at `radius_pc` along `distance_pc` of its circle."""
    return distance_pc / radius_pc * period_gy * 1e9 / (2 * math.pi)


MOVER_SPEED_KMS = 100.0
"""float: The speed the rows of a test that sets its own velocity are given along the tangent: the stepper carries a row by the velocity it stores (GEN.340), not by
a period written beside it."""


def _period_gy_at_speed(radius_pc):
    """The period, Gy, of a circle of `radius_pc` at `MOVER_SPEED_KMS`, so `_years_for` gives the years that speed takes to cover the distance."""
    return 2 * math.pi * radius_pc * physical_constants.PARSEC_M / (MOVER_SPEED_KMS * 1000.0) / physical_constants.SECONDS_PER_YEAR / 1e9


def _stored_galaxy_period_gy(mysql_config):
    """`radius_pc -> period_gy` of a circle in the potential of the galaxy the database stores, which the stepper and the circular velocity of a new point go by."""
    conn = store.get_connection(mysql_config)
    try:
        pot = store._galactic_potential(conn)
    finally:
        conn.close()
    return lambda radius: circular_orbit(pot, pc_to_ly(radius))[1]


def _tangent_velocity(point):
    return tuple(MOVER_SPEED_KMS * c for c in _tangent(point))


def _bits(value):
    return None if value is None else struct.pack("<d", value)


def _assert_followed(before_point, before_key, after_point, after_key, axis, sign, what):
    """The key is the Python key of the stored point, before and after, and it moved by exactly one cube on `axis`."""
    assert before_key == cellkey.cell_key(*before_point), (what, "before")
    assert after_key == cellkey.cell_key(*after_point), (what, "after")
    assert after_key != before_key, what
    moved = tuple(b - a for a, b in zip(_cube(before_point), _cube(after_point)))
    expected = tuple(sign if a == axis else 0 for a in range(3))
    assert moved == expected, (what, moved)


def _assert_stayed(before_point, before_key, after_point, after_key, what):
    assert before_key == cellkey.cell_key(*before_point), (what, "before")
    assert _cube(after_point) == _cube(before_point), what
    assert after_key == before_key == cellkey.cell_key(*after_point), what


def _place_at(points, name):
    """A sector holding one system at each galaxy-frame point of `points` (wedge 3's sector, so the points are near its center)."""
    center = sector_position_pc(SHELL, 0, 3, EDGE_PC)
    sector = SpaceSector(name, edge_ly=pc_to_ly(EDGE_PC))
    for point in points:
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = False
        cfg.BINARY_SYSTEM = False
        offset_ly = tuple(pc_to_ly(c) for c in galaxy_to_local_pc(center, point))
        sector.add_system(StarSystem(system_config=cfg), position=offset_ly, system_config=cfg)
    return sector, {"center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
                    "galactic_radius_pc": math.hypot(*center), "shell_index": SHELL, "band_index": 0, "wedge_index": 3}


def _seed_scatter_grid(mysql_config):
    shape = toy_shape(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
                               bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=2000.0,
                            config=mysql_config)


def _insert_scatter_point(conn, table, point_pc, kind=None, velocity=None):
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
        " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1, NULL)",
        (shell, band, wedge, kind, x, y, z, *v)).lastrowid


def _scatter_point(row):
    return tuple(row[f"position_{a}_mpc"] / 1000.0 for a in "xyz")


def _moved_by_orbit(mysql_config, make_row, read_row, period_gy, near, what):
    """
    Steps one orbiting row twice, once over a cube face and once inside a cube.

    `make_row(conn, point)` stores a row at `point` and returns its id; `read_row(conn, id)` returns `(point, key)`.
    """
    direction = _tangent(near)
    cross_start, axis = _just_short_of_a_face(near, direction)
    stay_start = _centered(near)
    conn = store.get_connection(mysql_config)
    try:
        for start, step, crossing in ((cross_start, STEP_PC, True), (stay_start, STAY_PC, False)):
            with conn:
                row_id = make_row(conn, start)
                before = read_row(conn, row_id)
                radius = math.hypot(start[0], start[1])
                store.advance_galactic_positions(conn, store.orbit_clock(conn, _years_for(step, radius, period_gy(radius))))
                after = read_row(conn, row_id)
                conn.execute("UPDATE orbit_simulation_state SET last_updated_at = NULL")  # each step starts from one clock
            moved = math.dist(before[0], after[0])
            assert moved == pytest.approx(step, rel=0.01), (what, "the step the test asked for")
            if crossing:
                _assert_followed(*before, *after, axis, 1 if direction[axis] > 0 else -1, what)
            else:
                _assert_stayed(*before, *after, what)
    finally:
        conn.close()


def test_a_system_that_steps_over_a_cube_face_gets_the_next_cubes_key(mysql_config):
    """`star_systems`: the sector center plus the offset is the point, and the key is the Python key of it, before and after."""
    center = sector_position_pc(SHELL, 0, 3, EDGE_PC)
    direction = _tangent(center)
    cross_start, axis = _just_short_of_a_face(center, direction)
    stay_start = _centered(center)
    sector, position = _place_at([cross_start, stay_start], "Keyed")
    store.save_sector(sector, config=mysql_config, galaxy_position=position)
    conn = store.get_connection(mysql_config)
    try:
        def read():
            rows = conn.execute("SELECT galaxy_x_pc, galaxy_y_pc, galaxy_z_pc, cell_key FROM star_systems"
                                " ORDER BY id").fetchall()
            return [((r["galaxy_x_pc"], r["galaxy_y_pc"], r["galaxy_z_pc"]), r["cell_key"]) for r in rows]

        (cross_before, stay_before) = read()
        assert math.dist(cross_before[0], cross_start) < 1e-6 and math.dist(stay_before[0], stay_start) < 1e-6
        with conn:
            for (_before, _key), start, row_id in zip((cross_before, stay_before), (cross_start, stay_start),
                                                       [r["id"] for r in conn.execute("SELECT id FROM star_systems ORDER BY id").fetchall()]):
                conn.execute("UPDATE star_systems SET velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ? WHERE id = ?",
                             (*_tangent_velocity(start), row_id))
            store.advance_galactic_positions(
                conn, store.orbit_clock(conn, _years_for(STEP_PC, math.hypot(*cross_start[:2]), _period_gy_at_speed(math.hypot(*cross_start[:2])))))
        cross_after, stay_after = read()
    finally:
        conn.close()
    assert math.dist(cross_before[0], cross_after[0]) == pytest.approx(STEP_PC, rel=0.01)
    _assert_followed(*cross_before, *cross_after, axis, 1 if direction[axis] > 0 else -1, "star_systems")
    # The second system started at the center of its cube and also moved STEP_PC (0.04 pc, a sixth of a cube): inside.
    assert math.dist(stay_before[0], stay_after[0]) == pytest.approx(STEP_PC, rel=0.01)
    _assert_stayed(*stay_before, *stay_after, "star_systems")


def test_a_nebula_and_a_facility_that_step_over_a_cube_face_get_the_next_cubes_key(mysql_config):
    """`nebulae` and `facilities` keep their point in `center_*_pc`; the orbit update turns it and the key follows."""
    from planetgen.generation.phenomena.nebula import Nebula
    sector, position = _place_at([sector_position_pc(SHELL, 0, 3, EDGE_PC)], "Hosts")
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position=position)
    center = (position["center_x_pc"], position["center_y_pc"], position["center_z_pc"])

    def nebula(conn, point):
        body = Nebula(SystemConfig())
        body.galactic_orbital_period_gy = 1.0
        row_id = store.insert_nebula(conn, body, sector_id=sector_id, placement={
            "center_x_pc": point[0], "center_y_pc": point[1], "center_z_pc": point[2],
            "galactic_radius_pc": math.hypot(*point)})
        conn.execute("UPDATE nebulae SET velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ? WHERE id = ?",
                     (*_tangent_velocity(point), row_id))
        return row_id

    def facility(conn, point):
        row_id = store.add_facility(conn, "Drifter", "station", "standalone", "space", sector_id)
        conn.execute("UPDATE facilities SET center_x_pc = ?, center_y_pc = ?, center_z_pc = ?, epoch_unix = NULL,"
                     " next_update_due = NULL WHERE id = ?", (*point, row_id))
        return row_id

    def read(table):
        def one(conn, row_id):
            row = conn.execute(f"SELECT center_x_pc, center_y_pc, center_z_pc, cell_key FROM {table} WHERE id = ?",
                               (row_id,)).fetchone()
            return (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]), row["cell_key"]
        return one

    _moved_by_orbit(mysql_config, nebula, read("nebulae"), _period_gy_at_speed, center, "nebulae")
    _moved_by_orbit(mysql_config, facility, read("facilities"), _stored_galaxy_period_gy(mysql_config), center, "facilities")


def test_scatter_points_that_step_over_a_cube_face_get_the_next_cubes_key(mysql_config):
    """`bright_stars` and the orbiting `phenomenon_scatter` points turn about the galaxy; the key follows the mpc columns."""
    _seed_scatter_grid(mysql_config)
    near = (RADIUS_PC, 0.1, 0.1)

    period = _stored_galaxy_period_gy(mysql_config)

    def bright_star(conn, point):
        return _insert_scatter_point(conn, "bright_stars", point)

    def black_hole(conn, point):
        return _insert_scatter_point(conn, "phenomenon_scatter", point, "black-hole")

    def read(table):
        def one(conn, row_id):
            row = conn.execute(f"SELECT * FROM {table} WHERE id = ?", (row_id,)).fetchone()
            return _scatter_point(row), row["cell_key"]
        return one

    _moved_by_orbit(mysql_config, bright_star, read("bright_stars"), period, near, "bright_stars")
    _moved_by_orbit(mysql_config, black_hole, read("phenomenon_scatter"), period, near, "phenomenon_scatter")


def test_a_hypervelocity_star_that_flies_over_a_cube_face_gets_the_next_cubes_key(mysql_config):
    """The straight-line mover: the same point columns, a different formula."""
    _seed_scatter_grid(mysql_config)
    elapsed_years = 1e5
    speed_kms = STEP_PC * physical_constants.PARSEC_M / 1000.0 / (elapsed_years * physical_constants.SECONDS_PER_YEAR)
    for direction, axis, sign in (((1.0, 0.0, 0.0), 0, 1), ((0.0, -1.0, 0.0), 1, -1), ((0.0, 0.0, 1.0), 2, 1)):
        start, _axis = _just_short_of_a_face((RADIUS_PC, 0.1, 0.1), direction)
        velocity = tuple(speed_kms * d for d in direction)
        conn = store.get_connection(mysql_config)
        try:
            with conn:
                row_id = _insert_scatter_point(conn, "phenomenon_scatter", start, "hypervelocity-star", velocity)
                before_row = conn.execute("SELECT * FROM phenomenon_scatter WHERE id = ?", (row_id,)).fetchone()
                store.advance_scatter_points(conn, store.orbit_clock(conn, elapsed_years))
                after_row = conn.execute("SELECT * FROM phenomenon_scatter WHERE id = ?", (row_id,)).fetchone()
        finally:
            conn.close()
        assert math.dist(_scatter_point(before_row), _scatter_point(after_row)) == pytest.approx(STEP_PC, rel=0.01)
        _assert_followed(_scatter_point(before_row), before_row["cell_key"], _scatter_point(after_row),
                         after_row["cell_key"], axis, sign, f"hypervelocity {direction}")


# ---------------------------------------------------------------------------------------------------------------
# Every table that holds a position
# ---------------------------------------------------------------------------------------------------------------

NOT_IN_THE_RICH_GALAXY = {
    "nebula_field": lambda conn: conn.execute(
        "INSERT INTO nebula_field (cell_i, cell_j, cell_k, n, nebula_class, radius_ly, center_x_pc, center_y_pc,"
        " center_z_pc) VALUES (0, 0, 0, 0, 'M', 10.0, 1.0, 2.0, 3.0)"),
    "phenomenon_scatter": lambda conn: _insert_scatter_point(conn, "phenomenon_scatter", (8000.0, 0.0, 0.0), "black-hole"),
}
"""dict: Tables with a cell key and no row in `rich_galaxy` (the whole-sphere fill and the scatter are not run there), each
with the statement that gives it one for the test; a table that gains a row in the rich galaxy must leave this list."""


def _positioned_tables(conn):
    """`{table: (the three point columns, divisor)}` for every table with a generated `cell_key`, read from its definition."""
    rows = conn.execute(
        "SELECT TABLE_NAME AS t, GENERATION_EXPRESSION AS e FROM information_schema.COLUMNS"
        " WHERE TABLE_SCHEMA = DATABASE() AND COLUMN_NAME = 'cell_key'").fetchall()
    tables = {}
    for row in rows:
        columns = re.findall(r"abs\(+`(\w+)`(?: / 1000\.0)?\)+ < 262144", row["e"], flags=re.I)
        assert len(columns) == 3, (row["t"], row["e"][:200])
        tables[row["t"]] = (columns, 1000.0 if all(c.endswith("_mpc") for c in columns) else 1.0)
    return tables


def _awkward_points():
    """Points on and around cube faces, with the doubles that do not survive a decimal round trip."""
    cell = cellkey.CELL_PC
    face = 8000.25
    below, above = math.nextafter(face, 0.0), math.nextafter(face, math.inf)
    points = [
        (8000.1 + 0.1 + 0.2, 0.1 + 0.2, -0.1 - 0.2),
        (face, 0.0, 0.0), (below, cell, -cell), (above, -cell, cell),
        (math.nextafter(cell, 0.0), math.nextafter(-cell, 0.0), 1.0 / 3.0),
        (-8000.25, math.nextafter(-8000.25, -math.inf), 8000.25),
        (1e-300, -1e-300, 5e-324),
        (cellkey.REACH_PC - cell, 123.456789012345678, -cellkey.REACH_PC + cell),
        (math.nextafter(cellkey.REACH_PC, 0.0), 0.0, 0.0), (cellkey.REACH_PC, 0.0, 0.0),  # the last cube; off the box
        (1e7, 0.0, 0.0),
    ]
    return points


def test_every_positioned_table_round_trips_its_point_and_key_bit_for_bit(rich_galaxy):  # noqa: F811
    conn = store.get_connection(rich_galaxy)
    try:
        tables = _positioned_tables(conn)
        assert len(tables) >= 14, sorted(tables)  # the systems, the phenomena, the facilities, the scatters, the components
        empty = sorted(t for t in tables if conn.execute(f"SELECT 1 FROM `{t}` LIMIT 1").fetchone() is None)
        assert empty == sorted(NOT_IN_THE_RICH_GALAXY), ("tables with no row, or rows that now exist", empty)
        points = _awkward_points()
        for table, (columns, divisor) in sorted(tables.items()):
            try:
                if table in NOT_IN_THE_RICH_GALAXY:
                    NOT_IN_THE_RICH_GALAXY[table](conn)
                key = "(cell_i, cell_j, cell_k, n)" if table == "nebula_field" else "id"
                where = conn.execute(f"SELECT * FROM `{table}` LIMIT 1").fetchone()
                row_id = tuple(where[c] for c in ("cell_i", "cell_j", "cell_k", "n")) if table == "nebula_field" \
                    else (where["id"],)
                marks = "(?, ?, ?, ?)" if table == "nebula_field" else "?"
                for point in points:
                    stored = tuple(c * divisor for c in point)
                    # A placed row also carries its radius from the galaxy's axis (a check ties the four together).
                    radius = ["galactic_radius_pc"] if "galactic_radius_pc" in where else []
                    conn.execute(f"UPDATE `{table}` SET " + ", ".join(f"`{c}` = ?" for c in [*columns, *radius])
                                 + f" WHERE {key} = {marks}",
                                 (*stored, *([math.hypot(stored[0], stored[1])] if radius else []), *row_id))
                    row = conn.execute("SELECT " + ", ".join(f"`{c}`" for c in columns) + f", cell_key FROM `{table}`"
                                       f" WHERE {key} = {marks}", row_id).fetchone()
                    assert [_bits(row[c]) for c in columns] == [_bits(v) for v in stored], (table, point)
                    read_back = tuple(row[c] / divisor for c in columns)
                    assert row["cell_key"] == cellkey.cell_key(*read_back), (table, point)
            finally:
                conn.rollback()
    finally:
        conn.close()
