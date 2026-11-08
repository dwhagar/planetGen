# tests/test_nav_dense_core.py
"""NAV.76: a course across the dense core of a stored galaxy goes round the objects that are made and says how many unresolved stars it crosses."""

import math

import pytest

from planetgen import tuning
from planetgen.db import obstacles, store
from planetgen.galaxy import shape_draw
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole

SEED = bytes(range(16))
START, END = (-60.0, 0.2, 0.0), (60.0, 0.2, 0.0)
FIELD_CENTER = (59.0, 0.25, 0.0)       # near the end, where the straight line (y 0.2) and the run the course takes round the hole's sphere (y 0.25 at x = 59) meet again


@pytest.fixture
def core(mysql_config):
    """A stored galaxy with its central supermassive hole and an asteroid field near the end of the course, both made; returns their ids."""
    store.save_galaxy_shape(shape_draw.default_shape(SEED), edge_pc=float(tuning.DEFAULT_SECTOR_EDGE_PC), outer_shell_index=1500,
                            expected_system_count_at_density_1=6.0, config=mysql_config, galaxy_seed=SEED, cell_star_ceiling=500.0)
    conn = store.get_connection(mysql_config)
    try:
        hole = store.insert_black_hole(
            conn, BlackHole(SystemConfig(), mass_class="supermassive", mass_solar=4.0e6),
            placement={"center_x_pc": 0.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 0.0})
        field = AsteroidField(SystemConfig())
        field.radius_ly = 0.9   # 0.28 pc: about the widest a field stands (past a light year it has spread and is retired, GEN.1107, and a course ignores it); it lies on both
        field_id = store.insert_asteroid_field(conn, field, placement={
            "center_x_pc": FIELD_CENTER[0], "center_y_pc": FIELD_CENTER[1], "center_z_pc": FIELD_CENTER[2],
            "galactic_radius_pc": math.hypot(FIELD_CENTER[0], FIELD_CENTER[1])})
        conn.commit()
        yield {"conn": conn, "hole": hole, "field": field_id}
    finally:
        conn.close()
        obstacles._UNRESOLVED_FIELDS.clear()


def test_a_course_across_the_core_goes_round_what_is_made_and_counts_what_is_not(core):
    course = obstacles.adjusted_course(core["conn"], START, END)
    assert sorted(course["avoided"]) == sorted([f"asteroid_field:{core['field']}", f"black_hole:{core['hole']}"])
    assert course["length"] > course["straight_length"] and course["warnings"] == []
    shown = course["passed_through_unresolved"]
    assert shown["cells"] > 20 and shown["crossings"] > 50.0
    assert 0.9 * course["length"] <= shown["length_pc"] <= 1.001 * course["length"]
    assert shown["unresolved_systems"] > 1.0e5


def test_going_round_the_hole_s_sphere_keeps_the_course_out_of_the_densest_stars(core):
    course = obstacles.adjusted_course(core["conn"], START, END)
    field = obstacles.unresolved_field(core["conn"])
    counted = course["passed_through_unresolved"]["crossings"]
    assert counted == pytest.approx(field.along(course["polyline"]).crossings, rel=1e-9)
    straight = field.crossings(START, END).crossings
    assert 0.5 * straight < counted < straight, "the bent course skirts the cusp the straight line runs through"


def test_a_course_past_the_core_is_not_counted(core):
    far = obstacles.unresolved_field(core["conn"]).reach_pc + 200.0
    course = obstacles.adjusted_course(core["conn"], (far, -80.0, 0.0), (far, 80.0, 0.0))
    assert "passed_through_unresolved" not in course
