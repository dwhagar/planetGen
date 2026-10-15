# tests/test_galaxy_core_ball.py

"""
The core ball's systems are built at a real distance from the galactic
center. The ball's center is the center itself, and a distance of zero gave
every star there a zero-size Hill sphere, so a batch fill of shell 0 made
systems without a single planet.
"""

import pytest

from planetgen import tuning
from planetgen.generation.run_galaxy import _systems_radius_pc
from planetgen.physics.units import ly_to_pc
from tests.bughunt_support import forced_system_config, mysql_argv, run_cli

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)


def test_core_ball_is_built_off_the_center():
    radius = _systems_radius_pc((0, 0, 0), 0.0, EDGE_PC)
    assert 0.0 < radius < EDGE_PC / 2


def test_other_shells_keep_their_centers_distance():
    assert _systems_radius_pc((2, 0, 0), 2 * EDGE_PC, EDGE_PC) == 2 * EDGE_PC


@pytest.mark.db
def test_core_ball_batch_fill_makes_planets(mysql_config):
    target = mysql_argv(mysql_config)
    run_cli("plan", ["--quiet", "--no-bright-stars"] + target)
    with forced_system_config(PLANETS=True):
        run_cli("galaxy", ["--shell", "0", "--band", "0", "--num-systems", "3", "--yes", "--quiet"] + target)
    from planetgen.db import store
    conn = store.get_connection(mysql_config)
    try:
        systems = conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
        planets = conn.execute("SELECT COUNT(*) AS n FROM planets").fetchone()["n"]
    finally:
        conn.close()
    assert systems == 3
    assert planets >= systems
