# tests/test_galaxy_core_ball.py

"""
A system at the galactic center keeps a Hill sphere. The core ball's center
is the center itself, and a distance of zero gave every star there a
zero-size Hill sphere, so a batch fill of shell 0 made systems without a
single planet.
"""

import pytest

from tests.bughunt_support import forced_system_config, mysql_argv, run_cli
from tests.test_galaxy_gen import REFERENCE_SHAPE_ARGV


def test_a_system_at_the_galactic_center_still_has_a_hill_sphere():
    from planetgen.generation.binary import BinaryStarProxy
    from planetgen.physics import constants

    mass = constants.SOLAR_MASS_TO_KG
    at_center = BinaryStarProxy._calculate_system_perimeter_static(mass, 0.0)
    assert at_center > 0.0
    assert at_center == BinaryStarProxy._calculate_system_perimeter_static(mass, constants.MIN_TIDAL_DISTANCE_LY)
    assert BinaryStarProxy._calculate_system_perimeter_static(mass, 1000.0) > at_center


@pytest.mark.db
def test_core_ball_batch_fill_makes_planets(mysql_config):
    target = mysql_argv(mysql_config)
    run_cli("plan", ["--quiet", "--nuclear-amplitude", "0"] + REFERENCE_SHAPE_ARGV + target)
    with forced_system_config(PLANETS=True):
        run_cli("galaxy", ["--shell", "0", "--band", "0", "--num-systems", "3", "--yes", "--quiet"] + target)
    from planetgen.db import store
    conn = store.get_connection(mysql_config)
    try:
        systems = conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
        # TEST.254: a fourth system was seen once and not again in 109 draws; if it comes back, the rows say which.
        rows = conn.execute("SELECT id, sector_id, name FROM star_systems ORDER BY id").fetchall() if systems != 3 else []
        seeds = conn.execute("SELECT * FROM galaxy_shape").fetchone() if systems != 3 else None
        planets = conn.execute("SELECT COUNT(*) AS n FROM planets").fetchone()["n"]
    finally:
        conn.close()
    assert systems == 3, (rows, {k: seeds[k] for k in seeds.keys() if 'seed' in k} if seeds else None)
    assert planets >= systems
