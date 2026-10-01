# tests/test_star_type_rendering.py

"""
White dwarf and giant systems render on the system page and System Map
(the star-type study of 2026-09-30: until the star population fix, no white dwarf or K/M giant had
ever been generated, so these paths were untested), and the Sector Map's
star dots keep them visibly apart in size.
"""

import pytest

from tests.bughunt_support import mysql_argv, run_cli
from tests.test_bughunt_end_to_end import _db_get_connection, _web_client

import web  # noqa: F401 -- puts src/html/lib on sys.path
import starmap  # noqa: E402


@pytest.mark.parametrize("star_type", ["G2VII", "K2III", "M1III", "B1IA"])
def test_white_dwarf_and_giant_systems_render(mysql_config, star_type):
    name = f"Render {star_type}"
    run_cli("system", ["--star-type", star_type, "--name", name, "--num-orbits", "4"] + mysql_argv(mysql_config))
    conn = _db_get_connection(mysql_config)
    try:
        row = conn.execute(
            "SELECT ss.id, s.star_type FROM star_systems ss JOIN stars s ON s.star_system_id = ss.id "
            "WHERE ss.name = ?", (name,),
        ).fetchone()
    finally:
        conn.close()
    assert row is not None
    html = _web_client(mysql_config).get(f"/system/{row['id']}").get_data(as_text=True)
    assert name in html
    assert row["star_type"] in html
    assert 'class="sysmap-body sysmap-star"' in html


def test_star_dots_grow_on_a_log_scale_from_white_dwarfs_to_supergiants():
    sun = starmap._SUN_RADIUS_KM
    sizes = [starmap._star_dot_radius(sun * r) for r in (0.01, 0.3, 1, 10, 100, 1000)]
    assert sizes == sorted(sizes) and len(set(sizes)) == len(sizes)
    assert sizes[0] == starmap._MIN_DOT_R
    assert sizes[2] == pytest.approx(starmap._SUN_DOT_R)
    assert sizes[-1] == pytest.approx(starmap._MAX_DOT_R, abs=0.05)
    # A giant (10 solar radii) is still larger than the Sun, but every core
    # stays small: a supergiant is a point, not a big disc.
    assert sizes[3] > sizes[2]
    assert starmap._MAX_DOT_R <= 8
    assert starmap._star_dot_radius(sun * 1e5) == starmap._MAX_DOT_R


def _glow_radius(luminosity_sol, radius_sol):
    dot_r = starmap._star_dot_radius(starmap._SUN_RADIUS_KM * radius_sol)
    glow = starmap._star_glow(luminosity_sol * starmap.SOLAR_LUMINOSITY, dot_r)
    return dot_r, glow, dot_r * glow["glowScale"]


def test_star_glow_shows_luminosity_while_the_core_stays_small():
    white_dwarf = _glow_radius(0.001, 0.01)
    sun = _glow_radius(1.0, 1.0)
    giant = _glow_radius(100.0, 15.0)
    supergiant = _glow_radius(1e5, 500.0)
    # The glow's outer radius grows with luminosity, and a supergiant's is
    # several times its core while a white dwarf's barely clears it.
    radii = [white_dwarf[2], sun[2], giant[2], supergiant[2]]
    assert radii == sorted(radii) and len(set(radii)) == 4
    assert supergiant[2] >= 5 * supergiant[0]
    assert white_dwarf[1]["glowScale"] == starmap._MIN_GLOW_SCALE
    # Stronger at the core's edge for a bright star, faint for a dim one.
    assert supergiant[1]["glowStrength"] > sun[1]["glowStrength"] > white_dwarf[1]["glowStrength"]
    # A hypergiant's halo is capped.
    hyper = _glow_radius(1e7, 1500.0)
    assert hyper[2] <= starmap._MAX_GLOW_R + 1e-6
