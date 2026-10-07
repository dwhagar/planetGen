# tests/test_star_type_rendering.py

"""
White dwarf and giant systems render on the system page and System Map
(the star-type study of 2026-09-30: until the star population fix, no white dwarf or K/M giant had
ever been generated, so these paths were untested), and the Sector Map's
points of light keep them visibly apart in size.
"""

import pytest

from tests.bughunt_support import mysql_argv, run_cli
from tests.test_bughunt_end_to_end import _db_get_connection, _web_client

from planetgen.web.maps import starmap  # noqa: E402


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


def _light(luminosity_sol, radius_sol, temperature_k=5778.0):
    return starmap._star_light(
        luminosity_sol * starmap.SOLAR_LUMINOSITY, starmap._SUN_RADIUS_KM * radius_sol, temperature_k,
    )


def test_star_light_shows_luminosity_as_a_point_with_a_halo():
    """MAP.15: every star is a point of light a few pixels across: the
    core grows with the star's radius, the halo's width and light with
    its luminosity, and both stay within the Galaxy Map's own ranges."""
    white_dwarf = _light(0.001, 0.01, 9000.0)
    sun = _light(1.0, 1.0)
    giant = _light(100.0, 15.0, 4500.0)
    supergiant = _light(1e5, 500.0, 3600.0)
    sizes = [star["sizePx"] for star in (white_dwarf, sun, giant, supergiant)]
    assert sizes == sorted(sizes) and len(set(sizes)) == 4
    assert supergiant["sizePx"] >= 3 * supergiant["corePx"]
    # The faint end is drawn brighter (MAP.87): the base halo widened by
    # the boost's fourth root.
    boost = starmap.star_light_boost(0.001)
    assert white_dwarf["sizePx"] == pytest.approx(starmap._LIGHT_SIZE_PX[0] * boost ** 0.25, abs=0.5)

    def halo_light(star):
        return star["glow"] * star["sizePx"] ** 2

    assert halo_light(supergiant) > halo_light(sun) > halo_light(white_dwarf)
    assert supergiant["corePx"] > sun["corePx"] > white_dwarf["corePx"]
    # A hypergiant's halo and core are capped.
    hyper = _light(1e7, 1500.0)
    assert hyper["sizePx"] <= starmap._LIGHT_SIZE_PX[1]
    assert hyper["corePx"] <= starmap._LIGHT_CORE_PX[1]
    # Color follows temperature: a hot white dwarf is bluer than a cool giant.
    assert int(white_dwarf["color"][5:7], 16) > int(supergiant["color"][5:7], 16)
