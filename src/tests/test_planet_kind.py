"""The formation kind without a class (GEN.681, physics/planet_kind.py)."""

from planetgen import tuning
from planetgen.physics import atmosphere_recipes, constants, planet_kind, planet_mass
from planetgen.util import draw

SUN_KG = constants.SOLAR_MASS_TO_KG
EARTH_KG = constants.EARTH_MASS_TO_KG
SNOW_LINE_AU = 2.7


def formations(distance_au, count=400, star_kg=SUN_KG, fe_h=0.0, seed=1):
    stream = draw.Stream(seed)
    return [planet_kind.draw_formation(distance_au, star_kg, SNOW_LINE_AU, fe_h, 10.0 ** fe_h, stream) for _ in range(count)]


def test_same_stream_same_formation():
    first = planet_kind.draw_formation(5.0, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, draw.Stream(7))
    second = planet_kind.draw_formation(5.0, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, draw.Stream(7))
    assert first == second


def test_it_draws_the_mass_law_unchanged():
    law = planet_mass.planet_mass(4.0, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, draw.Stream(11))
    got = planet_kind.draw_formation(4.0, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, draw.Stream(11))
    assert (got.mass_kg, got.core_kg, got.runaway) == law


def test_kinds_are_recipe_kinds():
    for distance in (0.1, 1.0, 5.0, 30.0):
        assert {f.kind for f in formations(distance, 100)} <= set(atmosphere_recipes.BODY_KINDS)


def test_a_runaway_planet_is_a_giant_and_formed_beyond_the_snow_line():
    runaways = [f for d in (3.0, 5.0, 10.0) for f in formations(d, 600) if f.runaway]
    assert runaways
    assert all(f.kind == "giant" and f.beyond_snow_line for f in runaways)


def test_inside_the_snow_line_nothing_is_a_giant_unless_it_ran_away():
    inner = formations(0.5, 600)
    assert all(f.kind == "terrestrial" for f in inner)
    assert not any(f.beyond_snow_line for f in inner)


def test_a_heavy_core_beyond_the_snow_line_is_an_ice_giant_without_a_runaway():
    critical = tuning.GIANT_CORE_CRITICAL_MASS_EARTH * EARTH_KG
    cores = [f for f in formations(15.0, 800, fe_h=0.3) if not f.runaway and f.core_kg >= critical]
    assert cores
    assert all(f.kind == "giant" for f in cores)
    assert all(f.kind == "terrestrial" for f in formations(15.0, 800) if not f.runaway and f.core_kg < critical)


def test_origin_flag_is_the_snow_line():
    assert not planet_kind.draw_formation(SNOW_LINE_AU - 0.01, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, draw.Stream(2)).beyond_snow_line
    assert planet_kind.draw_formation(SNOW_LINE_AU, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, draw.Stream(2)).beyond_snow_line


def test_a_given_mass_names_the_kind():
    assert planet_kind.kind_of_given(mass_kg=5 * EARTH_KG) == "terrestrial"
    assert planet_kind.kind_of_given(mass_kg=25 * EARTH_KG) == "giant"
    critical = tuning.GIANT_CORE_CRITICAL_MASS_EARTH * EARTH_KG
    assert planet_kind.kind_of_given(mass_kg=critical, beyond_snow_line=True) == "giant"
    assert planet_kind.kind_of_given(mass_kg=critical, beyond_snow_line=False) == "terrestrial"


def test_a_given_radius_names_the_kind_and_the_mass_wins():
    assert planet_kind.kind_of_given(radius_km=tuning.GIANT_MIN_RADIUS_KM) == "giant"
    assert planet_kind.kind_of_given(radius_km=tuning.GIANT_MIN_RADIUS_KM - 1) == "terrestrial"
    assert planet_kind.kind_of_given(mass_kg=EARTH_KG, radius_km=60000.0) == "terrestrial"
    assert planet_kind.kind_of_given() == "terrestrial"


def test_a_moon_ceiling_is_a_share_of_its_parent():
    assert planet_kind.moon_mass_ceiling_kg(1e27, True) == tuning.MOON_GIANT_PARENT_MAX_MASS_RATIO * 1e27
    assert planet_kind.moon_mass_ceiling_kg(6e24, False) == tuning.MOON_TERRESTRIAL_PARENT_MAX_MASS_RATIO * 6e24


def test_moon_masses_stay_between_the_floor_and_the_ceiling_and_favour_the_light():
    with draw.bound("moon-masses"):
        masses = [planet_kind.draw_moon_mass_kg(1e27, True) for _ in range(600)]
    ceiling = planet_kind.moon_mass_ceiling_kg(1e27, True)
    assert all(tuning.MOON_MIN_MASS_KG <= m <= ceiling for m in masses)
    middle = (tuning.MOON_MIN_MASS_KG * ceiling) ** 0.5
    assert sum(m < middle for m in masses) > 300


def test_a_limit_lowers_the_ceiling_and_a_tiny_parent_holds_no_moon():
    with draw.bound("moon-limit"):
        limit = 10 * tuning.MOON_MIN_MASS_KG
        assert all(planet_kind.draw_moon_mass_kg(1e27, True, limit) <= limit for _ in range(100))
        assert planet_kind.draw_moon_mass_kg(1e20, False) is None
