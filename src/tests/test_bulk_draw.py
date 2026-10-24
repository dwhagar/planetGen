"""
GEN.428: the draw of a planet's or moon's bulk composition and the mass solved from it (`planetgen.physics.bulk_draw`).
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import bulk_composition as bc, bulk_draw, bulk_structure, constants, hydrosphere
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

EARTH_RATIOS = bulk_draw.HostRatios(fe=0.81, si=0.91, al=0.07, ca=0.06, ni=0.045, c_over_o=0.55, y_envelope=0.27, z=0.014)


def _star(**extra):
    fields = dict(metallicity_fe_h=0.0, age=4.5, mass=2e30, galactic_center_dist_ly=26000.0, luminosity=3.8e26)
    fields.update(extra)
    return SimpleNamespace(**fields)


def _body(planet_class, radius_km, **extra):
    fields = dict(planet_class=planet_class, radius=radius_km, mass=None, star=_star(), is_moon=False, distance=1.0)
    fields.update(extra)
    return SimpleNamespace(**fields)


def test_an_earth_like_star_makes_an_earth_like_core():
    part = bulk_draw.dry_partition(EARTH_RATIOS, 0.081, 0.07, carbon_rock=False)
    assert part.metal + part.rock == pytest.approx(1.0)
    assert 0.27 < part.metal < 0.36
    assert part.rock_feo == pytest.approx(0.081)
    assert part.rock_mg_si == pytest.approx(1.0 / 0.91)


def test_an_iron_rich_star_makes_a_bigger_core():
    rich = bulk_draw.HostRatios(**{**EARTH_RATIOS.__dict__, "fe": 1.2})
    poor = bulk_draw.HostRatios(**{**EARTH_RATIOS.__dict__, "fe": 0.5})
    low = bulk_draw.dry_partition(poor, 0.08, 0.07, False).metal
    high = bulk_draw.dry_partition(rich, 0.08, 0.07, False).metal
    assert low < bulk_draw.dry_partition(EARTH_RATIOS, 0.08, 0.07, False).metal < high


def test_a_mantle_cannot_hold_more_iron_than_the_star_has():
    starved = bulk_draw.HostRatios(**{**EARTH_RATIOS.__dict__, "fe": 0.05})
    part = bulk_draw.dry_partition(starved, 0.18, 0.07, False)
    assert part.rock_feo < 0.18
    assert part.metal == pytest.approx(bulk_draw.dry_partition(starved, 0.18, 0.07, False).metal)
    assert 0.0 <= part.metal < 0.05


def test_carbon_rock_takes_silicon_carbide_and_graphite():
    silicate = bulk_draw.dry_partition(EARTH_RATIOS, 0.08, 0.07, False)
    carbon = bulk_draw.dry_partition(EARTH_RATIOS, 0.08, 0.07, True, graphite_over_rock=0.2)
    assert carbon.rock > silicate.rock - 0.2 and carbon.metal + carbon.rock == pytest.approx(1.0)


@pytest.mark.parametrize("planet_class, c_over_o, expected", [
    ("G", 0.5, "silicate"), ("G", 0.8, "carbide-graphite"), ("G", 1.1, "carbide-graphite"),
    ("C", 0.5, "carbide-graphite-by-class"), ("C", 0.9, "carbide-graphite"),
])
def test_the_rock_kind_follows_the_stars_c_over_o_and_the_class(planet_class, c_over_o, expected):
    for _ in range(3):
        assert bulk_draw.ROCK_KIND_TABLE.draw(planet_class, {"c_over_o": c_over_o}).name == expected


def test_every_terrestrial_class_has_exactly_one_rock_kind():
    from planetgen import tuning
    for letter, data in tuning.PLANET_CLASSES.items():
        for c_over_o in (0.1, 0.79, 0.8, 1.5):
            subject = bulk_draw.pt.Subject(letter, {"c_over_o": c_over_o})
            assert len(bulk_draw.ROCK_KIND_TABLE.passing(subject)) == 1


def test_a_class_radius_gets_the_mass_its_composition_needs():
    draw.set_run_seed(428)
    for planet_class, radius in (("G", 6371.0), ("K", 4000.0), ("E", 8000.0), ("C", 2000.0)):
        body = _body(planet_class, radius)
        bulk_draw.draw_terrestrial(body, True, False)
        bulk = bc.bulk_fractions(body)
        assert bulk.origin == "drawn" and sum(bulk.shares().values()) == pytest.approx(1.0)
        solved = bulk_structure.structure(bulk, body.mass, x_fe=bulk_draw._x_fe(bulk))
        assert solved.radius_km == pytest.approx(radius, rel=0.02)
        assert body.density == pytest.approx(body.mass / (4 / 3 * math.pi * (radius * 1000) ** 3) / 1000)


def test_an_earth_sized_body_around_a_sun_like_star_weighs_about_an_earth():
    draw.set_run_seed(1)
    masses = []
    for _ in range(40):
        body = _body("G", 6371.0)
        bulk_draw.draw_terrestrial(body, True, False)
        masses.append(body.mass / constants.EARTH_MASS_TO_KG)
    assert 0.7 < sum(masses) / len(masses) < 1.5


def test_a_given_mass_sets_the_radius():
    draw.set_run_seed(2)
    body = _body("G", 6000.0, mass=constants.EARTH_MASS_TO_KG)
    bulk_draw.draw_terrestrial(body, False, True)
    assert body.mass == constants.EARTH_MASS_TO_KG
    assert 5500.0 < body.radius < 7000.0


def test_class_s_and_a_body_given_both_are_worked_back():
    draw.set_run_seed(3)
    dense = _body("S", 9000.0)
    bulk_draw.draw_terrestrial(dense, True, False)
    assert dense.composition_origin == "inverse"
    low, high = 6.38, 9.45
    assert low <= dense.density <= high
    assert bc.bulk_fractions(dense).metal > 0.4

    both = _body("G", 6371.0, mass=constants.EARTH_MASS_TO_KG)
    bulk_draw.draw_terrestrial(both, True, True)
    assert both.composition_origin == "inverse" and both.mass == constants.EARTH_MASS_TO_KG and both.radius == 6371.0
    assert 0.15 < bc.bulk_fractions(both).dry_core_fraction < 0.6


def test_an_icy_small_body_is_a_porous_heap_of_ice_and_rock():
    draw.set_run_seed(4)
    body = _body("D", 150.0)
    bulk_draw.draw_terrestrial(body, True, False)
    bulk = bc.bulk_fractions(body)
    assert bulk.ice > 0.0 and bulk.water >= hydrosphere.WATER_MASS_FRACTION["D"][0] * 0.01 and 0.0 <= bulk.porosity <= 0.3
    assert body.mass < bulk_structure.SMALL_BODY_MASS_KG
    assert 0.5 < body.density < 3.0


def test_a_compact_remnant_makes_carbon_rock():
    assert bulk_draw.host_ratios(_star(metallicity_fe_h=None)) is bulk_draw.REMNANT_DISC_RATIOS
    body = _body("G", 6371.0, star=_star(metallicity_fe_h=None))
    draw.set_run_seed(6)
    bulk_draw.draw_terrestrial(body, True, False)
    assert body.carbon_rock


def test_the_same_seed_makes_the_same_body():
    results = []
    for _ in range(2):
        draw.set_run_seed(77)
        body = _body("G", 6371.0)
        bulk_draw.draw_terrestrial(body, True, False)
        results.append((body.mass, body.metal_mass_fraction, body.water_mass_fraction))
    assert results[0] == results[1]


def test_generated_bodies_store_a_composition_that_gives_their_density():
    seen_classes = set()
    for seed in range(8):
        with deterministic_entropy(seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets:
            if not hasattr(planet, "moons"):
                continue
            for body in [planet, *planet.moons]:
                bulk = bc.bulk_fractions(body)
                assert bulk is not None and bulk.origin in bc.ORIGINS
                assert sum(bulk.shares().values()) == pytest.approx(1.0, abs=1e-9)
                assert body.density == pytest.approx(body.mass / (4 / 3 * math.pi * (body.radius * 1000) ** 3) / 1000, rel=1e-9)
                if body.body_type == "g":
                    assert bulk.envelope > 0.0
                else:
                    assert bulk.envelope == 0.0 and body.water_mass_fraction == bulk.water
                seen_classes.add(body.planet_class)
    assert {"C", "D", "J"} <= seen_classes
