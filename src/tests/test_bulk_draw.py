"""
GEN.428: the draw of a planet's or moon's bulk composition and the mass solved from it (`planetgen.physics.bulk_draw`).
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen import tuning
from planetgen.physics import bulk_composition as bc, bulk_draw, bulk_structure, constants
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

EARTH_RATIOS = bulk_draw.HostRatios(fe=0.81, si=0.91, al=0.07, ca=0.06, ni=0.045, c_over_o=0.55, y_envelope=0.27, z=0.014)


def _star(**extra):
    fields = dict(metallicity_fe_h=0.0, age=4.5, mass=2e30, galactic_center_dist_ly=26000.0, luminosity=3.8e26, snow_line_au=2.7)
    fields.update(extra)
    return SimpleNamespace(**fields)


def _body(radius_km, icy=False, **extra):
    """A body of the given radius that formed inside the snow line, or (`icy`) well beyond it."""
    fields = dict(radius=radius_km, mass=None, star=_star(), is_moon=False, distance=20.0 if icy else 1.0)
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


@pytest.mark.parametrize("c_over_o, expected", [(0.5, "silicate"), (0.79, "silicate"), (0.8, "carbide-graphite"), (1.1, "carbide-graphite")])
def test_the_rock_kind_follows_the_stars_c_over_o(c_over_o, expected):
    for origin in (bulk_draw.ROCKY, bulk_draw.ICY):
        for _ in range(3):
            assert bulk_draw.ROCK_KIND_TABLE.draw(origin, {"c_over_o": c_over_o}).name == expected


def test_every_origin_has_exactly_one_rock_kind():
    for origin in (bulk_draw.ROCKY, bulk_draw.ICY):
        for c_over_o in (0.1, 0.79, 0.8, 1.5):
            subject = bulk_draw.pt.Subject(origin, {"c_over_o": c_over_o})
            assert len(bulk_draw.ROCK_KIND_TABLE.passing(subject)) == 1


def test_a_body_is_icy_by_where_it_formed():
    assert bulk_draw.origin_of(_body(6371.0)) == bulk_draw.ROCKY
    assert bulk_draw.origin_of(_body(6371.0, icy=True)) == bulk_draw.ICY
    assert bulk_draw.origin_of(_body(6371.0, formation_distance_au=5.0)) == bulk_draw.ICY
    assert bulk_draw.origin_of(_body(6371.0, icy=True, formation_distance_au=1.0)) == bulk_draw.ROCKY


def test_a_radius_gets_the_mass_its_composition_needs():
    draw.set_run_seed(428)
    for radius in (6371.0, 4000.0, 8000.0, 2000.0):
        body = _body(radius)
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
        body = _body(6371.0)
        bulk_draw.draw_terrestrial(body, True, False)
        masses.append(body.mass / constants.EARTH_MASS_TO_KG)
    assert 0.7 < sum(masses) / len(masses) < 1.5


def test_a_given_mass_sets_the_radius():
    draw.set_run_seed(2)
    body = _body(6000.0, mass=constants.EARTH_MASS_TO_KG)
    bulk_draw.draw_terrestrial(body, False, True)
    assert body.mass == constants.EARTH_MASS_TO_KG
    assert 5500.0 < body.radius < 7000.0


def test_a_body_given_both_is_worked_back_and_a_radius_alone_is_drawn():
    draw.set_run_seed(3)
    alone = _body(9000.0)
    bulk_draw.draw_terrestrial(alone, True, False)
    assert alone.composition_origin == "drawn"

    both = _body(6371.0, mass=constants.EARTH_MASS_TO_KG)
    bulk_draw.draw_terrestrial(both, True, True)
    assert both.composition_origin == "inverse" and both.mass == constants.EARTH_MASS_TO_KG and both.radius == 6371.0
    assert 0.15 < bc.bulk_fractions(both).dry_core_fraction < 0.6


def test_an_icy_small_body_is_a_porous_heap_of_ice_and_rock():
    draw.set_run_seed(4)
    body = _body(150.0, icy=True)
    bulk_draw.draw_terrestrial(body, True, False)
    bulk = bc.bulk_fractions(body)
    assert bulk.ice > 0.0 and bulk.water >= tuning.ICY_WATER_MASS_FRACTION_FLOOR * 0.01 and 0.0 <= bulk.porosity <= 0.3
    assert body.mass < bulk_structure.SMALL_BODY_MASS_KG
    assert 0.5 < body.density < 3.0


def test_a_compact_remnant_makes_carbon_rock():
    assert bulk_draw.host_ratios(_star(metallicity_fe_h=None)) is bulk_draw.REMNANT_DISC_RATIOS
    body = _body(6371.0, star=_star(metallicity_fe_h=None))
    draw.set_run_seed(6)
    bulk_draw.draw_terrestrial(body, True, False)
    assert body.carbon_rock


def test_the_same_seed_makes_the_same_body():
    results = []
    for _ in range(2):
        draw.set_run_seed(77)
        body = _body(6371.0)
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
    assert len(seen_classes) >= 3


def test_an_icy_bodys_water_tops_out_at_the_condensed_share_of_its_formation_temperature():
    """GEN.950: the water of an icy body is at most the water share of the solids that its disk temperature condenses (Lodders 2003 Table 11)."""
    from planetgen.physics import condensate_shares
    near = _body(1000.0, icy=True, distance=2.7)    # just beyond the snow line of 2.7 AU: 182 K, water alone
    _, top = bulk_draw._icy_water_range(near)
    assert top == pytest.approx(condensate_shares.water_share_of_solids(181.9), abs=1e-3) and top == pytest.approx(0.539, abs=0.002)
    far = _body(1000.0, icy=True, distance=2.7 * (182.0 / 100.0) ** 2)   # 100 K: ammonia is in as well
    _, top = bulk_draw._icy_water_range(far)
    assert top == pytest.approx(0.494, abs=0.002)
    farthest = _body(1000.0, icy=True, distance=2.7 * (182.0 / 30.0) ** 2)   # 30 K: methane too
    low, top = bulk_draw._icy_water_range(farthest)
    assert top == pytest.approx(0.384, abs=0.002) and low == tuning.ICY_WATER_MASS_FRACTION_FLOOR
