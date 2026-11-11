"""
GEN.1446: every body of a system raises a tide on every other (physics/tide_sources.py): the means of `1 / d^3`, the list of sources of a planet and a moon, the extra raisers of the
ocean tide, and the effect ledger entries.
"""

import math

import pytest

from planetgen.physics import constants, ocean_tides, tide_sources


def test_the_mean_of_inverse_cube_of_two_circles_is_the_quadrature():
    for a, b in ((1.0, 1.5), (1.0, 1.05), (2.0, 1.0), (1.0, 5.0)):
        n = 20000
        mean = sum((a * a + b * b - 2 * a * b * math.cos(2 * math.pi * (k + 0.5) / n)) ** -1.5 for k in range(n)) / n
        assert tide_sources.circles_mean_inverse_cube(a, b) == pytest.approx(mean, rel=1e-6)


def test_a_far_source_is_one_over_its_distance_cubed_and_the_kepler_mean_has_the_eccentricity_factor():
    assert tide_sources.circles_mean_inverse_cube(1.0, 100.0) == pytest.approx(1.0e-6, rel=1e-3)
    assert tide_sources.kepler_mean_inverse_cube(2.0, 0.0) == pytest.approx(1.0 / 8.0)
    assert tide_sources.kepler_mean_inverse_cube(2.0, 0.6) == pytest.approx(1.0 / (8.0 * 0.64 ** 1.5))
    assert tide_sources.orbit_mean_inverse_cube(1.0, 100.0, 0.5) > tide_sources.orbit_mean_inverse_cube(1.0, 100.0, 0.0)


class _Body:
    def __init__(self, mass, radius, distance=None, eccentricity=0.0, star=None, moons=(), name="b"):
        self.mass, self.radius, self.distance, self.eccentricity, self.star, self.moons, self.name = mass, radius, distance, eccentricity, star, list(moons), name


def _system():
    star = _Body(constants.SOLAR_MASS_TO_KG, 6.9e5, name="star")
    earth = _Body(constants.EARTH_MASS_TO_KG, 6371.0, 1.0, 0.0167, star, name="earth")
    moon = _Body(7.35e22, 1737.0, 0.00257, 0.05, star, name="moon")
    jupiter = _Body(317.8 * constants.EARTH_MASS_TO_KG, 69911.0, 5.2, 0.048, star, name="jupiter")
    earth.moons = [moon]
    return star, earth, moon, jupiter


def test_a_planet_is_raised_by_its_host_a_sibling_and_its_own_moon_and_the_sibling_is_tiny():
    star, earth, moon, jupiter = _system()
    sources = tide_sources.sources_of(earth, None, [earth, jupiter])
    by_name = {s.name: s for s in sources}
    assert set(by_name) == {"moon", "star", "jupiter"}
    accelerations = {name: tide_sources.acceleration_m_s2(earth, s) for name, s in by_name.items()}
    assert accelerations["moon"] > accelerations["star"] > accelerations["jupiter"] > 0.0
    assert accelerations["jupiter"] / accelerations["star"] < 1e-4


def test_only_the_sources_the_ocean_tide_did_not_have_become_raisers():
    star, earth, moon, jupiter = _system()
    sources = tide_sources.sources_of(earth, None, [earth, jupiter])
    extra = tide_sources.raisers_of(sources, tide_sources.counted_by_ocean_tides(earth, None))
    assert [r.name for r in extra] == ["jupiter"]
    assert isinstance(extra[0], ocean_tides.Raiser) and extra[0].distance_m > 3.0e11 and extra[0].eccentricity == 0.0
    assert extra[0].orbital_rate_rad_s > 0.0


def test_a_moon_is_raised_by_its_planet_the_host_the_siblings_and_the_other_moons():
    star, earth, moon, jupiter = _system()
    second = _Body(1.0e22, 1000.0, 0.0008, 0.01, star, name="inner moon")
    earth.moons.append(second)
    sources = tide_sources.sources_of(moon, earth, [earth, jupiter])
    names = [s.name for s in sources]
    assert names[0] == "earth" and {"star", "jupiter", "inner moon"} <= set(names)
    extra = tide_sources.raisers_of(sources, tide_sources.counted_by_ocean_tides(moon, earth))
    assert {r.name for r in extra} == {"jupiter", "inner moon"}


def test_a_lone_body_has_no_extra_raiser_and_the_ocean_tide_is_what_it_was():
    star, earth, moon, jupiter = _system()
    earth.moons = []
    sources = tide_sources.sources_of(earth, None, [earth])
    assert [s.name for s in sources] == ["star"]
    assert tide_sources.raisers_of(sources, tide_sources.counted_by_ocean_tides(earth, None)) == []


def test_extra_raisers_add_to_the_ranges_of_the_ocean_tide():
    base = ocean_tides.Raiser.of("star", constants.SOLAR_MASS_TO_KG, 1.496e11, constants.EARTH_MASS_TO_KG)
    extra = ocean_tides.Raiser("jupiter", 317.8 * constants.EARTH_MASS_TO_KG, 6.3e11, 0.0, 1.0e-8)
    args = (constants.EARTH_MASS_TO_KG, 6.371e6, 7.292e-5, [base], 0.3)
    without = ocean_tides.ocean_tide(*args)
    with_jupiter = ocean_tides.ocean_tide(constants.EARTH_MASS_TO_KG, 6.371e6, 7.292e-5, [base, extra], 0.3)
    assert with_jupiter.spring_range_m > without.spring_range_m
    assert with_jupiter.spring_range_m / without.spring_range_m - 1.0 < 1e-3


def test_register_tags_siblings_and_the_parent_and_the_spin_raisers_are_the_others_only():
    from planetgen.physics import spin_torques
    star, earth, moon, jupiter = _system()
    for body in (earth, jupiter):
        body.body_type = "t"
    moon.body_type = "t"
    tide_sources.register([earth, jupiter])
    assert earth._tide_siblings == (jupiter,) and moon._tide_parent is earth and moon._tide_siblings == (jupiter,)
    planet_raisers = tide_sources.spin_raisers(earth)
    assert len(planet_raisers) == 1 and isinstance(planet_raisers[0], spin_torques.Raiser)
    assert planet_raisers[0].mass_kg == pytest.approx(317.8 * constants.EARTH_MASS_TO_KG)
    moon_raisers = tide_sources.spin_raisers(moon)
    assert len(moon_raisers) == 2, "the host star at the planet's place and the sibling planet; the planet is the first raiser"
    tide_sources.register([earth])
    assert tide_sources.spin_raisers(earth) == ()


def test_an_inner_source_turns_the_line_at_the_receivers_rate_and_an_outer_one_at_its_own():
    star, earth, moon, jupiter = _system()
    inner = tide_sources.planet_sources(jupiter, [earth, jupiter])
    outer = tide_sources.planet_sources(earth, [earth, jupiter])
    n_earth = math.sqrt(constants.G * (star.mass + earth.mass) / (earth.distance * constants.AU_TO_M) ** 3)
    n_jupiter = math.sqrt(constants.G * (star.mass + jupiter.mass) / (jupiter.distance * constants.AU_TO_M) ** 3)
    assert {s.name: s.rate_rad_s for s in inner}["earth"] == pytest.approx(n_jupiter)
    assert {s.name: s.rate_rad_s for s in outer}["jupiter"] == pytest.approx(n_jupiter)
    assert {s.name: s.rate_rad_s for s in outer}["star"] == pytest.approx(n_earth)


def test_the_mean_inverse_sixth_of_two_circles_is_the_quadrature_and_a_close_pair_is_far_above_the_cube_squared():
    for a, b in ((1.0, 1.5), (1.0, 1.1), (2.0, 1.0)):
        n = 40000
        mean = sum((a * a + b * b - 2 * a * b * math.cos(2 * math.pi * (k + 0.5) / n)) ** -3 for k in range(n)) / n
        assert tide_sources.circles_mean_inverse_sixth(a, b) == pytest.approx(mean, rel=1e-6)
    assert tide_sources.circles_mean_inverse_sixth(1.0, 1.1) > 2.0 * tide_sources.circles_mean_inverse_cube(1.0, 1.1) ** 2


def test_heating_by_source_gives_the_orbit_to_the_primary_and_the_tide_to_the_others_and_nothing_without_a_day():
    star, earth, moon, jupiter = _system()
    earth.rotation_period_hours = 24.0
    earth.axial_tilt_deg = 23.0
    sources = tide_sources.sources_of(earth, None, [earth, jupiter])
    found = dict((getattr(b, "name", None), w) for b, w in tide_sources.heating_by_source(earth, sources, star))
    assert found["star"] > 0.0 and found["moon"] > 0.0 and found["jupiter"] > 0.0
    assert found["jupiter"] < found["moon"]
    earth.rotation_period_hours = None
    assert [b.name for b, w in tide_sources.heating_by_source(earth, sources, star)] == ["star"]


# GEN.1464 part B: belts and comets pull too

class _Belt:
    body_type = "a"

    def __init__(self, lower, upper, mass_kg, name="belt"):
        self.lower_limit, self.upper_limit, self.mass_kg, self.name = lower, upper, mass_kg, name


class _Comet:
    def __init__(self, q, e, mass_kg, name="comet"):
        self.perihelion_distance_au, self.eccentricity, self.mass_kg, self.name = q, e, mass_kg, name


def test_a_belt_is_one_source_whose_acceleration_is_the_sum_of_its_rings_and_whose_mass_is_the_belt_mass():
    from planetgen.physics import secular_chain
    star, earth, moon, jupiter = _system()
    belt = _Belt(2.0, 3.2, 3.0e21)
    source = tide_sources.belt_source(1.0, belt, star.mass, 1.9e-7)
    axes, masses = secular_chain.belt_rings(2.0, 3.2, 3.0e21)
    assert source.kind == "belt" and source.body is belt
    assert source.mass_kg == pytest.approx(3.0e21)       # no ring is within a tenth of an Earth orbit
    expected = sum(m * tide_sources.circles_mean_inverse_cube(1.0, a) for a, m in zip(axes, masses)) / 3.0e21 / constants.AU_TO_M ** 3
    assert source.g3 == pytest.approx(expected, rel=1e-12)
    # an independent read: 400 thin rings of one surface density (mass goes as a da), against the eight of the source
    edges = [2.0 + 1.2 * i / 400 for i in range(401)]
    fine = [(0.5 * (lo + hi), hi * hi - lo * lo) for lo, hi in zip(edges, edges[1:])]
    total = sum(w for _, w in fine)
    quadrature = sum(w * tide_sources.circles_mean_inverse_cube(1.0, a) for a, w in fine) / total
    assert source.g3 * constants.AU_TO_M ** 3 == pytest.approx(quadrature, rel=0.02)
    # rings raise their tides apart: M^2 g6 is the sum of m_i^2 g6_i, so it is under the point mass at the mean of 1/d^6
    assert 0.0 < source.g6 < source.g3 ** 2 * 1.0


def test_a_belt_outside_the_receiver_turns_the_tide_at_its_own_motion_and_one_inside_at_the_receivers():
    star, earth, moon, jupiter = _system()
    outside = tide_sources.belt_source(1.0, _Belt(2.0, 3.0, 1e21), star.mass, 123.0)
    inside = tide_sources.belt_source(5.0, _Belt(2.0, 3.0, 1e21), star.mass, 123.0)
    assert inside.rate_rad_s == pytest.approx(123.0)
    assert outside.rate_rad_s == pytest.approx(math.sqrt(constants.G * star.mass / (2.5 * constants.AU_TO_M) ** 3), rel=0.05)


def test_rings_next_to_the_receiver_are_left_out_and_a_belt_with_nothing_clear_gives_none():
    from planetgen.physics import secular_chain
    star, earth, moon, jupiter = _system()
    touching = tide_sources.belt_source(1.0, _Belt(0.5, 1.05, 1e21), star.mass, 1.9e-7)
    assert 0.0 < touching.mass_kg < 1e21
    assert tide_sources.belt_source(1.0, _Belt(1.0, 1.01, 1e21), star.mass, 1.9e-7) is None
    assert tide_sources.belt_source(1.0, _Belt(2.0, 3.0, 0.0), star.mass, 1.9e-7) is None
    assert secular_chain.BELT_ALPHA_LIMIT == 0.9      # the limit is the secular pull's own (reused)


def test_planet_and_moon_sources_list_the_belts_and_comets_among_the_siblings_and_nothing_else_changes():
    star, earth, moon, jupiter = _system()
    belt, comet = _Belt(2.0, 3.2, 3.0e21), _Comet(0.5, 0.9, 1.0e14)
    plain = tide_sources.sources_of(earth, None, [earth, jupiter])
    more = tide_sources.sources_of(earth, None, [earth, jupiter, belt, comet])
    assert [s.name for s in more[:len(plain)]] == [s.name for s in plain]
    assert [s.kind for s in more[len(plain):]] == ["belt", "comet"]
    assert [(s.mass_kg, s.g3) for s in more[:len(plain)]] == [(s.mass_kg, s.g3) for s in plain]
    moon_more = tide_sources.sources_of(moon, earth, [earth, jupiter, belt, comet])
    assert {"belt", "comet"} <= {s.kind for s in moon_more}
    cached = tide_sources.sources_of(moon, earth, [earth, jupiter, belt, comet], tide_sources.planet_sources(earth, [earth, jupiter, belt, comet]))
    assert [(s.name, s.g3) for s in cached] == [(s.name, s.g3) for s in moon_more]


def test_a_belt_pulls_far_less_than_the_host_and_far_more_than_nothing():
    star, earth, moon, jupiter = _system()
    belt = _Belt(2.0, 3.2, 3.0e21)
    sources = {s.kind: s for s in tide_sources.sources_of(earth, None, [belt])}
    ratio = tide_sources.acceleration_m_s2(earth, sources["belt"]) / tide_sources.acceleration_m_s2(earth, sources["star"])
    assert 1e-11 < ratio < 1e-8      # 3e21 kg against the Sun's 2e30 kg, at a mean separation of a few AU against one


def test_a_comet_is_a_source_on_its_own_orbit_and_an_unbound_or_massless_one_is_not():
    star, earth, moon, jupiter = _system()
    comet = _Comet(1.5, 0.8, 5.0e13)
    source = tide_sources.comet_source(1.0, comet, star.mass, 1.9e-7)
    axis = 1.5 / 0.2
    assert source.kind == "comet" and source.mass_kg == 5.0e13
    assert source.g3 == pytest.approx(tide_sources.orbit_mean_inverse_cube(1.0, axis, 0.8) / constants.AU_TO_M ** 3)
    assert source.rate_rad_s == pytest.approx(math.sqrt(constants.G * star.mass / (axis * constants.AU_TO_M) ** 3))
    assert tide_sources.comet_source(1.0, _Comet(1.5, 1.0, 5.0e13), star.mass, 1.9e-7) is None
    assert tide_sources.comet_source(1.0, _Comet(1.5, 1.4, 5.0e13), star.mass, 1.9e-7) is None
    assert tide_sources.comet_source(1.0, _Comet(1.5, 0.8, 0.0), star.mass, 1.9e-7) is None


def test_a_belt_becomes_one_raiser_at_the_distance_whose_inverse_cube_is_its_mean():
    star, earth, moon, jupiter = _system()
    sources = tide_sources.sources_of(earth, None, [_Belt(2.0, 3.2, 3.0e21)])
    extra = tide_sources.raisers_of(sources, tide_sources.counted_by_ocean_tides(earth, None))
    assert [r.name for r in extra] == ["belt"]
    belt = next(s for s in sources if s.kind == "belt")
    assert extra[0].mass_kg == belt.mass_kg and extra[0].distance_m ** -3 == pytest.approx(belt.g3)


def test_register_gives_every_planet_and_moon_the_belts_of_its_own_list_and_no_comets():
    star, earth, moon, jupiter = _system()
    belt = _Belt(2.0, 3.2, 3.0e21)
    tide_sources.register([earth, jupiter, belt])
    assert belt in earth._tide_siblings and belt in moon._tide_siblings and not hasattr(belt, "_tide_siblings")
    assert any(r for r in tide_sources.spin_raisers(earth))


def test_a_generated_system_with_a_belt_and_comets_notes_them_in_the_tide_ledger():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    config = SystemConfig()
    config.ASTEROID_BELT, config.COMETS = True, True
    with draw.bound(3003):
        system = StarSystem(system_config=config)
    belts = [p for p in system.planets if p.body_type == "a"]
    assert belts and system.comets
    noted = {id(source) for _, pairs in system.effect_ledger.contributions("tide") for source, _ in pairs}
    assert id(belts[0]) in noted and any(id(c) in noted for c in system.comets)


# GEN.1449: the sphere a planet's moons keep to counts every perturber

def test_the_perturber_ratio_is_the_other_sources_tide_over_the_hosts_and_a_lone_planet_has_none():
    star, earth, moon, jupiter = _system()
    assert tide_sources.perturber_ratio(earth, []) == 0.0
    ratio = tide_sources.perturber_ratio(earth, [jupiter])
    g3 = tide_sources.orbit_mean_inverse_cube(1.0, 5.2, 0.048) / constants.AU_TO_M ** 3
    expected = _mass_with_moons(jupiter) * g3 / (star.mass / constants.AU_TO_M ** 3)
    assert ratio == pytest.approx(expected, rel=1e-9)
    assert 1e-6 < ratio < 1e-4      # a Jupiter 4 AU out is a few parts in 10^5 of the Sun's tide on the Earth


def _mass_with_moons(body):
    return body.mass + sum(m.mass for m in body.moons)


def test_a_massive_close_neighbour_shrinks_the_sphere_of_the_moons_by_the_cube_root_of_one_plus_the_ratio():
    from planetgen.physics import planets as planet_physics
    star, earth, moon, jupiter = _system()
    earth.primary_mass_kg = star.mass
    earth.min_orbit_distance = None
    planet_physics.update_hill_sphere(earth)
    alone = earth.hill_radius
    near = _Body(10.0 * constants.EARTH_MASS_TO_KG * 300, 69911.0, 1.3, 0.0, star, name="near")
    ratio = tide_sources.perturber_ratio(earth, [near])
    earth._perturber_ratio = ratio
    planet_physics.update_hill_sphere(earth)
    assert ratio > 0.01
    assert earth.hill_radius == pytest.approx(alone * (1.0 + ratio) ** (-1.0 / 3.0), rel=1e-12)
    # the clearance the next body keeps is the star's Hill radius alone
    assert earth.min_orbit_distance == pytest.approx(5.0 * alone / constants.AU_TO_M * 1000.0, rel=1e-9)


def test_a_generated_system_gives_every_planet_its_perturber_ratio_and_a_sphere_no_bigger_than_the_stars():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics.orbits import calculate_hill_sphere
    from planetgen.util import draw
    config = SystemConfig()
    config.MAX_PLANETS = True
    with draw.bound(4001):
        system = StarSystem(system_config=config)
    planets = [p for p in system.planets if p.body_type != "a"]
    assert len(planets) > 2
    for planet in planets:
        assert planet._perturber_ratio >= 0.0
        host_hill_km = calculate_hill_sphere(planet.distance * constants.AU_TO_M, planet.mass, planet.star.mass) / 1000.0
        assert planet.hill_radius <= host_hill_km * (1.0 + 1e-12)
        assert planet.hill_radius == pytest.approx(host_hill_km * (1.0 + planet._perturber_ratio) ** (-1.0 / 3.0), rel=1e-9)
