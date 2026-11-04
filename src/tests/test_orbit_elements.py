# tests/test_orbit_elements.py

"""GEN.727: the eccentricity and inclination of planets and moons, and what follows from them (the light, the seasons, the spacing)."""

import math
import types

import pytest
import numpy as np
from scipy import special, stats

from planetgen import tuning
from planetgen.generation import validation
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.system import StarSystem
from planetgen.physics import icy_tide, interior_heat, tidal_heating
from planetgen.physics import constants, insolation, kepler, orbit_elements as oe
from planetgen.physics import planets as planet_physics
from planetgen.util import draw
from tests.tide_support import rocky

SAMPLES = 4000


def _draws(count, function, seed=1):
    draw.set_run_seed(seed)
    return [function() for _ in range(count)]


def _half_gaussian_cdf(sigma):
    return lambda x: special.erf(x / (sigma * math.sqrt(2.0)))


# --- The laws ---------------------------------------------------------------------------------------------------------


def test_the_eccentricities_of_a_multi_planet_system_follow_the_published_half_gaussian():
    values = _draws(SAMPLES, lambda: oe.draw_planet_eccentricity(True, 0.9))
    assert max(values) < 0.9
    assert stats.kstest(values, _half_gaussian_cdf(oe.MULTI_SIGMA)).pvalue > 0.001
    assert sum(values) / SAMPLES == pytest.approx(oe.MULTI_SIGMA * math.sqrt(2.0 / math.pi), rel=0.05)


def test_the_eccentricities_of_a_single_planet_follow_the_published_mixture():
    values = _draws(SAMPLES, lambda: oe.draw_planet_eccentricity(False, 0.9))
    rayleigh_mean = oe.SINGLE_RAYLEIGH_MODE * math.sqrt(math.pi / 2.0)
    multi_mean = oe.MULTI_SIGMA * math.sqrt(2.0 / math.pi)
    expected_mean = oe.SINGLE_MODERATE_FRACTION * rayleigh_mean + (1.0 - oe.SINGLE_MODERATE_FRACTION) * multi_mean
    assert sum(values) / SAMPLES == pytest.approx(expected_mean, rel=0.06)

    def cdf(x):
        return (oe.SINGLE_MODERATE_FRACTION * (1.0 - np.exp(-x * x / (2.0 * oe.SINGLE_RAYLEIGH_MODE ** 2)))
                + (1.0 - oe.SINGLE_MODERATE_FRACTION) * _half_gaussian_cdf(oe.MULTI_SIGMA)(x))

    # the laws are cut at 0.9, which they hardly reach: the cut only rescales the cumulative curve by the little they leave beyond it
    scale = cdf(0.9)
    assert stats.kstest(values, lambda x: cdf(x) / scale).pvalue > 0.001
    multi = _draws(SAMPLES, lambda: oe.draw_planet_eccentricity(True, 0.9))
    assert sum(values) > 2.5 * sum(multi)       # a single planet is the hotter population [R1307]


def test_a_cut_eccentricity_is_never_exceeded_and_a_cut_of_zero_is_a_circle():
    assert max(_draws(2000, lambda: oe.draw_planet_eccentricity(False, 0.12))) <= 0.12
    assert max(_draws(2000, lambda: oe.draw_planet_eccentricity(True, 0.05))) <= 0.05
    assert oe.draw_planet_eccentricity(False, 0.0) == 0.0
    assert oe.draw_planet_eccentricity(False, -1.0) == 0.0
    assert max(_draws(500, lambda: oe.draw_planet_eccentricity(False, 5.0))) <= oe.MAX_ECCENTRICITY


def test_a_planets_inclination_is_rayleigh_with_the_dispersion_of_the_kepler_multiples():
    values = _draws(SAMPLES, oe.draw_planet_inclination_deg)
    sigma = oe.PLANET_INCLINATION_SIGMA_DEG
    assert sigma * math.sqrt(2.0) == pytest.approx(oe.MUTUAL_INCLINATION_SIGMA_DEG)
    assert stats.kstest(values, lambda x: 1.0 - np.exp(-x * x / (2.0 * sigma * sigma))).pvalue > 0.001
    assert min(values) >= 0.0


def test_a_regular_moon_takes_a_small_forced_eccentricity_and_a_captured_one_the_thermal_law():
    low, high = oe.REGULAR_MOON_ECCENTRICITY_RANGE
    for origin in ("disk", "impact", "ring", "binary", None):
        values = _draws(500, lambda: oe.draw_moon_eccentricity(origin, 0.5))
        assert low <= min(values) and max(values) <= high
    captured = _draws(SAMPLES, lambda: oe.draw_moon_eccentricity("captured", 0.4))
    assert max(captured) <= 0.4
    assert stats.kstest(captured, lambda x: (x / 0.4) ** 2).pvalue > 0.001   # f(e) = 2 e, cut at e_max
    assert oe.draw_moon_eccentricity("captured", 0.0) == 0.0


def test_a_moon_at_the_edge_of_the_stable_sphere_has_no_room_for_an_eccentricity():
    # [R1138]: a_moon <= 0.4895 R_H (1 - 1.0305 e_planet - 0.2738 e_moon)
    hill = 0.001
    stable = 0.4895 * hill
    assert oe.max_moon_eccentricity(stable, hill) == 0.0
    assert oe.max_moon_eccentricity(stable * 1.2, hill) == 0.0
    assert oe.max_moon_eccentricity(stable / 1.2, hill) == pytest.approx((1.0 - 1.0 / 1.2) / 0.2738)
    assert oe.max_moon_eccentricity(stable / 4.0, hill) == oe.MAX_ECCENTRICITY
    assert oe.max_moon_eccentricity(stable, None) == oe.MAX_ECCENTRICITY
    # the planet's own eccentricity takes room from the moon's
    assert oe.max_moon_eccentricity(stable / 1.5, hill, 0.2) < oe.max_moon_eccentricity(stable / 1.5, hill, 0.0)
    assert oe.max_moon_eccentricity(stable / 1.5, hill, 0.2) == pytest.approx(
        (1.0 - 1.0305 * 0.2 - 1.0 / 1.5) / 0.2738)
    assert oe.max_moon_eccentricity(stable / 1.2, hill, 0.2) == 0.0   # the planet's eccentricity leaves no room at all
    assert oe.max_moon_eccentricity(stable * 0.9, hill, 0.5) == 0.0


def test_a_planet_cannot_pass_the_roche_limit_of_its_star():
    star = types.SimpleNamespace(radius=696_000.0, mass=1.989e30)
    planet = types.SimpleNamespace(radius=6371.0, mass=5.972e24)
    star_density = 1.989e30 / (4.0 / 3.0 * math.pi * 6.96e8 ** 3)
    planet_density = 5.972e24 / (4.0 / 3.0 * math.pi * 6.371e6 ** 3)
    limit = oe.roche_limit_au(696_000.0, star_density, planet_density)
    assert limit == pytest.approx(oe.ROCHE_COEFFICIENT * 696_000.0 / constants.AU_TO_KM * (star_density / planet_density) ** (1 / 3))
    assert 0.003 < limit < 0.02                        # a few hundredths of an AU round the Sun, for a planet as dense as the Earth
    assert oe.max_planet_eccentricity(limit, star, planet) == 0.0
    assert oe.max_planet_eccentricity(limit / 2.0, star, planet) == 0.0
    a = 10.0 * limit
    cut = oe.max_planet_eccentricity(a, star, planet)
    assert a * (1.0 - cut) == pytest.approx(limit)
    assert oe.max_planet_eccentricity(1e6, star, planet) == oe.MAX_ECCENTRICITY
    assert oe.max_planet_eccentricity(1.0, types.SimpleNamespace(radius=None, mass=None), planet) == oe.MAX_ECCENTRICITY


# --- The light ---------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("e", [0.0, 0.1, 0.3, 0.6, 0.9])
def test_the_mean_light_and_its_swing_are_those_of_an_eccentric_orbit(e):
    count = 20000
    flux = []
    for n in range(count):
        _nu, r = kepler.true_anomaly_and_distance_elliptical(2.0 * math.pi * (n + 0.5) / count, e, 1.0)
        flux.append(1.0 / (r * r))
    mean = sum(flux) / count                     # uniform in mean anomaly is uniform in time
    assert mean == pytest.approx(oe.mean_flux_factor(e), rel=1e-3)
    assert (max(flux) - min(flux)) / mean == pytest.approx(4.0 * e / (1.0 - e * e) ** 1.5, rel=1e-3)
    assert oe.flux_swing(e) == pytest.approx(min(4.0 * e / (1.0 - e * e) ** 1.5, 1.0))
    assert oe.flux_distance_au(2.0, e) ** -2 == pytest.approx(2.0 ** -2 * oe.mean_flux_factor(e))


def test_a_circle_has_no_swing_and_the_swing_never_passes_one():
    assert oe.flux_swing(0.0) == 0.0 and oe.flux_swing(None) == 0.0
    assert oe.mean_flux_factor(None) == 1.0
    assert oe.flux_swing(0.9) == 1.0


def _airless_planet():
    draw.set_run_seed(5)
    for _ in range(400):
        system = StarSystem(system_config=SystemConfig())
        for planet in system.planets:
            if isinstance(planet, Planet) and planet.atmosphere == "None" and planet.body_type == "t" and not planet.is_moon:
                return planet
    pytest.fail("no airless planet found")


def _airless_seed(planet, first, *runs):
    """The first seed from `first` whose draws leave the body with no air at all in every run (each run is the arguments of `calculate_atmospheric_conditions`):
    a thin air has a greenhouse that moves with the temperature, so the fourth-root law of the light holds only for a body with none."""
    for seed in range(first, first + 200):
        for arguments in runs:
            draw.set_run_seed(seed)
            planet_physics.calculate_atmospheric_conditions(planet, *arguments)
            if planet.atmospheric_pressure:
                break
        else:
            return seed
    pytest.fail("no draw leaves the body airless")


def _absorbed_light_ratio(planet, base_temperature, base_albedo):
    """How much more light the body absorbs than the base run, from its temperatures: an airless body is at the fourth root of the light it absorbs, `T^4 ~ L (1 - A)`.
    The albedo is the body's own and moves with the temperature (warmer ground loses its ice and frost), so it is read after each run and the law is held to the light, not to the temperature."""
    return (planet.surface_temperature / base_temperature) ** 4 * (1.0 - base_albedo) / (1.0 - planet.bond_albedo)


def test_an_eccentric_planet_is_as_warm_as_a_circle_with_the_same_mean_light():
    planet = _airless_planet()
    planet.eccentricity = 0.0
    draw.set_run_seed(77)
    planet_physics.calculate_atmospheric_conditions(planet)
    circle, circle_albedo = planet.surface_temperature, planet.bond_albedo
    planet.eccentricity = 0.6
    draw.set_run_seed(77)
    planet_physics.calculate_atmospheric_conditions(planet)
    # the light an orbit gives, averaged over it, is that of a circle of radius a (1 - e^2)^(1/4): 1 / sqrt(1 - e^2) times a circle's (GEN.727). The albedo is not the same
    # (the warmer eccentric planet has less ice), so the light is found from the temperature and the body's own albedo.
    assert _absorbed_light_ratio(planet, circle, circle_albedo) == pytest.approx((1.0 - 0.36) ** -0.5, rel=1e-6)


def test_a_moon_gets_the_light_of_its_parents_orbit():
    planet = _airless_planet()
    planet.is_moon = True            # an airless body, its light taken from the orbit round the star the parent has
    planet.eccentricity = 0.9        # its own orbit round its parent does not matter
    planet.host_eccentricity = 0.0
    seed = _airless_seed(planet, 78, (planet.distance,), (planet.distance, 0.5))
    draw.set_run_seed(seed)
    planet_physics.calculate_atmospheric_conditions(planet, planet.distance, 0.0)
    base, base_albedo = planet.surface_temperature, planet.bond_albedo
    draw.set_run_seed(seed)
    planet_physics.calculate_atmospheric_conditions(planet, planet.distance, 0.5)
    assert planet.host_eccentricity == 0.5
    assert _absorbed_light_ratio(planet, base, base_albedo) == pytest.approx((1.0 - 0.25) ** -0.5, rel=1e-6)


def test_the_orbits_swing_is_the_seasons_of_a_lone_star_and_adds_to_a_companions():
    system = _first_system(lambda s: s.binary_type is None and any(
        isinstance(p, Planet) and p.body_type == "t" and p.eccentricity > 0.02 for p in s.planets))
    for planet in system.planets:
        if isinstance(planet, Planet) and planet.body_type == "t":
            assert planet.flux_variation == pytest.approx(oe.flux_swing(planet.eccentricity))
            assert all(moon.flux_variation == planet.flux_variation for moon in planet.moons)
    circle = insolation.companion_swing(1.0, 1.0, 0.5, 20.0, 30.0)
    assert insolation.companion_swing(1.0, 1.0, 0.5, 20.0, 30.0, 0.0) == circle
    assert insolation.companion_swing(1.0, 1.0, 0.5, 20.0, 30.0, 0.1) > circle
    # the in-phase sum: peak over trough of both lights, against the mean of both
    host, companion = 1.0 / 1.0 ** 2, 0.5
    mean = host * oe.mean_flux_factor(0.1) + companion * insolation.mean_inverse_square(25.0, 0.2)
    peak = host / 0.9 ** 2 + companion / 20.0 ** 2
    trough = host / 1.1 ** 2 + companion / 30.0 ** 2
    assert insolation.companion_swing(1.0, 1.0, 0.5, 20.0, 30.0, 0.1) == pytest.approx((peak - trough) / mean)


# --- The tides ---------------------------------------------------------------------------------------------------------


def _body(mass_kg, radius_km, body_type="t"):
    if body_type == "g" or mass_kg is None:
        return types.SimpleNamespace(mass=mass_kg, radius=radius_km, body_type=body_type)
    return rocky(mass_kg, radius_km)


def test_a_hot_earth_is_circularised_in_a_billion_years_and_a_far_one_keeps_its_orbit():
    earth = _body(5.972e24, 6371.0)
    sun = 1.989e30
    near = oe.circularisation_time_yr(earth, sun, 0.02)
    far = oe.circularisation_time_yr(earth, sun, 1.0)
    assert 1e5 < near < 1e8
    assert far > 1e15
    assert oe.circularised(0.3, 5e9, near) < 1e-6
    assert oe.circularised(0.3, 5e9, far) == pytest.approx(0.3, rel=1e-5)
    # tau goes as a^(13/2) at a fixed Q (the fifth power of a/R and the 1/n of the period); the mantle's Q rises with the tide's frequency, so a closer
    # orbit is damped a little more slowly than that law says
    assert 2.0 ** 5.0 < oe.circularisation_time_yr(earth, sun, 0.04) / near < 2.0 ** 6.5
    assert oe.circularised(0.3, 1e9, 0.0) == 0.0
    assert oe.circularisation_time_yr(_body(None, 1.0), sun, 1.0) == math.inf


def test_a_gas_giant_dissipates_less_than_a_rock_of_its_size():
    giant = _body(1.9e27, 70_000.0, "g")
    n = 1.0e-7
    assert interior_heat.tidal_response(giant, n) == (tidal_heating.GIANT_K2, tidal_heating.GIANT_TIDAL_Q)
    k2, q = interior_heat.tidal_response(_body(5.972e24, 6371.0), n)
    assert q < tidal_heating.GIANT_TIDAL_Q and 0.1 < k2 < 1.5
    icy = types.SimpleNamespace(mass=1.0e22, radius=1000.0, body_type="t", water_mass_fraction=0.4)
    k2, q = interior_heat.tidal_response(icy, n)      # the ice shell's (GEN.960)
    assert k2 == pytest.approx(icy_tide.love_k2(1.0e22, 1.0e6, n)) and 1.0 < q < math.inf


# --- Generated systems -------------------------------------------------------------------------------------------------


def _first_system(predicate, seed=11, tries=600):
    draw.set_run_seed(seed)
    for _ in range(tries):
        system = StarSystem(system_config=SystemConfig())
        if predicate(system):
            return system
    pytest.fail("no system matching")


def _bodies(system):
    for planet in system.planets + getattr(system, "secondary_planets", []):
        if isinstance(planet, Planet):
            yield planet
            yield from planet.moons


def test_every_generated_orbit_is_the_ellipse_it_was_drawn_on():
    draw.set_run_seed(21)
    seen_eccentric = seen_moons = 0
    for _ in range(60):
        system = StarSystem(system_config=SystemConfig())
        for body in _bodies(system):
            assert 0.0 <= body.eccentricity <= oe.MAX_ECCENTRICITY
            orbit = body.orbit_from_vector()
            assert orbit["eccentricity"] == pytest.approx(body.eccentricity, abs=1e-7)
            assert orbit["semi_major_axis"] == pytest.approx(body.distance, rel=1e-9)
            assert orbit["periapsis_distance"] == pytest.approx(validation.periapsis_au(body), rel=1e-9)
            seen_eccentric += body.eccentricity > 0.05
            seen_moons += body.is_moon
        for planet in system.planets:
            if isinstance(planet, Planet):
                assert planet.orbital_inclination_deg >= 0.0
                assert validation.check_lunar_system(planet) == []
    assert seen_eccentric > 0 and seen_moons > 0


def test_the_orbits_of_neighbours_do_not_cross():
    draw.set_run_seed(22)
    for _ in range(60):
        system = StarSystem(system_config=SystemConfig())
        for planets in (system.planets, getattr(system, "secondary_planets", [])):
            assert validation.check_orbit_spacing(planets) == []
            for inner, outer in zip(planets, planets[1:]):
                if inner.body_type != "a" and outer.body_type != "a":
                    assert validation.apoapsis_au(inner) < validation.periapsis_au(outer)


def test_the_spacing_leaves_room_for_the_ellipses_and_is_the_old_one_for_circles():
    star = types.SimpleNamespace(mass=1.989e30)

    def pair(e_in, e_out):
        inner = types.SimpleNamespace(body_type="t", mass=6e24, star=star, distance=1.0, eccentricity=e_in)
        outer = types.SimpleNamespace(body_type="t", mass=6e24, star=star, distance=1.5, eccentricity=e_out)
        return outer, inner

    kappa = tuning.MUTUAL_HILL_RADII_SEPARATION * (2 * 6e24 / (3 * 1.989e30)) ** (1 / 3)
    outer, inner = pair(0.0, 0.0)
    assert validation.min_distance_after_au(outer, [inner]) == pytest.approx((1 + kappa / 2) / (1 - kappa / 2))
    outer, inner = pair(0.1, 0.05)
    needed = validation.min_distance_after_au(outer, [inner])
    # the gap from the inner apoapsis to the outer periapsis holds the separation
    assert needed * (1 - 0.05) - 1.0 * (1 + 0.1) == pytest.approx(kappa * (needed + 1.0) / 2.0)
    assert needed * (1 - 0.05) > 1.0 * (1 + 0.1)


def test_a_belt_clears_the_apoapsis_of_a_planet_inside_it_and_a_planet_the_belt_by_its_periapsis():
    star = types.SimpleNamespace(mass=1.989e30)
    planet = types.SimpleNamespace(body_type="t", distance=1.0, eccentricity=0.2, min_orbit_distance=0.05, hill_radius=3e5, mass=6e24,
                                   star=star)
    belt = types.SimpleNamespace(body_type="a", lower_limit=1.0, upper_limit=1.4, distance=1.0)
    assert validation.min_distance_after_au(belt, [planet]) == pytest.approx(1.2 + 0.05)
    outer = types.SimpleNamespace(body_type="t", distance=2.0, eccentricity=0.3, hill_radius=3e5, mass=6e24, star=star)
    needed = validation.min_distance_after_au(outer, [belt])
    assert needed * (1 - 0.3) > belt.upper_limit
    circle = types.SimpleNamespace(**{**outer.__dict__, "eccentricity": 0.0})
    assert validation.min_distance_after_au(circle, [belt]) < needed


def test_the_reach_of_a_planet_is_its_apoapsis_and_its_hill_sphere():
    planet = types.SimpleNamespace(body_type="t", distance=2.0, eccentricity=0.25, min_orbit_distance=0.1)
    assert validation.reach_au(planet) == pytest.approx(2.5 + 0.1)
    assert validation.exceeds_orbit_ceiling(planet, 2.55) and not validation.exceeds_orbit_ceiling(planet, 2.7)
