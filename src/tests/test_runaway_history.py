# tests/test_runaway_history.py

"""
GEN.639: the runaway greenhouse limit depends on the planet's mass (the three Kopparapu et al. 2014 sets), and a body that was
past the limit while its star was young and bright loses water to the energy-limited escape of its steam, deterministically.
"""

import types

import pytest

from planetgen import tuning
from planetgen.generation import plausibility
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.physics import constants, hydrosphere, orbits, runaway_history
from planetgen.util import draw


def _generate_hydrosphere(body, internal_flux_w_m2):
    """`generate_hydrosphere` after the surface water is drawn, as `planets.calculate_atmospheric_conditions` does."""
    return hydrosphere.generate_hydrosphere(body, internal_flux_w_m2, hydrosphere.draw_surface_water(body) if body.body_type != "g" else 0.0)


EARTH_OCEAN_KG = hydrosphere.EARTH_OCEAN_KG if hasattr(hydrosphere, "EARTH_OCEAN_KG") else 1.4e21
EARTH = (constants.EARTH_MASS_TO_KG, 6371.0)


def _star(mass_sol, luminosity_sol, temperature, age_gy, yerkes="V", log_lx_lbol=None):
    star = types.SimpleNamespace(
        mass=mass_sol * constants.SOLAR_MASS_TO_KG, luminosity=luminosity_sol * constants.SOLAR_LUMINOSITY,
        temperature=temperature, age=age_gy, yerkes_class=yerkes)
    if log_lx_lbol is not None:
        star.log_lx_lbol = log_lx_lbol
    return star


def test_the_one_earth_mass_set_is_the_old_edge():
    for teff in (2800.0, 4000.0, 5780.0, 7200.0):
        assert orbits.runaway_greenhouse_flux(teff, 1.0) == pytest.approx(
            orbits.kopparapu_flux(orbits.KOPPARAPU_RUNAWAY_GREENHOUSE, teff))
    inner, outer = orbits.hz_limits(constants.SOLAR_LUMINOSITY, 5772.0, 1.0)
    assert (inner, outer) == pytest.approx((0.9504, 1.676), abs=2e-3)


def test_a_heavier_planet_keeps_its_water_to_a_higher_flux():
    fluxes = [orbits.runaway_greenhouse_flux(5780.0, mass) for mass in (0.1, 0.3, 1.0, 3.0, 5.0)]
    assert fluxes == sorted(fluxes)
    assert fluxes[0] == pytest.approx(0.99 + 0.0, abs=0.02) and fluxes[-1] == pytest.approx(1.188, abs=0.02)
    # outside the span the nearest set is taken
    assert orbits.runaway_greenhouse_flux(5780.0, 0.01) == orbits.runaway_greenhouse_flux(5780.0, 0.1)
    assert orbits.runaway_greenhouse_flux(5780.0, 50.0) == orbits.runaway_greenhouse_flux(5780.0, 5.0)


def test_the_inner_edge_moves_with_mass_and_the_outer_does_not():
    light = orbits.hz_limits(constants.SOLAR_LUMINOSITY, 5772.0, 0.1)
    heavy = orbits.hz_limits(constants.SOLAR_LUMINOSITY, 5772.0, 5.0)
    assert heavy[0] < light[0]  # a heavier planet keeps its water to a higher flux, so its edge is nearer the star
    assert light[1] == heavy[1]


def test_a_sun_like_star_leaves_the_earth_its_ocean():
    assert runaway_history.loss_kg(_star(1.0, 1.0, 5772.0, 4.6), 1.0, *EARTH) == 0.0


def test_a_planet_well_inside_the_zone_of_a_sun_like_star_loses_oceans():
    # the median Sun leaves saturation at 23 Myr (GEN.493): a Venus-orbit Earth loses about two oceans, the 100 Myr of the earlier table gave seven
    assert runaway_history.loss_kg(_star(1.0, 1.0, 5772.0, 4.6), 0.72, *EARTH) > 1 * EARTH_OCEAN_KG


def test_a_cool_star_dries_the_planets_of_its_present_zone():
    star = _star(0.3, 0.011, 3416.0, 5.0)     # a 0.3 solar mass star on the MIST main sequence at 5 Gy (GEN.652)
    zone = orbits.hz_limits(star.luminosity, star.temperature)
    assert runaway_history.loss_kg(star, zone[1] * 0.5 + zone[0] * 0.5, *EARTH) > 1 * EARTH_OCEAN_KG
    assert runaway_history.loss_kg(star, 1.0, *EARTH) == 0.0


def test_the_loss_grows_as_the_planet_moves_in():
    star = _star(0.3, 0.013, 3400.0, 5.0)
    losses = [runaway_history.loss_kg(star, au, *EARTH) for au in (0.3, 0.2, 0.15, 0.1)]
    assert losses == sorted(losses) and losses[0] < losses[-1]


def test_a_heavier_planet_loses_less_of_its_water_fraction():
    star = _star(0.3, 0.013, 3400.0, 5.0)
    small = runaway_history.loss_kg(star, 0.2, 0.5 * constants.EARTH_MASS_TO_KG, 5000.0) / (0.5 * constants.EARTH_MASS_TO_KG)
    big = runaway_history.loss_kg(star, 0.2, 4.0 * constants.EARTH_MASS_TO_KG, 9500.0) / (4.0 * constants.EARTH_MASS_TO_KG)
    assert big < small


def test_no_loss_where_the_past_is_not_followed():
    young = _star(0.3, 0.013, 3400.0, tuning.RUNAWAY_HISTORY_START_YR / 1e9 / 2)
    assert runaway_history.loss_kg(young, 0.1, *EARTH) == 0.0
    giant = _star(1.0, 100.0, 4500.0, 10.0, yerkes="III")
    assert runaway_history.loss_kg(giant, 0.5, *EARTH) == 0.0
    assert runaway_history.loss_kg(_star(0.3, 0.013, 3400.0, 5.0), None, *EARTH) == 0.0
    assert runaway_history.loss_kg(_star(12.0, 1e4, 30000.0, 0.01), 0.1, *EARTH) == 0.0


def test_the_loss_is_deterministic_and_follows_the_stars_own_activity():
    import math
    star = _star(0.3, 0.013, 3400.0, 5.0)
    assert runaway_history.loss_kg(star, 0.1, *EARTH) == runaway_history.loss_kg(star, 0.1, *EARTH)
    activity = runaway_history.activity
    t_sat, decay = activity.saturation_age_gy(0.3), activity.decay_exponent(0.3)
    nominal = activity.lx_lbol_at(5.0, activity.saturated_lx_lbol(0.3), t_sat, decay)
    # at 0.1 AU the planet is past the limit the whole time, so the decaying activity after saturation counts
    quiet = _star(0.3, 0.013, 3400.0, 5.0, log_lx_lbol=math.log10(nominal) - 0.7)
    loud = _star(0.3, 0.013, 3400.0, 5.0, log_lx_lbol=math.log10(nominal) + 0.7)
    assert (runaway_history.loss_kg(quiet, 0.1, *EARTH) < runaway_history.loss_kg(star, 0.1, *EARTH)
            < runaway_history.loss_kg(loud, 0.1, *EARTH))


def _planet(star_type, mass_earth, seed):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    with draw.bound(1):
        star = Star(cfg)
    with draw.bound(seed):
        planet = Planet(cfg, star, star.habitable_zone, plausibility.distance_for_zone(star, "e"),
                        mass=mass_earth * constants.EARTH_MASS_TO_KG, zone_override="e", moon_count=0)
    return star, planet


def test_the_surface_water_a_planet_keeps_is_what_its_class_gave_less_what_it_lost(monkeypatch):
    lost_any = False
    for seed in range(6):
        _, sunlike = _planet("G2V", 1.0, seed)
        with monkeypatch.context() as patch:
            patch.setattr(hydrosphere.runaway_history, "loss_kg", lambda *args: 0.0)
            _, sunlike_untouched = _planet("G2V", 1.0, seed)
        assert sunlike.ocean_fraction == sunlike_untouched.ocean_fraction  # the Sun-like zone loses nothing
        _, cool = _planet("M3V", 1.0, seed)
        with monkeypatch.context() as patch:
            patch.setattr(hydrosphere.runaway_history, "loss_kg", lambda *args: 0.0)
            _, untouched = _planet("M3V", 1.0, seed)
        assert cool.water_mass_fraction == untouched.water_mass_fraction  # the bulk water of the composition stays
        assert cool.ocean_fraction <= untouched.ocean_fraction  # the escape takes surface water, never makes it
        # the cover is 1.0 on a world that holds many oceans, so the loss shows in the surface water itself (TEST.504)
        with draw.bound(seed):
            surface = hydrosphere.draw_surface_water(cool)
        with monkeypatch.context() as patch, draw.bound(seed):                    # the same buffer draw with and without the loss
            patch.setattr(hydrosphere.runaway_history, "loss_kg", lambda *args: 0.0)
            surface_untouched = hydrosphere.draw_surface_water(cool)
        assert surface <= surface_untouched
        lost_any = lost_any or surface < surface_untouched
    assert lost_any  # the zone of a 0.3 solar-mass star lost some of its surface water while the star was young


def test_the_loss_comes_off_the_surface_first_and_the_mantle_gives_part_back(monkeypatch):
    class Body:
        body_type = "r"
        planet_class = "M"
        mass = constants.EARTH_MASS_TO_KG
        radius = 6371.0
        gravity = 1.0
        surface_temperature = 288.0
        atmospheric_pressure = 101325.0
        water_mass_fraction = 1e-3
        star = None
        _star_distance_au = None

    monkeypatch.setattr(hydrosphere, "log_uniform", lambda *args, **kwargs: 2.0)  # the mantle buffer
    monkeypatch.setattr(hydrosphere.runaway_history, "loss_kg", lambda *args: 2e-4 * Body.mass)
    layers = []
    real = hydrosphere.water_layer_depth_km
    monkeypatch.setattr(hydrosphere, "water_layer_depth_km", lambda mass, radius, fraction: layers.append(fraction) or real(mass, radius, fraction))
    body = Body()
    _generate_hydrosphere(body, 0.08)
    # 1e-3 in all, half of it at the surface (buffer 2): the loss of 2e-4 is the surface's, and the mantle's 5e-4 gives back
    # its share (0.3, at most the loss), so 5e-4 - 2e-4 + 1.5e-4 = 4.5e-4 lies at the surface; the bulk water is the composition's and stays
    assert layers[0] == pytest.approx(5e-4 - 2e-4 + min(2e-4, tuning.MANTLE_OUTGAS_SHARE * 5e-4))
    assert body.water_mass_fraction == 1e-3
