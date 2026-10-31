# tests/test_final_age_refresh.py

"""
GEN.839: the ice shell and the volcanism gate follow the heat at the star's final age. A star's age only rises after its planets and
moons are made; the life pass runs the heat again at the final age (GEN.763) and then `planets.refresh_at_final_age` makes the hydrosphere,
the ground, the ground dose and the habitability score again from it, on the body's own stream, so no draw of the run's stream moves.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, hydrosphere, planets, surface_composition
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy


ICY_SEEDS = (5, 15, 19, 25, 39)
"""The seeds whose systems hold ice-covered oceans (the life pass moves many moons in, so a short range of seeds holds few)."""


def bodies(seeds):
    for seed in seeds:
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            for body in [planet, *getattr(planet, "moons", [])]:
                if body.body_type not in ("a", "g"):
                    yield body


def shell_km(body):
    """The ice shell the hydrosphere gives a body from the heat it stores (its own and its tide's)."""
    flux = (body.internal_heat_flux_w_m2 or 0.0) + (getattr(body, "tidal_flux_w_m2", None) or 0.0)
    return hydrosphere.pressure_melted_shell_km(max(flux, 1e-6), max(body.surface_temperature, 1.0),
                                                body.gravity * constants.EARTH_GRAVITY)


def icy(seeds):
    # a body with little water has a brine whose eutectic (drawn again at the refresh) is where its shell melts, not the 273.15 K `shell_km` takes
    return [body for body in bodies(seeds) if body.hydrosphere == "ice-covered ocean" and body.ocean_class != "chloride brine"
            and body._surface_w >= hydrosphere.BRINE_MAX_WATER_FRACTION and body.ice_shell_km and body.internal_heat_flux_w_m2]


def test_the_stored_ice_shell_is_the_one_the_stored_final_age_heat_gives():
    found = icy(ICY_SEEDS)
    assert len(found) > 10
    lifeless = [body for body in found if getattr(body, "life_class", None) is None]    # life changes the air after the refresh, by a few kelvin the shell does not follow
    assert len(lifeless) > 5
    for body in lifeless:
        assert body.ice_shell_km == pytest.approx(shell_km(body), rel=1e-9), (body.planet_class, body.internal_heat_flux_w_m2)


def test_a_body_made_again_at_an_older_age_has_the_shell_of_its_cooler_interior():
    # a thick shell on a heavy body makes `pressure_melted_shell_km`'s three passes alternate, and one that grows past the body's water (shell and ocean
    # together) freezes it through: the body keeps its ocean well below the thicker shell the cooler interior gives
    body = next(body for body in icy(ICY_SEEDS) if body.ice_shell_km < 15.0 and body.ocean_depth_km > 40.0)
    before_flux, before_shell = body.internal_heat_flux_w_m2, body.ice_shell_km
    body.star.age = body.star.age * 3.0
    planets.refresh_at_final_age(body)
    assert body.internal_heat_flux_w_m2 < before_flux
    assert body.ice_shell_km == pytest.approx(shell_km(body), rel=1e-9)
    assert body.ice_shell_km > before_shell


def test_the_ground_is_made_again_from_the_heat_at_the_final_age(monkeypatch):
    body = next(iter(bodies([12])))
    seen = []
    real = surface_composition.draw_surface
    monkeypatch.setattr(surface_composition, "draw_surface", lambda inputs: seen.append(inputs.internal_heat_flux_w_m2) or real(inputs))
    body.star.age = body.star.age * 2.0
    planets.refresh_at_final_age(body)
    assert seen == [body.internal_heat_flux_w_m2]


def test_no_draw_of_the_runs_stream_moves():
    body = next(iter(bodies([12])))
    with draw.bound(draw.Stream("untouched")):
        planets.refresh_at_final_age(body)
        after = draw.random()
    with draw.bound(draw.Stream("untouched")):
        assert draw.random() == after


def test_the_refresh_is_deterministic_and_repeatable():
    first, second = next(iter(bodies([17]))), next(iter(bodies([17])))
    planets.refresh_at_final_age(first)
    planets.refresh_at_final_age(second)
    planets.refresh_at_final_age(second)
    for name in (*hydrosphere.HYDROSPHERE_FIELDS, "ground_type", "dose_ground_msv_yr", "habitability_score"):
        assert getattr(first, name, None) == getattr(second, name, None), name


def test_a_giant_a_belt_and_a_body_with_no_water_drawn_are_left_alone():
    class Bare:
        body_type = "g"
    planets.refresh_at_final_age(Bare())
    body = next(iter(bodies([12])))
    body._surface_w = None
    before = body.ice_shell_km
    planets.refresh_at_final_age(body)
    assert body.ice_shell_km == before
