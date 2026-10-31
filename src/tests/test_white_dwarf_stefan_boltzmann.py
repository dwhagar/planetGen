# tests/test_white_dwarf_stefan_boltzmann.py

"""
GEN.365: the radius, luminosity and temperature stored for a white dwarf obey L = 4 pi R^2 sigma T^4. A white dwarf's radius is the degenerate star's (its mass sets it), and the
temperature stored is the rounded one or the forced type's, so `stellar_evolution.star_params` takes the luminosity from those two. Before, a forced type's temperature (the central star
of a planetary nebula, O3VII to B0VII) or the rounding of a cool one left the radius up to 10 percent off the radius that L and T give.
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.physics.stellar_evolution import sample_living_star, star_params
from planetgen.util import draw


def _ratio(radius_km, luminosity_w, temperature_k):
    """The stored radius over the one that L and T give."""
    return radius_km / (math.sqrt(luminosity_w / (constants.FOUR_PI * constants.STEFAN_BOLTZMANN_CONSTANT)) / temperature_k ** 2 / 1000.0)


def _star_ratio(star):
    return _ratio(star.radius, star.luminosity, star.temperature)


@pytest.mark.parametrize("star_type", ["O3VII", "O5VII", "O7VII", "O9VII", "B0VII", "A0VII", "F5VII", "K0VII"])
def test_a_forced_white_dwarf_obeys_stefan_boltzmann(star_type):
    for seed in range(25):
        config = SystemConfig()
        config.STAR_TYPE = star_type
        with draw.bound(seed):
            star = Star(config)
        assert star.yerkes_class == "VII"
        assert _star_ratio(star) == pytest.approx(1.0, abs=1e-6), (star_type, seed)


def test_a_hundred_hot_white_dwarfs_and_a_hundred_natural_ones_obey_it_within_one_percent():
    hot, natural = [], []
    for seed in range(100):
        config = SystemConfig()
        config.STAR_TYPE = ("O3VII", "O5VII", "O9VII", "B0VII")[seed % 4]
        with draw.bound(seed):
            hot.append(_star_ratio(Star(config)))
    draw.set_run_seed(4)
    while len(natural) < 100:
        mass, age, state = sample_living_star()
        if state["yerkes_class"] == "VII":
            params = star_params(mass, age, state)
            natural.append(_ratio(params["radius_km"], params["luminosity_w"], params["temperature_k"]))
    assert all(abs(r - 1.0) < 0.01 for r in hot + natural), (min(hot + natural), max(hot + natural))


def test_the_white_dwarfs_of_multi_star_systems_obey_it():
    seen = 0
    for stars, seeds in ((3, range(0, 60)), (4, range(100, 130))):
        for seed in seeds:
            config = SystemConfig()
            config.PLANETS = False
            config.STARS = stars
            with draw.bound(seed):
                system = StarSystem(system_config=config)
            for star in system.stars:
                if star.yerkes_class == "VII":
                    seen += 1
                    assert _star_ratio(star) == pytest.approx(1.0, abs=1e-6), (stars, seed)
    assert seen >= 5


def test_the_luminosity_is_from_the_radius_and_the_stored_temperature_not_the_cooling_sequences_unrounded_one():
    draw.set_run_seed(7)
    for _ in range(200):
        mass, age, state = sample_living_star()
        if state["yerkes_class"] != "VII":
            continue
        params = star_params(mass, age, state)
        expected = constants.FOUR_PI * (state["radius_km"] * 1000.0) ** 2 * constants.STEFAN_BOLTZMANN_CONSTANT * params["temperature_k"] ** 4
        assert params["luminosity_w"] == pytest.approx(expected, rel=1e-12)
        # the sequence's own luminosity is within the rounding of the temperature (50 K in the cool end) of it
        assert params["luminosity_w"] == pytest.approx(state["luminosity_sol"] * constants.SOLAR_LUMINOSITY, rel=4.0 * 50.0 / params["temperature_k"] + 1e-3)
