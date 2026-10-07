"""
A binary pair's two stars share one age, and a `--star-type` primary's
companion is a star whose type and luminosity follow from its own mass
(GEN.53, GEN.54).
"""

import random

import pytest

from planetgen.physics import constants
from planetgen.generation.config import SystemConfig
from planetgen.physics.stellar_evolution import main_sequence_luminosity_sol
from planetgen.generation.system import StarSystem

STAR_TYPES = [None, "G2V", "M2V", "O5V", "B3V", "A0V", "K1III", "M2VII", "M2IA"]
SEEDS = range(12)


def _binary(star_type, seed, **overrides):
    random.seed(seed)
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.BINARY_SYSTEM = True
    cfg.MAX_PLANETS = True
    for attr, value in overrides.items():
        setattr(cfg, attr, value)
    return StarSystem(system_config=cfg)


@pytest.mark.parametrize("wide", [False, True])
@pytest.mark.parametrize("star_type", STAR_TYPES)
def test_both_stars_of_a_binary_share_one_age(star_type, wide):
    for seed in SEEDS:
        system = _binary(star_type, seed, WIDE_BINARY=wide)
        primary, secondary = system.primary_star, system.secondary_star
        assert primary.age == secondary.age, (star_type, seed, system.binary_type)
        if system.binary_type == "close":
            assert system.star.age == primary.age
        # The shared age never pushes a population-model star out of the
        # phase it was generated in.
        for star in (primary, secondary):
            if star.phase_end_age_gy is not None and star.phase_end_age_gy != float("inf"):
                assert star.age <= star.phase_end_age_gy * 1.0001, (star_type, seed, star.type)


@pytest.mark.parametrize("star_type", ["G2V", "M2V", "O5V", "B3V", "A0V", "K1III", "M2IA"])
def test_a_specified_type_secondary_follows_from_its_own_mass(star_type):
    for seed in SEEDS:
        system = _binary(star_type, seed)
        secondary = system.secondary_star
        assert secondary.initial_mass_sol is not None, "companion comes from the population model"
        assert secondary.mass <= system.primary_star.mass
        if secondary.yerkes_class == "V":
            mass_sol = secondary.mass / constants.SOLAR_MASS_TO_KG
            lum_sol = secondary.luminosity / constants.SOLAR_LUMINOSITY
            # A main-sequence companion sits on (a little above, as it ages)
            # the mass-luminosity relation, not wherever its type's range
            # happened to put it.
            ratio = lum_sol / main_sequence_luminosity_sol(mass_sol)
            assert 0.5 < ratio < 3.0, (star_type, seed, secondary.type, mass_sol, lum_sol)
