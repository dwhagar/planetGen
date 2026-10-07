"""
The star population model (`planetgen.physics.stellar_evolution`, used by
`Star.generate_star` for every random star).

The shares are checked against the real census the model is built to
reproduce (see docs/design and the "Star population model" section of
program_constants): most stars are M dwarfs, about 5% are white dwarfs,
giants are a fraction of a percent and supergiants vanishingly rare.

Run with: pytest src/tests/test_star_population.py
"""
import collections
import math
import random

import pytest

from planetgen.physics import constants as phys_c
from planetgen import tuning as prog_c
from planetgen.physics import stellar_evolution as se
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import STAR_TYPE_PATTERN, Star
from planetgen.generation.system import StarSystem

SAMPLE_SIZE = 200_000


def _category(state, rng_letter):
    yerkes = state["yerkes_class"]
    if yerkes == "V":
        return f"{rng_letter} dwarf"
    if yerkes == "IV":
        return "subgiant"
    if yerkes in ("III", "II"):
        return "giant"
    if yerkes == "VII":
        return "white dwarf"
    return "supergiant"


@pytest.fixture(scope="module")
def population():
    rng = random.Random(20260930)
    counts = collections.Counter()
    states = []
    for _ in range(SAMPLE_SIZE):
        mass, age, state = se.sample_living_star(rng=rng)
        letter, _ = se.spectral_letter_and_subclass(state["temperature_k"])
        counts[_category(state, letter)] += 1
        states.append((mass, age, state))
    return counts, states


# (category, low %, high %): the model's expected share with room for
# sampling noise and later tuning of the letter boundaries.
EXPECTED_SHARES = [
    ("M dwarf", 68.0, 80.0),
    ("K dwarf", 11.0, 18.0),
    ("G dwarf", 2.0, 6.0),
    ("F dwarf", 1.2, 3.0),
    ("A dwarf", 0.4, 1.2),
    ("B dwarf", 0.1, 0.4),
    ("white dwarf", 4.0, 7.0),
    ("subgiant", 0.15, 0.6),
    ("giant", 0.12, 0.5),
]


@pytest.mark.parametrize("category, low, high", EXPECTED_SHARES)
def test_star_shares_match_the_census(population, category, low, high):
    counts, _ = population
    share = 100 * counts[category] / SAMPLE_SIZE
    assert low <= share <= high, f"{category}: {share:.3f}% (want {low}-{high}%)"


def test_supergiants_are_vanishingly_rare(population):
    counts, _ = population
    # Real share ~1e-6; allow a few in 200,000 for noise.
    assert counts["supergiant"] <= 5


def test_main_sequence_luminosity_follows_mass(population):
    _, states = population
    for mass, _, state in states:
        if state["yerkes_class"] == "V":
            assert state["luminosity_sol"] == pytest.approx(se.main_sequence_luminosity_sol(mass))


def test_every_state_is_in_its_own_phase(population):
    _, states = population
    for mass, age, state in states:
        t_ms = se.main_sequence_lifetime_gy(mass)
        if state["yerkes_class"] == "V":
            assert age < t_ms
        else:
            assert age >= t_ms
        assert age < state["phase_end_gy"]


def test_white_dwarfs_are_earth_sized_and_dim(population):
    _, states = population
    white_dwarfs = [state for _, _, state in states if state["yerkes_class"] == "VII"]
    assert white_dwarfs
    for state in white_dwarfs:
        low, high = prog_c.WD_MASS_RANGE_SOL
        assert low <= state["mass_sol"] <= high
        assert 3000 < state["radius_km"] < 10000
        assert state["luminosity_sol"] <= prog_c.WD_LUMINOSITY_RANGE_SOL[1]


def test_imf_sampler_respects_truncation():
    rng = random.Random(1)
    masses = [se.sample_imf_mass_sol(min_mass_sol=1.4, rng=rng) for _ in range(5000)]
    assert min(masses) >= 1.4
    assert max(masses) <= prog_c.IMF_BREAKS_SOL[-1]
    # Kroupa: about half of all stars are below 0.3 Msun.
    all_masses = [se.sample_imf_mass_sol(rng=rng) for _ in range(20000)]
    below = sum(m < 0.3 for m in all_masses) / len(all_masses)
    assert 0.4 < below < 0.65


def test_collapsed_massive_stars_return_none():
    assert se.evolve_star(20.0, 5.0) is None


@pytest.mark.parametrize("letter_temperature", [(35000, "O"), (5772, "G"), (3000, "M"), (100000, "O"), (1000, "M")])
def test_spectral_letter_from_temperature(letter_temperature):
    temperature, letter = letter_temperature
    assert se.spectral_letter_and_subclass(temperature)[0] == letter


def test_random_stars_are_self_consistent():
    for _ in range(300):
        star = Star(SystemConfig())
        assert STAR_TYPE_PATTERN.fullmatch(star.type.split()[0])
        assert star.initial_mass_sol is not None
        assert 0 <= star.age < star.phase_end_age_gy
        if star.yerkes_class == "V":
            mass_sol = star.mass / phys_c.SOLAR_MASS_TO_KG
            luminosity_sol = star.luminosity / phys_c.SOLAR_LUMINOSITY
            assert luminosity_sol == pytest.approx(se.main_sequence_luminosity_sol(mass_sol))


def test_large_stars_are_massive_and_alive():
    cfg = SystemConfig()
    cfg.LARGE_STAR = True
    for _ in range(200):
        star = Star(cfg)
        assert star.initial_mass_sol >= prog_c.LARGE_STAR_MIN_MASS_SOL
        assert star.yerkes_class != "VII"


def test_binary_companions_share_the_primarys_age_and_are_dimmer():
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = True
    cfg.PLANETS = False
    for _ in range(150):
        system = StarSystem(system_config=cfg)
        primary, secondary = system.primary_star, system.secondary_star
        assert secondary.mass <= primary.mass
        if primary.yerkes_class == "V" and secondary.yerkes_class == "V":
            assert secondary.luminosity <= primary.luminosity * (1 + 1e-9)
        assert math.isclose(primary.age, secondary.age, rel_tol=1e-9)


def test_planets_respect_their_stars_age_and_history():
    habitable = set(prog_c.HABITABLE_PLANET_CLASSES)
    for _ in range(400):
        system = StarSystem(system_config=SystemConfig())
        hosts = [(system.star, system.planets)]
        if system.binary_type == "wide":
            hosts.append((system.secondary_star, system.secondary_planets))
        for star, bodies in hosts:
            engulfed_au = system._engulfment_radius_au(star)
            for body in bodies:
                inner_edge = body.lower_limit if body.body_type == "a" else body.distance
                assert inner_edge >= engulfed_au
                if star.age < prog_c.PLANET_MIN_STAR_AGE_GY:
                    assert body.body_type == "a"
                if star.age < prog_c.LIFE_MIN_STAR_AGE_GY and body.body_type != "a":
                    assert body.planet_class not in habitable
                    assert all(moon.planet_class not in habitable for moon in body.moons)


def test_a_required_habitable_world_gets_an_old_enough_living_star():
    cfg = SystemConfig()
    cfg.HABITABLE_WORLD = True
    for _ in range(100):
        star = Star(cfg)
        assert star.age >= prog_c.LIFE_MIN_STAR_AGE_GY
        assert star.yerkes_class != "VII"


def test_from_params_rebuilds_a_star_without_rolling():
    rng = random.Random(7)
    mass, age, state = se.sample_living_star(rng=rng)
    params = se.star_params(mass, age, state)
    star = Star.from_params(params, SystemConfig(), name="Beacon")
    assert star.name == "Beacon"
    assert (star.type, star.mass, star.radius, star.temperature, star.luminosity, star.age) == (
        params["type"], params["mass_kg"], params["radius_km"], params["temperature_k"],
        params["luminosity_w"], params["age_gy"])
    assert star.habitable_zone and star.system_perimeter > 0

    system = StarSystem(system_config=SystemConfig(), primary_star_params=params)
    assert system.primary_star.mass == params["mass_kg"] or getattr(system, "secondary_star", system.primary_star).mass == params["mass_kg"]
