# tests/test_star_evolved_at_final_age.py

"""
GEN.1545: the star is the star of its final age. `Star.adjust_age_for_planets` raises the age after the planets are in; the star's light, temperature, radius, class, X-rays and
flares were worked out at the age it was born with, so a raised star lit its planets as a younger one. `Star.evolve_to_age` evolves it again at the final age (the same initial mass
and [Fe/H], the same tracks), `StarSystem._evolve_stars_to_final_age` does it for every star of the system, and the settle pass makes the planets again on the new light.
"""

import pytest

from planetgen.generation import star as star_module
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics.orbits import calculate_habitable_zone
from planetgen.physics.stellar_evolution import evolve_star, star_params
from planetgen.util import draw

SEEDS = range(0, 40)


@pytest.fixture(scope="module")
def systems():
    made, again = [], []
    real = Star.evolve_to_age

    def counted(self, age_gy):
        before = self.luminosity
        done = real(self, age_gy)
        if done:
            again.append((before, self.luminosity))
        return done

    Star.evolve_to_age = counted
    try:
        for seed in SEEDS:
            with draw.bound(seed):
                made.append(StarSystem(system_config=SystemConfig()))
    finally:
        Star.evolve_to_age = real
    return made, again


def _fresh(star):
    return evolve_star(star.initial_mass_sol, star.age, fe_h=star.metallicity_fe_h, population=getattr(star.system_config, "POPULATION", None))


def test_some_stars_of_the_sample_are_evolved_again(systems):
    assert len(systems[1]) >= 3


def test_some_of_them_moved_in_light(systems):
    assert any(abs(after / before - 1.0) > 0.05 for before, after in systems[1])


def test_a_main_sequence_star_has_the_light_of_its_final_age(systems):
    checked = 0
    for system in systems[0]:
        for star in system.stars:
            if star.yerkes_class != "V" or star.forced_type:
                continue
            state = _fresh(star)
            params = star_params(star.initial_mass_sol, star.age, state)
            assert params["luminosity_w"] == pytest.approx(star.luminosity, rel=1e-9)
            assert params["temperature_k"] == pytest.approx(star.temperature, rel=1e-9)
            assert params["type"] == star.type
            assert star.evolved_at_gy == star.age
            checked += 1
    assert checked >= 30


def test_the_habitable_zone_follows_the_light(systems):
    for system in systems[0]:
        for star in system.stars:
            if star.yerkes_class == "V" and not star.forced_type:
                assert tuple(star.habitable_zone) == pytest.approx(tuple(calculate_habitable_zone(star.luminosity, star.temperature)))


def test_the_bodies_of_a_star_hold_its_habitable_zone(systems):
    for system in systems[0]:
        if system.binary_type is not None:
            continue
        for planet in system.planets:
            if planet.body_type != "a":
                assert tuple(planet.habitable_zone) == pytest.approx(tuple(system.star.habitable_zone))


def test_a_close_pair_proxy_is_the_sum_of_its_stars_as_they_are(systems):
    close = [system for system in systems[0] if system.binary_type == "close"]
    for system in close:
        proxy = system.star
        assert proxy.luminosity == pytest.approx(sum(star.luminosity for star in proxy.stars))
        assert proxy.mass == pytest.approx(sum(star.mass for star in proxy.stars))
        assert proxy.age == proxy.stars[0].age


def _young_main_sequence_star():
    for seed in range(40):
        with draw.bound(seed):
            star = Star(SystemConfig())
        if star.yerkes_class == "V" and star.phase_end_age_gy and star.phase_end_age_gy > 1.5 * star.age and star.age > 0.3:
            return star
    pytest.fail("no young main-sequence star among the seeds")


def test_evolving_again_moves_the_light_with_the_age_and_draws_nothing_from_the_run():
    star = _young_main_sequence_star()
    before = star.luminosity
    age = star.age + 0.5 * (star.phase_end_age_gy - star.age)
    star.age = age      # what `adjust_age_for_planets` does: the age moves and the state is left behind
    assert star.evolved_at_gy != age
    with draw.bound(7):
        state = draw.getstate()
        assert star.evolve_to_age(age) is True
        assert draw.getstate() == state
    assert star.evolved_at_gy == age
    assert star.luminosity > before     # a main-sequence star brightens as it ages
    assert star.evolve_to_age(age) is False


def test_the_spin_and_the_activity_follow_the_age():
    star = _young_main_sequence_star()
    spin_before, xuv_before = star.rotation_period_hours, star.xuv_fluence_j
    age = star.age + 0.5 * (star.phase_end_age_gy - star.age)
    star.age = age
    star.evolve_to_age(age)
    assert star.rotation_period_hours != spin_before
    assert star.xuv_fluence_j != xuv_before


def test_a_star_whose_type_was_asked_for_is_left_as_asked():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    with draw.bound(3):
        star = Star(config)
    assert star.forced_type is True
    before = (star.luminosity, star.type, star.temperature)
    star.age = star.age + 1.0
    assert star.evolve_to_age(star.age) is False
    assert (star.luminosity, star.type, star.temperature) == before
