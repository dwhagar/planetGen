"""The FORCE_CLASS switch and the class sample built on it (GEN.791, generation/class_sample.py)."""

import pytest

from planetgen.generation import class_sample, condition_request
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


def forced_system(letter, seed, star_type=None):
    config = SystemConfig()
    config.FORCE_CLASS = letter
    if star_type:
        config.STAR_TYPE = star_type
    with draw.bound(draw.Stream(seed)):
        return StarSystem(system_config=config)


def test_the_switch_is_off_by_default():
    assert SystemConfig().FORCE_CLASS is None
    assert "FORCE_CLASS" not in SystemConfig().to_dict()


def _meets(letter, system):
    request = condition_request.from_letter(letter)
    return any(condition_request.meets(request, planet) for planet in system.planets)


@pytest.mark.parametrize("letter", ["E", "N"])
def test_a_forced_class_the_free_generator_gives_is_in_the_system(letter):
    # GEN.688: a request is judged on the finished planets' conditions, and the draws are the generator's own, so only a class a free draw reaches at some place is met.
    for seed in range(1, 25):
        system = forced_system(letter, seed)
        if not system.unmet_requirements:
            assert _meets(letter, system)
            return
    pytest.fail(f"no system with a planet that meets class {letter} in 24 tries")


@pytest.mark.parametrize("letter", ["P", "M", "H", "O", "R", "V", "G", "T", "K", "L", "F"])
def test_a_class_the_conditions_hardly_give_is_met_or_reported_unmet(letter):
    # GEN.646, GEN.1333 (Boss, 2026-11-08): the class is derived from the finished body and no draw is steered to it, so a class the free draw seldom makes (P needs oxygen,
    # which only oxygenic life gives, GEN.1339 and GEN.1340; V a dense bare core) stays requestable and the system says "unmet" and never claims what it does not hold.
    for seed in range(1, 4):
        system = forced_system(letter, seed)
        assert _meets(letter, system) != bool(system.unmet_requirements)
        if system.unmet_requirements:
            assert system.unmet_requirements == [f"a planet of class {letter}"]


def test_a_request_for_conditions_is_met_by_a_planet_that_meets_them_or_reported_unmet():
    for seed in range(1, 4):
        config = SystemConfig()
        config.FORCE_CONDITIONS = {"kind": "g", "mass_earth": [35.0, None]}
        with draw.bound(draw.Stream(seed)):
            system = StarSystem(system_config=config)
        request = condition_request.from_conditions(config.FORCE_CONDITIONS)
        met = any(condition_request.meets(request, planet) for planet in system.planets)
        assert met != bool(system.unmet_requirements)


def test_a_forced_request_never_draws_the_whole_system_again():
    # GEN.688: the class part of the system retry is gone: a request is judged once on the finished planets.
    with draw.counting() as counts:
        forced_system("P", 3)
    assert counts.get("system.planets_regenerated", 0) == 0 and counts.get("system.planets_regenerated_after_age", 0) == 0


def test_a_class_the_star_cannot_give_is_reported_unmet():
    system = forced_system("M", 1, star_type="O5V")
    assert system.unmet_requirements == ["a planet of class M"]


def test_the_switch_refuses_slots_and_unknown_letters():
    config = SystemConfig()
    config.FORCE_CLASS = "M"
    config.SLOTS = [{"type": "planet"}]
    with pytest.raises(ValueError, match="SLOTS"), draw.bound(draw.Stream(1)):
        StarSystem(system_config=config)
    config = SystemConfig()
    config.FORCE_CLASS = "Q"
    with pytest.raises(ValueError, match="not a planet class"), draw.bound(draw.Stream(1)):
        StarSystem(system_config=config)


def test_fill_reaches_the_count_and_is_repeatable():
    letters = ["E", "N"]
    systems, counts, short = class_sample.fill(per_class=2, seed=5, letters=letters)
    assert not short
    assert all(counts[letter] >= 2 for letter in letters)
    again, again_counts, _short = class_sample.fill(per_class=2, seed=5, letters=letters)
    assert [system.name for system in systems] == [system.name for system in again]
    assert counts == again_counts


def test_fill_counts_the_systems_it_is_given_and_reports_a_class_it_cannot_reach():
    systems, counts, _short = class_sample.fill(per_class=1, seed=2, letters=["J"])
    more, more_counts, short = class_sample.fill(per_class=counts["J"], seed=9, letters=["J"], systems=systems)
    assert len(more) == len(systems) and not short
    _none, _counts, short = class_sample.fill(per_class=1, seed=1, letters=["O"], tries_per_body=0)
    assert short == ["O"]
