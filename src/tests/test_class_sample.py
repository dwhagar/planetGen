"""The FORCE_CLASS switch and the class sample built on it (GEN.791, generation/class_sample.py)."""

import pytest

from planetgen.generation import class_sample
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


@pytest.mark.parametrize("letter", ["O", "P", "T", "V", "E", "K", "L", "N", "F", "G", "Q"])
def test_a_forced_rare_class_is_in_the_system(letter):
    for seed in range(1, 25):
        system = forced_system(letter, seed)
        if not system.unmet_requirements:
            assert any(getattr(planet, "planet_class", None) == letter for planet in system.planets)
            return
    pytest.fail(f"no system with a class {letter} planet in 24 tries")


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
    config.FORCE_CLASS = "Z"
    with pytest.raises(ValueError, match="not a planet class"), draw.bound(draw.Stream(1)):
        StarSystem(system_config=config)


def test_fill_reaches_the_count_and_is_repeatable():
    letters = ["O", "T", "P"]
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
