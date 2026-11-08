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


@pytest.mark.parametrize("letter", ["E", "V", "N"])
def test_a_forced_rare_class_is_in_the_system(letter):
    for seed in range(1, 25):
        system = forced_system(letter, seed)
        if not system.unmet_requirements:
            assert any(getattr(planet, "planet_class", None) == letter for planet in system.planets)
            return
    pytest.fail(f"no system with a class {letter} planet in 24 tries")


def test_a_forced_carbon_dioxide_world_is_met_at_a_sun_like_star():
    """GEN.1333: a request for class G draws its masses over the gravity G holds (0.5 to 2 g) and draws 80 planets at a place; at a Sun-like star it is met in
    18 of 24 seeds when measured (8 seeds here, half of them asked)."""
    met = 0
    for seed in range(1, 9):
        system = forced_system("G", seed, star_type="G2V")
        met += (not system.unmet_requirements) and any(getattr(planet, "planet_class", None) == "G" for planet in system.planets)
    assert met >= 4


def test_a_forced_glaciated_world_needs_the_oxygen_that_life_builds():
    """GEN.1333: class P is 16 to 35.5 kPa of oxygen at 240 to 270 K, and the draws of a request have no oxygen (only oxygenic life builds it, after the rise;
    GEN.1339 and GEN.1340 are about how often it does). The request never claims a P it does not hold: it is met or reported unmet. It is not asked to be met
    in 18 of 24 seeds until the oxygen chain gives a glaciated oxygen world in reach of a request."""
    for seed in range(1, 3):
        system = forced_system("P", seed, star_type="G2V")
        held = any(getattr(planet, "planet_class", None) == "P" for planet in system.planets)
        assert held != bool(system.unmet_requirements)


def test_a_class_request_draws_its_masses_over_the_gravity_of_its_class():
    from planetgen import tuning
    from planetgen.generation import system as system_module
    from planetgen.physics import class_fit

    assert class_fit.gravity_band("G") == (0.5, 2.0)
    assert class_fit.gravity_band("V")[1] is None and class_fit.gravity_band("J") == (None, None)
    low, high = system_module._request_mass_range_earth("G")
    exponent, widen = tuning.CLASS_REQUEST_GRAVITY_MASS_EXPONENT, tuning.CLASS_REQUEST_GRAVITY_WIDEN
    assert low == pytest.approx(0.5 ** (1.0 / exponent) / widen) and high == pytest.approx(2.0 ** (1.0 / exponent) * widen)
    assert 0.1 < low < 0.15 and 8.0 < high < 10.0
    assert system_module._request_mass_range_earth("J") == tuning.CLASS_REQUEST_MASS_EARTH       # a giant names no gravity: all masses
    assert system_module._request_mass_range_earth("V") == tuning.CLASS_REQUEST_MASS_EARTH       # V has no upper gravity edge


@pytest.mark.parametrize("letter", ["O", "T", "K", "L", "F"])
def test_a_class_the_conditions_hardly_give_is_met_or_reported_unmet(letter):
    # GEN.646: the class is derived from the finished body, so an ocean world (O), a gas dwarf (T: until the envelope law GEN.822), a temperate thin-air world (K), a marginal thick-air world (L) or a reducing world (F) (GEN.1142: the cold ones are E and P now, and
    # a body that keeps under 3 kPa at 273 to 320 K is rare) is not reachable by a request at most places; the system says so and never claims a class it does not hold.
    for seed in range(1, 4):
        system = forced_system(letter, seed)
        held = any(getattr(planet, "planet_class", None) == letter for planet in system.planets)
        assert held != bool(system.unmet_requirements)


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
    letters = ["E", "V", "N"]
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
