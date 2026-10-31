# tests/test_remnant_belts.py

"""GEN.1072: a belt through the giant branch and the mass loss of a white dwarf's progenitor (`generation/remnant_orbits.py`, `physics/belt_history.py`)."""

from types import SimpleNamespace

import pytest

from planetgen.generation import remnant_orbits, validation
from planetgen.physics import belt_history, constants

SOLAR_KG = constants.SOLAR_MASS_TO_KG


def _belt(lower=2.0, upper=6.0, mass_kg=1.0e21):
    return SimpleNamespace(body_type='a', distance=0.5 * (lower + upper), lower_limit=lower, upper_limit=upper, mass_kg=mass_kg)


def _white_dwarf(initial_mass_sol=2.0, mass_sol=0.6):
    return SimpleNamespace(yerkes_class="VII", mass=mass_sol * SOLAR_KG, initial_mass_sol=initial_mass_sol)


def test_a_belt_the_envelope_reaches_into_keeps_the_outer_part_and_its_share_of_the_mass():
    belt = _belt()
    remnant_orbits.engulf_belt(belt, 3.0)
    assert (belt.lower_limit, belt.upper_limit) == (3.0, 6.0)
    share = belt_history.BeltState.formed(2.0, 6.0)
    expected = belt_history.engulf(share, 3.0).mass_fraction
    assert belt.mass_kg == pytest.approx(1.0e21 * expected) and 0.0 < expected < 1.0
    assert belt.lower_limit <= belt.distance <= belt.upper_limit
    assert len(belt.history) == 1 and "swelled out to 3.00 au" in belt.history[0]


def test_a_belt_wholly_inside_or_wholly_outside_the_envelope_is_left_for_the_caller():
    inside, outside = _belt(1.0, 2.0), _belt(4.0, 6.0)
    remnant_orbits.engulf_belt(inside, 3.0)
    remnant_orbits.engulf_belt(outside, 3.0)
    assert (inside.lower_limit, inside.upper_limit, inside.mass_kg) == (1.0, 2.0, 1.0e21) and not hasattr(inside, "history")
    assert (outside.lower_limit, outside.upper_limit, outside.mass_kg) == (4.0, 6.0, 1.0e21) and not hasattr(outside, "history")


def _widen(planets, monkeypatch):
    monkeypatch.setattr(validation, "reconcile_moved_planet", lambda body: None)
    return remnant_orbits.widen_for_mass_loss(planets, _white_dwarf())


def test_the_mass_loss_widens_a_belt_and_with_a_giant_left_throws_a_share_of_it_inward(monkeypatch):
    giant = SimpleNamespace(body_type='g', distance=20.0, mass=1.9e27)
    belt = _belt(2.0, 6.0)
    assert _widen([belt, giant], monkeypatch)
    factor = _white_dwarf().mass and (remnant_orbits.formation_mass_kg(_white_dwarf()) / _white_dwarf().mass)
    assert belt.lower_limit == pytest.approx(2.0 * factor) and belt.upper_limit == pytest.approx(6.0 * factor)
    assert belt.lower_limit <= belt.distance <= belt.upper_limit
    thrown = belt_history.MASS_LOSS_SCATTER_SHARE * (1.0 - 1.0 / factor)
    assert belt.mass_kg == pytest.approx(1.0e21 * (1.0 - thrown)) and belt.scattered_mass_kg == pytest.approx(1.0e21 * thrown)
    assert "thrown inward" in belt.history[-1]
    rate = remnant_orbits.pollution_rate_g_s([belt, giant])
    assert rate == pytest.approx(belt_history.pollution_rate_g_s(belt.scattered_mass_kg)) and belt_history.is_polluted(rate)


def test_a_belt_with_no_giant_widens_with_nothing_thrown_in(monkeypatch):
    belt = _belt(2.0, 6.0)
    assert _widen([belt], monkeypatch)
    assert belt.scattered_mass_kg == 0.0 and belt.mass_kg == pytest.approx(1.0e21)
    assert remnant_orbits.pollution_rate_g_s([belt]) == 0.0


def test_a_belt_keeps_one_history_through_the_swelling_and_the_mass_loss(monkeypatch):
    giant = SimpleNamespace(body_type='g', distance=20.0, mass=1.9e27)
    belt = _belt(2.0, 6.0)
    remnant_orbits.engulf_belt(belt, 3.0)
    after_swelling = belt.mass_kg
    assert _widen([belt, giant], monkeypatch)
    assert len(belt.history) == 2 and "swelled" in belt.history[0] and "lost" in belt.history[1]
    assert belt.mass_kg < after_swelling and belt.formed_mass_kg == 1.0e21    # the mass it formed with is kept, not the mass left


def test_a_star_that_lost_no_mass_leaves_its_belt_alone(monkeypatch):
    monkeypatch.setattr(validation, "reconcile_moved_planet", lambda body: None)
    belt = _belt()
    star = SimpleNamespace(yerkes_class="V", mass=SOLAR_KG)
    assert not remnant_orbits.widen_for_mass_loss([belt], star)
    assert not hasattr(belt, "history") and belt.mass_kg == 1.0e21


def _white_dwarf_system(seed):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    config = SystemConfig()
    config.STAR_TYPE = "G2VII"
    config.PLANETS = True
    config.BINARY_SYSTEM = False
    with draw.bound(seed):
        return StarSystem(system_config=config)


def test_a_white_dwarfs_belt_has_been_through_the_swelling_and_the_mass_loss_of_its_progenitor():
    for seed in range(300):
        system = _white_dwarf_system(seed)
        belts = [body for body in system.planets if body.body_type == 'a']
        if belts:
            break
    else:
        pytest.fail("no white dwarf system with a belt in 300 seeds")
    from planetgen.generation.system import StarSystem
    engulfed_au = StarSystem._engulfment_radius_au(system.star)
    for belt in belts:
        assert belt.lower_limit >= engulfed_au - 1.0e-9 and belt.lower_limit <= belt.distance <= belt.upper_limit
        assert 0.0 < belt.mass_kg <= belt.formed_mass_kg
        assert any("lost" in line and "mass" in line for line in belt.history)    # the mass-loss line every white dwarf's belt has
        assert belt.scattered_mass_kg >= 0.0
