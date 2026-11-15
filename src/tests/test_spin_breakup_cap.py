"""A balance of the thermal and the gravitational tide cannot turn a planet past the fraction of its breakup rate no body exceeds (seen as the red
`test_orbital_motion` rotation-period test: a 3.5 Earth mass planet of 9 Earth radii at 6.2 hours, where it would leave the ground at 24)."""

import math

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import planets, spin, spin_torques
from planetgen.util import draw


def _rocky_planet_without_moons():
    for seed in range(1, 40):
        draw.set_run_seed(seed)
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        system = StarSystem(system_config=cfg)
        for planet in system.planets:
            if planet.body_type == "t" and not planet.moons and not getattr(planet, "_born_moons", None):
                return planet
    raise AssertionError("no rocky planet without moons in 39 seeds")


def test_an_asynchronous_balance_faster_than_the_breakup_fraction_is_held_at_it(monkeypatch):
    planet = _rocky_planet_without_moons()
    breakup = spin.BREAKUP_FRACTION * spin.critical_spin_rad_s(planet.mass, planet.radius)
    too_fast = spin_torques.SpinResult(3.0 * breakup, None, "asynchronous")
    monkeypatch.setattr(planets.body_spin, "evolve", lambda *args, **kwargs: too_fast)
    planets.generate_spin(planet, planet.primary_mass_kg, birth=planet._birth_spin)
    assert math.isclose(planet.rotation_period_hours, 2.0 * math.pi / breakup / 3600.0, rel_tol=1e-9)


def test_a_balance_slower_than_the_breakup_fraction_is_left_as_it_is(monkeypatch):
    planet = _rocky_planet_without_moons()
    breakup = spin.BREAKUP_FRACTION * spin.critical_spin_rad_s(planet.mass, planet.radius)
    slow = spin_torques.SpinResult(0.1 * breakup, None, "asynchronous")
    monkeypatch.setattr(planets.body_spin, "evolve", lambda *args, **kwargs: slow)
    planets.generate_spin(planet, planet.primary_mass_kg, birth=planet._birth_spin)
    assert math.isclose(planet.rotation_period_hours, 2.0 * math.pi / (0.1 * breakup) / 3600.0, rel_tol=1e-9)
