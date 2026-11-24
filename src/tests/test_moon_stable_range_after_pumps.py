"""
A moon a pump or a resonance leaves outside its planet's stable range is lost (found on the Train 61 gate: `test_wide_binaries_never_overlap` failed with a moon at 1.3e4 km and a stable
range of 6765 to 0 km). The stable range is the Domingos limit of the planet's Hill sphere (`clearance_ledger.moon_bounds_km`), and the planet's eccentricity and the moon's own both lower it
(to 0 for a planet of 0.8 and a moon of 0.9). Generation keeps a moon inside it when the moon is built, and the secular exchange only keeps an eccentricity that holds it; the Kozai pump of
the moons (`kozai_bodies.pump_moons`), the resonances and the pump of the planets set eccentricities nothing checked, so the settle pass now drops what the tide has not worn back into range.
"""

from planetgen.generation import validation
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import kozai_bodies, planets as planetPhysics
from planetgen.util import draw


def _system(seed):
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    with draw.bound(seed):
        return StarSystem(system_config=cfg)


def test_a_moon_whose_pumped_orbit_leaves_the_stable_range_is_lost(monkeypatch):
    pumped = []

    def pump(planet, ledger=None):
        """A pump no stable range constrains: a planet at 0.8 and each of its moons at 0.9, where the Domingos limit is 0 (the tide is held off below)."""
        planet.eccentricity = 0.8
        for moon in planet.moons:
            moon.eccentricity = 0.9
            pumped.append(moon)
        return []

    monkeypatch.setattr(kozai_bodies, "pump_moons", pump)
    monkeypatch.setattr(planetPhysics, "settle_orbit", lambda body: False)
    for seed in range(1, 7):
        system = _system(seed)
        assert not [problem for problem in validation.check_star_system(system) if "stable range" in problem.message], seed      # the planets at 0.8 may crowd each other; only the moons are judged
        kept = {id(moon) for planet in system.all_planets if planet.body_type != "a" for moon in planet.moons}
        assert not any(id(moon) in kept for moon in pumped), seed      # a planet made late (a belt member promoted) never met the pump and keeps its own moons
    assert pumped, "six systems of planets hold a moon somewhere"
