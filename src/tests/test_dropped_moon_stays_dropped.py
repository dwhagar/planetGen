"""
GEN.1856: a moon the settle pass drops for its size does not come back (a bare one, with no spin: the writer refused `rotation_period_hours` NULL on a moon).

`planets.drop_unstable_moons` takes a moon out of `planet.moons`, but the history of the moons (`moon_history.apply`, run again by the planet's `refinish_rocky`) is run from the moons the planet was
born with and put back every moon that stayed in its orbit. The orbit rules (the bounds and the stable range) were checked there; the size rule (`planets.moon_size_limits`) was not, so a moon dropped for being
too big for its planet returned as it was when the climate pass owed it its air, class and spin.
"""

from planetgen.generation import validation
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import planets as planetPhysics
from planetgen.util import draw

MOON_SEARCH_SEEDS = 40
"""int: The seeds searched for a system that holds a moon; the first is judged."""


def _system(seed):
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    with draw.bound(seed):
        return StarSystem(system_config=cfg)


def test_a_moon_dropped_for_its_size_is_not_put_back_bare(monkeypatch):
    real_limits = planetPhysics.moon_size_limits
    real_pumping = StarSystem._kozai_pumping
    real_remake = planetPhysics.needs_remake
    state = {"tight": False, "moons": 0}

    def limits(parent):
        radius, mass = real_limits(parent)
        return (radius, mass * 1.0e-9) if state["tight"] else (radius, mass)      # from the settle pass on, no moon is small enough for its planet

    def pumping(self):
        state["tight"] = True
        state["moons"] = sum(len(planet.moons) for planet in self.all_planets if planet.body_type != "a")
        return real_pumping(self)

    def remake(body):
        return True if state["tight"] and not body.is_moon else real_remake(body)      # a planet whose end stages are made again runs the history of its moons again (`refinish_rocky`)

    monkeypatch.setattr(planetPhysics, "moon_size_limits", limits)
    monkeypatch.setattr(planetPhysics, "needs_remake", remake)
    monkeypatch.setattr(StarSystem, "_kozai_pumping", pumping)
    for seed in range(1, MOON_SEARCH_SEEDS):
        state["tight"], state["moons"] = False, 0
        system = _system(seed)
        state["tight"] = False
        if not state["moons"]:
            continue
        bodies = [(planet, moon) for planet in system.all_planets if planet.body_type != "a" for moon in planet.moons]
        assert not bodies, f"seed {seed}: the moons {[moon.name for _planet, moon in bodies]} that the size rule dropped are back"
        return
    raise AssertionError(f"no system of the first {MOON_SEARCH_SEEDS - 1} seeds holds a moon at the settle pass")


def test_no_moon_is_left_bare_at_the_end_of_the_settle_pass(monkeypatch):
    """Whatever the history gives back (a moon dropped for its orbit may take a new one in the replay), every moon of a planet is finished: none keeps its air, class and spin owed (`planets.is_bare`)."""
    real_remake = planetPhysics.needs_remake
    real_pumping = StarSystem._kozai_pumping
    state = {"on": False, "moons": 0}

    def pumping(self):
        state["on"] = True
        state["moons"] = sum(len(planet.moons) for planet in self.all_planets if planet.body_type != "a")
        return real_pumping(self)

    monkeypatch.setattr(planetPhysics, "needs_remake", lambda body: True if state["on"] and not body.is_moon else real_remake(body))
    monkeypatch.setattr(StarSystem, "_kozai_pumping", pumping)
    for seed in range(1, MOON_SEARCH_SEEDS):
        state["on"], state["moons"] = False, 0
        system = _system(seed)
        state["on"] = False
        if not state["moons"]:
            continue
        bare = [moon.name for planet in system.all_planets if planet.body_type != "a" for moon in planet.moons if planetPhysics.is_bare(moon) or moon.rotation_period_hours is None]
        assert not bare, f"seed {seed}: the moons {bare} are left bare"
        return
    raise AssertionError(f"no system of the first {MOON_SEARCH_SEEDS - 1} seeds holds a moon at the settle pass")
