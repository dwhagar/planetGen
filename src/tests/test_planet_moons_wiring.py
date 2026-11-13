# tests/test_planet_moons_wiring.py

"""GEN.1157: a giant holds a captured moon of planet size, made as a planet is, with a captured orbit; the disk's moons stay under the disk's ceiling (planet_moons wired into generate_moons)."""

import pytest

from planetgen.physics import planet_kind
from planetgen.physics import moon_composition as mc
from planetgen.physics import planet_moons as pm
from planetgen.util import draw


def _giants(seed, systems, monkeypatch=None):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    found = []
    with draw.bound(seed):
        for _ in range(systems):
            config = SystemConfig()
            config.MOONS = True
            system = StarSystem(system_config=config)
            found += [planet for planet in system.planets + system.secondary_planets if planet.body_type == "g" and planet.moons]
    return found


def test_the_giant_cap_is_the_disk_ceiling():
    assert planet_kind.moon_mass_ceiling_kg(1e27, True) == pytest.approx(pm.DISK_TOTAL_RATIO_MAX * 1e27)
    assert "exchange" in mc.ORIGINS and "exchange" in mc.CAPTURE_ORBIT_ORIGINS


@pytest.fixture(scope="module")
def captured_giants():
    """Giants with moons of 25 seeded systems, at a capture chance of 1."""
    saved = pm.CAPTURE_CHANCE
    pm.CAPTURE_CHANCE = 1.0
    try:
        return _giants(31, 25)
    finally:
        pm.CAPTURE_CHANCE = saved


def test_a_capture_moon_is_a_planet_sized_moon_with_a_captured_orbit(captured_giants):
    low, high = mc.KOZAI_FORBIDDEN_DEG
    seen = [(giant, moon) for giant in captured_giants for moon in giant.moons if moon.moon_origin == "exchange"]
    assert seen, "no giant took a capture moon in 25 systems at a chance of 1"
    for giant, moon in seen:
        assert pm.PLANET_MOON_MIN_MASS_KG <= moon.mass <= pm.exchange_ceiling_kg(giant.mass)
        assert not low < moon.orbital_inclination_deg < high
        assert moon.radius > 0 and moon.density > 0
    for giant in captured_giants:
        for moon in giant.moons:
            if moon.moon_origin == "disk":
                assert moon.mass <= pm.DISK_TOTAL_RATIO_MAX * giant.mass * 1.0001


def test_no_capture_moon_at_a_chance_of_zero(monkeypatch):
    monkeypatch.setattr(pm, "CAPTURE_CHANCE", 0.0)
    assert not any(moon.moon_origin == "exchange" for giant in _giants(31, 10) for moon in giant.moons)
