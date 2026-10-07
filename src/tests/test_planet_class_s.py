"""
Class S, the barren rocky super-Earth (GEN.38, GEN.28's first class):
1.2-1.8 Earth radii, 2-10 Earth masses, every zone, rogue-eligible, never
a moon, and on the class reference page.
"""

import pytest

from planetgen.physics import constants as pc
from planetgen.generation import plausibility
from planetgen.physics import planets as planetPhysics
from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star

import web  # noqa: F401 -- puts src/html/lib on sys.path
import classref  # noqa: E402


@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


def test_class_s_ranges():
    data = tuning.PLANET_CLASSES["S"]
    low, high = data["radius_range"]
    assert low / pc.EARTH_RADIUS_KM == pytest.approx(1.2, rel=0.01)
    assert high / pc.EARTH_RADIUS_KM == pytest.approx(1.8, rel=0.01)
    low_kg, high_kg = planetPhysics.planet_mass_ranges["S"]
    assert low_kg / pc.EARTH_MASS_TO_KG == pytest.approx(2.0, rel=0.01)
    assert high_kg / pc.EARTH_MASS_TO_KG == pytest.approx(10.0, rel=0.01)
    assert data["h"] and data["e"] and data["c"] and data["r"]
    assert data["type"] == "t" and data["life_chemical"] is None
    assert "S" not in tuning.HABITABLE_PLANET_CLASSES
    assert "S" in tuning.MOON_BLACKLIST
    assert tuning.PLANET_CLASS_PROBABILITIES["S"] > 0


@pytest.mark.parametrize("zone", ["h", "e", "c"])
def test_class_s_planets_fit_the_class(host_star, zone):
    low, high = tuning.PLANET_CLASSES["S"]["radius_range"]
    low_kg, high_kg = planetPhysics.planet_mass_ranges["S"]
    for _ in range(50):
        planet = Planet(SystemConfig(), host_star, host_star.habitable_zone,
                        plausibility.distance_for_zone(host_star, zone),
                        planet_class="S", zone_override=zone, moon_count=0)
        assert planet.planet_class == "S"
        assert low <= planet.radius <= high
        assert low_kg * 0.999 <= planet.mass <= high_kg * 1.001


def test_class_s_is_on_the_class_reference_page():
    cat = classref.catalog()
    planets = next(kind for kind in cat.values() if "S" in kind.get("classes", {}) and "V" in kind["classes"])
    assert "super-earth" in str(planets["classes"]["S"]).lower()
