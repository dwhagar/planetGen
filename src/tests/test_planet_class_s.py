"""
Class S, the barren rocky super-Earth (GEN.38, GEN.28's first class): a rocky body of at least 1.5 Earth masses with almost no air
(GEN.645 rule), never a habitable class, and on the class reference page.
"""

import pytest

from planetgen.physics import constants as pc
from planetgen.physics import planet_classes
from planetgen import tuning

from planetgen.web.lib import classref  # noqa: E402


def test_class_s_is_a_rocky_class_with_no_life():
    assert tuning.PLANET_CLASSES["S"]["type"] == "t" and "life_class" not in tuning.PLANET_CLASSES["S"]
    assert "S" not in tuning.HABITABLE_PLANET_CLASSES


def test_the_rule_gives_s_to_a_heavy_airless_rock():
    body = {"body_type": "t", "mass_kg": 4 * pc.EARTH_MASS_TO_KG, "gravity_g": 1.5, "surface_temperature_k": 300.0,
            "atmospheric_pressure_pa": 0.0, "water_mass_fraction": 1e-6}
    assert planet_classes.classify(body) == "S"
    assert planet_classes.classify({**body, "mass_kg": 0.4 * pc.EARTH_MASS_TO_KG, "gravity_g": 0.5}) != "S"


def test_class_s_is_on_the_class_reference_page():
    cat = classref.catalog()
    planets = next(kind for kind in cat.values() if "S" in kind.get("classes", {}) and "V" in kind["classes"])
    assert "super-earth" in str(planets["classes"]["S"]).lower()
