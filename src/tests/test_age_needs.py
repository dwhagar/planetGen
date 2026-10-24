"""The age a star must have to host a body, from the body's conditions (GEN.686, physics/age_needs.py)."""

import copy

import pytest

from planetgen import tuning
from planetgen.physics import age_needs, constants


def body(mass_earth=1.0, gravity=1.0, temperature=288.0, pressure=101325.0, body_type="t", **gases):
    values = {"body_type": body_type, "mass_kg": mass_earth * constants.EARTH_MASS_TO_KG, "gravity_g": gravity,
              "surface_temperature_k": temperature, "atmospheric_pressure_pa": pressure}
    values.update({f"p_{gas}_kpa": share * pressure / 1000.0 for gas, share in gases.items()})
    return values


CASES = {
    "giant": (body(318.0, 2.4, 165.0, 1e5, "g"), "none"),
    "airless": (body(0.05, 0.3, 300.0, 0.0), "none"),
    "small": (body(0.04, 0.4, 300.0, 5000.0), "none"),
    "dense hot": (body(0.9, 0.9, 735.0, 9e6), "molten"),
    "reducing": (body(1.0, 1.0, 320.0, 80000.0, ch4=0.2, h2o=0.1), "molten"),
    "cold": (body(1.0, 1.0, 190.0, 60000.0), "frozen"),
    "heavy": (body(5.0, 1.6, 290.0, 150000.0), "settled"),
    "cold desert": (body(0.4, 0.55, 250.0, 1500.0), "settled"),
    "thin oxygen": (body(0.6, 0.7, 285.0, 2000.0, o2=0.2), "watered"),
    "carbon dioxide": (body(1.0, 1.0, 300.0, 90000.0, co2=0.5), "settled"),
    "temperate": (body(1.0, 1.0, 288.0, 101325.0, n2=0.78, o2=0.21), "watered"),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_each_rule_names_its_tier(name):
    values, expected = CASES[name]
    assert age_needs.tier_of_values(values) == expected


def test_rule_file_ends_with_a_catch_all_and_names_known_tiers():
    assert not age_needs.ROWS[-1]["when"]
    assert all(row["when"] for row in age_needs.ROWS[:-1])
    assert {row["tier"] for row in age_needs.ROWS} <= set(age_needs.TIERS)
    for ages in age_needs.TIERS.values():
        assert all(ages[scale] >= 0.0 for scale in age_needs.SCALES)


def test_slower_stars_need_older_ages_and_tiers_never_get_younger_with_settling():
    for ages in age_needs.TIERS.values():
        assert ages["fast"] <= ages["normal"] <= ages["slow"]
    order = ["none", "molten", "settled", "watered", "frozen"]
    for scale in age_needs.SCALES:
        values = [age_needs.TIERS[name][scale] for name in order]
        assert values == sorted(values)


def test_minimum_age_is_the_highest_over_the_stars_scales():
    class Stub:
        pass

    stub = Stub()
    values, _tier = CASES["temperate"]
    for attribute, column in (("body_type", "body_type"), ("mass", "mass_kg"), ("gravity", "gravity_g"),
                              ("surface_temperature", "surface_temperature_k"), ("atmospheric_pressure", "atmospheric_pressure_pa")):
        setattr(stub, attribute, values[column])
    assert age_needs.minimum_age_gy(stub, ["fast"]) == age_needs.TIERS["watered"]["fast"]
    assert age_needs.minimum_age_gy(stub, ["fast", "normal"]) == age_needs.TIERS["watered"]["normal"]
    assert age_needs.minimum_age_gy(stub, []) == 0.0


def test_class_table_no_longer_carries_age_ranges():
    assert not any("age_ranges" in data for data in tuning.PLANET_CLASSES.values())


def test_a_broken_rule_file_is_refused(monkeypatch):
    document = copy.deepcopy(age_needs.DOCUMENT)
    document["rows"][-1]["when"] = [{"quantity": "gravity_g", "min": 0.0}]
    monkeypatch.setattr(age_needs.json, "load", lambda handle: document)
    with pytest.raises(age_needs.AgeNeedsError):
        age_needs._load()
