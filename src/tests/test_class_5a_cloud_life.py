# tests/test_class_5a_cloud_life.py

"""GEN.1307: Class 5a (cloud-deck plankton) lives in the sulfuric-acid cloud layer, not on the ground: no animals, a clock at the layer's temperature, and the text says where."""

import copy

import pytest

from planetgen import tuning
from planetgen.describe import build, life_facts
from planetgen.physics import acid_medium, atmosphere, life_classes, life_complexity, life_timeline
from tests.test_life_classes import EARTH


def hot_acid_world(surface_k=1500.0, pressure_pa=2.0e8):
    """A dense, hot, dry air with sulfur dioxide and a trace of water: its dry adiabat cools through a liquid acid cloud (about 370 K) far above the ground."""
    values = copy.deepcopy(EARTH)
    values.update(surface_temperature_k=surface_k, atmospheric_pressure_pa=pressure_pa, hydrosphere="dry", ocean_ph=None, water_activity=None,
                  water_mass_fraction=0.0, gravity_g=1.0, energy_flux_w_m2=2000.0)
    values[atmosphere.GAS_FIELD] = {"CO2": pressure_pa * 0.95, "N2": pressure_pa * 0.02, "SO2": pressure_pa * 0.001, "H2O": pressure_pa * 0.001, "O2": pressure_pa * 0.01}
    return values


def test_the_rule_file_gives_class_5a_a_cloud_habitat_and_no_animals_and_no_other_class_a_cap():
    row = life_classes.BY_ID["5a"]
    assert row.habitat == "acid_cloud" and row.animal_cap["grade"] == 0 and row.max_animal_grade == 0
    assert row.animal_cap["tier"] in life_classes.TIERS and row.animal_cap["tweak"] in life_classes.TWEAKS and row.animal_cap["basis"]
    assert [c.id for c in life_classes.CLASSES if c.animal_cap or c.habitat] == ["5a"]
    assert all(c.max_animal_grade == life_classes.ANIMAL_GRADES for c in life_classes.CLASSES if c.id != "5a")


def test_a_cap_stops_the_grade_whatever_the_index_and_leaves_the_stage():
    plankton, plain = life_classes.BY_ID["5a"], life_classes.BY_ID["1a"]
    big = life_complexity.assess(plankton, 1.0, 21.0, 300.0)         # an index far over every threshold
    assert big.index > tuning.LIFE_COMPLEXITY["thresholds"]["A3"] and big.grade == 0 and big.top_stage == "multicellularity"
    assert life_complexity.assess(plankton, 1.0, 21.0, 300.0, anoxic_animals=True).grade == 0       # not even the anoxic route
    assert life_complexity.assess(plain, 1.0, 21.0, 300.0).grade == 3                              # the same inputs give an uncapped class its animals


def test_the_cloud_class_runs_its_clock_at_the_layer_and_the_other_classes_at_the_surface():
    values = hot_acid_world()
    layer = acid_medium.medium(values, cloud_only=True)
    assert layer.layer == "cloud" and 283.0 <= layer.temperature_k <= 373.0 < values["surface_temperature_k"]
    assert life_classes.clock_temperature_k(life_classes.BY_ID["5a"], values) == layer.temperature_k
    assert life_classes.clock_temperature_k(life_classes.BY_ID["1a"], values) == values["surface_temperature_k"]
    airless = dict(values, **{atmosphere.GAS_FIELD: {"N2": 1.0e3}})
    assert life_classes.clock_temperature_k(life_classes.BY_ID["5a"], airless) == values["surface_temperature_k"]      # no layer: the surface's
    assert life_timeline.rate_factor(layer.temperature_k) <= life_timeline.rate_factor(values["surface_temperature_k"])


def test_the_text_says_where_the_cloud_life_lives():
    values = hot_acid_world()
    said = life_classes.habitat_words(life_classes.BY_ID["5a"], values)
    assert said.startswith("aloft in the sulfuric-acid cloud layer, at ") and " K and " in said and said.endswith(" kPa")
    assert life_classes.habitat_words(life_classes.BY_ID["1a"], values) is None
    facts = life_facts.LifeFacts("planet", "5a", "Acid-world life", "fast", (), (), None, said)
    assert f"It lives {said}." in build._life_class_sentences(facts)
    assert len(build._life_class_sentences(life_facts.LifeFacts("planet", "1a", "Oxygenic", "normal", (), (), None))) == 1


@pytest.mark.parametrize("bad", [{"habitat": "orbit"}, {"animal_cap": {"grade": 4, "tier": "T9", "tweak": "L0", "basis": "x"}},
                                 {"animal_cap": {"grade": 0, "tier": "T9", "tweak": "L0"}}])
def test_a_bad_habitat_or_cap_is_refused(bad):
    import json
    row = json.loads(json.dumps(next(r for r in life_classes.DOCUMENT["classes"] if r["id"] == "5a")))
    row.update(bad)
    with pytest.raises(life_classes.LifeClassError):
        life_classes._check(row)
