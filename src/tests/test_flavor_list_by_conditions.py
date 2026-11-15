"""GEN.1399, GEN.1430: the flavor line is offered by the body's conditions, never by its class letter (the class is decided last)."""
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import life
from planetgen.generation.config import SystemConfig
from planetgen.physics import flavor_tree, planet_classes, surface_liquid_world


def _body(body_type, surface_liquid=None, life_stage=None, life_class=None, life_class_id=None, life_animal_grade=0):
    return SimpleNamespace(body_type=body_type, surface_liquid=surface_liquid, life_class=life_class, life_class_id=life_class_id, life_animal_grade=life_animal_grade, life_stage=life_stage, system_config=SystemConfig(),
                           flavor_key=None, flavor_text_count=0, mass=5.0e24, radius=6400.0, distance=1.0)


@pytest.fixture(autouse=True)
def _always_roll_and_never_read_the_class(monkeypatch):
    monkeypatch.setattr(tuning, "FLAVOR_CHANCE_PLANET", 1.0)

    def _forbidden(*_args, **_kwargs):
        raise AssertionError("decide_flavor_text read the class letter")
    monkeypatch.setattr(planet_classes, "classify", _forbidden)
    monkeypatch.setattr(planet_classes, "classify_detail", _forbidden)


def _family_of(body):
    life.decide_flavor_text(body)
    return body.flavor_key.split(":")[0] if body.flavor_key else None


def test_a_rock_with_a_ground_gets_the_ground_lines():
    assert _family_of(_body("t")) == "planet"
    assert _family_of(_body("t", "water")) == "planet"


@pytest.mark.parametrize("liquid", surface_liquid_world.LAVA_IDS)
def test_a_rock_that_is_molten_all_over_gets_the_lines_of_a_body_with_no_ground(liquid):
    assert _family_of(_body("t", liquid)) == "orbital"


def test_a_gas_giant_gets_the_lines_of_a_body_with_no_ground():
    assert _family_of(_body("g")) == "orbital"


def test_a_world_with_multicellular_life_is_offered_the_life_lines_as_well_as_the_ground_lines():
    options = flavor_tree.walk("planet", _body("t", life_stage="multicellularity", life_class=1, life_class_id="1a", life_animal_grade=2))
    families = {option.id.split(":")[0] for option in options}
    assert families == {"planet", "habitable", "life"}
    assert {option.id.split(":")[0] for option in flavor_tree.walk("planet", _body("t", life_stage=None, life_class=1, life_class_id="1a"))} == {"planet"}
