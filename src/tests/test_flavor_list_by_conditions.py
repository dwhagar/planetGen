"""GEN.1399: the flavor-text draw list is chosen by the body's conditions, never by its class letter (the class is decided last)."""
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import life
from planetgen.generation.config import SystemConfig
from planetgen.physics import planet_classes, surface_liquid_world


def _body(body_type, surface_liquid=None):
    return SimpleNamespace(body_type=body_type, surface_liquid=surface_liquid, life_class=None, life_stage=None, system_config=SystemConfig(),
                           flavor_key=None, flavor_text_count=0)


@pytest.fixture(autouse=True)
def _always_roll_and_never_read_the_class(monkeypatch):
    monkeypatch.setattr(tuning, "FLAVOR_CHANCE_PLANET", 1.0)

    def _forbidden(*_args, **_kwargs):
        raise AssertionError("decide_flavor_text read the class letter")
    monkeypatch.setattr(planet_classes, "classify", _forbidden)
    monkeypatch.setattr(planet_classes, "classify_detail", _forbidden)


def _list_of(body):
    life.decide_flavor_text(body)
    return body.flavor_key.split(":")[0] if body.flavor_key else None


def test_a_rock_with_a_ground_gets_the_planet_list():
    assert _list_of(_body("t")) == "planet"
    assert _list_of(_body("t", "water")) == "planet"


@pytest.mark.parametrize("liquid", surface_liquid_world.LAVA_IDS)
def test_a_rock_that_is_molten_all_over_gets_the_orbital_list(liquid):
    assert _list_of(_body("t", liquid)) == "orbital"


def test_a_gas_giant_gets_the_orbital_list():
    assert _list_of(_body("g")) == "orbital"
