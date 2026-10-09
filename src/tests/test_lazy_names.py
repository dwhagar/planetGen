"""
PERF.43: a phenomenon's word-salad name is drawn on first read, from its
own seed, so a placed phenomenon named by its object ID never pays for it.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.util import draw

CLASSES = [RoguePlanet, InterstellarComet, Nebula, AsteroidField, SupernovaRemnant, BlackHole, NeutronStar]


@pytest.mark.parametrize("cls", CLASSES, ids=lambda cls: cls.__name__)
def test_the_name_waits_until_read_and_does_not_move_other_draws(cls):
    draw.set_run_seed(43)
    unread = cls(SystemConfig())
    draw.set_run_seed(43)
    read = cls(SystemConfig())
    name = read.name
    assert name and name == read.name  # drawn once, then kept
    assert unread.to_dict() == read.to_dict()  # reading it changes nothing else
    draw.set_run_seed(43)
    assert cls(SystemConfig()).name == name  # the same seed names it the same


def test_a_rogue_planets_name_is_not_drawn_at_construction():
    draw.set_run_seed(1)
    planet = RoguePlanet(SystemConfig())
    assert planet.__dict__["_name"] is None
    planet.name = "Given"
    assert planet.name == "Given"
    assert RoguePlanet(SystemConfig(), name="Kept").name == "Kept"
