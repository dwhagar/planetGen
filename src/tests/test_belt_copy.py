"""
GEN.1610: a belt copied or sent to another process before its members were first read carries them (the star they are drawn around is held apart from the belt, in a weak map that is not copied).
"""

import copy
import pickle

import pytest

from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.test_systems import make_config


@pytest.fixture()
def system():
    return StarSystem(system_config=make_config("G2V", PLANETS=False))


def _unread_belt(system):
    with draw.bound(5):
        belt = AsteroidBelt(system.system_config, system.star, 2.5, 2.2, 3.4)
    assert belt._members is None
    return belt


@pytest.mark.parametrize("copier", [copy.copy, copy.deepcopy, lambda belt: pickle.loads(pickle.dumps(belt))], ids=["copy", "deepcopy", "pickle"])
def test_a_copy_of_a_belt_whose_members_were_never_read_has_them(system, copier):
    belt = _unread_belt(system)
    twin = copier(belt)
    assert twin.members and twin.members == belt.members
    assert twin.largest_body_km == belt.largest_body_km and twin.mass_kg == belt.mass_kg
