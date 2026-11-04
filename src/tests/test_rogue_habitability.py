# tests/test_rogue_habitability.py

"""A free-floating brown dwarf has no surface and scores nothing (`physics/rogue_habitability.py`); the rogue planets that systems throw out are scored by the
planet chain like any planet (GEN.1111, `test_rogue_planet_parity.py`)."""

import types

from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.physics import rogue_habitability
from planetgen.util import draw


def test_a_brown_dwarf_scores_nothing():
    dwarf = types.SimpleNamespace(planet_type="g")
    rogue_habitability.generate(dwarf)
    assert all(getattr(dwarf, name) is None for name in rogue_habitability.ROGUE_FIELDS)


def test_every_generated_brown_dwarf_has_every_score_field_empty():
    for seed in range(6):
        with draw.bound(seed):
            dwarf = RoguePlanet(SystemConfig())
        assert all(getattr(dwarf, name) is None for name in rogue_habitability.ROGUE_FIELDS)


def test_a_rebuilt_brown_dwarf_keeps_its_empty_score():
    with draw.bound(3):
        dwarf = RoguePlanet(SystemConfig())
    rebuilt = RoguePlanet.from_dict(dwarf.to_dict(), SystemConfig())
    assert all(getattr(rebuilt, name) is None for name in rogue_habitability.ROGUE_FIELDS)
