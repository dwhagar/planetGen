# tests/test_young_star_disk_kept.py

"""
GEN.1550: the belts of a young star are built after its class 0, I or II disk is drawn (`pre_main_sequence.draw_disk`), and `StarSystem._set_debris_disk` used to clear every disk field
before it asked whether the star could hold a cold debris disk, so a young star that kept a belt lost its own disk. Only a debris disk is that method's to set or clear.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import debris_dust
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw

FIELDS = ("disk_kind", "disk_mass_solar", "disk_inner_au", "disk_outer_au")


def _young(age_gy):
    draw.set_run_seed(5)
    return Star.from_params(evolution.star_params(1.0, age_gy, evolution.evolve_star(1.0, age_gy)), SystemConfig(), name="Young")


def _set(star):
    StarSystem.__new__(StarSystem)._set_debris_disk([], star)


@pytest.mark.parametrize("age_gy, kind", [(0.00005, "class-0"), (0.0003, "class-i"), (0.002, "class-ii")])
def test_a_young_stars_own_disk_survives_the_debris_step(age_gy, kind):
    star = _young(age_gy)
    assert star.disk_kind == kind
    before = tuple(getattr(star, name) for name in FIELDS)
    _set(star)
    assert tuple(getattr(star, name) for name in FIELDS) == before


def test_a_debris_disk_left_by_an_earlier_pass_is_still_cleared_when_its_belt_is_gone():
    star = _young(4.6)
    star.disk_kind, star.disk_mass_solar, star.disk_inner_au, star.disk_outer_au = debris_dust.DEBRIS_KIND, None, 30.0, 60.0
    _set(star)
    assert all(getattr(star, name) is None for name in FIELDS)
