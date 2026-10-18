"""
A supernova remnant in a superbubble (GEN.298): the remnant sized at build reads the cavity's interior density in its
Sedov-Taylor law where its sector lies in one, over a seeded bubble; outside, nothing changes.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.galaxy import profile, superbubble_field
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.generation import run_sector
from planetgen.generation.bright_stars import FillContext
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant, sedov_taylor_radius_ly
from planetgen.util import draw

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
SEED = bytes(range(16))
EDGE_PC = 10.0
E_VALUE = 1.0


def _bubble():
    """A seeded bubble of the galaxy and the profile at its centre."""
    bubbles = superbubble_field.all_bubbles(SEED, SHAPE, 6000.0)
    bubble = max(bubbles, key=lambda candidate: candidate.radius_pc)
    here = profile.profile_at(bubble.center_pc, SHAPE, EDGE_PC, E_VALUE, galaxy_seed=SEED)
    assert here.gas.superbubble == bubble
    return bubble, here


def _outside():
    """A profile far above the disk, in no bubble."""
    point = (8128.0, 0.0, 60000.0)
    here = profile.profile_at(point, SHAPE, EDGE_PC, E_VALUE, galaxy_seed=SEED)
    assert here.gas.superbubble is None
    return here


def _remnant(ambient, stream="remnant"):
    config = SystemConfig()
    config.AMBIENT_DENSITY_CM3 = ambient
    with draw.bound(stream):
        return SupernovaRemnant(config)


def test_a_remnant_in_a_bubble_is_sized_by_the_cavitys_density():
    bubble, _ = _bubble()
    inside = _remnant(bubble.ambient_density_cm3)
    assert inside.radius_ly == pytest.approx(sedov_taylor_radius_ly(inside.age_years, bubble.ambient_density_cm3))


def test_the_same_draws_give_a_larger_remnant_in_the_thin_cavity_and_the_same_contents():
    bubble, _ = _bubble()
    for stream in ("a", "b", "c", "d"):
        open_sky = _remnant(None, stream)
        inside = _remnant(bubble.ambient_density_cm3, stream)
        assert open_sky.radius_ly == pytest.approx(sedov_taylor_radius_ly(open_sky.age_years, open_sky.density_cm3))
        assert (inside.age_years, inside.density_cm3, inside.remnant_class) == \
            (open_sky.age_years, open_sky.density_cm3, open_sky.remnant_class)
        assert inside.density_cm3 > bubble.ambient_density_cm3
        assert inside.radius_ly >= open_sky.radius_ly


def test_the_fill_hands_the_sectors_bubble_density_to_its_remnants():
    bubble, here = _bubble()
    assert run_sector._bubble_density(FillContext(here)) == bubble.ambient_density_cm3
    assert run_sector._bubble_density(FillContext(_outside())) is None
    assert run_sector._bubble_density(SimpleNamespace()) is None
    assert run_sector._bubble_density(None) is None


def test_a_scattered_remnant_is_sized_by_the_density_it_is_given():
    class Sector:
        def add_phenomenon(self, phenomenon, kind, position=None):
            return phenomenon

    args = SimpleNamespace()
    row = {"kind": "supernova-remnant", "subtype": None, "mass_solar": None}
    bubble, _ = _bubble()
    with draw.bound("scattered"):
        inside = run_sector._build_scattered(Sector(), args, row, (0.0, 0.0, 0.0), 1000.0, None, None,
                                             bubble.ambient_density_cm3)
    assert inside.radius_ly == pytest.approx(sedov_taylor_radius_ly(inside.age_years, bubble.ambient_density_cm3))
    with draw.bound("scattered"):
        open_sky = run_sector._build_scattered(Sector(), args, row, (0.0, 0.0, 0.0), 1000.0, None, None)
    assert open_sky.radius_ly == pytest.approx(sedov_taylor_radius_ly(open_sky.age_years, open_sky.density_cm3))


def test_the_remnant_factor_still_ignores_the_bubble():
    bubble, in_bubble = _bubble()
    without = profile.profile_at(bubble.center_pc, SHAPE, EDGE_PC, E_VALUE)
    assert without.gas.superbubble is None
    assert in_bubble.kind_factors["supernova-remnant"] == without.kind_factors["supernova-remnant"] > 0.0
