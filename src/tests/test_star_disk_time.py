# tests/test_star_disk_time.py

"""GEN.1551 slice 2: a star's own disk time is read where the disk's end matters: the time of the architecture's losses, the start of the cluster's unbinding, the end of the pebble drift, and the disk of a close pair's young star."""

import types

import pytest

from planetgen.generation.system import StarSystem
from planetgen.physics import constants, disk_dispersal, disk_truncation, planet_architecture as pa, pre_main_sequence, radiation_field
from planetgen.util import draw

SUN = constants.SOLAR_MASS_TO_KG


def _star(mass_sol=1.0):
    return types.SimpleNamespace(mass=mass_sol * SUN, initial_mass_sol=mass_sol, snow_line_au=2.7, stars=None)


def test_a_stars_disk_time_is_its_own_draw_and_a_close_companion_makes_it_shorter():
    lone = disk_truncation.disk_time_yr(_star(0.8))
    assert lone == pytest.approx(pre_main_sequence.dispersal_yr(0.8))
    star = _star(0.8)
    radiation_field.attach(star, [radiation_field.Companion(types.SimpleNamespace(mass=0.5 * SUN), 1, 10.0, 0.0, None, None)])
    assert disk_truncation.has_close_companion(star)
    assert disk_truncation.disk_time_yr(star) <= lone
    wide = _star(0.8)
    radiation_field.attach(wide, [radiation_field.Companion(types.SimpleNamespace(mass=0.5 * SUN), 1, 400.0, 0.0, None, None)])
    assert disk_truncation.disk_time_yr(wide) == lone, "a wide companion changes nothing"


def test_the_stars_with_the_same_mass_have_one_time_and_different_masses_have_their_own():
    times = {disk_truncation.disk_time_yr(_star(1.0 + 0.001 * k)) for k in range(10)}
    assert len(times) > 1
    assert disk_truncation.disk_time_yr(_star(1.0)) == disk_truncation.disk_time_yr(_star(1.0))


def test_the_losses_of_a_build_are_recorded_at_the_end_of_the_stars_disk():
    stream = draw.Stream("star-disk-time-architecture")
    arch = pa.build(SUN, 2.7, 1.0, 0.0, 60.0, stream=stream, disk_time_yr=1.25e6)
    assert arch.lost.disk_time_yr == 1.25e6
    losses = pa.losses(arch)
    timed = [loss for loss in losses if loss.time_yr is not None]
    assert all(loss.time_yr == 1.25e6 for loss in timed if loss.kind in ("swept", "thrown", "sub_minimum"))


def _young(age_gy, mass_sol, kind="class-ii"):
    return types.SimpleNamespace(mass=mass_sol * SUN, initial_mass_sol=mass_sol, age=age_gy, disk_kind=kind, disk_mass_solar=0.01, disk_inner_au=0.05, disk_outer_au=60.0)


def test_a_young_star_with_a_close_companion_loses_its_class_ii_disk_once_the_shorter_time_has_passed(monkeypatch):
    """The star was built with the time of a lone star; the system takes the disk off it when its close companion's shorter time of the same draw is already past."""
    mass = next(0.8 + 0.001 * k for k in range(500) if 2.0e6 < pre_main_sequence.dispersal_yr(0.8 + 0.001 * k, True) < 5.0e6)       # a star whose close time falls after the embedded classes (0.5 Myr)
    close_yr, lone_yr = pre_main_sequence.dispersal_yr(mass, True), pre_main_sequence.dispersal_yr(mass)
    assert close_yr <= lone_yr
    monkeypatch.setattr(disk_truncation, "has_close_companion", lambda star: True)
    gone, kept, embedded = _young(close_yr * 1.5e-9 + 1e-9, mass), _young(close_yr * 0.5e-9 + 0.5e-3, mass), _young(0.1e-3, mass, kind="class-i")
    system = types.SimpleNamespace(stars=[gone, kept, embedded])
    StarSystem._lose_disks_to_close_companions(system)
    assert gone.disk_kind is None and gone.disk_mass_solar is None and gone.disk_inner_au is None and gone.disk_outer_au is None
    assert kept.disk_kind == "class-ii" and kept.disk_mass_solar == 0.01
    assert embedded.disk_kind == "class-i", "the embedded classes are by age"
    monkeypatch.setattr(disk_truncation, "has_close_companion", lambda star: False)
    lone = _young(min(lone_yr * 0.9, close_yr * 1.5) * 1e-9, mass)
    StarSystem._lose_disks_to_close_companions(types.SimpleNamespace(stars=[lone]))
    assert lone.disk_kind == "class-ii", "a star with no close companion keeps its own time"
