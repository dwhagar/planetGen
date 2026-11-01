# tests/test_belt_dose.py

"""
GEN.983: the dose a moon takes from the radiation belts of its parent's field (`radiation_belts.moon_exposure_msv_yr`, `radiation.apply_parent_belts`),
and that it reaches the stored surface dose, the habitability score's dose and the shielding plans.
"""

import math
import types

import pytest

from planetgen.physics import constants, magnetism, radiation, radiation_belts as rb
from tests.test_radiation_belts import JUPITER_MOMENT

AU_KM = constants.AU_M / 1000.0


def _parent(moment=JUPITER_MOMENT):
    """A Jupiter-like parent with just what `belts_of_body` reads."""
    return types.SimpleNamespace(
        magnetic_moment_a_m2=moment, radius=constants.JUPITER_RADIUS_KM, mass=constants.JUPITER_MASS_TO_KG, rotation_period_hours=9.925,
        star=types.SimpleNamespace(), _star_distance_au=5.2, surface_temperature=None, atmospheric_pressure=None)


@pytest.fixture(autouse=True)
def _wind(monkeypatch):
    monkeypatch.setattr(magnetism, "wind_pressure_pa", lambda star, distance: magnetism.EARTH_WIND_PA / distance ** 2)
    monkeypatch.setattr(rb, "air_scale_height_m", lambda body: None)


def _orbit_au(parent, radii):
    return radii * parent.radius / AU_KM


def test_a_moon_inside_the_belts_takes_their_dose_in_sieverts_from_half_the_sky():
    parent = _parent()
    belts = rb.belts_of_body(parent)
    radii = 9.4
    expected = sum(rb.BELT_WEIGHT[belt.species] * belt.dose_gy_yr(radii) for belt in belts) * 1000.0 * rb.SURFACE_VIEW_FRACTION
    assert expected > 0
    assert rb.moon_exposure_msv_yr(parent, _orbit_au(parent, radii)) == pytest.approx(expected)


def test_a_jupiter_strength_belt_gives_a_moon_thousands_of_sieverts_a_year():
    parent = _parent()
    assert rb.moon_exposure_msv_yr(parent, _orbit_au(parent, 9.4)) > 1e6


def test_a_moon_outside_the_belts_or_a_parent_without_a_field_takes_nothing():
    parent = _parent()
    assert rb.moon_exposure_msv_yr(parent, _orbit_au(parent, 500.0)) == 0.0
    assert rb.moon_exposure_msv_yr(_parent(moment=0.0), _orbit_au(parent, 5.0)) == 0.0
    assert rb.moon_exposure_msv_yr(parent, 0.0) == 0.0


def test_a_stronger_field_gives_a_moon_more():
    parent = _parent()
    weak = _parent(moment=JUPITER_MOMENT / 100.0)
    assert rb.moon_exposure_msv_yr(weak, _orbit_au(parent, 5.0)) < rb.moon_exposure_msv_yr(parent, _orbit_au(parent, 5.0))


def test_the_air_cuts_the_belts_by_e_for_each_stopping_column():
    assert rb.air_cut(0.0) == 1.0
    assert rb.air_cut(rb.STOP_COLUMN_G_CM2) == pytest.approx(1.0 / math.e)
    assert rb.air_cut(100.0) < 1e-20


def _moon(pressure_pa=0.0):
    return types.SimpleNamespace(
        distance=0.0, atmospheric_pressure=pressure_pa, gravity=0.1, surface_dose_msv_yr=10.0, dose_belt_msv_yr=0.0, _belt_exposure_msv_yr=None)


def test_apply_parent_belts_adds_the_dose_to_the_stored_surface_dose_once():
    parent = _parent()
    moon = _moon()
    moon.distance = _orbit_au(parent, 9.4)
    exposure = rb.moon_exposure_msv_yr(parent, moon.distance)
    radiation.apply_parent_belts(moon, parent)
    assert moon.dose_belt_msv_yr == pytest.approx(exposure)
    assert moon.surface_dose_msv_yr == pytest.approx(10.0 + exposure)
    radiation.apply_parent_belts(moon, parent)       # again: the same, not doubled
    assert moon.surface_dose_msv_yr == pytest.approx(10.0 + exposure)


def test_a_moons_air_cuts_its_belt_dose_by_the_column():
    parent = _parent()
    bare, aired = _moon(), _moon(pressure_pa=1e5)
    for moon in (bare, aired):
        moon.distance = _orbit_au(parent, 9.4)
        radiation.apply_parent_belts(moon, parent)
    column = radiation.column_g_cm2(1e5, 0.1)
    assert aired.dose_belt_msv_yr == pytest.approx(bare.dose_belt_msv_yr * rb.air_cut(column))
    assert aired.dose_belt_msv_yr < 0.01 * bare.dose_belt_msv_yr


def test_the_total_dose_has_the_belt_part_and_old_rows_without_it_still_total():
    assert radiation.total_dose_msv_yr(2.0, 1.5, 1.0, 0.5) == pytest.approx(4.5)
    assert radiation.total_dose_msv_yr(2.0, 1.5, 1.0, 0.5, 100.0) == pytest.approx(104.5)
    assert radiation.total_dose_msv_yr(2.0, 1.5, 1.0, 0.5, None) == pytest.approx(4.5)


def test_a_generated_system_stores_the_belt_dose_in_its_surface_dose():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw
    from tests.fuzz_support import deterministic_entropy
    moons_in_belts = 0
    for seed in range(5, 60):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            if planet.body_type == "a":
                continue
            if planet.body_type == "g":
                assert planet.dose_belt_msv_yr is None
            else:
                assert planet.dose_belt_msv_yr == 0.0
            for moon in planet.moons:
                assert moon.dose_belt_msv_yr >= 0.0
                assert moon.surface_dose_msv_yr >= moon.dose_belt_msv_yr + moon.dose_ground_msv_yr
                assert moon.dose_belt_msv_yr <= radiation_exposure(planet, moon) * (1.0 + 1e-9)
                moons_in_belts += moon.dose_belt_msv_yr > 0.0
        if moons_in_belts:
            break
    assert moons_in_belts > 0


def radiation_exposure(planet, moon):
    return rb.moon_exposure_msv_yr(planet, moon.distance)


def test_the_shielding_plan_cuts_the_belt_part_with_the_stopping_column():
    from planetgen.physics.equipment import _Shield
    from planetgen.physics.habitability import World

    def shield(belt):
        return _Shield(World(pressure_kpa=0.0, temperature_c=-150.0, gases={}, surface_dose_sv_yr=0.011 + belt, gravity_ms2=1.6,
                             dose_gcr_sv_yr=0.01, dose_sep_sv_yr=0.001, dose_ground_sv_yr=0.0, dose_belt_sv_yr=belt))

    bare = shield(100.0).dose(((1.0, 0.0),))
    assert bare == pytest.approx(100.011)
    behind = shield(100.0).dose(((1.0, 3 * rb.STOP_COLUMN_G_CM2),))
    assert behind < 0.06 * bare             # three stopping columns leave e^-3 of the belt
    assert shield(100.0).dose(((1.0, 3 * rb.STOP_COLUMN_G_CM2),)) > shield(0.0).dose(((1.0, 3 * rb.STOP_COLUMN_G_CM2),))
