# tests/test_gen_debris_disk.py

"""
GEN.1071: the cold debris disk of a main-sequence star is the dust of its outer-disk belt (`physics/debris_dust.py`), found when its light is bright enough for a survey,
and its edges are the belt's. The dust law is that of [R4785]; the shares of stars it gives are compared with the Herschel DEBRIS survey [R4081].
"""

import math

import pytest

from planetgen.describe import sheet
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import belt_mass, constants, debris_dust
from planetgen.util import draw

from tests.test_belt_birth import _belts, _config

SURVEY_SHARE_BY_MASS = ((0.45, 0.022), (0.8, 0.130), (1.04, 0.143), (1.4, 0.239), (math.inf, 0.24))
"""tuple: The share of stars with a cold disk by the mass band that ends at the first number ([R4081]: M 2.2, K 13.0, G 14.3, F 23.9 percent, A 24)."""


def test_a_belt_at_its_steady_state_has_dust_that_falls_as_one_over_the_age():
    lower, upper = 0.75, 1.25
    f = debris_dust.fractional_luminosity(lower, upper, 1.0, 1000.0, 2000.0, 1.0e-12 + belt_mass.collisional_constant(lower, upper, 1.0, 2000.0) / 1000.0, 1.0)
    assert f > 0.0
    # the same belt at half the age has twice the steady dust
    f2 = debris_dust.fractional_luminosity(lower, upper, 1.0, 500.0, 2000.0, belt_mass.collisional_constant(lower, upper, 1.0, 2000.0) / 500.0, 1.0)
    assert f2 == pytest.approx(2.0 * f, rel=1e-6)


def test_the_dust_goes_with_the_mass_the_collisions_have_left():
    lower, upper, age = 20.0, 60.0, 4000.0
    ceiling = belt_mass.collisional_constant(lower, upper, 1.0, 100.0) / age
    full = debris_dust.fractional_luminosity(lower, upper, 1.0, age, 100.0, ceiling, 1.0)
    half = debris_dust.fractional_luminosity(lower, upper, 1.0, age, 100.0, 0.5 * ceiling, 1.0)
    assert half == pytest.approx(0.5 * full, rel=1e-9)
    assert debris_dust.fractional_luminosity(lower, upper, 1.0, age, 100.0, 0.0, 1.0) == 0.0


def test_a_brighter_star_lights_the_same_belt_less():
    args = (20.0, 60.0)
    mass = 0.01
    near = debris_dust.fractional_luminosity(*args, 1.0, 4000.0, 100.0, mass, 1.0)
    far = debris_dust.fractional_luminosity(*args, 4.0, 4000.0, 100.0, mass, 1.0)
    assert far == pytest.approx(near / 2.0, rel=1e-9)


def test_a_largest_body_beyond_the_cascade_is_reckoned_at_the_cascade_ceiling():
    cap = debris_dust.CASCADE_MAX_KM
    args = (20.0, 60.0, 1.0, 4000.0)
    huge = debris_dust.fractional_luminosity(*args, 12000.0, 1.0e6, 1.0)
    assert huge > 0.0 and huge == pytest.approx(debris_dust.fractional_luminosity(*args, cap, 1.0e6, 1.0), rel=1e-12)


def test_a_belt_whose_cascade_cannot_start_makes_no_dust(monkeypatch):
    monkeypatch.setattr(belt_mass, "collisional_constant", lambda *a, **k: -1.0)
    assert debris_dust.fractional_luminosity(20.0, 60.0, 1.0, 4000.0, 100.0, 1.0, 1.0) == 0.0


def test_a_bad_belt_is_refused():
    for bad in ((60.0, 20.0, 1.0, 1000.0, 100.0, 1.0, 1.0), (20.0, 60.0, 0.0, 1000.0, 100.0, 1.0, 1.0), (20.0, 60.0, 1.0, 0.0, 100.0, 1.0, 1.0)):
        with pytest.raises(ValueError):
            debris_dust.fractional_luminosity(*bad)


def test_the_cut_separates_a_found_disk_from_a_faint_one():
    cut, sun = debris_dust.DETECTABLE_FRACTIONAL_LUMINOSITY, debris_dust.REFERENCE_TEMPERATURE_K
    assert debris_dust.is_detectable(cut, sun) and debris_dust.is_detectable(10.0 * cut, sun)
    assert not debris_dust.is_detectable(0.5 * cut, sun) and not debris_dust.is_detectable(0.0, sun)


def test_a_hotter_star_shows_a_fainter_disk_as_the_cube_of_its_temperature():
    sun = debris_dust.REFERENCE_TEMPERATURE_K
    assert debris_dust.detectable_fractional_luminosity(sun) == pytest.approx(debris_dust.DETECTABLE_FRACTIONAL_LUMINOSITY)
    assert debris_dust.detectable_fractional_luminosity(2.0 * sun) == pytest.approx(debris_dust.DETECTABLE_FRACTIONAL_LUMINOSITY / 8.0)
    assert debris_dust.is_detectable(1.5e-5, 9500.0) and not debris_dust.is_detectable(1.5e-5, 5772.0)
    with pytest.raises(ValueError):
        debris_dust.detectable_fractional_luminosity(0.0)


def test_a_belt_is_the_sum_of_its_rings_and_a_thin_one_is_the_ring():
    flat = lambda r: 1.0
    thin = (50.0, 55.0, 2.0, 300.0, 2000.0, 0.01, 1.0)
    assert debris_dust.belt_fractional_luminosity(*thin, flat) == pytest.approx(debris_dust.fractional_luminosity(*thin), rel=0.02)
    wide = (10.0, 100.0, 2.0, 300.0, 2000.0, 5.0, 1.0)
    assert debris_dust.belt_fractional_luminosity(*wide, flat) > 0.0
    assert debris_dust.belt_fractional_luminosity(10.0, 100.0, 2.0, 300.0, 2000.0, 0.0, 1.0, flat) == 0.0


def test_a_belt_whose_mass_lies_near_its_inner_edge_shines_more_than_one_whose_mass_lies_far_out():
    wide = (10.0, 200.0, 2.0, 300.0, 2000.0, 20.0, 1.0)
    inner = debris_dust.belt_fractional_luminosity(*wide, lambda r: r ** -1.5)
    outer = debris_dust.belt_fractional_luminosity(*wide, lambda r: r ** 0.5)
    assert inner > outer > 0.0                                        # the dust goes as the mass over the radius squared
    assert inner > debris_dust.fractional_luminosity(*wide)           # GEN.1278: the mass at the middle radius of a wide belt gives far too little


def test_a_belt_with_no_solids_has_no_dust_and_a_bad_belt_is_refused():
    assert debris_dust.belt_fractional_luminosity(10.0, 100.0, 2.0, 300.0, 2000.0, 5.0, 1.0, lambda r: 0.0) == 0.0
    with pytest.raises(ValueError):
        debris_dust.belt_fractional_luminosity(100.0, 10.0, 2.0, 300.0, 2000.0, 5.0, 1.0, lambda r: 1.0)


def _systems(count, **overrides):
    config = _config(**overrides)
    for seed in range(1, count + 1):
        draw.set_run_seed(seed)
        yield StarSystem(system_config=config)


def test_a_star_with_a_disk_has_the_edges_of_its_outer_belt():
    seen = 0
    for system in _systems(30):
        star = system.star
        outer = [belt for belt in _belts(system) if belt.origin == "outer-disk"]
        if star.disk_kind:
            seen += 1
            assert star.disk_kind == "debris" and star.disk_mass_solar is None and star.disk_inner_au < star.disk_outer_au
            assert any(belt.lower_limit == star.disk_inner_au and belt.upper_limit == star.disk_outer_au for belt in outer)
        else:
            assert star.disk_inner_au is None and star.disk_outer_au is None
    assert seen > 0


def test_no_giant_has_a_debris_disk():
    star = Star.from_params({"type": "K2III", "yerkes_class": "III", "mass_kg": 1.2 * constants.SOLAR_MASS_TO_KG, "radius_km": 2.0e7, "temperature_k": 4500.0,
                             "luminosity_w": 100.0 * constants.SOLAR_LUMINOSITY, "age_gy": 8.0, "lifespan_gy": 9.0, "initial_mass_sol": 1.2, "phase_end_age_gy": 9.0,
                             "metallicity_fe_h": 0.0}, SystemConfig(), name="Test")
    assert star.disk_kind is None


def test_the_sheet_rule_says_which_stars_have_disk_rows():
    star = {"disk_kind": "debris"}
    assert sheet.applies("star", "G2V", "disk_kind", star) and sheet.applies("star", "G2V", "disk_outer_au", star)
    assert not sheet.applies("star", "G2V", "disk_mass_solar", star)
    assert sheet.applies("star", "M5V", "disk_mass_solar", {"disk_kind": "class-ii"})
    assert not sheet.applies("star", "G2V", "disk_kind", {"disk_kind": None})
    assert ("star", "disk_kind") not in sheet.OPEN_GAPS and ("star", "disk_mass_solar") not in sheet.OPEN_GAPS


def test_the_share_of_sun_like_stars_with_a_disk_is_near_the_surveys():
    held = sum(bool(system.star.disk_kind) for system in _systems(60))
    survey = dict((upper, share) for upper, share in SURVEY_SHARE_BY_MASS)[1.04]
    assert 0.03 <= held / 60 <= 2.0 * survey, held         # the tuned cut (stand-in DU1) gives about 15 percent of G2V stars against the survey's 14 (GEN.1261: 1.3 percent at the earlier cut 1.5e-4, once the cold belts lost their unformed share)


def test_a_hot_young_star_holds_a_disk_as_often_as_the_surveys_find_one():
    # GEN.1278: A stars have no debris disk at the Sun's limit (0 of 60 at 3e-5); the limit that falls with the star's temperature and the dust of the rings of a wide belt
    # give them one in about six (stand-in DU1/DU4/DU5; the surveys find a quarter, 20 to 33 percent)
    held = sum(bool(system.star.disk_kind) for system in _systems(40, STAR_TYPE="A1V"))
    assert 1 <= held <= 20, held
