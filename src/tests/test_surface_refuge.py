"""
Surface refuges for oxygen-making life (GEN.1723, `planetgen.physics.surface_refuge`, reports/life-refuges-research-gen1722-2026-11-14.md).
"""

import copy
from types import SimpleNamespace

import pytest

from planetgen.physics import biogenic_oxygen as bo, habitability_world, life_classes, surface_refuge

from tests.test_life_classes import world

LIGHT_CLASSES = ("1a", "1b", "1c")


def _dry(depth_km, temperature_k=300.0, thickness_km=1.0):
    return world(hydrosphere="dry", aquifer_depth_km=depth_km, aquifer_thickness_km=thickness_km, surface_temperature_k=temperature_k)


def test_production_per_area_is_zero_under_the_limit_and_grows_with_the_hours():
    assert surface_refuge.production_share(surface_refuge.MIN_HOURS_YR - 1.0) == 0.0
    shares = [surface_refuge.production_share(hours) for hours in (75.0, 150.0, 500.0, 1000.0, 2000.0, 4380.0, 8760.0)]
    assert shares == sorted(shares) and shares[0] > 0.0
    assert shares[-1] == pytest.approx(surface_refuge.COVER_MAX)         # a year of light and water: the full cover, 0.1 of an ocean's production
    assert shares[0] < 1.0e-9           # at the limit it is a duty cycle of 75 h over a year, times a cover of 0.001, with the exponent at its largest


def test_an_open_sea_has_the_whole_year_and_leaves_its_source_to_the_ocean():
    refuge = surface_refuge.refuge_of(world())
    assert refuge.hours == surface_refuge.YEAR_HOURS and refuge.area == 0.0
    assert surface_refuge.refuge_of(world(hydrosphere="hycean")).hours == surface_refuge.YEAR_HOURS


def test_a_dry_world_with_its_water_table_at_the_surface_has_seeps_and_a_deeper_table_has_fewer_and_none():
    here = surface_refuge.refuge_of(_dry(0.0))
    assert here.hours == pytest.approx(0.5 * surface_refuge.YEAR_HOURS)          # lit half the time, in the range all the year
    assert here.area == pytest.approx(surface_refuge.SEEP_LAND_SHARE * surface_refuge.production_share(here.hours))
    assert here.temperature_k == pytest.approx(300.0, abs=0.5)
    assert 0.0 < surface_refuge.refuge_of(_dry(0.002)).area < here.area
    for deep in (0.05, 5.0):           # a table 50 m or 5 km down reaches no light
        assert surface_refuge.refuge_of(_dry(deep)) == surface_refuge.NONE


def test_a_dry_world_with_no_aquifer_has_no_refuge_and_scores_no_light_class():
    values = world(hydrosphere="dry")
    assert surface_refuge.refuge_of(values) == surface_refuge.NONE
    scores = life_classes.score_classes(values)
    assert all(scores[c].score == 0 and "lit_liquid_hours" in scores[c].failed for c in LIGHT_CLASSES)


def test_a_seep_world_is_scored_by_the_light_classes_on_the_light_it_gets_and_the_share_of_an_ocean_it_makes():
    scores = life_classes.score_classes(_dry(0.0))
    sea = life_classes.score_classes(world())
    for name in LIGHT_CLASSES:
        assert scores[name].score > 0 and scores[name].failed == []
        assert 0.0 < scores[name].npp < 0.05 * sea[name].npp          # a refuge makes a few hundredths of an ocean's production, not all of it
    assert sea["1a"].score == life_classes.score_classes(world(aquifer_depth_km=None))["1a"].score          # an open sea is scored as before


def test_a_world_under_ice_that_never_melts_has_no_light_class_and_one_that_melts_has_them():
    cold = life_classes.score_classes(world(hydrosphere="ice-covered ocean", surface_temperature_k=230.0))
    assert all(cold[c].score == 0 for c in (*LIGHT_CLASSES, "2"))
    assert cold["3a"].score > 0 and cold["4b"].score > 0          # the chemical classes still live under the ice
    melting = life_classes.score_classes(world(hydrosphere="ice-covered ocean", surface_temperature_k=288.0))
    assert all(melting[c].score > 0 for c in LIGHT_CLASSES)


def test_the_sea_under_a_thin_lid_counts_in_proportion_to_the_light_that_gets_through():
    def lid(km):
        return surface_refuge.refuge_of(world(hydrosphere="ice-covered ocean", surface_temperature_k=230.0, ice_shell_km=km, ocean_fraction=0.7))
    bare, thin, thick = lid(0.0), lid(0.003), lid(0.005)
    assert bare.area == pytest.approx(0.7 * (1.0 - surface_refuge.ICE_ALBEDO)) and bare.temperature_k == pytest.approx(surface_refuge.UNDER_ICE_K)
    assert thin.area == pytest.approx(0.7 * surface_refuge.under_ice_light_share({"hydrosphere": "ice-covered ocean", "ice_shell_km": 0.003}))        # 3 m: a few percent of the light, as the lake mats get
    assert 0.0 < thick.area < thin.area < bare.area
    assert lid(0.01) == surface_refuge.NONE          # 10 m lets through less than the light floor
    assert lid(1.0) == surface_refuge.NONE and lid(None) == surface_refuge.NONE        # a kilometre of ice, or a lid of unknown thickness, lets none through
    scores = life_classes.score_classes(world(hydrosphere="ice-covered ocean", surface_temperature_k=230.0, ice_shell_km=0.003, ocean_fraction=0.7))
    assert all(scores[c].score > 0 for c in LIGHT_CLASSES)


def test_rooted_plants_reach_a_deeper_water_table_and_more_of_the_land():
    microbes = surface_refuge.refuge_of(_dry(0.005))
    rooted = surface_refuge.refuge_of({**_dry(0.005), "native_life": "macroscopic"})
    assert rooted.area > 10.0 * microbes.area > 0.0           # 5 m is deep for a microbe's seep and shallow for a root (desert plants reach 9.5 m on the mean)
    assert surface_refuge.refuge_of({**_dry(0.2), "native_life": "biosphere"}) == surface_refuge.NONE          # 200 m is out of reach of the mean root too
    assert surface_refuge.refuge_of({**_dry(0.005), "native_life": "microbial"}).area == pytest.approx(microbes.area)


def test_a_brine_seep_is_a_refuge_of_the_size_of_its_liquid():
    values = world(hydrosphere="brine seep", ocean_class="chloride brine", ocean_fraction=0.05)
    refuge = surface_refuge.refuge_of(values)
    assert refuge.hours == pytest.approx(0.5 * surface_refuge.YEAR_HOURS)
    assert refuge.area == pytest.approx(0.05 * surface_refuge.production_share(refuge.hours))
    assert life_classes.score_classes(values)["1a"].score > 0


def test_too_hot_or_too_cold_a_seep_is_dry():
    assert surface_refuge.refuge_of(_dry(0.0, temperature_k=360.0)) == surface_refuge.NONE
    assert surface_refuge.refuge_of(_dry(0.0, temperature_k=240.0)) == surface_refuge.NONE


def test_a_wet_place_at_the_cold_edge_makes_next_to_nothing():
    warm, cool, edge = (surface_refuge.refuge_of(_dry(0.0, temperature_k=k)) for k in (290.0, 258.0, 254.0))
    assert warm.area > cool.area > edge.area > 0.0          # the rate falls linearly from 263 K to none at 253 K (LA1b), the gate stays open to 252 K
    assert cool.area == pytest.approx(warm.area * 0.5, rel=0.02) and edge.area < 0.2 * warm.area
    assert edge.hours == warm.hours


def _body(values, **over):
    """A planet of `values` with the attributes `biogenic_oxygen` reads (Earth's size, the Sun at 1 AU)."""
    star = SimpleNamespace(luminosity=3.828e26, temperature=5772.0, age=4.6)
    body = SimpleNamespace(star=star, _star_distance_au=1.0, radius=6371.0, gravity=1.0, surface_temperature=values["surface_temperature_k"], hydrosphere=values["hydrosphere"],
                           internal_heat_flux_w_m2=0.09, **{key: values.get(key) for key in ("aquifer_depth_km", "aquifer_thickness_km", "ocean_fraction", "ocean_class")})
    for key, value in over.items():
        setattr(body, key, value)
    return body


def test_the_refuge_source_is_a_hundredth_of_earths_at_the_most_and_under_the_reductants_of_an_earth_interior(monkeypatch):
    body = _body(_dry(0.0))
    monkeypatch.setattr(habitability_world, "values_from_body", lambda b: _dry(0.0))
    source = bo.productivity(body, bo.area_m2(body))
    assert 0.0 < source <= 0.02 * bo.EARTH_BURIAL_MOL_YR
    assert source < bo.reductant_flux(body, bo.area_m2(body))          # the dry world with Earth's volcanism stays at zero oxygen
    assert bo.profile_of(body).steady_pa == 0.0
    # a world whose interior has gone quiet (its reductant flux under the refuge source) holds some: the lever is R, not the refuge number
    quiet = _body(_dry(0.0), internal_heat_flux_w_m2=0.0)
    assert bo.reductant_flux(quiet, bo.area_m2(quiet)) < bo.productivity(quiet, bo.area_m2(quiet))


def test_an_open_sea_source_is_what_it_was():
    body = _body(world(), ocean_fraction=0.71)
    assert bo.productivity(body, bo.area_m2(body)) == pytest.approx(bo.EARTH_BURIAL_MOL_YR, rel=0.2)
    assert bo.refuge_productivity(_body(world(), ocean_fraction=0.71), bo.area_m2(body)) >= 0.0
