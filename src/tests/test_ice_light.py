# tests/test_ice_light.py

"""GEN.1726: the light that reaches the water under ice and snow (`physics/ice_light.py`): the law, its inverse, and the measurements it is checked against."""

import math

import pytest

from planetgen.physics import ice_light as il


def test_the_central_values_sit_inside_their_ranges():
    for kind in il.KINDS:
        assert il.K_ICE_RANGE_PER_M[kind][0] <= il.K_ICE_PER_M[kind] <= il.K_ICE_RANGE_PER_M[kind][1]
        assert il.ALBEDO_BARE_RANGE[kind][0] <= il.ALBEDO_BARE[kind] <= il.ALBEDO_BARE_RANGE[kind][1]
    assert il.K_SNOW_RANGE_PER_M[0] <= il.K_SNOW_PER_M <= il.K_SNOW_RANGE_PER_M[1]
    assert il.ALBEDO_SNOW_RANGE[0] <= il.ALBEDO_SNOW <= il.ALBEDO_SNOW_RANGE[1]
    assert il.ALBEDO_MELTING_RANGE[0] <= il.ALBEDO_MELTING <= il.ALBEDO_MELTING_RANGE[1]


def test_snow_takes_far_more_light_per_metre_than_ice():
    assert il.K_SNOW_PER_M > 20.0 * il.K_ICE_PER_M[il.SEA]
    assert il.K_SNOW_PER_M > 80.0 * il.K_ICE_PER_M[il.FRESH]
    assert il.K_ICE_PER_M[il.FRESH] < il.K_ICE_PER_M[il.SEA]


def test_no_lid_lets_through_what_the_surface_does_not_reflect():
    for kind in il.KINDS:
        assert il.transmission(0.0, 0.0, kind) == pytest.approx(1.0 - il.ALBEDO_BARE[kind])
    assert il.transmission(0.0, 0.0, il.SEA, melting=True) == pytest.approx(1.0 - il.ALBEDO_MELTING)


def test_each_metre_of_ice_takes_the_same_share():
    one, two = il.transmission(1.0), il.transmission(2.0)
    assert two / one == pytest.approx(math.exp(-il.K_ICE_PER_M[il.SEA]))
    assert il.transmission(3.0) / il.transmission(2.0) == pytest.approx(two / one)


def test_more_ice_or_more_snow_lets_less_light_through():
    for kind in il.KINDS:
        values = [il.transmission(h, 0.0, kind) for h in (0.0, 0.5, 1.0, 5.0, 20.0)]
        assert values == sorted(values, reverse=True) and values[-1] < values[0]
        snowy = [il.transmission(1.0, s, kind) for s in (0.0, 0.01, 0.05, 0.2, 1.0)]
        assert snowy == sorted(snowy, reverse=True) and snowy[-1] < snowy[0]


def test_a_thin_dusting_of_snow_does_not_make_the_albedo_jump():
    assert il.albedo(1e-9) == pytest.approx(il.ALBEDO_BARE[il.SEA], abs=1e-6)
    assert il.albedo(10.0) == pytest.approx(il.ALBEDO_SNOW, abs=1e-6)
    steps = [il.albedo(s) for s in (0.0, 0.005, 0.02, 0.05, 0.1, 0.3)]
    assert steps == sorted(steps)
    assert max(b - a for a, b in zip(steps, steps[1:])) < 0.2


def test_melting_ice_lets_more_light_through_than_cold_ice():
    assert il.transmission(1.0, melting=True) > il.transmission(1.0)


def test_lake_a_without_snow_is_within_a_third_of_the_measured_share():
    # [R5265]: 1.97 m of ice, snow removed: 5.7 to 6.6 percent of the light under the ice
    t = il.transmission(1.97, 0.0, il.FRESH)
    assert t == pytest.approx(0.085, abs=0.001)
    assert abs(t - 0.066) / 0.066 < 0.3


def test_lake_a_with_its_wind_packed_snow_is_the_measured_order():
    # [R5265]: 2 m of ice and 0.5 m of snow let 0.45 percent through; the snow there measured k 2.1 and albedo 0.97
    t = il.transmission(1.97, 0.5, il.FRESH, k_snow=2.1, albedo_snow=0.97)
    assert 0.003 < t < 0.007


def test_antarctic_pack_ice_in_spring_is_under_a_tenth_of_a_percent():
    # [R5264]: under 0.1 percent through a spring floe with 0.5 m of ice and 0.2 m of snow of albedo 0.8
    assert il.transmission(0.5, 0.2, il.SEA, albedo_snow=0.8) < 1.0e-3


def test_light_below_scales_the_surface_light():
    assert il.light_below(400.0, 1.0, 0.1) == pytest.approx(400.0 * il.transmission(1.0, 0.1))


def test_the_thickest_lid_that_passes_the_light_floor_of_an_earth_bright_sky():
    # GEN.1726 report 3.4: one part in ten thousand is about 6.6 m of bare sea ice and 23 to 30 m of bare clear lake ice (k 0.32 to 0.25)
    assert 6.0 < il.thickest_ice_m(1.0e-4, kind=il.SEA) < 7.0
    assert 22.0 < il.thickest_ice_m(1.0e-4, kind=il.FRESH) < 24.0
    assert 28.0 < il.thickest_ice_m(1.0e-4, kind=il.FRESH, k_ice=0.25) < 31.0


def test_thickest_ice_is_the_inverse_of_the_law_and_snow_costs_ice():
    for kind in il.KINDS:
        for snow in (0.0, 0.05, 0.1):
            h = il.thickest_ice_m(1.0e-6, snow, kind)
            if h > 0.0:
                assert il.transmission(h, snow, kind) == pytest.approx(1.0e-6, rel=1e-9)
        assert il.thickest_ice_m(1.0e-6, 0.1, kind) < il.thickest_ice_m(1.0e-6, 0.0, kind)
    assert il.thickest_ice_m(0.5, kind=il.SEA) == 0.0


def test_the_deepest_snow_that_still_passes_the_floor():
    s = il.deepest_snow_m(1.0e-4, ice_m=0.0)
    assert 0.15 < s < 0.35
    assert il.transmission(0.0, s) == pytest.approx(1.0e-4, rel=1e-6)
    assert il.deepest_snow_m(1.0e-4, ice_m=1.0) < s
    assert il.deepest_snow_m(0.99, ice_m=0.0) == 0.0
    assert il.deepest_snow_m(1.0e-4, ice_m=0.0, k_snow=1e-9, limit_m=5.0) == 5.0


@pytest.mark.parametrize("call", [
    lambda: il.transmission(-1.0),
    lambda: il.transmission(1.0, -0.1),
    lambda: il.transmission(float("nan")),
    lambda: il.transmission(1.0, kind="slush"),
    lambda: il.thickest_ice_m(0.0),
    lambda: il.thickest_ice_m(1.5),
    lambda: il.deepest_snow_m(0.0),
    lambda: il.albedo(-1.0),
])
def test_bad_input_is_refused(call):
    with pytest.raises(ValueError):
        call()
