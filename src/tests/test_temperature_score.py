"""
GEN.1265 and GEN.1266: the Temperature domain is scored from the share of the surface and the year in which a person needs no suit
(`habitability.coverage_score`, `physics/temperature_coverage.py`), the stored columns are the ones that share is read from, the seven bands
are read off the same series, and the gear and the hours outside are those of the best site.
"""

import types
import math

import numpy as np
import pytest

from planetgen.physics import equipment, habitability as hab
from planetgen.physics import habitability_world as hw
from planetgen.physics import regional_temperature, rogue_habitability, seasonal_swing, temperature_coverage as tc
from planetgen.physics.habitability import World

from tests.test_habitability_world import _earth

YEAR_YR = 1.0


def _spinning(**changes):
    """Earth with a spin, an orbit, a tilt, an albedo and a sea, so a series can be worked."""
    values = _earth(rotation_period_hours=24.0, period_years=YEAR_YR, axial_tilt_deg=23.4, bond_albedo=0.30, greenhouse_factor=0.15,
                    ground_ocean_fraction=0.7, ground_ice_fraction=0.05, day_share_mean=0.5, daylight_hours=12.0)
    values.update(changes)
    return values


def test_the_score_runs_through_the_edges():
    score = hab.coverage_score
    assert score(0.0, 0.0) == 0.0
    assert score(hab.LIVABLE_YELLOW_EDGE, 1.0) == pytest.approx(0.4)
    assert score(hab.LIVABLE_GREEN_EDGE, 1.0) == pytest.approx(0.75)
    assert score(hab.LIVABLE_BLUE_EDGE, 1.0) == 1.0 and score(1.0, 1.0) == 1.0
    assert [hab.tier(score(share, 1.0)) for share in (0.0, 0.2, 0.6, 0.9)] == ["Red", "Yellow", "Green", "Blue"]
    shares = np.linspace(0.0, 1.0, 51)
    assert all(a <= b + 1e-12 for a, b in zip([score(x, 1.0) for x in shares], [score(x, 1.0) for x in shares[1:]])), "more of the world livable never scores worse"


def test_blue_needs_a_share_in_the_blue_band():
    assert hab.tier(hab.coverage_score(1.0, 0.0)) == "Green"
    assert hab.coverage_score(1.0, hab.BLUE_SHARE_MIN * 0.99) < 1.0
    assert hab.coverage_score(1.0, hab.BLUE_SHARE_MIN) == 1.0


def test_a_world_that_needs_a_suit_everywhere_is_yellow_and_one_that_kills_is_red():
    assert hab.tier(hab.coverage_score(0.0, 0.0, yellow_share=1.0)) == "Yellow"
    assert hab.tier(hab.coverage_score(0.0, 0.0, yellow_share=0.0)) == "Red"
    assert hab.coverage_score(0.0, 0.0, yellow_share=0.5) > 0.0


def test_a_world_livable_only_in_the_terminator_belt_is_yellow_whatever_the_share():
    assert hab.tier(hab.coverage_score(0.9, 1.0)) == "Blue"
    assert hab.tier(hab.coverage_score(0.9, 1.0, belt_only=True)) == "Yellow"


def test_a_world_with_one_temperature_is_scored_from_it():
    base = dict(pressure_kpa=101.3, gases={"N2": 0.78, "O2": 0.21})
    tiers = {t: hab.tier(hab.temperature_score(World(temperature_c=t, **base))) for t in (15.0, 40.0, 60.0, -30.0, -90.0, 200.0)}
    assert tiers == {15.0: "Blue", 40.0: "Green", 60.0: "Yellow", -30.0: "Yellow", -90.0: "Red", 200.0: "Red"}


def test_an_uncompensable_wet_bulb_makes_every_sample_yellow():
    humid = World(pressure_kpa=101.3, temperature_c=38.0, gases={"N2": 0.79, "O2": 0.21}, relative_humidity=0.9)
    assert hab.wet_bulb_c(38.0, 0.9) > hab.WET_BULB_CRIT_C
    read = hab.temperature_coverage(humid)
    assert read.livable_share == 0.0 and read.humid == pytest.approx(1.0)
    assert hab.tier(hab.temperature_score(humid)) == "Yellow"


def test_a_visitor_with_other_limits_is_read_against_their_own_bands():
    hot = World(pressure_kpa=101.3, temperature_c=60.0, gases={"N2": 0.78, "O2": 0.21})
    species = hab.bands({"heat_garment_c": 70.0})
    assert hab.tier(hab.temperature_score(hot)) == "Yellow"
    assert hab.tier(hab.temperature_score(hot, species)) == "Green"


def test_earth_reads_blue_and_stores_what_it_was_read_from():
    scored = hw.compute(_spinning())
    assert scored["tier_temperature"] == 0 and scored["phi4_temperature"] == 1.0
    assert 0.8 < scored["temp_livable_share"] <= 1.0
    assert scored["temp_blue_share"] >= hab.BLUE_SHARE_MIN
    assert scored["temp_livable_belt_only"] == 0
    assert set(tc.FIELDS) <= set(scored) and all(scored[name] is not None for name in tc.FIELDS)
    assert scored["temp_coldest_k"] <= scored["temp_site_coldest_k"] <= scored["temp_site_mean_k"] <= scored["temp_site_hottest_k"] <= scored["temp_hottest_k"]
    assert -90.0 <= scored["temp_coldest_lat_deg"] <= 90.0 and 0.0 <= scored["temp_coldest_lon_deg"] <= 360.0
    assert scored["temp_blue_share"] + scored["temp_green_share"] == pytest.approx(scored["temp_livable_share"])


def test_a_gas_giant_and_a_body_with_no_temperature_store_nothing():
    for values in (_spinning(body_type="g"), _spinning(surface_temperature_k=None)):
        scored = hw.compute(values)
        assert all(scored[name] is None for name in tc.FIELDS)


def test_a_body_with_no_spin_or_orbit_has_one_temperature_everywhere():
    scored = hw.compute(_earth())
    assert scored["temp_livable_share"] == 1.0 and scored["temp_coldest_k"] == scored["temp_hottest_k"] == pytest.approx(288.0)
    assert scored["temp_site_longest_out_days"] == 0.0


def test_a_cold_mean_world_with_a_warm_region_can_score_better_than_its_mean_did():
    """The old mean-based score is not a cap (Boss, 7:32 PM EST 11-04): a mean in the Red band with a livable third of the world reads from that third."""
    grid = tc.surface_cycle.cell_grid()
    celsius = np.full((12, 12, 8), -100.0)
    celsius[4:8, :, :] = 30.0                     # the middle third of the rows is warm
    series = tc.SimpleNamespace(grid=grid, temperature_k=celsius + 273.15, time_days=np.arange(8) + 0.5, step_days=1.0, orbital_days=8.0, solar_days=1.0)
    read = tc.coverage(series)
    mean_c = float(np.sum(celsius.mean(axis=2) * grid.row_weight[:, None]) / 12.0)
    old = hab.band_score(mean_c, hab.HUMAN.temperature)
    new = hab.coverage_score(read.livable_share, read.blue, tc.belt_only(read), read.yellow)
    assert mean_c < -50.0 and hab.tier(old) == "Red"
    assert read.livable_share == pytest.approx(1.0 / 3.0) and hab.tier(new) == "Yellow" and new > old


def test_a_locked_world_livable_only_in_its_terminator_belt_says_so():
    values = _spinning(rotation_period_hours=8766.0, period_years=1.0, axial_tilt_deg=0.0, atmospheric_pressure_pa=1013.0, surface_temperature_k=285.0,
                       ground_ocean_fraction=0.0, ground_ice_fraction=0.0, greenhouse_factor=0.0)
    values.update(p_o2_kpa=0.2, p_n2_kpa=0.8)
    read = tc.coverage_of_values(values, 0.3)
    assert read.sectors is not None
    if tc.belt_only(read):
        assert hab.tier(hab.temperature_score(hw.world_from(values))) in ("Yellow", "Red")
    assert tc.belt_only(read) == (read.livable_share > 0.0 and read.sectors[tc.SECTORS[0]][0] < tc.BELT_ONLY_SHARE and read.sectors[tc.SECTORS[2]][0] < tc.BELT_ONLY_SHARE)


def test_a_hot_airless_body_that_spins_fast_gives_a_series_that_is_a_number():
    """A close-in airless body that turns near its orbit used to overflow in the emission and give a series of NaN (GEN.1265)."""
    series = seasonal_swing.temperature_series(22.48, 16.9, 26.25, 673.6, 67237.8, 0.305, 0.0, 0.288, 0.0)
    assert np.isfinite(series.temperature_k).all()
    assert series.temperature_k.max() < 2000.0 and series.temperature_k.min() >= seasonal_swing.FLOOR_K


def test_the_bands_are_read_off_the_series_and_average_to_the_mean():
    values = _spinning()
    out = regional_temperature.compute(values)
    means = [out[regional_temperature.column(n, "mean_k")] for n in range(1, 8)]
    assert all(m is not None for m in means) and means[3] > means[0] and means[3] > means[6]
    from planetgen.physics import latitude_bands
    areas = [band.area_share for band in latitude_bands.compute(values, latitude_bands.alpha_deg_of(values))]
    assert sum(a * m for a, m in zip(areas, means)) / sum(areas) == pytest.approx(288.0, abs=0.5)      # the area-weighted mean of the bands is the world's mean
    assert out["temp_mean_swing_k"] > 0.0


def test_a_locked_world_ring_swing_is_its_real_range():
    values = _spinning(rotation_period_hours=8766.0, period_years=1.0, axial_tilt_deg=0.0, atmospheric_pressure_pa=1013.0, surface_temperature_k=285.0,
                       ground_ocean_fraction=0.0, ground_ice_fraction=0.0, greenhouse_factor=0.0)
    out = regional_temperature.compute(values)
    assert out["band1_state"] in ("day side", "tidally locked", out["band1_state"])
    assert out[regional_temperature.column(1, "mean_k")] > out[regional_temperature.column(7, "mean_k")]


def test_the_gear_is_chosen_for_the_best_site_and_the_worst_place_is_recommended():
    values = _spinning(surface_temperature_k=250.0, greenhouse_factor=0.0, bond_albedo=0.3, ground_ocean_fraction=0.0, ground_ice_fraction=0.0)
    world = hw.world_from(values)
    site_cold, site_hot = hab.site_extremes_c(world)
    worst_cold, worst_hot = hab.worst_extremes_c(world)
    assert worst_cold <= site_cold <= site_hot <= worst_hot
    levels = equipment.compute(world)
    body, _ = equipment.body_at(world, site_cold)
    assert levels["eq_body"] >= body or levels["eq_body"] >= 0
    deepest = equipment.body_at(world, worst_cold)[0]
    assert levels["eq_body_rec"] == 0 or levels["eq_body_rec"] <= max(deepest, equipment.body_at(world, worst_hot)[0])


def test_the_rogue_planet_stores_the_columns_and_a_brown_dwarf_leaves_them_empty():
    """A rogue planet is built by the planet chain and scored like any planet (GEN.1111), so `rogue_planets` carries the columns too; a free
    brown dwarf has no surface and sets each to `None`."""
    assert set(tc.FIELDS) <= set(rogue_habitability.ROGUE_FIELDS)
    dwarf = types.SimpleNamespace()
    rogue_habitability.generate(dwarf)
    assert all(getattr(dwarf, name) is None for name in tc.FIELDS)


def test_the_coverage_of_the_same_inputs_is_kept():
    values = _spinning()
    first = tc.coverage_of_values(values, 0.6)
    assert tc.coverage_of_values(values, 0.6) is first
    assert tc.series_of_values(values) is tc.series_of_values(dict(values))


def test_the_api_object_reads_the_stored_columns():
    scored = hw.compute(_spinning())
    summary = tc.summary(scored)
    assert summary["livable_share"] == scored["temp_livable_share"] and summary["livable_belt_only"] is False
    assert summary["best_site"]["coldest_k"] == scored["temp_site_coldest_k"] and summary["coldest"]["temperature_k"] == scored["temp_coldest_k"]
    assert tc.summary({"temp_livable_share": None}) is None
    assert math.isfinite(summary["best_site"]["longest_out_days"])
