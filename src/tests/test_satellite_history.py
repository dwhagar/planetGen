# tests/test_satellite_history.py

"""
GEN.526: a satellite's young and intermediate stars come from its star formation history (`galaxy/satellite_history.py`),
not from a typed share: the history is set by the row's mean stellar age, the gas by the Schmidt-Kennicutt law, and a
satellite whose gas is under the star formation threshold forms no stars now.
"""

import math

import pytest

from planetgen import tuning
from planetgen.galaxy import component_field, satellite_history
from planetgen.generation import bright_stars as brightStars, component_pass
from tests.reference_satellites import MILKY_WAY_SATELLITES

SPAN = tuning.SATELLITE_HISTORY_SPAN_GY


def _row(name):
    return next(row for row in MILKY_WAY_SATELLITES if row.name == name)


def _history(name):
    row = _row(name)
    return satellite_history.history(row.stellar_mass_sol, row.half_light_radius_pc, row.age_gy)


@pytest.mark.parametrize("mean_age", [0.5, 3.0, 5.0, 6.5, 9.0, 12.0])
def test_the_history_has_the_mean_age_it_was_set_by(mean_age):
    rate = satellite_history.rate_for_mean_age(mean_age, SPAN)
    assert satellite_history.mean_age_gy(rate, SPAN) == pytest.approx(mean_age, abs=1e-6)


def test_a_flat_history_has_half_the_span_as_its_mean_age_and_its_shares_are_the_windows():
    assert satellite_history.mean_age_gy(0.0, SPAN) == pytest.approx(SPAN / 2)
    found = satellite_history.history(1.0e9, 1000.0, SPAN / 2)
    assert found.young_share == pytest.approx(0.1 / SPAN, rel=1e-6) or not found.forming
    assert found.young_share + found.intermediate_share + found.old_share == pytest.approx(1.0)
    assert found.intermediate_share == pytest.approx(2.9 / SPAN, rel=1e-4)


def test_the_shares_add_up_to_one_for_every_satellite():
    for row in MILKY_WAY_SATELLITES:
        found = _history(row.name)
        assert found.young_share >= 0.0 and found.intermediate_share >= 0.0 and found.old_share > 0.0
        assert found.young_share + found.intermediate_share + found.old_share == pytest.approx(1.0)


def test_a_younger_mean_age_puts_more_stars_in_the_recent_windows():
    shares = [satellite_history.history(1.0e9, 1500.0, age) for age in (4.0, 6.0, 8.0, 11.0)]
    assert [s.intermediate_share for s in shares] == sorted((s.intermediate_share for s in shares), reverse=True)


def test_the_large_cloud_forms_stars_at_about_the_rate_and_with_about_the_gas_observed():
    """The targets, not inputs: the Large Cloud's star formation has averaged roughly 0.2 solar masses a year since 5 Gy
    ago (Harris & Zaritsky 2009, AJ 138, 1243, arXiv:0908.1422, read) and holds some 4e8 solar masses of atomic gas
    [recalled, to verify]; the history gives the rate and the gas from the stars' mass, size and mean age alone."""
    found = _history("Large Magellanic Cloud")
    assert found.forming
    assert 0.1 < found.present_rate_msun_yr < 0.5
    assert 2.0e8 < found.gas_mass_solar < 8.0e8
    assert 0.005 < found.young_share < 0.03 and 0.2 < found.intermediate_share < 0.5


def test_the_small_cloud_forms_stars_too_and_the_dwarf_spheroidals_do_not():
    assert _history("Small Magellanic Cloud").forming and _history("Small Magellanic Cloud").young_share > 0.0
    for name in ("Sculptor dwarf", "Draco dwarf", "Ursa Minor dwarf", "Sextans dwarf", "Carina dwarf"):
        found = _history(name)
        assert not found.forming and found.young_share == 0.0, name
    for name in ("Sculptor dwarf", "Draco dwarf", "Ursa Minor dwarf", "Sextans dwarf"):
        assert _history(name).intermediate_share < 0.01       # old populations: almost nothing in the last 3 Gy


def test_the_gas_is_the_one_the_law_needs_and_the_threshold_decides():
    found = satellite_history.history(2.0e9, 2000.0, 5.0)
    area_kpc2 = math.pi * 2.0 ** 2
    rate_per_area = found.present_rate_msun_yr / area_kpc2
    assert rate_per_area == pytest.approx(
        tuning.STAR_FORMATION_LAW_COEFFICIENT * found.gas_surface_density ** tuning.STAR_FORMATION_LAW_EXPONENT)
    assert found.gas_mass_solar == pytest.approx(found.gas_surface_density * math.pi * 2000.0 ** 2)
    assert found.forming == (found.gas_surface_density >= tuning.SATELLITE_STAR_FORMATION_THRESHOLD_MSUN_PC2)
    spread = satellite_history.history(2.0e9, 20000.0, 5.0)         # the same stars over a hundred times the area
    assert spread.gas_surface_density < found.gas_surface_density
    assert not spread.forming and spread.young_share == 0.0 and spread.intermediate_share > 0.0


def test_a_field_splits_a_satellite_by_its_history_and_keeps_the_density_as_it_was():
    from tests.test_satellites import E1, EDGE, SHAPE, _by_name, _placed
    found = _placed()
    field = component_field.ComponentField(found, EDGE, E1)
    cloud = _by_name(found)["Large Magellanic Cloud"]
    densities = field.densities(cloud.center_pc)
    total = sum(densities[name] for name in component_field.SATELLITE_POPULATIONS)
    assert total == pytest.approx(field.relative_density(cloud.center_pc))
    shares = component_field.population_shares(cloud)
    for name in component_field.SATELLITE_POPULATIONS:
        assert densities[name] / total == pytest.approx(shares[name])
    assert component_field.population_shares(
        next(c for c in found if c.kind == "satellite"))["satellite_young"] >= 0.0
    del SHAPE


def test_the_clouds_place_massive_stars_at_the_default_floor_now():
    """The old satellite window (0.5 to 13 Gy) held no star this bright; the history's young share does."""
    fractions = brightStars.band_fractions(9000.0)
    assert fractions["satellite_young"] > 0.0 and fractions["satellite"] == 0.0
    from tests.test_satellites import E1, EDGE, _by_name, _placed
    cloud = _by_name(_placed())["Large Magellanic Cloud"]
    expected = component_pass.expected_stars([cloud], fractions)
    assert expected > 1.0e4
    dwarf = _by_name(_placed())["Draco dwarf"]
    assert component_pass.expected_stars([dwarf], fractions) < 0.01 * expected      # an old dwarf's last 3 Gy only
    del E1, EDGE
