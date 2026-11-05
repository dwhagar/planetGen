# tests/test_nitrogen_ice.py

"""Nitrogen ice from the disk model (GEN.1218a, GEN.1225, `physics/nitrogen_ice.py`) and the ethane an icy body can draw on (`ice_species.ethane_share_of_other`)."""

from types import SimpleNamespace

import pytest

from planetgen.physics import bulk_draw, ice_species, nitrogen_ice, surface_liquid_world as slw
from planetgen.util import draw


def _icy(temperature_k, snow_line_au=2.7):
    """An icy body whose formation zone is at `temperature_k` (the water line, 182 K, at the snow line)."""
    distance = snow_line_au * (182.0 / temperature_k) ** 2
    return SimpleNamespace(radius=1000.0, mass=None, star=SimpleNamespace(snow_line_au=snow_line_au), is_moon=False, distance=distance)


def test_nitrogen_condenses_at_about_29_k_in_the_disk():
    assert nitrogen_ice.condensation_temperature_k() == pytest.approx(28.7, abs=0.5)
    assert nitrogen_ice.DISK_N2_PER_H2 == pytest.approx(6.76e-5, rel=0.01)


def test_only_a_zone_colder_than_the_line_is_cold():
    assert nitrogen_ice.zone_is_cold(20.0) and not nitrogen_ice.zone_is_cold(45.0)


def test_the_trapping_range_gives_one_ten_thousandth_to_one_and_a_half_thousandths_of_the_water():
    low = nitrogen_ice.percent_of_water(nitrogen_ice.TRAPPING_EFFICIENCY_RANGE[0])
    high = nitrogen_ice.percent_of_water(nitrogen_ice.TRAPPING_EFFICIENCY_RANGE[1])
    assert low == pytest.approx(1.07e-2, rel=0.01) and high == pytest.approx(0.139, rel=0.01)


def test_a_share_over_one_is_refused():
    with pytest.raises(ValueError):
        nitrogen_ice.percent_of_water(1.5)


def test_a_cold_zone_adds_nitrogen_ice_in_the_range_and_a_warm_one_adds_none():
    base = {"CO": 1.4, "CO2": 16.5}
    draw.set_run_seed(3)
    warm = bulk_draw._with_nitrogen_ice(_icy(45.0), base)
    assert warm == base
    for seed in range(5):
        draw.set_run_seed(seed)
        cold = bulk_draw._with_nitrogen_ice(_icy(20.0), base)
        assert cold["CO"] == 1.4 and nitrogen_ice.percent_of_water(0.001) <= cold["N2"] <= nitrogen_ice.percent_of_water(0.013)


def test_nitrogen_percent_becomes_a_nitrogen_share_of_the_ice():
    shares = ice_species.shares_of_percent({"CO2": 16.5, "N2": 0.05})
    assert shares["N2"] == pytest.approx(0.05 / 16.55)


def test_a_share_of_the_unlisted_ices_is_ethane():
    assert 0.05 < ice_species.ethane_share_of_other() < 0.3


def test_ethane_in_the_unlisted_ices_supplies_an_ethane_reservoir():
    body = SimpleNamespace(ice_other_mass_fraction=0.02, water_mass_fraction=0.3)
    reservoir = slw._reservoir(body)
    assert reservoir[slw.ETHANE_LIQUID] == pytest.approx(0.02 * ice_species.ethane_share_of_other())
    assert slw.ETHANE_LIQUID not in slw._reservoir(SimpleNamespace(ice_other_mass_fraction=0.0))
