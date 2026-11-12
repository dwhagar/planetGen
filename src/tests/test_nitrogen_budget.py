# tests/test_nitrogen_budget.py

"""GEN.1584: where a rocky body's nitrogen goes before its air: the core's share by the redox, and the share of the degassed nitrogen that stays gas."""

import math

import pytest

from planetgen.physics import mixture_rocky, nitrogen_budget as nb, surface_pressure


def test_the_partition_coefficient_follows_the_table_it_was_taken_from():
    assert nb.metal_silicate_d(nb.REDUCED_DIW) == pytest.approx(1.0)
    assert nb.metal_silicate_d(nb.OXIDISED_DIW) == pytest.approx(39.0)
    assert nb.metal_silicate_d(-2.5) == pytest.approx(39.0 ** 0.5, rel=1e-9)                        # halfway between the ends, geometrically
    assert nb.metal_silicate_d(-20.0) == 1.0 and nb.metal_silicate_d(20.0) == nb.D_N_MAX            # held inside the table's range
    assert [nb.metal_silicate_d(x) for x in (-3.0, -2.0, -1.0, 0.0)] == sorted(nb.metal_silicate_d(x) for x in (-3.0, -2.0, -1.0, 0.0))     # more oxidised, more in the metal


def test_the_earth_twins_core_formation_sits_at_the_measured_partition():
    d = nb.metal_silicate_d(nb.REFERENCE_DIW + nb.CORE_FORMATION_OFFSET)
    assert 5.0 < d < 31.0                                   # 10 to 31 measured (Suer 2023), nominal 10 (Johansen 2023), 20 +/- 10 (Johnson and Goldblatt 2015)
    assert nb.mantle_share(nb.REFERENCE_DIW, 0.47) == pytest.approx(1.0 / (1.0 + d * 0.47))


def test_the_earth_twin_and_every_more_reduced_body_keep_the_inventory_and_an_oxidised_one_loses_to_the_metal():
    assert nb.mantle_retention(nb.REFERENCE_DIW, 0.47) == 1.0
    assert nb.mantle_retention(-1.5, 0.47) == 1.0
    assert nb.mantle_retention(1.5, 0.47) < nb.mantle_retention(1.0, 0.47) < 1.0
    assert nb.mantle_retention(1.5, 0.9) < nb.mantle_retention(1.5, 0.2)                    # a bigger core takes more
    assert nb.mantle_retention(None, 0.47) == 1.0 and nb.mantle_retention(1.5, None) == 1.0 and nb.mantle_retention(1.5, 0.0) == 1.0


def test_the_inventory_of_a_body_with_a_core_is_cut_by_the_retention_and_one_without_is_not(monkeypatch):
    subject = {"age_gy": 4.5, "metallicity_fe_h": 0.0}
    assert mixture_rocky._metal_retention(subject, 2.0) == 1.0                              # no metal or rock share: unchanged
    subject.update(metal_mass_fraction=0.32, rock_mass_fraction=0.68)
    assert mixture_rocky._metal_retention(subject, 2.0) == pytest.approx(nb.mantle_retention(2.0, 0.32 / 0.68))
    assert mixture_rocky._metal_retention(subject, 2.0) < 1.0
    assert mixture_rocky._metal_retention(subject, nb.REFERENCE_DIW) == 1.0


def test_nitrogen_is_drawn_between_the_published_atmospheric_shares_and_carbon_is_not():
    low, high, shape = surface_pressure.DEGAS_EFFICIENCY_N
    assert shape == "log" and 0.3 <= low < high <= 0.55
    assert 0.11 / 0.36 <= math.sqrt(low * high) <= 0.19 / 0.36 * 1.1                         # the middle sits in the range 4 of 35 and 4 of 21 give over the Earth twin's 0.36
    assert "N" not in surface_pressure.DEGAS_EFFICIENCY and surface_pressure.DEGAS_EFFICIENCY["C"] == 1.0
    efficiencies = surface_pressure.draw_efficiencies({})
    assert low <= efficiencies["N"] <= high


def test_the_nitrogen_the_delivered_ratio_and_the_degassing_give_at_the_middle_is_the_old_surface_nitrogen():
    low, high, _ = mixture_rocky.INVENTORY["nitrogen_per_carbon"]
    elow, ehigh, _ = surface_pressure.DEGAS_EFFICIENCY_N
    assert math.sqrt(low * high) * math.sqrt(elow * ehigh) == pytest.approx(0.022, rel=0.1)  # the old log-middle of the nitrogen per carbon, so the Earth twin keeps its 78 kPa


def _recipe(row):
    class Group:
        name = row
    class Recipe:
        groups = {"b": Group}
        values = {}
    return Recipe


def test_the_nitrogen_of_a_wet_class_follows_the_dry_carbon_and_only_a_share_of_the_rest():
    """GEN.1584 B (stand-in NB5): up to the dry class's log-middle all the carbon brings nitrogen, above it a drawn share."""
    low, high, _ = mixture_rocky.INVENTORY["co2_dry_wt_percent"]
    reference = math.sqrt(low * high) / mixture_rocky.AMU_G_MOL["CO2"]
    share_low, share_high, shape = mixture_rocky.INVENTORY["nitrogen_with_wet_carbon"]
    assert shape == "log" and 0.0 < share_low < share_high <= 1.0
    dry = mixture_rocky._nitrogen_carbon(_recipe("b1"), 100.0 * reference, {})
    assert dry == 100.0 * reference                                                      # the dry class is not touched
    for row in ("b2", "b3", "b4"):
        assert mixture_rocky._nitrogen_carbon(_recipe(row), 0.5 * reference, {}) == 0.5 * reference
        constants = {}
        counted = mixture_rocky._nitrogen_carbon(_recipe(row), 100.0 * reference, constants)
        assert counted == pytest.approx(reference + constants["nitrogen_with_wet_carbon"] * 99.0 * reference)
        assert share_low <= constants["nitrogen_with_wet_carbon"] <= share_high
        assert counted < 100.0 * reference
