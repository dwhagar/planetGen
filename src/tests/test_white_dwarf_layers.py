"""GEN.743: the layers of a helium-core and an oxygen-neon white dwarf (`physics.white_dwarf_layers`)."""

import pytest

from planetgen import tuning
from planetgen.physics import white_dwarf_layers as layers


@pytest.mark.parametrize("mass", [0.15, 0.2, 0.3, 0.45])
@pytest.mark.parametrize("envelope", ["thick", "thin"])
def test_the_helium_core_layers_sum_to_one_and_hold_no_carbon_oxygen_or_neon(mass, envelope):
    found = layers.helium_core(mass, envelope)
    assert sum(found.values()) == pytest.approx(1.0)
    assert found["carbon"] == found["oxygen"] == found["neon"] == found["trace"] == 0.0


def test_a_lighter_helium_core_star_has_the_thicker_hydrogen_envelope():
    assert layers.helium_core(0.17)["hydrogen"] > layers.helium_core(0.3)["hydrogen"] > layers.helium_core(0.44)["hydrogen"]


def test_the_helium_core_envelope_runs_between_its_two_points_and_is_held_outside_them():
    (low, high_fraction), (high, low_fraction) = tuning.WD_HELIUM_CORE_HYDROGEN_FRACTIONS
    assert layers.helium_core(low)["hydrogen"] == pytest.approx(high_fraction)
    assert layers.helium_core(high)["hydrogen"] == pytest.approx(low_fraction)
    assert layers.helium_core(0.1)["hydrogen"] == pytest.approx(high_fraction)
    assert layers.helium_core(0.5)["hydrogen"] == pytest.approx(low_fraction)


def test_the_thin_helium_core_star_holds_the_sequences_thin_hydrogen_case():
    assert layers.helium_core(0.3, "thin")["hydrogen"] == pytest.approx(tuning.WD_THIN_HYDROGEN_FRACTION)


@pytest.mark.parametrize("mass", [1.06, 1.10, 1.16, 1.2, 1.29, 1.33])
@pytest.mark.parametrize("envelope", ["thick", "thin"])
def test_the_oxygen_neon_layers_sum_to_one(mass, envelope):
    found = layers.oxygen_neon(mass, envelope)
    assert sum(found.values()) == pytest.approx(1.0)
    assert found["carbon"] == 0.0


def test_the_oxygen_neon_core_is_a_little_over_half_oxygen_and_about_a_third_neon():
    found = layers.oxygen_neon(1.16)
    core = found["oxygen"] + found["neon"] + found["trace"]
    assert found["oxygen"] / core == pytest.approx(0.55) and found["neon"] / core == pytest.approx(0.30)


def test_the_oxygen_neon_helium_layer_is_the_papers_mass_at_its_models_and_thinner_for_the_heavier_ones():
    for mass, helium_mass in tuning.WD_ONE_HELIUM_MASS_BY_WD_MASS_SOL:
        assert layers.oxygen_neon(mass)["helium"] * mass == pytest.approx(helium_mass)
    assert layers.oxygen_neon(1.12)["helium"] > layers.oxygen_neon(1.2)["helium"] > layers.oxygen_neon(1.3)["helium"]


def test_the_oxygen_neon_hydrogen_envelope_is_the_papers_maximum_for_a_hydrogen_atmosphere():
    assert layers.oxygen_neon(1.2)["hydrogen"] * 1.2 == pytest.approx(tuning.WD_ONE_HYDROGEN_MASS_SOL)
    assert layers.oxygen_neon(1.2, "thin")["hydrogen"] == pytest.approx(tuning.WD_THIN_HYDROGEN_FRACTION)
