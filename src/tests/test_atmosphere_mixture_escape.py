# tests/test_atmosphere_mixture_escape.py

"""Thermal escape of the light gases from a small body's drawn mixture (GEN.693)."""

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import mixture_draw as md, mixture_escape as escape
from planetgen.physics.mixture_base import MixtureDraw, MixtureError
from planetgen.util import draw

EARTH = dict(mass_earth=1.0, ts_k=288.0, pressure_bar=1.0, age_gy=4.5, long_lived_geology=True, ocean_present=True,
             snow_ratio=0.5, zone="e", host_type="G", is_moon=False, insolation=1.0, luminosity_solar=1.0, water_lost=False,
             runaway_ratio=0.5, teq_k=255.0, xuv_index=1.0, v_esc_kms=11.2, radius_km=6371.0, host_mass_solar=1.0)

SMALL_WARM = dict(EARTH, mass_earth=0.0327, ts_k=310.0, pressure_bar=0.198, v_esc_kms=3.54, radius_km=2100.0, teq_k=270.0, xuv_index=3.0,
                  long_lived_geology=False, ocean_present=False)


def made(fractions, family="rocky"):
    return MixtureDraw(mole_fractions=fractions, family=family, drawn={}, redox_diw=0.0, temperature_k=300.0, pressure_bar=1.0, inventory={})


def test_hydrogen_escapes_from_an_earth_like_body_and_the_heavy_gases_stay():
    air = {"H2": 0.002, "He": 0.001, "N2": 0.78, "O2": 0.21, "Ar": 0.007, "CO2": 4e-4}
    result = escape.apply(Subject("terrestrial", dict(EARTH, xuv_flux_earth=1.0)), made(air))
    assert "H2" not in result.mole_fractions
    assert result.constants["escape_kept"]["H2"] < 1e-3
    assert result.mole_fractions["N2"] == pytest.approx(0.78 / 0.9974, rel=0.02)
    assert sum(result.mole_fractions.values()) == pytest.approx(1.0)
    assert result.constants["exobase_k"] == pytest.approx(1155.0, rel=0.08)       # the neutral exobase of [R2184]


def test_an_air_with_no_carbon_dioxide_has_a_hotter_thermosphere_and_still_keeps_its_nitrogen():
    """Without the CO2 that cools it the thermosphere of the Earth is 1,300 K and more; [R2184] find the same atmosphere
    never hydrodynamic and its heavy species lost only by the very small Jeans rates of their Fig. 13."""
    result = escape.apply(Subject("terrestrial", dict(EARTH, xuv_flux_earth=1.0)), made({"N2": 0.78, "O2": 0.21, "Ar": 0.01}))
    assert result.constants["exobase_k"] > 1250.0
    assert "escape_kept" not in result.constants or result.constants["escape_kept"].get("N2", 1.0) > 0.9


def test_a_hydrogen_rich_air_drags_the_heavy_gases_only_while_it_has_the_hydrogen():
    """The history follows the air that is left: the hydrogen leaves first, and what is left of the air is as heavy as the Earth's."""
    result = escape.apply(Subject("terrestrial", dict(EARTH, xuv_flux_earth=1.0)), made({"H2": 0.5, "H2O": 0.3, "N2": 0.2}))
    assert "H2" not in result.mole_fractions
    assert result.mole_fractions["N2"] + result.mole_fractions["H2O"] == pytest.approx(1.0)
    assert result.constants["escape_kept"]["H2"] < 1e-3


def test_a_cold_carbon_dioxide_air_keeps_the_hydrogen_a_hot_nitrogen_air_loses():
    """The carbon dioxide cools the thermosphere ([R2184: Fig. 12]): the old filter fixed the ratio of the exobase to the
    equilibrium temperature by a threshold at half carbon dioxide; the model has the rise fall smoothly with the carbon dioxide."""
    cold = escape.apply(Subject("terrestrial", EARTH), made({"H2": 0.1, "N2": 0.4, "CO2": 0.5}))
    hot = escape.apply(Subject("terrestrial", EARTH), made({"H2": 0.1, "N2": 0.9}))
    assert cold.constants["exobase_k"] < 400.0 < hot.constants["exobase_k"]
    assert cold.constants["escape_kept"].get("H2", 1.0) > hot.constants["escape_kept"].get("H2", 1.0)


def test_a_carbon_dioxide_air_has_a_cold_exobase_and_keeps_its_gas_on_a_mars_sized_body():
    mars = dict(EARTH, v_esc_kms=5.03, teq_k=210.0, mass_earth=0.107, radius_km=3390.0, pressure_bar=0.00636)
    result = escape.apply(Subject("terrestrial", mars), made({"CO2": 0.95, "N2": 0.03, "H2": 0.02}, "mars"))
    assert result.mole_fractions["CO2"] > 0.95 and "H2" not in result.mole_fractions
    assert 180.0 <= result.constants["exobase_k"] <= 260.0                       # MAVEN


def test_the_radius_is_taken_from_the_escape_speed_and_the_mass_when_the_subject_has_no_radius():
    mass, radius, pressure, age = escape.body_inputs(Subject("terrestrial", {k: v for k, v in EARTH.items() if k != "radius_km"}))
    assert radius == pytest.approx(6.371e6, rel=0.02)
    assert pressure == pytest.approx(1e5) and age == 4.5


def test_a_small_warm_body_loses_the_hydrogen_it_would_outgas():
    for seed in range(40):
        with draw.bound(seed):
            result = md.draw_mixture(Subject("terrestrial", SMALL_WARM))
        assert "H2" not in result.mole_fractions
        assert result.family == "airless" or sum(result.mole_fractions.values()) == pytest.approx(1.0)


def test_families_without_a_loss_filter_and_airless_bodies_are_untouched():
    giant = made({"H2": 0.9, "He": 0.1}, "giant")
    assert escape.apply(Subject("giant", EARTH), giant) is giant
    lava = made({"Na": 1.0}, "lava")
    assert escape.apply(Subject("terrestrial", EARTH), lava) is lava
    airless = made({}, "airless")
    assert escape.apply(Subject("terrestrial", EARTH), airless) is airless


def test_a_missing_escape_speed_and_radius_is_refused_by_name():
    quantities = {k: v for k, v in EARTH.items() if k not in ("v_esc_kms", "radius_km")}
    with pytest.raises(MixtureError, match="v_esc_kms"):
        escape.apply(Subject("terrestrial", quantities), made({"N2": 1.0}))


def test_the_escape_inputs_are_required_inputs_of_a_rocky_class():
    assert {"v_esc_kms", "teq_k", "age_gy"} <= set(md.required_inputs("terrestrial"))


def test_a_lifeless_oxygen_world_is_capped_at_95_percent_with_a_nitrogen_argon_floor():
    from planetgen.physics import mixture_oxygen as oxygen
    capped = oxygen.cap_oxygen({"O2": 0.99, "CO2": 0.01}, oxygen.ABIOTIC_O2_MAX)
    assert capped["O2"] == pytest.approx(0.945)
    assert capped["N2"] + capped.get("Ar", 0.0) == pytest.approx(0.005)
    assert sum(capped.values()) == pytest.approx(1.0)
    pure = oxygen.cap_oxygen({"O2": 1.0}, oxygen.ABIOTIC_O2_MAX)
    assert pure["O2"] == pytest.approx(0.95) and pure["N2"] == pytest.approx(0.05)
    inside = {"O2": 0.5, "N2": 0.4, "Ar": 0.1}
    assert oxygen.cap_oxygen(inside, oxygen.ABIOTIC_O2_MAX) is inside


def test_a_living_world_is_capped_at_35_percent_and_keeps_the_ratios_of_the_other_gases():
    from planetgen.physics import mixture_oxygen as oxygen
    capped = oxygen.cap_oxygen({"O2": 0.6, "N2": 0.3, "CO2": 0.1}, oxygen.BIOTIC_O2_MAX)
    assert capped["O2"] == pytest.approx(0.35)
    assert capped["N2"] / capped["CO2"] == pytest.approx(3.0)


def test_no_abiotic_oxygen_world_is_drawn_above_95_percent_oxygen():
    from tests.test_atmosphere_mixture_evolved import ABIOTIC, BASE
    top, airless = 0.0, 0
    for seed in range(300):
        with draw.bound(seed):
            result = md.draw_mixture(Subject("terrestrial", {**BASE, **ABIOTIC}))
        if result.family == "airless":    # TEST.401: the thermosphere escape of GEN.738 strips this M-dwarf Earth in rare draws (seed 34)
            airless += 1
            assert not result.mole_fractions
            continue
        top = max(top, result.mole_fractions.get("O2", 0.0))
        assert sum(result.mole_fractions.values()) == pytest.approx(1.0)
    assert airless <= 15
    assert 0.5 < top <= 0.95 + 1e-12
