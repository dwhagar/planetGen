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


def test_hydrogen_and_helium_escape_from_an_earth_like_body_and_the_heavy_gases_stay():
    result = escape.apply(Subject("M", EARTH), made({"H2": 0.2, "He": 0.01, "N2": 0.5, "CO2": 0.1, "H2O": 0.19}))
    assert "H2" not in result.mole_fractions
    assert result.mole_fractions["He"] < 1e-4                      # Earth keeps a few ppm of helium
    assert result.mole_fractions["N2"] == pytest.approx(0.5 / 0.79, rel=1e-3)
    assert sum(result.mole_fractions.values()) == pytest.approx(1.0)
    assert result.constants["escape_kept"]["H2"] < 1e-3


def test_the_jeans_parameter_of_hydrogen_on_earth_is_near_fifteen():
    exobase = escape.exobase_temperature_k(Subject("M", EARTH), {"N2": 1.0})
    assert exobase == pytest.approx(1020.0, rel=0.01)
    assert escape.jeans_parameter("H2", 11.2, exobase) == pytest.approx(14.9, abs=0.5)
    assert escape.jeans_parameter("N2", 11.2, exobase) > 200


def test_a_carbon_dioxide_air_has_a_cold_exobase_and_keeps_its_gas_on_a_mars_sized_body():
    mars = dict(EARTH, v_esc_kms=5.03, teq_k=210.0, mass_earth=0.107)
    result = escape.apply(Subject("K", mars), made({"CO2": 0.95, "N2": 0.03, "H2": 0.02}, "mars"))
    assert result.mole_fractions["CO2"] > 0.95 and "H2" not in result.mole_fractions


def test_a_small_warm_body_loses_everything_it_would_outgas_and_is_airless():
    for seed in range(40):
        with draw.bound(seed):
            result = md.draw_mixture(Subject("B", SMALL_WARM))
        for formula in result.mole_fractions:
            assert escape.jeans_parameter(formula, 3.54, result.constants["exobase_k"]) > 20.0
        assert "H2" not in result.mole_fractions and "CO" not in result.mole_fractions
        assert result.family == "airless" or sum(result.mole_fractions.values()) == pytest.approx(1.0)


def test_the_share_kept_is_smooth_falls_with_the_parameter_and_older_bodies_need_more():
    shares = [escape.kept_share(value, 4.5) for value in (10.0, 30.0, 54.0, 80.0, 150.0)]
    assert shares == sorted(shares) and shares[2] == pytest.approx(0.5) and shares[0] < 1e-3 and shares[-1] > 0.999
    assert escape.kept_share(54.0, 9.0) < 0.5 < escape.kept_share(54.0, 2.0)


def test_families_without_a_loss_filter_and_airless_bodies_are_untouched():
    giant = made({"H2": 0.9, "He": 0.1}, "giant")
    assert escape.apply(Subject("J", EARTH), giant) is giant
    lava = made({"Na": 1.0}, "lava")
    assert escape.apply(Subject("B", EARTH), lava) is lava
    airless = made({}, "airless")
    assert escape.apply(Subject("A", EARTH), airless) is airless


def test_a_missing_escape_speed_is_refused_by_name():
    quantities = {k: v for k, v in EARTH.items() if k != "v_esc_kms"}
    with pytest.raises(MixtureError, match="v_esc_kms"):
        escape.apply(Subject("M", quantities), made({"N2": 1.0}))


def test_the_escape_inputs_are_required_inputs_of_a_rocky_class():
    assert {"v_esc_kms", "teq_k", "age_gy"} <= set(md.required_inputs("M"))
