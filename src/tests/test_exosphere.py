"""The exosphere of an airless body (GEN.455, part B): its species follow the source composition, its columns the sources and losses of part 04 section 5, and a body with any air gets none."""

import math
import random

import pytest

from planetgen.physics import exosphere

MOON_ROCK = {"Na2O": 0.5, "K2O": 0.05, "CaO": 12.0, "MgO": 8.0, "FeO": 15.0, "SiO2": 45.0, "Al2O3": 15.0}
MERCURY_ROCK = {"Na2O": 3.0, "K2O": 0.2, "CaO": 5.0, "MgO": 30.0, "FeO": 2.0, "SiO2": 55.0, "Al2O3": 6.0}
SUN_K = 5772.0


def _inputs(gravity, radius_km, temperature, rock, distance_au, teff=SUN_K, ice=None, cover=0.0, wind=1.0, shielded=False, luminosity=1.0):
    flux = luminosity / distance_au ** 2
    return exosphere.Inputs(gravity, radius_km, temperature, rock, ice or {}, cover, teff, flux, distance_au, wind, shielded)


def _columns(inputs, seed=1):
    return {entry["species"]: entry for entry in exosphere.columns(inputs, random.Random(seed))}


def test_the_sodium_of_a_moon_like_and_a_mercury_like_body_sits_in_the_range_of_the_test_points():
    """The test points of part 04 section 5.2 are targets, never inputs: sodium 1e9 to 1e12 cm-2 over the Moon, Mercury, Io, Europa and the comets."""
    moon = [_columns(_inputs(1.62, 1737.0, 250.0, MOON_ROCK, 1.0), seed)["Na"]["column_cm2"] for seed in range(30)]
    mercury = [_columns(_inputs(3.7, 2440.0, 440.0, MERCURY_ROCK, 0.387, luminosity=1.0), seed)["Na"]["column_cm2"] for seed in range(30)]
    assert all(1e9 <= column <= 1e12 for column in moon + mercury)
    assert max(moon) < min(mercury)                        # more sodium in the rock and nearer the star: more in the air
    assert all(entry for entry in _columns(_inputs(1.62, 1737.0, 250.0, MOON_ROCK, 1.0)).values())


def test_every_species_is_one_a_source_supplies_from_the_surface_the_body_has():
    """No species outside the rock and the ices: no K from a rock with none, no water from ground with no ice, no ion-sputtered species when there is no wind."""
    no_potassium = dict(MOON_ROCK, K2O=0.0)
    found = _columns(_inputs(1.62, 1737.0, 250.0, no_potassium, 1.0))
    assert "K" not in found and "Na" in found
    assert set(found) <= {"Na", "K", "Ca", "Mg", "Fe", "SiO", "H2O", "CO2", "CH4", "NH3", "N2", "SO2", "CO", "HCN"}
    assert "H2O" not in found                                  # no ice on the ground
    icy = _columns(_inputs(1.31, 1560.0, 100.0, MOON_ROCK, 5.2, ice={"H2O": 1.0}, cover=0.9, luminosity=1.0))
    assert "H2O" in icy and icy["H2O"]["source"] == exosphere.SOURCES[0]
    only_rock = {entry["source"] for entry in _columns(_inputs(1.62, 1737.0, 250.0, MOON_ROCK, 1.0, wind=0.0)).values()}
    assert exosphere.SOURCES[3] not in only_rock and exosphere.SOURCES[0] not in only_rock


def test_a_rock_with_no_oxides_to_read_and_a_dark_star_give_nothing():
    assert exosphere.atom_shares({}) == {}
    assert exosphere.columns(_inputs(1.62, 1737.0, 250.0, {}, 1.0), random.Random(1)) == []


def test_the_atoms_share_of_an_element_is_its_moles_over_the_atoms_of_the_rock():
    shares = exosphere.atom_shares({"Na2O": 1.0, "SiO2": 99.0})
    moles_na2o, moles_sio2 = 1.0 / 61.979, 99.0 / 60.084
    assert shares["Na"] == pytest.approx(2.0 * moles_na2o / (3.0 * moles_na2o + 3.0 * moles_sio2))
    assert shares["SiO"] == pytest.approx(moles_sio2 / (3.0 * moles_na2o + 3.0 * moles_sio2))


def test_sodium_follows_the_rate_the_star_ionises_it_at():
    for temperature, rate in exosphere.NA_IONISATION_RATE_1AU:
        assert exosphere.na_ionization_rate_1au(temperature) == pytest.approx(rate)
    assert exosphere.na_ionization_rate_1au(5000.0) < exosphere.na_ionization_rate_1au(6000.0)
    # the lifetime is 1 / rate and scales as the inverse of the light: halving the light doubles it
    near = exosphere._ionization_lifetime_s("Na", _inputs(1.62, 1737.0, 250.0, MOON_ROCK, 1.0))
    far = exosphere._ionization_lifetime_s("Na", _inputs(1.62, 1737.0, 250.0, MOON_ROCK, math.sqrt(2.0)))
    assert far == pytest.approx(2.0 * near)


def test_a_body_that_cannot_hold_an_atom_records_none_and_a_wind_cut_by_a_field_gives_less():
    assert exosphere.columns(_inputs(0.001, 10.0, 200.0, MOON_ROCK, 2.5), random.Random(1)) == []        # the scale height is over the radius: the atoms leave
    open_wind = _columns(_inputs(3.7, 2440.0, 440.0, MERCURY_ROCK, 0.387))
    shielded = _columns(_inputs(3.7, 2440.0, 440.0, MERCURY_ROCK, 0.387, shielded=True))
    assert shielded["Mg"]["column_cm2"] < open_wind["Mg"]["column_cm2"]


def test_the_escape_shares_follow_the_speed_laws():
    na = exosphere.mass_per_atom_kg("Na")
    assert exosphere.weibull_escape_share(2380.0) < 1e-2 and exosphere.weibull_escape_share(4250.0) < 1e-3     # the Moon, Mercury (part 04)
    assert exosphere.maxwell_escape_share(na, 3500.0, 100.0) == pytest.approx(1.0, abs=0.01)
    assert exosphere.maxwell_escape_share(na, 3500.0, 20000.0) < 1e-60
    assert 0.9 < exosphere.sputter_escape_share(na, 5.0, 2380.0) <= 1.0                  # 90 percent to all at the Moon (part 04)
    assert exosphere.sputter_escape_share(na, 5.0, 4250.0) < exosphere.sputter_escape_share(na, 5.0, 2380.0)    # held better by Mercury
    assert exosphere.time_up_s(100.0, None) == 100.0 and exosphere.time_up_s(100.0, 1.0e9) == pytest.approx(100.0, rel=1e-6)
    assert exosphere.time_up_s(1.0e6, 10.0) == pytest.approx(10.0)


def test_a_thermal_column_never_passes_the_edge_where_the_gas_is_an_atmosphere():
    """A frost of carbon dioxide at 130 K would sublime as an atmosphere: the column is held at the collisional one."""
    found = _columns(_inputs(1.24, 2410.0, 130.0, MOON_ROCK, 5.2, ice={"CO2": 0.01, "H2O": 0.99}, cover=0.6))
    assert found["CO2"]["column_cm2"] <= exosphere.COLLISIONAL_COLUMN_CM2 * 1.0000001


def test_the_draws_are_the_bodys_own_and_a_body_with_any_air_has_no_exosphere():
    inputs = _inputs(1.62, 1737.0, 250.0, MOON_ROCK, 1.0)
    assert _columns(inputs, 4) == _columns(inputs, 4) and _columns(inputs, 4) != _columns(inputs, 5)

    class Body:
        body_type = "t"
        atmospheric_pressure = 0.0
        ground_ice_fraction = 0.0

    body = Body()
    assert exosphere.has_exosphere(body)
    for pressure in (1.0e-12, 1.0, 1.0e5):
        body.atmospheric_pressure = pressure
        assert not exosphere.has_exosphere(body)
    body.atmospheric_pressure = None
    assert not exosphere.has_exosphere(body)
    body.atmospheric_pressure, body.body_type = 0.0, "g"
    assert not exosphere.has_exosphere(body)
    body.body_type, body.ground_ice_fraction = "t", None
    assert not exosphere.has_exosphere(body)               # no solid ground

