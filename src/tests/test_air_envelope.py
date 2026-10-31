"""The envelope share of a hydrogen-rich air (GEN.1226)."""

import math

import pytest

from planetgen.physics import air_envelope as ae, constants, planet_type

EARTH_KG, EARTH_KM = 5.972e24, 6371.0


def test_earth_air_is_not_an_envelope():
    assert ae.envelope_share({"N2": 78000.0, "O2": 21000.0}, 101325.0, EARTH_KG, EARTH_KM) == 0.0


def test_a_hydrogen_air_holds_its_true_mass_share():
    share = ae.envelope_share({"H2": 9e5, "He": 1e5}, 1e6, EARTH_KG, EARTH_KM)
    assert share == pytest.approx(ae.air_mass_kg(1e6, EARTH_KG, EARTH_KM) / EARTH_KG)
    assert 1e-6 < share < 1e-5                   # 10 bar on an Earth-sized body: 8.7e-6


def test_the_air_mass_is_the_weight_of_the_column():
    gravity = constants.G * EARTH_KG / (EARTH_KM * 1000.0) ** 2
    assert ae.air_mass_kg(101325.0, EARTH_KG, EARTH_KM) == pytest.approx(101325.0 * 4 * math.pi * (EARTH_KM * 1000.0) ** 2 / gravity)


def test_no_air_or_no_pressure_is_zero():
    assert ae.envelope_share({}, 0.0, EARTH_KG, EARTH_KM) == 0.0 and ae.envelope_share({"H2": 1.0}, None, EARTH_KG, EARTH_KM) == 0.0


def test_a_thin_hydrogen_air_stays_terrestrial_by_the_type():
    share = ae.envelope_share({"H2": 1e6}, 1e6, 5 * EARTH_KG, 8000.0)
    assert planet_type.type_of(share) == ("terrestrial", True)
