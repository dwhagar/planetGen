# tests/test_ionization.py

"""The ionization ceiling of an air (GEN.1580): the Saha half-ionization temperature of the main gas at the air's pressure."""

import json

import pytest

from planetgen.physics import elements, gases, ionization


def test_the_ceiling_of_an_air_is_tens_of_thousands_of_kelvin_and_rises_with_the_pressure():
    nitrogen = [ionization.ceiling_k({"N2": pressure}) for pressure in (1.0e3, 1.0e5, 1.0e7, 1.0e9)]
    assert nitrogen == sorted(nitrogen)
    assert nitrogen[1] == pytest.approx(15365.0, rel=2e-3)             # 1 bar of nitrogen, the Saha value for 14.534 eV
    assert ionization.ceiling_k({"CO2": 1.0e5}) == pytest.approx(12492.0, rel=2e-3)
    assert ionization.ceiling_k({"N2": 1.0e9}) > 40000.0


def test_the_main_gas_names_the_atom_that_ionizes_first():
    assert ionization.main_gas_energy_ev({"CO2": 9.0e4, "N2": 1.0e4}) == pytest.approx(ionization.first_ionization_energy_ev("C"))
    assert ionization.main_gas_energy_ev({"N2": 9.0e4, "CO2": 1.0e4}) == pytest.approx(ionization.first_ionization_energy_ev("N"))
    assert ionization.main_gas_energy_ev({"HD": 1.0}) == pytest.approx(ionization.first_ionization_energy_ev("H"))


def test_every_gas_of_the_table_has_an_ionization_ceiling():
    for item in json.load(open(gases.__file__.replace("gases.py", "data/gases.json"), encoding="utf-8"))["gases"]:
        formula = item["formula"] if isinstance(item, dict) else item
        assert ionization.ceiling_k({formula: 1.0e5}) > 5000.0


def test_the_half_ionized_share_is_a_half_at_the_ceiling():
    ceiling = ionization.half_ionization_temperature_k(1.0e5, 14.534)
    assert ionization.ionized_share(ceiling, 1.0e5, 14.534) == pytest.approx(0.5, rel=1e-4)
    assert ionization.ionized_share(0.5 * ceiling, 1.0e5, 14.534) < 1e-2


def test_the_energies_are_nist_first_ionization_energies():
    assert ionization.first_ionization_energy_ev("H") == pytest.approx(13.5984346, rel=1e-7)
    assert ionization.first_ionization_energy_ev("D") == ionization.first_ionization_energy_ev("H")
    assert ionization.first_ionization_energy_ev("He") == pytest.approx(24.5873890, rel=1e-7)


def test_bad_inputs_are_refused():
    with pytest.raises(ValueError):
        ionization.ceiling_k({})
    with pytest.raises(ValueError):
        ionization.half_ionization_temperature_k(0.0, 14.5)
    with pytest.raises(KeyError):
        ionization.first_ionization_energy_ev("Uue")
