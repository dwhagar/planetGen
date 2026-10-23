"""
Proton range in a column (GEN.453, `planetgen.physics.proton_range`): the NIST PSTAR table against the values
docs/design/atmosphere-composition-parts/07-radiation-shielding-and-ground-sources.md section 2.5 prints, the inverse,
the Bragg sum for mixtures and the file's checks.
"""

import json

import pytest

from planetgen.physics import proton_range as pr
from planetgen.physics.elements import DataError

AIR_COLUMN_TO_MEV = ((1, 30), (3, 55), (10, 108), (22.5, 171), (50, 274), (100, 420), (300, 868), (1033, 2241))
"""The kinetic energy, MeV, that crosses a column of air (part 07 table 2.5, inverted from PSTAR)."""


@pytest.mark.parametrize("column, energy", AIR_COLUMN_TO_MEV)
def test_air_columns_need_the_pstar_energies(column, energy):
    assert pr.energy_mev(column) == pytest.approx(energy, rel=0.01)


@pytest.mark.parametrize("material, column, energy", [
    ("CO2", 21.0, 165.0), ("H", 22.5, 284.0), ("He", 22.5, 181.0), ("H2O liquid", 22.5, 184.0), ("SiO2", 22.5, 165.0),
    ("CH4", 100.0, 508.0)])
def test_other_materials_follow_their_own_tables(material, column, energy):
    """Mars (CO2, 165 MeV at 21 g/cm2) and the light gases that need more energy per gram (part 07 F5)."""
    assert pr.energy_mev(column, material) == pytest.approx(energy, rel=0.02)


def test_the_range_and_the_energy_are_inverses():
    for material in pr.MATERIALS:
        for energy in (1.0, 12.0, 250.0, 2500.0, 10000.0, 30000.0):
            assert pr.energy_mev(pr.csda_range_g_cm2(energy, material), material) == pytest.approx(energy, rel=1e-9)


def test_the_table_values_are_returned_at_its_nodes():
    assert pr.csda_range_g_cm2(300.0, "air") == pytest.approx(58.16)
    assert pr.csda_range_g_cm2(10.0, "H") == pytest.approx(0.05346)
    assert pr.csda_range_g_cm2(10000.0, "Fe") == pytest.approx(6389.0)


def test_past_the_last_energy_the_power_law_continues():
    assert pr.energy_mev(10000.0) > 10000.0
    assert pr.csda_range_g_cm2(20000.0, "air") == pytest.approx(10100.0, rel=0.05)
    assert pr.energy_mev(0.0) == 0.0 and pr.csda_range_g_cm2(0.0) == 0.0


def test_a_mixture_of_the_elements_of_air_has_the_range_of_air():
    air = {"N": 0.7552, "O": 0.2318, "Ar": 0.0128, "C": 0.0001}
    for column, _ in AIR_COLUMN_TO_MEV:
        assert pr.energy_mev(column, air) == pytest.approx(pr.energy_mev(column, "air"), rel=0.01)
    assert pr.energy_mev(100.0, {"N": 2.0}) == pytest.approx(pr.energy_mev(100.0, "N"), rel=0.005)
    assert pr.energy_mev(100.0, {"N": 2, "O": 1}) == pytest.approx(pr.energy_mev(100.0, {"N": 0.6667, "O": 0.3333}))


def test_a_light_gas_in_the_mix_needs_more_energy_per_gram():
    nitrogen = pr.energy_mev(22.5, "N")
    half = pr.energy_mev(22.5, {"N": 0.9, "H": 0.1})
    assert nitrogen < half < pr.energy_mev(22.5, "H")


def test_unknown_materials_and_bad_fractions_are_refused():
    with pytest.raises(DataError):
        pr.energy_mev(10.0, "NH3")
    with pytest.raises(DataError):
        pr.energy_mev(10.0, {"N": -1.0, "O": 2.0})
    with pytest.raises(DataError):
        pr.energy_mev(10.0, {"N": 0.0})


def test_the_rigidity_of_a_proton():
    assert pr.rigidity_gv(0.0) == 0.0
    for energy, gv in ((10, 0.137), (30, 0.239), (100, 0.445), (300, 0.808), (1000, 1.696), (3000, 3.825), (10000, 10.90)):
        assert pr.rigidity_gv(energy) == pytest.approx(gv, rel=0.005)


def test_the_file_is_the_scripts_build_with_every_material_complete():
    with open(pr.DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    assert data["state"] == "verified" and data["source"].startswith("NIST")
    assert set(data["materials"]) == set(pr.MATERIALS) >= {"H", "He", "N", "O", "air", "CO2", "CH4", "SiO2", "Fe"}
    for key, row in data["materials"].items():
        assert len(row["csda_range_g_cm2"]) == len(row["stopping_mev_cm2_g"]) == len(data["energies_mev"]), key
