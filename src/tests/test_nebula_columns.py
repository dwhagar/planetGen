# tests/test_nebula_columns.py

"""
GEN.1276: the stored composition of a nebula and of a supernova remnant (`physics/nebula_columns.py`).

* every class gets a full set of columns whose fractions add up, the same numbers whenever its key is the same, other numbers for another key;
* a planetary nebula and a remnant store their progenitor's mass, the other nebulae none; a remnant stores its ejecta and swept-up parts;
* a row inserted through `store.insert_nebula` / `insert_supernova_remnant` carries its composition and the plain-words sentences written from it;
* migration 0205 fills the rows of a v204 galaxy with the same numbers, and a second run changes nothing.
"""

import math

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.physics import nebula_columns as nc
from tests.db_schema_support import load_old_schema

REMNANT_CLASSES = [letter for letter, data in tuning.NEBULA_CLASSES.items() if data["family"] == "supernova-remnant"]
NEBULA_CLASSES = [letter for letter, data in tuning.NEBULA_CLASSES.items() if data["family"] != "supernova-remnant"]


def nebula_row(letter, uid=b"\x02" * 12, density=100.0, radius=5.0, galactic_radius=8000.0):
    return {"nebula_class": letter, "density_cm3": density, "radius_ly": radius, "galactic_radius_pc": galactic_radius, "uid": uid, "name": "n", "hosts": "WN6 x1"}


def remnant_row(letter, uid=b"\x01" * 12, compact="neutron_star"):
    return {"remnant_class": letter, "progenitor_type": "Type Ia" if letter == "W" else "core-collapse", "compact_remnant_kind": None if letter == "W" else compact,
            "density_cm3": 1.0, "radius_ly": 20.0, "galactic_radius_pc": 8000.0, "uid": uid, "name": "r"}


@pytest.mark.parametrize("letter", NEBULA_CLASSES)
def test_a_nebula_of_every_class_has_a_full_composition_that_adds_up(letter):
    values = nc.nebula_values(nebula_row(letter))
    assert set(values) == set(nc.COMMON_COLUMNS)
    assert values["mass_fraction_h"] + values["mass_fraction_he"] + values["mass_fraction_z"] == pytest.approx(1.0)
    assert values["ionized_share"] + values["atomic_share"] + values["molecular_share"] == pytest.approx(1.0)
    assert sum(values[column] for column in nc.GRAIN_COLUMNS.values()) == pytest.approx(1.0)
    assert all(0.0 <= values[f"mass_fraction_{e.lower()}"] <= values["mass_fraction_z"] + 1e-12 for e in nc.ABUNDANCE_ELEMENTS)
    assert (values["progenitor_mass_sol"] is not None) == (letter in nc.PLANETARY)


@pytest.mark.parametrize("letter", REMNANT_CLASSES)
def test_a_remnant_of_every_class_stores_its_progenitor_and_its_two_parts(letter):
    values = nc.remnant_values(remnant_row(letter))
    assert set(values) == set(nc.COLUMNS["supernova_remnants"])
    assert values["mass_fraction_h"] + values["mass_fraction_he"] + values["mass_fraction_z"] == pytest.approx(1.0)
    assert values["progenitor_mass_sol"] > 0 and values["ejecta_mass_sol"] > 0 and values["swept_mass_sol"] > 0
    for part in ("ejecta", "swept"):
        assert values[f"{part}_mass_fraction_h"] + values[f"{part}_mass_fraction_he"] + values[f"{part}_mass_fraction_z"] == pytest.approx(1.0)
    assert values["ejecta_mass_fraction_h"] < values["swept_mass_fraction_h"]       # the ejecta are metal-rich, the swept-up gas is interstellar


def test_the_same_key_gives_the_same_numbers_and_another_key_other_ones():
    first = nc.nebula_values(nebula_row("K"))
    assert nc.nebula_values(nebula_row("K")) == first
    other = nc.nebula_values(nebula_row("K", uid=b"\x03" * 12))
    assert other["progenitor_mass_sol"] != first["progenitor_mass_sol"] and other["gas_oh_dex"] != first["gas_oh_dex"]
    assert nc.remnant_values(remnant_row("R")) == nc.remnant_values(remnant_row("R"))


def test_the_gas_is_poorer_further_out_and_a_carbon_rich_nebula_makes_carbon():
    near = [nc.nebula_values(nebula_row("D", uid=bytes([n]) * 12, galactic_radius=4000.0))["gas_oh_dex"] for n in range(1, 30)]
    far = [nc.nebula_values(nebula_row("D", uid=bytes([n]) * 12, galactic_radius=14000.0))["gas_oh_dex"] for n in range(1, 30)]
    assert sum(near) / len(near) > sum(far) / len(far)
    carbon = nc.nebula_values(nebula_row("K"))
    ordinary = nc.nebula_values(nebula_row("D"))
    assert carbon["mass_fraction_c"] > 2.0 * ordinary["mass_fraction_c"]


def test_a_row_with_no_uid_is_keyed_on_its_name():
    named = {**nebula_row("D"), "uid": None}
    assert nc.row_key("nebulae", named) == "nebulae:n"
    assert nc.nebula_values(named) == nc.nebula_values({**named})
    assert nc.nebula_values(named)["gas_oh_dex"] != nc.nebula_values({**named, "name": "other"})["gas_oh_dex"]


def test_the_sentences_are_written_from_the_numbers():
    values = nc.nebula_values(nebula_row("D"))
    sentence, species = nc.sentences_for(values)
    assert f"{100.0 * values['mass_fraction_he']:.1f}% helium" in sentence and "ionized" in sentence
    assert species.startswith("ionized hydrogen (H+)") and len(species) <= 255
    dark = nc.sentences_for(nc.nebula_values(nebula_row("N", density=1e4, radius=2.0)))[1]
    assert "molecular hydrogen" in dark or "atomic hydrogen" in dark and "ice mantles" in dark


def test_no_nan_is_stored():
    for letter in NEBULA_CLASSES:
        assert not any(isinstance(v, float) and math.isnan(v) for v in nc.nebula_values(nebula_row(letter)).values())
    for letter in REMNANT_CLASSES:
        assert not any(isinstance(v, float) and math.isnan(v) for v in nc.remnant_values(remnant_row(letter)).values())


# --- the database ------------------------------------------------------------------------------------------------------------------------------------

def test_an_inserted_nebula_and_remnant_carry_their_composition(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        nebula = Nebula(SystemConfig(), nebula_class="D")
        nebula_id = store.insert_nebula(conn, nebula)
        remnant = SupernovaRemnant(SystemConfig())
        remnant_id = store.insert_supernova_remnant(conn, remnant)
        conn.commit()
        row = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        assert row["mass_fraction_h"] is not None and row["dust_to_gas"] > 0 and row["progenitor_mass_sol"] is None
        assert f"{100.0 * row['mass_fraction_he']:.1f}% helium" in row["composition"] and row["dominant_species"].startswith(("ionized", "atomic", "molecular"))
        again = nc.values_for("nebulae", dict(conn.execute(f"SELECT {store._COMPOSITION_READ['nebulae']} FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()))
        assert {column: row[column] for column in nc.COMMON_COLUMNS} == pytest.approx(again, nan_ok=True)
        remnant_row_ = conn.execute("SELECT * FROM supernova_remnants WHERE id = ?", (remnant_id,)).fetchone()
        assert remnant_row_["progenitor_mass_sol"] > 0 and remnant_row_["ejecta_mass_sol"] > 0 and remnant_row_["swept_mass_sol"] > 0
    finally:
        conn.close()


@pytest.mark.fresh_database
def test_migration_0205_fills_the_rows_of_a_v204_galaxy_with_the_numbers_of_the_model(mysql_config):
    load_old_schema(mysql_config, 204)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for column in ("mass_fraction_h", "progenitor_mass_sol"):
            assert not conn.execute("SHOW COLUMNS FROM nebulae LIKE '" + column + "'").fetchall()
        conn.execute(
            "INSERT INTO nebulae (name, nebula_class, nebula_type, radius_ly, composition, formation_cause, dominant_species, density_cm3, temperature_k, extinction_av,"
            " galactic_orbital_period_gy, galactic_orbital_phase_deg, galactic_min_update_interval_years, galactic_radius_pc, center_x_pc, center_y_pc, center_z_pc, uid)"
            " VALUES ('Old', 'K', 'planetary', 1.0, 'old text', 'cause', 'old species', 100.0, 10000.0, 0.5, 0.2, 10.0, 1000.0, 8000.0, 8000.0, 0.0, 0.0, ?)", (b"\x05" * 12,))
        conn.execute(
            "INSERT INTO supernova_remnants (name, remnant_class, morphology, age_years, radius_ly, progenitor_type, dominant_species, density_cm3, temperature_k, extinction_av,"
            " galactic_orbital_period_gy, galactic_orbital_phase_deg, galactic_min_update_interval_years, galactic_radius_pc, center_x_pc, center_y_pc, center_z_pc, uid)"
            " VALUES ('Old remnant', 'W', 'shell', 500.0, 10.0, 'Type Ia', 'old species', 1.0, 1e6, 0.1, 0.2, 10.0, 1000.0, 8000.0, 8000.0, 0.0, 0.0, ?)", (b"\x06" * 12,))
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT * FROM nebulae WHERE name = 'Old'").fetchone()
        expected = nc.nebula_values({"nebula_class": "K", "density_cm3": 100.0, "radius_ly": 1.0, "galactic_radius_pc": 8000.0, "uid": b"\x05" * 12, "name": "Old", "hosts": None})
        assert {column: row[column] for column in nc.COMMON_COLUMNS} == pytest.approx(expected, nan_ok=True)
        assert row["composition"] != "old text" and row["dominant_species"] != "old species"
        remnant = conn.execute("SELECT * FROM supernova_remnants WHERE name = 'Old remnant'").fetchone()
        assert remnant["progenitor_mass_sol"] > 0 and remnant["ejecta_mass_sol"] > 0 and remnant["dominant_species"] != "old species"
    finally:
        conn.close()
