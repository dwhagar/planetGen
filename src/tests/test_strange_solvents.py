"""GEN.894: the ammonia of a body's water, its silicon-rich crust and the sulfur supply of the acid eaters."""

import json

import pytest

from planetgen.physics import atmosphere, life_supply, strange_solvents


def test_the_ammonia_is_a_tenth_of_the_other_ice_over_the_water():
    assert strange_solvents.ammonia_water_mass_fraction({"water_mass_fraction": 0.2, "ice_mass_fraction": 0.4}) == pytest.approx(0.2)
    assert strange_solvents.ammonia_water_mass_fraction({"water_mass_fraction": 0.2, "ice_mass_fraction": 0.0}) == 0.0


def test_a_body_without_water_has_no_ammonia_water_and_one_without_the_shares_has_none_to_read():
    assert strange_solvents.ammonia_water_mass_fraction({"water_mass_fraction": 0.0, "ice_mass_fraction": 0.4}) == 0.0
    assert strange_solvents.ammonia_water_mass_fraction({"water_mass_fraction": 0.2}) is None
    assert strange_solvents.ammonia_water_mass_fraction({}) is None


def test_the_crust_is_silicon_rich_at_a_molar_mg_si_of_one_or_under():
    assert strange_solvents.silicon_rich_crust({"rock_mg_si_ratio": 1.0}) == 1.0
    assert strange_solvents.silicon_rich_crust({"rock_mg_si_ratio": 1.07}) == 0.0
    assert strange_solvents.silicon_rich_crust({}) is None


def test_every_number_has_a_tier_a_tweak_and_a_basis_and_a_bad_one_is_refused(tmp_path, monkeypatch):
    for entry in strange_solvents.NUMBERS.values():
        assert entry["tier"] and entry["tweak"] and entry["basis"]
    with open(strange_solvents.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    del document["numbers"]["nh3_share_of_ice"]["basis"]
    path = tmp_path / "strange_solvents.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    monkeypatch.setattr(strange_solvents, "DATA_FILE", str(path))
    with pytest.raises(strange_solvents.StrangeSolventError):
        strange_solvents._load()


def _air(so2_share):
    total = 9.2e6
    return {"atmospheric_pressure_pa": total, "surface_temperature_k": 737.0,
            atmosphere.GAS_FIELD: {"CO2": total * (1.0 - so2_share), "SO2": total * so2_share}}


def test_the_sulfur_supply_is_the_sulfur_dioxide_share_over_venus_s_and_needs_no_water():
    reference = life_supply.number("sulfur_so2_ref_share")
    assert life_supply.supplies(_air(reference / 4.0))["sulfur_supply"] == pytest.approx(0.25)
    assert life_supply.supplies(_air(reference * 10.0))["sulfur_supply"] == 1.0
    assert life_supply.supplies(_air(0.0))["sulfur_supply"] == 0.0
    assert "sulfur_supply" in life_supply.QUANTITIES
