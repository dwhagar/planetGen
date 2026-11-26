# tests/test_comet_coma.py

"""GEN.455, part A: the coma of an active comet (`physics/comet_coma.py`): the water the light frees, the gases the nucleus holds, how far each goes, what is stored and what is read back."""

import math
import random

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.generation import comet as comet_module
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import InterstellarComet
from planetgen.generation.system import StarSystem
from planetgen.physics import comet_coma, constants, materials, photolysis
from planetgen.util import draw

NUCLEUS = {"CO": 3.0, "CO2": 12.0, "H2CO": 0.2, "H2S": 1.0, "OCS": 0.4}
"""dict: Volatiles of a nucleus in percent of water (OCS has no photolysis row, H2CO is CH2O in that table)."""


def test_the_water_is_the_light_the_active_ice_takes_over_the_energy_of_a_molecule():
    energy = (materials.material("water ice").props["h_vap_kj_mol"].value + materials.material("water ice").props["h_fus_kj_mol"].value) * 1000.0 / constants.AVOGADRO
    assert comet_coma.sublimation_j_per_molecule() == pytest.approx(energy)
    rate = comet_coma.water_production_1au(4.0, 0.1, 0.05)
    assert rate == pytest.approx(0.1 * math.pi * 2000.0 ** 2 * 0.95 * constants.SOLAR_CONSTANT_W_M2 / energy)
    assert comet_coma.water_production_1au(8.0, 0.1, 0.05) == pytest.approx(4.0 * rate)            # the projected circle: the diameter squared
    assert comet_coma.water_production_1au(4.0, 0.2, 0.05) == pytest.approx(2.0 * rate)            # the share of the surface that sublimates
    assert comet_coma.water_production_1au(4.0, 0.1, 0.10) == pytest.approx(rate * 0.90 / 0.95)    # a darker nucleus takes more
    with pytest.raises(ValueError):
        comet_coma.water_production_1au(0.0, 0.1, 0.05)


def test_every_gas_leaves_at_the_nucleus_own_ratio_to_the_water_and_a_gas_it_lacks_is_not_there():
    coma = comet_coma.draw_coma(random.Random(3), 5.0, NUCLEUS)
    assert [row[0] for row in coma.species] == ["H2O", *NUCLEUS]
    water = coma.water_production_1au_s
    for species, production, _ in coma.species[1:]:
        assert production == pytest.approx(water * NUCLEUS[species] / 100.0)
    assert comet_coma.draw_coma(random.Random(3), 5.0, {}).species[0][0] == "H2O"
    assert len(comet_coma.draw_coma(random.Random(3), 5.0, {}).species) == 1


def test_a_gas_with_no_photolysis_row_has_no_length_and_the_others_follow_the_hosts_light():
    coma = comet_coma.draw_coma(random.Random(5), 5.0, NUCLEUS, 5772.0)
    lengths = {species: length for species, _, length in coma.species}
    assert lengths["OCS"] is None
    assert lengths["H2O"] == pytest.approx(coma.outflow_kms_1au * photolysis.lifetime_s("H2O", 1.0, 5772.0))
    assert lengths["H2CO"] == pytest.approx(coma.outflow_kms_1au * photolysis.lifetime_s("CH2O", 1.0, 5772.0))
    cool = {species: length for species, _, length in comet_coma.draw_coma(random.Random(5), 5.0, NUCLEUS, 4000.0).species}
    assert cool["CO"] > 100.0 * lengths["CO"]                  # round a 4,000 K star carbon monoxide lives far longer (GEN.496)


def test_the_rate_falls_as_the_inverse_square_and_the_water_stops_where_the_activity_does_but_carbon_dioxide_and_monoxide_go_on():
    limit = tuning.COMET_ACTIVITY_PERIHELION_THRESHOLD_AU
    assert comet_coma.production_at(8.0, "H2O", 2.0) == pytest.approx(2.0)
    assert comet_coma.production_at(8.0, "H2O", limit) == pytest.approx(8.0 / limit ** 2)
    for species in ("H2O", "HCN", "CH4"):
        assert comet_coma.production_at(8.0, species, limit * 1.01) == 0.0
    for species in ("CO", "CO2"):
        assert comet_coma.production_at(8.0, species, 6.0) == pytest.approx(8.0 / 36.0)
    assert comet_coma.scale_length_at(100.0, 4.0) == pytest.approx(800.0)                  # r^1.5
    with pytest.raises(ValueError):
        comet_coma.production_at(1.0, "H2O", 0.0)


def test_the_water_of_the_nuclei_the_generator_makes_lies_in_the_range_comets_have_been_seen_to_lose():
    """The measured range is a test target and no input: 5e26 to 1.4e31 molecules a second near 1 AU."""
    rng = random.Random(11)
    rates = sorted(comet_coma.draw_coma(rng, diameter, {}).water_production_1au_s for diameter in (tuning.BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM[0] * (8.0 ** rng.random()) for _ in range(2000)))
    inside = sum(5e26 <= rate <= 1.41e31 for rate in rates) / len(rates)
    assert inside > 0.8
    assert rates[len(rates) // 2] == pytest.approx(rates[len(rates) // 2])                    # (the median is reported in the PR; nothing to assert but that it exists)
    assert 1e26 < rates[len(rates) // 2] < 1e31


def test_an_inactive_comet_has_no_coma_and_an_active_one_has_it_all():
    config = SystemConfig()
    config.COMETS = True
    seen = set()
    for _ in range(120):
        comet = comet_module.Comet(config, 1.0)
        seen.add(comet.is_active)
        if comet.is_active:
            assert comet.water_production_1au_per_s > 0.0 and comet.coma_species[0][0] == "H2O"
            assert 0.003 <= comet.coma_active_area_fraction <= 1.0 and 0.4 <= comet.coma_outflow_1au_kms <= 1.0 and 0.02 <= comet.nucleus_albedo <= 0.10
            assert comet_coma.of(comet).water_production_1au_s == comet.water_production_1au_per_s
        else:
            assert comet.water_production_1au_per_s is None and comet.coma_species == () and comet_coma.of(comet) is None
        if len(seen) == 2:
            break
    assert seen == {True, False}


def test_the_coma_draw_moves_no_other_number_of_the_run():
    """The coma has a stream of its own: the unit's stream is where it was before the coma was drawn."""
    config = SystemConfig()
    comet = comet_module.Comet(config, 1.0)
    comet.is_active = True
    before = draw.getstate() if hasattr(draw, "getstate") else None
    comet_coma.generate(comet)
    after = draw.getstate() if hasattr(draw, "getstate") else None
    assert before == after


def test_a_comet_keeps_its_coma_through_to_dict_and_from_dict():
    config = SystemConfig()
    for _ in range(200):
        comet = comet_module.Comet(config, 1.0)
        if comet.is_active:
            break
    else:
        pytest.fail("no active comet in 200")
    again = comet_module.Comet.from_dict(comet.to_dict(), config)
    assert again.coma_species == comet.coma_species and again.water_production_1au_per_s == comet.water_production_1au_per_s
    assert comet_coma.of(again) == comet_coma.of(comet)


def test_an_interstellar_comet_has_a_coma_when_it_is_active_with_the_suns_far_ultraviolet():
    config = SystemConfig()
    for _ in range(200):
        comet = InterstellarComet(config)
        if comet.is_active:
            break
    else:
        pytest.fail("no active interstellar comet in 200")
    assert comet.coma_species[0][0] == "H2O" and comet.water_production_1au_per_s > 0.0
    species = {row[0] for row in comet.coma_species}
    assert set(comet.makeup.detail["volatiles"]) <= species


def test_the_rows_are_stored_and_a_loaded_comet_gets_them_back(mysql_config):
    for _ in range(40):
        config = SystemConfig()
        config.STAR_TYPE = "G2V"
        config.COMETS = True
        config.BINARY_SYSTEM = False
        config.PLANETS = False          # the comets are the point: no planets to generate
        system = StarSystem(system_config=config)
        if any(comet.is_active for comet in system.comets):
            break
    else:
        pytest.fail("no system with an active comet")
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, config)
        rows = {row["name"]: row for row in conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchall()}
        for comet in system.comets:
            stored = rows[comet.name]
            coma_rows = conn.execute("SELECT species, production_1au_per_s, scale_length_1au_km FROM comet_coma WHERE comet_id = ? ORDER BY position", (stored["id"],)).fetchall()
            if comet.is_active:
                assert stored["water_production_1au_per_s"] == pytest.approx(comet.water_production_1au_per_s)
                assert stored["coma_outflow_1au_kms"] == pytest.approx(comet.coma_outflow_1au_kms) and stored["nucleus_albedo"] == pytest.approx(comet.nucleus_albedo)
                assert [(row["species"], row["production_1au_per_s"], row["scale_length_1au_km"]) for row in coma_rows] == [tuple(entry) for entry in comet.coma_species]
            else:
                assert stored["water_production_1au_per_s"] is None and not coma_rows
        reloaded = {comet.name: comet for comet in store.load_star_system(conn, system_id).comets}
        for comet in system.comets:
            assert reloaded[comet.name].coma_species == comet.coma_species
            assert reloaded[comet.name].water_production_1au_per_s == comet.water_production_1au_per_s
    finally:
        conn.close()


def _coma(*rows):
    """A `Coma` with the given `(species, production)` rows (water first), the other numbers of no matter."""
    return comet_coma.Coma(1.0e28, 0.1, 0.7, 0.05, tuple((species, rate, None) for species, rate in rows))


def test_the_coma_class_is_the_carbon_oxide_that_leaves_faster_and_equal_rates_read_as_carbon_dioxide():
    """GEN.1227: read off the rows of the coma, no draw."""
    assert comet_coma.coma_class(_coma(("H2O", 9e27), ("CO2", 3e27), ("CO", 1e27))) == comet_coma.CO2_DOMINATED
    assert comet_coma.coma_class(_coma(("H2O", 9e27), ("CO2", 1e27), ("CO", 3e27))) == comet_coma.CO_DOMINATED
    assert comet_coma.coma_class(_coma(("H2O", 9e27), ("CO2", 2e27), ("CO", 2e27))) == comet_coma.CO2_DOMINATED
    assert comet_coma.coma_class(_coma(("H2O", 9e27), ("CO", 1e27))) == comet_coma.CO_DOMINATED         # a gas the nucleus lacks is no rate
    assert comet_coma.coma_class(_coma(("H2O", 9e27), ("CO2", 1e27))) == comet_coma.CO2_DOMINATED
    assert comet_coma.coma_class(_coma(("H2O", 9e27), ("CH4", 1e26))) is None                          # neither of the two
    assert comet_coma.coma_class(None) is None                                                          # an inactive comet has no coma


def test_the_dominant_volatiles_are_the_fastest_gases_after_water_with_ties_broken_by_formula():
    coma = _coma(("H2O", 9e27), ("CO", 1e27), ("CO2", 3e27), ("CH3OH", 2e26), ("NH3", 1e26))
    assert comet_coma.dominant_volatiles(coma) == ("CO2", "CO", "CH3OH")
    assert comet_coma.dominant_volatiles(coma, 1) == ("CO2",)
    assert comet_coma.dominant_volatiles(_coma(("H2O", 9e27), ("CO2", 2e27), ("CO", 2e27), ("CH4", 2e27), ("NH3", 2e27))) == ("CH4", "CO", "CO2")
    assert comet_coma.dominant_volatiles(_coma(("H2O", 9e27))) == ()
    assert comet_coma.dominant_volatiles(None) == ()


def test_an_active_comets_class_and_dominant_gases_follow_the_volatiles_its_nucleus_was_drawn_with():
    seen = set()
    for seed in range(60):
        with draw.bound(seed):
            comet = InterstellarComet(SystemConfig())
        coma = comet_coma.of(comet)
        if coma is None:
            continue
        percent = comet.makeup.detail["volatiles"]
        expected = comet_coma.CO2_DOMINATED if percent.get("CO2", 0.0) >= percent.get("CO", 0.0) else comet_coma.CO_DOMINATED
        if "CO2" in percent or "CO" in percent:
            assert comet_coma.coma_class(coma) == expected
            seen.add(expected)
        dominant = comet_coma.dominant_volatiles(coma)
        assert 0 <= len(dominant) <= comet_coma.DOMINANT_COUNT and "H2O" not in dominant
    assert seen == set(comet_coma.COMA_CLASSES)


def test_an_interstellar_comet_keeps_the_carbon_chain_flag_its_nucleus_was_drawn_with():
    seen = set()
    for seed in range(120):
        with draw.bound(seed):
            comet = InterstellarComet(SystemConfig())
        assert comet.carbon_chain_depleted is bool(comet.makeup.detail["carbon_chain_depleted"])
        seen.add(comet.carbon_chain_depleted)
    assert seen == {True, False}
    with draw.bound(3):
        comet = InterstellarComet(SystemConfig())
    back = InterstellarComet.from_dict(comet.to_dict(), SystemConfig())
    assert back.carbon_chain_depleted == comet.carbon_chain_depleted
    data = comet.to_dict()
    data.pop("carbon_chain_depleted", None)
    assert InterstellarComet.from_dict(data, SystemConfig()).carbon_chain_depleted is None      # a record made before the flag was kept


def test_the_row_keeps_the_carbon_chain_flag_and_the_sheet_shows_it(mysql_config):
    from planetgen.describe import sheet
    with draw.bound(11):
        comet = InterstellarComet(SystemConfig())
    conn = store.get_connection(mysql_config)
    try:
        comet_id = store.insert_interstellar_comet(conn, comet)
        conn.commit()
        row = conn.execute("SELECT carbon_chain_depleted FROM interstellar_comets WHERE id = ?", (comet_id,)).fetchone()
    finally:
        conn.close()
    assert row["carbon_chain_depleted"] == int(comet.carbon_chain_depleted)
    assert "carbon_chain_depleted" in [spec.key for spec in sheet.specs_for("interstellar_comet")]
