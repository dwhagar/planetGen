# tests/test_ice_species.py

"""
GEN.834: the non-water ice of a body is stored by species. The eight amounts add up to `ice_mass_fraction`, an icy body's come from the comet
volatiles it drew, a body with no ice holds zeros, a body saved before holds none and reads back as the medians, and a Pluto-like and a
CO2-frost Mars-like body keep what they were given through the database.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import bulk_composition, ice_species, surface_composition
from planetgen.util import draw


def _bulk(ice, shares=None):
    return bulk_composition.BulkFractions(metal=0.2, rock=0.5 - ice, water=0.3, ice=ice, envelope=0.0, ice_shares=shares)


def test_the_amounts_add_up_to_the_ice_and_a_given_split_is_kept():
    shares = {"N2": 0.9, "CO": 0.06, "CH4": 0.04}
    values = _bulk(0.1, shares).to_attributes()
    assert sum(values[field] for field in ice_species.BODY_FIELDS) == pytest.approx(0.1)
    assert values["ice_n2_mass_fraction"] == pytest.approx(0.09)
    assert values["ice_co2_mass_fraction"] == 0.0


def test_an_ice_with_no_given_split_takes_the_medians_and_no_ice_holds_zeros():
    values = _bulk(0.1).to_attributes()
    assert values["ice_co2_mass_fraction"] == pytest.approx(0.1 * ice_species.median_shares()["CO2"])
    assert sum(values[field] for field in ice_species.BODY_FIELDS) == pytest.approx(0.1)
    assert all(_bulk(0.0).to_attributes()[field] == 0.0 for field in ice_species.BODY_FIELDS)


def test_the_stored_amounts_read_back_as_shares_of_the_ice():
    class Body:
        pass
    body = Body()
    for field, value in _bulk(0.2, {"CO2": 0.75, "other": 0.25}).to_attributes().items():
        setattr(body, field, value)
    assert ice_species.shares_of(body) == pytest.approx({"CO2": 0.75, "other": 0.25})
    assert ice_species.shares_of(Body()) is None


def test_the_ground_draw_reads_the_stored_split_and_a_body_saved_before_reads_the_medians():
    class Body:
        ice_mass_fraction = 0.2
    body = Body()
    assert surface_composition._volatile_shares(body).keys() == {"co2_ice", "ch4_ice", "nh3_ice", "co_ice", "hcn_ice"}
    for field, value in _bulk(0.2, {"N2": 0.98, "CH4": 0.02}).to_attributes().items():
        setattr(body, field, value)
    assert surface_composition._volatile_shares(body) == pytest.approx({"n2_ice": 0.98, "ch4_ice": 0.02})



def _system():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    for seed in range(60):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(hasattr(body, "moons") and body.body_type != "g" for body in system.planets):
            return system, cfg
    pytest.fail("no system with a rocky planet in 60 seeds")


@pytest.mark.fresh_database
def test_a_pluto_like_and_a_co2_frost_mars_like_body_keep_their_species_through_the_database(mysql_config):
    system, cfg = _system()
    rocky = [planet for planet in system.planets if hasattr(planet, "moons") and planet.body_type != "g"]
    pluto, mars = rocky[0], rocky[-1] if len(rocky) > 1 else rocky[0]
    split = {pluto: _bulk(0.15, {"N2": 0.97, "CH4": 0.02, "CO": 0.01}), mars: _bulk(0.01, {"CO2": 1.0})}
    for planet, bulk in split.items():
        for field, value in bulk.to_attributes().items():
            setattr(planet, field, value)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        index = {planet.name: planet for planet in loaded.planets if hasattr(planet, "moons")}
        for before in rocky:
            for field in ice_species.BODY_FIELDS:
                assert getattr(index[before.name], field) == pytest.approx(getattr(before, field)), (before.name, field)
        assert index[pluto.name].ice_n2_mass_fraction == pytest.approx(0.15 * 0.97)
        assert index[mars.name].ice_co2_mass_fraction == pytest.approx(0.01) or mars is pluto
    finally:
        conn.rollback()
        conn.close()
