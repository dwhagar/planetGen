# tests/test_mantle_reservoirs.py

"""
GEN.1365: the water and carbon a rocky mantle stores (`physics/mantle_reservoirs.py`, `data/mantle_reservoirs.json`).

The zone shares are held to the Earth (the core-mantle boundary near 136 GPa and a transition zone near a tenth of the mantle), Mars (a mantle that ends
inside the transition zone) and the Moon (no transition zone); the capacities to the ranges the sources print; the carbon form to its one limit.
"""

import json
import types

import pytest

from planetgen.describe import describe, model
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, hydrosphere as hy, mantle_reservoirs as mr
from planetgen.physics.elements import DataError
from planetgen.util import draw

EARTH_KG = constants.EARTH_MASS_TO_KG
EARTH_M = constants.EARTH_RADIUS_KM * 1000.0
OCEAN_KG = hy.EARTH_OCEAN_FRACTION * EARTH_KG


def test_the_earth_has_the_zones_it_has_and_the_core_where_it_is():
    zones = mr.zone_shares(EARTH_KG, EARTH_M, 0.675)
    upper, transition, lower = zones.shares
    assert zones.cmb_pressure_pa / 1e9 == pytest.approx(136.0, rel=0.06)             # PREM: 135.8 GPa
    assert 0.09 < transition < 0.13                                                    # PREM: the 410 to 660 km shell is about a tenth of the mantle
    assert 0.13 < upper < 0.19 and 0.70 < lower < 0.78
    assert sum(zones.shares) == pytest.approx(1.0)


def test_a_small_mantle_ends_inside_the_transition_zone_and_the_moon_has_none():
    mars = mr.zone_shares(0.107 * EARTH_KG, 0.532 * EARTH_M, 0.78)
    assert mars.shares[2] == 0.0 and mars.shares[1] > 0.1 and mars.cmb_pressure_pa / 1e9 == pytest.approx(20.0, rel=0.25)
    moon = mr.zone_shares(0.0123 * EARTH_KG, 0.273 * EARTH_M, 0.7)
    assert moon.shares == (1.0, 0.0, 0.0) and moon.cmb_pressure_pa / 1e9 < 5.0


def test_a_bigger_planet_has_a_thinner_transition_zone_and_a_deeper_core():
    earth = mr.zone_shares(EARTH_KG, EARTH_M, 0.675)
    big = mr.zone_shares(5.0 * EARTH_KG, 1.6 * EARTH_M, 0.67)
    assert big.shares[1] < earth.shares[1] / 2 and big.shares[2] > earth.shares[2] and big.cmb_pressure_pa > 2.5 * earth.cmb_pressure_pa


def test_the_earth_stores_a_few_oceans_and_the_transition_zone_most_of_it():
    capacity = mr.water_capacity_fraction(EARTH_KG, EARTH_M, 0.675) * EARTH_KG
    assert 3.0 * OCEAN_KG < capacity < 8.0 * OCEAN_KG                                  # the data note: three oceans at the primitive value, more than ten at saturation
    zones = mr.zone_shares(EARTH_KG, EARTH_M, 0.675).shares
    from_transition = zones[1] * mr.CAPACITY_FRACTION[1]
    assert from_transition > 0.9 * sum(share * cap for share, cap in zip(zones, mr.CAPACITY_FRACTION))
    # [R1717]: 2 to 8 oceans for 0.5 to 5 Earth masses, in true oceans
    assert 2.0 * OCEAN_KG < mr.water_capacity_fraction(0.5 * EARTH_KG, 0.8 * EARTH_M, 0.675) * 0.5 * EARTH_KG < 8.0 * OCEAN_KG
    assert 2.0 * OCEAN_KG < mr.water_capacity_fraction(5.0 * EARTH_KG, 1.6 * EARTH_M, 0.67) * 5.0 * EARTH_KG < 12.0 * OCEAN_KG


def test_a_body_without_a_transition_zone_stores_little():
    moon = mr.water_capacity_fraction(0.0123 * EARTH_KG, 0.273 * EARTH_M, 0.7) * 0.0123 * EARTH_KG
    assert moon < 0.05 * OCEAN_KG


def test_the_capacities_are_inside_what_the_sources_print():
    upper, transition, lower = (mr.CAPACITY_WT_PERCENT[zone] for zone in mr.ZONES)
    assert 0.005 <= upper <= 0.9 and 1.0 <= transition <= 3.0 and 0.0002 <= lower <= 0.1
    assert transition > upper > lower


def test_carbon_is_carbonate_above_its_limit_and_diamond_or_graphite_below():
    assert mr.carbon_form(mr.CARBON_REDOX_DIW + 0.5) == "carbonate"
    assert mr.carbon_form(mr.CARBON_REDOX_DIW) == "reduced" and mr.carbon_form(-2.0) == "reduced"
    for bad in (float("nan"), None, "3"):
        with pytest.raises(ValueError):
            mr.carbon_form(bad)


def test_refusals():
    for args in ((0.0, EARTH_M, 0.5), (EARTH_KG, 0.0, 0.5), (EARTH_KG, EARTH_M, 0.0), (EARTH_KG, EARTH_M, 1.5)):
        with pytest.raises(ValueError):
            mr.zone_shares(*args)


def test_a_malformed_file_is_refused(tmp_path):
    good = json.load(open(mr.DATA_FILE, encoding="utf-8"))

    def refused(change):
        document = json.loads(json.dumps(good))
        change(document)
        path = tmp_path / "bad.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        with pytest.raises(DataError):
            mr.load(str(path))

    refused(lambda d: d["boundaries_gpa"].update(transition_top=30.0))
    refused(lambda d: d["capacity_wt_percent"].update(transition_zone=0.001))
    refused(lambda d: d["capacity_wt_percent"].pop("lower_mantle"))
    refused(lambda d: d.update(carbon_redox_diw="three"))
    with pytest.raises(DataError):
        mr.load(str(tmp_path / "missing.json"))


def test_the_hydrosphere_keeps_what_the_mantle_cannot():
    """A water world keeps the rest at the surface: 40 Earth oceans of water on a 5 Earth mass body, the mantle keeps its capacity of them."""
    mass, radius, rock = 5.0 * EARTH_KG, 1.6 * EARTH_M, 0.67
    capacity = hy.mantle_max_water_fraction(mass, radius, rock)
    assert capacity == pytest.approx(mr.water_capacity_fraction(mass, radius, rock))
    highest_buffer = types.SimpleNamespace(uniform=lambda low, high: high)                  # the buffer at its top (the mantle keeps an eighth)
    water = 40.0 * hy.EARTH_OCEAN_FRACTION / 5.0
    assert hy.surface_water_fraction(water, mass, radius, rock, rng=highest_buffer) == pytest.approx(water - capacity)
    assert hy.mantle_max_water_fraction(mass, radius, 0.0) == 0.0


def test_the_page_of_a_rocky_planet_shows_its_mantle():
    seen = False
    for seed in range(8):
        draw.set_run_seed(seed)
        config = SystemConfig()
        document = describe(StarSystem(system_config=config))
        model.check(document)
        for section in model.sections(document):
            for block in section["blocks"]:
                if block["type"] == "facts" and block["caption"] == "Mantle":
                    keys = [row["key"] for row in block["rows"]]
                    assert keys[0] == "mantle_water_capacity" and set(keys) <= {"mantle_water_capacity", "mantle_carbon"}
                    seen = True
    assert seen
