# tests/test_phosphorus_budget.py

"""GEN.1367: the phosphorus chain from the mantle's P2O5 to the ocean's phosphate (`physics/phosphorus_budget.py`)."""

import math
from types import SimpleNamespace

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import interior_heat, phosphorus_budget as pb
from planetgen.util import draw
from tests import class_support

EARTH = 5.972e24
EARTH_SHARES = (1.0, 0.0, 0.0)


def _body(**over):
    """An Earth-like stored state: a mobile lid, a surface ocean over 71 percent of the surface, the Earth's mantle P2O5 and an oxidized mantle."""
    base = dict(body_type="t", rock_p2o5_fraction=2.1e-4, hydrosphere="surface ocean", mantle_delta_iw=3.5, ocean_class="neutral", land_fraction=0.29, ocean_fraction=0.71,
                ocean_depth_km=3.7, radius=6371.0)
    base.update(over)
    return SimpleNamespace(**base)


@pytest.fixture
def earth_heat(monkeypatch):
    monkeypatch.setattr(pb, "_flux_of", lambda body: 0.08)
    monkeypatch.setattr(pb.interior_heat, "lid_of", lambda body: getattr(body, "lid", interior_heat.MOBILE))


# --- the links -----------------------------------------------------------------------------------------------------------
def test_the_earths_mantle_makes_its_crust():
    assert pb.crust_p2o5_wt(2.1e-4, 0.15) == pytest.approx(pb.EARTH_CRUST_P2O5_WT, rel=0.05)
    assert pb.crust_p2o5_wt(2.1e-4, 0.03) > pb.crust_p2o5_wt(2.1e-4, 0.30)          # a smaller melt concentrates the incompatible element more


def test_a_reduced_mantle_puts_the_crusts_phosphorus_in_schreibersite_and_a_dry_one_in_merrillite():
    for delta_iw, wet in ((3.5, True), (3.5, False), (-3.0, True), (-6.0, False), (None, True)):
        assert math.fsum(pb.host_shares(delta_iw, wet)) == pytest.approx(1.0)
    assert pb.host_shares(3.5, True) == pytest.approx(EARTH_SHARES, abs=1e-3)
    assert pb.host_shares(None, True) == EARTH_SHARES
    assert pb.host_shares(pb.PHOSPHIDE_MIDPOINT_DIW, True)[2] == pytest.approx(0.5)
    assert pb.host_shares(-6.0, True)[2] > 0.99 > 0.01 > pb.host_shares(0.0, True)[2]
    dry = pb.host_shares(3.5, False)
    assert dry[0] == pytest.approx(pb.APATITE_DRY_SHARE, abs=1e-3) and dry[1] > dry[0]


def test_what_the_water_can_take_follows_the_hosts():
    assert pb.reactive_share(EARTH_SHARES) == 1.0
    assert pb.reactive_share((0.0, 1.0, 0.0)) == pb.MERRILLITE_FACTOR > 1.0
    assert pb.reactive_share((0.0, 0.0, 1.0)) == pb.SCHREIBERSITE_FACTOR < 0.1


def test_the_earth_reproduces_its_own_ocean():
    value = pb.ocean_phosphate(pb.EARTH_CRUST_P2O5_WT, EARTH_SHARES, 1.0, 1.0, 1.0, "neutral")
    assert value == pytest.approx(pb.EARTH_OCEAN_UMOL_L)


def test_each_ratio_scales_the_phosphate_as_its_link_does():
    base = pb.EARTH_OCEAN_UMOL_L
    half = pb.EARTH_CRUST_P2O5_WT / 2
    assert pb.ocean_phosphate(half, EARTH_SHARES, 1.0, 1.0, 1.0, "neutral") == pytest.approx(base / 2)
    assert pb.ocean_phosphate(pb.EARTH_CRUST_P2O5_WT, EARTH_SHARES, pb.STAGNANT_EROSION, 1.0, 1.0, "neutral") == pytest.approx(base * pb.STAGNANT_EROSION)
    assert pb.ocean_phosphate(pb.EARTH_CRUST_P2O5_WT, EARTH_SHARES, 1.0, 1.0, 10.0, "neutral") == pytest.approx(base / 10)            # ten oceans dilute it
    assert pb.ocean_phosphate(pb.EARTH_CRUST_P2O5_WT, (0.0, 0.0, 1.0), 1.0, 1.0, 1.0, "neutral") == pytest.approx(base * pb.SCHREIBERSITE_FACTOR)


def test_a_soda_lake_climbs_to_its_ceiling_and_a_sealed_ocean_holds_the_floor():
    soda = pb.ocean_phosphate(pb.EARTH_CRUST_P2O5_WT, EARTH_SHARES, 1.0, 1.0, 1.0, "soda")
    assert soda == pytest.approx(pb.EARTH_OCEAN_UMOL_L * pb.SODA_RESIDENCE) and soda < pb.CEILING_UMOL_L      # an Earth-like crust under soda water: 920
    assert pb.ocean_phosphate(2 * pb.EARTH_CRUST_P2O5_WT, EARTH_SHARES, 1.0, 1.0, 1.0, "soda") == pb.CEILING_UMOL_L
    assert pb.ocean_phosphate(pb.EARTH_CRUST_P2O5_WT, EARTH_SHARES, 1.0, 1.0, 1.0, "ice-sealed") == pb.FLOOR_UMOL_L
    assert pb.ocean_phosphate(0.0, EARTH_SHARES, 1.0, 1.0, 1.0, "neutral") == pb.FLOOR_UMOL_L


def test_the_bands_of_the_number():
    assert [pb.phosphorus_class(value) for value in (0.01, 0.099, 0.1, 2.3, 99.0, 100.0, 1000.0)] == ["starved", "starved", "limited", "limited", "limited", "high", "high"]


# --- a body ---------------------------------------------------------------------------------------------------------------
def test_an_earth_like_body_gets_the_earths_crust_hosts_and_ocean(earth_heat):
    body = _body()
    pb.generate(body)
    assert body.crust_p2o5_wt == pytest.approx(pb.crust_p2o5_wt(2.1e-4, pb.surface_composition.melt_fraction(0.08)))
    assert body.crust_apatite_share == pytest.approx(1.0, abs=1e-3)
    assert math.fsum((body.crust_apatite_share, body.crust_merrillite_share, body.crust_schreibersite_share)) == pytest.approx(1.0)
    assert 0.5 < body.ocean_phosphate_umol_l / pb.EARTH_OCEAN_UMOL_L < 2.0 and body.phosphorus == "limited"


def test_a_stagnant_lid_and_a_deep_ocean_both_lower_the_phosphate(earth_heat):
    earth, stagnant, deep = _body(), _body(lid=interior_heat.STAGNANT), _body(ocean_depth_km=37.0)
    for body in (earth, stagnant, deep):
        pb.generate(body)
    assert stagnant.ocean_phosphate_umol_l == pytest.approx(earth.ocean_phosphate_umol_l * pb.STAGNANT_EROSION)
    assert deep.ocean_phosphate_umol_l == pytest.approx(earth.ocean_phosphate_umol_l / 10.0)


def test_an_ocean_under_ice_is_reached_by_its_sea_floor_only(earth_heat):
    open_sea, covered = _body(), _body(hydrosphere="ice-covered ocean", land_fraction=0.0, ocean_fraction=1.0)
    pb.generate(open_sea)
    pb.generate(covered)
    assert covered.ocean_phosphate_umol_l < open_sea.ocean_phosphate_umol_l


def test_a_reduced_mantle_starves_the_ocean_of_weathering(earth_heat):
    oxidized, reduced = _body(), _body(mantle_delta_iw=-5.0)
    pb.generate(oxidized)
    pb.generate(reduced)
    assert reduced.crust_schreibersite_share > 0.95
    assert reduced.ocean_phosphate_umol_l < oxidized.ocean_phosphate_umol_l / 10.0


def test_a_dry_body_has_a_crust_and_no_ocean(earth_heat):
    body = _body(hydrosphere="dry", ocean_class=None, ocean_depth_km=None, ocean_fraction=0.0, land_fraction=1.0)
    pb.generate(body)
    assert body.crust_p2o5_wt > 0.0 and body.crust_merrillite_share > body.crust_apatite_share
    assert body.ocean_phosphate_umol_l is None and body.phosphorus is None


def test_a_body_with_no_stored_phosphorus_or_a_giant_gets_nothing(earth_heat):
    for body in (_body(rock_p2o5_fraction=None), _body(body_type="g")):
        pb.generate(body)
        assert all(getattr(body, name) is None for name in pb.COLUMN_FIELDS)


def test_the_stored_columns_are_the_fields_and_each_has_a_unit():
    from planetgen.db import column_units
    for table in ("planets", "moons", "rogue_planets"):
        assert all((table, name) in column_units.COLUMNS or name == "phosphorus" for name in pb.COLUMN_FIELDS)
    assert column_units.unit_of("planets", "ocean_phosphate_umol_l") == "umol/L" and column_units.unit_of("moons", "crust_p2o5_wt") == "pct"


# --- generated bodies -----------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def sun():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("phosphorus-budget:star"):
        star = Star(cfg)
    star.age = 4.5
    return star


def test_a_generated_rocky_planet_carries_the_chain_and_makes_no_draw_of_the_run(sun):
    planets = [class_support.draw_planet(sun, 1.0, seed, SystemConfig(), mass=EARTH, zone_override="c", moon_count=0) for seed in range(4)]
    for planet in planets:
        assert planet.crust_p2o5_wt > 0.0
        assert math.fsum((planet.crust_apatite_share, planet.crust_merrillite_share, planet.crust_schreibersite_share)) == pytest.approx(1.0)
        if planet.ocean_class is None:
            assert planet.ocean_phosphate_umol_l is None
        else:
            assert pb.FLOOR_UMOL_L <= planet.ocean_phosphate_umol_l <= pb.CEILING_UMOL_L and planet.phosphorus == pb.phosphorus_class(planet.ocean_phosphate_umol_l)
    again = planets[0]
    before = {name: getattr(again, name) for name in pb.COLUMN_FIELDS}
    pb.generate(again)
    assert {name: getattr(again, name) for name in pb.COLUMN_FIELDS} == before
