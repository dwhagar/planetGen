# tests/test_ore_deposits.py

"""GEN.1380: the three derived crust fields and the three magmatic deposit types of the ore design (docs/design/composition/ore-processes.md)."""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import interior_heat, ore_deposits as od
from planetgen.util import draw
from tests import class_support

EARTH = 5.972e24
OXIDES = {"rock_nio_fraction": 0.0025, "rock_cr2o3_fraction": 0.0038, "rock_tio2_fraction": 0.0020}


def _state(**over):
    """An Earth-like state: mobile lid, wet, 4.5 Gyr, an intermediate mantle with the Earth's bulk-silicate oxides."""
    base = dict(gravity_g=1.0, lid=interior_heat.MOBILE, mantle_k=od.REFERENCE_MANTLE_K, age_gy=4.5, magma_ratio=1.0, redox="intermediate", wet=True, sulfur_ratio=1.0,
                felsic_share=0.3, oxides=dict(OXIDES), icy=False)
    base.update(over)
    return od.OreState(**base)


def _rows(state, seed="a"):
    sulfur = od.crust_sulfur_index(state)
    return {kind.id: od.assess(kind, state, sulfur, draw.Stream(f"{seed}:{kind.id}")) for kind in od.TYPES}


# --- the derived fields ------------------------------------------------------------------------------------------------
def test_the_crust_is_thick_on_a_small_body_and_thin_on_a_hot_mantle():
    earth, mars_like = od.crust_thickness_km(_state()), od.crust_thickness_km(_state(gravity_g=0.38, lid=interior_heat.STAGNANT))
    assert earth == pytest.approx(17.0) and 40.0 < mars_like < 70.0                    # Mars about 50 km
    venus_like = od.crust_thickness_km(_state(gravity_g=0.90, lid=interior_heat.STAGNANT))
    assert 20.0 < venus_like < 36.0                                                    # Venus about 30 km
    assert od.crust_thickness_km(_state(mantle_k=2800.0)) < earth < od.crust_thickness_km(_state(mantle_k=1800.0))
    assert od.crust_thickness_km(_state(gravity_g=0.0, mantle_k=1.0)) == od.CRUST_KM_RANGE[1]


def test_a_craton_needs_a_mobile_lid_age_and_granitoid_ground():
    assert od.craton_fraction(_state(), 1.0) == pytest.approx(od.CRATON_SHARE)
    assert od.craton_fraction(_state(lid=interior_heat.STAGNANT), 1.0) == 0.0
    assert od.craton_fraction(_state(age_gy=1.9), 1.0) == 0.0
    assert od.craton_fraction(_state(felsic_share=0.0), 1.0) == 0.0
    assert 0.0 < od.craton_fraction(_state(age_gy=3.0), 1.0) < od.craton_fraction(_state(), 1.0)
    assert od.craton_fraction(_state(), 1.5) > od.craton_fraction(_state(), 0.5)


def test_the_sulfur_index_follows_water_age_and_the_stars_sulfur():
    wet, dry = od.crust_sulfur_index(_state()), od.crust_sulfur_index(_state(wet=False))
    assert wet == pytest.approx(0.8) and dry == pytest.approx(od.SULFUR_PRIMORDIAL)
    assert od.crust_sulfur_index(_state(age_gy=0.5)) < wet
    assert od.crust_sulfur_index(_state(wet=False, sulfur_ratio=2.0)) == pytest.approx(2 * od.SULFUR_PRIMORDIAL)
    assert 0.0 <= od.crust_sulfur_index(_state(sulfur_ratio=50.0)) <= 1.0


# --- the gates ---------------------------------------------------------------------------------------------------------
def test_an_earth_like_body_admits_all_three_magmatic_types_and_stores_each_row():
    seen = {kind: set() for kind in od.TYPE_BY_ID}
    for index in range(60):
        for kind, row in _rows(_state(), seed=str(index)).items():
            seen[kind].add(row["status"])
            assert set(row) == set(od.ROW_FIELDS)
            if row["status"] == od.PRESENT:
                low, high = od.TYPE_BY_ID[kind].grade_range
                assert low <= row["grade"] <= high and row["endowment_kg"] > 0.0
                assert 0.0 < row["age_window_gyr_start"] < row["age_window_gyr_end"] <= 4.5
            else:
                assert row["grade"] is None and row["endowment_kg"] is None and row["capped_by_host"] == 0
    assert all(od.PRESENT in statuses for statuses in seen.values())


def test_no_magmatic_ore_without_a_large_magma_volume_or_on_an_icy_body():
    for state in (_state(magma_ratio=od.MAGMA_MIN_RATIO / 2), _state(icy=True)):
        assert {row["status"] for row in _rows(state).values()} == {od.ABSENT_BY_STATE}


def test_the_sulfide_needs_sulfur_and_a_mantle_that_is_not_oxidised_and_the_chromite_a_reduced_one():
    dry = _rows(_state(wet=False))
    assert dry[od.NICKEL_COPPER_PGE]["status"] == od.ABSENT_BY_STATE                    # primordial sulfur only
    oxidised = _rows(_state(redox="oxidized"))
    assert oxidised[od.NICKEL_COPPER_PGE]["status"] == oxidised[od.CHROMITITE]["status"] == od.ABSENT_BY_STATE
    assert any(_rows(_state(redox="oxidized"), seed=str(i))[od.IRON_TITANIUM_VANADIUM]["status"] == od.PRESENT for i in range(20))   # the oxide type takes any mantle
    assert any(_rows(_state(wet=False, sulfur_ratio=3.0), seed=str(i))[od.NICKEL_COPPER_PGE]["status"] == od.PRESENT for i in range(40))   # a sulfur-rich star can make up for a dry body


def test_a_poor_host_gives_a_lean_ore_or_none_a_metal_poor_star_gives_poorer_ore():
    kind = od.TYPE_BY_ID[od.CHROMITITE]
    rich, lean, poor = (_state(oxides=dict(OXIDES, rock_cr2o3_fraction=value)) for value in (0.0060, 0.0030, 0.0010))
    assert od.grade_cap(kind, rich) > od.grade_cap(kind, lean) > od.grade_cap(kind, poor)
    assert od.grade_cap(kind, poor) < kind.grade_range[0]
    assert all(row["status"] == od.ABSENT_BY_STATE for row in (_rows(poor, str(i))[od.CHROMITITE] for i in range(10)))
    capped = [row for row in (_rows(lean, str(i))[od.CHROMITITE] for i in range(80)) if row["status"] == od.PRESENT]
    assert capped and all(row["grade"] <= od.grade_cap(kind, lean) + 1e-12 for row in capped) and any(row["capped_by_host"] for row in capped)


def test_the_earths_bulk_silicate_caps_the_grades_inside_the_ranges():
    for kind in od.TYPES:
        cap = od.grade_cap(kind, _state())
        assert kind.grade_range[0] <= cap <= 1.0, kind.id


def test_the_endowment_scales_with_the_magma_volume():
    mean = lambda ratio: sum((_rows(_state(magma_ratio=ratio), str(i))[od.IRON_TITANIUM_VANADIUM]["endowment_kg"] or 0.0) for i in range(120)) / 120
    assert mean(2.0) > 1.4 * mean(1.0) > 1.4 * 1.4 * mean(0.25) * 0.9


def test_a_type_has_its_own_stream():
    state = _state()
    first = _rows(state, "same")
    again = _rows(state, "same")
    assert first == again
    only = od.assess(od.TYPE_BY_ID[od.CHROMITITE], state, od.crust_sulfur_index(state), draw.Stream("same:chromitite"))
    assert only == first[od.CHROMITITE]                                                   # the other types' draws do not move it


def test_the_poisson_draw_has_the_mean_it_is_given():
    stream = draw.Stream("poisson")
    for mean in (0.5, 4.0, 60.0):
        values = [od.poisson(stream, mean) for _ in range(4000)]
        assert sum(values) / len(values) == pytest.approx(mean, rel=0.08)
    assert od.poisson(stream, 0.0) == 0


# --- a generated body --------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def sun():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("ore-deposits:star"):
        star = Star(cfg)
    star.age = 4.5
    return star


def _planet(star, mass_earth=1.0, seed=0, distance=1.0):
    return class_support.draw_planet(star, distance, seed, SystemConfig(), mass=mass_earth * EARTH, zone_override="c", moon_count=0)


def test_a_generated_planet_carries_the_fields_and_the_three_rows(sun):
    planet = _planet(sun)
    assert planet.crust_thickness_km is not None and planet.craton_fraction is not None and planet.crust_sulfur_index is not None
    assert [row["deposit_type"] for row in planet.deposit_types] == [kind.id for kind in od.TYPES]
    assert all(row["status"] in od.STATUSES for row in planet.deposit_types)
    assert 3.0 <= planet.crust_thickness_km <= 150.0 and 0.0 <= planet.craton_fraction <= 1.0 and 0.0 <= planet.crust_sulfur_index <= 1.0


def test_the_same_body_gets_the_same_record_and_the_run_stream_does_not_move(sun):
    planet = _planet(sun, seed=3)
    before = (planet.crust_thickness_km, planet.craton_fraction, planet.crust_sulfur_index, [dict(row) for row in planet.deposit_types])
    with draw.bound("ore-deposits:outer-a"):
        od.generate(planet)
        mark = draw.random()
    with draw.bound("ore-deposits:outer-a"):
        mark_without = draw.random()
    assert mark == mark_without
    assert before == (planet.crust_thickness_km, planet.craton_fraction, planet.crust_sulfur_index, planet.deposit_types)


def test_a_giant_and_a_body_without_content_are_not_assessed(sun):
    planet = _planet(sun)
    planet.body_type = "g"
    od.generate(planet)
    assert planet.deposit_types is None and planet.crust_thickness_km is None and planet.craton_fraction is None and planet.crust_sulfur_index is None


def test_the_magma_reference_is_an_earth_like_bodys_processed_mantle(sun):
    planet = _planet(sun, 1.0, seed=0)
    state = od.state_of(planet)
    assert 0.3 < state.magma_ratio < 3.0                                                 # an Earth-mass body at 4.5 Gyr is of the order of the reference
    small = od.state_of(_planet(sun, 0.03, seed=0))
    assert small.magma_ratio < od.MAGMA_MIN_RATIO < state.magma_ratio and math.isfinite(small.magma_ratio)
