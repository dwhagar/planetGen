# tests/test_ice_types.py

"""
GEN.1177: the strength of the other ices of an icy body next to ice Ih, and the naming of a body's ice type from its components. The numbers are the
stand-ins of the strength research (sections 2b, 9, 11c, 14); the tests pin the shapes the research gives (ordering, ranges), not a measurement.
"""

import copy

import pytest

from planetgen.physics import ice_types as it
from planetgen.physics import roundness, small_body_composition as sbc


def test_every_type_is_inside_its_range():
    for name, entry in it.TYPES.items():
        if "strength_ratio" in entry:
            low, high = entry["range"]
            assert low <= entry["strength_ratio"] <= high, name
        if "ratio_by_temperature_k" in entry:
            low, high = entry["range"]
            assert all(low <= ratio <= high for _, ratio in entry["ratio_by_temperature_k"]), name
        if "strength_pa" in entry:
            low, high = entry["range_pa"]
            assert low <= entry["strength_pa"] <= high, name
        assert entry["source"] and entry["tier"].startswith("T") and entry["tweak"].startswith("L"), name


def test_the_water_phases_follow_the_calculation_from_the_1997_laws():
    """GEN.1183: the ratio to ice Ih at the same temperature, from [R1190: Table 2] at 1e-15 per second."""
    assert it.strength_ratio("ice_ih") == 1.0
    assert it.strength_ratio("ice_ii", 150.0) == pytest.approx(1.7)
    assert it.strength_ratio("ice_ii", 180.0) == pytest.approx(1.0)
    assert all(5.0 <= it.strength_ratio("ice_ii", temperature) <= 6.5 for temperature in (200.0, 230.0, 250.0))
    assert all(0.45 <= it.strength_ratio("ice_iii", temperature) <= 0.6 for temperature in (240.0, 250.0, 256.0))     # weaker than Ih, not near zero
    assert all(7.0 <= it.strength_ratio("ice_v", temperature) <= 10.0 for temperature in (230.0, 250.0))
    assert all(4.0 <= it.strength_ratio("ice_vi", temperature) <= 5.5 for temperature in (150.0, 200.0, 230.0))
    assert it.strength_ratio("ice_vi", 180.0) == pytest.approx(1.3)


def test_ice_iii_is_the_weakest_of_the_high_pressure_phases_where_it_is_stable():
    for temperature in (240.0, 250.0, 256.0):
        assert it.strength_ratio("ice_iii", temperature) < min(it.strength_ratio(name, temperature) for name in ("ice_ii", "ice_v", "ice_vi"))
        assert it.strength_ratio("ice_iii", temperature) < 1.0


def test_a_water_phase_ratio_is_held_at_the_ends_of_its_table_and_is_log_linear_inside_it():
    for name in ("ice_ii", "ice_iii", "ice_v", "ice_vi"):
        points = it.TYPES[name]["ratio_by_temperature_k"]
        assert it.strength_ratio(name, points[0][0] - 50.0) == points[0][1]
        assert it.strength_ratio(name, points[-1][0] + 50.0) == points[-1][1]
        (t0, r0), (t1, r1) = points[0], points[1]
        assert it.strength_ratio(name, (t0 + t1) / 2.0) == pytest.approx((r0 * r1) ** 0.5)
    assert it.strength_pa("ice_v", 1.0e6, 250.0) == pytest.approx(7.2e6)


def test_methane_clathrate_is_twenty_to_thirty_times_ice():
    assert 20.0 <= it.strength_ratio("ch4_clathrate") <= 30.0
    assert it.strength_pa("ch4_clathrate", 2.0e6) == pytest.approx(it.strength_ratio("ch4_clathrate") * 2.0e6)


def test_solid_carbon_dioxide_has_its_own_strength_between_a_tenth_and_half_a_megapascal():
    assert 1.0e5 <= it.strength_pa("co2_ice", 5.0e6) <= 5.0e5       # a stronger ice Ih does not move it
    assert it.strength_pa("co2_ice", 5.0e6) == it.strength_pa("co2_ice", 1.0e5)
    assert it.strength_ratio("co2_ice") < 1.0


def test_ammonia_dihydrate_crosses_ice_at_135_k():
    assert it.strength_ratio("nh3_dihydrate", 135.0) == pytest.approx(1.0)
    assert it.strength_ratio("nh3_dihydrate", 110.0) == pytest.approx(2.0)
    assert it.strength_ratio("nh3_dihydrate", 170.0) == pytest.approx(0.3)
    assert it.strength_ratio("nh3_dihydrate", 100.0) == 2.0 and it.strength_ratio("nh3_dihydrate", 250.0) == 0.3     # held at the ends
    assert 1.0 < it.strength_ratio("nh3_dihydrate", 120.0) < 2.0 and 0.3 < it.strength_ratio("nh3_dihydrate", 150.0) < 1.0


def test_the_default_ice_is_the_ice_of_the_engine():
    assert it.ICE_IH_REFERENCE_PA == roundness.ICE_STRENGTH_PA


def test_an_unknown_type_is_refused():
    with pytest.raises(it.IceTypeError):
        it.strength_ratio("ice_ix")


@pytest.mark.parametrize("pressure,temperature,phase", [
    (None, None, "ice_ih"), (1.0e6, 150.0, "ice_ih"), (2.0e8, 150.0, "ice_ih"),
    (3.0e8, 150.0, "ice_ii"), (3.0e8, 250.0, "ice_iii"), (5.0e8, 150.0, "ice_ii"), (5.0e8, 260.0, "ice_v"),
    (1.0e9, 200.0, "ice_vi"), (5.0e9, 200.0, "ice_vi"),
])
def test_water_phase_from_the_central_pressure(pressure, temperature, phase):
    assert it.water_phase(pressure, temperature) == phase


def test_a_body_with_no_ice_has_no_ice_type():
    assert it.ice_type_of({"forsterite": 0.6, "iron": 0.4}) is None


def test_water_ice_is_named_by_phase_and_other_rows_by_their_own_name():
    assert it.ice_type_of({"water_ice": 0.5, "forsterite": 0.5}) == "ice_ih"
    assert it.ice_type_of({"water_ice": 0.5, "forsterite": 0.5}, central_pressure_pa=1.0e9) == "ice_vi"
    assert it.ice_type_of({"ice_ii": 0.5, "forsterite": 0.5}) == "ice_ii"
    assert it.ice_type_of({"ice_vii": 0.5, "forsterite": 0.5}) == "ice_vi"
    assert it.ice_type_of({"co2_ice": 0.6, "water_ice": 0.1, "forsterite": 0.3}) == "co2_ice"


def test_methane_is_a_clathrate_when_water_cages_enough_of_it():
    assert it.ice_type_of({"water_ice": 0.5, "ch4_ice": 0.01, "forsterite": 0.49}) == "ice_ih"     # a comet's trace of methane
    assert it.ice_type_of({"water_ice": 0.4, "ch4_ice": 0.1, "forsterite": 0.5}) == "ch4_clathrate"
    assert it.ice_type_of({"ch4_ice": 0.1, "forsterite": 0.9}) == "ch4_ice"                         # no water to cage it: solid methane


def test_ammonia_with_water_is_the_dihydrate():
    assert it.ice_type_of({"water_ice": 0.4, "nh3_ice": 0.1, "forsterite": 0.5}) == "nh3_dihydrate"
    assert it.ice_type_of({"water_ice": 0.5, "nh3_ice": 0.01, "forsterite": 0.49}) == "ice_ih"


def test_a_makeup_is_accepted_and_a_drawn_comet_has_an_ice_type():
    from planetgen.util import draw
    with draw.bound(draw.Stream(3)):
        comet = sbc.draw_comet("jupiter_family")
    assert it.ice_type_of(comet) in set(it.ICE_TYPES) | {"co_ice", "n2_ice"}


def test_a_malformed_table_is_refused(tmp_path):
    import json
    document = json.load(open(it.DATA_FILE, encoding="utf-8"))
    broken = copy.deepcopy(document)
    broken["ice_types"]["types"]["ice_ii"]["strength_pa"] = 1.0           # two strengths at once
    path = tmp_path / "broken.json"
    path.write_text(json.dumps(broken), encoding="utf-8")
    with pytest.raises(it.IceTypeError):
        it._load(str(path))
    del broken["ice_types"]["types"]["ch4_clathrate"]
    path.write_text(json.dumps(broken), encoding="utf-8")
    with pytest.raises(it.IceTypeError):
        it._load(str(path))
