# tests/test_superionic.py

"""
GEN.1422: the superionic boundaries of water and ammonia (`physics/superionic.py`, `physics/data/superionic_boundary.json`) and the deep ice layer of an ice giant
that reads them from the temperature of the model's adiabat (`giant_interior`, `giant_structure.Profile.temperature_at_pressure_k`).
"""

import pytest

from planetgen.physics import giant_interior as gi
from planetgen.physics import superionic
from planetgen.physics.elements import DataError

from tests.test_giant_interior import NEPTUNE, interior


def test_water_is_an_ice_under_its_floor_superionic_between_and_an_ionic_fluid_over_its_ceiling():
    assert superionic.state("water", 10.0, 1500.0) == "solid_or_molecular"          # under 20 GPa there is no superionic phase
    assert superionic.state("water", 25.0, 600.0) == "solid_or_molecular"
    assert superionic.state("water", 80.0, 3000.0) == "superionic"
    assert superionic.state("water", 50.0, 4500.0) == "ionic_fluid"
    assert superionic.state("water", 150.0, 8000.0) == "metallic_fluid"             # [R5172]: water turns metallic over 7000 K


def test_the_floor_follows_the_measured_points_and_the_widest_range_is_used():
    assert superionic.floor_k("water", 20.0) == pytest.approx(900.0)
    assert superionic.floor_k("water", 29.0) == pytest.approx(1300.0)
    assert superionic.floor_k("water", 60.0) == pytest.approx(1900.0)
    assert superionic.floor_k("water", 300.0) == pytest.approx(1900.0)               # held past the last point (STAND-IN SI1)
    assert superionic.floor_k("water", 19.9) is None
    assert superionic.ceiling_k("water", 105.0) == pytest.approx(4000.0)


def test_ammonia_has_its_own_boundary_and_turns_metallic_lower():
    assert superionic.state("ammonia", 100.0, 2500.0) == "superionic"
    assert superionic.state("ammonia", 100.0, 6000.0) == "metallic_fluid"
    assert superionic.state("ammonia", 10.0, 2500.0) == "solid_or_molecular"


def test_a_layer_takes_the_state_of_the_species_that_holds_more_of_it():
    assert superionic.layer_state({"water_liquid": 3.0, "nh3_ice": 1.0}, 100.0, 2500.0) == "superionic"
    assert superionic.layer_state({"water_liquid": 1.0, "nh3_ice": 3.0}, 100.0, 6000.0) == "metallic_fluid"
    assert superionic.layer_state({"rock": 1.0}, 100.0, 2500.0) is None


def test_a_bad_table_is_refused(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text('{"species": {"water": {"floor": [[20, 900], [10, 800]], "ceiling": [[20, 900]], "metallic_k": 7000}, "ammonia": {}}}')
    with pytest.raises(DataError):
        superionic.load(str(path))


def test_a_profile_without_an_adiabat_has_no_temperature_and_the_layers_fall_back_to_the_stand_in():
    from planetgen.physics import giant_structure

    assert giant_structure.Profile(1.0, 1.0, 100.0, True, 1e12, 5000.0).temperature_at_pressure_k(1e10) is None
    assert gi._bottom_temperature_k(lambda pressure_pa: None, 1e11) == superionic.FALLBACK_TEMPERATURE_K


def test_an_ice_giant_mantle_is_read_from_the_boundary_table_at_its_bottom():
    for seed in range(1, 9):
        body = interior(NEPTUNE, seed)
        mantle = body.layer("ice_mantle")
        assert mantle.state in ("superionic", "ionic_fluid")
        if mantle.state == "superionic":
            assert "ice_x" in mantle.components or "nh3_ice" in mantle.components or "water_liquid" not in mantle.components
        else:
            assert "ice_x" not in mantle.components


def test_the_stand_in_constants_agree_with_the_table():
    # STAND-INS SI1 to SI3 (section 6.233) name the points of the data file that are not measured
    water, ammonia = superionic.table()["water"], superionic.table()["ammonia"]
    assert water["floor"][-1][0] == superionic.FALLBACK_FLOOR_BEYOND_GPA
    assert (water["ceiling"][0][0], water["ceiling"][1][0]) == superionic.CEILING_INTERPOLATION_GPA
    assert ammonia["floor"][0] == superionic.AMMONIA_FLOOR_GPA_K
