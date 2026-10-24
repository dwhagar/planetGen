# tests/test_equipment.py

"""
GEN.655: the survival equipment line (`planetgen.physics.equipment`, docs/design/survival-equipment-fields.md). Named
worlds, the severity bands, the radiation ladder, the name and the code, and the species fit.
"""

import pytest

from planetgen import tuning
from planetgen.physics import equipment
from planetgen.physics.habitability import World
from planetgen.physics.protection_limits import RADIATION_STEPS

AIR = {"O2": 0.21}


def _world(pressure=101.3, temperature=15.0, gases=None, dose=0.0024, gravity=1.0, parts=None, uv=0.0, ph=None):
    extra = {}
    if parts:
        extra = dict(dose_gcr_sv_yr=parts[0], dose_sep_sv_yr=parts[1], dose_ground_sv_yr=parts[2])
    return World(pressure_kpa=pressure, temperature_c=temperature, gases=AIR if gases is None else gases,
                 surface_dose_sv_yr=dose, gravity_ms2=9.80665 * gravity, relative_humidity=0.5, uv_index=uv, ph=ph,
                 **extra)


def _code(**changes):
    return equipment.code(equipment.compute(_world(**changes)))


def _line(**changes):
    return equipment.line(equipment.compute(_world(**changes)))


NAMED = {
    "Earth": (dict(gases={"O2": 0.21, "CO2": 0.0004}), "", ""),
    "Arctic": (dict(temperature=-35.0), "B2c", "Required: insulated suit."),
    "Everest": (dict(pressure=33.7, temperature=-30.0), "B2c S1", "Required: insulated suit + gas-mix mask."),
    "Mars": (dict(pressure=0.6, temperature=-60.0, gases={"CO2": 0.95}, dose=0.236, parts=(0.236, 0.0, 0.001)),
             "B4p F100 S2 H3 [R2]",
             "Required: hardened pressure suit + dust-tight seals. Recommended: hardened habitat."),
    "Venus": (dict(pressure=9200.0, temperature=462.0, gases={"CO2": 0.965}, dose=0.0001, parts=(0.0001, 0.0, 0.0)),
              "B5h S4", "Required: armoured suit."),
    "Titan": (dict(pressure=146.7, temperature=-179.0, gases={"CH4": 0.0}), "B3c S2",
              "Required: heated sealed suit + self-contained air."),
    "10 bar air": (dict(pressure=1000.0, gases={"O2": 0.5}), "B2d S3", "Required: heated diving suit + mixed-gas rig."),
    "UV": (dict(uv=200.0), "B2u", "Required: full-cover suit."),
}


@pytest.mark.parametrize("name", list(NAMED))
def test_the_named_worlds_get_the_line_physics_gives(name):
    changes, code, line = NAMED[name]
    assert _code(**changes) == code and _line(**changes) == line


def test_a_world_that_needs_nothing_lists_nothing():
    levels = equipment.compute(_world())
    assert equipment.items(levels) == [] and equipment.line(levels) == "" and equipment.code(levels) == ""
    assert equipment.required_names(levels) == [] and equipment.recommended_names(levels) == []


def test_a_pressure_below_the_suit_limit_needs_a_pressure_suit_and_a_thin_air_only_recommends_one():
    thin = equipment.compute(_world(pressure=18.0, temperature=10.0))
    assert equipment.required_names(thin) == ["gas-mix mask"] and equipment.recommended_names(thin) == ["pressure suit"]
    vacuum = equipment.compute(_world(pressure=0.0, gases={}))
    assert "pressure suit" in equipment.required_names(vacuum)
    assert vacuum["eq_body"] == 4 and vacuum["eq_supply"] >= 2  # a suit that seals carries its own air


def test_a_toxic_gas_past_its_idlh_asks_for_its_own_air():
    assert _code(gases={"CO": 0.05, "O2": 0.2}, pressure=100.0, temperature=5.0) == "S2"
    assert _code(gases={"H2S": 0.0003, "SO2": 0.0005, "O2": 0.21}, pressure=100.0) == "S2"


def test_gravity_asks_for_a_frame_and_then_a_vehicle():
    assert equipment.compute(_world(gravity=1.1))["eq_mobility"] == 0
    assert equipment.compute(_world(gravity=1.3))["eq_mobility"] == 1
    assert equipment.compute(_world(gravity=2.0))["eq_mobility"] == 2
    assert "powered frame" in equipment.required_names(equipment.compute(_world(gravity=2.0)))
    assert equipment.compute(_world(gravity=3.5))["eq_body"] == 6


def test_the_radiation_ladder_runs_garment_vest_hardened_suit_shelter_habitat():
    assert [step.key for step in RADIATION_STEPS] == ["garment", "vest", "hardened", "shelter", "habitat"]
    assert equipment.compute(_world(dose=5.0, parts=(0.2, 4.8, 0.001), pressure=1e-9, gases={}))["eq_rad_suit"] == 3
    heavy = equipment.compute(_world(dose=2000.0, parts=(500.0, 1500.0, 0.001), pressure=1e-9, gases={}))
    assert heavy["eq_shelter"] == 2 and heavy["eq_rad_suit"] == 3
    assert heavy["outside_hours"] is not None and heavy["outside_hours"] < 2000.0
    assert heavy["shelter_column_g_cm2"] > 300.0


def test_a_dose_past_a_hundred_sieverts_needs_a_hardened_habitat():
    assert equipment.compute(_world(dose=50000.0))["eq_shelter"] == 2


def test_a_hardened_suit_is_a_hardened_version_of_the_suit_the_world_asks_for():
    required = equipment.required_names(equipment.compute(
        _world(dose=5.0, parts=(0.2, 4.8, 0.001), pressure=1e-9, gases={})))
    assert "hardened pressure suit" in required and "pressure suit" not in required


def test_the_summary_is_none_for_an_unscored_body_and_carries_the_code_and_items():
    assert equipment.summary({}) is None
    summary = equipment.summary(equipment.compute(_world(pressure=33.7, temperature=-30.0)))
    assert summary["code"] == "B2c S1" and summary["required"] == ["insulated suit", "gas-mix mask"]
    assert {item["severity"] for item in summary["items"]} == {"required"}


def test_every_item_the_tables_can_name_has_a_tech_level():
    names = {*equipment.BODY_NAMES.values(), *equipment.BODY_DEFAULT.values(), *equipment.SUPPLY_NAMES.values(),
             *equipment.GAS_FILTER_NAMES.values(), *equipment.SHELTER_NAMES.values(), *equipment.RAD_SUIT_NAMES.values(),
             "powered frame", "dust filter", "dust-tight seals", "CO2 scrubber"}
    for name in names:
        assert name in tuning.SPECIES_EQUIPMENT_TECH, name


def test_the_tech_a_species_needs_is_the_highest_of_its_required_items():
    levels = equipment.compute(_world(pressure=0.6, temperature=-60.0, gases={"CO2": 0.95}, dose=0.236,
                                      parts=(0.236, 0.0, 0.001)))
    assert equipment.tech_needed(levels, tuning.SPECIES_EQUIPMENT_TECH) == 3.5  # the hardened suit
    assert equipment.tech_needed(equipment.compute(_world()), tuning.SPECIES_EQUIPMENT_TECH) == 0.0
    vehicle = equipment.compute(_world(pressure=9200.0, temperature=600.0, gases={"CO2": 0.9}))
    assert equipment.tech_needed(vehicle, tuning.SPECIES_EQUIPMENT_TECH) == 99.0


def test_the_search_ticks_select_what_the_line_lists():
    for value, _label, clause in equipment.TICKS:
        assert equipment.tick_condition("p", [value]).count("p.") >= 1
        assert "{a}" not in equipment.tick_condition("p", [value])
        assert clause.format(a="p") in equipment.tick_condition("p", [value])


def test_where_no_garment_stops_cosmic_rays_the_hardened_suit_is_required_and_the_habitat_recommended():
    """GEN.671 (Boss, 2026-10-23): a hardened suit in the game works against cosmic rays, so Mars and the Moon need it."""
    for changes in (NAMED["Mars"][0], dict(pressure=1e-9, temperature=-20.0, gases={}, dose=0.5, parts=(0.5, 0.0, 0.0005))):
        levels = equipment.compute(_world(**changes))
        assert "hardened pressure suit" in equipment.required_names(levels), changes
        assert equipment.recommended_names(levels) == ["hardened habitat"], changes


def test_the_required_and_recommended_lines_leave_an_empty_list_off():
    nothing = equipment.compute(_world())
    assert equipment.lines(nothing) == ("", "")  # a list with nothing in it is left off, never "none"
    assert equipment.summary(nothing)["required_line"] == "" and equipment.summary(nothing)["recommended_line"] == ""
    mars = equipment.summary(equipment.compute(_world(**NAMED["Mars"][0])))
    assert mars["required_line"] == "Required: hardened pressure suit + dust-tight seals."
    assert mars["recommended_line"] == "Recommended: hardened habitat."
