# tests/test_evaporite_sequence.py

"""
GEN.1366: the salts a brine leaves, in the order of its chemical divides, and where the surface draw reads them.
"""

import math

import pytest

from planetgen.physics import evaporite_sequence as ev, hydrosphere, materials, rock_types, surface_composition as sc
from planetgen.util import draw

BASALT = {"CaO": 11.0, "MgO": 8.0, "Na2O": 2.5, "K2O": 0.3}
"""dict: A basaltic crust's four oxides, weight percent (a test input, never a body's)."""


def _sequence(cls, temperature_k=290.0, **kwargs):
    return ev.allocate(ev.brine_ions(BASALT, cls), temperature_k, acid=cls == "acid sulfate", **kwargs)


def test_every_row_the_sequence_can_make_is_a_row_of_the_table():
    assert not [row for row in ev.ROW_IDS if not ev.row_exists(row)]
    for row in ev.ROW_IDS:
        assert materials.material(row).formula, row


def test_the_gypsum_anhydrite_line_is_hardies_table():
    """[R5143: p. 181]: 12, 18, 23, 30, 39, 45, 50, 55 and 58 degrees C."""
    for celsius, activity in ((12, 0.725), (18, 0.750), (23, 0.770), (30, 0.800), (39, 0.845), (45, 0.880), (50, 0.915), (55, 0.960), (58, 1.000)):
        assert ev.gypsum_anhydrite_aw(celsius + 273.15) == pytest.approx(activity, abs=1e-9), celsius
    assert ev.gypsum_anhydrite_aw(333.15) == 1.0                                       # past 58 degrees C anhydrite is stable in water of any activity
    assert ev.gypsum_anhydrite_aw(273.15) < ev.gypsum_anhydrite_aw(285.15)             # colder, the line falls further


def test_a_halite_saturated_brine_makes_anhydrite_above_eighteen_degrees_and_gypsum_below():
    """[R5143]: the gypsum in a brine saturated with halite dehydrates above about 18 degrees C."""
    assert ev.calcium_sulfate_row(18.5 + 273.15, ev.HALITE_SATURATED_AW) == "anhydrite"
    assert ev.calcium_sulfate_row(17.0 + 273.15, ev.HALITE_SATURATED_AW) == "gypsum"
    assert ev.calcium_sulfate_row(25.0 + 273.15, 1.0) == "gypsum"                      # pure water: gypsum to 58 degrees C
    assert ev.calcium_sulfate_row(60.0 + 273.15, 1.0) == "anhydrite"


def test_the_hydrates_follow_the_temperature():
    assert ev.sodium_sulfate_row(300.0) == "mirabilite" and ev.sodium_sulfate_row(310.0) == "thenardite"
    assert [ev.magnesium_sulfate_row(t) for t in (300.0, 330.0, 360.0)] == ["epsomite", "hexahydrite", "kieserite"]
    assert [ev.sodium_carbonate_row(t, 1.0e3) for t in (300.0, 340.0, 400.0)] == ["trona", "thermonatrite", "natrite"]
    assert ev.sodium_carbonate_row(300.0, 1.0e5) == "nahcolite"
    assert ev.magnesium_carbonate_row(300.0) == "hydromagnesite" and ev.magnesium_carbonate_row(350.0) == "magnesite"


@pytest.mark.parametrize("cls", ev.ANION_PARAMS)
def test_the_ions_are_charge_balanced_and_every_step_conserves_them(cls):
    ions = ev.brine_ions(BASALT, cls)
    cation_charge = 2.0 * ions["Ca"] + 2.0 * ions["Mg"] + ions["Na"] + ions["K"]
    anion_charge = ions["Cl"] + 2.0 * ions["SO4"] + 2.0 * ions["CO3"]
    assert anion_charge == pytest.approx(cation_charge)
    result = _sequence(cls, pressure_pa=1.0e4)
    assert math.fsum(result.shares.values()) == pytest.approx(1.0)
    assert all(row in ev.FORMULAS for row in result.shares)
    assert all(moles >= 0.0 for moles in result.moles.values())


def test_a_neutral_water_takes_its_calcium_as_calcite_and_leaves_no_calcium_sulfate():
    result = _sequence("neutral")
    assert "calcite" in result.shares and "gypsum" not in result.shares and "anhydrite" not in result.shares
    assert {"halite", "epsomite"} <= set(result.shares)                                 # the divide of a calcium-poor water: the magnesium sulfate then the chloride


def test_a_chloride_brine_left_with_calcium_ends_in_calcium_chloride():
    result = _sequence("chloride brine", pressure_pa=0.0)
    assert "antarcticite" in result.shares and result.shares["antarcticite"] > 0.3
    assert result.calcium_sulfate in ("gypsum", "anhydrite")


def test_an_acid_sulfate_water_makes_jarosite_and_alunite_and_much_calcium_sulfate():
    result = _sequence("acid sulfate")
    assert {"jarosite", "alunite"} <= set(result.shares)
    assert result.shares.get("gypsum", 0.0) + result.shares.get("anhydrite", 0.0) > 0.3
    assert "calcite" not in result.shares and "nahcolite" not in result.shares


def test_a_soda_water_takes_the_alkaline_path_with_sodium_carbonate_and_magnesium_carbonate():
    cold = _sequence("soda", 290.0)
    assert {"calcite", "hydromagnesite", "trona"} <= set(cold.shares)
    assert "bischofite" not in cold.shares and "antarcticite" not in cold.shares
    hot = _sequence("soda", 400.0)
    assert {"magnesite", "natrite"} <= set(hot.shares) and "trona" not in hot.shares
    assert _sequence("soda", 290.0, pco2_pa=1.0e5).shares.get("nahcolite", 0.0) > 0.0


def test_the_calcium_sulfate_is_anhydrite_when_halite_forms_and_the_water_is_warm():
    warm, cold = _sequence("acid sulfate", 300.0), _sequence("acid sulfate", 275.0)
    assert warm.calcium_sulfate == "anhydrite" and "anhydrite" in warm.shares
    assert cold.calcium_sulfate == "gypsum" and "gypsum" in cold.shares


def test_the_magnesium_sulfate_changes_with_the_temperature():
    assert "epsomite" in _sequence("acid sulfate", 290.0).shares
    assert "kieserite" in _sequence("acid sulfate", 360.0).shares


def test_perchlorate_needs_air_and_no_ozone_shield():
    assert "mg_perchlorate" in _sequence("neutral", pressure_pa=600.0).shares
    assert "mg_perchlorate" not in _sequence("neutral", pressure_pa=0.0).shares
    assert "mg_perchlorate" not in _sequence("neutral", pressure_pa=1.0e5, o2_shield=True).shares


def test_a_rock_without_the_four_oxides_leaves_nothing_and_a_body_without_rock_leaves_halite():
    assert ev.allocate(ev.brine_ions({"SiO2": 100.0}, "neutral"), 290.0).shares == {}
    assert ev.salts_of(None, "neutral", 290.0) == {"halite": 1.0}


def test_an_unknown_class_takes_the_neutral_one():
    assert ev.brine_ions(BASALT, None) == ev.brine_ions(BASALT, "neutral") == ev.brine_ions(BASALT, "ice-sealed")


def test_the_classes_are_the_hydrospheres():
    assert set(ev.ANION_PARAMS) | {"ice-sealed"} == set(hydrosphere.OCEAN_CLASSES)


def test_a_mole_of_each_solid_has_its_formulas_mass():
    result = ev.allocate({"Ca": 1.0, "SO4": 1.0}, 280.0)
    assert result.moles == {"gypsum": pytest.approx(1.0)}
    assert result.shares == {"gypsum": pytest.approx(1.0)}
    assert result.calcium_sulfate == "gypsum"


# --- the surface draw -----------------------------------------------------------------------------------------------------------

def _inputs(**kwargs):
    base = dict(surface_temperature_k=288.0, pressure_pa=1.0e5, o2_partial_pa=2.1e4, partial_pressures_pa={"CO2": 40.0}, water_mass_fraction=2e-4, rock_feo_fraction=0.08,
                hydrosphere="surface ocean", ocean_fraction=0.71, ocean_class="neutral", internal_heat_flux_w_m2=0.08, mobile_lid=True,
                rock_oxides=rock_types.rock_oxides(0.08, 0.05, 1.2, 0.1))
    base.update(kwargs)
    return sc.SurfaceInputs(**base)


def test_the_ocean_holds_the_salts_of_its_own_brine_without_the_seafloors_calcite():
    draw.set_run_seed(11)
    result = sc.draw_surface(_inputs())
    salts = {row: share for row, share in result.components.items() if row in ev.ROW_IDS}
    assert "calcite" not in salts or salts["calcite"] < 0.71 * sc.RANGES["ocean_salt"]
    assert math.fsum(share for row, share in salts.items() if row != "calcite") >= 0.71 * sc.RANGES["ocean_salt"] - 1e-12
    assert math.fsum(result.components.values()) == pytest.approx(1.0)


def test_a_soda_ocean_and_a_chloride_brine_dry_to_different_ground():
    soda, chloride = sc.ocean_salts(_inputs(ocean_class="soda")), sc.ocean_salts(_inputs(ocean_class="chloride brine"))
    assert soda != chloride
    assert math.fsum(soda.values()) == pytest.approx(1.0) and math.fsum(chloride.values()) == pytest.approx(1.0)


def test_the_evaporite_ground_is_the_dried_brine_not_the_old_recipe():
    inputs = _inputs()
    dried = sc.dried_salts(inputs)
    assert dried and set(dried) != {"halite", "gypsum", "calcite"}
    wet = _inputs(ocean_class="acid sulfate")
    assert "jarosite" in sc.dried_salts(wet)


def test_a_body_with_no_silicate_rock_dries_to_halite():
    inputs = _inputs(rock_oxides=None, carbon_rock=True)
    assert sc.dried_salts(inputs) == {"halite": 1.0}
    assert sc.ocean_salts(inputs) == {"halite": 1.0}
