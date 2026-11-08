# tests/test_isotope_record.py

"""GEN.1377, GEN.1420 (part 1): the isotope record of a body: the elements it holds, the natural mix by default and the physical laws of D/H, 13C, 15N, 40Ar and the decaying parents."""

import math
from types import SimpleNamespace

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import elements, isotope_record as ir
from planetgen.util import draw
from tests import class_support

EARTH = 5.972e24


def _star(fe_h=0.0, age=4.5, snow_line_au=2.7):
    return SimpleNamespace(metallicity_fe_h=fe_h, age=age, snow_line_au=snow_line_au)


def _body(**over):
    """A bare body with the fields the record reads."""
    base = dict(star=_star(), is_moon=False, mass=EARTH, radius=6371.0, distance=1.0, formation_distance_au=1.0, gravity=1.0,
                rock_mass_fraction=0.67, metal_mass_fraction=0.32, water_mass_fraction=0.0, ice_mass_fraction=0.0, envelope_mass_fraction=0.0,
                envelope_helium_fraction=None, envelope_z_fraction=None, carbon_rock=False, gas_partials_pa={}, core_ni_fraction=0.06, core_s_fraction=0.02)
    base.update(over)
    return SimpleNamespace(**base)


# --- the laws ----------------------------------------------------------------------------------------------------------
def test_the_astration_law_gives_the_primordial_value_for_no_metals_and_the_protosolar_one_at_the_suns():
    assert ir.gas_d_h(_star(fe_h=0.0)) == pytest.approx(ir.PROTOSOLAR_D_H, rel=1e-9)
    assert ir.gas_d_h(_star(fe_h=-3.0)) == pytest.approx(ir.PRIMORDIAL_D_H, rel=0.01)
    assert ir.gas_d_h(_star(fe_h=-1.0)) > ir.gas_d_h(_star(fe_h=0.0)) > ir.gas_d_h(_star(fe_h=0.4))


def test_the_ice_law_passes_through_the_two_published_points_and_grows_colder():
    ocean = ir.ice_factor(ir.ICE_REFERENCE_K, ir.D_H_ENRICHMENT_REFERENCE, ir.D_H_T0_K) * ir.PROTOSOLAR_D_H
    comet = ir.ice_factor(ir.COLD_ICE_K, ir.D_H_ENRICHMENT_REFERENCE, ir.D_H_T0_K) * ir.PROTOSOLAR_D_H
    assert ocean == pytest.approx(1.56e-4, rel=1e-9) and comet == pytest.approx(5.3e-4, rel=1e-9)
    assert 20.0 < ir.D_H_T0_K < 80.0 and 10.0 < ir.N_T0_K < 40.0
    temperatures = (182.0, 120.0, 80.0, 50.0, 30.0)
    factors = [ir.ice_factor(t, ir.D_H_ENRICHMENT_REFERENCE, ir.D_H_T0_K) for t in temperatures]
    assert factors == sorted(factors)
    nitrogen = ir.ice_factor(ir.COLD_ICE_K, ir.NITROGEN_ENRICHMENT_REFERENCE, ir.N_T0_K)
    assert 441.0 / nitrogen == pytest.approx(140.0, rel=1e-9)                     # comets: 14N/15N about 140


def test_the_carbon_and_nitrogen_ratios_move_between_the_local_interstellar_and_the_birth_values_with_age_and_metals():
    sun_age = ir.tuning.STELLAR_SOLAR_SYSTEM_AGE_GY
    assert 1.0 / ir.gas_c13_over_c12(_star(0.0, sun_age)) == pytest.approx(89.0, rel=1e-6)
    assert 1.0 / ir.gas_n15_over_n14(_star(0.0, sun_age)) == pytest.approx(441.0, rel=1e-6)
    assert 1.0 / ir.gas_c13_over_c12(_star(0.0, 0.0)) == pytest.approx(68.0, rel=1e-6)
    assert 1.0 / ir.gas_n15_over_n14(_star(0.0, 0.0)) == pytest.approx(400.0, rel=1e-6)
    # a metal-rich star formed nearer the center, where 13C and 15N are more common
    assert ir.gas_c13_over_c12(_star(0.3, 4.5)) > ir.gas_c13_over_c12(_star(0.0, 4.5)) > ir.gas_c13_over_c12(_star(-0.5, 4.5))
    assert ir.gas_n15_over_n14(_star(0.3, 4.5)) > ir.gas_n15_over_n14(_star(-0.5, 4.5))


def test_a_warm_body_takes_its_ice_from_the_snow_line_and_a_cold_one_is_enriched():
    warm, cold = _body(water_mass_fraction=0.001, formation_distance_au=1.0), _body(water_mass_fraction=0.3, formation_distance_au=2.7 * (182.0 / 30.0) ** 2)
    assert ir._ice_temperature_k(warm) == ir.ICE_REFERENCE_K
    assert ir._ice_temperature_k(cold) == pytest.approx(30.0, rel=1e-6)
    stream = draw.Stream("t")
    assert ir.deuterium_fraction(cold, stream) > 2.0 * ir.deuterium_fraction(warm, draw.Stream("t")) > 0.0


def test_a_hydrogen_envelope_keeps_the_gas_value_and_water_is_enriched():
    giant = _body(envelope_mass_fraction=0.8, envelope_helium_fraction=0.25, envelope_z_fraction=0.02, water_mass_fraction=0.0)
    wet = _body(water_mass_fraction=0.05)
    d_giant, d_wet = ir.deuterium_fraction(giant, draw.Stream("g")), ir.deuterium_fraction(wet, draw.Stream("w"))
    assert d_giant == pytest.approx(ir.gas_d_h(giant.star), rel=0.001)
    assert d_wet > 3.0 * d_giant
    assert ir.deuterium_fraction(_body(), draw.Stream("n")) is None                # no hydrogen reservoir at all


def test_the_nitrogen_of_the_envelope_is_the_gas_value_and_the_solids_are_enriched():
    solids, gas = _body(), _body(envelope_mass_fraction=1.0)
    assert ir.nitrogen_15_fraction(solids) > ir.nitrogen_15_fraction(gas)
    assert 1.0 / ir.nitrogen_15_fraction(solids) - 1.0 == pytest.approx(274.0, rel=0.02)          # the terrestrial 14N/15N at the snow line's temperature


# --- the elements a body holds -----------------------------------------------------------------------------------------
def test_a_rocky_body_holds_its_rock_core_and_air_elements():
    body = _body(rock_na2o_fraction=0.003, rock_k2o_fraction=0.0003, rock_p2o5_fraction=0.0, gas_partials_pa={"N2": 7.8e4, "O2": 2.1e4, "Ar": 930.0, "H2O": 100.0}, water_mass_fraction=1e-4)
    held = ir.elements_held(body)
    assert {"O", "Mg", "Si", "Al", "Ca", "Fe", "Na", "K", "Ni", "S", "N", "Ar", "H"} <= set(held)
    assert "P" not in held and "He" not in held                                    # P2O5 is zero; no envelope
    assert list(held) == sorted(held, key=lambda symbol: elements.element(symbol).z)


def test_an_envelope_gives_hydrogen_and_helium_and_a_carbon_rock_gives_carbon_and_silicon():
    giant = _body(rock_mass_fraction=0.0, metal_mass_fraction=0.0, envelope_mass_fraction=1.0, envelope_z_fraction=0.03, core_ni_fraction=None, core_s_fraction=None)
    assert set(ir.elements_held(giant)) == {"H", "He", "C", "N", "O"}
    carbon = _body(carbon_rock=True, metal_mass_fraction=0.0, core_ni_fraction=None, core_s_fraction=None)
    assert set(ir.elements_held(carbon)) == {"C", "Si"}
    assert ir.elements_held(SimpleNamespace()) == ()


# --- a generated body --------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def sun():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("isotope-record:star"):
        star = Star(cfg)
    star.age = 4.5
    return star


def _planet(star, mass_earth=1.0, seed=0, distance=1.0):
    return class_support.draw_planet(star, distance, seed, SystemConfig(), mass=mass_earth * EARTH, zone_override="c", moon_count=0)


def test_a_generated_planet_carries_a_record_that_adds_up_per_element(sun):
    planet = _planet(sun)
    rows = planet.isotopes
    assert rows and all(set(row) == set(ir.ROW_FIELDS) and row["origin"] in ir.ORIGINS for row in rows)
    assert {row["isotope"] for row in rows if row["initial_ratio"] is not None} == {"K-40", "U-235"}
    by_element = {}
    for row in rows:
        by_element.setdefault(row["element"], []).append(row["fraction"])
    assert set(by_element) == set(ir.elements_held(planet))
    for symbol, fractions in by_element.items():
        assert math.fsum(fractions) == pytest.approx(1.0, abs=1e-9), symbol
        assert all(0.0 <= fraction <= 1.0 for fraction in fractions)
    names = {row["isotope"] for row in rows if row["element"] == "H"}
    assert names == {"H-1", "H-2"}                                                 # deuterium exists on a body with water


def test_every_isotope_of_the_table_with_an_abundance_can_be_drawn_for_the_elements_held(sun):
    planet = _planet(sun, 1.0, seed=2)
    recorded = {row["isotope"] for row in planet.isotopes}
    for symbol in ir.elements_held(planet):
        if symbol in ir.RADIOACTIVE_PARENTS or symbol == "Ar":
            continue
        assert {iso.name for iso in elements.isotopes_of(symbol) if iso.abundance} <= recorded, symbol


def test_an_earth_like_planets_deuterium_nitrogen_and_argon_are_of_the_order_of_the_real_ones(sun):
    rows = {row["isotope"]: row for row in _planet(sun, 1.0, seed=1).isotopes}
    d_h = rows["H-2"]["fraction"] / rows["H-1"]["fraction"]
    assert 0.8e-4 < d_h < 4.0e-4                                                    # the ocean 1.56e-4, scattered by the draw
    assert 1.0 / (rows["N-15"]["fraction"] / rows["N-14"]["fraction"]) == pytest.approx(274.0, rel=0.1)
    assert rows["Ar-40"]["origin"] == ir.RADIOGENIC and rows["Ar-40"]["fraction"] > 0.9 and rows["Ar-36"]["fraction"] / rows["Ar-38"]["fraction"] == pytest.approx(5.501, rel=1e-6)
    assert rows["C-13"]["fraction"] / rows["C-12"]["fraction"] == pytest.approx(1.0 / 89.0, rel=0.05)


def test_the_decaying_parents_fall_with_the_age_of_the_body(sun):
    young, old = _body(star=_star(0.0, 1.0)), _body(star=_star(0.0, 9.0))
    assert ir._decayed("K", 9.0)[0]["K-40"] < ir._decayed("K", 1.0)[0]["K-40"]
    assert ir._decayed("U", 9.0)[0]["U-235"] < ir._decayed("U", 1.0)[0]["U-235"]
    assert math.fsum(ir._decayed("U", 4.5)[0].values()) == pytest.approx(1.0)
    assert ir._decayed("U", 4.568)[0]["U-235"] == pytest.approx(0.00725, rel=0.2)    # today's 235U share of the natural mix, near the Earth's age
    initial = ir._decayed("U", 4.5)[1]
    assert set(initial) == {"U-235"} and 0.2 < initial["U-235"] < 0.4                 # 235U/238U at the Solar System's birth, about 0.3
    assert ir._decayed("K", 4.5)[1]["K-40"] > 1e-3
    assert young is not old


def test_a_gas_giant_with_no_interior_content_has_primordial_argon_only(sun):
    giant = _planet(sun, 300.0, seed=3, distance=5.2)
    argon = {row["isotope"]: row["fraction"] for row in (giant.isotopes or []) if row["element"] == "Ar"}
    if argon:
        assert argon["Ar-40"] == 0.0 and argon["Ar-36"] / argon["Ar-38"] == pytest.approx(5.501, rel=1e-6)


def test_the_same_body_gets_the_same_record_and_the_run_stream_does_not_move(sun):
    planet = _planet(sun, seed=3)
    before = [dict(row) for row in planet.isotopes]
    with draw.bound("isotope-record:outer-a"):
        ir.generate(planet)
        mark = draw.random()
    with draw.bound("isotope-record:outer-a"):
        mark_without = draw.random()
    assert mark == mark_without and before == planet.isotopes


def test_a_body_with_no_stored_composition_has_no_record():
    body = SimpleNamespace(star=None)
    ir.generate(body)
    assert body.isotopes is None


def test_the_sheet_lists_the_ratios_the_record_moves_and_counts_the_rest():
    rows = [{"isotope": "H-1", "fraction": 0.9998}, {"isotope": "H-2", "fraction": 0.0002}, {"isotope": "O-16", "fraction": 0.9976}, {"isotope": "Ar-40", "fraction": 0.99}, {"isotope": "Ar-36", "fraction": 0.01}]
    shown = {key: (text, ratio) for key, _label, text, ratio in ir.sheet_rows(rows)}
    assert shown["isotope_D/H"][1] == pytest.approx(0.0002 / 0.9998) and shown["isotope_40Ar/36Ar"][1] == pytest.approx(99.0)
    assert "isotope_13C/12C" not in shown and shown["isotope_count"][0] == "5 isotopes of 3 elements"
    assert ir.sheet_rows([]) == []


def test_a_pulsar_planet_with_no_orbit_distance_takes_the_snow_lines_ice_and_still_gets_a_record():
    for distance in (0.0, None):
        body = _body(distance=0.0, formation_distance_au=distance, water_mass_fraction=0.01)
        assert ir._ice_temperature_k(body) == ir.ICE_REFERENCE_K
        rows = {row["isotope"]: row for row in ir.record_of(body)}
        assert rows["H-2"]["fraction"] > 0.0
