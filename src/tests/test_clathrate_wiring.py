"""GEN.762: the generator reads the clathrate and brine limits (`physics/clathrate_brine.py`) in the air (`surface_air`, whose hydrate cages are tested in test_surface_air.py) and the hydrosphere."""

import ast
import pathlib
import types

import pytest

from planetgen.physics import clathrate_brine as cb, constants, hydrosphere as hy
from planetgen.util import draw


def _generate_hydrosphere(body, internal_flux_w_m2):
    """`generate_hydrosphere` after the surface water is drawn, as `planets.calculate_atmospheric_conditions` does."""
    return hy.generate_hydrosphere(body, internal_flux_w_m2, hy.draw_surface_water(body) if body.body_type != "g" else 0.0)


EARTH = constants.EARTH_MASS_TO_KG
PHYSICS = pathlib.Path(__file__).resolve().parents[1] / "planetgen" / "physics"


@pytest.mark.parametrize("module", ["surface_air.py", "hydrosphere.py"])
def test_the_generator_calls_the_clathrate_module(module):
    """The wiring guard for this module: it is called, not only imported (GEN.754 will make the check general)."""
    tree = ast.parse((PHYSICS / module).read_text(encoding="utf-8"))
    calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
             and isinstance(node.func.value, ast.Name) and node.func.value.id == "clathrate_brine"]
    assert calls, f"{module} never calls clathrate_brine"


def test_a_guest_above_its_curve_is_held_down_to_the_curve_when_the_cages_are_plentiful():
    floor = cb.dissociation_pressure_pa("CH4", 100.0)
    held = cb.hydrate_hold_pa({"ch4": 10.0 * floor, "n2": 1.0, "ar": 5e4}, 100.0, 1.35, 1e6)
    assert held["ch4"] == pytest.approx(9.0 * floor)
    assert "ar" not in held  # argon forms no hydrate here


def test_nothing_is_held_below_the_curve_above_the_ice_point_or_with_no_gravity():
    floor = cb.dissociation_pressure_pa("CH4", 100.0)
    assert cb.hydrate_hold_pa({"ch4": 0.5 * floor}, 100.0, 1.35, 1e6) == {}
    assert cb.hydrate_hold_pa({"co2": 5e5}, 275.0, 1.35, 1e6) == {}  # past CO2's quadruple point
    assert cb.hydrate_hold_pa({"ch4": 10.0 * floor}, 100.0, 0.0, 1e6) == {}


def test_a_thin_skin_holds_less_and_the_guests_share_the_cages():
    pressures = {"ch4": 1e4, "co2": 1e4, "ar": 1e4}
    plenty = cb.hydrate_hold_pa(pressures, 100.0, 1.35, 1e6)
    thin = cb.hydrate_hold_pa(pressures, 100.0, 1.35, 1.0)
    assert sum(thin.values()) < sum(plenty.values())
    shares = {gas: thin[gas] / plenty[gas] for gas in thin}
    assert shares["ch4"] == pytest.approx(shares["co2"])  # the same share of every guest's excess
    cages_mol = cb.ICE_DENSITY_KG_M3 * 1.0 / (6 * 0.018015)
    mean_molar_mass = sum(p * m for p, m in ((1e4, 0.016043), (1e4, 0.04401), (1e4, 0.039948))) / 3e4
    assert sum(thin.values()) / (mean_molar_mass * 1.35) == pytest.approx(cages_mol, rel=0.01)


def _dry_cold_body(temperature=210.0):
    return types.SimpleNamespace(body_type="t", planet_class="C", mass=0.3 * EARTH, radius=4000.0, water_mass_fraction=1e-5, gravity=0.5,
                                 surface_temperature=temperature, atmospheric_pressure=600.0, is_moon=False, mantle_redox="intermediate",
                                 p_co2_kpa=0.6, p_so2_kpa=0.0, p_h2_kpa=0.0)


def test_a_cold_body_with_little_water_keeps_a_brine_ocean_below_the_ice_point():
    classes, activities = set(), set()
    for seed in range(12):
        draw.set_run_seed(seed)
        body = _dry_cold_body()
        _generate_hydrosphere(body, 0.03)
        if body.ocean_depth_km:
            classes.add(body.ocean_class)
            activities.add(body.water_activity)
            assert body.hydrosphere in ("ice-covered ocean", "surface ocean")
    assert classes == {"chloride brine"}
    assert activities & {cb.water_activity("Ca(ClO4)2"), cb.water_activity("Mg(ClO4)2")}  # the drawn salt's own activity


def test_the_brine_lid_is_thinner_than_the_pure_water_lid():
    pure = hy.pressure_melted_shell_km(0.03, 150.0, 3.7)
    brine = hy.pressure_melted_shell_km(0.03, 150.0, 3.7, melting_k=cb.eutectic_k("Mg(ClO4)2"))
    assert brine < pure


def test_a_body_with_more_water_keeps_pure_water_melting():
    draw.set_run_seed(1)
    body = _dry_cold_body(200.0)
    body.water_mass_fraction = 1e-3
    _generate_hydrosphere(body, 0.03)
    assert body.ocean_class in (None, "neutral", "ice-sealed", "acid sulfate", "soda")
