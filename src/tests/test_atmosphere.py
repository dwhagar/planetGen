"""
Mantle redox and atmosphere species (GEN.85, `planetgen.physics.atmosphere`).
"""

import math
import statistics

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere, constants
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy


def _bodies(seeds):
    for seed in seeds:
        with deterministic_entropy(seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets:
            if hasattr(planet, "moons"):
                yield planet
                yield from planet.moons


def test_redox_rises_with_mass():
    draw.set_run_seed(85)
    moon = statistics.fmean(atmosphere.draw_mantle_delta_iw(0.0123 * constants.EARTH_MASS_TO_KG) for _ in range(2000))
    earth = statistics.fmean(atmosphere.draw_mantle_delta_iw(constants.EARTH_MASS_TO_KG) for _ in range(2000))
    assert moon == pytest.approx(-0.32, abs=0.1)
    assert earth == pytest.approx(atmosphere.REDOX_EARTH_DELTA_IW, abs=0.1)
    assert atmosphere.redox_class(-0.5) == "reduced"
    assert atmosphere.redox_class(1.0) == "intermediate"
    assert atmosphere.redox_class(3.5) == "oxidized"


def test_vapour_pressure_matches_the_design_table():
    # atmospheres-retention-and-classes.md 3.3: CO2 at 0.6 kPa freezes near
    # 148 K (Mars's polar frost point); water boils at 373 K.
    assert atmosphere.vapour_pressure_kpa("co2", 148.0) == pytest.approx(0.6, rel=0.25)
    assert atmosphere.vapour_pressure_kpa("h2o", 373.0) == pytest.approx(100.0)
    assert atmosphere.vapour_pressure_kpa("n2", 37.0) == pytest.approx(0.01, rel=0.2)


def test_class_mixes_sum_to_one():
    for planet_class, (template, balance) in atmosphere.CLASS_MIXES.items():
        assert set(template) <= set(atmosphere.SPECIES), planet_class
        total = sum(template.values())
        assert total <= 1.0 + 1e-3, planet_class
        if balance is None:
            assert total == pytest.approx(1.0, abs=0.01), planet_class


def test_a_reduced_mantle_outgasses_reduced_species_under_anoxic_air():
    draw.set_run_seed(1)
    oxidized, _, _ = atmosphere.draw_mixing_ratios("A", "oxidized")
    draw.set_run_seed(1)
    reduced, _, _ = atmosphere.draw_mixing_ratios("A", "reduced")
    assert oxidized["h2s"] == 0.0 and reduced["h2s"] > 0.0
    assert reduced["co"] > oxidized["co"]
    assert reduced["co2"] < oxidized["co2"]
    # Oxygenated air (life) keeps its mix whatever the mantle.
    draw.set_run_seed(2)
    earth_ox, _, _ = atmosphere.draw_mixing_ratios("M", "oxidized")
    draw.set_run_seed(2)
    earth_red, _, _ = atmosphere.draw_mixing_ratios("M", "reduced")
    assert earth_ox == earth_red


def test_describe_names_majors_and_traces():
    text = atmosphere.describe({"n2": 78.0, "o2": 21.0, "ar": 0.93, "co2": 0.04, "h2o": 0.0,
                                "co": 0.0, "h2": 0.0, "ch4": 0.0, "h2s": 0.0, "so2": 0.0}, 0.0, None)
    assert text == "a mix of nitrogen and oxygen, with traces of argon and carbon dioxide"
    assert atmosphere.describe({gas: 0.0 for gas in atmosphere.SPECIES}, 0.0, None) == "None"
    giant = atmosphere.describe({"h2": 86.0, "ch4": 0.3}, 13.7, "helium")
    assert giant.startswith("a mix of hydrogen and helium") and "reducing" not in giant


def test_generated_bodies_store_consistent_air():
    count = 0
    for body in _bodies(range(12)):
        partials = atmosphere.partial_pressures_kpa(body)
        if body.body_type == "g":
            assert body.mantle_redox is None
        else:
            assert body.mantle_redox == atmosphere.redox_class(body.mantle_delta_iw)
        if body.atmosphere == "None":
            assert all(value == 0.0 for value in partials.values())
            continue
        count += 1
        total_kpa = body.atmospheric_pressure / 1000.0
        assert sum(partials.values()) <= total_kpa * (1 + 1e-9)
        if body.body_type == "t" and atmosphere.CLASS_MIXES.get(body.planet_class, (None, None))[1] is None:
            assert sum(partials.values()) == pytest.approx(total_kpa)
            for gas, value in partials.items():
                if value > 0.0:
                    assert value <= atmosphere.vapour_pressure_kpa(gas, body.surface_temperature) * (1 + 1e-9)
    assert count > 0


def test_hydrogen_inventory_counts_every_hydrogen_carrier():
    body = type("Body", (), {f"p_{gas}_kpa": 0.0 for gas in atmosphere.SPECIES})()
    body.p_h2o_kpa, body.p_ch4_kpa, body.p_h2_kpa, body.p_h2s_kpa = 1.0, 0.5, 0.25, 0.1
    assert atmosphere.hydrogen_kpa(body) == pytest.approx(1.0 + 1.0 + 0.25 + 0.1)
    assert math.isclose(atmosphere.hydrogen_kpa(object()), 0.0)
