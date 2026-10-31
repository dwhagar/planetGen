"""
Mantle redox and atmosphere species (GEN.85, `planetgen.physics.atmosphere`).
"""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere, gases, saturation

from tests.air_support import air
from tests.fuzz_support import deterministic_entropy


def _bodies(seeds):
    for seed in seeds:
        with deterministic_entropy(seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets:
            if hasattr(planet, "moons"):
                yield planet
                yield from planet.moons


def test_redox_classes():
    assert atmosphere.redox_class(-0.5) == "reduced"
    assert atmosphere.redox_class(1.0) == "intermediate"
    assert atmosphere.redox_class(3.5) == "oxidized"


def test_vapour_pressure_is_the_saturation_curve():
    # CO2 at 0.6 kPa freezes near 148 K (Mars's polar frost point); water boils at 373 K; N2 at 37 K is Pluto's air
    # (about 1 Pa; the old two-point table gave 10 Pa)
    assert atmosphere.vapour_pressure_kpa("co2", 148.0) == pytest.approx(0.6, rel=0.25)
    assert atmosphere.vapour_pressure_kpa("h2o", 373.124) == pytest.approx(101.325, rel=0.003)
    assert atmosphere.vapour_pressure_kpa("n2", 37.0) == pytest.approx(0.0012, rel=0.1)
    assert atmosphere.vapour_pressure_kpa("n2", 300.0) == math.inf
    for gas in atmosphere.SPECIES:
        assert atmosphere.vapour_pressure_kpa(gas, 100.0) == saturation.pressure_pa(gases.stored(gas).formula, 100.0) / 1000.0


def test_the_two_point_condensation_table_is_gone_and_hydrogen_sulphide_has_its_own_curve():
    assert not hasattr(atmosphere, "CONDENSATION_K")
    # H2S boils near 213 K at 100 kPa, SO2 near 263 K (the old table gave H2S the SO2 row)
    assert atmosphere.vapour_pressure_kpa("h2s", 213.0) == pytest.approx(100.0, rel=0.03)
    assert atmosphere.vapour_pressure_kpa("so2", 263.0) == pytest.approx(100.0, rel=0.03)
    assert atmosphere.vapour_pressure_kpa("h2s", 150.0) > 5 * atmosphere.vapour_pressure_kpa("so2", 150.0)


def test_describe_names_majors_and_traces():
    text = atmosphere.describe({"N2": 78000.0, "O2": 21000.0, "Ar": 930.0, "CO2": 40.0})
    assert text == "a mix of nitrogen and oxygen, with traces of argon and carbon dioxide"
    assert atmosphere.describe({}) == "None"
    giant = atmosphere.describe({"H2": 86.0e3, "He": 13.7e3, "CH4": 0.3e3})
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
        # GEN.424: every gas is stored, so the rows add up to the surface pressure.
        air = atmosphere.gas_partials_pa(body)
        assert sum(air.values()) == pytest.approx(body.atmospheric_pressure, rel=1e-9)
        assert all(gases.has(formula) and value > 0.0 for formula, value in air.items())
    assert count > 0


def test_every_generated_body_says_its_stored_mixture_and_never_a_class_phrase():
    """GEN.415: the atmosphere text is `describe` of the stored partial pressures, "None" for an airless body, and the placeholder a
    class with air starts with never survives generation."""
    seen = set()
    for body in _bodies(range(12)):
        assert body.atmosphere != atmosphere.PENDING, body.name
        assert body.atmosphere == atmosphere.describe(atmosphere.gas_partials_pa(body)), body.name
        seen.add(body.atmosphere == "None")
    assert seen == {True, False}


