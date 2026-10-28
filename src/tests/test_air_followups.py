"""
What the air fixed point (GEN.411, GEN.422) left to wire: the ozone column reads the oxygen through the mixture accessor (GEN.667) and
every giant carries the photosphere pressure the hot-Jupiter rows gate on (GEN.673).
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.planet import Planet
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere, mixture, mixture_subject
from tests.fuzz_support import deterministic_entropy


def _giant(distance_au, seed):
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    with deterministic_entropy(seed):
        star = Star(config)
        return Planet(config, star, star.habitable_zone, distance_au, planet_class="J", moon_count=0)


def test_the_oxygen_the_ozone_reads_is_the_mixture_accessors():
    checked = 0
    for seed in range(12):
        with deterministic_entropy(seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets:
            if getattr(planet, "body_type", "t") == "g" or getattr(planet, "ozone_column_du", None) is None:
                continue
            via_accessor = mixture.partial_pressures_pa(planet).get("O2", 0.0)
            assert via_accessor == pytest.approx(atmosphere.partial_kpa(planet, "o2") * 1000.0, rel=1e-9)
            checked += 1
    assert checked


@pytest.mark.parametrize("distance_au", [0.04, 0.1, 1.0, 5.0])
def test_every_giant_has_the_photosphere_pressure_the_hot_rows_gate_on(distance_au):
    giant = _giant(distance_au, 1)
    subject = mixture_subject.subject_of_body(giant, ts_k=giant.surface_temperature, pressure_bar=1.0, ocean_present=False,
                                              long_lived_geology=False, water_lost=False, star_distance_au=distance_au)
    assert subject["photosphere_bar"] == mixture_subject.PHOTOSPHERE_BAR


def test_a_hot_jupiter_takes_the_hot_rows_and_a_cold_giant_the_cold_ones():
    hot = atmosphere.gas_partials_pa(_giant(0.04, 0))
    cold = atmosphere.gas_partials_pa(_giant(5.0, 0))
    assert hot["CO"] > hot.get("CH4", 0.0)      # carbon sits in CO where it is hot, in CH4 where it is cold
    assert cold["CH4"] > cold.get("CO", 0.0)
