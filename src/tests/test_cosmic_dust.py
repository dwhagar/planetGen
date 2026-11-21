# tests/test_cosmic_dust.py

"""GEN.1542 slice 1: the cosmic dust a rocky planet takes in decides the carbon V8 of its sulfur cloud sees (`generation/cosmic_dust.py`)."""

import pytest

from planetgen.generation import cosmic_dust
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import dust_supply, influx, sulfur_cloud_set
from planetgen.util import draw

ACID_CLOUD_SEED = 137    # a hot rocky planet of this default system holds an acid cloud (V1) with the hypotheses on or off; tried first (TEST.809: seed 32 held one until the sulfur sinks left almost no air with sulfur dioxide)

SEEDS_TRIED = 200
"""int: How many seeds the fixture tries, after `ACID_CLOUD_SEED`, for a default system with a planet that holds an acid cloud. The sulfur dioxide an acid cloud needs is in about 1 system in 20, so the first seed
that has it ends the search; the test fails, naming this count, when none of them does (it never passes without the cloud)."""


def _acid_planets(system):
    return [planet for planet in system.all_planets if planet.body_type != "a" and getattr(planet, "cloud_chemistry", None) and "V1" in planet.cloud_chemistry]


def _default_system(seed):
    draw.set_run_seed(seed)
    config = SystemConfig()
    config.PLANETS = True
    return StarSystem(system_config=config)


@pytest.fixture(scope="module")
def system():
    for seed in (ACID_CLOUD_SEED, *range(SEEDS_TRIED)):
        found = _default_system(seed)
        if _acid_planets(found):
            return found
    pytest.fail(f"none of the {SEEDS_TRIED + 1} default systems tried (seed {ACID_CLOUD_SEED} first, then 0 to {SEEDS_TRIED - 1}) holds a planet with an acid cloud (V1)")


def _comets(system):
    return [(system.star, system.comets), (getattr(system, "secondary_star", None), system.secondary_comets)]


def _with_cloud(system):
    return _acid_planets(system)


def test_the_pass_does_nothing_while_the_hypotheses_are_off(system, monkeypatch):
    assert not sulfur_cloud_set.INCLUDE_HYPOTHESES
    monkeypatch.setattr(sulfur_cloud_set, "fill", lambda *a, **k: pytest.fail("the cloud chemistry was made again with the hypotheses off"))
    assert cosmic_dust.apply(system.all_planets, _comets(system)) == 0


def test_the_pass_gives_each_cloudy_planet_the_carbon_its_own_dust_supplies(system, monkeypatch):
    monkeypatch.setattr(sulfur_cloud_set, "INCLUDE_HYPOTHESES", True)
    calls = []
    monkeypatch.setattr(sulfur_cloud_set, "fill", lambda body, values, carbon_dust=False: calls.append((body, carbon_dust)))
    cloudy = _with_cloud(system)
    assert cloudy, "the system holds a planet with an acid cloud"
    made = cosmic_dust.apply(system.all_planets, _comets(system))
    assert made == len(calls) >= len(cloudy)
    for body, carbon_dust in calls:
        supply = cosmic_dust.influx_of(body, system.all_planets, _comets(system))
        assert carbon_dust is dust_supply.v8_supplied(supply, body.radius)


def test_a_planet_reads_the_belts_and_comets_of_its_own_star_and_the_giants_in_the_way(system):
    planet = _with_cloud(system)[0]
    reservoirs = cosmic_dust.reservoirs_of(planet, system.all_planets, _comets(system))
    belts = [r for r in reservoirs if r.kind == "belt"]
    assert len(belts) == sum(1 for b in system.all_planets if b.body_type == "a" and b.mass_kg and cosmic_dust.late_accretion._same_star(b, planet))
    assert len([r for r in reservoirs if r.kind == "comet"]) == sum(1 for c in system.comets if c.mass_kg and c.perihelion_distance_au)
    assert cosmic_dust.influx_of(planet, system.all_planets, _comets(system)) > 0.0
    assert cosmic_dust.influx_of(planet, [], []) == 0.0, "a system with nothing to send in gives nothing"


def test_a_cloudless_system_is_left_alone(monkeypatch):
    monkeypatch.setattr(sulfur_cloud_set, "INCLUDE_HYPOTHESES", True)
    assert cosmic_dust.apply([], []) == 0
