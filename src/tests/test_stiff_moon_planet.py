"""
A planet whose moon history is stiff (GEN.1815): seed 9101 of `tests/test_aquifers` has a rocky planet of 8.9e24 kg and a moon of 5.0e23 kg born at 3.8 planet radii, which the three BDF runs of
`moon_orbits.STRATEGIES` could not integrate. The moons then stayed as born over a spin that had not felt them, and their tide put 2.5e7 W/m2 (4e8 times the Earth's heat flow) into the planet. Radau finishes the history.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import moon_orbits
from tests.fuzz_support import deterministic_entropy

IO_FLUX_W_M2 = 2.0
"""float: The surface heat flow of Io, the most tidally heated body there is (about 2 W/m2)."""


@pytest.fixture(scope="module")
def system():
    config = SystemConfig()
    config.MOONS = True
    with deterministic_entropy(9101):
        return StarSystem(system_config=config)


def _bodies(system):
    for planet in system.primary.orbiters + (system.secondary.orbiters if system.secondary else []):
        if planet.body_type == "a":
            continue
        yield planet
        yield from planet.moons


def test_the_last_integrator_of_the_ladder_is_the_implicit_one_that_finishes_a_stiff_history():
    assert moon_orbits.STRATEGIES[-1][0] == "Radau"


def test_no_body_of_the_system_is_heated_by_a_tide_more_than_io_is(system):
    flux = {body.name: (getattr(body, "tidal_flux_w_m2", None) or 0.0) for body in _bodies(system)}
    assert max(flux.values()) < 10.0 * IO_FLUX_W_M2, {name: value for name, value in flux.items() if value >= 10.0 * IO_FLUX_W_M2}


def test_the_planet_with_the_stiff_history_lost_its_close_moon_to_the_tide(system):
    planet = next(p for p in system.primary.orbiters if getattr(p, "name", "").endswith("Miorutor III"))
    assert not planet.moons and planet.fallen_moons, "the 5e23 kg moon born at 3.8 radii falls to its planet under the tide (1.4 Gyr)"
    assert abs(planet.fallen_moons[0].mass_kg - 5.0e23) < 1.0e22


def _forced_failure(monkeypatch, times_yr):
    """The seed 9101 system with every integrator failing and `inspiral_times_s` giving `times_yr` (years) for the moons of each planet."""
    def failing(*args, **kwargs):
        raise moon_orbits.IntegrationFailed("forced")
    monkeypatch.setattr(moon_orbits, "evolve", failing)
    monkeypatch.setattr(moon_orbits, "inspiral_times_s", lambda planet, star, dissipation, moons, *args: [times_yr(m) * 3.15576e7 for m in moons])
    config = SystemConfig()
    config.MOONS = True
    with deterministic_entropy(9101):
        forced = StarSystem(system_config=config)
    return next(p for p in forced.primary.orbiters if getattr(p, "name", "").endswith("Miorutor III"))


def test_a_history_no_integrator_finishes_loses_the_moon_the_tide_brings_in_within_the_age(monkeypatch):
    """GEN.1816: the moon whose inspiral time at birth is under the age is lost to its planet (recorded as fallen at that time); the planet is marked."""
    planet = _forced_failure(monkeypatch, lambda moon: 1.0e9)
    assert planet.build.moon_history_unfinished
    assert not planet.moons
    assert planet.fallen_moons and abs(planet.fallen_moons[0].mass_kg - 5.0e23) < 1.0e22
    assert planet.fallen_moons[0].time_yr == pytest.approx(1.0e9)


def test_a_history_no_integrator_finishes_keeps_a_moon_the_tide_does_not_bring_in_as_born(monkeypatch):
    planet = _forced_failure(monkeypatch, lambda moon: float("inf"))
    assert planet.build.moon_history_unfinished
    assert len(planet.moons) == 1 and not planet.fallen_moons
    assert planet.moons[0].distance == pytest.approx(planet.moons[0]._birth_distance_au)


def test_a_planet_whose_history_finished_is_not_marked(system):
    assert not any(body.build.moon_history_unfinished for body in _bodies(system) if getattr(body, "build", None) is not None)

