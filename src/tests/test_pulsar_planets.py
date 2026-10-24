"""
Planets of a millisecond pulsar (GEN.317, docs/design/systems/exotic-environments-planets-and-compact-binaries.md section 2.1):
the origin of the ones that hold any, the disc origin's planets, and the pulsar wind they sit in.
"""

import math

import pytest

from planetgen import tuning
from planetgen.generation import pulsar_planets
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import NeutronStar
from planetgen.generation.system import StarSystem
from planetgen.physics import activity, constants, magnetism
from planetgen.util import draw


@pytest.fixture
def recycled(monkeypatch):
    """A maker of old neutron stars that are always recycled millisecond pulsars."""
    monkeypatch.setattr(tuning, "PULSAR_RECYCLED_CHANCE", 1.0)

    def make(seed):
        with draw.bound(draw.Stream(f"pulsar-planets:{seed}")):
            return NeutronStar(SystemConfig(), age_gy=8.0)

    return make


def _planets(pulsar, seed, **config):
    """The planets of `pulsar` for a config with `config` set, drawn on the stream of `seed`."""
    cfg = SystemConfig()
    for name, value in config.items():
        setattr(cfg, name, value)
    with draw.bound(draw.Stream(f"pulsar-planets-draw:{seed}")):
        return pulsar_planets.generate(cfg, pulsar, pulsar.habitable_zone, 1.0e5)


def test_the_origin_weights_are_the_designs_quarter_quarter_half():
    assert tuning.PULSAR_PLANET_ORIGIN_WEIGHTS == {"disc": 25.0, "captured-giant": 25.0, "ablated-companion": 50.0}
    assert tuning.PULSAR_PLANET_BUILT_ORIGINS == ("disc",)


def test_a_pulsars_origin_follows_the_weights():
    config = SystemConfig()
    with draw.bound(draw.Stream("pulsar-origins")):
        origins = [pulsar_planets.draw_origin(config) for _ in range(2000)]
    shares = {name: origins.count(name) / len(origins) for name in tuning.PULSAR_PLANET_ORIGIN_WEIGHTS}
    assert shares["disc"] == pytest.approx(0.25, abs=0.04)
    assert shares["captured-giant"] == pytest.approx(0.25, abs=0.04)
    assert shares["ablated-companion"] == pytest.approx(0.50, abs=0.04)


def test_a_run_that_asks_for_planets_gets_the_disc():
    with draw.bound(draw.Stream("pulsar-forced")):
        assert {pulsar_planets.draw_origin(_config(PLANETS=True)) for _ in range(40)} == {"disc"}


def _config(**values):
    config = SystemConfig()
    for name, value in values.items():
        setattr(config, name, value)
    return config


def test_only_the_disc_origin_gives_planets(recycled):
    pulsar = recycled(1)
    counts = [len(_planets(pulsar, seed)) for seed in range(240)]
    with_planets = sum(1 for count in counts if count)
    assert with_planets / len(counts) == pytest.approx(0.25, abs=0.07)
    assert max(counts) <= tuning.PULSAR_DISC_PLANET_COUNT[1]


def test_the_discs_planets_are_airless_carbon_rich_and_lowest_tier(recycled):
    pulsar = recycled(2)
    planets = [p for seed in range(30) for p in _planets(pulsar, seed, PLANETS=True)]
    assert planets
    for planet in planets:
        assert planet.planet_class == tuning.PULSAR_DISC_CLASS
        assert planet.atmosphere == "None" and planet.atmospheric_pressure == 0.0
        assert "anthracite" in planet.composition
        assert planet.moons == []
        assert planet.tier_radiation == 3  # the lowest tier, whatever the world


def test_the_discs_planets_keep_to_its_count_distance_and_mass(recycled):
    pulsar = recycled(3)
    near, far = tuning.PULSAR_DISC_DISTANCE_AU_RANGE
    low, high = pulsar_planets.disc_mass_range_earth()
    assert low == tuning.PULSAR_DISC_MASS_EARTH_RANGE[0]
    assert low < high <= tuning.PULSAR_DISC_MASS_EARTH_RANGE[1]
    for seed in range(40):
        planets = _planets(pulsar, seed, PLANETS=True)
        assert tuning.PULSAR_DISC_PLANET_COUNT[0] <= len(planets) <= tuning.PULSAR_DISC_PLANET_COUNT[1]
        distances = [p.distance for p in planets]
        assert distances == sorted(distances)
        assert all(near * 0.999 <= d <= far * 1.001 for d in distances)
        # The mass follows the radius at a density in the class's range, so it lands a little either side of the draw.
        assert all(0.5 * low <= p.mass / constants.EARTH_MASS_TO_KG <= 1.3 * high for p in planets)


def test_the_discs_planets_sit_near_a_low_order_period_ratio(recycled):
    pulsar = recycled(4)
    found = 0
    for seed in range(60):
        planets = _planets(pulsar, seed, PLANETS=True, NUM_ORBITS=3)
        for inner, outer in zip(planets, planets[1:]):
            ratio = outer.period / inner.period
            assert any(abs(ratio / wanted - 1.0) < 0.06 for wanted in tuning.PULSAR_DISC_PERIOD_RATIOS)
            found += 1
    assert found >= 60


def test_num_orbits_and_planets_false_are_obeyed(recycled):
    pulsar = recycled(5)
    assert len(_planets(pulsar, 1, PLANETS=True, NUM_ORBITS=2)) == 2
    assert _planets(pulsar, 1, PLANETS=True, NUM_ORBITS=0) == []
    assert _planets(pulsar, 1, PLANETS=False) == []


def test_a_system_around_a_millisecond_pulsar_takes_the_discs_planets(recycled):
    pulsar = recycled(6)
    config = _config(PLANETS=True)
    with draw.bound(draw.Stream("pulsar-system")):
        system = StarSystem(system_config=config, compact_remnant=pulsar)
    assert system.planets
    assert all(p.planet_class == tuning.PULSAR_DISC_CLASS for p in system.planets)


def test_a_young_pulsar_keeps_the_generic_planets(monkeypatch):
    monkeypatch.setattr(tuning, "PULSAR_RECYCLED_CHANCE", 0.0)
    with draw.bound(draw.Stream("pulsar-young")):
        young = NeutronStar(SystemConfig(), age_gy=1.0e-3)
    assert not pulsar_planets.applies(young)


def test_the_pulsar_wind_pressure_matches_the_designs_table():
    # section 2.1: PSR B1257+12 at 0.19 AU and a typical millisecond pulsar at 1 AU
    assert magnetism.pulsar_wind_pressure_pa(1.9e27, 0.19) == pytest.approx(6.2e-4, rel=0.03)
    assert magnetism.pulsar_wind_pressure_pa(1.5e27, 1.0) == pytest.approx(1.8e-5, rel=0.03)
    # It falls with the square of the distance.
    assert (magnetism.pulsar_wind_pressure_pa(1.0e27, 0.2) / magnetism.pulsar_wind_pressure_pa(1.0e27, 0.4)
            == pytest.approx(4.0))


def test_a_pulsing_stars_wind_is_its_pulsar_wind_and_a_silent_ones_is_not(recycled):
    pulsar = recycled(7)
    assert pulsar.spin_down_power_w == pytest.approx(activity.neutron_star_spin_down_w(pulsar))
    assert magnetism.wind_pressure_pa(pulsar, 0.3) == pytest.approx(
        magnetism.pulsar_wind_pressure_pa(pulsar.spin_down_power_w, 0.3))
    pulsar.pulsar_type = "non-pulsing"
    assert magnetism.wind_pressure_pa(pulsar, 0.3) == pytest.approx(magnetism.EARTH_WIND_PA / 0.3 ** 2)


def test_the_xuv_still_reads_the_larger_of_the_thermal_and_the_spin_down_output(recycled):
    pulsar = recycled(8)
    expected = max(pulsar.luminosity, activity.DIPOLE_EDOT_TO_LX * pulsar.spin_down_power_w)
    assert activity.neutron_star_xuv_w(pulsar) == pytest.approx(expected)
    assert not math.isnan(pulsar.spin_down_power_w) and pulsar.spin_down_power_w > 0.0
