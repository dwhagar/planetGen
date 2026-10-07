# tests/test_prevalence.py

"""
Prevalence (GEN.52, TEST.75): a sector or galaxy run's percentage
deviation from a feature's normal chance. These generate systems with no
database and count how many come out with each feature, so the shares
are checked with a tolerance a few standard errors wide.
"""

import random

import pytest

from planetgen import tuning
from planetgen.cli import generate as generate_cli
from planetgen.generation import prevalence, run_system
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.population.model import has_civilization

SYSTEMS = 2000


def _features(system):
    planets, belts, _moons = system.count_objects()
    return {
        "habitable_world": system.count_habitable()[0] > 0,
        "asteroid_belt": belts > 0,
        "large_star": (system.primary_star.initial_mass_sol or 0.0) >= tuning.LARGE_STAR_MIN_MASS_SOL,
        "planets": planets + belts > 0,
        "comets": bool(system.comets or system.secondary_comets),
        "binary_system": system.binary_type is not None,
    }


def _shares(settings, seed, count=SYSTEMS):
    random.seed(seed)
    totals = {}
    for _ in range(count):
        config = SystemConfig()
        config.PREVALENCE = dict(settings)
        for feature, present in _features(StarSystem(system_config=config)).items():
            totals[feature] = totals.get(feature, 0) + present
    return {feature: n / count for feature, n in totals.items()}


def _tolerance(share, count=SYSTEMS):
    return 4 * (max(share * (1 - share), 0.01) / count) ** 0.5


def test_measured_base_shares_still_hold():
    """tuning.PREVALENCE_BASE_SHARES are what generation gives today; a
    change that moves one needs it remeasured."""
    shares = _shares({}, seed=11)
    for feature, base in tuning.PREVALENCE_BASE_SHARES.items():
        assert shares[feature] == pytest.approx(base, abs=_tolerance(base)), feature


@pytest.mark.parametrize("feature, percent", [
    ("habitable_world", 50), ("habitable_world", -50), ("asteroid_belt", -60), ("asteroid_belt", 30),
    ("large_star", 100), ("large_star", -100), ("planets", -40),
])
def test_a_measured_feature_moves_by_the_percentage(feature, percent):
    expected = min(1.0, tuning.PREVALENCE_BASE_SHARES[feature] * (1 + percent / 100))
    share = _shares({feature: percent}, seed=hash((feature, percent)) % 1000)[feature]
    assert share == pytest.approx(expected, abs=_tolerance(expected)), (feature, percent, share)


@pytest.mark.parametrize("feature, base, percent", [
    ("comets", tuning.SYSTEM_COMET_CHANCE, 100), ("comets", tuning.SYSTEM_COMET_CHANCE, -100),
])
def test_a_drawn_feature_moves_by_the_percentage(feature, base, percent):
    expected = min(1.0, base * (1 + percent / 100))
    share = _shares({feature: percent}, seed=7)[feature]
    assert share == pytest.approx(expected, abs=_tolerance(expected))


def test_no_binaries_at_minus_100_and_more_at_plus_50():
    assert _shares({"binary_system": -100}, seed=3, count=500)["binary_system"] == 0
    assert _shares({"binary_system": 50}, seed=3, count=1000)["binary_system"] > _shares({}, seed=3, count=1000)["binary_system"] + 0.05


def test_moons_and_max_planets_scale_their_draws():
    random.seed(5)
    for percent, check in ((-100, lambda planets: all(not p.moons for p in planets)),):
        for _ in range(200):
            config = SystemConfig()
            config.PREVALENCE = {"moons": percent}
            system = StarSystem(system_config=config)
            assert check([p for p in system.planets + system.secondary_planets if p.body_type != "a"])


def test_scaled_chance_is_capped_and_floored():
    config = SystemConfig()
    config.PREVALENCE = {"comets": 500}
    assert prevalence.scaled_chance(config, "comets", 0.3) == 1.0
    config.PREVALENCE = {"comets": -100}
    assert prevalence.scaled_chance(config, "comets", 0.3) == 0.0
    assert prevalence.scaled_chance(SystemConfig(), "comets", 0.3) == 0.3


def test_resolve_keeps_forced_flags_and_the_forcing_rules():
    random.seed(1)
    for _ in range(300):
        config = SystemConfig()
        config.COMETS = True
        config.PREVALENCE = {"habitable_world": 1000, "asteroid_belt": 1000, "planets": -100}
        prevalence.resolve(config)
        assert config.HABITABLE_WORLD is True and config.ASTEROID_BELT is True
        assert config.LARGE_STAR is True  # both need a large star
        assert config.PLANETS is None  # and planets, so -planets gives way
    config = SystemConfig()
    config.LARGE_STAR = False
    config.HABITABLE_WORLD = True
    config.PREVALENCE = {"asteroid_belt": 1000}
    prevalence.resolve(config)
    assert config.ASTEROID_BELT is None and config.LARGE_STAR is False


def test_civilizations_scale_with_intelligent_life_prevalence():
    class _Timeline:
        life_stage = "technological_civilization"
    rng = random.Random(2)
    assert not any(has_civilization(_Timeline(), None, rng, -100) for _ in range(2000))
    hits = sum(has_civilization(_Timeline(), None, rng, 99_900) for _ in range(2000))
    assert hits == pytest.approx(2000 * tuning.CIVILIZATION_CHANCE * 1000, rel=0.15)


def test_large_star_minus_forbids_one():
    random.seed(9)
    for _ in range(300):
        config = SystemConfig()
        config.LARGE_STAR = False
        assert StarSystem(system_config=config).primary_star.initial_mass_sol < tuning.LARGE_STAR_MIN_MASS_SOL


def test_cli_prevalence_reaches_every_system_config():
    args = generate_cli.build_parser()[0].parse_args(
        ["sector", "--prevalence", "comets=+50", "--prevalence", "Habitable-World=-20%", "--prevalence", "comets=10"])
    args.system_file = args.num_orbits = args.name = None  # as validate_sector_args sets them
    config = run_system.build_system_config(args)
    assert config.PREVALENCE == {"comets": 10.0, "habitable_world": -20.0}


@pytest.mark.parametrize("bad", ["comets", "spaceships=10", "comets=-101", "comets=lots"])
def test_cli_rejects_a_bad_prevalence(bad):
    with pytest.raises(SystemExit):
        generate_cli.build_parser()[0].parse_args(["sector", "--prevalence", bad])


def test_prevalence_round_trips_through_the_stored_config(mysql_config):
    from planetgen.db import store

    conn = store.get_connection(mysql_config)
    try:
        config = SystemConfig()
        config.PREVALENCE = {"comets": 50.0, "intelligent_life": -100.0}
        config_id = store.insert_system_config(conn, config)
        conn.commit()
        assert store.load_system_config(conn, config_id).PREVALENCE == {"comets": 50.0, "intelligent_life": -100.0}
        assert store.load_system_config(conn, store.insert_system_config(conn, SystemConfig())).PREVALENCE == {}
    finally:
        conn.close()
