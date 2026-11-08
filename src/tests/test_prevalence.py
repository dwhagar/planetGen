# tests/test_prevalence.py

"""
Prevalence (GEN.52, TEST.75): a sector or galaxy run's percentage
deviation from a feature's normal chance. These generate systems with no
database and count how many come out with each feature, so the shares
are checked with a tolerance a few standard errors wide.
"""

import functools
import random
import zlib

import pytest

from planetgen.util import draw
from planetgen import tuning
from planetgen.cli import generate as generate_cli
from planetgen.generation import prevalence, run_system
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.population.model import has_civilization

SYSTEMS = 140
"""int: Systems drawn to measure a share (PERF.236). A full system costs about 3 seconds on a lane machine (a `StarSystem` with its air, climate and spin
histories), so the 2,000 this was costs 6,000 seconds against the gate's 900 s limit per test; 140 cost about 440 s, half the limit as a margin for a slower
machine. The tolerance (`_tolerance`) is four standard errors of the count, which at 140 is 0.12 to 0.17 for a share of 0.3 to 0.9: wide, but the changes
these tests check move a share by 0.15 to 0.3, and each draw is seeded so it repeats."""
COMET_SYSTEMS = 100
"""int: Systems drawn for the comet and binary cases, whose systems cost more (a pair or comets make the system larger): 830 s for 140 on the slowest cell, so 100 (about 600 s)."""
FORCED_ON_SYSTEMS = 50
"""int: Systems drawn when a feature is forced on or barred in many of them (`habitable_world` +50 and -50): the forcing or the bar repeats a draw until the condition holds, which costs
about 8.6 seconds a system forced on and about 12 seconds a system barred (GEN.1275 judged every body), so 50 cost 430 to 600 s. The tolerance at 50 is 0.17 to 0.25; the shares compared lie 0.3 apart."""
BASELINE_SEED = 11
"""int: The seed of the draw with no prevalence change, which the checks that compare against the usual share share."""


def _features(system):
    planets, belts, _moons = system.count_objects()
    return {
        "habitable_world": any(system._barred(body) for planet in system.all_planets if planet.body_type != "a" for body in (planet, *planet.moons)),      # GEN.1397: the conditions of the request or a Blue or Green rating, never a class letter
        "asteroid_belt": belts > 0,
        "large_star": (system.primary_star.initial_mass_sol or 0.0) >= tuning.LARGE_STAR_MIN_MASS_SOL,
        "planets": planets + belts > 0,
        "comets": bool(system.comets or system.secondary_comets),
        "binary_system": system.binary_type is not None,
    }


def _shares(settings, seed, count=SYSTEMS):
    """The share of `count` drawn systems with each feature; a draw already made in this process (the same settings, seed and count) is
    reused, since drawing a system is most of this module's time (TEST.392)."""
    return _drawn(tuple(sorted(settings.items())), seed, count)


@functools.lru_cache(maxsize=None)
def _drawn(settings, seed, count):
    settings = dict(settings)
    draw.set_run_seed(seed)
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
    shares = _shares({}, seed=BASELINE_SEED)
    for feature, base in tuning.PREVALENCE_BASE_SHARES.items():
        assert shares[feature] == pytest.approx(base, abs=_tolerance(base)), feature
    # ADM.16: the Generate page shows these too.
    for feature in ("comets", "binary_system"):
        base = tuning.PREVALENCE_DRAW_SHARES[feature]
        assert shares[feature] == pytest.approx(base, abs=_tolerance(base)), feature


def test_every_feature_has_a_usual_share():
    assert set(prevalence.USUAL_SHARES) == set(prevalence.FEATURES)
    assert all(0 < share < 1 for share in prevalence.USUAL_SHARES.values())


@pytest.mark.parametrize("factor, percent", [(1.0, 0), (0.0, -100), (2.0, 100), (0.5, -50)])
def test_percent_for_share(factor, percent):
    usual = tuning.PREVALENCE_BASE_SHARES["habitable_world"]
    assert prevalence.percent_for_share("habitable_world", usual * factor) == pytest.approx(percent)


@pytest.mark.parametrize("feature, percent", [
    ("habitable_world", 50), ("habitable_world", -50), ("asteroid_belt", -60), ("asteroid_belt", 30),
    ("large_star", 100), ("large_star", -100), ("planets", -40),
])
def test_a_measured_feature_moves_by_the_percentage(feature, percent):
    expected = min(1.0, tuning.PREVALENCE_BASE_SHARES[feature] * (1 + percent / 100), prevalence.reachable_share(feature))    # GEN.1322: a habitable world is asked for, not placed
    # TEST.233: a seed that is the same in every process (hash() of a string differs per process), so a failure repeats
    seed = zlib.crc32(f"{feature}:{percent}".encode()) % 1000
    count = FORCED_ON_SYSTEMS if feature == "habitable_world" else SYSTEMS
    share = _shares({feature: percent}, seed=seed, count=count)[feature]
    assert share == pytest.approx(expected, abs=_tolerance(expected, count)), (feature, percent, seed, count, share)


@pytest.mark.parametrize("feature, base, percent", [
    ("comets", tuning.SYSTEM_COMET_CHANCE, 100), ("comets", tuning.SYSTEM_COMET_CHANCE, -100),
])
def test_a_drawn_feature_moves_by_the_percentage(feature, base, percent):
    # Each star of a wide pair draws its own comets, so those systems have
    # two chances, but a pair a few AU apart leaves each star a stable limit
    # under most comets' orbits (GEN.448): the second star's comets are mostly
    # dropped. The share lies between one chance and two for every wide pair.
    chance = min(1.0, base * (1 + percent / 100))
    wide = tuning.PREVALENCE_DRAW_SHARES["binary_system"] * tuning.PREVALENCE_DRAW_SHARES["wide_binary"]
    two = (1 - wide) * chance + wide * (1 - (1 - chance) ** 2)
    share = _shares({feature: percent}, seed=7, count=COMET_SYSTEMS)[feature]
    assert chance - _tolerance(chance, COMET_SYSTEMS) <= share <= two + _tolerance(two, COMET_SYSTEMS)
    if percent == 100:
        assert share > chance - _tolerance(chance, COMET_SYSTEMS) and share > _shares({}, seed=BASELINE_SEED)[feature]


def test_no_binaries_at_minus_100_and_more_at_plus_50():
    assert _shares({"binary_system": -100}, seed=3, count=40)["binary_system"] == 0    # none at all, so a few systems say it
    assert _shares({"binary_system": 50}, seed=3, count=COMET_SYSTEMS)["binary_system"] > _shares({}, seed=BASELINE_SEED)["binary_system"] + 0.02


def test_moons_and_max_planets_scale_their_draws():
    draw.set_run_seed(5)
    for percent, check in ((-100, lambda planets: all(not p.moons for p in planets)),):
        for _ in range(40):    # none at all, so a few systems say it (PERF.236)
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
    draw.set_run_seed(1)
    for _ in range(300):
        config = SystemConfig()
        config.COMETS = True
        # "all of them past 100%" needs p/100 * share / (1 - share) of 1 or more: a belt's usual share is 0.039, so +1000% forces it in about
        # 41% of systems and only +2464% and above in all (TEST.614: the test asked 1000 when the share was larger)
        config.PREVALENCE = {"habitable_world": 1000, "asteroid_belt": 10_000, "planets": -100}
        prevalence.resolve(config)
        assert config.HABITABLE_WORLD is True and config.ASTEROID_BELT is True
        assert config.LARGE_STAR is True  # both need a large star
        assert config.PLANETS is None  # and planets, so -planets gives way
    config = SystemConfig()
    config.LARGE_STAR = False
    config.HABITABLE_WORLD = True
    config.PREVALENCE = {"asteroid_belt": 10_000}    # forced on in every system, so the undo is what leaves it unset
    prevalence.resolve(config)
    assert config.ASTEROID_BELT is None and config.LARGE_STAR is False


def test_a_feature_a_forced_system_only_sometimes_has_is_forced_more_often_and_stops_at_its_hit_rate():
    """GEN.1322: +habitable_world forces just enough systems that the share with a habitable-class world comes to the wanted one when a forced system has one in
    `hit` of cases, and past `reachable_share` forces every system."""
    base, hit = tuning.PREVALENCE_BASE_SHARES["habitable_world"], tuning.PREVALENCE_FORCED_HIT_RATES["habitable_world"]
    assert prevalence.reachable_share("habitable_world") == hit and prevalence.reachable_share("asteroid_belt") == 1.0
    assert hit < 1.0 and base < hit

    def forced_share(percent, draws=4000):
        rng = random.Random(7)
        forced = 0
        for _ in range(draws):
            config = SystemConfig()
            config.PREVALENCE = {"habitable_world": percent}
            prevalence.resolve(config, rng)
            forced += config.HABITABLE_WORLD is True
        return forced / draws

    wanted = base * 1.2
    assert forced_share(20) == pytest.approx((wanted - base) / (hit - base), abs=0.03)
    assert forced_share(50) == pytest.approx((base * 1.5 - base) / (hit - base), abs=0.03)
    assert forced_share(200) == 1.0                 # +200 asks for 0.075, past the 0.069 a forced system gives
    assert forced_share(-100) == 0.0


def test_the_generate_page_says_a_habitable_world_is_only_asked_for():
    from planetgen.web import generate_page
    assert generate_page.PREVALENCE_LIMITS == {"prevalence_habitable_world": f"{prevalence.reachable_share('habitable_world') * 100:.3g}"}


def test_civilizations_scale_with_intelligent_life_prevalence():
    class _Timeline:
        life_stage = "technological_civilization"
    rng = random.Random(2)
    assert not any(has_civilization(_Timeline(), None, rng, -100) for _ in range(2000))
    hits = sum(has_civilization(_Timeline(), None, rng, 99_900) for _ in range(2000))
    assert hits == pytest.approx(2000 * tuning.CIVILIZATION_CHANCE * 1000, rel=0.15)


def test_large_star_minus_forbids_one():
    draw.set_run_seed(9)
    for _ in range(60):    # none at all, so a few systems say it (PERF.236)
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
