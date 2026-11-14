# tests/test_class_m_request_dial.py

"""
The Class M abundance tier's request dial (GEN.1704): the share of systems given the request for a habitable world.
`science` makes no draw and changes no system; a tier sets `HABITABLE_WORLD` (and `REQUEST_FROM_TIER`) in its share of the systems left to
chance, the star is drawn as without a request, and an unmet tier request is not reported as the user's.
"""

import random

import pytest

from planetgen import tuning
from planetgen.generation import abundance_tier, prevalence
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


class _Counting:
    """A random source that counts its draws."""

    def __init__(self, seed=3):
        self.rng, self.draws = random.Random(seed), 0

    def random(self):
        self.draws += 1
        return self.rng.random()


def _config(tier, **flags):
    config = SystemConfig()
    config.ABUNDANCE_TIER = tier
    for name, value in flags.items():
        setattr(config, name, value)
    return config


@pytest.mark.parametrize("tier", [None, "science"])
def test_science_makes_no_draw_and_sets_nothing(tier):
    rng = _Counting()
    for _ in range(50):
        config = _config(tier)
        prevalence.resolve(config, rng)
        assert config.HABITABLE_WORLD is None and config.REQUEST_FROM_TIER is False
    assert rng.draws == 0


@pytest.mark.parametrize("tier", abundance_tier.TIERS[1:])
def test_a_tier_makes_one_draw_only_when_its_share_is_over_zero(tier):
    rng = _Counting()
    share = abundance_tier.request_share(_config(tier))
    prevalence.resolve(_config(tier), rng)
    assert rng.draws == (1 if share > 0 else 0)


def test_share_one_requests_every_system_left_to_chance_and_share_half_about_half(monkeypatch):
    monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 1.0))
    for _ in range(100):
        config = _config("abundant")
        prevalence.resolve(config, random.Random(5))
        assert config.HABITABLE_WORLD is True and config.REQUEST_FROM_TIER is True
    monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 0.5))
    rng = random.Random(9)
    hits = 0
    for _ in range(2000):
        config = _config("abundant")
        prevalence.resolve(config, rng)
        hits += config.HABITABLE_WORLD is True
    assert hits / 2000 == pytest.approx(0.5, abs=0.04)


def test_a_config_that_already_says_true_or_false_or_cannot_hold_the_request_is_left_alone(monkeypatch):
    monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 1.0))
    rng = _Counting()
    for flags, expected in (({"HABITABLE_WORLD": False}, False), ({"HABITABLE_WORLD": True}, True),
                            ({"PLANETS": False}, None), ({"ASTEROID_BELT": True}, None)):
        config = _config("abundant", **flags)
        prevalence.resolve(config, rng)
        assert config.HABITABLE_WORLD is expected and config.REQUEST_FROM_TIER is False
    assert rng.draws == 0


def test_a_users_prevalence_resolves_first_and_stays_the_users_request(monkeypatch):
    monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 1.0))
    config = _config("abundant")
    config.PREVALENCE = {"habitable_world": 10_000}
    prevalence.resolve(config, random.Random(1))
    assert config.HABITABLE_WORLD is True and config.REQUEST_FROM_TIER is False


def _signature(system):
    return (system.star.mass, system.star.age, system.star.type, [(p.body_type, getattr(p, "mass", None), p.distance) for p in system.planets])


def test_science_builds_the_same_systems_as_no_tier():
    for seed in (101, 102, 103):
        signatures = []
        for tier in (None, "science"):
            config = _config(tier)
            config.PLANETS = True
            with draw.bound(seed):
                signatures.append(_signature(StarSystem(system_config=config)))
        assert signatures[0] == signatures[1]


def test_the_tier_request_leaves_the_host_star_as_drawn_without_a_request(monkeypatch):
    """The star sampler is told no habitable host for a tier's request (it would thin out white dwarfs, 36 to 5 of 500 systems) and told so for a user's."""
    from planetgen.generation import star as star_module
    monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 1.0))
    seen = []
    real = star_module.sample_living_star

    def spy(*args, **kwargs):
        seen.append(kwargs.get("habitable_host"))
        return real(*args, **kwargs)

    monkeypatch.setattr(star_module, "sample_living_star", spy)
    for tier, user, expected in ((None, True, True), ("abundant", None, False), (None, None, False)):
        seen.clear()
        config = _config(tier, HABITABLE_WORLD=user)
        config.PLANETS = False if user is None and tier is None else True
        with draw.bound(201):
            StarSystem(system_config=config)
        assert seen and seen[0] is expected, (tier, user, seen)
        assert config.REQUEST_FROM_TIER is (tier is not None and user is None)


def test_an_unmet_tier_request_says_nothing_and_an_unmet_user_request_says_unmet():
    holder = type("Holder", (), {})()
    holder.all_planets = []
    holder.system_config = _config("abundant", HABITABLE_WORLD=True)
    assert StarSystem.habitable_request_status(holder) == "unmet"
    holder.system_config.REQUEST_FROM_TIER = True
    assert StarSystem.habitable_request_status(holder) is None
