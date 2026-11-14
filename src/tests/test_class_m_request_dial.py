# tests/test_class_m_request_dial.py

"""
The Class M abundance tier's request dial (GEN.1704): the share of systems given the request for a habitable world, drawn by the `StarSystem` itself
after its stars and before its planets and kept on the system only (`StarSystem._tier_request`). `science` makes no draw and changes no system; the
config, the stars and the reported request status are never touched.
"""

import pytest

from planetgen import tuning
from planetgen.generation import abundance_tier
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


def _config(tier, **flags):
    config = SystemConfig()
    config.ABUNDANCE_TIER = tier
    config.PLANETS = True
    for name, value in flags.items():
        setattr(config, name, value)
    return config


def _stars_only(monkeypatch):
    """Planets are not the point of the star and draw checks: leave them out (the tier's draw is made before them)."""
    monkeypatch.setattr(StarSystem, "_generate_planets", lambda self, *args, **kwargs: [])


def _build(seed, tier, **flags):
    config = _config(tier, **flags)
    with draw.bound(seed):
        system = StarSystem(system_config=config)
    return system, config


def _star(system):
    star = system.stars[0]
    return (star.type, star.mass, star.age)


def test_the_request_is_a_system_attribute_and_never_enters_the_config(monkeypatch):
    _stars_only(monkeypatch)
    monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 1.0))
    system, config = _build(5, "abundant")
    assert system._tier_request is True and config.HABITABLE_WORLD is None
    assert system.habitable_request_status() is None and StarSystem._tier_request is False


def test_science_and_a_tier_without_a_request_share_make_the_same_draws(monkeypatch):
    _stars_only(monkeypatch)
    counts = {}
    real = draw.random
    for tier in (None, "science", "sparse"):         # sparse has request_share 0
        calls = []
        monkeypatch.setattr(draw, "random", lambda: calls.append(1) or real())
        system, _ = _build(7, tier)
        counts[tier] = len(calls)
        assert system._tier_request is False
    assert counts[None] == counts["science"] == counts["sparse"]


def test_a_tier_with_no_request_share_makes_the_same_system_as_science(monkeypatch):
    for seed in (101, 102, 103):
        signatures = []
        for tier in (None, "science", "sparse"):       # sparse has request_share 0
            system, _ = _build(seed, tier)
            signatures.append((_star(system), [(p.body_type, getattr(p, "mass", None), p.distance) for p in system.planets]))
        assert signatures[0] == signatures[1] == signatures[2]


def test_the_primary_star_is_identical_on_science_and_abundant_for_200_seeds(monkeypatch):
    _stars_only(monkeypatch)
    monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 1.0))
    for seed in range(5000, 5200):
        science, _ = _build(seed, "science")
        abundant, _ = _build(seed, "abundant")
        assert _star(science) == _star(abundant), seed
        assert abundant._tier_request is True and science._tier_request is False


def test_share_half_requests_about_half_of_the_systems(monkeypatch):
    _stars_only(monkeypatch)
    monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 0.5))
    hits = sum(_build(seed, "abundant")[0]._tier_request for seed in range(6000, 6200))
    assert hits / 200 == pytest.approx(0.5, abs=0.14)


@pytest.mark.parametrize("flags", [{"HABITABLE_WORLD": False}, {"HABITABLE_WORLD": True}, {"PLANETS": False}, {"ASTEROID_BELT": True}])
def test_a_config_the_user_set_is_left_alone(flags, monkeypatch):
    _stars_only(monkeypatch)
    monkeypatch.setitem(tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 1.0))
    system, config = _build(9, "abundant", **flags)
    assert system._tier_request is False
    assert config.HABITABLE_WORLD is flags.get("HABITABLE_WORLD")


def test_the_planets_are_drawn_under_the_request_a_tier_gave(monkeypatch):
    """Every planet-stage reader asks `_habitable_requested`, which counts the tier's request and the user's."""
    holder = type("Holder", (), {})()
    holder.system_config = _config("abundant")
    holder._tier_request = False
    assert StarSystem._habitable_requested(holder) is False
    holder._tier_request = True
    assert StarSystem._habitable_requested(holder) is True
    holder._tier_request, holder.system_config.HABITABLE_WORLD = False, True
    assert StarSystem._habitable_requested(holder) is True


def test_the_request_stands_check_runs_only_under_a_request():
    """`_habitable_request_stands` redraws planets and never stars; without a request it draws nothing."""
    holder = type("Holder", (), {})()
    holder.system_config, holder.planets = _config("abundant"), []
    holder._tier_request = False
    holder._habitable_requested = lambda: StarSystem._habitable_requested(holder)
    assert StarSystem._habitable_request_stands(holder) is True
