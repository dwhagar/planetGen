# tests/test_class_m_tier_water.py

"""
Class M tiers, thread 2: the water dial. Under a tier the share `abundance_tier.water_share(body)` of the candidates (rocky planets in the ecosphere)
take Earth's water class (`low`) in place of the Bern law's class, and the air of such a body takes the row b2. `science` changes nothing.
"""

import random

import pytest

from planetgen.generation import abundance_tier
from planetgen.generation.config import SystemConfig
from planetgen.generation.possibility_table import Subject
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import bulk_draw, delivered_water, mixture_rocky as md
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw
from tests.class_support import draw_planet
from tests.test_atmosphere_mixture_rocky import EARTH

LOW = delivered_water.bracket("low")


class _Star:
    snow_line_au = 2.0


class _Body:
    """A rocky planet at 1 AU of a star whose snow line is at 2 AU (a snow ratio of 0.5), under `tier`."""

    def __init__(self, tier=None, zone="e", is_moon=False, body_type="t"):
        self.system_config = SystemConfig()
        self.system_config.ABUNDANCE_TIER = tier
        self.star, self.distance, self.mass = _Star(), 1.0, 5.97e24
        self.zone, self.is_moon, self.body_type = zone, is_moon, body_type


def _waters(body, count, seed=11):
    out = []
    for index in range(count):
        with draw.bound(f"tier-water:{seed}:{index}"):
            out.append(bulk_draw._rocky_water(body))
    return out


def test_science_draws_exactly_as_before_the_dial():
    plain = random.Random(5)
    science = random.Random(5)
    before = [delivered_water.sample_fraction(1.0, 0.5, rng=plain) for _ in range(50)]
    after = [delivered_water.sample_fraction(1.0, 0.5, rng=science, earth=False) for _ in range(50)]
    assert before == after and plain.random() == science.random()
    assert _waters(_Body(None), 60) == _waters(_Body("science"), 60) == _waters(_Body(None, zone="h"), 60)


def test_the_dial_makes_one_extra_draw_only_under_a_tier_and_a_body_not_dialed_runs_the_law_untouched():
    class Counter(random.Random):
        calls = 0

        def random(self):
            type(self).calls += 1
            return super().random()

    assert not delivered_water.drawn_for_earth(0.0, rng=Counter(1)) and Counter.calls == 0
    assert delivered_water.drawn_for_earth(1.0, rng=Counter(1)) and Counter.calls == 1
    Counter.calls = 0
    delivered_water.sample_fraction(1.0, 0.4, rng=Counter(3), earth=True)
    assert Counter.calls == 1           # no class draw: the share inside Earth's class
    Counter.calls = 0
    delivered_water.sample_fraction(1.0, 0.4, rng=Counter(3), earth=False)
    assert Counter.calls == 2           # the law's own two
    assert delivered_water.sample_fraction(1.0, 0.4, rng=random.Random(8), earth=False) == delivered_water.sample_fraction(1.0, 0.4, rng=random.Random(8))


def test_at_share_one_every_candidate_holds_earths_class_and_nothing_else_moves(monkeypatch):
    monkeypatch.setitem(abundance_tier.tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 1.0))
    waters = _waters(_Body("abundant"), 200)
    assert all(LOW[0] <= water <= LOW[1] for water in waters)
    for other in (_Body("abundant", zone="h"), _Body("abundant", zone="c"), _Body("abundant", is_moon=True), _Body("abundant", zone=None)):
        assert _waters(other, 60) == _waters(_Body(None), 60)


def test_at_half_share_about_half_of_the_candidates_are_dialed(monkeypatch):
    monkeypatch.setitem(abundance_tier.tuning.ABUNDANCE_TIER_DIALS, "sparse", (0.5, 0.0, 0.0))
    natural = delivered_water.class_shares(bulk_draw._kg_to_mearth(5.97e24), 0.5)[delivered_water.EARTH_CLASS]
    count = 600
    low = sum(1 for water in _waters(_Body("sparse"), count) if LOW[0] <= water <= LOW[1])
    expected = 0.5 + 0.5 * natural          # the dialed half, and the law's own `low` among the rest
    assert abs(low / count - expected) < 0.06


def test_the_air_of_a_dialed_body_takes_the_row_of_its_class(monkeypatch):
    base = {**EARTH, "pinned_row_b": "b2"}
    for seed in range(40):
        with draw.bound(seed):
            recipe = RECIPES.draw_all(Subject("terrestrial", base), categories=md.ROCKY_CATEGORIES)
        assert recipe.groups["b"].name == "b2" and "b5" in recipe.values
    rows = set()
    for seed in range(60):
        with draw.bound(seed):
            rows.add(RECIPES.draw_all(Subject("terrestrial", {**EARTH, "snow_ratio": 0.3}), categories=md.ROCKY_CATEGORIES).groups["b"].name)
    assert len(rows) > 1            # without the pin the row is still the law's


def test_a_seeded_system_is_the_same_under_science_and_unset():
    def build(tier, seed):
        config = SystemConfig()
        config.STAR_TYPE = "G2V"
        config.ABUNDANCE_TIER = tier
        with draw.bound(f"tier-water-system:{seed}"):
            return StarSystem(system_config=config)

    def fingerprint(system):
        return [tuple(getattr(p, name, None) for name in ("body_type", "zone", "mass", "radius", "water_mass_fraction", "surface_pressure", "atmosphere")) for p in system.planets]

    for seed in range(3):
        assert fingerprint(build(None, seed)) == fingerprint(build("science", seed))


def test_whole_planets_in_the_ecosphere_of_a_sun_like_star_are_low_under_the_top_tier_and_their_air_agrees():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.ABUNDANCE_TIER = "abundant"
    with draw.bound("tier-water-star"):
        star = Star(config)
    found = 0
    for seed in range(40):
        planet = draw_planet(star, star.habitable_zone[0] * 1.1, seed, config, moon_count=0)
        if abundance_tier.candidate(planet) and bulk_draw.origin_of(planet) == bulk_draw.ROCKY:
            assert planet._water_dialed and planet.water_mass_fraction <= LOW[1] * 1.0001    # later water loss can only lower the stored share
            found += 1
    assert found >= 10
