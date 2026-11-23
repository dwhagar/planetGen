# tests/test_class_m_tier_water.py

"""
Class M tiers, thread 2: the water dial. Under a tier the share `abundance_tier.water_share(body)` of the candidates (rocky planets in the ecosphere)
take a water class from the tier in place of the Bern law's class: the part `abundance_tier.dry_part(body)` of them the dry class (GEN.1721, the water of
a desert world) and the rest Earth's (`low`); the air of such a body takes the row of its class (b1 or b2). `science` changes nothing.
"""

import random

import pytest

from planetgen.generation import abundance_tier
from planetgen.generation.config import SystemConfig
from planetgen.generation.possibility_table import Subject
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import bulk_draw, delivered_water, formation_snow_line, mixture_rocky as md
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw
from tests.class_support import draw_planet
from tests.test_atmosphere_mixture_rocky import EARTH

LOW = delivered_water.bracket("low")
DRY = delivered_water.bracket("dry")


class _Star:
    """A one-solar-mass star whose stored line is set so that the water snow line the draw reads (`formation_snow_line.water_snow_line_au`, GEN.1678) is 2 AU."""
    initial_mass_sol = 1.0
    snow_line_au = 2.0 * formation_snow_line.formation_snow_line_au(1.0, 0.0) / formation_snow_line.viscous_snow_line_au(1.0)


class _Body:
    """A rocky planet at 1 AU of a star whose snow line is at 2 AU (a snow ratio of 0.5), under `tier`."""

    def __init__(self, tier=None, zone="e", is_moon=False, body_type="t"):
        self.system_config = SystemConfig()
        self.system_config.ABUNDANCE_TIER = tier
        self.star, self.distance, self.mass = _Star(), 1.0, 5.97e24
        # the body stays at 1 AU (its water draw reads the distance) and the zone moves: the habitable zone widened by half its width is the effort zone (GEN.1807)
        self.habitable_zone = {"e": (0.95, 1.37), "h": (1.7, 2.4), "c": (0.2, 0.5), None: (1.7, 2.4)}[zone]
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
    after = [delivered_water.sample_fraction(1.0, 0.5, rng=science, dialed=None) for _ in range(50)]
    assert before == after and plain.random() == science.random()
    assert _waters(_Body(None), 60) == _waters(_Body("science"), 60) == _waters(_Body(None, zone="h"), 60)


def test_the_dial_makes_one_extra_draw_only_under_a_tier_and_a_body_not_dialed_runs_the_law_untouched():
    class Counter(random.Random):
        calls = 0

        def random(self):
            type(self).calls += 1
            return super().random()

    assert delivered_water.dialed_class(0.0, 0.5, rng=Counter(1)) is None and Counter.calls == 0
    assert delivered_water.dialed_class(1.0, rng=Counter(1)) == "low" and Counter.calls == 1
    Counter.calls = 0
    assert delivered_water.dialed_class(1.0, 1.0, rng=Counter(1)) == "dry" and Counter.calls == 1     # the dry part splits the same draw
    for dialed in ("low", "dry"):
        Counter.calls = 0
        water = delivered_water.sample_fraction(1.0, 0.4, rng=Counter(3), dialed=dialed)
        assert Counter.calls == 1           # no class draw: the share inside the dialed class
        low, high = abundance_tier.tuning.ABUNDANCE_DESERT_WATER_RANGE if dialed == "dry" else delivered_water.bracket(dialed)
        assert low <= water <= high          # a desert world draws inside the part of the dry class that leaves small seas (GEN.1721)
    Counter.calls = 0
    delivered_water.sample_fraction(1.0, 0.4, rng=Counter(3), dialed=None)
    assert Counter.calls == 2           # the law's own two
    assert delivered_water.sample_fraction(1.0, 0.4, rng=random.Random(8), dialed=None) == delivered_water.sample_fraction(1.0, 0.4, rng=random.Random(8))


def test_the_dry_part_splits_the_dials_one_draw_and_a_part_of_zero_draws_as_the_dial_did_alone():
    picks = [random.Random(seed).random() for seed in range(400)]
    for seed, pick in enumerate(picks):
        alone = delivered_water.dialed_class(0.5, rng=random.Random(seed))
        split = delivered_water.dialed_class(0.5, 0.2, rng=random.Random(seed))
        assert (alone is None) == (split is None) == (pick >= 0.5)       # the same bodies are reached, whatever the part
        assert alone in (None, "low") and split == (None if pick >= 0.5 else "dry" if pick < 0.1 else "low")
    dry = sum(1 for seed in range(400) if delivered_water.dialed_class(1.0, 0.2, rng=random.Random(seed)) == "dry")
    assert 50 <= dry <= 110 and dry < 200          # about a fifth, and Earth's class stays the larger part


def test_at_share_one_every_candidate_holds_earths_class_and_nothing_else_moves(monkeypatch):
    monkeypatch.setitem(abundance_tier.tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 1.0, 0.0, 0.0, 0.0))
    waters = _waters(_Body("abundant"), 200)
    assert all(LOW[0] <= water <= LOW[1] for water in waters)
    for other in (_Body("abundant", zone="h"), _Body("abundant", zone="c"), _Body("abundant", is_moon=True), _Body("abundant", zone=None)):
        assert _waters(other, 60) == _waters(_Body(None), 60)


def test_with_a_dry_part_every_candidate_holds_a_dialed_class_the_dry_one_the_smaller_part_and_remembers_which(monkeypatch):
    monkeypatch.setitem(abundance_tier.tuning.ABUNDANCE_TIER_DIALS, "abundant", (1.0, 1.0, 1.0, 0.25, 0.0, 0.0))
    body = _Body("abundant")
    dry = 0
    for index in range(300):
        with draw.bound(f"tier-dry:{index}"):
            water = bulk_draw._rocky_water(body)
        assert body._water_dialed in ("dry", "low")
        low, high = abundance_tier.tuning.ABUNDANCE_DESERT_WATER_RANGE if body._water_dialed == "dry" else delivered_water.bracket(body._water_dialed)
        assert low <= water <= high
        dry += body._water_dialed == "dry"
    assert abs(dry / 300 - 0.25) < 0.08
    for tier in (None, "science"):
        plain = _Body(tier)
        _waters(plain, 5)
        assert not hasattr(plain, "_water_dialed")           # science leaves the body untouched (the golden digests list every field)


def test_at_half_share_about_half_of_the_candidates_are_dialed(monkeypatch):
    monkeypatch.setitem(abundance_tier.tuning.ABUNDANCE_TIER_DIALS, "abundant", (0.5, 0.0, 0.0, 0.0, 0.0, 0.0))
    natural = delivered_water.class_shares(bulk_draw._kg_to_mearth(5.97e24), 0.5)[delivered_water.EARTH_CLASS]
    count = 600
    low = sum(1 for water in _waters(_Body("abundant"), count) if LOW[0] <= water <= LOW[1])
    expected = 0.5 + 0.5 * natural          # the dialed half, and the law's own `low` among the rest
    assert abs(low / count - expected) < 0.06


def test_the_air_of_a_dialed_body_takes_the_row_of_its_class(monkeypatch):
    base = {**EARTH, "pinned_row_b": "b2"}
    for seed in range(40):
        with draw.bound(seed):
            recipe = RECIPES.draw_all(Subject("terrestrial", base), categories=md.ROCKY_CATEGORIES)
        assert recipe.groups["b"].name == "b2" and "b5" in recipe.values
    assert delivered_water.ROWS["b1"] == "dry" and delivered_water.ROWS["b2"] == "low"       # the rows `mixture_subject` pins for the two dialed classes
    for seed in range(40):
        with draw.bound(seed):
            recipe = RECIPES.draw_all(Subject("terrestrial", {**EARTH, "pinned_row_b": "b1"}), categories=md.ROCKY_CATEGORIES)
        assert recipe.groups["b"].name == "b1"
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


def test_whole_planets_in_the_ecosphere_of_a_sun_like_star_hold_a_dialed_class_under_the_top_tier_and_their_air_agrees():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.ABUNDANCE_TIER = "abundant"
    with draw.bound("tier-water-star"):
        star = Star(config)
    found = 0
    for seed in range(40):
        planet = draw_planet(star, star.habitable_zone[0] * 1.1, seed, config, moon_count=0)
        if abundance_tier.candidate(planet) and bulk_draw.origin_of(planet) == bulk_draw.ROCKY:
            assert planet._water_dialed in ("dry", "low")
            assert planet.water_mass_fraction <= delivered_water.bracket(planet._water_dialed)[1] * 1.0001    # later water loss can only lower the stored share
            found += 1
    assert found >= 10
