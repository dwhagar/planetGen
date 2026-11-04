# tests/test_delivered_water.py

"""
GEN.849: a rocky body's water comes from the Bern delivery law (class from mass and snow-line distance, then a log-uniform share in the class), the
same law the atmosphere recipes b1 to b4 read.
"""

import math
import random

import pytest

from planetgen import tuning
from planetgen.physics import delivered_water, bulk_draw
from planetgen.physics.bulk_draw import ICY, ROCKY

DRAWS = 4000


def _fractions(mass_earth, snow_ratio, seed=7):
    rng = random.Random(seed)
    return [delivered_water.sample_fraction(mass_earth, snow_ratio, rng=rng) for _ in range(DRAWS)]


def _class_of(fraction):
    return delivered_water.CLASSES[sum(1 for edge in delivered_water.EDGES if fraction >= edge)]


@pytest.mark.parametrize("mass,ratio", [(0.5, 0.03), (1.0, 0.37), (4.0, 0.8), (30.0, 3.9)])
def test_class_shares_add_to_one_and_are_not_negative(mass, ratio):
    shares = delivered_water.class_shares(mass, ratio)
    assert set(shares) == set(delivered_water.CLASSES)
    assert all(value >= 0.0 for value in shares.values())
    assert sum(shares.values()) == pytest.approx(1.0)


def test_a_mass_or_distance_outside_the_model_is_held_to_its_range():
    assert delivered_water.class_shares(0.001, 0.37) == delivered_water.class_shares(delivered_water.MASS_RANGE_EARTH[0], 0.37)
    assert delivered_water.class_shares(1.0, 100.0) == delivered_water.class_shares(1.0, delivered_water.X_RANGE[1])


def test_the_brackets_tile_the_range_from_the_dry_floor_to_the_rich_cap():
    brackets = [delivered_water.bracket(name) for name in delivered_water.CLASSES]
    assert brackets[0][0] == tuning.DELIVERED_WATER_DRY_MIN
    assert brackets[-1][1] == tuning.DELIVERED_WATER_RICH_MAX
    for (_, high), (low, _) in zip(brackets, brackets[1:]):
        assert high == low
    assert all(low < high for low, high in brackets)


def test_the_draws_land_in_the_classes_the_law_gives():
    shares = delivered_water.class_shares(1.0, 0.37)
    fractions = _fractions(1.0, 0.37)
    assert all(tuning.DELIVERED_WATER_DRY_MIN <= f <= tuning.DELIVERED_WATER_RICH_MAX for f in fractions)
    for name in delivered_water.CLASSES:
        seen = sum(1 for f in fractions if _class_of(f) == name) / DRAWS
        assert seen == pytest.approx(shares[name], abs=0.03), name


def test_a_world_formed_near_the_star_is_drier_than_one_formed_at_the_snow_line():
    near = sum(1 for f in _fractions(1.0, 0.05) if _class_of(f) in ("high", "rich")) / DRAWS
    far = sum(1 for f in _fractions(1.0, 1.5) if _class_of(f) in ("high", "rich")) / DRAWS
    assert far > near


def test_two_draws_are_used_per_body():
    class Counter(random.Random):
        calls = 0

        def random(self):
            type(self).calls += 1
            return super().random()

    delivered_water.sample_fraction(1.0, 0.4, rng=Counter(3))
    assert Counter.calls == 2


class _Star:
    snow_line_au = 2.0


class _Body:
    def __init__(self, distance, mass_kg=5.97e24):
        self.star = _Star()
        self.distance = distance
        self.mass = mass_kg


def test_a_rocky_body_takes_the_delivery_law_and_an_icy_one_keeps_its_range(monkeypatch):
    rocky = _Body(0.8)
    assert bulk_draw.origin_of(rocky) == ROCKY
    water, ice, shares = bulk_draw._water_and_ice(rocky)
    assert delivered_water.bracket("dry")[0] <= water <= tuning.DELIVERED_WATER_RICH_MAX and ice == 0.0 and shares is None
    icy = _Body(5.0)
    assert bulk_draw.origin_of(icy) == ICY
    water, _ice, shares = bulk_draw._water_and_ice(icy)
    assert tuning.ICY_WATER_MASS_FRACTION_FLOOR * 0.5 <= water <= 0.95
    assert sum(shares.values()) == pytest.approx(1.0)


def test_a_body_with_no_mass_is_drawn_at_the_middle_of_the_models_masses():
    body = _Body(0.8)
    body.mass = None
    water, _ice, _shares = bulk_draw._water_and_ice(body)
    assert 0.0 < water <= tuning.DELIVERED_WATER_RICH_MAX
