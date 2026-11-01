# tests/test_potassium_depletion.py

"""
GEN.492: the potassium depletion of a rocky body against CI chondrites (`planetgen.physics.potassium_depletion`, [R1003: Table 1]):
the table and its refusals, the two classes, and the draw.
"""

import json
import random

import pytest

from planetgen.physics import potassium_depletion as kd
from planetgen.physics.elements import DataError


def written(tmp_path, mutate):
    with open(kd.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    mutate(document)
    path = tmp_path / "table.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return str(path)


def test_the_table_gives_the_published_depletions():
    for name, depletion in (("bulk silicate Earth", 0.136), ("bulk silicate Mars", 0.243), ("Mercury (surface)", 0.200),
                            ("bulk silicate Moon", 0.0343), ("bulk silicate Vesta", 0.0580), ("angrite parent body", 0.0018)):
        assert kd.depletion_of(name) == pytest.approx(depletion, rel=0.01)


def test_ci_chondrites_have_the_k_to_u_of_barrat_2012():
    assert kd.ci_k_over_u() == pytest.approx(550.0 / 0.0077, rel=1e-3)


def test_the_classes_split_between_the_moon_and_mercury():
    assert kd.is_planet(0.0553) and kd.is_planet(1.0) and not kd.is_planet(0.0123) and not kd.is_planet(1e-5)
    assert kd.statistics_of(True)[0] == 3 and kd.statistics_of(False)[0] == 3


def test_a_planet_is_centred_on_a_fifth_and_a_small_body_on_a_twelfth_of_that():
    assert kd.median_depletion(1.0) == pytest.approx(0.19, abs=0.01)
    assert kd.median_depletion(0.01) == pytest.approx(0.0153, abs=0.001)
    assert kd.median_depletion(1.0) / kd.median_depletion(0.01) > 10


def test_the_predictive_spread_is_wider_than_the_spread_of_the_bodies():
    _, _, planet = kd.statistics_of(True)
    _, _, small = kd.statistics_of(False)
    assert 0.12 < planet < 0.2 and 0.8 < small < 1.1 and small > 5 * planet


def test_the_draw_follows_the_class_and_never_exceeds_one():
    rng = random.Random(7)
    planets = [kd.draw_depletion(1.0, rng) for _ in range(4000)]
    small = [kd.draw_depletion(0.001, rng) for _ in range(4000)]
    assert max(planets + small) <= kd.MAX_DEPLETION
    assert sorted(planets)[2000] == pytest.approx(kd.median_depletion(1.0), rel=0.05)
    assert sorted(small)[2000] == pytest.approx(kd.median_depletion(0.001), rel=0.2)
    assert sum(1 for d in small if d >= 0.136) > 0.02 * len(small)


def test_the_draw_uses_one_gaussian_whatever_the_class():
    class Counter:
        calls = 0

        def gauss(self, mu, sigma):
            Counter.calls += 1
            return mu
    for mass in (2.0, 0.0001):
        Counter.calls = 0
        kd.draw_depletion(mass, Counter())
        assert Counter.calls == 1


def test_a_name_that_is_not_in_the_table_is_refused():
    with pytest.raises(DataError):
        kd.depletion_of("Venus")


@pytest.mark.parametrize("mutate", [
    lambda d: d["reference"].update(k_over_u=-1.0),
    lambda d: d["bodies"].__setitem__(0, {"name": "x", "mass_earth": 1.0, "k_over_u": float("nan")}),
    lambda d: d.update(bodies=d["bodies"][:3]),
    lambda d: d.update(class_mass_earth=0),
])
def test_a_table_that_is_empty_or_not_positive_is_refused(tmp_path, mutate):
    with pytest.raises(DataError):
        kd.load(written(tmp_path, mutate))
