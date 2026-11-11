# tests/test_first_seed_with.py

"""
TEST.714: `fuzz_support.first_seed_with`, the bounded search a test makes for a case that a seed draws rarely instead of pinning the seed that once had it.
"""

import secrets

import pytest

from tests.fuzz_support import first_seed_with


def test_the_search_stops_at_the_first_seed_that_has_the_case_and_gives_it_back():
    built = []

    def build(seed):
        built.append(seed)
        return seed * 2

    seed, value = first_seed_with(build, lambda value: value >= 6, seeds=range(10))
    assert (seed, value) == (3, 6) and built == [0, 1, 2, 3]


def test_each_seed_is_built_under_the_same_draws_every_time():
    first = first_seed_with(lambda seed: secrets.randbits(16) / 65536.0, lambda value: value > 0.9, seeds=range(200))
    again = first_seed_with(lambda seed: secrets.randbits(16) / 65536.0, lambda value: value > 0.9, seeds=range(200))
    assert first == again and first[1] > 0.9


def test_a_case_no_seed_has_fails_the_test_and_says_how_many_were_tried():
    with pytest.raises(pytest.fail.Exception, match="none of 5 seeds had a moon"):
        first_seed_with(lambda seed: seed, lambda value: False, seeds=range(5), what="a moon")


def test_a_search_with_no_seeds_given_tries_the_default_bound():
    tried = []
    with pytest.raises(pytest.fail.Exception, match="none of 64 seeds"):
        first_seed_with(lambda seed: tried.append(seed), lambda value: False)
    assert len(tried) == 64
