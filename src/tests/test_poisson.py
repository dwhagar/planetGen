# tests/test_poisson.py

"""PERF.60: Poisson counts for small and huge means, from the project's draw wrapper."""

import math
import statistics

import pytest

from planetgen import tuning
from planetgen.util import draw, poisson


@pytest.mark.parametrize("mean", [0.3, 4.0, 9.9, 10.0, 25.0, 400.0, 5_000.0, 2_000_000.0])
def test_the_mean_and_variance_of_the_counts_match_the_mean(mean):
    rng = draw.Stream(f"poisson-test:{mean}")
    samples = [poisson.poisson_count(mean, rng) for _ in range(20_000)]
    assert all(isinstance(sample, int) and sample >= 0 for sample in samples)
    error = math.sqrt(mean / len(samples))
    assert abs(statistics.fmean(samples) - mean) < 5 * error + 1e-9
    assert statistics.pvariance(samples) == pytest.approx(mean, rel=0.06)


def test_the_distribution_matches_the_poisson_probabilities_just_above_the_threshold():
    mean = tuning.POISSON_REJECTION_MEAN * 1.5
    rng = draw.Stream("poisson-shape")
    n = 60_000
    seen = {}
    for _ in range(n):
        value = poisson.poisson_count(mean, rng)
        seen[value] = seen.get(value, 0) + 1
    for k in range(int(mean) - 4, int(mean) + 5):
        expected = math.exp(-mean + k * math.log(mean) - math.lgamma(k + 1))
        assert seen.get(k, 0) / n == pytest.approx(expected, abs=5 * math.sqrt(expected / n))


def test_a_large_mean_takes_a_few_draws_not_a_number_proportional_to_it():
    class Counting:
        def __init__(self):
            self.stream, self.calls = draw.Stream("count"), 0

        def random(self):
            self.calls += 1
            return self.stream.random()

    rng = Counting()
    for _ in range(100):
        poisson.poisson_count(5_000_000.0, rng)
    assert rng.calls < 100 * 8


def test_the_same_seed_gives_the_same_counts():
    first = [poisson.poisson_count(m, draw.Stream("same")) for m in (3.0, 50.0, 1e6)]
    again = [poisson.poisson_count(m, draw.Stream("same")) for m in (3.0, 50.0, 1e6)]
    assert first == again


def test_a_mean_of_zero_or_less_gives_zero_and_a_bad_mean_is_refused():
    assert poisson.poisson_count(0.0) == 0 and poisson.poisson_count(-4.0) == 0
    for bad in (math.nan, math.inf):
        with pytest.raises(ValueError):
            poisson.poisson_count(bad)
