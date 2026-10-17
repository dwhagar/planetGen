# tests/test_block_sampler.py

"""
GEN.308: the block-first sampler of the backfill (`bright_stars.block_candidates`): one Poisson count from a block's
density bound, thinned to each cell's own rate. It is exact (a cell ends with a Poisson count of its own mean, whatever
the bound), keyed by the seed, the block and the luminosity band alone, and its bound holds for every cell.
"""

import collections
import math

import pytest

from planetgen.galaxy.spherical import band_range, wedge_count
from planetgen.generation import bright_stars as brightStars
from tests.test_bright_star_scatter import EDGE_PC, E_VALUE, SHAPE

SHELL, BAND = 120, 0                  # a band of 756 wedges: a full block of 256 and then some
CELLS = brightStars.BLOCK_WEDGES
REPLICATES = 3000
CHI_SQUARE_CRITICAL = {1: 10.83, 2: 13.82, 3: 16.27, 4: 18.47, 5: 20.52}      # p = 0.001
LAMBDAS = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0]                  # half-decade strata, the lowest first
BANDS = brightStars.canonical_bands(100.0)
SHARES = {population: brightStars.bright_band_fraction(BANDS[0][1], BANDS[0][2], population, None)
          for population in brightStars.POPULATIONS}


def _poisson_probabilities(mean, top):
    probabilities = [math.exp(-mean) * mean ** n / math.factorial(n) for n in range(top)]
    return probabilities + [max(1.0 - sum(probabilities), 0.0)]


def _multiplicities(mean, bound):
    """A histogram of how many stars each of the cells of REPLICATES blocks holds, from the real draw."""
    counts = collections.Counter()
    for replicate in range(REPLICATES):
        found = collections.Counter(c.wedge for c in brightStars.block_candidates(
            replicate, SHELL, BAND, 0, 16, bound, lambda wedge: mean))
        counts.update(collections.Counter(found.values()))
        counts[0] += CELLS - len(found)
    return counts


@pytest.mark.parametrize("mean", LAMBDAS)
def test_a_cells_count_is_poisson_whatever_the_bound(mean):
    """Stratified by the cell's mean (half decades, 3,000 blocks a stratum, the lowest strata being GEN.43's floor
    test): the multiplicities match exact Poisson probabilities, which drop-and-boost does not (P(2) would be
    mean^2 / (2q))."""
    counts = _multiplicities(mean, bound=4.0 * mean)
    cells = CELLS * REPLICATES
    expected = [cells * p for p in _poisson_probabilities(mean, 4)]
    observed = [counts[n] for n in range(4)] + [sum(v for n, v in counts.items() if n >= 4)]
    while len(expected) > 2 and expected[-1] < 5.0:          # merge the thin tail bins
        expected[-2] += expected.pop()
        observed[-2] += observed.pop()
    statistic = sum((o - e) ** 2 / e for o, e in zip(observed, expected))
    assert statistic < CHI_SQUARE_CRITICAL[len(expected) - 1], (mean, observed, expected, statistic)
    total = sum(n * v for n, v in counts.items())
    assert abs(total - mean * cells) < 4.5 * math.sqrt(mean * cells), (mean, total, mean * cells)


def test_a_cell_with_no_chance_is_never_kept_and_a_bound_of_zero_draws_nothing():
    assert brightStars.block_candidates(1, SHELL, BAND, 0, 16, 2.0, lambda wedge: 0.0) == []
    assert brightStars.block_candidates(1, SHELL, BAND, 0, 16, 0.0, lambda wedge: pytest.fail("asked")) == []


def test_a_cell_above_its_bound_is_refused():
    with pytest.raises(ValueError, match="above its block's bound"):
        brightStars.block_candidates(1, SHELL, BAND, 0, 16, 0.5, lambda wedge: 3.0)


def test_the_last_block_of_a_band_holds_only_the_wedges_that_are_left():
    cells = wedge_count(SHELL, BAND)
    last = (cells - 1) // CELLS
    wedges = {c.wedge for seed in range(40) for c in brightStars.block_candidates(
        seed, SHELL, BAND, last, 16, 1.0, lambda wedge: 1.0)}
    assert wedges and all(last * CELLS <= w < cells for w in wedges)
    assert brightStars.block_candidates(1, SHELL, BAND, cells // CELLS + 1, 16, 1.0, lambda wedge: 1.0) == []


# --- Keyed by the seed, the block and the band alone --------------------------------------------------------------


def _rate(wedge):
    return 0.05 + 0.01 * (wedge % 7)


def _block(seed, chunk, k=16):
    return brightStars.block_candidates(seed, SHELL, BAND, chunk, k, 0.12, _rate)


def test_a_block_gives_the_same_stars_whatever_region_order_or_worker_reaches_it():
    chunks = list(range(3))
    together = {chunk: _block(7, chunk) for chunk in chunks}
    assert all(together.values())
    assert {chunk: _block(7, chunk) for chunk in reversed(chunks)} == together
    assert {chunk: _block(7, chunk) for chunk in chunks[1:]} == {chunk: together[chunk] for chunk in chunks[1:]}
    assert _block(7, 0) == together[0]                                  # asked again later, alone


def test_another_seed_block_or_band_draws_another_set():
    base = _block(7, 0)
    assert _block(8, 0) != base and _block(7, 1) != base and _block(7, 0, k=17) != base


def test_a_candidates_star_stream_depends_only_on_its_place_in_the_block():
    one = brightStars.candidate_stream(7, SHELL, BAND, 0, 16, 3)
    again = brightStars.candidate_stream(7, SHELL, BAND, 0, 16, 3)
    other = brightStars.candidate_stream(7, SHELL, BAND, 0, 16, 4)
    assert [one.random() for _ in range(5)] == [again.random() for _ in range(5)] != [other.random() for _ in range(5)]


# --- The bound ----------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("shell", [1, 5, 20, 50, 150, 500, 2000])
@pytest.mark.parametrize("band", [0, 1, 3, 10, 30])
def test_every_cells_rate_is_at_most_its_blocks_bound(shell, band):
    if band > band_range(shell)[1]:
        pytest.skip("the shell has no such band")
    count = wedge_count(shell, band)
    for k, low, high in (BANDS[0], BANDS[8], BANDS[20]):
        shares = {population: brightStars.bright_band_fraction(low, high, population, None)
                  for population in brightStars.POPULATIONS}
        bound = brightStars.block_bound(SHAPE, shell, band, EDGE_PC, E_VALUE, shares)
        assert bound > 0.0
        for wedge in range(0, count, max(1, count // 400)):
            rate = brightStars.cell_band_rate(SHAPE, (shell, band, wedge), EDGE_PC, E_VALUE, shares)
            assert rate <= bound * (1.0 + 1e-9), (shell, band, wedge, k, rate, bound)


def test_a_real_block_keeps_each_cell_at_its_own_rate():
    """The block-first draw over a real band of the galaxy puts as many stars in the block as the cells' own rates
    sum to, and favours the denser cells."""
    shell, band, k = 30, 0, BANDS[8][0]
    shares = {population: brightStars.bright_band_fraction(BANDS[8][1], BANDS[8][2], population, None)
              for population in brightStars.POPULATIONS}
    bound = brightStars.block_bound(SHAPE, shell, band, EDGE_PC, E_VALUE, shares)
    first = {wedge: brightStars.cell_band_rate(SHAPE, (shell, band, wedge), EDGE_PC, E_VALUE, shares)
             for wedge in range(min(CELLS, wedge_count(shell, band)))}
    runs = 4000
    found = collections.Counter()
    for seed in range(runs):
        found.update(c.wedge for c in brightStars.block_candidates(seed, shell, band, 0, k, bound, first.__getitem__))
    expected = sum(first.values()) * runs
    assert abs(sum(found.values()) - expected) < 4.5 * math.sqrt(expected)
    dense = sorted(first, key=first.get)[-len(first) // 4:]
    sparse = sorted(first, key=first.get)[:len(first) // 4]
    expected_dense, expected_sparse = (sum(first[w] for w in group) * runs for group in (dense, sparse))
    assert abs(sum(found[w] for w in dense) - expected_dense) < 4.5 * math.sqrt(expected_dense)
    assert abs(sum(found[w] for w in sparse) - expected_sparse) < 4.5 * math.sqrt(expected_sparse)
    assert expected_dense > 2 * expected_sparse                          # the galaxy's bar and arms make a real contrast
