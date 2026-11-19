# tests/test_block_sampler.py

"""
GEN.308, PERF.280: the row-first sampler of the backfill (`bright_stars.row_candidates`): one Poisson count from a row's
density bound, thinned to each cell's own rate. It is exact (a cell ends with a Poisson count of its own mean, whatever
the bound), keyed by the seed, the row and the luminosity band alone, and its bound holds for every cell.
"""

import collections
import math

import pytest

from planetgen.galaxy.spherical import band_range, wedge_count
from planetgen.generation import bright_stars as brightStars
from tests.test_bright_star_scatter import EDGE_PC, E_VALUE, SHAPE

SHELL, BAND = 120, 0                  # a band of 756 wedges
CELLS = wedge_count(SHELL, BAND)
REPLICATES = 3000
CHI_SQUARE_CRITICAL = {1: 10.83, 2: 13.82, 3: 16.27, 4: 18.47, 5: 20.52}      # p = 0.001
LAMBDAS = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0]                  # half-decade strata, the lowest first
BANDS = brightStars.canonical_bands(100.0)
SHARES = {population: brightStars.bright_band_fraction(BANDS[0][1], BANDS[0][2], population, None)
          for population in brightStars.POPULATIONS}


def _poisson_probabilities(mean, top):
    probabilities = [math.exp(-mean) * mean ** n / math.factorial(n) for n in range(top)]
    return probabilities + [max(1.0 - sum(probabilities), 0.0)]


def _kept(seed, k, bound, rate_of, shell=SHELL, band=BAND):
    """The candidates of a row that the thinning keeps, as `candidate_path` keeps them: each with the chance `rate / bound` of its own cell."""
    kept = []
    for wedge, found in brightStars.row_candidates(seed, shell, band, k, bound).items():
        rate = rate_of(wedge)
        assert rate <= bound * (1.0 + 1e-9)
        kept += [c for c in found if c.accept * bound < rate]
    return kept


def _multiplicities(mean, bound):
    """A histogram of how many stars each of the cells of REPLICATES rows holds, from the real draw."""
    counts = collections.Counter()
    for replicate in range(REPLICATES):
        found = collections.Counter(c.wedge for c in _kept(replicate, 16, bound, lambda wedge: mean))
        counts.update(collections.Counter(found.values()))
        counts[0] += CELLS - len(found)
    return counts


@pytest.mark.parametrize("mean", LAMBDAS)
def test_a_cells_count_is_poisson_whatever_the_bound(mean):
    """Stratified by the cell's mean (half decades, 3,000 rows a stratum, the lowest strata being GEN.43's floor
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
    assert _kept(1, 16, 2.0, lambda wedge: 0.0) == []
    assert brightStars.row_candidates(1, SHELL, BAND, 16, 0.0) == {}


def test_a_cell_above_its_bound_is_refused():
    shares = {population: brightStars.bright_band_fraction(BANDS[0][1], BANDS[0][2], population, None) for population in brightStars.POPULATIONS}
    memo = {("candidates", SHELL, BAND, 16, None): {5: [brightStars.BlockCandidate(5, 0, 0.5, 0.1)]}}
    with pytest.raises(ValueError, match="above its row's bound"):
        brightStars._kept_candidates(memo, SHAPE, (SHELL, BAND, 5), EDGE_PC, E_VALUE, 1, 16, shares, None, 1e-12)


def test_a_row_holds_only_its_own_wedges():
    wedges = {wedge for seed in range(40) for wedge in brightStars.row_candidates(seed, SHELL, BAND, 16, 1.0)}
    assert wedges and all(0 <= wedge < CELLS for wedge in wedges)


# --- Keyed by the seed, the block and the band alone --------------------------------------------------------------


def _rate(wedge):
    return 0.05 + 0.01 * (wedge % 7)


def _row(seed, band=BAND, k=16):
    return _kept(seed, k, 0.12, _rate, band=band)


def test_a_row_gives_the_same_stars_whatever_region_order_or_worker_reaches_it():
    rows = [0, 1, 2]
    together = {band: _row(7, band) for band in rows}
    assert all(together.values())
    assert {band: _row(7, band) for band in reversed(rows)} == together
    assert _row(7, 0) == together[0]                                  # asked again later, alone


def test_another_seed_row_or_band_draws_another_set():
    base = _row(7, 0)
    assert _row(8, 0) != base and _row(7, 1) != base and _row(7, 0, k=17) != base


def test_a_candidates_star_stream_depends_only_on_its_place_in_the_row():
    one = brightStars.candidate_stream(7, SHELL, BAND, brightStars.ROW_STREAM, 16, 3)
    again = brightStars.candidate_stream(7, SHELL, BAND, brightStars.ROW_STREAM, 16, 3)
    other = brightStars.candidate_stream(7, SHELL, BAND, brightStars.ROW_STREAM, 16, 4)
    assert [one.random() for _ in range(5)] == [again.random() for _ in range(5)] != [other.random() for _ in range(5)]


# --- The bound ----------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("shell", [1, 5, 20, 50, 150, 500, 2000])
@pytest.mark.parametrize("band", [0, 1, 3, 10, 30])
def test_every_cells_rate_is_at_most_its_rows_bound(shell, band):
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


def test_a_real_row_keeps_each_cell_at_its_own_rate():
    """The row-first draw over a real band of the galaxy puts as many stars in the row as the cells' own rates
    sum to, and favours the denser cells."""
    shell, band, k = 30, 0, BANDS[8][0]
    shares = {population: brightStars.bright_band_fraction(BANDS[8][1], BANDS[8][2], population, None)
              for population in brightStars.POPULATIONS}
    bound = brightStars.block_bound(SHAPE, shell, band, EDGE_PC, E_VALUE, shares)
    first = {wedge: brightStars.cell_band_rate(SHAPE, (shell, band, wedge), EDGE_PC, E_VALUE, shares)
             for wedge in range(wedge_count(shell, band))}
    runs = 4000
    found = collections.Counter()
    for seed in range(runs):
        found.update(c.wedge for c in _kept(seed, k, bound, first.__getitem__, shell, band))
    expected = sum(first.values()) * runs
    assert abs(sum(found.values()) - expected) < 4.5 * math.sqrt(expected)
    dense = sorted(first, key=first.get)[-len(first) // 4:]
    sparse = sorted(first, key=first.get)[:len(first) // 4]
    expected_dense, expected_sparse = (sum(first[w] for w in group) * runs for group in (dense, sparse))
    assert abs(sum(found[w] for w in dense) - expected_dense) < 4.5 * math.sqrt(expected_dense)
    assert abs(sum(found[w] for w in sparse) - expected_sparse) < 4.5 * math.sqrt(expected_sparse)
    assert expected_dense > 2 * expected_sparse                          # the galaxy's bar and arms make a real contrast
