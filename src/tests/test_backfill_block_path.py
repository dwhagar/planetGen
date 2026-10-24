# tests/test_backfill_block_path.py

"""
GEN.309: the backfill picks its draw path per block. A block of a band whose bound is under `CANDIDATE_PATH_BOUND`
draws its stars block first (`block_candidates`, each candidate on its own stream); a denser block draws cell by cell.
The choice depends on the block alone, so a cell gets the same stars whichever other cells a run reaches, in whatever
order, split between workers, or restricted to a sub-region (GEN.39).
"""

import functools
import math

import pytest

from planetgen.galaxy.spherical import wedge_count
from planetgen.generation import bright_stars as brightStars
from tests.test_bright_star_scatter import EDGE_PC, E_VALUE, SHAPE
from tests.test_mass_backfill import EDGES

SPARSE_SHELL = 120
"""int: A shell far enough out that its blocks are sparse at the test galaxy's rate (a bound of about 0.001 a cell)."""
SPARSE_RATE = E_VALUE * 20.0
"""float: The systems per sector that make `SPARSE_SHELL` bound a cell at about 0.03 stars in a band, under the limit."""
DENSE_SHELL = 40
"""int: A shell whose cells expect a star or more in a band at the test rate: the cell-by-cell path."""
SHELL = 60
"""int: A shell of two blocks whose cells expect a few stars a band at the test rate, with the candidate path forced
(`forced`): enough stars to compare, and exact whatever the bound."""
FLOOR = 100.0
SEED = 17
BAND = 0
CELLS = wedge_count(SHELL, BAND)
EDGE_OF_BLOCKS = brightStars.BLOCK_WEDGES
ACROSS_A_BLOCK = [(SHELL, BAND, wedge) for wedge in range(EDGE_OF_BLOCKS - 100, min(EDGE_OF_BLOCKS + 100, CELLS))]
"""list: Two hundred cells that straddle the first and the second block of the band."""
AGE_LAWS = brightStars.age_laws_of(SHAPE)


@functools.lru_cache(maxsize=None)
def _shares_of(k, low, high):
    return tuple((population, brightStars.bright_band_fraction(low, high, population, None, AGE_LAWS))
                 for population in brightStars.SCATTERED_POPULATIONS)


def _shares(k, low, high):
    return dict(_shares_of(k, low, high))


@pytest.fixture
def forced(monkeypatch):
    """Every block takes the candidate path: the draw is exact whatever the bound, and these cells hold stars."""
    monkeypatch.setattr(brightStars, "CANDIDATE_PATH_BOUND", float("inf"))


def _luminosity(addresses, rate=E_VALUE, memo=None, seed=SEED):
    return list(brightStars.backfill_cells(SHAPE, addresses, EDGE_PC, rate, FLOOR, None, seed, memo=memo))


def _masses(addresses, rate=E_VALUE, memo=None, seed=SEED, target=1.0, held=8.0):
    return list(brightStars.backfill_mass_cells(SHAPE, addresses, EDGE_PC, rate, target, held, None, seed, EDGES,
                                                memo=memo))


def _wedges(rows):
    return {row[2] for row in rows}


def test_a_sparse_block_takes_the_candidate_path_and_a_dense_one_does_not():
    bands = brightStars.canonical_bands(FLOOR)
    for k, low, high in bands:
        shares = _shares(k, low, high)
        assert brightStars.candidate_path({}, SHAPE, (SPARSE_SHELL, BAND, 3), EDGE_PC, SPARSE_RATE, SEED, k, shares,
                                          None) is not None
    k, low, high = bands[0]
    assert brightStars.candidate_path({}, SHAPE, (DENSE_SHELL, BAND, 3), EDGE_PC, E_VALUE, SEED, k, _shares(k, low, high),
                                      None) is None


def test_the_brightest_band_of_a_dense_shell_is_sparse_and_takes_the_candidate_path():
    """The path is chosen per block and band: where the dim bands draw cell by cell the bright ones are rare enough to
    be drawn block first."""
    k, low, high = brightStars.canonical_bands(FLOOR)[-1]
    assert brightStars.candidate_path({}, SHAPE, (DENSE_SHELL, BAND, 3), EDGE_PC, E_VALUE, SEED, k,
                                      _shares(k, low, high), None) is not None


def test_the_path_changes_where_the_blocks_bound_crosses_the_limit():
    k, low, high = brightStars.canonical_bands(FLOOR)[0]
    shares = _shares(k, low, high)
    bound = brightStars.block_bound(SHAPE, SPARSE_SHELL, BAND, EDGE_PC, SPARSE_RATE, shares)
    just_under = SPARSE_RATE * brightStars.CANDIDATE_PATH_BOUND * 0.99 / bound
    just_over = SPARSE_RATE * brightStars.CANDIDATE_PATH_BOUND * 1.01 / bound
    address = (SPARSE_SHELL, BAND, 3)
    assert brightStars.candidate_path({}, SHAPE, address, EDGE_PC, just_under, SEED, k, shares, None) is not None
    assert brightStars.candidate_path({}, SHAPE, address, EDGE_PC, just_over, SEED, k, shares, None) is None


def test_the_candidates_of_a_cell_are_the_blocks_candidates_in_that_cell():
    k, low, high = brightStars.canonical_bands(FLOOR)[0]
    shares = _shares(k, low, high)
    bound = brightStars.block_bound(SHAPE, SPARSE_SHELL, BAND, EDGE_PC, SPARSE_RATE, shares)
    block = brightStars.block_candidates(
        SEED, SPARSE_SHELL, BAND, 0, k, bound,
        lambda wedge: brightStars.cell_band_rate(SHAPE, (SPARSE_SHELL, BAND, wedge), EDGE_PC, SPARSE_RATE, shares))
    assert block
    memo = {}
    for wedge in {candidate.wedge for candidate in block} | {0, 5}:
        found = brightStars.candidate_path(memo, SHAPE, (SPARSE_SHELL, BAND, wedge), EDGE_PC, SPARSE_RATE, SEED, k,
                                           shares, None)
        assert found == [candidate for candidate in block if candidate.wedge == wedge]


def test_the_same_stars_come_out_in_whatever_order_or_split_the_cells_are_reached(forced):
    whole = sorted(_luminosity(ACROSS_A_BLOCK))
    assert len(whole) > 5 and len(_wedges(whole)) > 3
    assert sorted(_luminosity(ACROSS_A_BLOCK[::-1])) == whole
    assert sorted(row for address in ACROSS_A_BLOCK for row in _luminosity([address])) == whole        # one at a time
    shared = {}
    assert sorted(row for address in ACROSS_A_BLOCK[::-1] for row in _luminosity([address], memo=shared)) == whole
    halves = _luminosity(ACROSS_A_BLOCK[::2]) + _luminosity(ACROSS_A_BLOCK[1::2])                    # two workers
    assert sorted(halves) == whole


def test_a_sub_region_gets_the_stars_the_whole_region_gives_its_cells(forced):
    whole = _luminosity(ACROSS_A_BLOCK)
    inside = ACROSS_A_BLOCK[70:130]
    wedges = {address[2] for address in inside}
    assert sorted(_luminosity(inside)) == sorted(row for row in whole if row[2] in wedges)


def test_a_cell_taken_down_in_steps_holds_the_stars_one_draw_gives(forced):
    """GEN.44 on the candidate path: a band is drawn whole and cut afterwards."""
    cells = ACROSS_A_BLOCK
    deep = _luminosity(cells)
    ceiling = brightStars.canonical_bands(FLOOR)[4][1]
    top = list(brightStars.backfill_cells(SHAPE, cells, EDGE_PC, E_VALUE, ceiling, None, SEED))
    rest = list(brightStars.backfill_cells(SHAPE, cells, EDGE_PC, E_VALUE, FLOOR, ceiling, SEED))
    assert sorted(top + rest) == sorted(deep)


def test_every_star_lies_in_its_cell_and_in_its_band(forced):
    rows = _luminosity(ACROSS_A_BLOCK)
    for row in rows:
        shell_index, band_index, wedge = row[0:3]
        assert (shell_index, band_index) == (SHELL, BAND)
        assert (shell_index, band_index, wedge) in ACROSS_A_BLOCK
        assert brightStars._stored_in_cell(row[3:6], shell_index, band_index, wedge, EDGE_PC)


def test_the_candidate_path_gives_each_cell_its_own_poisson_mean(forced):
    """Over ten seeds the stars placed match the sum of the cells' means to 4.5 standard deviations (the counts are
    Poisson, so their total has the sum for its variance)."""
    bands = brightStars.canonical_bands(FLOOR)
    cells = [(SHELL, BAND, wedge) for wedge in range(CELLS)]
    per_seed = 0.0
    for k, low, high in bands:
        shares = _shares(k, low, high)
        per_seed += sum(brightStars.cell_band_rate(SHAPE, address, EDGE_PC, E_VALUE, shares) for address in cells)
    seeds = range(30, 40)
    found = sum(len(_luminosity(cells, seed=seed)) for seed in seeds)
    expected = per_seed * len(seeds)
    assert expected > 100
    assert abs(found - expected) < 4.5 * math.sqrt(expected), (found, expected)


def test_a_dense_block_still_draws_cell_by_cell():
    """Where a cell expects a star or more the draw is the stream of the cell, whatever else is reached."""
    dense = [(DENSE_SHELL, BAND, wedge) for wedge in range(12)]
    alone = [row for address in dense for row in _luminosity([address], rate=E_VALUE)]
    assert alone and sorted(alone) == sorted(_luminosity(dense[::-1], rate=E_VALUE))


def test_the_mass_backfill_takes_the_same_paths(forced):
    whole = sorted(_masses(ACROSS_A_BLOCK))
    assert len(whole) > 3
    assert sorted(_masses(ACROSS_A_BLOCK[::-1])) == whole
    shared = {}
    assert sorted(row for address in ACROSS_A_BLOCK for row in _masses([address], memo=shared)) == whole
    inside = ACROSS_A_BLOCK[60:140]
    wedges = {address[2] for address in inside}
    assert sorted(_masses(inside)) == sorted(row for row in whole if row[2] in wedges)


def test_the_memo_is_emptied_when_it_grows_past_its_limit_and_the_stars_do_not_change(monkeypatch, forced):
    whole = sorted(_luminosity(ACROSS_A_BLOCK))
    untrimmed = {}
    assert sorted(row for address in ACROSS_A_BLOCK for row in _luminosity([address], memo=untrimmed)) == whole
    monkeypatch.setattr(brightStars, "BLOCK_MEMO_LIMIT", 3)
    memo = {}
    assert sorted(row for address in ACROSS_A_BLOCK for row in _luminosity([address], memo=memo)) == whole
    assert len(memo) < len(untrimmed) / 2
