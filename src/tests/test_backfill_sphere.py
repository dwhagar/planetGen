# tests/test_backfill_sphere.py

"""
VIEW.41: the sphere draw (`backfill_sphere`, `backfill_mass_sphere`) draws a sphere's backfill by blocks: it gives exactly the
stars `backfill_cells` and `backfill_mass_cells` give for the cells `enumerate_sectors_within_radius` lists, on the candidate
path, the cell path and a mix of the two, for a sphere that wraps past longitude 0, one that holds the core cell, and one
on the pole.
"""

import math

import pytest

from planetgen.galaxy import spherical
from planetgen.generation import bright_stars as brightStars
from tests.test_bright_star_scatter import EDGE_PC, E_VALUE, SHAPE
from tests.test_mass_backfill import EDGES

SEED = 29
FLOOR = 50.0
RATE = E_VALUE
CENTERS = {
    "wrapped": ((240.0, 0.0, 0.0), 16.0),
    "off the plane": ((-150.0, 180.0, 20.0), 16.0),
    "core": ((0.0, 0.0, 0.0), 6.0),
    "pole": ((0.0, 0.0, 60.0), 16.0),
}


def _cells(center, radius):
    return [row[:3] for row in spherical.enumerate_sectors_within_radius(center, radius, EDGE_PC)]


def _luminosity_by_cells(center, radius, rate=RATE, mass_range=None):
    return sorted(brightStars.backfill_cells(SHAPE, _cells(center, radius), EDGE_PC, rate, FLOOR, None, SEED, mass_range=mass_range))


def _luminosity_by_sphere(center, radius, rate=RATE, mass_range=None):
    return sorted(brightStars.backfill_sphere(SHAPE, center, radius, EDGE_PC, rate, FLOOR, None, SEED, mass_range=mass_range))


def _mass_by_cells(center, radius, rate=RATE):
    return sorted(brightStars.backfill_mass_cells(SHAPE, _cells(center, radius), EDGE_PC, rate, 1.0, 8.0, None, SEED, EDGES))


def _mass_by_sphere(center, radius, rate=RATE):
    return sorted(brightStars.backfill_mass_sphere(SHAPE, center, radius, EDGE_PC, rate, 1.0, 8.0, None, SEED, EDGES))


@pytest.mark.parametrize("name", sorted(CENTERS))
def test_sphere_rows_list_the_cells_enumerate_lists(name):
    center, radius = CENTERS[name]
    listed = {(shell, band, (first + k) % count)
              for shell, band, count, first, length in spherical.sphere_rows(center, radius, EDGE_PC) for k in range(length)}
    kept = {cell for cell in listed if spherical.cell_in_sphere(cell, center, radius, EDGE_PC)}
    expected = set(_cells(center, radius)) - {(0, 0, 0)}
    assert kept == expected
    assert len(listed) >= len(kept)


def test_sphere_rows_refuse_a_bad_sphere():
    with pytest.raises(ValueError):
        list(spherical.sphere_rows((0.0, 0.0, 0.0), -1.0, EDGE_PC))
    with pytest.raises(ValueError):
        list(spherical.sphere_rows((0.0, 0.0, 0.0), 1.0, 0.0))


def test_a_bound_from_its_terms_is_the_block_bound_to_the_bit():
    shares = {population: brightStars.bright_band_fraction(100.0, 133.0, population, None, brightStars.age_laws_of(SHAPE))
              for population in brightStars.SCATTERED_POPULATIONS}
    for shell, band in ((40, 3), (120, 0), (200, 17)):
        terms = brightStars.bound_terms(SHAPE, shell, band, EDGE_PC, RATE)
        assert brightStars.bound_of(terms, shares) == brightStars.block_bound(SHAPE, shell, band, EDGE_PC, RATE, shares)


@pytest.mark.parametrize("name", sorted(CENTERS))
@pytest.mark.parametrize("limit", [0.0, 0.05, float("inf")], ids=["cell path", "default mix", "candidate path"])
def test_the_luminosity_sphere_is_the_cells_draw(monkeypatch, name, limit):
    center, radius = CENTERS[name]
    if limit != 0.05:
        monkeypatch.setattr(brightStars, "CANDIDATE_PATH_BOUND", limit)
    expected = _luminosity_by_cells(center, radius)
    assert expected, "the test sphere must hold stars"
    assert _luminosity_by_sphere(center, radius) == expected


@pytest.mark.parametrize("name", ["wrapped"])
@pytest.mark.parametrize("limit", [0.0, float("inf")], ids=["cell path", "candidate path"])
def test_the_mass_sphere_is_the_cells_draw(monkeypatch, name, limit):
    center, radius = CENTERS[name]
    monkeypatch.setattr(brightStars, "CANDIDATE_PATH_BOUND", limit)
    expected = _mass_by_cells(center, radius)
    assert expected
    assert _mass_by_sphere(center, radius) == expected


def test_a_mass_range_and_a_ceiling_pass_through(monkeypatch):
    monkeypatch.setattr(brightStars, "CANDIDATE_PATH_BOUND", float("inf"))
    center, radius = CENTERS["wrapped"]
    cells = _cells(center, radius)
    expected = sorted(brightStars.backfill_cells(SHAPE, cells, EDGE_PC, RATE, 50.0, 400.0, SEED, mass_range=(None, 14.0)))
    got = sorted(brightStars.backfill_sphere(SHAPE, center, radius, EDGE_PC, RATE, 50.0, 400.0, SEED, mass_range=(None, 14.0)))
    assert expected and got == expected


def test_batches_do_not_change_the_draw(monkeypatch):
    monkeypatch.setattr(brightStars, "CANDIDATE_PATH_BOUND", float("inf"))
    center, radius = CENTERS["wrapped"]
    expected = _luminosity_by_cells(center, radius)
    monkeypatch.setattr(brightStars, "SPHERE_BATCH_CELLS", 3)
    monkeypatch.setattr(brightStars, "SPHERE_BATCH_ROWS", 2)
    assert _luminosity_by_sphere(center, radius) == expected


def test_a_shared_memo_keeps_only_the_shares():
    center, radius = CENTERS["wrapped"]
    memo = {}
    first = sorted(brightStars.backfill_sphere(SHAPE, center, radius, EDGE_PC, RATE, FLOOR, None, SEED, memo=memo))
    assert set(key[0] for key in memo) == {"shares"}
    assert sorted(brightStars.backfill_sphere(SHAPE, center, radius, EDGE_PC, RATE, FLOOR, None, SEED, memo=memo)) == first


def test_an_empty_sphere_draws_nothing():
    assert list(brightStars.backfill_sphere(SHAPE, (1.0e6, 0.0, 0.0), 10.0, EDGE_PC, RATE, FLOOR, None, SEED)) == []


def test_the_sphere_visits_blocks_not_cells(monkeypatch):
    """A sparse sphere asks the cell draw about the cells that hold a candidate, not the sphere's cells."""
    monkeypatch.setattr(brightStars, "CANDIDATE_PATH_BOUND", float("inf"))
    center, radius = CENTERS["wrapped"]
    seen = []
    real = brightStars.backfill_cells

    def spy(shape, addresses, *args, **kwargs):
        seen.append(len(addresses))
        return real(shape, addresses, *args, **kwargs)

    monkeypatch.setattr(brightStars, "backfill_cells", spy)
    list(brightStars.backfill_sphere(SHAPE, center, radius, EDGE_PC, E_VALUE, FLOOR, None, SEED))
    assert sum(seen) < len(_cells(center, radius)) / 4
