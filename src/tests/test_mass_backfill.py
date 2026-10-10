# tests/test_mass_backfill.py

"""
GEN.187: the mass backfill's pieces without a database. The population
model's mass bands add up with the luminosity pass, the fixed bands a sector
is drawn in do not depend on the steps taken to reach it, the stars drawn lie
in their band and below the luminosity already placed, and a white dwarf
stores its infinities as NULL and reads them back.
"""

import math
import random

import pytest

from planetgen.generation import bright_stars as brightStars
from planetgen.generation import star_population as sp
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.geometry import ring_sector_count, sector_address_at, sector_position_pc
from planetgen.physics import constants
from planetgen import tuning
from tests.test_bright_star_scatter import EDGE_PC, E_VALUE, SHAPE, THRESHOLD, _row_dict

EDGES = [1.0, 2.0, 5.0, 8.0, 20.0]


@pytest.mark.parametrize("population", [None, "young", "intermediate", "old"])
def test_mass_bands_add_up_to_the_open_band(population):
    whole = sp.mass_band_fraction(1.0, None, None, population)
    parts = (sp.mass_band_fraction(1.0, 2.0, None, population) + sp.mass_band_fraction(2.0, 5.0, None, population)
             + sp.mass_band_fraction(5.0, None, None, population))
    assert parts == pytest.approx(whole, rel=1e-9)
    assert whole == pytest.approx(sp.massive_star_fraction(1.0, population), rel=1e-12)
    assert 0.0 < sp.massive_star_fraction(2.0, population) < whole < 1.0


@pytest.mark.parametrize("population", ["young", "intermediate"])
def test_the_backfill_plus_the_luminosity_pass_is_every_star_born_in_the_band(population):
    # The stars born from 1 solar mass up are the dim ones a backfill draws plus the bright ones the scatter has.
    for low in (1.0, 2.0, 5.0):
        dim = sp.mass_band_fraction(low, None, THRESHOLD, population)
        bright = sp.bright_star_fraction(THRESHOLD, population, (low, None))
        assert dim + bright == pytest.approx(sp.massive_star_fraction(low, population), rel=1e-6)
        assert dim >= 0.0


def test_a_mass_band_needs_a_positive_lower_mass_and_a_real_range():
    with pytest.raises(ValueError):
        sp.mass_band_fraction(0.0)
    with pytest.raises(ValueError):
        sp.mass_band_fraction(2.0, 2.0)
    with pytest.raises(ValueError):
        sp.sample_mass_band_stars(1, 5.0, 3.0)


def test_sampled_stars_lie_in_their_band_and_under_the_luminosity_cap():
    rng = random.Random(4)
    stars = sp.sample_mass_band_stars(400, 2.0, 5.0, "young", rng, max_luminosity_sol=THRESHOLD)
    assert len(stars) == 400
    for star in stars:
        assert 2.0 * 0.999999 <= star["initial_mass_sol"] < 5.0
        assert star["yerkes_class"] == "VII" or star["luminosity_w"] < THRESHOLD * constants.SOLAR_LUMINOSITY
    # Old stars of these masses have become white dwarfs: the draw keeps them (the sector's fill refuses them).
    old = sp.sample_mass_band_stars(300, 1.0, 2.0, "old", rng)
    assert any(star["yerkes_class"] == "VII" for star in old)
    young_heavy = sp.sample_mass_band_stars(50, 8.0, None, "young", rng)
    assert all(star["initial_mass_sol"] >= 8.0 * 0.999999 for star in young_heavy)


def test_the_edges_hold_the_ring_masses_and_the_galaxys_limit():
    assert brightStars.mass_band_edges(None) == [1.0, 2.0, 5.0, 8.0]
    assert brightStars.mass_band_edges(20.0) == EDGES
    assert brightStars.mass_band_edges(3.0) == [1.0, 2.0, 3.0, 5.0, 8.0]
    assert brightStars.mass_band_edges(8.0) == [1.0, 2.0, 5.0, 8.0]
    assert brightStars.mass_band_edges(20.0, ring_masses=(4.0,)) == [4.0, 20.0]


def test_the_bands_a_draw_is_made_of_do_not_depend_on_where_it_started():
    bands = brightStars.canonical_mass_bands
    assert bands(1.0, 20.0, EDGES) == [(1.0, 2.0), (2.0, 5.0), (5.0, 8.0), (8.0, 20.0)]
    assert bands(5.0, 20.0, EDGES) == [(5.0, 8.0), (8.0, 20.0)]
    assert bands(1.0, 5.0, EDGES) == [(1.0, 2.0), (2.0, 5.0)]
    assert bands(8.0, None, [1.0, 2.0, 5.0, 8.0]) == [(8.0, None)]
    assert bands(1.0, None, [1.0, 2.0, 5.0, 8.0]) == [(1.0, 2.0), (2.0, 5.0), (5.0, 8.0), (8.0, None)]
    # A step from a held mass down to a target is the tail of the one draw from the top.
    assert bands(1.0, 5.0, EDGES) + bands(5.0, 20.0, EDGES) == bands(1.0, 20.0, EDGES)
    with pytest.raises(ValueError):
        bands(3.0, 20.0, EDGES)


ADDRESSES = [(3, 0, slot) for slot in range(ring_sector_count(3))] + [(5, 0, 7), (5, 0, 8)]


def _cells(target, held, cap, seed=11, addresses=ADDRESSES):
    return list(brightStars.backfill_mass_cells(SHAPE, addresses, EDGE_PC, E_VALUE, target, held, cap, seed, EDGES))


def test_a_mass_backfill_places_stars_in_their_cell_between_its_masses():
    rows = _cells(2.0, 8.0, THRESHOLD)
    assert len(rows) > 50
    for row in map(_row_dict, rows):
        point = (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000, row["position_z_mpc"] / 1000)
        address = (row["ring_index"], row["layer_index"], row["ring_slot_index"])
        assert address in ADDRESSES and sector_address_at(point, EDGE_PC) == address
        assert 2.0 * 0.999999 <= row["initial_mass_sol"] < 8.0
        assert row["yerkes_class"] == "VII" or row["luminosity_w"] < THRESHOLD * constants.SOLAR_LUMINOSITY


def test_the_same_seed_gives_the_same_stars_and_a_step_gives_the_tail_of_one_draw():
    assert _cells(1.0, 20.0, THRESHOLD) == _cells(1.0, 20.0, THRESHOLD)
    assert _cells(1.0, 20.0, THRESHOLD) != _cells(1.0, 20.0, THRESHOLD, seed=12)
    whole = _cells(1.0, 20.0, THRESHOLD)
    stepped = _cells(5.0, 20.0, THRESHOLD) + _cells(1.0, 5.0, THRESHOLD)
    assert sorted(stepped) == sorted(whole)


def test_the_count_matches_the_bands_expected_share():
    rows = _cells(1.0, 8.0, THRESHOLD, seed=5)
    expected = 0.0
    for address in ADDRESSES:
        center = sector_position_pc(*address, EDGE_PC)
        densities = brightStars._densities(center, SHAPE)
        for low, high in brightStars.canonical_mass_bands(1.0, 8.0, EDGES):
            expected += E_VALUE * sum(density * sp.mass_band_fraction(low, high, THRESHOLD, population)
                                      for population, density in densities.items())
    assert abs(len(rows) - expected) < 5 * math.sqrt(expected) + 5


def test_a_white_dwarf_row_stores_null_infinities_and_reads_them_back():
    rows = [_row_dict(row) for row in _cells(1.0, 20.0, None)]
    dwarfs = [row for row in rows if row["yerkes_class"] == "VII"]
    assert dwarfs and len(dwarfs) < len(rows)
    for row in dwarfs:
        assert row["lifespan_gy"] is None and row["phase_end_age_gy"] is None
        params = brightStars.star_params(row)
        assert params["lifespan_gy"] == math.inf and params["phase_end_age_gy"] == math.inf
    for row in rows:
        if row["yerkes_class"] != "VII":
            assert row["lifespan_gy"] is not None and math.isfinite(row["lifespan_gy"])
            assert brightStars.star_params(row)["lifespan_gy"] == row["lifespan_gy"]


def test_a_fill_below_a_backfilled_mass_caps_both_the_mass_and_the_luminosity():
    center = sector_position_pc(2, 0, 3, EDGE_PC)
    only_mass = brightStars.FillContext(center, SHAPE, star_mass_limit_sol=2.0)
    assert 0.02 < only_mass.bright_share() < 0.05   # about the share of stars born at 2 solar masses or more
    config = only_mass.apply(SystemConfig(), random.Random(1))
    assert config.MAX_STAR_MASS_SOL == 2.0 and config.MAX_STAR_LUMINOSITY_SOL is None

    both = brightStars.FillContext(center, SHAPE, min_luminosity_sol=THRESHOLD, star_mass_limit_sol=2.0)
    config = both.apply(SystemConfig(), random.Random(1))
    assert config.MAX_STAR_MASS_SOL == 2.0 and config.MAX_STAR_LUMINOSITY_SOL == THRESHOLD
    assert both.bright_share() >= only_mass.bright_share() - 1e-9

    neither = brightStars.FillContext(center, SHAPE)
    assert neither.bright_share() == 0.0
    config = neither.apply(SystemConfig(), random.Random(1))
    assert config.MAX_STAR_MASS_SOL is None and config.MAX_STAR_LUMINOSITY_SOL is None


def test_the_ring_masses_are_ascending_and_start_at_one_solar_mass():
    masses = tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL
    assert list(masses) == sorted(masses) and masses[0] == 1.0 and len(masses) == 4
