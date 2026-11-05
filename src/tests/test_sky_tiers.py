"""VIEW.24: the sky tier table (`planetgen/generation/sky_tiers.py`): the floor against distance, the reach of a band through the dust, the cut by a budget
of cell visits, and the split between the luminosity and the mass ledger."""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import spherical
from planetgen.generation import sky_tiers as tiers_module
from planetgen.generation.sky_tiers import build, cells_in_sphere, extinction_curve, floor_luminosity_sol, floor_table, reach_pc
from planetgen.physics import photometry

EDGE = 4.0
SUN = (8000.0, 0.0, 0.0)


@pytest.fixture(scope="module")
def clear():
    """A dust-free sky from the Sun's place."""
    return extinction_curve(SUN, count=12, av_per_kpc_at_sun=0.0)


@pytest.fixture(scope="module")
def dusty():
    return extinction_curve(SUN, count=12)


# --- the floor against distance -------------------------------------------------------------------------------------

def test_the_best_colour_gives_a_bolometric_correction_near_plus_0_03():
    assert 0.0 < tiers_module.BEST_CORRECTION < 0.1


def test_a_star_of_the_best_colour_and_one_solar_luminosity_shows_at_10_pc_at_the_limit_that_its_magnitude_gives():
    limit = photometry.SUN_ABSOLUTE_BOLOMETRIC_MAG - tiers_module.BEST_CORRECTION
    assert floor_luminosity_sol(10.0, limit) == pytest.approx(1.0)


def test_ten_times_the_distance_raises_the_floor_a_hundredfold_and_dust_by_the_dust_factor():
    near = floor_luminosity_sol(30.0, 7.5)
    assert floor_luminosity_sol(300.0, 7.5) == pytest.approx(100.0 * near)
    assert floor_luminosity_sol(30.0, 7.5, a_v=2.5) == pytest.approx(10.0 * near)
    assert floor_luminosity_sol(30.0, 6.5) == pytest.approx(10 ** 0.4 * near)           # a magnitude brighter limit needs 2.5 times the light


@pytest.mark.parametrize("bad", [0.0, -3.0, math.nan, math.inf])
def test_a_distance_that_is_not_positive_and_finite_raises(bad):
    with pytest.raises(ValueError):
        floor_luminosity_sol(bad, 7.5)


# --- the reach of a band --------------------------------------------------------------------------------------------

def test_without_dust_the_reach_is_the_distance_the_inverse_square_law_gives(clear):
    for luminosity in (50.0, 237.0, 3000.0):
        expected = 10.0 * 10 ** ((7.5 - (photometry.SUN_ABSOLUTE_BOLOMETRIC_MAG - 2.5 * math.log10(luminosity) - tiers_module.BEST_CORRECTION)) / 5.0)
        assert reach_pc(luminosity, 7.5, clear) == pytest.approx(expected, rel=2e-3)


def test_the_reach_matches_the_one_the_backfill_measurement_found_for_the_237_lsun_band(clear):
    assert reach_pc(237.0, 7.5, clear) == pytest.approx(555.0, rel=0.02)           # VIEW.19: 553 pc from sampled stars


def test_dust_shortens_the_reach_and_a_brighter_star_reaches_farther(clear, dusty):
    assert reach_pc(1000.0, 7.5, dusty) < reach_pc(1000.0, 7.5, clear)
    assert reach_pc(100.0, 7.5, dusty) < reach_pc(1000.0, 7.5, dusty)
    assert reach_pc(1000.0, 6.5, dusty) < reach_pc(1000.0, 7.5, dusty)


def test_a_star_too_faint_to_show_even_at_the_first_grid_distance_reaches_nowhere(clear):
    assert reach_pc(1e-9, 7.5, clear) == 0.0


def test_a_bad_luminosity_raises(clear):
    with pytest.raises(ValueError):
        reach_pc(0.0, 7.5, clear)


# --- the cells ------------------------------------------------------------------------------------------------------

def test_the_cell_count_of_a_sphere_agrees_with_the_cells_listed():
    listed = sum(1 for _ in spherical.enumerate_sectors_within_radius((300.0, 200.0, 50.0), 60.0, EDGE))
    assert cells_in_sphere(60.0, EDGE) == pytest.approx(listed, rel=0.01)


# --- the table ------------------------------------------------------------------------------------------------------

def test_the_luminosity_bands_start_at_the_brightest_white_dwarf_and_end_at_the_scatters_floor():
    bands = tiers_module.luminosity_bands()
    assert bands[0][1] == pytest.approx(tiers_module.RING_TOP_SOL, rel=0.2)
    assert bands[-1][1] < tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL <= bands[-1][2]


def test_the_rings_are_the_default_mass_edges_and_every_ring_reaches_as_far_as_a_star_of_50_lsun(clear):
    table = build(SUN, 7.5, EDGE, curve=clear)
    assert [ring.key for ring in table.rings] == [1.0, 2.0, 5.0, 8.0]
    assert table.rings[-1].high == tuning.PHENOMENON_MIN_MASS_SOLAR
    assert all(ring.reach_pc == pytest.approx(reach_pc(50.0, 7.5, clear)) for ring in table.rings)
    assert all(ring.ledger == "mass" for ring in table.rings)


def test_the_bands_taken_fit_the_budget_and_the_rest_are_dropped_from_the_bright_end(clear):
    table = build(SUN, 7.5, EDGE, curve=clear, cell_budget=3.0e8)
    assert sum(tier.cells for tier in table.luminosity) <= 3.0e8
    assert table.dropped and table.luminosity
    assert [tier.key for tier in table.luminosity + table.dropped] == [k for k, _l, _h in tiers_module.luminosity_bands()]
    assert table.luminosity[-1].high <= table.dropped[0].low * 1.0001            # bands run on from each other
    assert table.dropped[0].cells + sum(tier.cells for tier in table.luminosity) > 3.0e8


def test_a_bigger_budget_takes_more_bands_and_a_huge_one_takes_them_all(clear):
    small = build(SUN, 7.5, EDGE, curve=clear, cell_budget=1.0e7)
    large = build(SUN, 7.5, EDGE, curve=clear, cell_budget=1.0e9)
    everything = build(SUN, 7.5, EDGE, curve=clear, cell_budget=1.0e12)
    assert len(small.luminosity) < len(large.luminosity) < len(everything.luminosity)
    assert not everything.dropped
    assert tiers_module.cell_visits(everything) > 5e9                                # VIEW.19: 8e9 uncut at 7.5 with no dust, a thinner galaxy aside


def test_a_6_5_sky_is_cheaper_than_a_7_5_sky_and_dust_is_cheaper_still(clear, dusty):
    uncut = 1.0e13
    clear_75 = tiers_module.cell_visits(build(SUN, 7.5, EDGE, curve=clear, cell_budget=uncut))
    clear_65 = tiers_module.cell_visits(build(SUN, 6.5, EDGE, curve=clear, cell_budget=uncut))
    dusty_75 = tiers_module.cell_visits(build(SUN, 7.5, EDGE, curve=dusty, cell_budget=uncut))
    assert dusty_75 < clear_75 and clear_65 < clear_75


# --- the floor against distance, as rows ------------------------------------------------------------------------------

def test_the_table_changes_ledger_where_the_floor_passes_50_lsun(clear):
    table = build(SUN, 7.5, EDGE, curve=clear)
    boundary = reach_pc(50.0, 7.5, clear)
    rows = floor_table(table, clear, (boundary * 0.97, boundary * 1.03))
    assert [row["ledger"] for row in rows] == ["mass", "luminosity"]
    assert rows[0]["floor_sol"] < 50.0 < rows[1]["floor_sol"]


def test_a_distance_inside_the_cut_is_built_and_one_past_it_is_not(clear):
    table = build(SUN, 7.5, EDGE, curve=clear, cell_budget=3.0e8)
    near, far = floor_table(table, clear, (300.0, 3000.0))
    assert near["built"] and not far["built"]
    assert far["band"] is not None and far["band"] > near["band"]


def test_the_floor_through_the_least_dusty_direction_is_the_one_listed(clear, dusty):
    table = build(SUN, 7.5, EDGE, curve=dusty)
    row = floor_table(table, dusty, (1000.0,))[0]
    j = int(np.searchsorted(dusty.grid_pc, 1000.0))
    assert row["floor_sol"] == pytest.approx(float(floor_luminosity_sol(1000.0, 7.5, float(dusty.a_v[:, j].min()))))
    assert row["floor_sol"] > floor_luminosity_sol(1000.0, 7.5)                    # more than the clear sky's floor
