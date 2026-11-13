"""VIEW.24, VIEW.41: the sky tier table (`planetgen/generation/sky_tiers.py`): the floor against distance, the reach of a band through the dust, the cut by a
budget of block visits, the rings under one solar mass, and the split between the luminosity and the mass ledger."""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import spherical
from planetgen.generation import bright_stars as brightStars
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


@pytest.fixture
def few_samples(monkeypatch):
    """A ring's visible share is found from a sample; a small one keeps the tests quick."""
    monkeypatch.setattr(tiers_module, "RING_SAMPLES", 40)


def test_the_rings_of_one_solar_mass_and_up_reach_as_far_as_their_yield_and_no_farther_than_a_star_of_50_lsun(clear, few_samples):
    """GEN.1486: the heavy rings follow the yield rule of the light ones. The ring of 8 solar masses and up holds no star under 50 Lsun, so it draws nothing and is left out."""
    table = build(SUN, 7.5, EDGE, curve=clear)
    heavy = [ring for ring in table.rings if ring.key >= 1.0]
    assert [ring.key for ring in heavy] == [1.0, 2.0, 5.0]
    ceiling = reach_pc(50.0, 7.5, clear)
    assert all(0.0 < ring.reach_pc <= ceiling for ring in heavy)
    assert heavy[0].reach_pc < ceiling and heavy[1].reach_pc < ceiling                  # a ring reaches as far as 0.5 percent of its stars show, which is short of a 50 Lsun star
    assert heavy[-1].reach_pc < heavy[1].reach_pc                                       # the 5 to 8 ring holds few dim stars and they are the brighter ones
    assert all(ring.ledger == "mass" and ring.visits == ring.cells for ring in table.rings)


def test_the_yield_is_one_half_of_one_percent_for_every_ring():
    assert tiers_module.RING_MIN_YIELD == 0.005


def test_a_ring_with_no_star_under_the_cap_has_an_empty_sample_and_no_reach(clear, few_samples):
    top = tiers_module.ring_edges()[-2:]
    samples = tiers_module.ring_luminosities(top[0], top[1], 40)
    assert all(len(values) == 0 for values in samples.values())
    assert tiers_module.visible_share(samples, 1.0) == 0.0
    assert tiers_module.ring_reach_pc(top[0], top[1], 7.5, clear) == 0.0


def test_the_rings_under_one_solar_mass_are_the_near_masses_and_reach_less_far(clear, few_samples):
    table = build(SUN, 7.5, EDGE, curve=clear)
    near = [ring for ring in table.rings if ring.key < 1.0]
    assert near and [ring.key for ring in near] == sorted(ring.key for ring in near)
    assert all(ring.key in tiers_module.NEAR_RING_MASSES_SOL for ring in near)
    assert near[-1].high == 1.0
    ceiling = reach_pc(50.0, 7.5, clear)
    assert all(0.0 < ring.reach_pc < ceiling for ring in near)
    assert near[-1].reach_pc < 0.5 * ceiling                        # tens of parsecs, not the 250 of a ring of 1 solar mass and up
    assert [ring.reach_pc for ring in near] == sorted(ring.reach_pc for ring in near)      # the heavier the ring, the farther it shows


def test_a_ring_that_shows_nowhere_is_left_out_and_without_near_masses_the_rings_are_the_scatters(clear, few_samples):
    table = build(SUN, 7.5, EDGE, curve=clear, edges=tiers_module.ring_edges(()))
    assert [ring.key for ring in table.rings] == [1.0, 2.0, 5.0]
    dim = build(SUN, 1.0, EDGE, curve=clear)
    assert all(ring.key >= 1.0 or ring.reach_pc < 30.0 for ring in dim.rings)           # a 1.0 sky shows almost no dwarf


def test_the_share_of_a_ring_that_shows_falls_with_the_floor(few_samples):
    luminosities = tiers_module.ring_luminosities(0.5, 0.7, 40)
    shares = [tiers_module.visible_share(luminosities, floor) for floor in (1e-4, 1e-2, 1e-1, 1.0, 50.0)]
    assert shares[0] == 1.0 and shares == sorted(shares, reverse=True) and shares[-1] < 0.5
    assert tiers_module.ring_luminosities(0.5, 0.7, 40) is luminosities                  # kept: the same ring gives the same sample


def test_the_weights_of_the_populations_change_the_share_and_must_name_one(few_samples):
    luminosities = tiers_module.ring_luminosities(0.7, 0.85, 40)
    halo = tiers_module.visible_share(luminosities, 2.0, {"halo": 1.0})
    old = tiers_module.visible_share(luminosities, 2.0, {"old": 1.0})
    assert halo > old                                                                   # the halo's old giants are bright, the disk's are not yet
    both = tiers_module.visible_share(luminosities, 2.0, {"halo": 1.0, "old": 1.0})
    assert both == pytest.approx((halo + old) / 2.0)
    with pytest.raises(ValueError):
        tiers_module.visible_share(luminosities, 2.0, {"nowhere": 1.0})


def test_a_stricter_yield_shortens_a_ring(clear, few_samples):
    loose = tiers_module.ring_reach_pc(0.7, 0.85, 7.5, clear, min_yield=0.01)
    strict = tiers_module.ring_reach_pc(0.7, 0.85, 7.5, clear, min_yield=0.5)
    assert 0.0 <= strict < loose


def test_the_rows_of_a_sphere_are_its_cross_section_in_edges():
    for radius in (60.0, 250.0):
        rows = sum(1 for _ in spherical.sphere_rows((4000.0, 0.0, 0.0), radius, EDGE))
        assert rows == pytest.approx(math.pi * (radius / EDGE) ** 2, rel=0.02)


def test_the_visits_of_a_sphere_are_its_rows_and_its_blocks():
    radius = 250.0
    rows = sum(1 for _ in spherical.sphere_rows((4000.0, 0.0, 0.0), radius, EDGE))
    blocks = cells_in_sphere(radius, EDGE) / brightStars.BLOCK_WEDGES
    assert tiers_module.block_visits(radius, EDGE) == pytest.approx(rows + blocks, rel=0.02)
    assert tiers_module.block_visits(3000.0, EDGE) > 100 * tiers_module.block_visits(250.0, EDGE)       # the far spheres cost the most


def _lum_visits(table):
    return sum(tier.visits for tier in table.luminosity)


def test_the_bands_taken_fit_the_budget_and_the_rest_are_dropped_from_the_bright_end(clear):
    total = _lum_visits(build(SUN, 7.5, EDGE, curve=clear, visit_budget=1.0e15))
    budget = total / 3.0
    table = build(SUN, 7.5, EDGE, curve=clear, visit_budget=budget)
    assert _lum_visits(table) <= budget
    assert table.dropped and table.luminosity
    assert [tier.key for tier in table.luminosity + table.dropped] == [k for k, _l, _h in tiers_module.luminosity_bands()]
    assert table.luminosity[-1].high <= table.dropped[0].low * 1.0001            # bands run on from each other
    assert table.dropped[0].visits + _lum_visits(table) > budget


def test_a_bigger_budget_takes_more_bands_and_a_huge_one_takes_them_all(clear):
    everything = build(SUN, 7.5, EDGE, curve=clear, visit_budget=1.0e15)
    total = _lum_visits(everything)
    small = build(SUN, 7.5, EDGE, curve=clear, visit_budget=total / 100.0)
    large = build(SUN, 7.5, EDGE, curve=clear, visit_budget=total / 3.0)
    assert len(small.luminosity) < len(large.luminosity) < len(everything.luminosity)
    assert not everything.dropped
    assert 1e6 < total < 1e9                  # the whole luminosity side of a 7.5 sky is millions of block visits, VIEW.19 had 8e9 cell visits


def test_the_default_budget_takes_a_dust_free_7_5_sky_whole(clear):
    table = build(SUN, 7.5, EDGE, curve=clear)
    assert not table.dropped and len(table.luminosity) == len(tiers_module.luminosity_bands())


def test_a_6_5_sky_is_cheaper_than_a_7_5_sky_and_dust_is_cheaper_still(clear, dusty):
    uncut = 1.0e15
    clear_75 = _lum_visits(build(SUN, 7.5, EDGE, curve=clear, visit_budget=uncut))
    clear_65 = _lum_visits(build(SUN, 6.5, EDGE, curve=clear, visit_budget=uncut))
    dusty_75 = _lum_visits(build(SUN, 7.5, EDGE, curve=dusty, visit_budget=uncut))
    assert dusty_75 < clear_75 and clear_65 < clear_75


# --- the floor against distance, as rows ------------------------------------------------------------------------------

def test_the_table_changes_ledger_where_the_floor_passes_50_lsun(clear):
    table = build(SUN, 7.5, EDGE, curve=clear)
    boundary = reach_pc(50.0, 7.5, clear)
    rows = floor_table(table, clear, (boundary * 0.97, boundary * 1.03))
    assert [row["ledger"] for row in rows] == ["mass", "luminosity"]
    assert rows[0]["floor_sol"] < 50.0 < rows[1]["floor_sol"]


def test_a_distance_inside_the_cut_is_built_and_one_past_it_is_not(clear):
    table = build(SUN, 7.5, EDGE, curve=clear, visit_budget=_lum_visits(build(SUN, 7.5, EDGE, curve=clear, visit_budget=1.0e15)) / 3.0)
    near, far = floor_table(table, clear, (300.0, 3000.0))
    assert near["built"] and not far["built"]
    assert far["band"] is not None and far["band"] > near["band"]


def test_the_floor_through_the_least_dusty_direction_is_the_one_listed(clear, dusty):
    table = build(SUN, 7.5, EDGE, curve=dusty)
    row = floor_table(table, dusty, (1000.0,))[0]
    j = int(np.searchsorted(dusty.grid_pc, 1000.0))
    assert row["floor_sol"] == pytest.approx(float(floor_luminosity_sol(1000.0, 7.5, float(dusty.a_v[:, j].min()))))
    assert row["floor_sol"] > floor_luminosity_sol(1000.0, 7.5)                    # more than the clear sky's floor
