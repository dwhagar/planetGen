# tests/test_gen_phenomenon_reach.py

"""
GEN.323: the phenomena scatter covers the grid's reach. The part of a shell beyond the outline is drawn object first in
groups of bands (`phenomenon_reach`), the faint shells far out together, and the rows inside the outline are the rows
they were. Each class's count agrees with the band walk's, the bounds are never exceeded, and a plan with reach on
places objects beyond its outline.
"""

import collections
import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import density, remnant_distribution, spherical
from planetgen.generation import object_first, phenomenon_reach as reach, phenomenon_scatter as scatter, run_plan
from tests import test_gen_scatter_tiers as tiers
from tests.test_bright_star_scatter import EXTENTS
from tests.test_phenomenon_scatter import _plan_args

SHAPE, EDGE_PC = tiers.SHAPE, tiers.EDGE_PC
E_VALUE = 4000.0
CUT = 1.0
"""A cut of 1 solar mass scatters every neutron star and black hole."""
PLAN_SEED = bytes(range(16))


def _rows(function, *args, seed=3, **kwargs):
    kwargs.setdefault("min_mass_solar", CUT)
    return list(function(SHAPE, *args, EDGE_PC, E_VALUE, seed, **kwargs))


def _labels(rows):
    return collections.Counter(scatter.class_label(row[3], row[4]) for row in rows)


def _point_pc(row):
    return tuple(value / scatter.MPC_PER_PC for value in row[5:8])


@pytest.mark.parametrize("shell, high_band", [(5, -1), (800, 60), (3000, -1), (20, 5)])
def test_a_shells_rows_beyond_the_outline_are_in_their_own_cells_and_never_inside_the_outline(shell, high_band):
    rows = _rows(reach.scatter_shell, shell, high_band, check=True)
    assert len(rows) > 20
    assert {row[3] for row in rows} <= set(scatter.SCATTERED_KINDS)
    for row in rows:
        assert row[0] == shell and abs(row[1]) > high_band
        assert spherical.sector_address_at(_point_pc(row), EDGE_PC) == tuple(row[:3])
    assert len({row[11] for row in rows}) == len(rows)           # every row has a seed of its own
    assert all(row[12] is None or row[12] >= CUT for row in rows)


def test_the_same_seed_gives_the_same_rows_and_another_seed_others():
    assert _rows(reach.scatter_shell, 20, 5, seed=4) == _rows(reach.scatter_shell, 20, 5, seed=4)
    assert _rows(reach.scatter_shell, 20, 5, seed=4) != _rows(reach.scatter_shell, 20, 5, seed=5)


def test_sectors_already_filled_get_no_row():
    rows = _rows(reach.scatter_shell, 20, 5)
    skipped = {tuple(row[:3]) for row in rows[:25]}
    kept = _rows(reach.scatter_shell, 20, 5, skip_addresses=skipped)
    assert [row for row in rows if tuple(row[:3]) not in skipped] == kept


def test_the_rows_inside_the_outline_are_the_ones_the_band_walk_gives():
    # The reach has streams of its own and `scatter_shell` is untouched by it.
    before = _rows(scatter.scatter_shell, 8, 12)
    _rows(reach.scatter_shell, 8, 5)
    assert before == _rows(scatter.scatter_shell, 8, 12)


def test_the_inside_and_beyond_together_agree_with_the_band_walk_of_the_whole_shell():
    shell, outline = 8, 5
    top = spherical.band_range(shell)[1]
    walk, split = collections.Counter(), collections.Counter()
    for seed in range(4):
        walk.update(_labels(_rows(scatter.scatter_shell, shell, top, seed=seed)))
        split.update(_labels(_rows(scatter.scatter_shell, shell, outline, seed=seed)))
        split.update(_labels(_rows(reach.scatter_shell, shell, outline, seed=seed)))
    big = "stellar black-hole"
    assert walk[big] > 3000 and split[big] == pytest.approx(walk[big], rel=0.05)
    others = sum(walk.values()) - walk[big]
    assert others > 20 and abs((sum(split.values()) - split[big]) - others) < 5 * math.sqrt(others) + 3


@pytest.mark.parametrize("first, last", [(1500, 1540), (30000, 30400)])
def test_the_far_pass_agrees_with_the_shells_one_by_one(first, last):
    far = sum(len(_rows(reach.scatter_far, first, last, seed=seed, check=True)) for seed in range(2))
    each = sum(len(_rows(reach.scatter_shell, shell, -1, seed=seed)) for seed in range(2)
               for shell in range(first, last + 1))
    assert far > 100
    assert far == pytest.approx(each, rel=0.12)


def test_the_far_pass_has_its_own_stream_and_is_reproducible():
    a = _rows(reach.scatter_far, 1500, 1540, seed=6)
    assert a == _rows(reach.scatter_far, 1500, 1540, seed=6) and a != _rows(reach.scatter_far, 1500, 1540, seed=7)
    assert all(1500 <= row[0] <= 1540 for row in a)
    assert list(reach.scatter_far(SHAPE, 100, 99, EDGE_PC, E_VALUE, 1)) == []


def test_every_class_bound_covers_its_weight_at_random_points_of_a_group():
    rng = np.random.default_rng(8)
    for shell, high_band in [(5, -1), (60, 10), (800, 60)]:
        groups = object_first.ReachGroups.for_shell(shell, high_band, EDGE_PC)
        for kind, subtype, _rate in scatter._scattered_rates(CUT):
            bound = reach.class_bounds(SHAPE, groups, kind, subtype)
            for index in range(len(groups)):
                radius = np.cbrt(rng.uniform(groups.r_in[index] ** 3, groups.r_out[index] ** 3, 200))
                cos_colat = rng.uniform(groups.cos_low[index], groups.cos_high[index], 200)
                lon = rng.uniform(0.0, 2.0 * math.pi, 200)
                plane = radius * np.sqrt(1.0 - cos_colat ** 2)
                weight = reach._weights(kind, subtype, plane * np.cos(lon), plane * np.sin(lon), radius * cos_colat,
                                        SHAPE)
                assert np.all(weight <= bound[index] * (1.0 + 1e-9)), (shell, kind, subtype, index)


def test_the_remnants_are_thicker_than_the_stars_and_the_young_class_follows_the_young_stars():
    plane, high = (3000.0, 0.0, 0.0), (3000.0, 0.0, 600.0)
    for point in (plane, high):
        x, y, z = (np.array([value]) for value in point)
        stars = density.total_density_array(x, y, z, SHAPE)
        assert reach._weights("black-hole", "stellar", x, y, z, SHAPE)[0] == pytest.approx(
            stars[0] * remnant_distribution.placement_factor("black-hole", point, SHAPE.disk_scale_height_pc))
        assert reach._weights("planetary-nebula", None, x, y, z, SHAPE)[0] == pytest.approx(stars[0])
        assert reach._weights("neutron-star", scatter.MAGNETAR_SUBTYPE, x, y, z, SHAPE)[0] == pytest.approx(
            max(density.population_densities_array(x, y, z, SHAPE)["young"][0], 0.0))


def test_the_expected_objects_and_the_far_screen_follow_the_rates():
    fractions = reach.scattered_fractions(CUT)
    assert len(fractions) == 1 and all(value >= 0.0 for value in fractions[0].values())
    estimates = reach.reach_estimates(SHAPE, [(5, 2), (6, -1)], EDGE_PC, E_VALUE, CUT)
    assert estimates[5][0] > 0.0 and estimates[6][0] > 0.0 and estimates[6][1] > 0
    assert reach.far_expected(SHAPE, 1500, 1540, EDGE_PC, E_VALUE, CUT) > 0.0
    assert reach.far_expected(SHAPE, 100, 99, EDGE_PC, E_VALUE, CUT) == 0.0
    # (the seed's drawn halo holds 8.8% of the stars, so the far run starts past shell 40,000)
    first, means = reach.far_shells(SHAPE, 5, 80000, EDGE_PC, E_VALUE, CUT)
    assert 5 < first <= 80000 and len(means) == first - 5 and means[-1] > 0.0
    # A thinner galaxy sends the far run nearer in.
    assert reach.far_shells(SHAPE, 5, 40000, EDGE_PC, E_VALUE / 1000.0, CUT)[0] < first


@pytest.mark.db
def test_a_plan_scatter_covers_the_reach_and_leaves_the_outline_as_it_was(mysql_config, monkeypatch):
    from tests.test_bright_star_scatter import E_VALUE as e_value, EDGE_PC as edge_pc, SHAPE as shape

    def run(covers):
        monkeypatch.setattr(tuning, "SCATTER_COVERS_REACH", covers)
        store.save_galaxy_shape(shape, edge_pc=edge_pc, outer_shell_index=8, expected_system_count_at_density_1=e_value,
                                config=mysql_config, cell_star_ceiling=1e12, reach_shell_index=12,
                                galaxy_seed=PLAN_SEED)
        store.replace_galaxy_shells(EXTENTS, config=mysql_config)
        summary = run_plan.scatter_phenomena(_plan_args(mysql_config, "--workers", "1"))
        conn = store.get_connection(mysql_config)
        try:
            return summary, conn.execute(
                "SELECT shell_index, band_index, wedge_index, kind, position_x_mpc, position_y_mpc, position_z_mpc,"
                " seed FROM phenomenon_scatter ORDER BY seed").fetchall()
        finally:
            conn.close()

    plain_summary, plain = run(False)
    full_summary, full = run(True)
    high = dict(EXTENTS)
    inside = lambda row: row["shell_index"] in high and abs(row["band_index"]) <= high[row["shell_index"]]
    assert plain and all(row["shell_index"] in high for row in plain)
    # The hypervelocity stars are drawn to the reach (GEN.238), so a wider reach draws more of them, each at a distance
    # uniform to it: every other kind keeps its rows.
    others = lambda rows: [row for row in rows if row["kind"] != "hypervelocity-star"]
    assert others(row for row in full if inside(row)) == others(row for row in plain if inside(row))
    beyond = [row for row in full if not inside(row)]
    assert beyond and max(row["shell_index"] for row in beyond) > 8
    assert all(row["shell_index"] <= 12 for row in full)
    hypervelocity = lambda rows: len([row for row in rows if row["kind"] == "hypervelocity-star"])
    assert full_summary["total"] - hypervelocity(full) == plain_summary["total"] - hypervelocity(plain) + len(others(beyond))
    assert hypervelocity(full) > hypervelocity(plain)
