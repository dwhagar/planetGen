# tests/test_gen_reach_groups.py

"""
GEN.322: the bright-star and mass passes cover the grid's reach. The part of a shell beyond the outline is cut into
groups of 64 bands (`ReachGroups.for_shell`), the faint shells far out are drawn together (`ReachGroups.far`,
`scatter_far`), and the stars inside the outline are the rows they were without any of it.
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import scatter_limits, spherical
from planetgen.galaxy.scatter_limits import ScatterLimits
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import object_first, run_plan
from tests import test_gen_scatter_tiers as tiers
from tests.test_bright_star_scatter import EXTENTS, _plan_args

SHAPE, EDGE_PC, E_VALUE, FLOOR, MASS = tiers.SHAPE, tiers.EDGE_PC, tiers.E_VALUE, tiers.FLOOR, tiers.MASS
LIMITS = ScatterLimits.of(tiers._skeleton(), FLOOR, MASS)
TIERED = scatter_limits.luminosity_pass(LIMITS)
BY_TIER = brightStars.tier_fractions_of(TIERED)
SHELL, OUTLINE = 800, 60
"""A shell the toy halo fills and an outline that reaches 60 bands either side of its equator."""


def _outline_volume(shell, high_band):
    """The volume, in edge cubes, of the bands `-high_band .. high_band` of a shell."""
    return sum(spherical.cell_volume(shell, band, 0, EDGE_PC) * spherical.wedge_count(shell, band)
               for band in range(-high_band, high_band + 1)) / EDGE_PC ** 3


def _shell_volume(shell):
    inner, outer = spherical.shell_bounds_pc(shell, EDGE_PC)
    return (outer ** 3 - inner ** 3) * 4.0 * math.pi / 3.0 / EDGE_PC ** 3


@pytest.mark.parametrize("shell, high_band", [(1, 0), (2, 0), (3, 1), (5, 2), (100, 40), (800, 60), (1000, 300), (4000, 0)])
def test_the_groups_and_the_outline_bands_tile_a_shell(shell, high_band):
    groups = object_first.ReachGroups.for_shell(shell, high_band, EDGE_PC)
    assert groups.volume.sum() + _outline_volume(shell, high_band) == pytest.approx(_shell_volume(shell), rel=1e-12)


@pytest.mark.parametrize("shell", [1, 2, 9, 800, 5000])
def test_a_shell_the_outline_does_not_reach_is_one_run_of_groups_from_the_equator(shell):
    groups = object_first.ReachGroups.for_shell(shell, -1, EDGE_PC)
    assert groups.volume.sum() == pytest.approx(_shell_volume(shell), rel=1e-12)
    layout = spherical.shell_layout(shell)
    count = len(groups) // 2
    assert count >= 1 + (layout.high_band + 1) // tuning.OBJECT_FIRST_GROUP_BANDS
    # North then south mirror images, no group wider than the group size and the last one reaching the pole.
    assert groups.colat_near[count - 1] == 0.0 and groups.colat_near[2 * count - 1] == pytest.approx(math.pi)


def test_a_shell_the_outline_reaches_the_pole_of_has_no_groups():
    top = spherical.band_range(100)[1]
    assert len(object_first.ReachGroups.for_shell(100, top, EDGE_PC)) == 0
    assert len(object_first.ReachGroups.for_shell(0, 0, EDGE_PC)) == 0


def test_a_group_holds_at_most_the_group_size_in_bands():
    groups = object_first.ReachGroups.for_shell(1000, 300, EDGE_PC)
    layout = spherical.shell_layout(1000)
    steps = (groups.colat_far[:len(groups) // 2] - groups.colat_near[:len(groups) // 2]) / layout.band_step
    assert steps.max() <= tuning.OBJECT_FIRST_GROUP_BANDS + 1.0 + 1e-9    # the last takes the cap too


def test_the_far_runs_tile_the_shells_asked_for_and_stay_thin():
    blocks = object_first.ReachGroups.far_blocks(9000, 65535)
    assert blocks[0][0] == 9000 and blocks[-1][1] == 65535
    assert all(second[0] == first[1] + 1 for first, second in zip(blocks, blocks[1:]))
    assert all(end - start + 1 <= max(1, start // tuning.OBJECT_FIRST_FAR_STEPS) for start, end in blocks)
    assert len(blocks) < 1000
    groups = object_first.ReachGroups.far(9000, 65535, EDGE_PC)
    assert groups.volume.sum() == pytest.approx(
        (spherical.shell_bounds_pc(65535, EDGE_PC)[1] ** 3 - spherical.shell_bounds_pc(9000, EDGE_PC)[0] ** 3)
        * 4.0 * math.pi / 3.0 / EDGE_PC ** 3, rel=1e-12)


def _points_in(groups, index, rng, count):
    """Points uniform in group `index`: `(plane, z)` arrays."""
    radius = rng.uniform(groups.r_in[index] ** 3, groups.r_out[index] ** 3, count) ** (1.0 / 3.0)
    cos_colat = rng.uniform(groups.cos_low[index], groups.cos_high[index], count)
    return radius * np.sqrt(1.0 - cos_colat ** 2), radius * cos_colat


@pytest.mark.parametrize("make", [
    lambda: object_first.ReachGroups.for_shell(SHELL, OUTLINE, EDGE_PC),
    lambda: object_first.ReachGroups.for_shell(SHELL, -1, EDGE_PC),
    lambda: object_first.ReachGroups.far(3000, 6000, EDGE_PC),
])
def test_a_groups_bound_is_never_below_the_density_inside_it(make):
    groups = make()
    fractions = brightStars.band_fractions(FLOOR)
    bound = object_first.extent_majorants(SHAPE, groups.inner, groups.z_near, groups.z_far, fractions)
    rng = np.random.default_rng(7)
    for index in rng.choice(len(groups), size=min(len(groups), 40), replace=False):
        plane, z = _points_in(groups, index, rng, 60)
        angle = rng.uniform(0.0, 2.0 * math.pi, 60)
        found = brightStars.population_densities_array(plane * np.cos(angle), plane * np.sin(angle), z, SHAPE)
        wanted = sum(np.maximum(found[name], 0.0) * fractions[name] for name in brightStars.POPULATIONS)
        assert (wanted <= bound[index] * (1.0 + 1e-9)).all(), index


def test_the_stars_inside_the_outline_are_the_rows_they_were_without_the_reach():
    kwargs = dict(tiered=TIERED, check=True)
    inside = list(brightStars.scatter_shell(SHAPE, SHELL, OUTLINE, EDGE_PC, E_VALUE, FLOOR, "reach", **kwargs))
    both = list(brightStars.scatter_shell(SHAPE, SHELL, OUTLINE, EDGE_PC, E_VALUE, FLOOR, "reach", reach=True, **kwargs))
    assert inside and len(both) > len(inside)
    assert both[:len(inside)] == inside
    column = {name: index for index, name in enumerate(store.BRIGHT_STAR_COLUMNS)}
    beyond = both[len(inside):]
    assert all(row[column["shell_index"]] == SHELL and abs(row[column["band_index"]]) > OUTLINE for row in beyond)


def test_a_shell_is_drawn_the_same_every_time_and_with_a_seed_of_its_own():
    first = list(brightStars.scatter_shell(SHAPE, 1500, -1, EDGE_PC, E_VALUE, FLOOR, "a", tiered=TIERED, reach=True))
    again = list(brightStars.scatter_shell(SHAPE, 1500, -1, EDGE_PC, E_VALUE, FLOOR, "a", tiered=TIERED, reach=True))
    other = list(brightStars.scatter_shell(SHAPE, 1500, -1, EDGE_PC, E_VALUE, FLOOR, "b", tiered=TIERED, reach=True))
    assert first and first == again and first != other


def test_a_star_is_never_drawn_in_a_skipped_sector_beyond_the_outline():
    rows = list(brightStars.scatter_shell(SHAPE, 1500, -1, EDGE_PC, E_VALUE, FLOOR, "a", tiered=TIERED, reach=True))
    column = {name: index for index, name in enumerate(store.BRIGHT_STAR_COLUMNS)}
    skip = {(row[column["shell_index"]], row[column["band_index"]], row[column["wedge_index"]]) for row in rows[::2]}
    kept = list(brightStars.scatter_shell(SHAPE, 1500, -1, EDGE_PC, E_VALUE, FLOOR, "a", tiered=TIERED, reach=True,
                                          skip_addresses=skip))
    addresses = {(row[column["shell_index"]], row[column["band_index"]], row[column["wedge_index"]]) for row in kept}
    assert kept and not addresses & skip


def test_the_far_pass_draws_only_the_shells_it_was_given_and_the_same_every_time():
    rows = list(brightStars.scatter_far(SHAPE, 3000, 3400, EDGE_PC, E_VALUE, FLOOR, "far", tiered=TIERED, check=True))
    again = list(brightStars.scatter_far(SHAPE, 3000, 3400, EDGE_PC, E_VALUE, FLOOR, "far", tiered=TIERED))
    column = store.BRIGHT_STAR_COLUMNS.index("shell_index")
    assert rows and rows == again
    assert all(3000 <= row[column] <= 3400 for row in rows)


def test_the_far_pass_expects_what_the_shell_by_shell_pass_does():
    first, last, seeds = 1500, 1507, range(4)
    far = sum(len(list(brightStars.scatter_far(SHAPE, first, last, EDGE_PC, E_VALUE, FLOOR, f"s{seed}", tiered=TIERED)))
              for seed in seeds)
    each = sum(len(list(brightStars.scatter_shell(SHAPE, shell, -1, EDGE_PC, E_VALUE, FLOOR, f"s{seed}", tiered=TIERED,
                                                  reach=True)))
               for seed in seeds for shell in range(first, last + 1))
    assert far > 500 and each > 500
    assert far == pytest.approx(each, rel=0.12)


def test_far_shells_start_where_no_later_shell_can_expect_a_few_stars():
    first_far, means = brightStars.far_shells(SHAPE, 50, 30000, EDGE_PC, E_VALUE, BY_TIER)
    assert 50 < first_far <= 30000
    assert len(means) == first_far - 50
    assert means[-1] > 0.0
    later, _means = brightStars.far_shells(SHAPE, first_far, 30000, EDGE_PC, E_VALUE, BY_TIER)
    assert later == first_far
    # Nothing to screen, and a screen that finds every shell far.
    assert brightStars.far_shells(SHAPE, 100, 99, EDGE_PC, E_VALUE, BY_TIER)[0] == 100
    assert brightStars.far_shells(SHAPE, 29000, 30000, EDGE_PC, E_VALUE, BY_TIER)[0] == 29000


def test_the_expected_count_of_a_shell_includes_the_part_beyond_the_outline():
    inside = brightStars.shell_weight(SHAPE, SHELL, OUTLINE, EDGE_PC, E_VALUE, BY_TIER)
    beyond = brightStars.reach_expected_stars(SHAPE, SHELL, OUTLINE, EDGE_PC, E_VALUE, BY_TIER)
    both = brightStars.shell_weight(SHAPE, SHELL, OUTLINE, EDGE_PC, E_VALUE, BY_TIER, beyond=beyond)
    assert beyond[0] > 0.0 and beyond[1] > 0
    assert both[1] == pytest.approx(inside[1] + beyond[0]) and both[0] > inside[0]
    estimates = brightStars.reach_estimates(SHAPE, [(SHELL + k, OUTLINE) for k in range(40)], EDGE_PC, E_VALUE, BY_TIER)
    assert len(estimates) == 40 and estimates[SHELL + 17] == estimates[SHELL + 16]


@pytest.mark.db
def test_a_plan_scatter_covers_the_reach_and_a_pre_placed_star_is_never_drawn_twice(mysql_config, monkeypatch):
    monkeypatch.setattr(tuning, "SCATTER_COVERS_REACH", True)
    from tests.test_bright_star_scatter import E_VALUE as e_value, EDGE_PC as edge_pc, SHAPE as shape
    store.save_galaxy_shape(shape, edge_pc=edge_pc, outer_shell_index=8, expected_system_count_at_density_1=e_value,
                            config=mysql_config, cell_star_ceiling=1e12, reach_shell_index=12)
    store.replace_galaxy_shells(EXTENTS, config=mysql_config)
    summary = run_plan.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only",
                                                       "--bright-star-min-luminosity", "30000"))
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT shell_index, band_index, COUNT(*) AS n FROM bright_stars"
                            " GROUP BY shell_index, band_index").fetchall()
    finally:
        conn.close()
    high = dict(EXTENTS)
    assert summary["total"] == sum(row["n"] for row in rows) > 0
    assert max(row["shell_index"] for row in rows) > 8
    assert all(row["shell_index"] <= 12 for row in rows)
    assert any(row["shell_index"] in high and abs(row["band_index"]) > high[row["shell_index"]] for row in rows)
