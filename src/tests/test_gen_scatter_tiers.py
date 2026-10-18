# tests/test_gen_scatter_tiers.py

"""
GEN.333: the scatter's limits by place. A band takes the core, baseline or halo tier from the dominant term of its
density, and the scatter's passes, the marks, the clears, the backfill and the fill read that tier's floor and mass
limit from the address.
"""

import collections
import math
from types import SimpleNamespace

import numpy as np
import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import scatter_limits, structure_draw
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.scatter_limits import BASELINE, CORE, HALO, ScatterLimits
from planetgen.galaxy.spherical import band_range
from planetgen.generation import bright_stars as brightStars
from planetgen.physics import constants
from planetgen.physics.units import ly_to_pc

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
E_VALUE = 1.0e5
FLOOR = 9000.0
MASS = 14.0


def _shape(nuclear_amplitude=1.0):
    return build_galaxy_shape(
        disk_scale_length_pc=300.0, disk_scale_height_pc=60.0, bulge_scale_radius_pc=40.0, bulge_amplitude=2.0,
        arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4, nuclear_amplitude=nuclear_amplitude,
    )


# A seed whose drawn structure keeps the toy galaxy's three shells in the tiers they are named for (core, mixed, disc).
SEED = bytes([39]) * 16
SHAPE = _shape()._replace(structure=structure_draw.draw_structure(SEED))   # a stored shape's seed draws its structure
CORE_SHELL, MIXED_SHELL, DISC_SHELL = 1, 800, 40
"""A shell in the nuclear cluster, one the halo reaches near its poles and the disc crosses, and one in the disc."""


def _skeleton(shape=SHAPE, core_floor=None, core_mass=None, halo_floor=tuning.HALO_MIN_LUMINOSITY_SOL):
    return SimpleNamespace(shape=shape, edge_pc=EDGE_PC, core_min_luminosity_sol=core_floor,
                           core_mass_limit_sol=core_mass, halo_min_luminosity_sol=halo_floor)


def _high(shell):
    return min(band_range(shell)[1], 80)


def test_a_band_takes_the_tier_of_the_term_that_dominates_it():
    assert set(scatter_limits.tier_indices(SHAPE, CORE_SHELL, _high(CORE_SHELL), EDGE_PC).tolist()) == {CORE}
    assert set(scatter_limits.tier_indices(SHAPE, DISC_SHELL, _high(DISC_SHELL), EDGE_PC).tolist()) == {BASELINE}
    mixed = scatter_limits.tier_indices(SHAPE, MIXED_SHELL, _high(MIXED_SHELL), EDGE_PC)
    assert set(mixed.tolist()) == {BASELINE, HALO}
    # The halo takes the bands far from the plane, the baseline the ones the disc still reaches.
    high = _high(MIXED_SHELL)
    assert mixed[0] == mixed[-1] == HALO and mixed[high] == BASELINE


def test_with_no_nuclear_structures_there_is_no_core_tier():
    flat = _shape(nuclear_amplitude=0.0)
    for shell in (0, 1, 2, 5):
        assert CORE not in scatter_limits.tier_indices(flat, shell, _high(shell), EDGE_PC).tolist()


def test_the_tier_is_a_pure_function_of_the_address_whichever_way_it_is_asked():
    high = _high(MIXED_SHELL)
    whole = scatter_limits.tier_indices(SHAPE, MIXED_SHELL, high, EDGE_PC)
    for band in (-high, -high // 2, 0, 3, high):
        assert scatter_limits.tier_index(SHAPE, MIXED_SHELL, band, EDGE_PC) == whole[band + high]
    # A shorter run of the same shell gives the same tiers for the bands it holds.
    shorter = scatter_limits.tier_indices(SHAPE, MIXED_SHELL, high // 2, EDGE_PC)
    assert shorter.tolist() == whole[high - high // 2: high + high // 2 + 1].tolist()


def test_the_weights_choose_the_same_tier_as_the_scatter_does():
    """`sampled_band_expectations` reads the tier from the densities it samples; it must agree with `tier_indices`."""
    for shell in (CORE_SHELL, DISC_SHELL, MIXED_SHELL):
        high = _high(shell)
        whole = scatter_limits.tier_indices(SHAPE, shell, high, EDGE_PC)
        for band in range(-high, high + 1, max(1, high // 7)):
            from planetgen.galaxy.spherical import band_center_pc, wedge_count
            bins = min(wedge_count(shell, band), brightStars.WEIGHT_ANGLE_BINS)
            plane, z = band_center_pc(shell, band, EDGE_PC)
            points = [(plane * math.cos(t), plane * math.sin(t))
                      for t in ((j + 0.5) * 2 * math.pi / bins for j in range(bins))]
            found = brightStars._densities_at(points, z, SHAPE)
            assert scatter_limits.tier_of_densities(found) == whole[band + high]


def test_the_tiers_limits_follow_the_baseline_and_the_settings():
    limits = ScatterLimits.of(_skeleton(), FLOOR, MASS)
    core, baseline, halo = limits.tiers
    assert core.min_luminosity_sol == pytest.approx(FLOOR * tuning.CORE_TIER_FACTOR)
    assert core.mass_limit_sol == pytest.approx(MASS * tuning.CORE_TIER_FACTOR)
    assert (baseline.min_luminosity_sol, baseline.mass_limit_sol) == (FLOOR, MASS)
    assert (halo.min_luminosity_sol, halo.mass_limit_sol) == (tuning.HALO_MIN_LUMINOSITY_SOL, MASS)
    set_core = ScatterLimits.of(_skeleton(core_floor=20000.0, core_mass=30.0, halo_floor=80.0), FLOOR, MASS)
    assert set_core.tiers[CORE] == scatter_limits.TierLimit(20000.0, 30.0)
    assert set_core.tiers[HALO].min_luminosity_sol == 80.0
    nothing = ScatterLimits.of(_skeleton(), FLOOR, None)
    assert [found.mass_limit_sol for found in nothing.tiers] == [None, None, None]


def test_before_any_scatter_no_band_has_a_limit():
    limits = ScatterLimits.of(_skeleton(), None, None)
    assert limits.tiers is None
    assert limits.level_at(CORE_SHELL, 0) is None and limits.mass_limit_at(CORE_SHELL, 0) is None


def test_a_band_reads_its_own_tiers_limits_from_the_address():
    limits = ScatterLimits.of(_skeleton(), FLOOR, MASS)
    assert limits.level_at(CORE_SHELL, 0) == pytest.approx(FLOOR * 1.1)
    assert limits.mass_limit_at(CORE_SHELL, 0) == pytest.approx(MASS * 1.1)
    assert limits.level_at(DISC_SHELL, 0) == FLOOR
    high = _high(MIXED_SHELL)
    assert limits.level_at(MIXED_SHELL, 0) == FLOOR
    assert limits.level_at(MIXED_SHELL, -high) == tuning.HALO_MIN_LUMINOSITY_SOL
    # A place with no band takes the baseline's.
    assert limits.level_at(None, None) == FLOOR and limits.mass_limit_at(None, None) == MASS


def _drawn(shell, tiered, seed="tiers", **kwargs):
    high = _high(shell)
    rows = list(brightStars.scatter_shell(SHAPE, shell, high, EDGE_PC, E_VALUE, FLOOR, seed, tiered=tiered,
                                          check=True, **kwargs))
    tiers = scatter_limits.tier_indices(SHAPE, shell, high, EDGE_PC)
    by_tier = collections.defaultdict(list)
    for row in rows:
        by_tier[int(tiers[row[1] + high])].append(row)
    return by_tier


COLUMN = {name: index for index, name in enumerate(store.BRIGHT_STAR_COLUMNS)}


def _watts(row):
    return row[COLUMN["luminosity_w"]]


def test_the_luminosity_pass_draws_each_band_down_to_its_tiers_floor():
    limits = ScatterLimits.of(_skeleton(), FLOOR, MASS)
    tiered = scatter_limits.luminosity_pass(limits)
    core = _drawn(CORE_SHELL, tiered)
    assert core[CORE] and set(core) == {CORE}
    lowest = min(_watts(row) for row in core[CORE]) / constants.SOLAR_LUMINOSITY
    assert lowest >= FLOOR * tuning.CORE_TIER_FACTOR * (1 - 1e-9)
    assert max(row[COLUMN["initial_mass_sol"]] for row in core[CORE]) < MASS * tuning.CORE_TIER_FACTOR
    mixed = _drawn(MIXED_SHELL, tiered)
    assert len(mixed[HALO]) > 3 * len(mixed[BASELINE])      # the halo's low floor takes many dimmer giants (red giants stop at 2,500 solar, GEN.286)
    halo_lowest = min(_watts(row) for row in mixed[HALO]) / constants.SOLAR_LUMINOSITY
    assert tuning.HALO_MIN_LUMINOSITY_SOL * (1 - 1e-9) <= halo_lowest < FLOOR
    # The baseline bands of this shell hold about one star that bright, so several draws are looked at.
    baseline = [row for seed in ("tiers", "tiers-b", "tiers-c", "tiers-d", "tiers-e", "tiers-f")
                for row in _drawn(MIXED_SHELL, tiered, seed=seed)[BASELINE]]
    assert baseline
    assert min(_watts(row) for row in baseline) / constants.SOLAR_LUMINOSITY >= FLOOR * (1 - 1e-9)


def test_a_halo_floor_of_fifty_draws_down_to_fifty():
    limits = ScatterLimits.of(_skeleton(halo_floor=50.0), FLOOR, MASS)
    mixed = _drawn(MIXED_SHELL, scatter_limits.luminosity_pass(limits))
    lowest = min(_watts(row) for row in mixed[HALO]) / constants.SOLAR_LUMINOSITY
    assert 50.0 * (1 - 1e-9) <= lowest < 500.0
    assert tuning.WD_LUMINOSITY_RANGE_SOL[1] <= 50.0


def test_the_mass_pass_draws_each_band_from_its_tiers_mass_limit():
    limits = ScatterLimits.of(_skeleton(), FLOOR, MASS)
    tiered = scatter_limits.mass_pass(limits, 100.0)
    core = _drawn(CORE_SHELL, tiered)
    assert core[CORE]
    assert min(row[COLUMN["initial_mass_sol"]] for row in core[CORE]) >= MASS * tuning.CORE_TIER_FACTOR
    disc = _drawn(DISC_SHELL, tiered)
    assert disc[BASELINE]
    assert min(row[COLUMN["initial_mass_sol"]] for row in disc[BASELINE]) >= MASS


def test_with_the_same_limits_in_every_tier_the_pass_is_the_untiered_one():
    same = scatter_limits.TierLimit(FLOOR, MASS)
    limits = ScatterLimits(SHAPE, EDGE_PC, (same, same, same))
    for shell in (CORE_SHELL, MIXED_SHELL, DISC_SHELL):
        high = _high(shell)
        untiered = list(brightStars.scatter_shell(SHAPE, shell, high, EDGE_PC, E_VALUE, FLOOR, "tiers",
                                                  mass_range=(None, MASS)))
        tiered = list(brightStars.scatter_shell(SHAPE, shell, high, EDGE_PC, E_VALUE, FLOOR, "tiers",
                                                tiered=scatter_limits.luminosity_pass(limits)))
        assert untiered and tiered == untiered


def test_a_tier_whose_floor_does_not_fall_adds_nothing_to_a_staged_band():
    older = ScatterLimits.of(_skeleton(), FLOOR, MASS)
    newer = ScatterLimits.of(_skeleton(), 3000.0, MASS)
    staged = scatter_limits.band_pass(newer, older)
    assert staged.bands[HALO] is None                      # the halo's own floor did not move with the baseline
    assert staged.bands[BASELINE] == scatter_limits.PassBand(3000.0, FLOOR, (None, MASS))
    assert staged.bands[CORE].min_luminosity_sol == pytest.approx(3000.0 * tuning.CORE_TIER_FACTOR)
    assert staged.bands[CORE].max_luminosity_sol == pytest.approx(FLOOR * tuning.CORE_TIER_FACTOR)
    # A core floor that was set does not move either.
    fixed = scatter_limits.band_pass(ScatterLimits.of(_skeleton(core_floor=20000.0), 3000.0, MASS),
                                     ScatterLimits.of(_skeleton(core_floor=20000.0), FLOOR, MASS))
    assert fixed.bands[CORE] is None


def test_the_weights_count_each_bands_own_tier():
    limits = ScatterLimits.of(_skeleton(), FLOOR, MASS)
    tiered = scatter_limits.luminosity_pass(limits)
    high = _high(MIXED_SHELL)
    weight, expected = brightStars.shell_weight(SHAPE, MIXED_SHELL, high, EDGE_PC, E_VALUE,
                                                brightStars.tier_fractions_of(tiered))
    _, flat = brightStars.shell_weight(SHAPE, MIXED_SHELL, high, EDGE_PC, E_VALUE,
                                       brightStars.band_fractions(FLOOR, None, (None, MASS)))
    assert expected > 5 * flat > 0.0         # the halo's low floor draws far more stars than the baseline's
    drawn = sum(len(rows) for rows in _drawn(MIXED_SHELL, tiered).values())
    assert drawn == pytest.approx(expected, rel=0.25)


@pytest.mark.db
def test_the_marks_the_fills_and_the_clears_read_each_bands_tier(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=MIXED_SHELL + 1,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config, galaxy_seed=SEED)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        assert store.scatter_limits_of(conn, skeleton).tiers is None
        assert store.bright_star_fill_level(conn, CORE_SHELL, 0, 0) is None
        store.record_bright_star_scatter(conn, FLOOR, 7, MASS)
        conn.commit()
        limits = store.scatter_limits_of(conn, skeleton)
        assert store.bright_star_fill_level(conn, CORE_SHELL, 0, 0) == pytest.approx(FLOOR * tuning.CORE_TIER_FACTOR)
        assert store.bright_star_fill_level(conn, DISC_SHELL, 0, 0) == FLOOR
        high = _high(MIXED_SHELL)
        assert store.bright_star_fill_level(conn, MIXED_SHELL, -high, 0) == tuning.HALO_MIN_LUMINOSITY_SOL
        assert store.bright_star_fill_mass_limit(conn, CORE_SHELL, 0, 0) == pytest.approx(MASS * 1.1)
        assert limits.level_at(DISC_SHELL, 0) == FLOOR
    finally:
        conn.close()


@pytest.mark.db
def test_the_marks_skip_a_sector_by_its_own_tiers_limits(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=MIXED_SHELL + 1,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config, galaxy_seed=SEED)
    conn = store.get_connection(mysql_config)
    try:
        limits = store.scatter_limits_of(conn, store.get_galaxy_shape(conn))
        limits = ScatterLimits.of(store.get_galaxy_shape(conn), FLOOR, MASS)
        mass_pass = scatter_limits.mass_pass(limits, 100.0)
        written = 0
        for shell in (CORE_SHELL, DISC_SHELL):
            for row in brightStars.scatter_shell(SHAPE, shell, _high(shell), EDGE_PC, E_VALUE, 100.0, "marks",
                                                 tiered=mass_pass):
                store.insert_bright_stars(conn, [row])
                written += 1
        conn.commit()
        assert written
        marked = store.bright_star_marked_addresses(conn, limits)
        rows = conn.execute("SELECT shell_index, band_index, wedge_index, initial_mass_sol, luminosity_w"
                            " FROM bright_stars").fetchall()
        for row in rows:
            address = (row["shell_index"], row["band_index"], row["wedge_index"])
            tier = limits.limit_at(address[0], address[1])
            bright = row["luminosity_w"] >= tier.min_luminosity_sol * constants.SOLAR_LUMINOSITY
            if bright and row["initial_mass_sol"] >= tier.mass_limit_sol:
                assert address in marked
        assert marked <= {(row["shell_index"], row["band_index"], row["wedge_index"]) for row in rows}
    finally:
        conn.close()


@pytest.mark.db
def test_clearing_one_pass_tells_apart_the_stars_between_the_tiers_mass_limits(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=MIXED_SHELL + 1,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config, galaxy_seed=SEED)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        limits = ScatterLimits.of(skeleton, FLOOR, MASS)
        window = (MASS + MASS * tuning.CORE_TIER_FACTOR) / 2.0    # between the baseline's limit and the core's
        template = list(brightStars.scatter_shell(SHAPE, DISC_SHELL, _high(DISC_SHELL), EDGE_PC, E_VALUE, 100.0, "clear",
                                                  tiered=scatter_limits.mass_pass(limits, 100.0)))[0]
        core_row = list(template)
        core_row[:3] = (CORE_SHELL, 0, 0)
        disc_row = list(template)
        disc_row[:3] = (DISC_SHELL, 0, 0)
        for row in (core_row, disc_row):
            row[COLUMN["initial_mass_sol"]] = window
        store.insert_bright_stars(conn, [tuple(core_row), tuple(disc_row)])
        conn.commit()
        # In the disc a star of that mass is the mass pass's; in the core it is the luminosity pass's.
        store.clear_bright_star_pass(conn, "mass", limits)
        left = [(row["shell_index"]) for row in conn.execute("SELECT shell_index FROM bright_stars").fetchall()]
        assert left == [CORE_SHELL]
        store.insert_bright_stars(conn, [tuple(disc_row)])
        conn.commit()
        store.clear_bright_star_pass(conn, "luminosity", limits)
        left = [(row["shell_index"]) for row in conn.execute("SELECT shell_index FROM bright_stars").fetchall()]
        assert left == [DISC_SHELL]
    finally:
        conn.close()


@pytest.mark.db
def test_the_api_names_a_sectors_tier_and_its_limits(mysql_config):
    from planetgen.api import routes
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=MIXED_SHELL + 1,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config, galaxy_seed=SEED)
    detail = {"placed": True, "shell_index": CORE_SHELL, "band_index": 0, "wedge_index": 0}
    conn = store.get_connection(mysql_config)
    try:
        assert routes.scatter_tier_of(conn, detail) is None                 # no scatter has run
        store.record_bright_star_scatter(conn, FLOOR, 7, MASS)
        conn.commit()
        assert routes.scatter_tier_of(conn, detail) == {
            "tier": "core", "min_luminosity_sol": pytest.approx(9900.0), "mass_limit_sol": pytest.approx(15.4)}
        disc = routes.scatter_tier_of(conn, {**detail, "shell_index": DISC_SHELL})
        assert disc == {"tier": "baseline", "min_luminosity_sol": FLOOR, "mass_limit_sol": MASS}
        assert routes.scatter_tier_of(conn, {**detail, "placed": False}) is None
    finally:
        conn.close()
