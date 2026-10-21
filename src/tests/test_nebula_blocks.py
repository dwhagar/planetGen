"""
MAP.172: the nebula aggregates of the Galaxy Map's block ladder (`galaxy/nebula_blocks.py`). A point is found in the same
block `drill.py` bounds describe, the blocks of a level tile a shell, the counts, volumes and cover add up through the
ladder, and the stored cover is within 5% of the true cover found by throwing random points at the block
(docs/design/nebula-map-visibility.md section 9).
"""

import math

import numpy as np
import pytest
from scipy.spatial import cKDTree

from planetgen import tuning
from planetgen.galaxy import drill, nebula_blocks as rules, spherical

EDGE_PC = float(tuning.DEFAULT_SECTOR_EDGE_PC)
TWO_PI = 2 * math.pi


def _spheres(seed, count, centre, spread, radius=(4.0, 12.0), letters="MNMNMNP"):
    rng = np.random.default_rng(seed)
    where = np.asarray(centre) + rng.uniform(-spread, spread, (count, 3))
    r = rng.uniform(*radius, count)
    codes = rules.class_codes(rng.choice(list(letters), count))
    return where[:, 0], where[:, 1], where[:, 2], r, codes


def _by_key(level):
    return {(row["shell_index"], row["band_index"], row["wedge_index"]): row for row in level.blocks}


def test_a_point_is_in_the_block_whose_bounds_hold_it():
    rng = np.random.default_rng(3)
    points = rng.uniform(-9000, 9000, (400, 3)) * rng.uniform(0.001, 1.0, (400, 1))
    for m in rules.LEVELS:
        shells, bands, wedges = rules.locate(points[:, 0], points[:, 1], points[:, 2], m, EDGE_PC)
        for point, shell, band, wedge in zip(points, shells, bands, wedges):
            (inner, outer), (north, south), (start, end) = drill.drill_block_bounds(
                drill.DrillBlock(m, int(shell), int(band), int(wedge)), EDGE_PC)
            radius = float(np.linalg.norm(point))
            colatitude = math.atan2(math.hypot(point[0], point[1]), point[2])
            longitude = math.atan2(point[1], point[0]) % TWO_PI
            assert inner - 1e-9 <= radius < outer + 1e-9
            assert north - 1e-12 <= colatitude <= south + 1e-12
            assert start - 1e-12 <= longitude <= end + 1e-12


def test_a_sectors_centre_is_in_the_block_chain_of_the_sector():
    checked = 0
    for shell in (3, 17, 90, 400):
        low, high = spherical.band_range(shell)
        for band in sorted({low + 1, low // 2, 0, high // 2, high - 1}):
            if spherical.is_cap(shell, band):
                continue
            count = spherical.wedge_count(shell, band)
            for wedge in sorted({0, count // 3, count - 1}):
                x, y, z = spherical.sector_position_pc(shell, band, wedge, EDGE_PC)
                chain = drill.drill_chain_of(shell, band, wedge)[:-1]
                for block in chain:
                    found = rules.locate(np.array([x]), np.array([y]), np.array([z]), block.m, EDGE_PC)
                    assert tuple(int(v[0]) for v in found) == (block.shell, block.band, block.wedge), (
                        shell, band, wedge, block)
                    checked += 1
    assert checked > 100


@pytest.mark.parametrize("m", rules.LEVELS)
def test_the_blocks_of_a_block_shell_fill_its_shell(m):
    shell = 2
    total = 0.0
    for row in range(drill.drill_band_count(m, shell)):
        band = (drill.drill_band_count(m, shell) - 1) // 2 - row
        for wedge in range(drill.drill_wedge_count(m, shell, band)):
            total += rules.block_volume_pc3(m, shell, band, wedge, EDGE_PC)
    inner, outer = (shell * m - 0.5) * EDGE_PC, ((shell + 1) * m - 0.5) * EDGE_PC
    assert total == pytest.approx(4 / 3 * math.pi * (outer ** 3 - inner ** 3), rel=1e-9)


def test_cover_volume_is_all_the_nebula_volume_at_every_level():
    x, y, z, r, codes = _spheres(1, 800, (700.0, 150.0, 40.0), 250.0)
    volume = float(rules.sphere_volumes(r).sum())
    for m in rules.LEVELS:
        level = rules.level_rows(x, y, z, r, codes, m, EDGE_PC)
        assert sum(row["cover_volume_pc3"] for row in level.blocks) == pytest.approx(volume, rel=1e-9)
        assert sum(row["volume_pc3"] for row in level.blocks) == pytest.approx(volume, rel=1e-9)
        assert sum(row["nebula_count"] for row in level.blocks) == len(x)


def test_counts_volumes_and_families_add_up_through_the_ladder():
    x, y, z, r, codes = _spheres(2, 900, (-500.0, 900.0, -100.0), 300.0, letters="CDFGHMNPRS")
    levels = {m: _by_key(rules.level_rows(x, y, z, r, codes, m, EDGE_PC)) for m in rules.LEVELS}
    for m in rules.LEVELS[:-1]:
        child_m = m // drill.DRILL_STEP
        sums = {}
        for (shell, band, wedge), row in levels[child_m].items():
            parent = drill.drill_parent(drill.DrillBlock(child_m, shell, band, wedge))
            key = (parent.shell, parent.band, parent.wedge)
            total = sums.setdefault(key, [0, 0.0, [0.0] * len(rules.FAMILIES), 0.0])
            total[0] += row["nebula_count"]
            total[1] += row["volume_pc3"]
            total[2] = [a + b for a, b in zip(total[2], row["families"])]
            total[3] += row["cover_volume_pc3"]
        for key, (count, volume, families, cover_volume) in sums.items():
            parent_row = levels[m][key]
            assert parent_row["nebula_count"] == count
            assert parent_row["volume_pc3"] == pytest.approx(volume, rel=1e-9)
            assert parent_row["families"] == pytest.approx(families, rel=1e-9)
            assert parent_row["cover_volume_pc3"] == pytest.approx(cover_volume, rel=1e-6)


def test_a_block_holds_its_families_dominant_class_centroid_and_bound():
    # Three nebulae at known places near the middle of one level-3 block: a big M, a small N and a small C.
    shell, band, wedge = (int(v[0]) for v in rules.locate(*(np.array([c]) for c in (800.0, 200.0, 50.0)), 3, EDGE_PC))
    (inner, outer), (north, south), (start, end) = drill.drill_block_bounds(drill.DrillBlock(3, shell, band, wedge), EDGE_PC)
    middle, colatitude, longitude = (inner + outer) / 2, (north + south) / 2, (start + end) / 2
    centre = middle * np.array([math.sin(colatitude) * math.cos(longitude), math.sin(colatitude) * math.sin(longitude),
                                math.cos(colatitude)])
    inside = [centre, centre + [1.0, 0.5, -0.5], centre + [-0.5, 1.0, 0.5]]
    radii = np.array([3.0, 1.0, 0.5])
    letters = ["M", "N", "C"]
    xs, ys, zs = (np.array([p[i] for p in inside]) for i in range(3))
    level = rules.level_rows(xs, ys, zs, radii, rules.class_codes(letters), 3, EDGE_PC)
    row = _by_key(level)[(shell, band, wedge)]
    assert row["nebula_count"] == 3
    volumes = rules.sphere_volumes(radii)
    assert row["volume_pc3"] == pytest.approx(volumes.sum())
    assert row["dominant_class"] == "M"
    shares = dict(zip(rules.FAMILIES, row["families"]))
    assert shares["dark"] == pytest.approx(volumes[0] + volumes[1])
    assert shares["emission"] == pytest.approx(volumes[2])
    expected = sum(v * np.array(p) for v, p in zip(volumes, inside)) / volumes.sum()
    assert [row["centroid_x_pc"], row["centroid_y_pc"], row["centroid_z_pc"]] == pytest.approx(list(expected))
    centroid = np.array([row["centroid_x_pc"], row["centroid_y_pc"], row["centroid_z_pc"]])
    reach = max(float(np.linalg.norm(np.array(p) - centroid)) + r for p, r in zip(inside, radii))
    assert row["bound_radius_pc"] == pytest.approx(reach)


def test_a_block_a_cloud_only_reaches_has_a_cover_and_no_nebulae():
    # A cloud wider than a level-3 block sits across two of them; the one that holds its centre counts it.
    x, y, z = (np.array([v]) for v in (800.0, 200.0, 50.0))
    level = rules.level_rows(x, y, z, np.array([20.0]), rules.class_codes(["M"]), 3, EDGE_PC)
    counts = sorted(row["nebula_count"] for row in level.blocks)
    assert counts.count(1) == 1 and counts.count(0) == len(counts) - 1 and len(counts) > 4
    assert all(0 < row["cover"] <= 1 for row in level.blocks)
    assert all(row["bound_radius_pc"] is None for row in level.blocks if row["nebula_count"] == 0)


def _block_centre(m, shell, band, wedge, edge_pc):
    (inner, outer), (north, south), (start, end) = drill.drill_block_bounds(drill.DrillBlock(m, shell, band, wedge), edge_pc)
    middle, colatitude, longitude = (inner + outer) / 2, (north + south) / 2, (start + end) / 2
    return middle * np.array([math.sin(colatitude) * math.cos(longitude), math.sin(colatitude) * math.sin(longitude),
                              math.cos(colatitude)])


def _true_cover(block, m, edge_pc, centres, radii, tree, rng, points=40000):
    """The share of the block inside some sphere, from random points thrown at the block."""
    shell, band, wedge = block
    (inner, outer), (north, south), (start, end) = drill.drill_block_bounds(drill.DrillBlock(m, shell, band, wedge), edge_pc)
    radius = np.cbrt(rng.uniform(inner ** 3, outer ** 3, points))
    cos_theta = rng.uniform(math.cos(south), math.cos(north), points)
    sin_theta = np.sqrt(1 - cos_theta ** 2)
    phi = rng.uniform(start, end, points)
    sample = np.column_stack((radius * sin_theta * np.cos(phi), radius * sin_theta * np.sin(phi), radius * cos_theta))
    inside = 0
    for point, near in zip(sample, tree.query_ball_point(sample, float(radii.max()))):
        if near and (np.linalg.norm(centres[near] - point, axis=1) <= radii[near]).any():
            inside += 1
    return inside / points


@pytest.mark.parametrize("m", [243, 81, 27])
@pytest.mark.parametrize("direction", [(1.0, 0.3, 0.2), (-0.5, 0.8, -0.3), (0.2, -0.9, 0.4)])
def test_the_stored_cover_is_within_five_percent_of_the_true_cover(m, direction):
    """The check the item asks for: random points thrown at a block, the share inside some nebula against the stored
    cover. The block is one of the levels the cover shades (the ladder's 256 pc cells and up, here on a 1 pc edge so
    the clouds are a quarter of their size against the same block), in a field filled to 11%, the fill of an arm."""
    edge = 1.0
    rng = np.random.default_rng(int(m * 31 + sum(direction) * 1000))
    unit = np.asarray(direction) / np.linalg.norm(direction)
    shell, band, wedge = (int(v[0]) for v in rules.locate(*(np.array([c]) for c in unit * m * 1.5 * edge), m, edge))
    centre = _block_centre(m, shell, band, wedge, edge)
    half = 0.75 * m * edge + 4.0
    radius_range = (1.0, 3.5)
    count = int(0.11 * (2 * half) ** 3 / (4 / 3 * math.pi * np.mean(np.linspace(*radius_range, 200) ** 3)))
    centres = centre + rng.uniform(-half, half, (count, 3))
    radii = rng.uniform(*radius_range, count)
    level = rules.level_rows(centres[:, 0], centres[:, 1], centres[:, 2], radii,
                             rules.class_codes(["M"] * count), m, edge)
    row = _by_key(level)[(shell, band, wedge)]
    true = _true_cover((shell, band, wedge), m, edge, centres, radii, cKDTree(centres), rng)
    assert true > 0.05
    assert abs(row["cover"] - true) / true < 0.05, (row["cover"], true)
