# tests/test_influence.py

"""GEN.341: the influence-set search. The set is the masses within the largest Hill radius of the galaxy, read from the
sectors that can hold part of that sphere, checked against brute force and a `cKDTree`, with a negative control for the
cells' reach margin, the cap of exact members and the order the result is given in."""

import math
import random

import numpy as np
import pytest
from scipy.spatial import cKDTree

from planetgen.galaxy import enclosed_mass, influence, keepout, shape_draw, spherical
from planetgen.physics import constants

EDGE_PC = 4.0
SEED = bytes(range(31, 47))
SITES = {
    "core": (3.0, -2.0, 1.0),
    "ring face": (4000.0, 0.0, 0.0),
    "mid radius": (8128.0, 0.0, 0.0),
    "rim": (14900.0, 0.0, 0.0),
    "slot-0 seam above": (8128.0, 0.01, 0.0),
    "slot-0 seam below": (8128.0, -0.01, 0.0),
    "top layer": (8128.0, 100.0, 600.0),
    "bottom layer": (8128.0, 100.0, -600.0),
}
SPREAD_PC = 70.0
PER_SITE = 6000


@pytest.fixture(scope="module")
def shape():
    return shape_draw.default_shape(SEED)


def _salpeter(rng, n):
    """Masses, kg, from a Salpeter law over 0.1 to 100 solar masses."""
    u = np.array([rng.random() for _ in range(n)])
    low, high, slope = 0.1, 100.0, 2.35
    mass = (low ** (1 - slope) + u * (high ** (1 - slope) - low ** (1 - slope))) ** (1 / (1 - slope))
    return mass * constants.SOLAR_MASS_TO_KG


def _catalogue(extra=()):
    """Positions (n, 3) pc, masses (n,) kg and keys of a seeded catalogue around every site, plus `extra` rows."""
    rng = random.Random(341)
    positions, masses = [], []
    for site in SITES.values():
        points = []
        while len(points) < PER_SITE:
            p = np.array([rng.uniform(-SPREAD_PC, SPREAD_PC) for _ in range(3)])
            if p @ p <= SPREAD_PC ** 2:
                points.append(np.array(site) + p)
        positions.extend(points)
        masses.extend(_salpeter(rng, PER_SITE))
    for position, mass in extra:
        positions.append(np.array(position, dtype=float))
        masses.append(mass)
    return np.array(positions), np.array(masses), tuple(range(len(masses)))


def _tables(positions, masses, keys, shuffle=None):
    """`{address: SectorMasses}` of the catalogue, its sectors and members stored in the order `shuffle` gives."""
    members = {}
    for i, position in enumerate(positions):
        members.setdefault(spherical.sector_address_at(tuple(position), EDGE_PC), []).append(i)
    tables = {}
    for address, indices in members.items():
        if shuffle is not None:
            shuffle.shuffle(indices)
        tables[address] = influence.SectorMasses.of(address, positions[indices], masses[indices],
                                                    [keys[i] for i in indices])
    if shuffle is not None:
        items = list(tables.items())
        shuffle.shuffle(items)
        tables = dict(items)
    return tables


@pytest.fixture(scope="module")
def catalogue():
    positions, masses, keys = _catalogue()
    return positions, masses, keys, _tables(positions, masses, keys)


def _targets(count=6):
    rng = random.Random(7)
    return {name: [tuple(np.array(site) + np.array([rng.uniform(-20.0, 20.0) for _ in range(3)])) for _ in range(count)]
            for name, site in SITES.items()}


def test_the_set_is_brute_force_and_the_kd_tree_at_every_site_and_radius(catalogue):
    positions, _masses, keys, tables = catalogue
    tree = cKDTree(positions)
    for name, targets in _targets().items():
        for radius in (2.0, 6.0, 20.0, 60.0):
            for target in targets:
                found = influence.influence_set(target, radius, EDGE_PC, tables.get, max_exact_members=10 ** 9)
                brute = np.flatnonzero(np.sqrt(((positions - np.array(target)) ** 2).sum(axis=1)) <= radius)
                assert set(found.keys) == {keys[i] for i in brute}, (name, radius)
                assert set(found.keys) == {keys[i] for i in tree.query_ball_point(target, radius)}, (name, radius)
                assert not found.truncated
                assert found.member_count == len(brute)


def test_without_the_reach_margin_neighbours_are_missed(catalogue):
    positions, _masses, keys, tables = catalogue
    wanted = missed = 0
    for targets in _targets(count=4).values():
        for target in targets:
            brute = {keys[i] for i in np.flatnonzero(np.sqrt(((positions - np.array(target)) ** 2).sum(axis=1)) <= 6.0)}
            bare = influence.influence_set(target, 6.0, EDGE_PC, tables.get, max_exact_members=10 ** 9, reach_edges=0.0)
            full = influence.influence_set(target, 6.0, EDGE_PC, tables.get, max_exact_members=10 ** 9)
            assert set(bare.keys) <= brute and set(full.keys) == brute
            wanted += len(brute)
            missed += len(brute) - len(bare.keys)
    assert wanted > 0 and missed / wanted > 0.05


def test_every_object_whose_own_hill_sphere_holds_the_target_is_in_the_set(shape):
    heavy = [((8128.0 + 9.0, 4.0, 0.0), 3.0e3 * constants.SOLAR_MASS_TO_KG),
             ((14900.0 - 5.0, 0.0, 6.0), 4.0e4 * constants.SOLAR_MASS_TO_KG)]
    positions, masses, keys = _catalogue(heavy)
    tables = _tables(positions, masses, keys)
    radii = np.sqrt((positions ** 2).sum(axis=1))
    hill = influence.galactic_hill_radius_pc(masses, radii, shape, SEED)
    radius = influence.largest_hill_radius_pc(masses, radii, shape, SEED)
    assert radius == pytest.approx(hill.max())
    assert radius > 20.0          # the injected heavy objects reach well past a field star's few parsecs
    rng = random.Random(5)
    checked = 0
    for index in rng.sample(range(len(keys)), 60) + [len(keys) - 1, len(keys) - 2]:
        target = tuple(positions[index])
        found = set(influence.influence_set(target, radius, EDGE_PC, tables.get, exclude_keys=[keys[index]],
                                            max_exact_members=10 ** 9).keys)
        distance = np.sqrt(((positions - positions[index]) ** 2).sum(axis=1))
        inside_hill = {keys[j] for j in np.flatnonzero(distance < hill) if j != index}
        assert inside_hill <= found
        checked += len(inside_hill)
    assert checked > 0


def test_the_cap_keeps_the_nearest_exact_and_folds_the_rest_by_sector(catalogue):
    positions, masses, keys, tables = catalogue
    target = SITES["mid radius"]
    radius = 30.0
    found = influence.influence_set(target, radius, EDGE_PC, tables.get)
    distance = np.sqrt(((positions - np.array(target)) ** 2).sum(axis=1))
    inside = np.flatnonzero(distance <= radius)
    assert len(inside) > influence.DEFAULT_MAX_EXACT_MEMBERS
    nearest = inside[np.argsort(distance[inside], kind="stable")][:influence.DEFAULT_MAX_EXACT_MEMBERS]
    assert found.keys == tuple(keys[i] for i in nearest)
    assert found.truncated and found.member_count == len(inside)
    assert found.total_mass_kg == pytest.approx(masses[inside].sum(), rel=1e-12)
    assert [m.address for m in found.monopoles] == sorted({m.address for m in found.monopoles})
    folded = len(inside) - influence.DEFAULT_MAX_EXACT_MEMBERS
    assert sum(m.count for m in found.monopoles) == folded

    def pull(points, weights):
        offset = points - np.array(target)
        return (weights[:, None] * offset / (offset ** 2).sum(axis=1)[:, None] ** 1.5).sum(axis=0)

    full = pull(positions[inside], masses[inside])
    capped = pull(found.positions_pc, found.masses_kg) + pull(
        np.array([m.position_pc for m in found.monopoles]), np.array([m.mass_kg for m in found.monopoles]))
    assert np.linalg.norm(capped - full) / np.linalg.norm(full) < 0.05


def test_the_same_catalogue_in_another_order_gives_the_same_set(catalogue):
    positions, masses, keys, tables = catalogue
    reference = influence.influence_set(SITES["rim"], 25.0, EDGE_PC, tables.get)
    for seed in (1, 2, 3):
        shuffled = _tables(positions, masses, keys, shuffle=random.Random(seed))
        again = influence.influence_set(SITES["rim"], 25.0, EDGE_PC, shuffled.get)
        assert again.keys == reference.keys
        assert np.array_equal(again.positions_pc, reference.positions_pc)
        assert again.monopoles == reference.monopoles


def test_a_lone_rogue_has_an_empty_set_and_a_hierarchy_counts_once():
    far = np.array([[8128.0, 0.0, 0.0]])
    tables = _tables(far, np.array([1.0e30]), (0,))
    assert influence.influence_set((8128.0, 50.0, 0.0), 10.0, EDGE_PC, tables.get).member_count == 0
    empty = influence.influence_set((8128.0, 50.0, 0.0), 10.0, EDGE_PC, tables.get)
    assert empty.keys == () and empty.total_mass_kg == 0.0 and not empty.truncated
    binary = np.array([[8128.0, 3.0, 0.0]])        # a pair given as one point mass at its barycentre
    tables = _tables(binary, np.array([2.0e30]), ("pair",))
    found = influence.influence_set((8128.0, 0.0, 0.0), 10.0, EDGE_PC, tables.get)
    assert found.keys == ("pair",) and found.masses_kg[0] == 2.0e30


def test_the_object_itself_and_its_own_system_are_left_out(catalogue):
    positions, _masses, keys, tables = catalogue
    index = 100
    everything = influence.influence_set(tuple(positions[index]), 20.0, EDGE_PC, tables.get, max_exact_members=10 ** 9)
    assert keys[index] in everything.keys
    others = influence.influence_set(tuple(positions[index]), 20.0, EDGE_PC, tables.get, exclude_keys=[keys[index]],
                                     max_exact_members=10 ** 9)
    assert set(others.keys) == set(everything.keys) - {keys[index]}


def test_a_search_that_cannot_be_made_is_refused():
    with pytest.raises(ValueError):
        influence.influence_set((0, 0, 0), -1.0, EDGE_PC, lambda a: None)
    with pytest.raises(ValueError):
        influence.influence_set((0, 0, 0), 1.0, 0.0, lambda a: None)
    with pytest.raises(ValueError):
        influence.influence_set((0, 0, 0), 1.0, EDGE_PC, lambda a: None, max_exact_members=0)
    with pytest.raises(ValueError):
        influence.SectorMasses.of((1, 1, 1), [[0, 0, 0]], [1.0, 2.0], ["a"])


def test_the_search_radius_is_the_galaxys_largest_hill_radius_without_the_central_black_hole(shape):
    solar = constants.SOLAR_MASS_TO_KG
    masses = np.array([0.5, 1.0, 40.0, 1.0e3, 4.0e6]) * solar
    radii = np.array([5000.0, 8000.0, 3000.0, 6000.0, 0.0])
    hill = influence.galactic_hill_radius_pc(masses, radii, shape, SEED)
    # the same Hill sphere a star's own perimeter uses, from the mass inside the object's own radius
    central = enclosed_mass.enclosed_mass_kg(shape, SEED, radii[1])
    assert hill[1] == pytest.approx(keepout.galactic_hill_radius_km(masses[1], radii[1], central) * 1000.0 / constants.PARSEC_M)
    with_black_hole = influence.largest_hill_radius_pc(masses, radii, shape, SEED)
    without = influence.largest_hill_radius_pc(masses, radii, shape, SEED, skip=[False, False, False, False, True])
    assert without == pytest.approx(hill[:4].max())
    assert with_black_hole >= without
    assert influence.largest_hill_radius_pc([], [], shape, SEED) == 0.0
    assert influence.largest_hill_radius_pc(masses[:1], radii[:1], shape, SEED, skip=[True]) == 0.0


def test_a_bigger_galaxy_has_a_wider_search_radius(shape):
    other = shape._replace(disk_scale_length_pc=2.0 * shape.disk_scale_length_pc,
                           disk_scale_height_pc=2.0 * shape.disk_scale_height_pc,
                           bulge_scale_radius_pc=2.0 * shape.bulge_scale_radius_pc,
                           spiral_reference_radius_pc=2.0 * shape.spiral_reference_radius_pc)
    mass, radius = [1.0e3 * constants.SOLAR_MASS_TO_KG], [3.0 * shape.disk_scale_length_pc]
    assert influence.largest_hill_radius_pc(mass, radius, other, SEED) > influence.largest_hill_radius_pc(
        mass, radius, shape, SEED)
