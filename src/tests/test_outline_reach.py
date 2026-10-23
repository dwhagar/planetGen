# tests/test_outline_reach.py

"""
GEN.627: the outline takes a component in by the radius where its density reaches one expected system in a cell, not by
its tidal radius. A diffuse dwarf galaxy whose tidal sphere covers the host but whose density never reaches one system adds
no cells; a cluster reaches as far as its density does and no farther than its tidal radius; the plan prints how many sectors
really expect a system beside the candidates; a random start lands in a sector that expects one.
"""

import math
import types

import numpy as np
import pytest

from planetgen.galaxy import component_field, components, skeleton, spherical
from planetgen.galaxy import density as galaxyDensity
from planetgen.generation import plan_outline, run_galaxy
from planetgen.physics import king
from tests.test_component_field import E1, EDGE, _cluster, _lattice, _shape

SEED = bytes.fromhex("2139DD5841E680B009A3D17A6319B2AA")
THRESHOLD = 1.0 / E1


def _diffuse_satellite(center):
    """A Sagittarius-like dwarf: 5e7 solar masses spread over a 2 kpc half-mass radius, 1.8e4 systems a kiloparsec cube
    at best, a tidal sphere of 12 kpc, and a central density far under one system in a cell."""
    return components.Component(
        kind="satellite", name="Diffuse", source="generated", center_pc=center, velocity_kms=(0.0, 0.0, 0.0), epoch_unix=0.0,
        mass_solar=5.0e7, half_mass_radius_pc=2000.0, profile="exponential", age_gy=8.0, metallicity_fe_h=-1.0, seed=2,
        tidal_radius_pc=12000.0)


def test_a_dwarf_whose_density_never_reaches_one_system_in_a_cell_has_no_reach():
    field = component_field.ComponentField([_diffuse_satellite((15000.0, 0.0, 5000.0))], EDGE, E1)
    assert field._density_reach[0] == 0.0 and field._reach[0] == 12000.0
    assert field.qualifying_volume_pc3() == 0.0
    assert field.high_band(3500, spherical.band_at, EDGE) == -1


def test_such_a_dwarf_adds_no_cells_to_the_outline():
    plain, outer, _ = skeleton.build_shell_extents(_shape(), EDGE, THRESHOLD, max_shell=4000)
    shape = _shape([_diffuse_satellite((6000.0, 0.0, 6000.0))])
    with_dwarf, outer_dwarf, _ = skeleton.build_shell_extents(shape, EDGE, THRESHOLD, max_shell=4000)
    assert outer_dwarf == outer                    # its faint density still adds a band or so to the bound where it overlaps
    assert skeleton.shell_candidate_sector_count(with_dwarf) == pytest.approx(
        skeleton.shell_candidate_sector_count(plain), rel=0.002)


def test_a_cluster_reaches_where_its_density_does_and_no_farther_than_its_tidal_radius():
    cluster = _cluster(2.0e5, 5.0, 5)
    field = component_field.ComponentField([cluster], EDGE, E1)
    reach, tidal = float(field._density_reach[0]), float(field._reach[0])
    assert 0.0 < reach <= tidal
    cells = _lattice(field, tidal)
    n = int(tidal / EDGE) + 3
    axis = np.arange(-n, n + 1) * EDGE
    distance = np.sqrt(axis[:, None, None] ** 2 + axis[None, :, None] ** 2 + axis[None, None, :] ** 2)
    assert (cells >= 1.0).sum() > 100
    assert distance[cells >= 1.0].max() <= reach


def test_a_component_holding_under_one_system_in_all_never_qualifies_a_cell():
    tiny = _cluster(0.3, 1.0, 5)
    assert component_field.ComponentField([tiny], EDGE, E1)._density_reach[0] == 0.0


def test_the_outline_walks_a_cluster_off_the_plane_by_its_density_reach():
    center = (3000.0, 0.0, 12000.0)
    shape = _shape([_cluster(2.0e5, 5.0, 5, center=center)])
    extents, _outer, _confirmed = skeleton.build_shell_extents(shape, EDGE, THRESHOLD, max_shell=3400)
    high = dict(extents)
    reach = float(shape.components._density_reach[0])
    assert reach < king.tidal_radius_pc(5, 5.0) + 0.001
    for cell in spherical.cells_touching_sphere(center, reach, EDGE):
        shell, band = cell[:2]
        position = spherical.sector_position_pc(*cell[:3], EDGE)
        if shape.components.relative_density(position) >= THRESHOLD:
            assert shell in high and abs(band) <= high[shell]


def test_the_plan_prints_how_many_sectors_expect_a_system_beside_the_candidates():
    built = plan_outline.build(plan_outline._parser().parse_args(["--mass", "1e8", "--disk-scale-length", "2000pc"]), SEED, 4.0)
    made = plan_outline.report(built)
    assert 0 < made.reaching <= made.sectors
    assert f"{made.sectors:,} candidate sectors, about {made.reaching:,} of them expect at least one system" in made.line()


def test_the_reaching_count_matches_a_count_of_the_sectors_that_expect_a_system():
    """The integral against the model's own density on a coarse lattice of the disc."""
    built = plan_outline.build(plan_outline._parser().parse_args(["--mass", "1e8", "--disk-scale-length", "400pc"]), SEED, 4.0)
    made = plan_outline.report(built)
    edge, reach = built.edge_pc, made.radius_pc
    step = 3                                   # every third sector on each axis
    axis = np.arange(-reach, reach + edge, step * edge)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    rho = galaxyDensity.total_density_array(x, y, z, built.shape)
    counted = float((rho * built.e_value >= 1.0).sum()) * step ** 3
    assert made.reaching == pytest.approx(counted, rel=0.25)


def test_a_random_start_lands_in_a_sector_that_expects_a_system(monkeypatch):
    built = plan_outline.build(plan_outline._parser().parse_args(["--mass", "1e8", "--disk-scale-length", "2000pc"]), SEED, 4.0)
    shape = built.shape
    batch = types.SimpleNamespace(
        bounds=skeleton.GalaxyBounds(built.extents, built.edge_pc),
        skeleton=types.SimpleNamespace(expected_system_count_at_density_1=built.e_value),
        resolve=lambda args, address, position: types.SimpleNamespace(
            density=galaxyDensity.relative_density(position, shape), num_systems=None))
    monkeypatch.setattr(run_galaxy.store, "get_sector_id_at", lambda conn, *address: None)
    args = types.SimpleNamespace(neighborhoods=1, neighborhood_gamma=0.0, max_shell=None, avoid_filled_space=False,
                                 min_start_density=None, density=None, num_systems=None)
    for _ in range(1000):
        picked = run_galaxy._pick_random_start(args, built.edge_pc, batch, None, 10.0, [])
        assert picked is not None
        assert picked[2].density * built.e_value >= 1.0


def test_a_start_with_a_density_given_is_not_held_to_the_threshold(monkeypatch):
    built = plan_outline.build(plan_outline._parser().parse_args(["--mass", "1e8", "--disk-scale-length", "2000pc"]), SEED, 4.0)
    args = types.SimpleNamespace(neighborhoods=1, neighborhood_gamma=0.0, max_shell=None, avoid_filled_space=False,
                                 min_start_density=None, density=1.0e-9, num_systems=None)
    batch = types.SimpleNamespace(
        bounds=skeleton.GalaxyBounds(built.extents, built.edge_pc),
        skeleton=types.SimpleNamespace(expected_system_count_at_density_1=built.e_value),
        resolve=lambda given, address, position: given)
    monkeypatch.setattr(run_galaxy.store, "get_sector_id_at", lambda conn, *address: None)
    assert run_galaxy._pick_random_start(args, built.edge_pc, batch, None, 10.0, [])[2] is args
