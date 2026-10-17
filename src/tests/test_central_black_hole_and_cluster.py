# tests/test_central_black_hole_and_cluster.py

"""GEN.402: the central black hole and the nuclear star cluster come from the galaxy's own stellar mass, bulge and
seed (`galaxy.structure`, `galaxy.structure_draw`), not from fixed Milky Way masses."""

import math
import re
import statistics
from pathlib import Path

import pytest

from planetgen.galaxy import density, shape_draw, structure as galaxyStructure
from planetgen.galaxy import structure_draw as sd


def _seed(number):
    return number.to_bytes(16, "big")


SEEDS = [_seed(n) for n in range(1, 1201)]
SHAPES = [shape_draw.draw_galaxy(s) for s in SEEDS]
SOURCE = Path(__file__).resolve().parent.parent / "planetgen"


def _hole(drawn):
    return galaxyStructure.black_hole_mass_solar(drawn.structure, drawn.bulge_to_total)


def test_the_same_seed_gives_the_same_hole_and_cluster():
    first, second = shape_draw.draw_galaxy(_seed(5)), shape_draw.draw_galaxy(_seed(5))
    assert _hole(first) == _hole(second)
    assert first.structure.nuclear_cluster_mass_solar == second.structure.nuclear_cluster_mass_solar


def test_a_classical_bulge_follows_the_bulge_relation_and_the_rest_the_stellar_mass_relation():
    ref = galaxyStructure.REFERENCE_STRUCTURE
    classical = ref._replace(classical_bulge=True, stellar_mass_solar=1.0e11)
    assert galaxyStructure.black_hole_mass_solar(classical, 1.0) == pytest.approx(0.49e9)
    plain = ref._replace(classical_bulge=False, stellar_mass_solar=1.0e11)
    assert galaxyStructure.black_hole_mass_solar(plain, 0.9) == pytest.approx(10 ** 7.45)
    # A deviation of one sigma scales it by the relation's scatter.
    high = plain._replace(black_hole_deviation=1.0)
    assert galaxyStructure.black_hole_mass_solar(high, 0.0) / 10 ** 7.45 == pytest.approx(
        10 ** galaxyStructure.BLACK_HOLE_STELLAR_SCATTER_DEX)


def test_a_galaxy_without_a_classical_bulge_does_not_scale_its_hole_with_its_bulge_share():
    plain = galaxyStructure.REFERENCE_STRUCTURE._replace(classical_bulge=False)
    assert galaxyStructure.black_hole_mass_solar(plain, 0.01) == galaxyStructure.black_hole_mass_solar(plain, 0.5)


def test_the_hole_grows_with_the_stellar_mass_across_the_seeds():
    masses = [(s.structure.stellar_mass_solar, _hole(s)) for s in SHAPES]
    low = statistics.median(h for m, h in masses if m < 1.0e10)
    high = statistics.median(h for m, h in masses if m > 5.0e10)
    assert high > 5 * low
    assert all(1.0e3 < h < 1.0e11 for m, h in masses)


def test_the_m_sigma_relation_gives_the_dispersion_of_the_hole():
    assert galaxyStructure.black_hole_dispersion_kms(1.0e9 * 10 ** -0.510) == pytest.approx(200.0)
    assert galaxyStructure.black_hole_dispersion_kms(4.0e6) < galaxyStructure.black_hole_dispersion_kms(4.0e8)


def test_the_cluster_follows_the_power_law_between_the_two_verified_slopes():
    assert sd.nuclear_cluster_mass_solar(1.0e9, 0.9, 0.0) == pytest.approx(10 ** 6.51)
    assert sd.nuclear_cluster_mass_solar(1.0e11, 0.48, 0.0) == pytest.approx(10 ** (6.51 + 0.96))
    assert sd.nuclear_cluster_mass_solar(1.0e11, 1.36, 0.0) == pytest.approx(10 ** (6.51 + 2.72))
    assert sd.NUCLEAR_CLUSTER_SLOPE_RANGE == (0.48, 1.36)


def test_some_galaxies_have_no_nuclear_cluster_and_the_share_falls_with_mass():
    assert sd.nucleation_fraction(1.0e8) == pytest.approx(0.85)
    assert sd.nucleation_fraction(1.0e12) == pytest.approx(0.2)
    small = [s.structure.nuclear_cluster_mass_solar for s in SHAPES if s.structure.stellar_mass_solar < 1.0e9]
    large = [s.structure.nuclear_cluster_mass_solar for s in SHAPES if s.structure.stellar_mass_solar > 1.0e11]
    if small and large:
        assert sum(m == 0 for m in small) / len(small) < sum(m == 0 for m in large) / len(large)
    assert any(s.structure.nuclear_cluster_mass_solar == 0 for s in SHAPES)
    assert any(s.structure.nuclear_cluster_mass_solar > 0 for s in SHAPES)


def test_a_galaxy_without_a_cluster_has_none_in_its_density_terms():
    base = density.build_galaxy_shape(
        disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
        arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4, nuclear_amplitude=1.0)
    bare = base._replace(structure=base.structure._replace(nuclear_cluster_mass_solar=0.0))
    assert density.model_terms(bare)["nuclear_cluster_mass"] == 0
    assert density.model_terms(bare)["nuclear_disc_mass"] == 0
    assert density.model_terms(base)["nuclear_cluster_mass"] > 0


def test_a_hole_made_on_its_own_is_the_hole_of_a_random_galaxy():
    masses = {round(shape_draw.random_nucleus_mass_solar(), 0) for _ in range(5)}
    assert len(masses) > 1
    assert all(1.0e3 < m < 1.0e11 for m in masses)


def test_the_fixed_milky_way_nucleus_masses_are_gone():
    sources = "".join(p.read_text() for p in SOURCE.rglob("*.py"))
    for name in ("NUCLEAR_CLUSTER_BULGE_RATIO", "SUPERMASSIVE_BLACK_HOLE_MASS_SOLAR", "BLACK_HOLE_SUPERMASSIVE_VELOCITY_DISPERSION_KMS", "QUASAR_BLACK_HOLE_MASS_RANGE_SOLAR"):
        assert name not in sources, name
