# tests/test_gen_sector_profile.py

"""
GEN.235: the sector profile (`planetgen/galaxy/profile.py`). It is a pure function of the galaxy shape, the sector edge
and the address: the same answer every call, no database, the figures of docs/design/galaxy/intergalactic-and-halo-contents.md
(5.3 and 5.7) at the Sun, in the halo and at the edge of the reach. The engine paths that read it (an asked cell in the
halo is generated and saved, the neighbourhood pull, the fill) are in `test_galaxy_gen.py`.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import bright_stars as brightStars
from planetgen.galaxy import profile, spherical
from planetgen.galaxy.density import build_galaxy_shape, mean_system_mass_solar, relative_density
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.physics.units import ly_to_pc

pytestmark = pytest.mark.usefixtures("large_spiral_fallback")

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
E_VALUE = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
SUN = (8128.0, 0.0, 0.0)


def _at(position, shape=SHAPE):
    return profile.profile_at(position, shape, EDGE_PC, E_VALUE)


def test_a_profile_is_a_pure_function_of_the_place():
    assert _at(SUN) == _at(SUN)
    assert _at((50000.0, 3000.0, 20000.0)) == _at((50000.0, 3000.0, 20000.0))
    skeleton = SimpleNamespace(shape=SHAPE, edge_pc=EDGE_PC, expected_system_count_at_density_1=E_VALUE, galaxy_seed=bytes(range(16)))
    assert profile.sector_profile((30, 2, 5), skeleton) == profile.sector_profile((30, 2, 5), skeleton)


def test_density_is_the_relative_density_and_lambda_is_e_times_it():
    for position in (SUN, (0.0, 0.0, 0.0), (10000.0, 0.0, 10000.0), (50000.0, 0.0, 0.0)):
        found = _at(position)
        assert found.density == max(relative_density(position, SHAPE), 0.0)
        assert found.expected_systems == pytest.approx(E_VALUE * found.density)


def test_the_mix_adds_up_to_the_density_and_the_shares_to_one():
    for position in (SUN, (0.0, 0.0, 0.0), (10000.0, 0.0, 10000.0), (262000.0, 0.0, 0.0)):
        found = _at(position)
        assert set(found.mix) == set(profile.POPULATIONS) == set(found.populations)
        assert sum(found.mix.values()) == pytest.approx(found.density, rel=1e-9)
        assert sum(p.share for p in found.populations.values()) == pytest.approx(1.0)
        assert all(found.populations[name].density == found.mix[name] for name in found.mix)


@pytest.mark.parametrize("position, region", [
    (SUN, "disk"), ((0.0, 0.0, 0.0), "bulge"), ((50000.0, 0.0, 0.0), "halo"), ((10000.0, 0.0, 10000.0), "halo"),
])
def test_the_region_is_the_term_that_dominates(position, region):
    assert _at(position).region == region


def test_the_halo_cells_expect_the_systems_the_design_note_gives():
    # section 5.7: 3.2e-5 at 50 kpc and 1.8e-8 at the 262 kpc edge, for E = 6.3 systems at density 1 [C], worked out at the
    # recalled 0.4 solar masses a system; the derived mean (GEN.545) scales the galaxy's systems and so these.
    assert E_VALUE == pytest.approx(6.3, rel=0.01)
    scale = 0.4 / mean_system_mass_solar(SHAPE)
    assert _at((50000.0, 0.0, 0.0)).expected_systems == pytest.approx(3.2e-5 * scale, rel=0.15)
    assert _at((262000.0, 0.0, 0.0)).expected_systems == pytest.approx(1.8e-8 * scale, rel=0.15)
    far = [_at((r, 0.0, 0.0)).expected_systems for r in (30000.0, 50000.0, 100000.0, 262000.0)]
    assert far == sorted(far, reverse=True) and all(value > 0.0 for value in far)


def test_the_populations_carry_the_ages_metallicities_and_motions_of_the_note():
    found = _at(SUN)
    young, old, bulge, halo = (found.populations[name] for name in ("young", "old", "bulge", "halo"))
    assert young.age_gy == (0.0, 0.1) and old.age_gy == (3.0, 13.0) and halo.age_gy == (10.0, 13.0)
    assert bulge.rotation == "bar" and halo.rotation == "none" and young.rotation == "disk"
    assert young.dispersion_kms < old.dispersion_kms < bulge.dispersion_kms
    (weight, mean, sigma), = young.fe_h
    assert weight == 1.0 and mean == pytest.approx(0.0, abs=0.2) and sigma > 0.0
    # the halo is -1.6 inside 25 kpc and -2.2 outside, wherever the cell is
    assert _at((10000.0, 0.0, 5000.0)).populations["halo"].fe_h[0][1] == -1.6
    assert _at((60000.0, 0.0, 0.0)).populations["halo"].fe_h[0][1] == -2.2
    # a cluster has two peaks, and no density until the component table (GEN.237)
    peaks = found.populations["cluster"].fe_h
    assert [peak[1] for peak in peaks] == [-1.5, -0.5] and sum(peak[0] for peak in peaks) == pytest.approx(1.0)
    assert all(_at(position).mix["cluster"] == 0.0 for position in (SUN, (50000.0, 0.0, 0.0)))


def test_the_dark_matter_follows_the_nfw_halo_of_the_potential_note():
    # section 2.6: 0.0091 solar masses per pc3 at the Sun, 2.0e-4 at 50 kpc, 3.2e-5 at 100 kpc, 2.1e-6 at 262 kpc
    for radius, expected in ((8128.0, 0.0091), (50000.0, 2.0e-4), (100000.0, 3.2e-5), (262000.0, 2.1e-6)):
        found = _at((radius, 0.0, 0.0))
        assert found.dark_matter.density_solar_pc3 == pytest.approx(expected, rel=0.1)
        assert found.dark_matter.mass_solar == pytest.approx(found.dark_matter.density_solar_pc3 * found.volume_pc3)
    assert math.isfinite(_at((0.0, 0.0, 0.0)).dark_matter.density_solar_pc3)


def test_a_galaxy_twice_the_size_has_a_dark_matter_halo_twice_the_radius():
    big = SHAPE._replace(disk_scale_length_pc=2 * SHAPE.disk_scale_length_pc)
    for radius in (4000.0, 30000.0):
        # the same place in units of the halo's radius: the mass doubles and the volume grows eightfold
        assert profile.dark_matter_density(2 * radius, big) == pytest.approx(
            profile.dark_matter_density(radius, SHAPE) / 4.0, rel=1e-9)


def test_the_gas_is_cold_in_the_disk_and_hot_or_intergalactic_outside_it():
    disk = _at(SUN).gas
    assert 0.0 < disk.cold_fraction <= 1.0 and disk.ambient == "disk"
    assert disk.ambient_density_cm3 == pytest.approx(tuning.DISK_ISM_DENSITY_CM3 * disk.cold_fraction)
    halo = _at((50000.0, 0.0, 0.0)).gas
    assert halo.cold_fraction < 1e-6 and halo.ambient == "hot-halo"
    assert halo.ambient_density_cm3 == pytest.approx(1.0e-4, rel=1e-6)        # the note's 1e-4 per cm3 at 50 kpc
    assert _at((262000.0, 0.0, 0.0)).gas.ambient == "hot-halo"
    # far past the reach the hot halo has thinned below the Local Group's own gas
    out = _at((3.0e7, 0.0, 0.0)).gas
    assert out.ambient == "intergalactic"
    assert out.ambient_density_cm3 == pytest.approx(tuning.COSMIC_MEAN_DENSITY_CM3 * tuning.LOCAL_GROUP_GAS_OVERDENSITY)
    assert halo.mass_solar == pytest.approx(
        halo.ambient_density_cm3 * tuning.GAS_SOLAR_MASS_PER_PC3_PER_CM3 * _at((50000.0, 0.0, 0.0)).volume_pc3)


def test_a_cell_profile_records_its_address_and_the_exact_cell_volume():
    skeleton = SimpleNamespace(shape=SHAPE, edge_pc=EDGE_PC, expected_system_count_at_density_1=E_VALUE, galaxy_seed=bytes(range(16)))
    address = (500, 3, 7)
    found = profile.sector_profile(address, skeleton)
    assert found.address == address
    assert found.position_pc == pytest.approx(spherical.sector_position_pc(*address, EDGE_PC))
    assert found.volume_pc3 == pytest.approx(spherical.cell_volume(*address, EDGE_PC))
    with pytest.raises(ValueError):
        profile.sector_profile((2, 99, 0), skeleton)


def test_a_fill_context_reads_the_mix_and_the_centre_from_the_profile():
    found = _at(SUN)
    fill = brightStars.FillContext(found, min_luminosity_sol=2500.0)
    assert fill.profile is found and fill.center_pc == found.position_pc and fill.shape is SHAPE
    assert fill.densities == found.mix
    assert 0.0 <= fill.bright_share() < 1.0


def test_a_neighbourhood_pull_stays_in_the_outline_from_inside_it_and_takes_the_reach_from_outside():
    from planetgen.galaxy.skeleton import GalaxyBounds
    from planetgen.generation import run_galaxy
    bounds = GalaxyBounds([(0, 0), (1, 1), (2, 2), (3, 3), (4, 3), (5, 2), (6, 1)], EDGE_PC, reach_shell_index=100)
    inside = spherical.sector_position_pc(3, 0, 0, EDGE_PC)
    outside = spherical.sector_position_pc(40, 0, 0, EDGE_PC)
    assert run_galaxy._is_far(inside, EDGE_PC, bounds) is False
    assert run_galaxy._is_far(outside, EDGE_PC, bounds) is True
    # a halo pull is bounded by the reach, not by the outline
    assert run_galaxy._ReachOnly(bounds).qualifies(40, 0) and not bounds.qualifies(40, 0)
    assert not run_galaxy._ReachOnly(bounds).qualifies(101, 0)
