# tests/test_nuclear_disc_law.py

"""
GEN.628: the nuclear stellar disc's mass is drawn from the galaxy's stellar mass by the two fits of [R1487], not as a fixed multiple of the nuclear cluster, and only a barred galaxy has one.
"""

import math
import statistics

import pytest

from planetgen.galaxy import density, structure_draw as sd


def _seeds(count):
    return [bytes([i % 256, i // 256]) + bytes(14) for i in range(count)]


def test_the_two_fits_are_the_papers_equations():
    photometric, spectroscopic = sd.nuclear_disc_mass_bounds_solar(1.0e9)
    assert photometric == pytest.approx(10 ** 7.62) and spectroscopic == pytest.approx(10 ** 7.99)
    photometric, spectroscopic = sd.nuclear_disc_mass_bounds_solar(1.0e11)
    assert photometric == pytest.approx(10 ** (7.62 + 2 * 1.32)) and spectroscopic == pytest.approx(10 ** (7.99 + 2 * 0.86))


def test_the_same_seed_draws_the_same_disc_and_two_seeds_draw_two():
    first = sd.draw_nuclear_disc_mass_solar(bytes(16), 5.0e10)
    assert first == sd.draw_nuclear_disc_mass_solar(bytes(16), 5.0e10)
    assert len({sd.draw_nuclear_disc_mass_solar(seed, 5.0e10) for seed in _seeds(20)}) == 20


@pytest.mark.parametrize("mass", [3.0e9, 5.0e10, 2.0e11])
def test_over_many_seeds_the_disc_follows_the_law_with_its_scatter(mass):
    drawn = [math.log10(sd.draw_nuclear_disc_mass_solar(seed, mass)) for seed in _seeds(600)]
    low, high = sorted(math.log10(m) for m in sd.nuclear_disc_mass_bounds_solar(mass))
    assert low - 3 * sd.NUCLEAR_DISC_SCATTER_DEX <= min(drawn) and max(drawn) <= high + 3 * sd.NUCLEAR_DISC_SCATTER_DEX
    assert statistics.mean(drawn) == pytest.approx((low + high) / 2.0, abs=0.1)
    assert 0.8 * sd.NUCLEAR_DISC_SCATTER_DEX < statistics.pstdev(drawn) < 1.4 * sd.NUCLEAR_DISC_SCATTER_DEX + (high - low) / 2.0


def test_a_milky_way_mass_galaxy_lands_inside_the_range_the_fits_and_the_milky_ways_own_disc_span():
    median = statistics.median(sd.draw_nuclear_disc_mass_solar(seed, 5.0e10) for seed in _seeds(600))
    assert 1.0e9 < median < 1.0e10          # the measured 1e9 to 1.5e9 below it, the photometric fit's 7e9 above


def test_only_a_barred_galaxy_has_a_disc_and_it_no_longer_follows_the_cluster():
    structures = [sd.draw_structure(seed, 5.0e10) for seed in _seeds(200)]
    assert all((s.nuclear_disc_mass_solar > 0) == s.bar_present for s in structures)
    assert any(s.bar_present for s in structures) and any(not s.bar_present for s in structures)
    assert any(s.nuclear_disc_mass_solar > 0 and s.nuclear_cluster_mass_solar == 0 for s in structures)
    assert any(s.nuclear_disc_mass_solar == 0 and s.nuclear_cluster_mass_solar > 0 for s in structures)


def test_a_small_galaxys_disc_is_a_small_share_of_its_stars_not_a_fifth():
    shares = [s.nuclear_disc_mass_solar / s.stellar_mass_solar for s in (sd.draw_structure(seed, 4.0e9) for seed in _seeds(300))
              if s.nuclear_disc_mass_solar]
    assert shares and statistics.median(shares) < 0.15


def test_the_density_model_reads_the_drawn_disc_mass():
    shape = density.build_galaxy_shape(
        disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
        arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4, nuclear_amplitude=1.0)
    doubled = shape._replace(structure=shape.structure._replace(nuclear_disc_mass_solar=2 * shape.structure.nuclear_disc_mass_solar))
    assert density.model_terms(doubled)["nuclear_disc_mass"] == pytest.approx(2 * density.model_terms(shape)["nuclear_disc_mass"])
