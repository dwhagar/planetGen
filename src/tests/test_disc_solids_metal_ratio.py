# tests/test_disc_solids_metal_ratio.py

"""GEN.464: the disc's solids scale with the star's metal mass fraction over the Sun-like reference's (Z / Z0), not with 10^[Fe/H]."""

import numpy as np
import pytest

from planetgen.physics import bulk_draw, constants, formation, stellar_abundances as sa


def _ratio(fe_h, branch):
    offsets = sa.abundance_offsets(fe_h, branch, np.zeros(5), 1.0, 0.0)
    return sa.metal_ratio_to_solar(sa.composition(fe_h, offsets).z)


@pytest.mark.parametrize("fe_h", [0.0, -0.5, -1.0, -2.0])
def test_a_thin_disc_star_has_about_ten_to_the_fe_h_and_a_high_alpha_star_has_one_and_a_half_to_three_times_more(fe_h):
    thin = _ratio(fe_h, sa.BRANCH_THIN)
    high = _ratio(fe_h, sa.BRANCH_HIGH_ALPHA)
    # the thin disc's oxygen and magnesium rise a little as [Fe/H] falls, the high-alpha branch's much more (GEN.464's table)
    assert 0.95 <= thin / 10.0 ** fe_h <= (1.02 if fe_h == 0.0 else 2.0)
    if fe_h < 0.0:
        assert 1.7 <= high / 10.0 ** fe_h <= 3.1
        assert high > thin


def test_the_solar_reference_is_one():
    assert sa.metal_ratio_to_solar(sa.composition(0.0, {}).z) == pytest.approx(1.0)


def test_a_remnant_without_a_metallicity_scales_by_one():
    class Remnant:
        metallicity_fe_h = None
        age = 5.0

    assert bulk_draw.metal_ratio_to_solar(Remnant()) == 1.0


def test_the_disc_scale_follows_the_ratio_and_is_deterministic():
    mass = constants.SOLAR_MASS_TO_KG
    assert formation.disk_surface_density_scale(mass, 0.3) == pytest.approx(0.3 * formation.disk_surface_density_scale(mass))
    assert formation.disk_surface_density_scale(mass, 0.3) == formation.disk_surface_density_scale(mass, 0.3)


def test_a_metal_poor_star_gets_fewer_planet_building_solids_than_a_solar_one():
    class Star:
        age = 10.0
        mass = constants.SOLAR_MASS_TO_KG

        def __init__(self, fe_h):
            self.metallicity_fe_h = fe_h

    assert bulk_draw.metal_ratio_to_solar(Star(-1.0)) < 0.5 * bulk_draw.metal_ratio_to_solar(Star(0.0))
