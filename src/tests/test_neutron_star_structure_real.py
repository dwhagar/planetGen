# tests/test_neutron_star_structure_real.py

"""
GEN.968 on the 21 real equations of state of [R2834], against what NICER and the pulsar timing
measured. Skipped on a machine whose cache folder does not hold them (`python3 -m planetgen.cli.datasets fetch
neutron_star_matter`); `test_neutron_star_structure.py` covers the integration on equations of state built there.
"""

import math

import numpy as np
import pytest

from planetgen.datasets import fetch, neutron_star_matter
from planetgen.physics import neutron_star_structure as ns

pytestmark = [pytest.mark.real_datasets, pytest.mark.skipif(
    fetch.locate(neutron_star_matter.DATASET).missing, reason="the neutron-star tables are not in this machine's cache folder")]

INDICES = range(1, neutron_star_matter.EOS_COUNT + 1)


def test_the_committed_table_is_what_the_integration_gives():
    """`scripts/build_neutron_star_table.py --check` for three equations of state: the stiffest, a middle one and the softest."""
    for i in (1, 12, 6):
        assert ns.integrated_maximum_mass_solar(i) == pytest.approx(ns.maximum_mass_solar(i), abs=2e-4), i
        masses, radii = ns.integrated_curve(i)
        for mass in (0.8, 1.4, 1.9):
            assert ns.radius_km(i, mass) == pytest.approx(float(np.interp(mass, masses, radii)), abs=0.03), (i, mass)


def test_every_curve_rises_in_mass_and_no_star_is_denser_than_a_black_hole_or_a_buchdahl_star():
    for i in INDICES:
        masses, radii = ns.integrated_curve(i)
        assert (masses[1:] > masses[:-1]).all(), i
        assert (masses * ns.SOLAR_MASS_KM / radii < 4.0 / 9.0).all(), i
        assert masses[-1] == pytest.approx(ns.integrated_maximum_mass_solar(i), rel=1e-3)


def test_a_one_point_four_solar_mass_star_has_a_crust_about_a_kilometre_thick():
    for i in INDICES:
        eos = ns._eos(i)
        # bisect the central density for 1.4 solar masses
        low, high = 0.2, eos.density_max_fm3
        for _ in range(30):
            density = math.sqrt(low * high)
            star = ns.structure(i, math.exp(eos.log_pressure(density)) * ns.MEV_FM3_PER_KM2)
            low, high = (density, high) if star.mass_solar < 1.4 else (low, density)
        assert 0.4 < star.crust_thickness_km < 2.0, (i, star.crust_thickness_km)
        assert 0.0 < star.outer_crust_thickness_km < star.crust_thickness_km, i
        assert 0.0 < star.crust_mass_solar < 0.1 * star.mass_solar, i


def test_the_table_crust_is_the_integrated_crust():
    for i in (1, 12):
        eos = ns._eos(i)
        low, high = 0.2, eos.density_max_fm3
        for _ in range(30):
            density = math.sqrt(low * high)
            star = ns.structure(i, math.exp(eos.log_pressure(density)) * ns.MEV_FM3_PER_KM2)
            low, high = (density, high) if star.mass_solar < 1.4 else (low, density)
        assert ns.crust_thickness_km(i, star.mass_solar) == pytest.approx(star.crust_thickness_km, abs=0.05), i
