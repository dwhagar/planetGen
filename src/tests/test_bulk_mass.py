# tests/test_bulk_mass.py

"""DB.55 / GEN.336: the bulk mass of a comet nucleus and of a gas cloud, worked out from size and density."""

import math

import pytest

from planetgen import tuning
from planetgen.physics import bulk_mass, constants


def test_a_comet_nucleus_is_its_sphere_at_the_density_of_67p():
    # 67P has the volume of a sphere 3.29 km across (18.7 km^3) and weighs 9.98e12 kg (Rosetta).
    assert bulk_mass.comet_nucleus_mass_kg(3.29) == pytest.approx(9.98e12, rel=0.02)
    assert bulk_mass.comet_nucleus_mass_kg(2.0) == pytest.approx(
        tuning.COMET_NUCLEUS_DENSITY_KG_M3 * 4.0 / 3.0 * math.pi * 1000.0 ** 3)


def test_the_mass_of_a_nucleus_grows_with_the_cube_of_its_size():
    assert bulk_mass.comet_nucleus_mass_kg(2.0) / bulk_mass.comet_nucleus_mass_kg(1.0) == pytest.approx(8.0)


def test_a_gas_cloud_is_its_sphere_at_its_density_with_the_helium_counted():
    # A molecular cloud of 100 hydrogen nuclei per cm^3 and 50 light-years radius holds some 1e4 to 1e5 solar masses.
    solar = bulk_mass.gas_cloud_mass_kg(100.0, 50.0) / constants.SOLAR_MASS_TO_KG
    assert 1.0e4 < solar < 1.0e5
    assert bulk_mass.gas_cloud_mass_kg(200.0, 50.0) == pytest.approx(2.0 * bulk_mass.gas_cloud_mass_kg(100.0, 50.0))
    assert bulk_mass.gas_cloud_mass_kg(100.0, 100.0) == pytest.approx(8.0 * bulk_mass.gas_cloud_mass_kg(100.0, 50.0))
