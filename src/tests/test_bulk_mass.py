# tests/test_bulk_mass.py

"""DB.55 / GEN.336: the bulk mass of a gas cloud, worked out from size and density (a comet nucleus takes its mass from its makeup, GEN.442)."""

import pytest

from planetgen.physics import bulk_mass, constants


def test_a_gas_cloud_is_its_sphere_at_its_density_with_the_helium_counted():
    # A molecular cloud of 100 hydrogen nuclei per cm^3 and 50 light-years radius holds some 1e4 to 1e5 solar masses.
    solar = bulk_mass.gas_cloud_mass_kg(100.0, 50.0) / constants.SOLAR_MASS_TO_KG
    assert 1.0e4 < solar < 1.0e5
    assert bulk_mass.gas_cloud_mass_kg(200.0, 50.0) == pytest.approx(2.0 * bulk_mass.gas_cloud_mass_kg(100.0, 50.0))
    assert bulk_mass.gas_cloud_mass_kg(100.0, 100.0) == pytest.approx(8.0 * bulk_mass.gas_cloud_mass_kg(100.0, 50.0))
