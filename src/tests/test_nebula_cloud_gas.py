# tests/test_nebula_cloud_gas.py

"""
GEN.1273: a nebula inside a Magellanic Cloud takes its gas [O/H] from the Cloud's own metallicity (`ComponentField.gas_fe_h`, `nebula_composition.SATELLITE_GAS_OFFSETS`), so it is more
metal-poor than one at the solar circle; one of the Milky Way keeps the gradient.
"""

import pytest

from planetgen.galaxy import seed as galaxy_seed
from planetgen.physics import nebula_columns as nc
from planetgen.physics import nebula_composition
from tests.test_gen_magellanic_classes import GALAXY_SEED, LMC, SHAPE, SMC, SUN, WITH_SATELLITES

DISK = (8128.0, 0.0, 0.0)


def _row(center, index, galactic_radius=None):
    x, y, z = center
    return {"nebula_class": "H", "density_cm3": 100.0, "radius_ly": 5.0, "uid": index.to_bytes(12, "big"), "name": "n", "hosts": None,
            "galactic_radius_pc": (x * x + y * y + z * z) ** 0.5 if galactic_radius is None else galactic_radius,
            "center_x_pc": x, "center_y_pc": y, "center_z_pc": z}


def _oh(center, index, shape=WITH_SATELLITES):
    with galaxy_seed.seeded(GALAXY_SEED, "test", "gen1273", shape=shape):
        return nc.nebula_values(_row(center, index))["gas_oh_dex"]


def test_a_cloud_reports_its_own_metallicity_and_the_solar_circle_none():
    assert WITH_SATELLITES.components.gas_fe_h(tuple(LMC.center_pc)) == pytest.approx(LMC.metallicity_fe_h)
    assert WITH_SATELLITES.components.gas_fe_h(tuple(SMC.center_pc)) == pytest.approx(SMC.metallicity_fe_h)
    assert WITH_SATELLITES.components.gas_fe_h(SUN) is None


def test_a_cloud_nebula_is_drawn_inside_the_widest_range_of_the_two_methods_above_the_clouds_metallicity():
    for cloud in (LMC, SMC):
        low, high = nebula_composition.satellite_offsets(cloud.metallicity_fe_h)
        found = [_oh(tuple(cloud.center_pc), index) for index in range(40)]
        assert all(cloud.metallicity_fe_h + low <= value <= cloud.metallicity_fe_h + high for value in found)
        assert max(found) - min(found) > 0.5 * (high - low)        # drawn, not one value


def test_a_cloud_nebula_is_more_metal_poor_than_one_at_the_solar_circle():
    """On average, as the scatter of the solar circle's gas (a few tenths of a dex) overlaps the Large Cloud's range."""
    solar = sum(_oh(DISK, index) for index in range(200)) / 200
    for cloud in (LMC, SMC):
        assert sum(_oh(tuple(cloud.center_pc), index) for index in range(200)) / 200 < solar
    assert sum(_oh(tuple(SMC.center_pc), index) for index in range(200)) / 200 < sum(_oh(tuple(LMC.center_pc), index) for index in range(200)) / 200


def test_an_object_with_no_satellite_around_it_keeps_the_milky_way_gradient():
    """Without the galaxy bound, with no satellites, or with no center, the number is the gradient's, as before GEN.1273."""
    row = _row(tuple(LMC.center_pc), 7)
    key = nc.row_key("nebulae", row)
    assert nc.nebula_values(row)["gas_oh_dex"] == nebula_composition.gas_oh(row["galactic_radius_pc"], key)
    for shape, made in ((SHAPE, row), (WITH_SATELLITES, {k: v for k, v in row.items() if not k.startswith("center_")})):
        with galaxy_seed.seeded(GALAXY_SEED, "test", "gen1273", shape=shape):
            assert nc.nebula_values(made)["gas_oh_dex"] == nebula_composition.gas_oh(row["galactic_radius_pc"], nc.row_key("nebulae", made))
