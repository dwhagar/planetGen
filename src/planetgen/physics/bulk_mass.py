# planetgen/physics/bulk_mass.py

"""
The mass of the free objects whose mass is the bulk of what they hold (GEN.336): a comet's nucleus and a cloud of gas.

Both are worked out from the object's own size and density, so the gravity sum has a mass for every object rather than a
stated zero. A cloud is taken as its bounding sphere at its mean particle density, which is an upper bound for a nebula
that fills part of its sphere; GEN.295 refines it when a nebula becomes a group of points.
"""

import math

from planetgen import tuning
from planetgen.physics import constants


def comet_nucleus_mass_kg(diameter_km):
    """kg: a sphere of `diameter_km` at `tuning.COMET_NUCLEUS_DENSITY_KG_M3`."""
    radius_m = diameter_km * 500.0
    return 4.0 / 3.0 * math.pi * radius_m ** 3 * tuning.COMET_NUCLEUS_DENSITY_KG_M3


def gas_cloud_mass_kg(density_cm3, radius_ly):
    """kg: a sphere of `radius_ly` holding `density_cm3` hydrogen nuclei per cm^3, with the helium and the rest of the
    interstellar gas counted by `constants.ISM_MASS_PER_HYDROGEN`."""
    volume_m3 = 4.0 / 3.0 * math.pi * (radius_ly * constants.LY_TO_M) ** 3
    return volume_m3 * density_cm3 * constants.CM3_TO_M3 * constants.HYDROGEN_ATOM_MASS_KG * constants.ISM_MASS_PER_HYDROGEN
