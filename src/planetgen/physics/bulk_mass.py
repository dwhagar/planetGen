# planetgen/physics/bulk_mass.py

"""
The mass of the free objects whose mass is the bulk of what they hold (GEN.336): a cloud of gas (a comet's nucleus takes its mass from
its makeup, `small_body_composition.SmallBodyMakeup.mass_kg`, GEN.442).

It is worked out from the object's own size and density, so the gravity sum has a mass for every object rather than a
stated zero. A cloud is taken as its bounding sphere at its mean particle density, which is an upper bound for a nebula
that fills part of its sphere; GEN.295 refines it when a nebula becomes a group of points.
"""

from planetgen.physics import constants, relations


def gas_cloud_mass_kg(density_cm3, radius_ly):
    """kg: a sphere of `radius_ly` holding `density_cm3` hydrogen nuclei per cm^3, with the helium and the rest of the
    interstellar gas counted by `constants.ISM_MASS_PER_HYDROGEN`."""
    volume_m3 = relations.sphere_volume(radius_ly * constants.LY_TO_M)
    return volume_m3 * density_cm3 * constants.CM3_TO_M3 * constants.HYDROGEN_ATOM_MASS_KG * constants.ISM_MASS_PER_HYDROGEN
