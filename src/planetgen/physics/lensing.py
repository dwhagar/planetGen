# planetgen/physics/lensing.py

"""
Gravitational Lensing
=====================

The Einstein radius of a compact object (GEN.154): the radius, at the lens,
of the ring a source straight behind it is smeared into,

    R_E = sqrt(4 G M / c^2 * D_L * D_LS / D_S)

for a lens at `D_L` from the observer and a source at `D_S`, `D_LS = D_S -
D_L` behind it. Derived from the mass and the distances; nothing is stored.
Pure functions, SI units.
"""

import math

from planetgen.physics import constants


def einstein_radius_m(mass_kg, lens_distance_m, source_distance_m):
    """
    The Einstein radius, metres, of a lens of `mass_kg` at `lens_distance_m`
    from the observer, against a source `source_distance_m` away on the same
    line. `None` when the lens is not between the observer and the source.
    """
    if mass_kg <= 0 or not 0 < lens_distance_m < source_distance_m:
        return None
    schwarzschild_scale_m = 4 * constants.G * mass_kg / constants.SPEED_OF_LIGHT_M_S ** 2
    return math.sqrt(schwarzschild_scale_m * lens_distance_m * (1 - lens_distance_m / source_distance_m))


def einstein_radius_against_galactic_centre_m(mass_kg, galactic_radius_pc):
    """
    The Einstein radius, metres, seen from Sol against a source at the
    galactic centre, for a lens `galactic_radius_pc` from the centre on the
    line between them (the lens-Sol distance is Sol's `GALACTIC_CENTER_DISTANCE_LY`
    less that radius). `None` for a lens farther from the centre than Sol is:
    nothing of it lies between us and the centre.
    """
    source_m = constants.GALACTIC_CENTER_DISTANCE_LY * constants.LY_TO_M
    lens_m = source_m - galactic_radius_pc * constants.PARSEC_M
    return einstein_radius_m(mass_kg, lens_m, source_m)
