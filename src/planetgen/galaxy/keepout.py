# planetgen/galaxy/keepout.py

"""
The keep-out radius of an object (NAV.24): how close a course should come
before it counts as entering the object's gravity. One rule per kind, the
data for it read by `queryDb.keep_out_radius`:

- planet, moon: its stored Hill radius;
- star, system: the stored `system_perimeter_km` (the wider of the pair
  for a binary);
- black hole, neutron star, quasar, rogue planet: the Hill radius in the
  galaxy's field (`galactic_hill_radius_km`), never smaller than the
  object itself (a quasar sits at the galactic center, where the Hill
  radius is zero);
- nebula, supernova remnant, asteroid field: none -- they have no mass
  stored, so a course passes through them, with a note;
- asteroid belt, comet, interstellar comet: none (too small or too
  diffuse to steer around).
"""

from collections import namedtuple

from planetgen.physics import constants
from planetgen.physics.orbits import calculate_hill_sphere

KeepOut = namedtuple("KeepOut", ["radius_km", "basis", "note"])
"""
A keep-out radius: `radius_km` (`None` for none), `basis` (`"hill"`,
`"perimeter"`, `"galactic_hill"`, `"radius"` or `"none"`) and `note` (what
the NAV page says about it, or `None`).
"""

PASS_THROUGH_NOTE = "No mass is stored for it, so the course passes through it."
"""str: The note on a nebula, remnant or asteroid field."""

NO_KEEP_OUT = KeepOut(None, "none", None)


def galactic_hill_radius_km(mass_kg, galactic_radius_pc):
    """
    The Hill radius, in km, of an object of `mass_kg` that orbits the
    galaxy `galactic_radius_pc` parsecs from its center
    (`MILKY_WAY_MASS` as the central mass, as for a star's own
    `calculate_hill_sphere`).
    """
    distance_m = galactic_radius_pc * constants.PARSEC_M
    return calculate_hill_sphere(distance_m, mass_kg, constants.MILKY_WAY_MASS) / 1000.0


def compact_keep_out(mass_kg, galactic_radius_pc, radius_km):
    """
    The keep-out of a black hole, neutron star, quasar or rogue planet:
    its galactic Hill radius, but no smaller than `radius_km` (the object
    itself).
    """
    hill = galactic_hill_radius_km(mass_kg, galactic_radius_pc)
    if hill >= radius_km:
        return KeepOut(hill, "galactic_hill", None)
    return KeepOut(radius_km, "radius", None)
