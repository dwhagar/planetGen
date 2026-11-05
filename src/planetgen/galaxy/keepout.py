# planetgen/galaxy/keepout.py

"""
The keep-out radius of an object (NAV.24, NAV.54): how close a course should come
before it counts as entering the object's gravity. One rule per kind, the
data for it read by `planetgen.db.query.keep_out_radius`
(`docs/design/orbits/course-avoidance.md` section 2):

- planet, moon: its stored Hill radius (a moon's is against its planet's mass, GEN.138);
- star, system: the stored `system_perimeter_km` (the wider of the pair
  for a binary), at galaxy and sector scale; inside its own system a star
  keeps out its radiation radius (`star_in_system_keep_out`);
- black hole (stellar and intermediate), neutron star, rogue planet: the
  Hill radius in the galaxy's field (`galactic_hill_radius_km`), never
  smaller than the object itself;
- supermassive black hole: its sphere of influence `G M / sigma^2`
  (`sphere_of_influence_km`), which the galactic Hill radius (zero at the
  centre) cannot give; a quasar the larger of that and its radiation radius;
- asteroid field: its stored bounding radius (NAV.51: courses avoid asteroid fields);
- nebula, supernova remnant: none, with a note: the gas is thin and the
  drawn cloud fills a fraction of its bounding sphere, so a course passes
  through it;
- asteroid belt, comet, interstellar comet: none (too small or too
  diffuse to steer around), a belt with a note.
"""

import math
from collections import namedtuple

from planetgen import tuning
from planetgen.physics import constants
from planetgen.physics.orbits import calculate_hill_sphere
from planetgen.physics import unit_registry as units

_pc_to_m = units.converter("pc", "m")
_m_to_km = units.converter("m", "km")
_ly_to_km = units.converter("ly", "km")

KeepOut = namedtuple("KeepOut", ["radius_km", "basis", "note"])
"""
A keep-out radius: `radius_km` (`None` for none), `basis` (`"hill"`,
`"perimeter"`, `"galactic_hill"`, `"influence"`, `"radiation"`, `"radius"` or `"none"`) and `note` (what
the NAV page says about it, or `None`).
"""

PASS_THROUGH_NOTE = "It is thin gas and dust, so a course passes through it."
"""str: The note on a nebula or a supernova remnant (NAV.54: the drawn cloud fills 13% of its bounding sphere and the gas is thin, so there is
no hard keep-out)."""

BELT_NOTE = "An asteroid belt is thin, so a course crosses it without a close approach."
"""str: The note on an asteroid belt."""

NO_KEEP_OUT = KeepOut(None, "none", None)


def galactic_hill_radius_km(mass_kg, galactic_radius_pc, central_mass_kg):
    """
    The Hill radius, in km, of an object of `mass_kg` that orbits the
    galaxy `galactic_radius_pc` parsecs from its center, with
    `central_mass_kg` (the galaxy's mass inside that radius,
    `galaxy.enclosed_mass`) as the central mass, as for a star's own
    `calculate_hill_sphere`.
    """
    distance_m = _pc_to_m(galactic_radius_pc)
    return _m_to_km(calculate_hill_sphere(distance_m, mass_kg, central_mass_kg))


def compact_keep_out(mass_kg, galactic_radius_pc, radius_km, central_mass_kg):
    """
    The keep-out of a black hole, neutron star, quasar or rogue planet:
    its galactic Hill radius, but no smaller than `radius_km` (the object
    itself).
    """
    hill = galactic_hill_radius_km(mass_kg, galactic_radius_pc, central_mass_kg)
    if hill >= radius_km:
        return KeepOut(hill, "galactic_hill", None)
    return KeepOut(radius_km, "radius", None)


def sphere_of_influence_km(mass_solar):
    """
    The sphere of influence `G M / sigma^2`, km, of a supermassive black hole of `mass_solar` solar masses, `sigma` the bulge dispersion its
    mass sets (`galaxy.structure.black_hole_dispersion_kms`, the M-sigma relation the generator drew the hole with): 3.1 pc at 4e6 solar
    masses. A course should not graze the event horizon (0.08 AU there) and the galactic Hill radius is zero at the centre.
    """
    from planetgen.galaxy import structure     # not at the top: the galaxy package imports this module first
    sigma_m_s = structure.black_hole_dispersion_kms(mass_solar) * 1000.0
    return _m_to_km(constants.G * mass_solar * constants.SOLAR_MASS_TO_KG / sigma_m_s ** 2)


def radiation_radius_km(luminosity_w, flux_w_m2):
    """The distance, km, at which a source of `luminosity_w` gives `flux_w_m2`: `sqrt(L / (4 pi F))`."""
    return _m_to_km(math.sqrt(luminosity_w / (4.0 * math.pi * flux_w_m2)))


def supermassive_keep_out(mass_solar):
    """The keep-out of a supermassive black hole: its sphere of influence."""
    return KeepOut(sphere_of_influence_km(mass_solar), "influence", None)


def quasar_keep_out(mass_solar, luminosity_w):
    """The keep-out of a quasar nucleus: the larger of its sphere of influence and the radius at which its light is
    `tuning.KEEP_OUT_QUASAR_FLUX_W_M2`."""
    influence = sphere_of_influence_km(mass_solar)
    light = radiation_radius_km(luminosity_w, tuning.KEEP_OUT_QUASAR_FLUX_W_M2)
    return KeepOut(influence, "influence", None) if influence >= light else KeepOut(light, "radiation", None)


def star_in_system_keep_out(radius_km, luminosity_w):
    """The keep-out of a star inside its own system: `max(10 R_star, sqrt(L / (4 pi F_lim)))` (`tuning.KEEP_OUT_STAR_FLOOR_RADII`,
    `KEEP_OUT_STAR_FLUX_W_M2`). Hard, and a course with an end inside it is exempt."""
    floor = tuning.KEEP_OUT_STAR_FLOOR_RADII * radius_km
    light = radiation_radius_km(luminosity_w, tuning.KEEP_OUT_STAR_FLUX_W_M2)
    return KeepOut(floor, "radius", None) if floor >= light else KeepOut(light, "radiation", None)


def asteroid_field_keep_out(radius_ly):
    """The keep-out of an asteroid field: its stored bounding radius (the whole sphere, no density-based shrinking)."""
    return KeepOut(_ly_to_km(radius_ly), "radius", None)
