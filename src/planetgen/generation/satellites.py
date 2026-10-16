# planetgen/generation/satellites.py

"""
The Milky Way's satellite galaxies as galaxy components (GEN.237; docs/design/intergalactic-and-halo-contents.md 5.5).

The Large and Small Magellanic Clouds, the Sagittarius dwarf and eight of the classical dwarf spheroidals are real
objects at real places, so they are catalogue rows, not draws: `MILKY_WAY_SATELLITES` lists each with its galactic
longitude and latitude, its distance from the Sun, its stellar mass, half-light radius, mean [Fe/H] and age.
`place_satellites` turns the table into `Component` rows in the galaxy frame, putting the Sun where the density model
does (`galaxy.density.model_terms`: the solar radius at the inter-arm azimuth), which is the one anchor that makes a
direction on the sky a direction in the galaxy.

Every value in the table is recalled from McConnachie (2012, AJ 144, 4) and not checked against it (the sources were out
of reach when GEN.237 was built), so each is good to a few tens of percent in size and distance and a factor of two in
mass until Boss supplies the catalogue; they are filed with TODO as needing verification. The half-light radius stands
in for the half-mass radius (the stars' mass follows their light to this accuracy). The velocities are zero: the sector
gate reads only the position, and the satellites' proper motions (which move them by a few hundred parsecs over the
age of a game) belong with the encounter model of the orbital update.
"""

import collections
import math

from planetgen.galaxy import components
from planetgen.physics import component_profiles
from planetgen.util import draw

Satellite = collections.namedtuple(
    "Satellite", "name longitude_deg latitude_deg distance_kpc stellar_mass_sol half_light_radius_pc profile fe_h age_gy")
"""One real satellite: where it lies on the sky and how far, how heavy and how big, and its stars' [Fe/H] and age."""

MILKY_WAY_SATELLITES = (
    Satellite("Large Magellanic Cloud", 280.47, -32.89, 49.6, 2.5e9, 2700.0, "exponential", -0.4, 5.0),
    Satellite("Small Magellanic Cloud", 302.8, -44.3, 62.4, 5.2e8, 1100.0, "exponential", -0.9, 6.0),
    Satellite("Sagittarius dwarf", 5.6, -14.2, 24.0, 1.0e8, 2600.0, "exponential", -0.6, 8.0),
    Satellite("Fornax dwarf", 237.1, -65.7, 143.0, 2.0e7, 710.0, "plummer", -1.0, 7.0),
    Satellite("Sculptor dwarf", 287.5, -83.2, 86.0, 2.3e6, 283.0, "plummer", -1.7, 11.0),
    Satellite("Carina dwarf", 260.1, -22.2, 105.0, 3.8e5, 250.0, "plummer", -1.7, 7.0),
    Satellite("Draco dwarf", 86.4, 34.7, 76.0, 2.9e5, 221.0, "plummer", -2.0, 11.0),
    Satellite("Ursa Minor dwarf", 105.0, 44.8, 76.0, 2.9e5, 405.0, "plummer", -2.1, 12.0),
    Satellite("Sextans dwarf", 243.5, 42.3, 86.0, 4.4e5, 695.0, "plummer", -1.9, 11.0),
    Satellite("Leo I", 226.0, 49.1, 254.0, 5.5e6, 251.0, "plummer", -1.4, 5.0),
    Satellite("Leo II", 220.2, 67.2, 233.0, 7.4e5, 176.0, "plummer", -1.6, 9.0),
)
"""tuple: The Milky Way's satellites (recalled from McConnachie 2012, unchecked: see the module docstring)."""


def sun_frame(solar_radius_pc, solar_angle_rad):
    """`(sun, ex, ey, ez)`: the Sun's position in the galaxy frame and the unit vectors of the galactic coordinates there.

    `ex` points from the Sun to the galactic centre (longitude 0, latitude 0), `ey` towards longitude 90 degrees (the
    direction the disk turns) and `ez` along the galaxy frame's +z. Those three form a left-handed set, so a proper
    galactic frame has its north pole towards -z: `satellite_position_pc` takes a latitude `b` as a height of
    `-sin(b)` in the galaxy frame.
    """
    cos_a, sin_a = math.cos(solar_angle_rad), math.sin(solar_angle_rad)
    sun = (solar_radius_pc * cos_a, solar_radius_pc * sin_a, 0.0)
    return sun, (-cos_a, -sin_a, 0.0), (-sin_a, cos_a, 0.0), (0.0, 0.0, 1.0)


def satellite_position_pc(longitude_deg, latitude_deg, distance_pc, solar_radius_pc, solar_angle_rad):
    """The position in the galaxy frame (pc, the galactic centre at the origin) of an object seen from the Sun at this
    galactic longitude and latitude and this distance."""
    sun, ex, ey, ez = sun_frame(solar_radius_pc, solar_angle_rad)
    longitude, latitude = math.radians(longitude_deg), -math.radians(latitude_deg)
    along = (math.cos(latitude) * math.cos(longitude), math.cos(latitude) * math.sin(longitude), math.sin(latitude))
    return tuple(sun[i] + distance_pc * (along[0] * ex[i] + along[1] * ey[i] + along[2] * ez[i]) for i in range(3))


def _seed(galaxy_seed, name):
    text = galaxy_seed.hex() if isinstance(galaxy_seed, (bytes, bytearray)) else str(galaxy_seed)
    return draw.Stream(f"satellite:{text}:{name}").getrandbits(63)


def place_satellites(galaxy_seed, solar_radius_pc, solar_angle_rad, epoch_unix=0.0, scale=1.0):
    """
    The Milky Way's satellites as components of a galaxy.

    Args:
        galaxy_seed (bytes or int): The galaxy's seed, from which each satellite's own seed (the one its stars are drawn
            from) comes.
        solar_radius_pc (float): The Sun's distance from the galactic centre in this galaxy's frame.
        solar_angle_rad (float): The Sun's azimuth in that frame.
        epoch_unix (float): When the positions hold.
        scale (float): How much larger than the Milky Way the galaxy is, by its disk scale length: the satellites'
            distances from the Sun are multiplied by it, as the clusters' are, so a toy galaxy keeps its satellites
            inside itself. Their masses and sizes do not scale.

    Returns:
        list[Component]: One `satellite` component per row of `MILKY_WAY_SATELLITES`, in table order.
    """
    found = []
    for row in MILKY_WAY_SATELLITES:
        position = satellite_position_pc(
            row.longitude_deg, row.latitude_deg, row.distance_kpc * 1000.0 * scale, solar_radius_pc, solar_angle_rad)
        found.append(components.Component(
            kind="satellite", name=row.name, source="catalogue", center_pc=position, velocity_kms=(0.0, 0.0, 0.0),
            epoch_unix=epoch_unix, mass_solar=row.stellar_mass_sol,
            half_mass_radius_pc=row.half_light_radius_pc, profile=row.profile, age_gy=row.age_gy,
            metallicity_fe_h=row.fe_h, seed=_seed(galaxy_seed, row.name),
            tidal_radius_pc=component_profiles.cut_radius_pc(row.profile, row.half_light_radius_pc),
            catalogue_id=row.name))
    return found
