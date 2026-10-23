# tests/reference_satellites.py

"""
The Milky Way's satellite galaxies as reference data (GEN.369): the table the generator used to place (GEN.237), kept
for the plausibility checks of the drawn population (`test_satellite_law.py`) and as a fixture for the tests of what
the sector gate, the profile and the store do with `catalogue` satellites. The generator reads none of it: its
satellites are drawn from the galaxy's seed and halo (`generation/satellites.py`).

Every value in the table is recalled from McConnachie (2012, AJ 144, 4) and not checked against it, so each is good to a
few tens of percent in size and distance and a factor of two in mass. The half-light radius stands in for the half-mass
radius. The velocities are zero.
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


def sun_frame(frame_radius_pc, frame_angle_rad):
    """`(sun, ex, ey, ez)`: the Sun's position in the galaxy frame and the unit vectors of the galactic coordinates there.

    `ex` points from the Sun to the galactic centre (longitude 0, latitude 0), `ey` towards longitude 90 degrees (the
    direction the disk turns) and `ez` along the galaxy frame's +z. Those three form a left-handed set, so a proper
    galactic frame has its north pole towards -z: `satellite_position_pc` takes a latitude `b` as a height of
    `-sin(b)` in the galaxy frame.
    """
    cos_a, sin_a = math.cos(frame_angle_rad), math.sin(frame_angle_rad)
    sun = (frame_radius_pc * cos_a, frame_radius_pc * sin_a, 0.0)
    return sun, (-cos_a, -sin_a, 0.0), (-sin_a, cos_a, 0.0), (0.0, 0.0, 1.0)


def satellite_position_pc(longitude_deg, latitude_deg, distance_pc, frame_radius_pc, frame_angle_rad):
    """The position in the galaxy frame (pc, the galactic centre at the origin) of an object seen from the Sun at this
    galactic longitude and latitude and this distance."""
    sun, ex, ey, ez = sun_frame(frame_radius_pc, frame_angle_rad)
    longitude, latitude = math.radians(longitude_deg), -math.radians(latitude_deg)
    along = (math.cos(latitude) * math.cos(longitude), math.cos(latitude) * math.sin(longitude), math.sin(latitude))
    return tuple(sun[i] + distance_pc * (along[0] * ex[i] + along[1] * ey[i] + along[2] * ez[i]) for i in range(3))


def _seed(galaxy_seed, name):
    text = galaxy_seed.hex() if isinstance(galaxy_seed, (bytes, bytearray)) else str(galaxy_seed)
    return draw.Stream(f"satellite:{text}:{name}").getrandbits(63)


def place_reference_satellites(galaxy_seed, frame_radius_pc, frame_angle_rad, epoch_unix=0.0, scale=1.0):
    """
    The Milky Way's satellites as `catalogue` components of a galaxy (a fixture: no generator reads the table).

    Args:
        galaxy_seed (bytes or int): The galaxy's seed, from which each satellite's own seed (the one its stars are drawn
            from) comes.
        frame_radius_pc (float): The catalogue's observer's distance from the galactic centre in this galaxy's frame (the reference radius).
        frame_angle_rad (float): The Sun's azimuth in that frame.
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
            row.longitude_deg, row.latitude_deg, row.distance_kpc * 1000.0 * scale, frame_radius_pc, frame_angle_rad)
        found.append(components.Component(
            kind="satellite", name=row.name, source="catalogue", center_pc=position, velocity_kms=(0.0, 0.0, 0.0),
            epoch_unix=epoch_unix, mass_solar=row.stellar_mass_sol,
            half_mass_radius_pc=row.half_light_radius_pc, profile=row.profile, age_gy=row.age_gy,
            metallicity_fe_h=row.fe_h, seed=_seed(galaxy_seed, row.name),
            tidal_radius_pc=component_profiles.cut_radius_pc(row.profile, row.half_light_radius_pc),
            catalogue_id=row.name))
    return found
