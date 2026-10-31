# tests/reference_satellites.py

"""
The Milky Way's satellite galaxies as reference data (GEN.369): the table the generator used to place (GEN.237), kept
for the plausibility checks of the drawn population (`test_satellite_law.py`) and as a fixture for the tests of what
the sector gate, the profile and the store do with `catalogue` satellites. The generator reads none of it: its
satellites are drawn from the galaxy's seed and halo (`generation/satellites.py`).

Checked against McConnachie (2012, AJ 144, 4; arXiv:1204.1562) on 2026-10-31 (GEN.366): the position and distance of each row are Table 2
(l, b, D), the half-light radius is the semi-major axis radius of Table 3 in parsecs, the stellar mass is Table 4 (a stellar mass-to-light ratio of
1) and [Fe/H] is Table 5. The paper gives no half-light radius for the two Magellanic Clouds (they are not fitted by a single profile): 2700 and 1100 pc are
not from it. The stellar ages are not in the paper either (the ages are the table's own estimates, unchecked; Leo I's was 5 Gy, which asks for a rising history and so for gas, and Leo I has none: Grcevich and Putman 2009 (arXiv:0901.4975) find the Milky Way's dwarfs within about 270 kpc under 1e5 solar masses of HI, and a limit of under 3e4 for Leo I, TEST.508; it is 9 Gy here so the row forms no stars, and the age is a choice, not a measurement). Canes Venatici I (TEST.507) is Tables 2 to 5 as well, its age is a choice (an old dwarf: 11 Gy, as Draco and Sculptor). The table holds the 12 satellites of the paper above 1e5 solar masses; Canis Major (4.9e7 solar masses) is left out: the paper lists its type as unknown (Table 1, "MW disk substructure?") and says there is "considerable debate as to whether it is a dwarf galaxy or an unrelated disk substructure", so it is not a settled satellite for a table that checks a drawn population.
`test_reference_satellites.py` pins
the checked values. The half-light radius stands in for the half-mass radius. The velocities are zero.
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
    Satellite("Large Magellanic Cloud", 280.5, -32.9, 51.0, 1.5e9, 2700.0, "exponential", -0.5, 5.0),
    Satellite("Small Magellanic Cloud", 302.8, -44.3, 64.0, 4.6e8, 1100.0, "exponential", -1.0, 6.0),
    Satellite("Sagittarius dwarf", 5.6, -14.2, 26.0, 2.1e7, 2587.0, "exponential", -0.4, 8.0),
    Satellite("Fornax dwarf", 237.1, -65.7, 147.0, 2.0e7, 710.0, "plummer", -0.99, 7.0),
    Satellite("Sculptor dwarf", 287.5, -83.2, 86.0, 2.3e6, 283.0, "plummer", -1.68, 11.0),
    Satellite("Carina dwarf", 260.1, -22.2, 105.0, 3.8e5, 250.0, "plummer", -1.72, 7.0),
    Satellite("Draco dwarf", 86.4, 34.7, 76.0, 2.9e5, 221.0, "plummer", -1.93, 11.0),
    Satellite("Ursa Minor dwarf", 105.0, 44.8, 76.0, 2.9e5, 181.0, "plummer", -2.13, 12.0),
    Satellite("Sextans dwarf", 243.5, 42.3, 86.0, 4.4e5, 695.0, "plummer", -1.93, 11.0),
    Satellite("Leo I", 226.0, 49.1, 254.0, 5.5e6, 251.0, "plummer", -1.43, 9.0),
    Satellite("Leo II", 220.2, 67.2, 233.0, 7.4e5, 176.0, "plummer", -1.62, 9.0),
    Satellite("Canes Venatici I", 74.3, 79.8, 218.0, 2.3e5, 564.0, "plummer", -1.98, 11.0),
)
"""tuple: The Milky Way's satellites (checked against McConnachie 2012: see the module docstring)."""


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
