# planetgen/generation/sector_context.py

"""
The stellar context a system made inside a stored sector is drawn under (extracted from the API's system routes for GEN.229).

A system made for a galaxy-placed sector after the sector was filled (an admin's regeneration of one system, the re-roll of an unstable wide binary) is given what a sector
fill gives it: the sector's stellar population and, once a bright-star scatter ran, its dim-star cap (its own backfill level and mass, GEN.44, GEN.187), as
`brightStars.FillContext` applies them, and the distance of the sector from the galactic center for the star's metallicity.
"""

from planetgen.db import store
from planetgen.galaxy.profile import profile_at
from planetgen.generation import bright_stars as brightStars
from planetgen.physics.units import pc_to_ly


def apply_sector_context(conn, sector_id, system_config):
    """
    Gives `system_config` the stellar population of a galaxy-placed sector, as a sector fill would.

    Args:
        conn (store.Connection): An open connection to the galaxy.
        sector_id (int): The sector the system is (to be) filed in.
        system_config (SystemConfig): The recipe, changed in place.

    Returns:
        float or None: The sector's distance from the galactic center, light-years; `None` for a sector outside the galaxy.

    Raises:
        ValueError: There is no such sector.
    """
    placement = store.get_sector_galaxy_position(conn, sector_id)
    if placement is None:
        return None
    skeleton = store.get_galaxy_shape(conn)
    if skeleton is not None:
        address = (placement["shell_index"], placement["band_index"], placement["wedge_index"])
        level = store.bright_star_fill_level(conn, *address)
        mass_limit = store.bright_star_fill_mass_limit(conn, *address)
        center = (placement["center_x_pc"], placement["center_y_pc"], placement["center_z_pc"])
        profile = profile_at(center, skeleton.shape, skeleton.edge_pc, skeleton.expected_system_count_at_density_1,
                             address=address)
        brightStars.FillContext(profile, min_luminosity_sol=level, star_mass_limit_sol=mass_limit).apply(system_config)
    return pc_to_ly(placement["galactic_radius_pc"])
