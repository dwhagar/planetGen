# planetgen/physics/trojan_overlay.py

"""
Where a trojan region lies at a time (GEN.1150 slice 2, the libration overlay)
==============================================================================

A region is stored with its giant as parent and the libration of its bodies (`asteroid_belts.lagrange_point`, `libration_amplitude_deg`, `libration_phase_deg`), never with a position of its
own: the orbit updater moves the giant, and the region is the giant's mean longitude plus `+-60` degrees plus the slow libration `A sin(2 pi t / P_lib + phase)` (`trojan_swarm.longitude_offset_deg`,
`trojans.libration_period`). `t` is the time since Unix time zero, so the phase is a number of the region and no epoch column is needed. A reader (the system map, a sheet) asks for the offset
from the giant at the moment it wants, so a stored region is where it should be whenever it is read, with nothing for the updater to write.
"""

from planetgen.physics import constants, trojan_swarm, trojans

SECONDS_PER_YEAR = constants.SECONDS_PER_YEAR
"""float: Seconds in a Julian year, the unit of the libration period."""


def mass_parameter(secondary_mass, host_mass):
    """`mu = m2 / (m1 + m2)` of a giant of `secondary_mass` about a star of `host_mass` (one unit), or `None` when either is missing or not positive."""
    if not secondary_mass or not host_mass or secondary_mass <= 0.0 or host_mass <= 0.0:
        return None
    return trojans.mass_parameter(host_mass, secondary_mass)


def offset_deg(point, amplitude_deg, phase_deg, unix_s, orbital_period_years, mu):
    """
    The mean longitude of a region's center from its giant's, degrees, at Unix time `unix_s`: `+60` for `L4`, `-60` for `L5`, plus the libration. A giant with no known period or mass
    parameter, or a pair that holds no tadpole (`trojans.libration_period` is `None`), gives the point alone.
    """
    base = trojans.longitude_offset_deg(point)
    if not orbital_period_years or not mu or orbital_period_years <= 0.0:
        return base
    region = trojan_swarm.Region(point, 0.0, amplitude_deg or 0.0, phase_deg or 0.0, 0.0, 0.0)
    try:
        return trojan_swarm.longitude_offset_deg(region, unix_s / SECONDS_PER_YEAR, orbital_period_years, mu)
    except trojans.TrojanError:
        return base
