# planetgen/galaxy/plan_figures.py

"""
The two density figures the plan prints about its galaxy (GEN.626), read from the built shape and not from a table:

- the **crest**, how many systems the thin disc expects in a sector on an arm crest at its centre, the number the
  outline walks on (`k_norm x arm_level x (1 + arm_amplitude) x E`, `E` the systems a sector holds at the real local
  stellar density), which is also the thin disc's peak density in units of the one-star density;
- the **planned share**, the share of the galaxy's own systems that sit in sectors expecting at least one system,
  integrated from the model's density over radius, height and azimuth.

Both follow the galaxy's own `k_norm`, so the Census dial (richness) and the galaxy's own system mass are in them. The same
integral gives the **reaching sectors** (GEN.627), the sectors expecting at least one system, to set beside the outline's
candidate sectors (a superset by design: the outline takes the best longitude and the run of bands up to a component's).
"""

import collections
import math

import numpy as np

from planetgen.galaxy import density
from planetgen.physics import constants as physical_constants
from planetgen.physics.units import ly_to_pc

SHARE_RADIAL_STEPS = 240
"""int: Steps of the share's integral in radius and in height, geometric from `SHARE_INNER_PC` out."""

SHARE_INNER_PC = 0.5
"""float: The first edge, parsecs, of the share's grid in radius and height (the first step runs from the axis or the
plane); the nuclear parts are a few parsecs across."""

SHARE_AZIMUTHS = 24
"""int: Azimuths the share's integral samples around each ring, so the arms and the bar are in it."""


def crest(shape, expected_systems_at_density_1):
    """The thin disc's crest: the systems per sector it expects on an arm crest at its centre, the density in units of
    the one-star density (`expected_systems_at_density_1` is `E`, a sector's expected systems at density 1)."""
    return (shape.k_norm * shape.arm_level * (1.0 + shape.arm_amplitude) * expected_systems_at_density_1)


Planned = collections.namedtuple("Planned", ["share", "volume_pc3"])
"""The planned share of the galaxy's systems, and the volume (cubic parsecs) of the disc and bulge's sectors expecting one."""


def planned_share(shape, expected_systems_at_density_1, radius_pc):
    """The planned share alone (`planned(...).share`)."""
    return planned(shape, expected_systems_at_density_1, radius_pc).share


def reaching_sectors(shape, expected_systems_at_density_1, radius_pc, edge_pc):
    """
    How many sectors expect at least one system: the volume where the model's density reaches it (`planned`) over a
    sector's volume, plus each component's sphere out to where its own density reaches one system in a cell
    (`ComponentField.qualifying_volume_pc3`). Approximate to a few percent; a component inside the disc is counted twice.
    """
    volume = planned(shape, expected_systems_at_density_1, radius_pc).volume_pc3
    if shape.components is not None:
        volume += shape.components.qualifying_volume_pc3()
    return float(volume / edge_pc ** 3)


def planned(shape, expected_systems_at_density_1, radius_pc):
    """
    The share of the galaxy's own systems (its stars, not its clusters or satellites) that sit in sectors expecting at
    least one system, `0` to `1`, and the volume of those sectors, as a `Planned`.

    The model's density (`density.total_density_array`) is integrated over cylindrical radius and height out to
    `radius_pc` (the outline's outer edge) on a geometric grid and `SHARE_AZIMUTHS` azimuths; the systems in a cell of
    the grid count when `E x density` there reaches 1, and the total is the galaxy's systems (its stellar mass over
    its mean system mass). Approximate to a few percent.

    Args:
        shape (galaxy.density.GalaxyShape): The galaxy's shape, `k_norm` set.
        expected_systems_at_density_1 (float): `E`.
        radius_pc (float): How far out to integrate; nothing farther holds a sector that reaches one system.

    Returns:
        Planned: The share (`0.0` for a galaxy with no radius) and the volume.
    """
    if radius_pc <= SHARE_INNER_PC:
        return Planned(0.0, 0.0)
    edges = np.concatenate(([0.0], np.geomspace(SHARE_INNER_PC, radius_pc, SHARE_RADIAL_STEPS)))
    middles = np.sqrt(edges[:-1] * edges[1:])
    middles[0] = 0.5 * edges[1]
    ring_area = math.pi * np.diff(edges ** 2)
    slab = 2.0 * np.diff(edges)                 # both sides of the plane
    azimuths = (np.arange(SHARE_AZIMUTHS) + 0.5) * (2.0 * math.pi / SHARE_AZIMUTHS)
    r_grid = middles[:, None, None]
    z_grid = middles[None, :, None]
    x = r_grid * np.cos(azimuths)[None, None, :]
    y = r_grid * np.sin(azimuths)[None, None, :]
    rho = density.total_density_array(x, y, z_grid, shape)
    cell_volume = ring_area[:, None, None] * slab[None, :, None] / SHARE_AZIMUTHS
    reaching = rho * expected_systems_at_density_1 >= 1.0
    planned_systems = (rho * cell_volume)[reaching].sum()
    total = shape.structure.stellar_mass_solar / density.mean_system_mass_solar(shape)
    local_per_pc3 = physical_constants.LOCAL_STELLAR_DENSITY_LY3 / ly_to_pc(1.0) ** 3
    return Planned(float(min(1.0, planned_systems * local_per_pc3 / total)),
                   float((np.broadcast_to(cell_volume, rho.shape))[reaching].sum()))
