# planetgen/galaxy/superbubble_field.py

"""
The galaxy's superbubble field (GEN.254)
========================================

A superbubble is the cavity the supernovae of an OB association blow in the interstellar medium: 100 to 1,000 pc
across, filled with ionised gas of about a million kelvin and a thousandth of the disk's density. The Sun sits in one,
the Local Bubble, 100 to 200 pc across. Molecular clouds do not survive inside one, so the cloud field
(`nebula_field`) keeps none whose centre lies in a cavity.

The field is seeded and stores nothing, like the cloud field. The plane is cut into `SUPERBUBBLE_CELL_PC` columns;
each column's bubbles come from its own seed (`galaxySeed.seeded`, kind `"superbubble-cell"`, address `i/j`), so every
worker and every later run agrees where each one is. A column's bubble count is a Poisson draw whose mean puts
`SUPERBUBBLE_FILL_FRACTION` of the plane at the solar circle inside a bubble and follows the molecular gas
(`nebula_field.gas_factor`) elsewhere: more in the arms and the inner disk. Each bubble has a centre in its column,
at a height about the young stars' layer, a log-uniform radius, and the gas inside it.
"""

import functools
import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.galaxy import nebula_field, seed as galaxySeed
from planetgen.galaxy.sector import _sample_poisson_count
from planetgen.util import draw
from planetgen.util import detmath

FALLBACK_SEED = nebula_field.FALLBACK_SEED
"""bytes: The field's seed for a galaxy planned before it had one."""


@dataclass(frozen=True)
class Superbubble:
    """
    One cavity.

    Attributes:
        index (tuple): The `(i, j)` column that drew it.
        n (int): Its place in the column's draw.
        center_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        radius_pc (float): The cavity's radius.
        ambient_density_cm3 (float): The gas inside, H atoms per cm3.
        temperature_k (float): Its temperature.
    """
    index: tuple
    n: int
    center_pc: tuple
    radius_pc: float
    ambient_density_cm3: float
    temperature_k: float

    def depth(self, point_pc):
        """How far `point_pc` is from the centre as a share of the radius: under 1 inside."""
        return detmath.dist(self.center_pc, point_pc) / self.radius_pc


def _mean_cross_section_pc2():
    """The mean area `pi r^2` of a bubble cut through its centre, for a log-uniform radius."""
    low, high = tuning.SUPERBUBBLE_RADIUS_RANGE_PC
    return math.pi * (high * high - low * low) / (2.0 * math.log(high / low))


def centres_per_pc2():
    """
    The surface density of bubble centres at the solar circle that puts `SUPERBUBBLE_FILL_FRACTION` of the plane
    inside a bubble: coverage `1 - exp(-n * area)` of overlapping discs.
    """
    return -math.log(1.0 - tuning.SUPERBUBBLE_FILL_FRACTION) / _mean_cross_section_pc2()


def cell_index(position_pc):
    """The `(i, j)` column holding `position_pc`."""
    edge = tuning.SUPERBUBBLE_CELL_PC
    return (math.floor(position_pc[0] / edge), math.floor(position_pc[1] / edge))


@functools.lru_cache(maxsize=8192)
def cell_bubbles(galaxy_seed, shape, index):
    """The bubbles of column `index`, a tuple of `Superbubble`: the same every time for the same galaxy."""
    edge = tuning.SUPERBUBBLE_CELL_PC
    low = (index[0] * edge, index[1] * edge)
    middle = (low[0] + edge / 2, low[1] + edge / 2, 0.0)
    mean = centres_per_pc2() * edge * edge * nebula_field.gas_factor(middle, shape)
    bubbles = []
    with galaxySeed.seeded(galaxy_seed or FALLBACK_SEED, "superbubble-cell", index, shape=shape):
        for n in range(_sample_poisson_count(mean)):
            x, y = low[0] + draw.uniform(0.0, edge), low[1] + draw.uniform(0.0, edge)
            z = draw.gauss(0.0, tuning.SUPERBUBBLE_CENTRE_HEIGHT_SIGMA_PC)
            radius_low, radius_high = tuning.SUPERBUBBLE_RADIUS_RANGE_PC
            radius = math.exp(draw.uniform(math.log(radius_low), math.log(radius_high)))
            density_low, density_high = tuning.SUPERBUBBLE_AMBIENT_DENSITY_RANGE_CM3
            density = math.exp(draw.uniform(math.log(density_low), math.log(density_high)))
            bubbles.append(Superbubble(tuple(index), n, (x, y, z), radius, density, tuning.SUPERBUBBLE_TEMPERATURE_K))
    return tuple(bubbles)


def bubbles_reaching(galaxy_seed, shape, point_pc, reach_pc=0.0):
    """Every bubble whose sphere comes within `reach_pc` of `point_pc`, nearest column first."""
    edge = tuning.SUPERBUBBLE_CELL_PC
    span = tuning.SUPERBUBBLE_RADIUS_RANGE_PC[1] + reach_pc
    low = cell_index((point_pc[0] - span, point_pc[1] - span))
    high = cell_index((point_pc[0] + span, point_pc[1] + span))
    found = []
    for i in range(low[0], high[0] + 1):
        for j in range(low[1], high[1] + 1):
            for bubble in cell_bubbles(galaxy_seed, shape, (i, j)):
                if detmath.dist(bubble.center_pc, point_pc) <= bubble.radius_pc + reach_pc:
                    found.append(bubble)
    return found


def containing(galaxy_seed, shape, point_pc):
    """
    The bubble `point_pc` is inside, or `None`. Where several overlap, the one it is deepest in (the smallest share of
    its radius from the centre), so the answer does not depend on the order the columns are read.
    """
    inside = [bubble for bubble in bubbles_reaching(galaxy_seed, shape, point_pc) if bubble.depth(point_pc) < 1.0]
    return min(inside, key=lambda bubble: (bubble.depth(point_pc), bubble.index, bubble.n)) if inside else None


def all_bubbles(galaxy_seed, shape, radius_pc, half_height_pc=None):
    """
    Every bubble whose centre lies inside the cylinder of `radius_pc` around the galactic axis (and `half_height_pc`
    above and below the plane when given), by column then draw order: what the Galaxy Map draws.
    """
    edge = tuning.SUPERBUBBLE_CELL_PC
    reach = math.floor(radius_pc / edge) + 1
    found = []
    for i in range(-reach, reach + 1):
        for j in range(-reach, reach + 1):
            for bubble in cell_bubbles(galaxy_seed, shape, (i, j)):
                x, y, z = bubble.center_pc
                if detmath.hypot(x, y) <= radius_pc and (half_height_pc is None or abs(z) <= half_height_pc):
                    found.append(bubble)
    return found
