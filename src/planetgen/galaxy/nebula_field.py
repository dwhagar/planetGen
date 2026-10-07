# planetgen/galaxy/nebula_field.py

"""
The galaxy's molecular cloud field (GEN.47)
===========================================

Dark-family nebulae (giant molecular clouds, dark clouds, Bok globules and
star-forming cores, classes M-Q) are tens of light-years across, far
bigger than a 13 ly sector. Rolled once per sector and placed inside it,
as they used to be, a sector got about 3e-4 of them, so a generated
region almost never held one even when it sat inside a cloud.

Here the clouds belong to the galaxy instead. The galaxy is cut into
`NEBULA_FIELD_CELL_PC` cubes; each cube's clouds come from its own seed
(`galaxySeed.seeded`, kind `"nebula-cell"`, address `i/j/k`), so every
worker and every later run agrees where each cloud is, whatever order the
sectors run in. A cube's cloud count is a Poisson draw: the
`"molecular-cloud"` density (`PHENOMENON_DENSITY_PC3`) times the cube's
volume times the gas factor at its center (`gas_factor`, more in the arms
and near the plane). A sector takes every cloud whose sphere reaches it
(`clouds_reaching`), and `_db.insert_sector` stores each one once: the
first sector saved that it reaches is its home, and later sectors find it
already stored.
"""

import functools
import math
import random

from planetgen.galaxy import density as galaxyDensity, seed as galaxySeed
from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula, NEBULA_CLASS_LETTERS
from planetgen.galaxy.sector import _sample_poisson_count
from stellarObjects.utils import ly_to_pc

FIELD_FAMILY = "dark"
"""str: The nebula family the field places (classes M-Q). Emission and
reflection nebulae grow around a sector's own hot stars and planetary
nebulae around a new white dwarf (`generate.add_star_hosted_nebulae`,
`generate._add_planetary_nebula`); diffuse gas (A, B) isn't generated."""

FALLBACK_SEED = bytes(16)
"""bytes: The field's seed for a galaxy planned before it had one (schema
v51): still fixed, so its sectors agree with each other."""

FIELD_CLASSES = tuple(letter for letter in NEBULA_CLASS_LETTERS
                      if tuning.NEBULA_CLASSES[letter]["family"] == FIELD_FAMILY)
"""tuple: The field's class letters."""

MAX_CLOUD_RADIUS_PC = max(ly_to_pc(tuning.NEBULA_CLASSES[letter]["radius_range_ly"][1])
                          for letter in FIELD_CLASSES)
"""float: The largest cloud the field can place, parsecs: how far around
a sector its cells are read."""


REFERENCE_RADIUS_SCALE_LENGTHS = 2.82
"""float: The ring the `"molecular-cloud"` density is quoted on, in disk
scale lengths: the solar circle's ratio, `build_galaxy_shape`'s default
calibration radius."""


def _gas_tracer(position_pc, shape):
    """
    How much molecular gas sits at `position_pc`, in arbitrary units: the
    young stars' density (they form from it, so they share its thin layer
    and its radial profile), with the young stars' own arm contrast swapped
    for the gas's (`NEBULA_FIELD_ARM_AMPLITUDE`): molecular gas crowds the
    arms less than the stars born there do.
    """
    young = galaxyDensity.population_densities(position_pc, shape)["young"]
    arm_cos = galaxyDensity._arm_cosine(position_pc[0], position_pc[1], shape)
    star_contrast = 1 + tuning.STELLAR_POPULATION_ARM_AMPLITUDE["young"] * arm_cos
    if star_contrast <= 0.0:
        return 0.0
    return max(young, 0.0) / star_contrast * (1 + tuning.NEBULA_FIELD_ARM_AMPLITUDE * arm_cos)


@functools.lru_cache(maxsize=8)
def _reference_gas(shape):
    """`_gas_tracer` averaged around the reference ring on the plane
    (`REFERENCE_RADIUS_SCALE_LENGTHS`), arms and gaps alike."""
    radius = REFERENCE_RADIUS_SCALE_LENGTHS * shape.disk_scale_length_pc
    samples = [_gas_tracer(
        (radius * math.cos(2 * math.pi * step / 360), radius * math.sin(2 * math.pi * step / 360), 0.0),
        shape) for step in range(360)]
    return sum(samples) / len(samples)


def gas_factor(position_pc, shape):
    """
    Molecular gas at `position_pc` relative to the average around the
    solar circle (where the `"molecular-cloud"` density is quoted;
    `_gas_tracer`), raised to `GMC_GAS_DENSITY_EXPONENT`
    (Schmidt-Kennicutt): about 2 on an arm's crest and 0.3 between the
    arms at the solar circle, falling off within about 100 pc of the plane,
    and capped at `NEBULA_FIELD_MAX_GAS_FACTOR` toward the center.
    """
    reference = _reference_gas(shape)
    if reference <= 0.0:
        return 0.0
    return min(tuning.NEBULA_FIELD_MAX_GAS_FACTOR,
               (_gas_tracer(position_pc, shape) / reference) ** tuning.GMC_GAS_DENSITY_EXPONENT)


def cell_index(position_pc):
    """The `(i, j, k)` field cell holding `position_pc`."""
    edge = tuning.NEBULA_FIELD_CELL_PC
    return tuple(math.floor(coordinate / edge) for coordinate in position_pc)


def cell_clouds(galaxy_seed, shape, index):
    """
    The clouds of field cell `index`, from its own seed: a list of
    `(Nebula, center_pc)`, the same every time for the same galaxy.
    """
    edge = tuning.NEBULA_FIELD_CELL_PC
    low = tuple(i * edge for i in index)
    middle = tuple(corner + edge / 2 for corner in low)
    mean = (tuning.PHENOMENON_DENSITY_PC3["molecular-cloud"]
            * tuning.PHENOMENON_RATE_SCALE.get("molecular-cloud", 1.0)
            * edge ** 3 * gas_factor(middle, shape))
    clouds = []
    with galaxySeed.seeded(galaxy_seed or FALLBACK_SEED, "nebula-cell", index):
        count = _sample_poisson_count(mean)
        for _ in range(count):
            center = tuple(corner + random.uniform(0.0, edge) for corner in low)
            nebula = Nebula(SystemConfig(), nebula_type=FIELD_FAMILY)
            clouds.append((nebula, center))
    return clouds


def clouds_reaching(galaxy_seed, shape, center_pc, reach_pc):
    """
    Every field cloud whose sphere comes within `reach_pc` of `center_pc`
    (a sector's center and half diagonal, as `_db.sectors_reached_by`
    tests it), as `(Nebula, center_pc)` pairs, nearest cell first.
    """
    edge = tuning.NEBULA_FIELD_CELL_PC
    span = MAX_CLOUD_RADIUS_PC + reach_pc
    low = cell_index(tuple(c - span for c in center_pc))
    high = cell_index(tuple(c + span for c in center_pc))
    found = []
    for i in range(low[0], high[0] + 1):
        for j in range(low[1], high[1] + 1):
            for k in range(low[2], high[2] + 1):
                index = (i, j, k)
                nearest = tuple(min(max(c, n * edge), (n + 1) * edge) for c, n in zip(center_pc, index))
                if math.dist(nearest, center_pc) > span:
                    continue
                for nebula, cloud_center in cell_clouds(galaxy_seed, shape, index):
                    if math.dist(cloud_center, center_pc) <= ly_to_pc(nebula.radius_ly) + reach_pc:
                        found.append((nebula, cloud_center))
    return found
