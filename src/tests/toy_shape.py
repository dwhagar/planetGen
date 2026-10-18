# tests/toy_shape.py

"""
Toy galaxies for the tests (GEN.404).

A galaxy holds the systems its stellar mass makes (`galaxy.density.build_galaxy_shape`), so a toy one with a disk of a
few hundred parsecs and the Milky Way's stellar mass would hold them in a hundred-thousandth of the volume: a density
no test means to read. A toy galaxy here is read in relative units instead: its density is 1 at its reference radius
between two arms, the convention the tests were written to before the density followed the mass, and the sector
counts follow from the test's own `expected_system_count_at_density_1`.
"""

import math

from planetgen import tuning
from planetgen.galaxy import density


_BUILD = density.build_galaxy_shape


def at_unit_density(shape):
    """`shape` with its normalization set so the relative density is 1 at its reference radius, between two arms."""
    radius = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * shape.disk_scale_length_pc
    angle = density._interarm_angle(radius, shape)
    position = (radius * math.cos(angle), radius * math.sin(angle), 0.0)
    return shape._replace(k_norm=shape.k_norm / density.relative_density(position, shape))


def toy_shape(*args, **kwargs):
    """`build_galaxy_shape(*args, **kwargs)` read in relative units: 1.0 at its reference radius between two arms."""
    return at_unit_density(_BUILD(*args, **kwargs))


def plan_at_unit_density(plan_shape):
    """`shape_draw.plan_shape` wrapped so the shape it plans is read in relative units (see `unit_density_plans`)."""
    def planned(args, galaxy_seed):
        shape, drawn = plan_shape(args, galaxy_seed)
        return at_unit_density(shape), drawn
    return planned
