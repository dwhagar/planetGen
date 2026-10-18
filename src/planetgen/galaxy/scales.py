# planetgen/galaxy/scales.py

"""
The galaxy's own scales for the laws that used to be quoted at the Sun (GEN.404)
==================================================================================

A galaxy stores no home or Sun position. What a law once read as "the Sun's radius" is a radius in the galaxy's own
unit, its thin disk's scale length: the reference radius (`tuning.REFERENCE_RADIUS_SCALE_LENGTHS`), the radius the
hypervelocity stars' density is quoted at, the radius the clusters' sizes start to grow from, and the radius a system
with no place is taken at. The scale length is the running unit's galaxy's (`galaxy.seed.current_disk_scale_length_pc`),
else the one its seed draws (`shape_draw`), so work with no galaxy uses the fallback seed's galaxy.
"""

import functools

from planetgen import tuning
from planetgen.galaxy import mass_chain, seed as galaxySeed
from planetgen.physics.units import pc_to_ly


@functools.lru_cache(maxsize=64)
def _drawn_scale_length_pc(galaxy_seed):
    from planetgen.galaxy import shape_draw   # shape_draw imports structure_draw, which imports this package's draws
    return shape_draw.draw_galaxy(galaxy_seed).disk_scale_length_pc


def disk_scale_length_pc():
    """The thin disk scale length, parsecs, of the galaxy the running unit belongs to."""
    bound = galaxySeed.current_disk_scale_length_pc()
    if bound is not None:
        return bound
    return _drawn_scale_length_pc(bytes(galaxySeed.current_galaxy_seed() or mass_chain.FALLBACK_SEED))


def reference_radius_pc():
    """The reference radius, parsecs: `tuning.REFERENCE_RADIUS_SCALE_LENGTHS` scale lengths."""
    return tuning.REFERENCE_RADIUS_SCALE_LENGTHS * disk_scale_length_pc()


def hypervelocity_reference_radius_pc():
    """The radius the hypervelocity stars' density is quoted at, parsecs (`tuning.HVS_REFERENCE_SCALE_LENGTHS`)."""
    return tuning.HVS_REFERENCE_SCALE_LENGTHS * disk_scale_length_pc()


def default_galactocentric_pc():
    """Where a system with no place in the galaxy is taken to be, parsecs from the centre: the radius that holds
    half of an exponential disc's stars, `tuning.DISK_HALF_MASS_SCALE_LENGTHS` scale lengths (1.678)."""
    return tuning.DISK_HALF_MASS_SCALE_LENGTHS * disk_scale_length_pc()


def default_galactocentric_ly():
    """`default_galactocentric_pc` in light-years."""
    return pc_to_ly(default_galactocentric_pc())
