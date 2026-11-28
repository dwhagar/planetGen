# planetgen/galaxy/bubble_field.py

"""
The superbubble handle on the galaxy shape (GEN.392; docs/design/galaxy/superbubble-density.md)
================================================================================================

`galaxy.superbubble_field` (GEN.254) puts a seeded cavity set in every galaxy and stores nothing. This module is what a density
read holds to ask that set what a bubble does to the young stars: `GalaxyShape.bubbles` carries a `BubbleField`, the galaxy's seed
and nothing else, as `GalaxyShape.components` carries the clusters. It answers three things:

* **A box cut.** `BubbleField.bubbles_in_box` names the bubbles whose smooth edge reaches an axis-aligned box (a sector, a band's
  bin), so a build reads the few bubbles that touch it and not the field.
* **The factor.** `young_factor` is what a bubble does to the young population at a point, 1 far outside it (GEN.391, Boss
  2026-11-26): the stars that formed after the bubble opened are missing inside, so a young star is `interior_factor` as likely
  there as outside (0.6 for a bubble of the full 40 Myr), and the gas the bubble swept out is a cooled shell about 0.06 of the
  radius thick, in which young stars form faster by the shell's density to the Schmidt law's power. Overlaps multiply. The
  intermediate and old populations, the bulge, the halo and the components are not changed (Boss, 2026-10-16).
* **A certified maximum.** `max_young_factor` is the most the factor can be anywhere in a box: the product of the shell's peak for
  every bubble that reaches the box, so a scatter's bound is a bound.

The factor is smooth in the depth `d` (the distance from the center over the radius): `interior_factor` for `d` under the shell,
rising to the shell's peak across the inner edge of the shell, flat over its thickness, and falling to 1 across its outer edge. At
every depth it is a mix of 1, the interior factor and the peak, so it never goes beyond the peak: that is the certificate.

A bubble's age comes from its radius by Weaver's [R5374] expansion law at one driving power for all (Boss, 2026-11-26): the
law's radius `0.76 * (power * age ** 3 / density) ** (1 / 5)` solved for the age, no more than the 40 Myr the cluster's winds and
supernovae drive it for [R5373]. No real bubble is the anchor. The bubbles' own placement reads the
density with no bubbles in it (`bubble_free`, `superbubble_field.bubble_free`), so the factor never feeds back into where the
bubbles are.
"""

import math

from planetgen import tuning
from planetgen.galaxy import density as galaxyDensity, superbubble_field
from planetgen.galaxy.nodecut import smoothstep
from planetgen.physics import constants as physical_constants
from planetgen.util import detmath


def age_myr(bubble):
    """
    A bubble's age, Myr, from its radius by Weaver's [R5374] law at the one driving power and ambient density of every bubble
    (`radius = coefficient * (power * age ** 3 / density) ** (1 / 5)`, solved for the age), no more than the driving power's lifetime.
    A small bubble is young: about 2 Myr at 75 pc, 40 Myr from about 400 pc.
    """
    density_g_cm3 = (tuning.SUPERBUBBLE_DRIVING_DENSITY_CM3 * tuning.SUPERBUBBLE_MEAN_MASS_PER_PARTICLE
                     * physical_constants.HYDROGEN_ATOM_MASS_KG * 1000.0)
    radius_cm = bubble.radius_pc * physical_constants.PARSEC_M * 100.0
    age_s = (radius_cm / tuning.SUPERBUBBLE_WEAVER_COEFFICIENT) ** (5.0 / 3.0) * (
        density_g_cm3 / tuning.SUPERBUBBLE_DRIVING_POWER_ERG_S) ** (1.0 / 3.0)
    return min(age_s / (physical_constants.SECONDS_PER_YEAR * 1.0e6), tuning.SUPERBUBBLE_EXPANSION_LIFETIME_MYR)


def new_star_floor(bubble):
    """
    How many stars form inside a cavity over how many form in the disk's gas: the Schmidt law (`tuning.SUPERBUBBLE_SCHMIDT_EXPONENT`)
    of the cavity's gas over the disk's (`tuning.DISK_ISM_DENSITY_CM3`), 0.0003 to 0.002 for the field's 0.001 to 0.01 atoms per cm3.
    """
    return (bubble.ambient_density_cm3 / tuning.DISK_ISM_DENSITY_CM3) ** tuning.SUPERBUBBLE_SCHMIDT_EXPONENT


def interior_factor(bubble):
    """
    How likely a young star is inside the cavity over outside it: the stars younger than the bubble formed in the gas it swept
    out, so a share `age / young window` of them is down to `new_star_floor` and the rest is as outside. 0.6 for a bubble of
    40 Myr (a young window of 100 Myr); nearly 1 for a small young one.
    """
    missing = bubble_missing_share(bubble)
    return 1.0 - missing * (1.0 - new_star_floor(bubble))


def bubble_missing_share(bubble):
    """The share of the young window's ages the bubble is older than: the share of its young stars that never formed inside."""
    return min(1.0, age_myr(bubble) / tuning.SUPERBUBBLE_YOUNG_AGE_WINDOW_MYR)


def shell_factor():
    """
    How much faster young stars form in the shell than in the disk, the same for every bubble: the gas of the ball piled into a
    shell of thickness `t` of the radius is `1 / (3 t)` times denser, and star formation goes as the gas density to the Schmidt
    law's power (`tuning.SUPERBUBBLE_SCHMIDT_EXPONENT`, the disk's own).
    """
    return (1.0 / (3.0 * tuning.SUPERBUBBLE_SHELL_THICKNESS)) ** tuning.SUPERBUBBLE_SCHMIDT_EXPONENT


def reach_depth():
    """How far from a bubble's center its factor differs from 1, in radii: the shell's outer edge, `1 + edge` (the edge's half-width is a share of the shell's thickness)."""
    return 1.0 + tuning.SUPERBUBBLE_SHELL_THICKNESS * tuning.SUPERBUBBLE_EDGE_SHARE


def young_factor(bubble, depth):
    """
    What `bubble` does to the young population at `depth` (distance from its center over its radius): `interior_factor` inside
    the shell, `shell_factor` over the shell, 1 beyond it, smooth across both edges (a smoothstep over twice the edge half-width).
    """
    half = tuning.SUPERBUBBLE_SHELL_THICKNESS * tuning.SUPERBUBBLE_EDGE_SHARE
    inner = 1.0 - tuning.SUPERBUBBLE_SHELL_THICKNESS
    rise = smoothstep((depth - (inner - half)) / (2.0 * half))
    fall = smoothstep((depth - (1.0 - half)) / (2.0 * half))
    inside = 1.0 - (1.0 - interior_factor(bubble)) * (1.0 - rise)
    on_shell = rise * (1.0 - fall)
    return (1.0 - on_shell) * inside + on_shell * shell_factor()


def bubble_free(shape):
    """`shape` without its bubble handle: the shape the bubbles' own placement reads, so no factor feeds back into where the bubbles are."""
    return superbubble_field.bubble_free(shape)


class BubbleField:
    """
    The galaxy's superbubbles as a density modulation (GEN.392): the galaxy's seed, carried by `GalaxyShape.bubbles`. It holds no
    bubble: they come from the seed (`superbubble_field.cell_bubbles`), the same on every worker and every run. It compares and
    hashes by its seed, so two shapes of one galaxy are one cache key.

    Args:
        galaxy_seed (bytes or None): The galaxy's seed; `None` takes the field's fallback seed, as the cloud field does.
    """

    __slots__ = ("galaxy_seed",)

    def __init__(self, galaxy_seed):
        self.galaxy_seed = None if galaxy_seed is None else bytes(galaxy_seed)

    def __eq__(self, other):
        return isinstance(other, BubbleField) and self.galaxy_seed == other.galaxy_seed

    def __hash__(self):
        return hash((BubbleField, self.galaxy_seed))

    def __repr__(self):
        return f"BubbleField({'None' if self.galaxy_seed is None else self.galaxy_seed.hex()})"

    def bubbles_in_box(self, shape, low_pc, high_pc):
        """
        Every bubble whose smooth edge (`reach_depth` radii) comes within the axis-aligned box `low_pc`..`high_pc`, by column and
        draw order: the same list for the same galaxy and box on every run.
        """
        bare = superbubble_field.bubble_free(shape)
        span = tuning.SUPERBUBBLE_RADIUS_RANGE_PC[1] * reach_depth()
        first = superbubble_field.cell_index((low_pc[0] - span, low_pc[1] - span))
        last = superbubble_field.cell_index((high_pc[0] + span, high_pc[1] + span))
        found = []
        for i in range(first[0], last[0] + 1):
            for j in range(first[1], last[1] + 1):
                for bubble in superbubble_field.cell_bubbles(self.galaxy_seed, bare, (i, j)):
                    if _box_distance(bubble.center_pc, low_pc, high_pc) <= bubble.radius_pc * reach_depth():
                        found.append(bubble)
        return found

    def young_factor(self, shape, point_pc):
        """The product of every bubble's `young_factor` at `point_pc`: 1 where none reaches it, below 1 inside a cavity, above 1 on a shell."""
        factor = 1.0
        for bubble in self.bubbles_in_box(shape, point_pc, point_pc):
            factor *= young_factor(bubble, bubble.depth(point_pc))
        return factor

    def max_young_factor(self, shape, low_pc, high_pc):
        """
        A bound the young factor never exceeds anywhere in the box `low_pc`..`high_pc` (certified): the shell's peak for every
        bubble that reaches it, multiplied (the factor of one bubble never goes beyond its peak, and overlaps multiply); 1 where no
        bubble reaches the box. A scatter thins by the factor over this.
        """
        return shell_factor() ** len(self.bubbles_in_box(shape, low_pc, high_pc))

    def density_without_bubbles(self, shape, position_pc):
        """
        The stellar density at `position_pc` by population with no bubble in it (`galaxy.density.population_densities` on
        `bubble_free(shape)`): the read the bubbles' placement and the factor's own bound are made against.
        """
        return galaxyDensity.population_densities(position_pc, bubble_free(shape))


def _box_distance(point_pc, low_pc, high_pc):
    """The distance from `point_pc` to the nearest point of the axis-aligned box `low_pc`..`high_pc` (0 inside it)."""
    return detmath.hypot(*(point_pc[a] - min(max(point_pc[a], low_pc[a]), high_pc[a]) for a in range(3)))
