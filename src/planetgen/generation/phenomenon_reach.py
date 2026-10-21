"""
The phenomena scatter beyond the galaxy's outline (GEN.323): black holes, neutron stars, planetary nebulae and
supernova remnants out to the grid's reach, placed object first like the bright stars' (GEN.322).

`phenomenon_scatter.scatter_shell` walks every band of a shell inside the outline, which cannot reach 262 kpc. The
part of a shell beyond the outline is cut into groups of up to 64 bands (`object_first.ReachGroups`), the faint
shells far out are drawn together as a few hundred pieces, and each class draws one Poisson count from the groups'
bounds on its weight: the class's rate per star times the density of the stars there, times its own factors (the
thicker layer of the remnants and the black holes' excess at the core, the gas of the remnants, the young stars of
the magnetars). A candidate lands in a group in proportion to `E * volume * bound` at a point uniform in the group,
and is kept with probability true weight over bound. The point's own address then gives its shell, band and wedge.

A class draws from a stream of its own per reach (`{seed}:phenomena:reach:{shell}:{kind}:{subtype}`, or `far:{first shell}`),
so the rows inside the outline are the rows they were with or without it, and a change to one class moves no other's. Rates, the mass cut, the masses and the row layout
are `phenomenon_scatter`'s. The stars a row carries with it (members) are written by the same writer as the rest.
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import remnant_density, density, profile, remnant_distribution
from planetgen.galaxy.spherical import sector_address_at
from planetgen.generation import bright_stars, object_first, phenomenon_scatter
from planetgen.generation.bright_stars import MPC_PER_PC, POPULATIONS
from planetgen.generation.phenomena import registry
from planetgen.util.poisson import poisson_count

CHUNK = 200_000
"""int: Candidates drawn at a time, so a reach with millions of them holds a few arrays, not all."""


def _streams(prefix):
    """`streams(kind, subtype)` for a reach: the random stream of each class, named by `prefix` (GEN.563)."""
    return lambda kind, subtype: phenomenon_scatter.class_stream(prefix, kind, subtype)


def _follows_young_stars(kind, subtype):
    return registry.scatter_class_gate(kind, subtype) == "young stars"


def _remnant_kind(kind, subtype):
    return not _follows_young_stars(kind, subtype) and kind in tuning.REMNANT_SCALE_HEIGHT_RATIO


def _follows_metals(kind, subtype):
    return remnant_density.follows_metals(kind, subtype)


def class_fractions(kind, subtype, structure):
    """
    The fractions of each population (`bright_stars.POPULATIONS`) that bound a class's weight over the whole galaxy
    per unit of density: the young stars alone for a class that follows them; the largest the remnant factors reach
    for a neutron star or a black hole, and for a stellar-mass black hole or a neutron star also the most its metallicity
    factor reaches (`remnant_density.bound_factor` of the galaxy's `structure`); 1 for the rest.
    """
    if _follows_young_stars(kind, subtype):
        return {population: 1.0 if population == "young" else 0.0 for population in POPULATIONS}
    factor = 1.0
    if _remnant_kind(kind, subtype):
        factor = tuning.REMNANT_VERTICAL_FACTOR_MAX
        if kind == "black-hole":
            factor *= 1.0 + tuning.BLACK_HOLE_CORE_EXCESS
    if _follows_metals(kind, subtype):
        factor *= remnant_density.bound_factor(structure, kind)
    return {population: factor for population in POPULATIONS}


def scattered_fractions(min_mass_solar, structure):
    """`class_fractions` summed over the classes, each times its rate per star above the cut: the bound on the
    reach's rows per unit of stellar density, which the expected-count estimates of `bright_stars` take as the
    fractions of one tier."""
    total = {population: 0.0 for population in POPULATIONS}
    for kind, subtype, rate in phenomenon_scatter._scattered_rates(min_mass_solar):
        for population, fraction in class_fractions(kind, subtype, structure).items():
            total[population] += rate * fraction
    return [total]


def class_bounds(shape, groups, kind, subtype):
    """
    An upper bound per group on a class's weight per unit of stellar density anywhere in it (the rate left out):
    the stellar density's bound at the group's nearest point (`object_first.extent_majorants`), times the remnant
    factors at the group's most favourable end: the vertical factor grows with the height, so its largest is at the
    group's greatest `|z|`, and the black holes' radial factor falls with the radius, so its largest is at the group's
    smallest. A class that follows the young stars is bounded by their density alone.
    """
    if _follows_young_stars(kind, subtype):
        fractions = class_fractions(kind, subtype, shape.structure)
        return object_first.extent_majorants(shape, groups.inner, groups.z_near, groups.z_far, fractions, groups.r_out,
                                             r_near=groups.r_in)
    metals = remnant_density.bound_factor(shape.structure, kind) if _follows_metals(kind, subtype) else 1.0
    ones = {population: metals for population in POPULATIONS}
    bound = object_first.extent_majorants(shape, groups.inner, groups.z_near, groups.z_far, ones, groups.r_out,
                                          r_near=groups.r_in)
    if _remnant_kind(kind, subtype):
        bound = (bound * remnant_distribution.vertical_factor_array(kind, groups.z_far, shape.disk_scale_height_pc)
                 * remnant_distribution.radial_factor_array(kind, groups.inner))
    return bound


def _weights(kind, subtype, x, y, z, shape):
    """A class's weight per unit of stellar density's rate at arrays of points (parsecs): what
    `phenomenon_scatter._band_rows` weighs a longitude bin by, at the points themselves."""
    if _follows_young_stars(kind, subtype):
        return np.maximum(density.population_densities_array(x, y, z, shape)["young"], 0.0)
    weight = (remnant_density.weights_array(x, y, z, shape, kind) if _follows_metals(kind, subtype)
              else density.total_density_array(x, y, z, shape))
    if _remnant_kind(kind, subtype):
        weight = (weight * remnant_distribution.vertical_factor_array(kind, z, shape.disk_scale_height_pc)
                  * remnant_distribution.radial_factor_array(kind, np.hypot(x, y)))
    if kind in phenomenon_scatter.GAS_WEIGHTED_KINDS:
        weight = weight * phenomenon_scatter._gas_weight(x, y, z, shape)
    if kind == phenomenon_scatter.OLD_STAR_KIND:
        weight = weight * profile.planetary_nebula_factor_array(x, y, z, shape)
    return weight


def _group_rows(shape, groups, high_band, edge_pc, expected_at_density_1, streams, rates, skip_addresses, min_mass_solar,
                check):
    """The rows of the regions `groups`, class by class in `SCATTER_CLASSES` order, each from the stream
    `streams(kind, subtype)` gives it."""
    if not len(groups):
        return
    cube_low, cube_high = groups.r_in ** 3, groups.r_out ** 3
    cos_low, cos_high = groups.cos_low, groups.cos_high
    shells_low, shells_high = groups.shell_low.tolist(), groups.shell_high.tolist()
    for kind, subtype, rate in rates:
        rng = streams(kind, subtype)
        bound = class_bounds(shape, groups, kind, subtype)
        cumulative = np.cumsum(expected_at_density_1 * rate * groups.volume * bound)
        total = float(cumulative[-1])
        candidates = poisson_count(total, rng) if total > 0.0 else 0
        for start in range(0, candidates, CHUNK):
            size = min(CHUNK, candidates - start)
            uniform = rng.random_array(6 * size).reshape(6, size)
            pick = np.minimum(np.searchsorted(cumulative, uniform[0] * total, side="right"), len(bound) - 1)
            radius = np.cbrt(cube_low[pick] + uniform[1] * (cube_high[pick] - cube_low[pick]))
            cos_colat = cos_low[pick] + uniform[2] * (cos_high[pick] - cos_low[pick])
            lon = uniform[3] * (2.0 * math.pi)
            in_plane = radius * np.sqrt(np.maximum(0.0, 1.0 - cos_colat * cos_colat))
            mpc = np.stack([np.rint(in_plane * np.cos(lon) * MPC_PER_PC), np.rint(in_plane * np.sin(lon) * MPC_PER_PC),
                            np.rint(radius * cos_colat * MPC_PER_PC)]).astype(np.int64)
            weight = _weights(kind, subtype, mpc[0] / MPC_PER_PC, mpc[1] / MPC_PER_PC, mpc[2] / MPC_PER_PC, shape)
            if check and np.any(weight > bound[pick] * (1.0 + 1e-9)):
                worst = int(np.argmax(weight / bound[pick]))
                raise AssertionError(f"majorant violated for {kind} {subtype} in group {int(pick[worst])}: "
                                     f"{weight[worst]} > {bound[pick][worst]}")
            kept = np.nonzero(uniform[4] * bound[pick] < weight)[0]
            masses = phenomenon_scatter._masses(kind, subtype, min_mass_solar, uniform[5][kept])
            seeds = rng.randbits_list(63, len(kept))
            for position, index in enumerate(kept.tolist()):
                point = mpc[:, index]
                address = sector_address_at(tuple(float(value) / MPC_PER_PC for value in point), edge_pc)
                shell_index, band_index, wedge = address
                group = int(pick[index])
                if (not shells_low[group] <= shell_index <= shells_high[group] or abs(band_index) <= high_band
                        or address in skip_addresses):
                    continue
                mass = None if masses is None else float(masses[position])
                yield (shell_index, band_index, wedge, kind, subtype, int(point[0]), int(point[1]), int(point[2]),
                       None, None, None, seeds[position], mass, phenomenon_scatter.map_luminosity(mass))


def scatter_shell(shape, shell_index, high_band, edge_pc, expected_at_density_1, seed, skip_addresses=None,
                  min_mass_solar=None, check=False):
    """
    The part of shell `shell_index` beyond the outline (the bands past `high_band`, all of them when it is `-1`),
    the phenomena of `phenomenon_scatter.scatter_shell` placed object first in groups of bands. Arguments as there;
    `check` asserts that no class's weight exceeds its bound (tests).

    Yields:
        tuple: One row per object, in `PHENOMENON_SCATTER_COLUMNS` order.
    """
    if min_mass_solar is None:
        min_mass_solar = tuning.PHENOMENON_MIN_MASS_SOLAR
    rates = phenomenon_scatter._scattered_rates(min_mass_solar)
    if not rates:
        return
    groups = object_first.ReachGroups.for_shell(shell_index, high_band, edge_pc)
    yield from _group_rows(shape, groups, high_band, edge_pc, expected_at_density_1,
                           _streams(f"{seed}:phenomena:reach:{shell_index}"), rates, skip_addresses or set(),
                           min_mass_solar, check)


def scatter_far(shape, first_shell, last_shell, edge_pc, expected_at_density_1, seed, skip_addresses=None,
                min_mass_solar=None, check=False):
    """The phenomena of shells `first_shell .. last_shell`, all beyond the outline and so far out that each holds
    almost none (`far_shells`), drawn together from the stream `{seed}:phenomena:far:{first_shell}`."""
    if min_mass_solar is None:
        min_mass_solar = tuning.PHENOMENON_MIN_MASS_SOLAR
    rates = phenomenon_scatter._scattered_rates(min_mass_solar)
    if not rates or last_shell < first_shell:
        return
    groups = object_first.ReachGroups.far(first_shell, last_shell, edge_pc)
    yield from _group_rows(shape, groups, -1, edge_pc, expected_at_density_1,
                           _streams(f"{seed}:phenomena:far:{first_shell}"), rates, skip_addresses or set(),
                           min_mass_solar, check)


def far_shells(shape, first_shell, last_shell, edge_pc, expected_at_density_1, min_mass_solar):
    """`bright_stars.far_shells` for the phenomena: the first shell of the run drawn together, and the expected
    objects of each shell before it."""
    return bright_stars.far_shells(shape, first_shell, last_shell, edge_pc, expected_at_density_1,
                                   scattered_fractions(min_mass_solar, shape.structure))


def reach_estimates(shape, shells, edge_pc, expected_at_density_1, min_mass_solar):
    """`bright_stars.reach_estimates` for the phenomena: `{shell: (objects, groups)}` beyond the outline."""
    return bright_stars.reach_estimates(shape, shells, edge_pc, expected_at_density_1,
                                        scattered_fractions(min_mass_solar, shape.structure))


def far_expected(shape, first_shell, last_shell, edge_pc, expected_at_density_1, min_mass_solar):
    """About how many objects `scatter_far` places."""
    return bright_stars.far_expected_stars(shape, first_shell, last_shell, edge_pc, expected_at_density_1,
                                           scattered_fractions(min_mass_solar, shape.structure))
