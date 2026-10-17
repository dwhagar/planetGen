"""
The components' own pass (GEN.324): the bright and massive stars and the phenomena of the globular clusters and the
satellite galaxies (`galaxy.components`), placed where the shell walk and the object-first sampler cannot reach them.

A component is a density, not a list of stars: its mass over `tuning.CLUSTER_MEAN_SYSTEM_MASS_SOL` is its systems, and
its profile (`physics/component_profiles`) says where they lie. The scatters leave the clusters' and satellites' share of
the density out of their weights (`bright_stars._densities`), so this pass draws it: one Poisson count per component from
its systems times the share of a star of its population the pass keeps, each candidate a position by inverse transform of
the profile's mass table (a radius with that share of the mass inside it, a direction uniform on the sphere), kept with
probability its tier's share over the largest tier's. The point's own address then gives its shell, band and wedge, so
the tier (`ScatterLimits.tier_at`), the sectors already filled and the cell all follow from where the star is. A sector's
fill counts the same stars as already placed (`FillContext.bright_share`, which weighs each population's density by the
share the pass kept), so no pre-placed star is drawn twice.

Every component draws from a stream of its own (the scatter's seed and the component's kind and name), so one
component's stars do not depend on another's, on the order, or on the worker that draws them.

The phenomena follow the same way: each class's rate per star (`phenomenon_scatter._scattered_rates`) times the
component's systems, at the same positions. A component's remnants are not the disk's: the thin disk's thicker remnant
layer does not apply, a component has no gas to speak of (no supernova remnants) and no young stars (no magnetars).
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.generation import bright_stars, phenomenon_reach, phenomenon_scatter
from planetgen.generation.bright_stars import MPC_PER_PC, POPULATIONS
from planetgen.galaxy.spherical import sector_address_at
from planetgen.physics import component_profiles
from planetgen.util import draw
from planetgen.util.poisson import poisson_count

COMPONENT_KINDS = ("cluster", "satellite")
"""tuple: The kinds of component the pass draws; each is also the population of its stars."""


def systems_of(component):
    """The systems (stars, counting a pair as one) a component holds: its mass over the mean system mass."""
    return component.mass_solar / tuning.CLUSTER_MEAN_SYSTEM_MASS_SOL


def positions_mpc(component, uniform):
    """
    The positions of `uniform.shape[1]` of a component's systems, whole milliparsecs in the galaxy frame (an int64
    array of three rows): the radius is the profile's, by inverse transform of its table of the mass inside each radius
    (`uniform[0]`), the direction is uniform on the sphere (`uniform[1]`, `uniform[2]`).
    """
    table = component_profiles.table_for(component)
    radius = component_profiles.scale_radius_pc(component) * np.interp(uniform[0], table.mass_fraction, table.x)
    cos_colat = 2.0 * uniform[1] - 1.0
    in_plane = radius * np.sqrt(np.maximum(0.0, 1.0 - cos_colat * cos_colat))
    lon = uniform[2] * (2.0 * math.pi)
    centre = np.asarray(component.center_pc, dtype=float)
    offsets = np.stack([in_plane * np.cos(lon), in_plane * np.sin(lon), radius * cos_colat])
    return np.rint((centre[:, None] + offsets) * MPC_PER_PC).astype(np.int64)


def _candidates(component, mean, rng):
    """`(positions, uniform)` of the candidates of a Poisson `mean`: three rows of whole milliparsecs, and the uniform
    draws left for the caller's own use (a keep test, a mass), one row each."""
    count = poisson_count(mean, rng) if mean > 0.0 else 0
    uniform = phenomenon_scatter._uniforms(rng, 5 * count).reshape(5, count)
    return positions_mpc(component, uniform[:3]), uniform[3:]


def _address(point, edge_pc):
    return sector_address_at(tuple(float(value) / MPC_PER_PC for value in point), edge_pc)


def component_spots(component, edge_pc, by_tier, tier_of, rng, skip_addresses):
    """
    The stars of one component the pass keeps: Poisson in its systems times the largest share of its population any
    tier keeps, each kept with probability its tier's share over that.

    Args:
        by_tier (list): `band_fractions` per tier, one for an untiered pass.
        tier_of (callable or None): `tier_of(shell, band)` the tier of a band, `None` for one tier.

    Yields:
        tuple: `(shell, band, wedge, point, population, tier)`.
    """
    population = component.kind
    shares = [item[population] for item in by_tier]
    best = max(shares)
    if best <= 0.0:
        return
    points, (keep, _spare) = _candidates(component, systems_of(component) * best, rng)
    for index in range(points.shape[1]):
        point = tuple(int(value) for value in points[:, index])
        address = _address(point, edge_pc)
        if address in skip_addresses:
            continue
        tier = 0 if tier_of is None else tier_of(address[0], address[1])
        if keep[index] * best < shares[tier]:
            yield (*address, point, population, tier)


def scatter_stars(components, edge_pc, seed, fractions, pass_bands, tiered=None, skip_addresses=None):
    """
    The pass's stars in the components (`COMPONENT_KINDS`) among `components`, one stream each, in
    `store.BRIGHT_STAR_COLUMNS` order. Arguments as `bright_stars.scatter_far`; `fractions` and `pass_bands` are
    `bright_stars.pass_of`'s.
    """
    by_tier = fractions if tiered is not None else [fractions]
    tier_of = None if tiered is None else tiered.limits.tier_at
    skip_addresses = skip_addresses or set()
    for component in components:
        if component.kind not in COMPONENT_KINDS:
            continue
        rng = draw.Stream(f"{seed}:objects:component:{component.kind}:{component.name}")
        placed = {(tier, population): [] for population in POPULATIONS for tier in range(len(pass_bands))}
        for shell_index, band_index, wedge, point, population, tier in component_spots(
                component, edge_pc, by_tier, tier_of, rng, skip_addresses):
            placed[(tier, population)].append((shell_index, band_index, wedge, point))
        yield from bright_stars.draw_placed(placed, pass_bands, rng)


def expected_stars(components, fractions, tiered=None):
    """About how many stars `scatter_stars` places, the progress bar's weight: each component's systems times the
    largest share of its population any tier keeps (an upper bound, the tiers differing by place)."""
    by_tier = fractions if tiered is not None else [fractions]
    return sum(systems_of(component) * max(item[component.kind] for item in by_tier)
               for component in components if component.kind in COMPONENT_KINDS)


def phenomenon_rates(min_mass_solar):
    """`(kind, subtype, rate per system)` of the classes a component has: the black holes, neutron stars and planetary
    nebulae, not the supernova remnants (no gas) or the magnetars (no young stars)."""
    return [(kind, subtype, rate) for kind, subtype, rate in phenomenon_scatter._scattered_rates(min_mass_solar)
            if not phenomenon_reach._follows_young_stars(kind, subtype)
            and kind not in phenomenon_scatter.GAS_WEIGHTED_KINDS]


def scatter_phenomena(components, edge_pc, seed, skip_addresses=None, min_mass_solar=None):
    """
    The phenomena in the components, in `PHENOMENON_SCATTER_COLUMNS` order: per component and class, Poisson in its
    systems times the class's rate per system, at positions from the profile.
    """
    if min_mass_solar is None:
        min_mass_solar = tuning.PHENOMENON_MIN_MASS_SOLAR
    rates = phenomenon_rates(min_mass_solar)
    skip_addresses = skip_addresses or set()
    for component in components:
        if component.kind not in COMPONENT_KINDS:
            continue
        rng = draw.Stream(f"{seed}:phenomena:component:{component.kind}:{component.name}")
        for kind, subtype, rate in rates:
            points, (mass_draw, _spare) = _candidates(component, systems_of(component) * rate, rng)
            masses = phenomenon_scatter._masses(kind, subtype, min_mass_solar, mass_draw)
            for index in range(points.shape[1]):
                seed_here = rng.getrandbits(63)
                point = tuple(int(value) for value in points[:, index])
                address = _address(point, edge_pc)
                if address in skip_addresses:
                    continue
                mass = None if masses is None else float(masses[index])
                yield (*address, kind, subtype, *point, None, None, None, seed_here, mass,
                       phenomenon_scatter.map_luminosity(mass))
