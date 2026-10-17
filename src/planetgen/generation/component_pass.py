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
layer does not apply. A cluster has no gas to speak of (no supernova remnants) and no young stars (no magnetars). A
satellite's classes follow its history (GEN.528, `galaxy.satellite_history`): its planetary nebulae scale with the share
of its stars that are not young against the solar circle's (as a sector's do, `profile.kind_factors_at`), its magnetars
with the young share of its systems, and each candidate supernova remnant is kept with the gas fraction at its radius
(`nebula_field.gas_factor`, which counts a satellite's young stars), so a Cloud shows every class its profile allows.
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.generation import bright_stars, phenomenon_reach, phenomenon_scatter
from planetgen.generation.bright_stars import MPC_PER_PC, POPULATIONS
from planetgen.galaxy import component_field, nebula_field, profile
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
    uniform = rng.random_array(5 * count).reshape(5, count)
    return positions_mpc(component, uniform[:3]), uniform[3:]


def _address(point, edge_pc):
    return sector_address_at(tuple(float(value) / MPC_PER_PC for value in point), edge_pc)


def component_spots(component, edge_pc, by_tier, tier_of, rng, skip_addresses):
    """
    The stars of one component the pass keeps, population by population (a satellite's three, `galaxy.satellite_history`):
    Poisson in its systems times the population's share of its stars times the largest share of that population any
    tier keeps, each kept with probability its tier's share over that.

    Args:
        by_tier (list): `band_fractions` per tier, one for an untiered pass.
        tier_of (callable or None): `tier_of(shell, band)` the tier of a band, `None` for one tier.

    Yields:
        tuple: `(shell, band, wedge, point, population, tier)`; `population` is the cluster's or one of the satellite's
            three.
    """
    for population, share_of_stars in component_field.population_shares(component).items():
        if share_of_stars <= 0.0:
            continue
        shares = [item[population] for item in by_tier]
        best = max(shares)
        if best <= 0.0:
            continue
        points, (keep, _spare) = _candidates(component, systems_of(component) * share_of_stars * best, rng)
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
    return sum(systems_of(component) * share * max(item[population] for item in by_tier)
               for component in components if component.kind in COMPONENT_KINDS
               for population, share in component_field.population_shares(component).items())


def phenomenon_rates(min_mass_solar):
    """`(kind, subtype, rate per system)` of the classes every component has: the black holes, neutron stars and
    planetary nebulae, not the supernova remnants (gas) or the magnetars (young stars), which only a satellite can have
    (`satellite_rates`)."""
    return [(kind, subtype, rate) for kind, subtype, rate in phenomenon_scatter._scattered_rates(min_mass_solar)
            if not phenomenon_reach._follows_young_stars(kind, subtype)
            and kind not in phenomenon_scatter.GAS_WEIGHTED_KINDS]


def satellite_rates(min_mass_solar):
    """`(kind, subtype, rate per system)` of the classes only a satellite has (GEN.528): the supernova remnants, kept
    where its gas is, and the magnetars, born of its young stars."""
    return [(kind, subtype, rate) for kind, subtype, rate in phenomenon_scatter._scattered_rates(min_mass_solar)
            if phenomenon_reach._follows_young_stars(kind, subtype) or kind in phenomenon_scatter.GAS_WEIGHTED_KINDS]


def expected_phenomena(component, min_mass_solar):
    """About how many phenomena the pass places for one component, the progress bar's weight: its systems times the
    rates it has (a satellite's gas-kept and young-born classes at their full rate, an upper bound)."""
    rates = phenomenon_rates(min_mass_solar)
    if component.kind == "satellite":
        rates = rates + satellite_rates(min_mass_solar)
    return systems_of(component) * sum(rate for _kind, _subtype, rate in rates)


def planetary_nebula_factor(component, shape):
    """
    The factor on a component's planetary nebulae (GEN.528), the one `profile.kind_factors_at` gives a sector: the
    share of its stars that are not young against the solar circle's (`profile._reference_young_share`), capped at 1.
    A cluster has no young stars, so its factor is 1.
    """
    young = component_field.population_shares(component).get("satellite_young", 0.0)
    if young <= 0.0:
        return 1.0
    reference = profile._reference_young_share(shape)
    if reference >= 1.0:
        return 0.0
    return min(1.0, (1.0 - young) / (1.0 - reference))


GAS_PROFILE_POINTS = 96
"""int: The radii a satellite's gas fraction is tabulated at, evenly from its centre to its tidal radius."""


def gas_by_radius(component, shape):
    """
    `(radii_pc, factor)`: the gas fraction (`nebula_field.gas_factor`, capped at 1 as `SectorProfile.kind_factors` has
    it) along a ray from a satellite's centre. The satellite is spherical and lies away from the disk, so one ray gives
    the fraction at any radius, and a table of it lets a whole component's candidates be kept at once.
    """
    radii = np.linspace(0.0, component.tidal_radius_pc, GAS_PROFILE_POINTS)
    centre = component.center_pc
    factor = np.array([min(1.0, max(nebula_field.gas_factor((centre[0] + r, centre[1], centre[2]), shape), 0.0))
                       for r in radii])
    return radii, factor


def _class_mean(component, kind, subtype, rate, shape):
    """The Poisson mean of one class in one component: its systems times the rate per system, times the component's
    factor for the class (GEN.528): the planetary nebulae's old-star share, the magnetars' young share of the
    satellite's systems. A satellite's supernova remnants take the full rate and are kept by the gas fraction
    (`scatter_phenomena`)."""
    mean = systems_of(component) * rate
    if kind == phenomenon_scatter.OLD_STAR_KIND:
        mean *= planetary_nebula_factor(component, shape)
    elif phenomenon_reach._follows_young_stars(kind, subtype):
        mean *= component_field.population_shares(component).get("satellite_young", 0.0)
    return mean


def scatter_phenomena(components, edge_pc, seed, skip_addresses=None, min_mass_solar=None, shape=None):
    """
    The phenomena in the components, in `PHENOMENON_SCATTER_COLUMNS` order: per component and class, Poisson in its
    systems times the class's rate per system (and its factor, `_class_mean`), at positions from the profile. A
    satellite's supernova remnants and magnetars follow the others in its stream, so the draws of the classes every
    component has do not move; each remnant candidate is kept with the gas fraction at its radius.

    Args:
        shape (GalaxyShape): The galaxy's shape, with its components: the solar circle's young share and a satellite's
            gas (`gas_by_radius`) come from it. Needed for a satellite, and for the planetary nebulae of any component.
    """
    if min_mass_solar is None:
        min_mass_solar = tuning.PHENOMENON_MIN_MASS_SOLAR
    rates = phenomenon_rates(min_mass_solar)
    extra = satellite_rates(min_mass_solar)
    skip_addresses = skip_addresses or set()
    for component in components:
        if component.kind not in COMPONENT_KINDS:
            continue
        rng = draw.Stream(f"{seed}:phenomena:component:{component.kind}:{component.name}")
        gas = gas_by_radius(component, shape) if component.kind == "satellite" else None
        for kind, subtype, rate in rates + (extra if gas is not None else []):
            points, (mass_draw, keep) = _candidates(component, _class_mean(component, kind, subtype, rate, shape), rng)
            masses = phenomenon_scatter._masses(kind, subtype, min_mass_solar, mass_draw)
            if kind in phenomenon_scatter.GAS_WEIGHTED_KINDS and points.shape[1]:
                offsets = points / MPC_PER_PC - np.asarray(component.center_pc, dtype=float)[:, None]
                radius = np.sqrt((offsets ** 2).sum(axis=0))
                kept = keep < np.interp(radius, *gas)
            else:
                kept = np.ones(points.shape[1], dtype=bool)
            for index in range(points.shape[1]):
                seed_here = rng.getrandbits(63)
                if not kept[index]:
                    continue
                point = tuple(int(value) for value in points[:, index])
                address = _address(point, edge_pc)
                if address in skip_addresses:
                    continue
                mass = None if masses is None else float(masses[index])
                yield (*address, kind, subtype, *point, None, None, None, seed_here, mass,
                       phenomenon_scatter.map_luminosity(mass))
