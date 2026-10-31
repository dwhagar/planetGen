"""
The components' own pass (GEN.324): the bright and massive stars and the phenomena of the globular clusters and the
satellite galaxies (`galaxy.components`), placed where the shell walk and the object-first sampler cannot reach them.

A component is a density, not a list of stars: its mass over the mean mass of its systems (`Component.mean_system_mass_solar`) is its systems, and
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

A component with a great many stars (a Magellanic Cloud, at a low luminosity floor, holds a million) is drawn in slices,
each a unit of the queue's own: a slice has a stream of its own (the component's, then `:slice:` and its number) and
the share `1 / slices` of the component's Poisson mean, so the slices together are a draw of the whole component. How
many slices a component has comes from its expected stars alone (`parts`), never from the worker count, so a seed
draws the same stars on one worker or twenty. A component that needs no slices keeps its stream's name, and its stars.

The phenomena follow the same way: each class's rate per star (`phenomenon_scatter._scattered_rates`) times the
component's systems, at the same positions. A component's remnants are not the disk's: the thin disk's thicker remnant
layer does not apply. A cluster has no gas to speak of (no supernova remnants) and no young stars (no magnetars). A
satellite's classes follow its history (GEN.528, `galaxy.satellite_history`): its planetary nebulae scale with the share
of its stars that are not young against the solar circle's (as a sector's do, `profile.kind_factors_at`), its magnetars
with the young share of its systems, and its supernova remnants with its star formation against the solar circle's (GEN.1191, `supernova_remnant_factor`: none where it forms no stars,
more than the disk's per system where it forms them faster per star), so a Cloud shows every class its profile allows.
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.generation import bright_stars, phenomenon_reach, phenomenon_scatter
from planetgen.generation.bright_stars import MPC_PER_PC, POPULATIONS
from planetgen.galaxy import remnant_density, component_field, components as component_table, nebula_field, profile
from planetgen.galaxy.spherical import sector_address_at
from planetgen.physics import component_profiles
from planetgen.util import draw
from planetgen.util.poisson import poisson_count

STAR_SLICE_STARS = 25000.0
"""float: The most stars a slice of a component is expected to draw: about two minutes of a worker's time (the draw of
each star's properties, its members and the write of its row cost a few milliseconds)."""

PHENOMENON_SLICE_ROWS = 1000000.0
"""float: The most candidate phenomena a slice of a component draws; a phenomenon is cheap to draw (about 10 us) and
about a fifth of the candidates are kept."""

COMPONENT_KINDS = ("cluster", "satellite")
"""tuple: The kinds of component the pass draws; each is also the population of its stars (a satellite's three, a young
open cluster's one by its age)."""


def systems_of(component):
    """The systems (stars, counting a pair as one) a component holds: its mass over the mean system mass."""
    return component.expected_systems


def part_count(expected, size):
    """How many slices a component expected to draw `expected` rows is cut into: one for a component of `size` rows
    or fewer, the rows over `size` rounded up for a larger one."""
    return max(1, math.ceil(expected / size))


def parts(weighted, size):
    """
    The components cut into tasks, from their expected rows alone (the same at any worker count): a component that
    holds more than `size` is cut into slices, one task each, and the others are bundled in order, a task of them
    closed before the component that would take it past `size` rows, so no task holds more than `size`. A component expected to draw nothing is left out.

    Args:
        weighted (list): `(component, expected rows)`, in the order the tasks are to hold them.
        size (float): `STAR_SLICE_STARS` or `PHENOMENON_SLICE_ROWS`.

    Returns:
        list[tuple]: `(expected rows, parts)` per task, `parts` the `[kind, name, slice, slices]` lists that
            `chosen` turns back into components and their slices.
    """
    tasks, bundle, held = [], [], 0.0
    for component, expected in weighted:
        if expected <= 0.0:
            continue
        slices = part_count(expected, size)
        if slices > 1:
            tasks.extend((expected / slices, [[component.kind, component.name, index, slices]])
                         for index in range(slices))
            continue
        if bundle and held + expected > size:
            tasks.append((held, bundle))
            bundle, held = [], 0.0
        bundle.append([component.kind, component.name, 0, 1])
        held += expected
        if held >= size:
            tasks.append((held, bundle))
            bundle, held = [], 0.0
    if bundle:
        tasks.append((held, bundle))
    return tasks


def chosen(components, task_parts):
    """
    `(components, slices)` of one task: the stored `components` that `task_parts` (`parts`) names, in its order, and
    their slices as `scatter_stars` and `scatter_phenomena` take them.

    Raises:
        ValueError: A named component is not among `components` (it was removed after the run was cut).
    """
    by_name = {(component.kind, component.name): component for component in components}
    missing = [f"{kind} {name}" for kind, name, _index, _slices in task_parts if (kind, name) not in by_name]
    if missing:
        raise ValueError(f"The galaxy no longer has the component {missing[0]!r} this run was cut for; plan again.")
    return ([by_name[(kind, name)] for kind, name, _index, _slices in task_parts],
            {(kind, name): (index, slices) for kind, name, index, slices in task_parts})


def _stream_name(name, slices):
    """The stream of a component's slice: the component's own for a whole one."""
    index, count = slices
    return name if count == 1 else f"{name}:slice:{index}"


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


def _tiers_by_place(tiers_of_places, addresses):
    """`{(shell, band): tier}` for the bands of `addresses`, all looked up together (`tiers_of_places(places)`, an array)."""
    places = sorted({address[:2] for address in addresses})
    return dict(zip(places, tiers_of_places(places).tolist()))


def component_spots(component, edge_pc, by_tier, tiers_of_places, rng, skip_addresses, slices=1):
    """
    The stars of one component the pass keeps, population by population (a satellite's three, `galaxy.satellite_history`):
    Poisson in its systems times the population's share of its stars times the largest share of that population any
    tier keeps, each kept with probability its tier's share over that.

    Args:
        by_tier (list): `band_fractions` per tier, one for an untiered pass.
        tiers_of_places (callable or None): `tiers_of_places(places)` the tier of each `(shell, band)` (an array),
            `None` for one tier.
        slices (int): How many slices the component is drawn in; this one draws its share of them.

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
        points, (keep, _spare) = _candidates(component, systems_of(component) * share_of_stars * best / slices, rng)
        spots = [(tuple(int(value) for value in points[:, index]), index) for index in range(points.shape[1])]
        addresses = [_address(point, edge_pc) for point, _index in spots]
        tiers = None if tiers_of_places is None else _tiers_by_place(
            tiers_of_places, [address for address in addresses if address not in skip_addresses])
        for (point, index), address in zip(spots, addresses):
            if address in skip_addresses:
                continue
            tier = 0 if tiers is None else tiers[address[:2]]
            if keep[index] * best < shares[tier]:
                yield (*address, point, population, tier)


def scatter_stars(components, edge_pc, seed, fractions, pass_bands, tiered=None, skip_addresses=None, slices=None):
    """
    The pass's stars in the components (`COMPONENT_KINDS`) among `components`, one stream each, in
    `store.BRIGHT_STAR_COLUMNS` order. Arguments as `bright_stars.scatter_far`; `fractions` and `pass_bands` are
    `bright_stars.pass_of`'s. `slices` says which slice of a component to draw, `{(kind, name): (slice, slices)}`
    (`chosen`); a component it leaves out is drawn whole.
    """
    by_tier = fractions if tiered is not None else [fractions]
    tiers_of_places = None if tiered is None else tiered.limits.tiers_of_places
    skip_addresses = skip_addresses or set()
    for component in components:
        if component.kind not in COMPONENT_KINDS:
            continue
        here = (slices or {}).get((component.kind, component.name), (0, 1))
        rng = draw.Stream(f"{seed}:objects:component:{component.kind}:{_stream_name(component.name, here)}")
        placed = {(tier, population): [] for population in POPULATIONS for tier in range(len(pass_bands))}
        for shell_index, band_index, wedge, point, population, tier in component_spots(
                component, edge_pc, by_tier, tiers_of_places, rng, skip_addresses, here[1]):
            placed[(tier, population)].append((shell_index, band_index, wedge, point))
        yield from bright_stars.draw_placed(placed, pass_bands, rng, fe_h_peaks=component_field.fe_h_peaks_of(component))


def expected_stars(components, fractions, tiered=None):
    """The candidates `scatter_stars` draws and the most it places: each component's systems times the largest share of
    its population any tier keeps (an upper bound, the tiers differing by place; `kept_stars` is the stars placed)."""
    return sum(each for _component, each in star_weights(components, fractions, tiered))


def star_weights(components, fractions, tiered=None):
    """`(component, expected stars)` for each component of the pass's kinds, as `expected_stars` counts them."""
    by_tier = fractions if tiered is not None else [fractions]
    return [(component, sum(systems_of(component) * share * max(item[population] for item in by_tier)
                            for population, share in component_field.population_shares(component).items()))
            for component in components if component.kind in COMPONENT_KINDS]


KEPT_SAMPLE_POINTS = 256
"""int: The positions of a component's profile a component's kept stars are averaged over (`kept_stars`); a cluster lies
in one or a few cells and a satellite in a few hundred, so a few hundred evenly spread points give the tier shares to a
few percent."""


def _spread(count):
    """`(3, count)` uniforms spread evenly over the unit cube (the R2 sequence on the last two), the same at every call."""
    index = np.arange(count) + 0.5
    return np.stack([index / count, (index * 0.7548776662466927) % 1.0, (index * 0.5698402909980532) % 1.0])


def kept_stars(component, fractions, tiered, edge_pc):
    """
    About how many stars `scatter_stars` places for one component (GEN.678): its systems times, for each population, the
    share of its stars a pass keeps at the places the component's stars lie, `sum over populations of share * mean of
    the tier's share at KEPT_SAMPLE_POINTS positions of the profile`. `star_weights` is the number of candidates instead,
    every one at the share of the largest tier: the Poisson mean the candidates are drawn from, and the most a slice can
    hold. A component far from the galaxy lies in the faintest tier, which keeps a few percent of what the nearest does, so
    the candidates were 10 to 20 times the stars placed.

    Args:
        fractions (list or dict): `bright_stars.pass_of`'s, a list per tier for a tiered pass (`tiered`), else one dict.
        tiered: The pass's `TieredPass`, or `None`.
    """
    shares = component_field.population_shares(component)
    systems = systems_of(component)
    if tiered is None:
        return sum(systems * share * fractions[population] for population, share in shares.items() if share > 0.0)
    points = positions_mpc(component, _spread(KEPT_SAMPLE_POINTS))
    addresses = [_address(tuple(int(value) for value in points[:, index]), edge_pc)
                 for index in range(points.shape[1])]
    tier_of = _tiers_by_place(tiered.limits.tiers_of_places, addresses)
    tiers_here = [tier_of[address[:2]] for address in addresses]
    return sum(systems * share * sum(fractions[tier][population] for tier in tiers_here) / len(tiers_here)
               for population, share in shares.items() if share > 0.0)


def kept_weights(components, fractions, tiered, edge_pc):
    """`(component, expected stars kept)` for each component of the pass's kinds, as `kept_stars` counts them."""
    return [(component, kept_stars(component, fractions, tiered, edge_pc))
            for component in components if component.kind in COMPONENT_KINDS]


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
    rates it has (a satellite's gas-kept and young-born classes at their full rate, an upper bound). A young open cluster
    (GEN.244) places none: its remnants are the disk's, and its stars are a fraction of a percent of the disk's."""
    if component_table.is_open_cluster(component):
        return 0.0
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


def supernova_remnant_factor(component, shape):
    """
    The factor on a satellite's supernova remnants against the rate per system at the solar circle (GEN.1191): remnants follow the star formation. A
    core-collapse remnant comes from a young star and a Type Ia from one that is not, so the factor is the solar circle's odds of each
    (`tuning.SUPERNOVA_PROGENITOR_TYPE_IA_CHANCE`, 3 to 1) times the satellite's young share and old share against the solar circle's
    (`profile._reference_young_share`): 1 for the solar circle's mix, above 1 for a Cloud that forms stars faster per star than the disk, and
    0 for a satellite with no young stars, which forms none now and so holds none of the gas a remnant needs to be seen in (the old gas factor's zero). It is
    not capped, as the old gas factor was at 1, and not cut by a gas profile: the young stars of a satellite share its stars' profile, so the remnants are placed with them.
    """
    return profile.supernova_progenitor_factor(component_field.population_shares(component).get("satellite_young", 0.0), shape)


def _class_mean(component, kind, subtype, rate, shape):
    """The Poisson mean of one class in one component: its systems times the rate per system, times the component's
    factor for the class (GEN.528): the planetary nebulae's old-star share, the magnetars' young share of the
    satellite's systems, the stellar-mass black holes' metallicity factor of the component's populations (GEN.563). A satellite's supernova remnants take the rate times
    their star-formation factor (`supernova_remnant_factor`, GEN.1191)."""
    mean = systems_of(component) * rate
    if kind == "black-hole" and subtype == "stellar":
        mean *= remnant_density.component_factor(component_field.population_shares(component),
                                                    None if shape is None else shape.structure)
    if kind == phenomenon_scatter.OLD_STAR_KIND:
        mean *= planetary_nebula_factor(component, shape)
    elif phenomenon_reach._follows_young_stars(kind, subtype):
        mean *= component_field.population_shares(component).get("satellite_young", 0.0)
    elif kind in phenomenon_scatter.GAS_WEIGHTED_KINDS and component.kind == "satellite":
        mean *= supernova_remnant_factor(component, shape)
    return mean


def scatter_phenomena(components, edge_pc, seed, skip_addresses=None, min_mass_solar=None, shape=None, slices=None):
    """
    The phenomena in the components, in `PHENOMENON_SCATTER_COLUMNS` order: per component and class, Poisson in its
    systems times the class's rate per system (and its factor, `_class_mean`), at positions from the profile. A
    satellite's supernova remnants and magnetars follow the others in its stream, so the draws of the classes every
    component has do not move; the remnants' count scales with the satellite's star formation (GEN.1191).

    Args:
        shape (GalaxyShape): The galaxy's shape, with its components: the solar circle's young share comes from it. Needed for a satellite, and for the planetary nebulae of any component.
        slices (dict): Which slice of a component to draw, `{(kind, name): (slice, slices)}` (`chosen`), as
            `scatter_stars` takes them; a component it leaves out is drawn whole.
    """
    if min_mass_solar is None:
        min_mass_solar = tuning.PHENOMENON_MIN_MASS_SOLAR
    rates = phenomenon_rates(min_mass_solar)
    extra = satellite_rates(min_mass_solar)
    skip_addresses = skip_addresses or set()
    for component in components:
        if component.kind not in COMPONENT_KINDS or component_table.is_open_cluster(component):
            continue
        is_satellite = component.kind == "satellite"
        here = (slices or {}).get((component.kind, component.name), (0, 1))
        for kind, subtype, rate in rates + (extra if is_satellite else []):
            rng = phenomenon_scatter.class_stream(
                f"{seed}:phenomena:component:{component.kind}:{_stream_name(component.name, here)}", kind, subtype)
            points, (mass_draw, keep) = _candidates(
                component, _class_mean(component, kind, subtype, rate, shape) / here[1], rng)
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
