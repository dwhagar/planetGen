"""
Bright-star pre-placement: every star at least
`tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL` bright is drawn and
placed galaxy-wide right after `planetgen plan`, into `bright_stars`,
while every sector stays unfilled. A sector's fill later builds a full
system around each of its pre-placed stars (`fill_context`, used by
`generate.generate_sector`) and draws the rest of its systems from dimmer
stars only, so its expected total is unchanged.

The scatter is object first, never sector by sector (a full galaxy has
billions of cells): a layer draws one Poisson count of candidates from a
certified density bound, gives each a ring, a slot and a point, and keeps it
with probability true density over bound (`generation/object_first.py`,
PERF.58). Every cell inside the outline can hold one, however sparse
(GEN.78): the density never falls below the halo floor
(`tuning.MIN_RELATIVE_DENSITY`).

The backfill goes the other way: around the generated sectors it adds the
stars the scatter left out, sector by sector and by initial mass
(GEN.187, `backfill_mass_cells`): the sectors a face away from the filled
region down to 1 solar mass, the next ring down to 2, then 5, then 8
(`tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL`), each sector keeping the
stars below the luminosity the scatter already placed. A staged scatter's
luminosity bands (`backfill_cells`, GEN.44) still top up what a sector
lacks.

The scatter can go down in stages (`planetgen plan
--bright-stars-down-to`): a galaxy scattered at 500 Lsun can later add
only the band from, say, 100 up to (not including) 500, keeping every
star already placed. `galaxy_shape.bright_star_min_luminosity_sol` holds
the level reached so far, and `sector_stats.bright_level_sol` a sector's
own level where a backfill took it deeper (GEN.44). A sector already
filled gets none of the new band: its own systems were drawn below the
old level, so they already include stars that bright.

A sector's own draws (the backfill, and a band's top-up of a sector a
backfill took deeper) are split into fixed luminosity bands
(`canonical_bands`, eight to a decade), each drawn whole from its own
random stream and then cut to the range asked for (GEN.44). So a
sector's stars between two levels never depend on the steps taken to get
there: down to 1000 Lsun and later down to 500 gives exactly the stars
one draw down to 500 gives. The galaxy-wide layer scatter keeps one
stream per layer and band run: splitting it the same way would cost
about four times the draws (a narrow band is drawn by redrawing the
stars that overshoot it).
"""

import bisect
import collections
import math

import numpy as np

from planetgen.physics import constants
from planetgen import tuning
from planetgen.galaxy.density import population_densities
from planetgen.galaxy.geometry import (
    galaxy_to_local_pc,
    layer_bounds_pc,
    layer_center_z_pc,
    ring_bounds_pc,
    ring_radius_pc,
    ring_sector_count,
    sector_address_at,
    sector_position_pc,
)
from planetgen.galaxy.sector import _sample_poisson_count
from planetgen.generation import object_first
from planetgen.generation.star_population import (
    bright_band_fraction, mass_band_fraction, pick_population, placed_star_fraction, sample_bright_stars,
    sample_mass_band_stars,
)
from planetgen.physics.units import pc_to_ly
from planetgen.util import draw
from planetgen.util.poisson import poisson_count

POPULATIONS = ("young", "intermediate", "old", "bulge")
"""tuple: The stellar populations `galaxyDensity.population_densities`
splits a position's density into, each scattered on its own."""

ANGLE_BINS = 32
"""int: How many angle bins a ring's density is averaged over (fewer when
the ring has fewer slots). The spiral arms' density varies smoothly over
a bin this size."""

CANDIDATE_REPORT = 512
"""int: How many candidates a layer tries between progress reports (PERF.4)."""

SLOT_REDRAWS = 8
"""int: How many points a star tries before giving up when rounding to
whole milliparsecs keeps carrying it over its cell's edge."""

MPC_PER_PC = 1000

BANDS_PER_DECADE = 8
"""int: How many fixed luminosity bands (`canonical_bands`) a decade of
luminosity is split into: band `k` runs from `10 ** (k / 8)` to
`10 ** ((k + 1) / 8)` Lsun. A draw cut at a level between two edges
draws the band it falls in whole and drops the stars outside the range,
so finer bands waste less; 8 keeps that under a third of one band."""

TOP_BAND = 56
"""int: The last fixed band (`canonical_bands`), open-ended from
10^7 Lsun up: no star model reaches that, so it is almost always empty."""


def canonical_bands(min_luminosity_sol, max_luminosity_sol=None):
    """
    The fixed luminosity bands (GEN.44) a draw from `min_luminosity_sol`
    up to (not including) `max_luminosity_sol` (`None`: no limit) is made
    of, dimmest first, as `(k, low_sol, high_sol)`: band `k` covers
    `10 ** (k / BANDS_PER_DECADE)` to the next edge (`high_sol` `None` for
    the open top band, `TOP_BAND`); a band's low edge is raised to the
    brightest white dwarf, the dimmest threshold the star model takes.
    Each band is drawn whole from its own stream and cut to the range
    afterwards, so stars never depend on where earlier draws stopped.

    Raises:
        ValueError: For a threshold below the brightest white dwarf (as
            `band_fractions` does).
    """
    band_fractions(min_luminosity_sol, max_luminosity_sol)  # refuses a bad threshold first
    lowest = tuning.WD_LUMINOSITY_RANGE_SOL[1]
    k = min(math.floor(BANDS_PER_DECADE * math.log10(min_luminosity_sol) + 1e-9), TOP_BAND)
    while k > -10 ** 6 and 10 ** (k / BANDS_PER_DECADE) > min_luminosity_sol:
        k -= 1
    bands = []
    while True:
        low = max(10 ** (k / BANDS_PER_DECADE), lowest)
        if max_luminosity_sol is not None and low >= max_luminosity_sol:
            break
        high = None if k >= TOP_BAND else 10 ** ((k + 1) / BANDS_PER_DECADE)
        if high is None or high > low:
            bands.append((k, low, high))
        if high is None:
            break
        k += 1
    return bands


def _in_range(row, min_luminosity_sol, max_luminosity_sol):
    """Whether a drawn row's luminosity (column 12, watts) is in `[min,
    max)` Lsun."""
    watts = row[12]
    return (watts >= min_luminosity_sol * constants.SOLAR_LUMINOSITY
            and (max_luminosity_sol is None or watts < max_luminosity_sol * constants.SOLAR_LUMINOSITY))


def _densities(position_pc, shape):
    """`population_densities`, with a negative term (a shape with a
    negative bulge amplitude) counted as none."""
    return {population: max(density, 0.0) for population, density in population_densities(position_pc, shape).items()}


def _ring_bins(ring_index, layer_index, shape, expected_at_density_1, edge_pc):
    """Per angle bin, the population densities at the ring's centerline.
    No bin is zeroed for being sparse (GEN.78): every cell in the outline
    keeps its chance of a bright star."""
    slots = ring_sector_count(ring_index)
    count = min(slots, ANGLE_BINS)
    radius = ring_radius_pc(ring_index, edge_pc)
    z = layer_center_z_pc(layer_index, edge_pc)
    bins = []
    for k in range(count):
        theta = (k + 0.5) * 2 * math.pi / count
        point = (radius * math.cos(theta), radius * math.sin(theta), z)
        bins.append(_densities(point, shape))
    return slots, bins


def _place_one(rng, weights, ring_index, layer_index, slots, edge_pc):
    """A uniform point in a slot of a bin picked by `weights`, as `(slot,
    (x, y, z))` in whole milliparsecs (as stored), or `None` if rounding
    carried every try over its cell's edge (a point that does is redrawn)."""
    total = sum(weights)
    bin_width = 2 * math.pi / len(weights)
    slot_width = 2 * math.pi / slots
    for _ in range(SLOT_REDRAWS):
        pick = rng.random() * total
        for k, weight in enumerate(weights):
            if pick < weight:
                break
            pick -= weight
        else:
            # Rounding left `pick` past every bin: the last one with any
            # weight, never a trailing zero-weight one (TEST.24).
            k = max(index for index, weight in enumerate(weights) if weight > 0)
        theta = (k + rng.random()) * bin_width
        slot = min(int(theta / slot_width), slots - 1)
        stored = _point_in_cell(rng, ring_index, layer_index, slot, slots, edge_pc)
        if stored is not None:
            return slot, stored
    return None


def _point_in_cell(rng, ring_index, layer_index, slot, slots, edge_pc):
    """A uniform point in one cell in whole milliparsecs, or `None` when
    rounding carried it over the cell's edge."""
    stored = _draw_point_in_cell(rng, ring_index, layer_index, slot, slots, edge_pc)
    return stored if _stored_in_cell(stored, ring_index, layer_index, slot, edge_pc) else None


def _draw_point_in_cell(rng, ring_index, layer_index, slot, slots, edge_pc):
    """A uniform point in one cell in whole milliparsecs, not yet checked to have stayed in it after rounding
    (`_stored_in_cell`): the object-first scatter checks only the candidates it keeps."""
    slot_width = 2 * math.pi / slots
    r_low, r_high = ring_bounds_pc(ring_index, edge_pc)
    z_low, z_high = layer_bounds_pc(layer_index, edge_pc)
    theta = (slot + rng.random()) * slot_width
    # Uniform over the cell's area, which grows with radius.
    radius = math.sqrt(rng.uniform(r_low * r_low, r_high * r_high))
    point = (radius * math.cos(theta), radius * math.sin(theta), rng.uniform(z_low, z_high))
    return tuple(round(value * MPC_PER_PC) for value in point)


def _stored_in_cell(stored, ring_index, layer_index, slot, edge_pc):
    """Whether a stored point (whole milliparsecs) is still in the cell after rounding."""
    return sector_address_at(tuple(value / MPC_PER_PC for value in stored), edge_pc) == (ring_index, layer_index, slot)


def _finite_or_none(value):
    """A lifespan or phase end as stored: `NULL` for the white dwarf's
    infinity (a DOUBLE column cannot hold one; `star_params` reads it back)."""
    return value if value is None or math.isfinite(value) else None


def _row(ring_index, layer_index, slot, point, population, params, rng):
    """One star in `_db.BRIGHT_STAR_COLUMNS` order."""
    x, y, z = point
    return (
        ring_index, layer_index, slot,
        x, y, z,
        population, params["type"], params["yerkes_class"], params["mass_kg"],
        params["radius_km"], params["temperature_k"], params["luminosity_w"],
        params["age_gy"], _finite_or_none(params["lifespan_gy"]), params["initial_mass_sol"],
        _finite_or_none(params["phase_end_age_gy"]), rng.getrandbits(63),
    )


def scatter(shape, extents, edge_pc, expected_at_density_1, min_luminosity_sol, seed,
            skip_addresses=None, on_layer=None, max_luminosity_sol=None, mass_range=None):
    """
    Draws and places every bright star in the galaxy's outline, one
    layer after another (`scatter_layer`), or, with `max_luminosity_sol`,
    only the band below it (a staged scatter going one layer dimmer keeps
    the brighter stars already placed).

    Args:
        shape (GalaxyShape): The galaxy's shape.
        extents (list): `(layer_index, outer_ring_index)` per layer, as
            `galaxySkeleton.build_layer_extents` returns.
        edge_pc (float): The sector edge.
        expected_at_density_1 (float): Systems per sector at density 1.
        min_luminosity_sol (float): The threshold, in solar luminosities.
        seed (int): The scatter's seed; the same seed and outline give the
            same stars.
        skip_addresses (set, optional): `(ring, layer, slot)` cells to leave
            out (sectors already filled).
        on_layer (callable, optional): Called as `on_layer(done, total)`
            after each layer.
        max_luminosity_sol (float, optional): The band's upper limit
            (exclusive), the level already scattered; `None` for every
            star at or above the threshold.
        mass_range (tuple, optional): Only stars born in this `(low, high)`
            range of solar masses (an end `None` for no limit): the mass
            pass takes `(mass limit, None)`, the luminosity pass under it
            `(None, mass limit)`.

    Yields:
        tuple: One row per star, in `_db.BRIGHT_STAR_COLUMNS` order.
    """
    for done, (layer_index, outer_ring) in enumerate(extents, start=1):
        yield from scatter_layer(shape, layer_index, outer_ring, edge_pc, expected_at_density_1,
                                 min_luminosity_sol, seed, skip_addresses, max_luminosity_sol=max_luminosity_sol,
                                 mass_range=mass_range)
        if on_layer is not None:
            on_layer(done, len(extents))


def band_fractions(min_luminosity_sol, max_luminosity_sol=None, mass_range=None):
    """Per population, the share of its stars in the scatter's band
    (`bright_band_fraction`)."""
    return {population: bright_band_fraction(min_luminosity_sol, max_luminosity_sol, population, mass_range)
            for population in POPULATIONS}


RING_WEIGHT_STARS = 5.0
"""float: What walking one ring of a layer costs, in stars drawn: a
ring's density bins take about as long as placing and drawing five
stars (measured 2026-10-01: 0.2 ms a ring against about 40 us a star),
so an empty edge layer still counts for something (PERF.9). The phenomena
pass still walks its rings; the star passes do not (`OBJECT_FIRST_RING_WEIGHT_STARS`)."""

OBJECT_FIRST_RING_WEIGHT_STARS = 0.02
"""float: What one ring of a star layer costs, in stars drawn, with the
object-first sampler (PERF.58): its density bound, about 1.2 us a ring
against about 60 us a star."""

WEIGHT_RING_SAMPLES = 48
"""int: Rings sampled per layer by `layer_weight` (evenly spaced)."""

WEIGHT_ANGLE_BINS = 8
"""int: Angle bins per sampled ring in `layer_weight`."""


def layer_expected_stars(shape, layer_index, outer_ring, edge_pc, expected_at_density_1, fractions,
                         ring_samples=WEIGHT_RING_SAMPLES, angle_bins=WEIGHT_ANGLE_BINS):
    """
    About how many bright stars `scatter_layer` will place in one layer:
    the same per-ring expected counts it draws its Poisson counts from,
    on `ring_samples` evenly spaced rings with `angle_bins` bins each
    instead of every ring with `ANGLE_BINS`, so the whole galaxy takes a
    second or two rather than as long as a scatter (PERF.9).

    Args:
        fractions (dict): `band_fractions`.
    """
    rings = outer_ring + 1
    if rings <= 0:
        return 0.0
    count = min(rings, ring_samples)
    step = rings / count
    total_fraction = sum(fractions.values())
    if total_fraction <= 0.0:
        return 0.0
    z = layer_center_z_pc(layer_index, edge_pc)
    expected = 0.0
    for k in range(count):
        ring_index = min(int((k + 0.5) * step), rings - 1)
        slots = ring_sector_count(ring_index)
        bins = min(slots, angle_bins)
        radius = ring_radius_pc(ring_index, edge_pc)
        ring_expected = 0.0
        for j in range(bins):
            theta = (j + 0.5) * 2 * math.pi / bins
            point = (radius * math.cos(theta), radius * math.sin(theta), z)
            densities = _densities(point, shape)
            ring_expected += sum(densities[population] * fractions[population] for population in POPULATIONS)
        expected += expected_at_density_1 * slots / bins * ring_expected * step
    return expected


def layer_weight(shape, layer_index, outer_ring, edge_pc, expected_at_density_1, fractions):
    """
    The work one layer of the scatter is expected to take, in stars: its
    expected stars (`layer_expected_stars`) plus `OBJECT_FIRST_RING_WEIGHT_STARS`
    per ring bounded. The bright-star progress bar and its ETA count these, so
    the near-empty layers at the top and bottom of the disk no longer
    count as much as the dense ones in the middle (PERF.9).

    Returns:
        tuple: `(weight, expected_stars)`.
    """
    expected = layer_expected_stars(shape, layer_index, outer_ring, edge_pc, expected_at_density_1, fractions)
    return expected + OBJECT_FIRST_RING_WEIGHT_STARS * (outer_ring + 1), expected


def _spots(shape, layer_index, outer_ring, edge_pc, expected_at_density_1, fractions, rng, skip_addresses,
           check=False, on_candidates=None):
    """
    The stars of one layer, object first (PERF.58, `generation/object_first.py`):
    one Poisson count of candidates from the layer's density majorant, each
    given a ring, a slot and a point, and kept with probability true density
    over majorant. A star in a skipped sector is dropped (thinning); a sector
    takes up to its capacity (`object_first.sector_capacity`).

    Args:
        fractions (dict): `band_fractions`: the share of each population's stars in the band drawn.
        check (bool): Assert that the true density never exceeds the majorant.
        on_candidates (callable, optional): Called as `on_candidates(done, total, kept)` every few hundred
            candidates.

    Yields:
        tuple: `(ring_index, slot, point, population)`, `point` in whole milliparsecs.
    """
    if sum(fractions.values()) <= 0.0 or outer_ring < 0:
        return
    majorants = object_first.ring_majorants(shape, layer_index, outer_ring, edge_pc, fractions)
    slot_counts = [ring_sector_count(ring_index) for ring_index in range(outer_ring + 1)]
    cumulative = np.cumsum(expected_at_density_1 * np.array(slot_counts, dtype=float) * majorants)
    total_mean = float(cumulative[-1])
    if not total_mean > 0.0:
        return
    bounds = majorants.tolist()
    cumulative = cumulative.tolist()
    candidates = poisson_count(total_mean, rng)
    kept = {}
    taken = 0
    for done in range(candidates):
        if on_candidates is not None and done % CANDIDATE_REPORT == 0:
            on_candidates(done, candidates, taken)
        ring_index = min(bisect.bisect_right(cumulative, rng.random() * total_mean), outer_ring)
        slots = slot_counts[ring_index]
        slot = min(int(rng.random() * slots), slots - 1)
        point = _draw_point_in_cell(rng, ring_index, layer_index, slot, slots, edge_pc)
        densities = _densities(tuple(value / MPC_PER_PC for value in point), shape)
        weights = [densities[population] * fractions[population] for population in POPULATIONS]
        wanted = sum(weights)
        if check and wanted > bounds[ring_index] * (1.0 + 1e-9):
            raise AssertionError(f"majorant violated in ring {ring_index} of layer {layer_index}: "
                                 f"{wanted} > {bounds[ring_index]}")
        if rng.random() * bounds[ring_index] >= wanted:
            continue
        address = (ring_index, layer_index, slot)
        if address in skip_addresses or not _stored_in_cell(point, ring_index, layer_index, slot, edge_pc):
            continue
        held = kept.get(address, 0)
        if held:
            centre = _densities(sector_position_pc(ring_index, layer_index, slot, edge_pc), shape)
            expected = expected_at_density_1 * sum(centre[population] * fractions[population]
                                                    for population in POPULATIONS)
            if held >= object_first.sector_capacity(expected):
                continue
        kept[address] = held + 1
        taken += 1
        pick = rng.random() * wanted
        for population, weight in zip(POPULATIONS, weights):
            pick -= weight
            if pick < 0.0:
                break
        else:
            population = max(zip(POPULATIONS, weights), key=lambda item: item[1])[0]
        yield ring_index, slot, point, population
    if on_candidates is not None:
        on_candidates(candidates, candidates, taken)


def scatter_layer(shape, layer_index, outer_ring, edge_pc, expected_at_density_1, min_luminosity_sol, seed,
                  skip_addresses=None, max_luminosity_sol=None, on_progress=None, mass_range=None, check=False):
    """
    One layer of `scatter`: every bright star from ring 0 out to
    `outer_ring` at `layer_index`, placed object first (`_spots`, PERF.58).
    Each layer draws from its own random stream (the scatter's seed and the
    layer index), so layers can be drawn in any order, or side by side in
    worker processes (PERF.7), and still give the same stars.

    Args:
        on_progress (callable, optional): PERF.4: called as `on_progress(done,
            estimate)`, in stars. Stars are placed, then drawn (their type and
            luminosity), and count half at each step, so `done` reaches the
            layer's count when it is finished. `estimate` is the stars placed
            so far scaled up by the share of candidates tried, then the placed
            total once all are.
        check (bool): Assert that the density never exceeds the majorant (tests).

    Yields:
        tuple: One row per star, in `_db.BRIGHT_STAR_COLUMNS` order.
    """
    rng = draw.Stream(f"{seed}:objects:{layer_index}")
    skip_addresses = skip_addresses or set()
    fractions = band_fractions(min_luminosity_sol, max_luminosity_sol, mass_range)
    placed = {population: [] for population in POPULATIONS}

    def candidates_done(done, total, taken):
        if on_progress is not None and total:
            estimate = taken * total / done if done else float(total)
            on_progress(taken / 2, max(estimate, float(taken)))

    for ring_index, slot, point, population in _spots(
            shape, layer_index, outer_ring, edge_pc, expected_at_density_1, fractions, rng, skip_addresses,
            check=check, on_candidates=candidates_done):
        placed[population].append((ring_index, slot, point))
    total_placed = sum(len(spots) for spots in placed.values())
    drawn = 0
    for population, spots in placed.items():
        if not spots:
            continue
        stars = sample_bright_stars(len(spots), min_luminosity_sol, population, rng,
                                    max_luminosity_sol=max_luminosity_sol, mass_range=mass_range)
        for (ring_index, slot, point), params in zip(spots, stars):
            drawn += 1
            if on_progress is not None:
                on_progress((total_placed + drawn) / 2, total_placed)
            yield _row(ring_index, layer_index, slot, point, population, params, rng)


def backfill_cells(shape, addresses, edge_pc, expected_at_density_1, min_luminosity_sol, max_luminosity_sol, seed,
                   mass_range=None):
    """
    Draws and places every star in a luminosity band for a few cells (the
    sectors a backfill reaches, GEN.23, one by one since GEN.44): the band
    below what was already placed there, so no star is drawn twice.

    Per cell (every one, however sparse: GEN.78), fixed band (`canonical_bands`) and population, a
    Poisson count with mean `expected_at_density_1 * density * band
    share` at the cell's center (the per-cell rate `scatter` averages over
    a bin), each star uniform in the cell. Each cell and band draws from
    its own stream (the seed, the address and the band), cut to the band
    asked for afterwards, so a cell taken down to 1000 Lsun and later to
    500 holds exactly the stars one draw down to 500 gives (GEN.44).

    Args:
        shape (GalaxyShape): The galaxy's shape.
        addresses (iterable): `(ring, layer, slot)` cells to fill.
        edge_pc (float): The sector edge.
        expected_at_density_1 (float): Systems per sector at density 1.
        min_luminosity_sol (float): The band's floor (Lsun).
        max_luminosity_sol (float or None): Its ceiling, the level the
            cells were already filled to; `None` when nothing was placed.
        seed (int): The galaxy's bright-star seed.
        mass_range (tuple, optional): As in `scatter`: a galaxy whose mass
            pass placed every heavy star draws the lighter ones only,
            `(None, mass limit)`.

    Yields:
        tuple: One row per star, in `_db.BRIGHT_STAR_COLUMNS` order.
    """
    bands = canonical_bands(min_luminosity_sol, max_luminosity_sol)
    shares = {(k, population): bright_band_fraction(low, high, population, mass_range)
              for k, low, high in bands for population in POPULATIONS}
    for ring_index, layer_index, slot in addresses:
        center = sector_position_pc(ring_index, layer_index, slot, edge_pc)
        densities = _densities(center, shape)
        slots = ring_sector_count(ring_index)
        for k, low, high in bands:
            means = {population: expected_at_density_1 * densities[population] * shares[(k, population)]
                     for population in POPULATIONS}
            if not any(mean > 0.0 for mean in means.values()):
                continue
            rng = draw.Stream(f"{seed}:{ring_index}:{layer_index}:{slot}:{k}")
            for population in POPULATIONS:
                if means[population] <= 0.0:
                    continue
                count = _sample_poisson_count(means[population], rng=rng)
                if not count:
                    continue
                for params in sample_bright_stars(count, low, population, rng, max_luminosity_sol=high,
                                                mass_range=mass_range):
                    for _ in range(SLOT_REDRAWS):
                        point = _point_in_cell(rng, ring_index, layer_index, slot, slots, edge_pc)
                        if point is not None:
                            row = _row(ring_index, layer_index, slot, point, population, params, rng)
                            if _in_range(row, min_luminosity_sol, max_luminosity_sol):
                                yield row
                            break


def mass_band_edges(mass_limit_sol, ring_masses=None):
    """
    The fixed mass edges (GEN.187) the mass backfill's bands lie between,
    ascending: the ring masses (`tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL`)
    and the galaxy's mass limit (the mass pass placed every star from it up).
    """
    edges = set(ring_masses if ring_masses is not None else tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL)
    if mass_limit_sol is not None:
        edges.add(float(mass_limit_sol))
    return sorted(float(edge) for edge in edges)


def canonical_mass_bands(target_mass_sol, held_mass_sol, edges):
    """
    The fixed mass bands (GEN.187) a draw from `target_mass_sol` up to (not
    including) `held_mass_sol` (`None`: no limit, no mass pass placed any)
    is made of, lightest first, as `(low, high)` pairs (`high` `None` for the
    open top band). Each band is drawn whole from its own stream, so a
    sector's stars between two masses never depend on the steps taken to get
    there. `target_mass_sol` and `held_mass_sol` are among `edges`.

    Raises:
        ValueError: For a target that is not an edge.
    """
    if target_mass_sol not in edges:
        raise ValueError(f"{target_mass_sol:g} solar masses is not a mass edge ({edges})")
    bands = []
    for index, low in enumerate(edges):
        if low < target_mass_sol:
            continue
        if held_mass_sol is not None and low >= held_mass_sol:
            break
        high = edges[index + 1] if index + 1 < len(edges) else None
        if held_mass_sol is not None and (high is None or high > held_mass_sol):
            high = held_mass_sol
        bands.append((low, high))
    return bands


def backfill_mass_cells(shape, addresses, edge_pc, expected_at_density_1, target_mass_sol, held_mass_sol,
                        max_luminosity_sol, seed, edges):
    """
    Draws and places every star born with `target_mass_sol <= mass <
    held_mass_sol` for some cells (the sectors a mass backfill reaches,
    GEN.187, one by one): the stars below what was already placed there,
    so none is drawn twice.

    Per cell, fixed mass band (`canonical_mass_bands`) and population, a
    Poisson count with mean `expected_at_density_1 * density * band share`
    at the cell's center (`mass_band_fraction`: the living stars born in the
    band that are dimmer than `max_luminosity_sol`, the luminosity pass's
    floor, since the brighter ones were placed already), each star uniform
    in the cell. A band draws from its own stream (the seed, the address and
    the band's lower mass).

    Args:
        shape (GalaxyShape): The galaxy's shape.
        addresses (iterable): `(ring, layer, slot)` cells to fill.
        edge_pc (float): The sector edge.
        expected_at_density_1 (float): Systems per sector at density 1.
        target_mass_sol (float): The lightest mass to place.
        held_mass_sol (float or None): The mass the cells hold stars from
            already (the galaxy's mass limit, or an earlier backfill's);
            `None` when no mass pass placed any.
        max_luminosity_sol (float or None): The luminosity the scatter
            placed everything above; `None` for no luminosity pass.
        seed (int): The galaxy's bright-star seed.
        edges (list): `mass_band_edges`.

    Yields:
        tuple: One row per star, in `_db.BRIGHT_STAR_COLUMNS` order.
    """
    bands = canonical_mass_bands(target_mass_sol, held_mass_sol, edges)
    shares = {(low, population): mass_band_fraction(low, high, max_luminosity_sol, population)
              for low, high in bands for population in POPULATIONS}
    for ring_index, layer_index, slot in addresses:
        center = sector_position_pc(ring_index, layer_index, slot, edge_pc)
        densities = _densities(center, shape)
        slots = ring_sector_count(ring_index)
        for low, high in bands:
            means = {population: expected_at_density_1 * densities[population] * shares[(low, population)]
                     for population in POPULATIONS}
            if not any(mean > 0.0 for mean in means.values()):
                continue
            rng = draw.Stream(f"{seed}:{ring_index}:{layer_index}:{slot}:m{low:g}")
            for population in POPULATIONS:
                if means[population] <= 0.0:
                    continue
                count = _sample_poisson_count(means[population], rng=rng)
                if not count:
                    continue
                for params in sample_mass_band_stars(count, low, high, population, rng,
                                                     max_luminosity_sol=max_luminosity_sol):
                    for _ in range(SLOT_REDRAWS):
                        point = _point_in_cell(rng, ring_index, layer_index, slot, slots, edge_pc)
                        if point is not None:
                            yield _row(ring_index, layer_index, slot, point, population, params, rng)
                            break


def star_params(row):
    """A `bright_stars` row back as the `stellarEvolution.star_params`
    dict `Star.from_params` takes. A white dwarf (a mass backfill places
    them, GEN.187) stores its infinite lifespan and phase end as `NULL`."""
    white_dwarf = row["yerkes_class"] == "VII"
    lifespan, phase_end = row["lifespan_gy"], row["phase_end_age_gy"]
    if white_dwarf:
        lifespan = math.inf if lifespan is None else lifespan
        phase_end = math.inf if phase_end is None else phase_end
    return {
        "type": row["star_type"], "yerkes_class": row["yerkes_class"], "mass_kg": row["mass_kg"],
        "radius_km": row["radius_km"], "temperature_k": row["temperature_k"],
        "luminosity_w": row["luminosity_w"], "age_gy": row["age_gy"], "lifespan_gy": lifespan,
        "initial_mass_sol": row["initial_mass_sol"], "phase_end_age_gy": phase_end,
    }


def local_position_ly(row, center_pc):
    """A `bright_stars` row's stored point as a sector-local `(x, y, z)`
    in light-years (the frame `SpaceSector` places systems in)."""
    point = (row["position_x_mpc"] / MPC_PER_PC, row["position_y_mpc"] / MPC_PER_PC,
             row["position_z_mpc"] / MPC_PER_PC)
    return tuple(pc_to_ly(value) for value in galaxy_to_local_pc(center_pc, point))


class BelowCut(collections.namedtuple("BelowCut", "address min_mass_solar seed expected_stars")):
    """A galaxy sector's share of the phenomenon scatter's mass cut
    (GEN.168): its `(ring, layer, slot)` address, the cut in solar masses,
    the scatter's seed, and the stars the plan expects there (relative
    density times the calibration constant)."""

    __slots__ = ()


class FillContext:
    """
    What a galaxy-placed sector's fill needs beyond its arguments: the
    stellar population mix at its center, and its pre-placed bright stars
    with the threshold they were scattered at.

    Attributes:
        center_pc (tuple): The sector's center, galaxy-frame parsecs.
        densities (dict): `population_densities` at the center.
        bright_rows (list): The sector's unfilled `bright_stars` rows.
        min_luminosity_sol (float or None): The scatter's threshold;
            `None` when no scatter ran (no star is capped then).
        shape (GalaxyShape): The galaxy's shape (the molecular cloud
            field reads its gas from it, GEN.47).
        phenomena_scattered (bool): Whether the galaxy's phenomenon scatter
            ran (GEN.100): the sector then builds its black holes, neutron
            stars, planetary nebulae, supernova remnants and hypervelocity
            stars (and the nucleus) from `phenomenon_rows` and rolls none.
        phenomenon_rows (list): The sector's unbuilt `phenomenon_scatter`
            rows.
        below_cut (BelowCut or None): The scatter's mass cut, when it had
            one (GEN.167): the sector draws the neutron stars and black
            holes below it (GEN.168), and builds its scattered ones above it.
        star_mass_limit_sol (float or None): The mass every star at or above is
            already placed from (the scatter's mass pass, or lower where a
            mass backfill took this sector, GEN.187); the sector draws only
            lighter stars. `None` when nothing was placed by mass.
    """

    def __init__(self, center_pc, shape, bright_rows=(), min_luminosity_sol=None, phenomenon_rows=None,
                 below_cut=None, star_mass_limit_sol=None):
        self.star_mass_limit_sol = star_mass_limit_sol
        self.center_pc = center_pc
        self.shape = shape
        self.densities = _densities(center_pc, shape)
        self.bright_rows = list(bright_rows)
        self.min_luminosity_sol = min_luminosity_sol
        self.phenomena_scattered = phenomenon_rows is not None
        self.phenomenon_rows = list(phenomenon_rows or ())
        self.below_cut = below_cut

    def bright_share(self):
        """The share of this position's stars at or above the threshold
        (0 when no scatter ran): the part of the sector's expected count
        its pre-placed stars already stand for."""
        if self.min_luminosity_sol is None and self.star_mass_limit_sol is None:
            return 0.0
        total = sum(self.densities.values())
        if total <= 0.0:
            return 0.0
        return sum(density * placed_star_fraction(self.min_luminosity_sol, self.star_mass_limit_sol, population)
                   for population, density in self.densities.items()) / total

    def apply(self, system_config, rng=draw):
        """Gives one of the sector's own (dim) systems its population and,
        after a scatter, the luminosity cap."""
        system_config.POPULATION = pick_population(self.densities, rng)
        if self.min_luminosity_sol is not None:
            system_config.MAX_STAR_LUMINOSITY_SOL = self.min_luminosity_sol
        if self.star_mass_limit_sol is not None:
            system_config.MAX_STAR_MASS_SOL = self.star_mass_limit_sol
        return system_config
