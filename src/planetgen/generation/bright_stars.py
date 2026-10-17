"""
Bright-star pre-placement: every star at least
`tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL` bright is drawn and
placed galaxy-wide right after `planetgen plan`, into `bright_stars`,
while every sector stays unfilled. A sector's fill later builds a full
system around each of its pre-placed stars (`fill_context`, used by
`generate.generate_sector`) and draws the rest of its systems from dimmer
stars only, so its expected total is unchanged.

The scatter is object first, never sector by sector (a full galaxy has
billions of cells): a shell draws one Poisson count of candidates from a
certified density bound, gives each a band, a wedge and a point, and keeps it
with probability true density over bound (`generation/object_first.py`,
PERF.58). Every cell inside the outline can hold one, however sparse
(GEN.78): the outline keeps only cells that expect at least one system, and
outside it the stellar halo (GEN.234, `galaxy.density.halo_terms`) is the
density, a measured spheroid and not a flat floor.

The backfill goes the other way: around the generated sectors it adds the
stars the scatter left out, sector by sector and by initial mass
(GEN.187, `backfill_mass_cells`): the sectors a face away from the filled
region down to 1 solar mass, the next sectors out down to 2, then 5, then 8
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
one draw down to 500 gives. The galaxy-wide shell scatter keeps one
stream per shell and band run: splitting it the same way would cost
about four times the draws (a narrow band is drawn by redrawing the
stars that overshoot it).
"""

import bisect
import collections
import functools
import math

import numpy as np

from planetgen.physics import constants, metallicity
from planetgen import tuning
from planetgen.galaxy import scatter_limits as scatterLimits
from planetgen.galaxy.component_field import KINDS_IN_THE_GATE
from planetgen.galaxy.density import population_densities, population_densities_array
from planetgen.galaxy.spherical import (
    band_center_pc,
    band_colatitude_bounds,
    cell_volume,
    galaxy_to_local_pc,
    sector_address_at,
    sector_position_pc,
    shell_bounds_pc,
    wedge_count,
    wedge_longitude_bounds,
)
from planetgen.galaxy.sector import _sample_poisson_count
from planetgen.generation import object_first, place_terms
from planetgen.generation.star_population import (
    bright_band_fraction, mass_band_fraction, pick_population, placed_star_fraction, sample_bright_stars,
    sample_mass_band_stars,
)
from planetgen.physics.units import pc_to_ly
from planetgen.util import draw
from planetgen.util.poisson import poisson_count

POPULATIONS = tuple(tuning.STELLAR_POPULATION_AGE_RANGES_GY)
"""tuple: The stellar populations `galaxy.density.population_densities`
splits a position's density into, each scattered on its own."""

ANGLE_BINS = 32
"""int: How many longitude bins a band's density is averaged over (fewer when
the band has fewer wedges). The spiral arms' density varies smoothly over
a bin this size."""

CANDIDATE_REPORT = 512
"""int: How many candidates a shell tries between progress reports (PERF.4)."""

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
    negative bulge amplitude) counted as none, and the clusters' and satellites' systems counted as none: the scatter
    pre-places no component star (they are filled, GEN.161, or scattered by their own pass, GEN.240), which is what
    `_densities_at` and the band majorants already assume."""
    found = {population: max(density, 0.0) for population, density in population_densities(position_pc, shape).items()}
    for kind in KINDS_IN_THE_GATE:
        found[kind] = 0.0
    return found


def _densities_at(points, z, shape):
    """`_densities` at each of `points` (`(x, y)` pairs at height `z`), as a list of dicts,
    from one array evaluation (PERF.89)."""
    if not points:
        return []
    xs, ys = zip(*points)
    arrays = population_densities_array(np.array(xs), np.array(ys), z, shape)
    columns = {population: np.maximum(array, 0.0).tolist() for population, array in arrays.items()}
    return [{population: column[i] for population, column in columns.items()} for i in range(len(points))]


def _band_bins(shell_index, band_index, high_band, shape, edge_pc, bins=ANGLE_BINS):
    """
    Per band `-high_band .. high_band` of a shell, the population densities at the band's centreline
    in `bins` longitude bins (fewer where the band has fewer wedges): `[(wedges, [densities per bin])]`.
    No bin is zeroed for being sparse (GEN.78): every cell in the outline keeps its chance of a bright star.
    """
    out = []
    for band in range(-high_band, high_band + 1):
        wedges = wedge_count(shell_index, band)
        count = min(wedges, bins)
        plane, z = band_center_pc(shell_index, band, edge_pc)
        points = [(plane * math.cos(theta), plane * math.sin(theta))
                  for theta in ((k + 0.5) * 2 * math.pi / count for k in range(count))]
        out.append((wedges, _densities_at(points, z, shape)))
    return out


@functools.lru_cache(maxsize=65536)
def volume_ratio(shell_index, band_index, edge_pc):
    """A cell of the band's volume in edge cubes (the cells of a shell differ a little, and near the poles a lot);
    every expected count scales with it, as a sector's own does."""
    return cell_volume(shell_index, band_index, 0, edge_pc) / edge_pc ** 3


def _point_in_cell(rng, shell_index, band_index, wedge, edge_pc):
    """A uniform point in one cell in whole milliparsecs, or `None` when
    rounding carried it over the cell's edge."""
    stored = _draw_point_in_cell(rng, shell_index, band_index, wedge, edge_pc)
    return stored if _stored_in_cell(stored, shell_index, band_index, wedge, edge_pc) else None


def _draw_point_in_cell(rng, shell_index, band_index, wedge, edge_pc):
    """A uniform point in one cell in whole milliparsecs, not yet checked to have stayed in it after rounding
    (`_stored_in_cell`): the object-first scatter checks only the candidates it keeps."""
    r_low, r_high = shell_bounds_pc(shell_index, edge_pc)
    north, south = band_colatitude_bounds(shell_index, band_index)
    start, end = wedge_longitude_bounds(shell_index, band_index, wedge)
    # Uniform over the cell's volume: the cube of the radius and the cosine of the colatitude are uniform.
    radius = rng.uniform(r_low ** 3, r_high ** 3) ** (1.0 / 3.0)
    cos_colat = rng.uniform(math.cos(south), math.cos(north))
    lon = rng.uniform(start, end)
    plane = radius * math.sqrt(max(0.0, 1.0 - cos_colat * cos_colat))
    point = (plane * math.cos(lon), plane * math.sin(lon), radius * cos_colat)
    return tuple(round(value * MPC_PER_PC) for value in point)


def _stored_in_cell(stored, shell_index, band_index, wedge, edge_pc):
    """Whether a stored point (whole milliparsecs) is still in the cell after rounding."""
    return sector_address_at(tuple(value / MPC_PER_PC for value in stored), edge_pc) == (shell_index, band_index, wedge)


def _finite_or_none(value):
    """A lifespan or phase end as stored: `NULL` for the white dwarf's
    infinity (a DOUBLE column cannot hold one; `star_params` reads it back)."""
    return value if value is None or math.isfinite(value) else None


def star_fe_h(population, point, seed, params=None):
    """A bright star's [Fe/H] (GEN.158): drawn from its population and its distance from the galactic centre (its
    stored point, whole milliparsecs), on a stream of its own seed so no scatter draw moves. `None` for a cluster
    or satellite star, which carries its component's value (GEN.159, GEN.237). A horizontal-branch star keeps the
    [Fe/H] its colour was drawn from (`params["metallicity_fe_h"]`, GEN.286)."""
    if population in ("cluster", "satellite"):
        return None
    if params is not None and params.get("metallicity_fe_h") is not None:
        return params["metallicity_fe_h"]
    distance_ly = pc_to_ly(math.sqrt(sum((value / MPC_PER_PC) ** 2 for value in point)))
    return metallicity.sample_fe_h(population, distance_ly, draw.Stream(f"metallicity:bright:{seed}"))


def _row(shell_index, band_index, wedge, point, population, params, rng):
    """One star in `store.BRIGHT_STAR_COLUMNS` order."""
    x, y, z = point
    seed = rng.getrandbits(63)
    return (
        shell_index, band_index, wedge,
        x, y, z,
        population, params["type"], params["yerkes_class"], params["mass_kg"],
        params["radius_km"], params["temperature_k"], params["luminosity_w"],
        params["age_gy"], _finite_or_none(params["lifespan_gy"]), params["initial_mass_sol"],
        _finite_or_none(params["phase_end_age_gy"]), seed, star_fe_h(population, point, seed, params),
    )


def scatter(shape, extents, edge_pc, expected_at_density_1, min_luminosity_sol, seed,
            skip_addresses=None, on_shell=None, max_luminosity_sol=None, mass_range=None):
    """
    Draws and places every bright star in the galaxy's outline, one
    shell after another (`scatter_shell`), or, with `max_luminosity_sol`,
    only the band below it (a staged scatter going one level dimmer keeps
    the brighter stars already placed).

    Args:
        shape (GalaxyShape): The galaxy's shape.
        extents (list): `(shell_index, high_band)` per shell, as
            `galaxy.skeleton.build_shell_extents` returns them.
        edge_pc (float): The sector edge.
        expected_at_density_1 (float): Systems per sector at density 1.
        min_luminosity_sol (float): The threshold, in solar luminosities.
        seed (int): The scatter's seed; the same seed and outline give the
            same stars.
        skip_addresses (set, optional): `(shell, band, wedge)` cells to leave
            out (sectors already filled).
        on_shell (callable, optional): Called as `on_shell(done, total)`
            after each shell.
        max_luminosity_sol (float, optional): The band's upper limit
            (exclusive), the level already scattered; `None` for every
            star at or above the threshold.
        mass_range (tuple, optional): Only stars born in this `(low, high)`
            range of solar masses (an end `None` for no limit): the mass
            pass takes `(mass limit, None)`, the luminosity pass under it
            `(None, mass limit)`.

    Yields:
        tuple: One row per star, in `store.BRIGHT_STAR_COLUMNS` order.
    """
    for done, (shell_index, high_band) in enumerate(extents, start=1):
        yield from scatter_shell(shape, shell_index, high_band, edge_pc, expected_at_density_1,
                                 min_luminosity_sol, seed, skip_addresses, max_luminosity_sol=max_luminosity_sol,
                                 mass_range=mass_range)
        if on_shell is not None:
            on_shell(done, len(extents))


def band_fractions(min_luminosity_sol, max_luminosity_sol=None, mass_range=None):
    """Per population, the share of its stars in the scatter's band
    (`bright_band_fraction`)."""
    return {population: bright_band_fraction(min_luminosity_sol, max_luminosity_sol, population, mass_range)
            for population in POPULATIONS}


BAND_WEIGHT_STARS = 5.0
"""float: What walking one band of a shell costs, in stars drawn: a
band's density bins take about as long as placing and drawing five
stars (measured 2026-10-01: 0.2 ms a band against about 40 us a star),
so an empty edge shell still counts for something (PERF.9). The phenomena
pass still walks its bands; the star passes do not (`OBJECT_FIRST_BAND_WEIGHT_STARS`)."""

OBJECT_FIRST_BAND_WEIGHT_STARS = 0.02
"""float: What one band of a star shell costs, in stars drawn, with the
object-first sampler (PERF.58): its density bound, about 1.2 us a band
against about 60 us a star."""

WEIGHT_BAND_SAMPLES = 48
"""int: Bands sampled per shell by `shell_weight` (evenly spaced)."""

WEIGHT_ANGLE_BINS = 8
"""int: Angle bins per sampled band in `shell_weight`."""


def tier_fractions_of(tiered):
    """The three `band_fractions` of a `TieredPass` (the core's, the baseline's and the halo's), a tier that
    draws nothing counting as none."""
    return [{population: 0.0 for population in POPULATIONS} if band is None
            else band_fractions(band.min_luminosity_sol, band.max_luminosity_sol, band.mass_range)
            for band in tiered.bands]


def sampled_band_expectations(shape, shell_index, high_band, edge_pc, expected_at_density_1, fractions,
                              band_samples=WEIGHT_BAND_SAMPLES, angle_bins=WEIGHT_ANGLE_BINS):
    """
    The expected bright stars of a shell on `band_samples` evenly spaced bands with `angle_bins` bins each:
    `[(band, plane_radius_pc, z_pc, expected)]`, each `expected` scaled up by the bands it stands for, so they
    add up to the shell's.

    `fractions` is one `band_fractions` dict, or (GEN.333) the three of a `TieredPass` (`tier_fractions_of`): then
    each sampled band takes the fractions of its own tier, chosen from the very densities the sample reads.
    """
    bands = 2 * high_band + 1
    by_tier = isinstance(fractions, (list, tuple))
    if high_band < 0 or sum(sum(item.values()) for item in (fractions if by_tier else [fractions])) <= 0.0:
        return []
    count = min(bands, band_samples)
    step = bands / count
    out = []
    for k in range(count):
        band = -high_band + min(int((k + 0.5) * step), bands - 1)
        wedges = wedge_count(shell_index, band)
        bins = min(wedges, angle_bins)
        plane, z = band_center_pc(shell_index, band, edge_pc)
        points = [(plane * math.cos(theta), plane * math.sin(theta))
                  for theta in ((j + 0.5) * 2 * math.pi / bins for j in range(bins))]
        band_expected = 0.0
        found = _densities_at(points, z, shape)
        here = fractions[scatterLimits.tier_of_densities(found)] if by_tier else fractions
        for densities in found:
            band_expected += sum(densities[population] * here[population] for population in POPULATIONS)
        out.append((band, plane, z, expected_at_density_1 * wedges * volume_ratio(shell_index, band, edge_pc) / bins
                    * band_expected * step))
    return out


def shell_expected_stars(shape, shell_index, high_band, edge_pc, expected_at_density_1, fractions,
                         band_samples=WEIGHT_BAND_SAMPLES, angle_bins=WEIGHT_ANGLE_BINS):
    """
    About how many bright stars `scatter_shell` will place in one shell:
    the same per-band expected counts it draws its Poisson counts from,
    on `band_samples` evenly spaced bands with `angle_bins` bins each, so the
    whole galaxy takes a second or two rather than as long as a scatter (PERF.9).

    Args:
        fractions (dict): `band_fractions`.
    """
    return sum(item[3] for item in sampled_band_expectations(
        shape, shell_index, high_band, edge_pc, expected_at_density_1, fractions, band_samples, angle_bins))


def shell_weight(shape, shell_index, high_band, edge_pc, expected_at_density_1, fractions, beyond=None):
    """
    The work one shell of the scatter is expected to take, in stars: its
    expected stars (`shell_expected_stars`) plus `OBJECT_FIRST_BAND_WEIGHT_STARS`
    per band bounded. The bright-star progress bar and its ETA count these, so
    the near-empty shells at the edge of the disk no longer
    count as much as the dense ones in the middle (PERF.9).

    With `beyond` (the shell's `reach_expected_stars`, or `reach_estimates`' entry), the part of the shell beyond the
    outline counts too: its expected stars plus `OBJECT_FIRST_BAND_WEIGHT_STARS` per group bounded.

    Returns:
        tuple: `(weight, expected_stars)`.
    """
    expected = shell_expected_stars(shape, shell_index, high_band, edge_pc, expected_at_density_1, fractions)
    weight = expected + OBJECT_FIRST_BAND_WEIGHT_STARS * max(2 * high_band + 1, 0)
    if beyond is not None:
        stars, groups = beyond
        expected += stars
        weight += stars + OBJECT_FIRST_BAND_WEIGHT_STARS * groups
    return weight, expected


def _capacity(shape, address, fractions, expected_at_density_1, edge_pc):
    """How many objects the sector at `address` may take (`object_first.sector_capacity`), from the density at its
    centre and the `fractions` of its band's tier."""
    shell_index, band_index, wedge = address
    centre = _densities(sector_position_pc(shell_index, band_index, wedge, edge_pc), shape)
    expected = (expected_at_density_1 * volume_ratio(shell_index, band_index, edge_pc)
                * sum(centre[population] * fractions[population] for population in POPULATIONS))
    return object_first.sector_capacity(expected)


def _pick_population(rng, weights, wanted):
    """The population of a kept candidate, in proportion to `weights` (which add up to `wanted`)."""
    pick_weight = rng.random() * wanted
    for population, weight in zip(POPULATIONS, weights):
        pick_weight -= weight
        if pick_weight < 0.0:
            return population
    return max(zip(POPULATIONS, weights), key=lambda item: item[1])[0]


def _spots(shape, shell_index, high_band, edge_pc, expected_at_density_1, fractions, rng, skip_addresses,
           check=False, on_candidates=None, tiers=None):
    """
    The stars of one shell, object first (PERF.58, `generation/object_first.py`):
    one Poisson count of candidates from the shell's density majorant, each
    given a band, a wedge and a point, and kept with probability true density
    over majorant. A star in a skipped sector is dropped (thinning); a sector
    takes up to its capacity (`object_first.sector_capacity`).

    Args:
        fractions (dict or list): `band_fractions`: the share of each population's stars in the band drawn. With
            `tiers` (GEN.333), a list of three, one per tier, and each band takes its own tier's.
        check (bool): Assert that the true density never exceeds the majorant.
        on_candidates (callable, optional): Called as `on_candidates(done, total, kept)` every few hundred
            candidates.
        tiers (sequence, optional): The tier index of each band `-high_band .. high_band`
            (`galaxy.scatter_limits.tier_indices`); `None` for one `fractions` everywhere.

    Yields:
        tuple: `(band_index, wedge, point, population, tier)`, `point` in whole milliparsecs, `tier` the band's
            tier index (0 without `tiers`).
    """
    by_tier = fractions if tiers is not None else [fractions]
    if high_band < 0 or sum(sum(item.values()) for item in by_tier) <= 0.0:
        return
    if tiers is None:
        majorants = object_first.band_majorants(shape, shell_index, high_band, edge_pc, fractions)
        tier_of_band = [0] * (2 * high_band + 1)
    else:
        tier_of_band = [int(tier) for tier in tiers]
        majorants = np.zeros(2 * high_band + 1)
        for tier in sorted(set(tier_of_band)):
            if sum(by_tier[tier].values()) > 0.0:
                bound = object_first.band_majorants(shape, shell_index, high_band, edge_pc, by_tier[tier])
                majorants = np.where(np.array(tier_of_band) == tier, bound, majorants)
    bands = list(range(-high_band, high_band + 1))
    wedge_counts = [wedge_count(shell_index, band) for band in bands]
    volumes = np.array([volume_ratio(shell_index, band, edge_pc) for band in bands])
    cumulative = np.cumsum(expected_at_density_1 * np.array(wedge_counts, dtype=float) * volumes * majorants)
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
        pick = min(bisect.bisect_right(cumulative, rng.random() * total_mean), len(bands) - 1)
        band_index = bands[pick]
        wedges = wedge_counts[pick]
        wedge = min(int(rng.random() * wedges), wedges - 1)
        point = _draw_point_in_cell(rng, shell_index, band_index, wedge, edge_pc)
        here = by_tier[tier_of_band[pick]]
        densities = _densities(tuple(value / MPC_PER_PC for value in point), shape)
        weights = [densities[population] * here[population] for population in POPULATIONS]
        wanted = sum(weights)
        if check and wanted > bounds[pick] * (1.0 + 1e-9):
            raise AssertionError(f"majorant violated in band {band_index} of shell {shell_index}: "
                                 f"{wanted} > {bounds[pick]}")
        if rng.random() * bounds[pick] >= wanted:
            continue
        address = (shell_index, band_index, wedge)
        if address in skip_addresses or not _stored_in_cell(point, shell_index, band_index, wedge, edge_pc):
            continue
        held = kept.get(address, 0)
        if held and held >= _capacity(shape, address, here, expected_at_density_1, edge_pc):
            continue
        kept[address] = held + 1
        taken += 1
        yield band_index, wedge, point, _pick_population(rng, weights, wanted), tier_of_band[pick]
    if on_candidates is not None:
        on_candidates(candidates, candidates, taken)


def reach_bounds(shape, groups, by_tier, coarse=False):
    """
    An upper bound per group of `groups` (`object_first.ReachGroups`) on `sum_p density_p * fractions[p]` anywhere
    in it (GEN.322): `object_first.extent_majorants` for the group's nearest point. `by_tier` is the `fractions` of
    each tier (one for an untiered pass): a group takes the fractions of the tiers its nearest and farthest places
    are judged to be in, the larger of each population's, since the tier changes along a group at most once and a
    band's own tier is known only when a candidate lands in it. A `coarse` bound (for the far shells, which hold
    almost nothing) skips the judging and takes the larger of each population's fractions over all the tiers.
    """
    if not len(groups):
        return np.zeros(0)
    if len(by_tier) == 1 or coarse:
        merged = {population: max(item[population] for item in by_tier) for population in POPULATIONS}
        return object_first.extent_majorants(shape, groups.inner, groups.z_near, groups.z_far, merged)
    near, far = groups.probes()
    near_tier = scatterLimits.probe_densities(shape, *near)[1]
    far_tier = scatterLimits.probe_densities(shape, *far)[1]
    bounds = np.zeros(len(groups))
    for pair in sorted(set(zip(near_tier.tolist(), far_tier.tolist()))):
        merged = {population: max(by_tier[pair[0]][population], by_tier[pair[1]][population])
                  for population in POPULATIONS}
        if sum(merged.values()) <= 0.0:
            continue
        mask = (near_tier == pair[0]) & (far_tier == pair[1])
        bounds[mask] = object_first.extent_majorants(
            shape, groups.inner[mask], groups.z_near[mask], groups.z_far[mask], merged)
    return bounds


def _expected_in_groups(shape, groups, expected_at_density_1, fractions):
    """About how many stars `groups` hold: the density at each group's middle times its volume, with the fractions of
    the tier it is in."""
    by_tier = fractions if isinstance(fractions, (list, tuple)) else [fractions]
    if not len(groups) or sum(sum(item.values()) for item in by_tier) <= 0.0:
        return 0.0
    r_mid = (groups.r_in + groups.r_out) / 2.0
    colat = (groups.colat_near + groups.colat_far) / 2.0
    means, tiers = scatterLimits.probe_densities(shape, r_mid * np.sin(colat), r_mid * np.cos(colat))
    if len(by_tier) == 1:
        tiers = np.zeros(len(groups), dtype=int)
    wanted = np.zeros(len(groups))
    for tier, here in enumerate(by_tier):
        wanted += np.where(tiers == tier, sum(means[population] * here[population] for population in POPULATIONS), 0.0)
    return float(expected_at_density_1 * np.sum(groups.volume * wanted))


def reach_expected_stars(shape, shell_index, high_band, edge_pc, expected_at_density_1, fractions):
    """
    About how many stars a shell holds beyond the outline (`ReachGroups`), and how many groups that is:
    `(stars, groups)`. What the progress bar and the expected-stars line read, so it is cheap, not exact.

    Args:
        fractions (dict or list): `band_fractions`, or the three of a `TieredPass`.
    """
    groups = object_first.ReachGroups.for_shell(shell_index, high_band, edge_pc)
    return _expected_in_groups(shape, groups, expected_at_density_1, fractions), len(groups)


REACH_ESTIMATE_STRIDE = 16
"""int: `reach_estimates` works out every this many shells' beyond-outline stars and takes them for the ones between."""


def reach_estimates(shape, shells, edge_pc, expected_at_density_1, fractions, stride=REACH_ESTIMATE_STRIDE):
    """
    `reach_expected_stars` for each of `shells` (`(shell, high_band)` pairs, consecutive shells), worked out on every
    `stride`th one and taken for the ones between, since the estimate varies slowly with the shell.

    Returns:
        dict: `{shell: (stars, groups)}`.
    """
    found, last = {}, None
    for index, (shell_index, high_band) in enumerate(shells):
        if index % stride == 0 or last is None:
            last = reach_expected_stars(shape, shell_index, high_band, edge_pc, expected_at_density_1, fractions)
        found[shell_index] = last
    return found


def far_expected_stars(shape, first_shell, last_shell, edge_pc, expected_at_density_1, fractions):
    """About how many stars shells `first_shell .. last_shell` hold, drawn together (`scatter_far`)."""
    if last_shell < first_shell:
        return 0.0
    return _expected_in_groups(shape, object_first.ReachGroups.far(first_shell, last_shell, edge_pc),
                               expected_at_density_1, fractions)


def far_shells(shape, first_shell, last_shell, edge_pc, expected_at_density_1, fractions):
    """
    Where the shells beyond the outline stop being worth bounding one by one (GEN.322): the coarse bound of each
    run of shells `first_shell .. last_shell` (`ReachGroups.far_blocks`, the larger fractions of every tier) gives the
    stars a shell of it can expect at most, and the shells from the first one after the last that can expect
    `tuning.OBJECT_FIRST_FAR_MEAN` or more are drawn together (`scatter_far`).

    Returns:
        tuple: `(first_far_shell, means)`: the first shell of the far run (`last_shell + 1` if none is far enough
            out, `first_shell` if every one is) and the coarse expected stars of each shell `first_shell ..
            first_far_shell - 1`.
    """
    if last_shell < first_shell:
        return first_shell, np.zeros(0)
    by_tier = fractions if isinstance(fractions, (list, tuple)) else [fractions]
    blocks = object_first.ReachGroups.far_blocks(first_shell, last_shell)
    groups = object_first.ReachGroups.far(first_shell, last_shell, edge_pc)
    weight = expected_at_density_1 * groups.volume * reach_bounds(shape, groups, by_tier, coarse=True)
    half, pieces = len(weight) // 2, tuning.OBJECT_FIRST_FAR_PIECES
    per_block = weight[:half].reshape(-1, pieces).sum(axis=1) + weight[half:].reshape(-1, pieces).sum(axis=1)
    widths = np.array([end - start + 1 for start, end in blocks])
    per_shell = per_block / widths
    near = np.nonzero(per_shell >= tuning.OBJECT_FIRST_FAR_MEAN)[0]
    first_far = first_shell if not len(near) else blocks[int(near[-1])][1] + 1
    return first_far, np.repeat(per_shell, widths)[:first_far - first_shell]


def _reach_spots(shape, groups, high_band, edge_pc, expected_at_density_1, by_tier, tier_of, rng, skip_addresses,
                 check=False, on_candidates=None, coarse=False):
    """
    The stars of the regions beyond the galaxy's outline (`ReachGroups`), object first (GEN.322): the same
    process as `_spots`, one Poisson count from the groups' bounds (`reach_bounds`), each candidate given a group in
    proportion to `E * volume * bound`, a point uniform in the group, and kept with probability true density over
    bound. The point's own address then decides its band, its tier (`tier_of(shell, band)`) and whether it is a
    sector's first or later object. Nothing here draws from the stream `_spots` does, so a shell's stars inside
    the outline are the same with or without it.

    Args:
        high_band (int): The highest band inside the outline of the groups' shell (`-1` for none): a point that
            falls in the outline's own bands (the milliparsec rounding of a point on the edge) is dropped.
        by_tier (list): `band_fractions` per tier, one for an untiered pass.
        tier_of (callable or None): `tier_of(shell, band)` the tier index of a band, `None` for one tier.
        rng: The stream of the regions' draws.
        coarse (bool): Bound the groups with `reach_bounds`' coarse bound.

    Yields:
        tuple: `(shell, band_index, wedge, point, population, tier)`.
    """
    if sum(sum(item.values()) for item in by_tier) <= 0.0:
        return
    bounds = reach_bounds(shape, groups, by_tier, coarse=coarse)
    cumulative = np.cumsum(expected_at_density_1 * groups.volume * bounds)
    total_mean = float(cumulative[-1]) if len(cumulative) else 0.0
    if not total_mean > 0.0:
        return
    bounds, cumulative = bounds.tolist(), cumulative.tolist()
    cube_low, cube_high = (groups.r_in ** 3).tolist(), (groups.r_out ** 3).tolist()
    cos_low, cos_high = groups.cos_low.tolist(), groups.cos_high.tolist()
    shells_low, shells_high = groups.shell_low.tolist(), groups.shell_high.tolist()
    candidates = poisson_count(total_mean, rng)
    kept = {}
    taken = 0
    for done in range(candidates):
        if on_candidates is not None and done % CANDIDATE_REPORT == 0:
            on_candidates(done, candidates, taken)
        pick = min(bisect.bisect_right(cumulative, rng.random() * total_mean), len(bounds) - 1)
        radius = rng.uniform(cube_low[pick], cube_high[pick]) ** (1.0 / 3.0)
        cos_colat = rng.uniform(cos_low[pick], cos_high[pick])
        lon = rng.uniform(0.0, 2.0 * math.pi)
        plane = radius * math.sqrt(max(0.0, 1.0 - cos_colat * cos_colat))
        point = tuple(round(value * MPC_PER_PC) for value in (plane * math.cos(lon), plane * math.sin(lon),
                                                                radius * cos_colat))
        address = sector_address_at(tuple(value / MPC_PER_PC for value in point), edge_pc)
        shell_index, band_index, wedge = address
        if not shells_low[pick] <= shell_index <= shells_high[pick] or abs(band_index) <= high_band:
            continue
        tier = 0 if tier_of is None else tier_of(shell_index, band_index)
        here = by_tier[tier]
        densities = _densities(tuple(value / MPC_PER_PC for value in point), shape)
        weights = [densities[population] * here[population] for population in POPULATIONS]
        wanted = sum(weights)
        if check and wanted > bounds[pick] * (1.0 + 1e-9):
            raise AssertionError(f"majorant violated in group {pick} of shell {shell_index}: "
                                 f"{wanted} > {bounds[pick]}")
        if rng.random() * bounds[pick] >= wanted or address in skip_addresses:
            continue
        held = kept.get(address, 0)
        if held and held >= _capacity(shape, address, here, expected_at_density_1, edge_pc):
            continue
        kept[address] = held + 1
        taken += 1
        yield shell_index, band_index, wedge, point, _pick_population(rng, weights, wanted), tier
    if on_candidates is not None:
        on_candidates(candidates, candidates, taken)


def _draw_groups(shape, groups, high_band, edge_pc, expected_at_density_1, rng, skip_addresses, fractions, pass_bands,
                 tiered, on_progress, base, check, coarse=False):
    """The rows of regions beyond the outline (`ReachGroups`): their stars placed (`_reach_spots`), then each star's
    type and luminosity drawn, from the one stream `rng`. `base` is the progress already reported, so `on_progress`
    carries on from it."""
    by_tier = fractions if tiered is not None else [fractions]
    tier_of = None if tiered is None else tiered.limits.tier_at
    placed = {(tier, population): [] for population in POPULATIONS for tier in range(len(pass_bands))}

    def candidates_done(done, total, taken):
        if on_progress is not None and total:
            estimate = taken * total / done if done else float(total)
            on_progress(base + taken / 2, base + max(estimate, float(taken)))

    for shell_index, band_index, wedge, point, population, tier in _reach_spots(
            shape, groups, high_band, edge_pc, expected_at_density_1, by_tier, tier_of, rng, skip_addresses,
            check=check, on_candidates=candidates_done, coarse=coarse):
        placed[(tier, population)].append((shell_index, band_index, wedge, point))
    total_placed = sum(len(spots) for spots in placed.values())
    drawn = 0
    for population in POPULATIONS:
        for tier, band in enumerate(pass_bands):
            spots = placed[(tier, population)]
            if not spots:
                continue
            stars = sample_bright_stars(len(spots), band.min_luminosity_sol, population, rng,
                                        max_luminosity_sol=band.max_luminosity_sol, mass_range=band.mass_range)
            for (shell_index, band_index, wedge, point), params in zip(spots, stars):
                drawn += 1
                if on_progress is not None:
                    on_progress(base + (total_placed + drawn) / 2, base + total_placed)
                yield _row(shell_index, band_index, wedge, point, population, params, rng)


def _pass_of(tiered, min_luminosity_sol, max_luminosity_sol, mass_range):
    """`(fractions, pass_bands)` of a pass: the three tiers' `band_fractions` and bands for a `TieredPass`, else the
    one pass's."""
    if tiered is None:
        return (band_fractions(min_luminosity_sol, max_luminosity_sol, mass_range),
                [scatterLimits.PassBand(min_luminosity_sol, max_luminosity_sol, mass_range)])
    return tier_fractions_of(tiered), list(tiered.bands)


def scatter_far(shape, first_shell, last_shell, edge_pc, expected_at_density_1, min_luminosity_sol, seed,
                skip_addresses=None, max_luminosity_sol=None, on_progress=None, mass_range=None, check=False,
                tiered=None):
    """
    The stars of shells `first_shell .. last_shell`, all of them beyond the galaxy's outline and far enough out that
    each holds almost nothing (`far_shells`), drawn together (GEN.322): one Poisson count from the coarse bounds of
    all their groups (`ReachGroups.far`), from the stream `{seed}:objects:far:{first}`. Arguments as `scatter_shell`.

    Yields:
        tuple: One row per star, in `store.BRIGHT_STAR_COLUMNS` order.
    """
    fractions, pass_bands = _pass_of(tiered, min_luminosity_sol, max_luminosity_sol, mass_range)
    groups = object_first.ReachGroups.far(first_shell, last_shell, edge_pc)
    rng = draw.Stream(f"{seed}:objects:far:{first_shell}")
    yield from _draw_groups(shape, groups, -1, edge_pc, expected_at_density_1, rng, skip_addresses or set(),
                            fractions, pass_bands, tiered, on_progress, 0.0, check, coarse=True)


def scatter_shell(shape, shell_index, high_band, edge_pc, expected_at_density_1, min_luminosity_sol, seed,
                  skip_addresses=None, max_luminosity_sol=None, on_progress=None, mass_range=None, check=False,
                  tiered=None, reach=False):
    """
    One shell of `scatter`: every bright star of shell `shell_index` in the
    bands `-high_band .. high_band`, placed object first (`_spots`, PERF.58).
    Each shell draws from its own random stream (the scatter's seed and the
    shell index), so shells can be drawn in any order, or side by side in
    worker processes (PERF.7), and still give the same stars.

    Args:
        on_progress (callable, optional): PERF.4: called as `on_progress(done,
            estimate)`, in stars. Stars are placed, then drawn (their type and
            luminosity), and count half at each step, so `done` reaches the
            shell's count when it is finished. `estimate` is the stars placed
            so far scaled up by the share of candidates tried, then the placed
            total once all are.
        check (bool): Assert that the density never exceeds the majorant (tests).
        tiered (galaxy.scatter_limits.TieredPass, optional): The pass by place (GEN.333): each band is drawn at
            its tier's limits and `min_luminosity_sol`, `max_luminosity_sol` and `mass_range` are not read. Without it
            one limit stands everywhere.
        reach (bool): Also draw the part of the shell beyond the outline (bands past `high_band`, all of them when
            `high_band` is `-1`) in groups of bands (`_reach_spots`, GEN.322). It has a stream of its own and is drawn
            after the outline's stars, which are the same either way.

    Yields:
        tuple: One row per star, in `store.BRIGHT_STAR_COLUMNS` order.
    """
    rng = draw.Stream(f"{seed}:objects:shell:{shell_index}")
    skip_addresses = skip_addresses or set()
    if tiered is None:
        fractions = band_fractions(min_luminosity_sol, max_luminosity_sol, mass_range)
        tiers = None
        pass_bands = [scatterLimits.PassBand(min_luminosity_sol, max_luminosity_sol, mass_range)]
    else:
        fractions = tier_fractions_of(tiered)
        tiers = tiered.tiers_in(shell_index, high_band)
        pass_bands = list(tiered.bands)
    placed = {(tier, population): [] for population in POPULATIONS for tier in range(len(pass_bands))}

    def candidates_done(done, total, taken):
        if on_progress is not None and total:
            estimate = taken * total / done if done else float(total)
            on_progress(taken / 2, max(estimate, float(taken)))

    for band_index, wedge, point, population, tier in _spots(
            shape, shell_index, high_band, edge_pc, expected_at_density_1, fractions, rng, skip_addresses,
            check=check, on_candidates=candidates_done, tiers=tiers):
        placed[(tier, population)].append((band_index, wedge, point))
    total_placed = sum(len(spots) for spots in placed.values())
    drawn = 0
    for population in POPULATIONS:
        for tier, band in enumerate(pass_bands):
            spots = placed[(tier, population)]
            if not spots:
                continue
            stars = sample_bright_stars(len(spots), band.min_luminosity_sol, population, rng,
                                        max_luminosity_sol=band.max_luminosity_sol, mass_range=band.mass_range)
            for (band_index, wedge, point), params in zip(spots, stars):
                drawn += 1
                if on_progress is not None:
                    on_progress((total_placed + drawn) / 2, total_placed)
                yield _row(shell_index, band_index, wedge, point, population, params, rng)
    if reach:
        groups = object_first.ReachGroups.for_shell(shell_index, high_band, edge_pc)
        yield from _draw_groups(shape, groups, high_band, edge_pc, expected_at_density_1,
                                draw.Stream(f"{seed}:objects:reach:{shell_index}"), skip_addresses, fractions,
                                pass_bands, tiered, on_progress, total_placed, check)


def backfill_cells(shape, addresses, edge_pc, expected_at_density_1, min_luminosity_sol, max_luminosity_sol, seed,
                   mass_range=None):
    """
    Draws and places every star in a luminosity band for a few cells (the
    sectors a backfill reaches, GEN.23, one by one since GEN.44): the band
    below what was already placed there, so no star is drawn twice.

    Per cell (every one, however sparse: GEN.78), fixed band (`canonical_bands`) and population, a
    Poisson count with mean `expected_at_density_1 * density * band
    share` at the cell's center, times the cell's volume in edge cubes (the per-cell rate `scatter` averages over
    a bin), each star uniform in the cell. Each cell and band draws from
    its own stream (the seed, the address and the band), cut to the band
    asked for afterwards, so a cell taken down to 1000 Lsun and later to
    500 holds exactly the stars one draw down to 500 gives (GEN.44).

    Args:
        shape (GalaxyShape): The galaxy's shape.
        addresses (iterable): `(shell, band, wedge)` cells to fill.
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
        tuple: One row per star, in `store.BRIGHT_STAR_COLUMNS` order.
    """
    bands = canonical_bands(min_luminosity_sol, max_luminosity_sol)
    shares = {(k, population): bright_band_fraction(low, high, population, mass_range)
              for k, low, high in bands for population in POPULATIONS}
    for shell_index, band_index, wedge in addresses:
        center = sector_position_pc(shell_index, band_index, wedge, edge_pc)
        densities = _densities(center, shape)
        volume = volume_ratio(shell_index, band_index, edge_pc)
        for k, low, high in bands:
            means = {population: expected_at_density_1 * volume * densities[population] * shares[(k, population)]
                     for population in POPULATIONS}
            if not any(mean > 0.0 for mean in means.values()):
                continue
            rng = draw.Stream(f"{seed}:{shell_index}:{band_index}:{wedge}:{k}")
            for population in POPULATIONS:
                if means[population] <= 0.0:
                    continue
                count = _sample_poisson_count(means[population], rng=rng)
                if not count:
                    continue
                for params in sample_bright_stars(count, low, population, rng, max_luminosity_sol=high,
                                                mass_range=mass_range):
                    for _ in range(SLOT_REDRAWS):
                        point = _point_in_cell(rng, shell_index, band_index, wedge, edge_pc)
                        if point is not None:
                            row = _row(shell_index, band_index, wedge, point, population, params, rng)
                            if _in_range(row, min_luminosity_sol, max_luminosity_sol):
                                yield row
                            break


def mass_band_edges(mass_limit_sol, ring_masses=None):
    """
    The fixed mass edges (GEN.187) the mass backfill's bands lie between,
    ascending: the ring masses, one per step away from the filled region (`tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL`)
    and the mass limit (the mass pass placed every star from it up): a number, or (GEN.333) the limits of each tier
    (`None` for a tier with none), so no band straddles one.
    """
    edges = set(ring_masses if ring_masses is not None else tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL)
    limits = mass_limit_sol if isinstance(mass_limit_sol, (list, tuple)) else [mass_limit_sol]
    edges.update(float(limit) for limit in limits if limit is not None)
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
        addresses (iterable): `(shell, band, wedge)` cells to fill.
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
        tuple: One row per star, in `store.BRIGHT_STAR_COLUMNS` order.
    """
    bands = canonical_mass_bands(target_mass_sol, held_mass_sol, edges)
    shares = {(low, population): mass_band_fraction(low, high, max_luminosity_sol, population)
              for low, high in bands for population in POPULATIONS}
    for shell_index, band_index, wedge in addresses:
        center = sector_position_pc(shell_index, band_index, wedge, edge_pc)
        densities = _densities(center, shape)
        volume = volume_ratio(shell_index, band_index, edge_pc)
        for low, high in bands:
            means = {population: expected_at_density_1 * volume * densities[population] * shares[(low, population)]
                     for population in POPULATIONS}
            if not any(mean > 0.0 for mean in means.values()):
                continue
            rng = draw.Stream(f"{seed}:{shell_index}:{band_index}:{wedge}:m{low:g}")
            for population in POPULATIONS:
                if means[population] <= 0.0:
                    continue
                count = _sample_poisson_count(means[population], rng=rng)
                if not count:
                    continue
                for params in sample_mass_band_stars(count, low, high, population, rng,
                                                     max_luminosity_sol=max_luminosity_sol):
                    for _ in range(SLOT_REDRAWS):
                        point = _point_in_cell(rng, shell_index, band_index, wedge, edge_pc)
                        if point is not None:
                            yield _row(shell_index, band_index, wedge, point, population, params, rng)
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
        "metallicity_fe_h": row.get("metallicity_fe_h"),
    }


def local_position_ly(row, center_pc):
    """A `bright_stars` row's stored point as a sector-local `(x, y, z)`
    in light-years (the frame `SpaceSector` places systems in)."""
    point = (row["position_x_mpc"] / MPC_PER_PC, row["position_y_mpc"] / MPC_PER_PC,
             row["position_z_mpc"] / MPC_PER_PC)
    return tuple(pc_to_ly(value) for value in galaxy_to_local_pc(center_pc, point))


class BelowCut(collections.namedtuple("BelowCut", "address min_mass_solar seed expected_stars")):
    """A galaxy sector's share of the phenomenon scatter's mass cut
    (GEN.168): its `(shell, band, wedge)` address, the cut in solar masses,
    the scatter's seed, and the stars the plan expects there (relative
    density times the calibration constant)."""

    __slots__ = ()


class FillContext:
    """
    What a galaxy-placed sector's fill needs beyond its arguments: its sector
    profile (`galaxy.profile`, GEN.235: the stellar population mix at its
    center and the rest), and its pre-placed bright stars with the threshold
    they were scattered at.

    Attributes:
        profile (SectorProfile): The sector's profile.
        center_pc (tuple): The sector's center, galaxy-frame parsecs.
        densities (dict): The profile's `mix`: `population_densities` at the center.
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
        resolution (CellResolution or None): What the cell resolves (GEN.281), `None` when no ceiling applies. A
            capped cell (`capped`) holds its pre-placed stars down to the fill luminosity cap and no others.
    """

    def __init__(self, profile, bright_rows=(), min_luminosity_sol=None, phenomenon_rows=None,
                 below_cut=None, star_mass_limit_sol=None, resolution=None):
        self.resolution = resolution
        self.star_mass_limit_sol = star_mass_limit_sol
        self.profile = profile
        self.center_pc = profile.position_pc
        self.shape = profile.shape
        self.densities = dict(profile.mix)
        self.bright_rows = list(bright_rows)
        self.min_luminosity_sol = min_luminosity_sol
        self.phenomena_scattered = phenomenon_rows is not None
        self.phenomenon_rows = list(phenomenon_rows or ())
        self.below_cut = below_cut

    @functools.cached_property
    def place_terms(self):
        """The profile's density and velocity dispersion for the pair draws (GEN.103, `place_terms.from_profile`)."""
        return place_terms.from_profile(self.profile)

    @property
    def capped(self):
        """Whether the cell expects more systems than the cell star ceiling (GEN.281): it is built from its
        pre-placed stars alone, and draws no dimmer system of its own."""
        return self.resolution is not None and self.resolution.capped

    def bright_share(self):
        """The share of this position's stars at or above the threshold
        (0 when no scatter ran): the part of the sector's expected count
        its pre-placed stars already stand for. A capped cell's pre-placed stars stand for all of what it makes."""
        if self.capped:
            return 1.0
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
        system_config.PLACE_TERMS = self.place_terms
        if system_config.POPULATION in ("cluster", "satellite"):
            system_config.FE_H_PEAKS = self.profile.populations[system_config.POPULATION].fe_h
        if self.min_luminosity_sol is not None:
            system_config.MAX_STAR_LUMINOSITY_SOL = self.min_luminosity_sol
        if self.star_mass_limit_sol is not None:
            system_config.MAX_STAR_MASS_SOL = self.star_mass_limit_sol
        return system_config
