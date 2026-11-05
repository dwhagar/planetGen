# generation/sky_tiers.py

"""
The sky tier table (VIEW.24): how deep the backfill must reach, in cells, so that every star a system's sky shows at a limiting magnitude is stored.

A star shows from distance d when its apparent magnitude `M_V + 5 log10(d / 10 pc) + A_V(d)` is at or under the limit. That makes a *floor* for each distance:
the dimmest star that still shows there, `floor_luminosity_sol`. The floor is worked out for the star whose colour gives the most visual light for its
luminosity (the bolometric correction's maximum, `BEST_CORRECTION`, +0.03 near 7,000 K), so no star that could show is left out, and for the direction
of least dust, so no direction is short.

The two ledgers of the backfill (VIEW.19, reports/system-sky-backfill-cost-2026-11-05.md section 5) split the floor:

* **Luminosity ledger, a floor of 50 Lsun and up.** The fixed luminosity bands of `bright_stars.canonical_bands` from `RING_TOP_SOL` (50 Lsun, the brightest
  white dwarf, under which the luminosity backfill refuses a floor) to the scatter's own floor. A band shows out to the reach of its brightest member
  (`reach_pc`), and is drawn over every cell inside that sphere.
* **Mass ledger, a floor under 50 Lsun.** The mass rings (`bright_stars.mass_band_edges`: 1, 2, 5, 8 and the mass limit) draw the stars born at 1 solar mass and up
  that are dimmer than 50 Lsun, so every ring reaches as far as a star of 50 Lsun does. The rings under 1 solar mass (`NEAR_RING_MASSES_SOL`, VIEW.41) draw
  the stars born lighter than the Sun, which are faint: a ring reaches only as far as `RING_MIN_YIELD` of its stars still show, tens of parsecs.

The bright far bands hold few stars in an enormous number of cells (a 7.5 sky is 3e9 cells in its farthest band), and since VIEW.41 they are drawn by
blocks (`bright_stars.backfill_sphere`), so they cost the blocks and rows of their sphere, not its cells (`block_visits`). The table is cut by a budget of
those visits (`BLOCK_VISIT_BUDGET`): bands are taken from the 50 Lsun band up while the visits fit, and the rest are listed as `dropped`. The stars of a
dropped band that the scatter did not pre-place are not stored; the glow band (VIEW.27) carries their light.

    curve = extinction_curve(observer_pc)                       # the dust along 48 directions, once per observer
    tiers = build(observer_pc, 7.5, edge_pc, curve=curve)
    tiers.luminosity, tiers.rings, tiers.dropped, visits(tiers)
    floor_table(tiers, curve, (10, 50, 250, 1000))              # the floor against distance
"""

import collections
import functools
import math

import numpy as np

from planetgen import tuning
from planetgen.generation import star_population
from planetgen.generation import bright_stars
from planetgen.physics import constants, photometry
from planetgen.physics import sky_photometry as photo
from planetgen.util import draw

BLOCK_VISIT_BUDGET = 1.0e8
"""float: The visits the luminosity bands may take for one sky (stand-in SKT1): a visit is a block or a row of the sphere draw (`block_visits`), about
`BLOCK_VISIT_SECONDS` of one worker, so this is about an hour of one worker for the bands, more than a dust-free 7.5 sky takes (4e7)."""

BLOCK_VISIT_SECONDS = 4.0e-5
"""float: What the sphere draw spends on a block or a row, one worker, warm, in the far bands that hold most of the visits (stand-in SKT3;
reports/system-sky-sphere-draw-2026-11-05.md): the near bands, whose cells hold more candidates, spend two to six times as much a visit (21 microseconds a visit at 3,500 pc, 246 at 250 pc)."""

RING_TOP_SOL = float(tuning.WD_LUMINOSITY_RANGE_SOL[1])
"""float: The luminosity the rings stop at and the luminosity bands start at (the brightest white dwarf, the lowest floor the luminosity backfill takes;
stand-in SKT2). A ring reaches as far as a star of this luminosity does."""

DISTANCE_GRID_PC = np.geomspace(0.5, 60000.0, 360)
"""numpy.ndarray: The distances the dust is worked out at; a reach between them is interpolated."""

DIRECTIONS = 48
"""int: Lines of sight the dust is sampled along (Fibonacci spread over the sphere)."""


def _best_correction():
    """The largest bolometric correction over the temperatures a star can have: the star of a given luminosity that shows brightest."""
    return float(np.max(photo.bolometric_correction(np.geomspace(1500.0, 300000.0, 2000))))


BEST_CORRECTION = _best_correction()
"""float: `max BC_V`, magnitudes: the colour of the star that gives the most V light for its luminosity."""

Tier = collections.namedtuple("Tier", "ledger key low high reach_pc cells visits")
"""One tier: `ledger` is "luminosity" or "mass", `key` the band number (luminosity) or the ring's low mass, `low` and `high` its limits (Lsun, or solar
masses; `high` may be `None`), `reach_pc` the radius of the sphere it is drawn over, `cells` the cells in that sphere and `visits` what drawing it takes
(`block_visits` for a luminosity band, a visit to every cell for a ring, whose band is on the cell path)."""

SkyTiers = collections.namedtuple("SkyTiers", "limit_mag observer_pc edge_pc budget luminosity dropped rings")
"""The table for one observer: the luminosity bands to draw, the bands the budget cut off, the mass rings, and the inputs."""


def cells_in_sphere(radius_pc, edge_pc):
    """The cells a sphere holds: its volume over the sector edge cubed (the cells tile space at that mean volume; 14,137 against 14,134 listed in a 60 pc
    sphere at a 4 pc edge)."""
    return 4.0 / 3.0 * math.pi * radius_pc ** 3 / edge_pc ** 3


def block_visits(radius_pc, edge_pc):
    """
    The visits the sphere draw makes over a sphere of `radius_pc`: one to each row it crosses and one to each block of `BLOCK_WEDGES` cells. The rows are
    the sphere's cross section in edges, pi times the radius in edges squared (counted in `spherical.sphere_rows` for spheres of 60 to 1,000 pc: 701 to
    196,345 rows, within 1 percent), and the blocks are the sphere's cells over `BLOCK_WEDGES`, since a row's window is long once the sphere is wide.
    """
    return math.pi * (radius_pc / edge_pc) ** 2 + cells_in_sphere(radius_pc, edge_pc) / bright_stars.BLOCK_WEDGES


def visits(tiers):
    """The visits of the whole table: the luminosity bands' block visits and the rings' cell visits."""
    return sum(tier.visits for tier in tiers.luminosity) + sum(tier.visits for tier in tiers.rings)


def floor_luminosity_sol(distance_pc, limit_mag, a_v=0.0):
    """
    The dimmest luminosity, solar units, a star can have and still show at `limit_mag` from `distance_pc` through `a_v` magnitudes of dust, for the best
    colour: `L = 10 ** ((M_bol,sun - BC_best + 5 log10(d / 10) + A_V - limit) / 2.5)`. Arrays or numbers.

    Raises:
        ValueError: A distance that is not positive and finite.
    """
    d = np.asarray(distance_pc, dtype=float)
    if np.any(d <= 0) or not np.all(np.isfinite(d)):
        raise ValueError("a distance must be positive and finite")
    return 10.0 ** ((photometry.SUN_ABSOLUTE_BOLOMETRIC_MAG - BEST_CORRECTION + 5.0 * np.log10(d / 10.0) + a_v - limit_mag) / 2.5)


def directions(count=DIRECTIONS):
    """`count` unit vectors spread evenly over the sphere."""
    i = np.arange(count) + 0.5
    phi = math.pi * (1.0 + 5.0 ** 0.5) * i
    z = 1.0 - 2.0 * i / count
    s = np.sqrt(1.0 - z * z)
    return np.stack([s * np.cos(phi), s * np.sin(phi), z], axis=1)


Extinction = collections.namedtuple("Extinction", "observer_pc directions grid_pc a_v")
"""The dust seen from one observer: `a_v[i, j]` magnitudes along direction `i` to distance `grid_pc[j]`."""


def extinction_curve(observer_pc, *, count=DIRECTIONS, grid_pc=DISTANCE_GRID_PC, **dust):
    """
    The V extinction from `observer_pc` (galactocentric `(x, y, z)`, parsecs) along `count` directions to each distance of `grid_pc`. Keyword arguments go to
    `sky_photometry.galactic_extinction_av` (`av_per_kpc_at_sun=0` is a dust-free sky), to scale the dust for a galaxy of another size.
    """
    units = directions(count)
    grid = np.asarray(grid_pc, dtype=float)
    a_v = np.empty((len(units), len(grid)))
    for i, unit in enumerate(units):
        a_v[i] = photo.galactic_extinction_av(observer_pc, unit, grid, **dust)
    return Extinction(tuple(float(x) for x in observer_pc), units, grid, a_v)


def reach_pc(luminosity_sol, limit_mag, curve):
    """
    The farthest a star of `luminosity_sol` (best colour) shows at `limit_mag`: the largest over the lines of sight, so the sphere is long enough for the
    clearest one. Zero when it shows nowhere; the end of the grid when it shows to the end.
    """
    if not (math.isfinite(luminosity_sol) and luminosity_sol > 0):
        raise ValueError("a luminosity must be positive and finite")
    base = photometry.SUN_ABSOLUTE_BOLOMETRIC_MAG - 2.5 * math.log10(luminosity_sol) - BEST_CORRECTION
    apparent = base + 5.0 * np.log10(curve.grid_pc / 10.0)[None, :] + curve.a_v
    best = 0.0
    for row in apparent:
        if row[0] > limit_mag:
            continue
        reach = curve.grid_pc[-1] if row[-1] <= limit_mag else float(np.interp(limit_mag, row, curve.grid_pc))
        best = max(best, reach)
    return best


def luminosity_bands(floor_sol=RING_TOP_SOL, ceiling_sol=None):
    """The fixed luminosity bands `(k, low, high)` from `floor_sol` up to the scatter's own floor (the default `BRIGHT_STAR_MIN_LUMINOSITY_SOL`)."""
    top = float(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL if ceiling_sol is None else ceiling_sol)
    return bright_stars.canonical_bands(floor_sol, top)


NEAR_RING_MASSES_SOL = (0.1, 0.3, 0.5, 0.7, 0.85)
"""tuple: The mass edges under 1 solar mass the sky's rings add (stand-in SKT4), so the stars born lighter than the Sun that show near the observer are
drawn. The sampler takes 0.1 as the lightest mass (VIEW.19 section 5)."""

RING_MIN_YIELD = 0.01
"""float: A ring under 1 solar mass reaches as far as this share of its stars (the largest of any population) is still bright enough to show (stand-in
SKT5): beyond, the ring draws 100 stars or more for every one that shows, as the rings over 1 solar mass do (VIEW.19: 140)."""

RING_SAMPLES = 400
"""int: The stars sampled per population to find a ring's share of visible stars."""

LOCAL_MIX = {"young": 1.0, "intermediate": 1.0, "old": 1.0}
"""dict: The weight each population has in a ring's share of visible stars when the caller gives none (stand-in SKT6): the three disk populations alike,
the mix near the Sun, where the halo and the bulge hold a thousandth of the stars. A caller that knows the observer's place passes the densities there."""


def ring_edges(near_masses=NEAR_RING_MASSES_SOL):
    """The mass edges of the sky's mass backfill, ascending: the near rings under 1 solar mass, the default scatter's ring masses and the mass limit."""
    return bright_stars.mass_band_edges([float(tuning.PHENOMENON_MIN_MASS_SOLAR)], sorted({*tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL, *near_masses}))


@functools.lru_cache(maxsize=None)
def ring_luminosities(low, high, count):
    """
    `{population: sorted luminosities, Lsun}` of `count` living stars sampled per population from the ring born with `low <= mass < high` and dimmer than
    `RING_TOP_SOL` (the brighter ones belong to the luminosity ledger), on streams named by the ring alone, so the same ring gives the same sample.
    """
    out = {}
    for population in bright_stars.SCATTERED_POPULATIONS:
        rng = draw.Stream(f"sky-tiers:ring:{low:g}:{high:g}:{population}")
        stars = star_population.sample_mass_band_stars(count, low, high, population, rng, max_luminosity_sol=RING_TOP_SOL)
        out[population] = np.sort(np.array([star["luminosity_w"] / constants.SOLAR_LUMINOSITY for star in stars]))
    return out


def visible_share(luminosities, floor_sol, weights=None):
    """The share of stars at or over `floor_sol`, averaged over the populations of `luminosities` (`ring_luminosities`) with `weights` (`LOCAL_MIX` by default)."""
    weights = LOCAL_MIX if weights is None else weights
    total = sum(weights.get(population, 0.0) for population in luminosities)
    if not total > 0.0:
        raise ValueError("no population of the weights is among the ring's")
    return sum(weights.get(population, 0.0) * (len(values) - int(np.searchsorted(values, floor_sol, side="left"))) / len(values)
               for population, values in luminosities.items()) / total


def ring_reach_pc(low, high, limit_mag, curve, *, min_yield=None, count=None, weights=None):
    """
    How far a ring under 1 solar mass is drawn: the largest distance of the curve's grid at which `min_yield` of its stars (`visible_share`, the mix of `weights`) still show
    at `limit_mag` through the least dust of any direction; never more than a star of `RING_TOP_SOL` reaches. Zero when the ring shows nowhere.
    """
    min_yield = RING_MIN_YIELD if min_yield is None else min_yield
    luminosities = ring_luminosities(low, high, RING_SAMPLES if count is None else count)
    floors = floor_luminosity_sol(curve.grid_pc, limit_mag, np.min(curve.a_v, axis=0))
    shown = [j for j, floor in enumerate(floors) if visible_share(luminosities, float(floor), weights) >= min_yield]
    reach = float(curve.grid_pc[shown[-1]]) if shown else 0.0
    return min(reach, reach_pc(RING_TOP_SOL, limit_mag, curve))


def build(observer_pc, limit_mag, edge_pc, *, visit_budget=BLOCK_VISIT_BUDGET, curve=None, bands=None, edges=None, weights=None, **dust):
    """
    The sky tier table for an observer at `observer_pc` (galactocentric parsecs) and a limiting magnitude.

    Args:
        observer_pc (tuple): The observer's place.
        limit_mag (float): The limiting magnitude of the sky (7.5 by default in the product, lower on a phone).
        edge_pc (float): The sector edge.
        visit_budget (float): The visits (`block_visits`) the luminosity bands may take.
        curve (Extinction, optional): The dust from `extinction_curve`; made from `dust` when left out.
        bands (list, optional): `luminosity_bands()` by default.
        edges (list, optional): `ring_edges()` by default; leave out the masses under 1 solar mass (`ring_edges(())`) for the rings of the default scatter alone.
        weights (dict, optional): The populations' weights in the share of a ring's stars that show (`LOCAL_MIX` by default).
        **dust: To `extinction_curve`.

    Returns:
        SkyTiers: The luminosity bands taken from the 50 Lsun band up while their visits fit `visit_budget` (a band that does not fit is dropped with every
            brighter one), the mass rings of 1 solar mass and up, each over the sphere of a star of `RING_TOP_SOL`, and the rings under it, each over the
            sphere `ring_reach_pc` gives (a ring that shows nowhere is left out).
    """
    if curve is None:
        curve = extinction_curve(observer_pc, **dust)
    bands = luminosity_bands() if bands is None else bands
    edges = ring_edges() if edges is None else edges
    taken, dropped, used, over = [], [], 0.0, False
    for k, low, high in bands:
        top = float(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL)
        radius = reach_pc(min(high, top) if high is not None else top, limit_mag, curve)
        tier = Tier("luminosity", k, low, high, radius, cells_in_sphere(radius, edge_pc), block_visits(radius, edge_pc))
        if over or used + tier.visits > visit_budget:
            over = True
            dropped.append(tier)
        else:
            used += tier.visits
            taken.append(tier)
    rings = []
    for low, high in zip(edges[:-1], edges[1:]):
        radius = reach_pc(RING_TOP_SOL, limit_mag, curve) if low >= 1.0 else ring_reach_pc(low, high, limit_mag, curve, weights=weights)
        if radius > 0.0:
            cells = cells_in_sphere(radius, edge_pc)
            rings.append(Tier("mass", low, low, high, radius, cells, cells))
    return SkyTiers(limit_mag, tuple(observer_pc), edge_pc, visit_budget, tuple(taken), tuple(dropped), tuple(rings))


def floor_table(tiers, curve, distances_pc):
    """
    The floor against distance: for each distance a row `{distance_pc, floor_sol, ledger, band, built}`. `floor_sol` is the dimmest star that shows there
    through the least dust of any direction; `ledger` is "mass" under `RING_TOP_SOL` (a ring draws it) and "luminosity" over it; `band` is the fixed band
    that holds the floor (its number, or the low mass of the lightest ring that reaches the distance); `built` is whether the table draws it (a dropped band
    is not).
    """
    rows = []
    for distance in distances_pc:
        j = int(np.searchsorted(curve.grid_pc, distance))
        j = min(max(j, 0), len(curve.grid_pc) - 1)
        a_v = float(np.min(curve.a_v[:, j]))
        floor = float(floor_luminosity_sol(distance, tiers.limit_mag, a_v))
        if floor < RING_TOP_SOL:
            reaching = [ring for ring in tiers.rings if ring.reach_pc >= distance]
            rows.append({"distance_pc": float(distance), "floor_sol": floor, "ledger": "mass", "band": reaching[0].key if reaching else None,
                         "built": bool(reaching)})
            continue
        held = [tier for tier in tiers.luminosity + tiers.dropped if tier.low <= floor]
        band = held[-1] if held else None
        rows.append({"distance_pc": float(distance), "floor_sol": floor, "ledger": "luminosity", "band": None if band is None else band.key,
                     "built": band is not None and band in tiers.luminosity})
    return rows
