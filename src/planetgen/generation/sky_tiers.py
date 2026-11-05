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
  that are dimmer than 50 Lsun, so every ring reaches as far as a star of 50 Lsun does. Stars born under 1 solar mass are drawn by neither (VIEW.41).

The bright far bands hold few stars in an enormous number of cells (a 7.5 sky is 8e9 cell visits uncut), so the table is cut by a budget of cell visits
(`CELL_VISIT_BUDGET`): bands are taken from the 50 Lsun band up while the visits fit, and the rest are listed as `dropped`. The stars of a dropped band
that the scatter did not pre-place are not stored; the glow band (VIEW.27) carries their light.

    curve = extinction_curve(observer_pc)                       # the dust along 48 directions, once per observer
    tiers = build(observer_pc, 7.5, edge_pc, curve=curve)
    tiers.luminosity, tiers.rings, tiers.dropped, tiers.cell_visits
    floor_table(tiers, curve, (10, 50, 250, 1000))              # the floor against distance
"""

import collections
import math

import numpy as np

from planetgen import tuning
from planetgen.generation import bright_stars
from planetgen.physics import photometry
from planetgen.physics import sky_photometry as photo

CELL_VISIT_BUDGET = 3.0e8
"""float: The cell visits the luminosity bands may take for one sky (stand-in SKT1). VIEW.19 measured 2 to 30 microseconds a visit, so this is 10 minutes
to 2.5 hours one worker, and at a 7.5 limit with no dust it reaches the bands under about 1,000 Lsun, half of the luminosity ledger's stars."""

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

Tier = collections.namedtuple("Tier", "ledger key low high reach_pc cells")
"""One tier: `ledger` is "luminosity" or "mass", `key` the band number (luminosity) or the ring's low mass, `low` and `high` its limits (Lsun, or solar
masses; `high` may be `None`), `reach_pc` the radius of the sphere it is drawn over and `cells` the cells in that sphere."""

SkyTiers = collections.namedtuple("SkyTiers", "limit_mag observer_pc edge_pc budget luminosity dropped rings")
"""The table for one observer: the luminosity bands to draw, the bands the budget cut off, the mass rings, and the inputs."""


def cells_in_sphere(radius_pc, edge_pc):
    """The cells a sphere holds: its volume over the sector edge cubed (the cells tile space at that mean volume; 14,137 against 14,134 listed in a 60 pc
    sphere at a 4 pc edge)."""
    return 4.0 / 3.0 * math.pi * radius_pc ** 3 / edge_pc ** 3


def cell_visits(tiers):
    """The cell visits of the whole table: the luminosity bands and the rings."""
    return sum(tier.cells for tier in tiers.luminosity) + sum(tier.cells for tier in tiers.rings)


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


def ring_edges():
    """The mass edges of the default scatter's mass backfill: the ring masses and the mass limit, ascending."""
    return bright_stars.mass_band_edges([float(tuning.PHENOMENON_MIN_MASS_SOLAR)])


def build(observer_pc, limit_mag, edge_pc, *, cell_budget=CELL_VISIT_BUDGET, curve=None, bands=None, edges=None, **dust):
    """
    The sky tier table for an observer at `observer_pc` (galactocentric parsecs) and a limiting magnitude.

    Args:
        observer_pc (tuple): The observer's place.
        limit_mag (float): The limiting magnitude of the sky (7.5 by default in the product, lower on a phone).
        edge_pc (float): The sector edge.
        cell_budget (float): The cell visits the luminosity bands may take.
        curve (Extinction, optional): The dust from `extinction_curve`; made from `dust` when left out.
        bands (list, optional): `luminosity_bands()` by default.
        edges (list, optional): `ring_edges()` by default.
        **dust: To `extinction_curve`.

    Returns:
        SkyTiers: The luminosity bands taken from the 50 Lsun band up while their visits fit `cell_budget` (a band that does not fit is dropped with every
            brighter one), the mass rings, each over the sphere of a star of `RING_TOP_SOL`.
    """
    if curve is None:
        curve = extinction_curve(observer_pc, **dust)
    bands = luminosity_bands() if bands is None else bands
    edges = ring_edges() if edges is None else edges
    taken, dropped, used, over = [], [], 0.0, False
    for k, low, high in bands:
        top = float(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL)
        radius = reach_pc(min(high, top) if high is not None else top, limit_mag, curve)
        tier = Tier("luminosity", k, low, high, radius, cells_in_sphere(radius, edge_pc))
        if over or used + tier.cells > cell_budget:
            over = True
            dropped.append(tier)
        else:
            used += tier.cells
            taken.append(tier)
    ring_radius = reach_pc(RING_TOP_SOL, limit_mag, curve)
    rings = tuple(Tier("mass", low, low, high, ring_radius, cells_in_sphere(ring_radius, edge_pc)) for low, high in zip(edges[:-1], edges[1:]))
    return SkyTiers(limit_mag, tuple(observer_pc), edge_pc, cell_budget, tuple(taken), tuple(dropped), rings)


def floor_table(tiers, curve, distances_pc):
    """
    The floor against distance: for each distance a row `{distance_pc, floor_sol, ledger, band, built}`. `floor_sol` is the dimmest star that shows there
    through the least dust of any direction; `ledger` is "mass" under `RING_TOP_SOL` (a ring draws it) and "luminosity" over it; `band` is the fixed band
    that holds the floor (its number, or the ring's low mass for the lowest ring); `built` is whether the table draws it (a dropped band is not).
    """
    rows = []
    for distance in distances_pc:
        j = int(np.searchsorted(curve.grid_pc, distance))
        j = min(max(j, 0), len(curve.grid_pc) - 1)
        a_v = float(np.min(curve.a_v[:, j]))
        floor = float(floor_luminosity_sol(distance, tiers.limit_mag, a_v))
        if floor < RING_TOP_SOL:
            rows.append({"distance_pc": float(distance), "floor_sol": floor, "ledger": "mass", "band": tiers.rings[0].key if tiers.rings else None, "built": True})
            continue
        held = [tier for tier in tiers.luminosity + tiers.dropped if tier.low <= floor]
        band = held[-1] if held else None
        rows.append({"distance_pc": float(distance), "floor_sol": floor, "ledger": "luminosity", "band": None if band is None else band.key,
                     "built": band is not None and band in tiers.luminosity})
    return rows
