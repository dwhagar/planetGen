# planetgen/generation/sky_charting.py

"""
Charting builds the sky tiers (VIEW.25, System sky S7)
======================================================

The sky a system sees (`galaxy/sky.py`, the route of `api/system_sky.py`) is made of the stars around it that are bright enough to show, and only the stars the scatter and the
backfill have stored can be listed. The tier table (`generation/sky_tiers.py`) says how far out each luminosity band and each mass ring must be drawn for a limiting magnitude;
this module draws them, when a system is charted, with the star backfill, and records what it drew (`db/sky_tier_ledger.py`).

* **Cells, not systems.** The tiers are built for the center of the observer's 50 parsec cell (`galaxy/sky_band.py`, the cell the glow band is kept for), with every radius padded by
  the cell's half diagonal, so every system in the cell is covered by the one set of draws. A system charted in a cell that is done adds nothing.
* **Spheres, not sectors.** A tier is drawn over a sphere (`bright_stars.luminosity_sphere_batches`, `mass_sphere_batches`: the cells that hold a candidate, a few hundred at a time) and
  recorded as a sphere. No cell of it gets a `sector_stats` row; the backfill and the fill settle a sector's level from the ledger when they touch it (`db/sky_tier_ledger.py`).
* **Each cell as it stands.** A cell's draw is cut by what it already holds: its own backfill level (`sector_stats`), else the level of the scatter at its place
  (`galaxy/scatter_limits.py`), and by the spheres already drawn; a charted cell is left alone (its systems were drawn below its level), and so is a cell the scatter's levels leave
  nothing under. The luminosity bands are drawn first and the rings after, so a ring draws the stars born in its mass range that are dimmer than the cell's luminosity level
  by then, the order the sector backfill itself takes.
* **Stopped draws go on.** A sphere is drawn in batches, each committed with the count of batches done (`sky_tier_draws.batches_done`); the walk of a sphere is the same each time, so a
  stopped draw is found by its kind, band, center and radius and goes on from there. A sphere counts only when its last batch is in.

    center, tiers = plan(observer_pc, 7.5, edge_pc)
    chart(conn, observer_pc, 7.5)             # draws and records every tier the ledger does not hold; returns the summary
"""

import collections
import dataclasses
import math
import time

from planetgen import tuning
from planetgen.db import chart_ledger, contention, sky_tier_ledger, store
from planetgen.galaxy import sky_band
from planetgen.galaxy.density import relative_density
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.generation import bright_stars, scatter_companions, sky_tiers
from planetgen.util import detmath

CELL_PAD_PC = sky_band.CELL_PC * math.sqrt(3.0) / 2.0
"""float: What every radius is padded by: the half diagonal of the observer's cell, so a system anywhere in it is inside the sphere (stand-in SKT7 in 6.176)."""

BATCH_COMMIT_STARS = 20000
"""int: A batch that has written this many stars is committed at once (a batch is a few hundred cells; this only bounds a transaction)."""


@dataclasses.dataclass(frozen=True)
class Drawn:
    """What one tier cost: the tier, its padded radius, the batches walked, the cells drawn in, the stars written, the seconds, and whether it was already in the ledger."""

    tier: object
    radius_pc: float
    batches: int
    cells: int
    stars: int
    seconds: float
    skipped: bool = False


MIN_SKY_LIMIT = 2.0
"""float: The brightest limiting magnitude a sky is built for; a setting under it builds none (`queue.sky_limit` 0 is off)."""


def limit_of(args):
    """
    The limiting magnitude a run's charting builds the sky tiers for, or `None` for none: the command's `--sky-limit` when given (0 is none), else the `queue.sky_limit` setting
    (7.5 unless the deployment changed it; 0 is none).
    """
    from planetgen.util.settings import get_settings
    given = getattr(args, "sky_limit", None)
    limit = get_settings().queue.sky_limit if given is None else float(given)
    return float(limit) if limit >= MIN_SKY_LIMIT else None


def observer_cell(observer_pc):
    """The 50 parsec cell of an observer and its center, galactocentric parsecs."""
    cell = sky_band.cell_of(observer_pc)
    return cell, sky_band.cell_center(cell)


RING_COST_ROWS = 2000
"""int: About how many rows of a ring's sphere are looked at to cost it (`cost_rings`); the sphere of a ring holds some hundred thousand rows at a 7.5 sky."""


def walk_edges_and_ceiling(limits):
    """
    What a ring's walk is made of, from the galaxy's `ScatterLimits`: the mass edges of the sector backfill (1, 2, 5, 8 and the galaxy's mass limits) with the sky's near rings under
    them, and the luminosity ceiling a ring leaves to the bright-star scatter (`_Context` draws with both; `cost_rings` costs with both).
    """
    found = list(limits.tiers or ())
    edges = bright_stars.mass_band_edges([tier.mass_limit_sol for tier in found],
                                         sorted({*sky_tiers.NEAR_RING_MASSES_SOL, *tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL}))
    return edges, max([tier.min_luminosity_sol for tier in found] or [0.0])


def cost_rings(center, tiers, skeleton, limits, rows=RING_COST_ROWS):
    """
    `tiers` with each mass ring costed by what its walk does (GEN.1485): `Tier.visits` of a ring is the visits `bright_stars.mass_sphere_walk_visits` counts over its drawn sphere (a
    visit to every block of a row whose bands are on the candidate path, to every cell of a row with a band on the cell path), looked at on about `rows` rows of the sphere. Before, every
    ring was costed one visit to each cell of its sphere (`sky_tiers.build` still does, with no galaxy to ask), which overstates a sparse ring by the cells a block holds.

    Args:
        skeleton (GalaxyShape): The stored galaxy.
        limits (ScatterLimits): The galaxy's limits, which say the ring's edges and ceiling (`walk_edges_and_ceiling`).
    """
    if not tiers.rings:
        return tiers
    edges, ceiling = walk_edges_and_ceiling(limits)
    costed = []
    for ring in tiers.rings:
        radius = padded_radius(ring)
        stride = max(1, int(math.ceil(math.pi * (radius / tiers.edge_pc) ** 2 / rows)))
        visits = bright_stars.mass_sphere_walk_visits(skeleton.shape, center, radius, tiers.edge_pc, skeleton.expected_system_count_at_density_1,
                                                      float(ring.low), float(ring.high), ceiling, edges, stride)
        costed.append(ring._replace(visits=visits))
    return tiers._replace(rings=tuple(costed))


def plan(observer_pc, limit_mag, edge_pc, *, limits=None, skeleton=None, **options):
    """
    The tier table for the cell of `observer_pc` (`sky_tiers.build` at the cell's center) and that center. With `limits` (the galaxy's `ScatterLimits`) the top ring ends at the
    galaxy's own mass limit instead of the default one, so no mass between the ring and the scatter is left out. With `skeleton` as well, each ring is costed by what its walk does
    (`cost_rings`) and not by the cells of its sphere.

    Returns:
        tuple: `(center, SkyTiers)`; the radii are not padded (see `padded_radius`).
    """
    _cell, center = observer_cell(observer_pc)
    tiers = sky_tiers.build(center, limit_mag, edge_pc, **options)
    if limits is not None and tiers.rings:
        top = [found.mass_limit_sol for found in (limits.tiers or ()) if found.mass_limit_sol is not None]
        if top:
            rings = list(tiers.rings)
            rings[-1] = rings[-1]._replace(high=float(max(top)))
            tiers = tiers._replace(rings=tuple(rings))
    if skeleton is not None and limits is not None:
        tiers = cost_rings(center, tiers, skeleton, limits)
    return center, tiers


def padded_radius(tier):
    """The radius a tier is drawn to: its reach plus `CELL_PAD_PC`."""
    return float(tier.reach_pc) + CELL_PAD_PC


class _Context:
    """What a draw reads once: the connection, the stored galaxy, the scatter's limits and seed, the mass edges and the widest caps the walk must allow for."""

    def __init__(self, conn):
        from planetgen.generation import run_plan
        self.conn = conn
        self.skeleton = store.get_galaxy_shape(conn)
        if self.skeleton is None:
            raise ValueError("the galaxy has no skeleton: nothing to chart a sky in")
        self.limits, self.seed = run_plan._scatter_limits_and_seed(conn, self.skeleton)
        self.edge_pc = self.skeleton.edge_pc
        self.expected = self.skeleton.expected_system_count_at_density_1
        tiers = list(self.limits.tiers or ())
        caps = [found.mass_limit_sol for found in tiers]
        self.walk_mass_range = None if (not caps or any(cap is None for cap in caps)) else (None, float(max(caps)))
        # The edges of the sector backfill (1, 2, 5, 8 and the galaxy's mass limits) and the sky's near rings under them: the bands from 1 solar mass up are the same bands, so
        # a star is the same star whichever backfill draws it.
        self.edges, self.walk_ceiling = walk_edges_and_ceiling(self.limits)
        self.blocks = {}

    def own_level(self, address, locked):
        """The luminosity level to store with a cell's mass (`store.set_sector_bright_masses`): the scatter's, when the cell had no level of its own, so it counts as having its own state from then on."""
        galaxy = self.limits.level_at(address[0], address[1])
        stored = locked.get(address)
        return galaxy if galaxy is not None and (stored is None or stored < 0) else None


def cell_states(context, ledger, addresses, locked, masses):
    """
    What each of `addresses` holds, as `{address: (luminosity level, held mass)}` for the cells a tier may draw in. `locked` and `masses` are the cells' stored rows
    (`store.lock_sector_stats`, `store.sector_bright_masses`). Left out: a charted cell, a filled one, and a cell with no scatter level (nothing was placed above anything). The level is the
    cell's own backfill level, else the scatter's at its place, then taken lower by the spheres drawn; the mass is the cell's backfill mass, else the scatter's mass limit at its place
    (`None`: no mass pass), then taken lower by the rings drawn.
    """
    conn, edge_pc = context.conn, context.edge_pc
    charted = chart_ledger.charted_among(conn, addresses)
    lum_start, mass_start = {}, {}
    for address in addresses:
        own = locked.get(address)
        if address in charted or own == store.FILLED_LEVEL:
            continue
        galaxy = context.limits.level_at(address[0], address[1])
        level = own if own is not None and own > 0 else galaxy
        if level is None:
            continue
        lum_start[address] = float(level)
        limit = context.limits.mass_limit_at(address[0], address[1])
        held = masses.get(address)
        held = limit if held is None else held if limit is None else min(held, limit)
        if held is not None:
            mass_start[address] = float(held)
    lum = ledger.luminosity_reached(list(lum_start), edge_pc, lum_start)
    mass = ledger.mass_reached(list(mass_start), edge_pc, mass_start)
    return {address: (lum.get(address, level), mass.get(address, mass_start.get(address))) for address, level in lum_start.items()}


def _luminosity_rows(context, tier, addresses, work, states):
    """
    The rows of one batch of a luminosity tier: each cell drawn from its tier's floor to the lower of the tier's top and the cell's level, with the mass cap its place and rings put on it.

    Returns:
        tuple: `(rows, levels)`; `levels` is `{address: tier floor}` for each cell drawn whose level the tier touched (the tier reaches up to it), to be stored in the cell's row.
    """
    groups = collections.defaultdict(list)
    levels = {}
    for address in addresses:
        state = states.get(address)
        if state is None:
            continue
        level, held = state
        top = level if tier.high is None else min(float(tier.high), level)
        if top <= tier.low + sky_tier_ledger.EDGE_TOLERANCE:
            continue
        groups[(top, held)].append(address)
        if tier.high is None or float(tier.high) >= level - sky_tier_ledger.EDGE_TOLERANCE:
            levels[address] = float(tier.low)
    rows = []
    for (top, held), group in sorted(groups.items(), key=lambda item: (item[0][0], -1.0 if item[0][1] is None else item[0][1])):
        rows.extend(bright_stars.backfill_cells(context.skeleton.shape, group, context.edge_pc, context.expected, float(tier.low), top, context.seed,
                                                mass_range=None if held is None else (None, held), memo=work))
    return rows, levels


def _mass_rows(context, tier, addresses, work, states):
    """
    The rows of one batch of a mass ring: each cell drawn from the ring's low mass to the lower of its high mass and the mass the cell is already held to, dimmer than its luminosity level.

    Returns:
        tuple: `(rows, masses)`; `masses` is `{address: tier low}` for each cell drawn that the ring reaches up to (its held mass is at or over the ring's high), to be stored in the cell's row.
    """
    groups = collections.defaultdict(list)
    masses = {}
    for address in addresses:
        state = states.get(address)
        if state is None:
            continue
        level, held = state
        top = float(tier.high) if held is None else min(float(tier.high), held)
        if top <= tier.low + sky_tier_ledger.EDGE_TOLERANCE:
            continue
        groups[(top, level)].append(address)
        if held is None or held >= float(tier.high) - sky_tier_ledger.EDGE_TOLERANCE:
            masses[address] = float(tier.low)
    rows = []
    for (top, level), group in sorted(groups.items()):
        rows.extend(bright_stars.backfill_mass_cells(context.skeleton.shape, group, context.edge_pc, context.expected, float(tier.low), top, level, context.seed,
                                                     context.edges, memo=work))
    return rows, masses


class DrawLock:
    """
    One sky draw at a time in the galaxy (a MySQL named lock held on the connection): two spheres that overlap would each draw a cell the other has not recorded yet. The sector
    backfill and the fill are not held back; they meet the draw at the cells' row locks.
    """

    NAME = "planetgen-sky-tier-draw"
    WAIT_SECONDS = 6 * 3600

    def __init__(self, conn):
        self.conn = conn

    def __enter__(self):
        got = self.conn.execute("SELECT GET_LOCK(?, ?) AS got", (self.NAME, self.WAIT_SECONDS)).fetchone()["got"]
        if got != 1:
            raise RuntimeError("another sky draw has held the lock for hours")
        return self

    def __exit__(self, *exc):
        self.conn.execute("SELECT RELEASE_LOCK(?)", (self.NAME,))


def draw_tier(context, tier, center, *, generator_epoch=None, on_batch=None, ledger=None):
    """
    Draws one tier over the sphere of its padded reach around `center`, in batches, each committed with the ledger's count, and records the sphere when the last is in. A tier the ledger
    already holds a sphere for is not drawn (`Drawn.skipped`).

    Args:
        context (_Context): `_Context(conn)`.
        tier (sky_tiers.Tier): A luminosity band or mass ring of `plan`.
        center (tuple): The cell's center, galactocentric parsecs.
        generator_epoch (int, optional): Recorded with the sphere.
        on_batch (callable, optional): Called as `on_batch(batches, cells, stars)` after each batch commits (a progress tick).
        ledger (Ledger, optional): The ledger as it stands; read when left out.

    Returns:
        Drawn: What it cost.
    """
    conn = context.conn
    radius = padded_radius(tier)
    ledger = sky_tier_ledger.read(conn) if ledger is None else ledger
    if ledger.drawn(tier.ledger, float(tier.low), tier.high, center, radius):
        return Drawn(tier, radius, 0, 0, 0, 0.0, skipped=True)
    started = time.perf_counter()
    draw = sky_tier_ledger.Draw(tier.ledger, float(tier.low), None if tier.high is None else float(tier.high), tuple(center), radius, 0.0)
    draw_id, done = sky_tier_ledger.begin(conn, draw, generator_epoch)
    conn.commit()
    if tier.ledger == sky_tier_ledger.LUMINOSITY:
        batches = bright_stars.luminosity_sphere_batches(context.skeleton.shape, center, radius, context.edge_pc, context.expected, float(tier.low),
                                                         None if tier.high is None else float(tier.high), context.seed, mass_range=context.walk_mass_range,
                                                         memo=context.blocks)
        make_rows = _luminosity_rows
    else:
        batches = bright_stars.mass_sphere_batches(context.skeleton.shape, center, radius, context.edge_pc, context.expected, float(tier.low),
                                                   float(tier.high), context.walk_ceiling, context.seed, context.edges, memo=context.blocks)
        make_rows = _mass_rows
    row = conn.execute("SELECT stars, seconds FROM sky_tier_draws WHERE id = ?", (draw_id,)).fetchone()
    stars, spent, cells, walked = int(row["stars"]), float(row["seconds"]), 0, 0
    for addresses, work in batches:
        walked += 1
        if walked <= done:
            continue                                   # written by the draw this one goes on from
        # The cells that hold a star here are locked, as a sector backfill locks the cells it draws, so the two never draw one at once; each gets the depth its tier reached.

        def batch(c):
            """One batch as one transaction (DB.114): locks, draws, writes the stars with their IDs and the levels, advances the ledger and commits; it keeps nothing between tries."""
            contention.begin_read_committed(c)
            entries = []
            for address in addresses:
                density = max(relative_density(sector_position_pc(*address, context.edge_pc), context.skeleton.shape), 0.0)
                entries.append((address, density, density * context.expected))
            locked = store.lock_sector_stats(c, entries)
            masses = store.sector_bright_masses(c, addresses, current=True)
            states = cell_states(context, ledger, addresses, locked, masses)
            rows, reached = make_rows(context, tier, addresses, work, states)
            scatter_companions.insert_bright_stars(c, rows)
            if tier.ledger == sky_tier_ledger.LUMINOSITY:
                store.set_sector_bright_levels(c, reached)
            else:
                store.set_sector_bright_masses(c, {address: (mass, context.own_level(address, locked)) for address, mass in reached.items()})
            sky_tier_ledger.advance(c, draw_id, walked, stars + len(rows), spent + time.perf_counter() - started)
            c.commit()
            return rows, states

        rows, states = contention.run_with_retries(conn, batch, gate="scatter-ids", label="Sky charting batch")
        stars += len(rows)
        cells += len(states)
        if on_batch is not None:
            on_batch(walked, len(states), len(rows))
    seconds = spent + time.perf_counter() - started
    sky_tier_ledger.finish(conn, draw_id, stars, seconds)
    conn.commit()
    return Drawn(tier, radius, walked, cells, stars, seconds)


def built_to_pc(conn, observer_pc):
    """
    How far the stored stars are complete around `observer_pc` because the sky tiers were drawn (VIEW.42): the ledger's spheres that hold the observer, followed from the scatter's own luminosity
    level down through the luminosity bands to `sky_tiers.RING_TOP_SOL` and from the mass limit down through the mass rings to the lightest ring edge. When both chains reach that far every tier
    that shows a star to the observer is drawn, the largest of their spheres' margins (radius less the observer's distance from its center) is how far the sky is built, and the faint floors that
    the rings carry (under 50 Lsun) count, which the cells' luminosity levels alone cannot say. 0.0 when a chain stops short (a tier not drawn yet, or dropped by the visit budget), when no
    sphere holds the observer, or when the galaxy has no scatter level.
    """
    scatter = store.bright_star_scatter_settings(conn)
    ledger = sky_tier_ledger.read(conn)
    if scatter is None or not ledger:
        return 0.0
    observer = tuple(float(c) for c in observer_pc)
    holding = sky_tier_ledger.Ledger(d for d in ledger.draws if detmath.dist(d.center, observer) < d.radius_pc)
    if not holding:
        return 0.0
    if holding.sphere_reached(sky_tier_ledger.LUMINOSITY, observer, 0.0, float(scatter[0])) > sky_tiers.RING_TOP_SOL + sky_tier_ledger.EDGE_TOLERANCE:
        return 0.0
    mass_limit = store.bright_star_mass_limit(conn)
    if mass_limit is not None and holding.sphere_reached(sky_tier_ledger.MASS, observer, 0.0, mass_limit) > min(sky_tiers.ring_edges()) + sky_tier_ledger.EDGE_TOLERANCE:
        return 0.0
    return max(d.radius_pc - detmath.dist(d.center, observer) for d in holding.draws)


def chart(conn, observer_pc, limit_mag, *, generator_epoch=None, on_batch=None, tiers=None, center=None):
    """
    Draws the sky tiers of the cell of `observer_pc` for `limit_mag`: the luminosity bands from the top band down, then the mass rings from the top ring down, each as `draw_tier`, one draw at a time in the galaxy (`DrawLock`). Tiers the ledger holds are skipped. `tiers` (with
    its `center`) is a table already built by `plan`, so a caller that costed it first does not build it twice.

    Returns:
        dict: `center`, `cell`, `limit_mag`, `tiers` (a `Drawn` each, in the order drawn), `stars`, `cells`, `seconds`, `dropped` (the luminosity bands the visit budget cut off) and `floors` (`sky_tiers.floor_table` at `FLOOR_DISTANCES_PC`: the dimmest star the stored sky shows at each distance and whether the table draws it).
    """
    context = _Context(conn)
    if tiers is None:
        center, tiers = plan(observer_pc, limit_mag, context.edge_pc, limits=context.limits)
    started = time.perf_counter()
    drawn = []
    with DrawLock(conn):
        ledger = sky_tier_ledger.read(conn)
        # From the top down: a cell is inside every sphere of the bands above the one being drawn, so each band touches the level the one before it left.
        for tier in (*sorted(tiers.luminosity, key=lambda tier: -tier.low), *sorted(tiers.rings, key=lambda tier: -tier.low)):
            result = draw_tier(context, tier, center, generator_epoch=generator_epoch, on_batch=on_batch, ledger=ledger)
            drawn.append(result)
            if not result.skipped:
                ledger = sky_tier_ledger.read(conn)        # the next tier's cells see what this one drew
    return {"center": tuple(center), "cell": sky_band.cell_of(observer_pc), "limit_mag": float(limit_mag), "tiers": drawn,
            "stars": sum(d.stars for d in drawn), "cells": sum(d.cells for d in drawn), "seconds": time.perf_counter() - started,
            "dropped": [tier.key for tier in tiers.dropped],
            "floors": [] if tiers.curve is None else sky_tiers.floor_table(tiers, tiers.curve, sky_tiers.FLOOR_DISTANCES_PC)}
