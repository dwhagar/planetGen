# planetgen/generation/run_galaxy.py

"""
Run Galaxy
==========

The `galaxy` command: many sectors placed as one galaxy, each at its
ring/layer/slot address in the cylindrical grid, with its systems inside
the sector's own cell. Four modes: `--ring` (batch), `--ring --slot`
(one address), `--center-sector`/`--radius-pc` (local neighborhood),
or neither (random start), plus blocks, columns and shells.
`ensure_sector_generated` and `generate_sector_neighborhood` are the
same per-sector logic for the web, run when a visitor reaches a sector,
against the skeleton `run_plan` builds.
"""

import argparse
import copy
import math
import threading
import time

import pymysql

from planetgen.db import sector_paths, store
from planetgen.generation import bright_stars as brightStars
from planetgen.galaxy import seed as galaxySeed, version_check
from planetgen.physics import mathcheck
from planetgen.queue import redisqueue
from planetgen import tuning as program_constants
from planetgen.util import log
from planetgen.galaxy.density import relative_density
from planetgen.galaxy.drill import (
    DRILL_LEVELS, drill_block_sectors, drill_children, drill_slabs, format_drill_key,
)
from planetgen.galaxy.geometry import (
    SectorCell, enumerate_sectors_within_radius, galactic_radius_pc, provisional_sector_designation,
    ring_sector_count, sector_address_at, sector_position_pc,
)
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.generation import run_common
from planetgen.generation import run_plan
from planetgen.generation import run_population
from planetgen.generation import run_sector
from planetgen.generation import steps
from planetgen.util import draw


LARGE_RING_WARNING_THRESHOLD = 2000
"""int: `--ring I` requires `--limit` or `--yes` when ring `I` holds more
slots than this (see `ring_sector_count`) -- from ring 321, ~4,200 ly
out. Anything larger takes a real, unbounded amount of time and disk, so
it needs an explicit choice."""


BACKFILL_FROM_CHOICES = ("edge", "none")
"""tuple: `--backfill-from`'s choices (GEN.30, GEN.98)."""


def _format_address(address):
    """`(ring, layer, slot)` as the words the CLI prints."""
    ring_index, layer_index, slot_index = address
    return f"ring {ring_index} layer {layer_index} slot {slot_index}"


class _BatchDensity:
    """
    Checks each address against the galaxy's stored outline and resolves
    each sector's own `--density` from the skeleton's real position-based
    `relative_density`, for `galaxy` mode's batch/local-neighborhood/
    random-start generation -- the same mechanism `ensure_sector_generated`
    uses for a single lazily-generated sector, applied across a whole run
    instead of every sector sharing one flat CLI value.

    The outline check (`galaxySkeleton.GalaxyBounds.contains`) always
    applies: nothing is ever generated outside the galaxy, even with an
    explicit `--density`/`--num-systems`. Past it, an explicit flag is a
    deliberate, uniform override for the whole run (`resolve` returns
    `args` unchanged); otherwise `resolve` sets each sector's density
    from its position. A sparse sector inside the outline is generated
    like any other (GEN.76).

    The skeleton and outline are fetched once, on first use -- neither
    changes mid-run.
    """

    def __init__(self, config):
        self._config = config
        self._skeleton = None
        self._bounds = None

    def _load(self):
        if self._skeleton is None:
            conn = store.get_connection(self._config)
            try:
                self._skeleton = store.get_galaxy_shape(conn)
                self._bounds = store.get_galaxy_bounds(conn)
            finally:
                conn.close()
            if self._skeleton is None:
                raise RuntimeError(
                    "The galaxy's skeleton has never been built (no galaxy_shape row) -- run "
                    "'planetgen plan' first, so every sector can be checked against the galaxy's "
                    "bounds before it is generated."
                )

    @property
    def skeleton(self):
        """The stored `store.GalaxySkeletonInfo`."""
        self._load()
        return self._skeleton

    @property
    def bounds(self):
        """The stored outline, a `galaxySkeleton.GalaxyBounds`."""
        self._load()
        return self._bounds

    def resolve(self, args, address, position_pc, outside_ok=False):
        """
        Args:
            args (argparse.Namespace): The `galaxy` subcommand's own parsed
                (and validated) arguments.
            address (tuple): This sector's `(ring, layer, slot)`.
            position_pc (tuple): This sector's `(x, y, z)` center, parsecs.
            outside_ok (bool): Resolve an address outside the outline too
                (one the user named outright, GEN.81) instead of `None`.

        Returns:
            argparse.Namespace or None: `None` if this address is outside
                the galaxy's stored outline (and `outside_ok` is false). Otherwise `args` itself when a
                density/count was given explicitly, else a fresh copy with
                `.density` set to this position's own `relative_density`
                and `.num_systems` cleared. A sparse address is never
                skipped (GEN.76): the predicted density only sets how many
                systems the draw expects, and the draw always runs.
        """
        if not outside_ok and not self.bounds.contains(address[0], address[1]):
            return None
        if args.density is not None or args.num_systems is not None:
            return args

        resolved = copy.copy(args)
        resolved.density = relative_density(position_pc, self.skeleton.shape)
        resolved.num_systems = None
        return resolved


def _fill_context(args, address, position_pc):
    """The `brightStars.FillContext` for one galaxy sector: its population
    mix, its unfilled pre-placed bright stars down to its own level
    (`store.bright_star_fill_level`: its backfill's, GEN.44, else the galaxy
    scatter's), and its unbuilt scattered phenomena when the galaxy's
    phenomenon scatter ran (GEN.100). `None` without a stored skeleton
    (nothing to take the mix from)."""
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            return None
        phenomena = below_cut = None
        settings = store.phenomenon_scatter_settings(conn)
        if settings is not None:
            phenomena = store.phenomena_for_sector(conn, *address)
            seed, min_mass_solar = settings
            if min_mass_solar is not None:
                # GEN.168: what the scatter left below its cut, at the
                # density it used (the plan's, not this sector's own roll).
                expected = (max(relative_density(position_pc, skeleton.shape), 0.0)
                            * skeleton.expected_system_count_at_density_1)
                below_cut = brightStars.BelowCut(tuple(address), min_mass_solar, seed, expected)
        level = store.bright_star_fill_level(conn, *address)
        if level is None:
            return brightStars.FillContext(position_pc, skeleton.shape, phenomenon_rows=phenomena,
                                           below_cut=below_cut)
        rows = store.bright_stars_for_sector(conn, *address)
    finally:
        conn.close()
    return brightStars.FillContext(position_pc, skeleton.shape, rows, min_luminosity_sol=level,
                                   phenomenon_rows=phenomena, below_cut=below_cut)


def backfill_tiers(radius_ly=None, min_luminosity_sol=None, tiers=None):
    """
    The backfill's distance tiers (GEN.30), as `(out_to_ly,
    min_luminosity_sol)` pairs, nearest first: `tiers` when given, else
    one tier when either `radius_ly` or `min_luminosity_sol` is (a single
    floor out to a single distance, GEN.23's original form; the other
    defaults to the nearest tier's floor or the farthest tier's
    distance), else `program_constants.BRIGHT_STAR_BACKFILL_TIERS`.

    Returns:
        tuple: `(out_to_ly, min_luminosity_sol)` float pairs, sorted by
            distance.
    """
    default = program_constants.BRIGHT_STAR_BACKFILL_TIERS
    if tiers is None and (radius_ly is not None or min_luminosity_sol is not None):
        tiers = ((default[-1][0] if radius_ly is None else radius_ly,
                  default[0][1] if min_luminosity_sol is None else min_luminosity_sol),)
    return tuple(sorted((float(out_to), float(floor)) for out_to, floor in (tiers or default)))


def format_backfill_tiers(tiers):
    """`tiers` as text for the log ("10:100,25:250,...", ly:L_sun)."""
    return ",".join(f"{out_to:g}:{floor:g}" for out_to, floor in tiers)


def _tier_floor(tiers, distance_ly):
    """The floor of the first tier reaching `distance_ly`, or None past the last."""
    for out_to, floor in tiers:
        if distance_ly < out_to or (out_to, floor) == tiers[-1] and distance_ly <= out_to:
            return floor
    return None


def backfill_bright_stars(config, center_pc, radius_ly=None, min_luminosity_sol=None, tiers=None):
    """`backfill_bright_stars_around` one center (see there)."""
    return backfill_bright_stars_around(config, [center_pc], radius_ly, min_luminosity_sol, tiers)


def backfill_bright_stars_around(config, centers_pc, radius_ly=None, min_luminosity_sol=None, tiers=None,
                                 progress=None):
    """
    The bright-star backfill around generated sectors (GEN.23, tiered by
    GEN.30, per sector since GEN.44): every unfilled sector within the
    farthest tier of any of `centers_pc` gets every star from its tier's
    floor (by its own distance from the nearest center) up to the level it
    already holds: its own `sector_stats` level, else the galaxy scatter's
    threshold, else no ceiling when no scatter ran. A sector at level 0
    (generated) or already that deep is skipped (one query for the whole
    sphere), and one a nearer sector reaches later is topped up with only
    the band it lacks, so no star is ever drawn twice. A sector with no
    level anywhere that still holds stars is left over from a run that
    failed part way: they are wiped and it is drawn whole.

    Sectors are drawn a few hundred at a time, under their rows' locks,
    and each one's new level is written in the transaction that writes
    its stars. A sector's draw is seeded from the galaxy's scatter seed,
    its address and fixed luminosity bands (`brightStars.backfill_cells`),
    so it is the same whichever sector reached it first, and whatever
    steps took it there.

    Args:
        config (MySQLConfig): Connection parameters.
        centers_pc (list): Generated sectors' `(x, y, z)` centers,
            galaxy-frame parsecs.
        radius_ly (float, optional): With `min_luminosity_sol`, one tier
            instead of the defaults (`backfill_tiers`).
        min_luminosity_sol (float, optional): See `radius_ly`.
        tiers (tuple, optional): `(out_to_ly, min_luminosity_sol)` pairs;
            defaults to `program_constants.BRIGHT_STAR_BACKFILL_TIERS`
            (100 L_sun under 10 ly, 250 to 25 ly, 500 to 50 ly, 750 to
            100 ly).
        progress (Progress, optional): A live `_generation_progress`
            display: the backfill adds its own bar, counting the sectors
            it visits, with its ETA (PERF.28).

    Returns:
        dict: `sectors` (drawn now) and `stars` (placed now), both int;
            zeros without a stored skeleton or when the galaxy scatter
            already went that deep.
    """
    tiers = backfill_tiers(radius_ly, min_luminosity_sol, tiers)
    radius_pc = ly_to_pc(tiers[-1][0])
    summary = {"sectors": 0, "stars": 0}
    layers = {}
    # Finding the sectors to draw can take a while on a large radius; the
    # bar is there from the start (ADM.26), unmeasured until the count is known.
    bar = steps.Step("Bright-star backfill (finding sectors)", "backfill", None, progress=progress,
                     stats=run_common._stats_for_config(config)).__enter__()
    failed = True
    conn = store.get_connection(config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            return summary
        galaxy_level, seed = run_plan._scatter_level_and_seed(conn, skeleton)
        if galaxy_level is not None and galaxy_level <= min(floor for _out_to, floor in tiers):
            return summary
        bounds = store.get_galaxy_bounds(conn)
        floors = {}
        for center_pc in centers_pc:
            for ring, layer, slot, *_xyz, distance_pc in enumerate_sectors_within_radius(
                    center_pc, radius_pc, skeleton.edge_pc):
                if not bounds.contains(ring, layer):
                    continue
                floor = _tier_floor(tiers, pc_to_ly(distance_pc))
                if floor is not None:
                    floors[(ring, layer, slot)] = min(floor, floors.get((ring, layer, slot), math.inf))
        if galaxy_level is not None:
            floors = {address: floor for address, floor in floors.items() if floor < galaxy_level}
        levels = store.sector_bright_levels(conn, floors)
        filled = store.get_occupied_addresses(conn, {address[0] for address in floors})
        todo = sorted(address for address, floor in floors.items()
                      if address not in filled and run_plan._needs_band(levels.get(address), galaxy_level, floor))
        on_sector = None
        if todo:
            bar.update(description="Bright-star backfill (sectors)", total=len(todo))

            def on_sector():
                bar.update(advance=1)
        drawn = run_plan._draw_sector_bands(conn, skeleton, todo, floors, galaxy_level, seed, on_sector=on_sector,
                                            layers=layers)
        summary["sectors"], summary["stars"] = drawn
        failed = False
    except BaseException:
        conn.rollback()
        raise
    finally:
        conn.close()
        bar.close(success=not failed)     # nothing to draw: no bar left unmeasured
    run_plan.log_layers("Bright-star backfill", layers)
    if summary["sectors"]:
        log.debug(f"bright-star backfill: {summary['stars']} stars in {summary['sectors']} sector(s), tiers "
                  f"{format_backfill_tiers(tiers)}")
    return summary


def backfill_after_run(args, edge_pc, started_at):
    """
    The bright-star backfill once a `galaxy` run has generated every
    sector it was asked for (GEN.30, GEN.98): out from the run's edge
    (`--backfill-from edge`, the default), the backfill distance past the
    farthest generated sector in every direction (the union of the tier
    radii around every sector the run generated), or not at all (`none`).
    Nothing when the run generated no sector.

    Args:
        args (argparse.Namespace): The run's arguments.
        edge_pc (float): The sector edge, parsecs.
        started_at (datetime): The database's clock when the run started
            (`store.database_now`); sectors created since are the run's.

    Returns:
        dict: `backfill_bright_stars_around`'s summary (zeros when skipped).
    """
    mode = getattr(args, "backfill_from", "edge") or "edge"
    summary = {"sectors": 0, "stars": 0}
    if mode == "none":
        return summary
    config = store.mysql_config_from_args(args)
    conn = store.get_connection(config)
    try:
        generated = store.sector_centers_since(conn, started_at)
        if not generated:
            return summary
        centers = generated
    finally:
        conn.close()
    log.normal("Backfilling the bright stars out from the edge of the run...")
    # Its own bar and ETA (PERF.28): a backfill can take minutes.
    with run_common._generation_progress() as progress:
        log.set_console(progress.console)
        try:
            summary = backfill_bright_stars_around(config, centers, progress=progress)
        finally:
            log.reset_console()
    log.normal(f"Backfilled {summary['stars']:,} bright stars in {summary['sectors']:,} sectors.")
    return summary


def link_after_run(args, started_at):
    """
    PERF.45: links the sectors a `galaxy` run created to their neighbours
    (`store.link_sector_neighbors`): containment, nearest systems and
    quadrants, and the new systems merged into the lists of the sectors
    around. Done once for the whole run, in place of one locked step per
    saved sector, so the workers no longer queue for it; the stored result is
    the same. Runs even when the run is cancelled or fails, for the sectors it
    saved; a failure here only warns (the correlative update links
    everything again).

    Args:
        args (argparse.Namespace): The run's arguments.
        started_at (datetime): The database's clock when the run started.

    Returns:
        int: How many sectors were linked.
    """
    if not getattr(args, "link_later", False):
        return 0
    config = store.mysql_config_from_args(args)
    try:
        conn = store.get_connection(config)
        try:
            created = sector_paths.sector_ids_since(conn, started_at)
        finally:
            conn.close()
        if not created:
            return 0
        log.normal("Linking the new sectors to their neighbours...")
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step("Neighbours", "link", len(created) * store.LINK_PHASES, args=args,
                                progress=progress) as bar:
                    def on_progress(done, total, step=""):
                        # Never back: a batch retried after a deadlock repeats its steps.
                        bar.update(completed=max(done, bar.done), total=total,
                                   description="Neighbours" + (f" ({step})" if step else ""))

                    linked = store.link_sector_neighbors(config, created, on_progress)
            finally:
                log.reset_console()
    except Exception as exc:  # noqa: BLE001 -- the sectors are saved; the links can be redone
        log.normal(f"Warning: linking the new sectors to their neighbours failed: {exc}")
        return 0
    log.normal(f"Linked {linked:,} sectors to their neighbours.")
    return linked


def settle_after_run(args, started_at):
    """
    The last step of a `galaxy` run (GEN.126): once every sector, bright star and population is
    final, saves the path of every body in the sectors the run created and in the sectors around
    them (`sector_paths.settle_sectors`), so each path sees the final set of neighbours whatever
    order the workers filled sectors in. Skipped with `--no-settle` (and for an estimate).

    Args:
        args (argparse.Namespace): The run's arguments.
        started_at (datetime): The database's clock when the run started.

    Returns:
        int: How many paths were saved (0 when skipped or nothing was generated).
    """
    if getattr(args, "no_settle", False):
        return 0
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        created = sector_paths.sector_ids_since(conn, started_at)
        if not created:
            return 0
        log.normal("Saving the sector paths...")
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step("Sector paths", "paths", len(created), args=args, progress=progress) as bar:
                    def on_progress(done, total):
                        bar.update(completed=done, total=total)

                    saved = sector_paths.settle_sectors(conn, created, on_progress)
            finally:
                log.reset_console()
    finally:
        conn.close()
    log.normal(f"Saved {saved:,} sector paths.")
    return saved


def generate_and_save_sector_at(args, address, position_pc, edge_pc, channel=None):
    """
    Generates one sector via `generate_sector` inside its real grid cell
    and saves it -- the per-sector unit of work every `galaxy` mode
    repeats. The bright-star backfill runs once the whole run is done
    (`backfill_after_run`, GEN.30), not per sector.

    Args:
        args (argparse.Namespace): Parsed arguments (see `add_galaxy_arguments`).
        address (tuple): This sector's `(ring, layer, slot)`.
        position_pc (tuple): Its `(x, y, z)` center, parsecs.
        edge_pc (float): The sector edge length, parsecs.
        channel (optional): Where a worker reports the save's progress
            (`steps.worker_step`, PERF.50); `None` reports nowhere.

    Returns:
        tuple: `(sector_id, sector_name, sector)` of the newly saved sector
            -- `sector_name` is read back after saving, since name
            uniqueness (v22) may rename it on save.
    """
    ring_index, layer_index, slot_index = address
    x, y, z = position_pc
    radius_pc = galactic_radius_pc(position_pc)
    cell = SectorCell.for_ring(ring_index, pc_to_ly(edge_pc))

    fill = _fill_context(args, address, position_pc)
    galaxy_position = {
        "center_x_pc": x, "center_y_pc": y, "center_z_pc": z,
        "galactic_radius_pc": radius_pc,
        "ring_index": ring_index, "layer_index": layer_index, "ring_slot_index": slot_index,
    }
    galaxy_seed = _galaxy_seed(args)
    # GEN.47: the galaxy's molecular clouds that reach this sector.
    cloud_field = None
    if fill is not None:
        cloud_field = (galaxy_seed, fill.shape, tuple(position_pc), edge_pc * math.sqrt(3) / 2)
    # GEN.39: every draw for this sector, its save included, comes from
    # its own seed, so the run's order and worker count don't matter.
    with galaxySeed.seeded(galaxy_seed, "sector", address):
        _sector_name, sector = run_sector.generate_directed_sector(
            args, galaxy_seed=galaxy_seed, address=address,
            galactic_center_dist_ly=pc_to_ly(radius_pc), cell=cell, fill=fill, cloud_field=cloud_field)
        if address == run_sector.NUCLEUS_ADDRESS and not (fill is not None and fill.phenomena_scattered):
            run_sector.add_galactic_nucleus(sector, args, pc_to_ly(radius_pc))
        sector.place_in_galaxy(tuple(pc_to_ly(c) for c in position_pc))
        # PERF.45: a run links its sectors to their neighbours once, at its end
        # (`link_after_run`), instead of one at a time under a lock.
        units = len(sector.entries) + len(sector.phenomena)
        with steps.worker_step(channel, "save:" + ",".join(str(part) for part in address),
                               f"Saving sector {provisional_sector_designation(*address)}", "save", units,
                               workers=run_common._worker_count(args)) as save_step:
            sector_id = store.save_sector(sector, config=store.mysql_config_from_args(args),
                                        galaxy_position=galaxy_position,
                                        link_neighbors=not getattr(args, "link_later", False),
                                        progress=save_step)
    run_common._count_sector(sector)
    return sector_id, sector.name, sector


def _galaxy_seed(args):
    """The galaxy's stored 16-byte seed (GEN.39), or `None` when it has
    none (never planned, or planned before schema v51)."""
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        return store.get_galaxy_seed(conn)
    finally:
        conn.close()


def _default_generation_args(config=None):
    """
    Builds a minimal `argparse.Namespace`, shaped exactly like the
    `sector` subcommand's own parsed output, for a caller with no actual
    command line to parse -- `ensure_sector_generated` in particular,
    which runs at sector-visit time.

    Built by feeding an empty argument list through
    `add_shared_generation_options`/`validate_shared_generation_args`
    (rather than hand-listing every field here) so this can never
    silently drift from whatever options those functions actually
    declare.

    Args:
        config (MySQLConfig, optional): Connection parameters, copied
            onto the returned namespace's `mysql_*` attributes. Defaults
            to `DEFAULT_MYSQL_CONFIG`.

    Returns:
        argparse.Namespace: Every shared generation option at its
            documented default (`num_systems` resolved to `10`) -- a caller
            driving generation by density should set `args.density` and
            clear `args.num_systems` back to `None` first.
    """
    # The command line imports this module, so its options are imported
    # here, when a visit first needs them.
    from planetgen.cli.generate import add_shared_generation_options, validate_shared_generation_args

    parser = argparse.ArgumentParser(prefix_chars='-+')
    add_shared_generation_options(parser)
    args = parser.parse_args([])
    validate_shared_generation_args(args, parser)

    args.sector_name = None
    args.system_file = None
    args.num_orbits = None
    args.name = None

    config = config or store.MySQLConfig()
    args.mysql_host = config.host
    args.mysql_port = config.port
    args.mysql_user = config.user
    args.mysql_password = config.password
    args.mysql_database = config.database
    return args


def queue_settle(config, sector_ids):
    """
    Saves the sector paths of `sector_ids` and the sectors around them
    (GEN.126) as a queued job, so a caller that generated a sector on the
    spot doesn't wait for up to 125 sectors of integration. With no Redis
    server to queue it on, it runs here instead.

    Returns:
        str | None: The job's id, or `None` when it ran inline.
    """
    from planetgen.queue import api_jobs
    sector_ids = sorted(sector_ids)
    try:
        return api_jobs.submit(api_jobs.settle_sectors, sector_ids, config)
    except api_jobs.NoQueue:
        api_jobs.settle_sectors(sector_ids, config)
        return None


def ensure_sector_generated(ring_index, layer_index, ring_slot_index, config=None, backfill=True,
                            outside_ok=False, settle=True, link=True):
    """
    The galaxy map's "recalculate on visit" entry point: returns the
    sector already generated at this address if one exists; otherwise
    uses the stored skeleton's outline (see section 4, `build_skeleton`)
    to check the address is inside the galaxy, and if so generates and
    saves it on the spot, using its own position's `relative_density` as
    the `--density` multiplier. A sparse address generates too (GEN.76):
    its draw may come out with no systems, and the sector is still saved
    and marked generated.

    A concurrent visit to the same never-generated address is handled by
    `sectors`'s `UNIQUE (ring_index, layer_index, ring_slot_index)`: the
    losing `INSERT` raises `pymysql.err.IntegrityError`, caught here and
    turned into "return what the other call just created".

    A new sector's bright-star backfill runs here too, unless `backfill`
    is false: `planetgen galaxy --slot` leaves it to the end of the run
    (`backfill_after_run`), where it has its own progress bar (PERF.28).

    An address outside the galaxy's stored outline is generated only with
    `outside_ok` (`planetgen galaxy --slot` naming it outright, GEN.81),
    at the halo floor's density.

    A new sector's paths (and its neighbours') are saved by a queued job
    unless `settle` is false (GEN.126): the caller that settles at the end
    of its own run, or after its own edit, passes false. Likewise the sector is
    linked to its neighbours (containment, nearest systems) as it is saved
    unless `link` is false (PERF.45): the run that links its sectors at its
    end passes false.

    Returns:
        dict: `created` (bool), `qualifies` (bool: inside the galaxy's
              stored outline; a sector made with `outside_ok` is not), `sector_id` (int or `None`), `sector_name`
              (str or `None`, only when `created`), and `summary`
              (`sector_generation_summary_lines`, only when `created`).

    Raises:
        ValueError: If `ring_slot_index` is out of range for the ring.
        RuntimeError: If the galaxy's skeleton has never been built --
                     run `planetgen plan` first.
    """
    address = (ring_index, layer_index, ring_slot_index)
    conn = store.get_connection(config)
    try:
        existing_id = store.get_sector_id_at(conn, *address)
        if existing_id is not None:
            return {"created": False, "qualifies": True, "sector_id": existing_id, "sector_name": None}

        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError(
                "The galaxy's skeleton has never been built (no galaxy_shape row) -- run "
                "'planetgen plan' first."
            )

        bounds = store.get_galaxy_bounds(conn)
    finally:
        conn.close()

    position_pc = sector_position_pc(ring_index, layer_index, ring_slot_index, skeleton.edge_pc)
    qualifies = bounds.contains(ring_index, layer_index)
    if not qualifies and not outside_ok:
        # Past the layer's stored outer ring -- galaxySkeleton's bound is
        # exact, so this is a certain "no".
        return {"created": False, "qualifies": False, "sector_id": None, "sector_name": None}

    # Every address inside the outline generates, however sparse
    # (GEN.76): the density sets the draw's expected count, never whether
    # the draw happens.
    args = _default_generation_args(config=config)
    args.density = relative_density(position_pc, skeleton.shape)
    args.num_systems = None
    args.link_later = not link

    try:
        sector_id, sector_name, sector = generate_and_save_sector_at(
            args, address, position_pc, skeleton.edge_pc,
        )
    except pymysql.err.IntegrityError:
        conn = store.get_connection(config)
        try:
            existing_id = store.get_sector_id_at(conn, *address)
        finally:
            conn.close()
        if existing_id is None:
            raise
        return {"created": False, "qualifies": True, "sector_id": existing_id, "sector_name": None}

    if backfill:
        backfill_bright_stars(config, position_pc)  # GEN.30: around the requested sector
    if settle:
        queue_settle(config, [sector_id])
    return {"created": True, "qualifies": qualifies, "sector_id": sector_id, "sector_name": sector_name,
            "summary": run_sector.sector_generation_summary_lines(sector, args)}


def _log_saved(saved, address, suffix=""):
    designation = provisional_sector_designation(*address)
    log.normal(
        f"Saved sector '{saved['name']}' [{designation}] at {_format_address(address)}{suffix} "
        f"(sector_id={saved['sector_id']})."
    )
    run_sector._log_summary(saved["summary"])


def _submit_batch(args, batch, title, edge_pc, progress, bar):
    """Queues every `(address, position_pc, sector_args, suffix)` of
    `batch` (`_submit_sector`) and waits for them."""
    relay = steps.Relay(args, progress)
    stop = threading.Event()
    channel = drain = None
    try:
        with run_common._work_queue(args, title) as queue:
            if queue.parallel:
                channel = queue.channel("sector-progress")
                drain = threading.Thread(target=relay.drain, args=(channel, stop), name="sector-progress", daemon=True)
                drain.start()
            else:
                channel = relay   # one worker saves in this process: its reports go straight to the relay
            queue.expect(len(batch))
            log.normal(f"Generating {len(batch):,} sector(s) with {queue.workers} worker(s); each is reported below as "
                       f"it is saved, so the first report can take a while in a dense region.")
            for address, position_pc, sector_args, suffix in batch:
                _submit_sector(queue, sector_args, address, position_pc, edge_pc, bar, suffix=suffix, channel=channel)
    finally:
        stop.set()
        if drain is not None:
            drain.join(timeout=5)
        if isinstance(channel, redisqueue.Channel):
            channel.close()
        relay.close()


def _fill_sector_task(payload):
    """
    One galaxy sector, start to finish -- a work queue task (PERF.8): runs
    in a worker process (or in this one, with one worker), generates the
    sector in its grid cell and saves it in one transaction, and returns
    what the run reports for it.

    Returns:
        dict: `sector_id`, `name` (as saved), `systems`, `stars`,
            `density`, `phenomena` and `summary`
            (`sector_generation_summary_lines`).
    """
    sector_args = payload["args"]
    sector_id, sector_name, sector = generate_and_save_sector_at(
        sector_args, payload["address"], payload["position_pc"], payload["edge_pc"],
        channel=payload.get("channel"),
    )
    return {
        "sector_id": sector_id, "name": sector_name, "systems": len(sector.entries),
        "stars": run_sector.sector_star_count(sector), "density": run_common._sector_density(sector_args),
        "phenomena": len(sector.phenomena), "summary": run_sector.sector_generation_summary_lines(sector, sector_args),
    }


def _submit_sector(queue, sector_args, address, position_pc, edge_pc, bar, suffix="", channel=None):
    """Queues one galaxy sector (`_fill_sector_task`); when it's saved,
    advances `bar`, logs it and adds its time to the speed stats."""
    def saved(result, seconds, _weight):
        if queue.parallel:
            # A worker's own RUN_COUNTS die with it; the run's are here.
            run_common.RUN_COUNTS["sectors"] += 1
            run_common.RUN_COUNTS["systems"] += result["systems"]
            run_common.RUN_COUNTS["phenomena"] += result["phenomena"]
        run_common._record_sector(sector_args, result, seconds)
        bar.update(advance=1)
        _log_saved(result, address, suffix=suffix)

    payload = {"args": sector_args, "address": address, "position_pc": position_pc, "edge_pc": edge_pc,
               "channel": channel}
    queue.submit("sector", ",".join(str(part) for part in address), _fill_sector_task, payload, on_done=saved)


def _require_inside(args, bounds, ring_index, layer_index, what):
    """
    Warns (GEN.81), or under `--strict` stops the run before anything is
    generated, when `(ring_index, layer_index)` lies outside the galaxy's
    stored outline.

    Returns:
        bool: Whether it lies inside.
    """
    if bounds.contains(ring_index, layer_index):
        return True
    run_common._refuse_or_warn(args, f"{what} -- {bounds.describe_miss(ring_index, layer_index)}.")
    return False


def run_ring_batch(args, edge_pc, progress):
    """
    Batch mode: generates every not-yet-generated sector in ring
    `args.ring` at layer `args.layer` (up to `args.limit`, if given) --
    when neither `--density` nor `--num-systems` was given, each takes its
    own position's density (see `_BatchDensity.resolve`); sparse ones are
    generated too (GEN.76).

    Args:
        args (argparse.Namespace): Parsed arguments; `args.ring` set.
        edge_pc (float): The sector edge length, in parsecs (`_edge_pc`).
        progress (rich.progress.Progress): `run_galaxy`'s shared progress
            display -- an outer "Sectors" task is added here and advanced
            once per slot visited.

    Raises:
        SystemExit: Under `--strict`, if the ring and layer lie outside the
                   galaxy's stored outline, or the ring's slot count
                   exceeds `LARGE_RING_WARNING_THRESHOLD` and neither
                   `--limit` nor `--yes` was given. Without it each is a
                   warning and the ring is generated (GEN.81).
    """
    ring_index, layer_index = args.ring, args.layer
    total_slots = ring_sector_count(ring_index)

    if total_slots > LARGE_RING_WARNING_THRESHOLD and args.limit is None and not args.yes:
        run_common._refuse_or_warn(args, f"Ring {ring_index} holds {total_slots} sector slots, a very large run "
                              f"(--limit N generates only the first N not-yet-generated slots).")

    mysql_config = store.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    # Named outright, so generated even outside the outline (GEN.81).
    outside_ok = not _require_inside(args, batch_density.bounds, ring_index, layer_index,
                                     f"ring {ring_index} layer {layer_index}")

    conn = store.get_connection(mysql_config)
    try:
        occupied = {a for a in store.get_occupied_addresses(conn, [ring_index]) if a[1] == layer_index}
    finally:
        conn.close()

    what = f"ring {ring_index} layer {layer_index}"
    batch = []
    skipped = 0
    for slot_index in range(total_slots):
        if args.limit is not None and len(batch) >= args.limit:
            break
        address = (ring_index, layer_index, slot_index)
        if address in occupied:
            continue

        position_pc = sector_position_pc(ring_index, layer_index, slot_index, edge_pc)
        sector_args = batch_density.resolve(args, address, position_pc, outside_ok=outside_ok)
        if sector_args is None:
            log.debug(f"{_format_address(address)}: skipped (outside its layer's "
                      f"stored extent)")
            skipped += 1
            continue
        batch.append((address, position_pc, sector_args, ""))

    run_common._check_estimate(args, [item[2] for item in batch], what, progress)
    with run_common.sector_bar(progress, args, f"Sectors ({what})", len(batch)) as bar:
        _submit_batch(args, batch, f"Sectors ({what})", edge_pc, progress, bar)
    generated = len(batch)
    skip_note = f", {skipped} skipped (outside the outline)" if skipped else ""
    log.normal(
        f"Generated {generated} new sector(s) in ring {ring_index} layer {layer_index} "
        f"({total_slots} total slots, {len(occupied)} already existed{skip_note})."
    )


def _neighborhood_candidates(center, radius_pc, edge_pc, config, bounds):
    """
    Every address within `radius_pc` of `center` that lies inside the
    galaxy's outline, the set of those already occupied, and how many
    addresses in the sphere were left out for lying outside it -- the
    sphere is trimmed to the galaxy before anything is generated.
    """
    if bounds:
        # Nothing in the galaxy lies farther from `center` than this, so a
        # larger radius only enumerates (and discards) empty space --
        # `--radius-pc 1e300` used to walk ~1e300 rings before trimming.
        top_layer = max(abs(layer) for layer in bounds.outer_ring)
        reach_pc = (math.hypot(center[0], center[1]) + abs(center[2])
                    + (bounds.outer_ring_index + top_layer + 2) * edge_pc)
        radius_pc = min(radius_pc, reach_pc)
    candidates = []
    outside = 0
    for candidate in enumerate_sectors_within_radius(center, radius_pc, edge_pc):
        if bounds.contains(candidate[0], candidate[1]):
            candidates.append(candidate)
        else:
            outside += 1
    conn = store.get_connection(config)
    try:
        occupied = store.get_occupied_addresses(conn, {c[0] for c in candidates})
    finally:
        conn.close()
    return candidates, occupied, outside


def _neighborhood_batch(args, candidates, occupied, batch_density, suffix="", skip=()):
    """
    The not-yet-generated sectors among `candidates`
    (`_neighborhood_candidates`), as `_submit_batch` items, in
    the order `candidates` came (ring by ring, not nearest first: sort them
    by distance before applying any limit, GEN.101); `suffix` may hold
    `{distance}` (pc).

    Returns:
        tuple: `(batch, already_existed, skipped)`.
    """
    batch = []
    skipped = 0
    already_existed = 0
    for ring_index, layer_index, slot_index, x, y, z, distance_pc in candidates:
        address = (ring_index, layer_index, slot_index)
        if address in occupied or address in skip:
            already_existed += 1
            continue
        sector_args = batch_density.resolve(args, address, (x, y, z))
        if sector_args is None:
            log.debug(f"{_format_address(address)}: skipped (outside its layer's "
                      f"stored extent)")
            skipped += 1
            continue
        batch.append((address, (x, y, z), sector_args, suffix.format(distance=distance_pc)))
    return batch, already_existed, skipped


def run_local_neighborhood(args, edge_pc, progress):
    """
    Local-neighborhood mode: generates every not-yet-generated
    sector within `args.radius_pc` parsecs of `args.center_sector`'s own
    stored center -- see `run_ring_batch` for how each sector's density is set.

    Args:
        args (argparse.Namespace): Parsed arguments;
            `args.center_sector`/`args.radius_pc` must not be `None`.
        edge_pc (float): The sector edge length, in parsecs (`_edge_pc`).
        progress (rich.progress.Progress): `run_galaxy`'s shared progress
            display (see `run_ring_batch`). `run_random_start` also calls
            this directly after generating its own seed sector.

    Raises:
        SystemExit: If `args.center_sector` doesn't exist, or has never
                   been placed in a galaxy.
    """
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        try:
            center_position = store.get_sector_galaxy_position(conn, args.center_sector)
        except ValueError as exc:
            log.error(str(exc))
            raise SystemExit(1) from exc
    finally:
        conn.close()

    if center_position is None:
        log.error(
            f"sector_id={args.center_sector} has never been placed in a galaxy (its galaxy-position "
            f"columns are NULL) -- --center-sector requires an already galaxy-placed sector (one "
            f"generated via 'planetgen galaxy' itself, not 'planetgen sector'). Use --ring to "
            f"generate placed sectors from scratch instead."
        )
        raise SystemExit(1)

    center = (
        center_position["center_x_pc"], center_position["center_y_pc"], center_position["center_z_pc"],
    )
    mysql_config = store.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    center_ring, center_layer, _slot = sector_address_at(center, edge_pc)
    _require_inside(args, batch_density.bounds, center_ring, center_layer,
                    f"sector_id={args.center_sector} sits at {_format_address(sector_address_at(center, edge_pc))}")
    candidates, occupied, outside = _neighborhood_candidates(
        center, args.radius_pc, edge_pc, mysql_config, batch_density.bounds,
    )

    batch, already_existed, skipped = _neighborhood_batch(
        args, candidates, occupied, batch_density, f", {{distance:.2f}} pc from sector_id={args.center_sector}",
    )
    run_common._check_estimate(args, [item[2] for item in batch],
                    f"the sectors within {args.radius_pc:g} pc of sector {args.center_sector}", progress)
    with run_common.sector_bar(progress, args, "Sectors (local neighborhood)", len(batch)) as bar:
        _submit_batch(args, batch, f"Sectors (within {args.radius_pc:g} pc of sector {args.center_sector})",
                      edge_pc, progress, bar)
    generated = len(batch)

    skip_note = f", {skipped} skipped (outside the outline)" if skipped else ""
    outside_note = f", {outside} beyond the galaxy's edge left out" if outside else ""
    log.normal(
        f"Generated {generated} new sector(s) within {args.radius_pc} pc of sector_id={args.center_sector} "
        f"({len(candidates)} candidate slot(s) inside the galaxy{outside_note}, {already_existed} already "
        f"existed{skip_note})."
    )


class GenerationRefused(RuntimeError):
    """PERF.3: a bulk generation the database disk can't hold; the
    message says why and how much it needs."""


class MathCheckFailed(RuntimeError):
    """TEST.68: the math check (`planetgen.physics.mathcheck`) failed, so a
    bulk generation was refused before writing anything; the message
    names the failed checks."""


def require_math_check():
    """
    The gate in front of every bulk generation (TEST.68): runs the math
    check once per process (`mathcheck.startup_failures`, cached after
    that) and raises if any check failed.

    Raises:
        MathCheckFailed: Naming the failed checks.
    """
    failed = mathcheck.startup_failures()
    if failed:
        raise MathCheckFailed(
            f"the math check failed ({', '.join(r.name for r in failed)}), so nothing was generated. "
            f"Run 'planetgen check-math' for details.")


def generate_sector_neighborhood(center_sector_id, radius_ly=None, config=None, estimate_only=False):
    """
    Non-CLI counterpart to `run_local_neighborhood` -- for the admin web
    UI's "generate more sectors around this one" action
    (`planetgen/api/routes.py`'s `generate_sector_neighborhood_route`). Same
    work, a plain result dict instead of prints, and a catchable
    `ValueError` instead of `SystemExit` for an invalid/unplaced sector.

    The default radius is `program_constants.DEFAULT_GENERATE_RADIUS_PC`
    (12 pc, about 100 candidate addresses; GEN.23). Once they are all
    generated, the 100 ly around the center sector gets its bright stars
    (`backfill_bright_stars`, GEN.30).
    A larger `radius_ly` (100 ly holds roughly 1,500-2,000 addresses) can
    take minutes to hours.

    Args:
        center_sector_id (int): The already galaxy-placed sector to
            generate a neighborhood around.
        radius_ly (float, optional): Defaults to
            `program_constants.DEFAULT_GENERATE_RADIUS_PC` (12 pc).
        config (MySQLConfig, optional): Connection parameters.
        estimate_only (bool): Only work out the size and time (PERF.3);
            nothing is written.

    Returns:
        dict: `generated`, `already_existed`, `skipped`, `candidates`
            (addresses inside the galaxy), `outside_galaxy` (addresses in
            the sphere past the galaxy's edge, left out) -- all int --
            and `estimate` (`generationStats.Estimate.as_dict`).

    Raises:
        ValueError: If `center_sector_id` doesn't exist, has never been
                   placed in a galaxy, or lies outside its outline.
        GenerationRefused: The database disk can't hold it (nothing was
                   written).
        RuntimeError: If the galaxy's skeleton has never been built.
        MathCheckFailed: If the math check failed (a real run only);
            nothing is written.
    """
    edge_pc = run_common._edge_pc()
    radius_pc = (
        ly_to_pc(radius_ly) if radius_ly is not None
        else program_constants.DEFAULT_GENERATE_RADIUS_PC
    )
    config = config or store.DEFAULT_MYSQL_CONFIG
    if not estimate_only:
        require_math_check()

    conn = store.get_connection(config)
    try:
        center_position = store.get_sector_galaxy_position(conn, center_sector_id)
    finally:
        conn.close()

    if center_position is None:
        raise ValueError(
            f"sector_id={center_sector_id} has never been placed in a galaxy (its galaxy-position "
            f"columns are NULL) -- generating a neighborhood requires an already galaxy-placed sector "
            f"(one generated via 'planetgen galaxy', not 'planetgen sector')."
        )

    center = (center_position["center_x_pc"], center_position["center_y_pc"], center_position["center_z_pc"])
    batch_density = _BatchDensity(config)
    center_ring, center_layer, _slot = sector_address_at(center, edge_pc)
    if not batch_density.bounds.contains(center_ring, center_layer):
        raise ValueError(
            f"sector_id={center_sector_id} lies outside the galaxy: "
            f"{batch_density.bounds.describe_miss(center_ring, center_layer)}."
        )
    candidates, occupied, outside = _neighborhood_candidates(
        center, radius_pc, edge_pc, config, batch_density.bounds,
    )

    args = _default_generation_args(config=config)
    # Density-driven from the skeleton (_BatchDensity), not the flat
    # num_systems=10 default.
    args.density = None
    args.num_systems = None
    args.workers = 1   # one sector at a time, in whichever process called this

    batch, already_existed, skipped = _neighborhood_batch(args, candidates, occupied, batch_density)
    counts = {
        "already_existed": already_existed,
        "skipped": skipped,
        "candidates": len(candidates),
        "outside_galaxy": outside,
    }
    estimate = run_common._estimate_sectors(args, [item[2] for item in batch])
    if estimate_only:
        run_common._RUN_STATS.clear()
        return {"generated": 0, "estimate": estimate.as_dict(), **counts}
    if estimate.refusal:
        run_common._RUN_STATS.clear()
        raise GenerationRefused(estimate.refusal)

    generated = 0
    created_ids = []
    try:
        for address, position_pc, sector_args, _suffix in batch:
            started = time.monotonic()
            sector_id, _name, sector = generate_and_save_sector_at(sector_args, address, position_pc, edge_pc)
            created_ids.append(sector_id)
            run_common._record_sector(args, {"density": run_common._sector_density(sector_args), "systems": len(sector.entries),
                                  "stars": run_sector.sector_star_count(sector)}, time.monotonic() - started)
            generated += 1
    finally:
        run_common._finish_stats(args)
    if generated:
        backfill_bright_stars(config, center)  # GEN.30: once, around the requested sector
        # GEN.126: once the neighbours are final; this runs as the API's queued job.
        conn = store.get_connection(config)
        try:
            sector_paths.settle_sectors(conn, created_ids)
        finally:
            conn.close()
    return {"generated": generated, "estimate": estimate.as_dict(), **counts}


def run_random_start(args, edge_pc, progress):
    """
    Random-start mode (no `--ring`/`--center-sector` given): picks a
    random, not-yet-occupied sector address -- drawn only from
    inside the galaxy's stored outline (`GalaxyBounds.random_address`,
    every sector equally likely, optionally only out to `--max-ring`), so
    the start and the neighborhood around it are always in the galaxy --
    retried up to `program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS`
    times, generates it, then generates every not-yet-generated sector
    within `args.radius_pc` of it (default
    `program_constants.DEFAULT_GENERATE_RADIUS_PC`, 12 pc) via
    `run_local_neighborhood`. `--min-start-density` tightens the retry:
    an address below it is retried too.

    Raises:
        SystemExit: If no suitable address was found within the attempt
                   budget.
    """
    radius_pc = (
        args.radius_pc if args.radius_pc is not None
        else program_constants.DEFAULT_GENERATE_RADIUS_PC
    )

    mysql_config = store.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        sector_args = None
        for _ in range(program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS):
            address = batch_density.bounds.random_address(draw, max_ring=args.max_ring)
            if store.get_sector_id_at(conn, *address) is not None:
                continue
            position_pc = sector_position_pc(*address, edge_pc)
            sector_args = batch_density.resolve(args, address, position_pc)
            if sector_args is None:
                continue
            if args.min_start_density is not None and sector_args.density < args.min_start_density:
                continue
            break
        else:
            density_note = (
                f", meeting --min-start-density {args.min_start_density} (try lowering it or --max-ring "
                f"-- most volume-weighted draws land in the galaxy's own sparser outskirts)"
                if args.min_start_density is not None else ""
            )
            within = f"within {args.max_ring} rings" if args.max_ring is not None else "in the galaxy"
            log.error(
                f"Could not find an unoccupied sector address {within} after {program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS} attempts{density_note} "
                f"-- this galaxy may already be almost entirely generated within that range, or that range "
                f"may hold too little real stellar density; try a larger --max-ring."
            )
            raise SystemExit(1)
    finally:
        conn.close()

    # The estimate covers the start and its whole neighborhood, before
    # either is written.
    candidates, occupied, _outside = _neighborhood_candidates(
        position_pc, radius_pc, edge_pc, mysql_config, batch_density.bounds,
    )
    batch, _existed, _skipped = _neighborhood_batch(args, candidates, occupied, batch_density, skip={address})
    run_common._check_estimate(args, [sector_args] + [item[2] for item in batch],
                    f"a random start at {_format_address(address)} and the sectors within {radius_pc:g} pc of it",
                    progress)

    started = time.monotonic()
    sector_id, sector_name, sector = generate_and_save_sector_at(sector_args, address, position_pc, edge_pc)
    run_common._record_sector(args, {"density": run_common._sector_density(sector_args), "systems": len(sector.entries),
                          "stars": run_sector.sector_star_count(sector)}, time.monotonic() - started)
    designation = provisional_sector_designation(*address)
    log.normal(
        f"Saved random starting sector '{sector_name}' [{designation}] at {_format_address(address)} "
        f"(sector_id={sector_id})."
    )
    run_sector._log_summary(run_sector.sector_generation_summary_lines(sector, sector_args))

    args.center_sector = sector_id
    args.radius_pc = radius_pc
    run_local_neighborhood(args, edge_pc, progress)


def run_single_slot(args, edge_pc, progress):
    """
    Single-address mode: generates exactly the one sector at `(--ring I,
    --layer J, --slot K)` via the same `ensure_sector_generated` the
    galaxy map's "recalculate on visit" path uses -- the direct route from
    an address copied out of the interactive 3D Galaxy Map. With
    `--radius-pc`, then generates that sector's neighborhood
    (`run_local_neighborhood`), the Sector Map's "Generate neighborhood".

    Raises:
        SystemExit: If `--slot` is out of range for the ring, or the
                   address is outside the galaxy's stored outline.
    """
    address = (args.ring, args.layer, args.slot)
    total_slots = ring_sector_count(args.ring)
    if not (0 <= args.slot < total_slots):
        log.error(
            f"--slot {args.slot} is out of range for ring {args.ring} (holds {total_slots} slots, "
            f"0..{total_slots - 1})."
        )
        raise SystemExit(1)

    mysql_config = store.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    # Named outright, so generated even outside the outline (GEN.81).
    outside_ok = not _require_inside(args, batch_density.bounds, args.ring, args.layer,
                                     _format_address(address))
    position_pc = sector_position_pc(*address, edge_pc)
    conn = store.get_connection(mysql_config)
    try:
        exists = store.get_sector_id_at(conn, *address) is not None
    finally:
        conn.close()
    sector_args = None
    if not exists:
        # As `ensure_sector_generated` does it: density-driven.
        single = _default_generation_args(config=mysql_config)
        single.density = single.num_systems = None
        sector_args = batch_density.resolve(single, address, position_pc, outside_ok=outside_ok)
    sectors = [sector_args] if sector_args is not None else []
    what = _format_address(address)
    if args.radius_pc is not None:
        candidates, occupied, _outside = _neighborhood_candidates(
            position_pc, args.radius_pc, edge_pc, mysql_config, batch_density.bounds,
        )
        batch, _existed, _skipped = _neighborhood_batch(args, candidates, occupied, batch_density, skip={address})
        sectors += [item[2] for item in batch]
        what += f" and the sectors within {args.radius_pc:g} pc of it"
    run_common._check_estimate(args, sectors, what, progress)
    # The backfill waits for the end of the run (backfill_after_run), with
    # its own bar, instead of stalling this one at 0 of 1 (PERF.28).
    with run_common.sector_bar(progress, args, f"Sector ({_format_address(address)})", 1) as bar:
        result = ensure_sector_generated(*address, config=mysql_config, backfill=False, outside_ok=outside_ok,
                                         settle=False, link=False)  # the run settles and links at its end
        bar.update(advance=1)

    designation = provisional_sector_designation(*address)
    if result["created"]:
        log.normal(
            f"Saved sector '{result['sector_name']}' [{designation}] at {_format_address(address)} "
            f"(sector_id={result['sector_id']})."
        )
        run_sector._log_summary(result["summary"])
    else:
        log.normal(
            f"Sector [{designation}] at {_format_address(address)} already existed "
            f"(sector_id={result['sector_id']})."
        )

    if args.radius_pc is not None:
        # Then its neighborhood, the same way random-start mode follows
        # its seed sector.
        args.center_sector = result["sector_id"]
        run_local_neighborhood(args, edge_pc, progress)


def _layers_reaching(bounds, ring_index):
    """Every layer the galaxy's stored outline reaches at `ring_index`,
    lowest first (the column `galaxy_column` stores for that ring)."""
    return sorted(layer for layer in bounds.outer_ring if bounds.contains(ring_index, layer))


def _generate_addresses(args, addresses, what, edge_pc, progress, batch_density):
    """
    Generates every not-yet-generated address of `addresses`
    in order (up to `args.limit`, when set) -- the loop behind column and
    shell modes; see `run_ring_batch` for how each sector's density is set.
    """
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        occupied = set(store.get_occupied_addresses(conn, {a[0] for a in addresses}))
    finally:
        conn.close()

    pending = [a for a in addresses if a not in occupied]
    limit = getattr(args, "limit", None)
    batch = []
    skipped = 0
    for address in pending:
        if limit is not None and len(batch) >= limit:
            break
        position_pc = sector_position_pc(*address, edge_pc)
        sector_args = batch_density.resolve(args, address, position_pc)
        if sector_args is None:
            skipped += 1
            continue
        batch.append((address, position_pc, sector_args, ""))

    run_common._check_estimate(args, [item[2] for item in batch], what, progress)
    with run_common.sector_bar(progress, args, f"Sectors ({what})", len(batch)) as bar:
        _submit_batch(args, batch, f"Sectors ({what})", edge_pc, progress, bar)
    generated = len(batch)
    skip_note = f", {skipped} skipped (outside the outline)" if skipped else ""
    log.normal(
        f"Generated {generated} new sector(s) in {what} ({len(addresses)} total, "
        f"{len(addresses) - len(pending)} already existed{skip_note})."
    )


def block_layers(block):
    """
    Every sector layer a drill-down block spans, lowest first: a size-3
    block's three (`drill_slabs`), and for a bigger one the layers of
    each of its child slabs in turn.
    """
    m, _ring, _wedge, slab = block
    if m == 3:
        return drill_slabs(block)
    child_m = DRILL_LEVELS[DRILL_LEVELS.index(m) + 1]
    layers = []
    for child_slab in drill_slabs(block):
        layers.extend(block_layers(type(block)(child_m, 0, 0, child_slab)))
    return layers


def block_addresses(block, layer=None):
    """
    Every sector address `(ring, layer, slot)` inside drill-down block
    `block` (only on sector layer `layer` when given), ring by ring --
    a size-3 block's sectors (`drill_block_sectors`, design doc section
    3.5), and a bigger block's children's, recursively.
    """
    if block.m == 3:
        layers = drill_slabs(block) if layer is None else [layer]
        for sector_layer in layers:
            for sector in drill_block_sectors(block, sector_layer):
                yield (sector.ring, sector_layer, sector.wedge)
        return
    for _slab, children in drill_children(block):
        for child in children:
            if layer is not None and layer not in block_layers(child):
                continue
            yield from block_addresses(child, layer)


def run_block(args, edge_pc, progress):
    """
    Block mode: every sector the galaxy's outline allows inside one
    drill-down block (`--block`, optionally one `--block-layer`). Needs
    `--limit` or `--yes` when that is more than
    `LARGE_RING_WARNING_THRESHOLD` sectors.

    Raises:
        SystemExit: If the outline allows nothing in the block, or it is
                   too large and neither `--limit` nor `--yes` was given.
    """
    batch_density = _BatchDensity(store.mysql_config_from_args(args))
    bounds = batch_density.bounds
    key = format_drill_key(args.block)
    where = f"block {key}" + (f" layer {args.block_layer}" if args.block_layer is not None else "")
    confirmed = args.limit is not None or args.yes
    addresses = []
    for address in block_addresses(args.block, args.block_layer):
        if not bounds.contains(address[0], address[1]):
            continue
        addresses.append(address)
        if not confirmed and len(addresses) == LARGE_RING_WARNING_THRESHOLD + 1:
            run_common._refuse_or_warn(args, f"{where.capitalize()} holds more than {LARGE_RING_WARNING_THRESHOLD} "
                                  f"sector slots, a very large run (--limit N generates only the first N).")
    if not addresses:
        run_common._refuse_or_warn(args, f"The galaxy's outline allows no sector in {where}, so there is nothing "
                              f"to generate.")
        return
    _generate_addresses(args, addresses, where, edge_pc, progress, batch_density)


def run_column(args, edge_pc, progress):
    """
    Column mode: every sector at `(--ring I, --slot K)` through every
    layer the galaxy's outline reaches at ring I.

    Raises:
        SystemExit: If `--slot` is out of range for the ring, or the
                   outline doesn't reach the ring at all.
    """
    total_slots = ring_sector_count(args.ring)
    if not (0 <= args.slot < total_slots):
        log.error(
            f"--slot {args.slot} is out of range for ring {args.ring} (holds {total_slots} slots, "
            f"0..{total_slots - 1})."
        )
        raise SystemExit(1)
    batch_density = _BatchDensity(store.mysql_config_from_args(args))
    layers = _layers_reaching(batch_density.bounds, args.ring)
    if not layers:
        _require_inside(args, batch_density.bounds, args.ring, 0, f"ring {args.ring}")
        log.normal(f"No layer of the outline reaches ring {args.ring}: nothing to generate.")
        return
    addresses = [(args.ring, layer, args.slot) for layer in layers]
    _generate_addresses(args, addresses, f"the column at ring {args.ring} slot {args.slot}",
                        edge_pc, progress, batch_density)


def run_shell(args, edge_pc, progress):
    """
    Shell mode: every sector of ring `--ring I` through every layer the
    galaxy's outline reaches there. Needs `--limit` or `--yes` when that
    is more than `LARGE_RING_WARNING_THRESHOLD` sectors.

    Raises:
        SystemExit: If the outline doesn't reach the ring, or the shell is
                   too large and neither `--limit` nor `--yes` was given.
    """
    batch_density = _BatchDensity(store.mysql_config_from_args(args))
    layers = _layers_reaching(batch_density.bounds, args.ring)
    if not layers:
        _require_inside(args, batch_density.bounds, args.ring, 0, f"ring {args.ring}")
        log.normal(f"No layer of the outline reaches ring {args.ring}: nothing to generate.")
        return
    total_slots = ring_sector_count(args.ring)
    total = total_slots * len(layers)
    if total > LARGE_RING_WARNING_THRESHOLD and args.limit is None and not args.yes:
        run_common._refuse_or_warn(args, f"The shell at ring {args.ring} holds {total} sector slots ({total_slots} "
                              f"slots across {len(layers)} layers), a very large run (--limit N generates "
                              f"only the first N).")
    addresses = [(args.ring, layer, slot) for layer in layers for slot in range(total_slots)]
    _generate_addresses(args, addresses, f"the shell at ring {args.ring}", edge_pc, progress, batch_density)


def run_galaxy(args):
    """
    Dispatches to block, column, shell, single-address, ring-batch,
    local-neighborhood, or random-start mode, owning the one `rich.progress.Progress` display
    they share. First checks the galaxy has been planned at the standard
    sector edge, since every mode validates its addresses against that
    plan's stored outline before generating anything.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "galaxy"`).
    """
    edge_pc = run_common._edge_pc()
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        bounds = store.get_galaxy_bounds(conn)
        warning = version_check.galaxy_warning(conn)
    finally:
        conn.close()
    if warning:
        log.normal(f"Warning: {warning}")  # DB.7: a heads-up before it extends a mixed-version galaxy
    if bounds is None:
        log.error("The galaxy has never been planned -- run 'planetgen plan' first, so every sector "
                  "can be checked against the galaxy's bounds before it is generated.")
        raise SystemExit(1)
    if not math.isclose(bounds.edge_pc, edge_pc):
        log.error(f"The stored plan was built with {bounds.edge_pc:g} pc sectors, not the standard "
                  f"{edge_pc:g} pc -- re-run 'planetgen plan'.")
        raise SystemExit(1)
    if not bounds:
        log.error("The stored plan has no layers at all (nothing in this galaxy shape expects a star "
                  "per sector) -- re-run 'planetgen plan' with a different shape.")
        raise SystemExit(1)

    estimate_only = getattr(args, "estimate_only", False)
    args.link_later = not estimate_only
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        started_at = store.database_now(conn)
    finally:
        conn.close()
    try:
        with run_common._generation_progress(disable=estimate_only) as progress:
            log.set_console(progress.console)
            try:
                _run_galaxy_mode(args, edge_pc, progress)
            finally:
                log.reset_console()
    except run_common._EstimateOnly as exc:
        run_common._print_estimate(exc)
        return
    finally:
        run_common._finish_stats(args)
        link_after_run(args, started_at)
    # GEN.30: the bright stars come after the sectors, so neither the
    # scatter nor the backfill draws stars for a sector the run filled.
    if getattr(args, "then_scatter", False):
        run_plan.scatter_bright_stars(args)
    backfill_after_run(args, edge_pc, started_at)
    run_population.run_population_after(args)
    settle_after_run(args, started_at)
    run_common._finish_stats(args)    # the steps after the sectors recorded their speeds too


def _run_galaxy_mode(args, edge_pc, progress):
    """`run_galaxy`'s dispatch to the mode its arguments ask for."""
    if getattr(args, "block", None) is not None:
        run_block(args, edge_pc, progress)
    elif args.column:
        run_column(args, edge_pc, progress)
    elif args.shell:
        run_shell(args, edge_pc, progress)
    elif args.slot is not None:
        run_single_slot(args, edge_pc, progress)
    elif args.ring is not None:
        run_ring_batch(args, edge_pc, progress)
    elif args.center_sector is not None:
        run_local_neighborhood(args, edge_pc, progress)
    else:
        run_random_start(args, edge_pc, progress)
