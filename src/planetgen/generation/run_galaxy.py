# planetgen/generation/run_galaxy.py

"""
Run Galaxy
==========

The `galaxy` command: many sectors placed as one galaxy, each at its
shell/band/wedge address in the spherical grid, with its systems inside
the sector's own cell. Modes: `--shell` (a shell, or with `--band` one band
of it, or with `--wedge` too one sector), `--shells/--bands/--wedges` (a
span), `--center-sector`/`--radius-pc` (local neighborhood), `--course`,
or none (random start).
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

from planetgen.datasets import catalogue as datasetCatalogue
from planetgen.db import region, sector_paths, store
from planetgen.generation import bright_stars as brightStars
from planetgen.galaxy import density, resolution as galaxyResolution, seed as galaxySeed, version_check
from planetgen.physics import mathcheck
from planetgen.queue import redisqueue
from planetgen import tuning as program_constants
from planetgen.util import log
from planetgen.galaxy.profile import profile_at, sector_profile
from planetgen.galaxy.region import BALL_EDGE_ALLOWANCE, Region
from planetgen.galaxy import spherical
from planetgen.galaxy.frame import galactic_radius_pc
from planetgen.galaxy.span import Span
from planetgen.galaxy.spherical import (
    SectorCell, enumerate_sectors_within_radius, neighbor_addresses,
    provisional_sector_designation, sector_address_at, sector_position_pc,
)
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.generation import fill_order, object_timing
from planetgen.generation import nebula_contents as nebulaContents
from planetgen.generation import run_common
from planetgen.generation import run_plan
from planetgen.generation import run_population
from planetgen.generation import run_sector
from planetgen.generation import stages, steps
from planetgen.generation import stats as generationStats
from planetgen.util import draw
from planetgen.util import detmath


LARGE_RUN_WARNING_THRESHOLD = 2000
"""int: `--shell I` (and a span) requires `--limit` or `--yes` when it holds
more sectors than this -- a band of shell 200 holds a few hundred, a whole
shell from about shell 12 on. Anything larger takes a real, unbounded amount
of time and disk, so it needs an explicit choice."""


BACKFILL_FROM_CHOICES = ("edge", "none")
"""tuple: `--backfill-from`'s choices (GEN.30, GEN.98)."""


def _format_address(address):
    """`(shell, band, wedge)` as the words the CLI prints."""
    shell_index, band_index, wedge_index = address
    return f"shell {shell_index} band {band_index} wedge {wedge_index}"


class _BatchDensity:
    """
    Checks each address against the galaxy's reach and resolves each
    sector's own `--density` from its sector profile (`galaxy.profile`, the
    skeleton's real position-based `relative_density`), for `galaxy` mode's
    batch/local-neighborhood/random-start generation -- the same mechanism
    `ensure_sector_generated` uses for a single lazily-generated sector,
    applied across a whole run instead of every sector sharing one flat CLI
    value.

    The reach check (`galaxy.skeleton.GalaxyBounds.contains`) always
    applies and is the only refusal (GEN.235): a cell the grid has no
    address for is never generated, even with an explicit
    `--density`/`--num-systems`. Past it, an explicit flag is a deliberate,
    uniform override for the whole run (`resolve` returns `args`
    unchanged); otherwise `resolve` sets each sector's density from its
    profile. A sparse sector, in the outline or far out in the halo, is
    generated like any other (GEN.76).

    The skeleton and bounds are fetched once, on first use -- neither
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
        """The stored outline, a `galaxy.skeleton.GalaxyBounds`."""
        self._load()
        return self._bounds

    def profile(self, address):
        """The address's `galaxy.profile.SectorProfile`."""
        return sector_profile(address, self.skeleton)

    def resolve(self, args, address, position_pc):
        """
        Args:
            args (argparse.Namespace): The `galaxy` subcommand's own parsed
                (and validated) arguments.
            address (tuple): This sector's `(shell, band, wedge)`.
            position_pc (tuple): This sector's `(x, y, z)` center, parsecs.

        Returns:
            argparse.Namespace or None: `None` if this address is beyond the
                grid's reach. Otherwise `args` itself when a
                density/count was given explicitly, else a fresh copy with
                `.density` set to this address's profile density (its
                `relative_density`) and `.num_systems` cleared. A sparse
                address is never skipped (GEN.76): the predicted density only
                sets how many systems the draw expects, and the draw always
                runs.
        """
        if not self.bounds.contains(address[0], address[1]):
            return None
        if args.density is not None or args.num_systems is not None:
            return args

        resolved = copy.copy(args)
        resolved.density = profile_at(position_pc, self.skeleton.shape, self.skeleton.edge_pc,
                                      self.skeleton.expected_system_count_at_density_1).density
        resolved.num_systems = None
        return resolved


class _Pending:
    """
    The sectors a run still has to generate, walked and not listed (PERF.36): a batch of a million sectors would hold a
    million resolved copies of the arguments, so a run counts its sectors in a first pass (`survey`), works its estimate
    out from a sample of them (`sample_args`, over `generationStats.SAMPLE_ABOVE`; below that, `items` as a list) and
    resolves each sector only as it is handed to the work queue (`items`, one more pass). The same occupancy,
    reach and `--limit` rules decide every pass.

    Attributes:
        count (int): The sectors to generate; set by `survey`.
        existed (int): The addresses left out because they were already generated, or named in `skip`.
        skipped (int): The addresses left out because they lie beyond the grid's reach.
    """

    def __init__(self, args, walk, place, occupied, batch_density, skip=(), limit=None):
        """
        Args:
            args (argparse.Namespace): The run's arguments.
            walk (callable): `walk()` -> a new iterable of `(address, extra)` each time, in the order to generate.
            place (callable): `place(address, extra)` -> `(position_pc, suffix)` of the sector.
            occupied (set): The addresses already generated.
            batch_density (_BatchDensity): Resolves each sector's own arguments.
            skip (iterable): Addresses to leave out too.
            limit (int or None): Generate at most this many.
        """
        self.args = args
        self._walk = walk
        self._place = place
        self._occupied = occupied
        self._density = batch_density
        self._skip = skip
        self._limit = limit
        self.count = self.existed = self.skipped = 0

    @classmethod
    def of_addresses(cls, args, addresses, occupied, batch_density, edge_pc, limit=None):
        """The addresses `addresses()` yields (a callable giving a fresh iterable each time)."""
        return cls(args, lambda: ((address, None) for address in addresses()),
                   lambda address, _extra: (sector_position_pc(*address, edge_pc), ""),
                   occupied, batch_density, limit=limit)

    @classmethod
    def of_candidates(cls, args, candidates, occupied, batch_density, suffix="", skip=()):
        """The `_neighborhood_candidates` list, with `suffix` (it may hold `{distance}`, pc) on each."""
        return cls(args, lambda: ((candidate[:3], candidate) for candidate in candidates),
                   lambda _address, candidate: (candidate[3:6], suffix.format(distance=candidate[6])),
                   occupied, batch_density, skip=skip)

    def _eligible(self):
        """The `(address, extra)` to generate, counting `existed` and `skipped` as it goes."""
        existed = skipped = taken = 0
        bounds = self._density.bounds
        for address, extra in self._walk():
            if address in self._occupied or address in self._skip:
                existed += 1
                continue
            if self._limit is not None and taken >= self._limit:
                break
            if not bounds.contains(address[0], address[1]):
                log.debug(f"{_format_address(address)}: skipped (beyond the grid's reach)")
                skipped += 1
                continue
            taken += 1
            yield address, extra
        self.count, self.existed, self.skipped = taken, existed, skipped

    def survey(self):
        """Counts the sectors to generate (`count`, `existed`, `skipped`)."""
        for _entry in self._eligible():
            pass
        return self.count

    def items(self):
        """Each sector to generate as a `_submit_batch` item `(address, position_pc, sector_args, suffix)`, resolved as
        it is reached."""
        for address, extra in self._eligible():
            position_pc, suffix = self._place(address, extra)
            yield address, position_pc, self._density.resolve(self.args, address, position_pc), suffix

    def sample_args(self, size=None):
        """The resolved arguments of `size` sectors spread evenly through the run (`generationStats.sample_positions`);
        `survey` first."""
        wanted = generationStats.sample_positions(self.count, size)
        sample = []
        if not wanted:
            return sample
        for place, (address, extra) in enumerate(self._eligible()):
            if place == wanted[len(sample)]:
                position_pc, _suffix = self._place(address, extra)
                sample.append(self._density.resolve(self.args, address, position_pc))
                if len(sample) == len(wanted):
                    break
        return sample


def _plan_pending(pending):
    """
    Counts what `pending` (a `_Pending`) has to generate and chooses how the run goes on (PERF.36): up to
    `generationStats.SAMPLE_ABOVE` sectors, a list with the estimate worked out from all of them; above it, the sectors
    streamed as they are submitted and the estimate worked out from a sample.

    Returns:
        tuple: `(batch, estimate_args, total)`: the items to submit (a list or a generator), the arguments the estimate
            reads and the size of the run when they are a sample of it (else `None`).
    """
    if pending.survey() > generationStats.SAMPLE_ABOVE:
        return pending.items(), pending.sample_args(), pending.count
    batch = list(pending.items())
    return batch, [item[2] for item in batch], None


HOT_HOST_PATTERN = "^(O[0-9]|B[012])V"
"""str: The spectral types that can light an emission region (`tuning.NEBULA_HOST_RULES`: O and B0 to B2 rows)."""


def _emission_hosts(conn, skeleton, address, profile):
    """
    The emission regions (D and E) lit by the pre-placed hot stars of the sectors around `address` (GEN.321): every
    `bright_stars` row, built or not, within the reach of a region's members of the sector's cube, rolled as its own
    sector rolled it (`nebula_contents.row_host`, with that sector's `"nebula"` factor), as `(nebula, host_name,
    host_position_pc)`. The sector's own stars are left out: its fill makes their regions itself. A hot star the
    scatter did not place is no row, so its region reaches no neighbour.
    """
    reach = profile.edge_pc * math.sqrt(3) / 2 + nebulaContents.emission_reach_pc()
    rows = region.region_rows(
        conn, "bright_stars",
        "p.shell_index, p.band_index, p.wedge_index, p.star_type, p.age_gy, p.mass_kg, p.temperature_k, "
        "p.position_x_mpc, p.position_y_mpc, p.position_z_mpc",
        ("p.position_x_mpc / 1000.0", "p.position_y_mpc / 1000.0", "p.position_z_mpc / 1000.0"),
        profile.position_pc, reach, where=f"p.star_type REGEXP '{HOT_HOST_PATTERN}'")
    factors, found = {}, []
    for row in sorted((dict(r) for r in rows),
                      key=lambda r: (r["shell_index"], r["band_index"], r["wedge_index"], r["position_x_mpc"],
                                     r["position_y_mpc"], r["position_z_mpc"])):
        host_address = (row["shell_index"], row["band_index"], row["wedge_index"])
        if host_address == tuple(address) or nebulaContents.row_host(row) is None:
            continue
        if host_address not in factors:
            factors[host_address] = sector_profile(host_address, skeleton).kind_factors.get("nebula", 1.0)
        host = nebulaContents.row_host(row, factors[host_address])
        if host is not None:
            found.append(host)
    return found


def _fill_context(args, address):
    """The `brightStars.FillContext` for one galaxy sector: its sector
    profile (`galaxy.profile`: the population mix and the rest), its unfilled pre-placed bright stars down to its own level and
    mass (`store.bright_star_fill_level`, `store.bright_star_fill_mass_limit`:
    its backfill's, GEN.44, GEN.187, else the galaxy scatter's), and its unbuilt scattered phenomena when the galaxy's
    phenomenon scatter ran (GEN.100). `None` without a stored skeleton
    (nothing to take the mix from)."""
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            return None
        profile = sector_profile(address, skeleton)
        limits, backfill_seed = run_plan._scatter_limits_and_seed(conn, skeleton)
        resolution = galaxyResolution.cell_resolution(
            profile, skeleton.cell_star_ceiling, skeleton.fill_luminosity_cap_sol,
            limits.mass_limit_at(address[0], address[1]))
        if resolution.capped:
            # GEN.281: a cell above the ceiling is taken down to the cap and no further.
            run_plan._draw_sector_bands(conn, skeleton, [tuple(address)], resolution.luminosity_cap_sol,
                                        limits, backfill_seed)
        phenomena = below_cut = None
        settings = store.phenomenon_scatter_settings(conn)
        if settings is not None:
            phenomena = store.phenomena_for_sector(conn, *address)
            seed, min_mass_solar = settings
            if min_mass_solar is not None and not resolution.capped:
                # GEN.168: what the scatter left below its cut, at the
                # density it used (the plan's, not this sector's own roll). A capped cell makes no remnant of its
                # own: they are dark, below the cap, and part of what it leaves unresolved.
                expected = profile.expected_systems
                below_cut = brightStars.BelowCut(tuple(address), min_mass_solar, seed, expected)
        level = store.bright_star_fill_level(conn, *address)
        mass_limit = store.bright_star_fill_mass_limit(conn, *address)
        emission_hosts = _emission_hosts(conn, skeleton, address, profile)
        if level is None and mass_limit is None:
            return brightStars.FillContext(profile, phenomenon_rows=phenomena, below_cut=below_cut,
                                           resolution=resolution, emission_hosts=emission_hosts)
        rows = store.bright_stars_for_sector(conn, *address)
    finally:
        conn.close()
    return brightStars.FillContext(profile, rows, min_luminosity_sol=level,
                                   phenomenon_rows=phenomena, below_cut=below_cut, star_mass_limit_sol=mass_limit,
                                   resolution=resolution, emission_hosts=emission_hosts)


def backfill_step_targets(addresses, bounds, step_masses=None):
    """
    The mass backfill's steps (GEN.187) outward from some sectors: step 1 is
    the sectors a face away from `addresses` (no diagonals,
    `spherical.neighbor_addresses`), step 2 the ones a face away from step 1,
    and so on, one step per entry of `step_masses` (default
    `program_constants.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL`). Cells outside the
    galaxy's outline (`bounds`) are left out.

    Returns:
        dict: `(shell, band, wedge)` -> the lightest initial mass (solar
            masses) its step takes stars down to; the nearest step has the
            lowest. `addresses` themselves are not in it.
    """
    masses = program_constants.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL if step_masses is None else step_masses
    seen = set(addresses)
    frontier = sorted(seen)
    targets = {}
    for mass in masses:
        following = set()
        for address in frontier:
            for neighbor in neighbor_addresses(*address):
                if neighbor in seen or not bounds.qualifies(neighbor[0], neighbor[1]):
                    continue
                seen.add(neighbor)
                following.add(neighbor)
        for address in following:
            targets[address] = float(mass)
        frontier = sorted(following)
        if not frontier:
            break
    return targets


def backfill_bright_stars(config, center_pc):
    """`backfill_bright_stars_around` one center (see there)."""
    return backfill_bright_stars_around(config, [center_pc])


def backfill_bright_stars_around(config, centers_pc, progress=None):
    """
    The bright-star backfill around generated sectors (GEN.23; by initial
    mass in steps since GEN.187, per sector since GEN.44): the sectors a face
    away from `centers_pc`' sectors get every star born from 1 solar mass
    up, the step after that from 2, then 5, then 8
    (`backfill_step_targets`); past the fourth step the scatter's own mass
    limit stands. A sector filled, or already taken that deep (by an earlier
    backfill, or by the scatter's mass pass), is skipped, and one a nearer
    sector reaches later is topped up with only the masses it lacks, so no
    star is ever drawn twice. A sector with no level anywhere that still
    holds stars is left over from a run that failed part way: they are wiped
    and it is drawn whole.

    Sectors are drawn a few hundred at a time, under their rows' locks,
    and each one's new mass is written in the transaction that writes its
    stars. A sector's draw is seeded from the galaxy's scatter seed, its
    address and fixed mass bands (`brightStars.backfill_mass_cells`), so it
    is the same whichever sector reached it first, and whatever steps took
    it there.

    Args:
        config (MySQLConfig): Connection parameters.
        centers_pc (list): Generated sectors' `(x, y, z)` centers,
            galaxy-frame parsecs.
        progress (Progress, optional): A live `_generation_progress`
            display: the backfill adds its own bar, counting the sectors
            it visits, with its ETA (PERF.28).

    Returns:
        dict: `sectors` (drawn now) and `stars` (placed now), both int;
            zeros without a stored skeleton.
    """
    summary = {"sectors": 0, "stars": 0}
    shells = {}
    # Finding the sectors to draw can take a while on a large run; the
    # bar is there from the start (ADM.26), unmeasured until the count is known.
    bar = steps.Step("Bright-star backfill (finding sectors)", "backfill", None, progress=progress,
                     stats=run_common._stats_for_config(config)).__enter__()
    failed = True
    conn = store.get_connection(config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            return summary
        limits, seed = run_plan._scatter_limits_and_seed(conn, skeleton)
        bounds = store.get_galaxy_bounds(conn)
        centers = {sector_address_at(center, skeleton.edge_pc) for center in centers_pc}
        targets = backfill_step_targets(centers, bounds)
        # GEN.333: a sector backfills the masses under its own band's mass limit.
        targets = {address: mass for address, mass in targets.items()
                   if limits.mass_limit_at(address[0], address[1]) is None
                   or mass < limits.mass_limit_at(address[0], address[1])}
        held = store.sector_bright_masses(conn, targets)
        filled = store.get_occupied_addresses(conn, {address[0] for address in targets})
        todo = {address: mass for address, mass in targets.items()
                if address not in filled and mass < held.get(address, math.inf)}
        on_sector = None
        if todo:
            bar.update(description="Bright-star backfill (sectors)", total=len(todo))

            def on_sector():
                bar.update(advance=1)
        drawn = run_plan._draw_sector_masses(conn, skeleton, todo, limits, seed,
                                             on_sector=on_sector, shells=shells)
        summary["sectors"], summary["stars"] = drawn
        failed = False
    except BaseException:
        conn.rollback()
        raise
    finally:
        conn.close()
        bar.close(success=not failed)     # nothing to draw: no bar left unmeasured
    run_plan.log_shells("Bright-star backfill", shells)
    if summary["sectors"]:
        log.debug(f"bright-star backfill: {summary['stars']} stars in {summary['sectors']} sector(s), masses "
                  f"{','.join(f'{mass:g}' for mass in program_constants.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL)} Msun")
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
            stages.skip("backfill", "the run generated no sector")
            return summary
        centers = generated
    finally:
        conn.close()
    stages.enter("backfill")
    log.normal("Backfilling the bright stars out from the edge of the run...")
    # Its own bar and ETA (PERF.28): a backfill can take minutes.
    with run_common._generation_progress() as progress:
        log.set_console(progress.console)
        try:
            summary = backfill_bright_stars_around(config, centers, progress=progress)
        finally:
            log.reset_console()
    stages.note(sectors=summary["sectors"], objects=summary["stars"])
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
            stages.skip("link", "the run created no sectors")
            return 0
        stages.enter("link")
        stages.note(sectors=len(created))
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
            stages.skip("settle", "the run created no sectors")
            return 0
        stages.enter("settle")
        stages.note(sectors=len(created))
        log.normal("Saving the sector paths...")
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step("Sector paths", "paths", len(created), args=args, progress=progress) as bar:
                    def on_progress(done, total):
                        bar.update(completed=done, total=total)

                    saved = sector_paths.settle_sectors(conn, created, on_progress,
                                                       workers=run_common._worker_count(args))
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
        address (tuple): This sector's `(shell, band, wedge)`.
        position_pc (tuple): Its `(x, y, z)` center, parsecs.
        edge_pc (float): The sector edge length, parsecs.
        channel (optional): Where a worker reports the save's progress
            (`steps.worker_step`, PERF.50); `None` reports nowhere.

    Returns:
        tuple: `(sector_id, sector_name, sector)` of the newly saved sector
            -- `sector_name` is read back after saving, since name
            uniqueness (v22) may rename it on save.
    """
    shell_index, band_index, wedge_index = address
    x, y, z = position_pc
    radius_pc = galactic_radius_pc(position_pc)
    cell = SectorCell.for_address(shell_index, band_index, wedge_index, pc_to_ly(edge_pc))

    fill = _fill_context(args, address)
    galaxy_position = {
        "center_x_pc": x, "center_y_pc": y, "center_z_pc": z,
        "galactic_radius_pc": radius_pc,
        "shell_index": shell_index, "band_index": band_index, "wedge_index": wedge_index,
    }
    galaxy_seed = _galaxy_seed(args)
    # GEN.47: the galaxy's molecular clouds that reach this sector.
    cloud_field = None
    if fill is not None:
        cloud_field = (galaxy_seed, fill.shape, tuple(position_pc), edge_pc * math.sqrt(3) / 2)
    # GEN.39: every draw for this sector, its save included, comes from
    # its own seed, so the run's order and worker count don't matter.
    with galaxySeed.seeded(galaxy_seed, "sector", address,
                           disk_scale_length_pc=fill.shape.disk_scale_length_pc if fill is not None else None
                           ), object_timing.unit(
            ",".join(str(part) for part in address)) as timings:
        _sector_name, sector = run_sector.generate_directed_sector(
            args, galaxy_seed=galaxy_seed, address=address,
            galactic_center_dist_ly=pc_to_ly(radius_pc), cell=cell, fill=fill, cloud_field=cloud_field)
        if address == run_sector.NUCLEUS_ADDRESS and not (fill is not None and fill.phenomena_scattered):
            run_sector.add_galactic_nucleus(
                sector, args, pc_to_ly(radius_pc),
                mass_solar=None if fill is None else density.central_black_hole_mass_solar(fill.shape))
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
    _save_object_timing(args, timings, sector_id)
    return sector_id, sector.name, sector


def _save_object_timing(args, timings, sector_id):
    """Writes the sector's per-object timings (PERF.92) under the run that
    made it; never fails the run."""
    run_id = getattr(args, "generation_run_id", None)
    if run_id is None or not object_timing.stats_enabled():
        return
    try:
        conn = store.get_connection(store.mysql_config_from_args(args))
        try:
            object_timing.save(conn, run_id, timings, sector_id)
            conn.commit()
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- the timings never fail a run
        log.debug(f"Object timings: can't save them for sector {sector_id} ({exc}).")


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


def ensure_sector_generated(shell_index, band_index, wedge_index, config=None, backfill=True,
                            settle=True, link=True):
    """
    The galaxy map's "recalculate on visit" entry point: returns the
    sector already generated at this address if one exists; otherwise
    checks the address is a cell of the grid's reach (`GalaxyBounds.contains`,
    the only refusal, GEN.235) and generates and saves it on the spot from
    its sector profile (`galaxy.profile`), using the profile's density (its
    position's `relative_density`) as the `--density` multiplier. A sparse
    address generates too (GEN.76): its draw may come out with no systems,
    and the sector is still saved and marked generated, in the disk or out
    in the halo.

    A concurrent visit to the same never-generated address is handled by
    `sectors`'s `UNIQUE (shell_index, band_index, wedge_index)`: the
    losing `INSERT` raises `pymysql.err.IntegrityError`, caught here and
    turned into "return what the other call just created".

    A new sector's bright-star backfill runs here too, unless `backfill`
    is false: `planetgen galaxy --shell --band --wedge` leaves it to the end of the run
    (`backfill_after_run`), where it has its own progress bar (PERF.28).

    A new sector's paths (and its neighbours') are saved by a queued job
    unless `settle` is false (GEN.126): the caller that settles at the end
    of its own run, or after its own edit, passes false. Likewise the sector is
    linked to its neighbours (containment, nearest systems) as it is saved
    unless `link` is false (PERF.45): the run that links its sectors at its
    end passes false.

    Returns:
        dict: `created` (bool), `qualifies` (bool: inside the galaxy's
              stored outline, where the disk and bulge expect a system; a halo
              sector is generated too; `True` for a sector that already existed), `in_reach` (bool: a cell of the grid's
              reach; `False` means nothing was generated), `sector_id` (int or `None`), `sector_name`
              (str or `None`, only when `created`), and `summary`
              (`sector_generation_summary_lines`, only when `created`).

    Raises:
        ValueError: If `wedge_index` is out of range for the band.
        RuntimeError: If the galaxy's skeleton has never been built --
                     run `planetgen plan` first.
    """
    address = (shell_index, band_index, wedge_index)
    conn = store.get_connection(config)
    try:
        existing_id = store.get_sector_id_at(conn, *address)
        if existing_id is not None:
            return {"created": False, "qualifies": True, "in_reach": True, "sector_id": existing_id, "sector_name": None}

        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError(
                "The galaxy's skeleton has never been built (no galaxy_shape row) -- run "
                "'planetgen plan' first."
            )

        bounds = store.get_galaxy_bounds(conn)
    finally:
        conn.close()

    if not bounds.contains(shell_index, band_index):
        # Past the grid's reach: the cell has no address (the object ID cannot name it).
        return {"created": False, "qualifies": False, "in_reach": False, "sector_id": None, "sector_name": None}

    position_pc = sector_position_pc(shell_index, band_index, wedge_index, skeleton.edge_pc)
    qualifies = bounds.qualifies(shell_index, band_index)

    # Every address in the reach generates, however sparse (GEN.76): the
    # profile's density sets the draw's expected count, never whether the
    # draw happens.
    args = _default_generation_args(config=config)
    args.density = sector_profile(address, skeleton).density
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
        return {"created": False, "qualifies": True, "in_reach": True, "sector_id": existing_id, "sector_name": None}

    if backfill:
        backfill_bright_stars(config, position_pc)  # GEN.30: around the requested sector
    if settle:
        queue_settle(config, [sector_id])
    return {"created": True, "qualifies": qualifies, "in_reach": True, "sector_id": sector_id, "sector_name": sector_name,
            "summary": run_sector.sector_generation_summary_lines(sector, args)}


def _log_saved(saved, address, suffix=""):
    designation = provisional_sector_designation(*address)
    log.normal(
        f"Saved sector '{saved['name']}' [{designation}] at {_format_address(address)}{suffix} "
        f"(sector_id={saved['sector_id']})."
    )
    run_sector._log_summary(saved["summary"])


def _submit_batch(args, batch, title, edge_pc, progress, bar, count=None):
    """Queues every `(address, position_pc, sector_args, suffix)` of
    `batch` (`_submit_sector`) and waits for them. `batch` may be a generator (PERF.36): then `count` says how many
    it holds; the queue takes them as its workers are free, so a run of millions is never held at once."""
    count = len(batch) if count is None else count
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
            queue.expect(count)
            log.normal(f"Generating {count:,} sector(s) with {queue.workers} worker(s); each is reported below as "
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


def _require_in_reach(bounds, shell_index, band_index, what):
    """
    Stops the run before anything is generated when `(shell_index, band_index)`
    is not a cell of the grid's reach (GEN.235): the only refusal. A cell outside
    the outline but inside the reach is generated from its profile like any other.

    Raises:
        SystemExit: Beyond the reach.
    """
    if bounds.contains(shell_index, band_index):
        return
    log.error(f"{what} -- {bounds.describe_beyond_reach(shell_index, band_index)}.")
    raise SystemExit(1)


def ball_extent(sectors, edge_pc):
    """
    The radius, in parsecs, of a radial fill of `sectors` sectors across
    (ADM.30): `sectors` edges plus `BALL_EDGE_ALLOWANCE`, so one sector holds
    the centre and its face neighbours and a bigger one a round shape joined
    by face adjacency. The fill is a ball (`--radius-pc`), so the Generate
    page turns its "sectors" field into the radius with this.
    """
    return (sectors + BALL_EDGE_ALLOWANCE) * edge_pc


def _region_text(args):
    """The fill's size in words, for logs and the estimate."""
    return f"within {args.radius_pc:g} pc"


class _ReachOnly:
    """Stands for a `GalaxyBounds` where only the grid's reach limits a request (`Region` asks `qualifies`)."""

    def __init__(self, bounds):
        self.qualifies = bounds.contains


def _is_far(center, edge_pc, bounds):
    """
    Whether the cell at `center` lies outside the galaxy's outline: the pull from a
    cell inside it stays in the galactic region (`qualifies`), and one from a cell
    outside, a halo sector someone asked for, takes every cell of the reach
    (`contains`), each generated from its own profile (GEN.235).
    """
    shell, band, _wedge = sector_address_at(center, edge_pc)
    return not bounds.qualifies(shell, band)


def _report_far_neighborhood(args, center, radius_pc, edge_pc, bounds, skeleton):
    """
    A neighbourhood request centred outside the outline says, before anything is
    enumerated, how many sectors it holds and how many systems they expect at
    the centre's profile (a cubic kiloparsec of empty halo is 1.6e7 cells), and
    warns (`--strict` stops) over `LARGE_RUN_WARNING_THRESHOLD`.
    """
    address = sector_address_at(center, edge_pc)
    count = Region(center, radius_pc, edge_pc).count(_ReachOnly(bounds))
    expected = count * profile_at(center, skeleton.shape, edge_pc, skeleton.expected_system_count_at_density_1).expected_systems
    log.normal(f"{count:,} sector(s) lie within {radius_pc:g} pc of {_format_address(address)}, outside the "
               f"galaxy's outline; at its profile they expect about {expected:,.3g} systems in all.")
    if count > LARGE_RUN_WARNING_THRESHOLD and args is not None and not getattr(args, "yes", False):
        run_common._refuse_or_warn(args, f"That is {count:,} sectors, a very large run (--yes goes ahead without asking).")


class EnumerationTooLarge(Exception):
    """PERF.36: a web request asked for a ball of more sectors than it may list (`run_common.WEB_ENUMERATION_LIMIT`);
    the message says how many and what to do."""


def _neighborhood_candidates(center, radius_pc, edge_pc, config, bounds, args=None, skeleton=None, ordered=True,
                             limit=None):
    """
    Every address within `radius_pc` of `center` that the neighbourhood
    takes (`_neighborhood_filter`: inside the galaxy's outline for a centre
    inside it, anywhere in the reach for one outside), the set of those
    already occupied, and how many addresses in the sphere were left out --
    the sphere is trimmed to the galaxy before anything is generated.
    A centre outside the outline reports the size of the request first
    (`_report_far_neighborhood`, once per run) when `args` and `skeleton` are given.

    The candidates come in the order the fill generates them (GEN.101,
    `fill_order.fill_order`): the centre's sector first, then outward in blocks of
    four sectors a side along a Hilbert curve; a caller that only wants the set
    passes `ordered=False` and gets the enumeration's order. With `args` given, the
    steps that are not between neighbours are counted in a line of the log, once
    per run.

    Raises:
        EnumerationTooLarge: The ball holds more than `limit` sectors (PERF.36: a web request is not to enumerate more than
            `run_common.WEB_ENUMERATION_LIMIT`); nothing was enumerated.
    """
    far = _is_far(center, edge_pc, bounds)
    keep = bounds.contains if far else bounds.qualifies
    if bounds:
        # Nothing in the galaxy lies farther from `center` than this, so a
        # larger radius only enumerates (and discards) empty space --
        # `--radius-pc 1e300` used to walk ~1e300 shells before trimming.
        farthest = bounds.reach_shell_index if far else bounds.outer_shell_index
        reach_pc = math.sqrt(sum(c * c for c in center)) + (farthest + 2) * edge_pc
        radius_pc = min(radius_pc, reach_pc)
    if limit is not None:
        held = Region(center, radius_pc, edge_pc).count(_ReachOnly(bounds) if far else bounds)
        if held > limit:
            raise EnumerationTooLarge(f"That ball holds {held:,} sectors, more than the {limit:,} a web request may list. "
                             f"Choose a smaller radius, or run the fill from the command line.")
    if far and skeleton is not None and not getattr(args, "_far_reported", False):
        if args is not None:
            args._far_reported = True
        _report_far_neighborhood(args, center, radius_pc, edge_pc, bounds, skeleton)
    candidates = []
    outside = 0
    for candidate in enumerate_sectors_within_radius(center, radius_pc, edge_pc):
        if keep(candidate[0], candidate[1]):
            candidates.append(candidate)
        else:
            outside += 1
    if ordered:
        candidates = fill_order.fill_order(candidates, center, edge_pc)
        if args is not None and not getattr(args, "_order_reported", False):
            args._order_reported = True
            _report_fill_order(candidates, edge_pc)
    conn = store.get_connection(config)
    try:
        occupied = store.get_occupied_addresses(conn, {c[0] for c in candidates})
    finally:
        conn.close()
    return candidates, occupied, outside


def _report_fill_order(candidates, edge_pc):
    """One line on how well the fill order walks (GEN.101): the sectors, the steps between neighbours and the jumps
    where the ball cuts the curve, with their mean and longest length in sector edges; nothing for a run of one."""
    steps, jumps, mean, longest = fill_order.step_statistics(candidates, edge_pc)
    if not steps:
        return
    log.normal(f"Fill order: {len(candidates):,} sector(s), {steps - jumps:,} of {steps:,} steps ({(steps - jumps) / steps:.1%}) "
               f"between neighbouring sectors; {jumps:,} jump(s) where the ball cuts the curve"
               + (f", {mean:.1f} sector edges on average and {longest:.1f} at most." if jumps else "."))


def _neighborhood_batch(args, candidates, occupied, batch_density, suffix="", skip=()):
    """
    The not-yet-generated sectors among `candidates`
    (`_neighborhood_candidates`), as `_submit_batch` items, in
    the order `candidates` came (the fill order of GEN.101: the centre first, then
    outward); `suffix` may hold `{distance}` (pc).

    Returns:
        tuple: `(batch, already_existed, skipped)`.
    """
    pending = _Pending.of_candidates(args, candidates, occupied, batch_density, suffix, skip)
    batch = list(pending.items())
    return batch, pending.existed, pending.skipped


def run_local_neighborhood(args, edge_pc, progress):
    """
    Local-neighborhood mode: generates every not-yet-generated
    sector within `args.radius_pc` parsecs of `args.center_sector`'s own
    stored center -- see `_BatchDensity.resolve` for how each sector's density is set.

    Args:
        args (argparse.Namespace): Parsed arguments;
            `args.center_sector`/`args.radius_pc` must not be `None`.
        edge_pc (float): The sector edge length, in parsecs (`_edge_pc`).
        progress (rich.progress.Progress): `run_galaxy`'s shared progress
            display (see `run_shell`). `run_random_start` also calls
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
            f"generated via 'planetgen galaxy' itself, not 'planetgen sector'). Use --shell to "
            f"generate placed sectors from scratch instead."
        )
        raise SystemExit(1)

    center = (
        center_position["center_x_pc"], center_position["center_y_pc"], center_position["center_z_pc"],
    )
    mysql_config = store.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    center_shell, center_band, _wedge = sector_address_at(center, edge_pc)
    _require_in_reach(batch_density.bounds, center_shell, center_band,
                      f"sector_id={args.center_sector} sits at {_format_address(sector_address_at(center, edge_pc))}")
    estimating = getattr(args, "estimate_only", False)   # the Generate page's web request (PERF.36)
    try:
        candidates, occupied, outside = _neighborhood_candidates(
            center, args.radius_pc, edge_pc, mysql_config, batch_density.bounds, args, batch_density.skeleton,
            ordered=not estimating, limit=run_common.WEB_ENUMERATION_LIMIT if estimating else None,
        )
    except EnumerationTooLarge as exc:
        log.error(str(exc))
        raise SystemExit(1) from exc
    region = _region_text(args)

    pending = _Pending.of_candidates(args, candidates, occupied, batch_density,
                                     f", {{distance:.2f}} pc from sector_id={args.center_sector}")
    batch, estimate_args, total = _plan_pending(pending)
    already_existed, skipped = pending.existed, pending.skipped
    run_common._check_estimate(args, estimate_args,
                    f"the sectors {region} of sector {args.center_sector}", progress, total=total)
    with run_common.sector_bar(progress, args, "Sectors (local neighborhood)", pending.count) as bar:
        _submit_batch(args, batch, f"Sectors ({region} of sector {args.center_sector})",
                      edge_pc, progress, bar, count=pending.count)
    generated = pending.count

    skip_note = f", {skipped} skipped (beyond the grid's reach)" if skipped else ""
    outside_note = f", {outside} beyond the galaxy's edge left out" if outside else ""
    log.normal(
        f"Generated {generated} new sector(s) {region} of sector_id={args.center_sector} "
        f"({len(candidates)} candidate sector(s) inside the galaxy{outside_note}, {already_existed} already "
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


class DatasetsMissing(RuntimeError):
    """GEN.590: a file the generator reads was never downloaded (`planetgen.datasets`), so a bulk generation was refused
    before writing anything; the message names the datasets."""


def require_datasets():
    """
    The second gate in front of every bulk generation (GEN.590, Boss's "Require the data"): the data files star composition
    reads (isochrones, white-dwarf sequences, neutron-star tables) must be in the cache folder that setup fills. Nothing is
    downloaded here (setup does that); the files are looked at and hashed once per process.

    Raises:
        DatasetsMissing: Naming each dataset that is missing and why.
    """
    missing = datasetCatalogue.missing()
    if missing:
        named = "; ".join(f"{located.dataset.key} ({located.dataset.fill}): {located.reason}" for located in missing)
        raise DatasetsMissing(
            f"the downloaded data is missing ({named}), so nothing was generated. Setup downloads it "
            f"('python3 -m planetgen.cli.datasets status' says which files are missing).")


def estimate_cell(density, systems, config=None):
    """
    The size and time of filling one cell (UX.113), from the server's own stored speeds and sizes
    (`generationStats.estimate`; the defaults before any run is recorded). Reads, writes nothing, and does not check the
    disk: the Generate page's estimate does that when a run is asked for.

    Args:
        density (float): The cell's `relative_density` (picks the speed bucket).
        systems (float): The systems the fill makes (`CellResolution.resolved_systems`).
        config (MySQLConfig, optional): The galaxy's database.

    Returns:
        dict: `systems`, `stars`, `bytes`, `seconds` and `measured` (whether the speed is this server's own).
    """
    config = config or store.DEFAULT_MYSQL_CONFIG
    try:
        result = run_common.generationStats.estimate([(density, systems)], run_common._stats_for_config(config),
                                                     config.database, workers=1)
    finally:
        run_common._RUN_STATS.clear()
    return {"systems": round(result.systems), "stars": round(result.stars), "bytes": result.bytes,
            "seconds": round(result.seconds, 1), "measured": result.measured}


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
                   placed in a galaxy, or lies beyond the grid's reach.
        EnumerationTooLarge: The radius holds more sectors than a web request may list (PERF.36).
        GenerationRefused: The database disk can't hold it (nothing was
                   written).
        RuntimeError: If the galaxy's skeleton has never been built.
        MathCheckFailed: If the math check failed (a real run only);
            nothing is written.
        DatasetsMissing: If a downloaded data file is missing (a real run only);
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
        require_datasets()

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
    center_shell, center_band, _wedge = sector_address_at(center, edge_pc)
    if not batch_density.bounds.contains(center_shell, center_band):
        raise ValueError(
            f"sector_id={center_sector_id} lies outside the grid: "
            f"{batch_density.bounds.describe_beyond_reach(center_shell, center_band)}."
        )
    candidates, occupied, outside = _neighborhood_candidates(
        center, radius_pc, edge_pc, config, batch_density.bounds,
        ordered=not estimate_only, limit=run_common.WEB_ENUMERATION_LIMIT,
    )

    args = _default_generation_args(config=config)
    # Density-driven from the skeleton (_BatchDensity), not the flat
    # num_systems=10 default.
    args.density = None
    args.num_systems = None
    args.workers = 1   # one sector at a time, in whichever process called this

    pending = _Pending.of_candidates(args, candidates, occupied, batch_density)
    batch, estimate_args, total = _plan_pending(pending)
    counts = {
        "already_existed": pending.existed,
        "skipped": pending.skipped,
        "candidates": len(candidates),
        "outside_galaxy": outside,
    }
    estimate = run_common._estimate_sectors(args, estimate_args, total)
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


def _whole_neighborhood_inside(bounds, position_pc, radius_pc, edge_pc):
    """Whether the centre and the six points `radius_pc` out along the axes
    all lie inside the outline: a cheap test that the neighborhood is not
    cut by the galaxy's edge (GEN.97)."""
    x, y, z = position_pc
    for dx, dy, dz in ((0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
        shell, band, _wedge = sector_address_at((x + dx * radius_pc, y + dy * radius_pc, z + dz * radius_pc), edge_pc)
        if not bounds.qualifies(shell, band):
            return False
    return True


def _touches_filled_space(args, position_pc, radius_pc, edge_pc, batch_density):
    """Whether any sector within `radius_pc` of `position_pc` is already generated (GEN.186)."""
    candidates, occupied, _outside = _neighborhood_candidates(
        position_pc, radius_pc, edge_pc, store.mysql_config_from_args(args), batch_density.bounds, ordered=False)
    return bool(occupied)


def _pick_random_start(args, edge_pc, batch_density, conn, radius_pc, taken, drawn=None):
    """
    One random start for random-start mode: an unoccupied address drawn
    inside the outline by volume (`GalaxyBounds.random_address`), retried up
    to `RANDOM_START_MAX_PLACEMENT_ATTEMPTS` times until it meets
    `--min-start-density`, and, for `--neighborhoods` (GEN.97), until its
    whole neighborhood lies inside the outline, it is at least twice the
    radius from every start in `taken`, and it passes the
    `--neighborhood-gamma` draw (accepted with probability
    `min(1, density) ** gamma`; gamma 0, the default, accepts every one).
    With `--avoid-filled-space` (GEN.186) a centre whose neighborhood
    touches an already generated sector is retried too; `drawn`, when
    given, is a dict that gets `centres` (examined) and `filled` (dropped
    for touching filled space) counted up.

    Returns:
        tuple or None: `(address, position_pc, sector_args)`; `None` when
            the attempt budget ran out.
    """
    many = args.neighborhoods > 1
    gamma = args.neighborhood_gamma
    for _ in range(program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS):
        address = batch_density.bounds.random_address(draw, max_shell=args.max_shell)
        if store.get_sector_id_at(conn, *address) is not None:
            continue
        position_pc = sector_position_pc(*address, edge_pc)
        if many and (
                any(detmath.dist(position_pc, other) < 2 * radius_pc for other in taken)
                or not _whole_neighborhood_inside(batch_density.bounds, position_pc, radius_pc, edge_pc)):
            continue
        if args.avoid_filled_space:
            if drawn is not None:
                drawn["centres"] = drawn.get("centres", 0) + 1
            if _touches_filled_space(args, position_pc, radius_pc, edge_pc, batch_density):
                if drawn is not None:
                    drawn["filled"] = drawn.get("filled", 0) + 1
                continue
        sector_args = batch_density.resolve(args, address, position_pc)
        if sector_args is None:
            continue
        if args.min_start_density is not None and sector_args.density < args.min_start_density:
            continue
        if gamma and many and draw.random() >= min(1.0, sector_args.density) ** gamma:
            continue
        return address, position_pc, sector_args
    return None


def run_random_start(args, edge_pc, progress):
    """
    Random-start mode (no `--shell`/`--center-sector` given): picks a
    random, not-yet-occupied sector address -- drawn only from
    inside the galaxy's stored outline (`GalaxyBounds.random_address`,
    every sector equally likely, optionally only out to `--max-shell`), so
    the start and the neighborhood around it are always in the galaxy --
    retried up to `program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS`
    times, generates it, then generates every not-yet-generated sector
    within `args.radius_pc` of it (default
    `program_constants.DEFAULT_GENERATE_RADIUS_PC`, 12 pc) via
    `run_local_neighborhood`. `--min-start-density` tightens the retry:
    an address below it is retried too.

    With `--neighborhoods N` (GEN.97) it makes N of them: all N starts are
    chosen first (see `_pick_random_start`), the estimate covers all of
    them, then each is generated in turn. Too few places for N is a
    warning when at least one was found.

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
    starts = []
    drawn = {}
    try:
        for _ in range(args.neighborhoods):
            start = _pick_random_start(args, edge_pc, batch_density, conn, radius_pc,
                                       [position for _address, position, _args in starts], drawn)
            if start is None:
                break
            starts.append(start)
    finally:
        conn.close()
    if args.avoid_filled_space:
        log.normal(f"Keeping away from filled space: {drawn.get('filled', 0):,} of {drawn.get('centres', 0):,} "
                   f"candidate centres touched sectors that are already generated and were dropped, "
                   f"{drawn.get('centres', 0) - drawn.get('filled', 0):,} qualified.")
    if not starts:
        density_note = (
            f", meeting --min-start-density {args.min_start_density} (try lowering it or --max-shell "
            f"-- most volume-weighted draws land in the galaxy's own sparser outskirts)"
            if args.min_start_density is not None else ""
        )
        within = f"within {args.max_shell} shells" if args.max_shell is not None else "in the galaxy"
        log.error(
            f"Could not find an unoccupied sector address {within} after {program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS} attempts{density_note} "
            f"-- this galaxy may already be almost entirely generated within that range, or that range "
            f"may hold too little real stellar density; try a larger --max-shell."
        )
        raise SystemExit(1)
    if len(starts) < args.neighborhoods:
        log.normal(f"Found room for only {len(starts)} of {args.neighborhoods} neighborhoods "
                   f"(each whole inside the galaxy, at least {2 * radius_pc:g} pc apart).")

    # The estimate covers every start and its whole neighborhood, before
    # any is written.
    estimate_args = []
    for address, position_pc, sector_args in starts:
        candidates, occupied, _outside = _neighborhood_candidates(
            position_pc, radius_pc, edge_pc, mysql_config, batch_density.bounds,
        )
        batch, _existed, _skipped = _neighborhood_batch(args, candidates, occupied, batch_density, skip={address})
        estimate_args += [sector_args] + [item[2] for item in batch]
    first_address = starts[0][0]
    what = (f"a random start at {_format_address(first_address)} and the sectors within {radius_pc:g} pc of it"
            if len(starts) == 1 else f"{len(starts)} random neighborhoods of {radius_pc:g} pc")
    run_common._check_estimate(args, estimate_args, what, progress)

    for address, position_pc, sector_args in starts:
        stages.enter("start")
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
        stages.enter("neighborhood")
        run_local_neighborhood(args, edge_pc, progress)


def run_single_sector(args, edge_pc, progress):
    """
    Single-address mode: generates exactly the one sector at `(--shell I,
    --band J, --wedge K)` via the same `ensure_sector_generated` the
    galaxy map's "recalculate on visit" path uses -- the direct route from
    an address copied out of the interactive 3D Galaxy Map. With
    `--radius-pc`, then generates that sector's neighborhood
    (`run_local_neighborhood`), the Sector Map's "Generate neighborhood".

    Raises:
        SystemExit: If the address is out of range for the band, or is beyond
                   the grid's reach.
    """
    address = (args.shell, args.band, args.wedge)
    try:
        total_wedges = spherical.wedge_count(args.shell, args.band)
    except ValueError as error:
        log.error(f"{_format_address(address)} -- {error}.")
        raise SystemExit(1) from error
    if not (0 <= args.wedge < total_wedges):
        log.error(
            f"--wedge {args.wedge} is out of range for shell {args.shell} band {args.band} (holds {total_wedges} "
            f"wedges, 0..{total_wedges - 1})."
        )
        raise SystemExit(1)

    mysql_config = store.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    # Named outright, so generated from its profile anywhere in the reach (GEN.235).
    _require_in_reach(batch_density.bounds, args.shell, args.band, _format_address(address))
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
        sector_args = batch_density.resolve(single, address, position_pc)
    sectors = [sector_args] if sector_args is not None else []
    what = _format_address(address)
    if args.radius_pc is not None:
        candidates, occupied, _outside = _neighborhood_candidates(
            position_pc, args.radius_pc, edge_pc, mysql_config, batch_density.bounds, args, batch_density.skeleton,
        )
        batch, _existed, _skipped = _neighborhood_batch(args, candidates, occupied, batch_density, skip={address})
        sectors += [item[2] for item in batch]
        what += f" and the sectors {_region_text(args)} of it"
    run_common._check_estimate(args, sectors, what, progress)
    # The backfill waits for the end of the run (backfill_after_run), with
    # its own bar, instead of stalling this one at 0 of 1 (PERF.28).
    with run_common.sector_bar(progress, args, f"Sector ({_format_address(address)})", 1) as bar:
        result = ensure_sector_generated(*address, config=mysql_config, backfill=False,
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


def _generate_addresses(args, addresses, shells, total, what, edge_pc, progress, batch_density):
    """
    Generates every not-yet-generated address of `addresses` (a callable giving a new iterable of
    `(shell, band, wedge)` each time it is called, walked lazily, in order) up to `args.limit` when
    set -- the loop behind the shell, span and course modes; see
    `_BatchDensity.resolve` for how each sector's density is set. A run of more than
    `generationStats.SAMPLE_ABOVE` sectors is estimated from a sample and its sectors are resolved as they are
    submitted, not held (`_Pending`, PERF.36).

    Args:
        shells (iterable): The shells `addresses` lie in, to look up what is
            already generated there.
        total (int): How many addresses there are, for the closing line.
    """
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        occupied = set(store.get_occupied_addresses(conn, set(shells)))
    finally:
        conn.close()

    limit = getattr(args, "limit", None)
    pending = _Pending.of_addresses(args, addresses, occupied, batch_density, edge_pc, limit=limit)
    batch, estimate_args, sampled_total = _plan_pending(pending)
    skipped, existed = pending.skipped, pending.existed

    run_common._check_estimate(args, estimate_args, what, progress, total=sampled_total)
    with run_common.sector_bar(progress, args, f"Sectors ({what})", pending.count) as bar:
        _submit_batch(args, batch, f"Sectors ({what})", edge_pc, progress, bar, count=pending.count)
    generated = pending.count
    skip_note = f", {skipped} skipped (beyond the grid's reach)" if skipped else ""
    already = f"{total - generated - skipped} already existed" if limit is None else f"{existed} already existed"
    log.normal(f"Generated {generated} new sector(s) in {what} ({total} total, {already}{skip_note}).")


def run_shell(args, edge_pc, progress):
    """
    Shell mode: every not-yet-generated sector of shell `--shell I` the
    galaxy's outline holds -- or, with `--band J`, of that one band (named
    outright, so generated even outside the outline, GEN.81). When the shell
    holds more than `LARGE_RUN_WARNING_THRESHOLD` sectors it needs `--limit`
    or `--yes`.

    Raises:
        SystemExit: If the shell (or the named band) is beyond the grid's reach,
                   or the run is too large and neither `--limit` nor `--yes`
                   was given.
    """
    shell_index = args.shell
    batch_density = _BatchDensity(store.mysql_config_from_args(args))
    bounds = batch_density.bounds
    if args.band is not None:
        band = args.band
        what = f"shell {shell_index} band {band}"
        _require_in_reach(bounds, shell_index, band, what)
        total = spherical.wedge_count(shell_index, band)
        def addresses():
            return ((shell_index, band, wedge) for wedge in range(total))
    else:
        what = f"shell {shell_index}"
        _require_in_reach(bounds, shell_index, 0, what)
        if shell_index not in bounds.high_band:
            log.normal(f"No band of the outline reaches shell {shell_index}: nothing to generate "
                       f"(name a band with --band to generate its sectors from their profiles).")
            return
        span = Span(shells=(shell_index, shell_index))
        total = span.count(bounds)
        def addresses():
            return span.addresses(bounds)
    if total > LARGE_RUN_WARNING_THRESHOLD and args.limit is None and not args.yes:
        run_common._refuse_or_warn(args, f"{what.capitalize()} holds {total:,} sectors, a very large run "
                                         f"(--limit N generates only the first N not-yet-generated sectors).")
    _generate_addresses(args, addresses, [shell_index], total, what, edge_pc, progress, batch_density)


def run_span(args, edge_pc, progress):
    """
    Span mode (ADM.29): every not-yet-generated sector in `args.span`, a
    range of shells, bands and wedges, inside the galaxy's outline. The
    count comes from the outline's stored cell counts, so the size warning is
    instant. Needs `--limit` or `--yes` when it is more than
    `LARGE_RUN_WARNING_THRESHOLD` sectors.
    """
    span = args.span
    batch_density = _BatchDensity(store.mysql_config_from_args(args))
    bounds = batch_density.bounds
    total = span.count(bounds)
    what = span.describe()
    if total == 0:
        log.normal(f"{what.capitalize()} holds no sector inside the galaxy's outline: nothing to generate.")
        return
    if total > LARGE_RUN_WARNING_THRESHOLD and args.limit is None and not args.yes:
        run_common._refuse_or_warn(args, f"{what.capitalize()} holds {total:,} sectors, a very large run "
                                         f"(--limit N generates only the first N).")
    _generate_addresses(args, lambda: span.addresses(bounds), span._shells(bounds), total, what, edge_pc, progress,
                        batch_density)


def course_cells(config, from_ref, to_ref, border=False):
    """
    The uncharted cells that block the NAV course between two systems (NAV.48), the plan
    `db.query.nav_chart_plan` works out from the stored galaxy as it is now.

    Args:
        config (MySQLConfig): The database.
        from_ref (str): The first object: a system's object ID, or `<kind>:<ID>`.
        to_ref (str): The last object, same forms.
        border (bool): Also the uncharted cells sharing a face with those on the course.

    Returns:
        dict: `nav_chart_plan`'s answer.

    Raises:
        ValueError: A system that does not exist, or a course NAV cannot plot.
    """
    from planetgen.api import ids
    from planetgen.db import query

    conn = store.get_connection(config)
    try:
        refs = []
        for text in (from_ref, to_ref):
            resolved = ids.resolve_ref(conn, text)
            if resolved is None:
                raise ValueError(f"nothing has the ID {text}")
            refs.append(resolved)
        try:
            return query.nav_chart_plan(conn, refs[0], refs[1], border=border)
        except query.NavUnavailable as exc:
            raise ValueError(str(exc)) from exc
    finally:
        conn.close()


def run_course(args, edge_pc, progress):
    """
    Course mode (NAV.48): every uncharted cell on the hops through unknown space of the NAV course
    between two objects (`course_cells`), inside the galaxy's outline. Past
    `program_constants.NAV_CHART_CONFIRM_SECTORS` sectors it needs `--yes`; the size and disk checks are
    the ones every run has. The bright-star backfill and the settle that follow every run cover the
    new sectors.
    """
    config = store.mysql_config_from_args(args)
    try:
        plan = course_cells(config, args.course[0], args.course[1], border=getattr(args, "course_border", False))
    except ValueError as exc:
        log.error(str(exc))
        raise SystemExit(1) from exc
    cells = plan["cells"]
    what = f"the uncharted sectors on the course {args.course[0]} to {args.course[1]}"
    if not cells:
        why = ("the course crosses no uncharted sector" if not plan["unknown_hops"]
               else "the uncharted sectors on the course all lie outside the galaxy")
        log.normal(f"Nothing to chart: {why}.")
        return
    bypass = plan.get("bypass")
    if bypass and bypass["found"]:
        log.normal(f"Note: a route through charted space exists ({bypass['distance_ly']:,.1f} ly against "
                   f"{plan['route_distance_ly']:,.1f} ly); charting is optional.")
    if len(cells) > program_constants.NAV_CHART_CONFIRM_SECTORS and not args.yes:
        log.error(f"Charting this course takes {len(cells):,} sectors, more than the {program_constants.NAV_CHART_CONFIRM_SECTORS:,} "
                  f"that need confirming. Run it again with --yes.")
        raise SystemExit(1)
    if plan["outside_galaxy"]:
        log.normal(f"{plan['outside_galaxy']} uncharted cell(s) lie outside the galaxy and are left out.")
    _generate_addresses(args, lambda: cells, {cell[0] for cell in cells}, len(cells), what, edge_pc, progress,
                        _BatchDensity(config))


def run_galaxy(args):
    """
    Dispatches to span, course, single-address, shell,
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
        log.error("The stored plan has no shells at all (nothing in this galaxy shape expects a star "
                  "per sector) -- re-run 'planetgen plan' with a different shape.")
        raise SystemExit(1)

    estimate_only = getattr(args, "estimate_only", False)
    if not estimate_only:
        stages.begin(stages.galaxy_stages(args), args, "galaxy")
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
        run_plan.scatter_phenomena(args)    # the last pass of GEN.185: the other phenomena, after the stars
    backfill_after_run(args, edge_pc, started_at)
    run_population.run_population_after(args)
    settle_after_run(args, started_at)
    stages.finish()
    run_common._finish_stats(args)    # the steps after the sectors recorded their speeds too


def _random_start_mode(args):
    """Whether `_run_galaxy_mode` will pick a random start (no mode argument was given)."""
    return stages.galaxy_stages(args)[0].key == "start"


def _run_galaxy_mode(args, edge_pc, progress):
    """`run_galaxy`'s dispatch to the mode its arguments ask for."""
    if not _random_start_mode(args):
        stages.enter("sectors")
    if getattr(args, "span", None) is not None:
        run_span(args, edge_pc, progress)
    elif getattr(args, "course", None) is not None:
        run_course(args, edge_pc, progress)
    elif args.wedge is not None:
        run_single_sector(args, edge_pc, progress)
    elif args.shell is not None:
        run_shell(args, edge_pc, progress)
    elif args.center_sector is not None:
        run_local_neighborhood(args, edge_pc, progress)
    else:
        run_random_start(args, edge_pc, progress)
