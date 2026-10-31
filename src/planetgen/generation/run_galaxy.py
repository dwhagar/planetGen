# planetgen/generation/run_galaxy.py

"""
Run Galaxy
==========

The `galaxy` command: many sectors placed as one galaxy, each at its
shell/band/wedge address in the spherical grid, with its systems inside
the sector's own cell. Modes: `--shell` (a shell, or with `--band` one band
of it, or with `--wedge` too one sector), `--shells/--bands/--wedges` (a
span), `--block` (a block of the Galaxy Map's ladder), `--center-sector`/`--radius-pc` (local neighborhood), `--course`,
or none (random start).
The sector logic is for runs only: no page or API route generates a sector while it answers (PERF.149), so a visit to an
uncharted cell shows its profile and generates nothing. `regenerate_sector_cell` makes one sector as a run, and an admin's
visit queuing a one-sector run is PERF.204.
"""

import argparse
import collections
import copy
import math
import os
import time

import pymysql

from planetgen.datasets import catalogue as datasetCatalogue, fetch
from planetgen.db import block_lights, block_stats, nebula_blocks, query, region, region_stats, scatter_systems, sector_paths, store
from planetgen.generation import bright_stars as brightStars
from planetgen.galaxy import density, drill, named_regions, region_run, resolution as galaxyResolution, seed as galaxySeed, version_check
from planetgen.queue import commands, wrapper
from planetgen import tuning as program_constants
from planetgen.util import log
from planetgen.galaxy.profile import profile_at, sector_profile
from planetgen.galaxy.region import Region
from planetgen.galaxy import spherical
from planetgen.galaxy.frame import galactic_radius_pc
from planetgen.galaxy.span import Span, parse_range
from planetgen.galaxy.spherical import (
    SectorCell, enumerate_sectors_within_radius, neighbor_addresses,
    provisional_sector_designation, sector_address_at, sector_position_pc,
)
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.generation import fill_order, followed, object_timing
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
    a single-address run uses, applied across a whole run instead of every
    sector sharing one flat CLI value.

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

    def located(self):
        """Each address to generate with its suffix (`place`'s second answer), in order, its sector unresolved."""
        for address, extra in self._eligible():
            yield address, self._place(address, extra)[1]

    def addresses(self):
        """The addresses to generate, in order, unresolved (a run's cut needs no more than where its sectors are)."""
        for address, _extra in self._eligible():
            yield address

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
        # The resolution reads the reference metallicity, which depends on the galaxy's own disk scale length: bind the
        # galaxy as the fill will, so the cell is aimed at the number its stored stats then record.
        with galaxySeed.seeded(skeleton.galaxy_seed, "sector", address,
                               disk_scale_length_pc=skeleton.shape.disk_scale_length_pc, shape=skeleton.shape):
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
        if not resolution.capped:
            _take_to_the_deepest_floor(store.mysql_config_from_args(args), skeleton, address, limits, backfill_seed)
            conn.commit()   # the draw ran on a connection of its own: this one's snapshot must see its stars
        level = store.bright_star_fill_level(conn, *address)
        mass_limit = store.bright_star_fill_mass_limit(conn, *address)
        emission_hosts = _emission_hosts(conn, skeleton, address, profile)
        if level is None and mass_limit is None:
            return brightStars.FillContext(profile, phenomenon_rows=phenomena, below_cut=below_cut,
                                           resolution=resolution, emission_hosts=emission_hosts)
        rows = store.bright_stars_for_sector(conn, *address)
        # GEN.315: a row with members is built from the tree it stored
        member_trees = scatter_systems.read_stored_trees(conn, "bright_stars", [r["id"] for r in rows if r.get("member_count")])
        orbit_epoch_unix = store.get_orbit_epoch_unix(conn)
    finally:
        conn.close()
    return brightStars.FillContext(profile, rows, min_luminosity_sol=level,
                                   phenomenon_rows=phenomena, below_cut=below_cut, star_mass_limit_sol=mass_limit,
                                   resolution=resolution, emission_hosts=emission_hosts, member_trees=member_trees,
                                   orbit_epoch_unix=orbit_epoch_unix)


def _take_to_the_deepest_floor(config, skeleton, address, limits, seed):
    """
    GEN.633: before a sector is filled, takes its pre-placed bright stars down to the deepest floor the mass backfill
    ever uses (the first of `BRIGHT_STAR_BACKFILL_RING_MASSES_SOL`), whatever the backfill has reached so far. The
    backfill skips a sector once it is filled, so without this a sector the backfill reached first held the dimmer
    stars of its ring and one generated before the backfill reached it did not, and its contents depended on the order
    the sectors were generated in. A sector's draw is seeded from the scatter's seed, its address and fixed mass
    bands, so the stars drawn here are the ones the backfill would have drawn, and a sector the backfill already took
    to a heavier floor is topped up with only the masses it lacks. A sector at or below the floor (a band whose mass
    limit already reaches it) is left alone.

    The draw runs on a connection of its own at READ COMMITTED and starts over after a deadlock or lock wait timeout
    (`store.RETRYABLE_ERRORS`): the fills of several workers write `bright_stars` rows at once.
    """
    floor = float(min(program_constants.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL))
    mass_limit = limits.mass_limit_at(address[0], address[1])
    if mass_limit is not None and mass_limit <= floor:
        return
    for attempt in range(1, store.SECTOR_SAVE_ATTEMPTS + 1):
        conn = store.get_connection(config)
        try:
            conn.execute("SET TRANSACTION ISOLATION LEVEL READ COMMITTED")
            run_plan._draw_sector_masses(conn, skeleton, {tuple(address): floor}, limits, seed)
            return
        except pymysql.err.OperationalError as exc:
            if not exc.args or exc.args[0] not in store.RETRYABLE_ERRORS or attempt == store.SECTOR_SAVE_ATTEMPTS:
                raise
            conn.rollback()
            log.debug(f"The first fill of {tuple(address)} hit MySQL error {exc.args[0]}; trying again "
                      f"({attempt}/{store.SECTOR_SAVE_ATTEMPTS - 1}).")
            time.sleep((0.05 + os.getpid() % 8 * 0.03) * attempt)   # not a draw: the process id spreads the workers
        finally:
            conn.close()


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
        todo = backfill_targets(conn, skeleton, limits, centers_pc)
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


def backfill_targets(conn, skeleton, limits, centers_pc):
    """
    The sectors a bright-star backfill around `centers_pc` has stars to draw in, found once from the stored galaxy
    (PERF.108): the sectors out from the centers' sectors (`backfill_step_targets`) that a band's mass limit does not
    already cover (GEN.333), that are not filled and that are not already taken that deep.

    Returns:
        dict: `(shell, band, wedge)` -> the lightest initial mass to place there, solar masses.
    """
    bounds = store.get_galaxy_bounds(conn)
    centers = {sector_address_at(center, skeleton.edge_pc) for center in centers_pc}
    targets = backfill_step_targets(centers, bounds)
    # GEN.333: a sector backfills the masses under its own band's mass limit.
    targets = {address: mass for address, mass in targets.items()
               if limits.mass_limit_at(address[0], address[1]) is None
               or mass < limits.mass_limit_at(address[0], address[1])}
    held = store.sector_bright_masses(conn, targets)
    filled = store.get_occupied_addresses(conn, {address[0] for address in targets})
    return {address: mass for address, mass in targets.items()
            if address not in filled and mass < held.get(address, math.inf)}


BACKFILL_UNIT_CELLS = 10_000
"""int: The fewest sectors a backfill unit draws (PERF.108): the draw is about 0.06 ms a sector, so a smaller unit is
mostly overhead."""

BACKFILL_UNITS_MAX = 2000
"""int: The most units a backfill is cut into, however many sectors it has."""


def expand_backfill_bright_stars(params):
    """
    The `scatter-backfill` command's stage (PERF.108): the sectors to draw (`backfill_targets` around the sectors
    created since `params["since"]`, or around `params["center_pc"]` for an API job's one requested sector) are found once, here, and cut in order into static chunks of at least
    `BACKFILL_UNIT_CELLS` sectors (and at most `BACKFILL_UNITS_MAX` chunks), the same at any worker count. A stage with
    nothing to draw is one empty unit, so the stage still stores its row.
    """
    conn = store.get_connection(_mysql_config_of(params))
    try:
        skeleton = store.get_galaxy_shape(conn)
        todo = {}
        skip = None
        if skeleton is not None:
            limits, _seed = run_plan._scatter_limits_and_seed(conn, skeleton)
            centers = ([tuple(params["center_pc"])] if params.get("center_pc") is not None
                       else store.sector_centers_since(conn, params["since"]))
            todo = backfill_targets(conn, skeleton, limits, centers)
            skip = None if centers else "the run generated no sector"
    finally:
        conn.close()
    addresses = sorted(todo)
    size = max(BACKFILL_UNIT_CELLS, -(-len(addresses) // BACKFILL_UNITS_MAX))
    units = [({"database": params["database"], "items": [[*address, todo[address]] for address in addresses[start:start + size]]},
              float(len(addresses[start:start + size]))) for start in range(0, len(addresses), size)]
    return commands.Expansion(units or [({"database": params["database"], "items": []}, 1.0)],
                              {"sectors": len(addresses), **({} if skip is None else {"skip": skip})})


def backfill_bright_stars_unit(unit_args, report):
    """
    The `scatter-backfill` command's unit (PERF.108): draws the chunk's sectors down to their target masses
    (`run_plan._draw_sector_masses`, a few hundred a transaction, each sector rechecked under its lock), every sector
    from its own seed, so the stars are the same whichever sector reached it first and at any worker count.

    Returns:
        dict: `sectors` (drawn), `stars` (placed), `visited` (the chunk's sectors) and `shells`, the stars placed per
            shell and kind (`shell index as text` -> `[[label, plural, count], ...]`, for the starter's lines).
    """
    items = unit_args["items"]
    targets = {tuple(item[:3]): item[3] for item in items}
    conn = store.get_connection(_mysql_config_of(unit_args))
    shells = {}
    visited = [0]

    def on_sector():
        visited[0] += 1
        report(visited[0], len(items))

    try:
        skeleton = store.get_galaxy_shape(conn)
        sectors = stars = 0
        if skeleton is not None and targets:
            limits, seed = run_plan._scatter_limits_and_seed(conn, skeleton)
            sectors, stars = run_plan._draw_sector_masses(conn, skeleton, targets, limits, seed, on_sector=on_sector,
                                                          shells=shells)
    except BaseException:
        conn.rollback()
        raise
    finally:
        conn.close()
    return {"sectors": sectors, "stars": stars, "visited": len(items),
            "shells": {str(shell): [[*label, count] for label, count in types.items()]
                       for shell, types in shells.items()}}


class BackfillFollower(followed.Follower):
    """
    The bright-star backfill of a run on the queue (PERF.108): the wrapper selects the `scatter-backfill` command from
    this run's steps, finds the sectors to draw once and cuts them into units (`expand_backfill_bright_stars`), this
    process follows the units with a bar of the sectors they visit (PERF.28) and says how many stars each shell was
    given, and the unit that completes the stage stores the stage's row.

    Args:
        started_at (datetime): The database's clock when the run started (or the deferred clock of the sectors stage,
            `wrapper.opened_by`): the backfill goes out from the sectors created since.
        center_pc (tuple, optional): Backfill around this one point instead (an API job's requested sector, GEN.30).

    Attributes:
        summary (dict): `sectors` (drawn) and `stars` (placed), summed over the units.
    """

    key = "backfill"
    label = "Bright-star backfill"

    def __init__(self, args, started_at, center_pc=None):
        super().__init__(args)
        self.started_at, self.center_pc = started_at, center_pc
        self.summary = {"sectors": 0, "stars": 0}
        self.shells = {}
        self.bar = None

    def params(self):
        return {"since": since_param(self.started_at),
                **({} if self.center_pc is None else {"center_pc": list(self.center_pc)})}

    def on_plan(self, plan):
        if plan["detail"].get("skip"):
            log.normal(f"Nothing to backfill: {plan['detail']['skip']}.")
            return
        log.normal("Backfilling the bright stars out from the edge of the run...")
        self.bar = self.hold(steps.Step("Bright-star backfill (sectors)", "backfill", 1, args=self.args,
                                        progress=self.progress, workers=self.workers))
        self.bar.update(total=max(plan["detail"]["sectors"], 1))

    def on_result(self, unit, result):
        self.summary["sectors"] += result["sectors"]
        self.summary["stars"] += result["stars"]
        for shell, types in result["shells"].items():
            counts = self.shells.setdefault(int(shell), collections.Counter())
            for label, plural, count in types:
                counts[(label, plural)] += count
        if self.bar is not None:
            self.bar.update(advance=result["visited"])

    def on_done(self, conclusion):
        if self.bar is None:
            return
        self._held.close()
        run_plan.log_shells("Bright-star backfill", self.shells)
        stages.note(sectors=self.summary["sectors"], objects=self.summary["stars"])
        log.normal(f"Backfilled {self.summary['stars']:,} bright stars in {self.summary['sectors']:,} sectors.")


def _backfill_on_queue(args, started_at, center_pc=None):
    """The backfill as a run of its own (`BackfillFollower`). Returns its `summary`: `sectors` (drawn) and `stars`
    (placed)."""
    return followed.run_planned(args, planned_backfill(args, started_at, center_pc))


def planned_backfill(args, started_at, center_pc=None):
    """The bright-star backfill as a `followed.Planned`; its result is the `summary`."""
    follower = BackfillFollower(args, started_at, center_pc)
    return followed.Planned([follower], "Bright-star backfill", lambda: follower.summary)


def since_param(started_at):
    """The `since` parameter of a stage that works on the sectors a run saved: the clock as text, or (PERF.213) the deferred
    parameter `wrapper.opened_by` that a run takes the clock of the sectors stage from when the stage opens."""
    return started_at if wrapper.is_deferred(started_at) else started_at.isoformat(sep=" ")


def _sectors_since(params):
    """The ids of the galaxy sectors a run created: those saved at or after `params["since"]` (the database's clock when
    the run started, as text), the set a unit of the link or the settle stage works on; or, for a run that names them
    (`params["sector_ids"]`), exactly those."""
    if params.get("sector_ids") is not None:      # an API job names the sectors it made or changed (PERF.119)
        return list(params["sector_ids"])
    conn = store.get_connection(_mysql_config_of(params))
    try:
        return sector_paths.sector_ids_since(conn, params["since"])
    finally:
        conn.close()


def _mysql_config_of(params):
    """The `MySQLConfig` a command's params name (`wrapper.database_args`)."""
    return run_plan._mysql_config_of(params)


SAVED_NO_SECTOR = "no sector was saved since the run started"
"""str: The reason a stage that works on the run's sectors gives for having nothing to do."""

LINK_UNIT_SECTORS = 500
"""int: Sectors in a unit of the `link-neighbours` stage (PERF.113). A sector links in about 70 ms at the default galaxy size,
so a unit is about 35 s and a long link is many short units (PERF.155)."""

LINK_UNITS_MAX = 1000
"""int: Most units the `link-neighbours` stage is cut into; a larger stage has bigger units."""


def expand_link_neighbours(params):
    """
    The `link-neighbours` command's stage (PERF.107, PERF.113): the sectors saved since `params["since"]` (or named in
    `params["sector_ids"]`) are found once, here, put in neighbourhood order (`store.link_order`) and cut into static
    chunks of `LINK_UNIT_SECTORS` (at most `LINK_UNITS_MAX` chunks), the same at any worker count. Each unit links its
    chunk and leaves the run's other sectors alone, as the batches of one call do: a sector's lists are computed whole from
    the placed systems, and a system merged into the lists of the sectors around is the same whichever batch merges it, so
    the stored links do not depend on the order the units run in. A stage with nothing to link is one empty unit.
    The units take the neighbour lock one batch at a time (`store.link_sector_neighbors`), so they queue for it rather than
    run side by side; the cut bounds a unit's length, it does not speed the stage up.
    """
    created = _sectors_since(params)
    ordered = store.link_order(_mysql_config_of(params), created)
    size = max(LINK_UNIT_SECTORS, -(-len(ordered) // LINK_UNITS_MAX))
    run = {key: params[key] for key in ("since", "sector_ids") if params.get(key) is not None}
    units = [({"database": params["database"], "chunk": ordered[start:start + size], **run},
              float(len(ordered[start:start + size]) * store.LINK_PHASES)) for start in range(0, len(ordered), size)]
    return commands.Expansion(units or [({"database": params["database"], "chunk": [], **run}, 1.0)],
                              {"sectors": len(ordered), **({} if ordered else {"skip": SAVED_NO_SECTOR})})


def link_neighbours_unit(params, report):
    """
    The `link-neighbours` command's unit (PERF.45, PERF.107, PERF.113): links a chunk of the sectors a `galaxy` run created
    to their neighbours (`store.link_sector_neighbors`): containment, nearest systems and quadrants, and the new systems
    merged into the lists of the sectors around. Done once for the whole run, in place of one locked step per saved sector,
    so the workers no longer queue for it; the stored result is the same. The run's other sectors (found again from
    `params["since"]` or `params["sector_ids"]`) are left to their own units.

    Returns:
        dict: `linked`, how many sectors were linked.
    """
    def on_progress(done, total, step=""):
        report(done, total)

    return {"linked": store.link_sector_neighbors(_mysql_config_of(params), params["chunk"], on_progress,
                                                  skip=_sectors_since(params))}


class LinkFollower(followed.Follower):
    """
    The link stage of a run on the new queue (PERF.113): the wrapper selects the `link-neighbours` command from this run's
    steps, finds the sectors to link once and cuts them into units (`expand_link_neighbours`), this process follows the
    units with a bar of the sectors they link, and the unit that completes the stage stores the stage's row.

    Attributes:
        linked (int): How many sectors were linked, summed over the units.
    """

    key = "link"
    label = "Neighbours"

    def __init__(self, args, started_at):
        super().__init__(args)
        self.started_at = started_at
        self.linked = 0
        self.bar = None

    def params(self):
        return {"since": since_param(self.started_at)}

    def on_plan(self, plan):
        if plan["detail"].get("skip"):
            log.normal(f"Nothing to link: {plan['detail']['skip']}.")
            return
        stages.note(sectors=plan["detail"]["sectors"])
        log.normal("Linking the new sectors to their neighbours...")
        self.bar = self.hold(steps.Step("Neighbours", "link", 1, args=self.args, progress=self.progress,
                                        workers=self.workers))
        self.bar.update(total=max(plan["detail"]["sectors"], 1))

    def on_result(self, unit, result):
        self.linked += result["linked"]
        if self.bar is not None:
            self.bar.update(advance=result["linked"])

    def on_done(self, conclusion):
        if self.bar is not None:
            log.normal(f"Linked {self.linked:,} sectors to their neighbours.")


def _link_on_queue(args, started_at):
    """The link stage as a run of its own (`LinkFollower`). Returns how many sectors were linked."""
    return followed.run_planned(args, planned_link(args, started_at))


def planned_link(args, started_at):
    """The link stage as a `followed.Planned`; its result is how many sectors were linked."""
    follower = LinkFollower(args, started_at)
    return followed.Planned([follower], "Neighbours", lambda: follower.linked)


def expand_refresh_run_blocks(params):
    """The `refresh-run-blocks` command's stage: one unit."""
    return commands.Expansion([(params, 1.0)], {})


def refresh_run_blocks_unit(params, report):
    """
    The `refresh-run-blocks` command's unit (MAP.172, PERF.172): brings the Galaxy Map's block totals up to what a `galaxy`
    run stored, once at the end of the run so its workers never queue on the same top-level blocks: the nebulae the
    `nebula_field` table does not hold (emission and reflection nebulae around a star, a class P or Q cloud;
    `db.nebula_blocks.catch_up`), then the generated totals of the sectors saved since the run started
    (`db.block_stats.refresh_sectors`, which a sector save does not do) and the light of the pre-placed objects, which a saved
    sector builds into its own sums (`db.block_lights.refresh_sectors`). A failure only warns: `planetgen fill-nebula-blocks` and
    `planetgen fill-block-totals` make the totals again.

    Returns:
        dict: `nebulae` (added), `sectors` (brought in) and `warnings` (the failures, as text).
    """
    config = _mysql_config_of(params)
    nebulae, sectors, warnings = 0, 0, []
    try:
        conn = store.get_connection(config)
        try:
            nebulae = nebula_blocks.catch_up(conn, query.galaxy_edge_pc(conn)) or 0
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- the nebulae are saved; the totals can be made again
        warnings.append(f"adding the new nebulae to the block totals failed: {exc}")
    report(1, 2)
    try:
        conn = store.get_connection(config)
        try:
            addresses = sector_paths.sector_addresses_since(conn, params["since"])
            if addresses:
                block_stats.refresh_sectors(conn, addresses)
                block_lights.refresh_sectors(conn, addresses)
                conn.commit()
            sectors = len(addresses)
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- the sectors are saved; the totals can be made again
        warnings.append(f"adding the new sectors to the block totals failed: {exc}")
    return {"nebulae": nebulae, "sectors": sectors, "warnings": warnings}


class RunBlocksFollower(followed.Follower):
    """The Galaxy Map's block totals brought up to date at the end of a run (`refresh_run_blocks_unit`)."""

    key = "run-blocks"
    label = "Update the map's block totals"

    def __init__(self, args, started_at):
        super().__init__(args)
        self.started_at = started_at

    def params(self):
        return {"since": since_param(self.started_at)}

    def on_result(self, unit, result):
        if result["nebulae"]:
            log.normal(f"Added {result['nebulae']:,} nebulae to the block totals.")
        for warning in result["warnings"]:
            log.normal(f"Warning: {warning}")


SETTLE_UNIT_SECTORS = 100
"""int: Sectors in a unit of the `settle-sector-paths` stage (PERF.115). A sector takes about 0.18 s to settle at the default
galaxy size, so a unit is about 18 s and a long settle is many short units (PERF.155)."""

SETTLE_UNITS_MAX = 1000
"""int: Most units the `settle-sector-paths` stage is cut into; a larger stage has bigger units."""


def expand_save_sector_paths(params):
    """
    The `settle-sector-paths` command's stage (PERF.107, PERF.115): the sectors to settle (those saved since `params["since"]`,
    or `params["sector_ids"]`, and with `params["expand"]` (the default) the sectors around them,
    `sector_paths.sectors_to_settle`) are found once, here, and cut in order into static chunks of `SETTLE_UNIT_SECTORS`
    (at most `SETTLE_UNITS_MAX` chunks), the same at any worker count. Each sector's paths depend only on the final
    sectors around it and are saved with the sector, so the chunks cannot collide. A stage with nothing to settle is one
    empty unit, so the stage still stores its row.
    """
    created = _sectors_since(params)
    if params.get("expand", True):
        conn = store.get_connection(_mysql_config_of(params))
        try:
            sectors = sector_paths.sectors_to_settle(conn, created)
        finally:
            conn.close()
    else:
        sectors = sorted(set(created))
    size = max(SETTLE_UNIT_SECTORS, -(-len(sectors) // SETTLE_UNITS_MAX))
    # One unit may spread a crowded sector's bodies over processes; a stage of many units already fills the workers.
    workers = params["workers"] if len(sectors) <= size else 1
    units = [({"database": params["database"], "sector_ids": sectors[start:start + size], "workers": workers},
              float(len(sectors[start:start + size]))) for start in range(0, len(sectors), size)]
    return commands.Expansion(units or [({"database": params["database"], "sector_ids": [], "workers": 1}, 1.0)],
                              {"sectors": len(sectors), **({"skip": SAVED_NO_SECTOR} if not created and params.get("sector_ids") is None else {})})


def save_sector_paths_unit(params, report):
    """
    The `settle-sector-paths` command's unit (GEN.126, PERF.107, PERF.115): saves the path of every body in the chunk of
    sectors it was given (`sector_paths.settle_sectors`, `expand=False`: the stage's expansion already added the sectors
    around the run's), so each path sees the final set of neighbours whatever order the workers filled sectors in. A
    crowded sector's bodies are integrated in `params["workers"]` processes.

    Returns:
        dict: `saved`, how many paths were saved, and `visited`, the chunk's sectors.
    """
    conn = store.get_connection(_mysql_config_of(params))
    try:
        return {"saved": sector_paths.settle_sectors(conn, params["sector_ids"], report, workers=params["workers"],
                                                     expand=False),
                "visited": len(params["sector_ids"])}
    finally:
        conn.close()


class SettleFollower(followed.Follower):
    """
    The settle stage of a run on the new queue (PERF.115): the wrapper selects the `settle-sector-paths` command from this
    run's steps, finds the sectors to settle once and cuts them into units (`expand_save_sector_paths`), this process
    follows the units with a bar of the sectors they settle, and the unit that completes the stage stores the stage's row.

    Args:
        started_at: The sectors saved since this (the database's clock when the run started) are the run's, unless
            `extra["sector_ids"]` names them.
        extra (dict, optional): More params for the stage: `sector_ids` and `expand`.

    Attributes:
        saved (int): How many paths were saved, summed over the units.
    """

    key = "settle"
    label = "Sector paths"

    def __init__(self, args, started_at, extra=None):
        super().__init__(args)
        self.started_at, self.extra = started_at, extra or {}
        self.saved = 0
        self.bar = None

    def params(self):
        return {"since": since_param(self.started_at), **self.extra}

    def on_plan(self, plan):
        if plan["detail"].get("skip"):
            log.normal(f"Nothing to settle: {plan['detail']['skip']}.")
            return
        stages.note(sectors=plan["detail"]["sectors"])
        log.normal("Saving the sector paths...")
        self.bar = self.hold(steps.Step("Sector paths", "paths", 1, args=self.args, progress=self.progress,
                                        workers=self.workers))
        self.bar.update(total=max(plan["detail"]["sectors"], 1))

    def on_result(self, unit, result):
        self.saved += result["saved"]
        if self.bar is not None:
            self.bar.update(advance=result["visited"])

    def on_done(self, conclusion):
        if self.bar is not None:
            log.normal(f"Saved {self.saved:,} sector paths.")


def _settle_on_queue(args, started_at, extra=None):
    """The settle stage as a run of its own (`SettleFollower`). Returns how many paths were saved."""
    return followed.run_planned(args, planned_settle(args, started_at, extra))


def planned_settle(args, started_at, extra=None):
    """The settle stage as a `followed.Planned`; its result is how many paths were saved."""
    follower = SettleFollower(args, started_at, extra)
    return followed.Planned([follower], "Sector paths", lambda: follower.saved)


def generate_and_save_sector_at(args, address, position_pc, edge_pc, channel=None):
    """
    Generates one sector via `generate_sector` inside its real grid cell
    and saves it -- the per-sector unit of work every `galaxy` mode
    repeats. The bright-star backfill runs once the whole run is done
    (the `backfill` stage, GEN.30), not per sector.

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
                           disk_scale_length_pc=fill.shape.disk_scale_length_pc if fill is not None else None,
                           shape=fill.shape if fill is not None else None
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
        # (the `link` stage), instead of one at a time under a lock.
        units = len(sector.entries) + len(sector.phenomena)
        with steps.worker_step(channel, "save:" + ",".join(str(part) for part in address),
                               f"Saving sector {provisional_sector_designation(*address)}", "save", units,
                               workers=run_common._worker_count(args)) as save_step:
            sector_id = store.save_sector(sector, config=store.mysql_config_from_args(args),
                                        galaxy_position=galaxy_position,
                                        link_neighbors=not getattr(args, "link_later", False),
                                        progress=save_step, run_id=getattr(args, "generation_run_id", None))
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
    command line to parse -- `regenerate_sector_cell` in particular,
    which runs for an API job.

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


def _since_text(moment):
    """The text of a start time as `--since` takes it (the one place that spells it, so a message never words a value of its own)."""
    return moment.strftime("%Y-%m-%d %H:%M:%S")


def _resolve_since(args):
    """Fills in `args.since` when `--since` was not given: the start of the newest `sectors` or `galaxy` run.

    Raises:
        SystemExit: No such run is recorded, so there is no start to work from.
    """
    if args.since is not None:
        return
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        args.since = store.newest_sectors_run_start(conn)
    finally:
        conn.close()
    if args.since is None:
        log.error("No 'sectors' or 'galaxy' run is recorded, so there is no start time: give --since TIME.")
        raise SystemExit(1)
    log.normal(f"Working on the sectors saved since {_since_text(args.since)}, the start of the newest "
               f"'sectors' or 'galaxy' run.")


def _sectors_made_since(args):
    """The ids of the galaxy sectors saved at or after `args.since` (the database's clock, `--since`)."""
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        return sector_paths.sector_ids_since(conn, args.since)
    finally:
        conn.close()


def run_link(args):
    """
    `planetgen link --since TIME` (PERF.123): links the sectors saved since `TIME` to their neighbours (the
    `link-neighbours` command, one unit): containment, nearest systems and quadrants, and the new systems merged into the
    lists of the sectors around them. Nothing else; a failure ends the run.
    """
    _resolve_since(args)
    stages.begin(stages.for_command("link", args), args, "link")
    created = _sectors_made_since(args)
    if not created:
        stages.skip("link", "no sector was saved since then")
        log.normal(f"No sector was saved since {args.since}: nothing to link.")
        stages.finish()
        return
    stages.enter("link")
    _link_on_queue(args, args.since)
    stages.finish()


def run_backfill(args):
    """
    `planetgen backfill --since TIME` (PERF.123): the bright-star backfill out from the sectors saved since `TIME`
    (`scatter-backfill`, GEN.30): the sectors a face away from them get the stars born from 1, 2, 5 and 8 solar masses
    up, each only where it is not drawn yet. Nothing else.
    """
    _resolve_since(args)
    stages.begin(stages.for_command("backfill", args), args, "backfill")
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        any_made = bool(store.sector_centers_since(conn, args.since))
    finally:
        conn.close()
    if not any_made:
        stages.skip("backfill", "no sector was saved since then")
        log.normal(f"No sector was saved since {args.since}: nothing to backfill.")
        stages.finish()
        return
    stages.enter("backfill")
    _backfill_on_queue(args, args.since)
    stages.finish()


def run_settle(args):
    """
    `planetgen settle --since TIME` (PERF.123, GEN.126): saves the path of every body in the sectors saved since `TIME` and
    in the sectors around them (`settle-sector-paths`), so each path sees the final set of neighbours. Nothing else.
    """
    if args.sector_id:
        _settle_named(args)
        return
    _resolve_since(args)
    stages.begin(stages.for_command("settle", args), args, "settle")
    created = _sectors_made_since(args)
    if not created:
        stages.skip("settle", "no sector was saved since then")
        log.normal(f"No sector was saved since {args.since}: nothing to settle.")
        stages.finish()
        return
    stages.enter("settle")
    _settle_on_queue(args, args.since)
    stages.finish()


def _settle_named(args):
    """`planetgen settle --sector-id ID...`: the settle stage over the sectors named (and the sectors around them), for a
    sector a `phenomenon --sector-id` run or an edit changed, which no `--since` time finds (PERF.124)."""
    ids = sorted(set(args.sector_id))
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        known = {row["id"] for row in conn.execute(
            f"SELECT id FROM sectors WHERE id IN ({', '.join('?' * len(ids))})", ids).fetchall()}
        started_at = store.database_now(conn)
    finally:
        conn.close()
    unknown = [sector_id for sector_id in ids if sector_id not in known]
    if unknown:
        log.error(f"Error: --sector-id: no such sector: {', '.join(str(sector_id) for sector_id in unknown)}")
        raise SystemExit(1)
    stages.begin(stages.for_command("settle", args), args, "settle")
    stages.enter("settle")
    _settle_on_queue(args, started_at, extra={"sector_ids": ids})
    stages.finish()


def settle_named_sectors(sector_ids, config, expand=True):
    """
    Saves the sector paths of `sector_ids` (and, with `expand`, of the sectors around them) as a run of
    `settle-sector-paths` (GEN.126; PERF.119: the settle of a regenerate or a filled sector, and of `planetgen settle`).

    Returns:
        int: How many paths were saved.
    """
    sector_ids = sorted(sector_ids)
    if not sector_ids:
        return 0
    args = _default_generation_args(config=config)
    args.workers = 1
    conn = store.get_connection(config)
    try:
        started_at = store.database_now(conn)
    finally:
        conn.close()
    return _settle_on_queue(args, started_at, extra={"sector_ids": sector_ids, "expand": expand})


def _cell_to_generate(address, config, link=True):
    """
    What generating the sector of grid cell `address` takes (`regenerate_sector_cell`): the
    answer already when there is nothing to generate, else the arguments, the center and the sector edge to generate with.

    Returns:
        tuple: `(early, args, position_pc, qualifies, edge_pc)`. `early` is the answer
            for a cell that already holds a sector or lies past the grid's reach, else `None` (and the rest is set).

    Raises:
        RuntimeError: The galaxy's skeleton has never been built.
    """
    conn = store.get_connection(config)
    try:
        existing_id = store.get_sector_id_at(conn, *address)
        if existing_id is not None:
            return ({"created": False, "qualifies": True, "in_reach": True, "sector_id": existing_id,
                     "sector_name": None}, None, None, None, None)

        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError(
                "The galaxy's skeleton has never been built (no galaxy_shape row) -- run "
                "'planetgen plan' first."
            )

        bounds = store.get_galaxy_bounds(conn)
    finally:
        conn.close()

    if not bounds.contains(address[0], address[1]):
        # Past the grid's reach: the cell has no address (the object ID cannot name it).
        return ({"created": False, "qualifies": False, "in_reach": False, "sector_id": None, "sector_name": None},
                None, None, None, None)

    position_pc = sector_position_pc(*address, skeleton.edge_pc)
    qualifies = bounds.qualifies(address[0], address[1])

    # Every address in the reach generates, however sparse (GEN.76): the
    # profile's density sets the draw's expected count, never whether the
    # draw happens.
    args = _default_generation_args(config=config)
    args.density = sector_profile(address, skeleton).density
    args.num_systems = None
    args.link_later = not link
    return None, args, position_pc, qualifies, skeleton.edge_pc


def _refresh_block_totals(config, address):
    """
    Brings the Galaxy Map's block totals up to the sector just saved at `address` (DB.83): the generated totals
    (`db.block_stats.refresh_chain`, which a sector save does not do) and the light of the pre-placed objects
    (`db.block_lights.refresh_chain`). A failure only warns: `planetgen fill-block-totals` makes them again.
    """
    try:
        conn = store.get_connection(config)
        try:
            block_stats.refresh_chain(conn, address)
            block_lights.refresh_chain(conn, address)
            conn.commit()
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- the sector is saved; the totals can be made again
        log.normal(f"Warning: adding the sector at {_format_address(address)} to the block totals failed: {exc}")


def regenerate_sector_cell(address, config):
    """
    Generates the sector of grid cell `address` again as the runs of an API job (PERF.119; the new sector of
    `api_jobs.regenerate_sector`, whose old one is deleted): the sector as a run of `generate-sectors` and the bright-star
    backfill around it as a run of `scatter-backfill`. The caller settles the paths.

    Returns:
        dict: `sector_id` and `sector_name` of the new sector (`None` when the cell is beyond the grid's reach).
    """
    early, args, position_pc, _qualifies, edge_pc = _cell_to_generate(address, config)
    if early is not None:
        return {"sector_id": early["sector_id"], "sector_name": early["sector_name"]}
    args.workers = 1
    conn = store.get_connection(config)
    try:
        started_at = store.database_now(conn)
    finally:
        conn.close()
    [saved] = _submit_batch(args, [(address, position_pc, args, "")], f"Sector ({_format_address(address)})", edge_pc,
                            count=1, collect=True)
    _backfill_on_queue(args, started_at, center_pc=position_pc)
    _refresh_block_totals(config, address)
    return {"sector_id": saved["sector_id"], "sector_name": saved["name"]}


SECTOR_UNITS_MAX = 4000
"""int: The most units a stage of sectors is cut into (PERF.107): about that many static chunks, however many workers."""

SECTOR_STAGE_KEYS = ("start", "neighborhood", "sectors")
"""tuple: The stages of a `galaxy` run that generate sectors; each is a command of its own (`commands.KIND_OF_STAGE`)."""


def _sector_args_json(args):
    """
    The arguments a sector's generation reads, as plain JSON: every value of the run's namespace that is a number,
    a string, a boolean, `None` or a list of those, without the password (a worker reads that from its own settings)
    or the connection (`database` carries it) or a value a stage worked out for itself (names that start with `_`).
    """
    def plain(value):
        if isinstance(value, (list, tuple)):
            return all(plain(item) for item in value)
        return value is None or isinstance(value, (str, int, float, bool))

    skipped = {"mysql_password", "mysql_host", "mysql_port", "mysql_user", "mysql_database"}
    return {name: list(value) if isinstance(value, tuple) else value for name, value in vars(args).items()
            if name not in skipped and not name.startswith("_") and plain(value)}


def _sector_args_of(unit_args):
    """The namespace a unit generates with: the run's arguments and the database the unit names, the password
    from this process's settings (`store.MySQLConfig`)."""
    named = unit_args["database"]
    return argparse.Namespace(**unit_args["args"], mysql_host=named["host"], mysql_port=named["port"],
                              mysql_user=named["user"], mysql_password=None, mysql_database=named["database"])


def _planned_around_centers(args, centers, edge_pc, config, batch_density):
    """
    The sectors of a "neighborhood" stage that follows the stage that makes the sectors at `centers` (addresses): those
    within `args.radius_pc` of each, in the order the centers are listed, as `(items, suffixes)` (`planned_items`). The
    stage before has saved the centers by the time this is asked, so each one's id is read here for the words after its
    sectors; its place is the one the plan was made for (`sector_position_pc`).
    """
    conn = store.get_connection(config)
    try:
        ids = [store.get_sector_id_at(conn, *address) for address in centers]
    finally:
        conn.close()
    located = []
    for address, sector_id in zip(centers, ids):
        if sector_id is None:
            raise RuntimeError(f"the sector at {_format_address(address)} was not saved, so there is no neighborhood to "
                               f"generate around it")
        candidates, occupied, _outside = _neighborhood_candidates(
            sector_position_pc(*address, edge_pc), args.radius_pc, edge_pc, config, batch_density.bounds)
        located.extend(_Pending.of_candidates(args, candidates, occupied, batch_density,
                                              f", {{distance:.2f}} pc from sector_id={sector_id}").located())
    return [index for address, _suffix in located for index in address], [suffix for _address, suffix in located]


def _planned_single(address, config, batch_density):
    """The sector of a "single" stage (PERF.214), `[shell, band, wedge]`, or `[]` when it exists: the `address` `--shell --band --wedge` names
    (the stage's own arguments are the defaults of `sector`, not the run's), anywhere in the grid's reach (GEN.235). The part of
    `_plan_single_sector` that reads the database, at the time the stage opens."""
    _require_in_reach(batch_density.bounds, address[0], address[1], _format_address(address))
    conn = store.get_connection(config)
    try:
        exists = store.get_sector_id_at(conn, *address) is not None
    finally:
        conn.close()
    return [] if exists else list(address)


def _planned_starts(args, edge_pc, config, batch_density, run_seed):
    """The starts of a "random-start" stage (PERF.214), as flat `shell, band, wedge` triples: the places `_plan_random_start` picks, from the
    run seed the starter drew, so the same seed gives the same starts wherever it is asked. Fewer than `--neighborhoods` is said; none
    fails the stage."""
    radius_pc = _random_start_radius(args)
    conn = store.get_connection(config)
    starts, drawn = [], {}
    try:
        for index in range(args.neighborhoods):
            start = _pick_random_start(args, edge_pc, batch_density, conn, radius_pc, [position for _address, position, _args in starts],
                                       drawn, run_seed=run_seed, index=index)
            if start is None:
                break
            starts.append(start)
    finally:
        conn.close()
    if not starts:
        raise RuntimeError(f"no unoccupied sector address was found after {program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS} attempts "
                           f"(this galaxy may already be almost entirely generated within the range asked, or hold too little stellar density)")
    if len(starts) < args.neighborhoods:
        log.normal(f"Found room for only {len(starts)} of {args.neighborhoods} neighborhoods (each whole inside the galaxy, at least "
                   f"{2 * radius_pc:g} pc apart).")
    return [index for address, _position, _args in starts for index in address]


def planned_items(params):
    """
    The sectors a `planned` sector stage generates (PERF.123), as `(items, suffixes)`: `items` are flat `shell, band,
    wedge` triples in the order to generate them, `suffixes` what the log says after each address (none for a stage
    whose sectors need none). The stage's mode says where they are: "addresses", the shell, band, span or course the
    arguments name (`_address_plan`); "neighborhood", the sectors within `--radius-pc` of `--center-sector` in the fill
    order (`_neighborhood_candidates`), or around the `centers` the stage before it made (`_planned_around_centers`). Either way, less
    the sectors that exist now, up to `--limit`. The starter has
    counted, estimated and asked about them already; a run planned with `_plan_sectors(deferred=True)` (PERF.214) has not, and two more
    modes are for it: "single", the one sector `--shell --band --wedge` names (none when it exists), and "random-start", the starts
    drawn from `params["run_seed"]` (the neighborhood after them reads their addresses from this stage's plan, `centers`).
    """
    args = _sector_args_of(params)
    edge_pc = params["edge_pc"]
    config = store.mysql_config_from_args(args)
    batch_density = _BatchDensity(config)
    if params["planned"] == "single":
        return _planned_single(tuple(params["address"]), config, batch_density), []
    if params["planned"] == "random-start":
        return _planned_starts(args, edge_pc, config, batch_density, params["run_seed"]), []
    if params["planned"] == "neighborhood" and params.get("centers") is not None:
        return _planned_around_centers(args, params["centers"], edge_pc, config, batch_density)
    if params["planned"] == "neighborhood":
        candidates, occupied, _outside = _neighborhood_candidates(
            _neighborhood_center(args), args.radius_pc, edge_pc, config, batch_density.bounds)
        pending = _Pending.of_candidates(args, candidates, occupied, batch_density, _neighborhood_suffix(args))
        located = list(pending.located())
        return [index for address, _suffix in located for index in address], [suffix for _address, suffix in located]
    plan = _address_plan(args, edge_pc, batch_density, say=_quiet)
    if plan is None:
        return [], []
    pending = _pending_of_plan(args, plan, edge_pc, batch_density)
    return [index for address in pending.addresses() for index in address], []


def expand_generate_sectors(params):
    """
    The sector stages' cut (PERF.107; `generate-start-sector`, `generate-neighborhood` and `generate-sectors`): the
    addresses of `params["items"]` (the run's sectors to generate, flat `shell, band, wedge` triples, in the order to
    generate them; for a stage that is `planned`, the ones `planned_items` works out when the stage opens) in about `SECTOR_UNITS_MAX` static chunks of consecutive sectors, the same at any worker count.
    A unit's cost is its sectors.
    """
    if getattr(_sector_args_of(params), "around", None) is not None:
        return _expand_region(params)
    warning = None
    if params.get("check_plan"):        # the stage of a run that plans the galaxy first (PERF.203) finds the plan when it opens
        problem, warning = _plan_check(store.mysql_config_from_args(_sector_args_of(params)), params["edge_pc"])
        if problem:
            raise RuntimeError(problem)
    items, suffixes = planned_items(params) if params.get("planned") else (params["items"], [])
    count = len(items) // 3
    size = max(1, -(-count // SECTOR_UNITS_MAX))
    units = []
    for start in range(0, count, size):
        chunk = items[3 * start:3 * (start + size)]
        unit = {"database": params["database"], "args": params["args"], "edge_pc": params["edge_pc"], "items": chunk}
        if params.get("again"):
            unit["again"] = True
        if suffixes:
            unit["suffixes"] = suffixes[start:start + size]
        units.append((unit, float(len(chunk) // 3)))
    detail = {"sectors": count}
    if warning:
        detail["warning"] = warning
    if params.get("planned") == "random-start":
        detail["centers"] = [list(items[place:place + 3]) for place in range(0, len(items), 3)]       # the neighborhood stage reads them
    return commands.Expansion(units, detail)


def _saved_addresses(args):
    """The `(shell, band, wedge)` of every galaxy sector saved now: a stage run again after it stopped (PERF.182) leaves them
    alone, for each was saved whole in a transaction of its own."""
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        return store.filled_sector_addresses(conn)
    finally:
        conn.close()


def generate_sectors_unit(unit_args, report):
    """
    The unit of the sector stages (PERF.107): generates and saves each sector of the chunk (`generate_and_save_sector_at`,
    one transaction each) in its own grid cell, its density resolved from its sector profile as `_BatchDensity` does for
    the run (a density or count the run was given is used for every sector). Every draw for a sector comes from its
    own seed (GEN.39), so the stars are the same at any worker count. A unit of a stage run again after it stopped (`again`,
    PERF.182) skips the sectors already saved.

    Returns:
        dict: `pid` (this process: a unit run by the process that started the run has counted its sectors there already),
            `sectors`, one entry per sector in order: `address`, `sector_id`, `name`, `systems`, `stars`, `density`,
            `phenomena`, `summary` (`sector_generation_summary_lines`), `seconds` and `suffix` (what the log says after the address).
    """
    args = _sector_args_of(unit_args)
    edge_pc = unit_args["edge_pc"]
    batch_density = _BatchDensity(store.mysql_config_from_args(args))
    items = unit_args["items"]
    suffixes = unit_args.get("suffixes", [])
    saved = _saved_addresses(args) if unit_args.get("again") else ()
    done = []
    for place in range(0, len(items), 3):
        address = tuple(items[place:place + 3])
        if address in saved:
            continue
        position_pc = sector_position_pc(*address, edge_pc)
        sector_args = batch_density.resolve(args, address, position_pc)
        started = time.monotonic()
        sector_id, sector_name, sector = generate_and_save_sector_at(sector_args, address, position_pc, edge_pc)
        done.append({
            "address": list(address), "sector_id": sector_id, "name": sector_name, "systems": len(sector.entries),
            "stars": run_sector.sector_star_count(sector), "density": run_common._sector_density(sector_args),
            "phenomena": len(sector.phenomena), "summary": run_sector.sector_generation_summary_lines(sector, sector_args),
            "seconds": time.monotonic() - started, "suffix": suffixes[len(done)] if suffixes else ""})
        report(len(done), len(items) // 3)
    return {"sectors": done, "pid": os.getpid()}


def _log_saved(saved, address, suffix=""):
    designation = provisional_sector_designation(*address)
    log.normal(
        f"Saved sector '{saved['name']}' [{designation}] at {_format_address(address)}{suffix} "
        f"(sector_id={saved['sector_id']})."
    )
    run_sector._log_summary(saved["summary"])


class SectorsFollower(followed.Follower):
    """
    The sectors of a stage of a run on the new queue (PERF.107): the wrapper selects the stage's command
    (`generate-start-sector`, `generate-neighborhood` or `generate-sectors`, by the stage's key) from the run's steps, cuts
    the sectors into units (`expand_generate_sectors`) and queues them, and this process follows the units' results: it
    logs each sector as it is saved, counts it and adds its time to the speed stats, and advances the stage's bar. The unit
    that completes the stage stores the stage's row.

    Args:
        key (str): The stage: `start`, `neighborhood` or `sectors`.
        label (str): What the stage is called.
        edge_pc (float): The sector edge, parsecs.
        batch (iterable): The sectors as `(address, position_pc, sector_args, suffix)` items, when the stage is given them
            (a stage that is `planned` works them out when it opens, and `batch` is empty).
        count (int, optional): How many sectors the stage generates (default `len(batch)`).
        bar_description (str, optional): The words of the stage's bar of sectors; none draws no bar.
        base_args (argparse.Namespace, optional): What the units generate with when it is not the run's own `args`
            (a single sector of `galaxy --wedge` is drawn from its profile with the defaults of `sector`).
        collect (bool): Whether to keep what each sector reported (`reported`) instead of logging it, for a stage of few
            sectors whose caller says it in its own words (`announce`).
        planned (str): The way the stage's cut works the sectors out from the run's arguments (`planned_items`:
            "addresses" or "neighborhood") instead of taking them from `batch`; empty when `batch` holds them.
        extra (dict, optional): More for the stage's parameters (what the "neighborhood" cut needs: `centers`).
        announce (callable, optional): `announce(saved, address, suffix)` says one sector as it is saved (default
            `_log_saved`); not called for a stage that `collect`s.
        finish (callable, optional): `finish(follower)` when the stage is done, for the lines that sum it up.

    Attributes:
        reported (list): What each sector reported, in the order it was saved (`generate_sectors_unit`), when `collect`.
    """

    def __init__(self, args, key, label, edge_pc, batch=(), count=None, bar_description=None, base_args=None,
                 collect=False, planned="", extra=None, announce=None, finish=None):
        super().__init__(args)
        self.key, self.label, self.edge_pc = key, label, edge_pc
        self.batch = batch
        self.count = len(batch) if count is None else count
        self.bar_description, self.base_args, self.collect = bar_description, base_args, collect
        self.planned, self.extra = planned, extra or {}
        self.announce, self.finish = announce or _log_saved, finish
        self.reported = []
        self.generated = 0
        self.bar = None
        self.items = []
        self.suffixes = {}
        for address, _position_pc, _sector_args, suffix in batch:
            self.items.extend(address)
            if suffix:
                self.suffixes[tuple(address)] = suffix

    def params(self):
        return {"args": _sector_args_json(self.base_args or self.args), "edge_pc": self.edge_pc, "items": self.items,
                "planned": self.planned, **self.extra}

    def on_plan(self, plan):
        sectors = plan["detail"]["sectors"]
        if plan["detail"].get("warning"):
            log.normal(f"Warning: {plan['detail']['warning']}")
        if not self.collect:
            log.normal(f"Generating {sectors:,} sector(s) with {self.workers} worker(s); each is reported below as "
                       f"it is saved, so the first report can take a while in a dense region.")
        if self.bar_description is not None:
            self.bar = self.hold(run_common.sector_bar(self.progress, self.args, self.bar_description, sectors))

    def on_result(self, unit, result):
        elsewhere = result["pid"] != os.getpid()     # a worker's own RUN_COUNTS die with it; the run's are here
        sector_args = self.base_args or self.args
        for saved in result["sectors"]:
            address = tuple(saved["address"])
            if elsewhere:
                run_common.RUN_COUNTS["sectors"] += 1
                run_common.RUN_COUNTS["systems"] += saved["systems"]
                run_common.RUN_COUNTS["phenomena"] += saved["phenomena"]
            run_common._record_sector(sector_args, saved, saved["seconds"])
            self.generated += 1
            if self.bar is not None:
                self.bar.update(advance=1)
            if self.collect:
                self.reported.append(saved)
            else:
                self.announce(saved, address, saved.get("suffix") or self.suffixes.get(address, ""))

    def on_done(self, conclusion):
        if self.finish is not None:
            self.finish(self)


def _submit_batch(args, batch, title, edge_pc, count=None, base_args=None, collect=False, planned="", extra=None,
                  bar_description=None, key="sectors"):
    """
    Generates every `(address, position_pc, sector_args, suffix)` of `batch` as a run of its own (`SectorsFollower`);
    `batch` may be a generator (PERF.36), and `count` says how many it holds.

    Returns:
        list: What each sector reported, in the order it was saved (`generate_sectors_unit`), when `collect`; else `[]`.
    """
    follower = SectorsFollower(args, key, title, edge_pc, batch, count, bar_description, base_args, collect, planned,
                               extra)
    followed.run_followed(args, [follower], title)
    return follower.reported


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


def _neighborhood_center(args):
    """The stored center, in parsecs, of `args.center_sector`.

    Raises:
        SystemExit: The sector does not exist, or has never been placed in a galaxy.
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
    return (center_position["center_x_pc"], center_position["center_y_pc"], center_position["center_z_pc"])


def _neighborhood_suffix(args):
    """What the log says after a neighborhood sector's address; `{distance}` is its distance from the center, in pc."""
    return f", {{distance:.2f}} pc from sector_id={args.center_sector}"


def _plan_neighborhood(args, edge_pc, progress):
    """
    Local-neighborhood mode: generates every not-yet-generated
    sector within `args.radius_pc` parsecs of `args.center_sector`'s own
    stored center -- see `_BatchDensity.resolve` for how each sector's density is set. The plan is made here (counted,
    estimated, confirmed); the stage that generates the sectors cuts them again when it opens (`planned_items`).

    Args:
        args (argparse.Namespace): Parsed arguments;
            `args.center_sector`/`args.radius_pc` must not be `None`.
        edge_pc (float): The sector edge length, in parsecs (`_edge_pc`).
        progress (rich.progress.Progress): The display the plan's questions are asked on (see `_plan_sectors`).

    Returns:
        list[SectorsFollower]: The stage.

    Raises:
        SystemExit: If `args.center_sector` doesn't exist, or has never
                   been placed in a galaxy.
    """
    center = _neighborhood_center(args)
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

    pending = _Pending.of_candidates(args, candidates, occupied, batch_density, _neighborhood_suffix(args))
    batch, estimate_args, total = _plan_pending(pending)
    already_existed, skipped = pending.existed, pending.skipped
    run_common._check_estimate(args, estimate_args,
                    f"the sectors {region} of sector {args.center_sector}", progress, total=total)

    def finished(_follower):
        skip_note = f", {skipped} skipped (beyond the grid's reach)" if skipped else ""
        outside_note = f", {outside} beyond the galaxy's edge left out" if outside else ""
        log.normal(
            f"Generated {pending.count} new sector(s) {region} of sector_id={args.center_sector} "
            f"({len(candidates)} candidate sector(s) inside the galaxy{outside_note}, {already_existed} already "
            f"existed{skip_note})."
        )

    return [SectorsFollower(args, "sectors", f"Sectors ({region} of sector {args.center_sector})", edge_pc,
                            count=pending.count, bar_description="Sectors (local neighborhood)",
                            planned="neighborhood", finish=finished)]


def _neighborhood_after(args, edge_pc, config, batch_density, places, radius_pc):
    """
    The plan of the neighborhoods that follow stages that make sectors at `places` (`(address, position_pc)` each): the
    sectors within `radius_pc` of each, less the center itself, which the stage before makes. Returns `(estimate_args,
    count, existed, skipped)`: the arguments of every sector of them, for the estimate, how many there are, and how many
    were left out as already generated or beyond the grid's reach.
    """
    estimate_args, count, existed, skipped = [], 0, 0, 0
    for address, position_pc in places:
        candidates, occupied, _outside = _neighborhood_candidates(position_pc, radius_pc, edge_pc, config,
                                                                  batch_density.bounds)
        pending = _Pending.of_candidates(args, candidates, occupied, batch_density, skip={address})
        batch = list(pending.items())
        estimate_args += [item[2] for item in batch]
        count += len(batch)
        existed, skipped = existed + pending.existed, skipped + pending.skipped
    return estimate_args, count, existed, skipped


def _neighborhood_follower(args, edge_pc, centers, count, existed, skipped):
    """The "neighborhood" stage that follows the stage that makes the sectors at `centers`: their `args.radius_pc`
    neighborhoods, cut when the stage opens (`_planned_around_centers`)."""
    def finished(_follower):
        skip_note = f", {skipped} skipped (beyond the grid's reach)" if skipped else ""
        log.normal(f"Generated {count} new sector(s) within {args.radius_pc:g} pc of "
                   f"{'the starting sector' if len(centers) == 1 else f'{len(centers)} starting sectors'} "
                   f"({existed} already existed{skip_note}).")

    return SectorsFollower(args, "neighborhood", "Neighborhood", edge_pc, count=count,
                           bar_description="Sectors (local neighborhood)", planned="neighborhood",
                           extra={"centers": [list(address) for address in centers]}, finish=finished)


class GenerationRefused(RuntimeError):
    """PERF.3: a bulk generation the database disk can't hold; the
    message says why and how much it needs."""


class DatasetsMissing(RuntimeError):
    """GEN.590: a file the generator reads was never downloaded (`planetgen.datasets`), so a bulk generation was refused
    before writing anything; the message names the datasets."""


def require_datasets(disregard=None):
    """
    The second gate in front of every bulk generation (GEN.590, Boss's "Require the data"): the data files star composition
    reads (isochrones, white-dwarf sequences, neutron-star tables) must be in the cache folder that setup fills. Nothing is
    downloaded here (setup does that); the files are looked at and hashed once per process.

    A run told to disregard missing data (GEN.599) is not refused: the datasets that are missing are returned, so the run
    can say that it takes each one's labelled fallback.

    Args:
        disregard (bool, optional): Whether the run was told to disregard missing data. Default: what this process was told
            (`fetch.disregarding_missing_data`).

    Returns:
        list: The `fetch.Located` of each dataset that is missing and will be replaced by its fallback (empty when every
            file is there).

    Raises:
        DatasetsMissing: Naming each dataset that is missing and why, each file's URL and the setup command.
    """
    disregard = fetch.disregarding_missing_data() if disregard is None else disregard
    missing = datasetCatalogue.missing()
    if not missing or (disregard and all(located.dataset.fallback for located in missing)):
        return missing
    raise DatasetsMissing(
        "the downloaded data is missing, so nothing was generated: " + "; ".join(located.describe() for located in missing)
        + ". Setup downloads it ('python3 -m planetgen.cli.datasets fetch'); to generate anyway with a labelled fallback "
        "for each missing dataset, say so with --disregard-missing-data.")


def disregard_notice(missing):
    """The lines a run that disregards missing data prints at its start: that it does, and each missing dataset with
    the fallback it takes (GEN.599). Empty when nothing is missing."""
    if not missing:
        return []
    return (["Missing data is disregarded (--disregard-missing-data): a labelled fallback stands in for each missing dataset."]
            + [f"  {located.dataset.key} ({located.dataset.fill}) is missing ({located.reason}); fallback: {located.dataset.fallback}"
               for located in missing])


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


def _generate_batch_as_runs(args, batch, title, edge_pc, config, center_pc):
    """
    The runs of an API job that generates sectors (PERF.119): the sectors of `batch` as one run of `generate-sectors`,
    then, when any were made, the bright-star backfill around `center_pc` (GEN.30, once, around the requested sector) as a
    run of `scatter-backfill` and the paths of the new sectors and their neighbours as a run of `settle-sector-paths`
    (GEN.126), each through the wrapper like the stages of `planetgen galaxy`. A run that finds the machine taken waits in
    line for it (PERF.153).

    Returns:
        tuple: `(created_ids, generated)`: the ids of the sectors made, in the order they were saved, and how many.
    """
    conn = store.get_connection(config)
    try:
        started_at = store.database_now(conn)
    finally:
        conn.close()
    saved = _submit_batch(args, batch, title, edge_pc, count=len(batch), collect=True)
    created_ids = [item["sector_id"] for item in saved]
    if created_ids:
        _backfill_on_queue(args, started_at, center_pc=center_pc)
        _settle_on_queue(args, started_at, extra={"sector_ids": sorted(created_ids)})
    return created_ids, len(created_ids)


def generate_sector_neighborhood(center_sector_id, radius_ly=None, config=None, estimate_only=False):
    """
    Non-CLI counterpart to `_plan_neighborhood` -- for the admin web
    UI's "generate more sectors around this one" action
    (`planetgen/api/routes.py`'s `generate_sector_neighborhood_route`). Same
    work, a plain result dict instead of prints, and a catchable
    `ValueError` instead of `SystemExit` for an invalid/unplaced sector.

    The default radius is `program_constants.DEFAULT_GENERATE_RADIUS_PC`
    (12 pc, about 100 candidate addresses; GEN.23). Once they are all
    generated, the 100 ly around the center sector gets its bright stars
    (the `scatter-backfill` stage, GEN.30).
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
    args.workers = 1   # one worker generates the run's sectors, one at a time

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

    _created_ids, generated = _generate_batch_as_runs(args, batch, f"Sectors (around sector {center_sector_id})", edge_pc,
                                                     config, center_pc=center)
    return {"generated": generated, "estimate": estimate.as_dict(), **counts}


def generate_block(block_key, config=None, estimate_only=False):
    """
    Non-CLI counterpart to `planetgen sectors --block` (GEN.310) -- for the Galaxy Map's "generate this block" action
    (`planetgen/api/routes.py`'s `generate_block_route`). The same work as the command: every not-yet-generated
    sector the stored outline holds inside the block `block_key` names (`drill.block_addresses`), shell by shell
    outward, as runs through the wrapper (`_generate_batch_as_runs`), with the bright-star backfill and the settle
    after them for the sectors made. A plain result dict instead of prints, and a catchable `ValueError` for a key that
    names no block.

    Args:
        block_key (str): `m.shell.band.wedge` of a block of a stage level (`drill.parse_stage_key`).
        config (MySQLConfig, optional): Connection parameters.
        estimate_only (bool): Only work out the size and time (PERF.3); nothing is written.

    Returns:
        dict: `generated`, `already_existed`, `skipped`, `candidates` (the block's addresses inside the outline) and
            `outside_galaxy` (its addresses past the outline, left out) -- all int -- and `estimate`
            (`generationStats.Estimate.as_dict`).

    Raises:
        ValueError: `block_key` names no block of a stage level.
        EnumerationTooLarge: The block holds more sectors than a web request may list (PERF.36).
        GenerationRefused: The database disk can't hold it (nothing was written).
        RuntimeError: The galaxy's skeleton has never been built.
        DatasetsMissing: A downloaded data file is missing (a real run only); nothing is written.
    """
    block = drill.parse_stage_key(block_key)
    edge_pc = run_common._edge_pc()
    config = config or store.DEFAULT_MYSQL_CONFIG
    if not estimate_only:
        require_datasets()
    batch_density = _BatchDensity(config)
    bounds = batch_density.bounds

    inside, outside, shells = 0, 0, set()
    for shell, band, _wedge in drill.block_addresses(block):
        if bounds.qualifies(shell, band):
            inside += 1
            shells.add(shell)
            if inside > run_common.WEB_ENUMERATION_LIMIT:
                raise EnumerationTooLarge(
                    f"Block {drill.format_drill_key(block)} holds more than {run_common.WEB_ENUMERATION_LIMIT:,} sectors, "
                    f"too many to list in one request: fill a smaller block, or run 'planetgen sectors --block' with --limit.")
        else:
            outside += 1

    def addresses():
        return ((shell, band, wedge) for shell, band, wedge in drill.block_addresses(block)
                if bounds.qualifies(shell, band))

    conn = store.get_connection(config)
    try:
        occupied = set(store.get_occupied_addresses(conn, shells))
    finally:
        conn.close()

    args = _default_generation_args(config=config)
    args.density = None            # density-driven from the skeleton (_BatchDensity), not the flat default
    args.num_systems = None
    args.workers = 1               # one worker generates the run's sectors, one at a time

    pending = _Pending.of_addresses(args, addresses, occupied, batch_density, edge_pc)
    batch, estimate_args, total = _plan_pending(pending)
    counts = {"already_existed": pending.existed, "skipped": pending.skipped, "candidates": inside,
              "outside_galaxy": outside}
    estimate = run_common._estimate_sectors(args, estimate_args, total)
    if estimate_only:
        run_common._RUN_STATS.clear()
        return {"generated": 0, "estimate": estimate.as_dict(), **counts}
    if estimate.refusal:
        run_common._RUN_STATS.clear()
        raise GenerationRefused(estimate.refusal)

    _created_ids, generated = _generate_batch_as_runs(
        args, batch, f"Sectors (block {drill.format_drill_key(block)})", edge_pc, config, center_pc=None)
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


def _pick_random_start(args, edge_pc, batch_density, conn, radius_pc, taken, drawn=None, run_seed=None, index=0):
    """
    One random start for random-start mode: the anchor sector of an uncharted region (GEN.826: a 9-block with no generated
    sector), the region found by drawing an address
    inside the outline by volume (`GalaxyBounds.random_address`); where the anchor is outside the outline's reach, or beyond
    `--max-shell`, the drawn address itself starts the fill. An address is retried up
    to `RANDOM_START_MAX_PLACEMENT_ATTEMPTS` times until its density expects at least one system (the outline is a
    superset of the sectors that do, GEN.627, and a start in an empty sector is no start; skipped when the density is given
    explicitly), until it meets `--min-start-density`, and, for `--neighborhoods` (GEN.97), until its
    whole neighborhood lies inside the outline, it is at least twice the
    radius from every start in `taken`, and it passes the
    `--neighborhood-gamma` draw (accepted with probability
    `min(1, density) ** gamma`; gamma 0, the default, accepts every one).
    With `--avoid-filled-space` (GEN.186) a centre whose neighborhood
    touches an already generated sector is retried too; `drawn`, when
    given, is a dict that gets `centres` (examined) and `filled` (dropped
    for touching filled space) counted up. With `run_seed` the start is drawn from the stream of that seed and
    `index` (the same seed and index give the same start at any worker count, PERF.123); without it, from the
    run's own stream.

    Returns:
        tuple or None: `(address, position_pc, sector_args)`; `None` when
            the attempt budget ran out.
    """
    if run_seed is not None:
        from planetgen.queue import work
        with draw.bound(work.task_seed(run_seed, f"random-start-{index}")):
            return _pick_random_start(args, edge_pc, batch_density, conn, radius_pc, taken, drawn)
    many = args.neighborhoods > 1
    gamma = args.neighborhood_gamma
    for _ in range(program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS):
        address = batch_density.bounds.random_address(draw, max_shell=args.max_shell)
        region = named_regions.region_of(address)
        if region_stats.region_row(conn, region) is not None:
            continue                # GEN.826: a start is in a region nothing is generated in yet ...
        anchor = named_regions.anchor_address(region, edge_pc)
        if batch_density.bounds.contains(anchor[0], anchor[1]) and (args.max_shell is None or anchor[0] <= args.max_shell):
            address = anchor        # ... at its anchor, so the fill of the start's sphere charts the whole region
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
        if (sector_args is not args
                and sector_args.density * batch_density.skeleton.expected_system_count_at_density_1 < 1.0):
            continue
        if args.min_start_density is not None and sector_args.density < args.min_start_density:
            continue
        if gamma and many and draw.random() >= min(1.0, sector_args.density) ** gamma:
            continue
        return address, position_pc, sector_args
    return None


def _random_start_radius(args):
    """The radius of a random start's neighborhood: `--radius-pc`, or the default."""
    return args.radius_pc if args.radius_pc is not None else program_constants.DEFAULT_GENERATE_RADIUS_PC


def _log_start(saved, address, _suffix):
    """Says a random starting sector as it is saved."""
    designation = provisional_sector_designation(*address)
    log.normal(
        f"Saved random starting sector '{saved['name']}' [{designation}] at {_format_address(address)} "
        f"(sector_id={saved['sector_id']})."
    )
    run_sector._log_summary(saved["summary"])


def _plan_random_start(args, edge_pc, progress):
    """
    Random-start mode (no `--shell`/`--center-sector` given): picks a
    random, not-yet-occupied sector address -- drawn only from
    inside the galaxy's stored outline (`GalaxyBounds.random_address`,
    every sector equally likely, optionally only out to `--max-shell`), so
    the start and the neighborhood around it are always in the galaxy --
    retried up to `program_constants.RANDOM_START_MAX_PLACEMENT_ATTEMPTS`
    times, generates it, then generates every not-yet-generated sector
    within `args.radius_pc` of it (default
    `program_constants.DEFAULT_GENERATE_RADIUS_PC`, 12 pc). `--min-start-density` tightens the retry:
    an address below it is retried too.

    With `--neighborhoods N` (GEN.97) it makes N of them: all N starts are
    chosen first (see `_pick_random_start`) and the estimate covers all of
    them; the starts are one stage and their neighborhoods the next (PERF.172). Too few places for N is a
    warning when at least one was found.

    Returns:
        list[SectorsFollower]: The stage of the starts and the stage of their neighborhoods.

    Raises:
        SystemExit: If no suitable address was found within the attempt
                   budget.
    """
    radius_pc = _random_start_radius(args)
    run_seed = draw.getrandbits(128)
    log.debug(f"Random starts are drawn from the run seed {run_seed}.")

    mysql_config = store.mysql_config_from_args(args)
    batch_density = _BatchDensity(mysql_config)
    conn = store.get_connection(mysql_config)
    starts = []
    drawn = {}
    try:
        for index in range(args.neighborhoods):
            start = _pick_random_start(args, edge_pc, batch_density, conn, radius_pc,
                                       [position for _address, position, _args in starts], drawn,
                                       run_seed=run_seed, index=index)
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

    # The estimate covers every start and its whole neighborhood, before any is written.
    around, count, existed, skipped = _neighborhood_after(
        args, edge_pc, mysql_config, batch_density, [(address, position) for address, position, _args in starts],
        radius_pc)
    estimate_args = [sector_args for _address, _position, sector_args in starts] + around
    first_address = starts[0][0]
    what = (f"a random start at {_format_address(first_address)} and the sectors within {radius_pc:g} pc of it"
            if len(starts) == 1 else f"{len(starts)} random neighborhoods of {radius_pc:g} pc")
    run_common._check_estimate(args, estimate_args, what, progress)

    args.radius_pc = radius_pc
    addresses = [address for address, _position, _args in starts]
    return [SectorsFollower(args, "start", "Starting sector", edge_pc,
                            [(address, position, sector_args, "") for address, position, sector_args in starts],
                            announce=_log_start),
            _neighborhood_follower(args, edge_pc, addresses, count, existed, skipped)]


def _plan_single_sector(args, edge_pc, progress):
    """
    Single-address mode: generates exactly the one sector at `(--shell I,
    --band J, --wedge K)` -- the direct route from
    an address copied out of the interactive 3D Galaxy Map. With
    `--radius-pc`, then generates that sector's neighborhood
    (the Sector Map's "Generate neighborhood"), as the stage after it.

    Returns:
        list[SectorsFollower]: The stage of the sector (none when it exists) and, with `--radius-pc`, of its neighborhood.

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
        existing_id = store.get_sector_id_at(conn, *address)
    finally:
        conn.close()
    sector_args = single = None
    if existing_id is None:
        # Density-driven, from the address's own profile.
        single = _default_generation_args(config=mysql_config)
        single.density = single.num_systems = None
        single.link_later = True            # the run links its sectors at its end
        sector_args = batch_density.resolve(single, address, position_pc)
    sectors = [sector_args] if sector_args is not None else []
    what = _format_address(address)
    around = None
    if args.radius_pc is not None:
        around = _neighborhood_after(args, edge_pc, mysql_config, batch_density, [(address, position_pc)], args.radius_pc)
        sectors += around[0]
        what += f" and the sectors {_region_text(args)} of it"
    run_common._check_estimate(args, sectors, what, progress)
    # The backfill waits for the end of the run, with its own bar, instead of stalling this one at 0 of 1 (PERF.28).
    designation = provisional_sector_designation(*address)

    def saved_line(follower):
        if not follower.reported:
            log.normal(f"Sector [{designation}] at {_format_address(address)} already existed (sector_id={existing_id}).")
            return
        [saved] = follower.reported
        log.normal(
            f"Saved sector '{saved['name']}' [{designation}] at {_format_address(address)} "
            f"(sector_id={saved['sector_id']})."
        )
        run_sector._log_summary(saved["summary"])

    label = f"Sector ({_format_address(address)})"
    batch = [(address, position_pc, sector_args, "")] if single is not None else []
    followers = [SectorsFollower(args, "sectors", label, edge_pc, batch, bar_description=label if batch else None,
                                 base_args=single, collect=True, finish=saved_line)]
    if around is not None:
        followers.append(_neighborhood_follower(args, edge_pc, [address], *around[1:]))
    return followers


AddressPlan = collections.namedtuple("AddressPlan", "what total shells addresses")
"""The sectors of a shell, band, span or course run: `what` (the words for them), `total` (how many addresses),
`shells` (the shells they lie in, to look up what is generated there) and `addresses` (a callable giving a fresh
iterable of `(shell, band, wedge)` each time, in order)."""


def _quiet(_text):
    """The `say` of a cut that leaves the telling to the run's starter."""


def _span_of(args):
    """The `Span` of a span run: the one the command line built, or one rebuilt from `--shells/--bands/--wedges` (a
    cut reads its arguments back from plain JSON)."""
    span = getattr(args, "span", None)
    if span is None:
        ranges = {name: None if getattr(args, name, None) is None else parse_range(getattr(args, name), f"--{name}")
                  for name in ("shells", "bands", "wedges")}
        span = Span(**ranges)
    return span


def _address_plan(args, edge_pc, batch_density, say=log.normal):
    """
    The sectors of the shell, band, span or course mode the arguments ask for (GEN.81, ADM.29, NAV.48), checked as the run
    checks them (`_require_in_reach`, the size and confirmation rules) and described: the part of `run_addresses`
    that decides where the sectors are, which the run's cut (`expand_generate_sectors`)
    shares with the starter.

    Args:
        say (callable): Where the notes go (`log.normal`, or `_quiet` for a cut whose starter has said them).

    Returns:
        AddressPlan or None: None when there is nothing to generate (said).

    Raises:
        SystemExit: A shell beyond the grid's reach, a course that cannot be plotted, or a run too large that was not
            confirmed.
    """
    if getattr(args, "course", None) is not None:
        return _course_plan(args, say)
    bounds = batch_density.bounds
    if getattr(args, "around", None) is not None:
        return _region_plan(args, edge_pc, bounds, say)
    if getattr(args, "block", None) is not None:
        return _block_plan(args, bounds, say)
    if getattr(args, "span", None) is not None or any(
            getattr(args, name, None) is not None for name in ("shells", "bands", "wedges")):
        span = _span_of(args)
        total = span.count(bounds)
        what = span.describe()
        if total == 0:
            say(f"{what.capitalize()} holds no sector inside the galaxy's outline: nothing to generate.")
            return None
        _refuse_if_very_large(args, what, total, "(--limit N generates only the first N).")
        return AddressPlan(what, total, span._shells(bounds), lambda: span.addresses(bounds))
    shell_index = args.shell
    if args.band is not None:
        band = args.band
        what = f"shell {shell_index} band {band}"
        _require_in_reach(bounds, shell_index, band, what)
        total = spherical.wedge_count(shell_index, band)
        addresses = lambda: ((shell_index, band, wedge) for wedge in range(total))
    else:
        what = f"shell {shell_index}"
        _require_in_reach(bounds, shell_index, 0, what)
        if shell_index not in bounds.high_band:
            say(f"No band of the outline reaches shell {shell_index}: nothing to generate "
                f"(name a band with --band to generate its sectors from their profiles).")
            return None
        span = Span(shells=(shell_index, shell_index))
        total = span.count(bounds)
        addresses = lambda: span.addresses(bounds)
    _refuse_if_very_large(args, what, total, "(--limit N generates only the first N not-yet-generated sectors).")
    return AddressPlan(what, total, [shell_index], addresses)


def _block_plan(args, bounds, say):
    """
    `_address_plan` of block mode (GEN.310): the sectors the galaxy's outline holds inside the block `args.block` names
    (`drill.block_addresses`), shell by shell outward, so the inner sectors are generated first. A block the outline
    does not reach holds none. The key was checked by the command line; a key read back from a run's arguments is
    checked again.
    """
    try:
        block = drill.parse_stage_key(args.block)
    except ValueError as exc:
        log.error(f"--block: {exc}")
        raise SystemExit(1) from exc
    what = f"block {drill.format_drill_key(block)}"

    def addresses():
        return ((shell, band, wedge) for shell, band, wedge in drill.block_addresses(block)
                if bounds.qualifies(shell, band))

    total = sum(1 for _address in addresses())
    if total == 0:
        say(f"{what.capitalize()} holds no sector inside the galaxy's outline: nothing to generate.")
        return None
    _refuse_if_very_large(args, what, total, "(--limit N generates only the first N not-yet-generated sectors).")
    shells = sorted({shell for shell, _band, _wedge in addresses()})
    return AddressPlan(what, total, shells, addresses)


def _refuse_if_very_large(args, what, total, hint):
    """Refuses (`--strict`) or warns about a run of more than `LARGE_RUN_WARNING_THRESHOLD` sectors that has neither
    `--limit` nor `--yes`."""
    if total > LARGE_RUN_WARNING_THRESHOLD and args.limit is None and not args.yes:
        run_common._refuse_or_warn(args, f"{what.capitalize()} holds {total:,} sectors, a very large run {hint}")


def _region_run_plan(args, edge_pc, bounds, skip=()):
    """The `region_run.Plan` of region mode (GEN.1128): the sphere of `--radius-pc` around `--around`, `--grow-shells` more shells, in units of
    `--unit-sectors` (default `region_run.BLOCK_SECTORS`), less the sectors in `skip`."""
    return region_run.plan(tuple(args.around), args.radius_pc, edge_pc, block_sectors=args.unit_sectors or region_run.BLOCK_SECTORS,
                           grow_shells=args.grow_shells or 0, skip=skip, bounds=bounds)


def _region_plan(args, edge_pc, bounds, say):
    """`_address_plan` of region mode (GEN.1128): the sectors of the region in the order of its units, nearest first, inside the galaxy's outline."""
    plan = _region_run_plan(args, edge_pc, bounds)
    x, y, z = args.around
    what = (f"the region of {args.radius_pc:g} pc around ({x:g}, {y:g}, {z:g})"
            + (f" grown by {args.grow_shells} shell(s)" if args.grow_shells else ""))
    if not plan.sectors:
        say(f"{what.capitalize()} holds no sector inside the galaxy's outline: nothing to generate.")
        return None
    _refuse_if_very_large(args, what, plan.sectors, "(--yes goes ahead).")
    def addresses():
        for unit in plan.units:
            items = region_run.unit_items(unit, spherical.wedge_count)
            for place in range(0, len(items), 3):
                yield tuple(items[place:place + 3])

    return AddressPlan(what, plan.sectors, sorted({unit.shell for unit in plan.units}), addresses)


def _expand_region(params):
    """The cut of a region run (GEN.1128): one unit of work per `region_run` unit, nearest first, less the sectors saved now (a rerun after a stop
    builds only what is missing). A unit's cost is its sectors."""
    args = _sector_args_of(params)
    edge_pc = params["edge_pc"]
    bounds = _BatchDensity(store.mysql_config_from_args(args)).bounds
    plan = _region_run_plan(args, edge_pc, bounds, skip=_saved_addresses(args))
    units = []
    for unit in plan.units:
        units.append(({"database": params["database"], "args": params["args"], "edge_pc": edge_pc,
                       "items": region_run.unit_items(unit, spherical.wedge_count), "again": True}, float(unit.count)))
    return commands.Expansion(units, {"sectors": plan.sectors, "units": len(units), "skipped": plan.skipped})


def _course_plan(args, say):
    """`_address_plan` of course mode (NAV.48): the uncharted cells on the hops through unknown space of the NAV course
    between two objects (`course_cells`), inside the galaxy's outline. Past `program_constants.NAV_CHART_CONFIRM_SECTORS`
    sectors it needs `--yes`."""
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
        say(f"Nothing to chart: {why}.")
        return None
    bypass = plan.get("bypass")
    if bypass and bypass["found"]:
        say(f"Note: a route through charted space exists ({bypass['distance_ly']:,.1f} ly against "
            f"{plan['route_distance_ly']:,.1f} ly); charting is optional.")
    if len(cells) > program_constants.NAV_CHART_CONFIRM_SECTORS and not args.yes:
        log.error(f"Charting this course takes {len(cells):,} sectors, more than the {program_constants.NAV_CHART_CONFIRM_SECTORS:,} "
                  f"that need confirming. Run it again with --yes.")
        raise SystemExit(1)
    if plan["outside_galaxy"]:
        say(f"{plan['outside_galaxy']} uncharted cell(s) lie outside the galaxy and are left out.")
    return AddressPlan(what, len(cells), {cell[0] for cell in cells}, lambda: cells)


def _pending_of_plan(args, plan, edge_pc, batch_density):
    """The `_Pending` of an `AddressPlan`: its addresses that are not generated yet, up to `--limit`."""
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        occupied = set(store.get_occupied_addresses(conn, set(plan.shells)))
    finally:
        conn.close()
    return _Pending.of_addresses(args, plan.addresses, occupied, batch_density, edge_pc,
                                 limit=getattr(args, "limit", None))


def _plan_addresses(args, edge_pc, progress):
    """
    The shell, band, span and course modes (GEN.81, ADM.29, NAV.48): every not-yet-generated sector of the addresses
    `_address_plan` names -- a shell of the galaxy's outline, one band (named outright, so generated even outside the
    outline), a range of shells, bands and wedges, or the uncharted cells of a NAV course, up to `args.limit` when set; see
    `_BatchDensity.resolve` for how each sector's density is set. A run of more than `generationStats.SAMPLE_ABOVE`
    sectors is estimated from a sample (`_Pending`, PERF.36). A run of more than `LARGE_RUN_WARNING_THRESHOLD` sectors needs
    `--limit` or `--yes` (a course, `NAV_CHART_CONFIRM_SECTORS` and `--yes`). The stage's cut enumerates the same sectors
    again (`expand_generate_sectors`) when it opens, so they are handed on as the arguments that name them, not as a list.

    Returns:
        list[SectorsFollower]: The stage.

    Raises:
        SystemExit: A shell beyond the grid's reach, a course NAV cannot plot, or a run too large that was neither
            limited nor confirmed.
    """
    batch_density = _BatchDensity(store.mysql_config_from_args(args))
    plan = _address_plan(args, edge_pc, batch_density)
    if plan is None:
        return [SectorsFollower(args, "sectors", "Sectors", edge_pc, count=0)]
    what, total = plan.what, plan.total
    limit = getattr(args, "limit", None)
    pending = _pending_of_plan(args, plan, edge_pc, batch_density)
    batch, estimate_args, sampled_total = _plan_pending(pending)
    skipped, existed = pending.skipped, pending.existed
    run_common._check_estimate(args, estimate_args, what, progress, total=sampled_total)

    def finished(_follower):
        skip_note = f", {skipped} skipped (beyond the grid's reach)" if skipped else ""
        already = f"{total - pending.count - skipped} already existed" if limit is None else f"{existed} already existed"
        log.normal(f"Generated {pending.count} new sector(s) in {what} ({total} total, {already}{skip_note}).")

    return [SectorsFollower(args, "sectors", f"Sectors ({what})", edge_pc, count=pending.count,
                            bar_description=f"Sectors ({what})", planned="addresses", finish=finished)]


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


def _plan_check(config, edge_pc):
    """What the stored plan says to a run of sectors: `(problem, warning)`, the words that refuse a galaxy that was never planned (at
    the standard sector edge, with shells) or `None`, and the heads-up about extending a mixed-version galaxy (DB.7) or `None`."""
    conn = store.get_connection(config)
    try:
        bounds = store.get_galaxy_bounds(conn)
        warning = version_check.galaxy_warning(conn, compute_battery=True)
    finally:
        conn.close()
    if bounds is None:
        return ("The galaxy has never been planned -- run 'planetgen plan' first, so every sector "
                "can be checked against the galaxy's bounds before it is generated."), warning
    if not math.isclose(bounds.edge_pc, edge_pc):
        return (f"The stored plan was built with {bounds.edge_pc:g} pc sectors, not the standard "
                f"{edge_pc:g} pc -- re-run 'planetgen plan'."), warning
    if not bounds:
        return ("The stored plan has no shells at all (nothing in this galaxy shape expects a star "
                "per sector) -- re-run 'planetgen plan' with a different shape."), warning
    return None, warning


def _check_plan(args, edge_pc):
    """Refuses a run of a galaxy that was never planned (at the standard sector edge, with shells); warns when it would
    extend a mixed-version galaxy."""
    problem, warning = _plan_check(store.mysql_config_from_args(args), edge_pc)
    if warning:
        log.normal(f"Warning: {warning}")  # DB.7: a heads-up before it extends a mixed-version galaxy
    if problem:
        log.error(problem)
        raise SystemExit(1)


def _start_clock(args):
    """The database's clock now (`store.database_now`), the start time of a run: the sectors saved since are the run's.
    Says it for the stages that take it as `--since`, unless the run only estimates."""
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        started_at = store.database_now(conn)
    finally:
        conn.close()
    if not getattr(args, "estimate_only", False):
        log.normal(f"Run started at {_since_text(started_at)} on the database's clock (the --since of "
                   f"'planetgen link', 'backfill' and 'settle' for this run's sectors).")
    return started_at


def _deferred_sectors(args, edge_pc):
    """
    The sector stages of a run whose plan cannot be made yet (PERF.214): the same stages `_plan_sectors` gives (`stages.sectors_stages`), each
    told only the arguments that name its sectors. Nothing is read from the database, counted, estimated or asked: the stage works its
    sectors out when it opens, in a worker (`planned_items`), and a stage after it reads what it planned (`wrapper.plan_detail`). The sectors
    are the ones the starter's plan makes. A run built this way is for a caller that has its own estimate (the Generate page's job).

    Returns:
        list[SectorsFollower]: One for each of `stages.sectors_stages`.

    Raises:
        SystemExit: A wedge out of range for its band (said).
    """
    def generated_in(what):
        def finished(follower):
            log.normal(f"Generated {follower.generated} new sector(s) in {what}.")
        return finished

    if any(getattr(args, name, None) is not None for name in ("span", "course", "block", "around")) or (
            args.wedge is None and args.shell is not None):
        return [SectorsFollower(args, "sectors", "Sectors", edge_pc, count=0, bar_description="Sectors", planned="addresses",
                                finish=generated_in("the sectors asked for"))]
    if args.wedge is not None:
        address = (args.shell, args.band, args.wedge)
        try:
            total_wedges = spherical.wedge_count(args.shell, args.band)
        except ValueError as error:
            log.error(f"{_format_address(address)} -- {error}.")
            raise SystemExit(1) from error
        if not (0 <= args.wedge < total_wedges):
            log.error(f"--wedge {args.wedge} is out of range for shell {args.shell} band {args.band} (holds {total_wedges} "
                      f"wedges, 0..{total_wedges - 1}).")
            raise SystemExit(1)
        single = _default_generation_args(config=store.mysql_config_from_args(args))
        single.density = single.num_systems = None
        single.link_later = True            # the run links its sectors at its end
        label = f"Sector ({_format_address(address)})"

        def saved_line(follower):
            designation = provisional_sector_designation(*address)
            if not follower.reported:
                log.normal(f"Sector [{designation}] at {_format_address(address)} already existed.")
                return
            [saved] = follower.reported
            log.normal(f"Saved sector '{saved['name']}' [{designation}] at {_format_address(address)} (sector_id={saved['sector_id']}).")
            run_sector._log_summary(saved["summary"])

        followers = [SectorsFollower(args, "sectors", label, edge_pc, count=0, bar_description=label, base_args=single, collect=True,
                                     planned="single", extra={"address": list(address)}, finish=saved_line)]
        if args.radius_pc is not None:
            followers.append(_deferred_neighborhood(args, edge_pc, [list(address)]))
        return followers
    if args.center_sector is not None:
        region = _region_text(args)
        return [SectorsFollower(args, "sectors", f"Sectors ({region} of sector {args.center_sector})", edge_pc, count=0,
                                bar_description="Sectors (local neighborhood)", planned="neighborhood",
                                finish=generated_in(f"{region} of sector_id={args.center_sector}"))]
    radius_pc = _random_start_radius(args)
    args.radius_pc = radius_pc
    run_seed = draw.getrandbits(128)        # the one draw the starter's plan makes, so the draws after it are the same
    log.debug(f"Random starts are drawn from the run seed {run_seed}.")
    return [SectorsFollower(args, "start", "Starting sector", edge_pc, count=0, planned="random-start", extra={"run_seed": run_seed},
                            announce=_log_start),
            _deferred_neighborhood(args, edge_pc, wrapper.plan_detail("start", "centers"))]


def _deferred_neighborhood(args, edge_pc, centers):
    """The "neighborhood" stage after the stage that makes the sectors at `centers` (a list of addresses, or the deferred parameter that
    reads them from the plan of the stage before it): their `args.radius_pc` neighborhoods, cut when the stage opens."""
    def finished(follower):
        log.normal(f"Generated {follower.generated} new sector(s) within {args.radius_pc:g} pc of the starting sector(s).")

    return SectorsFollower(args, "neighborhood", "Neighborhood", edge_pc, count=0, bar_description="Sectors (local neighborhood)",
                           planned="neighborhood", extra={"centers": centers}, finish=finished)


def _plan_sectors(args, edge_pc, deferred=False):
    """
    The sector stages of a `galaxy` or `sectors` run, in the mode the arguments ask for (span, course, block, single
    address, shell, local neighborhood or random start), as followers: the plan is made, counted, estimated and asked about
    here, under a short display of its own, and nothing is generated.

    Args:
        deferred (bool): Plan nothing here (PERF.214): the stages work their sectors out when they open (`_deferred_sectors`), for a run
            that is built before the galaxy it reads is, or that has its own estimate. It cannot be only estimated.

    Returns:
        list[SectorsFollower]: One for each of `stages.sectors_stages`.

    Raises:
        run_common._EstimateOnly: The run was asked only for its estimate (`--estimate-only`).
    """
    if deferred:
        if getattr(args, "estimate_only", False):
            raise ValueError("a run whose sectors are planned when its stage opens cannot be estimated beforehand")
        return _deferred_sectors(args, edge_pc)
    with run_common._generation_progress(disable=getattr(args, "estimate_only", False)) as progress:
        log.set_console(progress.console)
        try:
            if any(getattr(args, name, None) is not None for name in ("span", "course", "block", "around")):
                return _plan_addresses(args, edge_pc, progress)
            if args.wedge is not None:
                return _plan_single_sector(args, edge_pc, progress)
            if args.shell is not None:
                return _plan_addresses(args, edge_pc, progress)
            if args.center_sector is not None:
                return _plan_neighborhood(args, edge_pc, progress)
            return _plan_random_start(args, edge_pc, progress)
        finally:
            log.reset_console()


def _snap_to_region(args, edge_pc):
    """
    `--region` (GEN.826): moves the center of a single-address or local-neighborhood run to the anchor sector of the region
    (9-block) that holds it and sets the radius to the standard fill sphere (`named_regions.REGION_RADIUS_PC`), so the run is the
    fill of one anchor sphere and charts the whole region. `--center-sector` becomes the anchor's `--shell --band --wedge`
    (generated first when it is empty, as that mode does), so both ways run the same.

    Raises:
        SystemExit: The sector of `--center-sector` does not exist or has never been placed in a galaxy.
    """
    if not getattr(args, "region", False):
        return
    if args.center_sector is not None:
        conn = store.get_connection(store.mysql_config_from_args(args))
        try:
            row = conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors WHERE id = ?",
                               (args.center_sector,)).fetchone()
        finally:
            conn.close()
        if row is None or row["shell_index"] is None:
            log.error(f"sector_id={args.center_sector} " + ("does not exist." if row is None else
                      "has never been placed in a galaxy -- --region needs a sector with a grid address."))
            raise SystemExit(1)
        address = (row["shell_index"], row["band_index"], row["wedge_index"])
    else:
        address = (args.shell, args.band, args.wedge)
    anchor = named_regions.snap_to_anchor(address, edge_pc)
    args.center_sector = None
    args.shell, args.band, args.wedge = anchor
    args.radius_pc = named_regions.REGION_RADIUS_PC
    log.normal(f"The region of {_format_address(address)} is filled about its anchor sector {_format_address(anchor)}, "
               f"{args.radius_pc:g} pc around it.")


def _run_stages(args, title, tail):
    """
    One `galaxy` or `sectors` run: the sector stages (`_plan_sectors`) and the stages of `tail(since)` (followers),
    all as one run in the queue (PERF.172), numbered as `stages.begin` is told.

    Args:
        args (argparse.Namespace): Validated arguments.
        title (str): What the run is called.
        tail (callable): `tail(since)` gives the followers that come after the sectors; `since` is the deferred
            `wrapper.opened_by` parameter of the first sectors stage (PERF.213): the stages that work on the run's sectors
            take the database's clock when it opens, not the moment the run was built.

    Returns:
        datetime: The database's clock when the run started, or `None` when the run only estimated.
    """
    edge_pc = run_common._edge_pc()
    _check_plan(args, edge_pc)
    _snap_to_region(args, edge_pc)
    estimate_only = getattr(args, "estimate_only", False)
    args.link_later = not estimate_only
    started_at = _start_clock(args)
    try:
        sectors = _plan_sectors(args, edge_pc)
    except run_common._EstimateOnly as exc:
        run_common._print_estimate(exc)
        run_common._finish_stats(args)
        return None
    followers = sectors + tail(wrapper.opened_by(sectors[0].key))
    listed = stages.galaxy_stages(args) if title == "galaxy" else stages.sectors_stages(args)
    assert [f.key for f in followers] == [stage.key for stage in listed], "the followers are the stage list"
    stages.begin(listed, args, title)
    followed.run_followed(args, followers, f"{title.capitalize()}")
    stages.finish()
    return started_at


def run_sectors(args):
    """
    `planetgen sectors` (PERF.123): generates sectors in the mode the arguments ask for (span, course, single address,
    shell, local neighborhood or random start) and nothing else. The sectors are saved unlinked: `planetgen link`,
    `backfill` and `settle` do the stages after them, for the sectors saved since the start time the run says.

    Args:
        args (argparse.Namespace): Validated arguments (`command == "sectors"`).
    """
    started_at = _run_stages(args, "sectors", lambda _started_at: [])
    if started_at is not None:
        log.normal(f"The sectors are saved. 'planetgen link --since \"{_since_text(started_at)}\"' links them to "
                   f"their neighbours.")


def run_galaxy(args):
    """
    `planetgen galaxy`: the sectors in the mode the arguments ask for (span, course, single-address, shell,
    local-neighborhood, or random-start), then the link, the scatters of the bright stars and phenomena (GEN.30: they come
    after the sectors, so neither the scatter nor the backfill draws stars for a sector the run filled), the backfill, the
    population and territories, the sector paths (GEN.126) and the block totals -- every one a stage of the one run
    (PERF.172), so one estimate covers it and it is one entry in the queue. First checks the galaxy has been planned at the
    standard sector edge, since every mode validates its addresses against that plan's stored outline before generating
    anything.

    Args:
        args (argparse.Namespace): Validated arguments (`command == "galaxy"`).
    """
    def tail(since):
        stars, _book = run_plan.bright_star_followers(args)
        phenomena, _book = run_plan.phenomena_followers(args)
        return [LinkFollower(args, since), *stars, *phenomena, run_plan.BlockLightsFollower(args),
                BackfillFollower(args, since), run_population.population_follower(args),
                run_population.territories_follower(args), SettleFollower(args, since),
                RunBlocksFollower(args, since)]

    _run_stages(args, "galaxy", tail)
