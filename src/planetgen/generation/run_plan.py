# planetgen/generation/run_plan.py

"""
Run Plan
========

The `plan` command: the galaxy's density skeleton (`galaxy_shape`/
`galaxy_shell`, the outline with one row per shell) that
`ensure_sector_generated` consults to decide, cheaply and exactly,
whether an address is worth generating, and the galaxy-wide bright-star
scatter and its per-sector bands.
"""

import collections
import dataclasses
import math
import queue as queue_module
import threading
import time

from planetgen.queue import commands, progress_rate, redisqueue, units, work as workQueue, wrapper
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import component_pass as componentPass
from planetgen.generation import clusters, satellites
from planetgen.generation import phenomenon_reach, phenomenon_scatter
from planetgen.generation import scatter_companions as scatterCompanions
from planetgen.generation import star_population as starPopulation
from planetgen.generation.star_labels import describe_types, star_label
from planetgen.galaxy import component_field as componentField, nebula_field as nebulaField, seed as galaxySeed
from planetgen.galaxy import shape_draw as shapeDraw
from planetgen.names import naming_key
from planetgen.physics import constants
from planetgen import tuning as program_constants
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw
from planetgen.util import log
from planetgen.galaxy.density import model_terms, relative_density
from planetgen.galaxy import scatter_limits, spherical
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.galaxy.skeleton import (
    DEFAULT_REACH_SHELL, build_shell_extents, expected_system_count_at_density_1, shell_candidate_sector_count,
)
from planetgen.generation import run_common
from planetgen.generation import stages, steps


BACKFILL_CHUNK_SECTORS = 200
"""int: How many sectors one backfill transaction locks, draws and
records (GEN.44)."""


def _outline(conn):
    """The stored outline as `(shell_index, high_band)` pairs, the core first."""
    return [(shell, high) for shell, high, _cells in store.get_galaxy_shells(conn)]


def _scatter_limits_and_seed(conn, skeleton):
    """The galaxy scatter's limits by place (`store.scatter_limits_of`; none before any scatter) and the seed
    every sector's own draws use: the scatter's, else the galaxy seed's
    (GEN.39), else 0 for a galaxy planned before v51."""
    limits = store.scatter_limits_of(conn, skeleton)
    settings = store.bright_star_scatter_settings(conn)
    if settings is not None:
        return limits, settings[1]
    if skeleton.galaxy_seed is None:
        return limits, 0
    return limits, galaxySeed.short_seed(skeleton.galaxy_seed, "bright-stars", "scatter")


def _limits_text(limits):
    """One line naming the limits by place a scatter runs at (the core, baseline and halo tiers)."""
    parts = []
    for name, found in zip(scatter_limits.TIER_NAMES, limits.tiers):
        mass = "no mass pass" if found.mass_limit_sol is None else f"{found.mass_limit_sol:g} solar masses"
        parts.append(f"{name} {found.min_luminosity_sol:g} L_sun and {mass}")
    return "Scatter limits by place: " + "; ".join(parts) + "."


def _effective_level(level, galaxy_level):
    """A sector's bright-star level (GEN.44): its own when it has one (0
    filled, or a backfill's floor), else the galaxy scatter's (`None`
    when no scatter ran: untouched)."""
    if level is not None and level >= 0:
        return level
    return galaxy_level


def _needs_band(level, galaxy_level, floor):
    """Whether a sector at stored `level` lacks stars down to `floor`."""
    current = _effective_level(level, galaxy_level)
    return current is None or current > floor


def _draw_sector_bands(conn, skeleton, addresses, floors, limits, seed, ceiling_cap=None, counts=None,
                       on_sector=None, shells=None):
    """
    Draws each of `addresses` down to its floor (`floors`, a dict or one
    number) from the level it holds (else its band's scatter level, `limits`, a
    `galaxy.scatter_limits.ScatterLimits`: GEN.333), `BACKFILL_CHUNK_SECTORS` at a time:
    locks the chunk's `sector_stats` rows, rechecks each level under the
    lock (another run may have got there first, or filled it), wipes the
    stray stars of a sector with no level at all (a failed run, GEN.44),
    writes the stars and the new levels, and commits. `ceiling_cap` caps
    the ceiling (a band run never draws above the galaxy's old level; a dict gives each address its own).
    `counts`, when given, adds up the stars written per population.
    `on_sector`, when given, is called once per address visited (a
    progress bar's tick, PERF.28). `shells`, when given, adds up the stars
    written per shell and kind: `shells[shell_index]` is a `Counter`
    of `star_label` pairs (GEN.131).

    Returns:
        tuple: `(sectors drawn, stars written)`.
    """
    sectors = stars = 0
    e_value = skeleton.expected_system_count_at_density_1
    for start in range(0, len(addresses), BACKFILL_CHUNK_SECTORS):
        chunk = addresses[start:start + BACKFILL_CHUNK_SECTORS]
        entries = []
        for address in chunk:
            density = max(relative_density(sector_position_pc(*address, skeleton.edge_pc), skeleton.shape), 0.0)
            entries.append((address, density, density * e_value))
        locked = store.lock_sector_stats(conn, entries)
        held_masses = store.sector_bright_masses(conn, chunk, current=True)
        filled = store.get_occupied_addresses(conn, {address[0] for address in chunk})
        rows, new_levels = [], {}
        for address in chunk:
            floor = floors if isinstance(floors, (int, float)) else floors[address]
            level = locked.get(address)
            galaxy_level = limits.level_at(address[0], address[1])
            if address in filled or not _needs_band(level, galaxy_level, floor):
                if on_sector is not None:
                    on_sector()
                continue
            ceiling = _effective_level(level, galaxy_level)
            if ceiling is None:
                # No level here and no finished scatter: any unbuilt star in
                # this cell is left over from a run that failed before it
                # recorded one (TEST.25, GEN.44). This draw has no ceiling,
                # so it is wiped and drawn again whole.
                store.delete_unbuilt_bright_stars(conn, address)
            if ceiling_cap is not None:
                cap = ceiling_cap if isinstance(ceiling_cap, (int, float)) else ceiling_cap[address]
                ceiling = cap if ceiling is None else min(ceiling, cap)
            held = held_masses.get(address)
            # GEN.185: a galaxy whose mass pass placed every heavy star draws the lighter ones here
            # (GEN.187: a sector a mass backfill took lower holds every star from its own mass up).
            mass_limit = limits.mass_limit_at(address[0], address[1])
            mass_cap = mass_limit if held is None else held if mass_limit is None else min(held, mass_limit)
            rows.extend(brightStars.backfill_cells(skeleton.shape, [address], skeleton.edge_pc, e_value, floor,
                                                   ceiling, seed,
                                                   mass_range=None if mass_cap is None else (None, mass_cap)))
            new_levels[address] = floor
            if on_sector is not None:
                on_sector()
        scatterCompanions.insert_bright_stars(conn, rows)
        store.set_sector_bright_levels(conn, new_levels)
        conn.commit()
        if counts is not None:
            for row in rows:
                counts[row[6]] += 1
        if shells is not None:
            for row in rows:
                shells.setdefault(row[0], collections.Counter())[star_label(row[7], row[8])] += 1
        sectors += len(new_levels)
        stars += len(rows)
    return sectors, stars


def _draw_sector_masses(conn, skeleton, targets, limits, seed, on_sector=None, shells=None):
    """
    The mass backfill's draw (GEN.187): takes each sector of `targets` (a
    dict, `(shell, band, wedge)` -> the lightest initial mass to place,
    solar masses) down to its target, `BACKFILL_CHUNK_SECTORS` at a time,
    the way `_draw_sector_bands` does for luminosity: locks the chunk's
    `sector_stats` rows, rechecks each sector under the lock (another run
    may have got there first, or filled it), writes the stars and the new
    mass (`store.set_sector_bright_masses`; the galaxy scatter's level too,
    so the sector counts as having its own state from then on), and
    commits.

    A sector holds its stars from `held` up: its own backfill's mass, else
    its band's mass limit (`limits`, GEN.333; `None`: none placed by mass). The stars drawn
    are those born in [target, held) that are dimmer than the sector's
    luminosity level (`_effective_level`): the brighter ones were placed
    already. A sector at or below its target, or filled, is skipped.

    Returns:
        tuple: `(sectors drawn, stars written)`.
    """
    sectors = stars = 0
    e_value = skeleton.expected_system_count_at_density_1
    edges = brightStars.mass_band_edges([found.mass_limit_sol for found in limits.tiers or ()])
    addresses = sorted(targets)
    for start in range(0, len(addresses), BACKFILL_CHUNK_SECTORS):
        chunk = addresses[start:start + BACKFILL_CHUNK_SECTORS]
        entries = []
        for address in chunk:
            density = max(relative_density(sector_position_pc(*address, skeleton.edge_pc), skeleton.shape), 0.0)
            entries.append((address, density, density * e_value))
        locked = store.lock_sector_stats(conn, entries)
        masses = store.sector_bright_masses(conn, chunk, current=True)
        filled = store.get_occupied_addresses(conn, {address[0] for address in chunk})
        rows, states = [], {}
        for address in chunk:
            target = targets[address]
            level = locked.get(address)
            galaxy_level = limits.level_at(address[0], address[1])
            held = masses.get(address, limits.mass_limit_at(address[0], address[1]))
            if address in filled or level == 0 or (held is not None and held <= target):
                if on_sector is not None:
                    on_sector()
                continue
            ceiling = _effective_level(level, galaxy_level)
            if ceiling is None and held is None:
                # Nothing was placed here by any pass: an unbuilt star in this
                # cell is left over from a run that failed before it recorded
                # a level (TEST.25, GEN.44), so it is wiped and drawn again.
                store.delete_unbuilt_bright_stars(conn, address)
            rows.extend(brightStars.backfill_mass_cells(skeleton.shape, [address], skeleton.edge_pc, e_value,
                                                        target, held, ceiling, seed, edges))
            own = galaxy_level if galaxy_level is not None and (level is None or level < 0) else None
            states[address] = (target, own)
            if on_sector is not None:
                on_sector()
        scatterCompanions.insert_bright_stars(conn, rows)
        store.set_sector_bright_masses(conn, states)
        conn.commit()
        if shells is not None:
            for row in rows:
                shells.setdefault(row[0], collections.Counter())[star_label(row[7], row[8])] += 1
        sectors += len(states)
        stars += len(rows)
    return sectors, stars


def build_skeleton(args):
    """
    Builds and persists the galaxy's skeleton: `galaxy_shape` (the shape
    parameters, calibration constant, edge length and outer shell -- one
    singleton row) and `galaxy_shell` (the galaxy's outline: one row per
    shell, from the core out, holding the highest band that shell reaches --
    see `galaxy.skeleton.build_shell_extents`). Together they are the bounds
    every generation path checks first. The edge is always the standard
    `program_constants.DEFAULT_SECTOR_EDGE_PC`.

    No sector content or individual addresses are stored -- that stays
    lazy (`ensure_sector_generated`). Re-running replaces the whole
    skeleton; a full Milky-Way-scale build takes a few milliseconds.

    Args:
        args (argparse.Namespace): Parsed arguments.

    Returns:
        dict: `outer_shell_index`, `shell_count`, `total_candidate_sectors`,
              `elapsed_s`, `edge_confirmed`.
    """
    edge_pc = run_common._edge_pc()
    e_value = expected_system_count_at_density_1(program_constants.DEFAULT_SECTOR_EDGE_LY)
    threshold_rho = 1.0 / e_value

    reach_shell_index = getattr(args, "reach_shells", DEFAULT_REACH_SHELL)
    try:
        spherical.check_reach(reach_shell_index)
    except ValueError as exc:
        log.error(f"--reach-shells {reach_shell_index}: {exc}")
        raise SystemExit(1)

    mysql_config = store.mysql_config_from_args(args)
    requested_seed = getattr(args, "seed", None)
    conn = store.get_connection(mysql_config)
    try:
        stored_seed = store.get_galaxy_seed(conn)
        if (requested_seed is not None and stored_seed is not None and requested_seed != stored_seed
                and store.has_galaxy_sectors(conn)):
            log.error(f"The galaxy already has sectors made from seed {galaxySeed.format_seed(stored_seed)}; "
                      f"a different --seed needs a wiped galaxy.")
            raise SystemExit(1)
        # A new outline invalidates every pre-placed bright star (schema v43).
        store.clear_bright_stars(conn)
        # The seed is settled first: the shape (GEN.401), the globular clusters (GEN.159, GEN.160) and the satellites
        # (GEN.237) all come from it, and the components go in the outline's bound, so they are placed before the
        # outline is walked.
        galaxy_seed = requested_seed or stored_seed or galaxySeed.new_seed()
        shape, drawn = shapeDraw.plan_shape(args, galaxy_seed)
        log.normal(
            f"Building skeleton: edge_pc={edge_pc:g} stellar_mass={drawn.chain.stellar_mass_solar:.3g} ({drawn.stage}-type) "
            f"disk_scale_length_pc={shape.disk_scale_length_pc:.0f} disk_scale_height_pc={shape.disk_scale_height_pc:.0f} "
            f"bulge_scale_radius_pc={shape.bulge_scale_radius_pc:.0f} bulge_to_total={drawn.bulge_to_total:.3f} "
            f"arm_count={shape.arm_count} pitch_angle_deg={drawn.pitch_angle_deg:.1f} "
            f"k_norm={shape.k_norm:.4f} threshold_rho={threshold_rho:.6f}"
        )
        scale = shape.disk_scale_length_pc / program_constants.HALO_SCALE_LENGTH_PC
        found = clusters.generate_clusters(
            galaxy_seed, epoch_unix=program_constants.SCATTER_PHASE_ORIGIN_UNIX, scale=scale)
        dwarfs = satellites.place_satellites(
            galaxy_seed, program_constants.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH * shape.disk_scale_length_pc,
            model_terms(shape)["solar_angle_rad"], epoch_unix=program_constants.SCATTER_PHASE_ORIGIN_UNIX, scale=scale)
        store.replace_galaxy_components(conn, found, kind="cluster")
        store.replace_galaxy_components(conn, dwarfs, kind="satellite")
        conn.commit()
        shape = shape._replace(components=componentField.ComponentField(found + dwarfs, edge_pc, e_value))
    finally:
        conn.close()
    log.normal(f"Drew {len(found):,} globular clusters holding {sum(c.expected_systems for c in found):,.0f} systems "
               f"and placed {len(dwarfs):,} satellite galaxies holding {sum(c.expected_systems for c in dwarfs):,.0f}.")

    t0 = time.perf_counter()
    extents, outer_shell_index, edge_confirmed = build_shell_extents(
        shape, edge_pc, threshold_rho, max_shell=min(args.max_shell, reach_shell_index),
    )
    elapsed = time.perf_counter() - t0

    store.replace_galaxy_shells(extents, config=mysql_config)
    galaxy_seed = store.save_galaxy_shape(
        shape, edge_pc=edge_pc, outer_shell_index=outer_shell_index,
        expected_system_count_at_density_1=e_value, config=mysql_config, galaxy_seed=galaxy_seed,
        reach_shell_index=reach_shell_index, cell_star_ceiling=getattr(args, "cell_star_ceiling", None),
        fill_luminosity_cap_sol=getattr(args, "fill_luminosity_cap", None),
        core_min_luminosity_sol=getattr(args, "core_min_luminosity", None),
        core_mass_limit_sol=getattr(args, "core_mass_limit", None),
        halo_min_luminosity_sol=getattr(args, "halo_min_luminosity", None),
    )
    core_floor = getattr(args, "core_min_luminosity", None)
    core_mass = getattr(args, "core_mass_limit", None)
    halo_floor = getattr(args, "halo_min_luminosity", None)
    limits = store.GalaxySkeletonInfo(
        shape, edge_pc, outer_shell_index, e_value, None, core_min_luminosity_sol=core_floor,
        core_mass_limit_sol=core_mass,
        halo_min_luminosity_sol=program_constants.HALO_MIN_LUMINOSITY_SOL if halo_floor is None else halo_floor)
    log.normal(
        f"Scatter limits by place: core luminosity floor {scatter_limits.core_min_luminosity_sol(limits, None):,.0f} "
        f"L_sun{'' if core_floor is not None else ' (the baseline plus a tenth)'}, core mass limit "
        f"{scatter_limits.core_mass_limit_sol(limits, None):g} M_sun{'' if core_mass is not None else ' (the baseline plus a tenth)'}, "
        f"halo luminosity floor {scatter_limits.halo_min_luminosity_sol(limits):,.0f} L_sun.")
    if requested_seed is not None:
        source = "from --seed"
    elif stored_seed is not None:
        source = "kept from the earlier plan"
    else:
        source = "drawn at random"
    log.normal(f"Galaxy seed {galaxySeed.format_seed(galaxy_seed)} ({source}).")
    _draw_naming_key(mysql_config, galaxy_seed, new_seed=stored_seed != galaxy_seed)

    return {
        "outer_shell_index": outer_shell_index,
        "reach_shell_index": reach_shell_index,
        "shell_count": len(extents),
        "total_candidate_sectors": shell_candidate_sector_count(extents),
        "elapsed_s": elapsed,
        "edge_confirmed": edge_confirmed,
    }


def _draw_naming_key(mysql_config, galaxy_seed, new_seed):
    """
    GEN.70: the galaxy's naming key, drawn from its seed into the control
    database. A new seed draws a new key; planning again over the same seed
    keeps the key an admin may have changed. Without a control database the
    plan still succeeds and the key is drawn the first time an admin asks.
    """
    try:
        conn = store.get_control_connection(store.control_mysql_config(mysql_config))
    except Exception as exc:  # noqa: BLE001 -- no control database yet
        log.normal(f"No naming key drawn: the control database can't be opened ({exc}). Run update.sh.")
        return None
    try:
        key = naming_key.draw(conn, mysql_config.database, galaxy_seed, replace=new_seed)
    except Exception as exc:  # noqa: BLE001 -- control schema older than v9
        log.normal(f"No naming key drawn: {exc}. Run update.sh.")
        return None
    finally:
        conn.close()
    log.normal(f"Naming key {key}.")
    return key


def scatter_bright_stars(args, passes=("mass", "luminosity")):
    """
    Pre-places the galaxy's bright stars on the stored plan into
    `bright_stars`, replacing any earlier scatter, in the three star passes
    of GEN.185 (pass 1, the phenomena above the mass limit, and pass 5, the
    other phenomena, are `scatter_phenomena`'s):

    2. the mass pass: every star born at or above the mass limit
       (`--phenomenon-min-mass`, the stellar mass limit, one of 8 to 20 solar masses), whatever its
       luminosity (`brightStars.scatter` with `mass_range=(limit, None)`);
    3. the marks: each sector holding a pass 2 star at least as bright as
       the luminosity floor is marked (`store.bright_star_marked_addresses`);
    4. the luminosity pass: every star born lighter than the mass limit
       and at least `--bright-star-min-luminosity` bright, skipping the
       marked sectors, which already hold a star that bright.

    Each run of shells is one task per pass, with a progress bar. Sectors already
    filled are always left out (GEN.30): the scatter never adds stars to a
    generated sector. `--force` is still accepted and does nothing more.

    `passes` (GEN.196) names the passes to redo. Both (the default) replace
    the whole scatter. One alone clears and rewrites only that pass's stars
    (`store.clear_bright_star_pass`) and keeps the other's: the luminosity
    pass alone works at the mass limit already stored, the mass pass alone at
    the luminosity floor already stored. With no scatter stored yet a single
    pass becomes both.

    Returns:
        dict: `counts` (per population), `total` and `elapsed_s`.
    """
    mysql_config = store.mysql_config_from_args(args)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError("The galaxy's skeleton has never been built -- run 'planetgen plan' first.")
        extents = _outline(conn)
        filled = store.filled_sector_addresses(conn)
        if filled:
            log.normal(f"Leaving out the {len(filled):,} sectors already filled.")
        min_luminosity_sol = float(args.bright_star_min_luminosity)
        mass_limit = _stellar_mass_limit(args, conn)
        stored_limit = store.bright_star_mass_limit(conn)
        stored_settings = store.bright_star_scatter_settings(conn)
        passes = set(passes)
        if passes != {"mass", "luminosity"} and (stored_limit is None or stored_settings is None):
            log.normal("No complete scatter is stored yet, so both passes run.")
            passes = {"mass", "luminosity"}
        elif passes == {"mass"} and not math.isclose(mass_limit, stored_limit):
            log.normal(f"The mass limit changed from {stored_limit:g} to {mass_limit:g} solar masses, so the "
                       "luminosity pass is redone with it.")
            passes = {"mass", "luminosity"}
        seed = _bright_star_seed(skeleton, "scatter")
        mass_seed = _bright_star_seed(skeleton, "mass-scatter")
        if passes == {"mass", "luminosity"}:
            store.clear_bright_stars(conn)
        elif passes == {"luminosity"}:
            mass_limit = stored_limit
            store.clear_bright_star_pass(conn, "luminosity", scatter_limits.ScatterLimits.of(
                skeleton, stored_settings[0], stored_limit))
        else:
            min_luminosity_sol = float(stored_settings[0])
            store.clear_bright_star_pass(conn, "mass", scatter_limits.ScatterLimits.of(
                skeleton, min_luminosity_sol, stored_limit))
        # GEN.333: the baseline floor and mass limit, and the core's and halo's settings, give each band its tier's.
        limits = scatter_limits.ScatterLimits.of(skeleton, min_luminosity_sol, mass_limit)
        log.normal(_limits_text(limits))

        t0 = time.perf_counter()
        counts = {population: 0 for population in brightStars.POPULATIONS}
        massive = 0
        if "mass" in passes:
            stages.enter("mass")
            counts = _scatter_mass(args, mysql_config, min_luminosity_sol, mass_limit, mass_seed)
            massive = sum(counts.values())
        else:
            stages.skip("mass", "not chosen to be redone")
        if "luminosity" in passes:
            marked = store.bright_star_marked_addresses(conn, limits)
            log.normal(f"{len(marked):,} sectors hold a star born at their band's mass limit or more and at "
                       "least its luminosity floor bright; the luminosity pass skips them.")
            stages.enter("luminosity")
            light = _scatter_shells(args, mysql_config, skeleton, extents, set(filled) | marked, min_luminosity_sol,
                                    seed, "Bright stars", tiered=scatter_limits.luminosity_pass(limits))
            for population, count in light.items():
                counts[population] += count
        else:
            stages.skip("luminosity", "not chosen to be redone")
        store.record_bright_star_scatter(conn, min_luminosity_sol, seed, mass_limit)
        conn.commit()
    finally:
        conn.close()
    elapsed = time.perf_counter() - t0
    total = sum(counts.values())
    log.normal(
        f"Placed {total:,} bright stars ({massive:,} born at {mass_limit:g} solar masses or more, the rest at least "
        f"{min_luminosity_sol:g} L_sun bright) in {elapsed:.1f}s: "
        + ", ".join(f"{count:,} {population}" for population, count in counts.items()) + "."
    )
    log.debug(f"bright-star scatter: {total} stars, seed {seed}, threshold {min_luminosity_sol:g} L_sun, "
              f"mass limit {mass_limit:g} solar masses")
    return {"counts": counts, "total": total, "elapsed_s": elapsed}


def redo_scatters(args):
    """
    Redoes the scatters `--redo-scatters` names, each at this run's settings
    (GEN.196), in the order of the plan's own run: the phenomena (neutron
    stars and black holes at the compact-object limit, everything else, and
    the nucleus), then the star passes. A star pass alone clears and rewrites
    only its own stars (`scatter_bright_stars(passes=...)`); the phenomena
    scatter clears and rewrites the phenomena. A scatter not named is left
    exactly as it was, and its stage is listed as skipped.
    """
    chosen = list(args.redo_scatters)
    if "phenomena" in chosen:
        with workQueue.job_node("phenomena", "Phenomena"):
            scatter_phenomena(args)
    star_passes = [name for name in ("mass", "luminosity") if name in chosen]
    if star_passes:
        with workQueue.job_node("bright-stars", "Bright stars"):
            scatter_bright_stars(args, passes=star_passes)


def _bright_star_seed(skeleton, address):
    """
    The 63-bit seed of a bright-star scatter (`"scatter"`) or band, which
    its shells and the backfill's blocks derive their own streams from:
    the top 63 bits of the unit seed `bright-stars:<address>` (GEN.39),
    since `galaxy_shape.bright_star_seed` is a `BIGINT UNSIGNED`. A galaxy
    with no seed (planned before schema v51) draws one from the run's
    `random` stream.
    """
    if skeleton.galaxy_seed is None:
        return draw.getrandbits(63)
    return galaxySeed.short_seed(skeleton.galaxy_seed, "bright-stars", address)


class _ShellTracker:
    """
    The bright-star scatter's progress (PERF.4, PERF.9): the main bar
    counts expected work (`brightStars.shell_weight`, in stars) rather
    than tasks, credited as each task in progress reports its stars,
    so its ETA follows the galaxy's shape; and while tasks finish
    slower than one per `SLOW_TASK_SECONDS`, a second bar under it
    shows the stars of the tasks being drawn, done of their estimate,
    with their own ETA. The second bar goes again once tasks finish
    faster than one per `FAST_TASK_SECONDS` (the gap keeps it from
    flashing on and off near the line). A task is a run of consecutive
    shells (`shell_chunks`); `names` says which ones, `sizes` how many.
    Called from the run's own thread and the progress channel's
    (`_drain_channel`), so every change holds `lock`.
    """

    SLOW_TASK_SECONDS = 30.0
    FAST_TASK_SECONDS = 20.0

    def __init__(self, progress, label, weights, clock=time.monotonic, prior=None, bar=None, sizes=None, names=None):
        self.progress = progress
        self.bar = bar
        self.label = label
        self.weights = weights
        self.sizes = sizes if sizes is not None else {key: 1 for key in weights}
        self.names = names if names is not None else {key: str(key) for key in weights}
        self.clock = clock
        self.lock = threading.Lock()
        self.in_flight = {}
        self.credited = {}
        self.finished = set()
        self.done_weight = 0.0
        self.done_shells = 0
        self.task_rate = progress_rate.DecayingRate(clock=clock)
        self.last_done = clock()
        self.detail = None
        self.task = None
        if bar is None:
            self.task = progress.add_task(self._description(), total=max(sum(weights.values()), 1.0), percent=True,
                                          prior=prior)
            progress.main_task = self.task
        else:
            bar.update(description=self._description())     # a step (UX.84): it draws itself when it is long

    def _description(self):
        return f"{self.label} ({self.done_shells:,} of {sum(self.sizes.values()):,} shells)"

    def task_progress(self, key, done, estimate):
        """A task in progress has drawn `done` of about `estimate` stars.
        A report that arrives after its task finished (reports come
        through the channel's thread) is ignored, or that task would be
        counted twice (PERF.23)."""
        with self.lock:
            if key in self.weights and key not in self.finished:
                self.in_flight[key] = (done, estimate)
                self._refresh()

    def task_done(self, key):
        with self.lock:
            if key in self.finished:
                return
            self.finished.add(key)
            self.in_flight.pop(key, None)
            self.credited.pop(key, None)
            self.done_weight += self.weights.get(key, 0.0)
            self.done_shells += self.sizes.get(key, 1)
            self.last_done = self.clock()
            self.task_rate.add(1)
            self._refresh()

    def slow(self):
        """Whether tasks are finishing slower than one per
        `SLOW_TASK_SECONDS` (or, once the second bar shows, not yet
        faster than one per `FAST_TASK_SECONDS`)."""
        limit = self.FAST_TASK_SECONDS if self.detail is not None else self.SLOW_TASK_SECONDS
        since = self.clock() - self.last_done
        rate = self.task_rate.rate
        return since > limit or (rate is not None and rate < 1.0 / limit)

    def _refresh(self):
        partial = 0.0
        for key, (done, estimate) in self.in_flight.items():
            share = min(done / estimate, 1.0) if estimate > 0 else 0.0
            # Never back: an estimate that grows doesn't take credit away.
            credit = max(self.credited.get(key, 0.0), share * self.weights[key])
            self.credited[key] = credit
            partial += credit
        if self.bar is not None:
            self.bar.update(completed=self.done_weight + partial, description=self._description())
        else:
            self.progress.update(self.task, completed=self.done_weight + partial, description=self._description())
        if self.in_flight and self.slow():
            done = sum(item[0] for item in self.in_flight.values())
            estimate = sum(max(item[0], item[1]) for item in self.in_flight.values())
            keys = sorted(self.in_flight)
            names = ", ".join(self.names[key] for key in keys[:4]) + (", ..." if len(keys) > 4 else "")
            description = f"  Shell{'s' if sum(self.sizes[key] for key in keys) > 1 else ''} {names}: stars"
            if self.detail is None:
                self.detail = self.progress.add_task(description, total=max(estimate, 1.0), completed=done)
                self.progress.detail_task = self.detail
            # A new total writes the progress file at once.
            self.progress.update(self.detail, completed=done, total=max(estimate, 1.0), description=description)
        elif self.detail is not None and not self.slow():
            self.progress.remove_task(self.detail)
            self.progress.detail_task = None
            self.detail = None


class _DirectChannel:
    """The progress channel of a scatter run with one worker: its tasks
    are drawn in this process, so reports go straight to the tracker."""

    def __init__(self, tracker):
        self.tracker = tracker

    def put(self, item):
        self.tracker.task_progress(*item)


def _drain_channel(channel, tracker, stop):
    """Hands the workers' task reports to the tracker until `stop`."""
    while True:
        try:
            item = channel.get(timeout=0.25)
        except queue_module.Empty:
            if stop.is_set():
                return
            continue
        except (EOFError, OSError):
            return
        tracker.task_progress(*item)


CHANNEL_INTERVAL_SECONDS = 0.25
"""float: How often a worker drawing a task reports its stars (PERF.4)."""

CHUNKS_PER_WORKER = 32
"""int: About how many tasks of equal expected work a scatter is cut into for each worker, so that no worker is left
with a big one at the end while the others sit idle."""

MAX_SHELLS_PER_CHUNK = 64
"""int: The most consecutive shells one task draws, so the progress moves steadily at the thin edge of the galaxy,
where thousands of shells hold almost nothing."""


def shell_chunks(rows, workers):
    """
    Cuts the outline into tasks of consecutive shells with about equal
    expected work. A shell draws from a stream of its own, so how the shells
    are grouped changes nothing in what is drawn; only how well the workers
    are kept busy.

    Args:
        rows (list): `(shell_index, high_band, weight)`, innermost shell first.
        workers (int): How many workers share the run.

    Returns:
        list[list]: The rows of each task, the heaviest task first.
    """
    total = sum(weight for _shell, _high, weight in rows)
    target = total / max(workers, 1) / CHUNKS_PER_WORKER
    chunks, current, weight = [], [], 0.0
    for row in rows:
        current.append(row)
        weight += row[2]
        if weight >= target or len(current) >= MAX_SHELLS_PER_CHUNK:
            chunks.append(current)
            current, weight = [], 0.0
    if current:
        chunks.append(current)
    return sorted(chunks, key=lambda chunk: -sum(row[2] for row in chunk))


def static_shell_chunks(rows):
    """
    The units of a scatter on the new queue (PERF.97): the outline cut into runs of consecutive shells of about
    equal expected work, from the rows alone (`units.static_chunks`), so the same galaxy is the same units at any
    worker count. The heaviest unit first, as `shell_chunks`.

    Args:
        rows (list): `(shell_index, high_band, weight)`, innermost shell first.
    """
    chunks = units.static_chunks(rows, max_rows=MAX_SHELLS_PER_CHUNK)
    return sorted(chunks, key=lambda chunk: -sum(row[2] for row in chunk))


def _shell_names(shells):
    """`shells` (a chunk's rows) as "7" or "100-131"."""
    first, last = shells[0][0], shells[-1][0]
    return str(first) if first == last else f"{first}-{last}"


def _cells_in(shells):
    """Sectors in `shells` (`(shell_index, high_band)` pairs) inside the outline (none in a shell the outline does not
    reach, `high_band` -1)."""
    return sum(spherical.cells_in_band_range(shell, high) for shell, high in shells if high >= 0)


def log_shells(label, shells):
    """`_log_shell` for each shell of a sector-by-sector draw (`_draw_sector_bands`'s `shells`), innermost first."""
    for shell_index in sorted(shells):
        _log_shell(label, shell_index, shells[shell_index])


def _log_shell(label, shell_index, types):
    """One line saying how many stars a shell was given, by kind (GEN.131);
    nothing for a shell that drew none (GEN.79's note says how many did)."""
    total = sum(types.values())
    if total:
        log.normal(f"{label}, shell {shell_index}: added {total:,} stars: {describe_types(types)}.")


def _log_shell_range(label, shells, types):
    """`_log_shell` for a task of several shells: the stars it added, by kind, over the shells it covers."""
    total = sum(types.values())
    if total:
        log.normal(f"{label}, shell{'s' if len(shells) > 1 else ''} {_shell_names(shells)}: added {total:,} stars: "
                   f"{describe_types(types)}.")


COMPONENT_TASK_SHELL = -1
"""int: The shell the components' task is listed under in a scatter's tasks (it has none of its own)."""


def _galaxy_components(mysql_config):
    """The galaxy's stored clusters and satellites (`component_pass.COMPONENT_KINDS`), none for a galaxy without."""
    conn = store.get_connection(mysql_config)
    try:
        return [component for component in store.get_galaxy_components(conn)
                if component.kind in componentPass.COMPONENT_KINDS]
    finally:
        conn.close()


@dataclasses.dataclass
class _ScatterShape:
    """
    How a shell scatter is cut (`_plan_shape`): the tasks (`tasks`, unit -> its `(shell_index, high_band, weight)`
    rows), each one's expected work (`weights`, in stars), shell count (`sizes`) and name (`names`), the stars expected
    by shell (`expected`; the components' task is listed under `COMPONENT_TASK_SHELL`) and the keys of the far shells'
    and the components' tasks. What a process following the run needs: it travels as JSON (`to_json`).
    """

    tasks: dict
    weights: dict
    sizes: dict
    names: dict
    expected: dict
    far_key: int = None
    component_key: int = None
    first_far: int = 0
    reach_shell: int = 0
    reach_on: bool = False
    e_value: float = 0.0
    band_share: float = 0.0
    component_count: int = 0

    def to_json(self):
        order = range(len(self.tasks))
        return {"tasks": [[list(row) for row in self.tasks[key]] for key in order],
                "weights": [self.weights[key] for key in order], "sizes": [self.sizes[key] for key in order],
                "names": [self.names[key] for key in order],
                "expected": [[shell, stars] for shell, stars in self.expected.items()],
                "far_key": self.far_key, "component_key": self.component_key, "first_far": self.first_far,
                "reach_shell": self.reach_shell, "reach_on": self.reach_on, "e_value": self.e_value,
                "band_share": self.band_share, "component_count": self.component_count}

    @classmethod
    def from_json(cls, data):
        order = range(len(data["tasks"]))
        return cls(
            tasks={key: [tuple(row) for row in data["tasks"][key]] for key in order},
            weights={key: data["weights"][key] for key in order}, sizes={key: data["sizes"][key] for key in order},
            names={key: data["names"][key] for key in order},
            expected={shell: stars for shell, stars in data["expected"]}, far_key=data["far_key"],
            component_key=data["component_key"], first_far=data["first_far"], reach_shell=data["reach_shell"],
            reach_on=data["reach_on"], e_value=data["e_value"], band_share=data["band_share"],
            component_count=data["component_count"])

    def log_line(self, label):
        return (f"{label}: about {round(sum(self.expected.values())):,} stars expected in "
                f"{sum(self.sizes.values()):,} shells, {self.reach_shell:,} shells being the grid's reach"
                + (f", and {self.component_count:,} components." if self.component_count else "."))


def _band_fractions(tiered, min_luminosity_sol, max_luminosity_sol, mass_range):
    """A scatter pass's expected share of the stars by population (`fractions`) and its sum for the baseline
    (`band_share`), for a tiered pass or a single band."""
    if tiered is None:
        fractions = brightStars.band_fractions(min_luminosity_sol, max_luminosity_sol, mass_range)
        return fractions, sum(fractions.values())
    fractions = brightStars.tier_fractions_of(tiered)
    shares = [sum(item.values()) for item in fractions]
    return fractions, shares[scatter_limits.BASELINE] or max(shares)


def _plan_shape(skeleton, extents, fractions, band_share, components, tiered, workers=None):
    """
    Cuts a shell scatter into tasks (`_ScatterShape`): the outline's shells and the shells past it (GEN.322: the
    pass covers the grid's whole reach; the shells past the outline are drawn one by one while they hold stars to
    speak of, and the faint ones further out together, as one task), then the components' task (GEN.324). With
    `workers` the cut follows the worker count (the old queue); without, it is the static cut of the new queue
    (`static_shell_chunks`), the same at any worker count.
    """
    e_value = skeleton.expected_system_count_at_density_1
    rows, expected = [], {}
    outer_shell = max((shell for shell, _high in extents), default=0)
    reach_on = program_constants.SCATTER_COVERS_REACH
    reach_shell = getattr(skeleton, "reach_shell_index", None) if reach_on else None
    reach_shell = outer_shell if reach_shell is None else max(reach_shell, outer_shell)
    first_beyond = outer_shell + 1
    first_far, _means = brightStars.far_shells(
        skeleton.shape, first_beyond, reach_shell, skeleton.edge_pc, e_value, fractions)
    beyond_shells = [(shell, -1) for shell in range(first_beyond, first_far)]
    estimates = brightStars.reach_estimates(
        skeleton.shape, list(extents) + beyond_shells, skeleton.edge_pc, e_value, fractions) if reach_on else {}
    for shell_index, high_band in list(extents) + beyond_shells:
        weight, expected[shell_index] = brightStars.shell_weight(
            skeleton.shape, shell_index, high_band, skeleton.edge_pc, e_value, fractions,
            beyond=estimates[shell_index] if reach_on else None)
        rows.append((shell_index, high_band, weight))
    tasks = dict(enumerate(static_shell_chunks(rows) if workers is None else shell_chunks(rows, workers)))
    weights = {key: sum(row[2] for row in chunk) for key, chunk in tasks.items()}
    sizes = {key: len(chunk) for key, chunk in tasks.items()}
    names = {key: _shell_names(chunk) for key, chunk in tasks.items()}
    far_key = component_key = None
    if first_far <= reach_shell:
        far_key = len(tasks)
        far_expected = brightStars.far_expected_stars(
            skeleton.shape, first_far, reach_shell, skeleton.edge_pc, e_value, fractions)
        expected[first_far] = far_expected
        tasks[far_key] = [(first_far, reach_shell, far_expected)]
        weights[far_key] = far_expected + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS
        sizes[far_key] = reach_shell - first_far + 1
        names[far_key] = f"{first_far}-{reach_shell}"
    # GEN.324: the globular clusters and satellite galaxies place their own stars, as one more task.
    if components:
        component_key = len(tasks)
        component_expected = componentPass.expected_stars(components, fractions, tiered)
        expected[COMPONENT_TASK_SHELL] = component_expected
        tasks[component_key] = [(COMPONENT_TASK_SHELL, COMPONENT_TASK_SHELL, component_expected)]
        weights[component_key] = component_expected + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS
        sizes[component_key] = 0
        names[component_key] = "components"
    return _ScatterShape(tasks, weights, sizes, names, expected, far_key, component_key, first_far, reach_shell,
                         reach_on, e_value, band_share, len(components or []))


@dataclasses.dataclass
class _ScatterPass:
    """What a shell scatter draws with besides its shape: the galaxy's skeleton, the sectors to leave out (`filled`),
    the pass's limits and seed and the components (`None` for a galaxy without)."""

    skeleton: object
    filled: set
    min_luminosity_sol: float
    max_luminosity_sol: float
    seed: int
    mass_range: tuple
    tiered: object
    components: list
    skip_by_shell: dict = dataclasses.field(init=False)

    def __post_init__(self):
        self.skip_by_shell = collections.defaultdict(set)
        for address in self.filled:
            self.skip_by_shell[address[0]].add(address)

    def payload(self, shape, key):
        """What one task of the pass is given besides the database, the shape and the progress channel."""
        chunk = shape.tasks[key]
        special = key in (shape.far_key, shape.component_key)
        return {
            "shells": [] if special else [(shell, high) for shell, high, _weight in chunk],
            "far": (shape.first_far, shape.reach_shell) if key == shape.far_key else None, "reach": shape.reach_on,
            "components": self.components if key == shape.component_key else None,
            "expected_by_shell": ({shape.first_far: shape.expected[shape.first_far]} if key == shape.far_key else
                                  {COMPONENT_TASK_SHELL: shape.expected[COMPONENT_TASK_SHELL]}
                                  if key == shape.component_key else
                                  {shell: shape.expected[shell] for shell, _high, _weight in chunk}),
            "edge_pc": self.skeleton.edge_pc, "expected": shape.e_value,
            "min_luminosity_sol": self.min_luminosity_sol, "max_luminosity_sol": self.max_luminosity_sol,
            "seed": self.seed, "skip": (set(self.filled) if key == shape.component_key else
                                        set().union(*(self.skip_by_shell.get(shell, set())
                                                      for shell, _high, _weight in chunk))),
            "mass_range": self.mass_range, "tiered": self.tiered,
        }

    def unit_args(self, shape, key, mysql_config):
        """The arguments of the unit that runs task `key` (`scatter_shells_unit`): `payload` as plain JSON, and the
        database's name, host, port and user -- never its password. The components are read by the unit from the
        database."""
        payload = self.payload(shape, key)
        return {
            "database": wrapper.database_args(mysql_config), "key": key,
            "shells": [list(pair) for pair in payload["shells"]],
            "far": None if payload["far"] is None else list(payload["far"]),
            "reach": payload["reach"], "components": payload["components"] is not None,
            "expected_by_shell": [[shell, expected] for shell, expected in payload["expected_by_shell"].items()],
            "skip": sorted(list(address) for address in payload["skip"]),
            "edge_pc": payload["edge_pc"], "expected": payload["expected"],
            "min_luminosity_sol": payload["min_luminosity_sol"], "max_luminosity_sol": payload["max_luminosity_sol"],
            "seed": payload["seed"],
            "mass_range": None if payload["mass_range"] is None else list(payload["mass_range"]),
            "tiered": scatter_limits.tiered_to_json(payload["tiered"]),
        }


class _ShellBars:
    """
    The progress of a shell scatter: the console's progress display, and, once the scatter's shape is known (`open`),
    the bar counting each task's expected work (PERF.9) with its tracker (`_ShellTracker`). Closes the bar, as a
    success only when the scatter did not fail.
    """

    def __init__(self, args, label, stats_kind, workers):
        self.args, self.label, self.stats_kind, self.workers = args, label, stats_kind, workers
        self.bar = self.tracker = self.progress = None
        self._display = run_common._generation_progress()

    def __enter__(self):
        self.progress = self._display.__enter__()
        log.set_console(self.progress.console)
        return self

    def open(self, shape):
        # PERF.33: the bar starts from the stars a second this server recorded for this many workers.
        prior = run_common._generation_stats(self.args).pool_rate(self.stats_kind, self.workers)
        self.bar = steps.Step(self.label, self.stats_kind, max(sum(shape.weights.values()), 1.0), args=self.args,
                              progress=self.progress, workers=self.workers, percent=True, record=False,
                              prior=prior).__enter__()
        self.tracker = _ShellTracker(self.progress, self.label, shape.weights, prior=prior, bar=self.bar,
                                     sizes=shape.sizes, names=shape.names)

    def __exit__(self, exc_type, exc, traceback):
        try:
            if self.bar is not None:
                self.bar.close(success=exc_type is None)
        finally:
            log.reset_console()
            suppressed = self._display.__exit__(exc_type, exc, traceback)
        return suppressed


class _ScatterTally:
    """What a shell scatter adds up as its tasks finish: the stars written per population, the shells that drew
    some, the stars with companions, and the speed stats of each task (PERF.10, PERF.53)."""

    def __init__(self, args, label, shape, stats_kind, workers, tracker):
        self.args, self.label, self.shape, self.stats_kind, self.workers = args, label, shape, stats_kind, workers
        self.tracker = tracker
        self.counts = {population: 0 for population in brightStars.POPULATIONS}
        self.drew = set()
        self.systems = 0

    def on_done_for(self, key):
        shape, chunk = self.shape, self.shape.tasks[key]

        def task_done(result, seconds, _weight):
            task_counts = result["counts"]
            types = collections.Counter({(name, plural): count for name, plural, count in result["types"]})
            if key == shape.component_key:
                if sum(types.values()):
                    log.normal(f"{self.label}, components: added {sum(types.values()):,} stars: "
                               f"{describe_types(types)}.")
            else:
                _log_shell_range(self.label, chunk, types)
            for population, count in task_counts.items():
                self.counts[population] += count
            self.systems += result.get("systems", 0)
            self.tracker.task_done(key)
            stars = sum(task_counts.values())
            for shell_index, shell_stars in result["shells"]:
                if shell_stars:
                    self.drew.add(shell_index)
            if key in (shape.far_key, shape.component_key):
                return
            cells = _cells_in([(shell, high) for shell, high, _weight in chunk])
            wanted = sum(shape.expected[shell] for shell, _high, _weight in chunk)
            density = wanted / (shape.e_value * shape.band_share * cells) if shape.band_share and cells else 0.0
            if stars:   # PERF.53: a task that drew nothing says nothing about how long a star takes
                run_common._generation_stats(self.args).record(self.stats_kind, density, seconds, systems=stars,
                                                                stars=stars, workers=self.workers)
        return task_done

    def finish(self):
        """Says where the stars landed and notes the stage; returns the stars written per population."""
        # Every shell is walked, but only the ones that drew a star count as holding any (GEN.79).
        shells = sum(self.shape.sizes.values())
        if self.drew:
            log.normal(f"{self.label}: stars landed in {len(self.drew):,} of {shells:,} shells "
                       f"(shells {min(self.drew)} to {max(self.drew)}); {self.systems:,} of them have companions.")
        else:
            log.normal(f"{self.label}: no stars landed in any of the {shells:,} shells.")
        stages.note(shells=shells, shells_modified=len(self.drew), objects=sum(self.counts.values()))
        return self.counts


def _scatter_shells(args, mysql_config, skeleton, extents, filled, min_luminosity_sol, seed, label,
                    max_luminosity_sol=None, mass_range=None, stats_kind="scatter-bright", task_kind="bright-stars",
                    tiered=None):
    """
    Draws and writes the bright stars of every shell through the work
    queue, a task per run of consecutive shells (`shell_chunks`), with a
    progress bar: every star at or above `min_luminosity_sol`, or only those
    below `max_luminosity_sol` too (one band of a staged scatter), and born
    in `mass_range` (initial solar masses `(low, high)`, an end `None` for no
    limit; GEN.185). Cells in `filled` (the sectors already filled, and the
    ones the mass pass marked) are left out. `tiered` (a `scatter_limits.TieredPass`, GEN.333) draws each band at its
    tier's limits instead, and `min_luminosity_sol`, `max_luminosity_sol` and `mass_range` are then not read.

    The bar counts each task's expected work (PERF.9,
    `brightStars.shell_weight`, worked out first in a second or two),
    credited star by star as the workers report (PERF.4, `_ShellTracker`,
    which also adds the second bar for slow tasks). Each task's time
    goes into the speed stats as a `stats_kind` task (PERF.10) and has its own work queue task kind `task_kind`
    (PERF.54). The massive-star pass is `_scatter_mass`'s, on the new queue.

    Returns:
        dict: Stars written per population.
    """
    fractions, band_share = _band_fractions(tiered, min_luminosity_sol, max_luminosity_sol, mass_range)
    workers = run_common._worker_count(args)
    components = _galaxy_components(mysql_config)
    shape = _plan_shape(skeleton, extents, fractions, band_share, components, tiered, workers)
    drawing = _ScatterPass(skeleton, filled, min_luminosity_sol, max_luminosity_sol, seed, mass_range, tiered,
                           components)
    log.normal(shape.log_line(label))
    with _ShellBars(args, label, stats_kind, workers) as bars:
        bars.open(shape)
        tally = _ScatterTally(args, label, shape, stats_kind, workers, bars.tracker)
        stop = threading.Event()
        channel = drain = None
        try:
            with run_common._work_queue(args, label) as queue:
                if queue.parallel:
                    channel = queue.channel("scatter-progress")
                    drain = threading.Thread(target=_drain_channel, args=(channel, bars.tracker, stop),
                                             name="scatter-progress", daemon=True)
                    drain.start()
                else:
                    channel = _DirectChannel(bars.tracker)
                queue.submit_each(task_kind, [
                    (f"shell {shape.names[key]}", {
                        "mysql_config": mysql_config, "shape": skeleton.shape, "key": key,
                        **drawing.payload(shape, key), "channel": channel,
                    }, shape.weights[key], tally.on_done_for(key))
                    for key in shape.tasks], _scatter_shells_task)
        finally:
            stop.set()
            if drain is not None:
                drain.join(timeout=5)
            if isinstance(channel, redisqueue.Channel):
                channel.close()
            run_common._finish_stats(args)
    return tally.finish()


def expand_scatter_mass(params):
    """
    The `scatter-mass` command's stage cut into its units (PERF.118): the massive-star pass over the stored galaxy,
    a unit per run of shells (`static_shell_chunks`, the same at any worker count), from `params` and the database
    alone, so the same params and galaxy give the same units.

    Args:
        params (dict): `database` (`wrapper.database_args`), `min_luminosity_sol` (the baseline luminosity floor),
            `mass_limit_sol` (the stellar mass limit) and `seed` (the pass's).

    Returns:
        commands.Expansion: The units, and the `_ScatterShape` as JSON.
    """
    named = params["database"]
    mysql_config = store.MySQLConfig(host=named["host"], port=named["port"], user=named["user"],
                                     database=named["database"])
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError("The galaxy's skeleton has never been built -- run 'planetgen plan' first.")
        extents = _outline(conn)
        filled = store.filled_sector_addresses(conn)
    finally:
        conn.close()
    limits = scatter_limits.ScatterLimits.of(skeleton, params["min_luminosity_sol"], params["mass_limit_sol"])
    floor = starPopulation.MASS_PASS_MIN_LUMINOSITY_SOL
    tiered = scatter_limits.mass_pass(limits, floor)
    fractions, band_share = _band_fractions(tiered, floor, None, None)
    components = _galaxy_components(mysql_config)
    shape = _plan_shape(skeleton, extents, fractions, band_share, components, tiered)
    drawing = _ScatterPass(skeleton, filled, floor, None, params["seed"], None, tiered, components)
    return commands.Expansion([(drawing.unit_args(shape, key, mysql_config), shape.weights[key]) for key in shape.tasks],
                              shape.to_json())


def _scatter_mass(args, mysql_config, min_luminosity_sol, mass_limit, seed, label="Massive stars"):
    """
    The massive-star pass (GEN.185 pass 2) on the new queue (PERF.97, PERF.118): the wrapper (`queue/wrapper.py`)
    selects the `scatter-mass` command from this run's steps, records the run, cuts the pass into units
    (`expand_scatter_mass`) and queues them, and this process follows the units' results and progress with the same
    bar as the other scatters until all are done. `--workers` (or `PLANETGEN_WORKERS`) is how many workers the run
    starts, one included.

    Returns:
        dict: Stars written per population.
    """
    workers = run_common._worker_count(args)
    step, listing = wrapper.step_for(args, "mass", label, {
        "min_luminosity_sol": float(min_luminosity_sol), "mass_limit_sol": float(mass_limit), "seed": seed,
        "database": wrapper.database_args(mysql_config)})
    held = {}
    with _ShellBars(args, label, "scatter-mass", workers) as bars:
        def on_plan(_index, plan):
            held["shape"] = _ScatterShape.from_json(plan["detail"])
            log.normal(held["shape"].log_line(label))
            bars.open(held["shape"])
            held["tally"] = _ScatterTally(args, label, held["shape"], "scatter-mass", workers, bars.tracker)

        def on_result(_index, unit, result):
            held["tally"].on_done_for(unit)(result, result["seconds"], held["shape"].weights[unit])

        try:
            wrapper.run(label, [step], on_plan, on_result, lambda _index, unit, done, estimate:
                        bars.tracker.task_progress(unit, done, estimate),
                        workers_count=workers, who=getattr(args, "who", None), database=mysql_config.database,
                        argv=getattr(args, "argv", None), seed=seed, listing=listing)
        finally:
            run_common._finish_stats(args)
    return held["tally"].finish()


def add_bright_star_band(args):
    """
    Lowers the galaxy's star-fill level to `--bright-stars-down-to`:
    keeps every bright star already placed and scatters only the band from
    the new level up to (not including) the stored one, then stores the new
    level. Sectors already filled are left out: their own systems were
    drawn below the old level, so they already hold stars that bright.

    Returns:
        dict or None: `counts` (per population), `total`, `elapsed_s`,
            `from_luminosity_sol` and `to_luminosity_sol`; `None` when there
            was nothing to do (no scatter yet, or the level asked for is not
            below the current one).
    """
    mysql_config = store.mysql_config_from_args(args)
    target = float(args.bright_stars_down_to)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError("The galaxy's skeleton has never been built -- run 'planetgen plan' first.")
        settings = store.bright_star_scatter_settings(conn)
        mass_limit = store.bright_star_mass_limit(conn)
        mass_range = None if mass_limit is None else (None, mass_limit)
        if settings is None:
            log.normal("No bright stars are scattered yet, so there is no level to go below. Run "
                       f"'planetgen plan --bright-stars-only --bright-star-min-luminosity {target:g}' instead.")
            return None
        current, first_seed = float(settings[0]), settings[1]
        # GEN.333: each tier's floor moves with the baseline unless it was set, so each band goes from its tier's new
        # floor up to its old one, and a tier whose floor does not fall adds nothing.
        older = scatter_limits.ScatterLimits.of(skeleton, current, mass_limit)
        newer = scatter_limits.ScatterLimits.of(skeleton, target, mass_limit)
        tiered = scatter_limits.band_pass(newer, older)
        if target >= current:
            log.normal(f"Nothing to add: every star of {current:g} L_sun or more is already placed, "
                       f"and {target:g} is not below that.")
            return None
        extents = _outline(conn)
        filled = store.filled_sector_addresses(conn)
        # Sectors a backfill took to their own level (GEN.44) already hold
        # part of the band: the shell scatter leaves them out and each gets
        # only what it lacks below its level, sector by sector.
        own = store.sector_bright_level_keys(conn)
        skip = set(filled) | set(own)
        # GEN.32: a band run that stopped part way left the band in the
        # shells it finished; this run draws the whole band again, so it
        # starts from none of it.
        stale = store.delete_unfinished_band(conn, older, skip)
        conn.commit()
        if stale:
            log.normal(f"Removed {stale:,} bright stars an unfinished earlier run left below {current:g} L_sun.")
        seed = _bright_star_seed(skeleton, f"band/{target:g}-{current:g}")
        log.normal(f"Adding the bright stars from {target:g} up to {current:g} L_sun"
                   + (f", leaving out {len(filled):,} filled sectors" if filled else "")
                   + (f" and {len(own):,} backfilled sectors" if own else "") + ".")
        t0 = time.perf_counter()
        counts = _scatter_shells(args, mysql_config, skeleton, extents, skip, target, seed,
                                 f"Bright stars {target:g}-{current:g} L_sun", max_luminosity_sol=current,
                                 mass_range=mass_range, tiered=tiered)
        floors = {address: newer.level_at(address[0], address[1]) for address in own}
        topped = sorted(address for address, level in own.items() if address not in filled and level > floors[address]
                        and floors[address] < older.level_at(address[0], address[1]))
        if topped:
            # Its own bar and ETA, like the backfill's (ADM.26): the web
            # page shows the progress file, and this phase can take minutes.
            with run_common._generation_progress() as progress:
                log.set_console(progress.console)
                try:
                    with steps.Step("Topping up backfilled sectors", "topup", len(topped), args=args,
                                    progress=progress) as bar:
                        topped_shells = {}
                        _sectors, stars = _draw_sector_bands(
                            conn, skeleton, topped, floors, older, first_seed,
                            ceiling_cap={address: older.level_at(address[0], address[1]) for address in topped},
                            counts=counts, shells=topped_shells, on_sector=bar.advance)
                finally:
                    log.reset_console()
            log_shells("Topping up backfilled sectors", topped_shells)
            log.normal(f"Topped up {len(topped):,} backfilled sectors with {stars:,} bright stars.")
        # The first scatter's seed stays: it names the galaxy's scatter.
        store.record_bright_star_scatter(conn, target, first_seed)
        conn.commit()
    finally:
        conn.close()
    elapsed = time.perf_counter() - t0
    total = sum(counts.values())
    log.normal(
        f"Added {total:,} bright stars ({target:g} to {current:g} L_sun) in {elapsed:.1f}s: "
        + ", ".join(f"{count:,} {population}" for population, count in counts.items())
        + f". The star-fill level is now {target:g} L_sun."
    )
    log.debug(f"bright-star band: {total} stars, seed {seed}, {target:g} to {current:g} L_sun")
    return {"counts": counts, "total": total, "elapsed_s": elapsed,
            "from_luminosity_sol": current, "to_luminosity_sol": target}


class _ReportChannel:
    """What `_scatter_shells_task` puts its progress on, for a unit on the new queue: `(key, done, estimate)`
    becomes the unit's `report(done, estimate)`."""

    def __init__(self, report):
        self.report = report

    def put(self, item):
        _key, done, estimate = item
        self.report(done, estimate)


def scatter_shells_unit(unit_args, report):
    """
    A `scatter-mass` unit of the new queue (PERF.97, `queue/units.py`; PERF.118, `queue/commands.py`): the shells `_scatter_shells_task` draws and
    writes for a run of shells, from plain JSON arguments. The galaxy's shape comes from the database, which this
    unit opens with the credentials of this worker's settings and the database, host, port and user the run named.

    Returns:
        dict: `_scatter_shells_task`'s, plus `seconds` the unit took.
    """
    started = time.monotonic()
    named = unit_args["database"]
    mysql_config = store.MySQLConfig(host=named["host"], port=named["port"], user=named["user"],
                                     database=named["database"])
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    mass_range = unit_args.get("mass_range")
    far = unit_args.get("far")
    result = _scatter_shells_task({
        "mysql_config": mysql_config, "shape": skeleton.shape, "key": unit_args["key"],
        "shells": [tuple(pair) for pair in unit_args["shells"]],
        "far": None if far is None else tuple(far), "reach": unit_args.get("reach", False),
        "components": _galaxy_components(mysql_config) if unit_args.get("components") else None,
        "expected_by_shell": {shell: expected for shell, expected in unit_args["expected_by_shell"]},
        "edge_pc": skeleton.edge_pc, "expected": skeleton.expected_system_count_at_density_1,
        "min_luminosity_sol": unit_args["min_luminosity_sol"], "max_luminosity_sol": unit_args.get("max_luminosity_sol"),
        "seed": unit_args["seed"], "skip": {tuple(address) for address in unit_args["skip"]},
        "mass_range": None if mass_range is None else tuple(mass_range),
        "tiered": scatter_limits.tiered_from_json(unit_args.get("tiered"), skeleton), "channel": _ReportChannel(report),
    })
    result["seconds"] = time.monotonic() - started
    return result


def _scatter_shells_task(payload):
    """
    A run of consecutive shells of the bright-star scatter -- a work queue
    task (PERF.7): draws each shell (`brightStars.scatter_shell`, its own
    random stream) and writes its stars, committing every 10,000.

    Returns:
        dict: `counts` (stars written per population), `types` (stars
            per kind, `[label, plural, count]` triples, GEN.131),
            `shells` (`[shell_index, stars]` per shell) and `systems` (the
            stars written with members, GEN.262).
    """
    counts = {population: 0 for population in brightStars.POPULATIONS}
    types = collections.Counter()
    per_shell = []
    channel = payload.get("channel")
    key = payload["key"]
    expected_by_shell = payload["expected_by_shell"]
    last = [0.0]
    before = [0.0, 0.0]        # the stars and the expectation of the shells already drawn
    remaining = [sum(expected_by_shell.values())]

    def report(done, estimate):
        # PERF.4: the stars drawn so far, a few times a second.
        now = time.monotonic()
        if channel is not None and now - last[0] >= CHANNEL_INTERVAL_SECONDS:
            last[0] = now
            try:
                channel.put((key, before[0] + done, before[0] + max(estimate, done) + remaining[0]))
            except (EOFError, OSError):
                pass

    conn = store.get_connection(payload["mysql_config"])
    try:
        epoch = scatterCompanions.default_epoch_unix(conn)
        systems = 0
        batch = []
        runs = [(shell_index, brightStars.scatter_shell(
            payload["shape"], shell_index, high_band, payload["edge_pc"],
            payload["expected"], payload["min_luminosity_sol"], payload["seed"], skip_addresses=payload["skip"],
            max_luminosity_sol=payload.get("max_luminosity_sol"), on_progress=report,
            mass_range=payload.get("mass_range"), tiered=payload.get("tiered"), reach=payload.get("reach", False)))
                for shell_index, high_band in payload["shells"]]
        if payload.get("far"):
            # GEN.322: the faint shells far out, drawn together (`brightStars.scatter_far`).
            first_far, last_far = payload["far"]
            runs.append((first_far, brightStars.scatter_far(
                payload["shape"], first_far, last_far, payload["edge_pc"], payload["expected"],
                payload["min_luminosity_sol"], payload["seed"], skip_addresses=payload["skip"],
                max_luminosity_sol=payload.get("max_luminosity_sol"), on_progress=report,
                mass_range=payload.get("mass_range"), tiered=payload.get("tiered"))))
        if payload.get("components"):
            # GEN.324: the clusters' and satellites' stars (`component_pass`).
            fractions, pass_bands = brightStars.pass_of(
                payload.get("tiered"), payload["min_luminosity_sol"], payload.get("max_luminosity_sol"),
                payload.get("mass_range"))
            runs.append((COMPONENT_TASK_SHELL, componentPass.scatter_stars(
                payload["components"], payload["edge_pc"], payload["seed"], fractions, pass_bands,
                tiered=payload.get("tiered"), skip_addresses=payload["skip"])))
        for shell_index, rows in runs:
            remaining[0] -= expected_by_shell[shell_index]
            drawn = 0
            landed = collections.Counter()
            for row in rows:
                counts[row[6]] += 1
                types[star_label(row[7], row[8])] += 1
                drawn += 1
                landed[row[0]] += 1
                batch.append(row)
                if len(batch) >= 10000:
                    systems += scatterCompanions.insert_bright_stars(conn, batch, epoch).systems
                    conn.commit()
                    batch = []
            if (payload.get("far") and shell_index == payload["far"][0]) or shell_index == COMPONENT_TASK_SHELL:
                per_shell.extend([shell, stars] for shell, stars in sorted(landed.items()))
            else:
                per_shell.append([shell_index, drawn])
            before[0] += drawn
        if batch:
            systems += scatterCompanions.insert_bright_stars(conn, batch, epoch).systems
        conn.commit()
    finally:
        conn.close()
    return {"counts": counts, "types": [[label, plural, count] for (label, plural), count in types.items()],
            "shells": per_shell, "systems": systems}


def _phenomenon_scatter_seed(skeleton):
    """
    The 63-bit seed of the phenomenon scatter (GEN.100): the top 63 bits of
    the unit seed `phenomenon-scatter:scatter` (GEN.39), since
    `galaxy_shape.phenomenon_scatter_seed` is a `BIGINT UNSIGNED`. A galaxy
    with no seed draws one from the run's `random` stream.
    """
    if skeleton.galaxy_seed is None:
        return draw.getrandbits(63)
    return galaxySeed.short_seed(skeleton.galaxy_seed, "phenomenon-scatter", "scatter")


def _stellar_mass_limit(args, conn):
    """The stellar mass limit (GEN.167, GEN.183, GEN.195): `--phenomenon-min-mass`,
    else the one already stored with the star scatter, else
    `tuning.PHENOMENON_MIN_MASS_SOLAR`. It governs stars only."""
    value = getattr(args, "phenomenon_min_mass", None)
    if value is None:
        value = store.bright_star_mass_limit(conn)
    return float(program_constants.PHENOMENON_MIN_MASS_SOLAR if value is None else value)


def _compact_mass_limit(args, conn):
    """The lowest mass of a neutron star or black hole the phenomena scatter
    places (GEN.195): `--compact-min-mass` when it is a number, the stellar
    mass limit for "star". Left out, it follows a `--phenomenon-min-mass` given
    in the same run, else keeps the cut already stored with the galaxy, else
    is the stellar limit."""
    value = getattr(args, "compact_min_mass", None)
    if value == program_constants.COMPACT_MIN_MASS_STAR:
        return _stellar_mass_limit(args, conn)
    if value is not None:
        return float(value)
    if getattr(args, "phenomenon_min_mass", None) is None:
        stored = store.phenomenon_scatter_settings(conn)
        if stored is not None and stored[1] is not None:
            return float(stored[1])
    return _stellar_mass_limit(args, conn)


def _log_phenomena_shells(shells, counts):
    """One line saying how many phenomena a run of shells was given, by kind, as `_log_shell_range` does for stars;
    nothing for a run that drew none."""
    total = sum(counts.values())
    if total:
        log.normal(f"Phenomena, shell{'s' if len(shells) > 1 else ''} {_shell_names(shells)}: placed {total:,}: "
                   + ", ".join(f"{count:,} {kind}" for kind, count in sorted(counts.items())) + ".")


def _scatter_phenomena_shells(args, bar, tasks, weights, task_done, mysql_config, skeleton, e_value, seed,
                              min_mass_solar, filled, far=None, reach=False, components=None):
    """The shells of the phenomenon scatter through the work queue, a task to a run of consecutive shells, the step
    `bar` credited with each task's expected work as it finishes (UX.83, PERF.51). `far` is `(key, first, last)` of
    the task that draws the faint shells far out together (GEN.323), whose chunk is that one entry; `reach` also draws
    each shell's part beyond the outline; `components` is `(key, components)` of the task that places the clusters' and
    satellites' own phenomena (GEN.324). Returns the rows written with members (GEN.262)."""
    shells_done = [0]
    systems = [0]
    far_key = None if far is None else far[0]
    component_key = None if components is None else components[0]
    total_shells = sum(len(chunk) if key not in (far_key, component_key) else
                       (far[2] - far[1] + 1 if key == far_key else 0) for key, chunk in tasks.items())

    def shells_of(key):
        """A task's shells as the `(shell, high, weight)` rows `_shell_names` reads."""
        if key == component_key:
            return []
        return [(shell, -1, 0.0) for shell in range(far[1], far[2] + 1)] if key == far_key else tasks[key]

    def done_with(key):
        def on_done(result, seconds, weight):
            task_counts = result["counts"]
            systems[0] += result.get("systems", 0)
            task_done(task_counts, result["write_seconds"], seconds, weight)
            if key == component_key:
                if sum(task_counts.values()):
                    log.normal(f"Phenomena, components: placed {sum(task_counts.values()):,}: " + ", ".join(
                        f"{count:,} {kind}" for kind, count in sorted(task_counts.items())) + ".")
            else:
                _log_phenomena_shells(shells_of(key), task_counts)
            shells_done[0] += len(shells_of(key))
            bar.update(advance=weights[key],
                       description=f"Phenomena ({shells_done[0]:,} of {total_shells:,} shells)")
        return on_done

    with run_common._work_queue(args, "Phenomena") as queue:
        queue.submit_each("phenomena", [
            ("components" if key == component_key else f"shell {_shell_names(shells_of(key))}", {
                "mysql_config": mysql_config, "shape": skeleton.shape,
                "shells": [] if key == far_key else [(shell, high) for shell, high, _w in chunk],
                "far": (far[1], far[2]) if key == far_key else None, "reach": reach,
                "components": components[1] if key == component_key else None,
                "edge_pc": skeleton.edge_pc, "expected": e_value, "seed": seed,
                "min_mass_solar": min_mass_solar,
                "skip": (set(filled) if key == component_key else
                         {address for address in filled if address[0] in {row[0] for row in shells_of(key)}}),
            }, weights[key], done_with(key))
            for key, chunk in tasks.items()], _phenomenon_shells_task)
    return systems[0]


def _nebula_field_task(payload):
    """
    One first-index slab of field cells of the nebula field table -- a work queue task: draws its kept classes M
    and N clouds from the seed and writes them (DB.24).

    Returns:
        dict: `written` (the rows written), `counts` (rows per class letter) and `shells` (`(lowest, highest)`
            shell any of the clouds reaches, or `None` for a slab without clouds, GEN.222).
    """
    edge_pc = payload["edge_pc"]
    counts = {}
    reach = [None, None]

    def watched(rows):
        for row in rows:
            letter = row[4]
            counts[letter] = counts.get(letter, 0) + 1
            radius_pc = ly_to_pc(row[5])
            distance = math.sqrt(row[6] ** 2 + row[7] ** 2 + row[8] ** 2)
            low = spherical.shell_index_at(max(distance - radius_pc, 0.0), edge_pc)
            high = spherical.shell_index_at(distance + radius_pc, edge_pc)
            reach[0] = low if reach[0] is None else min(reach[0], low)
            reach[1] = high if reach[1] is None else max(reach[1], high)
            yield row

    conn = store.get_connection(payload["mysql_config"])
    try:
        written = store.insert_nebula_field(conn, watched(nebulaField.table_rows(
            payload["galaxy_seed"], payload["shape"], payload["i"], payload["radius_pc"], payload["half_height_pc"])))
        conn.commit()
    finally:
        conn.close()
    outer = payload["outer_shell"]
    shells = None if reach[0] is None else (min(reach[0], outer), min(reach[1], outer))
    return {"written": written, "counts": counts, "shells": shells}


def _shell_range_text(shells):
    """`shells` (lowest, highest) as "shell 2" or "shells 1 to 3"."""
    low, high = shells
    return f"shell {low}" if low == high else f"shells {low} to {high}"


def _log_nebula_slab(i, edge_pc, result):
    """One line saying what a slab of the nebula field was given, by class, and which shells its clouds
    reach, as `_log_phenomena_shells` does for the phenomena; nothing for a slab that drew none (GEN.222)."""
    counts = result["counts"]
    total = sum(counts.values())
    if total:
        kinds = ", ".join(f"{count:,} {program_constants.NEBULA_CLASSES[letter]['name']} (class {letter})"
                          for letter, count in sorted(counts.items()))
        log.normal(f"Nebula field, slab {i} (x {i * edge_pc:,.0f} to {(i + 1) * edge_pc:,.0f} pc): {total:,} clouds"
                   f" in {_shell_range_text(result['shells'])}: {kinds}.")


def build_nebula_field(args, conn, skeleton, extents, progress):
    """
    Draws `nebula_field` (DB.24): every kept cloud of classes M and N of the galaxy's molecular cloud field, from
    the seed, one slab of field cells per task, replacing the earlier table. A nebula already stored from the field
    is linked back to its row afterwards. Returns the rows written.
    """
    galaxy_seed = store.get_galaxy_seed(conn)
    outer_shell = max(shell for shell, _high in extents)
    radius_pc, half_height_pc = nebulaField.field_extent_pc(skeleton.edge_pc, extents)
    edge = program_constants.NEBULA_FIELD_CELL_PC
    slabs = range(math.floor(-radius_pc / edge), math.floor(radius_pc / edge) + 1)
    store.clear_nebula_field(conn)
    written = [0]
    with steps.Step(f"Drawing the nebula field (0 of {len(slabs):,} slabs)", "nebula-field", len(slabs), args=args,
                    progress=progress, workers=run_common._worker_count(args)) as bar:
        done = [0]

        def on_done_slab(i):
            def on_done(result, _seconds, _weight):
                written[0] += result["written"]
                _log_nebula_slab(i, edge, result)
                done[0] += 1
                bar.update(advance=1, description=f"Drawing the nebula field ({done[0]:,} of {len(slabs):,} slabs)")
            return on_done

        mysql_config = store.mysql_config_from_args(args)
        with run_common._work_queue(args, "Nebula field") as queue:
            queue.submit_each("nebula-field", [
                (f"slab {i}", {"mysql_config": mysql_config, "galaxy_seed": galaxy_seed, "shape": skeleton.shape,
                               "i": i, "radius_pc": radius_pc, "half_height_pc": half_height_pc,
                               "edge_pc": skeleton.edge_pc, "outer_shell": outer_shell}, 1.0, on_done_slab(i))
                for i in slabs], _nebula_field_task)
    linked = store.link_nebula_field_stored(conn)
    log.normal(f"Nebula field: {written[0]:,} clouds of classes M and N"
               + (f", {linked:,} already stored." if linked else "."))
    return written[0]


MAP_LUMINOSITY_BATCH_IDS = 200_000
"""int: Ids per batch of the map luminosity fill: one short statement and commit each, so readers never wait."""


def map_luminosity_sql(column="mass_solar"):
    """The SQL expression for `stellar_evolution.map_luminosity_sol` of `column`, from the same pieces."""
    clauses = []
    for piece in program_constants.MAP_MASS_LUMINOSITY_PIECES:
        value = f"{piece['coeff']!r} * POW({column}, {piece['exponent']!r})"
        if piece["max_mass_sol"] is None:
            clauses.append(f"ELSE {value}")
        else:
            clauses.append(f"WHEN {column} < {piece['max_mass_sol']!r} THEN {value}")
    return "CASE " + " ".join(clauses) + " END"


def fill_map_luminosity(args):
    """
    Fills `phenomenon_scatter.map_luminosity_sol` (DB.23, GEN.201) for the rows that have a mass and no value yet,
    in id-range batches of `MAP_LUMINOSITY_BATCH_IDS`, each its own short statement and commit. It only touches
    rows still empty, so a stopped run carries on where it was. The values are `stellar_evolution.map_luminosity_sol`
    of the stored mass, worked out by the same pieces in SQL and stored as FLOAT, so a row agrees with the one a
    new scatter writes to float precision.

    Returns:
        int: The rows filled.
    """
    expression = map_luminosity_sql()
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        top = int(conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM phenomenon_scatter").fetchone()["n"])
        remaining = conn.execute("SELECT COUNT(*) AS n FROM phenomenon_scatter WHERE mass_solar IS NOT NULL"
                                 " AND map_luminosity_sol IS NULL").fetchone()["n"]
        log.normal(f"Map luminosity: {remaining:,} scattered rows have a mass and no value yet.")
        filled = 0
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step("Filling the map luminosity", "map-luminosity", max(top, 1), args=args,
                                progress=progress) as bar:
                    for low in range(1, top + 1, MAP_LUMINOSITY_BATCH_IDS):
                        high = min(low + MAP_LUMINOSITY_BATCH_IDS - 1, top)
                        cursor = conn.execute(
                            f"UPDATE phenomenon_scatter SET map_luminosity_sol = {expression}"
                            " WHERE id BETWEEN ? AND ? AND mass_solar IS NOT NULL AND map_luminosity_sol IS NULL",
                            (low, high))
                        conn.commit()
                        filled += cursor.rowcount
                        bar.update(advance=high - low + 1)
            finally:
                log.reset_console()
        log.normal(f"Filled the map luminosity of {filled:,} scattered rows.")
        return filled
    finally:
        conn.close()


HABITABILITY_BATCH_SYSTEMS = 500
"""int: Systems per batch of the habitability fill: one transaction each, so readers never wait."""


def fill_habitability(args):
    """
    Re-derives the stored habitability numbers of an existing galaxy (GEN.217) from each body's stored conditions
    (`store.refresh_habitability`): the flux swing of a body in a double-star system and the whole score (index,
    six domains, tiers, equipment, note). For a galaxy generated before a change to the rules (GEN.148, GEN.207,
    GEN.215), so the planet chip, the search page and the facet counts agree with the description, which is always
    worked out from the same conditions when the page is drawn. Nothing but those columns is touched. Each batch of
    `HABITABILITY_BATCH_SYSTEMS` systems is its own transaction, in id order, and the recomputation is the same
    whenever it is run, so a stopped run is started again and ends in the same place.

    Returns:
        int: The bodies updated.
    """
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        top = int(conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM star_systems").fetchone()["n"])
        total = conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
        log.normal(f"Habitability: {total:,} systems to re-score.")
        updated = 0
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step("Re-scoring habitability", "habitability", max(top, 1), args=args,
                                progress=progress) as bar:
                    low = 0
                    while True:
                        ids = [row["id"] for row in conn.execute(
                            "SELECT id FROM star_systems WHERE id > ? ORDER BY id LIMIT ?",
                            (low, HABITABILITY_BATCH_SYSTEMS)).fetchall()]
                        if not ids:
                            break
                        for system_id in ids:
                            updated += store.refresh_habitability(conn, system_id)
                        conn.commit()
                        bar.update(advance=ids[-1] - low)
                        low = ids[-1]
            finally:
                log.reset_console()
        log.normal(f"Re-scored the habitability of {updated:,} bodies.")
        return updated
    finally:
        conn.close()


def redo_nebula_field(args):
    """
    Brings an existing galaxy's stored field clouds in line with the current keep rates (GEN.205). Every stored dark
    nebula that is a field cloud (found in its cell's draw by its centre) is checked against `cloud_kept`; the ones
    the rule now rejects are counted per class (the dry run, always done first). With `--redo-nebula-field apply` they are deleted,
    which clears the `inside_nebula_id` of everything inside them (the foreign key sets it to NULL), containment is
    refreshed in the sectors they reached, and `nebula_field` is drawn again. Nothing else is touched: a nebula
    with no cell draw at its centre was not made by the field and stays.

    Returns:
        dict: `stored` and `rejected` (per class, field clouds only) and `deleted` (the count removed).
    """
    apply = args.redo_nebula_field == "apply"
    mysql_config = store.mysql_config_from_args(args)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError("The galaxy's skeleton has never been built -- run 'planetgen plan' first.")
        galaxy_seed = store.get_galaxy_seed(conn)
        draws = {}
        stored, rejected, doomed = collections.Counter(), collections.Counter(), []
        for row in conn.execute(
                "SELECT id, nebula_class, radius_ly, center_x_pc, center_y_pc, center_z_pc FROM nebulae"
                " WHERE nebula_type = ? AND center_x_pc IS NOT NULL ORDER BY id", (nebulaField.FIELD_FAMILY,)).fetchall():
            center = (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])
            index = nebulaField.cell_index(center)
            if index not in draws:
                draws[index] = nebulaField.cell_draw(galaxy_seed, skeleton.shape, index)
            match = next(((n, nebula) for n, nebula, drawn in draws[index]
                          if all(abs(a - b) <= store.FIELD_NEBULA_MATCH_PC for a, b in zip(drawn, center))), None)
            if match is None:
                continue
            stored[row["nebula_class"]] += 1
            if not nebulaField.cloud_kept(galaxy_seed, index, match[0], match[1], center, skeleton.shape):
                rejected[row["nebula_class"]] += 1
                doomed.append((row["id"], center, row["radius_ly"]))
        log.normal("Nebula field cleanup (dry run): " + (
            ", ".join(f"class {letter}: {rejected[letter]:,} of {stored[letter]:,} stored field clouds rejected"
                      for letter in sorted(stored)) or "no stored field clouds") + f"; {len(doomed):,} would be deleted.")
        if not apply:
            log.normal("Nothing was changed. Run again with --redo-nebula-field apply to delete them and draw "
                       "the nebula_field table again.")
            return {"stored": dict(stored), "rejected": dict(rejected), "deleted": 0}
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step("Deleting the rejected field clouds", "nebula-field-clean", max(len(doomed), 1),
                                args=args, progress=progress) as bar:
                    sector_ids = set()
                    for _nebula_id, center, radius_ly in doomed:
                        sector_ids.update(store.sectors_reached_by(conn, center, ly_to_pc(radius_ly)))
                    for start in range(0, len(doomed), 500):
                        ids = [nebula_id for nebula_id, _c, _r in doomed[start:start + 500]]
                        conn.execute(f"DELETE FROM nebulae WHERE id IN ({','.join('?' * len(ids))})", ids)
                        conn.commit()
                        bar.update(advance=len(ids))
                    store.refresh_containment(conn, sorted(sector_ids))
                    conn.commit()
                extents = _outline(conn)
                if extents:
                    build_nebula_field(args, conn, skeleton, extents, progress)
            finally:
                log.reset_console()
        log.normal(f"Deleted {len(doomed):,} stored field clouds the current keep rates reject.")
        return {"stored": dict(stored), "rejected": dict(rejected), "deleted": len(doomed)}
    finally:
        conn.close()


def scatter_phenomena(args):
    """
    Pre-places the galaxy's black holes, neutron stars, planetary nebulae
    and supernova remnants (`phenomenon_scatter.scatter_shell`; neutron
    stars and black holes only from `--compact-min-mass` up, GEN.195), its
    hypervelocity stars and its nucleus into `phenomenon_scatter`, replacing
    any earlier scatter, a run of shells per task. Sectors already filled are left
    out, as the bright-star scatter does (GEN.30), and the nucleus and the
    hypervelocity stars are left out of a cell already filled.

    Returns:
        dict: `counts` (per kind), `total` and `elapsed_s`.
    """
    mysql_config = store.mysql_config_from_args(args)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError("The galaxy's skeleton has never been built -- run 'planetgen plan' first.")
        extents = _outline(conn)
        filled = store.filled_sector_addresses(conn)
        seed = _phenomenon_scatter_seed(skeleton)
        min_mass_solar = _compact_mass_limit(args, conn)
        companions = 0
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step("Clearing the earlier phenomena scatter", "phenomena-clear",
                                max(store.estimated_rows(conn, "phenomenon_scatter"), 1), args=args,
                                progress=progress):
                    store.clear_phenomenon_scatter(conn)
                    # No index to keep up to date while tens of millions of rows go in; built once below (PERF.73).
                    store.drop_phenomenon_scatter_indexes(conn)
            finally:
                log.reset_console()

        t0 = time.perf_counter()
        stages.enter("phenomena")
        counts = {}
        e_value = skeleton.expected_system_count_at_density_1
        # GEN.323: the pass covers the grid's whole reach, as the stars' does (GEN.322): the outline's shells also draw
        # their bands beyond it, the shells past the outline one by one while they hold objects to speak of, and the
        # faint ones further out together, as one task.
        outer_shell = max((shell for shell, _high in extents), default=0)
        reach_on = program_constants.SCATTER_COVERS_REACH
        reach_shell = getattr(skeleton, "reach_shell_index", None) if reach_on else None
        reach_shell = outer_shell if reach_shell is None else max(reach_shell, outer_shell)
        first_far, beyond_shells, estimates = reach_shell + 1, [], {}
        if reach_on:
            first_far, _means = phenomenon_reach.far_shells(
                skeleton.shape, outer_shell + 1, reach_shell, skeleton.edge_pc, e_value, min_mass_solar)
            beyond_shells = [(shell, -1) for shell in range(outer_shell + 1, first_far)]
            estimates = phenomenon_reach.reach_estimates(
                skeleton.shape, list(extents) + beyond_shells, skeleton.edge_pc, e_value, min_mass_solar)
        rows = []
        for shell, high in list(extents) + beyond_shells:
            weight = 0.0
            if high >= 0:
                weight = (phenomenon_scatter.shell_expected(skeleton.shape, shell, high, skeleton.edge_pc, e_value,
                                                           min_mass_solar)
                          + brightStars.BAND_WEIGHT_STARS * (2 * high + 1))
            if reach_on:
                weight += estimates[shell][0] + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS * estimates[shell][1]
            rows.append((shell, high, weight))
        tasks = dict(enumerate(shell_chunks(rows, run_common._worker_count(args))))
        weights = {key: sum(row[2] for row in chunk) for key, chunk in tasks.items()}
        far = None
        if reach_on and first_far <= reach_shell:
            far = (len(tasks), first_far, reach_shell)
            far_weight = phenomenon_reach.far_expected(
                skeleton.shape, first_far, reach_shell, skeleton.edge_pc, e_value, min_mass_solar)
            tasks[far[0]] = [(first_far, reach_shell, far_weight)]
            weights[far[0]] = far_weight + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS
        # GEN.324: the clusters' and satellites' phenomena, one more task.
        galaxy_components = _galaxy_components(mysql_config)
        component_task = None
        if galaxy_components:
            component_task = (len(tasks), galaxy_components)
            component_weight = sum(componentPass.expected_phenomena(component, min_mass_solar)
                                   for component in galaxy_components)
            tasks[component_task[0]] = [(COMPONENT_TASK_SHELL, COMPONENT_TASK_SHELL, component_weight)]
            weights[component_task[0]] = component_weight + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS
        shell_total = len(extents) + len(beyond_shells) + (0 if far is None else reach_shell - first_far + 1)
        log.normal(f"Phenomena: about {round(sum(weights.values())):,} to place in {shell_total:,} shells "
                   f"(neutron stars and black holes from {min_mass_solar:g} solar masses).")
        if filled:
            log.normal(f"Leaving out the {len(filled):,} sectors already filled.")
        landed = [0]

        def task_done(task_counts, write_seconds, seconds, weight):
            for kind, count in task_counts.items():
                counts[kind] = counts.get(kind, 0) + count
            if sum(task_counts.values()):
                landed[0] += 1
            # Recorded in the bar's own units (the task's weight), so the next run's bar can start from it;
            # not for a task that placed nothing (PERF.53). The task's whole time is the bar's pace; the part
            # of it spent writing the rows has a row of its own, the rest being the drawing (PERF.54).
            if sum(task_counts.values()):
                stats = run_common._generation_stats(args)
                stats.record("phenomena", 0.0, seconds, systems=weight, stars=weight,
                             workers=run_common._worker_count(args))
                stats.record("phenomena-write", 0.0, write_seconds, systems=weight, stars=weight,
                             workers=run_common._worker_count(args))

        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step(f"Phenomena (0 of {shell_total:,} shells)", "phenomena",
                                max(sum(weights.values()), 1.0), args=args, progress=progress,
                                workers=run_common._worker_count(args), percent=True, record=False) as bar:
                    companions = _scatter_phenomena_shells(args, bar, tasks, weights, task_done, mysql_config,
                                                           skeleton, e_value, seed, min_mass_solar, filled,
                                                           far=far, reach=reach_on, components=component_task)
            finally:
                log.reset_console()
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                build_nebula_field(args, conn, skeleton, extents, progress)
            finally:
                log.reset_console()
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step("Drawing the special phenomena", "phenomena-special", max(len(extents), 1), args=args,
                                progress=progress):
                    reach = getattr(skeleton, "reach_shell_index", None) if program_constants.SCATTER_COVERS_REACH else None
                    special = list(phenomenon_scatter.special_rows(extents, skeleton.edge_pc, seed, filled, reach))
                with steps.Step("Writing the special phenomena", "phenomena-insert", max(len(special), 1), args=args,
                                progress=progress):
                    store.insert_phenomenon_scatter(conn, special)
                # One index at a time, each worth the table's rows, so the bar moves after every index and its
                # time left comes from the rows per second earlier runs recorded (PERF.101).
                rows = max(store.estimated_rows(conn, "phenomenon_scatter"), 1)
                pending = store.missing_phenomenon_scatter_indexes(conn)
                log.normal(f"Building the phenomena indexes: {len(pending)} indexes over about {rows:,} rows, one "
                           f"at a time; a full galaxy takes a while.")
                with steps.Step("Building the phenomena indexes", "phenomena-index", rows * max(len(pending), 1),
                                args=args, progress=progress) as bar:
                    store.add_phenomenon_scatter_indexes(
                        conn, on_index=lambda name: (bar.advance(rows), log.normal(f"Built index {name}.")))
                with steps.Step("Stamping the hypervelocity stars' epoch", "phenomena-stamp",
                                max(store.estimated_rows(conn, "phenomenon_scatter"), 1), args=args,
                                progress=progress):
                    store.stamp_phenomenon_scatter_epoch(conn)
            finally:
                log.reset_console()
        special_counts = {}
        for row in special:
            label = phenomenon_scatter.class_label(row[3], row[4])
            special_counts[label] = special_counts.get(label, 0) + 1
        for kind, count in special_counts.items():
            counts[kind] = counts.get(kind, 0) + count
        if special_counts:
            log.normal("Special phenomena: " + ", ".join(f"{count:,} {kind}" for kind, count in sorted(special_counts.items())) + ".")
        store.record_phenomenon_scatter(conn, seed, min_mass_solar)
        store.record_phenomenon_scatter_classes(
            conn, {phenomenon_scatter.split_label(label): count for label, count in counts.items()})
        conn.commit()
    finally:
        conn.close()
    elapsed = time.perf_counter() - t0
    total = sum(counts.values())
    stages.note(shells=len(extents), shells_modified=landed[0], objects=total)
    if landed[0]:
        log.normal(f"Phenomena landed in {landed[0]:,} of {len(tasks):,} runs of shells.")
    else:
        log.normal(f"No phenomena landed in any of the {len(extents):,} shells.")
    # Like the star count: only what was created.
    labels = [label for label in list(phenomenon_scatter.EXPECTED_LABELS) + sorted(
        set(counts) - set(phenomenon_scatter.EXPECTED_LABELS)) if counts.get(label)]
    listed = ": " + ", ".join(f"{counts[label]:,} {label}" for label in labels) if labels else ""
    log.normal(f"Placed {total:,} phenomena in {elapsed:.1f}s{listed}; {companions:,} of them have companions.")
    return {"counts": counts, "total": total, "elapsed_s": elapsed}


def _phenomenon_shells_task(payload):
    """
    A run of consecutive shells of the phenomenon scatter -- a work queue
    task: draws each shell (`phenomenon_scatter.scatter_shell`, then its part beyond the outline with
    `phenomenon_reach.scatter_shell` when the scatter covers the reach, GEN.323; or the faint shells far out
    together, `phenomenon_reach.scatter_far`) and writes its rows with their members (`scatter_companions`),
    committing every 10,000.

    Returns:
        dict: `counts` (rows written per kind) and `write_seconds` (the time spent writing them, the rest of the
            task being the drawing, PERF.54).
    """
    counts = {}
    write_seconds = 0.0
    systems = 0
    conn = store.get_connection(payload["mysql_config"])
    try:
        epoch = scatterCompanions.default_epoch_unix(conn)
        batch = []
        shape, edge_pc, expected = payload["shape"], payload["edge_pc"], payload["expected"]
        seed, skip, min_mass = payload["seed"], payload["skip"], payload["min_mass_solar"]
        runs = []
        for shell_index, high_band in payload["shells"]:
            runs.append(phenomenon_scatter.scatter_shell(
                shape, shell_index, high_band, edge_pc, expected, seed, skip_addresses=skip, min_mass_solar=min_mass)
                if high_band >= 0 else iter(()))
            if payload.get("reach"):
                runs.append(phenomenon_reach.scatter_shell(
                    shape, shell_index, high_band, edge_pc, expected, seed, skip_addresses=skip,
                    min_mass_solar=min_mass))
        if payload.get("far"):
            first_far, last_far = payload["far"]
            runs.append(phenomenon_reach.scatter_far(
                shape, first_far, last_far, edge_pc, expected, seed, skip_addresses=skip, min_mass_solar=min_mass))
        if payload.get("components"):
            runs.append(componentPass.scatter_phenomena(
                payload["components"], edge_pc, seed, skip_addresses=skip, min_mass_solar=min_mass, shape=shape))
        for rows in runs:
            for row in rows:
                label = phenomenon_scatter.class_label(row[3], row[4])
                counts[label] = counts.get(label, 0) + 1
                batch.append(row)
                if len(batch) >= 10000:
                    started = time.perf_counter()
                    systems += scatterCompanions.insert_phenomenon_scatter(conn, batch, epoch).systems
                    conn.commit()
                    write_seconds += time.perf_counter() - started
                    batch = []
        started = time.perf_counter()
        if batch:
            systems += scatterCompanions.insert_phenomenon_scatter(conn, batch, epoch).systems
        conn.commit()
        write_seconds += time.perf_counter() - started
    finally:
        conn.close()
    return {"counts": counts, "write_seconds": write_seconds, "systems": systems}


def run_plan(args):
    """
    Builds and persists the galaxy's density skeleton, then scatters its
    bright stars (unless `--no-bright-stars`; `--bright-stars-only` skips
    the rebuild, and `--bright-stars-down-to` adds one dimmer band to the
    stored scatter instead).

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "plan"`).
    """
    stages.begin(stages.plan_stages(args), args, "plan")
    if getattr(args, "bright_stars_down_to", None) is not None:
        stages.enter("band")
        with workQueue.job_node("bright-stars", f"Bright stars down to {args.bright_stars_down_to:g} L_sun"):
            add_bright_star_band(args)
        stages.finish()
        return
    if getattr(args, "fill_map_luminosity", False):
        stages.enter("map-luminosity")
        with workQueue.job_node("map-luminosity", "Map luminosity"):
            fill_map_luminosity(args)
        stages.finish()
        return
    if getattr(args, "fill_habitability", False):
        stages.enter("habitability")
        with workQueue.job_node("habitability", "Habitability"):
            fill_habitability(args)
        stages.finish()
        return
    if getattr(args, "redo_nebula_field", None):
        stages.enter("nebula-field")
        with workQueue.job_node("nebula-field", "Nebula field"):
            redo_nebula_field(args)
        stages.finish()
        return
    if getattr(args, "redo_scatters", None):
        redo_scatters(args)
        stages.finish()
        return
    if getattr(args, "phenomena_only", False):
        with workQueue.job_node("phenomena", "Phenomena"):
            scatter_phenomena(args)
        stages.finish()
        return
    if getattr(args, "bright_stars_only", False):
        with workQueue.job_node("bright-stars", "Bright stars"):
            scatter_bright_stars(args)
        stages.finish()
        return
    stages.enter("skeleton")
    with workQueue.job_node("skeleton", "Galaxy skeleton"):
        summary = build_skeleton(args)
    if summary["shell_count"]:
        shells = f"{summary['shell_count']:,} shells"
    else:
        shells = "no shells (nothing clears the one-star-per-sector threshold)"
    log.normal(
        f"Skeleton built in {summary['elapsed_s']:.2f}s: {shells}, outer edge = shell "
        f"{summary['outer_shell_index']}, ~{summary['total_candidate_sectors']:,} candidate sectors."
    )
    if not summary["edge_confirmed"]:
        log.normal(
            f"WARNING: the galactic plane still qualified at --max-shell ({args.max_shell}) -- the "
            f"galaxy's true edge was not reached. Re-run with a larger --max-shell if these shape "
            f"parameters really do produce a galaxy this large."
        )
    if summary["shell_count"] and not getattr(args, "no_bright_stars", False):
        with workQueue.job_node("phenomena", "Phenomena"):
            scatter_phenomena(args)
        with workQueue.job_node("bright-stars", "Bright stars"):
            scatter_bright_stars(args)
    stages.finish("the galaxy has no shells" if not summary["shell_count"] else None)
