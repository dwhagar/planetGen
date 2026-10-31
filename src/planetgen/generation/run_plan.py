# planetgen/generation/run_plan.py

"""
Run Plan
========

The `plan` command: the galaxy's density skeleton (`galaxy_shape`/
`galaxy_shell`, the outline with one row per shell) that
a run consults to decide, cheaply and exactly,
whether an address is worth generating, and the galaxy-wide bright-star
scatter and its per-sector bands.
"""

import collections
import concurrent.futures
import dataclasses
import math
import multiprocessing
import os
import threading
import time
import types

from planetgen.queue import commands, progress_rate, units, wrapper
from planetgen.db import block_lights, block_stats, containers as containersDb, nebula_blocks, query, storage_checksum, store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import component_pass as componentPass
from planetgen.generation import generator_battery as generatorBattery
from planetgen.generation import containers as generationContainers
from planetgen.generation import phenomenon_reach, phenomenon_scatter
from planetgen.generation import scatter_companions as scatterCompanions
from planetgen.generation import star_population as starPopulation
from planetgen.generation import stats as generationStats
from planetgen.generation.star_labels import describe_types, star_label
from planetgen.galaxy import nebula_field as nebulaField, seed as galaxySeed
from planetgen.galaxy import dials
from planetgen.galaxy import shape_draw as shapeDraw
from planetgen.names import naming_key
from planetgen import tuning as program_constants
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw
from planetgen.util import log
from planetgen.util import step_timer as stepTimer
from planetgen.galaxy.density import relative_density
from planetgen.galaxy import scatter_limits, spherical
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.galaxy.skeleton import (
    DEFAULT_REACH_SHELL, expected_system_count_at_density_1, shell_candidate_sector_count,
)
from planetgen.generation import plan_outline as planOutline
from planetgen.generation import followed, run_common
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
    blocks = {}        # the blocks' candidates, drawn once for the whole draw (GEN.309)
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
                                                   mass_range=None if mass_cap is None else (None, mass_cap),
                                                   memo=blocks))
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
    blocks = {}        # the blocks' candidates, drawn once for the whole draw (GEN.309)
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
                                                        target, held, ceiling, seed, edges, memo=blocks))
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


SKELETON_OPTIONS = ("reach_shells", "cell_star_ceiling", "fill_luminosity_cap", "core_min_luminosity",
                    "core_mass_limit", "halo_min_luminosity", *(row.dest for row in dials.LEVERS))
"""tuple: The options of `planetgen plan` that shape the skeleton, the parameters of the `plan-skeleton` command: the
scan's own and every row of the lever table (`galaxy.dials`)."""


def skeleton_params(args, mysql_config):
    """The parameters of the `plan-skeleton` command for the parsed `args` (plain JSON): the database, the seed the
    user asked for (32 hex digits, or `None`), whether the unit prints its step table (`--timing` or `--debug`, UX.162)
    and each of `SKELETON_OPTIONS` the arguments carry."""
    seed = getattr(args, "seed", None)
    return {"database": wrapper.database_args(mysql_config),
            "galaxy_seed": None if seed is None else galaxySeed.format_seed(seed),
            "timing": bool(getattr(args, "timing", False) or getattr(args, "debug", None) is not None),
            **{name: getattr(args, name) for name in SKELETON_OPTIONS if hasattr(args, name)}}


def expand_plan_skeleton(params):
    """The `plan-skeleton` command's stage: one unit, the whole skeleton, given the command's parameters."""
    return commands.Expansion([(params, 1.0)], {})


def plan_skeleton_unit(params, report):
    """
    The `plan-skeleton` command (PERF.106): builds and persists the galaxy's skeleton: `galaxy_shape` (the shape
    parameters, calibration constant, edge length and outer shell -- one singleton row) and `galaxy_shell` (the
    galaxy's outline: one row per shell, from the core out, holding the highest band that shell reaches -- see
    `galaxy.skeleton.build_shell_extents`). Together they are the bounds every generation path checks first. The edge
    is always the standard `program_constants.DEFAULT_SECTOR_EDGE_PC`.

    No sector content or individual addresses are stored -- that stays lazy (a run makes it). Re-running
    replaces the whole skeleton; a full Milky-Way-scale build takes a few milliseconds.

    Args:
        params (dict): `skeleton_params`.
        report (callable): Unused: the unit is one step.

    Returns:
        dict: `summary` (`outer_shell_index`, `reach_shell_index`, `shell_count`, `total_candidate_sectors`,
            `elapsed_s` and `edge_confirmed`) and `lines`, what the run says about it, for the process that started the
            run to log (this unit may run in a worker); with `params["timing"]` they end with the table of where the
            unit's time went.

    Raises:
        ValueError: `reach_shells` is out of range, or a different seed is asked for in a galaxy that has sectors.
    """
    unit_started = time.perf_counter()
    lines = []
    say = lines.append
    timer = stepTimer.StepTimer()
    edge_pc = run_common._edge_pc()
    e_value = expected_system_count_at_density_1(program_constants.DEFAULT_SECTOR_EDGE_LY)
    threshold_rho = 1.0 / e_value

    reach_shell_index = params.get("reach_shells", DEFAULT_REACH_SHELL)
    try:
        spherical.check_reach(reach_shell_index)
    except ValueError as exc:
        raise ValueError(f"--reach-shells {reach_shell_index}: {exc}") from exc

    mysql_config = _mysql_config_of(params)
    requested_seed = None if params["galaxy_seed"] is None else galaxySeed.parse_seed(params["galaxy_seed"])
    conn = store.get_connection(mysql_config)
    try:
        stored_seed = store.get_galaxy_seed(conn)
        if (requested_seed is not None and stored_seed is not None and requested_seed != stored_seed
                and store.has_galaxy_sectors(conn)):
            raise ValueError(f"The galaxy already has sectors made from seed {galaxySeed.format_seed(stored_seed)}; "
                             f"a different --seed needs a wiped galaxy.")
        # A new outline invalidates every pre-placed bright star (schema v43) and with it the phenomena and the nebula
        # field drawn on the old one: a plan alone leaves them to the scatters, which draw them again (PERF.121).
        with timer.step("clearing the scatters the old outline drew"):
            store.clear_bright_stars(conn)
            store.clear_phenomenon_scatter(conn)
            store.clear_nebula_field(conn)
            block_lights.clear(conn)
            conn.commit()
        # The seed is settled first: the shape (GEN.401), the globular clusters (GEN.159, GEN.160) and the satellites
        # (GEN.237) all come from it, and the components go in the outline's bound, so they are placed before the
        # outline is walked.
        galaxy_seed = requested_seed or stored_seed or galaxySeed.new_seed()
        # A new seed records the generator battery's digest with the galaxy (store.save_galaxy_shape): worked out in a
        # process of its own while the rest of the plan runs, instead of waited for at the end (PERF.176).
        finish_battery = generatorBattery.start_running_digest() if galaxy_seed != stored_seed else None
        shape_args = types.SimpleNamespace(**{name: dials.restored(params.get(name)) for name in SKELETON_OPTIONS})
        built = planOutline.build(shape_args, galaxy_seed, edge_pc, reach_shell_index)
        shape, drawn, found, young, dwarfs = built.shape, built.drawn, built.clusters, built.young, built.satellites
        say(
            f"Building skeleton: edge {edge_pc:g} pc, stellar mass {drawn.chain.stellar_mass_solar:.3g} solar masses "
            f"({drawn.stage}-type), disk scale length {shape.disk_scale_length_pc:.0f} pc, "
            f"disk scale height {shape.disk_scale_height_pc:.0f} pc, "
            f"bulge scale radius {shape.bulge_scale_radius_pc:.0f} pc, bulge-to-total mass ratio {drawn.bulge_to_total:.3f}, "
            f"{shape.arm_count} spiral arms, pitch angle {drawn.pitch_angle_deg:.1f} degrees, "
            f"density normalisation {shape.k_norm:.4f} (a factor, no unit), "
            f"lowest density kept {threshold_rho:.6f} x the local average"
        )
        say(f"{dials.describe_pins(shapeDraw.pins_of(shape_args, shape.disk_scale_length_pc))}.")
        with timer.step("the report lines (pins, halo, outline)"):
            halo_line, halo_warning = shapeDraw.halo_report(drawn)
            made = planOutline.report(built)
        say(f"{halo_line}.")
        if halo_warning:
            say(f"WARNING: {halo_warning}")
        say(f"{made.line()}.")
        for note in made.warnings:
            say(f"WARNING: {note}")
        if made.empty:
            raise ValueError(planOutline.EMPTY_OUTLINE_REASON)
        timer.steps.extend(built.timer.steps)
        with timer.step("storing the clusters and satellites"):
            store.replace_galaxy_components(conn, found + young, kind="cluster")
            store.replace_galaxy_components(conn, dwarfs, kind="satellite")
            conn.commit()
    finally:
        conn.close()
    with timer.step("the system counts of the 'Drew' line"):
        say(f"Drew {len(found):,} globular clusters holding {sum(c.expected_systems for c in found):,.0f} systems, "
            f"{len(young):,} young open clusters holding {sum(c.expected_systems for c in young):,.0f} systems "
            f"and placed {len(dwarfs):,} satellite galaxies holding {sum(c.expected_systems for c in dwarfs):,.0f} systems.")

    extents, outer_shell_index, edge_confirmed = built.extents, built.outer_shell_index, built.edge_confirmed
    elapsed = built.elapsed_s

    with timer.step("storing the shells"):
        store.replace_galaxy_shells(extents, config=mysql_config)
    if finish_battery is not None:
        with timer.step("the generator battery's digest, as far as it was waited for"):
            finish_battery()
    with timer.step("storing the shape"):
        galaxy_seed = store.save_galaxy_shape(
            shape, edge_pc=edge_pc, outer_shell_index=outer_shell_index,
            expected_system_count_at_density_1=e_value, config=mysql_config, galaxy_seed=galaxy_seed,
            reach_shell_index=reach_shell_index, cell_star_ceiling=params.get("cell_star_ceiling"),
            fill_luminosity_cap_sol=params.get("fill_luminosity_cap"),
            core_min_luminosity_sol=params.get("core_min_luminosity"),
            core_mass_limit_sol=params.get("core_mass_limit"),
            halo_min_luminosity_sol=params.get("halo_min_luminosity"),
            pinned_stellar_mass_solar=shapeDraw.pinned_stellar_mass(shape_args, drawn),
            pins=shapeDraw.pins_of(shape_args, shape.disk_scale_length_pc),
        )
    with timer.step("storing the containers"):
        say(_store_containers(mysql_config, galaxy_seed, edge_pc, outer_shell_index))
    core_floor = params.get("core_min_luminosity")
    core_mass = params.get("core_mass_limit")
    halo_floor = params.get("halo_min_luminosity")
    limits = store.GalaxySkeletonInfo(
        shape, edge_pc, outer_shell_index, e_value, None, core_min_luminosity_sol=core_floor,
        core_mass_limit_sol=core_mass,
        halo_min_luminosity_sol=program_constants.HALO_MIN_LUMINOSITY_SOL if halo_floor is None else halo_floor)
    with timer.step("the scatter limits line"):
        say(
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
    say(f"Galaxy seed {galaxySeed.format_seed(galaxy_seed)} ({source}).")
    with timer.step("the naming key"):
        say(_draw_naming_key(mysql_config, galaxy_seed, new_seed=stored_seed != galaxy_seed))
    timer.steps.append(("the rest of the unit (the checks, the report, the limits)",
                        max(0.0, time.perf_counter() - unit_started - timer.total())))
    if params.get("timing"):
        lines.extend(timer.lines("Where the plan's time went"))

    return {"summary": {
        "outer_shell_index": outer_shell_index,
        "reach_shell_index": reach_shell_index,
        "shell_count": len(extents),
        "total_candidate_sectors": shell_candidate_sector_count(extents),
        "elapsed_s": elapsed,
        "edge_confirmed": edge_confirmed,
    }, "lines": lines}


def _store_containers(mysql_config, galaxy_seed, edge_pc, outer_shell_index):
    """
    Writes the galaxy's containers (DB.85): each globular cluster and each superbubble, drawn from the stored shape and
    seed (the shape as the API reads it back, so the stored bubbles are the ones it draws). Returns the line to report.
    """
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        components = store.get_galaxy_components(conn, kind="cluster", profile=generationContainers.GLOBULAR_PROFILE)
        made = generationContainers.build(galaxy_seed, skeleton.shape, edge_pc, (outer_shell_index + 1) * edge_pc,
                                          components)
        containersDb.replace(conn, made)
        conn.commit()
    finally:
        conn.close()
    clusters = sum(1 for c in made if c.kind == generationContainers.GLOBULAR_CLUSTER)
    return f"Stored {clusters:,} globular clusters and {len(made) - clusters:,} superbubbles as containers."


class SkeletonFollower(followed.Follower):
    """
    The galaxy skeleton as a stage of a run (PERF.106, PERF.202): the `plan-skeleton` command builds it in one unit
    (`plan_skeleton_unit`) and this side logs what the unit says. The skeleton stage is only that unit, so the unit that
    completes it stores the stage's row.

    Attributes:
        summary (dict): The build's summary (`plan_skeleton_unit`), once the unit is back.
    """

    key = "skeleton"
    label = "Galaxy skeleton"

    def __init__(self, args):
        super().__init__(args)
        self.summary = None

    def params(self):
        return skeleton_params(self.args, store.mysql_config_from_args(self.args))

    def on_result(self, unit, result):
        for line in result["lines"]:
            log.normal(line)
        self.summary = result["summary"]


def _skeleton_on_queue(args):
    """
    Builds the skeleton on the new queue as a run of the `SkeletonFollower` stage.

    Returns:
        dict: The build's summary (`plan_skeleton_unit`).

    Raises:
        SystemExit: The unit refused its parameters (the reach, or a seed that does not match a galaxy with sectors).
    """
    try:
        return followed.run_planned(args, planned_skeleton(args))
    except units.RunStopped as exc:
        marker = "ValueError: "
        if marker not in str(exc):
            raise
        log.error(str(exc).split(marker, 1)[1])
        raise SystemExit(1) from exc


def planned_skeleton(args):
    """The skeleton build as a `followed.Planned`; its result is the build's summary."""
    follower = SkeletonFollower(args)
    return followed.Planned([follower], "Galaxy skeleton", lambda: follower.summary)


def _draw_naming_key(mysql_config, galaxy_seed, new_seed):
    """
    GEN.70: the galaxy's naming key, drawn from its seed into the control
    database. A new seed draws a new key; planning again over the same seed
    keeps the key an admin may have changed. Without a control database the
    plan still succeeds and the key is drawn the first time an admin asks.

    Returns:
        str: What to say about it: the key, or why none was drawn.
    """
    try:
        conn = store.get_control_connection(store.control_mysql_config(mysql_config))
    except Exception as exc:  # noqa: BLE001 -- no control database yet
        return f"No naming key drawn: the control database can't be opened ({exc}). Run update.sh."
    try:
        key = naming_key.draw(conn, mysql_config.database, galaxy_seed, replace=new_seed)
    except Exception as exc:  # noqa: BLE001 -- control schema older than v9
        return f"No naming key drawn: {exc}. Run update.sh."
    finally:
        conn.close()
    return f"Naming key {key}."


class _BrightStarBook:
    """What the passes of one bright-star scatter add up to: the limits it says before the first, the sectors it leaves
    out, and the line that sums the stars up after the last. The mass limit, the seed and the limits line are worked out when
    the first pass opens (PERF.212) and come with its plan."""

    def __init__(self, min_luminosity_sol):
        self.min_luminosity_sol = min_luminosity_sol
        self.mass_limit = self.seed = self.limits_text = None
        self.followers = []
        self.started = None

    def begin(self, detail):
        """The first pass's plan has arrived: says the limits, and how many sectors the scatter leaves out."""
        if self.started is not None:
            return
        self.started = time.perf_counter()
        self.mass_limit, self.seed, self.limits_text = detail["mass_limit_sol"], detail["seed"], detail["limits_text"]
        log.normal(self.limits_text)
        if detail.get("filled"):
            log.normal(f"Leaving out the {detail['filled']:,} sectors already filled.")

    def summary(self):
        """`{"counts", "total", "elapsed_s"}` of the passes so far: the stars written per population."""
        counts = {population: 0 for population in brightStars.POPULATIONS}
        for follower in self.followers:
            for population, count in (follower.counts or {}).items():
                counts[population] += count
        return {"counts": counts, "total": sum(counts.values()),
                "elapsed_s": 0.0 if self.started is None else time.perf_counter() - self.started}

    def log(self):
        """Sums the scatter up, once its last pass is done."""
        summary = self.summary()
        massive = sum((self.followers[0].counts or {}).values()) if self.followers[0].key == "mass" else 0
        log.normal(
            f"Placed {summary['total']:,} bright stars ({massive:,} born at {self.mass_limit:g} solar masses or more, the rest at least "
            f"{self.min_luminosity_sol:g} L_sun bright) in {summary['elapsed_s']:.1f}s: "
            + ", ".join(f"{count:,} {population}" for population, count in summary["counts"].items()) + "."
        )
        log.debug(f"bright-star scatter: {summary['total']} stars, seed {self.seed}, threshold {self.min_luminosity_sol:g} L_sun, "
                  f"mass limit {self.mass_limit:g} solar masses")


def bright_star_followers(args, passes=("mass", "luminosity"), refuse=None, tag=""):
    """
    The stages that pre-place the galaxy's bright stars on the stored plan into `bright_stars`, in the three star
    passes of GEN.185 (pass 1, the phenomena above the mass limit, and pass 5, the other phenomena, are
    `phenomena_followers`'):

    2. the mass pass: every star born at or above the mass limit (`--phenomenon-min-mass`, the stellar mass limit, one
       of 8 to 20 solar masses), whatever its luminosity (`brightStars.scatter` with `mass_range=(limit, None)`);
    3. the marks: each sector holding a pass 2 star at least as bright as the luminosity floor is marked
       (`store.bright_star_marked_addresses`);
    4. the luminosity pass: every star born lighter than the mass limit and at least `--bright-star-min-luminosity`
       bright, skipping the marked sectors, which already hold a star that bright.

    Each pass is a stage with a unit per run of shells. Sectors already filled are always left out (GEN.30): the
    scatter never adds stars to a generated sector. `passes` names the passes to draw. Both (the default, what `galaxy`
    draws) replace the whole scatter: the first stage empties `bright_stars` when it opens (`prepare_scatter_pass`).
    One alone is what `planetgen scatter --mass` or `--luminosity` draws (PERF.121) after `run_scatter` has found none
    of it stored (Boss's F7: `planetgen clear-scatter` empties a kind first, so no draw replaces another): the
    luminosity pass works at the mass limit the mass pass stored, and needs one (`planetgen scatter --mass` first).
    The last pass stores what the scatter used (`conclude_scatter_pass`), which the stages after it read.

    Args:
        refuse (str, optional): The kind of scatter that is refused when its rows are stored (`_refuse_stored_scatter`, what
            `planetgen scatter` does), when the first pass opens.
        tag (str): Names this scatter among the others of the run that holds it (PERF.203), for what it reads when it opens.

    Returns:
        tuple: `(followers, book)`: a `ScatterFollower` for each pass, and the `_BrightStarBook` that adds them up.

    Raises:
        RuntimeError: The skeleton was never built, or a luminosity pass is asked for with no mass pass stored.
    """
    mysql_config = store.mysql_config_from_args(args)
    min_luminosity_sol = float(args.bright_star_min_luminosity)
    passes = [name for name in ("mass", "luminosity") if name in set(passes)]

    def settings(stage):
        """What the stage's seed, mass limit and clearing are, worked out from the database when it opens (`resolve_scatter_settings`)."""
        return wrapper.deferred("scatter-settings", merge=True, scope="bright-stars", stage=stage,
                                database=wrapper.database_args(mysql_config), min_luminosity_sol=min_luminosity_sol,
                                phenomenon_min_mass=getattr(args, "phenomenon_min_mass", None), passes=passes,
                                run_id=getattr(args, "generation_run_id", None), refuse=refuse, tag=tag)

    book = _BrightStarBook(min_luminosity_sol)
    if "mass" in passes:
        book.followers.append(ScatterFollower(args, "mass", "Massive stars", {
            "min_luminosity_sol": min_luminosity_sol, "settings": settings("mass")}, None, book, last=passes == ["mass"]))
    if "luminosity" in passes:
        book.followers.append(ScatterFollower(args, "luminosity", "Bright stars", {
            "min_luminosity_sol": min_luminosity_sol, "settings": settings("luminosity")}, None, book, last=True))
    return list(book.followers), book


def _bright_star_settings(marker):
    """
    What the passes of one bright-star scatter work from, read from the database when the first of them opens (PERF.212; the
    values `bright_star_followers` read when it built the followers): `{pass: params}` for each pass of the marker. The
    mass limit is `--phenomenon-min-mass`, else the one stored with the star scatter; the luminosity pass alone works at the mass
    limit the mass pass stored and clears only itself, split from the mass pass at the stored limits.

    Raises:
        RuntimeError: The skeleton was never built, or a luminosity pass is asked for with no mass pass stored.
    """
    conn = store.get_connection(_mysql_config_of(marker))
    try:
        if marker.get("refuse"):
            _refuse_stored_scatter(conn, marker["refuse"])
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError(SKELETON_NEVER_BUILT)
        mass_limit = _stellar_mass_limit(types.SimpleNamespace(phenomenon_min_mass=marker["phenomenon_min_mass"]), conn)
        stored_limit = store.bright_star_mass_limit(conn)
        stored_settings = store.bright_star_scatter_settings(conn)
    finally:
        conn.close()
    passes, min_luminosity_sol = marker["passes"], marker["min_luminosity_sol"]
    only_light = passes == ["luminosity"]
    if only_light and stored_limit is None:
        raise RuntimeError("No massive-star pass is stored yet, and the luminosity pass is split from it at its "
                           "mass limit -- run 'planetgen scatter --mass' first.")
    seed = _bright_star_seed(skeleton, "scatter")
    mass_seed = _bright_star_seed(skeleton, "mass-scatter")
    if only_light:
        mass_limit = stored_limit
    record = {"min_luminosity_sol": min_luminosity_sol, "seed": seed, "mass_limit_sol": float(mass_limit),
              "passes": passes, "run_id": marker["run_id"]}
    found = {}
    if "mass" in passes:
        found["mass"] = {"mass_limit_sol": float(mass_limit), "clear": "all", "seed": mass_seed,
                         **({"record": record} if passes == ["mass"] else {})}
    if "luminosity" in passes:
        clear = ({"clear": "luminosity", "clear_floor_sol": float(stored_settings[0] if stored_settings is not None
                                                                 else min_luminosity_sol),
                  "clear_limit_sol": float(stored_limit)} if only_light else {})
        found["luminosity"] = {"mass_limit_sol": float(mass_limit), **clear, "record": record, "seed": seed}
    return found


def resolve_scatter_settings(conn, run_id, marker):
    """
    The resolver of the `scatter-settings` deferred parameter (PERF.212): the params of one stage (`marker["stage"]`) of a
    bright-star scatter (`marker["scope"]` `bright-stars`) or of the phenomenon scatter (`phenomena`), read from the database when
    the first stage of the scatter opens and kept in the run record (`wrapper.cached`), so a stage that comes after another has
    emptied what they were read from, and the stage's conclude, take the same.
    """
    work = _bright_star_settings if marker["scope"] == "bright-stars" else _phenomena_settings
    return wrapper.cached(conn, run_id, wrapper.scatter_settings_name(marker), lambda: work(marker))[marker["stage"]]


def prepare_scatter_pass(params):
    """
    The `prepare` of the shell scatter commands (PERF.172), when the first pass opens: empties what the scatter is
    about to draw again: `params["clear"]` "all" empties `bright_stars` and the scatter's record
    (`store.clear_bright_stars`); "luminosity" empties the luminosity pass alone, split from the mass pass at the limits
    `clear_floor_sol` and `clear_limit_sol` name (`store.clear_bright_star_pass`). No `clear`, nothing is emptied.
    """
    clear = params.get("clear")
    if clear is None:
        return
    conn = store.get_connection(_mysql_config_of(params))
    try:
        if clear == "all":
            store.clear_bright_stars(conn)
        else:
            store.clear_bright_star_pass(conn, "luminosity", scatter_limits.ScatterLimits.of(
                store.get_galaxy_shape(conn), params["clear_floor_sol"], params["clear_limit_sol"]))
    finally:
        conn.close()


def conclude_scatter_pass(params, plan, results, seconds):
    """The `conclude` of the shell scatter commands (PERF.172): after the scatter's last pass, stores the threshold, seed
    and mass limit it used (`store.record_bright_star_scatter`), which the stages after it read (`params["record"]`), and
    the generator epoch the passes were drawn under (`store.record_scatter_epoch`)."""
    record = params.get("record")
    if record is None:
        return {}
    conn = store.get_connection(_mysql_config_of(params))
    try:
        store.record_bright_star_scatter(conn, record["min_luminosity_sol"], record["seed"], record["mass_limit_sol"])
        store.record_scatter_epoch(conn, record["passes"], record["run_id"])
        conn.commit()
    finally:
        conn.close()
    return {}


def scatter_bright_stars(args, passes=("mass", "luminosity"), refuse=None):
    """
    Pre-places the galaxy's bright stars (`bright_star_followers`) as one run: the passes, then the Galaxy Map's block
    light totals made again (`BlockLightsFollower`).

    Returns:
        dict: `counts` (per population), `total` and `elapsed_s`.

    Raises:
        RuntimeError: The skeleton was never built, or a luminosity pass is asked for with no mass pass stored.
    """
    planned = planned_bright_stars(args, passes, refuse)
    if len(passes) == 1:          # the seed the run record names is the one pass's own (a galaxy's two passes have two)
        planned.options["seed"] = _recorded_seed(args, lambda skeleton: _bright_star_seed(
            skeleton, "mass-scatter" if "mass" in passes else "scatter"))
    return _refusals_as_errors(lambda: followed.run_planned(args, planned))


def _recorded_seed(args, seed_of):
    """`seed_of(skeleton)` for the stored galaxy, to record with a run that starts from its own process; `None` with no skeleton or
    no seed (the stage says the skeleton is missing when it opens, and a galaxy with no seed draws one then)."""
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        skeleton = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    return None if skeleton is None or skeleton.galaxy_seed is None else seed_of(skeleton)


def _refusals_as_errors(run):
    """`run()`, with a refusal a stage made when it opened (a `RuntimeError` in a worker, PERF.212) raised again as the `RuntimeError`
    it was, the way the verbs raised it before their stages took their params when they opened."""
    try:
        return run()
    except units.RunStopped as exc:
        marker = "RuntimeError: "
        if marker not in str(exc):
            raise
        raise RuntimeError(str(exc).split(marker, 1)[1]) from exc


def planned_bright_stars(args, passes=("mass", "luminosity"), refuse=None, tag=""):
    """The bright-star passes of `scatter_bright_stars` as a `followed.Planned`; its result is `book.summary()`."""
    followers, book = bright_star_followers(args, passes, refuse, tag)
    return followed.Planned(followers + [BlockLightsFollower(args)], "Bright stars", book.summary)


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
    shells (`static_shell_chunks`); `names` says which ones, `sizes` how many.
    Called from the run's own thread and the one following the units' progress
    (`wrapper.run`), so every change holds `lock`.
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


CHANNEL_INTERVAL_SECONDS = 0.25
"""float: How often a worker drawing a task reports its stars (PERF.4)."""

MAX_SHELLS_PER_CHUNK = 64
"""int: The most consecutive shells one task draws, so the progress moves steadily at the thin edge of the galaxy,
where thousands of shells hold almost nothing."""


def static_shell_chunks(rows):
    """
    The units of a scatter on the new queue (PERF.97): the outline cut into runs of consecutive shells of about
    equal expected work, from the rows alone (`units.static_chunks`), so the same galaxy is the same units at any
    worker count. The heaviest unit first.

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


def _unit_note(position, count, expected, what="stars"):
    """` (unit 3 of 512, about 1,200 stars expected)`: which unit of its stage a log line is about and how big the cut
    made it, so the sizes of the units can be read off the log."""
    return f" (unit {position:,} of {count:,}, about {round(expected):,} {what} expected)"


def _log_shell_range(label, shells, types, unit=""):
    """`_log_shell` for a task of several shells: the stars it added, by kind, over the shells it covers. `unit` is
    `_unit_note`'s, put after the shells."""
    total = sum(types.values())
    if total:
        log.normal(f"{label}, shell{'s' if len(shells) > 1 else ''} {_shell_names(shells)}{unit}: added {total:,} stars: "
                   f"{describe_types(types)}.")


COMPONENT_TASK_SHELL = -1
"""int: The shell the components' tasks are listed under in a scatter's tasks (they have none of their own)."""


COMPONENT_REPORT_ROWS = 200
"""int: How many of a component task's stars are drawn between the reports of its progress (`report` itself lets a few
through a second); a star takes a few milliseconds."""


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
    by shell (`expected`; the components' tasks are listed under `COMPONENT_TASK_SHELL`, as one total), the key of the
    far shells' task and the components' tasks (`component_parts`: task -> the `[kind, name, slice, slices]` lists of
    `component_pass.parts`). What a process following the run needs: it travels as JSON (`to_json`).
    """

    tasks: dict
    weights: dict
    sizes: dict
    names: dict
    expected: dict
    far_key: int = None
    component_parts: dict = dataclasses.field(default_factory=dict)
    first_far: int = 0
    reach_shell: int = 0
    reach_on: bool = False
    e_value: float = 0.0
    band_share: float = 0.0
    component_count: int = 0

    @property
    def special_keys(self):
        """The keys of the tasks that are not a run of shells: the far shells' and the components'."""
        return {self.far_key, *self.component_parts}

    def to_json(self):
        order = range(len(self.tasks))
        return {"tasks": [[list(row) for row in self.tasks[key]] for key in order],
                "weights": [self.weights[key] for key in order], "sizes": [self.sizes[key] for key in order],
                "names": [self.names[key] for key in order],
                "expected": [[shell, stars] for shell, stars in self.expected.items()],
                "far_key": self.far_key, "component_parts": [[key, parts] for key, parts in self.component_parts.items()],
                "first_far": self.first_far, "reach_shell": self.reach_shell, "reach_on": self.reach_on, "e_value": self.e_value,
                "band_share": self.band_share, "component_count": self.component_count}

    @classmethod
    def from_json(cls, data):
        order = range(len(data["tasks"]))
        return cls(
            tasks={key: [tuple(row) for row in data["tasks"][key]] for key in order},
            weights={key: data["weights"][key] for key in order}, sizes={key: data["sizes"][key] for key in order},
            names={key: data["names"][key] for key in order},
            expected={shell: stars for shell, stars in data["expected"]}, far_key=data["far_key"],
            component_parts={key: parts for key, parts in data["component_parts"]}, first_far=data["first_far"],
            reach_shell=data["reach_shell"],
            reach_on=data["reach_on"], e_value=data["e_value"], band_share=data["band_share"],
            component_count=data["component_count"])

    def expected_of(self, key):
        """The stars task `key` is expected to place: its shells' expected stars (the outline's and the part of each
        shell beyond it), the far shells' together, or the components' share. The task's `weights` also count the
        bands it walks, so they are not a count of stars (GEN.398)."""
        if key == self.far_key:
            return self.expected[self.first_far]
        if key in self.component_parts:
            return self.tasks[key][0][2]
        return sum(self.expected[shell] for shell, _high, _weight in self.tasks[key])

    def log_line(self, label):
        return (f"{label}: about {round(sum(self.expected.values())):,} stars expected in "
                f"{sum(self.sizes.values()):,} shells, {self.reach_shell:,} shells being the grid's reach"
                + (f", and {self.component_count:,} components." if self.component_count else "."))


def _band_fractions(tiered, min_luminosity_sol, max_luminosity_sol, mass_range, galaxy_shape):
    """A scatter pass's expected share of the stars by population (`fractions`) and its sum for the baseline
    (`band_share`), for a tiered pass or a single band, with the ages spread as the draw of `galaxy_shape` spreads
    them (GEN.546)."""
    age_laws = brightStars.age_laws_of(galaxy_shape)
    if tiered is None:
        fractions = brightStars.band_fractions(min_luminosity_sol, max_luminosity_sol, mass_range, age_laws)
        return fractions, sum(fractions.values())
    fractions = brightStars.tier_fractions_of(tiered, age_laws)
    shares = [sum(item.values()) for item in fractions]
    return fractions, shares[scatter_limits.BASELINE] or max(shares)


SHELL_WEIGHT_PARALLEL_MIN = 64
"""int: Fewer shells than this are weighed in the process itself (a pool costs more than it saves)."""

SHELL_WEIGHT_PROCESSES = 8
"""int: Most processes that weigh a pass's shells (`_shell_weights`)."""

_weighing = None
"""tuple or None: What a forked weighing process reads (`_shell_weights` sets it before the fork)."""


def _weigh(shells, start, stop):
    skeleton, e_value, fractions, estimates = _weighing
    return [brightStars.shell_weight(skeleton.shape, shell_index, high_band, skeleton.edge_pc, e_value, fractions,
                                     beyond=None if estimates is None else estimates[shell_index])
            for shell_index, high_band in shells[start:stop]]


def _shell_weights(skeleton, shells, e_value, fractions, estimates):
    """
    `brightStars.shell_weight` of each `(shell, high_band)` of `shells`, in order: `(weight, expected stars)`. Each shell
    samples the density at a few dozen places, 10 ms a shell, and the passes weigh 2,600 shells before any unit can start
    (26 s at the default size, per pass, with the workers idle), so the shells go to forked processes in runs. Each shell
    is worked out by the same call on the same numbers wherever it runs, so the weights are those of a loop (PERF.159).
    """
    global _weighing
    count = min(os.cpu_count() or 1, SHELL_WEIGHT_PROCESSES)
    if count < 2 or len(shells) < SHELL_WEIGHT_PARALLEL_MIN or "fork" not in multiprocessing.get_all_start_methods():
        _weighing = (skeleton, e_value, fractions, estimates)
        try:
            return _weigh(shells, 0, len(shells))
        finally:
            _weighing = None
    runs = max(1, len(shells) // (count * 8))
    cuts = [(start, min(start + runs, len(shells))) for start in range(0, len(shells), runs)]
    _weighing = (skeleton, e_value, fractions, estimates)
    try:
        with concurrent.futures.ProcessPoolExecutor(count, mp_context=multiprocessing.get_context("fork")) as pool:
            parts = list(pool.map(_weigh, [shells] * len(cuts), [cut[0] for cut in cuts], [cut[1] for cut in cuts]))
    finally:
        _weighing = None
    return [weighed for part in parts for weighed in part]


def _plan_shape(skeleton, extents, fractions, band_share, components, tiered):
    """
    Cuts a shell scatter into tasks (`_ScatterShape`): the outline's shells and the shells past it (GEN.322: the
    pass covers the grid's whole reach; the shells past the outline are drawn one by one while they hold stars to
    speak of, and the faint ones further out together, as one task), then the components' tasks (GEN.324): the small
    components bundled, a large one in slices (`component_pass.parts`). The cut is the static one
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
    shells = list(extents) + beyond_shells
    for (shell_index, high_band), (weight, expected[shell_index]) in zip(
            shells, _shell_weights(skeleton, shells, e_value, fractions, estimates if reach_on else None)):
        rows.append((shell_index, high_band, weight))
    tasks = dict(enumerate(static_shell_chunks(rows)))
    weights = {key: sum(row[2] for row in chunk) for key, chunk in tasks.items()}
    sizes = {key: len(chunk) for key, chunk in tasks.items()}
    names = {key: _shell_names(chunk) for key, chunk in tasks.items()}
    far_key = None
    component_parts = {}
    if first_far <= reach_shell:
        far_key = len(tasks)
        far_expected = brightStars.far_expected_stars(
            skeleton.shape, first_far, reach_shell, skeleton.edge_pc, e_value, fractions)
        expected[first_far] = far_expected
        tasks[far_key] = [(first_far, reach_shell, far_expected)]
        weights[far_key] = far_expected + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS
        sizes[far_key] = reach_shell - first_far + 1
        names[far_key] = f"{first_far}-{reach_shell}"
    # GEN.324: the globular clusters and satellite galaxies place their own stars, in tasks of their own.
    if components:
        cut = componentPass.parts(componentPass.star_weights(components, fractions, tiered),
                                  componentPass.STAR_SLICE_STARS)
        # The cut is by the candidates (the most a slice holds, and what keeps the streams where they were); the stars
        # a task is expected to place are what its components keep at the places they lie (GEN.678).
        kept = {(component.kind, component.name): stars for component, stars in
                componentPass.kept_weights(components, fractions, tiered, skeleton.edge_pc)}
        cut = [(sum(kept[(kind, name)] / slices for kind, name, _index, slices in task_parts), task_parts)
               for _candidates, task_parts in cut]
        expected[COMPONENT_TASK_SHELL] = sum(task_expected for task_expected, _parts in cut)
        for number, (task_expected, task_parts) in enumerate(cut, 1):
            key = len(tasks)
            component_parts[key] = task_parts
            tasks[key] = [(COMPONENT_TASK_SHELL, COMPONENT_TASK_SHELL, task_expected)]
            weights[key] = task_expected + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS
            sizes[key] = 0
            names[key] = f"components {number} of {len(cut)}"
    return _ScatterShape(tasks, weights, sizes, names, expected, far_key, component_parts, first_far, reach_shell,
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
        special = key in shape.special_keys
        return {
            "shells": [] if special else [(shell, high) for shell, high, _weight in chunk],
            "far": (shape.first_far, shape.reach_shell) if key == shape.far_key else None, "reach": shape.reach_on,
            "components": self.components if key in shape.component_parts else None,
            "component_parts": shape.component_parts.get(key),
            "expected_by_shell": ({shape.first_far: shape.expected[shape.first_far]} if key == shape.far_key else
                                  {COMPONENT_TASK_SHELL: chunk[0][2]} if key in shape.component_parts else
                                  {shell: shape.expected[shell] for shell, _high, _weight in chunk}),
            "edge_pc": self.skeleton.edge_pc, "expected": shape.e_value,
            "min_luminosity_sol": self.min_luminosity_sol, "max_luminosity_sol": self.max_luminosity_sol,
            "seed": self.seed, "skip": (set(self.filled) if key in shape.component_parts else
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
            "reach": payload["reach"], "component_parts": payload["component_parts"],
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

    def __init__(self, args, label, stats_kind, workers, progress=None):
        self.args, self.label, self.stats_kind, self.workers = args, label, stats_kind, workers
        self.bar = self.tracker = None
        self.progress = progress        # a run's display, which is not this scatter's to open or close
        self._display = None if progress is not None else run_common._generation_progress()

    def __enter__(self):
        if self._display is not None:
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
        if self._display is None:
            if self.bar is not None:
                self.bar.close(success=exc_type is None)
            return False
        try:
            if self.bar is not None:
                self.bar.close(success=exc_type is None)
        finally:
            log.reset_console()
            suppressed = self._display.__exit__(exc_type, exc, traceback)
        return suppressed


class _ScatterTally:
    """What a shell scatter adds up as its tasks finish: the stars written per population, the shells that drew
    some and the stars with companions. The speed of each task is the stage's to record (`finish_scatter`)."""

    def __init__(self, label, shape, tracker):
        self.label, self.shape, self.tracker = label, shape, tracker
        self.counts = {population: 0 for population in brightStars.POPULATIONS}
        self.drew = set()
        self.systems = 0

    def on_done_for(self, key):
        shape, chunk = self.shape, self.shape.tasks[key]

        def task_done(result):
            task_counts = result["counts"]
            types = collections.Counter({(name, plural): count for name, plural, count in result["types"]})
            note = _unit_note(key + 1, len(shape.tasks), shape.expected_of(key))
            if key in shape.component_parts:
                if sum(types.values()):
                    log.normal(f"{self.label}, {shape.names[key]}{note}: added {sum(types.values()):,} stars: "
                               f"{describe_types(types)}.")
            else:
                _log_shell_range(self.label, chunk, types, note)
            for population, count in task_counts.items():
                self.counts[population] += count
            self.systems += result.get("systems", 0)
            self.tracker.task_done(key)
            for shell_index, shell_stars in result["shells"]:
                if shell_stars:
                    self.drew.add(shell_index)
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


STATS_KIND_OF_STAGE = {"mass": "scatter-mass", "luminosity": "scatter-bright", "band": "scatter-bright"}
"""dict: The speed stats kind of the tasks of each shell scatter (PERF.10)."""


def _task_density(shape, key):
    """The density a shell scatter task's speed is recorded under: the stars it was expected to draw over what a
    standard sector at density 1 holds in its cells, scaled to the pass's share of the stars (PERF.53)."""
    chunk = shape.tasks[key]
    cells = _cells_in([(shell, high) for shell, high, _weight in chunk])
    wanted = sum(shape.expected[shell] for shell, _high, _weight in chunk)
    return wanted / (shape.e_value * shape.band_share * cells) if shape.band_share and cells else 0.0


def finish_scatter(params, plan, results, seconds):
    """
    The `finish` of the shell scatter commands (PERF.106), called by the unit that completes the stage: records each
    unit's speed in the generation stats (PERF.10, PERF.53) and gives the metrics the stage's row is stored with.

    Args:
        params (dict): The step's parameters (`database`, and `workers` for the stats' bucket).
        plan (dict): The stage's plan (`_ScatterShape` as JSON, and `stats_kind`).
        results (list[dict]): Every unit's result (`scatter_shells_unit`), in unit order.
        seconds (float): How long the stage took.

    Returns:
        dict: `shells` (visited), `shells_modified` (that drew a star) and `objects` (stars written).
    """
    shape = _ScatterShape.from_json(plan)
    stats = _speed_stats(_mysql_config_of(params))
    drew, objects = set(), 0
    for key, result in enumerate(results):
        stars = sum(result["counts"].values())
        objects += stars
        drew.update(shell for shell, shell_stars in result["shells"] if shell_stars)
        if key not in shape.special_keys and stars:   # PERF.53: no star, no time per star
            stats.record(plan["stats_kind"], _task_density(shape, key), result["seconds"], systems=stars, stars=stars,
                         workers=params.get("workers", 1))
    stats.flush()
    return {"shells": sum(shape.sizes.values()), "shells_modified": len(drew), "objects": objects}


def _speed_stats(mysql_config):
    """The stored speeds, read fresh (a worker lives across runs, so a copy kept from an earlier one would write
    older buckets back), or nothing to record to when statistics are off."""
    if not stages.recording_enabled():
        return generationStats.GenerationStats()
    return generationStats.GenerationStats(store.control_mysql_config(mysql_config))


def _mysql_config_of(params):
    """The `MySQLConfig` a command's params name (`wrapper.database_args`; the password is this process's settings')."""
    named = params["database"]
    return store.MySQLConfig(host=named["host"], port=named["port"], user=named["user"], database=named["database"])


def _expansion(mysql_config, skeleton, extents, skip, min_luminosity_sol, max_luminosity_sol, mass_range, tiered, seed,
               stage, extra=None):
    """
    A shell scatter's stage cut into its units (`_ScatterShape`, `static_shell_chunks`): the same at any worker count,
    from the galaxy and the pass's limits alone. `skip` is the sectors the pass leaves out; `stage` names the pass
    (`STATS_KIND_OF_STAGE`); `extra` is added to the stage's JSON summary.

    Returns:
        commands.Expansion: A unit per task, its cost the task's expected stars, and the shape as JSON.
    """
    fractions, band_share = _band_fractions(tiered, min_luminosity_sol, max_luminosity_sol, mass_range, skeleton.shape)
    components = _galaxy_components(mysql_config)
    shape = _plan_shape(skeleton, extents, fractions, band_share, components, tiered)
    drawing = _ScatterPass(skeleton, skip, min_luminosity_sol, max_luminosity_sol, seed, mass_range, tiered, components)
    return commands.Expansion([(drawing.unit_args(shape, key, mysql_config), shape.weights[key]) for key in shape.tasks],
                              {**shape.to_json(), "stats_kind": STATS_KIND_OF_STAGE[stage], **(extra or {})})


def _stored_galaxy(conn):
    """`(skeleton, extents, filled)` of the stored galaxy, for a shell scatter.

    Raises:
        RuntimeError: The skeleton was never built.
    """
    skeleton = store.get_galaxy_shape(conn)
    if skeleton is None:
        raise RuntimeError("The galaxy's skeleton has never been built -- run 'planetgen plan' first.")
    return skeleton, _outline(conn), store.filled_sector_addresses(conn)


def _scatter_facts(params, limits):
    """What a bright-star pass's plan says about the scatter it belongs to, for the log the followers write (`_BrightStarBook.begin`):
    the mass limit, the seed and the limits line, which were worked out when the stage opened (PERF.212)."""
    return {"mass_limit_sol": params["mass_limit_sol"], "seed": params["seed"], "limits_text": _limits_text(limits)}


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
    mysql_config = _mysql_config_of(params)
    conn = store.get_connection(mysql_config)
    try:
        skeleton, extents, filled = _stored_galaxy(conn)
    finally:
        conn.close()
    limits = scatter_limits.ScatterLimits.of(skeleton, params["min_luminosity_sol"], params["mass_limit_sol"])
    floor = starPopulation.MASS_PASS_MIN_LUMINOSITY_SOL
    return _expansion(mysql_config, skeleton, extents, filled, floor, None, None,
                      scatter_limits.mass_pass(limits, floor), params["seed"], "mass",
                      {"filled": len(filled), **_scatter_facts(params, limits)})


def expand_scatter_luminosity(params):
    """
    The `scatter-luminosity` command's stage cut into its units (PERF.106): the luminosity pass (GEN.185 pass 4):
    every star born lighter than the mass limit and at least the luminosity floor bright, in the sectors that are
    neither filled nor marked by the mass pass (a sector holding a mass-pass star that bright already has one). The
    marks are read here, when the stage starts, between the passes, while the table holds only the mass pass's stars.

    Args:
        params (dict): `database` (`wrapper.database_args`), `min_luminosity_sol` (the baseline luminosity floor),
            `mass_limit_sol` (the stellar mass limit) and `seed` (the pass's).

    Returns:
        commands.Expansion: The units, and the `_ScatterShape` as JSON with `marked`, the sectors left out as marked.
    """
    mysql_config = _mysql_config_of(params)
    conn = store.get_connection(mysql_config)
    try:
        skeleton, extents, filled = _stored_galaxy(conn)
        limits = scatter_limits.ScatterLimits.of(skeleton, params["min_luminosity_sol"], params["mass_limit_sol"])
        marked = store.bright_star_marked_addresses(conn, limits)
    finally:
        conn.close()
    return _expansion(mysql_config, skeleton, extents, set(filled) | marked, params["min_luminosity_sol"], None, None,
                      scatter_limits.luminosity_pass(limits), params["seed"], "luminosity",
                      {"marked": len(marked), "filled": len(filled), **_scatter_facts(params, limits)})


def expand_scatter_band(params):
    """
    The `scatter-band` command's stage cut into its units (PERF.106): the band of bright stars from the new level up
    to the stored one (`planetgen scatter --band --down-to`), each tier from its new floor up to its old one, in the sectors that
    are neither filled nor backfilled to their own level.

    Args:
        params (dict): `database` (`wrapper.database_args`), `target_luminosity_sol` (the new level),
            `current_luminosity_sol` (the stored one), `mass_limit_sol` (the stored mass limit, or `None`) and `seed`.

    Returns:
        commands.Expansion: The units, and the `_ScatterShape` as JSON with `own`, the backfilled sectors left out.
    """
    mysql_config = _mysql_config_of(params)
    mass_limit = params["mass_limit_sol"]
    conn = store.get_connection(mysql_config)
    try:
        skeleton, extents, filled = _stored_galaxy(conn)
        own = store.sector_bright_level_keys(conn)
    finally:
        conn.close()
    older = scatter_limits.ScatterLimits.of(skeleton, params["current_luminosity_sol"], mass_limit)
    newer = scatter_limits.ScatterLimits.of(skeleton, params["target_luminosity_sol"], mass_limit)
    return _expansion(mysql_config, skeleton, extents, set(filled) | set(own), params["target_luminosity_sol"],
                      params["current_luminosity_sol"], None if mass_limit is None else (None, mass_limit),
                      scatter_limits.band_pass(newer, older), params["seed"], "band", {"own": len(own)})


class ScatterFollower(followed.Follower):
    """
    A shell scatter pass (PERF.97, PERF.106, PERF.118) as a stage of a run: the wrapper selects the stage's command
    (`scatter-mass`, `scatter-luminosity` or `scatter-band`) from the run's steps, cuts the pass into units (the
    command's `expand`) and queues them, and this process follows the units' results and progress with the same bar
    as the other scatters until all are done. `--workers` (or `PLANETGEN_WORKERS`) is how many workers the run starts,
    one included.

    Args:
        key (str): The stage (`mass`, `luminosity` or `band`).
        params (dict): What the command's expansion and hooks read besides the database (`expand_scatter_mass` and the
            others, `prepare_scatter_pass`, `conclude_scatter_pass`).
        seed (int or None): The pass's seed; `None` when it comes with the params' deferred settings (PERF.212).
        book (_BrightStarBook, optional): The scatter the pass is one of, which says its limits and adds it up.
        last (bool): Whether the pass is the scatter's last, which sums it up.
        own_row (bool): Whether the pass is the whole stage, so the unit that completes it stores the stage's row; false
            when the starter still does more inside the stage (the band tops up backfilled sectors), so the starter
            stores the row for all of it.

    Attributes:
        counts (dict): Stars written per population, once the pass is done.
    """

    def __init__(self, args, key, label, params, seed, book=None, last=False, own_row=True):
        super().__init__(args)
        self.key, self.label, self.own_row = key, label, own_row
        self._params, self._seed, self.book, self.last = params, seed, book, last
        self.counts = None
        self.shape = self.tally = self.bars = None

    def params(self):
        return {**self._params, **({} if self._seed is None else {"seed": self._seed})}

    def on_plan(self, plan):
        detail = plan["detail"]
        if self.book is not None:
            self.book.begin(detail)
        self.shape = _ScatterShape.from_json(detail)
        if detail.get("marked") is not None:
            log.normal(f"{detail['marked']:,} sectors hold a star born at their band's mass limit or more "
                       "and at least its luminosity floor bright; the luminosity pass skips them.")
        log.normal(self.shape.log_line(self.label))
        self.bars = self.hold(_ShellBars(self.args, self.label, STATS_KIND_OF_STAGE[self.key], self.workers,
                                         progress=self.progress))
        self.bars.open(self.shape)
        self.tally = _ScatterTally(self.label, self.shape, self.bars.tracker)

    def on_result(self, unit, result):
        self.tally.on_done_for(unit)(result)

    def on_progress(self, unit, done, estimate):
        self.bars.tracker.task_progress(unit, done, estimate)

    def on_done(self, conclusion):
        self.counts = self.tally.finish()
        if self.book is not None and self.last:
            self.book.log()


def _scatter_on_queue(args, key, label, params, seed, own_row=True):
    """A shell scatter pass as a run of its own (`ScatterFollower`). Returns the stars written per population."""
    return followed.run_planned(args, planned_scatter(args, key, label, params, seed, own_row))


def planned_scatter(args, key, label, params, seed, own_row=True):
    """A shell scatter pass as a `followed.Planned`; its result is the stars written per population."""
    follower = ScatterFollower(args, key, label, params, seed, own_row=own_row)
    return followed.Planned([follower], label, lambda: follower.counts, seed=seed)


def add_bright_star_band(args):
    """
    Lowers the galaxy's star-fill level to `--down-to`:
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
    target = float(args.down_to)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError("The galaxy's skeleton has never been built -- run 'planetgen plan' first.")
        settings = store.bright_star_scatter_settings(conn)
        mass_limit = store.bright_star_mass_limit(conn)
        if settings is None:
            log.normal("No bright stars are scattered yet, so there is no level to go below. Run "
                       f"'planetgen scatter --luminosity --bright-star-min-luminosity {target:g}' instead.")
            return None
        current, first_seed = float(settings[0]), settings[1]
        # GEN.333: each tier's floor moves with the baseline unless it was set, so each band goes from its tier's new
        # floor up to its old one, and a tier whose floor does not fall adds nothing.
        older = scatter_limits.ScatterLimits.of(skeleton, current, mass_limit)
        newer = scatter_limits.ScatterLimits.of(skeleton, target, mass_limit)
        if target >= current:
            log.normal(f"Nothing to add: every star of {current:g} L_sun or more is already placed, "
                       f"and {target:g} is not below that.")
            return None
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
        counts = _scatter_on_queue(args, "band", f"Bright stars {target:g}-{current:g} L_sun", {
            "target_luminosity_sol": target, "current_luminosity_sol": current,
            "mass_limit_sol": None if mass_limit is None else float(mass_limit)}, seed, own_row=False)
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
        store.record_scatter_epoch(conn, "band", getattr(args, "generation_run_id", None))
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


def _unit_skeleton(unit_args):
    """`(mysql_config, skeleton)` of a scatter unit: the database it names, with the credentials of this worker's
    settings, and the galaxy's shape read from it."""
    named = unit_args["database"]
    mysql_config = store.MySQLConfig(host=named["host"], port=named["port"], user=named["user"],
                                     database=named["database"])
    conn = store.get_connection(mysql_config)
    try:
        return mysql_config, store.get_galaxy_shape(conn)
    finally:
        conn.close()


def scatter_shells_warm(unit_args):
    """
    The `warm` function of the shell scatters (`units.warm_unit`, PERF.146): reads the shape a unit reads, which builds the
    component field, and works out the pass's sampling tables, the share of each population's stars in the band
    (`brightStars.pass_of`, the cached `star_population._bright_table`), which `scatter_shell` looks up first thing
    and which take seconds to build. The worker does it once and every horse it forks starts with both.
    """
    mysql_config, skeleton = _unit_skeleton(unit_args)
    try:
        mass_range = unit_args.get("mass_range")
        mass_range = None if mass_range is None else tuple(mass_range)
        tiered = scatter_limits.tiered_from_json(unit_args.get("tiered"), skeleton)
        age_laws = brightStars.age_laws_of(skeleton.shape)
        _, pass_bands = brightStars.pass_of(tiered, unit_args["min_luminosity_sol"], unit_args.get("max_luminosity_sol"),
                                            mass_range, age_laws)
        if unit_args.get("component_parts"):       # the clusters' and satellites' stars: an even spread of ages
            brightStars.pass_of(tiered, unit_args["min_luminosity_sol"], unit_args.get("max_luminosity_sol"),
                                mass_range)
        for band in pass_bands:                    # a star draws from the table of its own [Fe/H] bin (GEN.651)
            for population in starPopulation.POPULATIONS:
                starPopulation.warm_bins(band.min_luminosity_sol, population, band.mass_range, age_laws)
    finally:
        store.close_pool(mysql_config)


def scatter_shells_unit(unit_args, report):
    """
    A `scatter-mass` unit of the new queue (PERF.97, `queue/units.py`; PERF.118, `queue/commands.py`): the shells `_scatter_shells_task` draws and
    writes for a run of shells, from plain JSON arguments. The galaxy's shape comes from the database, which this
    unit opens with the credentials of this worker's settings and the database, host, port and user the run named.

    Returns:
        dict: `_scatter_shells_task`'s, plus `seconds` the unit took.
    """
    started = time.monotonic()
    mysql_config, skeleton = _unit_skeleton(unit_args)
    mass_range = unit_args.get("mass_range")
    far = unit_args.get("far")
    result = _scatter_shells_task({
        "mysql_config": mysql_config, "shape": skeleton.shape, "key": unit_args["key"],
        "shells": [tuple(pair) for pair in unit_args["shells"]],
        "far": None if far is None else tuple(far), "reach": unit_args.get("reach", False),
        "components": _galaxy_components(mysql_config) if unit_args.get("component_parts") else None,
        "component_parts": unit_args.get("component_parts"),
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
            these, slices = componentPass.chosen(payload["components"], payload["component_parts"])
            runs.append((COMPONENT_TASK_SHELL, componentPass.scatter_stars(
                these, payload["edge_pc"], payload["seed"], fractions, pass_bands,
                tiered=payload.get("tiered"), skip_addresses=payload["skip"], slices=slices)))
        for shell_index, rows in runs:
            remaining[0] -= expected_by_shell[shell_index]
            drawn = 0
            landed = collections.Counter()
            for row in rows:
                counts[row[6]] += 1
                types[star_label(row[7], row[8])] += 1
                drawn += 1
                if shell_index == COMPONENT_TASK_SHELL and drawn % COMPONENT_REPORT_ROWS == 0:
                    report(drawn, expected_by_shell[shell_index])     # a component's stars come one at a time
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


def _log_phenomena_shells(shells, counts, unit=""):
    """One line saying how many phenomena a run of shells was given, by kind, as `_log_shell_range` does for stars;
    nothing for a run that drew none."""
    total = sum(counts.values())
    if total:
        log.normal(f"Phenomena, shell{'s' if len(shells) > 1 else ''} {_shell_names(shells)}{unit}: placed {total:,} phenomena: "
                   + ", ".join(f"{count:,} {kind}" for kind, count in sorted(counts.items())) + ".")


@dataclasses.dataclass
class _PhenomenaShape:
    """
    How the phenomenon scatter is cut (`_plan_phenomena`): the tasks (`tasks`, unit -> its `(shell_index, high_band,
    weight)` rows), each one's expected work (`weights`), the `far` task that draws the faint shells far out together
    (`(key, first shell, last shell)`, GEN.323), the components' tasks (GEN.324; `component_parts`: task -> the
    `[kind, name, slice, slices]` lists of `component_pass.parts`), whether each shell also draws its part beyond the
    outline (`reach`) and the totals the log and the bar count. It travels as JSON.
    """

    tasks: dict
    weights: dict
    far: tuple = None
    component_parts: dict = dataclasses.field(default_factory=dict)
    reach: bool = False
    shell_total: int = 0
    expected_rows: float = 0.0
    expected: dict = dataclasses.field(default_factory=dict)

    def expected_of(self, key):
        """The phenomena a task is expected to place: its rows, not the bar's weight (which adds a cost per band)."""
        return self.expected[key]

    @property
    def far_key(self):
        return None if self.far is None else self.far[0]

    def shells_of(self, key):
        """A task's shells as the `(shell, high, weight)` rows `_shell_names` reads."""
        if key in self.component_parts:
            return []
        if key == self.far_key:
            return [(shell, -1, 0.0) for shell in range(self.far[1], self.far[2] + 1)]
        return self.tasks[key]

    def to_json(self):
        order = range(len(self.tasks))
        return {"tasks": [[list(row) for row in self.tasks[key]] for key in order],
                "weights": [self.weights[key] for key in order], "far": None if self.far is None else list(self.far),
                "component_parts": [[key, parts] for key, parts in self.component_parts.items()], "reach": self.reach, "shell_total": self.shell_total,
                "expected_rows": self.expected_rows, "expected": [self.expected[key] for key in order]}

    @classmethod
    def from_json(cls, data):
        order = range(len(data["tasks"]))
        return cls(tasks={key: [tuple(row) for row in data["tasks"][key]] for key in order},
                   weights={key: data["weights"][key] for key in order},
                   far=None if data["far"] is None else tuple(data["far"]),
                   component_parts={key: parts for key, parts in data["component_parts"]},
                   reach=data["reach"], shell_total=data["shell_total"], expected_rows=data["expected_rows"],
                   expected={key: data["expected"][key] for key in order})


def _plan_phenomena(skeleton, extents, galaxy_components, min_mass_solar):
    """
    Cuts the phenomenon scatter into tasks (`_PhenomenaShape`), the same at any worker count (`static_shell_chunks`).
    GEN.323: the pass covers the grid's whole reach, as the stars' does (GEN.322): the outline's shells also draw their
    bands beyond it, the shells past the outline one by one while they hold objects to speak of, and the faint ones
    further out together, as one task. GEN.324: the clusters' and satellites' phenomena are one more task.
    """
    e_value = skeleton.expected_system_count_at_density_1
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
    rows, rows_of_shell = [], {}
    expected_rows = 0.0     # the rows expected, not the bar's weight (that adds a cost per band, GEN.397)
    for shell, high in list(extents) + beyond_shells:
        weight = shell_rows = 0.0
        if high >= 0:
            shell_rows = phenomenon_scatter.shell_expected(skeleton.shape, shell, high, skeleton.edge_pc, e_value,
                                                           min_mass_solar)
            weight = shell_rows + brightStars.BAND_WEIGHT_STARS * (2 * high + 1)
        if reach_on:
            shell_rows += estimates[shell][0]
            weight += estimates[shell][0] + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS * estimates[shell][1]
        expected_rows += shell_rows
        rows_of_shell[shell] = shell_rows
        rows.append((shell, high, weight))
    tasks = dict(enumerate(static_shell_chunks(rows)))
    weights = {key: sum(row[2] for row in chunk) for key, chunk in tasks.items()}
    expected = {key: sum(rows_of_shell[row[0]] for row in chunk) for key, chunk in tasks.items()}
    far = None
    if reach_on and first_far <= reach_shell:
        far = (len(tasks), first_far, reach_shell)
        far_weight = phenomenon_reach.far_expected(
            skeleton.shape, first_far, reach_shell, skeleton.edge_pc, e_value, min_mass_solar)
        expected_rows += far_weight
        tasks[far[0]] = [(first_far, reach_shell, far_weight)]
        weights[far[0]] = far_weight + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS
        expected[far[0]] = far_weight
    component_parts = {}
    if galaxy_components:
        cut = componentPass.parts([(component, componentPass.expected_phenomena(component, min_mass_solar))
                                   for component in galaxy_components], componentPass.PHENOMENON_SLICE_ROWS)
        for task_weight, task_parts in cut:
            key = len(tasks)
            component_parts[key] = task_parts
            expected_rows += task_weight
            tasks[key] = [(COMPONENT_TASK_SHELL, COMPONENT_TASK_SHELL, task_weight)]
            weights[key] = task_weight + brightStars.OBJECT_FIRST_BAND_WEIGHT_STARS
            expected[key] = task_weight
    shell_total = len(extents) + len(beyond_shells) + (0 if far is None else reach_shell - first_far + 1)
    return _PhenomenaShape(tasks, weights, far, component_parts, reach_on, shell_total, expected_rows, expected)


def expand_scatter_phenomena(params):
    """
    The `scatter-phenomena` command's stage cut into its units (PERF.106): the shells of the phenomenon scatter
    (neutron stars, black holes, planetary nebulae, supernova remnants, hypervelocity stars, and the clusters' and
    satellites' own), a unit per run of shells, from `params` and the stored galaxy alone. The earlier scatter is
    cleared, and the special phenomena, the nebula field and the indexes are done, around the stage, not by it.

    Args:
        params (dict): `database` (`wrapper.database_args`), `min_mass_solar` (the lowest mass of a neutron star or
            black hole placed) and `seed` (the scatter's).

    Returns:
        commands.Expansion: The units, a unit's cost its task's expected work, and the `_PhenomenaShape` as JSON.
    """
    mysql_config = _mysql_config_of(params)
    conn = store.get_connection(mysql_config)
    try:
        skeleton, extents, filled = _stored_galaxy(conn)
    finally:
        conn.close()
    shape = _plan_phenomena(skeleton, extents, _galaxy_components(mysql_config), params["min_mass_solar"])
    filled_by_shell = collections.defaultdict(list)
    for address in filled:
        filled_by_shell[address[0]].append(list(address))
    units_ = []
    for key in range(len(shape.tasks)):
        special = key == shape.far_key or key in shape.component_parts
        shells = [] if special else [(shell, high) for shell, high, _weight in shape.tasks[key]]
        if key in shape.component_parts:
            skip = sorted(list(address) for address in filled)
        else:
            skip = sorted(address for shell, _high, _weight in shape.shells_of(key) for address in filled_by_shell[shell])
        units_.append(({"database": wrapper.database_args(mysql_config), "key": key, "shells": [list(pair) for pair in shells],
                        "far": None if key != shape.far_key else list(shape.far[1:]), "reach": shape.reach,
                        "component_parts": shape.component_parts.get(key), "seed": params["seed"],
                        "min_mass_solar": params["min_mass_solar"], "skip": skip}, shape.weights[key]))
    return commands.Expansion(units_, {**shape.to_json(), "filled": len(filled), "outline_shells": len(extents),
                                       "min_mass_solar": params["min_mass_solar"]})


def scatter_phenomena_unit(unit_args, report):
    """
    A `scatter-phenomena` unit (PERF.106): the shells `_phenomenon_shells_task` draws and writes for a run of shells,
    from plain JSON arguments; the galaxy's shape and components come from the database as `scatter_shells_unit`'s do.

    Returns:
        dict: `_phenomenon_shells_task`'s, plus `seconds` the unit took.
    """
    started = time.monotonic()
    mysql_config = store.MySQLConfig(**unit_args["database"])
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    far = unit_args.get("far")
    result = _phenomenon_shells_task({
        "mysql_config": mysql_config, "shape": skeleton.shape,
        "shells": [tuple(pair) for pair in unit_args["shells"]], "far": None if far is None else tuple(far),
        "reach": unit_args["reach"], "components": _galaxy_components(mysql_config) if unit_args["component_parts"] else None,
        "component_parts": unit_args["component_parts"],
        "edge_pc": skeleton.edge_pc, "expected": skeleton.expected_system_count_at_density_1, "seed": unit_args["seed"],
        "min_mass_solar": unit_args["min_mass_solar"], "skip": {tuple(address) for address in unit_args["skip"]}})
    result["seconds"] = time.monotonic() - started
    return result


def finish_scatter_phenomena(params, plan, results, seconds):
    """
    The `finish` of the `scatter-phenomena` command (PERF.106), called by the unit that completes the stage: records
    each unit's speed in the generation stats, as the bar's next run starts from it (PERF.53: not for a unit that
    placed nothing), the drawing and the writing of its rows each under a kind of their own (PERF.54).

    Returns:
        dict: `shells` (visited), `shells_modified` (units that placed something) and `objects` (rows written).
    """
    stats = _speed_stats(_mysql_config_of(params))
    objects = landed = 0
    for key, result in enumerate(results):
        placed = sum(result["counts"].values())
        objects += placed
        if placed:
            landed += 1
            weight = plan["weights"][key]
            stats.record("phenomena", 0.0, result["seconds"], systems=weight, stars=weight,
                         workers=params.get("workers", 1))
            stats.record("phenomena-write", 0.0, result["write_seconds"], systems=weight, stars=weight,
                         workers=params.get("workers", 1))
    stats.flush()
    return {"shells": plan["shell_total"], "shells_modified": landed, "objects": objects}


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


SKELETON_NEVER_BUILT = "The galaxy's skeleton has never been built -- run 'planetgen plan' first."
"""str: What a stage that needs the plan says when the galaxy has none (said when the stage opens, PERF.203)."""

NEBULA_DRAW_KEY = "nebula-draw"
"""str: The key of the step that draws the nebula field: a step of a run of its own, not a stage the plan lists."""


def nebula_slab_cut(first, last):
    """
    The units of the nebula field for the slabs `first` to `last` (inclusive): `(first slab, slabs)` for each, a unit a
    slab for a galaxy of up to `units.MAX_UNITS` slabs and `units.coarsen`'s whole multiple of that beyond it (PERF.192),
    so a large galaxy makes bigger units, not a stage `enqueue_units` refuses. The same slabs always give the same cut.
    """
    count = last - first + 1
    per = units.coarsen(count, 1)
    return [(start, min(per, last + 1 - start)) for start in range(first, last + 1, per)]


def expand_draw_nebula_field(params):
    """
    The `draw-nebula-field` command's stage (PERF.106, DB.24): a unit per first-index slab of the field's cells (a run of
    consecutive slabs per unit in a galaxy with more than `units.MAX_UNITS` of them, PERF.192), from the stored galaxy
    alone. `detail` holds the cell edge, so a follower can name a slab's extent. A galaxy with no outline has no unit.
    """
    conn = store.get_connection(_mysql_config_of(params))
    try:
        skeleton = store.get_galaxy_shape(conn)
        extents = _outline(conn)
    finally:
        conn.close()
    if skeleton is None:
        raise RuntimeError(SKELETON_NEVER_BUILT)
    if not extents:
        return commands.Expansion([], {"cell_pc": program_constants.NEBULA_FIELD_CELL_PC, "slabs": 0})
    outer_shell = max(shell for shell, _high in extents)
    radius_pc, half_height_pc = nebulaField.field_extent_pc(skeleton.edge_pc, extents)
    edge = program_constants.NEBULA_FIELD_CELL_PC
    first, last = math.floor(-radius_pc / edge), math.floor(radius_pc / edge)
    return commands.Expansion(
        [({"database": params["database"], "i": start, "count": count, "radius_pc": radius_pc,
           "half_height_pc": half_height_pc, "edge_pc": skeleton.edge_pc, "outer_shell": outer_shell}, float(count))
         for start, count in nebula_slab_cut(first, last)],
        {"cell_pc": edge, "slabs": last - first + 1})


def draw_nebula_field_unit(params, report):
    """
    The `draw-nebula-field` command's unit: the kept classes M and N clouds of its slabs (`count` consecutive ones from
    `i`, one in a galaxy that needs no coarser cut), from the galaxy's seed and shape, written to `nebula_field`
    (`_nebula_field_task`). The table is cleared before the run.

    Returns:
        dict: `written` (rows written), `slab` (the first slab's index), `count` (slabs drawn) and `drawn`, one
            `_nebula_field_task` result with its `slab` for each slab.
    """
    mysql_config = _mysql_config_of(params)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        galaxy_seed = store.get_galaxy_seed(conn)
    finally:
        conn.close()
    drawn = []
    for i in range(params["i"], params["i"] + params["count"]):
        drawn.append({**_nebula_field_task({
            "mysql_config": mysql_config, "galaxy_seed": galaxy_seed, "shape": skeleton.shape, "i": i,
            "radius_pc": params["radius_pc"], "half_height_pc": params["half_height_pc"], "edge_pc": params["edge_pc"],
            "outer_shell": params["outer_shell"]}), "slab": i})
    return {"written": sum(result["written"] for result in drawn), "slab": params["i"], "count": params["count"],
            "drawn": drawn}


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


class NebulaFieldFollower(followed.Follower):
    """
    The nebula field as a stage of a run (DB.24): the `draw-nebula-field` command draws every kept cloud of classes M and N
    of the galaxy's molecular cloud field from the seed, a unit per slab of field cells (`draw_nebula_field_unit`), after
    the earlier table was emptied (`prepare_draw_nebula_field`). A nebula already stored from the field is linked back to
    its row afterwards, and the Galaxy Map's nebula block totals are made again from the new field
    (`conclude_draw_nebula_field`).

    Attributes:
        written (int): The rows written.
    """

    key = NEBULA_DRAW_KEY
    label = "Nebula field"

    def __init__(self, args):
        super().__init__(args)
        self.written = 0
        self.slabs = self.done = 0
        self.cell_pc = None
        self.bar = None

    def params(self):
        return {"run_id": getattr(self.args, "generation_run_id", None)}

    def on_plan(self, plan):
        self.slabs, self.cell_pc = plan["detail"]["slabs"], plan["detail"]["cell_pc"]
        self.bar = self.hold(steps.Step("Drawing the nebula field", "nebula-field", max(self.slabs, 1), args=self.args,
                                        progress=self.progress, workers=self.workers))
        self.bar.update(description=f"Drawing the nebula field (0 of {self.slabs:,} slabs)")

    def on_result(self, unit, result):
        self.written += result["written"]
        for drawn in result["drawn"]:
            _log_nebula_slab(drawn["slab"], self.cell_pc, drawn)
        self.done += result["count"]
        self.bar.update(advance=result["count"], description=f"Drawing the nebula field ({self.done:,} of {self.slabs:,} slabs)")

    def on_done(self, conclusion):
        linked = conclusion["linked"]
        log.normal(f"Nebula field: {self.written:,} clouds of classes M and N"
                   + (f", {linked:,} already stored." if linked else "."))
        log.normal(f"Made {conclusion['blocks']:,} nebula block totals.")


def prepare_draw_nebula_field(params):
    """
    The `prepare` of the `draw-nebula-field` command (PERF.172), when its stage opens: empties the earlier table. A galaxy with
    no outline has nothing to draw and nothing is emptied (PERF.203, checked here and not when the run starts, because an
    earlier stage of the same run may be what builds the plan).

    Raises:
        RuntimeError: The galaxy's skeleton has never been built.
    """
    conn = store.get_connection(_mysql_config_of(params))
    try:
        if store.get_galaxy_shape(conn) is None:
            raise RuntimeError(SKELETON_NEVER_BUILT)
        if _outline(conn):
            store.clear_nebula_field(conn)
    finally:
        conn.close()


def conclude_draw_nebula_field(params, plan, results, seconds):
    """
    The `conclude` of the `draw-nebula-field` command (PERF.172): links the nebulae already stored from the field back to
    their rows (`store.link_nebula_field_stored`) and makes the nebula block totals of the Galaxy Map from the new field
    (MAP.172, `nebula_blocks.rebuild`).

    Returns:
        dict: `linked` (nebulae linked) and `blocks` (block totals made).
    """
    if not plan["slabs"]:
        return {"linked": 0, "blocks": 0}
    conn = store.get_connection(_mysql_config_of(params))
    try:
        store.record_scatter_epoch(conn, "nebula-field", params.get("run_id"))
        linked = store.link_nebula_field_stored(conn)
        blocks = nebula_blocks.rebuild(conn, store.get_galaxy_shape(conn).edge_pc)
        conn.commit()
    finally:
        conn.close()
    return {"linked": linked, "blocks": blocks}


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


def _id_batches(top, size):
    """`(low, high)` id ranges of `size` ids covering 1 to `top`, in order: the cut of a batch fill, the same at any
    worker count."""
    return [(low, min(low + size - 1, top)) for low in range(1, top + 1, size)]


def expand_fill_map_luminosity(params):
    """The `fill-map-luminosity` command's stage (PERF.106): a unit per id range of `MAP_LUMINOSITY_BATCH_IDS` (a whole
    multiple of it past `units.MAX_UNITS` ranges) of `phenomenon_scatter`, its cost the ids it covers; `detail` holds the
    highest id and the rows still empty."""
    conn = store.get_connection(_mysql_config_of(params))
    try:
        top = int(conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM phenomenon_scatter").fetchone()["n"])
        remaining = int(conn.execute("SELECT COUNT(*) AS n FROM phenomenon_scatter WHERE mass_solar IS NOT NULL"
                                     " AND map_luminosity_sol IS NULL").fetchone()["n"])
    finally:
        conn.close()
    return commands.Expansion(
        [({"database": params["database"], "low": low, "high": high}, float(high - low + 1))
         for low, high in _id_batches(top, units.coarsen(top, MAP_LUMINOSITY_BATCH_IDS))],
        {"top": top, "remaining": remaining})


def fill_map_luminosity_unit(params, report):
    """
    The `fill-map-luminosity` command's unit: `phenomenon_scatter.map_luminosity_sol` (DB.23, GEN.201) for the rows
    of ids `low` to `high` that have a mass and no value yet, one short statement and one commit, so readers never
    wait. It only touches rows still empty, so a stopped run carries on where it was. The values are
    `stellar_evolution.map_luminosity_sol` of the stored mass, worked out by the same pieces in SQL and stored as
    FLOAT, so a row agrees with the one a new scatter writes to float precision.

    Returns:
        dict: `filled`, the rows filled, and `ids`, how many ids the unit covered (for the bar).
    """
    conn = store.get_connection(_mysql_config_of(params))
    try:
        cursor = conn.execute(
            f"UPDATE phenomenon_scatter SET map_luminosity_sol = {map_luminosity_sql()}"
            " WHERE id BETWEEN ? AND ? AND mass_solar IS NOT NULL AND map_luminosity_sol IS NULL",
            (params["low"], params["high"]))
        conn.commit()
        return {"filled": cursor.rowcount, "ids": params["high"] - params["low"] + 1}
    finally:
        conn.close()


class BatchFollower(followed.Follower):
    """
    A stage of batch units over ids as a stage of a run (PERF.106, PERF.202): the wrapper selects the stage's command from
    this run's steps and cuts and queues the units, and this side follows their results with a bar of the ids they cover.
    The unit that completes the stage stores the stage's row.

    A subclass names the stage's `key`, `label`, `bar_label` and `noun` (the key of the count in a unit's result) and gives
    `plan_line(detail)`, the line to log when the stage is cut.

    Attributes:
        total (int): The sum of the units' counts.
    """

    bar_label = ""
    noun = ""

    def __init__(self, args):
        super().__init__(args)
        self.total = 0
        self.bar = None

    def plan_line(self, detail):
        raise NotImplementedError

    def on_plan(self, plan):
        log.normal(self.plan_line(plan["detail"]))
        self.bar = self.hold(steps.Step(self.bar_label, self.key, 1, args=self.args, progress=self.progress,
                                        workers=self.workers))
        self.bar.update(total=max(plan["detail"]["top"], 1))

    def on_result(self, unit, result):
        self.total += result[self.noun]
        self.bar.update(advance=result["ids"])


class MapLuminosityFollower(BatchFollower):
    """The map luminosity fill (`fill_map_luminosity_unit`) as a stage of a run."""

    key = "map-luminosity"
    label = "Map luminosity"
    bar_label = "Filling the map luminosity"
    noun = "filled"

    def plan_line(self, detail):
        return f"Map luminosity: {detail['remaining']:,} scattered rows have a mass and no value yet."

    def on_done(self, conclusion):
        log.normal(f"Filled the map luminosity of {self.total:,} scattered rows.")


def fill_map_luminosity(args):
    """
    Fills `phenomenon_scatter.map_luminosity_sol` for the rows that have a mass and no value yet, on the new queue:
    the `fill-map-luminosity` command, a unit per id range (`fill_map_luminosity_unit`), then the Galaxy Map's block light
    totals made again (`BlockLightsFollower`).

    Returns:
        int: The rows filled.
    """
    return followed.run_planned(args, planned_map_luminosity(args))


def planned_map_luminosity(args):
    """The map luminosity fill and the block lights after it as a `followed.Planned`; its result is the rows filled."""
    follower = MapLuminosityFollower(args)
    return followed.Planned([follower, BlockLightsFollower(args)], "Map luminosity", lambda: follower.total)


HABITABILITY_BATCH_SYSTEMS = 500
"""int: Systems per batch of the habitability fill: one transaction each, so readers never wait."""


def expand_rescore_habitability(params):
    """The `rescore-habitability` command's stage (PERF.106): a unit per id range of `HABITABILITY_BATCH_SYSTEMS`
    (a whole multiple of it past `units.MAX_UNITS` ranges) of `star_systems`, its cost the ids it covers; `detail` holds
    the highest id and the systems to re-score."""
    conn = store.get_connection(_mysql_config_of(params))
    try:
        top = int(conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM star_systems").fetchone()["n"])
        total = int(conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"])
    finally:
        conn.close()
    return commands.Expansion(
        [({"database": params["database"], "low": low, "high": high}, float(high - low + 1))
         for low, high in _id_batches(top, units.coarsen(top, HABITABILITY_BATCH_SYSTEMS))],
        {"top": top, "total": total})


def rescore_habitability_unit(params, report):
    """
    The `rescore-habitability` command's unit: re-derives the stored habitability numbers of the systems of ids `low`
    to `high` (GEN.217) from each body's stored conditions (`store.refresh_habitability`): the flux swing of a body in
    a double-star system and the whole score (index, six domains, tiers, equipment, note). For a galaxy generated
    before a change to the rules (GEN.148, GEN.207, GEN.215), so the planet chip, the search page and the facet counts
    agree with the description, which is always worked out from the same conditions when the page is drawn. Nothing but
    those columns is touched. The unit is one transaction, and the recomputation is the same whenever it is run, so a
    stopped run is started again and ends in the same place.

    Returns:
        dict: `updated`, the bodies updated, and `ids`, how many ids the unit covered (for the bar).
    """
    conn = store.get_connection(_mysql_config_of(params))
    try:
        updated = 0
        system_ids = [row["id"] for row in conn.execute(
            "SELECT id FROM star_systems WHERE id BETWEEN ? AND ? ORDER BY id", (params["low"], params["high"])).fetchall()]
        for system_id in system_ids:
            updated += store.refresh_habitability(conn, system_id)
        store.recount_habitable_worlds(conn, system_ids)    # MAP.143; the block totals follow once, in `fill_habitability`
        conn.commit()
        return {"updated": updated, "ids": params["high"] - params["low"] + 1}
    finally:
        conn.close()


class HabitabilityFollower(BatchFollower):
    """The habitability re-score (`rescore_habitability_unit`) as a stage of a run."""

    key = "habitability"
    label = "Habitability"
    bar_label = "Re-scoring habitability"
    noun = "updated"

    def plan_line(self, detail):
        return f"Habitability: {detail['total']:,} systems to re-score."

    def on_done(self, conclusion):
        log.normal(f"Re-scored the habitability of {self.total:,} bodies.")


def conclude_rescore_habitability(params, plan, results, seconds):
    """
    The `conclude` of the `rescore-habitability` command (PERF.202): makes the blocks' habitable-world totals follow the
    sectors' counts (MAP.143, `block_stats.rebuild`), once, when every range is done.
    """
    conn = store.get_connection(_mysql_config_of(params))
    try:
        block_stats.rebuild(conn)
        conn.commit()
    finally:
        conn.close()
    return {}


def fill_habitability(args):
    """
    Re-derives the stored habitability numbers of an existing galaxy (GEN.217) on the new queue: the
    `rescore-habitability` command, a unit per id range of systems (`rescore_habitability_unit`).

    Returns:
        int: The bodies updated.
    """
    return followed.run_planned(args, planned_habitability(args))


def planned_habitability(args):
    """The habitability rescore as a `followed.Planned`; its result is the bodies updated."""
    follower = HabitabilityFollower(args)
    return followed.Planned([follower], "Habitability", lambda: follower.total)


def refresh_block_lights(mysql_config, args=None):
    """
    MAP.259: makes the light totals of the Galaxy Map's blocks again (`db.block_lights.rebuild`) once a scatter, a
    clear or the map-luminosity fill has changed the pre-placed objects, so the nodes' glow follows them. It has its
    own step with a bar over the units, and a log line before the scan, after it and every `program_constants.PROGRESS_LOG_SECONDS`
    (PERF.189). A failure only warns: `planetgen fill-block-totals` makes the totals again.

    Returns:
        int: The block rows written (0 when it failed).
    """
    try:
        conn = store.get_connection(mysql_config)
        try:
            with run_common._generation_progress() as progress:
                log.set_console(progress.console)
                try:
                    log.normal("Making the block light totals for the Galaxy Map (this reads every pre-placed star and "
                               "phenomenon again; a full galaxy takes a while)...")
                    last_line = [time.monotonic()]
                    with steps.Step("Making the block light totals", "block-lights", None, args=args, progress=progress,
                                    force=True) as bar:
                        def on_progress(done, total):
                            bar.update(completed=done, total=total)
                            if time.monotonic() - last_line[0] >= program_constants.PROGRESS_LOG_SECONDS:
                                last_line[0] = time.monotonic()
                                log.normal(f"Making the block light totals: unit {done:,} of {total:,}.")

                        written = block_lights.rebuild(conn, on_progress=on_progress, on_phase=log.normal)
                finally:
                    log.reset_console()
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- the scatter is stored; the totals can be made again
        log.normal(f"Warning: making the block light totals failed: {exc}")
        return 0
    log.normal(f"Made {written:,} block light totals.")
    return written


def fill_block_totals(args):
    """
    Makes the totals of the Galaxy Map's blocks again, in the two tables the map's nodes read: what the generated
    sectors hold (MAP.258, `db.block_stats.rebuild`, which nothing but a sector's deletion kept up to date before) and
    the light of the pre-placed stars and scattered phenomena (MAP.259, `db.block_lights.rebuild`, a unit of shells at a
    time). The sectors' masses are made again with the first (DB.108); a sector holding a body with no stored mass is
    listed in the log, left out of the sums and never guessed.

    Returns:
        int: The block rows written, both tables together.
    """
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        log.normal("Making the block totals...")
        stats = block_stats.rebuild(conn, on_skipped=lambda addresses: log.normal(
            f"{len(addresses):,} sectors hold a body with no stored mass and are left out of the mass sums: "
            + ", ".join("-".join(str(part) for part in address) for address in addresses)))
        conn.commit()
        lights = block_lights.rebuild(conn, on_progress=lambda done, total: log.debug(f"units {done}/{total}"))
    finally:
        conn.close()
    log.normal(f"Made {stats:,} generated block totals and {lights:,} block light totals.")
    return stats + lights


def fill_nebula_blocks(args):
    """
    Makes the nebula totals of the Galaxy Map's blocks again (MAP.172, `db.nebula_blocks.rebuild`) from
    `nebula_field` and the stored nebulae, in one transaction.

    Returns:
        int: The block rows written.
    """
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        log.normal("Making the nebula block totals...")
        written = nebula_blocks.rebuild(conn, query.galaxy_edge_pc(conn))
        conn.commit()
    finally:
        conn.close()
    log.normal(f"Made {written:,} nebula block totals.")
    return written


def expand_clean_nebula_field(params):
    """The `clean-nebula-field` command's stage (PERF.106): one unit, the whole check and, with `apply`, the delete."""
    return commands.Expansion([(params, 1.0)], {})


def clean_nebula_field_unit(params, report):
    """
    The `clean-nebula-field` command's unit (GEN.205): brings an existing galaxy's stored field clouds in line with
    the current keep rates. Every stored dark nebula that is a field cloud (found in its cell's draw by its centre) is
    checked against `cloud_kept`; the ones the rule now rejects are counted per class (the dry run, always done first).
    With `apply` they are deleted, which clears the `inside_nebula_id` of everything inside them (the foreign key sets
    it to NULL), and containment is refreshed in the sectors they reached. Nothing else is touched: a nebula with no
    cell draw at its centre was not made by the field and stays. Drawing `nebula_field` again is the starter's next
    command (`draw-nebula-field`).

    Returns:
        dict: `stored` and `rejected` (per class, field clouds only), `deleted` (the count removed) and `lines`, what
            the run says about it, for the process that started the run to log (this unit may run in a worker).

    Raises:
        RuntimeError: The galaxy's skeleton has never been built.
    """
    conn = store.get_connection(_mysql_config_of(params))
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
        lines = ["Nebula field cleanup (dry run): " + (
            ", ".join(f"class {letter}: {rejected[letter]:,} of {stored[letter]:,} stored field clouds rejected"
                      for letter in sorted(stored)) or "no stored field clouds") + f"; {len(doomed):,} would be deleted."]
        if not params["apply"]:
            lines.append("Nothing was changed. Run again with --apply to delete them and draw "
                         "the nebula_field table again.")
            return {"stored": dict(stored), "rejected": dict(rejected), "deleted": 0, "lines": lines}
        sector_ids = set()
        for _nebula_id, center, radius_ly in doomed:
            sector_ids.update(store.sectors_reached_by(conn, center, ly_to_pc(radius_ly)))
        for start in range(0, len(doomed), 500):
            ids = [nebula_id for nebula_id, _c, _r in doomed[start:start + 500]]
            storage_checksum.mark_rows(conn, "nebulae", ids)
            conn.execute(f"DELETE FROM nebulae WHERE id IN ({','.join('?' * len(ids))})", ids)
            conn.commit()
            report(min(start + 500, len(doomed)), len(doomed))
        store.refresh_containment(conn, sorted(sector_ids))
        conn.commit()
        lines.append(f"Deleted {len(doomed):,} stored field clouds the current keep rates reject.")
        return {"stored": dict(stored), "rejected": dict(rejected), "deleted": len(doomed), "lines": lines}
    finally:
        conn.close()


class CleanNebulaFieldFollower(followed.Follower):
    """
    The clean-up of the stored nebula field as a stage of a run (GEN.205, PERF.122, PERF.202): the `clean-nebula-field`
    command checks the clouds (the dry run, always done and logged first) and, with `--apply`, deletes the rejected ones
    (`clean_nebula_field_unit`); this side logs what the unit says.

    Attributes:
        result (dict): The unit's result: `stored`, `rejected`, `deleted`, and the `lines` it said.
    """

    key = "nebula-field"
    label = "Nebula field"

    def __init__(self, args):
        super().__init__(args)
        self.result = None

    def params(self):
        return {"apply": bool(self.args.apply)}

    def on_plan(self, plan):
        self.hold(steps.Step("Cleaning up the nebula field", "nebula-field-clean", 1, args=self.args,
                             progress=self.progress))

    def on_result(self, unit, result):
        for line in result["lines"]:
            log.normal(line)
        self.result = result


def clean_nebula_field(args):
    """
    Brings an existing galaxy's stored field clouds in line with the current keep rates (GEN.205), on the new queue:
    the `clean-nebula-field` command checks them (the dry run, always done and logged first) and, with `--apply`,
    deletes the rejected ones (`clean_nebula_field_unit`). It draws nothing: the table is drawn again by
    `planetgen scatter --nebula-field` (PERF.121, PERF.122).

    Returns:
        dict: `stored` and `rejected` (per class, field clouds only) and `deleted` (the count removed).
    """
    try:
        return followed.run_planned(args, planned_clean_nebula_field(args))
    except units.RunStopped as exc:
        marker = "RuntimeError: "
        if marker not in str(exc):
            raise
        raise RuntimeError(str(exc).split(marker, 1)[1]) from exc


def planned_clean_nebula_field(args):
    """The nebula field check (and clean) as a `followed.Planned`; its result is `stored`, `rejected` and `deleted`."""
    follower = CleanNebulaFieldFollower(args)
    return followed.Planned([follower], "Nebula field", lambda: {key: follower.result[key] for key in ("stored", "rejected", "deleted")})


def redraw_nebula_field(args):
    """
    Draws the nebula field table again from the galaxy seed, the one thing `planetgen scatter --nebula-field` does
    (PERF.121): the earlier table is cleared and the `draw-nebula-field` command draws it (`NebulaFieldFollower`).

    Returns:
        int: The rows written (0 for a galaxy with no outline).
    """
    try:
        return followed.run_planned(args, planned_nebula_field(args))
    except units.RunStopped as exc:
        marker = "RuntimeError: "
        if marker not in str(exc):
            raise
        raise RuntimeError(str(exc).split(marker, 1)[1]) from exc


def planned_nebula_field(args):
    """The nebula field table drawn again as a `followed.Planned`; its result is the rows written."""
    follower = NebulaFieldFollower(args)
    return followed.Planned([follower], "Nebula field", lambda: follower.written)


def expand_refresh_block_lights(params):
    """The `refresh-block-lights` command's stage: one unit."""
    return commands.Expansion([(params, 1.0)], {})


def refresh_block_lights_unit(params, report):
    """
    The `refresh-block-lights` command's unit (MAP.259, PERF.172): makes the light totals of the Galaxy Map's blocks again
    (`db.block_lights.rebuild`) once a scatter, a clear or the map-luminosity fill has changed the pre-placed objects, so the
    nodes' glow follows them. It reports the units it has done (PERF.189). A failure only warns: `planetgen fill-block-totals`
    makes the totals again.

    Returns:
        dict: `written` (the block rows written, 0 when it failed), `lines` (what the run says about it) and `warning`
            (the failure, or `None`).
    """
    lines = []
    conn = store.get_connection(_mysql_config_of(params))
    try:
        written = block_lights.rebuild(conn, on_progress=report, on_phase=lines.append)
    except Exception as exc:  # noqa: BLE001 -- the scatter is stored; the totals can be made again
        return {"written": 0, "lines": lines, "warning": f"{type(exc).__name__}: {exc}"}
    finally:
        conn.close()
    return {"written": written, "lines": lines, "warning": None}


class BlockLightsFollower(followed.Follower):
    """The Galaxy Map's block light totals made again, as the last stage of a scatter (`refresh_block_lights_unit`)."""

    key = "block-lights"
    label = "Make the block light totals"
    repeatable = True

    def __init__(self, args):
        super().__init__(args)
        self.bar = None

    def on_plan(self, plan):
        log.normal("Making the block light totals for the Galaxy Map (this reads every pre-placed star and phenomenon again; "
                   "a full galaxy takes a while)...")
        self.bar = self.hold(steps.Step("Making the block light totals", "block-lights", None, args=self.args,
                                        progress=self.progress, force=True))

    def on_progress(self, unit, done, estimate):
        self.bar.update(completed=done, total=estimate)

    def on_result(self, unit, result):
        for line in result["lines"]:
            log.normal(line)
        if result["warning"]:
            log.normal(f"Warning: making the block light totals failed: {result['warning']}")
        else:
            log.normal(f"Made {result['written']:,} block light totals.")


class _PhenomenaBook:
    """What the stages of one phenomenon scatter add up to: the phenomena placed per class (the shells' and the special
    ones), the stars with companions, and the line that sums them up when the scatter's table is complete."""

    def __init__(self):
        self.shells = self.min_mass_solar = None        # the plan of the shells' stage says them (PERF.212)
        self.counts = collections.Counter()
        self.companions = 0
        self.landed = 0
        self.tasks = 0
        self.started = time.perf_counter()

    def summary(self):
        """`{"counts", "total", "elapsed_s"}`: the phenomena placed per class."""
        return {"counts": dict(self.counts), "total": sum(self.counts.values()),
                "elapsed_s": time.perf_counter() - self.started}

    def log(self):
        """Sums the scatter up, once its table is complete."""
        if self.landed:
            log.normal(f"Phenomena landed in {self.landed:,} of {self.tasks:,} runs of shells.")
        else:
            log.normal(f"No phenomena landed in any of the {self.shells:,} shells.")
        # Like the star count: only what was created.
        labels = [label for label in list(phenomenon_scatter.EXPECTED_LABELS) + sorted(
            set(self.counts) - set(phenomenon_scatter.EXPECTED_LABELS)) if self.counts.get(label)]
        listed = ": " + ", ".join(f"{self.counts[label]:,} {label}" for label in labels) if labels else ""
        summary = self.summary()
        log.normal(f"Placed {summary['total']:,} phenomena in {summary['elapsed_s']:.1f}s{listed}; "
                   f"{self.companions:,} of them have companions.")


class PhenomenaFollower(followed.Follower):
    """
    The shells of the phenomenon scatter as a stage of a run (PERF.106): the `scatter-phenomena` command's units, a
    run of shells each. The earlier scatter was emptied when the stage opened (`prepare_scatter_phenomena`).
    """

    key = "phenomena"
    label = "Phenomena"

    def __init__(self, args, book, settings):
        super().__init__(args)
        self._params, self.book = {"settings": settings}, book
        self.shape = self.bar = None
        self.shells_done = 0

    def params(self):
        return dict(self._params)

    def on_plan(self, plan):
        shape = self.shape = _PhenomenaShape.from_json(plan["detail"])
        self.book.tasks = len(shape.tasks)
        self.book.shells, self.book.min_mass_solar = plan["detail"]["outline_shells"], plan["detail"]["min_mass_solar"]
        log.normal(f"Phenomena: about {round(shape.expected_rows):,} phenomena to place in {shape.shell_total:,} shells "
                   f"({phenomenon_scatter.cut_summary(self.book.min_mass_solar)}).")
        if plan["detail"].get("filled"):
            log.normal(f"Leaving out the {plan['detail']['filled']:,} sectors already filled.")
        self.bar = self.hold(steps.Step(f"Phenomena (0 of {shape.shell_total:,} shells)", "phenomena",
                                        max(sum(shape.weights.values()), 1.0), args=self.args, progress=self.progress,
                                        workers=self.workers, percent=True, record=False))

    def on_result(self, unit, result):
        shape, book = self.shape, self.book
        task_counts = result["counts"]
        book.companions += result.get("systems", 0)
        for kind, count in task_counts.items():
            book.counts[kind] += count
        if sum(task_counts.values()):
            book.landed += 1
        note = _unit_note(unit + 1, len(shape.tasks), shape.expected_of(unit), "phenomena")
        if unit in shape.component_parts:
            if sum(task_counts.values()):
                log.normal(f"Phenomena, components{note}: placed {sum(task_counts.values()):,} phenomena: " + ", ".join(
                    f"{count:,} {kind}" for kind, count in sorted(task_counts.items())) + ".")
        else:
            _log_phenomena_shells(shape.shells_of(unit), task_counts, note)
        self.shells_done += len(shape.shells_of(unit))
        self.bar.update(advance=shape.weights[unit],
                        description=f"Phenomena ({self.shells_done:,} of {shape.shell_total:,} shells)")

    def on_done(self, conclusion):
        stages.note(shells=self.shape.shell_total, shells_modified=self.book.landed, objects=sum(self.book.counts.values()))


class PhenomenaIndexFollower(followed.Follower):
    """The end of the phenomenon scatter as a stage of a run: the special phenomena written, the indexes built, the
    hypervelocity stars stamped and the scatter recorded (`finish_phenomena_unit`)."""

    key = "phenomena-index"
    label = "Index the phenomena"

    def __init__(self, args, book, settings):
        super().__init__(args)
        self._params, self.book = {"settings": settings, "run_id": getattr(args, "generation_run_id", None)}, book
        self.bar = None

    def params(self):
        return dict(self._params)

    def on_plan(self, plan):
        self.bar = self.hold(steps.Step("Drawing the special phenomena", "phenomena-index", 1, args=self.args,
                                        progress=self.progress, force=True))

    def on_progress(self, unit, done, estimate):
        # The unit reports in quarters: the special rows drawn, written, the indexes built, the rest (UX.84).
        phase = ("Drawing the special phenomena", "Writing the special phenomena", "Indexing the phenomena",
                 "Stamping the hypervelocity stars")[min(int(done), 3)]
        self.bar.update(completed=done, total=estimate, description=phase)

    def on_result(self, unit, result):
        for line in result["lines"]:
            log.normal(line)
        special = result["special"]
        for label, count in special.items():
            self.book.counts[label] += count
        if special:
            log.normal("Special phenomena: " + ", ".join(f"{count:,} {kind}" for kind, count in sorted(special.items())) + ".")

    def on_done(self, conclusion):
        self.book.log()


def phenomena_followers(args, refuse=None, tag=""):
    """
    The stages that pre-place the galaxy's black holes, neutron stars, planetary nebulae and supernova remnants
    (`phenomenon_scatter.scatter_shell`; neutron stars and black holes only from `--compact-min-mass` up, GEN.195), its
    hypervelocity stars and its nucleus into `phenomenon_scatter`, replacing any earlier scatter, and the nebula field:
    the shells, a run of shells to a unit; the special phenomena and the indexes; and the nebula field, a slab of
    cells to a unit. Sectors already filled are left out, as the bright-star scatter does (GEN.30), and the nucleus and
    the hypervelocity stars are left out of a cell already filled.

    Args:
        refuse (str, optional): The kind of scatter that is refused when its rows are stored (`bright_star_followers`).
        tag (str): Names this scatter among the others of the run that holds it (PERF.203).

    Returns:
        tuple: `(followers, book)`: the stages' followers and the `_PhenomenaBook` that adds the scatter up.

    Raises:
        RuntimeError: The skeleton was never built.
    """
    book = _PhenomenaBook()
    followers = [PhenomenaFollower(args, book, _phenomena_settings_marker(args, "phenomena", refuse, tag)),
                 PhenomenaIndexFollower(args, book, _phenomena_settings_marker(args, "phenomena-index", refuse, tag)),
                 NebulaFieldFollower(args)]
    return followers, book


def _phenomena_settings_marker(args, stage, refuse=None, tag=""):
    """The deferred parameter that gives a stage of the phenomenon scatter its seed and lowest compact mass when it opens."""
    return wrapper.deferred("scatter-settings", merge=True, scope="phenomena", stage=stage,
                            database=wrapper.database_args(store.mysql_config_from_args(args)),
                            phenomenon_min_mass=getattr(args, "phenomenon_min_mass", None),
                            compact_min_mass=getattr(args, "compact_min_mass", None), refuse=refuse, tag=tag)


def _phenomena_settings(marker):
    """
    What the stages of the phenomenon scatter work from, read from the database when the first of them opens (PERF.212; the values
    `phenomena_followers` read when it built the followers): `{stage: params}` with the scatter's seed and the lowest mass of a neutron
    star or black hole it places.

    Raises:
        RuntimeError: The skeleton was never built.
    """
    args = types.SimpleNamespace(phenomenon_min_mass=marker["phenomenon_min_mass"], compact_min_mass=marker["compact_min_mass"])
    conn = store.get_connection(_mysql_config_of(marker))
    try:
        if marker.get("refuse"):
            _refuse_stored_scatter(conn, marker["refuse"])
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError(SKELETON_NEVER_BUILT)
        seed = _phenomenon_scatter_seed(skeleton)
        min_mass_solar = _compact_mass_limit(args, conn)
    finally:
        conn.close()
    params = {"min_mass_solar": min_mass_solar, "seed": seed}
    return {"phenomena": params, "phenomena-index": params}


def scatter_phenomena(args, refuse=None):
    """
    Pre-places the galaxy's phenomena and its nebula field (`phenomena_followers`) as one run, then makes the Galaxy
    Map's block light totals again (`BlockLightsFollower`).

    Returns:
        dict: `counts` (per kind), `total` and `elapsed_s`.
    """
    planned = planned_phenomena(args, refuse)
    planned.options["seed"] = _recorded_seed(args, _phenomenon_scatter_seed)
    return _refusals_as_errors(lambda: followed.run_planned(args, planned))


def planned_phenomena(args, refuse=None, tag=""):
    """The phenomena scatter and the block lights after it as a `followed.Planned`; its result is `book.summary()`."""
    followers, book = phenomena_followers(args, refuse, tag)
    return followed.Planned(followers + [BlockLightsFollower(args)], "Phenomena", book.summary)


def prepare_scatter_phenomena(params):
    """The `prepare` of the `scatter-phenomena` command (PERF.172), when its stage opens: empties the earlier scatter, and
    takes off the indexes: none to keep up to date while tens of millions of rows go in, built once at the end
    (`finish_phenomena_unit`, PERF.73)."""
    conn = store.get_connection(_mysql_config_of(params))
    try:
        store.clear_phenomenon_scatter(conn)
        store.drop_phenomenon_scatter_indexes(conn)
    finally:
        conn.close()


def conclude_scatter_phenomena(params, plan, results, seconds):
    """The `conclude` of the `scatter-phenomena` command (PERF.172): stores how many phenomena of each class the shells
    placed, for the Phenomena table (`store.record_phenomenon_scatter_classes`); the special phenomena add theirs."""
    counts = collections.Counter()
    for result in results:
        counts.update(result["counts"])
    conn = store.get_connection(_mysql_config_of(params))
    try:
        store.record_phenomenon_scatter_classes(
            conn, {phenomenon_scatter.split_label(label): count for label, count in counts.items()})
        conn.commit()
    finally:
        conn.close()
    return {}


def expand_finish_phenomena(params):
    """The `finish-phenomena` command's stage: one unit."""
    return commands.Expansion([(params, 1.0)], {})


def finish_phenomena_unit(params, report):
    """
    The `finish-phenomena` command's unit (PERF.172): ends the phenomenon scatter once its shells are written: draws the
    special phenomena (the nucleus and the hypervelocity stars, minus any in a cell already filled) and writes them, builds
    the table's indexes one at a time (PERF.101), stamps the hypervelocity stars' epoch and records the scatter's seed and
    mass cut, which tells a sector fill to build its objects from the table (`store.record_phenomenon_scatter`). The
    unit reports its progress in quarters: the special rows drawn, written, the indexes built, the rest.

    Returns:
        dict: `special`, how many special phenomena of each class were written, and `lines`, what the run says about it.
    """
    mysql_config = _mysql_config_of(params)
    conn = store.get_connection(mysql_config)
    lines = []
    try:
        skeleton = store.get_galaxy_shape(conn)
        extents = _outline(conn)
        filled = store.filled_sector_addresses(conn)
        reach = getattr(skeleton, "reach_shell_index", None) if program_constants.SCATTER_COVERS_REACH else None
        special = list(phenomenon_scatter.special_rows(extents, skeleton.edge_pc, params["seed"], skeleton.shape, filled,
                                                       reach))
        report(1, 4)
        store.insert_phenomenon_scatter(conn, special)
        report(2, 4)
        # One index at a time, each worth the table's rows, so the bar moves after every index and its time left comes
        # from the rows per second earlier runs recorded (PERF.101).
        rows = max(store.estimated_rows(conn, "phenomenon_scatter"), 1)
        pending = store.missing_phenomenon_scatter_indexes(conn)
        lines.append(f"Building the phenomena indexes: {len(pending)} indexes over about {rows:,} rows, one at a time; a "
                     f"full galaxy takes a while.")

        def built(name):
            lines.append(f"Built index {name}.")
            report(2 + len(lines) / (len(pending) + 1), 4)

        store.add_phenomenon_scatter_indexes(conn, on_index=built)
        store.stamp_phenomenon_scatter_epoch(conn)
        report(4, 4)
        store.record_phenomenon_scatter(conn, params["seed"], params["min_mass_solar"])
        store.record_scatter_epoch(conn, "phenomena", params.get("run_id"))
        special_counts = collections.Counter(phenomenon_scatter.class_label(row[3], row[4]) for row in special)
        store.add_phenomenon_scatter_classes(
            conn, {phenomenon_scatter.split_label(label): count for label, count in special_counts.items()})
        conn.commit()
    finally:
        conn.close()
    return {"special": dict(special_counts), "lines": lines}


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
            these, slices = componentPass.chosen(payload["components"], payload["component_parts"])
            runs.append(componentPass.scatter_phenomena(
                these, edge_pc, seed, skip_addresses=skip, min_mass_solar=min_mass, shape=shape, slices=slices))
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
    Builds and persists the galaxy's density skeleton and nothing else (PERF.121): the scatters, the fills and the
    nebula field clean-up are commands of their own (`run_scatter`, `run_fill_map_luminosity`,
    `run_rescore_habitability`, `run_clean_nebula_field`), asked for in the run's step list.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "plan"`).
    """
    stages.begin(stages.plan_stages(args), args, "plan")
    stages.enter("skeleton")
    summary = _skeleton_on_queue(args)
    log_skeleton_summary(args, summary)
    stages.finish("the galaxy has no shells" if not summary["shell_count"] else None)


def log_skeleton_summary(args, summary):
    """Says what the skeleton build (`plan_skeleton_unit`'s `summary`) made, and warns when the galaxy's edge was not reached."""
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
            f"WARNING: the galactic plane still qualified at the --edge ({dials.lever('edge').format(args.edge)}) -- the "
            f"galaxy's true edge was not reached. Re-run with a larger --edge, or none, if these shape "
            f"parameters really do produce a galaxy this large."
        )


SCATTER_NAMES = {"mass": "massive-star", "luminosity": "bright-star", "phenomena": "phenomena",
                 "nebula-field": "nebula field"}
"""dict: What each kind of scatter is called in a message."""


def _refuse_stored_scatter(conn, kind):
    """
    Stops a `planetgen scatter` of `kind` when that scatter already has rows (Boss's F7): a scatter draws where none is
    stored, and `planetgen clear-scatter` empties a kind first. The band adds to the bright stars and is not refused. It is
    made when the scatter's first stage opens (PERF.203), so a run that clears a kind and draws it again is not refused.

    Raises:
        RuntimeError: The kind's rows are stored.
    """
    if kind != "band" and store.scatter_stored(conn, kind):
        raise RuntimeError(f"The {SCATTER_NAMES[kind]} scatter is already stored -- run "
                           f"'planetgen clear-scatter --{kind}' first.")


def _refuse_stored_scatter_of(args, kind):
    """`_refuse_stored_scatter` for the database `args` name, for a scatter that is drawn from the starting process."""
    conn = store.get_connection(store.mysql_config_from_args(args))
    try:
        _refuse_stored_scatter(conn, kind)
    finally:
        conn.close()


def run_scatter(args):
    """
    `planetgen scatter` (PERF.121): draws the one kind of scatter `args` names (`stages.scatter_kind`) on the stored
    plan and nothing else. `--mass` and `--luminosity` each clear and rewrite only their own pass of the bright stars
    (`scatter_bright_stars`); `--phenomena` clears and rewrites the phenomena (`scatter_phenomena`); `--nebula-field`
    draws the nebula field table again (`redraw_nebula_field`); `--band` adds the band of dimmer stars down to
    `--down-to` (`add_bright_star_band`).
    """
    kind = stages.scatter_kind(args)
    if kind == "nebula-field":
        _refuse_stored_scatter_of(args, kind)
    stages.begin(stages.scatter_stages(args), args, "scatter")
    if kind == "band":
        stages.enter("band")
        add_bright_star_band(args)
    elif kind == "phenomena":
        scatter_phenomena(args, refuse=kind)
    elif kind == "nebula-field":
        stages.enter(NEBULA_DRAW_KEY)
        redraw_nebula_field(args)
    else:
        scatter_bright_stars(args, passes=(kind,), refuse=kind)
    stages.finish()


def expand_clear_scatter(params):
    """The `clear-scatter` command's stage: one unit."""
    return commands.Expansion([(params, 1.0)], {})


def clear_scatter_unit(params, report):
    """
    The `clear-scatter` command's unit (Boss's F7): empties the stored rows of one kind of scatter, `params["part"]`:
    "mass" (the massive stars, refused while the bright stars are stored, since they are split from these at the mass
    limit), "luminosity" (the bright stars, the band with them), "phenomena" or "nebula-field". Nothing else is touched.

    Returns:
        dict: `part` and `lines`, what the run says about it, for the process that started the run to log (this unit
            may run in a worker).

    Raises:
        RuntimeError: The massive stars are asked for while the bright stars are stored.
    """
    part = params["part"]
    name = SCATTER_NAMES[part]
    conn = store.get_connection(_mysql_config_of(params))
    try:
        if not store.scatter_stored(conn, part):
            store.forget_scatter_epoch(conn, ("luminosity", "band") if part == "luminosity" else part)   # drawn, and empty
            conn.commit()
            return {"part": part, "lines": [f"The {name} scatter has nothing stored."]}
        if part == "mass" and store.scatter_stored(conn, "luminosity"):
            raise RuntimeError("The bright stars are split from the massive stars at the mass limit -- run "
                               "'planetgen clear-scatter --luminosity' first.")
        if part in ("mass", "luminosity"):
            if store.bright_star_mass_limit(conn) is None:      # no mass pass to split the stars at: they all go
                store.clear_bright_stars(conn)
            else:
                store.clear_bright_star_pass(conn, part, store.scatter_limits_of(conn))
        elif part == "phenomena":
            store.clear_phenomenon_scatter(conn)
        else:
            store.clear_nebula_field(conn)
    finally:
        conn.close()
    if part != "nebula-field":
        refresh_block_lights(_mysql_config_of(params))
    return {"part": part, "lines": [f"Cleared the {name} scatter."]}


class ClearScatterFollower(followed.Follower):
    """
    The clearing of one scatter's stored rows as a stage of a run (Boss's F7, PERF.121, PERF.202): the `clear-scatter`
    command does it in one unit (`clear_scatter_unit`, which refuses when the bright stars are split from the massive
    ones) and this side logs what the unit says. The stage is the kind's (`stages.CLEAR_STAGE`).
    """

    def __init__(self, args):
        super().__init__(args)
        kind = stages.clear_kind(args)
        self.key, self.label = stages.CLEAR_STAGE[kind].key, stages.CLEAR_STAGE[kind].label

    def on_result(self, unit, result):
        for line in result["lines"]:
            log.normal(line)


def clear_scatter(args):
    """`planetgen clear-scatter`'s work on the new queue: the `clear-scatter` command, one unit, and its lines logged."""
    try:
        followed.run_planned(args, planned_clear_scatter(args))
    except units.RunStopped as exc:
        marker = "RuntimeError: "
        if marker not in str(exc):
            raise
        raise RuntimeError(str(exc).split(marker, 1)[1]) from exc


def planned_clear_scatter(args):
    """The scatter clear as a `followed.Planned`."""
    follower = ClearScatterFollower(args)
    return followed.Planned([follower], follower.label)


def run_clear_scatter(args):
    """`planetgen clear-scatter` (Boss's F7): `clear_scatter`, the one stage of its run."""
    kind = stages.clear_kind(args)
    stages.begin(stages.clear_stages(args), args, "clear-scatter")
    stages.enter(stages.CLEAR_STAGE[kind].key)
    clear_scatter(args)
    stages.finish()


def run_fill_map_luminosity(args):
    """`planetgen fill-map-luminosity` (PERF.122): `fill_map_luminosity`, its stage and the block light totals."""
    stages.begin(stages.for_command("fill-map-luminosity", args), args, "fill-map-luminosity")
    fill_map_luminosity(args)
    stages.finish()


def run_rescore_habitability(args):
    """`planetgen rescore-habitability` (PERF.122): `fill_habitability`, the one stage of its run."""
    stages.begin(stages.for_command("rescore-habitability", args), args, "rescore-habitability")
    fill_habitability(args)
    stages.finish()


def run_fill_nebula_blocks(args):
    """`planetgen fill-nebula-blocks` (MAP.172, PERF.122): `fill_nebula_blocks`, the one stage of its run."""
    stages.begin(stages.for_command("fill-nebula-blocks", args), args, "fill-nebula-blocks")
    stages.enter("nebula-blocks")
    fill_nebula_blocks(args)
    stages.finish()


def run_fill_block_totals(args):
    """`planetgen fill-block-totals` (MAP.259): `fill_block_totals`, the one stage of its run."""
    stages.begin(stages.for_command("fill-block-totals", args), args, "fill-block-totals")
    stages.enter("block-totals")
    fill_block_totals(args)
    stages.finish()


def run_clean_nebula_field(args):
    """`planetgen clean-nebula-field` (PERF.122): `clean_nebula_field`, the one stage of its run."""
    stages.begin(stages.for_command("clean-nebula-field", args), args, "clean-nebula-field")
    stages.enter("nebula-field")
    clean_nebula_field(args)
    stages.finish()
