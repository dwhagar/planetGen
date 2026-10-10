# planetgen/generation/run_plan.py

"""
Run Plan
========

The `plan` command: the galaxy's density skeleton (`galaxy_shape`/
`galaxy_layer`, the outline with one row per layer) that
`ensure_sector_generated` consults to decide, cheaply and exactly,
whether an address is worth generating, and the galaxy-wide bright-star
scatter and its per-sector bands.
"""

import collections
import math
import queue as queue_module
import threading
import time

from planetgen.queue import progress_rate, redisqueue, work as workQueue
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import phenomenon_scatter
from planetgen.generation import star_population as starPopulation
from planetgen.generation.star_labels import describe_types, star_label
from planetgen.galaxy import seed as galaxySeed, settings_file
from planetgen.names import naming_key
from planetgen.physics import constants
from planetgen import tuning as program_constants
from planetgen.util import draw
from planetgen.util import log
from planetgen.galaxy.density import build_galaxy_shape, relative_density
from planetgen.galaxy.geometry import ring_sector_count, sector_position_pc
from planetgen.galaxy.skeleton import (
    build_layer_extents, candidate_sector_count, expected_system_count_at_density_1,
)
from planetgen.generation import run_common
from planetgen.generation import steps


BACKFILL_CHUNK_SECTORS = 200
"""int: How many sectors one backfill transaction locks, draws and
records (GEN.44)."""


def _scatter_level_and_seed(conn, skeleton):
    """The galaxy scatter's level (`None` before any scatter) and the seed
    every sector's own draws use: the scatter's, else the galaxy seed's
    (GEN.39), else 0 for a galaxy planned before v51."""
    settings = store.bright_star_scatter_settings(conn)
    if settings is not None:
        return settings
    if skeleton.galaxy_seed is None:
        return None, 0
    return None, galaxySeed.short_seed(skeleton.galaxy_seed, "bright-stars", "scatter")


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


def _draw_sector_bands(conn, skeleton, addresses, floors, galaxy_level, seed, ceiling_cap=None, counts=None,
                       on_sector=None, layers=None):
    """
    Draws each of `addresses` down to its floor (`floors`, a dict or one
    number) from the level it holds, `BACKFILL_CHUNK_SECTORS` at a time:
    locks the chunk's `sector_stats` rows, rechecks each level under the
    lock (another run may have got there first, or filled it), wipes the
    stray stars of a sector with no level at all (a failed run, GEN.44),
    writes the stars and the new levels, and commits. `ceiling_cap` caps
    the ceiling (a band run never draws above the galaxy's old level).
    `counts`, when given, adds up the stars written per population.
    `on_sector`, when given, is called once per address visited (a
    progress bar's tick, PERF.28). `layers`, when given, adds up the stars
    written per sector layer and kind: `layers[layer_index]` is a `Counter`
    of `star_label` pairs (GEN.131).

    Returns:
        tuple: `(sectors drawn, stars written)`.
    """
    sectors = stars = 0
    e_value = skeleton.expected_system_count_at_density_1
    # GEN.185: a galaxy whose mass pass placed every heavy star draws the lighter ones here.
    mass_limit = store.bright_star_mass_limit(conn)
    mass_range = None if mass_limit is None else (None, mass_limit)
    for start in range(0, len(addresses), BACKFILL_CHUNK_SECTORS):
        chunk = addresses[start:start + BACKFILL_CHUNK_SECTORS]
        entries = []
        for address in chunk:
            density = max(relative_density(sector_position_pc(*address, skeleton.edge_pc), skeleton.shape), 0.0)
            entries.append((address, density, density * e_value))
        locked = store.lock_sector_stats(conn, entries)
        filled = store.get_occupied_addresses(conn, {address[0] for address in chunk})
        rows, new_levels = [], {}
        for address in chunk:
            floor = floors if isinstance(floors, (int, float)) else floors[address]
            level = locked.get(address)
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
                conn.execute("DELETE FROM bright_stars WHERE ring_index = ? AND layer_index = ?"
                             " AND ring_slot_index = ? AND star_system_id IS NULL", address)
            if ceiling_cap is not None:
                ceiling = ceiling_cap if ceiling is None else min(ceiling, ceiling_cap)
            rows.extend(brightStars.backfill_cells(skeleton.shape, [address], skeleton.edge_pc, e_value, floor,
                                                   ceiling, seed, mass_range=mass_range))
            new_levels[address] = floor
            if on_sector is not None:
                on_sector()
        store.insert_bright_stars(conn, rows)
        store.set_sector_bright_levels(conn, new_levels)
        conn.commit()
        if counts is not None:
            for row in rows:
                counts[row[6]] += 1
        if layers is not None:
            for row in rows:
                layers.setdefault(row[1], collections.Counter())[star_label(row[7], row[8])] += 1
        sectors += len(new_levels)
        stars += len(rows)
    return sectors, stars


def build_skeleton(args):
    """
    Builds and persists the galaxy's skeleton: `galaxy_shape` (the shape
    parameters, calibration constant, edge length and outer ring -- one
    singleton row) and `galaxy_layer` (the galaxy's outline: one row per
    layer, highest to lowest, holding the last ring that layer reaches --
    see `galaxySkeleton.build_layer_extents`), plus `galaxy_column` (the
    same outline per ring: the highest and lowest layer each ring
    reaches). Together they are the bounds every generation path checks
    first. The edge is always the standard
    `program_constants.DEFAULT_SECTOR_EDGE_PC`.

    No sector content or individual addresses are stored -- that stays
    lazy (`ensure_sector_generated`). Re-running replaces the whole
    skeleton; a full Milky-Way-scale build takes a few milliseconds.

    Args:
        args (argparse.Namespace): Parsed arguments.

    Returns:
        dict: `outer_ring_index`, `layer_count`, `top_layer_index`,
              `total_candidate_sectors`, `elapsed_s`, `edge_confirmed`.
    """
    edge_pc = run_common._edge_pc()
    e_value = expected_system_count_at_density_1(program_constants.DEFAULT_SECTOR_EDGE_LY)
    threshold_rho = 1.0 / e_value

    shape = build_galaxy_shape(
        disk_scale_length_pc=args.disk_scale_length_pc,
        disk_scale_height_pc=args.disk_scale_height_pc,
        bulge_scale_radius_pc=args.bulge_scale_radius_pc,
        bulge_amplitude=args.bulge_amplitude,
        arm_count=args.arm_count,
        pitch_angle_rad=math.radians(args.pitch_angle_deg),
        arm_amplitude=args.arm_amplitude,
        calibration_radius_pc=args.calibration_radius_pc,
    )

    log.normal(
        f"Building skeleton: edge_pc={edge_pc:g} disk_scale_length_pc={shape.disk_scale_length_pc} "
        f"disk_scale_height_pc={shape.disk_scale_height_pc} "
        f"bulge_scale_radius_pc={shape.bulge_scale_radius_pc} "
        f"bulge_amplitude={shape.bulge_amplitude} arm_count={shape.arm_count} "
        f"k_norm={shape.k_norm:.4f} threshold_rho={threshold_rho:.6f}"
    )

    t0 = time.perf_counter()
    extents, outer_ring_index, edge_confirmed = build_layer_extents(
        shape, edge_pc, threshold_rho, max_ring=args.max_ring,
    )
    elapsed = time.perf_counter() - t0

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
    finally:
        conn.close()
    store.replace_galaxy_layers(extents, config=mysql_config)
    galaxy_seed = store.save_galaxy_shape(
        shape, edge_pc=edge_pc, outer_ring_index=outer_ring_index,
        expected_system_count_at_density_1=e_value, config=mysql_config, galaxy_seed=requested_seed,
    )
    if requested_seed is not None:
        source = "from --seed"
    elif stored_seed is not None:
        source = "kept from the earlier plan"
    else:
        source = "drawn at random"
    log.normal(f"Galaxy seed {galaxySeed.format_seed(galaxy_seed)} ({source}).")
    key = _draw_naming_key(mysql_config, galaxy_seed, new_seed=stored_seed != galaxy_seed)
    _write_settings_file(args, galaxy_seed, key, edge_pc, outer_ring_index, e_value, shape)

    return {
        "outer_ring_index": outer_ring_index,
        "layer_count": len(extents),
        "top_layer_index": extents[0][0] if extents else None,
        "total_candidate_sectors": candidate_sector_count(extents),
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


def _write_settings_file(args, galaxy_seed, naming_key_value, edge_pc, outer_ring_index, e_value, shape):
    """
    ADM.18: the creation settings, seed, version and word lists as a JSON
    file for the Admin dashboard (`galaxy/settings_file.py`). A plan that
    changes nothing writes nothing; one that changes a setting keeps the old
    file as a dated backup. A file that can't be written only warns: the
    plan itself is done.
    """
    try:
        document = settings_file.build(
            vars(args), galaxy_seed, naming_key=naming_key_value,
            shape={"edge_pc": edge_pc, "outer_ring_index": outer_ring_index,
                   "expected_system_count_at_density_1": e_value, "k_norm": shape.k_norm})
        entry = settings_file.write(document, galaxy_seed)
    except Exception as exc:  # noqa: BLE001 -- no writable folder, or no word lists installed
        log.normal(f"Warning: the settings file was not written: {exc}")
        return
    log.normal(f"Settings file {entry['path']}.")


def scatter_bright_stars(args):
    """
    Pre-places the galaxy's bright stars on the stored plan into
    `bright_stars`, replacing any earlier scatter, in the three star passes
    of GEN.185 (pass 1, the phenomena above the mass limit, and pass 5, the
    other phenomena, are `scatter_phenomena`'s):

    2. the mass pass: every star born at or above the mass limit
       (`--phenomenon-min-mass`, one of 8 to 20 solar masses), whatever its
       luminosity (`brightStars.scatter` with `mass_range=(limit, None)`);
    3. the marks: each sector holding a pass 2 star at least as bright as
       the luminosity floor is marked (`store.bright_star_marked_addresses`);
    4. the luminosity pass: every star born lighter than the mass limit
       and at least `--bright-star-min-luminosity` bright, skipping the
       marked sectors, which already hold a star that bright.

    Each layer is one task per pass, with a progress bar. Sectors already
    filled are always left out (GEN.30): the scatter never adds stars to a
    generated sector. `--force` is still accepted and does nothing more.

    Returns:
        dict: `counts` (per population), `total` and `elapsed_s`.
    """
    mysql_config = store.mysql_config_from_args(args)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        if skeleton is None:
            raise RuntimeError("The galaxy's skeleton has never been built -- run 'planetgen plan' first.")
        extents = store.get_galaxy_layers(conn)
        filled = store.filled_sector_addresses(conn)
        if filled:
            log.normal(f"Leaving out the {len(filled):,} sectors already filled.")
        min_luminosity_sol = float(args.bright_star_min_luminosity)
        mass_limit = _phenomenon_min_mass(args, conn)
        seed = _bright_star_seed(skeleton, "scatter")
        mass_seed = _bright_star_seed(skeleton, "mass-scatter")
        store.clear_bright_stars(conn)

        t0 = time.perf_counter()
        counts = _scatter_layers(args, mysql_config, skeleton, extents, filled,
                                 starPopulation.MASS_PASS_MIN_LUMINOSITY_SOL, mass_seed, "Massive stars",
                                 mass_range=(mass_limit, None))
        massive = sum(counts.values())
        marked = store.bright_star_marked_addresses(conn, mass_limit, min_luminosity_sol * constants.SOLAR_LUMINOSITY)
        log.normal(f"{len(marked):,} sectors hold a star born at {mass_limit:g} solar masses or more and at least "
                   f"{min_luminosity_sol:g} L_sun bright; the luminosity pass skips them.")
        light = _scatter_layers(args, mysql_config, skeleton, extents, set(filled) | marked, min_luminosity_sol, seed,
                                "Bright stars", mass_range=(None, mass_limit))
        for population, count in light.items():
            counts[population] += count
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


def _bright_star_seed(skeleton, address):
    """
    The 63-bit seed of a bright-star scatter (`"scatter"`) or band, which
    its layers and the backfill's blocks derive their own streams from:
    the top 63 bits of the unit seed `bright-stars:<address>` (GEN.39),
    since `galaxy_shape.bright_star_seed` is a `BIGINT UNSIGNED`. A galaxy
    with no seed (planned before schema v51) draws one from the run's
    `random` stream.
    """
    if skeleton.galaxy_seed is None:
        return draw.getrandbits(63)
    return galaxySeed.short_seed(skeleton.galaxy_seed, "bright-stars", address)


class _LayerTracker:
    """
    The bright-star scatter's progress (PERF.4, PERF.9): the main bar
    counts expected work (`brightStars.layer_weight`, in stars) rather
    than layers, credited as each layer in progress reports its stars,
    so its ETA follows the galaxy's shape; and while layers finish
    slower than one per `SLOW_LAYER_SECONDS`, a second bar under it
    shows the stars of the layers being drawn, done of their estimate,
    with their own ETA. The second bar goes again once layers finish
    faster than one per `FAST_LAYER_SECONDS` (the gap keeps it from
    flashing on and off near the line). Called from the run's own
    thread and the progress channel's (`_drain_channel`), so every
    change holds `lock`.
    """

    SLOW_LAYER_SECONDS = 30.0
    FAST_LAYER_SECONDS = 20.0

    def __init__(self, progress, label, weights, clock=time.monotonic, prior=None, bar=None):
        self.progress = progress
        self.bar = bar
        self.label = label
        self.weights = weights
        self.clock = clock
        self.lock = threading.Lock()
        self.in_flight = {}
        self.credited = {}
        self.finished = set()
        self.done_weight = 0.0
        self.done_layers = 0
        self.layer_rate = progress_rate.DecayingRate(clock=clock)
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
        return f"{self.label} ({self.done_layers:,} of {len(self.weights):,} layers)"

    def layer_progress(self, layer_index, done, estimate):
        """A layer in progress has drawn `done` of about `estimate` stars.
        A report that arrives after its layer finished (reports come
        through the channel's thread) is ignored, or that layer would be
        counted twice (PERF.23)."""
        with self.lock:
            if layer_index in self.weights and layer_index not in self.finished:
                self.in_flight[layer_index] = (done, estimate)
                self._refresh()

    def layer_done(self, layer_index):
        with self.lock:
            if layer_index in self.finished:
                return
            self.finished.add(layer_index)
            self.in_flight.pop(layer_index, None)
            self.credited.pop(layer_index, None)
            self.done_weight += self.weights.get(layer_index, 0.0)
            self.done_layers += 1
            self.last_done = self.clock()
            self.layer_rate.add(1)
            self._refresh()

    def slow(self):
        """Whether layers are finishing slower than one per
        `SLOW_LAYER_SECONDS` (or, once the second bar shows, not yet
        faster than one per `FAST_LAYER_SECONDS`)."""
        limit = self.FAST_LAYER_SECONDS if self.detail is not None else self.SLOW_LAYER_SECONDS
        since = self.clock() - self.last_done
        rate = self.layer_rate.rate
        return since > limit or (rate is not None and rate < 1.0 / limit)

    def _refresh(self):
        partial = 0.0
        for layer_index, (done, estimate) in self.in_flight.items():
            share = min(done / estimate, 1.0) if estimate > 0 else 0.0
            # Never back: an estimate that grows doesn't take credit away.
            credit = max(self.credited.get(layer_index, 0.0), share * self.weights[layer_index])
            self.credited[layer_index] = credit
            partial += credit
        if self.bar is not None:
            self.bar.update(completed=self.done_weight + partial, description=self._description())
        else:
            self.progress.update(self.task, completed=self.done_weight + partial, description=self._description())
        if self.in_flight and self.slow():
            done = sum(item[0] for item in self.in_flight.values())
            estimate = sum(max(item[0], item[1]) for item in self.in_flight.values())
            layers = sorted(self.in_flight, key=lambda layer: (abs(layer), layer))
            names = ", ".join(str(layer) for layer in layers[:4]) + (", ..." if len(layers) > 4 else "")
            description = f"  Layer{'s' if len(layers) > 1 else ''} {names}: stars"
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
    """The progress channel of a scatter run with one worker: its layers
    are drawn in this process, so reports go straight to the tracker."""

    def __init__(self, tracker):
        self.tracker = tracker

    def put(self, item):
        self.tracker.layer_progress(*item)


def _drain_channel(channel, tracker, stop):
    """Hands the workers' layer reports to the tracker until `stop`."""
    while True:
        try:
            item = channel.get(timeout=0.25)
        except queue_module.Empty:
            if stop.is_set():
                return
            continue
        except (EOFError, OSError):
            return
        tracker.layer_progress(*item)


CHANNEL_INTERVAL_SECONDS = 0.25
"""float: How often a worker drawing a layer reports its stars (PERF.4)."""


def _layer_slots(outer_ring):
    """Sector slots in rings 0 to `outer_ring` of one layer."""
    return sum(ring_sector_count(ring_index) for ring_index in range(outer_ring + 1))


def log_layers(label, layers):
    """`_log_layer` for each layer of a sector-by-sector draw (`_draw_sector_bands`'s `layers`), top layer first."""
    for layer_index in sorted(layers):
        _log_layer(label, layer_index, layers[layer_index])


def _log_layer(label, layer_index, types):
    """One line saying how many stars a layer was given, by kind (GEN.131);
    nothing for a layer that drew none (GEN.79's note says how many did)."""
    total = sum(types.values())
    if total:
        log.normal(f"{label}, layer {layer_index}: added {total:,} stars: {describe_types(types)}.")


def _scatter_layers(args, mysql_config, skeleton, extents, filled, min_luminosity_sol, seed, label,
                    max_luminosity_sol=None, mass_range=None):
    """
    Draws and writes the bright stars of every layer through the work
    queue, one task per layer, with a progress bar: every star at or
    above `min_luminosity_sol`, or only those below `max_luminosity_sol`
    too (one band of a staged scatter), and born in `mass_range` (initial
    solar masses `(low, high)`, an end `None` for no limit; GEN.185). Cells
    in `filled` (the sectors already filled, and the ones the mass pass
    marked) are left out.

    The bar counts each layer's expected work (PERF.9,
    `brightStars.layer_weight`, worked out first in a second or two),
    credited star by star as the workers report (PERF.4, `_LayerTracker`,
    which also adds the second bar for slow layers). Each layer's time
    goes into the speed stats as a "scatter" task (PERF.10).

    Returns:
        dict: Stars written per population.
    """
    counts = {population: 0 for population in brightStars.POPULATIONS}
    drew = set()
    # Densest layers (nearest the plane) first, so no worker is left
    # with a big one at the end while the others sit idle.
    layers = sorted(extents, key=lambda extent: (abs(extent[0]), extent[0]))
    fractions = brightStars.band_fractions(min_luminosity_sol, max_luminosity_sol, mass_range)
    e_value = skeleton.expected_system_count_at_density_1
    skip_by_layer = collections.defaultdict(set)
    for address in filled:
        skip_by_layer[address[1]].add(address)
    weights, expected = {}, {}
    for layer_index, outer_ring in layers:
        weights[layer_index], expected[layer_index] = brightStars.layer_weight(
            skeleton.shape, layer_index, outer_ring, skeleton.edge_pc, e_value, fractions,
        )
    log.normal(f"{label}: about {round(sum(expected.values())):,} stars expected in {len(layers):,} layers.")
    band_share = sum(fractions.values())
    with run_common._generation_progress() as progress:
        log.set_console(progress.console)
        bar, finished = None, False
        try:
            # PERF.33: the bar starts from the stars a second this server recorded for this many workers.
            prior = run_common._generation_stats(args).pool_rate("scatter", run_common._worker_count(args))
            bar = steps.Step(label, "scatter", max(sum(weights.values()), 1.0), args=args, progress=progress,
                             workers=run_common._worker_count(args), percent=True, record=False,
                             prior=prior).__enter__()
            tracker = _LayerTracker(progress, label, weights, prior=prior, bar=bar)
            outer_rings = dict(layers)

            def on_done_for(layer_index):
                def layer_done(result, seconds, _weight):
                    layer_counts = result["counts"]
                    _log_layer(label, layer_index, collections.Counter(
                        {(name, plural): count for name, plural, count in result["types"]}))
                    for population, count in layer_counts.items():
                        counts[population] += count
                    tracker.layer_done(layer_index)
                    stars = sum(layer_counts.values())
                    if stars:
                        drew.add(layer_index)
                    slots = _layer_slots(outer_rings[layer_index])
                    density = expected[layer_index] / (e_value * band_share * slots) if band_share and slots else 0.0
                    run_common._generation_stats(args).record("scatter", density, seconds, systems=stars, stars=stars,
                                                                 workers=run_common._worker_count(args))
                return layer_done

            stop = threading.Event()
            channel = drain = None
            try:
                with run_common._work_queue(args, label) as queue:
                    if queue.parallel:
                        channel = queue.channel("scatter-progress")
                        drain = threading.Thread(target=_drain_channel, args=(channel, tracker, stop),
                                                 name="scatter-progress", daemon=True)
                        drain.start()
                    else:
                        channel = _DirectChannel(tracker)
                    queue.expect(len(layers))
                    for layer_index, outer_ring in layers:
                        payload = {
                            "mysql_config": mysql_config, "shape": skeleton.shape, "layer_index": layer_index,
                            "outer_ring": outer_ring, "edge_pc": skeleton.edge_pc, "expected": e_value,
                            "min_luminosity_sol": min_luminosity_sol, "max_luminosity_sol": max_luminosity_sol,
                            "seed": seed, "skip": skip_by_layer.get(layer_index, set()), "mass_range": mass_range,
                            "channel": channel,
                        }
                        queue.submit("bright-stars", f"layer {layer_index}", _scatter_layer_task, payload,
                                     weight=weights[layer_index], on_done=on_done_for(layer_index))
            finally:
                stop.set()
                if drain is not None:
                    drain.join(timeout=5)
                if isinstance(channel, redisqueue.Channel):
                    channel.close()
                run_common._finish_stats(args)
            finished = True
        finally:
            if bar is not None:
                bar.close(success=finished)
            log.reset_console()
    # Every layer is walked, but only the ones that drew a star count as
    # holding any (GEN.79).
    if drew:
        log.normal(f"{label}: stars landed in {len(drew):,} of {len(layers):,} layers "
                   f"(layers {min(drew)} to {max(drew)}).")
    else:
        log.normal(f"{label}: no stars landed in any of the {len(layers):,} layers.")
    return counts


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
            log.normal("No bright stars are scattered yet, so there is no layer to go below. Run "
                       f"'planetgen plan --bright-stars-only --bright-star-min-luminosity {target:g}' instead.")
            return None
        current, first_seed = float(settings[0]), settings[1]
        if target >= current:
            log.normal(f"Nothing to add: every star of {current:g} L_sun or more is already placed, "
                       f"and {target:g} is not below that.")
            return None
        extents = store.get_galaxy_layers(conn)
        filled = store.filled_sector_addresses(conn)
        # Sectors a backfill took to their own level (GEN.44) already hold
        # part of the band: the layer scatter leaves them out and each gets
        # only what it lacks below its level, sector by sector.
        own = store.sector_bright_level_keys(conn)
        skip = set(filled) | set(own)
        # GEN.32: a band run that stopped part way left the band in the
        # layers it finished; this run draws the whole band again, so it
        # starts from none of it.
        stale = store.delete_unfinished_band(conn, current * constants.SOLAR_LUMINOSITY, skip)
        conn.commit()
        if stale:
            log.normal(f"Removed {stale:,} bright stars an unfinished earlier run left below {current:g} L_sun.")
        seed = _bright_star_seed(skeleton, f"band/{target:g}-{current:g}")
        log.normal(f"Adding the bright stars from {target:g} up to {current:g} L_sun"
                   + (f", leaving out {len(filled):,} filled sectors" if filled else "")
                   + (f" and {len(own):,} backfilled sectors" if own else "") + ".")
        t0 = time.perf_counter()
        counts = _scatter_layers(args, mysql_config, skeleton, extents, skip, target, seed,
                                 f"Bright stars {target:g}-{current:g} L_sun", max_luminosity_sol=current,
                                 mass_range=mass_range)
        topped = sorted(address for address, level in own.items() if address not in filled and level > target)
        if topped:
            # Its own bar and ETA, like the backfill's (ADM.26): the web
            # page shows the progress file, and this phase can take minutes.
            with run_common._generation_progress() as progress:
                log.set_console(progress.console)
                try:
                    with steps.Step("Topping up backfilled sectors", "topup", len(topped), args=args,
                                    progress=progress) as bar:
                        topped_layers = {}
                        _sectors, stars = _draw_sector_bands(conn, skeleton, topped, target, current, first_seed,
                                                             ceiling_cap=current, counts=counts, layers=topped_layers,
                                                             on_sector=bar.advance)
                finally:
                    log.reset_console()
            log_layers("Topping up backfilled sectors", topped_layers)
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


def _scatter_layer_task(payload):
    """
    One layer of the bright-star scatter -- a work queue task (PERF.7):
    draws the layer (`brightStars.scatter_layer`, its own random stream)
    and writes its stars, committing every 10,000.

    Returns:
        dict: `counts` (stars written per population) and `types` (stars
            per kind, `[label, plural, count]` triples, GEN.131).
    """
    counts = {population: 0 for population in brightStars.POPULATIONS}
    types = collections.Counter()
    channel = payload.get("channel")
    layer_index = payload["layer_index"]
    last = [0.0]

    def report(done, estimate):
        # PERF.4: the stars drawn so far, a few times a second.
        now = time.monotonic()
        if channel is not None and now - last[0] >= CHANNEL_INTERVAL_SECONDS:
            last[0] = now
            try:
                channel.put((layer_index, done, estimate))
            except (EOFError, OSError):
                pass

    conn = store.get_connection(payload["mysql_config"])
    try:
        batch = []
        for row in brightStars.scatter_layer(
            payload["shape"], layer_index, payload["outer_ring"], payload["edge_pc"],
            payload["expected"], payload["min_luminosity_sol"], payload["seed"], skip_addresses=payload["skip"],
            max_luminosity_sol=payload.get("max_luminosity_sol"), on_progress=report,
            mass_range=payload.get("mass_range"),
        ):
            counts[row[6]] += 1
            types[star_label(row[7], row[8])] += 1
            batch.append(row)
            if len(batch) >= 10000:
                store.insert_bright_stars(conn, batch)
                conn.commit()
                batch = []
        if batch:
            store.insert_bright_stars(conn, batch)
        conn.commit()
    finally:
        conn.close()
    return {"counts": counts, "types": [[label, plural, count] for (label, plural), count in types.items()]}


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


def _phenomenon_min_mass(args, conn):
    """The mass limit (GEN.167, GEN.183): `--phenomenon-min-mass`, else the
    one already stored with the galaxy, else
    `tuning.PHENOMENON_MIN_MASS_SOLAR`."""
    value = getattr(args, "phenomenon_min_mass", None)
    if value is None:
        stored = store.phenomenon_scatter_settings(conn)
        value = stored[1] if stored is not None else None
    return float(program_constants.PHENOMENON_MIN_MASS_SOLAR if value is None else value)


def _log_phenomena_layer(layer_index, counts):
    """One line saying how many phenomena a layer was given, by kind, as `_log_layer` does for stars;
    nothing for a layer that drew none."""
    total = sum(counts.values())
    if total:
        log.normal(f"Phenomena, layer {layer_index}: placed {total:,}: "
                   + ", ".join(f"{count:,} {kind}" for kind, count in sorted(counts.items())) + ".")


def _scatter_phenomena_layers(args, bar, layers, weights, layer_done, mysql_config, skeleton, e_value, seed,
                              min_mass_solar, filled):
    """The layers of the phenomenon scatter through the work queue, the step `bar` credited with each
    layer's expected work as it finishes (UX.83, PERF.51)."""
    layers_done = [0]

    def done_with(layer_index):
        def on_done(layer_counts, seconds, weight):
            layer_done(layer_counts, seconds, weight)
            _log_phenomena_layer(layer_index, layer_counts)
            layers_done[0] += 1
            bar.update(advance=weights[layer_index],
                       description=f"Phenomena ({layers_done[0]:,} of {len(layers):,} layers)")
        return on_done

    with run_common._work_queue(args, "Phenomena") as queue:
        queue.expect(len(layers))
        for layer_index, outer_ring in layers:
            payload = {
                "mysql_config": mysql_config, "shape": skeleton.shape, "layer_index": layer_index,
                "outer_ring": outer_ring, "edge_pc": skeleton.edge_pc, "expected": e_value, "seed": seed,
                "min_mass_solar": min_mass_solar,
                "skip": {address for address in filled if address[1] == layer_index},
            }
            queue.submit("phenomena", f"layer {layer_index}", _phenomenon_layer_task, payload,
                         weight=weights[layer_index], on_done=done_with(layer_index))


def scatter_phenomena(args):
    """
    Pre-places the galaxy's black holes, neutron stars, planetary nebulae
    and supernova remnants (`phenomenon_scatter.scatter_layer`; neutron
    stars and black holes only from `--phenomenon-min-mass` up), its
    hypervelocity stars and its nucleus into `phenomenon_scatter`, replacing
    any earlier scatter, one layer per task. Sectors already filled are left
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
        extents = store.get_galaxy_layers(conn)
        filled = store.filled_sector_addresses(conn)
        seed = _phenomenon_scatter_seed(skeleton)
        min_mass_solar = _phenomenon_min_mass(args, conn)
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step("Clearing the earlier phenomena scatter", "phenomena-clear",
                                max(store.estimated_rows(conn, "phenomenon_scatter"), 1), args=args,
                                progress=progress):
                    store.clear_phenomenon_scatter(conn)
            finally:
                log.reset_console()

        t0 = time.perf_counter()
        counts = {}
        e_value = skeleton.expected_system_count_at_density_1
        layers = sorted(extents, key=lambda extent: (abs(extent[0]), extent[0]))
        weights = {layer_index: phenomenon_scatter.layer_expected(skeleton.shape, layer_index, outer_ring,
                                                                  skeleton.edge_pc, e_value, min_mass_solar)
                   + brightStars.RING_WEIGHT_STARS * (outer_ring + 1)
                   for layer_index, outer_ring in layers}
        log.normal(f"Phenomena: about {round(sum(weights.values())):,} to place in {len(layers):,} layers "
                   f"(neutron stars and black holes from {min_mass_solar:g} solar masses).")
        if filled:
            log.normal(f"Leaving out the {len(filled):,} sectors already filled.")
        landed = [0]

        def layer_done(layer_counts, seconds, weight):
            for kind, count in layer_counts.items():
                counts[kind] = counts.get(kind, 0) + count
            if sum(layer_counts.values()):
                landed[0] += 1
            # Recorded in the bar's own units (the layer's weight), so the next run's bar can start from it.
            run_common._generation_stats(args).record("phenomena", 0.0, seconds, systems=weight, stars=weight,
                                                      workers=run_common._worker_count(args))

        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step(f"Phenomena (0 of {len(layers):,} layers)", "phenomena",
                                max(sum(weights.values()), 1.0), args=args, progress=progress,
                                workers=run_common._worker_count(args), percent=True, record=False) as bar:
                    _scatter_phenomena_layers(args, bar, layers, weights, layer_done, mysql_config, skeleton,
                                              e_value, seed, min_mass_solar, filled)
            finally:
                log.reset_console()
        with run_common._generation_progress() as progress:
            log.set_console(progress.console)
            try:
                with steps.Step("Drawing the special phenomena", "phenomena-special", max(len(extents), 1), args=args,
                                progress=progress):
                    special = list(phenomenon_scatter.special_rows(extents, skeleton.edge_pc, seed, filled))
                with steps.Step("Writing the special phenomena", "phenomena-insert", max(len(special), 1), args=args,
                                progress=progress):
                    store.insert_phenomenon_scatter(conn, special)
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
        conn.commit()
    finally:
        conn.close()
    elapsed = time.perf_counter() - t0
    total = sum(counts.values())
    if landed[0]:
        log.normal(f"Phenomena landed in {landed[0]:,} of {len(layers):,} layers.")
    else:
        log.normal(f"No phenomena landed in any of the {len(layers):,} layers.")
    labels = list(phenomenon_scatter.EXPECTED_LABELS) + sorted(set(counts) - set(phenomenon_scatter.EXPECTED_LABELS))
    log.normal(f"Placed {total:,} phenomena in {elapsed:.1f}s: "
               + ", ".join(f"{counts.get(label, 0):,} {label}" for label in labels) + ".")
    return {"counts": counts, "total": total, "elapsed_s": elapsed}


def _phenomenon_layer_task(payload):
    """
    One layer of the phenomenon scatter -- a work queue task: draws the
    layer (`phenomenon_scatter.scatter_layer`) and writes its rows,
    committing every 10,000.

    Returns:
        dict: Rows written per kind.
    """
    counts = {}
    conn = store.get_connection(payload["mysql_config"])
    try:
        batch = []
        for row in phenomenon_scatter.scatter_layer(
            payload["shape"], payload["layer_index"], payload["outer_ring"], payload["edge_pc"],
            payload["expected"], payload["seed"], skip_addresses=payload["skip"],
            min_mass_solar=payload["min_mass_solar"],
        ):
            label = phenomenon_scatter.class_label(row[3], row[4])
            counts[label] = counts.get(label, 0) + 1
            batch.append(row)
            if len(batch) >= 10000:
                store.insert_phenomenon_scatter(conn, batch)
                conn.commit()
                batch = []
        if batch:
            store.insert_phenomenon_scatter(conn, batch)
        conn.commit()
    finally:
        conn.close()
    return counts


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
    if getattr(args, "bright_stars_down_to", None) is not None:
        with workQueue.job_node("bright-stars", f"Bright stars down to {args.bright_stars_down_to:g} L_sun"):
            add_bright_star_band(args)
        return
    if getattr(args, "phenomena_only", False):
        with workQueue.job_node("phenomena", "Phenomena"):
            scatter_phenomena(args)
        return
    if getattr(args, "bright_stars_only", False):
        with workQueue.job_node("bright-stars", "Bright stars"):
            scatter_bright_stars(args)
        return
    with workQueue.job_node("skeleton", "Galaxy skeleton"):
        summary = build_skeleton(args)
    if summary["layer_count"]:
        layers = f"{summary['layer_count']} layers ({summary['top_layer_index']} to {-summary['top_layer_index']})"
    else:
        layers = "no layers (nothing clears the one-star-per-sector threshold)"
    log.normal(
        f"Skeleton built in {summary['elapsed_s']:.2f}s: {layers}, outer edge = ring "
        f"{summary['outer_ring_index']}, ~{summary['total_candidate_sectors']:,} candidate sectors."
    )
    if not summary["edge_confirmed"]:
        log.normal(
            f"WARNING: the galactic plane still qualified at --max-ring ({args.max_ring}) -- the "
            f"galaxy's true edge was not reached. Re-run with a larger --max-ring if these shape "
            f"parameters really do produce a galaxy this large."
        )
    if summary["layer_count"] and not getattr(args, "no_bright_stars", False):
        with workQueue.job_node("phenomena", "Phenomena"):
            scatter_phenomena(args)
        with workQueue.job_node("bright-stars", "Bright stars"):
            scatter_bright_stars(args)
