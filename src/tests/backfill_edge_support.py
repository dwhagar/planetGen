# tests/backfill_edge_support.py

"""
The proof of the backfill's edge by sector address (OPS.242, the narrower edges, "backfill waits per cell for the sectors of that cell").

`scatter-backfill` draws the bright stars of the sectors around the sectors a run made (`run_galaxy.backfill_targets`): one stage, one list of sectors, each taken down to the mass its distance from the nearest
new sector gives it (`run_galaxy.backfill_step_targets`). A sector the run made itself is never in the list: it is charted, so the draw skips it (`store.get_occupied_addresses`, rechecked under the lock of its
`sector_stats` row). What the draw reads is the galaxy's shape and its mass limit (the scatter, a stage before), the sky tier ledger (the sky units, another), the sector's own stored mass and level, and whether the
sector is charted; what it writes is that sector's unbuilt bright stars and its `sector_stats` row. So the edge is one sector and no more: backfill of the sector A waits for the save of A when A is among the
sectors the run saves, and for nothing else of the run.

This module proves that on a small toy galaxy of three sectors, side by side:

* **Plan order** saves every sector, then runs one backfill unit with the stage's list (`backfill_targets` around the three).
* **By address** runs one backfill unit for each saved sector `S`: the sectors a face away from `S` out to the last step, each at the mass of that step as seen from `S` alone, which is all that is known when only `S`
  is saved. The unit of `S` is released when `S` is saved and when every sector of its list that the run saves is saved (`needs`). The saves keep their order (the ids the sectors take). Orders at random on one
  thread and on a pool of threads, every unit a real `engine.run` of the command.
* **The verdict** is as in `tests/address_edge_support.py`: every run is dumped and compared with plan order row by row (the ids no output depends on are left out), and each unit's SQL is read by the tap and must
  be what `queue/touches.py` declares for the command. Two sectors a step from two of the saved ones are reached from both: the draw to the lower mass after the draw to the higher is the draw to the lower.
* **A scope that is too small must fail.** `prove(..., scope="none")` releases each unit at once, before any save: the sectors of the run that lie in its list are filled with stars before they are saved, which is not
  what plan order stores.
"""

import concurrent.futures
import copy
import json
import random
import threading

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import spherical
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import run_galaxy, run_plan
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.queue import engine, sqltap, touches, wrapper
from tests import cut_neutral_support as support
from tests import test_bright_star_scatter as toy
from tests.address_edge_support import Verdict, _comparable
from tests.fuzz_support import deterministic_entropy

BACKFILL, SAVE = "scatter-backfill", "save"
SEEDS = (21, 22, 23)
"""tuple: The seeds of the random orders of the single-thread runs."""
WORKERS = 15
"""int: The size of the pool of the threaded runs."""

_SECTOR_CACHE = {}


class Recipe:
    """The comparison: the id of a bright star is the order the rows landed in, which no output depends on. The members and nodes of its stars refer to it by `source_id`, so they are compared through a join with the star's seed."""

    ignore_columns = (("bright_stars", "id"), ("scatter_members", "source_id"), ("scatter_nodes", "source_id"))
    views = (("members of the bright stars",
              "SELECT b.seed AS star_seed, m.source, m.member_no, m.parent_node_no, m.slot, m.kind, m.state, m.flags, m.offset_x_au, m.offset_y_au, m.offset_z_au, m.velocity_x_kms, m.velocity_y_kms, m.velocity_z_kms,"
              " m.decay_age_yr, m.mass_solar, m.initial_mass_solar, m.seed, m.epoch_unix, m.radius_au, m.radius_max_au, m.mass_max_solar"
              " FROM scatter_members m JOIN bright_stars b ON m.source = 1 AND b.id = m.source_id"),
             ("nodes of the bright stars",
              "SELECT b.seed AS star_seed, n.source, n.node_no, n.parent_node_no, n.slot, n.mass_solar, n.separation_km, n.eccentricity, n.inclination_deg, n.ascending_node_deg, n.arg_periapsis_deg,"
              " n.mean_anomaly_deg, n.epoch_unix, n.period_years, n.kind, n.interaction, n.flags"
              " FROM scatter_nodes n JOIN bright_stars b ON n.source = 1 AND b.id = n.source_id"))


def _addresses(config):
    """Three sectors in a row, each a face away from the one before, inside the toy galaxy's outline: `[(shell, band, wedge)]`, found with the code's own neighbor rule."""
    conn = store.get_connection(config)
    try:
        bounds = store.get_galaxy_bounds(conn)
    finally:
        conn.close()
    found = [(4, 0, 5)]
    while len(found) < 3:
        found.append(next(address for address in sorted(spherical.neighbor_addresses(*found[-1])) if address not in found and bounds.qualifies(address[0], address[1])))
    return found


def _build(address):
    """The `SpaceSector` of one address, made once (a pure function of the address) and copied for each run: one system near its center."""
    if address not in _SECTOR_CACHE:
        with deterministic_entropy(7000 + address[2]):
            center = spherical.sector_position_pc(*address, toy.EDGE_PC)
            sector = SpaceSector(f"Backfill Edge {address[0]}-{address[1]}-{address[2]}", edge_ly=tuning.DEFAULT_SECTOR_EDGE_LY)
            cfg = SystemConfig()
            cfg.STAR_TYPE = "G2V"
            cfg.BINARY_SYSTEM = False
            system = StarSystem(system_config=cfg)
            system.name = f"Backfill star {address[2]}"
            sector.add_system(system, position=(0.3, 0.2, 0.0), system_config=cfg)
            sector.place_in_galaxy(tuple(c * constants.PARSEC_M / constants.LIGHTYEAR_M for c in center))
            _SECTOR_CACHE[address] = (sector, center)
    sector, center = _SECTOR_CACHE[address]
    return copy.deepcopy(sector), center


def save_sector(config, address):
    """The generate stage's unit for the sector at `address`: saved, charted, its link left to the link stage; returns its `sectors.id`."""
    sector, center = _build(address)
    return store.save_sector(sector, config=config, link_neighbors=False, galaxy_position={
        "center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
        "galactic_radius_pc": (center[0] ** 2 + center[1] ** 2) ** 0.5,
        "shell_index": address[0], "band_index": address[1], "wedge_index": address[2]})


def _items(config, addresses):
    """`{source address: [[shell, band, wedge, mass], ...]}`: what the stage hands each unit, found as `run_galaxy.backfill_targets` does for one source, from the empty galaxy (the draw checks the rest under lock)."""
    conn = store.get_connection(config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        limits, _seed = run_plan._scatter_limits_and_seed(conn, skeleton)
        bounds = store.get_galaxy_bounds(conn)
        found = {}
        for source in addresses:
            targets = run_galaxy.backfill_step_targets({source}, bounds)
            targets = {address: mass for address, mass in targets.items() if (limit := limits.mass_limit_at(address[0], address[1])) is None or mass < limit}
            found[source] = [[*address, targets[address]] for address in sorted(targets)]
        return found
    finally:
        conn.close()


def _stage_items(config, addresses):
    """The one list plan order draws: `run_galaxy.backfill_targets` around the centers of all the sectors the run saved."""
    conn = store.get_connection(config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        limits, _seed = run_plan._scatter_limits_and_seed(conn, skeleton)
        centers = [spherical.sector_position_pc(*address, toy.EDGE_PC) for address in addresses]
        targets = run_galaxy.backfill_targets(conn, skeleton, limits, centers)
    finally:
        conn.close()
    return [[*address, targets[address]] for address in sorted(targets)]


def _unit(config, items):
    return {"database": wrapper.database_args(config), "items": items}


def _run_unit(config, items, tap):
    with tap.running(BACKFILL):
        json.loads(json.dumps(engine.run(BACKFILL, _unit(config, items), lambda *args, **kwargs: None)))


def _plan_order(config, addresses, tap):
    ids = {address: save_sector(config, address) for address in addresses}
    _run_unit(config, _stage_items(config, addresses), tap)
    return ids


def _actions(addresses, items, scope):
    """`[(name, depends on)]`: the saves in a chain (the same ids in every run), then the unit of each source after its own save and the saves of the sectors of its list that the run saves."""
    actions = {}
    for position, address in enumerate(addresses):
        actions[(SAVE, address)] = {(SAVE, addresses[position - 1])} if position else set()
    for address in addresses:
        listed = {tuple(item[:3]) for item in items[address]}
        wanted = set() if scope == "none" else {address} | (listed & set(addresses))
        actions[(BACKFILL, address)] = {(SAVE, other) for other in wanted}
    return actions


def _released(actions, done):
    return [name for name, deps in actions.items() if name not in done and deps <= done]


def _run_action(config, name, items, tap):
    kind, address = name
    if kind == SAVE:
        save_sector(config, address)
    else:
        _run_unit(config, items[address], tap)


def _by_address_one_thread(config, addresses, items, scope, seed, tap):
    rng = random.Random(seed)
    actions = _actions(addresses, items, scope)
    done = set()
    while len(done) < len(actions):
        choice = rng.choice(_released(actions, done))
        _run_action(config, choice, items, tap)
        done.add(choice)


def _by_address_pool(config, addresses, items, scope, tap):
    actions = _actions(addresses, items, scope)
    done, submitted, lock = set(), set(), threading.Lock()
    with concurrent.futures.ThreadPoolExecutor(WORKERS) as pool:
        pending = set()

        def run(name):
            _run_action(config, name, items, tap)
            return name

        def submit_ready():
            with lock:
                for name in _released(actions, done):
                    if name not in submitted:
                        submitted.add(name)
                        pending.add(pool.submit(run, name))

        submit_ready()
        while pending:
            finished, _ = concurrent.futures.wait(pending, return_when=concurrent.futures.FIRST_COMPLETED)
            for future in finished:
                pending.discard(future)
                with lock:
                    done.add(future.result())
            submit_ready()


def _undeclared(tap):
    slot = tap.touched.get(BACKFILL)
    if slot is None:
        return []
    reads, writes = set(touches.reads_of(BACKFILL)), set(touches.writes_of(BACKFILL))
    lines = []
    extra_reads = slot["reads"] - reads - writes - touches.ENGINE_TABLES if touches.ANY not in reads else set()
    extra_writes = slot["writes"] - writes - touches.ENGINE_TABLES if touches.ANY not in writes else set()
    if extra_reads:
        lines.append(f"{BACKFILL} reads {sorted(extra_reads)}, not declared")
    if extra_writes:
        lines.append(f"{BACKFILL} writes {sorted(extra_writes)}, not declared")
    return lines


def seed_galaxy(config):
    """The toy galaxy the backfill tests use, its mass scatter drawn (the stage before the backfill)."""
    toy._seed_galaxy(config, e_value=toy.SMALL_E_VALUE)
    run_plan.scatter_bright_stars(toy._scatter_args(config, "--mass", "--phenomenon-min-mass", "20"))


def prove(config, scope="declared"):
    """
    The verdict on the backfill's edge (module note). `scope` is `"declared"` (each unit waits for its source and the sectors of its list that the run saves) or `"none"` (released at once: the control, which must fail).
    The database is empty when it is called.
    """
    seed_galaxy(config)
    seeded = support.dump(config)
    addresses = _addresses(config)
    items = _items(config, addresses)
    tap = sqltap.Tap()
    with tap:
        ids = _plan_order(config, addresses, tap)
        reference = _comparable(config, Recipe)
        stars = bool(reference.get("bright_stars"))
        runs = [("plan order again", None)] + [(f"by address, one thread, order {seed}", seed) for seed in SEEDS] + [(f"by address, {WORKERS} threads", "pool"), (f"by address, {WORKERS} threads again", "pool")]
        found = {}
        for name, how in runs:
            support.restore(config, seeded)
            if how is None:
                assert _plan_order(config, addresses, tap) == ids
            elif how == "pool":
                _by_address_pool(config, addresses, items, scope, tap)
            else:
                _by_address_one_thread(config, addresses, items, scope, how, tap)
            found[name] = _comparable(config, Recipe)
    differences = []
    for name, rows in found.items():
        differences += [f"{name}: {line}" for line in support._differences(reference, rows)]
    listed = {tuple(item[:3]) for unit in items.values() for item in unit}
    return Verdict(scope, [name for name, _ in runs], differences, _undeclared(tap),
                   {"stars drawn": stars, "a saved sector is in a list": bool(listed & set(addresses))})
