# tests/address_edge_support.py

"""
The proof of an edge by address (OPS.242, the narrower edges): a command that works on one sector may start as soon as the sectors it reads are saved, while the others are still being saved.

`queue/independence.json` lists whole-command pairs that may overlap (`tests/pair_overlap_support.py`). The plan's narrow edges, "link of sector S waits for S and the sectors around it, settle of S for the
sectors around it", are by address, not by table: the unit of the link or the settle stage that holds sector S starts when the sectors it reads exist, whatever else the generate stage has still to save. This
module proves such an edge on a small galaxy of `SECTORS` sectors in a row, each next to the one before:

* **Plan order** saves every sector (the generate stage's unit: `store.save_sector` with the link left to the link stage), then links every sector, then settles every sector, one unit at a time.
* **By address** runs the same units in an order the edges allow and plan order does not: a link or a settle unit is released the moment the sectors it needs are saved (`needs`, found with the code's own
  neighbor rules: `sector_paths._neighbour_centers` for the settle, `store._sectors_near` at the nearest-systems reach for the link), in several random orders on one thread and on a pool of 15 threads.
  The saves keep their order (the ids the sectors take are the order they landed in, and the comparison needs the same ids); every link and settle unit is a real `engine.run` of the command, as a worker's.
* **The verdict.** Every run is dumped and compared with plan order, table by table, row by row (`cut_neutral_support`; the ids no output depends on are left out, and the paths with their knots, the
  links with their systems are compared through their joins). Each unit's SQL is read by the tap (`queue/sqltap.py`) and must be what `queue/touches.py` declares for the command.
* **A scope that is too small must fail.** `prove(..., scope="self")` releases each unit when its own sector is saved alone; on the toy that gives other bytes than plan order, which shows the toy sees what the
  neighbors are for. The edge a test lists is the one that passes with its declared scope; the control is the same edge with a narrower scope.

What the link unit is given: the run's own sectors, `sector_ids`, not `since` (the sectors saved since a time are the ones saved so far when the unit starts, and a link that is told the run's sectors by the clock
links a sector whole or merges it into the lists of the sectors around, whichever are saved by then); `skip` is that fixed set. The supervisor must give the same.
"""

import concurrent.futures
import copy
import json
import random
import threading
from dataclasses import dataclass, field

from planetgen.db import sector_paths, store
from planetgen.galaxy import spherical
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.queue import commands, engine, sqltap, touches, wrapper
from tests import cut_neutral_recipes as recipes
from tests import cut_neutral_support as support
from tests.fuzz_support import deterministic_entropy

SECTORS = 5
"""int: Sectors in the row. The black hole is in the one before last, near the border of the last one's neighbor."""
EDGE_LY = 11.5
EDGE_PC = EDGE_LY * constants.LIGHTYEAR_M / constants.PARSEC_M
SPACING_PC = 0.9 * EDGE_PC
"""float: Distance between two sectors' centers along the row: under one edge, so every sector has the next within the neighbor reach of both rules."""
SYSTEMS = 2
SEEDS = (11, 12, 13)
"""tuple: The seeds of the random orders of the single-thread runs."""
WORKERS = 15
"""int: The size of the pool of the threaded runs."""
LINK, SETTLE, SAVE = "link-neighbors", "settle-sector-paths", "save"

_SECTOR_CACHE = {}


def _build(index):
    """The `SpaceSector` of position `index` of the row, made once (a pure function of its seed) and copied for each run: two systems near its borders, a black hole at the border of the one before last."""
    if index not in _SECTOR_CACHE:
        with deterministic_entropy(5000 + index):
            base = spherical.sector_position_pc(1500, 2, 700, EDGE_PC)
            center = (base[0] + index * SPACING_PC, base[1], base[2])
            sector = SpaceSector(f"Edge Row {index}", edge_ly=EDGE_LY)
            for number in range(SYSTEMS):
                cfg = SystemConfig()
                cfg.STAR_TYPE = "G2V"
                cfg.BINARY_SYSTEM = False
                system = StarSystem(system_config=cfg)
                system.name = f"Rowstar {index}.{number}"
                side = -1.6 if number == 0 else 1.6
                sector.add_system(system, position=(side, 0.2 * number, 0.0), system_config=cfg)
            if index == SECTORS - 2:
                sector.add_phenomenon(BlackHole(SystemConfig()), "black-hole", position=(1.7, 0.0, 0.0))
            sector.place_in_galaxy(tuple(c * constants.PARSEC_M / constants.LIGHTYEAR_M for c in center))
            _SECTOR_CACHE[index] = (sector, center)
    sector, center = _SECTOR_CACHE[index]
    return copy.deepcopy(sector), center


def save_sector(config, index):
    """The generate stage's unit for the sector at `index`: saved, its link left to the link stage; returns its `sectors.id`."""
    sector, center = _build(index)
    sector_id = store.save_sector(sector, config=config, link_neighbors=False, galaxy_position={
        "center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
        "galactic_radius_pc": (center[0] ** 2 + center[1] ** 2 + center[2] ** 2) ** 0.5,
        "shell_index": 1500, "band_index": 2, "wedge_index": 700 + index})
    if index == SECTORS - 2:
        conn = store.get_connection(config)
        try:
            conn.execute("UPDATE black_holes SET mass_solar = 1.0e7 WHERE sector_id = ?", (sector_id,))
            conn.commit()
        finally:
            conn.close()
    return sector_id


@dataclass
class Plan:
    """The toy's run: `ids` of the sectors in the order saved, `needs[kind][sector id]` the ids the unit of that command for the sector waits for, the recipe's comparison."""

    ids: list
    needs: dict
    recipe: object
    base: dict = field(repr=False, default=None)


def _recipe():
    """The comparison of the link and the settle together: the ignore columns and the views of the two recipes of the cut-neutral check."""
    class Both:
        ignore_columns = recipes.LINK.ignore_columns + recipes.SETTLE.ignore_columns
        views = recipes.SETTLE.views
    return Both


def _needs(config, ids):
    """`{kind: {sector id: set of ids}}` found with the code's own rules on the fully saved galaxy: the settle's neighbor reach and the link's nearest-systems reach."""
    conn = store.get_connection(config)
    try:
        edge = store._edge_pc(conn)
        half_diagonal = edge * (3 ** 0.5) / 2
        found = {LINK: {}, SETTLE: {}}
        for sector_id in ids:
            sector = sector_paths._sector(conn, sector_id)
            found[SETTLE][sector_id] = set(sector_paths._neighbour_centers(conn, sector, sector["edge_mpc"] / 1000.0))
            centers = store._sector_centers(conn, [sector_id])
            found[LINK][sector_id] = set(store._sectors_near(conn, centers, 2 * half_diagonal + store.NEAREST_SYSTEMS_SEARCH_PC)) | {sector_id}
        return found
    finally:
        conn.close()


def _unit(kind, config, sector_id, ids):
    """The params of the unit of `kind` for one sector: what `expand_link_neighbors` and `expand_save_sector_paths` cut a run into, for a run that names its sectors."""
    database = wrapper.database_args(config)
    if kind == LINK:
        return {"database": database, "chunk": [sector_id], "sector_ids": list(ids)}
    return {"database": database, "sector_ids": [sector_id], "workers": 1}


def _run_unit(kind, config, sector_id, ids, tap):
    with tap.running(kind):
        json.loads(json.dumps(engine.run(kind, _unit(kind, config, sector_id, ids), lambda *args, **kwargs: None)))


def _plan_order(config, tap):
    ids = [save_sector(config, index) for index in range(SECTORS)]
    for kind in (LINK, SETTLE):
        for sector_id in ids:
            _run_unit(kind, config, sector_id, ids, tap)
    return ids


def _actions(ids, needs, scope):
    """`[(name, depends on)]`: the saves in a chain (so the ids are the same in every run), then the link and settle unit of each sector, each after the saves of `needs` (its own sector alone for `scope="self"`)."""
    actions = {}
    for position, sector_id in enumerate(ids):
        actions[(SAVE, sector_id)] = {(SAVE, ids[position - 1])} if position else set()
    for kind in (LINK, SETTLE):
        for sector_id in ids:
            wanted = {sector_id} if scope == "self" else needs[kind][sector_id]
            actions[(kind, sector_id)] = {(SAVE, other) for other in wanted}
    return actions


def _released(actions, done):
    return [name for name, deps in actions.items() if name not in done and deps <= done]


def _by_address_one_thread(config, ids, needs, scope, seed, tap):
    """The units one at a time, in a random order among those released: a unit is chosen at random whenever several are ready."""
    rng = random.Random(seed)
    actions = _actions(ids, needs, scope)
    index_of = {sector_id: position for position, sector_id in enumerate(ids)}
    done = set()
    while len(done) < len(actions):
        ready = _released(actions, done)
        choice = rng.choice(ready)
        kind, sector_id = choice
        if kind == SAVE:
            assert save_sector(config, index_of[sector_id]) == sector_id, "the sectors must take the same ids in every run"
        else:
            _run_unit(kind, config, sector_id, ids, tap)
        done.add(choice)


def _by_address_pool(config, ids, needs, scope, tap):
    """The same on a pool of `WORKERS` threads: a unit is submitted when its needs are done, and the saves wait for each other (the ids)."""
    actions = _actions(ids, needs, scope)
    index_of = {sector_id: position for position, sector_id in enumerate(ids)}
    done, submitted, lock = set(), set(), threading.Lock()
    with concurrent.futures.ThreadPoolExecutor(WORKERS) as pool:
        pending = set()

        def run(name):
            kind, sector_id = name
            if kind == SAVE:
                assert save_sector(config, index_of[sector_id]) == sector_id, "the sectors must take the same ids in every run"
            else:
                _run_unit(kind, config, sector_id, ids, tap)
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


@dataclass
class Verdict:
    """What `prove` found: the runs, the rows that differ, the commands that touch what they do not declare."""

    scope: str
    runs: list
    differences: list
    undeclared: list
    stored: dict

    @property
    def passed(self):
        return not self.differences and not self.undeclared and all(self.stored.values())

    def line(self):
        words = "EDGES HOLD" if self.passed else "FAILS"
        return f"{words} (scope {self.scope}): {len(self.runs)} runs compared with plan order, differences {len(self.differences)}, undeclared {len(self.undeclared)}, stored {self.stored}"


def _comparable(config, recipe):
    class Holder:
        views = recipe.views
    return {**support._comparable(support.dump(config), set(recipe.ignore_columns)), **support._view_rows(config, Holder)}


def _undeclared(tap):
    lines = []
    for kind in (LINK, SETTLE):
        slot = tap.touched.get(kind)
        if slot is None:
            continue
        reads, writes = set(touches.reads_of(kind)), set(touches.writes_of(kind))
        extra_reads = slot["reads"] - reads - writes - touches.ENGINE_TABLES if touches.ANY not in reads else set()
        extra_writes = slot["writes"] - writes - touches.ENGINE_TABLES if touches.ANY not in writes else set()
        if extra_reads:
            lines.append(f"{kind} reads {sorted(extra_reads)}, not declared")
        if extra_writes:
            lines.append(f"{kind} writes {sorted(extra_writes)}, not declared")
    return lines


def prove(config, scope="declared"):
    """
    The verdict on the edges (module note). `scope` is `"declared"` (each unit waits for the sectors the code's own rule says it reads) or `"self"` (only its own sector: the control, which must fail).
    The database is empty when it is called and is put back to that state before each run.
    """
    recipe = _recipe()
    empty = support.dump(config)
    tap = sqltap.Tap()
    with tap:
        support.restore(config, empty)
        ids = _plan_order(config, tap)
        reference = _comparable(config, recipe)
        needs = _needs(config, ids)
        stored = {LINK: bool(reference.get("nearest_systems")), SETTLE: bool(reference.get("sector_paths"))}
        runs = [("plan order again", None)] + [(f"by address, one thread, order {seed}", seed) for seed in SEEDS] + [(f"by address, {WORKERS} threads", "pool"), (f"by address, {WORKERS} threads again", "pool")]
        found = {}
        for name, how in runs:
            support.restore(config, empty)
            if how is None:
                assert _plan_order(config, tap) == ids
            elif how == "pool":
                _by_address_pool(config, ids, needs, scope, tap)
            else:
                _by_address_one_thread(config, ids, needs, scope, how, tap)
            found[name] = _comparable(config, recipe)
    differences = []
    for name, rows in found.items():
        differences += [f"{name}: {line}" for line in support._differences(reference, rows)]
    return Verdict(scope, [name for name, _ in runs], differences, _undeclared(tap), stored)
