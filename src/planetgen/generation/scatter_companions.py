# planetgen/generation/scatter_companions.py

"""
Companions in the scatter passes (GEN.262)
==========================================

The star passes (the mass pass, the luminosity pass and the backfills that go with them) and the phenomenon pass
write one row per placed object. This module gives the rows that have company a member tree
(`multistar.draw_member_tree`, docs/design/systems/multi-star-systems.md section 6) and stores it with the row
(`db.scatter_systems`): the row's point becomes the system's barycenter, and its `member_count` and
`system_mass_solar` say so.

* **A tree is a function of the row's seed.** It is drawn from `draw.Stream(f"{seed}:members")` in a fixed order
  (the system's age first, then the tree), so two scatters with the same galaxy seed store the same members, and a
  later stage that draws from another suffix (`planets`) never moves them (6.2). The rows themselves are untouched:
  the passes draw exactly the rows they drew before.
* **Bulk screening (6.3).** A black hole, a neutron star or a planetary nebula's central star has a companion at the
  odds of `multistar.MEMBER_STATE_BY_PLACED` (1.25% for a stellar black hole). Whether a row has one is one uniform
  derived from its seed by an integer hash (`screen_unit`): a row it does not name costs a few operations and no
  generator call, and only the hits (7.4e5 of 5.9e7 black holes) draw a tree. A star row always draws its tree: its
  number of stars is the table of 3.3 (most O and B stars are not alone), so a screen would pass most of them.
* **A row is a system (6.4).** The phenomenon scatter divides a class's object rate by `1 + p_placed`
  (`phenomenon_scatter.row_density_divisor`), so the objects of a class, rows and same-class companions together,
  keep the class's density.
* **Ids.** A bulk insert does not return the ids of its rows, so each batch is written first and its hits are then
  found by seed among the ids above the table's top at the start of the batch (`_ids_by_seed`), however many other
  workers write to the table meanwhile.

`insert_bright_stars` and `insert_phenomenon_scatter` are the writers the scatter, the backfills and the plan call in
place of the plain ones in `db.store`.
"""

import collections
import dataclasses
import functools
import itertools
import math

from planetgen import tuning
from planetgen.db import contention, scatter_systems, store
from planetgen.generation import cluster_compact_binaries, multistar, phenomenon_scatter
from planetgen.physics import imf, remnant_mass
from planetgen.physics import stellar_evolution as evolution, unit_registry as units
from planetgen.util import draw

_kg_to_msun = units.converter("kg", "Msun")

MEMBERS_STREAM = "members"
"""str: The suffix of a row's seed that names the stream its tree is drawn from (and, later, `planets`)."""

SCREEN_SALTS = {"companion": 0x9E3779B97F4A7C15, "recycled": 0xD1B54A32D192ED03,
                "runaway": 0xA24BAED4963EE407, "runaway-origin": 0x9FB21C651E98DF25, "runaway-bound": 0xC13FA9A902A6328F}
"""dict: What `screen_unit` mixes into the seed for each question it answers, so the questions are independent."""

_MASK = (1 << 64) - 1
_STAR_COLUMN = {name: i for i, name in enumerate(store.BRIGHT_STAR_COLUMNS)}
_PHENOMENON_COLUMN = {name: i for i, name in enumerate(store.PHENOMENON_SCATTER_COLUMNS)}
_LOOKUP_SEEDS = 1000

Written = collections.namedtuple("Written", "rows systems")
"""namedtuple: What an insert reports: the rows written and how many of them are systems with members."""


def screen_unit(seed, question="companion"):
    """
    A uniform in [0, 1) that is a pure function of a row's `seed` and the `question` (the SplitMix64 finalizer, 64-bit
    integer arithmetic only, so it is the same on every machine and version and costs about a microsecond).
    """
    z = (seed + SCREEN_SALTS[question]) & _MASK
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _MASK
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _MASK
    z ^= z >> 31
    return (z >> 11) / 9007199254740992.0


@functools.lru_cache(maxsize=None)
def _odds(placed_kind, recycled, planetary_nebula, cluster_key=None):
    """`(total, ((state, cumulative upper bound), ...))` of a placed class's companion table; `cluster_key` is `(mass, half-mass radius)` of the cluster whose table it is (GEN.330), `None` for the field's."""
    probe = multistar.Placed(placed_kind, 1.0, recycled=recycled, planetary_nebula=planetary_nebula)
    tables = multistar.MEMBER_STATE_BY_PLACED if cluster_key is None else cluster_compact_binaries.tables_of(*cluster_key)
    table = tables.get(multistar.placed_class(probe), {})
    steps, running = [], 0.0
    for state, (default, _low, _high) in table.items():
        running += default
        steps.append((state, running))
    return running, tuple(steps)


def _cluster_key(cluster):
    """The hashable key of a cluster for `_odds`: its mass and half-mass radius, the two numbers its companion odds depend on; `None` for no cluster."""
    return None if cluster is None else (float(cluster.mass_solar), float(cluster.half_mass_radius_pc))


def phenomenon_companion(kind, subtype, seed, cluster=None):
    """
    The screen of one phenomenon row: `(placed_kind, recycled, state)` when the row has a companion, else `None`.
    `placed_kind` is the `multistar.Placed` kind, `recycled` whether a neutron star is a millisecond pulsar, `state`
    the companion state (`star`, `black-hole`, ...). A class without a tree (`phenomenon_scatter.PLACED_KINDS`) is
    never a hit. The odds are `multistar.MEMBER_STATE_BY_PLACED`'s, walked with the row's `screen_unit`; a row of a globular cluster (`cluster`, its `galaxy.components.Component`) takes the cluster's
    (`cluster_compact_binaries.tables`, GEN.330).
    """
    if phenomenon_scatter.is_pair(kind, subtype):
        return _pair_companion(kind, subtype, seed, cluster)
    placed_kind = phenomenon_scatter.PLACED_KINDS.get((kind, subtype))
    if placed_kind is None:
        return None
    recycled = (placed_kind == "neutron-star"
                and screen_unit(seed, "recycled") < phenomenon_scatter.NEUTRON_STAR_RECYCLED_SHARE)
    total, steps = _odds(placed_kind, recycled, kind == "planetary-nebula", _cluster_key(cluster))
    pick = screen_unit(seed)
    if pick >= total:
        return None
    for state, upper in steps:
        if pick < upper:
            return placed_kind, recycled, state
    return None


def _pair_companion(kind, subtype, seed, cluster=None):
    """
    The companion of a promoted compact pair row (GEN.263): always a hit, whose state is picked among the compact companions of the placed object
    (`phenomenon_scatter.PAIR_COMPANIONS`) by the odds the class table gives each, with the row's own `screen_unit`. A neutron star of a pair row is a millisecond pulsar
    at the share the pair rows have (`phenomenon_scatter.pair_recycled_share`), not the share the neutron stars have.
    """
    placed_kind = phenomenon_scatter.PLACED_KINDS[(kind, phenomenon_scatter.base_subtype(kind, subtype))]
    recycled = (placed_kind == "neutron-star"
                and screen_unit(seed, "recycled") < phenomenon_scatter.pair_recycled_share(kind, subtype, cluster))
    odds = phenomenon_scatter._pair_table(placed_kind, recycled, cluster)
    pick = screen_unit(seed) * sum(odds.values())
    running = 0.0
    for state, share in odds.items():
        running += share
        if pick < running:
            return placed_kind, recycled, state
    return placed_kind, recycled, next(reversed(odds))


def _planetary_nebula_star(rng):
    """`(Placed, age_yr)` of a planetary nebula's central white dwarf: a progenitor from the white-dwarf range, the
    initial-to-final mass relation, and an age of the progenitor's main-sequence life (the nebula is the end of its
    giant phase)."""
    low, high = tuning.WHITE_DWARF_PROGENITOR_SOL
    initial = imf.sample(low, high, rng)
    mass = remnant_mass.white_dwarf_mass_sol(initial)
    age_yr = evolution.main_sequence_lifetime_gy(initial) * 1.0e9
    return multistar.Placed("white-dwarf", mass, initial, planetary_nebula=True), age_yr


def phenomenon_tree(row, cluster=None):
    """
    The member tree of one `phenomenon_scatter` row (`store.PHENOMENON_SCATTER_COLUMNS` order), or `None` when the
    row has no companion or no room for one. A compact object's system age is drawn from the compact remnant age
    range (`tuning.COMPACT_REMNANT_AGE_RANGE_GY`, the age the fill gives such a remnant); a planetary nebula's from its
    progenitor's life.

    Args:
        row (tuple): The row.
        cluster (galaxy.components.Component, optional): The globular cluster the row belongs to (its `cluster_id`), when it has one (GEN.330): the companion is drawn at the cluster's odds
            (`cluster_compact_binaries.tables`), the system is as old as the cluster, and a pair wider than the cluster keeps (`cluster_compact_binaries.axis_limit_au`) is not drawn. `None`: the field's.
    """
    kind, subtype = row[_PHENOMENON_COLUMN["kind"]], row[_PHENOMENON_COLUMN["subtype"]]
    seed = row[_PHENOMENON_COLUMN["seed"]]
    hit = phenomenon_companion(kind, subtype, seed, cluster)
    if hit is None:
        return None
    placed_kind, recycled, state = hit
    rng = draw.Stream(f"{seed}:{MEMBERS_STREAM}")
    if kind == "planetary-nebula":
        placed, age_yr = _planetary_nebula_star(rng)
    else:
        mass = row[_PHENOMENON_COLUMN["mass_solar"]]
        if mass is None:
            return None
        placed = multistar.Placed(placed_kind, float(mass), recycled=recycled)
        age_yr = rng.uniform(*tuning.COMPACT_REMNANT_AGE_RANGE_GY) * 1.0e9
    max_axis_au = math.inf
    if cluster is not None and kind != "planetary-nebula":
        age_yr = float(cluster.age_gy) * 1.0e9              # the cluster's stars were born together
        max_axis_au = cluster_compact_binaries.axis_limit_au(cluster, placed.mass_solar)
    tree = multistar.draw_member_tree(placed, age_yr, rng, companion=state, max_axis_au=max_axis_au)
    return tree if tree.nodes else None


def runaway_draw(seed, initial_mass_sol):
    """
    Whether the star of a row with this `seed` and birth mass is a runaway, and by which channel
    (`multistar.draw_runaway_origin`, GEN.329): `None` or a `multistar.RunawayDraw`. A pure function of the seed, like the
    companion screen, so the scatter that stores the star's tree and the fill that builds the star agree without either
    storing the answer for a single star.
    """
    return multistar.draw_runaway_origin(initial_mass_sol, screen_unit(seed, "runaway"),
                                         screen_unit(seed, "runaway-origin"), screen_unit(seed, "runaway-bound"))


def star_tree(row):
    """
    The member tree of one `bright_stars` row (`store.BRIGHT_STAR_COLUMNS` order), or `None` for a single object.
    The row is the placed object: a star of its stored mass, birth mass and age (a white dwarf, class VII, a white
    dwarf), and the system is as old as the star. A runaway star (`runaway_draw`) carries the flag of its channel on
    member 0, and a supernova-channel runaway that kept a compact companion has it in its tree. A massive star under
    `tuning.YOUNG_GROUP_AGE_MAX_YR` may be born in a young group (`multistar.draw_young_group`, GEN.327): a `group` node
    while it is younger than the group's decay age, the hierarchy the group left otherwise.
    """
    seed = row[_STAR_COLUMN["seed"]]
    white_dwarf = row[_STAR_COLUMN["yerkes_class"]] == "VII"
    placed = multistar.Placed("white-dwarf" if white_dwarf else "star",
                              _kg_to_msun(row[_STAR_COLUMN["mass_kg"]]),
                              row[_STAR_COLUMN["initial_mass_sol"]], fe_h=row[_STAR_COLUMN["metallicity_fe_h"]])
    age_yr = max(row[_STAR_COLUMN["age_gy"]], 1.0e-6) * 1.0e9
    runaway = None if white_dwarf else runaway_draw(seed, placed.birth_mass_solar)
    tree = multistar.draw_member_tree(placed, age_yr, draw.Stream(f"{seed}:{MEMBERS_STREAM}"),
                                      bound_runaway=runaway is not None and runaway.bound_companion, groups=True)
    if not tree.nodes:
        return None
    if runaway is not None:
        flagged = dataclasses.replace(tree.members[0], flags=tree.members[0].flags | multistar.RUNAWAY_ORIGIN_FLAGS[runaway.origin])
        tree = dataclasses.replace(tree, members=(flagged, *tree.members[1:]))
    return tree


def star_trees(rows, conn=None):
    """`[(index, tree), ...]` for the rows of `rows` (`store.BRIGHT_STAR_COLUMNS` order) that have members. `conn` is read by `phenomenon_trees` and not here (the writers pass it to both)."""
    return [(i, tree) for i, row in enumerate(rows) if (tree := star_tree(row)) is not None]


def phenomenon_trees(rows, conn=None):
    """
    `[(index, tree), ...]` for the rows of `rows` (`store.PHENOMENON_SCATTER_COLUMNS` order) that have members. With a connection, a row of a globular cluster is drawn at its cluster's odds
    (`phenomenon_tree`, GEN.330): the clusters of the batch are read once, by id; without one (or a cluster the table does not hold) every row is the field's.
    """
    clusters = _clusters_of(conn, rows)
    cluster_at = _PHENOMENON_COLUMN["cluster_id"]
    return [(i, tree) for i, row in enumerate(rows) if (tree := phenomenon_tree(row, clusters.get(row[cluster_at]))) is not None]


def _clusters_of(conn, rows):
    """`{id: Component}` of the globular clusters the `rows` (`store.PHENOMENON_SCATTER_COLUMNS` order) name in their `cluster_id`; empty without a connection or a cluster."""
    if conn is None:
        return {}
    ids = {row[_PHENOMENON_COLUMN["cluster_id"]] for row in rows} - {None}
    if not ids:
        return {}
    return {each.id: each for each in store.get_galaxy_components(conn, kind="cluster", ids=ids)}


def default_epoch_unix(conn):
    """The time to stamp new nodes with: the orbit epoch every stored position holds at, else the fixed
    `tuning.SCATTER_PHASE_ORIGIN_UNIX` (no orbit update has run, so no clock exists yet; a fixed origin keeps two scatters
    of one seed identical, the database's clock would not)."""
    epoch = store.get_orbit_epoch_unix(conn)
    return float(epoch) if epoch is not None else tuning.SCATTER_PHASE_ORIGIN_UNIX


def _ids_by_seed(conn, table, above_id, seeds):
    """`{seed: [id, ...]}` of the rows of `table` above `above_id` with one of `seeds`: the ids the batch just
    written got. Other workers' rows in the same range are skipped by their seeds."""
    found = collections.defaultdict(list)
    seeds = list(seeds)
    for start in range(0, len(seeds), _LOOKUP_SEEDS):
        chunk = seeds[start:start + _LOOKUP_SEEDS]
        for row in conn.execute(
                f"SELECT id, seed FROM {table} WHERE id > ? AND seed IN ({', '.join('?' * len(chunk))}) ORDER BY id",
                (above_id, *chunk)).fetchall():
            found[row["seed"]].append(row["id"])
    return found


def _write(conn, table, rows, trees_of, seed_index, insert, epoch_unix, batch_size):
    written = systems = 0
    rows = iter(rows)
    while True:
        batch = list(itertools.islice(rows, batch_size))
        if not batch:
            return Written(written, systems)
        hits = trees_of(batch, conn)
        top = conn.execute(f"SELECT COALESCE(MAX(id), 0) AS top FROM {table}").fetchone()["top"] if hits else 0
        written += insert(conn, batch)
        if not hits:
            continue
        ids = _ids_by_seed(conn, table, top, {batch[i][seed_index] for i, _tree in hits})
        taken = collections.Counter()
        items = []
        for index, tree in hits:
            seed = batch[index][seed_index]
            candidates = ids.get(seed, ())
            if taken[seed] < len(candidates):
                items.append((candidates[taken[seed]], tree))
                taken[seed] += 1
        systems += scatter_systems.store_member_trees(conn, table, items, epoch_unix)


def insert_bright_stars(conn, rows, epoch_unix=None, batch_size=10000):
    """
    `store.insert_bright_stars` that also draws and stores the members of the rows that have them
    (`star_tree`). The caller commits.

    Args:
        conn (Connection): An open, read-write connection.
        rows (iterable): `store.BRIGHT_STAR_COLUMNS` rows.
        epoch_unix (float, optional): The time the nodes' phases hold at (`default_epoch_unix` when
            not given).
        batch_size (int): Rows to a write.

    Returns:
        Written: The rows written and the systems among them.
    """
    return _insert(conn, "bright_stars", rows, star_trees, _STAR_COLUMN["seed"], store.insert_bright_stars, epoch_unix, batch_size)


def _insert(conn, table, rows, trees_of, seed_index, plain, epoch_unix, batch_size):
    """`insert_bright_stars` and `insert_phenomenon_scatter`: the epoch when none is given, then `_write` with `plain` (the
    `db.store` writer) as the insert of each batch."""
    if epoch_unix is None:
        epoch_unix = default_epoch_unix(conn)
    return _write(conn, table, rows, trees_of, seed_index, lambda c, batch: plain(c, batch, batch_size), epoch_unix, batch_size)


def write_committed(conn, write, rows, epoch_unix, label="Scatter write"):
    """
    Writes one batch of scatter rows and commits it; a deadlock or lock wait timeout rolls it back and the engine runs the unit again (the scatter writers' ID
    step, `scatter_uids.assign`, reads every row of the table that has no ID yet, so workers writing at once meet on each other's rows). The
    unit's retry is the engine's (`queue/failures.py`, OPS.220; the unit's own tries, waits and exclusive grant went with OPS.261, `db.contention.run_once`).
    The batch is READ COMMITTED, so its reads take no gap locks. A batch is the whole transaction, so a rollback takes back exactly it and the
    retry writes the same rows with the same IDs: they rank by luminosity and seed, not by the order the rows landed.

    Args:
        conn (Connection): Open, with no uncommitted work of the caller's (it is committed first).
        write (callable): `insert_bright_stars` or `insert_phenomenon_scatter`.
        rows (list): The batch, kept whole for a retry.
        epoch_unix (float): As the writers take it.
        label (str): What the log calls the batch.

    Returns:
        Written: What `write` reported for the batch that landed.

    Raises:
        pymysql.err.OperationalError: Any other error, or a deadlock that persists through the granted try.
    """
    def attempt(c):
        c.commit()          # a connection that just read has a transaction open, and MySQL refuses the next line inside one
        c.execute("SET TRANSACTION ISOLATION LEVEL READ COMMITTED")
        written = write(c, rows, epoch_unix)
        c.commit()
        return written

    return contention.run_once(conn, attempt)


_ONCE = {"insert_bright_stars": ("bright_stars", star_trees, _STAR_COLUMN["seed"]),
         "insert_phenomenon_scatter": ("phenomenon_scatter", phenomenon_trees, _PHENOMENON_COLUMN["seed"])}
"""dict: What `write_once` needs of each writer, by its function name: the table, the rows' member trees and the seed's place in a row."""


def _landed(conn, table, row, seed_index):
    """Whether a row of `table` at the row's cell holds the row's seed (the seed is a draw of 63 bits, the cell is the table's address index)."""
    return conn.execute(f"SELECT 1 AS found FROM {table} WHERE shell_index = ? AND band_index = ? AND wedge_index = ? AND seed = ? LIMIT 1",
                        (*row[:3], row[seed_index])).fetchone() is not None


def write_once(conn, write, rows, epoch_unix, label="Scatter write"):
    """
    `write_committed` for a piece that may run again (OPS.264): a batch that an earlier try of the piece already wrote is not written a second time. A piece commits its batches one by one and the
    engine marks the piece done after the last, so a unit run again after a stop or a lost connection draws the same rows (each piece has its own random stream) and would write every batch it
    reached again, under new Object IDs, with nothing to refuse them. A batch is one transaction, so its first row and its last are both in the table or neither is; the check is two reads
    by the cell's index.

    Returns:
        Written: What the batch wrote; for a batch found written, its rows and the systems among them (`star_trees` or `phenomenon_trees`, drawn from the seeds as the write drew them).

    Raises:
        RuntimeError: One of the two rows is in the table and the other is not (the batch was not written whole, or the seed repeats).
    """
    once = _ONCE.get(getattr(write, "__name__", ""))
    if once is None or not rows:
        return write_committed(conn, write, rows, epoch_unix, label)
    table, trees_of, seed_index = once
    conn.commit()           # a connection that just read has a transaction open: the check must see what the earlier try committed
    found = [_landed(conn, table, row, seed_index) for row in (rows[0], rows[-1])]
    if all(found):
        return Written(len(rows), len(trees_of(rows, conn)))
    if any(found):
        raise RuntimeError(f"{label}: a batch of {len(rows):,} rows is in {table} in part (its first row is {'there' if found[0] else 'not there'} and its last {'there' if found[1] else 'not there'})")
    return write_committed(conn, write, rows, epoch_unix, label)


def insert_phenomenon_scatter(conn, rows, epoch_unix=None, batch_size=5000):
    """
    `store.insert_phenomenon_scatter` that also draws and stores the members of the rows that have them
    (`phenomenon_tree`). The caller commits. Arguments and result as `insert_bright_stars`.
    """
    return _insert(conn, "phenomenon_scatter", rows, phenomenon_trees, _PHENOMENON_COLUMN["seed"], store.insert_phenomenon_scatter,
                   epoch_unix, batch_size)
