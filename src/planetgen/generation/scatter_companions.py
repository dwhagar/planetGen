# planetgen/generation/scatter_companions.py

"""
Companions in the scatter passes (GEN.262)
==========================================

The star passes (the mass pass, the luminosity pass and the backfills that go with them) and the phenomenon pass
write one row per placed object. This module gives the rows that have company a member tree
(`multistar.draw_member_tree`, docs/design/systems/multi-star-systems.md section 6) and stores it with the row
(`db.scatter_systems`): the row's point becomes the system's barycentre, and its `member_count` and
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

from planetgen import tuning
from planetgen.db import scatter_systems, store
from planetgen.generation import multistar, phenomenon_scatter
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
def _odds(placed_kind, recycled, planetary_nebula):
    """`(total, ((state, cumulative upper bound), ...))` of a placed class's companion table."""
    probe = multistar.Placed(placed_kind, 1.0, recycled=recycled, planetary_nebula=planetary_nebula)
    table = multistar.MEMBER_STATE_BY_PLACED.get(multistar.placed_class(probe), {})
    steps, running = [], 0.0
    for state, (default, _low, _high) in table.items():
        running += default
        steps.append((state, running))
    return running, tuple(steps)


def phenomenon_companion(kind, subtype, seed):
    """
    The screen of one phenomenon row: `(placed_kind, recycled, state)` when the row has a companion, else `None`.
    `placed_kind` is the `multistar.Placed` kind, `recycled` whether a neutron star is a millisecond pulsar, `state`
    the companion state (`star`, `black-hole`, ...). A class without a tree (`phenomenon_scatter.PLACED_KINDS`) is
    never a hit. The odds are `multistar.MEMBER_STATE_BY_PLACED`'s, walked with the row's `screen_unit`.
    """
    placed_kind = phenomenon_scatter.PLACED_KINDS.get((kind, subtype))
    if placed_kind is None:
        return None
    recycled = (placed_kind == "neutron-star"
                and screen_unit(seed, "recycled") < phenomenon_scatter.NEUTRON_STAR_RECYCLED_SHARE)
    total, steps = _odds(placed_kind, recycled, kind == "planetary-nebula")
    pick = screen_unit(seed)
    if pick >= total:
        return None
    for state, upper in steps:
        if pick < upper:
            return placed_kind, recycled, state
    return None


def _planetary_nebula_star(rng):
    """`(Placed, age_yr)` of a planetary nebula's central white dwarf: a progenitor from the white-dwarf range, the
    initial-to-final mass relation, and an age of the progenitor's main-sequence life (the nebula is the end of its
    giant phase)."""
    low, high = tuning.WHITE_DWARF_PROGENITOR_SOL
    initial = imf.sample(low, high, rng)
    mass = remnant_mass.white_dwarf_mass_sol(initial)
    age_yr = evolution.main_sequence_lifetime_gy(initial) * 1.0e9
    return multistar.Placed("white-dwarf", mass, initial, planetary_nebula=True), age_yr


def phenomenon_tree(row):
    """
    The member tree of one `phenomenon_scatter` row (`store.PHENOMENON_SCATTER_COLUMNS` order), or `None` when the
    row has no companion or no room for one. A compact object's system age is drawn from the compact remnant age
    range (`tuning.COMPACT_REMNANT_AGE_RANGE_GY`, the age the fill gives such a remnant); a planetary nebula's from its
    progenitor's life.
    """
    kind, subtype = row[_PHENOMENON_COLUMN["kind"]], row[_PHENOMENON_COLUMN["subtype"]]
    seed = row[_PHENOMENON_COLUMN["seed"]]
    hit = phenomenon_companion(kind, subtype, seed)
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
    tree = multistar.draw_member_tree(placed, age_yr, rng, companion=state)
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


def star_trees(rows):
    """`[(index, tree), ...]` for the rows of `rows` (`store.BRIGHT_STAR_COLUMNS` order) that have members."""
    return [(i, tree) for i, row in enumerate(rows) if (tree := star_tree(row)) is not None]


def phenomenon_trees(rows):
    """`[(index, tree), ...]` for the rows of `rows` (`store.PHENOMENON_SCATTER_COLUMNS` order) that have members."""
    return [(i, tree) for i, row in enumerate(rows) if (tree := phenomenon_tree(row)) is not None]


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
        hits = trees_of(batch)
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
    if epoch_unix is None:
        epoch_unix = default_epoch_unix(conn)
    return _write(conn, "bright_stars", rows, star_trees, _STAR_COLUMN["seed"],
                  lambda c, batch: store.insert_bright_stars(c, batch, batch_size), epoch_unix, batch_size)


def insert_phenomenon_scatter(conn, rows, epoch_unix=None, batch_size=5000):
    """
    `store.insert_phenomenon_scatter` that also draws and stores the members of the rows that have them
    (`phenomenon_tree`). The caller commits. Arguments and result as `insert_bright_stars`.
    """
    if epoch_unix is None:
        epoch_unix = default_epoch_unix(conn)
    return _write(conn, "phenomenon_scatter", rows, phenomenon_trees, _PHENOMENON_COLUMN["seed"],
                  lambda c, batch: store.insert_phenomenon_scatter(c, batch, batch_size), epoch_unix, batch_size)
