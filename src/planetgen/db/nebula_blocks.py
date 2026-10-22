# planetgen/db/nebula_blocks.py

"""
Stored nebula aggregates for the Galaxy Map's block ladder (MAP.172, `nebula_block_stats` and
`nebula_block_classes`; the rules are `galaxy/nebula_blocks.py`).

For every block of the five levels above a sector (`drill.DRILL_LEVELS[:-1]`) the table holds what the nebulae add
up to: count, volume, centroid, bounding radius, the volume of each family, the dominant class and the cover. The
nebulae are the rows of `nebula_field` (the molecular clouds the seed places, classes M and N) and the stored
`nebulae` that are not one of them (matched by `uid`: emission and reflection nebulae around a star, planetary
nebulae, remnants, a class P or Q cloud). No sector has to exist, so a block no sector was generated in has its cover
too. The table is:

* `rebuild`: made again from every nebula in one pass (`planetgen scatter --nebula-field`'s step, the migration, `fill-nebula-blocks`);
* `catch_up`: the nebulae stored since (not one at a time as they are stored, which would make every worker of a run
  wait on the same top-level blocks) are added to the blocks they are in and reach, once at the end of a `galaxy` or
  `sector` run.
"""


import numpy as np

from planetgen.galaxy import nebula_blocks as rules
from planetgen.physics.units import ly_to_pc
from planetgen.util import detmath

LEVELS = rules.LEVELS

FAMILY_COLUMNS = ("volume_diffuse_pc3", "volume_emission_pc3", "volume_reflection_pc3", "volume_planetary_pc3",
                  "volume_dark_pc3", "volume_remnant_pc3")
"""tuple: The table's family volume columns, in the order of `rules.FAMILIES`."""

BLOCK_COLUMNS = ("m", "shell_index", "band_index", "wedge_index", "nebula_count", "volume_pc3", "centroid_x_pc",
                 "centroid_y_pc", "centroid_z_pc", "bound_radius_pc", *FAMILY_COLUMNS, "dominant_class",
                 "cover_volume_pc3", "cover")
"""tuple: `nebula_block_stats`' columns in the order `block_params` writes them."""

CLASS_COLUMNS = ("m", "shell_index", "band_index", "wedge_index", "nebula_class", "volume_pc3")

INSERT_BLOCK = (f"INSERT INTO nebula_block_stats ({', '.join(BLOCK_COLUMNS)}) VALUES ({', '.join('?' for _ in BLOCK_COLUMNS)})")
INSERT_CLASS = (f"INSERT INTO nebula_block_classes ({', '.join(CLASS_COLUMNS)}) VALUES ({', '.join('?' for _ in CLASS_COLUMNS)})")

FIELD_SLAB_CELLS = 10
"""int: Cells along x of the nebula field read at a time."""

WRITE_BATCH = 2000


def block_params(row):
    """The `INSERT_BLOCK` values of one block dict (`rules.level_rows`)."""
    return (row["m"], row["shell_index"], row["band_index"], row["wedge_index"], row["nebula_count"],
            row["volume_pc3"], row["centroid_x_pc"], row["centroid_y_pc"], row["centroid_z_pc"],
            row["bound_radius_pc"], *row["families"], row["dominant_class"], row["cover_volume_pc3"], row["cover"])


def class_params(m, row):
    """The `INSERT_CLASS` values of one `rules.LevelRows.classes` entry of level `m`."""
    shell, band, wedge, letter, volume = row
    return (m, shell, band, wedge, letter, volume)


def _arrays(parts):
    if not parts:
        return tuple(np.zeros(0) for _ in range(5))
    letters = [letter for part in parts for letter in part[0]]
    columns = [np.concatenate([part[i] for part in parts]) for i in range(1, 5)]
    return (rules.class_codes(letters), *columns)


def _part(rows):
    """`(class letters, x, y, z, radius_pc)` of rows `(class, radius_ly, x, y, z)`."""
    return ([r[0] for r in rows], np.array([r[2] for r in rows], dtype=float), np.array([r[3] for r in rows], dtype=float),
            np.array([r[4] for r in rows], dtype=float), np.array([ly_to_pc(r[1]) for r in rows], dtype=float))


def sources(conn):
    """
    Every nebula the aggregates count, as `(class codes, x, y, z, radius in parsecs)` arrays: the rows of
    `nebula_field` and the stored placed `nebulae` that no `nebula_field` row holds.
    """
    parts = []
    span = conn.execute("SELECT MIN(cell_i) AS low, MAX(cell_i) AS high FROM nebula_field").fetchone()
    if span and span["low"] is not None:
        for low in range(int(span["low"]), int(span["high"]) + 1, FIELD_SLAB_CELLS):
            rows = conn.execute(
                "SELECT nebula_class, radius_ly, center_x_pc, center_y_pc, center_z_pc FROM nebula_field"
                " WHERE cell_i BETWEEN ? AND ?", (low, low + FIELD_SLAB_CELLS - 1)).fetchall()
            if rows:
                parts.append(_part([(r["nebula_class"], r["radius_ly"], r["center_x_pc"], r["center_y_pc"],
                                     r["center_z_pc"]) for r in rows]))
    rows = conn.execute(
        "SELECT n.nebula_class, n.radius_ly, n.center_x_pc, n.center_y_pc, n.center_z_pc FROM nebulae n"
        " WHERE n.center_x_pc IS NOT NULL AND NOT EXISTS (SELECT 1 FROM nebula_field f WHERE f.uid = n.uid)").fetchall()
    if rows:
        parts.append(_part([(r["nebula_class"], r["radius_ly"], r["center_x_pc"], r["center_y_pc"],
                             r["center_z_pc"]) for r in rows]))
    return _arrays(parts)


def clear(conn):
    """Empties both tables and the mark."""
    conn.execute("DELETE FROM nebula_block_stats")
    conn.execute("DELETE FROM nebula_block_classes")
    conn.execute("DELETE FROM nebula_block_state")


def write(conn, rows_by_level):
    """Writes `rules.all_rows`' result into the tables, a batch of rows at a statement. Returns the block rows."""
    written = 0
    for m, level in rows_by_level.items():
        blocks = [block_params(row) for row in level.blocks]
        for start in range(0, len(blocks), WRITE_BATCH):
            conn.executemany(INSERT_BLOCK, blocks[start:start + WRITE_BATCH])
        classes = [class_params(m, row) for row in level.classes]
        for start in range(0, len(classes), WRITE_BATCH):
            conn.executemany(INSERT_CLASS, classes[start:start + WRITE_BATCH])
        written += len(blocks)
    return written


def rebuild(conn, edge_pc, on_progress=None):
    """
    Computes every row from `nebula_field` and the stored nebulae and replaces the tables' rows with them (the caller
    commits, so readers see the old rows until then). `edge_pc` is the galaxy's sector edge. Returns the block rows
    written.
    """
    top = int(conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM nebulae").fetchone()["n"])
    codes, x, y, z, radius = sources(conn)
    rows = rules.all_rows(x, y, z, radius, codes, edge_pc, on_level=on_progress)
    clear(conn)
    written = write(conn, rows)
    _set_through(conn, top)
    return written


STATE_ID = 1

CATCH_UP_IDS = 50000
"""int: Nebula ids folded in at a time by `catch_up`: one transaction, so a stopped catch-up loses at most this many."""

KEY_BATCH = 400
"""int: Blocks read at a statement."""


def _set_through(conn, nebula_id):
    conn.execute("INSERT INTO nebula_block_state (id, nebulae_through_id) VALUES (?, ?)"
                 " ON DUPLICATE KEY UPDATE nebulae_through_id = ?", (STATE_ID, nebula_id, nebula_id))


def through_id(conn):
    """The highest `nebulae.id` the aggregates have folded in, or `None` when they were never built."""
    row = conn.execute("SELECT nebulae_through_id AS n FROM nebula_block_state WHERE id = ?", (STATE_ID,)).fetchone()
    return None if row is None else int(row["n"])


def _stored_blocks(conn, m, keys):
    """The stored `nebula_block_stats` rows of level `m` at `keys` (`(shell, band, wedge)` each), by key."""
    found = {}
    keys = list(keys)
    for start in range(0, len(keys), KEY_BATCH):
        part = keys[start:start + KEY_BATCH]
        for row in conn.execute(
                "SELECT * FROM nebula_block_stats WHERE m = ? AND (shell_index, band_index, wedge_index) IN ("
                + ", ".join("(?, ?, ?)" for _ in part) + ")", [m] + [v for key in part for v in key]).fetchall():
            found[(row["shell_index"], row["band_index"], row["wedge_index"])] = row
    return found


def _merge(old, new, volume):
    """The values of block `old` (a table row) with the batch's block dict `new` added: `(count, volume, centroid,
    bound, families, cover volume, cover)`. A bounding radius stays an upper bound: the old one carried by the
    centroid's shift, or the new members' reach from the merged centroid, whichever is larger; a `rebuild` tightens
    it."""
    old_volume = old["volume_pc3"] or 0.0
    total = old_volume + new["volume_pc3"]
    old_centre = (None if old["centroid_x_pc"] is None else
                  (old["centroid_x_pc"], old["centroid_y_pc"], old["centroid_z_pc"]))
    new_centre = (None if new["centroid_x_pc"] is None else
                  (new["centroid_x_pc"], new["centroid_y_pc"], new["centroid_z_pc"]))
    if new_centre is None:
        centre, bound = old_centre, old["bound_radius_pc"]
    elif old_centre is None or old_volume <= 0:
        centre, bound = new_centre, new["bound_radius_pc"]
    else:
        centre = tuple((o * old_volume + n * new["volume_pc3"]) / total for o, n in zip(old_centre, new_centre))
        bound = max((old["bound_radius_pc"] or 0.0) + detmath.dist(old_centre, centre),
                    detmath.dist(new_centre, centre) + new["bound_radius_pc"])
    families = [o + n for o, n in zip((old[c] for c in FAMILY_COLUMNS), new["families"])]
    cover_volume = old["cover_volume_pc3"] + new["cover_volume_pc3"]
    return (old["nebula_count"] + new["nebula_count"], total, centre, bound, families, cover_volume,
            rules.cover_of(cover_volume, volume))


UPDATE_BLOCK = ("UPDATE nebula_block_stats SET nebula_count = ?, volume_pc3 = ?, centroid_x_pc = ?, centroid_y_pc = ?,"
                " centroid_z_pc = ?, bound_radius_pc = ?, " + ", ".join(f"{column} = ?" for column in FAMILY_COLUMNS)
                + ", cover_volume_pc3 = ?, cover = ? WHERE m = ? AND shell_index = ? AND band_index = ?"
                " AND wedge_index = ?")


def _stored_classes(conn, m, keys):
    """The stored `nebula_block_classes` volumes of level `m` at `keys`: `{(shell, band, wedge): {class: volume}}`."""
    found = {}
    keys = list(keys)
    for start in range(0, len(keys), KEY_BATCH):
        part = keys[start:start + KEY_BATCH]
        for row in conn.execute(
                "SELECT shell_index, band_index, wedge_index, nebula_class, volume_pc3 FROM nebula_block_classes"
                " WHERE m = ? AND (shell_index, band_index, wedge_index) IN ("
                + ", ".join("(?, ?, ?)" for _ in part) + ")", [m] + [v for key in part for v in key]).fetchall():
            found.setdefault((row["shell_index"], row["band_index"], row["wedge_index"]), {})[
                row["nebula_class"]] = row["volume_pc3"]
    return found


def _add_classes(conn, m, added):
    """Adds `added` (`rules.LevelRows.classes` entries of level `m`) to the stored class volumes and sets the dominant
    class of each block they touch."""
    new = {}
    for shell, band, wedge, letter, volume in added:
        new.setdefault((shell, band, wedge), {})[letter] = volume
    stored = _stored_classes(conn, m, new)
    inserts, updates, dominant = [], [], []
    for key, volumes in new.items():
        held = stored.get(key, {})
        for letter, volume in volumes.items():
            if letter in held:
                updates.append((held[letter] + volume, m, *key, letter))
                held[letter] += volume
            else:
                inserts.append((m, *key, letter, volume))
                held[letter] = volume
        dominant.append((max((volume, letter) for letter, volume in held.items())[1], m, *key))
    for start in range(0, len(inserts), WRITE_BATCH):
        conn.executemany(INSERT_CLASS, inserts[start:start + WRITE_BATCH])
    for start in range(0, len(updates), WRITE_BATCH):
        conn.executemany("UPDATE nebula_block_classes SET volume_pc3 = ? WHERE m = ? AND shell_index = ?"
                         " AND band_index = ? AND wedge_index = ? AND nebula_class = ?",
                         updates[start:start + WRITE_BATCH])
    for start in range(0, len(dominant), WRITE_BATCH):
        conn.executemany("UPDATE nebula_block_stats SET dominant_class = ? WHERE m = ? AND shell_index = ?"
                         " AND band_index = ? AND wedge_index = ?", dominant[start:start + WRITE_BATCH])


def add_nebulae(conn, edge_pc, codes, x, y, z, radius):
    """
    Adds nebulae (class codes, and centres and radii in parsecs, as `sources` gives them) to the blocks they are in
    and reach, at every level: count, volume, family and class volume to the block holding the centre, the part of
    the sphere in each block to that block's cover. Returns the nebulae added.
    """
    if not len(x):
        return 0
    for m, level in rules.all_rows(x, y, z, radius, codes, edge_pc).items():
        keys = [(row["shell_index"], row["band_index"], row["wedge_index"]) for row in level.blocks]
        stored = _stored_blocks(conn, m, keys)
        inserts, updates = [], []
        for row, key in zip(level.blocks, keys):
            old = stored.get(key)
            if old is None:
                inserts.append(block_params(row))
                continue
            count, total, centre, bound, families, cover_volume, cover = _merge(
                old, row, rules.block_volume_pc3(m, *key, edge_pc))
            centre = centre or (None, None, None)
            updates.append((count, total, *centre, bound, *families, cover_volume, cover, m, *key))
        for start in range(0, len(inserts), WRITE_BATCH):
            conn.executemany(INSERT_BLOCK, inserts[start:start + WRITE_BATCH])
        for start in range(0, len(updates), WRITE_BATCH):
            conn.executemany(UPDATE_BLOCK, updates[start:start + WRITE_BATCH])
        _add_classes(conn, m, level.classes)
    return int(len(x))


NOT_IN_FIELD = " AND NOT EXISTS (SELECT 1 FROM nebula_field f WHERE f.uid = n.uid)"
"""str: Leaves out a stored nebula that `nebula_field` holds (counted from the field)."""


def catch_up(conn, edge_pc, on_progress=None):
    """
    Folds into the aggregates every nebula stored since the last time (the highest `nebulae.id` folded in is kept in
    `nebula_block_state`, the rows and the mark in one transaction, so a stopped run starts again where it stopped).
    A `galaxy` or `sector` run calls it once at its end, so its workers never queue on the shared top-level blocks.
    A nebula the field holds is left out (the field counted it). Returns the nebulae added, or `None` when the
    aggregates were never built (`rebuild` makes them).
    """
    done = through_id(conn)
    if done is None:
        return None
    top = int(conn.execute("SELECT COALESCE(MAX(id), 0) AS n FROM nebulae").fetchone()["n"])
    added = 0
    while done < top:
        upper = min(done + CATCH_UP_IDS, top)
        rows = conn.execute(
            "SELECT n.nebula_class, n.radius_ly, n.center_x_pc, n.center_y_pc, n.center_z_pc FROM nebulae n"
            " WHERE n.id > ? AND n.id <= ? AND n.center_x_pc IS NOT NULL" + NOT_IN_FIELD, (done, upper)).fetchall()
        if rows:
            codes, x, y, z, radius = _arrays([_part([(r["nebula_class"], r["radius_ly"], r["center_x_pc"],
                                                      r["center_y_pc"], r["center_z_pc"]) for r in rows])])
            added += add_nebulae(conn, edge_pc, codes, x, y, z, radius)
        _set_through(conn, upper)
        conn.commit()
        done = upper
        if on_progress is not None:
            on_progress(done, top)
    return added


def read(conn, m, cells):
    """
    The stored aggregates of the blocks `cells` (`(shell, band, wedge)` each) of level `m`, in one statement:
    `{(shell, band, wedge): {...}}` with `count`, `volume_pc3`, `centroid` (`(x, y, z)` or `None`),
    `bound_radius_pc`, `shares` (each family's share of the volume, by family name), `dominant_class` and `cover`.
    A block no nebula's sphere reaches is left out.
    """
    cells = list(cells)
    if not cells:
        return {}
    rows = conn.execute(
        f"SELECT {', '.join(BLOCK_COLUMNS)} FROM nebula_block_stats WHERE m = ? AND "
        "(shell_index, band_index, wedge_index) IN (" + ", ".join("(?, ?, ?)" for _ in cells) + ")",
        [m] + [value for cell in cells for value in cell]).fetchall()
    found = {}
    for r in rows:
        total = r["volume_pc3"]
        found[(r["shell_index"], r["band_index"], r["wedge_index"])] = {
            "count": r["nebula_count"], "volume_pc3": total,
            "centroid": None if r["centroid_x_pc"] is None else (
                r["centroid_x_pc"], r["centroid_y_pc"], r["centroid_z_pc"]),
            "bound_radius_pc": r["bound_radius_pc"],
            "shares": {family: (r[column] / total if total else 0.0)
                       for family, column in zip(rules.FAMILIES, FAMILY_COLUMNS)},
            "dominant_class": r["dominant_class"], "cover": r["cover"]}
    return found
