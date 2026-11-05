# planetgen/db/fingerprint_tree.py

"""
The stored fingerprint tree (GEN.305)
=====================================

A sector's leaf digest is taken from the sector in memory when it is saved (`leaf_digest`, over `SpaceSector.to_dict()`
and the sector's address, through the canonical encoder of GEN.304) and stored in `sectors.leaf_digest`. Above the
leaves are the nodes of `sector_tree_nodes` (DB.49): a band node (level 1: a shell and a band, over its wedges'
leaves in wedge order) and a shell node (level 2: a shell, over its band nodes in band order), each with the count of
sectors under it; the galaxy root is worked out from the shell nodes on demand (`galaxy_root`).

Saving a sector stores its leaf and marks its band and shell nodes dirty (`store_leaf`), which is one small upsert a
sector and takes no lock a neighbour waits on, so sectors saved at once by several workers do not race over a node.
`refresh` rehashes the dirty nodes from what is stored (one band node reads its up to a thousand leaves, not the
galaxy), band nodes before shell nodes; the roots and `planetgen fingerprint --stored` call it first. A node is a pure
function of the leaves under it, so the tree does not depend on the order sectors were saved in or on the number of
workers. `tree` works the whole tree out in memory from the leaves (nothing written): the stored nodes must equal it.

A rogue planet is filed in the sector it landed in, which is often not its donor's and may never have been charted, and it is written after
that sector's leaf was taken. So the leaf the tree reads (`sectors.leaf_digest`) is the content part (`sectors.content_leaf`, over the sector in
memory) joined to `resident_digest`, a digest of the rogue planets filed in the sector with their child rows and moons, whenever it holds any
(`combined`); every writer of a rogue planet calls `store_residents` for the sectors it changed, which makes the leaf again from the content part.
A sector with no content part, one made only for a landing, has a leaf of its residents alone. A sector with no rogues keeps the content part as its leaf.

A region's digest (`region_digest`) is over its sectors' leaves in address order together with their count, so adding
a sector changes it.

What is hashed is the generated content with the clocks and the collision decoration left out of every name (`plain`): which of two
colliding names keeps the plain one depends on what else is filled (GEN.57, phase 1), and the digest must not. That
stopgap goes when GEN.57 makes a sector's content independent of the order sectors are saved in. An edit to a stored
system does not rehash its sector's leaf yet: every writer marks its sector stale with the storage checksum (DB.45).
"""

import re

from planetgen.names.wordlists import DIMINUTIVE_PREFIXES, GREEK_LETTERS
from planetgen.util import fingerprint_encoding as encoding

_DECORATION = re.compile(r"\b(?:%s) (?=[A-Z])" % "|".join(GREEK_LETTERS + DIMINUTIVE_PREFIXES))
"""Pattern: A name-collision decoration word (`names.uniqueness`) before a name."""

CLOCK_KEYS = frozenset({"epoch_unix", "next_update_due", "binary_epoch_unix", "binary_next_update_due"})
"""frozenset: Keys left out of every dict: the clocks of GEN.106 record when an orbit was last moved, which for a pair
saved before the database has an orbit epoch is the moment of the save (GEN.311), not what was made (the same columns
`fingerprint.LEFT_OUT_COLUMNS` leaves out)."""

BAND, SHELL = 1, 2
"""int: The levels of `sector_tree_nodes`."""


def plain(value):
    """`value` (a tree of the encoder's types) with the collision decoration left out of every string and the clocks
    (`CLOCK_KEYS`) left out of every dict."""
    if isinstance(value, str):
        return _DECORATION.sub("", value)
    if isinstance(value, dict):
        return {key: plain(item) for key, item in value.items() if key not in CLOCK_KEYS}
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    return value


def leaf_digest(sector, address):
    """
    The 32-byte leaf digest of a saved sector at `address` (`(shell, band, wedge)`), from the sector in memory: its
    `to_dict()` with the decoration left out of the names, and the address.
    """
    value = {"address": list(address), "sector": plain(sector.to_dict())}
    return bytes.fromhex(encoding.digest(value))


def _node_digest(kind, key, count, children):
    """The 32-byte digest of a node: its kind and key, the count under it and its children's `[key, digest]`s."""
    value = {"node": kind, "key": list(key), "count": count, "children": [[list(child), digest.hex()]
                                                                           for child, digest in children]}
    return bytes.fromhex(encoding.digest(value))


def _band(shell, band, leaves):
    """`(digest, count)` of a band node over `leaves`, `[(wedge, digest)]` in wedge order."""
    return _node_digest("band", (shell, band), len(leaves), [((wedge,), digest) for wedge, digest in leaves]), len(leaves)


def _shell(shell, bands):
    """`(digest, count)` of a shell node over `bands`, `[(band, digest, count)]` in band order."""
    count = sum(n for _band_index, _digest, n in bands)
    return _node_digest("shell", (shell,), count, [((band,), digest) for band, digest, _n in bands]), count


def _root(shells):
    """`(digest, count)` of the galaxy root over `shells`, `[(shell, digest, count)]` in shell order."""
    count = sum(n for _shell_index, _digest, n in shells)
    return _node_digest("galaxy", (), count, [((shell,), digest) for shell, digest, _n in shells]), count


def _mark_dirty(conn, address):
    shell, band, _wedge = address
    for level, node_band in ((BAND, band), (SHELL, 0)):
        conn.execute(
            "INSERT INTO sector_tree_nodes (level, shell_index, band_index, dirty) VALUES (?, ?, ?, 1)"
            " ON DUPLICATE KEY UPDATE dirty = 1", (level, shell, node_band))


def combined(content, residents):
    """
    The leaf of a sector from its content part (bytes or `None`) and its residents' digest (bytes or `None`): the content part alone when it
    has no residents, `None` when it has neither, else a digest of both.
    """
    if residents is None:
        return content
    return bytes.fromhex(encoding.digest({"content": None if content is None else content.hex(), "residents": residents.hex()}))


_RUN_TIME_COLUMNS = frozenset({"uid", "parent_id", "rogue_planet_id", "moon_id"})
"""frozenset: Columns left out of a rogue planet's rows in `resident_digest`: its run-time object ID (`uid`) comes from a counter that moves each time
the donor is saved again, and the rows that point at it (its child rows, its moons' `parent_id`) carry it, so none of them says what was made."""


def resident_digest(conn, sector_id):
    """
    The 32-byte digest of the rogue planets filed in a sector, or `None` when it holds none. Each rogue is taken whole: its `rogue_planets` row,
    its child rows, its moons and their child rows, as `fingerprint` writes a row (a foreign key by the key of what it points at, names without
    the collision decoration), less the run-time object IDs and the columns that carry them (`_RUN_TIME_COLUMNS`); the rogues' digests are put in
    canonical order, so the order they landed in does not count.
    """
    from planetgen.db import fingerprint
    reader = fingerprint._Reader(conn, plain_names=True)
    rogues = reader.select("rogue_planets", "sector_id = ?", (sector_id,))
    if not rogues:
        return None

    def bare(table, row):
        value = reader.row_value(table, row)
        value["row"] = {column: item for column, item in value["row"].items() if column not in _RUN_TIME_COLUMNS}
        return value

    moons = reader.select_in("moons", "rogue_planet_id", [row["id"] for row in rogues])
    moon_children = {}
    for table, key, owner in fingerprint.CHILD_TABLES:
        if owner == "moons":
            for row in reader.select_in(table, key, [moon["id"] for moon in moons]):
                moon_children.setdefault(row[key], []).append(bare(table, row))
    own_children = {}
    for table, key, owner in fingerprint.CHILD_TABLES:
        if owner == "rogue_planets":
            for row in reader.select_in(table, key, [rogue["id"] for rogue in rogues]):
                own_children.setdefault(row[key], []).append(bare(table, row))
    moons_of = {}
    for moon in moons:
        moons_of.setdefault(moon["rogue_planet_id"], []).append(
            {"moon": bare("moons", moon), "rows": encoding.sorted_rows(moon_children.get(moon["id"], []))})
    digests = []
    for rogue in rogues:
        value = {"rogue": bare("rogue_planets", rogue), "rows": encoding.sorted_rows(own_children.get(rogue["id"], [])),
                 "moons": encoding.sorted_rows(moons_of.get(rogue["id"], []))}
        digests.append({"rogue": encoding.digest(value)})
    return bytes.fromhex(encoding.row_set_digests(digests).full)


def store_leaf(conn, sector_id, address, content):
    """
    Stores a sector's content part and the leaf made from it and the rogue planets filed in the sector (`combined`), and marks its band and
    shell nodes dirty. `content` is `leaf_digest(sector, address)` of the sector as it is saved.
    """
    conn.execute("UPDATE sectors SET content_leaf = ?, leaf_digest = ? WHERE id = ?",
                 (content, combined(content, resident_digest(conn, sector_id)), sector_id))
    _mark_dirty(conn, address)


def drop_content(conn, sector_id):
    """Forgets the content part of a sector that was emptied but stays for the objects filed in it, and makes its leaf again from them (`store_residents`)."""
    conn.execute("UPDATE sectors SET content_leaf = NULL WHERE id = ?", (sector_id,))
    store_residents(conn, sector_id)


def store_residents(conn, sector_id):
    """
    Makes the leaf of a sector again after a rogue planet was written in it or deleted from it, from the content part stored with it (none for
    a sector made only for a landing), and marks its band and shell nodes dirty. A no-op for a sector with no address. Run it after the rogues'
    object IDs are assigned, which their rows carry.
    """
    row = conn.execute("SELECT shell_index, band_index, wedge_index, content_leaf FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if row is None or row["shell_index"] is None:
        return
    content = None if row["content_leaf"] is None else bytes(row["content_leaf"])
    conn.execute("UPDATE sectors SET leaf_digest = ? WHERE id = ?", (combined(content, resident_digest(conn, sector_id)), sector_id))
    _mark_dirty(conn, (row["shell_index"], row["band_index"], row["wedge_index"]))


def _leaves(conn, shell, band):
    return [(row["wedge_index"], bytes(row["leaf_digest"])) for row in conn.execute(
        "SELECT wedge_index, leaf_digest FROM sectors WHERE shell_index = ? AND band_index = ?"
        " AND leaf_digest IS NOT NULL ORDER BY wedge_index", (shell, band)).fetchall()]


def refresh(conn):
    """
    Rehashes every dirty node from what is stored: the band nodes from their leaves, then the shell nodes from their
    band nodes. A node with nothing under it is deleted. Returns the number of nodes written.
    """
    written = 0
    for row in conn.execute("SELECT shell_index, band_index FROM sector_tree_nodes WHERE level = ? AND dirty = 1",
                            (BAND,)).fetchall():
        shell, band = row["shell_index"], row["band_index"]
        leaves = _leaves(conn, shell, band)
        if not leaves:
            conn.execute("DELETE FROM sector_tree_nodes WHERE level = ? AND shell_index = ? AND band_index = ?",
                         (BAND, shell, band))
        else:
            digest, count = _band(shell, band, leaves)
            conn.execute("UPDATE sector_tree_nodes SET digest = ?, child_count = ?, dirty = 0"
                         " WHERE level = ? AND shell_index = ? AND band_index = ?", (digest, count, BAND, shell, band))
        written += 1
    for row in conn.execute("SELECT shell_index FROM sector_tree_nodes WHERE level = ? AND dirty = 1",
                            (SHELL,)).fetchall():
        shell = row["shell_index"]
        bands = [(r["band_index"], bytes(r["digest"]), r["child_count"]) for r in conn.execute(
            "SELECT band_index, digest, child_count FROM sector_tree_nodes WHERE level = ? AND shell_index = ?"
            " AND dirty = 0 ORDER BY band_index", (BAND, shell)).fetchall()]
        if not bands:
            conn.execute("DELETE FROM sector_tree_nodes WHERE level = ? AND shell_index = ?", (SHELL, shell))
        else:
            digest, count = _shell(shell, bands)
            conn.execute("UPDATE sector_tree_nodes SET digest = ?, child_count = ?, dirty = 0"
                         " WHERE level = ? AND shell_index = ? AND band_index = 0", (digest, count, SHELL, shell))
        written += 1
    return written


def galaxy_root(conn):
    """`(hex digest, sector count)` of the galaxy root over the shell nodes, after `refresh`; `(None, 0)` for a
    galaxy with no stored leaf."""
    refresh(conn)
    shells = [(row["shell_index"], bytes(row["digest"]), row["child_count"]) for row in conn.execute(
        "SELECT shell_index, digest, child_count FROM sector_tree_nodes WHERE level = ? AND digest IS NOT NULL"
        " ORDER BY shell_index", (SHELL,)).fetchall()]
    if not shells:
        return None, 0
    digest, count = _root(shells)
    return digest.hex(), count


def tree(conn):
    """
    The whole tree worked out in memory from the stored leaves, nothing written: `{"bands": {(shell, band): (hex,
    count)}, "shells": {shell: (hex, count)}, "root": (hex, count) or None}`.
    """
    leaves = {}
    for row in conn.execute("SELECT shell_index, band_index, wedge_index, leaf_digest FROM sectors"
                            " WHERE leaf_digest IS NOT NULL AND shell_index IS NOT NULL"
                            " ORDER BY shell_index, band_index, wedge_index").fetchall():
        leaves.setdefault((row["shell_index"], row["band_index"]), []).append((row["wedge_index"],
                                                                              bytes(row["leaf_digest"])))
    bands, by_shell = {}, {}
    for (shell, band), found in leaves.items():
        digest, count = _band(shell, band, found)
        bands[shell, band] = (digest.hex(), count)
        by_shell.setdefault(shell, []).append((band, digest, count))
    shells = {}
    for shell, found in by_shell.items():
        digest, count = _shell(shell, sorted(found))
        shells[shell] = (digest.hex(), count)
    root = None
    if shells:
        digest, count = _root([(shell, bytes.fromhex(hexed), n) for shell, (hexed, n) in sorted(shells.items())])
        root = (digest.hex(), count)
    return {"bands": bands, "shells": shells, "root": root}


def region_digest(conn, shells=None, addresses=None):
    """
    `(hex digest, sector count)` of a region's stored leaves in address order, with the count in the digest: the
    sectors on `shells` and at `addresses` (`(shell, band, wedge)`), or with neither the whole galaxy.
    """
    shells = None if shells is None else set(shells)
    addresses = None if addresses is None else {tuple(a) for a in addresses}
    leaves = []
    for row in conn.execute("SELECT shell_index, band_index, wedge_index, leaf_digest FROM sectors"
                            " WHERE leaf_digest IS NOT NULL AND shell_index IS NOT NULL"
                            " ORDER BY shell_index, band_index, wedge_index").fetchall():
        address = (row["shell_index"], row["band_index"], row["wedge_index"])
        if shells is None and addresses is None or address[0] in (shells or ()) or address in (addresses or ()):
            leaves.append((address, bytes(row["leaf_digest"])))
    value = {"region": "sectors", "count": len(leaves), "leaves": [[list(a), d.hex()] for a, d in leaves]}
    return encoding.digest(value), len(leaves)


def lines(conn, shells=None, addresses=None):
    """The printed form of the stored tree (`planetgen fingerprint --stored`): a line a band, a line a shell, the root
    and the region's digest; computed in memory, nothing written."""
    found = tree(conn)
    out = [f"band {shell} {band} {digest} ({count} sector(s))" for (shell, band), (digest, count) in found["bands"].items()]
    out += [f"shell {shell} {digest} ({count} sector(s))" for shell, (digest, count) in found["shells"].items()]
    if found["root"] is not None:
        out.append(f"root {found['root'][0]} ({found['root'][1]} sector(s))")
    digest, count = region_digest(conn, shells, addresses)
    out.append(f"stored-region {digest} ({count} sector(s))")
    return out
