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

A region's digest (`region_digest`) is over its sectors' leaves in address order together with their count, so adding
a sector changes it.

What is hashed is the generated content with the collision decoration left out of every name (`plain`): which of two
colliding names keeps the plain one depends on what else is filled (GEN.57, phase 1), and the digest must not. That
stopgap goes when GEN.57 makes a sector's content independent of the order sectors are saved in. An edit to a stored
system does not rehash its sector's leaf yet: every writer marks its sector stale with the storage checksum (DB.45).
"""

import re

from planetgen.names.wordlists import DIMINUTIVE_PREFIXES, GREEK_LETTERS
from planetgen.util import fingerprint_encoding as encoding

_DECORATION = re.compile(r"\b(?:%s) (?=[A-Z])" % "|".join(GREEK_LETTERS + DIMINUTIVE_PREFIXES))
"""Pattern: A name-collision decoration word (`names.uniqueness`) before a name."""

BAND, SHELL = 1, 2
"""int: The levels of `sector_tree_nodes`."""


def plain(value):
    """`value` (a tree of the encoder's types) with the collision decoration left out of every string."""
    if isinstance(value, str):
        return _DECORATION.sub("", value)
    if isinstance(value, dict):
        return {key: plain(item) for key, item in value.items()}
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


def store_leaf(conn, sector_id, address, digest):
    """Stores a sector's leaf digest and marks its band and shell nodes dirty."""
    shell, band, _wedge = address
    conn.execute("UPDATE sectors SET leaf_digest = ? WHERE id = ?", (digest, sector_id))
    for level, node_band in ((BAND, band), (SHELL, 0)):
        conn.execute(
            "INSERT INTO sector_tree_nodes (level, shell_index, band_index, dirty) VALUES (?, ?, ?, 1)"
            " ON DUPLICATE KEY UPDATE dirty = 1", (level, shell, node_band))


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
