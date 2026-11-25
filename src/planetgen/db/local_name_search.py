# planetgen/db/local_name_search.py

"""
Search by a bare local name (API.38).

A local name (GEN.1255, `names/naming_key.py`) spells the serial and body number of an Object ID and nothing of its birth sector, and the spelling of the serial is
scrambled by a keyed bijection of its own for every sector. So a name typed without its sector stands for one Object ID per sector, and the only way to find what has it
is to undo the scramble of each sector in turn and look the Object ID up. Nothing is stored in words (a key change renames every object without a rewrite), so there is
no index to read the name from; this module is that pass.

* **The sectors to search** (`cells`) are every sector row (charted cells and the placeholders of the cells a drifting free planet landed in), every cell of the scatter tables
  (`bright_stars`, `phenomenon_scatter`: uncharted space has no sector row) and every birth sector of the standalone phenomena (a field cloud or a container is born in the
  sector that holds its centre, which may have no row). Each is counted once.
* **Nearest first** (`ordered`): the sectors are ordered by the distance of their centre from the sector the person is in (`near`), so the answer for the neighbourhood is
  there at once and the rest follows.
* **The pass** (`search`) takes the sectors in batches of `BATCH`: the Object ID the name stands for in each sector, then one indexed read of each table that holds ID-named
  objects (`TABLES`) for the whole batch. It reports `searched` of `total` after each batch and can stop at a deadline, saying how many it did, so a caller can carry it on
  (`start`). It never stops without saying so.
"""

import time

from planetgen.db import kind_tables
from planetgen.galaxy import object_uid, spherical, uid as sector_uid_module
from planetgen.names import naming_key

SCATTER_TABLES = ("bright_stars", "phenomenon_scatter")
"""tuple: The scatter tables: a row has the cell it was placed in (`shell_index`, `band_index`, `wedge_index`) and an Object ID of its own."""

TABLES = tuple(kind_tables.tables("phenomenon")) + SCATTER_TABLES
"""tuple: The tables whose objects are named by their Object ID (GEN.1256, GEN.1254): the standalone phenomena and the scatter. Systems, stars, planets, moons and belts keep their
generated names, so a local name never names them."""

BATCH = 2000
"""int: Sectors looked at between two reports of the pass. One read of each table per batch."""


def cells(conn):
    """
    The sectors that can hold an object named by its Object ID, each once, as the 56-bit integers of their designations (`galaxy.uid.sector_uid`).

    Returns:
        list: The designations, in no particular order.
    """
    found = set()
    for row in conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors WHERE shell_index IS NOT NULL").fetchall():
        found.add(_designation(row["shell_index"], row["band_index"], row["wedge_index"]))
    for table in SCATTER_TABLES:
        for row in conn.execute(f"SELECT shell_index, band_index, wedge_index FROM {table} GROUP BY shell_index, band_index, wedge_index").fetchall():
            found.add(_designation(row["shell_index"], row["band_index"], row["wedge_index"]))
    for table in kind_tables.tables("phenomenon"):
        for row in conn.execute(f"SELECT DISTINCT SUBSTRING(uid, 1, {object_uid.LAYOUT.sector_bits // 8}) AS sector FROM {table} WHERE uid IS NOT NULL").fetchall():
            found.add(int.from_bytes(bytes(row["sector"]), "big"))
    found.discard(None)
    return list(found)


def _designation(shell, band, wedge):
    """The designation of a cell as an integer, or `None` for an address the grid does not have (a row of an older grid)."""
    try:
        return sector_uid_module.sector_uid(shell, band, wedge)
    except ValueError:
        return None


def ordered(designations, near, edge_pc):
    """
    `designations` nearest first.

    Args:
        designations (iterable): Sector designations (integers).
        near (int or None): The designation of the sector the person is in; `None` for the galaxy's core.
        edge_pc (float): The edge of a sector in parsecs (`galaxy_shape.edge_pc`), for the centres.

    Returns:
        list: The same designations, the nearest centre first, a sector whose address the grid does not hold last.
    """
    if near is None:
        origin = (0.0, 0.0, 0.0)
    else:
        origin = spherical.sector_position_pc(*sector_uid_module.sector_address(near), edge_pc)
    keyed = []
    for designation in designations:
        try:
            x, y, z = spherical.sector_position_pc(*sector_uid_module.sector_address(designation), edge_pc)
            distance = (x - origin[0]) ** 2 + (y - origin[1]) ** 2 + (z - origin[2]) ** 2
        except ValueError:
            distance = float("inf")
        keyed.append((distance, designation))
    keyed.sort()
    return [designation for _distance, designation in keyed]


def parts_of(text, key):
    """What `text` spells as a local name (`naming_key.local_parts`), or `None` when it is not one: not two words the codec reads, or no naming key."""
    if key is None or not isinstance(text, str) or not text.strip():
        return None
    try:
        return naming_key.local_parts(text, key)
    except Exception:  # noqa: BLE001 -- the codec rejects what it cannot read as a local name; ordinary text is the common case
        return None


def search(conn, parts, key, designations, start=0, deadline=None, progress=None, batch=BATCH):
    """
    Looks for the objects whose local name `parts` spell, sector by sector from `start`.

    Args:
        conn (Connection): An open connection.
        parts (tuple): `naming_key.local_parts` of the typed name.
        key (str): The naming key in force.
        designations (list): The sectors, in the order to take them (`ordered`).
        start (int): How many of them an earlier pass already took.
        deadline (float, optional): A `time.monotonic()` after which no new batch is begun.
        progress (callable, optional): `progress(searched, total)` after each batch.
        batch (int): Sectors per read.

    Returns:
        dict: `matches` (`[{"table", "id", "uid", "sector"}]`, `uid` an integer and `sector` the 14 hex digits of its birth sector), `searched` (how many of `designations` have
            been taken, counting `start`), `total` and `complete` (every sector taken).
    """
    matches = []
    total = len(designations)
    position = start
    while position < total:
        if deadline is not None and time.monotonic() >= deadline:
            break
        chunk = designations[position:position + batch]
        wanted = {}
        for designation in chunk:
            wanted[object_uid.to_bytes(naming_key.uid_in_sector(parts, designation, key, cached=False))] = designation
        raws = list(wanted)
        marks = ", ".join(["?"] * len(raws))
        for table in TABLES:
            for row in conn.execute(f"SELECT id, uid FROM {table} WHERE uid IN ({marks})", raws).fetchall():
                uid = object_uid.from_bytes(bytes(row["uid"]))
                matches.append({"table": table, "id": row["id"], "uid": uid, "sector": naming_key.sector_text(uid)})
        position += len(chunk)
        if progress is not None:
            progress(position, total)
    return {"matches": matches, "searched": position, "total": total, "complete": position >= total}
