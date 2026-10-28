# planetgen/db/table_rebuild.py

"""
Rebuilding a table whose clustered index has a page the server cannot read (DB.48; docs/design/operations/db-check-and-parity-repair.md
sections 4 and 5.11). `planetgen repair-db --salvage TABLE` and `--reload TABLE` are `salvage` and `reload`.

**Why two steps.** SQL can neither read nor delete the rows on a bad clustered-index page, and a rebuild of one sector in place does
not help. The server has to be told to skip the page (`innodb_force_recovery = 1` in the `[mysqld]` section of its configuration, then
a restart), and while it runs that way InnoDB refuses `INSERT`, `UPDATE` and `DELETE` (only reading and DDL work). So:

1. **Salvage** (server in force recovery, website and generation stopped): every readable row of the table is written to a dump file
   (`<dir>/<table>.salvage.jsonl`, each value in `storage_checksum.encode_value`'s tagged form, so nothing is rounded), the file is
   synced and renamed into place, and only then is the table dropped (foreign key checks off) and made again empty from `schema.sql`.
2. **Reload** (after a normal restart, `innodb_force_recovery` removed): the rows are inserted again with their original ids, foreign
   key checks off, and the count is checked. Then `repair-db` runs as usual: a sector that lost rows with the page no longer hashes to
   its checksum, and the parity file rebuilds it, lost rows included (`repair.apply_export` inserts them).

`reload` of `sectors` itself (`restore_sectors=True`, `--restore-sectors`) also puts back the sector rows the parity file holds and the
table no longer has: the one loss `check-db` cannot see, since a sector with no row is not in any list. It is opt-in because a
sector deleted while `parity.enabled` was off is still a member of the file and would come back.

Only the tables a sector's checksum covers (`storage_checksum.CONTENT_TABLES`) can be rebuilt this way, because the parity file is what
brings back the rows the page took. A damaged secondary index needs no parity: `ALTER TABLE ... FORCE` or drop and re-add the index.
"""

import json
import os
import time

from planetgen.db import repair, sector_export, storage_checksum, store
from planetgen.util import log

CHUNK = 5000
BATCH = 500
SUFFIX = ".salvage.jsonl"
LOADED = ".loaded"


class RebuildRefused(Exception):
    """A step that was not done, and why (nothing was changed)."""


def dump_path(directory, table):
    return os.path.join(directory, table + SUFFIX)


def force_recovery(conn):
    """The server's `innodb_force_recovery` (0 when it is off)."""
    row = conn.execute("SELECT @@innodb_force_recovery AS level").fetchone()
    conn.commit()
    return int(row["level"])


def _table_columns(conn, table):
    """`(stored columns in table order, primary key columns)`: generated columns are left out, they are worked out from the rest."""
    rows = conn.execute(
        "SELECT column_name AS name, column_key AS k, extra AS x FROM information_schema.columns"
        " WHERE table_schema = DATABASE() AND table_name = ? ORDER BY ordinal_position", (table,)).fetchall()
    conn.commit()
    if not rows:
        raise RebuildRefused(f"there is no table {table}")
    columns = [row["name"] for row in rows if "GENERATED" not in (row["x"] or "").upper()]
    key = [row["name"] for row in rows if row["k"] == "PRI"]
    return columns, key


def _check_table_name(table):
    if table not in storage_checksum.CONTENT_TABLES:
        raise RebuildRefused(f"{table} is not a table the parity file can bring rows back for; the tables are "
                             + ", ".join(storage_checksum.CONTENT_TABLES))


def probe_unreadable(conn, say=log.normal):
    """The content tables whose clustered index cannot be read through (a count forced onto the primary key fails). Reads every
    table once."""
    bad = []
    for table in storage_checksum.CONTENT_TABLES:
        try:
            conn.execute(f"SELECT COUNT(*) AS n FROM `{table}` FORCE INDEX (PRIMARY)").fetchone()
        except repair.READ_ERRORS as exc:
            conn.rollback()
            say(f"Repair: {table} cannot be read: {type(exc).__name__}: {exc}")
            bad.append(table)
        conn.commit()
    return bad


def _scan(conn, table, columns, key, last=None):
    """The next `CHUNK` rows of `table` in key order after `last` (a key tuple), as dicts."""
    names = ", ".join(f"`{c}`" for c in columns)
    order = ", ".join(f"`{c}`" for c in key)
    if last is None:
        sql, params = f"SELECT {names} FROM `{table}` ORDER BY {order} LIMIT {CHUNK}", ()
    else:
        after = f"({', '.join(f'`{c}`' for c in key)}) > ({', '.join('?' * len(key))})"
        sql, params = f"SELECT {names} FROM `{table}` WHERE {after} ORDER BY {order} LIMIT {CHUNK}", tuple(last)
    rows = conn.execute(sql, params).fetchall()
    conn.commit()
    return rows


def _line(row, columns):
    return json.dumps([storage_checksum.encode_value(row[column]).decode("utf-8") for column in columns])


def salvage(conn, table, directory, say=log.normal):
    """
    Step 1: dumps every readable row of `table` to `directory`, then drops the table and makes it again empty. Returns the
    number of rows dumped.

    Raises:
        RebuildRefused: The table cannot be rebuilt this way, a dump of it is already waiting for a reload, the scan fails and the
            server is not in force recovery, or the dumped count differs from a count the server gives.
    """
    _check_table_name(table)
    columns, key = _table_columns(conn, table)
    if not key:
        raise RebuildRefused(f"{table} has no primary key to read it in order by")
    path = dump_path(directory, table)
    if os.path.exists(path):
        raise RebuildRefused(f"{path} is already waiting for `--reload {table}`; reload it first (or remove the file if it is "
                             "a dump you no longer want)")
    level = force_recovery(conn)
    os.makedirs(directory, exist_ok=True)
    partial = path + ".partial"
    count, last = 0, None
    try:
        with open(partial, "w", encoding="utf-8") as out:
            out.write(json.dumps({"table": table, "columns": columns, "key": key, "created": time.time(), "force_recovery": level}) + "\n")
            while True:
                try:
                    rows = _scan(conn, table, columns, key, last)
                except repair.READ_ERRORS as exc:
                    conn.rollback()
                    if level == 0:
                        raise RebuildRefused(
                            f"{table} cannot be read ({type(exc).__name__}: {exc}). Stop the website and any generation, add "
                            "`innodb_force_recovery = 1` to the [mysqld] section of the server's configuration, restart the "
                            f"server and run `planetgen repair-db --salvage {table}` again.") from exc
                    raise RebuildRefused(
                        f"{table} still cannot be read after {count:,} rows with innodb_force_recovery = {level} "
                        f"({type(exc).__name__}: {exc}); a higher level (2 or 3) may skip more; nothing was dropped.") from exc
                if not rows:
                    break
                for row in rows:
                    out.write(_line(row, columns) + "\n")
                count += len(rows)
                last = tuple(rows[-1][c] for c in key)
                if count % (CHUNK * 20) < CHUNK:
                    say(f"Salvage: {count:,} rows of {table}.")
            out.flush()
            os.fsync(out.fileno())
        if level == 0:                           # a server that reads normally can say how many rows there should be
            expected = conn.execute(f"SELECT COUNT(*) AS n FROM `{table}`").fetchone()["n"]
            conn.commit()
            if expected != count:
                raise RebuildRefused(f"{count:,} rows were dumped but {table} holds {expected:,}; nothing was dropped")
        os.replace(partial, path)
    finally:
        if os.path.exists(partial):
            os.unlink(partial)
    say(f"Salvage: {count:,} row(s) of {table} dumped to {path}.")
    conn.execute("SET foreign_key_checks = 0")
    try:
        conn.execute(f"DROP TABLE `{table}`")
        conn.commit()
        store._ensure_schema(conn)
    finally:
        conn.execute("SET foreign_key_checks = 1")
    columns_after, _ = _table_columns(conn, table)
    if columns_after != columns:
        raise RebuildRefused(f"{table} was made again with other columns than it had ({columns_after} against {columns}); "
                             f"the dump is kept at {path}")
    say(f"Salvage: {table} made again, empty. Restart the server without innodb_force_recovery, then run "
        f"`planetgen repair-db --reload {table}`.")
    return count


def _read_dump(path):
    with open(path, "r", encoding="utf-8") as handle:
        header = json.loads(handle.readline())
        yield header
        for line in handle:
            yield line


def reload(conn, table, directory, say=log.normal):
    """
    Step 2: inserts the rows of the dump back into the empty table, with their ids and foreign key checks off, and checks the
    count. A reload that was interrupted is finished by running it again (rows already in the table are skipped). The dump is
    renamed `<table>.salvage.jsonl.loaded` when every row is in. Returns the number of rows loaded.

    Raises:
        RebuildRefused: The server is still in force recovery, there is no dump, the table differs from the dump's, or the
            rows in the table are not a part of the dump's.
    """
    _check_table_name(table)
    path = dump_path(directory, table)
    if not os.path.exists(path):
        raise RebuildRefused(f"there is no dump of {table} at {path}; `--salvage {table}` makes one")
    level = force_recovery(conn)
    if level:
        raise RebuildRefused(f"the server is still running with innodb_force_recovery = {level}, which refuses inserts; restart it "
                             "without that setting first")
    columns, key = _table_columns(conn, table)
    reader = _read_dump(path)
    header = next(reader)
    if header["table"] != table or header["columns"] != columns:
        raise RebuildRefused(f"the dump's columns are not the table's now ({header['columns']} against {columns})")
    present = conn.execute(f"SELECT COUNT(*) AS n FROM `{table}`").fetchone()["n"]
    conn.commit()
    ignore = "IGNORE " if present else ""
    if present:
        say(f"Reload: {table} already holds {present:,} row(s); finishing an interrupted reload.")
    sql = (f"INSERT {ignore}INTO `{table}` ({', '.join(f'`{c}`' for c in columns)}) "
           f"VALUES ({', '.join('?' * len(columns))})")
    total, batch = 0, []
    conn.execute("SET foreign_key_checks = 0")
    try:
        for line in reader:
            batch.append(tuple(sector_export.decode_value(text.encode("utf-8")) for text in json.loads(line)))
            if len(batch) == BATCH:
                conn.executemany(sql, batch)
                conn.commit()
                total += len(batch)
                batch = []
                if total % (BATCH * 100) == 0:
                    say(f"Reload: {total:,} rows of {table}.")
        if batch:
            conn.executemany(sql, batch)
            conn.commit()
            total += len(batch)
    finally:
        conn.execute("SET foreign_key_checks = 1")
    held = conn.execute(f"SELECT COUNT(*) AS n FROM `{table}`").fetchone()["n"]
    conn.commit()
    if held != total:
        raise RebuildRefused(f"{total:,} rows were read from the dump but {table} holds {held:,}; the dump is kept at {path}")
    os.replace(path, path + LOADED)
    say(f"Reload: {total:,} row(s) of {table} loaded.")
    return total


def restore_sectors(conn, dry_run=False, say=log.normal):
    """
    Puts back the sector rows the parity file holds and the `sectors` table lacks, each rebuilt from its group like a damaged sector
    (`repair.apply_export` inserts the `sectors` row first). Returns the `repair.Result`s.
    """
    from planetgen.db import parity, parity_file as pfile
    parity_store = parity.store()
    present = {row["id"] for row in conn.execute("SELECT id FROM sectors").fetchall()}
    conn.commit()
    missing = {}
    with parity_store.file.locked():
        for sector_id, (group_id, _slot) in parity_store.file.sector_group.items():
            if sector_id != pfile.EMPTY and sector_id not in present:
                missing.setdefault(group_id, []).append(sector_id)
    say(f"Restore: {sum(len(ids) for ids in missing.values()):,} sector row(s) are in the parity file and not in the table.")
    results = []
    for group_id, ids in sorted(missing.items()):
        results.extend(repair._repair_group(conn, parity_store, group_id, sorted(ids), dry_run))
    for result in results:
        say("Restore: " + result.line())
    return results
