# planetgen/db/repair.py

"""
Repair of damaged sectors from the parity file (DB.47; docs/design/operations/db-check-and-parity-repair.md sections 5.1, 5.2 and
5.6). `planetgen repair-db` is `repair`.

**What is damaged.** A sector whose `parity_dirty` flag is clear and whose rows no longer hash to the storage checksum sealed
for it (`storage_checksum.check_sector`, the check DB.8's `storage` check makes). A sector that is merely marked or never sealed
is not damage and is left to `refresh-parity`.

**What is done about it**, by the decision table of design 5.6 (`parity.classify`):

- the rows hash to what the parity record holds for the sector, but not to the stored checksum: the checksum column is the
  damaged part, two independent copies agree on the rows, and the checksum is sealed again (`RESEALED`);
- the record holds the checksum stored for the sector and the rows differ: real damage with current parity. The group's other
  members are read and each must still hash to what the record holds for it; every block that does not (damaged, changed since
  the record was made, or gone) counts as an erasure, and up to `m` erasures are mended together (`util/rscodec.RS.decode`).
  The rebuilt block must hash to the stored checksum before anything is written (`REPAIRED`, `WOULD_REPAIR` under `--dry-run`);
- the record holds another checksum: damage and a stale record at once, the parity cannot be trusted for the sector, and it
  is reported (`UNREPAIRABLE`); the way out is the seed regeneration or a backup. So are a sector in no group, a group with more
  erasures than parity blocks, and a record that fails its CRC;
- a table that cannot be read at all (a clustered-index page the server refuses) is `UNREADABLE`, DB.48's case.

**How a rebuilt block is written.** The block is read back into rows (`sector_export.parse`) and made the rows of the sector by
primary key, in one transaction under the sector's row lock and the parity file's lock: a row that differs is updated in its
covered columns, a row that is missing is inserted, and a row of the sector the block does not hold is deleted. Only the columns
the checksum covers are written, so the columns the orbit update rewrites (design 5.2) are left exactly as they are on every row
that still exists; a row that has to be inserted gets them as the column's default, zero where there is none (and its
`sector_id` the sector's), and the next orbit update carries on from there (the design's option (a): a repaired row carries
no history of its orbit). The sector is checked against the checksum again before the transaction commits, and rolled back if
it differs.

**The parity file itself.** A record that fails its CRC, or a file damaged behind its last record, is made again from the
database: the file is compacted (the unreadable records are dropped) and `parity.refresh` protects every sector the file no
longer holds, which it does only for a sector that still hashes to its checksum.
"""

import dataclasses
import datetime
import os

import sqlalchemy.exc

from planetgen.db import parity, parity_file as pfile, sector_export, storage_checksum, store
from planetgen.db.check import Scope
from planetgen.util import log

READ_ERRORS = (store.pymysql.err.MySQLError, sqlalchemy.exc.SQLAlchemyError)
"""tuple: What reading a damaged table raises."""

REPAIRED = "repaired"
RESEALED = "resealed"
WOULD_REPAIR = "would repair"
UNREPAIRABLE = "unrepairable"
UNREADABLE = "unreadable"
PROGRESS_EVERY = 1000


class RepairFailed(Exception):
    """A rebuilt block that could not be made the rows of its sector (nothing was changed)."""


@dataclasses.dataclass
class Result:
    """What became of one damaged sector: `outcome` is one of the constants above, `reason` says why not for the last two."""
    sector_id: int
    outcome: str
    reason: str = ""

    def line(self):
        return f"sector id {self.sector_id}: {self.outcome}" + (f" ({self.reason})" if self.reason else "")


@dataclasses.dataclass
class Report:
    """The whole run: a `Result` for each damaged sector, the number of sectors that were only unprotected, and what was done
    about the parity file (`records_damaged` records that failed their CRC, `file_damaged_at` the byte the file stops
    parsing at, `file_remade` whether it was compacted and filled again)."""
    results: list = dataclasses.field(default_factory=list)
    unprotected: int = 0
    records_damaged: int = 0
    file_damaged_at: int = None
    file_remade: bool = False

    def count(self, outcome):
        return sum(1 for result in self.results if result.outcome == outcome)

    @property
    def clean(self):
        """True when nothing is left damaged (a dry run's `WOULD_REPAIR` is still damage)."""
        return not any(result.outcome in (WOULD_REPAIR, UNREPAIRABLE, UNREADABLE) for result in self.results)


def find_damaged(conn, scope=None, say=log.normal):
    """
    The sectors in `scope` that are damaged, as `storage_checksum.check_sector` tells them.

    Returns:
        tuple: `(damaged ids, unreadable ids, unprotected count)`; an unreadable sector is one whose tables the server refuses to read.
    """
    scope_sql, params = (scope or Scope()).where("s.id")
    rows = conn.execute("SELECT s.id AS id, s.parity_dirty AS dirty, s.storage_checksum AS sum FROM sectors s"
                        f" WHERE 1 = 1{scope_sql} ORDER BY s.id", params).fetchall()
    conn.commit()
    damaged, unreadable, unprotected = [], [], 0
    say(f"Repair: checking {len(rows):,} sector(s).")
    for done, row in enumerate(rows, 1):
        if row["dirty"] or row["sum"] is None:
            unprotected += 1
        else:
            try:
                if storage_checksum.check_sector(conn, row["id"]) == storage_checksum.DAMAGED:
                    damaged.append(row["id"])
            except READ_ERRORS as exc:
                conn.rollback()
                unreadable.append((row["id"], f"{type(exc).__name__}: {exc}"))
        if done % PROGRESS_EVERY == 0:
            say(f"Repair: checked {done:,} of {len(rows):,}.")
    conn.commit()
    return damaged, unreadable, unprotected


# ---------------------------------------------------------------------
# Writing a rebuilt block back
# ---------------------------------------------------------------------

_NUMERIC = ("tinyint", "smallint", "mediumint", "int", "bigint", "decimal", "numeric", "float", "double", "bit", "year")
_BLOB = ("tinyblob", "blob", "mediumblob", "longblob", "binary", "varbinary")


def _zero(data_type):
    """A value to write to a NOT NULL column that has no default, for a row the repair has to insert."""
    data_type = data_type.lower()
    if data_type in _NUMERIC:
        return 0
    if data_type in _BLOB:
        return b""
    if data_type in ("datetime", "timestamp"):
        return datetime.datetime(1970, 1, 2)
    if data_type == "date":
        return datetime.date(1970, 1, 2)
    if data_type == "time":
        return datetime.time(0, 0)
    if data_type == "json":
        return "{}"
    return ""


def _filler(conn, table, covered, sector_id, cache):
    """`[(column, value)]` for the columns of `table` the export does not carry that an inserted row must be given: its `sector_id`
    when the table files by one, and a zero for any other NOT NULL column with no default."""
    if table not in cache:
        cache[table] = conn.execute(
            "SELECT column_name AS name, is_nullable AS nullable, column_default AS dflt, extra AS x, data_type AS t"
            " FROM information_schema.columns WHERE table_schema = DATABASE() AND table_name = ? ORDER BY ordinal_position",
            (table,)).fetchall()
    fill = []
    for column in cache[table]:
        name = column["name"]
        extra = (column["x"] or "").upper()
        if name in covered or "GENERATED" in extra or "AUTO_INCREMENT" in extra:
            continue
        if name == "sector_id":
            fill.append((name, sector_id))
        elif column["nullable"] == "NO" and column["dflt"] is None:
            fill.append((name, _zero(column["t"])))
    return fill


def _key_of(row, key):
    return tuple(storage_checksum.encode_value(row[column]) for column in key)


def _values_of(row, columns):
    return tuple(storage_checksum.encode_value(row[column]) for column in columns)


def _insert_row(conn, table, covered, key, row, sector_id, cache):
    """Inserts one exported row: any row already holding its key is deleted first (it was filed elsewhere or is a stray copy),
    and the columns the export leaves out take `_filler`'s values."""
    fill = _filler(conn, table, covered, sector_id, cache)
    conn.execute(f"DELETE FROM `{table}` WHERE " + " AND ".join(f"`{c}` = ?" for c in key), tuple(row[c] for c in key))
    names = list(covered) + [name for name, _value in fill]
    conn.execute(f"INSERT INTO `{table}` ({', '.join(f'`{n}`' for n in names)}) VALUES ({', '.join('?' * len(names))})",
                 tuple(row[c] for c in covered) + tuple(value for _name, value in fill))


def apply_export(conn, sector_id, tables, expected):
    """
    Makes the rows of a sector the rows of an export (`sector_export.parse`), in one transaction, and checks the result hashes to
    `expected` before committing (design: this module's docstring). Returns `{"updated": n, "inserted": n, "deleted": n}`.

    Raises:
        RepairFailed: The rows could not be written, or do not hash to `expected` afterwards (everything is rolled back).
    """
    reader = storage_checksum.reader(conn)
    cache = {}
    done = {"updated": 0, "inserted": 0, "deleted": 0}
    try:
        locked = conn.execute("SELECT id FROM sectors WHERE id = ? FOR UPDATE", (sector_id,)).fetchall()
        if not locked:          # the sector's own row is gone (a lost page of `sectors`): put it back first, so the rest can be read by its id
            conn.execute("SET foreign_key_checks = 0")
            try:
                for table, covered, rows in tables:
                    if table == "sectors":
                        for row in rows:
                            _insert_row(conn, table, covered, storage_checksum.reader(conn).columns(table)[1], row, sector_id, cache)
                            done["inserted"] += 1
            finally:
                conn.execute("SET foreign_key_checks = 1")
        current = {table: rows for table, _columns, rows in storage_checksum.sector_rows(conn, sector_id)}
        conn.execute("SET foreign_key_checks = 0")
        try:
            for table, covered, rows in tables:
                _, key = reader.columns(table)
                have = {_key_of(row, key): row for row in current.get(table, [])}
                want = {_key_of(row, key): row for row in rows}
                body = [column for column in covered if column not in key]
                for ident, row in want.items():
                    held = have.get(ident)
                    if held is None:
                        _insert_row(conn, table, covered, key, row, sector_id, cache)
                        done["inserted"] += 1
                    elif body and _values_of(held, body) != _values_of(row, body):
                        conn.execute(f"UPDATE `{table}` SET " + ", ".join(f"`{c}` = ?" for c in body)
                                     + " WHERE " + " AND ".join(f"`{c}` = ?" for c in key),
                                     tuple(row[c] for c in body) + tuple(row[c] for c in key))
                        done["updated"] += 1
            for table, _covered, rows in reversed(tables):
                _, key = reader.columns(table)
                want = {_key_of(row, key) for row in rows}
                for row in current.get(table, []):
                    if _key_of(row, key) not in want:
                        conn.execute(f"DELETE FROM `{table}` WHERE " + " AND ".join(f"`{c}` = ?" for c in key),
                                     tuple(row[c] for c in key))
                        done["deleted"] += 1
        finally:
            conn.execute("SET foreign_key_checks = 1")
        value = storage_checksum.checksum(conn, sector_id)
        if value != expected:
            raise RepairFailed("the rows written do not hash to the checksum stored for the sector")
        storage_checksum.seal(conn, sector_id, value)
        conn.commit()
    except RepairFailed:
        conn.rollback()
        raise
    except READ_ERRORS as exc:
        conn.rollback()
        raise RepairFailed(f"{type(exc).__name__}: {exc}") from exc
    return done


# ---------------------------------------------------------------------
# The decision per sector and per group
# ---------------------------------------------------------------------

def _current_block(conn, member):
    """The export of a group member as it is stored, when it is still what the record holds for it; `None` when it is not."""
    try:
        data, value = sector_export.export_sector(conn, member.sector_id)
    except (storage_checksum.StorageChecksumError, *READ_ERRORS):
        conn.rollback()
        return None
    conn.commit()
    return data if value == member.checksum and len(data) == member.length else None


def _repair_group(conn, parity_store, group_id, sector_ids, dry_run):
    """Rebuilds the damaged `sector_ids` of one group at once. Returns their `Result`s."""
    results = []

    def refuse(ids, reason):
        results.extend(Result(sector_id, UNREPAIRABLE, reason) for sector_id in ids)
        return results

    with parity_store.file.locked():
        try:
            record = parity_store.file.read(group_id)
        except pfile.RecordDamaged:
            return refuse(sector_ids, "its group's parity record is damaged as well")
        code = parity._code(record.k, record.m)
        have, erased = {}, []
        for slot, member in enumerate(record.members):
            if member.sector_id == pfile.EMPTY:
                have[slot] = b""
                continue
            block = None if member.sector_id in sector_ids else _current_block(conn, member)
            if block is None:
                erased.append(slot)
            else:
                have[slot] = block
        if len(erased) > record.m:
            return refuse(sector_ids, f"{len(erased)} blocks of its group are damaged or have changed since the parity was made;"
                                      f" {record.m} can be mended")
        for index, block in enumerate(record.parity):
            have[record.k + index] = block
        rebuilt = code.decode(have)
        for sector_id in sector_ids:
            slot = [member.sector_id for member in record.members].index(sector_id)
            member = record.members[slot]
            block = rebuilt[slot][:member.length]
            try:
                tables = sector_export.parse(block)
            except sector_export.ExportError as exc:
                refuse([sector_id], f"the rebuilt block is not a sector export ({exc})")
                continue
            if storage_checksum.checksum_of(tables) != bytes(member.checksum):
                refuse([sector_id], "the rebuilt block does not hash to the stored checksum; the parity or another member is damaged")
                continue
            if dry_run:
                results.append(Result(sector_id, WOULD_REPAIR))
                continue
            try:
                done = apply_export(conn, sector_id, tables, bytes(member.checksum))
            except RepairFailed as exc:
                refuse([sector_id], f"the rebuilt rows could not be written: {exc}")
                continue
            results.append(Result(sector_id, REPAIRED, f"{done['updated']} row(s) updated, {done['inserted']} inserted, "
                                                       f"{done['deleted']} deleted"))
    return results


def _damaged_records(parity_store):
    """`(group ids whose latest record fails its CRC, the byte the file stops parsing at or None)`."""
    bad = []
    with parity_store.file.locked():
        for group_id in sorted(parity_store.file.entries):
            try:
                parity_store.file.read(group_id)
            except pfile.RecordDamaged:
                bad.append(group_id)
        return bad, parity_store.file.damaged_at


def repair(conn, scope=None, dry_run=False, say=log.normal):
    """
    Finds the damaged sectors in `scope` and repairs what the parity file allows (this module's docstring). With `dry_run`
    nothing is written: each sector that could be rebuilt is `WOULD_REPAIR`.

    Returns:
        Report: What became of each damaged sector and of the parity file.

    Raises:
        FileNotFoundError: There is no parity file to repair from.
    """
    path = parity.parity_path()
    if not os.path.exists(path):
        raise FileNotFoundError(f"no parity file at {path}; `planetgen refresh-parity` makes one (setting parity.enabled keeps it)")
    parity_store = parity.store()
    report = Report()
    damaged, unreadable, report.unprotected = find_damaged(conn, scope, say)
    report.results.extend(Result(sector_id, UNREADABLE, f"{why}; DB.48 rebuilds the table") for sector_id, why in unreadable)
    say(f"Repair: {len(damaged):,} damaged sector(s), {len(unreadable):,} unreadable, {report.unprotected:,} unprotected.")
    plan = {}
    for sector_id in damaged:
        with parity_store.file.locked():
            held = parity_store.member(sector_id)
        if held is None:
            report.results.append(Result(sector_id, UNREPAIRABLE, "it is in no group of the parity file"))
            continue
        group_id, _slot, member = held
        stored = bytes(conn.execute("SELECT storage_checksum AS c FROM sectors WHERE id = ?", (sector_id,)).fetchone()["c"])
        conn.commit()
        try:
            current = storage_checksum.checksum(conn, sector_id)
        except READ_ERRORS as exc:
            conn.rollback()
            report.results.append(Result(sector_id, UNREADABLE, f"{type(exc).__name__}: {exc}; DB.48 rebuilds the table"))
            continue
        conn.commit()
        if current == bytes(member.checksum):
            if not dry_run:
                storage_checksum.seal(conn, sector_id, current)
                conn.commit()
            report.results.append(Result(sector_id, WOULD_REPAIR if dry_run else RESEALED,
                                         "only its stored checksum was wrong" if dry_run else ""))
        elif bytes(member.checksum) == stored:
            plan.setdefault(group_id, []).append(sector_id)
        else:
            report.results.append(Result(sector_id, UNREPAIRABLE, "the parity record is behind the database as well as the "
                                                                    "rows being damaged; regenerate it from its seed or restore a backup"))
    for group_id, ids in sorted(plan.items()):
        report.results.extend(_repair_group(conn, parity_store, group_id, ids, dry_run))
    for result in report.results:
        say("Repair: " + result.line())
    bad, report.file_damaged_at = _damaged_records(parity_store)
    report.records_damaged = len(bad)
    if (bad or report.file_damaged_at is not None) and not dry_run:
        kept, dropped = parity.compact()
        counts = parity.refresh(conn, say=say)
        report.file_remade = True
        say(f"Repair: the parity file was compacted ({kept:,} kept, {dropped:,} unreadable dropped) and "
            + ", ".join(f"{count:,} {result}" for result, count in sorted(counts.items())) + ".")
    elif bad or report.file_damaged_at is not None:
        say(f"Repair: the parity file has {len(bad):,} damaged record(s)" +
            (f" and stops parsing at byte {report.file_damaged_at:,}" if report.file_damaged_at is not None else "") +
            "; a run without --dry-run makes it again from the database.")
    return report
