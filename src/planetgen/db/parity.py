# planetgen/db/parity.py

"""
Parity over groups of sectors (DB.46; docs/design/operations/db-check-and-parity-repair.md sections 5.1 to 5.7). A parity file
(`db/parity_file.py`) outside the database holds, for each group of up to `G` sectors, `m` Reed-Solomon parity blocks
(`util/rscodec.py`) over the sectors' exports (`db/sector_export.py`), so that up to `m` damaged sectors of a group can be
rebuilt from the others (DB.47). This module forms the groups, keeps the parity current as sectors are saved, edited and
deleted, and holds the decision table that tells a stale record from damage.

**Groups** are made at save time. A sector joins an open group of its **size class** (`ceil(log(size / 2000) / log 1.25)`,
about 20 classes for exports of 2 to 200 KB, so a group's blocks are about as long as each other and the padding is small)
that has no member within `W` sector ids (a damaged stretch of a table file takes neighbouring rows, so the members of a
group are spread out; `W` is 256, less for a small galaxy: about the sector count over 100). Of the groups that qualify the
fullest wins; with none it opens a new one. An export is zero-padded to the longest member of its group, and each record keeps
every member's true length, so a rebuilt block is cut back.

**Updates.** A new member of a group has no before-image (it was zeros), so joining reads no other sector. A member that
changes is updated from its before-image (`P_j' = P_j xor c_ji x (old xor new)`, `RS.update`) when the caller has one whose
checksum is the one the record holds, and otherwise the group's parity is made again from the current rows of all its
members (`G - 1` reads, only an edit or a deletion pays it). A member whose rows no longer match the checksum stored for it
is damage, not a change: nothing is recomputed over it (that would bake the damage into the parity) and the group is left for
the repair.

**Order and crashes.** The writers only mark a sector (`storage_checksum.mark_*`: `parity_dirty = 1`) and commit. After the
commit `protect_sector` takes the sector's row lock, exports it, writes the parity record and then seals the sector (stores
the checksum, clears the flag). A crash before the record leaves the flag set; between the record and the seal it leaves the
flag set with a current record; either way the next `protect_sector` finishes it. The record's member checksum against the
database's is the decision table (`classify`): a record behind the database is stale, not damage.

The connection reports the sectors a transaction marked (`Connection.parity_touch`) and `Connection.after_commit` runs
`protect_sectors` for them once the transaction has committed, so every save, edit and deletion keeps the parity current
without a writer calling it. `refresh` makes a whole database current: sectors that were never sealed, sealed sectors the
file does not hold, and members of the file whose sector is gone. Nothing here runs unless `parity.enabled` is on.
"""

import math

from planetgen.db import parity_file as pfile, sector_export, storage_checksum
from planetgen.util import log, rscodec
from planetgen.util.settings import get_settings

DEFAULT_DIR = "/var/lib/planetGen/parity"
"""str: The folder of the parity file when `parity.dir` is empty."""
FILE_NAME = "sectors.parity"

SIZE_CLASS_BASE = 2000
SIZE_CLASS_GROWTH = 1.25
WINDOW_MAX = 256
"""int: The widest spread rule: a new member is at least this many sector ids from every member of its group."""
_BUCKET = WINDOW_MAX


class MemberUnreadable(Exception):
    """A group's member whose rows no longer match its stored checksum, or that cannot be read: the group cannot be made
    again from the database. `sector_id` names it."""

    def __init__(self, sector_id, reason):
        super().__init__(f"sector {sector_id}: {reason}")
        self.sector_id = sector_id


class NeedBefore(Exception):
    """An existing member was changed with neither a before-image nor a way to read the other members."""


def size_class(length):
    """The size class of an export of `length` bytes (0 for 2,000 or fewer)."""
    if length <= SIZE_CLASS_BASE:
        return 0
    return math.ceil(math.log(length / SIZE_CLASS_BASE) / math.log(SIZE_CLASS_GROWTH))


def window(members):
    """`W`: how far apart in sector ids the members of a group must be, for a file holding `members` sectors."""
    return max(1, min(WINDOW_MAX, members // 100))


class GroupIndex:
    """
    The group layer's view of the file (kept current by `ParityFile.observer`): per size class the groups with a free slot,
    and per class and range of ids the groups with a member there, so the groups a new sector may join are found without
    looking at every group.
    """

    def __init__(self):
        self.reset()

    def reset(self):
        self.open = {}
        self.blocks = {}
        self.info = {}
        self.max_group = 0

    def changed(self, group_id, old, new):
        if old is not None:
            self._drop(group_id)
        self.max_group = max(self.max_group, group_id)
        ids = [member.sector_id for member in new.members if member.sector_id != pfile.EMPTY]
        lengths = [member.length for member in new.members if member.sector_id != pfile.EMPTY]
        cls = size_class(max(lengths)) if lengths else 0
        self.info[group_id] = (cls, ids, new.k)
        if len(ids) < new.k:
            self.open.setdefault(cls, set()).add(group_id)
        for sector_id in ids:
            self.blocks.setdefault((cls, sector_id // _BUCKET), set()).add(group_id)

    def _drop(self, group_id):
        cls, ids, _k = self.info.pop(group_id)
        self.open.get(cls, set()).discard(group_id)
        for sector_id in ids:
            self.blocks.get((cls, sector_id // _BUCKET), set()).discard(group_id)

    def choose(self, sector_id, length, total_members):
        """The open group of the right size class, with no member within `W` ids of `sector_id`, that has the most members
        (the lowest id on a tie); `None` when there is none."""
        cls = size_class(length)
        candidates = self.open.get(cls)
        if not candidates:
            return None
        width = window(total_members)
        blocked = set()
        bucket = sector_id // _BUCKET
        for near in (bucket - 1, bucket, bucket + 1):
            for group_id in self.blocks.get((cls, near), ()):
                if group_id in candidates and any(abs(other - sector_id) < width for other in self.info[group_id][1]):
                    blocked.add(group_id)
        free = [group_id for group_id in candidates if group_id not in blocked]
        return max(free, key=lambda group_id: (len(self.info[group_id][1]), -group_id), default=None)


_CODES = {}


def _code(k, m):
    if (k, m) not in _CODES:
        _CODES[(k, m)] = rscodec.RS(k, m)
    return _CODES[(k, m)]


class ParityStore:
    """
    The groups of one parity file, without the database: `put` makes a sector's block the given bytes and `release` takes
    it out, each as one appended record.

    Args:
        path (str): The parity file.
        group_size (int): `G`, the member slots of a new group.
        parity_blocks (int): `m`, the parity blocks of a new group.
    """

    def __init__(self, path, group_size=32, parity_blocks=2):
        self.file = pfile.open_file(path)
        self.k = group_size
        self.m = parity_blocks
        _code(group_size, parity_blocks)        # refuses a code that cannot exist
        with self.file.locked():
            if not isinstance(self.file.observer, GroupIndex):
                self.file.observer = GroupIndex()
                self.file.observer.reset()
                for group_id, entry in self.file.entries.items():
                    self.file.observer.changed(group_id, None, entry)
        self.index = self.file.observer

    # -- reading ---------------------------------------------------------------------------------------------

    def member(self, sector_id):
        """`(group id, slot, Member)` of a sector in the file, or `None`. Call it inside `file.locked()`."""
        located = self.file.sector_group.get(sector_id)
        if located is None:
            return None
        return located[0], located[1], self.file.entries[located[0]].members[located[1]]

    def members(self):
        """The sector ids the file holds."""
        with self.file.locked():
            return set(self.file.sector_group)

    # -- writing ---------------------------------------------------------------------------------------------

    def put(self, sector_id, data, checksum, before=None, read_member=None):
        """
        Makes `data` (a sector's export, with its storage `checksum`) the block of `sector_id`.

        Args:
            before (tuple, optional): `(old export, its checksum)` of a sector that is already in the file; used for the
                delta only when the checksum is the one the record holds.
            read_member (callable, optional): `read_member(sector_id)` returns `(export, checksum)` of another member as it
                is stored now (or `(b"", None)` for a sector that is gone), and raises `MemberUnreadable` for one that
                cannot be trusted. Used to make a group's parity again when there is no usable before-image.

        Returns:
            str: `"joined"` (a new member), `"delta"` or `"recomputed"` (an existing one).

        Raises:
            NeedBefore: The sector is in the file and neither `before` nor `read_member` can bring its old block.
            MemberUnreadable: Raised by `read_member`.
            pfile.RecordDamaged: The group's record fails its CRC.
        """
        with self.file.locked():
            located = self.member(sector_id)
            if located is None:
                return self._join(sector_id, data, checksum)
            group_id, slot, held = located
            record = self.file.read(group_id)
            code = _code(record.k, record.m)
            if before is not None and before[1] == held.checksum and len(before[0]) == held.length:
                parity = code.update(record.parity, slot, before[0], data)
                record.members[slot] = pfile.Member(sector_id, len(data), checksum)
                mode = "delta"
            elif read_member is not None:
                blocks = []
                for index, other in enumerate(record.members):
                    if index == slot:
                        blocks.append(data)
                    elif other.sector_id == pfile.EMPTY:
                        blocks.append(b"")
                    else:
                        export, other_sum = read_member(other.sector_id)
                        blocks.append(export)
                        record.members[index] = (pfile.Member(other.sector_id, len(export), other_sum) if export
                                                 else pfile.Member())
                record.members[slot] = pfile.Member(sector_id, len(data), checksum)
                parity = code.encode(blocks)
                mode = "recomputed"
            else:
                raise NeedBefore(f"sector {sector_id} is in the file; its old block is needed")
            self._append(record, parity)
            return mode

    def release(self, sector_id, before=None, read_member=None):
        """
        Takes a sector out of its group (a deleted sector, a sector that moved). Same arguments as `put`'s; with neither
        the sector is only forgotten and the group's parity is no longer right until the group is made again.

        Returns:
            str: `"released"`, `"delta"`, `"recomputed"`, or `"absent"` (it was not in the file).
        """
        with self.file.locked():
            located = self.member(sector_id)
            if located is None:
                return "absent"
            group_id, slot, held = located
            record = self.file.read(group_id)
            code = _code(record.k, record.m)
            record.members[slot] = pfile.Member()
            if before is not None and before[1] == held.checksum and len(before[0]) == held.length:
                self._append(record, code.update(record.parity, slot, before[0], b""))
                return "delta"
            if read_member is None:
                raise NeedBefore(f"sector {sector_id} is in the file; its old block is needed")
            blocks = []
            for index, other in enumerate(record.members):
                if other.sector_id == pfile.EMPTY:
                    blocks.append(b"")
                    continue
                export, other_sum = read_member(other.sector_id)
                blocks.append(export)
                record.members[index] = pfile.Member(other.sector_id, len(export), other_sum) if export else pfile.Member()
            self._append(record, code.encode(blocks))
            return "recomputed"

    def remake(self, group_id, read_member):
        """Makes a group's parity again from the current rows of its members (a stale record)."""
        with self.file.locked():
            record = self.file.read(group_id)
            code = _code(record.k, record.m)
            blocks = []
            for index, other in enumerate(record.members):
                if other.sector_id == pfile.EMPTY:
                    blocks.append(b"")
                    continue
                export, other_sum = read_member(other.sector_id)
                blocks.append(export)
                record.members[index] = pfile.Member(other.sector_id, len(export), other_sum) if export else pfile.Member()
            self._append(record, code.encode(blocks))

    def _join(self, sector_id, data, checksum):
        total = len(self.file.sector_group)
        group_id = self.index.choose(sector_id, len(data), total)
        if group_id is None:
            record = pfile.Record(self.index.max_group + 1, 0, self.k, self.m, [pfile.Member()] * self.k, [b""] * self.m)
        else:
            record = self.file.read(group_id)
        slot = next(index for index, member in enumerate(record.members) if member.sector_id == pfile.EMPTY)
        code = _code(record.k, record.m)
        parity = code.update(record.parity, slot, b"", data)
        record.members[slot] = pfile.Member(sector_id, len(data), checksum)
        self._append(record, parity)
        return "joined"

    def _append(self, record, parity):
        """Writes `record` with `parity`, cut to the longest member (the rest is zeros), as the group's next state."""
        longest = max((member.length for member in record.members if member.sector_id != pfile.EMPTY), default=0)
        record.parity = [bytes(block[:longest]).ljust(longest, b"\0") for block in parity]
        record.seq += 1
        self.file.append(record)


# ---------------------------------------------------------------------
# The database side
# ---------------------------------------------------------------------

def enabled():
    """Whether parity is kept (`parity.enabled`)."""
    return get_settings().parity.enabled


_STORES = {}


def parity_path(settings=None):
    import os
    settings = settings or get_settings().parity
    return os.path.join(settings.dir or DEFAULT_DIR, FILE_NAME)


def store(settings=None):
    """The `ParityStore` the settings name, one per file in this process."""
    settings = settings or get_settings().parity
    key = (parity_path(settings), settings.group_size, settings.parity_blocks)
    if key not in _STORES:
        _STORES[key] = ParityStore(key[0], settings.group_size, settings.parity_blocks)
    return _STORES[key]


HEALTHY = "healthy"
RECORD_STALE = "record_stale"
REBUILDABLE = "rebuildable"
LOST = "lost"
PARITY_DAMAGED = "parity_damaged"
NOT_PROTECTED = "not_protected"
"""str: The states `classify` returns (design 5.6): everything matches; the record is behind the database (recompute it, it
is not damage); the rows are damaged and the record is current (rebuild from parity); damaged and the record is behind
(parity cannot be trusted for it: seed regeneration or a backup); the record's CRC fails (recompute it from the database if
every member verifies); the sector is in no group yet."""


def classify(rows_match, record_checksum, stored_checksum, record_ok=True):
    """
    The stale-record decision table of design 5.6 as code.

    Args:
        rows_match (bool): The sector's rows still hash to the checksum stored for it.
        record_checksum (bytes or None): The checksum the parity record holds for it, `None` when it is in no group.
        stored_checksum (bytes): The checksum stored in the database.
        record_ok (bool): The group's record passes its CRC.

    Returns:
        str: One of `HEALTHY`, `RECORD_STALE`, `REBUILDABLE`, `LOST`, `PARITY_DAMAGED`, `NOT_PROTECTED`.
    """
    if record_checksum is None:
        return NOT_PROTECTED
    if not record_ok:
        return PARITY_DAMAGED
    current = bytes(record_checksum) == bytes(stored_checksum)
    if rows_match:
        return HEALTHY if current else RECORD_STALE
    return REBUILDABLE if current else LOST


def member_state(conn, sector_id, parity_store=None):
    """
    `classify`'s answer for one sector of the database (reads its rows; changes nothing). A sector that is `unprotected`
    (flag set, or no checksum yet) is not damaged and is reported as `NOT_PROTECTED` when the file lacks it and
    `RECORD_STALE` when the file has it.
    """
    parity_store = parity_store or store()
    row = conn.execute("SELECT storage_checksum AS c, parity_dirty AS d FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if row is None:
        raise storage_checksum.StorageChecksumError(f"sector {sector_id} has no row")
    with parity_store.file.locked():
        held = parity_store.member(sector_id)
        record_ok = True
        if held is not None:
            try:
                parity_store.file.read(held[0])
            except pfile.RecordDamaged:
                record_ok = False
    record_sum = None if held is None else held[2].checksum
    if row["c"] is None or row["d"]:
        return NOT_PROTECTED if held is None else RECORD_STALE
    stored = bytes(row["c"])
    return classify(storage_checksum.checksum(conn, sector_id) == stored, record_sum, stored, record_ok)


def _read_member(conn):
    def read(other_id):
        row = conn.execute("SELECT parity_dirty AS d, storage_checksum AS c FROM sectors WHERE id = ?", (other_id,)).fetchone()
        if row is None:
            return b"", None
        export, checksum = sector_export.export_sector(conn, other_id)
        if not row["d"] and row["c"] is not None and bytes(row["c"]) != checksum:
            raise MemberUnreadable(other_id, "its rows do not match the checksum stored for it")
        return export, checksum
    return read


CURRENT = "current"
PROTECTED = "protected"
RELEASED = "released"
UNPROTECTABLE = "unprotectable"


def protect_sector(conn, sector_id, parity_store=None, before=None):
    """
    Makes the parity record of one sector current and seals it, in a transaction of its own (call it with no transaction open,
    after the writer committed). The sector's row is locked, so a writer marking it waits and the sector is marked again.

    Args:
        before (tuple, optional): `(old export, its checksum)`, as `ParityStore.put`.

    Returns:
        str: `"protected"`, `"current"` (the record already holds these rows; only the flag is cleared), `"released"` (the
        sector is gone and was taken out of its group), or `"unprotectable"` (a member of the group is damaged or its record
        is unreadable: nothing was changed, the sector stays marked).
    """
    parity_store = parity_store or store()
    row = conn.execute("SELECT id FROM sectors WHERE id = ? FOR UPDATE", (sector_id,)).fetchone()
    try:
        if row is None:
            result = _release(conn, sector_id, parity_store)
            conn.commit()
            return result
        data, checksum = sector_export.export_sector(conn, sector_id)
        with parity_store.file.locked():
            held = parity_store.member(sector_id)
            if held is not None and held[2].checksum == checksum and held[2].length == len(data):
                result = CURRENT
            else:
                try:
                    parity_store.put(sector_id, data, checksum, before=before, read_member=_read_member(conn))
                except (MemberUnreadable, pfile.RecordDamaged) as exc:
                    conn.rollback()
                    log.error(f"Parity: sector {sector_id} left unprotected: {exc}")
                    return UNPROTECTABLE
                result = PROTECTED
        storage_checksum.seal(conn, sector_id, checksum)
        conn.commit()
        return result
    except Exception:
        conn.rollback()
        raise


def _release(conn, sector_id, parity_store):
    try:
        parity_store.release(sector_id, read_member=_read_member(conn))
    except (MemberUnreadable, pfile.RecordDamaged) as exc:
        log.error(f"Parity: sector {sector_id} is gone but its group could not be made again: {exc}")
        return UNPROTECTABLE
    return RELEASED


def protect_sectors(config, sector_ids):
    """`protect_sector` for each of `sector_ids`, on a connection of its own; a failure is logged and never raised (the sector
    stays marked and a later `refresh` finishes it). Returns `{result: count}`."""
    from planetgen.db import store as dbStore
    counts = {}
    conn = dbStore.get_connection(config)
    try:
        parity_store = store()
        for sector_id in sorted(set(sector_ids)):
            try:
                result = protect_sector(conn, sector_id, parity_store)
            except Exception as exc:  # noqa: BLE001 -- the writer has committed; the parity catches up at the next refresh
                log.error(f"Parity: sector {sector_id} could not be protected: {type(exc).__name__}: {exc}")
                result = "failed"
            counts[result] = counts.get(result, 0) + 1
    finally:
        conn.close()
    return counts


def compact():
    """Rewrites the parity file with only the latest record of each group. Returns `(kept, dropped)` groups."""
    parity_store = store()
    with parity_store.file.locked():
        return parity_store.file.compact()


def refresh(conn, say=log.normal, limit=None):
    """
    Makes the parity of the whole database current: every sector that is marked or was never sealed, every sealed sector the
    file does not hold, and every member of the file whose sector is gone. Returns `{result: count}`.
    """
    parity_store = store()
    counts = {}

    def note(result):
        counts[result] = counts.get(result, 0) + 1

    ids = [row["id"] for row in conn.execute("SELECT id FROM sectors ORDER BY id").fetchall()]
    conn.commit()
    present = set(ids)
    held = parity_store.members()
    todo = []
    for row in conn.execute("SELECT id, parity_dirty AS d, storage_checksum AS c FROM sectors ORDER BY id").fetchall():
        if row["d"] or row["c"] is None or row["id"] not in held:
            todo.append(row["id"])
    conn.commit()
    todo.extend(sorted(held - present))
    if limit is not None:
        todo = todo[:limit]
    say(f"Parity: {len(todo):,} sector(s) to make current.")
    for done, sector_id in enumerate(todo, 1):
        note(protect_sector(conn, sector_id, parity_store))
        if done % 1000 == 0:
            say(f"Parity: {done:,} of {len(todo):,}.")
    return counts
