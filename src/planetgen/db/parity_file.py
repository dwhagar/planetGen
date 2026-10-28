# planetgen/db/parity_file.py

"""
The parity file: an append-only log of parity records, one record per group of sectors (DB.46;
docs/design/operations/db-check-and-parity-repair.md section 5.6). The groups and the codec are `db/parity.py` and
`util/rscodec.py`; this module is only the file.

**A record** is self-describing, little endian: a 32 byte header (`PGPR`, a format byte, the record's total length, the group
id, `seq`, `k` (member slots), `m` (parity blocks) and `P` (the length of a parity block)), `k` member slots of 44 bytes (the
sector id, 0 for an empty slot; the length of the member's export; its storage checksum), the `m` parity blocks of `P` bytes,
and a CRC32 of all that, which is there to tell a torn write from a whole record. The last record of a group is the group's
state (a later `seq` replaces an earlier one); a replaced record is dead space until `compact`.

**The index** (group id to its latest record's place and members, sector id to its group and slot) is rebuilt by a header scan
(`_scan`): it reads each record's header and slots and skips its parity bytes, and it checks the CRC of the last record only,
since only a tail can be torn. A record read to be used (`read`) has its CRC checked. A header that does not parse ends the
scan: a short tail is a torn write and nothing is lost (the next append overwrites it), anything longer is damage
(`damaged_at`), and `append` refuses to write behind it until `compact` rewrites the live records.

**Many writers.** Every operation runs inside `locked()`: a thread lock and an `flock` on the file, then a catch-up scan of what
other processes appended since this one last looked (and a reopen if `compact` replaced the file). So generation workers saving
different sectors at once see each other's groups and append one record at a time.
"""

import contextlib
import fcntl
import os
import struct
import threading
import zlib

MAGIC = b"PGPR"
FORMAT = 1
_HEADER = struct.Struct("<4sBIQQHBI")
_SLOT = struct.Struct("<QI32s")
HEADER_SIZE = _HEADER.size
SLOT_SIZE = _SLOT.size
_CRC = struct.Struct("<I")
EMPTY = 0
"""int: The sector id of an empty slot."""
MAX_RECORD_BYTES = 1 << 31
"""int: A record longer than this is not believed (a damaged length field)."""


class ParityFileError(OSError):
    """The file cannot be used (damaged behind the last record, or not writable)."""


class RecordDamaged(ParityFileError):
    """A record whose CRC does not match what was read."""


class Member:
    """One slot of a record: `sector_id` (`EMPTY` for a free slot), `length` of its export and its storage `checksum`."""

    __slots__ = ("sector_id", "length", "checksum")

    def __init__(self, sector_id=EMPTY, length=0, checksum=b"\0" * 32):
        self.sector_id = sector_id
        self.length = length
        self.checksum = bytes(checksum)

    def __eq__(self, other):
        return (isinstance(other, Member) and (self.sector_id, self.length, self.checksum)
                == (other.sector_id, other.length, other.checksum))

    def __repr__(self):
        return f"Member({self.sector_id}, {self.length}, {self.checksum[:4].hex()}..)"


class Record:
    """A group's state: `group_id`, `seq`, `k`, `m`, `members` (`k` of `Member`) and `parity` (`m` blocks of equal length)."""

    __slots__ = ("group_id", "seq", "k", "m", "members", "parity")

    def __init__(self, group_id, seq, k, m, members, parity):
        if len(members) != k or len(parity) != m or len({len(block) for block in parity}) > 1:
            raise ValueError("a record has k members and m parity blocks of one length")
        self.group_id = group_id
        self.seq = seq
        self.k = k
        self.m = m
        self.members = list(members)
        self.parity = [bytes(block) for block in parity]

    @property
    def parity_length(self):
        return len(self.parity[0]) if self.parity else 0

    def encode(self):
        length = self.parity_length
        total = HEADER_SIZE + self.k * SLOT_SIZE + self.m * length + _CRC.size
        body = [_HEADER.pack(MAGIC, FORMAT, total, self.group_id, self.seq, self.k, self.m, length)]
        body.extend(_SLOT.pack(member.sector_id, member.length, member.checksum) for member in self.members)
        body.extend(self.parity)
        data = b"".join(body)
        return data + _CRC.pack(zlib.crc32(data))


def _decode_header(raw):
    """`(total, group_id, seq, k, m, parity_length)` of a header, or `None` if it is not one."""
    if len(raw) < HEADER_SIZE:
        return None
    magic, fmt, total, group_id, seq, k, m, length = _HEADER.unpack(raw[:HEADER_SIZE])
    if magic != MAGIC or fmt != FORMAT or not k or not m:
        return None
    if total != HEADER_SIZE + k * SLOT_SIZE + m * length + _CRC.size or total > MAX_RECORD_BYTES:
        return None
    return total, group_id, seq, k, m, length


def decode(data):
    """A `Record` from its bytes.

    Raises:
        RecordDamaged: The bytes are not a whole record, or the CRC fails.
    """
    head = _decode_header(data)
    if head is None or head[0] != len(data):
        raise RecordDamaged("not a parity record")
    total, group_id, seq, k, m, length = head
    if _CRC.unpack(data[-_CRC.size:])[0] != zlib.crc32(data[:-_CRC.size]):
        raise RecordDamaged(f"the CRC of group {group_id} record {seq} does not match")
    at = HEADER_SIZE
    members = []
    for _ in range(k):
        sector_id, member_length, checksum = _SLOT.unpack(data[at:at + SLOT_SIZE])
        members.append(Member(sector_id, member_length, checksum))
        at += SLOT_SIZE
    parity = [data[at + j * length:at + (j + 1) * length] for j in range(m)]
    return Record(group_id, seq, k, m, members, parity)


class Entry:
    """What the index keeps of a group's latest record: where it is and what it says about its members."""

    __slots__ = ("offset", "size", "seq", "k", "m", "length", "members")

    def __init__(self, offset, size, seq, k, m, length, members):
        self.offset = offset
        self.size = size
        self.seq = seq
        self.k = k
        self.m = m
        self.length = length
        self.members = members


class ParityFile:
    """
    One parity file. Use `locked()` around every read or write; the methods that need it say so.

    Args:
        path (str): The file. It is created, with its folder, by the first `append`.
    """

    def __init__(self, path):
        self.path = path
        self.entries = {}
        """dict: Group id to the `Entry` of its latest record."""
        self.sector_group = {}
        """dict: Sector id to `(group id, slot)`."""
        self.valid_end = 0
        """int: The offset just past the last whole record."""
        self.damaged_at = None
        """int or None: Where the scan met a header that is not a torn tail (so records may be hidden behind it)."""
        self.records_seen = 0
        self.sync = True
        """bool: Whether `append` waits for the disk (`fsync`); tests that write thousands of records turn it off."""
        self.observer = None
        """object or None: Told of each change to the index, `changed(group_id, old entry or None, new entry)`, and of a
        reset, `reset()` (the group layer keeps its own view of the groups that way)."""
        self._thread_lock = threading.RLock()
        self._fd = None
        self._inode = None
        self._depth = 0

    # -- the lock and the index ------------------------------------------------------------------------------

    @contextlib.contextmanager
    def locked(self):
        """Holds the file against other threads and processes, with the index current (re-entrant)."""
        with self._thread_lock:
            if self._depth == 0:
                self._open_and_lock()
            self._depth += 1
            try:
                yield self
            finally:
                self._depth -= 1
                if self._depth == 0:
                    self._unlock()

    def _open_and_lock(self):
        os.makedirs(os.path.dirname(os.path.abspath(self.path)), exist_ok=True)
        while True:
            fd = os.open(self.path, os.O_RDWR | os.O_CREAT, 0o640)
            fcntl.flock(fd, fcntl.LOCK_EX)
            try:
                same = os.stat(self.path).st_ino == os.fstat(fd).st_ino
            except FileNotFoundError:
                same = False
            if same:
                break
            os.close(fd)       # `compact` replaced the file while this waited for the lock
        if self._inode != os.fstat(fd).st_ino:
            self._reset()
            self._inode = os.fstat(fd).st_ino
        self._fd = fd
        self._scan()

    def _unlock(self):
        fd, self._fd = self._fd, None
        if fd is not None:
            with contextlib.suppress(OSError):
                fcntl.flock(fd, fcntl.LOCK_UN)
            os.close(fd)

    def _reset(self):
        self.entries = {}
        self.sector_group = {}
        self.valid_end = 0
        self.damaged_at = None
        self.records_seen = 0
        if self.observer is not None:
            self.observer.reset()

    def _read_at(self, offset, count):
        return os.pread(self._fd, count, offset)

    def _scan(self):
        """Reads the records past `valid_end`, adding them to the index."""
        size = os.fstat(self._fd).st_size
        if size < self.valid_end:       # the file is shorter than what was seen: start again
            self._reset()
        offset = self.valid_end
        self.damaged_at = None
        while offset < size:
            head = _decode_header(self._read_at(offset, HEADER_SIZE))
            if head is None or offset + head[0] > size:
                self.damaged_at = None if size - offset < HEADER_SIZE or head is not None else offset
                break
            total, group_id, seq, k, m, length = head
            slots = self._read_at(offset + HEADER_SIZE, k * SLOT_SIZE)
            if offset + total == size:       # the last record: the only one a torn write can have hit
                try:
                    decode(self._read_at(offset, total))
                except RecordDamaged:
                    break
            members = tuple(Member(*_SLOT.unpack(slots[i * SLOT_SIZE:(i + 1) * SLOT_SIZE])) for i in range(k))
            self._index(Entry(offset, total, seq, k, m, length, members), group_id)
            offset += total
        self.valid_end = offset

    def _index(self, entry, group_id):
        old = self.entries.get(group_id)
        if old is not None and old.seq > entry.seq:
            return
        if old is not None:
            for member in old.members:
                if member.sector_id != EMPTY and self.sector_group.get(member.sector_id, (None,))[0] == group_id:
                    del self.sector_group[member.sector_id]
        self.entries[group_id] = entry
        for slot, member in enumerate(entry.members):
            if member.sector_id != EMPTY:
                self.sector_group[member.sector_id] = (group_id, slot)
        self.records_seen += 1
        if self.observer is not None:
            self.observer.changed(group_id, old, entry)

    # -- reading and writing (inside `locked()`) -------------------------------------------------------------

    def read(self, group_id):
        """The latest `Record` of a group, its CRC checked.

        Raises:
            KeyError: No such group.
            RecordDamaged: The record's bytes fail their CRC.
        """
        entry = self.entries[group_id]
        return decode(self._read_at(entry.offset, entry.size))

    def append(self, record):
        """Writes `record` as the latest state of its group and syncs it. Returns the offset it was written at.

        Raises:
            ParityFileError: The file is damaged behind its last whole record (run `compact`).
        """
        if self.damaged_at is not None:
            raise ParityFileError(f"{self.path} is damaged at byte {self.damaged_at}; compact it first")
        data = record.encode()
        if os.fstat(self._fd).st_size > self.valid_end:        # a torn tail: written over
            os.ftruncate(self._fd, self.valid_end)
        offset = self.valid_end
        written = 0
        while written < len(data):
            written += os.pwrite(self._fd, data[written:], offset + written)
        if self.sync:
            os.fsync(self._fd)
        head = _decode_header(data[:HEADER_SIZE])
        self._index(Entry(offset, len(data), record.seq, record.k, record.m, head[5], tuple(record.members)), record.group_id)
        self.valid_end = offset + len(data)
        return offset

    def compact(self):
        """Rewrites the file with only each group's latest record (and drops a group with no member left), keeping
        what is readable; a record that fails its CRC is dropped and counted. Returns `(kept, dropped)` groups."""
        kept = dropped = 0
        temp = self.path + ".compact"
        with open(temp, "wb") as out:
            for group_id in sorted(self.entries):
                entry = self.entries[group_id]
                if all(member.sector_id == EMPTY for member in entry.members):
                    continue
                data = self._read_at(entry.offset, entry.size)
                try:
                    decode(data)
                except RecordDamaged:
                    dropped += 1
                    continue
                out.write(data)
                kept += 1
            out.flush()
            os.fsync(out.fileno())
        os.replace(temp, self.path)
        self._unlock_after_replace()
        return kept, dropped

    def _unlock_after_replace(self):
        # The lock is on the old file's inode; the next `locked()` opens the new one and scans it.
        self._reset()
        self._inode = None

    def live_bytes(self):
        """Bytes of the latest records of the groups that still have a member."""
        return sum(entry.size for entry in self.entries.values()
                   if any(member.sector_id != EMPTY for member in entry.members))

    def size(self):
        return os.fstat(self._fd).st_size if self._fd is not None else (os.path.getsize(self.path) if os.path.exists(self.path) else 0)


_FILES = {}
_FILES_LOCK = threading.Lock()


def open_file(path):
    """The one `ParityFile` of `path` in this process."""
    path = os.path.abspath(path)
    with _FILES_LOCK:
        if path not in _FILES:
            _FILES[path] = ParityFile(path)
        return _FILES[path]
