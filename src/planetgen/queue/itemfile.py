# planetgen/queue/itemfile.py

"""
The items of a windowed stage (OPS.240): what a stage with more units than the ledger should hold at once is made of, written once when the stage is cut and read a window at a time when the next
window of unit rows is written (`planrows.spawn_window`).

A sector stage is one unit for each sector, however many sectors the galaxy has; the cut writes the first window of unit rows and the sectors of the rest wait here, as fixed records in the order the
stage generates them. A stage that is cut again after a restart writes the same file again; a run that carries on after it was held reads the file its stage wrote (`Reader`), and a file that is
missing, cut short or different from what its header says fails loudly (`ItemFileError`), never as a stage that quietly has fewer units.

**The file.** A header (`HEADER`: the magic, the format version, the number of records, whether suffixes follow, and the CRC-32 of everything after the header), then `count` records of four 32-bit
integers (`RECORD`: the three numbers of the item's address and its weight, the expected work the stage's bar and ETA count), then, only when the stage has suffixes (the words a log line says after an
address), `count` 64-bit offsets (the end of each suffix in the blob that follows them) and the UTF-8 blob. It is written to a file of its own name and renamed into place, so a reader sees all of it or
none of it. It lives in the state folder (`settings.state_dir`), under `stage-items`, named by the run and the stage, and is deleted when the stage ends (`remove`).
"""

import os
import struct
import zlib

from planetgen.util import settings

MAGIC = b"PGSI"
"""bytes: What the first four bytes of an item file are."""

VERSION = 1
"""int: The format of the file; a reader of another version fails (no file of an older format is read)."""

HEADER = struct.Struct("<4sIQII")
"""struct.Struct: The header: the magic, the version, the number of records, the flags (`HAS_SUFFIXES`) and the CRC-32 of the rest of the file."""

RECORD = struct.Struct("<iiiI")
"""struct.Struct: One record: `shell, band, wedge` and the weight."""

HAS_SUFFIXES = 1
"""int: The flag that says a table of suffix offsets and the suffix blob follow the records."""

FOLDER = "stage-items"
"""str: The folder of the state folder the files are in."""

SUFFIX_OFFSET = struct.Struct("<Q")
"""struct.Struct: One entry of the suffix offsets."""


class ItemFileError(RuntimeError):
    """An item file is missing, torn or not the one its stage wrote."""


def path_of(run_id, stage):
    """The file of the stage `stage` of the run `run_id`."""
    return os.path.join(settings.state_dir(), FOLDER, f"{run_id}-{stage}.items")


def write(path, records):
    """
    Writes the item file `path`: `records` is a list of `(shell, band, wedge, weight, suffix)` in stage order (`suffix` is `""` for an item with none). The file is made whole under a name of its own and renamed, and
    the folder is made when it is not there; a file already there (a stage cut again) is replaced by the same content. Returns the number of records.
    """
    records = list(records)
    body = bytearray()
    for shell, band, wedge, weight, _suffix in records:
        body += RECORD.pack(shell, band, wedge, weight)
    flags = 0
    if any(record[4] for record in records):
        flags |= HAS_SUFFIXES
        blob, ends = bytearray(), bytearray()
        for record in records:
            blob += record[4].encode("utf-8")
            ends += SUFFIX_OFFSET.pack(len(blob))
        body += ends
        body += blob
    made = HEADER.pack(MAGIC, VERSION, len(records), flags, zlib.crc32(body)) + bytes(body)
    folder = os.path.dirname(path)
    os.makedirs(folder, exist_ok=True)
    scratch = f"{path}.part{os.getpid()}"
    try:
        with open(scratch, "wb") as handle:
            handle.write(made)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(scratch, path)
    finally:
        if os.path.exists(scratch):
            os.remove(scratch)
    return len(records)


def remove(path):
    """Deletes the file (none is not an error: the stage may have ended already)."""
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


class Reader:
    """
    An item file opened for reading windows. Opening it checks the header and that the file is as long as the header says; `verify` reads the whole of it against its CRC (a run that carries on or a spawner that starts
    again does so once).

    Attributes:
        count (int): The number of records.
        has_suffixes (bool): Whether the items have suffixes.
    """

    def __init__(self, path):
        self.path = path
        try:
            self._handle = open(path, "rb")
        except FileNotFoundError as exc:
            raise ItemFileError(f"the item file {path} of the stage is missing") from exc
        try:
            head = self._handle.read(HEADER.size)
            if len(head) < HEADER.size:
                raise ItemFileError(f"the item file {path} is cut short: {len(head)} bytes is not a header")
            magic, version, count, flags, crc = HEADER.unpack(head)
            if magic != MAGIC or version != VERSION:
                raise ItemFileError(f"the item file {path} is not a version {VERSION} file of this program (magic {magic!r}, version {version})")
            self.count, self.has_suffixes, self._crc = count, bool(flags & HAS_SUFFIXES), crc
            size = os.fstat(self._handle.fileno()).st_size
            table = count * SUFFIX_OFFSET.size if self.has_suffixes else 0
            floor = HEADER.size + count * RECORD.size + table
            if size < floor or (not self.has_suffixes and size != floor):
                raise ItemFileError(f"the item file {path} is {size:,} bytes, not the {floor:,} its header says: it is torn")
            self._blob_at = floor
            self._size = size
        except BaseException:
            self._handle.close()
            raise

    def close(self):
        self._handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def verify(self):
        """Reads everything after the header and compares its CRC-32 with the header's; a file that does not match fails."""
        self._handle.seek(HEADER.size)
        crc = 0
        while True:
            block = self._handle.read(1 << 22)
            if not block:
                break
            crc = zlib.crc32(block, crc)
        if crc != self._crc:
            raise ItemFileError(f"the item file {self.path} does not match its checksum: it is torn or was changed")
        if self.has_suffixes and self.count and self._suffix_end(self.count - 1) != self._size - self._blob_at:
            raise ItemFileError(f"the item file {self.path} has a suffix blob of the wrong length")

    def _suffix_end(self, index):
        self._handle.seek(HEADER.size + self.count * RECORD.size + index * SUFFIX_OFFSET.size)
        return SUFFIX_OFFSET.unpack(self._handle.read(SUFFIX_OFFSET.size))[0]

    def window(self, start, count):
        """The records `start` up to `start + count` (fewer at the end): `(shell, band, wedge, weight, suffix)` each."""
        if start < 0 or start > self.count:
            raise ItemFileError(f"the item file {self.path} has {self.count:,} records, not a record {start:,}")
        count = max(0, min(count, self.count - start))
        self._handle.seek(HEADER.size + start * RECORD.size)
        raw = self._handle.read(count * RECORD.size)
        if len(raw) != count * RECORD.size:
            raise ItemFileError(f"the item file {self.path} ended inside records {start:,} to {start + count:,}")
        rows = [RECORD.unpack_from(raw, place * RECORD.size) for place in range(count)]
        if not self.has_suffixes:
            return [(*row, "") for row in rows]
        first = 0 if start == 0 else self._suffix_end(start - 1)
        self._handle.seek(HEADER.size + self.count * RECORD.size + start * SUFFIX_OFFSET.size)
        ends = struct.unpack(f"<{count}Q", self._handle.read(count * SUFFIX_OFFSET.size)) if count else ()
        self._handle.seek(self._blob_at + first)
        blob = self._handle.read((ends[-1] - first) if count else 0)
        suffixes, place = [], 0
        for end in ends:
            suffixes.append(blob[place:end - first].decode("utf-8"))
            place = end - first
        return [(*row, suffix) for row, suffix in zip(rows, suffixes)]
