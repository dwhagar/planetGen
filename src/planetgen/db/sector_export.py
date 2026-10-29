# planetgen/db/sector_export.py

"""
A sector as bytes, and back (DB.46; docs/design/operations/db-check-and-parity-repair.md sections 5.2, 5.3 and 5.6).

The parity file (`db/parity.py`) protects each sector as one block of bytes, its **export**: the rows the storage checksum
covers (`storage_checksum.sector_rows`: the sector, everything filed in it and the child rows, without the columns the
orbit update rewrites), every value written exactly as the checksum writes it, so the export of a sector that has not
changed is the same bytes every time and `storage_checksum.checksum_of(parse(export))` is the checksum stored for it. The
export is also re-importable: it names its tables and columns and keeps every id, so a repair (DB.47) can read a rebuilt
block back into rows and insert them with their original ids.

Layout, little endian: `PGSX`, a format byte, then per table: its name, its column names, its row count and each row's
values; a name is a u16 length and its ASCII bytes, a count a u32, a value a u32 length and `storage_checksum.encode_value`'s
bytes (a tag letter, then the value exactly: an integer's digits, a float's 64 bits as hex, a string's UTF-8).
"""

import datetime
import struct
from decimal import Decimal

from planetgen.db import storage_checksum

MAGIC = b"PGSX"
FORMAT = 1


class ExportError(ValueError):
    """Bytes that are not a sector export (a torn or damaged block)."""


def _name(text):
    raw = text.encode("ascii")
    return struct.pack("<H", len(raw)) + raw


def dump(tables):
    """The export of `tables`, `storage_checksum.sector_rows`'s list `[(table, covered columns, [row, ...])]`."""
    out = [MAGIC, bytes([FORMAT])]
    for table, covered, rows in tables:
        out.append(_name(table))
        out.append(struct.pack("<H", len(covered)))
        out.extend(_name(column) for column in covered)
        out.append(struct.pack("<I", len(rows)))
        for row in rows:
            for column in covered:
                value = storage_checksum.encode_value(row[column])
                out.append(struct.pack("<I", len(value)))
                out.append(value)
    return b"".join(out)


def export_sector(conn, sector_id):
    """
    The export of a sector as it is stored now, with its storage checksum, both from one read of its rows.

    Returns:
        tuple: `(export bytes, checksum bytes)`.

    Raises:
        storage_checksum.StorageChecksumError: The sector has no row.
    """
    tables = storage_checksum.sector_rows(conn, sector_id)
    return dump(tables), storage_checksum.checksum_of(tables)


def decode_value(raw):
    """One value back from `storage_checksum.encode_value`'s bytes (a tag letter, then the value exactly)."""
    tag, body = raw[:1], raw[1:]
    if tag == b"N":
        return None
    text = body.decode("utf-8") if tag in (b"I", b"D", b"S", b"T", b"A", b"M", b"L") else None
    if tag == b"I":
        return int(text)
    if tag == b"F":
        return struct.unpack("<d", bytes.fromhex(body.decode("ascii")))[0]
    if tag == b"D":
        return Decimal(text)
    if tag == b"B":
        return bytes.fromhex(body.decode("ascii"))
    if tag == b"S":
        return text
    if tag == b"T":
        return datetime.datetime.fromisoformat(text)
    if tag == b"A":
        return datetime.date.fromisoformat(text)
    if tag == b"M":
        return datetime.time.fromisoformat(text)
    if tag == b"L":
        return datetime.timedelta(seconds=float(text))
    raise ExportError(f"unknown value tag {tag!r}")


class _Cursor:
    def __init__(self, data):
        self.data = memoryview(data)
        self.at = 0

    def take(self, count):
        if self.at + count > len(self.data):
            raise ExportError("the export ends inside a value")
        chunk = bytes(self.data[self.at:self.at + count])
        self.at += count
        return chunk

    def number(self, fmt):
        return struct.unpack(fmt, self.take(struct.calcsize(fmt)))[0]

    def name(self):
        return self.take(self.number("<H")).decode("ascii")


def parse(data):
    """
    The tables of an export, as `storage_checksum.sector_rows` gives them: `[(table, columns, [row, ...])]`, a row a dict
    from column name to value (so `storage_checksum.checksum_of` reads it).

    Raises:
        ExportError: The bytes are not an export, or end early.
    """
    cursor = _Cursor(data)
    if cursor.take(len(MAGIC)) != MAGIC:
        raise ExportError("not a sector export")
    if cursor.take(1)[0] != FORMAT:
        raise ExportError("an export of another format")
    tables = []
    while cursor.at < len(cursor.data):
        table = cursor.name()
        columns = [cursor.name() for _ in range(cursor.number("<H"))]
        rows = []
        for _ in range(cursor.number("<I")):
            rows.append({column: decode_value(cursor.take(cursor.number("<I"))) for column in columns})
        tables.append((table, columns, rows))
    return tables
