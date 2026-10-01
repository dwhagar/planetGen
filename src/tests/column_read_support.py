# tests/column_read_support.py

"""
Records which columns of which tables the code reads (TEST.11): while
`tracking_reads()` is active, every row `Connection` returns remembers the
keys asked of it, and a key counts as read from each table its statement
names. Handing a whole row on (`dict(row)`, `**row`, iterating it) reads
every column it has.
"""

import contextlib
import re

from stellarObjects import _db


class _TrackedRow(dict):
    def __init__(self, row, used):
        super().__init__(row)
        self._used = used

    def __getitem__(self, key):
        self._used.add(key)
        return super().__getitem__(key)

    def get(self, key, default=None):
        self._used.add(key)
        return super().get(key, default)

    def __contains__(self, key):
        self._used.add(key)
        return super().__contains__(key)

    def _all(self):
        self._used.update(super().keys())

    def __iter__(self):
        self._all()
        return super().__iter__()

    def keys(self):
        self._all()
        return super().keys()

    def items(self):
        self._all()
        return super().items()

    def values(self):
        self._all()
        return super().values()

    def copy(self):
        self._all()
        return dict(super().items())


class _TrackedCursor:
    def __init__(self, cursor, used):
        self._cursor = cursor
        self._used = used

    def _wrap(self, row):
        return _TrackedRow(row, self._used) if isinstance(row, dict) else row

    def fetchone(self):
        return self._wrap(self._cursor.fetchone())

    def fetchall(self):
        return [self._wrap(row) for row in self._cursor.fetchall()]

    def fetchmany(self, *args):
        return [self._wrap(row) for row in self._cursor.fetchmany(*args)]

    def __iter__(self):
        return iter(self.fetchall())

    def __getattr__(self, name):
        return getattr(self._cursor, name)


@contextlib.contextmanager
def tracking_reads(monkeypatch):
    """Yields `reads`: `{sql: set of keys read from its rows}`."""
    reads = {}
    real_run = _db.Connection._run

    def run(self, sql, params):
        return _TrackedCursor(real_run(self, sql, params), reads.setdefault(sql, set()))

    monkeypatch.setattr(_db.Connection, "_run", run)
    yield reads


def columns_read(reads, tables):
    """`{table: columns read}` for `tables`, from `tracking_reads`' record."""
    found = {table: set() for table in tables}
    for sql, keys in reads.items():
        if not keys:
            continue
        for table in tables:
            if re.search(rf"\b(?:FROM|JOIN)\s+`?{table}`?\b", sql, re.IGNORECASE):
                found[table] |= keys
    return found
