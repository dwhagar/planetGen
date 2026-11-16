# planetgen/queue/sqltap.py

"""
The SQL tap of the independence harness (OPS.241): what tables a command touches, read from the statements it sends.

`tables_of(sql)` names the tables a statement reads and the tables it writes (a table in the `FROM` or `JOIN` of a `SELECT`, or of the subquery of a write, is a read; the target of an
`INSERT`, `REPLACE`, `UPDATE`, `DELETE`, `TRUNCATE`, `LOAD DATA`, `ALTER`, `DROP`, `CREATE`, `RENAME`, `OPTIMIZE` or `ANALYZE` is a write; a `LOCK TABLES` names a lock of each table;
`GET_LOCK('name')` is a lock of the named lock). Schema reads (`information_schema`, `SHOW`, `DESCRIBE`, `SET`) touch no table. `Tap` is the context manager that records, for the command that
is running (`Tap.running`), what its statements touch by wrapping the driver's cursors (`pymysql`), so every connection of the process is read, whoever opened it.

The declarations the tap is held to are in `queue/touches.py`.
"""

import contextlib
import re
import threading

_COMMENT = re.compile(r"/\*.*?\*/|--[^\n]*|#[^\n]*", re.S)
_STRING = re.compile(r"'(?:[^'\\]|\\.|'')*'|\"(?:[^\"\\]|\\.)*\"")
_NAME = r"`?([A-Za-z_][\w$]*)`?(?:\.`?([A-Za-z_][\w$]*)`?)?"
_WRITE_TARGET = re.compile(
    r"\b(?:INSERT\s+(?:IGNORE\s+)?INTO|REPLACE\s+INTO|INSERT\s+IGNORE|REPLACE|UPDATE\s+(?:IGNORE\s+)?|DELETE\s+(?:QUICK\s+)?(?:IGNORE\s+)?FROM|TRUNCATE\s+(?:TABLE\s+)?"
    r"|LOAD\s+DATA\s+(?:LOCAL\s+)?INFILE\s+\S+\s+(?:IGNORE\s+|REPLACE\s+)?INTO\s+TABLE|ALTER\s+(?:IGNORE\s+)?TABLE|DROP\s+TABLE\s+(?:IF\s+EXISTS\s+)?"
    r"|CREATE\s+(?:TEMPORARY\s+)?TABLE\s+(?:IF\s+NOT\s+EXISTS\s+)?|RENAME\s+TABLE|OPTIMIZE\s+(?:NO_WRITE_TO_BINLOG\s+)?TABLE|ANALYZE\s+(?:NO_WRITE_TO_BINLOG\s+)?TABLE"
    r"|CREATE\s+(?:UNIQUE\s+)?INDEX\s+\S+\s+ON|DROP\s+INDEX\s+\S+\s+ON)\s*" + _NAME, re.I)
_READ_SOURCE = re.compile(r"\b(?:FROM|JOIN)\s*" + _NAME, re.I)
_LOCK_TABLES = re.compile(r"\bLOCK\s+TABLES\b(.*)", re.I | re.S)
_GET_LOCK = re.compile(r"\bGET_LOCK\s*\(\s*'([^']*)'", re.I)
_SCHEMA = {"information_schema", "performance_schema", "mysql", "sys"}
_NOT_TABLES = {"dual", "select", "values", "lateral", "unnest"}
_TEMPORARY = re.compile(r"\bCREATE\s+TEMPORARY\s+TABLE\s+(?:IF\s+NOT\s+EXISTS\s+)?`?([A-Za-z_][\w$]*)`?", re.I)


_ON_DUPLICATE = re.compile(r"\bON\s+DUPLICATE\s+KEY\s+UPDATE\b.*", re.I | re.S)


def _clean(sql):
    """The statement with its comments removed and every string literal reduced to `''` (a table name in a string is not a table), except the name `GET_LOCK` is given."""
    locks = _GET_LOCK.findall(_COMMENT.sub(" ", sql))
    text = _ON_DUPLICATE.sub(" ", _STRING.sub("''", _COMMENT.sub(" ", sql)))      # `ON DUPLICATE KEY UPDATE col = ...` names columns, not a table
    return text, locks


def _table(match):
    schema, name = match.group(1), match.group(2)
    if name is not None:
        return None if schema.lower() in _SCHEMA else name
    return None if schema.lower() in _NOT_TABLES else schema


def tables_of(sql):
    """
    `(reads, writes, locks)` of one statement: three sorted tuples of names. A table that is both read and written by the statement is in both. A name that a `WITH` clause gives to a subquery is not
    a table and is left out.
    """
    text, locks = _clean(sql)
    ctes = {name.lower() for name in re.findall(r"\bWITH\s+(?:RECURSIVE\s+)?`?(\w+)`?\s+AS\s*\(|\)\s*,\s*`?(\w+)`?\s+AS\s*\(", text, re.I) for name in name if name}
    reads, writes, lockset = set(), set(), set(f"lock:{name}" for name in locks)
    head = text.lstrip().split(None, 1)[0].upper() if text.strip() else ""
    if head in ("SET", "SHOW", "DESCRIBE", "DESC", "EXPLAIN", "USE", "BEGIN", "COMMIT", "ROLLBACK", "START", "SAVEPOINT", "RELEASE", "KILL", "FLUSH"):
        return (), (), tuple(sorted(lockset))
    lock = _LOCK_TABLES.search(text)
    if lock:
        for item in lock.group(1).split(","):
            found = re.match(r"\s*" + _NAME, item)
            name = _table(found) if found else None
            if name:
                lockset.add(f"lock:{name}")
                (writes if re.search(r"\bWRITE\b", item, re.I) else reads).add(name)
        return tuple(sorted(reads)), tuple(sorted(writes)), tuple(sorted(lockset))
    for match in _WRITE_TARGET.finditer(text):
        name = _table(match)
        if name:
            writes.add(name)
    for match in _READ_SOURCE.finditer(text):
        name = _table(match)
        if name and name.lower() not in ctes:
            reads.add(name)
    # `INSERT ... SELECT` and `UPDATE ... JOIN` name a table once as a target and, in a join, as a source; a table in `FROM` of a `DELETE` is a target already counted
    if head == "DELETE":
        reads.discard(next(iter(writes), None))
    return tuple(sorted(reads)), tuple(sorted(writes)), tuple(sorted(lockset))


class Tap:
    """
    Records, per command, the tables the statements of this process touch while that command runs. `with Tap() as tap:` wraps `pymysql`'s cursors; `tap.running(kind)` is the context a command's
    unit (or cut, or end hook) runs in; `tap.touched[kind]` is `{"reads": set, "writes": set, "locks": set}`; a table the process creates `TEMPORARY` is left out (`tap.temporary`). A statement sent while no command is running is kept under `None`.
    """

    def __init__(self):
        self.touched = {}
        self.statements = []
        self.temporary = set()
        self._local = threading.local()
        self._originals = []

    def _stack(self):
        stack = getattr(self._local, "stack", None)
        if stack is None:
            stack = self._local.stack = []
        return stack

    @contextlib.contextmanager
    def running(self, kind):
        """The statements sent in this block (by this thread) are `kind`'s."""
        stack = self._stack()
        stack.append(kind)
        try:
            yield
        finally:
            stack.pop()

    def record(self, sql):
        """Adds one statement to the command that is running."""
        stack = self._stack()
        kind = stack[-1] if stack else None
        sql = sql if isinstance(sql, str) else sql.decode("utf-8", "replace")
        self.temporary.update(_TEMPORARY.findall(sql))      # a temporary table is one connection's own: no other command can touch it
        reads, writes, locks = tables_of(sql)
        reads, writes = [name for name in reads if name not in self.temporary], [name for name in writes if name not in self.temporary]
        slot = self.touched.setdefault(kind, {"reads": set(), "writes": set(), "locks": set()})
        slot["reads"].update(reads)
        slot["writes"].update(writes)
        slot["locks"].update(locks)

    def __enter__(self):
        import pymysql.cursors
        tap = self
        for cls in (pymysql.cursors.Cursor,):
            for method in ("execute", "executemany"):
                original = getattr(cls, method)

                def wrapped(self_, query, *args, _original=original, **kwargs):
                    if isinstance(query, (str, bytes)):
                        tap.record(query)
                    return _original(self_, query, *args, **kwargs)

                wrapped.__name__ = method
                self._originals.append((cls, method, original))
                setattr(cls, method, wrapped)
        return self

    def __exit__(self, *_exc):
        for cls, method, original in reversed(self._originals):
            setattr(cls, method, original)
        self._originals.clear()
        return False
