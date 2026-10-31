# planetgen/db/migrations/helpers.py

"""
Helpers for revisions that touch a big table (DB.18, docs/design/operations/db-check-and-parity-repair.md section 7).

Alembic's `op.add_column`, `op.create_index` and `op.alter_column` cannot name an `ALGORITHM`, and without one
the server may copy the whole table and block writers for the length of the copy; `op.batch_alter_table` is a
table copy built for SQLite and is never right here. These helpers write the statement out with `ALGORITHM=` (and
`LOCK=NONE`) so that the server refuses when it cannot do the work online, instead of copying silently. Every
helper skips work that is already done, so a revision that was stopped half way can be started again (the
`exists_*` checks), and works on the SQLAlchemy connection a revision gets from `op.get_bind()`.

    from planetgen.db.migrations import helpers

    def upgrade():
        connection = op.get_bind()
        helpers.add_column_instant(connection, "sectors", "storage_checksum", "BIGINT UNSIGNED NULL")
        helpers.add_index_online(connection, "stars", "idx_stars_type", ("star_type",))
        helpers.batched_update(connection, "stars", "mass_solar = mass_solar * 1", "mass_solar IS NULL",
                               label="Star masses")

Measured costs (a 5 million row table, one thread writing meanwhile): an instant column or a dropped one takes a
fraction of a second; an online index takes a minute and writers carry on; `ALGORITHM=COPY` or a new CHECK
blocks writes for the whole run. Scaled to 1.6 x 10^8 rows an index is 11 to 31 minutes.

`sectors` is the exception to "instant" (DB.84): its FULLTEXT index on `name` (`ft_sectors_name`) makes MySQL refuse
`ALGORITHM=INSTANT` and rebuild the whole table, so a revision that changes a `sectors` column pays a full copy of the
table and needs free disk for a second copy. The rule, decided with Boss's defaults for the DB items: the FULLTEXT index
stays; every remaining change to `sectors` is batched into as few revisions as possible (one revision for all the
columns a version adds); and such a revision says in its docstring that it rebuilds `sectors` and why, so the cost is
read before the revision is run. Moving the index to a side table (so `sectors` can use INSTANT again) is the other
choice and is not taken.
"""

import contextlib
import re
import time

import sqlalchemy as sa
from sqlalchemy import exc

from planetgen.db import alembic_runner

BATCH_ROWS = 20000
"""int: Rows an update step changes before it commits (about 1.5 s a batch measured at 50,000)."""

LOCK_RETRIES = 3
"""int: How often a statement that timed out waiting for a metadata lock is tried again."""

LONG_TRANSACTION_S = 60
"""int: A transaction open this long makes `require_no_long_transactions` refuse: an online ALTER needs the table's
metadata lock briefly at its start and end, and one long reader leaves everything behind it waiting."""

_LOCK_ERRORS = (1205, 1213)
"""1205 lock wait timeout exceeded; 1213 deadlock."""

_NOT_IN_PLACE = 1846
"""The server cannot do an ALTER with the algorithm or lock asked for."""

_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _name(value):
    """`value` as an SQL identifier; refuses anything but letters, digits and underscores."""
    if not _NAME.match(value):
        raise ValueError(f"not a plain identifier: {value!r}")
    return value


# --- what exists ----------------------------------------------------------------------------------------------------

def exists_table(connection, table):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = :t"),
        {"t": table}).scalar())


def exists_column(connection, table, column):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND COLUMN_NAME = :c"), {"t": table, "c": column}).scalar())


def exists_index(connection, table, index):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND INDEX_NAME = :i"), {"t": table, "i": index}).scalar())


def index_columns(connection, table, index):
    """The columns of an index in order (empty when there is none)."""
    rows = connection.execute(sa.text(
        "SELECT COLUMN_NAME FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND INDEX_NAME = :i ORDER BY SEQ_IN_INDEX"), {"t": table, "i": index}).fetchall()
    return [row[0] for row in rows]


def column_type(connection, table, column):
    """A column's type as the server writes it (`bigint unsigned`), or None when there is no such column."""
    return connection.execute(sa.text(
        "SELECT COLUMN_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND COLUMN_NAME = :c"), {"t": table, "c": column}).scalar()


# --- the store's connection ------------------------------------------------------------------------------------------

class StoreConnection:
    """The store connection's `execute(sql, params)` (`?` marks, dict rows), `executemany` and `commit` on a revision's
    SQLAlchemy connection, for the `db` modules that a revision runs (`db.nebula_blocks.rebuild`). The revision commits
    as a whole, so `commit` does nothing."""

    def __init__(self, connection):
        self.connection = connection

    @staticmethod
    def _named(sql, params):
        names = iter(f"p{i}" for i in range(len(params)))
        return re.sub(r"\?", lambda _m: ":" + next(names), sql), {f"p{i}": v for i, v in enumerate(params)}

    def execute(self, sql, params=()):
        sql, bound = self._named(sql, tuple(params))
        return _StoreResult(self.connection.execute(sa.text(sql), bound))

    def executemany(self, sql, rows):
        rows = list(rows)
        if rows:
            sql, _ = self._named(sql, rows[0])
            self.connection.execute(sa.text(sql), [{f"p{i}": v for i, v in enumerate(row)} for row in rows])

    def execute_literal(self, sql):
        """A statement that already holds its values (the store connection's `execute_literal`)."""
        self.connection.execute(sa.text(sql))

    def commit(self):
        """The revision commits as a whole."""


class _StoreResult:
    def __init__(self, result):
        self.result = result

    def fetchone(self):
        row = self.result.mappings().fetchone()
        return None if row is None else dict(row)

    def fetchall(self):
        return [dict(row) for row in self.result.mappings().fetchall()]


# --- lock waits -----------------------------------------------------------------------------------------------------

@contextlib.contextmanager
def lock_wait_timeout(connection, seconds):
    """Sets the session's `lock_wait_timeout` for the block and puts the old value back after it."""
    old = connection.execute(sa.text("SELECT @@SESSION.lock_wait_timeout")).scalar()
    connection.execute(sa.text("SET SESSION lock_wait_timeout = :s"), {"s": int(seconds)})
    try:
        yield
    finally:
        connection.execute(sa.text("SET SESSION lock_wait_timeout = :s"), {"s": int(old)})


def long_transactions(connection, seconds=LONG_TRANSACTION_S):
    """
    The other connections' transactions that have been open `seconds` or more, as `(thread id, seconds open)`.
    Empty when the server does not let this account see them (no PROCESS privilege), as nothing can be said.
    """
    try:
        rows = connection.execute(sa.text(
            "SELECT trx_mysql_thread_id, TIMESTAMPDIFF(SECOND, trx_started, NOW()) AS age"
            " FROM information_schema.INNODB_TRX"
            " WHERE trx_mysql_thread_id <> CONNECTION_ID() AND TIMESTAMPDIFF(SECOND, trx_started, NOW()) >= :s"
            " ORDER BY age DESC"), {"s": int(seconds)}).fetchall()
    except exc.DBAPIError:
        return []
    return [(row[0], row[1]) for row in rows]


def require_no_long_transactions(connection, seconds=LONG_TRANSACTION_S):
    """
    Raises RuntimeError naming the transactions when any is `seconds` old or more: an index build started behind
    one waits for its metadata lock and holds every later query on the table back (docs/design/operations/
    db-check-and-parity-repair.md section 7.3).
    """
    found = long_transactions(connection, seconds)
    if found:
        detail = ", ".join(f"connection {thread} for {age} s" for thread, age in found)
        raise RuntimeError(f"Transactions have been open too long to change the schema safely ({detail}); "
                           f"stop what is using the database and run the update again.")


def _with_retries(connection, statement):
    """Runs `statement`, trying again after a lock wait timeout or deadlock, up to `LOCK_RETRIES` times."""
    for attempt in range(LOCK_RETRIES):
        try:
            connection.execute(sa.text(statement))
            return
        except exc.OperationalError as error:
            code = error.orig.args[0] if getattr(error.orig, "args", None) else None
            if code not in _LOCK_ERRORS or attempt == LOCK_RETRIES - 1:
                raise
            time.sleep(1 + attempt)


# --- schema changes -------------------------------------------------------------------------------------------------

def add_column_instant(connection, table, column, definition, after=None):
    """
    `ALTER TABLE table ADD COLUMN column definition [AFTER after], ALGORITHM=INSTANT`: only the table's definition
    changes, however many rows it has. `definition` is the type and everything after it (`INT NOT NULL DEFAULT 0`).
    Does nothing when the column is there. The server refuses (error 1845) when it cannot do it instantly.
    """
    _name(table), _name(column)
    if exists_column(connection, table, column):
        return False
    position = f" AFTER {_name(after)}" if after else ""
    _with_retries(connection, f"ALTER TABLE {table} ADD COLUMN {column} {definition}{position}, ALGORITHM=INSTANT")
    return True


def drop_column_instant(connection, table, column):
    """`DROP COLUMN ... ALGORITHM=INSTANT`; does nothing when the column is already gone."""
    _name(table), _name(column)
    if not exists_column(connection, table, column):
        return False
    _with_retries(connection, f"ALTER TABLE {table} DROP COLUMN {column}, ALGORITHM=INSTANT")
    return True


def add_index_online(connection, table, name, columns, unique=False):
    """
    `ALTER TABLE table ADD [UNIQUE] KEY name (columns), ALGORITHM=INPLACE, LOCK=NONE`: readers and writers carry
    on while the index builds (minutes on a table of 10^8 rows, with no progress to report). An index already
    there with the same columns is left alone; one of that name with other columns is an error.

    Args:
        columns (sequence[str]): Column names, each optionally `"name DESC"`.
    """
    _name(table), _name(name)
    wanted = [_name(column.split()[0]) for column in columns]
    parts = []
    for column in columns:
        words = column.split()
        if len(words) > 2 or (len(words) == 2 and words[1].upper() not in ("ASC", "DESC")):
            raise ValueError(f"not a column with an optional direction: {column!r}")
        parts.append(" ".join([words[0]] + [word.upper() for word in words[1:]]))
    present = index_columns(connection, table, name)
    if present:
        if present != wanted:
            raise RuntimeError(f"{table}.{name} exists with columns {present}, not {wanted}")
        return False
    kind = "UNIQUE KEY" if unique else "KEY"
    _with_retries(connection, f"ALTER TABLE {table} ADD {kind} {name} ({', '.join(parts)}), "
                              f"ALGORITHM=INPLACE, LOCK=NONE")
    return True


def add_expression_indexes(connection, table, indexes):
    """
    Adds indexes with expression key parts (functional indexes) to `table` in one `ALTER TABLE ... ALGORITHM=INPLACE, LOCK=SHARED`. The
    server keeps a hidden virtual column for each expression and refuses `LOCK=NONE` for them, so readers carry on while the indexes build
    and writers to the table wait. It also refuses to add a second one in place to a table that already has one in some orders (error
    1846), which is why several go in one statement; when it still refuses (a run that stopped halfway left one of the set), the
    statement is asked again with the lock only and the server chooses how.

    Args:
        indexes (dict[str, sequence[str]]): Index name -> its key parts, each a column name or an expression written in brackets
            (`"(x IS NOT NULL)"`; the brackets of the key part are added). The parts are SQL of the caller's own, not of a request.
            An index already there is left alone.

    Returns:
        bool: Whether any index was added.
    """
    _name(table)
    clauses = []
    for name, parts in indexes.items():
        _name(name)
        if exists_index(connection, table, name):
            continue
        keys = []
        for part in parts:
            if part.startswith("("):
                if ";" in part:
                    raise ValueError(f"not an expression: {part!r}")
                keys.append(f"({part})")
            else:
                keys.append(_name(part))
        clauses.append(f"ADD KEY {name} ({', '.join(keys)})")
    if not clauses:
        return False
    try:
        _with_retries(connection, f"ALTER TABLE {table} {', '.join(clauses)}, ALGORITHM=INPLACE, LOCK=SHARED")
    except exc.OperationalError as error:
        if (error.orig.args[0] if getattr(error.orig, "args", None) else None) != _NOT_IN_PLACE:
            raise
        _with_retries(connection, f"ALTER TABLE {table} {', '.join(clauses)}, LOCK=SHARED")
    return True


def drop_index_online(connection, table, name):
    """`DROP KEY ... ALGORITHM=INPLACE, LOCK=NONE`; does nothing when the index is gone."""
    _name(table), _name(name)
    if not exists_index(connection, table, name):
        return False
    _with_retries(connection, f"ALTER TABLE {table} DROP KEY {name}, ALGORITHM=INPLACE, LOCK=NONE")
    return True


# --- updating rows --------------------------------------------------------------------------------------------------

def batched_update(connection, table, assignments, where, key="id", batch=BATCH_ROWS, label=None, params=None):
    """
    `UPDATE table SET assignments WHERE where` in pieces: pages of `batch` rows by `key` (the primary key, an
    increasing integer), one commit after each, so there is no giant undo log and the purge keeps up. `where`
    must stop matching a row once it is updated (`new_col IS NULL`), which is also what lets a stopped run start
    again where it stopped. Progress goes to the migration's bar (`alembic_runner.report_progress`) under `label`.

    Args:
        assignments (str): The SET list, SQL (`b = a * 2`).
        where (str): The condition that picks the rows still to do.
        params (dict, optional): Values for `:name` markers in `assignments` and `where`.

    Returns:
        int: Rows looked at (the rows `where` matched when the call began).
    """
    _name(table), _name(key)
    params = dict(params or {})
    total = connection.execute(sa.text(f"SELECT COUNT(*) FROM {table} WHERE {where}"), params).scalar()
    done = 0
    last = -1
    while True:
        ids = connection.execute(sa.text(
            f"SELECT {key} FROM {table} WHERE {key} > :_last AND ({where}) ORDER BY {key} LIMIT :_n"),
            {**params, "_last": last, "_n": batch}).scalars().all()
        if not ids:
            break
        connection.execute(sa.text(
            f"UPDATE {table} SET {assignments} WHERE {key} BETWEEN :_lo AND :_hi AND ({where})"),
            {**params, "_lo": ids[0], "_hi": ids[-1]})
        connection.commit()
        last = ids[-1]
        done += len(ids)
        if label:
            alembic_runner.report_progress(label, done, total)
    if label and total == 0:
        alembic_runner.report_progress(label, 0, 0)
    return done
