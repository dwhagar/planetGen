# tests/test_sql_portability.py

"""
SQL portability lint (TEST.6): the schema files and every SQL string in
the non-test code under src/ must run unchanged on MySQL 8.0, MySQL 8.4
and MariaDB 10.x/11.x. A query once used `AS generated` -- fine on
MariaDB, reserved on MySQL 8 -- and broke production.

Checked, without a database:

- No unquoted table, column, index, constraint name or alias is a
  reserved word on any of those engines (`sql_reserved_words.py`);
  in backticks it's fine.
- No new `VALUES(col)` in `ON DUPLICATE KEY UPDATE` (deprecated since
  MySQL 8.0.20). See `VALUES_FUNCTION_ALLOWED` for why the existing ones
  stay.

The SQL in .py files is every string constant that looks like a statement
or a fragment of one (docstrings aside), implicitly concatenated strings
joined and f-strings' `{...}` parts read as a plain identifier.
"""

import ast
import re
from pathlib import Path

import pytest

from tests import sql_reserved_words

SRC = Path(__file__).resolve().parent.parent

SCHEMA_FILES = ("planetgen/db/schema.sql", "planetgen/db/control_schema.sql")

# MySQL 8.0's replacement for VALUES(col), `INSERT ... VALUES (...) AS new
# ON DUPLICATE KEY UPDATE x = new.x`, is a syntax error on MariaDB (10.11
# and 11.x), and VALUES(col) still works on MySQL 8.0 and 8.4 (a
# deprecation warning, 1287) -- tried on both servers, 2026-10-01. So
# VALUES() is the only portable form for now and these uses stay; a new
# one fails the lint until it's added here, so they stay counted. Swap
# them all when MySQL removes VALUES() or MariaDB gains the row alias.
# {file: {column, ...}} -- "__expr__" is a column an f-string names.
VALUES_FUNCTION_ALLOWED = {
    "planetgen/db/store.py": {
        "occurrence_count", "first_object_table", "first_star_system_id", "first_object_id",
        "diminutive_index", "disk_scale_length_pc", "disk_scale_height_pc", "bulge_scale_radius_pc",
        "bulge_amplitude", "arm_count", "pitch_angle_rad", "arm_amplitude", "spiral_reference_radius_pc",
        "spiral_reference_angle_rad", "k_norm", "edge_pc", "expected_system_count_at_density_1",
        "outer_ring_index",
    },
    "planetgen/population/model.py": {"scanned_planet_id"},
    "planetgen/generation/stats.py": {"__expr__", "bytes_per_system", "systems", "total_bytes"},
}

# --- Tokens ----------------------------------------------------------------

_TOKEN = re.compile(r"""
    (?P<ws>\s+)
  | (?P<comment>--[^\n]*|/\*.*?\*/)
  | (?P<string>'(?:[^'\\]|\\.|'')*'|"(?:[^"\\]|\\.|"")*")
  | (?P<quoted>`(?:[^`]|``)*`)
  | (?P<number>\d+(?:\.\d+)?(?:[eE][-+]?\d+)?\b)
  | (?P<word>[A-Za-z_$][A-Za-z0-9_$]*)
  | (?P<op><=>|<=|>=|<>|!=|:=|\S)
""", re.X | re.S)

PLACEHOLDER = "__expr__"

# Python-side placeholders, read as a plain value or identifier.
_PY_PLACEHOLDERS = re.compile(r"\{[^{}'\"]*\}|%\(\w+\)s|%s")


class Token:
    def __init__(self, kind, text, pos):
        self.kind, self.text, self.pos = kind, text, pos
        self.upper = text.upper()

    def is_word(self, *words):
        return self.kind == "word" and (not words or self.upper in words)

    def __repr__(self):
        return f"Token({self.kind}, {self.text!r})"


def tokenize(sql):
    tokens = []
    for match in _TOKEN.finditer(sql):
        kind = match.lastgroup
        if kind not in ("ws", "comment"):
            tokens.append(Token(kind, match.group(), match.start()))
    return tokens


# --- Identifier positions ----------------------------------------------------

# Reserved words this code uses as SQL syntax. Any other reserved word,
# bare, is an identifier wherever it is (`SELECT rank FROM`), unless it
# calls a function (`RANK() OVER`, `LEFT(name, 1)`) or follows a dot
# (`t.rank` -- MySQL never needs quotes after a dot). In the positions
# `_identifier_positions` finds, even these are flagged (`AS order`).
SYNTAX_KEYWORDS = frozenset("""
    ADD ALL ALTER AND AS ASC BETWEEN BIGINT BINARY BLOB BY CASCADE CASE CHANGE CHAR CHARACTER CHECK
    COLLATE COLUMN CONSTRAINT CONVERT CREATE CROSS CURRENT_DATE CURRENT_TIME CURRENT_TIMESTAMP
    DATABASE DECIMAL DEFAULT DELETE DESC DISTINCT DIV DOUBLE DROP DUAL ELSE ELSEIF EXISTS FALSE FLOAT
    FOR FOREIGN FROM FULLTEXT GROUP HAVING IF IGNORE IN INDEX INNER INSERT INT INTEGER INTERVAL INTO IS
    JOIN KEY KEYS LEFT LIKE LIMIT LOCK LONGBLOB LONGTEXT MEDIUMBLOB MEDIUMINT MEDIUMTEXT MOD NATURAL NOT
    NULL ON OR ORDER OUTER PRIMARY REFERENCES REGEXP RENAME REPLACE RESTRICT RIGHT SCHEMA SELECT SET
    SHOW SMALLINT SPATIAL STRAIGHT_JOIN TABLE THEN TINYINT TINYTEXT TO TRUE UNION UNIQUE UNSIGNED UPDATE
    FORCE USE USING UTC_DATE UTC_TIME UTC_TIMESTAMP VALUES VARBINARY VARCHAR WHEN WHERE WITH XOR ZEROFILL
    LOCALTIME LOCALTIMESTAMP DUPLICATE OVER PARTITION WINDOW ROWS RANGE SQL_CALC_FOUND_ROWS LINES
    SQL_NO_CACHE HIGH_PRIORITY LOW_PRIORITY DELAYED READ WRITE OPTIMIZE ANALYZE EXPLAIN DESCRIBE
    TRIGGER EACH ROW BEFORE AFTER GRANT REVOKE USAGE LEADING TRAILING BOTH SEPARATOR
""".split())


def _keyword_in_context(tokens, i):
    """A reserved word outside `SYNTAX_KEYWORDS` that is syntax here."""
    tok = tokens[i]
    prev = tokens[i - 1] if i else None
    nxt = tokens[i + 1] if i + 1 < len(tokens) else None
    if tok.upper == "GENERATED":
        return nxt is not None and nxt.is_word("ALWAYS")
    if tok.upper in ("STORED", "VIRTUAL"):
        return prev is not None and prev.text == ")"
    if tok.upper == "OFFSET":
        return any(t.is_word("LIMIT") for t in tokens[:i]) and nxt is not None and (
            nxt.text == "?" or nxt.kind in ("number", "word"))
    return False


def _paren_end(tokens, i):
    """Index of the `)` closing the `(` at `i` (or the end)."""
    depth = 0
    for j in range(i, len(tokens)):
        if tokens[j].text == "(":
            depth += 1
        elif tokens[j].text == ")":
            depth -= 1
            if depth == 0:
                return j
    return len(tokens)


def _list_items(tokens, start, end):
    """Index of the first token of each top-level comma-separated item
    between `start` and `end` (exclusive)."""
    items, depth, first = [], 0, True
    for j in range(start, end):
        text = tokens[j].text
        if first:
            items.append(j)
            first = False
        if text == "(":
            depth += 1
        elif text == ")":
            depth -= 1
        elif text == "," and depth == 0:
            first = True
    return items


def _column_list(tokens, open_paren):
    """The column names of the `(a, b(10), c DESC)` list at `open_paren`."""
    end = _paren_end(tokens, open_paren)
    return [j for j in _list_items(tokens, open_paren + 1, end) if tokens[j].kind == "word"]


_TABLE_KEYWORDS = ("FROM", "JOIN", "INTO", "UPDATE", "TABLE", "REFERENCES")
_INDEX_WORDS = ("KEY", "INDEX")
_NOT_COLUMN = ("PRIMARY", "KEY", "INDEX", "UNIQUE", "CONSTRAINT", "FOREIGN", "FULLTEXT", "SPATIAL", "CHECK")
_AFTER_AS_SYNTAX = ("SELECT", "WITH", "VALUES", "TABLE")


def _skip_if_exists(tokens, j):
    """Steps over `IF [NOT] EXISTS` at `j`."""
    if j < len(tokens) and tokens[j].is_word("IF"):
        j += 1
        if j < len(tokens) and tokens[j].is_word("NOT"):
            j += 1
        if j < len(tokens) and tokens[j].is_word("EXISTS"):
            j += 1
    return j


def _index_definition(tokens, j, found):
    """`[UNIQUE|FULLTEXT|SPATIAL] [KEY|INDEX] [name] (cols)` at `j`."""
    if j < len(tokens) and tokens[j].is_word("CHECK"):
        return
    while j < len(tokens) and tokens[j].is_word("UNIQUE", "FULLTEXT", "SPATIAL", "PRIMARY", "FOREIGN"):
        j += 1
    if j < len(tokens) and tokens[j].is_word(*_INDEX_WORDS):
        j += 1
    if j < len(tokens) and tokens[j].kind == "word" and not tokens[j].is_word("USING"):
        found.add(j)
        j += 1
    if j < len(tokens) and tokens[j].text == "(":
        found.update(_column_list(tokens, j))
        j = _paren_end(tokens, j) + 1
    if j < len(tokens) and tokens[j].is_word("REFERENCES"):
        j += 2  # the table name is a `_TABLE_KEYWORDS` position
        if j < len(tokens) and tokens[j].text == "(":
            found.update(_column_list(tokens, j))


def _table_element(tokens, j, found):
    """One item of `CREATE TABLE t (...)`: a column or an index."""
    tok = tokens[j]
    if tok.kind != "word":
        return
    if tok.is_word("CONSTRAINT"):
        if j + 1 < len(tokens) and tokens[j + 1].kind == "word" and not tokens[j + 1].is_word(*_NOT_COLUMN):
            found.add(j + 1)
            j += 1
        _index_definition(tokens, j + 1, found)
    elif tok.is_word(*_NOT_COLUMN):
        _index_definition(tokens, j, found)
    else:
        found.add(j)


def _alter_item(tokens, j, found):
    """One item of `ALTER TABLE t ...`."""
    if j >= len(tokens) or tokens[j].kind != "word":
        return
    action = tokens[j].upper
    j += 1
    if action in ("ADD", "DROP", "MODIFY", "CHANGE", "RENAME", "ALTER"):
        if j < len(tokens) and tokens[j].is_word("COLUMN"):
            j += 1
        j = _skip_if_exists(tokens, j)
        if j >= len(tokens):
            return
        if action in ("ADD", "DROP") and tokens[j].is_word(*_NOT_COLUMN):
            if action == "ADD":
                _table_element(tokens, j, found)
            else:
                while j < len(tokens) and tokens[j].is_word("PRIMARY", "FOREIGN", "CONSTRAINT", *_INDEX_WORDS):
                    j += 1
                if j < len(tokens) and tokens[j].kind == "word":
                    found.add(j)
            return
        if action == "RENAME" and tokens[j].is_word("TO", "AS"):
            return
        if action == "RENAME" and tokens[j].is_word(*_INDEX_WORDS):
            j += 1
        if action == "ALTER" and tokens[j].is_word(*_INDEX_WORDS):
            j += 1
        if tokens[j].kind == "word":
            found.add(j)
        if action == "CHANGE" and j + 1 < len(tokens) and tokens[j + 1].kind == "word":
            found.add(j + 1)
        if action == "RENAME" and j + 2 < len(tokens) and tokens[j + 1].is_word("TO"):
            found.add(j + 2)


def _identifier_positions(tokens):
    """Indexes of the tokens that can only be an identifier here."""
    found = set()
    parens = []  # for each open "(", the word before it
    for i, tok in enumerate(tokens):
        prev = tokens[i - 1] if i else None
        nxt = tokens[i + 1] if i + 1 < len(tokens) else None
        if tok.text == "(":
            parens.append(prev.upper if prev is not None else "")
        elif tok.text == ")" and parens:
            parens.pop()
        if tok.kind != "word":
            continue
        # A table or alias qualifying a column: `rank.x`.
        if nxt is not None and nxt.text == "." and (prev is None or prev.text != "."):
            found.add(i)
        # `AS alias`, but not `CAST(x AS CHAR)` or `CREATE TABLE t AS SELECT`.
        if tok.upper == "AS" and nxt is not None and nxt.kind == "word" \
                and not nxt.is_word(*_AFTER_AS_SYNTAX) and not (parens and parens[-1] in ("CAST", "CONVERT")):
            found.add(i + 1)
        # A table name.
        if tok.upper in _TABLE_KEYWORDS and not (
            tok.upper == "UPDATE" and prev is not None and prev.is_word("KEY", "ON", "FOR")
        ):
            j = i + 1
            while j < len(tokens) and tokens[j].is_word("TEMPORARY", "IGNORE", "LOW_PRIORITY"):
                j += 1
            j = _skip_if_exists(tokens, j)
            if j < len(tokens) and tokens[j].kind == "word" and not (
                tok.upper == "FROM" and tokens[j].is_word("DUAL")
            ) and not tokens[j].is_word("SELECT"):
                table = j
                if j + 2 < len(tokens) and tokens[j + 1].text == ".":
                    table = j + 2  # `db.table`: `db` is flagged above
                found.add(table)
                after = table + 1
                # CREATE TABLE t (...) and INSERT INTO t (cols).
                if after < len(tokens) and tokens[after].text == "(":
                    if tok.upper == "TABLE" and prev is not None and prev.is_word("CREATE", "TEMPORARY"):
                        end = _paren_end(tokens, after)
                        for item in _list_items(tokens, after + 1, end):
                            _table_element(tokens, item, found)
                    elif tok.upper == "INTO":
                        found.update(_column_list(tokens, after))
                # ALTER TABLE t ADD ..., DROP ..., ...
                if tok.upper == "TABLE" and prev is not None and prev.is_word("ALTER"):
                    for item in _list_items(tokens, after, len(tokens)):
                        _alter_item(tokens, item, found)
        # CREATE [UNIQUE|FULLTEXT] INDEX name ON t (cols).
        if tok.upper == "INDEX" and prev is not None and prev.is_word("CREATE", "UNIQUE", "FULLTEXT", "SPATIAL"):
            if nxt is not None and nxt.kind == "word":
                found.add(i + 1)
            if i + 4 < len(tokens) and tokens[i + 2].is_word("ON") and tokens[i + 4].text == "(":
                found.add(i + 3)
                found.update(_column_list(tokens, i + 4))
        # `AFTER col` in ALTER TABLE ... ADD COLUMN.
        if tok.upper == "AFTER" and nxt is not None and nxt.kind == "word" and any(
            t.is_word("ALTER") for t in tokens[:i]
        ):
            found.add(i + 1)
        # `col = ...`, `col < ...` and the like, but not `CHARACTER SET = x`.
        if nxt is not None and nxt.text in ("=", "<", ">", "<=", ">=", "<>", "!=", "<=>") \
                and not tok.is_word("SET", "COLLATE", "CHARSET"):
            found.add(i)
    return found


def reserved_identifiers(sql):
    """`[(word, snippet), ...]`: reserved words used unquoted as an
    identifier in `sql`."""
    tokens = tokenize(sql)
    strict = _identifier_positions(tokens)
    problems = []
    for i, tok in enumerate(tokens):
        if tok.kind != "word" or tok.upper not in sql_reserved_words.RESERVED:
            continue
        prev = tokens[i - 1] if i else None
        nxt = tokens[i + 1] if i + 1 < len(tokens) else None
        if prev is not None and prev.text == ".":
            continue
        if i not in strict:
            if tok.upper in SYNTAX_KEYWORDS or (nxt is not None and nxt.text == "("):
                continue
            if _keyword_in_context(tokens, i):
                continue
        problems.append((tok.text, _snippet(sql, tok.pos)))
    return problems


_VALUES_FUNCTION = re.compile(r"\bVALUES\s*\(\s*(`?)([A-Za-z_$][\w$]*)\1\s*\)", re.I)
_ON_DUPLICATE = re.compile(r"\bON\s+DUPLICATE\s+KEY\s+UPDATE\b", re.I)


def values_function_uses(sql):
    """`[(column, snippet), ...]`: every `VALUES(col)` after `ON
    DUPLICATE KEY UPDATE` in `sql`. A fragment of an upsert (built
    without its `ON DUPLICATE KEY UPDATE`) that reads `col = VALUES(col)`
    counts too."""
    stripped = " ".join(t.text for t in tokenize(sql))
    duplicate = _ON_DUPLICATE.search(stripped)
    found = []
    for match in _VALUES_FUNCTION.finditer(stripped):
        before = stripped[:match.start()]
        in_upsert = (duplicate is not None and duplicate.start() < match.start()) or \
            re.search(r"(=|,|\+|-|\(|\bIF\s*\()\s*$", before, re.I) is not None and duplicate is None \
            and not re.search(r"\bINSERT\b", stripped[:match.start()], re.I)
        if in_upsert:
            found.append((match.group(2), _snippet(stripped, match.start())))
    return found


def _snippet(sql, pos, width=40):
    start, end = max(0, pos - width), min(len(sql), pos + width)
    text = " ".join(sql[start:end].split())
    return ("..." if start else "") + text + ("..." if end < len(sql) else "")


# --- Where the SQL is ---------------------------------------------------------

_STATEMENT = re.compile(
    r"^\s*\(?\s*(SELECT|INSERT|UPDATE|DELETE|REPLACE|CREATE|ALTER|DROP|TRUNCATE|RENAME|SHOW|WITH|"
    r"GRANT|REVOKE|LOCK TABLES|UNLOCK TABLES|ANALYZE|OPTIMIZE|EXPLAIN|START TRANSACTION|SET (?:SESSION|@@|NAMES))\b"
)
_FRAGMENT = re.compile(
    r"\b(FROM \w|WHERE\b|(?:LEFT |INNER )?JOIN \w|ORDER BY\b|GROUP BY\b|ON DUPLICATE KEY UPDATE\b|"
    r"VALUES ?\(|LIMIT [?\d{%]|AS [a-z_`]|IS (?:NOT )?NULL\b|(?:AND|OR) \w+ (?:=|<|>|IN|LIKE|IS)|"
    r"\w+ (?:=|<=?|>=?|LIKE|REGEXP|IN) (?:\?|%s|\(\?)|BETWEEN \? AND \?|SET \w+ = |"
    r"ADD (?:COLUMN|CONSTRAINT|INDEX|KEY|UNIQUE)\b|NOT NULL\b|REFERENCES \w|\w IN \(|MATCH ?\()"
)


def looks_like_sql(text):
    """A statement, or a fragment of one: SQL in this project is written
    with upper-case keywords, prose isn't."""
    return bool(_STATEMENT.match(text) or _FRAGMENT.search(text))


def _string_value(node):
    """A constant or f-string's text, with each `{...}` as `PLACEHOLDER`."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.JoinedStr):
        parts = []
        for value in node.values:
            if isinstance(value, ast.Constant):
                parts.append(value.value)
            else:
                parts.append(f" {PLACEHOLDER} ")
        return "".join(parts)
    return None


def python_sql_strings(path):
    """`[(line, sql), ...]` for every SQL-looking string in a .py file."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    docstrings, inner = set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
            docstrings.add(id(node.value))  # docstrings and bare string statements
        if isinstance(node, ast.JoinedStr):
            inner.update(id(v) for v in node.values)
            inner.update(id(v.format_spec) for v in node.values
                         if isinstance(v, ast.FormattedValue) and v.format_spec is not None)
    found = []
    for node in ast.walk(tree):
        if id(node) in docstrings or id(node) in inner:
            continue
        text = _string_value(node)
        if text is None or not looks_like_sql(text):
            continue
        found.append((node.lineno, _PY_PLACEHOLDERS.sub(PLACEHOLDER, text)))
    return found


def sql_sources():
    """`[(where, sql), ...]` for the schema files' statements and every
    SQL string in the non-test code."""
    sources = []
    for name in SCHEMA_FILES:
        text = (SRC / name).read_text(encoding="utf-8")
        # Split at the `;` tokens, so one in a comment or string doesn't.
        start = None
        for tok in tokenize(text):
            if start is None:
                start = tok.pos
            if tok.text == ";":
                sources.append((f"{name}:{text.count(chr(10), 0, start) + 1}", text[start:tok.pos]))
                start = None
        if start is not None:
            sources.append((f"{name}:{text.count(chr(10), 0, start) + 1}", text[start:]))
    for path in sorted(SRC.rglob("*.py")):
        rel = path.relative_to(SRC).as_posix()
        if rel.startswith("tests/") or "/vendor/" in rel or "__pycache__" in rel:
            continue
        for line, sql in python_sql_strings(path):
            sources.append((f"{rel}:{line}", sql))
    return sources


@pytest.fixture(scope="module")
def sources():
    return sql_sources()


# --- The lint itself -------------------------------------------------------------

def test_the_lint_finds_the_sql_it_should_check(sources):
    where = {w.split(":")[0] for w, _sql in sources}
    assert set(SCHEMA_FILES) <= where
    assert {"planetgen/db/store.py", "planetgen/db/query.py"} <= where
    assert sum(1 for w, _sql in sources if w.startswith("planetgen/db/query.py")) > 100
    assert any("CREATE TABLE" in sql for w, sql in sources if w.startswith(SCHEMA_FILES[0]))


def test_no_reserved_word_is_an_unquoted_identifier(sources):
    problems = [
        f"{where}: {word!r} (reserved on {sql_reserved_words.engines(word)}) in: {snippet}"
        for where, sql in sources for word, snippet in reserved_identifiers(sql)
    ]
    assert not problems, "Reserved words used as unquoted identifiers -- quote them in backticks:\n" + \
        "\n".join(problems)


def test_no_new_values_function_in_on_duplicate_key_update(sources):
    problems = [
        f"{where}: VALUES({column}) in: {snippet}"
        for where, sql in sources for column, snippet in values_function_uses(sql)
        if column not in VALUES_FUNCTION_ALLOWED.get(where.split(":")[0], ())
    ]
    assert not problems, (
        "VALUES(col) in ON DUPLICATE KEY UPDATE is deprecated (MySQL 8.0.20+); see VALUES_FUNCTION_ALLOWED "
        "before adding one:\n" + "\n".join(problems)
    )


# --- The lint catches what it should --------------------------------------------------

@pytest.mark.parametrize("sql, word", [
    ("SELECT COUNT(*) AS generated FROM sectors", "generated"),
    ("SELECT COUNT(*) AS `generated` FROM sectors, (SELECT 1 AS rank) r", "rank"),
    ("CREATE TABLE t (id INT PRIMARY KEY, rank INT NOT NULL)", "rank"),
    ("CREATE TABLE t (id INT, KEY offset (id))", "offset"),
    ("CREATE TABLE t (id INT, CONSTRAINT window FOREIGN KEY (id) REFERENCES u (id))", "window"),
    ("CREATE INDEX idx ON t (groups)", "groups"),
    ("ALTER TABLE t ADD COLUMN lead INT AFTER id", "lead"),
    ("INSERT INTO t (id, row) VALUES (?, ?)", "row"),
    ("UPDATE t SET system = ? WHERE id = ?", "system"),
    ("SELECT id, rank FROM t", "rank"),
    ("SELECT x FROM returning WHERE id = ?", "returning"),
    ("SELECT x FROM t WHERE offset > 3", "offset"),
    ("SELECT x FROM t AS order", "order"),
    ("SELECT key.x FROM t", "key"),
])
def test_lint_flags_reserved_identifiers(sql, word):
    assert [found for found, _snippet in reserved_identifiers(sql)] == [word]


@pytest.mark.parametrize("sql", [
    "SELECT COUNT(*) AS `generated` FROM sectors",
    "CREATE TABLE t (id INT, `rank` INT, KEY `offset` (`rank`)) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4",
    "SELECT t.rank, RANK() OVER (ORDER BY id) AS r, CAST(x AS CHAR) AS c, LEFT(name, 1) FROM t LIMIT ? OFFSET ?",
    "CREATE TABLE t (id INT, g INT GENERATED ALWAYS AS (id + 1) STORED, v INT AS (id) VIRTUAL)",
    "INSERT INTO t (id, n) VALUES (?, ?) ON DUPLICATE KEY UPDATE n = n + 1, updated_at = CURRENT_TIMESTAMP",
    "SELECT 'rank' AS label, \"offset\" AS other FROM DUAL -- generated, rank",
    "UPDATE t SET name = ? WHERE id = ? AND x IS NOT NULL ORDER BY id DESC FOR UPDATE",
])
def test_lint_passes_keywords_and_quoted_identifiers(sql):
    assert reserved_identifiers(sql) == []


def test_lint_flags_values_function_in_upserts():
    sql = "INSERT INTO t (id, x) VALUES (?, ?) ON DUPLICATE KEY UPDATE x = VALUES(x), y = IF(a, VALUES(`y`), y)"
    assert [column for column, _snippet in values_function_uses(sql)] == ["x", "y"]
    # A fragment built on its own, and VALUES (...) as the row list, which is fine.
    assert [c for c, _s in values_function_uses(" first_id = IF(v, VALUES(first_id), first_id),")] == ["first_id"]
    assert values_function_uses("INSERT INTO t (id, x) VALUES (?, ?) AS new ON DUPLICATE KEY UPDATE x = new.x") == []
    assert values_function_uses("INSERT INTO t (x) VALUES (x)") == []


def test_python_strings_are_read_with_their_concatenated_and_f_string_parts(tmp_path):
    path = tmp_path / "module.py"
    path.write_text(
        '"""SELECT rank FROM docstrings is prose."""\n'
        'col = "x"\n'
        'sql = ("SELECT COUNT(*) "\n'
        '       "AS generated FROM t")\n'
        'other = f"SELECT {col} AS rank FROM t WHERE id = ?"\n'
        'prose = "Generated 3 sectors"\n'
    )
    found = python_sql_strings(path)
    assert [line for line, _sql in found] == [3, 5]
    assert [[w for w, _s in reserved_identifiers(sql)] for _line, sql in found] == [["generated"], ["rank"]]
