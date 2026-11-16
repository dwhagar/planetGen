# tests/test_command_touches.py

"""
OPS.241, the independence harness: every command declares the tables it reads and writes, and a SQL tap on a toy galaxy holds the declarations to what the statements really touch.

Reuses `queue/needs.py`'s idea (a declared row per command kind), `commands.REGISTRY` and `commands.END_COMMANDS` for the list of commands, `tests/inline_queue.py` through `bughunt_support.run_cli`
for in-process runs, and `test_galaxy_gen._seed_skeleton` for the toy galaxy.
"""

import json
import time

import pytest

from planetgen.queue import commands, engine, sqltap, touches, wrapper
from tests.bughunt_support import mysql_argv, run_cli
from tests.test_galaxy_gen import _seed_skeleton
from tests.test_gen_resume import SHELLS

UNTAPPED = {
    "check-datasets": "reads the raw dataset files and the setup tables; a toy galaxy has no datasets",
    "check-math": "a pure arithmetic check; it opens no connection (declared empty)",
    "end-map-rebuild-maintenance": "the maintenance of the map rebuild (analyze the map tables), declared from its code",
    "plan-skeleton": "writes the four galaxy-shape tables; the scenarios start from a seeded skeleton",
    "reset-galaxy": "wipes every content table; declared ANY and not run against the toy galaxy",
    "empty-nebula-field": "empties nebula_field; declared from its code",
    "empty-phenomenon-scatter": "empties the phenomenon scatter tables; declared from its code",
    "scatter-band": "one band of the scatter; its tables are those of scatter-mass and scatter-luminosity (declared as the union)",
    "conclude-scatter-backfill": "the end hook of scatter-backfill; declared from its code",
    "conclude-warm-cluster-field": "writes the cluster-field file only",
    "warm-cluster-field": "writes the cluster-field file only",
    "warm-galaxy-deep-nodes": "a warm-up answers through the site's pages; declared ANY",
    "warm-galaxy-map": "a warm-up answers through the site's pages; declared ANY",
    "warm-galaxy-nodes": "a warm-up answers through the site's pages; declared ANY",
    "warm-lists": "a warm-up answers through the site's pages; declared ANY",
    "warm-object-index": "a warm-up answers through the site's pages; declared ANY",
    "warm-phenomena": "a warm-up answers through the site's pages; declared ANY",
    "warm-sectors": "a warm-up answers through the site's pages; declared ANY",
    "warm-search": "a warm-up answers through the site's pages; its tables were read once by hand (reads only) and are declared as they were seen",
    "warm-systems": "a warm-up answers through the site's pages; declared ANY",
    "finish-warm-galaxy-deep-nodes": "end hook of a warm-up; declared ANY",
    "finish-warm-galaxy-nodes": "end hook of a warm-up; declared ANY",
    "finish-warm-phenomena": "end hook of a warm-up; declared ANY",
    "finish-warm-sectors": "end hook of a warm-up; declared ANY",
    "finish-warm-systems": "end hook of a warm-up; declared ANY",
}
"""dict: Commands the tap run does not exercise, with the reason; their declaration is read from their code and made wide where in doubt. Anything not here that the run does not exercise fails the test."""

SCENARIOS = [
    ("galaxy", ["--shell", "1", "--band", "0", "--limit", "1", "--num-systems", "1"]),
    ("galaxy", ["--radius-pc", "4", "--num-systems", "1"]),
    ("clear-scatter", ["--luminosity"]),
    ("clear-scatter", ["--mass"]),
    ("clear-scatter", ["--phenomena"]),
    ("clear-scatter", ["--nebula-field"]),
    ("scatter", ["--mass"]),
    ("scatter", ["--luminosity"]),
    ("scatter", ["--phenomena"]),
    ("clear-scatter", ["--nebula-field"]),
    ("scatter", ["--nebula-field"]),
    ("check-db", ["--yes"]),
    ("rescore-habitability", []),
    ("fill-map-luminosity", []),
    ("refill-effect-ledger", []),
    ("clean-nebula-field", []),
    ("link", []),
    ("settle", []),
    ("population", []),
    ("territories", []),
    ("backfill", []),
    ("fill-block-totals", []),
    ("fill-nebula-blocks", []),
    ("sector", ["--num-sectors", "1", "--num-systems", "1"]),
    ("system", ["--name", "Tapped"]),
]
"""list: The command line runs of the toy galaxy that exercise the commands (name, arguments); a few at a time, because flags that combine in one run share one plan."""


def _all_kinds():
    return set(commands.REGISTRY) | {commands.end_name(kind, role) for role, kind in commands.END_COMMANDS}


def _held(observed, declared, writes_cover=()):
    return set(observed) - set(declared) - set(writes_cover) - touches.ENGINE_TABLES


class TestDeclarations:
    """The declarations themselves: complete, well formed, and the conflict rule."""

    def test_every_registered_command_declares(self):
        kinds = _all_kinds()
        assert not kinds - set(touches.READS), "commands with no reads row (add one to queue/touches.py)"
        assert not kinds - set(touches.WRITES), "commands with no writes row (add one to queue/touches.py)"

    def test_no_declaration_names_an_unregistered_command(self):
        kinds = _all_kinds()
        assert not set(touches.READS) - kinds
        assert not set(touches.WRITES) - kinds
        assert not set(touches.OTHER) - kinds

    def test_declarations_are_sorted_tuples_of_names(self):
        for table in (touches.READS, touches.WRITES):
            for kind, names in table.items():
                assert isinstance(names, tuple), kind
                assert list(names) == sorted(set(names)), f"{kind}: sorted and without repeats"
                assert touches.ANY not in names or names == (touches.ANY,), f"{kind}: ANY stands alone"
                assert not set(names) & touches.ENGINE_TABLES, f"{kind}: the engine's own tables are not declared"

    def test_unknown_commands_conflict_with_everything(self):
        assert touches.conflict("not-a-command", "check-counts")
        assert touches.conflict("check-counts", "not-a-command")

    def test_a_write_conflicts_with_a_read_and_a_write_of_the_same_table(self):
        writer = next(kind for kind in touches.WRITES if "sectors" in touches.WRITES[kind] and kind != "reset-galaxy")
        reader = next(kind for kind in touches.READS if "sectors" in touches.READS[kind] and not touches.WRITES[kind] and touches.ANY not in touches.READS[kind])
        assert touches.conflict(writer, reader) and touches.conflict(reader, writer)
        assert touches.conflict(writer, writer)

    def test_two_readers_do_not_conflict(self):
        readers = [kind for kind in touches.READS if not touches.WRITES[kind] and touches.ANY not in touches.READS[kind]]
        assert len(readers) >= 2
        assert not touches.conflict(readers[0], readers[1])

    def test_any_conflicts_with_a_writer(self):
        assert touches.conflict("reset-galaxy", "check-counts")
        assert touches.conflict("warm-sectors", "generate-sectors")

    def test_files_are_conflicts(self):
        assert touches.conflict("scatter-backfill", "warm-cluster-field")

    def test_shared_row_ids_conflict(self):
        for kinds in (("generate-sectors", "generate-system"), ("generate-start-sector", "generate-sectors")):
            assert touches.conflict(*kinds), kinds


class TestIndependenceFile:
    """`independence.json`: the proven pairs, and that nothing overlaps that is not listed."""

    def test_the_shipped_file_loads(self):
        pairs = touches.load_pairs()
        for pair, entry in pairs.items():
            assert all(key in entry for key in touches.PAIR_KEYS)

    def test_an_unlisted_pair_never_overlaps(self):
        assert not touches.may_overlap("scatter-mass", "scatter-luminosity", pairs={})
        assert touches.must_wait_for(["scatter-mass", "scatter-luminosity", "check-counts"], 2, pairs={}) == [0, 1]

    def test_a_listed_pair_overlaps_and_is_symmetric(self):
        pairs = {frozenset(("check-counts", "check-ids")): {"a": "check-counts", "b": "check-ids", "grain": "table", "shared": [], "proof": "x"}}
        assert touches.may_overlap("check-ids", "check-counts", pairs=pairs)
        assert touches.must_wait_for(["check-counts", "check-ids"], 1, pairs=pairs) == []

    def test_a_bad_file_is_refused(self, tmp_path):
        def load(data):
            path = tmp_path / "i.json"
            path.write_text(json.dumps(data))
            return touches.load_pairs(str(path))
        with pytest.raises(ValueError):
            load({"format": 2, "pairs": []})
        with pytest.raises(ValueError):
            load({"format": 1, "pairs": [{"a": "check-ids", "b": "check-ids", "grain": "table", "shared": [], "proof": "x"}]})
        with pytest.raises(ValueError):
            load({"format": 1, "pairs": [{"a": "check-ids", "b": "no-such-command", "grain": "table", "shared": [], "proof": "x"}]})
        with pytest.raises(ValueError):
            load({"format": 1, "pairs": [{"a": "check-ids", "b": "check-counts"}]})
        entry = {"a": "check-ids", "b": "check-counts", "grain": "table", "shared": [], "proof": "x"}
        with pytest.raises(ValueError):
            load({"format": 1, "pairs": [entry, dict(entry, a="check-counts", b="check-ids")]})
        assert len(load({"format": 1, "pairs": [entry]})) == 1

    def test_a_listed_pair_names_what_it_shares(self):
        for pair, entry in touches.load_pairs().items():
            kinds = sorted(pair)
            assert entry["shared"] == touches.shared_between(*kinds), f"{kinds}: what the commands share is not what the entry says"


class TestSqlTap:
    """`tables_of`, the statement reader."""

    def test_select_insert_update_delete(self):
        assert sqltap.tables_of("SELECT a FROM stars s JOIN planets p ON p.id = s.id") == (("planets", "stars"), (), ())
        assert sqltap.tables_of("INSERT INTO sectors (a) VALUES (%s)") == ((), ("sectors",), ())
        assert sqltap.tables_of("INSERT IGNORE INTO sectors (a) SELECT a FROM stars") == (("stars",), ("sectors",), ())
        assert sqltap.tables_of("UPDATE planets SET a = 1 WHERE id IN (SELECT id FROM moons)") == (("moons",), ("planets",), ())
        assert sqltap.tables_of("DELETE FROM moons WHERE id = 1") == ((), ("moons",), ())
        assert sqltap.tables_of("TRUNCATE TABLE nebula_field") == ((), ("nebula_field",), ())

    def test_on_duplicate_key_update_names_columns_not_tables(self):
        assert sqltap.tables_of("INSERT INTO engine_jobs (a) VALUES (1) ON DUPLICATE KEY UPDATE attempt = attempt + 1, result = VALUES(result)") == ((), ("engine_jobs",), ())

    def test_with_names_and_strings_and_comments_are_not_tables(self):
        sql = "/* FROM hidden */ WITH t AS (SELECT a FROM stars) SELECT a FROM t WHERE b = 'FROM nothing'"
        assert sqltap.tables_of(sql) == (("stars",), (), ())

    def test_schema_and_session_statements_touch_nothing(self):
        assert sqltap.tables_of("SELECT x FROM information_schema.tables") == ((), (), ())
        assert sqltap.tables_of("SET SESSION foo = 1") == ((), (), ())
        assert sqltap.tables_of("SHOW TABLES") == ((), (), ())

    def test_locks(self):
        assert sqltap.tables_of("SELECT GET_LOCK('link_neighbors', 5)")[2] == ("lock:link_neighbors",)
        reads, writes, locks = sqltap.tables_of("LOCK TABLES stars WRITE, planets READ")
        assert (reads, writes, locks) == (("planets",), ("stars",), ("lock:planets", "lock:stars"))

    def test_the_running_command_gets_the_statement(self):
        tap = sqltap.Tap()
        with tap.running("a"):
            tap.record("SELECT 1 FROM stars")
            with tap.running("b"):
                tap.record("DELETE FROM moons")
            tap.record("INSERT INTO planets (x) VALUES (1)")
        tap.record("SELECT * FROM sectors")
        assert tap.touched["a"]["reads"] == {"stars"} and tap.touched["a"]["writes"] == {"planets"}
        assert tap.touched["b"]["writes"] == {"moons"}
        assert tap.touched[None]["reads"] == {"sectors"}


@pytest.mark.db
@pytest.mark.usefixtures("unit_density_plans")
def test_commands_touch_only_what_they_declare(mysql_config, monkeypatch):
    """
    Runs the toy galaxy's commands under the tap and fails any command that reads a table it does not declare (a write covers a read), writes a table it does not declare, or is neither
    exercised nor in `UNTAPPED`. The engine's own tables are left out.
    """
    _seed_skeleton(mysql_config, shells=SHELLS)
    tap = sqltap.Tap()
    real_run = engine.run
    real_expand = commands.expand_function
    real_end = commands.end_function

    def run(kind, *args, **kwargs):
        with tap.running(kind):
            return real_run(kind, *args, **kwargs)

    def expand(kind):
        function = real_expand(kind)

        def cut(*args, **kwargs):
            with tap.running(kind):
                return function(*args, **kwargs)
        return cut

    def end(kind, role):
        function = real_end(kind, role)
        if function is None:
            return function
        name = commands.end_name(kind, role)

        def hook(*args, **kwargs):
            with tap.running(name):
                return function(*args, **kwargs)
        return hook

    monkeypatch.setattr(engine, "run", run)
    monkeypatch.setattr(commands, "expand_function", expand)
    monkeypatch.setattr(commands, "end_function", end)
    failed = []
    with tap:
        for name, argv in SCENARIOS:
            try:
                run_cli(name, [*argv, *mysql_argv(mysql_config)])
            except SystemExit as exit_:
                if exit_.code:
                    failed.append(f"{name} {argv} exited {exit_.code}")
            except Exception as error:
                failed.append(f"{name} {argv} raised {error!r}")
        params = {"database": wrapper.database_args(mysql_config)}
        for kind in sorted(commands.REGISTRY):
            if kind.startswith("check-") and kind not in UNTAPPED:
                commands.run_alone(kind, params)
    problems = []
    for kind, slot in sorted((k, v) for k, v in tap.touched.items() if k is not None):
        if kind not in touches.READS:
            problems.append(f"{kind}: touched tables but has no declaration")
            continue
        extra_reads = _held(slot["reads"], touches.reads_of(kind), touches.writes_of(kind)) if touches.ANY not in touches.reads_of(kind) else set()
        extra_writes = _held(slot["writes"], touches.writes_of(kind)) if touches.ANY not in touches.writes_of(kind) else set()
        if extra_reads:
            problems.append(f"{kind}: reads undeclared {sorted(extra_reads)}")
        if extra_writes:
            problems.append(f"{kind}: writes undeclared {sorted(extra_writes)}")
    assert not failed, "\n".join(failed)
    assert not problems, "\n".join(problems)
    unexercised = _all_kinds() - set(tap.touched) - set(UNTAPPED)
    assert not unexercised, f"commands the tap run never exercised (run them, or list them in UNTAPPED with the reason): {sorted(unexercised)}"
    stale = set(UNTAPPED) & set(tap.touched)
    assert not stale, f"in UNTAPPED but exercised now (remove the row): {sorted(stale)}"


class TestPlanEdges:
    """The planner's edges (`touches.plan_edges`): the chain of plan order until a pair is proven."""

    KINDS = ["generate-sectors", "scatter-mass", "scatter-luminosity", "check-counts"]

    def test_with_no_pair_listed_each_entry_waits_for_the_one_before(self):
        assert touches.plan_edges(self.KINDS, pairs={}) == [[], [0], [1], [2]]

    def test_the_shipped_plan_is_the_chain_until_a_pair_is_listed(self):
        if not touches.load_pairs():
            assert touches.plan_edges(self.KINDS) == [[], [0], [1], [2]]

    def test_a_proven_pair_leaves_its_two_entries_free_of_each_other(self):
        pairs = {frozenset(("scatter-mass", "scatter-luminosity")): {"a": "scatter-mass", "b": "scatter-luminosity", "grain": "table", "shared": [], "proof": "x"}}
        assert touches.plan_edges(self.KINDS, pairs=pairs) == [[], [0], [0], [1, 2]]

    def test_a_proven_pair_does_not_free_a_third_entry(self):
        pairs = {frozenset(("scatter-mass", "scatter-luminosity")): {"a": "scatter-mass", "b": "scatter-luminosity", "grain": "table", "shared": [], "proof": "x"}}
        edges = touches.plan_edges(["generate-sectors", "scatter-mass", "scatter-luminosity", "generate-system"], pairs=pairs)
        assert edges[3] == [1, 2]
