# tests/test_job_tree.py

"""
Tests for the job tree (ADM.12, `stellarObjects/workQueue.py`'s
`open_node`/`job_node`/`load_tree`, control schema v7): nodes nest under
the one a process has open, under the job runner's step node across
processes, record their own start, end and duration, and add their
totals and ETA up to their parents; an older control schema gains the
v7 columns on the next deploy.

Tests that take the `mysql_config` fixture (see `conftest.py`) are
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import os
import sys
import time

import pytest

import generate
import jobRunner
from stellarObjects import workQueue
from planetgen.db import store as _db

from tests.test_galaxy_gen import _mysql_argv, _plan_wide_galaxy


def _square(payload):
    return payload * payload


def _slow(payload):
    time.sleep(0.05)
    return payload


@pytest.fixture
def control_config(mysql_config, monkeypatch):
    conn = _db.get_control_connection(mysql_config, ensure_schema=True)
    conn.close()
    monkeypatch.setenv(_db.CONTROL_DB_ENV_VAR, mysql_config.database)
    monkeypatch.delenv(workQueue.PARENT_ENV_VAR, raising=False)
    return mysql_config


@pytest.fixture(autouse=True)
def _no_open_nodes():
    yield
    assert workQueue.current_node() is None, "a test left a job tree node open"


def _conn(config):
    return _db.get_control_connection(config)


def _rows(config, sql, params=()):
    conn = _conn(config)
    try:
        return conn.execute(sql, params).fetchall()
    finally:
        conn.close()


def _tree(config, root_id):
    conn = _conn(config)
    try:
        return workQueue.load_tree(conn, root_id)
    finally:
        conn.close()


def _walk(node, depth=0):
    yield depth, node
    for child in node["children"]:
        yield from _walk(child, depth + 1)


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

_V6_WORK_JOBS = """CREATE TABLE work_jobs (
    id VARCHAR(32) NOT NULL PRIMARY KEY, title VARCHAR(255) NOT NULL, holder VARCHAR(255) NOT NULL,
    state VARCHAR(16) NOT NULL, workers INT UNSIGNED NOT NULL,
    tasks_queued INT UNSIGNED NOT NULL DEFAULT 0, tasks_done INT UNSIGNED NOT NULL DEFAULT 0,
    tasks_failed INT UNSIGNED NOT NULL DEFAULT 0, created_at DATETIME(6) NOT NULL,
    started_at DATETIME(6) NULL, finished_at DATETIME(6) NULL, heartbeat_at DATETIME(6) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def test_an_older_control_schema_gains_the_tree_columns(mysql_config):
    conn = _db.get_control_connection(mysql_config)
    try:
        conn.execute(_V6_WORK_JOBS)
        conn.execute("CREATE TABLE work_lease (id TINYINT UNSIGNED NOT NULL PRIMARY KEY, holder VARCHAR(255) NULL,"
                     " job_id VARCHAR(32) NULL, heartbeat_at DATETIME(6) NULL)"
                     " ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci")
        conn.execute("INSERT INTO work_jobs (id, title, holder, state, workers, created_at, heartbeat_at)"
                     " VALUES ('old-run', 'Before v7', 'h:1:a', 'done', 2, NOW(6), NOW(6))")
        conn.commit()
        _db._ensure_control_schema(conn)
        _db._ensure_control_schema(conn)  # a second deploy changes nothing
        columns = {row["name"] for row in conn.execute(
            "SELECT COLUMN_NAME AS name FROM information_schema.COLUMNS"
            " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME IN ('work_jobs', 'work_lease')").fetchall()}
        version = conn.execute("SELECT MAX(version) AS v FROM control_schema_migrations").fetchone()["v"]
    finally:
        conn.close()
    for table, added in _db._CONTROL_COLUMNS.items():
        assert {name for name, _definition in added} <= columns, table
    assert version == _db.CONTROL_SCHEMA_VERSION == 7
    # The old run is still there, now a root of its own tree.
    tree = _tree(mysql_config, "old-run")
    assert tree["title"] == "Before v7" and tree["kind"] == "queue" and tree["children"] == []


def test_a_fresh_control_schema_matches_the_added_columns(control_config):
    """`control_schema.sql`'s CREATE TABLEs and `_CONTROL_COLUMNS` agree,
    so a new install and an upgraded one end up the same."""
    for table, added in _db._CONTROL_COLUMNS.items():
        rows = _rows(control_config, "SELECT COLUMN_NAME AS name FROM information_schema.COLUMNS"
                                     " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,))
        assert {name for name, _definition in added} <= {row["name"] for row in rows}


# ---------------------------------------------------------------------------
# Nodes
# ---------------------------------------------------------------------------

def test_nodes_nest_and_every_node_is_timed(control_config):
    with workQueue.job_node("plan", "Plan the galaxy", control_config, argv=["plan"], database="g1") as root:
        with workQueue.job_node("skeleton", "Galaxy skeleton") as skeleton:
            time.sleep(0.02)
        with workQueue.job_node("bright-stars", "Bright stars"):
            with workQueue.WorkQueue("Layers", workers=2, control_config=control_config) as queue:
                queue.expect(4)
                for n in range(4):
                    queue.submit("bright-stars", f"layer {n}", _slow, n)
    tree = _tree(control_config, root.id)
    shape = [(depth, node["kind"], node["title"]) for depth, node in _walk(tree)]
    assert shape == [
        (0, "plan", "Plan the galaxy"),
        (1, "skeleton", "Galaxy skeleton"),
        (1, "bright-stars", "Bright stars"),
        (2, "queue", "Layers"),
    ]
    for _depth, node in _walk(tree):
        assert node["state"] == node["status"] == "done"
        assert node["started_at"] and node["finished_at"] and node["seconds"] is not None
        assert node["root_id"] == root.id and node["database_name"] == "g1"
    assert tree["argv"] == ["plan"]
    assert tree["children"][0]["id"] == skeleton.id and tree["children"][0]["seconds"] >= 0.02
    layers = tree["children"][1]["children"][0]
    assert layers["workers"] == 2 and [task["task_key"] for task in layers["tasks"]] == [f"layer {n}" for n in range(4)]
    # Totals roll up from the queue to the root.
    for node in (tree, tree["children"][1], layers):
        assert (node["totals"]["tasks"], node["totals"]["done"], node["totals"]["queued"]) == (4, 4, 0)
        assert node["totals"]["work_seconds"] >= 4 * 0.05
        assert node["totals"]["eta_seconds"] is None
    assert tree["totals"]["started_at"] <= layers["started_at"]
    assert tree["totals"]["finished_at"] >= layers["finished_at"]


def test_one_worker_records_its_tasks_without_the_lease(control_config):
    with workQueue.job_node("galaxy", "Serial run", control_config) as root:
        with workQueue.WorkQueue("Sectors", workers=1, control_config=control_config) as queue:
            for n in range(3):
                queue.submit("sector", f"0,0,{n}", _square, n)
    queue_node = _tree(control_config, root.id)["children"][0]
    assert queue_node["totals"]["done"] == 3 and queue_node["workers"] == 1
    assert all(task["seconds"] is not None for task in queue_node["tasks"])
    assert _rows(control_config, "SELECT COUNT(*) AS n FROM work_lease WHERE holder IS NOT NULL")[0]["n"] == 0


def test_a_failing_block_fails_its_node_and_its_parents(control_config):
    with pytest.raises(ValueError):
        with workQueue.job_node("galaxy", "Breaks", control_config) as root:
            with workQueue.job_node("population", "Population pass"):
                raise ValueError("broken")
    tree = _tree(control_config, root.id)
    assert tree["state"] == "failed" and tree["children"][0]["state"] == "failed"


@pytest.mark.parametrize("raised", [KeyboardInterrupt, SystemExit])
def test_an_interrupted_block_cancels_its_node(control_config, raised):
    with pytest.raises(raised):
        with workQueue.job_node("galaxy", "Stopped", control_config) as root:
            raise raised()
    assert _tree(control_config, root.id)["state"] == "cancelled"


def test_a_process_hangs_its_root_under_the_parent_it_is_given(control_config, monkeypatch):
    step = workQueue.open_node("step", "Step 1", control_config)
    try:
        # What a step's own process sees: its node is open only there.
        workQueue._open_nodes.remove(step)
        monkeypatch.setenv(workQueue.PARENT_ENV_VAR, step.id)
        with workQueue.job_node("galaxy", "Child run", control_config) as child:
            pass
        assert child.parent_id == step.id and child.root_id == step.root_id
        monkeypatch.setenv(workQueue.PARENT_ENV_VAR, "20260101-000000-deadbeef")
        with workQueue.job_node("galaxy", "Orphan", control_config) as orphan:
            pass
        assert orphan.parent_id is None and orphan.root_id == orphan.id
    finally:
        workQueue.close_node(step, "done")
    assert [node["title"] for node in _tree(control_config, step.id)["children"]] == ["Child run"]


def test_without_a_control_database_nodes_cost_nothing(mysql_config):
    with workQueue.job_node("galaxy", "Nowhere", mysql_config) as root:
        with workQueue.WorkQueue("Sectors", workers=1, control_config=mysql_config) as queue:
            queue.submit("sector", "k", _square, 3)
    assert not root.store.available
    assert queue.finished == 1


def test_a_dead_run_reads_as_interrupted_and_a_live_one_has_an_eta(control_config):
    conn = _conn(control_config)
    try:
        with conn:
            conn.execute("INSERT INTO work_jobs (id, root_id, kind, title, holder, state, workers, tasks_total,"
                         " created_at, started_at, heartbeat_at) VALUES ('dead', 'dead', 'galaxy', 'Died', 'h:1',"
                         " 'running', 0, NULL, NOW(6) - INTERVAL 1 HOUR, NOW(6) - INTERVAL 1 HOUR,"
                         " NOW(6) - INTERVAL 10 MINUTE)")
            conn.execute("INSERT INTO work_jobs (id, root_id, kind, title, holder, state, workers, tasks_total,"
                         " created_at, started_at, heartbeat_at) VALUES ('live', 'live', 'queue', 'Running', 'h:2',"
                         " 'running', 2, 10, NOW(6), NOW(6), NOW(6))")
            for n in range(4):
                conn.execute("INSERT INTO work_tasks (job_id, kind, task_key, state, created_at, seconds)"
                             " VALUES ('live', 'sector', ?, 'done', NOW(6), 3)", (str(n),))
            conn.execute("INSERT INTO work_tasks (job_id, kind, task_key, state, created_at)"
                         " VALUES ('live', 'sector', '4', 'running', NOW(6))")
    finally:
        conn.close()
    dead = _tree(control_config, "dead")
    assert (dead["state"], dead["status"], dead["live"]) == ("running", "interrupted", False)
    assert dead["totals"]["eta_seconds"] is None
    live = _tree(control_config, "live")
    totals = live["totals"]
    assert (totals["tasks"], totals["done"], totals["running"], totals["queued"]) == (10, 4, 1, 5)
    # 6 left at 3 s each on 2 workers.
    assert totals["eta_seconds"] == pytest.approx(9.0)


def test_old_finished_trees_are_pruned_with_their_subtrees(control_config):
    conn = _conn(control_config)
    try:
        with conn:
            for node_id, parent, state in (("old", None, "done"), ("old-child", "old", "done"),
                                           ("old-live", None, "running")):
                conn.execute("INSERT INTO work_jobs (id, parent_id, root_id, title, holder, state, workers,"
                             " created_at, heartbeat_at) VALUES (?, ?, 'old', 'x', 'h', ?, 0,"
                             " NOW(6) - INTERVAL 30 DAY, NOW(6))", (node_id, parent, state))
            conn.execute("INSERT INTO work_tasks (job_id, kind, task_key, state, created_at)"
                         " VALUES ('old-child', 'sector', 'k', 'done', NOW(6) - INTERVAL 30 DAY)")
    finally:
        conn.close()
    with workQueue.job_node("galaxy", "New run", control_config):
        pass
    left = {row["id"] for row in _rows(control_config, "SELECT id FROM work_jobs")}
    assert "old" not in left and "old-child" not in left and "old-live" in left
    assert _rows(control_config, "SELECT COUNT(*) AS n FROM work_tasks WHERE job_id = 'old-child'")[0]["n"] == 0


def test_timing_by_kind(control_config):
    with workQueue.job_node("galaxy", "Timed", control_config):
        with workQueue.WorkQueue("Sectors", workers=1, control_config=control_config) as queue:
            for n in range(3):
                queue.submit("sector", str(n), _slow, n)
    conn = _conn(control_config)
    try:
        timing = workQueue.timing_by_kind(conn)
    finally:
        conn.close()
    assert timing["tasks"]["sector"]["count"] == 3
    assert timing["tasks"]["sector"]["mean_seconds"] >= 0.05
    assert {"galaxy", "queue"} <= set(timing["nodes"])


# ---------------------------------------------------------------------------
# generate.py and the job runner
# ---------------------------------------------------------------------------

def test_the_run_argv_leaves_out_the_database_and_debug_options():
    argv = ["galaxy", "--ring", "3", "--mysql-password", "secret", "--debug", "f.log", "--layer=1",
            "--mysql-user=u", "--mysql-database", "db", "--debug", "--yes"]
    assert generate._run_argv(argv) == ["galaxy", "--ring", "3", "--layer=1", "--yes"]


def _run(argv):
    old_argv = sys.argv
    try:
        sys.argv = ["generate.py"] + argv
        generate.main()
    finally:
        sys.argv = old_argv


def test_a_galaxy_run_is_a_tree_down_to_its_sectors(control_config):
    _plan_wide_galaxy(control_config)
    _run(["galaxy", "--ring", "1", "--num-systems", "2"] + _mysql_argv(control_config))
    root_id = _rows(control_config, "SELECT id FROM work_jobs WHERE parent_id IS NULL")[0]["id"]
    tree = _tree(control_config, root_id)
    assert tree["kind"] == "galaxy" and tree["state"] == "done"
    assert tree["argv"] == ["galaxy", "--ring", "1", "--num-systems", "2"]
    assert "secret" not in tree["title"] and control_config.database not in tree["title"]
    assert tree["database_name"] == control_config.database
    queue = tree["children"][0]
    expected = generate.ring_sector_count(1)
    assert queue["kind"] == "queue" and queue["tasks_total"] == expected
    assert tree["totals"]["done"] == expected and {task["kind"] for task in queue["tasks"]} == {"sector"}


def test_a_web_job_is_the_root_of_its_steps_runs(control_config, tmp_path, monkeypatch):
    jobs_dir = tmp_path / "jobs"
    job_dir = jobs_dir / "20261001-120000-abcd"
    job_dir.mkdir(parents=True)
    src = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    child = (
        "import sys; sys.path.insert(0, %r)\n"
        "from stellarObjects import _db, workQueue\n"
        "with workQueue.job_node('galaxy', 'Inside the step', _db.control_mysql_config()):\n"
        "    pass\n"
    ) % src
    env = {
        "PLANETGEN_MYSQL_HOST": control_config.host, "PLANETGEN_MYSQL_PORT": str(control_config.port),
        "PLANETGEN_MYSQL_USER": control_config.user, "PLANETGEN_MYSQL_PASSWORD": control_config.password or "",
        "PLANETGEN_MYSQL_DATABASE": control_config.database,
    }
    import json
    (job_dir / "job.json").write_text(json.dumps({
        "id": job_dir.name, "kind": "galaxy", "title": "Generate sectors", "database": control_config.database,
        "created_at": time.time(), "cwd": src, "env": env,
        "steps": [{"label": "Step one", "argv": [sys.executable, "-c", child]},
                  {"label": "Step two", "argv": [sys.executable, "-c", "pass"]}],
    }))
    assert jobRunner.run(str(job_dir)) == 0
    root = _rows(control_config, "SELECT id FROM work_jobs WHERE kind = 'web-job'")[0]["id"]
    tree = _tree(control_config, root)
    assert (tree["title"], tree["web_job_id"], tree["state"]) == ("Generate sectors", job_dir.name, "done")
    assert [(step["kind"], step["title"], step["state"]) for step in tree["children"]] == [
        ("step", "Step one", "done"), ("step", "Step two", "done")]
    inside = tree["children"][0]["children"]
    assert [(node["title"], node["web_job_id"]) for node in inside] == [("Inside the step", None)]
    assert all(node["seconds"] is not None for _depth, node in _walk(tree))
