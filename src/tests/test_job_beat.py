"""OPS.142: the API jobs, the web jobs and the tile prefetch chunks give a beat while they run, and the Queue page lists them with a state."""

import json
import sys

import pytest

from planetgen.queue import api_jobs, allowlist, beatjudge, jobbeat, redisqueue, unitbeat
from planetgen.web import job_runner
from planetgen.web.lib import tileprefetch, tilecache


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


def _listed(conn):
    return {(row["kind"], row["job"]): row for row in jobbeat.states(conn)}


def test_a_job_has_a_record_while_it_runs_and_none_after(conn):
    with jobbeat.beating("api", "abc", label="regenerate_sector", conn=conn) as beat:
        beat.phase("step 2 of 4")
        beat.write()
        row = _listed(conn)[("api", "abc")]
        assert row["state"] == beatjudge.WORKING and row["label"] == "regenerate_sector" and row["phase"] == "step 2 of 4"
    assert jobbeat.states(conn) == []


def test_a_job_that_died_reads_dead_and_is_dropped_after_an_hour(conn):
    with jobbeat.beating("web", "w1", conn=conn) as beat:
        beat.write()
        raw = unitbeat.parse(conn.hget(jobbeat.hash_key(), "web:w1"))
        now = unitbeat.server_time(conn)
        assert jobbeat.states(conn, now=now + 120.0)[0]["state"] == beatjudge.DEAD, "no beat for 2 minutes"
        assert jobbeat.states(conn, now=now + jobbeat.FORGET_SECONDS + 60.0) == []
        assert conn.hget(jobbeat.hash_key(), "web:w1") is None, "the dead record is dropped from the hash"
        assert raw["kind"] == "web"


def test_a_job_runs_without_a_beat_when_redis_does_not_answer():
    class _Broken:
        def __getattr__(self, name):
            raise ConnectionError("down")

    with jobbeat.beating("api", "x", conn=_Broken()) as beat:
        beat.phase("anything")
        with beat.waiting("mysql"):
            beat(1, 2)




def test_an_api_job_beats_while_its_function_runs(conn, monkeypatch):
    seen = {}

    def function(value):
        beat = unitbeat._current
        beat.write()
        seen["rows"] = jobbeat.states(conn)
        return value * 2

    monkeypatch.setattr(allowlist, "check_reference", lambda reference: "tests.test_job_beat:_unused")
    monkeypatch.setitem(sys.modules, "tests.test_job_beat", type("M", (), {"_unused": staticmethod(function)}))
    assert api_jobs.execute("ignored", [api_jobs.encode_argument(21)]) == {"result": 42}
    assert [(row["kind"], row["label"], row["state"]) for row in seen["rows"]] == [("api", "_unused", "working")]
    assert jobbeat.states(conn) == []


def test_a_tile_prefetch_chunk_beats_while_it_builds(conn, monkeypatch):
    seen = {}
    from planetgen.web import warmup

    class _App:
        def test_request_context(self, path):
            import contextlib
            return contextlib.nullcontext()

    def fetch_nodes(db, keys, ahead=False):
        unitbeat._current.write()
        seen["rows"] = jobbeat.states(conn)
        return {"cached": 1}

    monkeypatch.setattr(warmup, "build_app", lambda: _App())
    monkeypatch.setattr(tilecache, "fetch_nodes", fetch_nodes)
    assert tileprefetch.build_nodes("db", ["a", "b", "c"]) == {"nodes": 3, "built": 2}
    assert [(row["kind"], row["label"]) for row in seen["rows"]] == [("tile-prefetch", "3 nodes of db")]
    assert jobbeat.states(conn) == []


def test_a_web_job_beats_for_the_whole_run_and_names_each_step(tmp_path, monkeypatch):
    job_dir = tmp_path / "jobs" / "job-9"
    job_dir.mkdir(parents=True)
    job = {"id": "job-9", "kind": "galaxy", "title": "Two steps", "cwd": str(tmp_path), "env": {},
           "steps": [{"label": "One", "argv": [sys.executable, "-c", "pass"]}, {"label": "Two", "argv": [sys.executable, "-c", "pass"]}]}
    (job_dir / "job.json").write_text(json.dumps(job), encoding="utf-8")
    spy = []

    class _Beat:
        def phase(self, text):
            spy.append(("phase", text))

    from contextlib import contextmanager

    @contextmanager
    def beating(kind, job_id, label="", **kwargs):
        spy.append(("beating", kind, job_id, label))
        yield _Beat()

    monkeypatch.setattr(jobbeat, "beating", beating)
    assert job_runner.run(str(job_dir)) == 0
    assert spy == [("beating", "web", "job-9", "Two steps"), ("phase", "step 1 of 2: One"), ("phase", "step 2 of 2: Two")]
