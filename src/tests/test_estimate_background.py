# tests/test_estimate_background.py

"""
ADM.141: the Generate page never works its size estimate out inside a web request. A big galaxy takes minutes, Apache
restarts the daemon of a request that outlives `request-timeout` (60 s), and the live size line asked again on every
change, so the command ran over and over. A request now answers from the cache or queues the work and answers
"working" at once; one command runs for a question, and for a galaxy database, until it ends.
"""

import json
import sys
import textwrap
import uuid

import pytest
import redis

from planetgen.web import generate_page
from planetgen.web.lib import estimatecache, estimatejob, singleflight

from tests.test_web_generate import (  # noqa: F401 -- fixtures
    DB, _ESTIMATE, _post, client, jobs_root, no_spawn, site,
)

REAL_START = estimatejob.start
"""The real `start`, which the suite's autouse fixture (`conftest._estimates_run_in_the_request`) replaces."""

ARGV = ["planetgen", "sectors", "--center-sector", "5", "--radius-pc", "100"]
VERSIONS = {"generator_epoch": 10, "data": ["a", {"_n": 1, "epoch": 0}], "speeds": "2@t/1@t"}
ESTIMATE = {"sectors": 9, "systems": 120, "what": "the sectors"}


@pytest.fixture
def database(redis_server):
    name = f"estimate_background_{uuid.uuid4().hex[:10]}"
    yield name
    client = redis.Redis.from_url(redis_server)
    for pattern in (f"planetgen:estimate:{name}:*", f"planetgen:build:planetgen:estimate:{name}:*",
                    f"planetgen:build:{estimatecache.GATE_PREFIX}{name}"):
        for key in client.scan_iter(match=pattern):
            client.delete(key)


class _Queue:
    """Stands for the work queue: records what was queued and does not run it."""

    def __init__(self, accepts=True):
        self.accepts = accepts
        self.started = []

    def __call__(self, key, versions):
        self.started.append((key, versions))
        return self.accepts


def _ask(database, queue, argv=ARGV, versions=VERSIONS, compute=None):
    return estimatecache.estimate_in_background(
        database, "galaxy", argv, lambda: versions, compute or (lambda: pytest.fail("worked out in the request")), queue)


def test_a_missing_estimate_is_queued_and_the_request_returns_at_once(database):
    queue = _Queue()
    with pytest.raises(estimatecache.Working):
        _ask(database, queue)
    assert len(queue.started) == 1


def test_asking_again_while_it_runs_does_not_start_the_command_again(database):
    queue = _Queue()
    for _ in range(5):
        with pytest.raises(estimatecache.Working):
            _ask(database, queue)
    assert len(queue.started) == 1


def test_another_form_of_the_same_galaxy_waits_for_the_running_command(database):
    """A form typed into asks for every value on the way; one command runs at a time for the galaxy."""
    queue = _Queue()
    with pytest.raises(estimatecache.Working):
        _ask(database, queue)
    with pytest.raises(estimatecache.Working):
        _ask(database, queue, argv=ARGV[:-1] + ["1000"])
    assert len(queue.started) == 1


def test_the_estimate_of_another_galaxy_is_not_held_back(database):
    queue = _Queue()
    with pytest.raises(estimatecache.Working):
        _ask(database, queue)
    with pytest.raises(estimatecache.Working):
        _ask(database + "_other", queue)
    assert len(queue.started) == 2
    other = database + "_other"
    singleflight.release([estimatecache.cache_key(other, "galaxy", ARGV), estimatecache.GATE_PREFIX + other])


def test_a_stored_answer_is_served_and_gives_up_the_claims(database):
    queue = _Queue()
    with pytest.raises(estimatecache.Working):
        _ask(database, queue)
    key, versions = queue.started[0]
    estimatecache.store(key, versions, estimate=ESTIMATE)
    assert _ask(database, queue) == ESTIMATE
    # The claims are gone: a change of the versions queues the command again.
    with pytest.raises(estimatecache.Working):
        _ask(database, queue, versions={**VERSIONS, "speeds": "3@t/1@t"})
    assert len(queue.started) == 2


def test_a_failed_estimate_says_why_for_a_moment_and_then_asks_again(database):
    queue = _Queue()
    with pytest.raises(estimatecache.Working):
        _ask(database, queue)
    key, versions = queue.started[0]
    estimatecache.store(key, versions, error="Nothing was generated: the shell is too large.")
    with pytest.raises(estimatecache.Failed, match="shell is too large"):
        _ask(database, queue)
    assert len(queue.started) == 1
    redis.Redis.from_url(redis_url()).delete(key)          # its time is up
    with pytest.raises(estimatecache.Working):
        _ask(database, queue)
    assert len(queue.started) == 2


def redis_url():
    from planetgen.queue import redisqueue
    return redisqueue.redis_url()


def test_with_no_queue_the_estimate_is_worked_out_in_the_request_and_cached(database):
    queue = _Queue(accepts=False)
    runs = []
    first = _ask(database, queue, compute=lambda: runs.append(1) or dict(ESTIMATE))
    second = _ask(database, queue, compute=lambda: pytest.fail("cached"))
    assert first == second == ESTIMATE and runs == [1]


def test_a_queue_that_refuses_falls_back_to_the_request(database):
    def broken(key, versions):
        raise RuntimeError("no Redis server")

    assert _ask(database, broken, compute=lambda: dict(ESTIMATE)) == ESTIMATE


def test_versions_that_cannot_be_read_run_the_estimate_without_queueing(database):
    def unreadable():
        raise OSError("control database is down")

    queue = _Queue()
    got = estimatecache.estimate_in_background(database, "galaxy", ARGV, unreadable, lambda: dict(ESTIMATE), queue)
    assert got == ESTIMATE and queue.started == []


def test_the_claim_outlives_the_command():
    assert estimatecache.CLAIM_SECONDS > estimatejob.TIMEOUT_SECONDS


# ---- the queued body -------------------------------------------------------------------------------------------------

def _script(tmp_path, body):
    path = tmp_path / "estimate.py"
    path.write_text(textwrap.dedent(body))
    return [sys.executable, str(path)]


def test_the_worker_stores_the_estimate_the_command_prints(database, tmp_path):
    queue = _Queue()
    argv = _script(tmp_path, '''
        import json, sys
        assert sys.argv[-1] == "--estimate-only"
        print("Estimate for the sectors: 9 sectors")
        print("ESTIMATE " + json.dumps({"sectors": 9, "what": "the sectors"}))
    ''')
    with pytest.raises(estimatecache.Working):
        _ask(database, queue, argv=argv)
    key, versions = queue.started[0]
    assert estimatejob.work(key, versions, argv, {}, str(tmp_path)) == {"stored": "estimate"}
    assert _ask(database, queue, argv=argv) == {"sectors": 9, "what": "the sectors"}


def test_the_worker_stores_why_the_command_made_no_estimate(database, tmp_path):
    queue = _Queue()
    argv = _script(tmp_path, '''
        import sys
        print("shell 900 is too large to generate without a limit", file=sys.stderr)
        sys.exit(1)
    ''')
    with pytest.raises(estimatecache.Working):
        _ask(database, queue, argv=argv)
    key, versions = queue.started[0]
    assert estimatejob.work(key, versions, argv, {}, str(tmp_path)) == {"stored": "error"}
    with pytest.raises(estimatecache.Failed, match="too large to generate without a limit"):
        _ask(database, queue, argv=argv)


def test_the_worker_stores_a_timeout_and_gives_up_the_claims(database, tmp_path):
    queue = _Queue()
    argv = _script(tmp_path, "import time\ntime.sleep(30)\n")
    with pytest.raises(estimatecache.Working):
        _ask(database, queue, argv=argv)
    key, versions = queue.started[0]
    assert estimatejob.work(key, versions, argv, {}, str(tmp_path), timeout=1) == {"stored": "error"}
    with pytest.raises(estimatecache.Failed, match="took too long"):
        _ask(database, queue, argv=argv)
    owned, held = singleflight.claim([key, estimatecache._gate_of(key)])
    assert held == [] and len(owned) == 2
    singleflight.release(owned)


def test_a_worker_that_crashes_stores_the_failure_and_gives_up_the_claims(database, monkeypatch):
    queue = _Queue()
    with pytest.raises(estimatecache.Working):
        _ask(database, queue)
    key, versions = queue.started[0]

    def boom(*args, **kwargs):
        raise MemoryError("out of memory")

    monkeypatch.setattr(estimatejob.api_jobs, "run_command", boom)
    with pytest.raises(MemoryError):
        estimatejob.work(key, versions, ARGV, {}, ".")
    with pytest.raises(estimatecache.Failed, match="out of memory"):
        _ask(database, queue)


@pytest.mark.parametrize("done, message", [
    ({"timed_out": True, "error": None, "stdout": "", "stderr": "", "returncode": None}, "took too long"),
    ({"timed_out": False, "error": "no such file", "stdout": "", "stderr": "", "returncode": None},
     "couldn't be worked out: no such file"),
    ({"timed_out": False, "error": None, "stdout": "", "stderr": "", "returncode": 3}, "exited with 3"),
    ({"timed_out": False, "error": None, "stdout": "ESTIMATE {not json", "stderr": "bad", "returncode": 0},
     "Nothing was generated: bad"),
])
def test_read_says_why_there_is_no_estimate(done, message):
    with pytest.raises(estimatejob.EstimateError, match=message):
        estimatejob.read(done)


def test_start_queues_the_worker_by_name(monkeypatch):
    monkeypatch.setattr(estimatejob, "start", REAL_START)
    queued = []
    monkeypatch.setattr(estimatejob.api_jobs, "submit", lambda function, *args, **kw: queued.append((function, args)))
    assert estimatejob.start("k", {"v": 1}, ["a", "b"], {"X": "1"}, "/repo") is True
    function, args = queued[0]
    assert function is estimatejob.work and args == ("k", {"v": 1}, ["a", "b"], {"X": "1"}, "/repo")

    def refuse(*args, **kwargs):
        raise estimatejob.api_jobs.NoQueue("no Redis server")

    monkeypatch.setattr(estimatejob.api_jobs, "submit", refuse)
    assert estimatejob.start("k", {}, [], {}, ".") is False


# ---- the Generate page -------------------------------------------------------------------------------------------------

SIZE = "/admin/generate/task/add-sectors/size"
SHELL_40 = {"what": "grid", "mode": "shell", "shell_index": "40", "shell_limit": "200"}


@pytest.fixture
def queued(site, monkeypatch, redis_server):
    """The page with a work queue that records what it is given, a readable cache and an estimate that must not run
    in a request."""
    queue = _Queue()
    monkeypatch.setattr(estimatejob, "start", lambda key, versions, argv, env, cwd: queue(key, versions))
    monkeypatch.setattr(generate_page, "_estimate_versions", lambda name: dict(VERSIONS))
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: pytest.fail("worked out in the request"))
    server = redis.Redis.from_url(redis_server)

    def clean():
        for pattern in (f"planetgen:estimate:{DB}:*", f"planetgen:build:planetgen:estimate:{DB}:*",
                        f"planetgen:build:{estimatecache.GATE_PREFIX}{DB}"):
            for key in server.scan_iter(match=pattern):
                server.delete(key)

    clean()
    yield queue
    clean()


def test_the_size_line_answers_working_at_once_and_ready_when_the_worker_is_done(client, queued, no_spawn):
    answer = client.get(SIZE, query_string=SHELL_40).get_json()
    assert answer["state"] == "working" and answer["retry_seconds"] == generate_page.ESTIMATE_RETRY_SECONDS
    assert "Working out how big" in answer["message"]
    for _ in range(3):
        assert client.get(SIZE, query_string=SHELL_40).get_json()["state"] == "working"
    assert len(queued.started) == 1
    key, versions = queued.started[0]
    estimatecache.store(key, versions, estimate=dict(_ESTIMATE))
    ready = client.get(SIZE, query_string=SHELL_40).get_json()
    assert ready["state"] == "ready" and ready["estimate"]["sectors"] == 9
    assert no_spawn == []


def test_the_review_step_says_it_is_working_and_keeps_the_form(client, queued, no_spawn):
    resp = _post(client, action="galaxy", mode="shell", shell_index="40")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 202 and "Working out how big this would be" in html
    assert 'name="shell_index"' in html and 'data-size-report' in html
    assert no_spawn == [] and len(queued.started) == 1


def test_a_script_posting_the_form_is_told_to_ask_again(client, queued, no_spawn):
    resp = _post(client, headers={"Accept": "application/json"}, action="galaxy", mode="shell", shell_index="40")
    assert resp.status_code == 202 and "Working out how big" in resp.get_json()["error"]
    assert no_spawn == []


def test_the_review_step_shows_the_estimate_once_it_is_stored(client, queued, no_spawn):
    _post(client, action="galaxy", mode="shell", shell_index="40")
    key, versions = queued.started[0]
    estimatecache.store(key, versions, estimate=dict(_ESTIMATE))
    resp = _post(client, action="galaxy", mode="shell", shell_index="40")
    assert resp.status_code == 200 and "About 8.4 MB" in resp.get_data(as_text=True)


def test_a_failed_estimate_is_shown_on_the_form(client, queued, no_spawn):
    _post(client, action="galaxy", mode="shell", shell_index="40")
    key, versions = queued.started[0]
    estimatecache.store(key, versions, error="Nothing was generated: shell 40 is beyond the plan.")
    resp = _post(client, action="galaxy", mode="shell", shell_index="40")
    assert resp.status_code == 400 and "beyond the plan" in resp.get_data(as_text=True)
    size = client.get(SIZE, query_string={"action": "galaxy", "mode": "shell", "shell_index": "40"}).get_json()
    assert size == {"state": "invalid", "message": "Nothing was generated: shell 40 is beyond the plan."}


def test_a_confirmed_form_starts_the_job_without_asking_for_the_estimate(client, queued, no_spawn):
    """Pressing Run on the Review step must not wait for, or start, the command."""
    resp = _post(client, action="galaxy", estimate_ok="1", mode="shell", shell_index="40")
    assert resp.status_code == 303 and len(no_spawn) == 1 and queued.started == []


def test_the_answers_keep_their_headers_small(client, queued):
    """The daemon's answer to a request that took too long was 'truncated or oversized response headers'; the page's
    own answers while it works carry only the usual few headers, nothing of the estimate."""
    for resp in (client.get(SIZE, query_string=SHELL_40),
                 _post(client, action="galaxy", mode="shell", shell_index="40")):
        size = sum(len(name) + len(value) + 4 for name, value in resp.headers)
        assert size < 2048, dict(resp.headers)


# --- Through a real Redis and the guarded worker (the allowlist let nothing else catch ADM.141's missing entry) ----------

def test_the_estimate_runs_end_to_end_through_the_queue_and_the_guarded_worker(database, redis_server, monkeypatch):
    """The request queues `estimatejob.work` with `api_jobs.submit`; a worker built as the site builds it (the allowlist's `GuardedJob`) runs it, and the next ask reads the estimate."""
    import rq

    from planetgen.queue import allowlist, api_jobs, redisqueue, workers

    name = f"planetgen-estimate-e2e-{uuid.uuid4().hex[:8]}"
    monkeypatch.setattr(api_jobs, "lane", lambda expected_seconds=None: (name, 1))
    monkeypatch.setattr(workers, "ensure_workers", lambda *args, **kwargs: None)       # this test's worker is the one below
    ran = []

    def command(argv, env, cwd, timeout, merge_stderr=False):
        ran.append(argv)
        return {"returncode": 0, "stdout": "noise\nESTIMATE " + json.dumps(ESTIMATE) + "\n", "stderr": "", "timed_out": False, "error": None}

    monkeypatch.setattr(api_jobs, "run_command", command)
    with pytest.raises(estimatecache.Working):
        _ask(database, lambda key, versions: REAL_START(key, versions, ARGV, {"PLANETGEN_MYSQL_DATABASE": database}, "."))
    connection = redisqueue.connect(redis_server)
    try:
        queue = redisqueue.queue(name, connection)
        assert queue.count == 1
        rq.SimpleWorker([queue], connection=connection, serializer=redisqueue.serializer_for(name), job_class=allowlist.GuardedJob).work(burst=True, logging_level="CRITICAL")
        assert ran == [ARGV + ["--estimate-only"]]
        assert _ask(database, _Queue()) == ESTIMATE
    finally:
        queue.empty()
        for job_id in queue.failed_job_registry.get_job_ids():
            queue.failed_job_registry.remove(job_id, delete_job=True)
