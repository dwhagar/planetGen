# tests/test_web_jobs_damage.py

"""
The Generate page's background jobs (`web/jobs.py`) when things go
wrong around them: two admins starting a job at the same moment
(TEST.40), damaged job files, a lock holding garbage, a job id drawn
twice, an unwritable jobs directory, pruning while a job runs and
cancelling ids that aren't running jobs (TEST.41), and a job outliving
the web server process that started it (ADM.11).

No database: the steps are tiny Python one-liners, like
`test_web_generate.py`'s job tests.
"""

import json
import os
import signal
import subprocess
import sys
import threading
import time

import pytest

from planetgen.web import jobs

PY = sys.executable
REPO = jobs.REPO_DIR


@pytest.fixture
def jobs_root(tmp_path, monkeypatch):
    root = tmp_path / "jobs"
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(root))
    return str(root)


def _step(label, code):
    return {"label": label, "argv": [PY, "-c", code]}


def _start(title="Job", spawn=False, root=None):
    return jobs.start_job("reset", title, [_step("x", "pass")], spawn=spawn, root=root)


def _wait_finished(job_id, root, timeout=30):
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = jobs.get_job(job_id, root)
        if job and job["finished"]:
            return job
        time.sleep(0.05)
    raise AssertionError(f"job {job_id} did not finish: {jobs.get_job(job_id, root)}")


def _age(path, seconds):
    then = time.time() - seconds
    os.utime(path, (then, then))


def _backdate_job(root, job_id, seconds=3600):
    """Pretend the job was spawned long ago by a runner that's gone."""
    path = os.path.join(root, job_id, "job.json")
    with open(path) as f:
        body = json.load(f)
    body["created_at"] -= seconds
    with open(path, "w") as f:
        json.dump(body, f)


# --- TEST.40 Two admins start a job at once ----------------------------------------

def test_two_threads_starting_at_once_start_one_job(jobs_root):
    for _round in range(25):
        barrier = threading.Barrier(2)
        started, busy = [], []

        def start():
            barrier.wait()
            try:
                started.append(_start(root=jobs_root))
            except jobs.JobBusy:
                busy.append(True)

        threads = [threading.Thread(target=start) for _ in range(2)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        assert len(started) == 1 and len(busy) == 1
        assert jobs.active_job(jobs_root)["id"] == started[0]
        _backdate_job(jobs_root, started[0])  # free the lock for the next round


STARTER = r"""
import os, sys, time
from planetgen.web import jobs
go = sys.argv[1]
while not os.path.exists(go):
    time.sleep(0.001)
try:
    print(jobs.start_job("reset", "Race", [{"label": "x", "argv": ["true"]}], spawn=False))
except jobs.JobBusy:
    print("busy")
"""


def test_many_processes_starting_at_once_start_one_job(jobs_root, tmp_path):
    """Separate processes, as two Apache workers would be."""
    os.makedirs(jobs_root)
    go = str(tmp_path / "go")
    env = dict(os.environ, PLANETGEN_JOBS_DIR=jobs_root)
    procs = [subprocess.Popen([PY, "-c", STARTER, go], env=env, cwd=REPO,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
             for _ in range(6)]
    time.sleep(1.0)  # let them all import and wait
    open(go, "w").close()
    outputs = [proc.communicate(timeout=60) for proc in procs]
    lines = [out.strip() for out, _err in outputs]
    assert all(proc.returncode == 0 for proc in procs), [err for _out, err in outputs]
    started = [line for line in lines if line != "busy"]
    assert len(started) == 1, lines
    assert jobs.active_job(jobs_root)["id"] == started[0]
    # The losers left no job directories behind.
    assert [name for name in os.listdir(jobs_root) if jobs.JOB_ID_RE.match(name)] == started


def test_the_lock_is_never_seen_empty(jobs_root, monkeypatch):
    """The lock appears with the job id already in it: a caller reading it
    between creating and writing used to see an empty lock, judge it
    stale and start a second job."""
    seen = []
    real_link = os.link

    def watching_link(src, dst):
        real_link(src, dst)
        with open(dst) as f:
            seen.append(f.read())

    monkeypatch.setattr(jobs.os, "link", watching_link)
    job_id = _start(root=jobs_root)
    assert seen == [job_id]


def test_lock_works_without_hard_links(jobs_root, monkeypatch):
    """Some filesystems (FAT, some network shares) have no hard links."""
    def no_link(src, dst):
        raise OSError(1, "Operation not permitted")

    monkeypatch.setattr(jobs.os, "link", no_link)
    first = _start(root=jobs_root)
    with pytest.raises(jobs.JobBusy):
        _start(root=jobs_root)
    assert jobs.active_job(jobs_root)["id"] == first
    assert not [name for name in os.listdir(jobs_root) if name.startswith(".active-")]


def test_a_stale_lock_is_cleared_by_one_caller_only(jobs_root):
    first = _start(root=jobs_root)
    _backdate_job(jobs_root, first)
    clearing = os.path.join(jobs_root, jobs.CLEARING_NAME)
    open(clearing, "w").close()  # someone else is clearing it right now
    with pytest.raises(jobs.JobBusy):
        _start(root=jobs_root)
    os.remove(clearing)
    second = _start(root=jobs_root)
    assert jobs.active_job(jobs_root)["id"] == second


def test_a_clearing_mark_left_by_a_dead_process_expires(jobs_root):
    first = _start(root=jobs_root)
    _backdate_job(jobs_root, first)
    clearing = os.path.join(jobs_root, jobs.CLEARING_NAME)
    open(clearing, "w").close()
    _age(clearing, jobs.STARTING_GRACE_SECONDS + 5)
    with pytest.raises(jobs.JobBusy):
        _start(root=jobs_root)  # removes the old mark
    assert not os.path.exists(clearing)
    _start(root=jobs_root)


# --- TEST.41 Job files damaged -------------------------------------------------------

@pytest.mark.parametrize("content", ["", "{", "[1, 2]x", "\x00\xff garbage"])
def test_corrupt_job_json_is_an_unknown_job(jobs_root, content):
    job_id = _start(root=jobs_root)
    with open(os.path.join(jobs_root, job_id, "job.json"), "w", encoding="utf-8") as f:
        f.write(content)
    assert jobs.get_job(job_id, jobs_root) is None
    assert jobs.list_jobs(root=jobs_root) == []
    assert jobs.active_job(jobs_root) is None
    assert jobs.cancel_job(job_id, jobs_root) is False
    # Its lock still blocks a new job while it might be a job being
    # written, then counts as stale.
    with pytest.raises(jobs.JobBusy) as busy:
        _start(root=jobs_root)
    assert busy.value.job is None
    _age(os.path.join(jobs_root, jobs.LOCK_NAME), jobs.STARTING_GRACE_SECONDS + 5)
    second = _start(root=jobs_root)
    assert jobs.active_job(jobs_root)["id"] == second


@pytest.mark.parametrize("content", ["", "{\"status\": \"runn", "null", "not json"])
def test_corrupt_state_json_reads_as_not_started(jobs_root, content):
    job_id = _start(root=jobs_root)
    with open(os.path.join(jobs_root, job_id, "state.json"), "w", encoding="utf-8") as f:
        f.write(content)
    assert jobs.get_job(job_id, jobs_root)["status"] == "starting"
    _backdate_job(jobs_root, job_id)
    job = jobs.get_job(job_id, jobs_root)
    assert job["status"] == "interrupted" and job["finished"]
    assert "stopped without finishing" in job["error"]


def test_corrupt_progress_json_shows_no_progress(jobs_root, monkeypatch):
    monkeypatch.setattr(jobs, "_runner_alive", lambda *args: True)
    job_id = _start(root=jobs_root)
    path = os.path.join(jobs_root, job_id)
    with open(os.path.join(path, "state.json"), "w") as f:
        json.dump({"status": "running", "step": 1, "pid": os.getpid(), "started_at": time.time()}, f)
    with open(os.path.join(path, "progress.json"), "w") as f:
        f.write('{"completed": 3, "tot')
    job = jobs.get_job(job_id, jobs_root)
    assert job["status"] == "running" and job["progress"] is None


@pytest.mark.parametrize("garbage", ["", "garbage", "../../etc/passwd", "20990101-000000-abcd", "\x00" * 50])
def test_a_lock_holding_garbage(jobs_root, garbage):
    os.makedirs(jobs_root)
    lock = os.path.join(jobs_root, jobs.LOCK_NAME)
    with open(lock, "w", encoding="utf-8") as f:
        f.write(garbage)
    assert jobs.active_job(jobs_root) is None
    # Fresh: maybe a job being written right now, so wait.
    with pytest.raises(jobs.JobBusy):
        _start(root=jobs_root)
    # Old: nothing will ever finish writing it.
    _age(lock, jobs.STARTING_GRACE_SECONDS + 5)
    job_id = _start(root=jobs_root)
    assert jobs.active_job(jobs_root)["id"] == job_id
    assert not os.path.exists(os.path.join(os.path.dirname(jobs_root), "etc"))


def test_a_job_id_drawn_twice_gets_a_fresh_one(jobs_root, monkeypatch):
    first = _start(root=jobs_root)
    _backdate_job(jobs_root, first)
    ids = iter([first, first, "20260101-000000-beef"])
    monkeypatch.setattr(jobs, "new_job_id", lambda: next(ids))
    second = _start(root=jobs_root)
    assert second == "20260101-000000-beef"
    # The first job's files are untouched.
    assert jobs.get_job(first, jobs_root)["status"] == "interrupted"


def test_an_unwritable_jobs_directory(tmp_path, monkeypatch):
    blocker = tmp_path / "not-a-dir"
    blocker.write_text("x")
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(blocker / "jobs"))
    with pytest.raises(OSError):
        jobs.jobs_dir()
    with pytest.raises(OSError):
        _start()
    # Given straight to the functions, nothing breaks.
    root = str(blocker / "jobs")
    with pytest.raises(OSError):
        _start(root=root)
    assert jobs.list_jobs(root=root) == []
    assert jobs.active_job(root) is None
    assert jobs.get_job("20260101-000000-abcd", root) is None
    assert jobs.cancel_job("20260101-000000-abcd", root) is False


def test_a_lock_refused_leaves_no_job_directory(jobs_root):
    first = _start(root=jobs_root)
    with pytest.raises(jobs.JobBusy):
        _start(root=jobs_root)
    assert [name for name in os.listdir(jobs_root) if jobs.JOB_ID_RE.match(name)] == [first]


def test_prune_never_removes_the_running_job(jobs_root):
    """Even when it's older than every job kept, and even when its own
    files can't be read: the lock names it."""
    running = _start(root=jobs_root)
    _age(os.path.join(jobs_root, running), jobs.STARTING_GRACE_SECONDS + 5)
    jobs._prune(jobs_root, 0)
    assert jobs.get_job(running, jobs_root)["status"] == "starting"
    os.remove(os.path.join(jobs_root, running, "job.json"))
    jobs._prune(jobs_root, 0)
    assert os.path.isdir(os.path.join(jobs_root, running))


def test_prune_spares_a_directory_still_being_written(jobs_root):
    os.makedirs(os.path.join(jobs_root, "20000101-000000-abcd"))
    jobs._prune(jobs_root, 0)
    assert os.path.isdir(os.path.join(jobs_root, "20000101-000000-abcd"))
    _age(os.path.join(jobs_root, "20000101-000000-abcd"), jobs.STARTING_GRACE_SECONDS + 5)
    jobs._prune(jobs_root, 0)
    assert not os.path.exists(os.path.join(jobs_root, "20000101-000000-abcd"))


@pytest.mark.parametrize("job_id", [
    "20260101-000000-abcd",          # well formed, never existed
    "../../etc", "", None, 42, "20260101-000000-ABCD", "20260101-000000-abcd/../x",
])
def test_cancel_with_ids_that_are_not_jobs(jobs_root, job_id):
    os.makedirs(jobs_root)
    assert jobs.cancel_job(job_id, jobs_root) is False
    assert os.listdir(jobs_root) == []


def test_cancel_of_a_finished_job(jobs_root, redis_server):
    job_id = jobs.start_job("reset", "Quick", [_step("x", "pass")], root=jobs_root)
    _wait_finished(job_id, jobs_root)
    assert jobs.cancel_job(job_id, jobs_root) is False
    assert not os.path.exists(os.path.join(jobs_root, job_id, jobs.CANCEL_NAME))


# --- ADM.11 Jobs outlive the page and the server process -----------------------------

SERVER = r"""
import os, sys
from planetgen.web import jobs
job_id = jobs.start_job("reset", "Outlives", [
    {"label": "Slow", "argv": [sys.executable, "-c", "import time; time.sleep(2); print('still here')"]},
])
print(job_id, flush=True)
"""


@pytest.mark.skipif(jobs.WINDOWS, reason="kills a POSIX process group")
def test_a_job_outlives_the_server_process_that_started_it(jobs_root, redis_server):
    """The web server process (here a stand-in in its own process group)
    starts a job and is then killed with its whole group, as closing the
    browser can't but a server stop might: the runner, detached into its
    own session, carries on and finishes."""
    env = dict(os.environ, PLANETGEN_JOBS_DIR=jobs_root)
    server = subprocess.Popen([PY, "-c", SERVER + "import time; time.sleep(60)"], env=env, cwd=REPO,
                              stdout=subprocess.PIPE, text=True, start_new_session=True)
    try:
        job_id = server.stdout.readline().strip()
        assert jobs.JOB_ID_RE.match(job_id), job_id
    finally:
        os.killpg(server.pid, signal.SIGKILL)
        server.wait(timeout=10)
    job = _wait_finished(job_id, jobs_root)
    assert job["status"] == "succeeded"
    assert "still here" in jobs.log_tail(job_id, root=jobs_root)
