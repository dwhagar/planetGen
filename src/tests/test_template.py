"""OPS.226: the warm template and the jobs forked from it (`queue/template.py`, `queue/jobhost.py`): real processes, no database. Each job's exit is known at once and exactly, a closed pipe
makes a job finish its piece and exit, and a job outlives neither a template's death nor a supervisor's without being noticed."""

import os
import signal
import time

import pytest

from planetgen.queue import allowlist, jobhost, template

R = "tests.template_runners"


@pytest.fixture
def host(monkeypatch):
    monkeypatch.setenv(allowlist.EXTRA_ENV, f"{R}:*")
    made = jobhost.JobHost(warm=[f"{R}.warm"])
    made.start(ready_seconds=60)
    yield made
    made.stop()


def _exits(host, count, seconds=20.0):
    found, deadline = [], time.monotonic() + seconds
    while len(found) < count and time.monotonic() < deadline:
        found += host.poll(0.2)
    return found


def test_each_job_ends_with_the_code_or_the_signal_the_template_read(host):
    ok, three, killed, boom = (host.fork(f"{R}.quick", {"code": 0}), host.fork(f"{R}.quick", {"code": 3}), host.fork(f"{R}.kills_itself"), host.fork(f"{R}.explodes"))
    by_id = {exit_.id: exit_ for exit_ in _exits(host, 4)}
    assert [(by_id[job.id].code, by_id[job.id].signal) for job in (ok, three, killed, boom)] == [(0, None), (3, None), (None, signal.SIGKILL), (template.EXIT_FAILED, None)]
    assert all(by_id[job.id].source == "template" for job in (ok, three, killed, boom)) and not host.jobs


def test_a_job_that_raised_says_so_on_its_pipe(host):
    job = host.fork(f"{R}.explodes")
    _exits(host, 1)
    assert any("ValueError: the piece broke" in message.get("error", "") for message in job.messages())
    job.close()


def test_every_job_is_a_fork_of_the_one_warm_template(host):
    job = host.fork(f"{R}.uses_the_warm_set")
    assert _exits(host, 1)[0].code == 0 and job.pid != host.template_pid


def test_the_supervisor_and_a_job_talk_over_the_pipe(host):
    job = host.fork(f"{R}.chats")
    heard, deadline = [], time.monotonic() + 10
    while not heard and time.monotonic() < deadline:
        heard = job.messages()
        time.sleep(0.02)
    assert heard == [{"hello": job.pid}]
    job.send({"op": "stop"})
    assert _exits(host, 1)[0].code == 0


def test_a_job_whose_supervisor_is_gone_finishes_its_piece_and_exits(host, tmp_path):
    marks = tmp_path / "marks.txt"
    job = host.fork(f"{R}.pieces", {"file": str(marks), "count": 50, "piece_seconds": 0.1})
    deadline = time.monotonic() + 10
    while not marks.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    job.close()                    # the supervisor lets go of the pipe
    exit_ = _exits(host, 1)[0]
    lines = marks.read_text().splitlines()
    assert exit_.code == template.EXIT_SUPERVISOR_GONE
    assert lines[-1].startswith("done") and len(lines) < 100, "it finished the piece it was in and then stopped, long before its fifty"
    assert len([line for line in lines if line.startswith("start")]) == len([line for line in lines if line.startswith("done")])


def test_a_job_that_outlives_its_template_is_still_noticed_by_the_kernel_handle(host):
    job = host.fork(f"{R}.lingers")
    os.kill(host.template_pid, signal.SIGKILL)
    os.kill(job.pid, signal.SIGKILL)
    if not hasattr(os, "pidfd_open"):
        pytest.skip("no pidfd on this system")
    exit_ = _exits(host, 1)[0]
    assert (exit_.id, exit_.code, exit_.signal, exit_.source) == (job.id, None, None, "pidfd")


def test_a_job_asked_of_a_template_that_has_died_is_not_ready_and_is_not_failed(host):
    """OPS.236: the launcher marked the row `failed` for the BrokenPipeError of a fork on a template that was just killed; it is `NotReady` now, the host is dead, and a `TemplateSet` starts a new template."""
    os.kill(host.template_pid, signal.SIGKILL)
    host.process.wait()
    with pytest.raises(jobhost.NotReady, match="the template ended"):
        host.fork(f"{R}.quick", {"code": 0})
    assert host._dead and not host.jobs


def test_a_template_with_a_second_thread_does_not_start(monkeypatch):
    monkeypatch.setenv(allowlist.EXTRA_ENV, f"{R}:*")
    with pytest.raises(jobhost.HostError):
        jobhost.JobHost(warm=[f"{R}.warm_with_a_thread"]).start(ready_seconds=30)


def test_a_runner_that_is_not_on_the_allow_list_is_refused_and_nothing_runs(host, monkeypatch):
    with pytest.raises(jobhost.HostError, match="not a function this program queues"):
        host.fork("os.system", {"command": "true"})
    assert not host.jobs


def test_the_star_tables_are_a_warm_set_a_template_may_build():
    """The warm set the bulk lane's template builds (OPS.228) is on the allow-list and is the function workers build the tables with today (building them takes about a minute, so this does not run it)."""
    from planetgen.queue import units
    allowlist.check_function("planetgen.queue.units.warm_star_tables", "a warm set")
    assert template._function("planetgen.queue.units.warm_star_tables") is units.warm_star_tables


def test_without_kernel_handles_the_template_reports_every_exit_and_a_dead_template_leaves_the_pids_to_look_for(monkeypatch):
    """The path macOS takes (OPS.226): no `pidfd`, so the template's report is the way, and when the template is the one that died each job's pid is looked for."""
    monkeypatch.setenv(allowlist.EXTRA_ENV, f"{R}:*")
    made = jobhost.JobHost(pidfd=False)
    made.start(ready_seconds=60)
    try:
        quick, lingering = made.fork(f"{R}.quick", {"code": 4}), made.fork(f"{R}.lingers")
        assert quick.pidfd is None and lingering.pidfd is None
        assert [(exit_.id, exit_.code, exit_.source) for exit_ in _exits(made, 1)] == [(quick.id, 4, "template")]
        os.kill(made.template_pid, signal.SIGKILL)
        time.sleep(0.2)
        os.kill(lingering.pid, signal.SIGKILL)
        exit_ = _exits(made, 1)[0]
        assert (exit_.id, exit_.code, exit_.signal, exit_.source) == (lingering.id, None, None, "pid")
    finally:
        made.stop()


def test_each_kind_is_forked_from_the_template_of_its_lane_and_a_template_starts_when_first_asked_for(monkeypatch):
    monkeypatch.setenv(allowlist.EXTRA_ENV, f"{R}:*")
    assert jobhost.template_key("scatter-mass") == "tables" and jobhost.template_key("generate-system") == "plain"
    templates = jobhost.TemplateSet(warm_sets={"tables": [f"{R}.warm"], "plain": []}, ready_seconds=60)
    try:
        assert templates.hosts == {}, "nothing is started before a job asks"
        key, plain = templates.fork("generate-system", f"{R}.quick", {"code": 0})
        assert key == "plain" and list(templates.hosts) == ["plain"], "a host with no star-drawing work never builds the tables"
        key, tables = templates.fork("scatter-mass", f"{R}.uses_the_warm_set")
        assert key == "tables" and set(templates.hosts) == {"plain", "tables"}
        found, deadline = {}, time.monotonic() + 20
        while len(found) < 2 and time.monotonic() < deadline:
            for lane, exit_ in templates.poll(0.2):
                found[(lane, exit_.id)] = exit_.code
        assert sorted(found.values()) == [0, 0] and {lane for lane, _ in found} == {"plain", "tables"}
    finally:
        templates.stop()


def test_a_template_set_that_builds_in_the_background_answers_not_ready_until_the_warm_set_is_built(monkeypatch):
    monkeypatch.setenv(allowlist.EXTRA_ENV, f"{R}:*")
    templates = jobhost.TemplateSet(warm_sets={"plain": [f"{R}.warm"]}, ready_seconds=60, background=True)
    try:
        started = time.monotonic()
        with pytest.raises(jobhost.NotReady):
            templates.fork("generate-system", f"{R}.quick", {"code": 0})
        assert time.monotonic() - started < 5, "the ask returns at once, the build goes on in the background"
        deadline, forked = time.monotonic() + 60, None
        while forked is None and time.monotonic() < deadline:
            try:
                forked = templates.fork("generate-system", f"{R}.quick", {"code": 0})
            except jobhost.NotReady:
                time.sleep(0.2)
        assert forked is not None and forked[0] == "plain"
    finally:
        templates.stop()


# --- OPS.283: arguments bigger than the socket's limit -------------------------------------------------------------------------------------------

BIG = 300 * 1024                                            # over the 213 KB a control-socket message can carry (the socket's send buffer)


@pytest.fixture
def state_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))          # the template process inherits it, so both ends name the same folder
    return tmp_path


def _spilled(folder):
    return sorted(os.listdir(folder)) if os.path.isdir(folder) else []


def test_a_job_whose_arguments_are_over_the_socket_limit_gets_them_by_a_file_that_is_gone_when_it_ends(state_dir, host):
    blob = "x" * BIG
    from socket import SOL_SOCKET, SO_SNDBUF
    assert len(template.encode({"op": "fork", "id": 1, "runner": "r", "job": {"blob": blob}})) > host.control.getsockopt(SOL_SOCKET, SO_SNDBUF), "this test must give more than the socket takes in one message"
    job = host.fork(f"{R}.measures_its_arguments", {"blob": blob, "size": BIG})
    assert _exits(host, 1)[0].code == 0, "the job saw all of its arguments"
    assert _spilled(template.spill_folder()) == [], "the file the arguments came in is deleted"
    job.close()


def test_the_file_is_deleted_when_the_job_could_not_be_handed_out(state_dir, host):
    class Refusing:
        def __init__(self, real):
            self.real = real

        def send(self, data):
            raise BrokenPipeError

        def __getattr__(self, name):
            return getattr(self.real, name)

    host.control = Refusing(host.control)
    with pytest.raises(jobhost.NotReady):
        host.fork(f"{R}.measures_its_arguments", {"blob": "x" * BIG, "size": BIG})
    assert _spilled(template.spill_folder()) == []


def test_a_job_with_small_arguments_gets_no_file(state_dir, host):
    host.fork(f"{R}.measures_its_arguments", {"blob": "x" * 100, "size": 100})
    assert _exits(host, 1)[0].code == 0 and _spilled(template.spill_folder()) == []


def test_the_template_reads_only_files_in_the_arguments_folder(state_dir, tmp_path_factory):
    outside = tmp_path_factory.mktemp("elsewhere") / "job.json"
    outside.write_bytes(template.encode({"blob": "x"}))
    with pytest.raises(ValueError):
        template.unspill(str(outside))
    assert outside.exists(), "a file that is not ours is left where it was"
    path = template.spill({"blob": "y"}, 7)
    assert template.unspill(path) == {"blob": "y"} and not os.path.exists(path)


def test_argument_files_an_earlier_supervisor_left_are_swept_after_a_day(state_dir):
    old, new = template.spill({"a": 1}, 1), template.spill({"a": 2}, 2)
    aged = time.time() - template.SPILL_KEPT_SECONDS - 60
    os.utime(old, (aged, aged))
    assert template.sweep() == 1 and not os.path.exists(old) and os.path.exists(new)
