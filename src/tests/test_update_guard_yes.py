# tests/test_update_guard_yes.py

"""
DB.29: `update.sh` and every migration refuse to start while jobs are active unless `--yes` was given or the operator
typed a yes. Nobody to ask is a no: a run started by systemd or cron, one with no terminal, and 30 seconds of silence.
The five paths (no jobs, jobs and `--yes`, jobs and a typed yes, jobs and silence, jobs and a run from systemd or
cron), the live runs of the run record as a source, and the flag in `update.sh`.
"""

import os
import subprocess

import pytest

from planetgen.admin import updateguard
from planetgen.queue import commands, redisqueue, runrecord, units
from planetgen.queue.runrecord import Stage
from planetgen.admin import maintenance
from tests.test_update_guard import FakeWork, REPO, begin, lines, say  # noqa: F401 -- fixtures


def sure_nothing_changed(fake=None):
    assert not maintenance.is_on() and updateguard.read_lock() is None
    if fake is not None:
        assert fake.cancelled == []


# --- the five paths ---------------------------------------------------------------------------------------------

def test_with_no_jobs_nobody_is_asked_and_the_update_goes_on(say):
    asked = []
    assert begin(sources=[FakeWork()], interactive=True, answer=lambda text: asked.append(text) or "n", write=say)
    assert asked == []
    assert begin(sources=[FakeWork()], unattended="cron", write=say)       # a timer with nothing running goes on too


def test_with_jobs_the_flag_cancels_them_and_asks_nothing(say, lines):
    fake = FakeWork("Generate sectors", "Warm-up")
    assert begin(sources=[fake], yes=True, interactive=True, answer=lambda text: pytest.fail(text), write=say)
    assert sorted(fake.cancelled) == ["job-0", "job-1"]


def test_with_jobs_a_typed_yes_cancels_them_and_a_typed_no_changes_nothing(say):
    fake = FakeWork("Generate sectors")
    assert not begin(sources=[fake], interactive=True, answer=lambda text: "n\n", write=say)
    assert not begin(sources=[fake], interactive=True, answer=lambda text: "\n", write=say)
    sure_nothing_changed(fake)
    assert begin(sources=[fake], interactive=True, answer=lambda text: "y\n", write=say)
    assert fake.cancelled == ["job-0"]


def test_with_jobs_thirty_seconds_of_silence_is_a_no(say, lines):
    fake = FakeWork("Generate sectors")
    assert not begin(sources=[fake], interactive=True, answer=lambda text: None, write=say)
    text = "\n".join(lines)
    assert f"No answer in {updateguard.PROMPT_SECONDS} seconds: taking that as no." in text
    assert "nothing was changed" in text
    sure_nothing_changed(fake)
    assert updateguard.PROMPT_SECONDS == 30


@pytest.mark.parametrize("starter", ["systemd", "cron"])
def test_with_jobs_a_run_from_systemd_or_cron_asks_nothing_and_is_a_no_unless_the_flag_is_given(say, lines, starter):
    fake = FakeWork("Generate sectors")
    asked = []
    assert not begin(sources=[fake], interactive=True, unattended=starter,
                     answer=lambda text: asked.append(text) or "y\n", write=say)       # even a terminal does not count
    assert asked == []
    text = "\n".join(lines)
    assert f"started by {starter}" in text and "add --yes" in text
    sure_nothing_changed(fake)
    assert begin(sources=[fake], unattended=starter, yes=True, write=say)
    assert fake.cancelled == ["job-0"]


# --- who started the run ----------------------------------------------------------------------------------------

def test_a_service_with_its_output_in_the_journal_is_systemd_and_cron_is_an_ancestor_by_name():
    unattended = updateguard.started_unattended
    assert unattended({"JOURNAL_STREAM": "8:12345"}, parents=[], stdin_is_terminal=True) == "systemd"
    assert unattended({"INVOCATION_ID": "abc"}, parents=[], stdin_is_terminal=False) == "systemd"
    assert unattended({}, parents=["sh", "sudo", "cron", "systemd"], stdin_is_terminal=False) == "cron"
    assert unattended({}, parents=["crond"], stdin_is_terminal=False) == "cron"


def test_a_person_at_a_terminal_is_not_taken_for_a_timer():
    unattended = updateguard.started_unattended
    assert unattended({"INVOCATION_ID": "abc"}, parents=["bash", "sshd"], stdin_is_terminal=True) is None
    assert unattended({}, parents=["bash", "sudo", "sshd", "systemd"], stdin_is_terminal=True) is None


def test_the_question_waits_for_a_line_and_gives_up_after_the_seconds(capsys):
    read, write = os.pipe()
    with os.fdopen(read) as stdin, os.fdopen(write, "w") as feed:
        assert updateguard.answer_within("Cancel? ", seconds=0.1, stdin=stdin) is None
        feed.write("y\n")
        feed.flush()
        assert updateguard.answer_within("Cancel? ", seconds=5, stdin=stdin) == "y\n"
    assert capsys.readouterr().out.count("Cancel? ") == 2


# --- the runs in the run record ---------------------------------------------------------------------------------

KIND = "guard-yes-test"


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture(autouse=True)
def _kind(monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, KIND, commands.CommandSpec(KIND, "A test", unit=(__name__, "_noop"),
                                                                     expand=(__name__, "_cut")))


def _noop(args, report):
    return {}


def _cut(params):
    return commands.Expansion([({"n": n}, 1.0) for n in range(params["count"])], {})


def _live_run(conn, count=20):
    run_id = runrecord.start_run(conn, "Guard test", [Stage("s0", "Stage 0", units=count)], who="boss", now=100.0)
    runrecord.begin_stage(conn, run_id, 0, now=101.0)
    units.enqueue_units(conn, run_id, 0, KIND, [({"n": n}, 1.0) for n in range(count)])
    return run_id


def test_a_live_run_is_listed_once_without_its_units_and_a_finished_one_is_not(conn, redis_server):
    run_id = _live_run(conn)
    finished = runrecord.start_run(conn, "Done", [Stage("s0", "Stage 0", units=1)], kind=runrecord.CHECK, now=100.0)
    runrecord.finish_run(conn, finished, "done", now=110.0)
    found = updateguard.RunRecordSource(redis_server).items()
    assert [(item["id"], item["title"], item["kind"], item["queued"]) for item in found] == \
        [(run_id, "Guard test", runrecord.GENERATION, False)]
    assert found[0]["started"] and found[0]["holds"] == "the generation slot"
    assert updateguard.QueuedJobSource(redis_server).items() == []        # its 20 units are the run's, not jobs
    assert "Guard test (generation " in "\n".join(updateguard.describe(found))
    assert any(isinstance(source, updateguard.RunRecordSource) for source in updateguard.default_sources())


def test_a_live_run_stops_the_update_without_the_word_and_is_cancelled_with_it(conn, redis_server, say, lines):
    run_id = _live_run(conn)
    sources = [updateguard.RunRecordSource(redis_server), updateguard.QueuedJobSource(redis_server)]
    assert not begin(sources=sources, write=say)
    assert "Guard test" in "\n".join(lines) and "--yes" in "\n".join(lines)
    assert runrecord.read_run(conn, run_id)["state"] == "running"
    sure_nothing_changed()
    assert begin(sources=sources, yes=True, write=say, timeout=30, poll=0.1)
    run = runrecord.read_run(conn, run_id)
    assert run["state"] == "cancelled" and run["error"] == "cancelled by update"
    assert runrecord.active_run(conn) is None
    assert redisqueue.queue(redisqueue.BULK_QUEUE, conn).count == 0


def test_without_a_redis_server_the_source_lists_nothing():
    assert updateguard.RunRecordSource("redis://127.0.0.1:1/0").items() == []


# --- update.sh --------------------------------------------------------------------------------------------------

def test_update_sh_takes_yes_and_no_longer_the_old_flag():
    done = subprocess.run(["bash", os.path.join(REPO, "update.sh"), "--cancel-running"], capture_output=True, text=True,
                          timeout=60)
    assert done.returncode == 2
    assert "unknown argument '--cancel-running'" in done.stdout + done.stderr
    assert "the options are --yes, --reset-jobs and --kill-all" in done.stdout + done.stderr
    with open(os.path.join(REPO, "update.sh"), encoding="utf-8") as handle:
        text = handle.read()
    assert "--yes) ANSWER_YES=(--yes) ;;" in text and "cancel-running" not in text
