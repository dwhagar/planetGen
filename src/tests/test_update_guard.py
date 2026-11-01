# tests/test_update_guard.py

"""
The update guard (OPS.41, OPS.53): `admin/updateguard.py`, its command `cli/update_guard.py`, the manual
`cli/migrate.py` run, and the wiring in `update.sh`. The running work is fake (`FakeWork`) except where a test
says it is a real process or a real Generate job directory.
"""

import os
import re
import subprocess
import sys

import pytest

from planetgen.admin import maintenance, updateguard
from planetgen.cli import migrate


REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))


class FakeWork:
    """A source of running work: `items` are the jobs it holds; `stops` says whether cancelling ends them."""

    name = "fake"

    def __init__(self, *titles, stops=True, queued=()):
        self.live = [{"source": self.name, "kind": "generate", "id": f"job-{index}", "title": title,
                      "started": "2026-10-16T03:00:00Z", "holds": "the work lease", "queued": title in queued}
                     for index, title in enumerate(titles)]
        self.stops = stops
        self.cancelled = []

    def items(self):
        return [dict(item) for item in self.live]

    def cancel(self, item):
        self.cancelled.append(item["id"])
        if self.stops:
            self.live = [other for other in self.live if other["id"] != item["id"]]
        return True


def begin(pid=None, **options):
    options.setdefault("write", lambda text: None)
    options.setdefault("poll", 0.0)
    options.setdefault("timeout", 0.2)
    return updateguard.begin(pid or os.getpid(), **options)


@pytest.fixture
def lines():
    return []


@pytest.fixture
def say(lines):
    return lines.append


# --- OPS.41: the site is closed while the update runs ------------------------------------------------------

def test_an_update_closes_the_site_and_opens_it_when_it_has_finished(say, lines):
    assert begin(sources=[], write=say)
    state = maintenance.read()
    assert state["enabled"] and state["cause"] == updateguard.UPDATE_CAUSE and state["by"] == "update.sh"
    assert state["reason"] == updateguard.UPDATE_REASON
    assert updateguard.read_lock()["pid"] == os.getpid() and updateguard.is_updating()
    assert updateguard.end(os.getpid(), ok=True, write=say)
    assert not maintenance.is_on() and updateguard.read_lock() is None and not updateguard.is_updating()
    assert lines == []


def test_a_failed_update_leaves_the_site_closed_and_says_how_to_open_it(say, lines):
    assert begin(sources=[])
    assert not updateguard.end(os.getpid(), ok=False, step="step 5/9, migrating the database", write=say)
    assert maintenance.is_on() and maintenance.read()["cause"] == updateguard.UPDATE_CAUSE
    message = " ".join(lines)
    assert "step 5/9, migrating the database" in message and "half migrated" in message
    assert "python -m planetgen.cli.maintenance off" in message
    # the lock is released so the update can be run again, but nothing else may start while the site is closed
    assert updateguard.read_lock() is None and updateguard.is_updating()
    with pytest.raises(updateguard.UpdateRunning):
        updateguard.refuse_while_updating("A Generate job")
    # the admin opens the site: work may start, and the next update closes and opens it as usual
    maintenance.turn_off(by="boss")
    updateguard.refuse_while_updating("A Generate job")
    assert begin(sources=[]) and updateguard.end(os.getpid(), ok=True)
    assert not maintenance.is_on()


def test_a_site_an_admin_closed_stays_closed_after_the_update():
    maintenance.turn_on("until noon", by="boss")
    assert begin(sources=[])
    assert maintenance.read()["by"] == "boss" and maintenance.read()["eta"] == "until noon"
    updateguard.end(os.getpid(), ok=True)
    assert maintenance.is_on() and maintenance.read()["by"] == "boss"


def test_the_lock_of_a_dead_process_is_ignored():
    child = subprocess.Popen([sys.executable, "-c", "pass"])
    child.wait()
    updateguard._write_lock(child.pid, "update.sh", kept=False)
    assert updateguard.read_lock() is None and not updateguard.is_updating()
    assert begin(sources=[]) and updateguard.read_lock()["pid"] == os.getpid()


def test_a_second_update_is_refused_while_one_runs(say, lines):
    sleeper = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    try:
        updateguard._write_lock(sleeper.pid, "update.sh", kept=False)
        assert not begin(sources=[], write=say)
        assert f"process {sleeper.pid}" in " ".join(lines)
        assert not maintenance.is_on()
    finally:
        sleeper.kill()
        sleeper.wait()


def test_the_updates_own_commands_pass_the_guard(monkeypatch):
    assert begin(sources=[])
    with pytest.raises(updateguard.UpdateRunning, match="the site is being updated"):
        updateguard.refuse_while_updating("A Generate job")
    monkeypatch.setenv(updateguard.OWNER_ENV, str(os.getpid()))
    updateguard.refuse_while_updating("A Generate job")
    monkeypatch.setenv(updateguard.OWNER_ENV, "1")
    with pytest.raises(updateguard.UpdateRunning):
        updateguard.refuse_while_updating("A Generate job")


def test_the_sector_switch_takes_over_the_closing_and_the_update_leaves_it_to_the_warm_up(monkeypatch):
    monkeypatch.setattr(migrate, "schema_status", lambda config: (80, 12))
    assert begin(sources=[])
    assert migrate.close_site_for_the_sector_switch(object(), write=lambda text: None)
    assert maintenance.read()["cause"] == maintenance.SECTOR_SWITCH
    updateguard.end(os.getpid(), ok=True)
    assert maintenance.is_on() and maintenance.read()["cause"] == maintenance.SECTOR_SWITCH
    assert maintenance.turn_off_for(maintenance.SECTOR_SWITCH)


def test_the_sector_switch_leaves_an_admins_own_closing_alone(monkeypatch):
    monkeypatch.setattr(migrate, "schema_status", lambda config: (80, 12))
    maintenance.turn_on("noon", by="boss")
    assert not migrate.close_site_for_the_sector_switch(object(), write=lambda text: None)
    assert maintenance.read()["by"] == "boss"


# --- the command -------------------------------------------------------------------------------------------

def run_cli(*args, env=None):
    return subprocess.run([sys.executable, "-m", "planetgen.cli.update_guard", *args], capture_output=True,
                          text=True, env={**os.environ, **(env or {})}, stdin=subprocess.DEVNULL, timeout=120)


def test_the_command_begins_ends_and_reports(tmp_path):
    shell = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    try:
        status = run_cli("status")
        assert status.returncode == 0 and "No update is running." in status.stdout
        done = run_cli("begin", "--pid", str(shell.pid))
        assert done.returncode == 0, done.stderr
        assert "An update is running" in run_cli("status").stdout
        assert "Maintenance mode is ON" in run_cli("status").stdout
        failed = run_cli("end", "--pid", str(shell.pid), "--failed", "--step", "step 3/9, checking things")
        assert failed.returncode == 0 and "step 3/9, checking things" in failed.stdout
        assert "left the site closed" in run_cli("status").stdout
        assert run_cli("begin", "--pid", str(shell.pid)).returncode == 0
        assert run_cli("end", "--pid", str(shell.pid)).returncode == 0
        assert "No update is running." in run_cli("status").stdout
    finally:
        shell.kill()
        shell.wait()


def test_the_command_needs_the_process_to_hold_the_lock_for():
    result = run_cli("begin")
    assert result.returncode == 2 and "needs --pid" in result.stderr


# --- a migration run by hand --------------------------------------------------------------------------------

@pytest.fixture
def fake_machine(monkeypatch):
    """A machine with `FakeWork` running, in place of the real sources."""
    work = FakeWork()
    monkeypatch.setattr(updateguard, "default_sources", lambda config=None, url=None: [work])
    monkeypatch.setattr(updateguard, "STOP_SECONDS", 0.2)
    monkeypatch.setattr(updateguard, "POLL_SECONDS", 0.0)
    monkeypatch.setattr(sys.stdin, "isatty", lambda: False)
    return work


def test_a_manual_migration_with_something_to_migrate_closes_the_site_and_opens_it(monkeypatch, fake_machine):
    monkeypatch.setattr(migrate, "schema_status", lambda config: (80, 3))
    seen = []
    with migrate._update_guard(object(), False):
        seen.append(maintenance.read()["cause"])
    assert seen == [updateguard.UPDATE_CAUSE] and not maintenance.is_on() and updateguard.read_lock() is None


def test_a_manual_migration_that_fails_leaves_the_site_closed(monkeypatch, fake_machine):
    monkeypatch.setattr(migrate, "schema_status", lambda config: (80, 3))
    with pytest.raises(SystemExit):
        with migrate._update_guard(object(), False):
            raise SystemExit(1)
    assert maintenance.is_on() and maintenance.read()["cause"] == updateguard.UPDATE_CAUSE
    assert updateguard.read_lock() is None


def test_a_migration_with_nothing_to_migrate_closes_nothing(monkeypatch, fake_machine):
    monkeypatch.setattr(migrate, "schema_status", lambda config: (91, 0))
    with migrate._update_guard(object(), False):
        assert not maintenance.is_on()


def test_a_migration_that_update_sh_runs_leaves_the_guard_to_it(monkeypatch, fake_machine):
    monkeypatch.setattr(migrate, "schema_status", lambda config: (80, 3))
    assert begin(sources=[])
    monkeypatch.setenv(updateguard.OWNER_ENV, str(os.getpid()))
    with migrate._update_guard(object(), False):
        pass
    assert maintenance.is_on() and updateguard.read_lock() is not None      # still update.sh's


def test_a_manual_migration_refuses_running_work_without_a_terminal_and_changes_nothing(monkeypatch, fake_machine):
    monkeypatch.setattr(migrate, "schema_status", lambda config: (80, 3))
    fake_machine.live = FakeWork("Generate sectors").live
    with pytest.raises(SystemExit) as stopped:
        with migrate._update_guard(object(), False):
            raise AssertionError("the migration must not start")
    assert stopped.value.code == 1 and fake_machine.cancelled == []
    assert not maintenance.is_on() and updateguard.read_lock() is None


# --- update.sh -----------------------------------------------------------------------------------------------

def update_sh():
    with open(os.path.join(REPO, "update.sh"), encoding="utf-8") as handle:
        return handle.read()


def test_update_sh_takes_the_guard_before_its_first_change_and_requires_it_by_step_two():
    text = update_sh()
    first = text.index("\nensure_update_guard\n")
    assert text.index("trap finish_update_guard EXIT") < first < text.index("git fetch origin")
    assert text.index("export PLANETGEN_UPDATE_PID=$$") < first
    # taken again after the pull (a checkout that did not have it) and required once the venv exists
    assert text.index("git reset --hard") < text.index("\nensure_update_guard\n", first + 1) < text.index("log_stage 2 9")
    assert text.index("use_site_python\nensure_update_guard required") < text.index("log_stage 3 9")
    assert text.index("ensure_update_guard required") < text.index("\nmigrate_or_reset_db\n")
    assert "--yes" in text


def test_update_sh_names_every_step_for_the_failure_message():
    text = update_sh()
    for step in range(1, 10):
        assert re.search(rf'^UPDATE_STEP="step {step}/9, [^"]+"\n(?:.*\n)?log_stage {step} 9', text, re.M), step


def test_update_sh_never_edits_the_redis_config_or_the_apache_site():
    text = update_sh() + open(os.path.join(REPO, "scripts", "deploy-common.sh"), encoding="utf-8").read()
    guard = text[text.index("UPDATE_GUARD_HELD=0"):text.index("run_warm_stage() {")]
    assert "redis.conf" not in guard and "sites-available" not in guard and "sed -i" not in guard


def run_shell(script, *, status_out=None):
    """Runs `script` in bash with the guard functions loaded and `update_guard` replaced by a recorder."""
    prelude = f'''
        set -euo pipefail
        SCRIPT_DIR={REPO!r}
        export PLANETGEN_VENV=$(mktemp -d); VENV=$PLANETGEN_VENV
        mkdir -p "$VENV/bin" "$VENV/src"
        printf '#!/bin/sh\\nexit 0\\n' > "$VENV/bin/python"; chmod +x "$VENV/bin/python"
        ANSWER_YES=()
        source {REPO}/scripts/deploy-common.sh
        update_guard() {{ echo "update_guard $*" >> "$LOG"; return "${{GUARD_STATUS:-0}}"; }}
    '''
    return subprocess.run(["bash", "-c", prelude + script], capture_output=True, text=True,
                          env={**os.environ, "LOG": status_out}, timeout=60)


def test_the_shell_functions_hold_the_guard_once_and_end_it_with_the_step(tmp_path):
    log = tmp_path / "calls"
    done = run_shell('''
        ANSWER_YES=(--yes)
        ensure_update_guard; ensure_update_guard
        UPDATE_STEP="step 5/9, migrating the database"
        (exit 3) || finish_update_guard
    ''', status_out=str(log))
    assert done.returncode == 0, done.stderr
    calls = log.read_text().splitlines()
    assert calls[0].startswith("update_guard begin --pid ") and calls[0].endswith("--yes")
    assert len(calls) == 2
    assert calls[1].startswith("update_guard end --pid ") and "--failed --step step 5/9, migrating the database" in calls[1]


def test_a_refused_guard_stops_the_update_and_ends_nothing(tmp_path):
    log = tmp_path / "calls"
    done = run_shell('''
        trap finish_update_guard EXIT
        GUARD_STATUS=1 ensure_update_guard
        echo "unreachable"
    ''', status_out=str(log))
    assert done.returncode == 1 and "unreachable" not in done.stdout
    assert "the update did not start" in done.stderr
    assert log.read_text().splitlines()[0].startswith("update_guard begin")
    assert len(log.read_text().splitlines()) == 1       # no end: the guard was never held


def test_a_good_update_ends_the_guard_without_the_failed_flag(tmp_path):
    log = tmp_path / "calls"
    done = run_shell('''
        ensure_update_guard
        finish_update_guard
    ''', status_out=str(log))
    assert done.returncode == 0, done.stderr
    assert "--failed" not in log.read_text().splitlines()[1]
