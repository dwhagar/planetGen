# tests/test_lane_services.py

"""
The check before a test run and the watchdog during it (TEST.376, `tests/lane_services.py`, `scripts/test.py`). The
servers are fakes here: no test stops or starts the MySQL or Redis the other tests are using.
"""

import json
import threading

import pytest

from tests import lane_services as services
from tests.test_test_runner import runner

ENVIRON = {"PLANETGEN_TEST_MYSQL_USER": "pgtest", "PLANETGEN_TEST_MYSQL_PASSWORD": "x",
           "PLANETGEN_TEST_REDIS_URL": "redis://127.0.0.1:6379/0"}


@pytest.fixture(autouse=True)
def _runner_in_a_throwaway_footing(tmp_path, monkeypatch):
    """`runner.main` writes a progress log, clears the test Redis and, here, would look at the real servers: it does not."""
    monkeypatch.setattr(runner, "PROGRESS_DIR", str(tmp_path / "progress"))
    monkeypatch.setattr(runner, "clear_test_redis", lambda: "")
    monkeypatch.setattr(runner, "_services", lambda: _Quiet)


class _Quiet:
    """Services that are in order and never trip."""

    class Watchdog:
        def __init__(self, environ, on_trip, **options):
            self.events, self.what, self.failed, self.tripped = [], "", False, threading.Event()

        def start(self):
            pass

        def stop(self):
            pass

    @staticmethod
    def preflight(environ, browser=True, say=print):
        return []


@pytest.fixture
def ordinary(monkeypatch):
    """Everything in order: the servers answer, the disk has room, the packages are there."""
    monkeypatch.setattr(services, "mysql_state", lambda environ: "up")
    monkeypatch.setattr(services, "redis_state", lambda environ: "up")
    monkeypatch.setattr(services, "free_mb", lambda path="/": 50_000)
    monkeypatch.setattr(services, "package_problems", lambda browser: [])
    monkeypatch.setattr(services, "load_profile", lambda environ, path=services.PROFILE: [])


# --- the variables ---------------------------------------------------------------------

def test_the_profile_supplies_the_variables_that_are_missing_and_only_those(tmp_path):
    profile = tmp_path / "planetgen-test.sh"
    profile.write_text("# written by setup\nexport PLANETGEN_TEST_MYSQL_HOST=127.0.0.1\n"
                       "export PLANETGEN_TEST_MYSQL_PASSWORD='s3 cret'\nexport PLANETGEN_TEST_MYSQL_USER=root\n"
                       "export PIP_BREAK_SYSTEM_PACKAGES=1\nnot an export\n")
    environ = {"PLANETGEN_TEST_MYSQL_USER": "pgtest"}
    loaded = services.load_profile(environ, str(profile))
    assert environ["PLANETGEN_TEST_MYSQL_PASSWORD"] == "s3 cret" and environ["PLANETGEN_TEST_MYSQL_HOST"] == "127.0.0.1"
    assert environ["PLANETGEN_TEST_MYSQL_USER"] == "pgtest"                  # what is set stays
    assert sorted(loaded) == ["PIP_BREAK_SYSTEM_PACKAGES", "PLANETGEN_TEST_MYSQL_HOST", "PLANETGEN_TEST_MYSQL_PASSWORD"]
    assert services.load_profile({}, str(tmp_path / "missing")) == []


def test_without_the_test_variables_the_run_does_not_start(monkeypatch):
    monkeypatch.setattr(services, "load_profile", lambda environ, path=services.PROFILE: [])
    problems = services.preflight({}, say=lambda text: None)
    assert len(problems) == 1 and "PLANETGEN_TEST_MYSQL_USER" in problems[0] and "wrong account" in problems[0]


# --- the check ---------------------------------------------------------------------------

def test_servers_that_answer_leave_nothing_to_say(ordinary):
    said = []
    assert services.preflight(dict(ENVIRON), say=said.append) == [] and said == []


def test_a_dead_mysql_and_a_dead_redis_are_started(ordinary, monkeypatch):
    states = {"mysql": ["down"], "redis": ["down"]}
    started = []
    monkeypatch.setattr(services, "mysql_state", lambda environ: states["mysql"][0])
    monkeypatch.setattr(services, "redis_state", lambda environ: states["redis"][0])
    monkeypatch.setattr(services, "start_mysql", lambda environ: started.append("mysql") or states["mysql"].__setitem__(0, "up") or True)
    monkeypatch.setattr(services, "start_redis", lambda environ: started.append("redis") or states["redis"].__setitem__(0, "up") or True)
    said = []
    assert services.preflight(dict(ENVIRON), say=said.append) == []
    assert started == ["mysql", "redis"]
    assert any("MySQL is up again" in line for line in said) and any("Redis is up again" in line for line in said)


def test_a_server_that_cannot_be_started_is_a_problem(ordinary, monkeypatch):
    monkeypatch.setattr(services, "mysql_state", lambda environ: "down")
    monkeypatch.setattr(services, "start_mysql", lambda environ: False)
    monkeypatch.setattr(services, "redis_state", lambda environ: "down")
    monkeypatch.setattr(services, "start_redis", lambda environ: False)
    problems = services.preflight(dict(ENVIRON), say=lambda text: None)
    assert len(problems) == 2 and "MySQL does not answer" in problems[0] and "Redis does not answer" in problems[1]


def test_a_refused_login_is_reported_and_the_server_is_left_alone(ordinary, monkeypatch):
    monkeypatch.setattr(services, "mysql_state", lambda environ: "denied")
    monkeypatch.setattr(services, "start_mysql", lambda environ: pytest.fail("restarted a server that is up"))
    problems = services.preflight(dict(ENVIRON), say=lambda text: None)
    assert len(problems) == 1 and "refused" in problems[0]


def test_a_nearly_full_disk_has_the_binlogs_purged(ordinary, monkeypatch):
    free = [800]
    purged = []
    monkeypatch.setattr(services, "free_mb", lambda path="/": free[0])
    monkeypatch.setattr(services, "purge_binlogs", lambda environ: purged.append(1) or free.__setitem__(0, 20_000) or True)
    assert services.preflight(dict(ENVIRON), say=lambda text: None) == [] and purged == [1]
    free[0] = 800
    monkeypatch.setattr(services, "purge_binlogs", lambda environ: True)
    problems = services.preflight(dict(ENVIRON), say=lambda text: None)
    assert len(problems) == 1 and "free space" in problems[0]


def test_packages_are_checked_by_version(monkeypatch):
    versions = {"argon2-cffi": "25.1.0", "playwright": services.PLAYWRIGHT_VERSION}

    def version(name):
        if name not in versions:
            raise services.importlib.metadata.PackageNotFoundError(name)
        return versions[name]

    monkeypatch.setattr(services.importlib.metadata, "version", version)
    assert services.package_problems(browser=True) == []
    versions["playwright"] = "1.57.0"
    assert services.package_problems(browser=False) == []                    # not needed without browser tests
    problem, = services.package_problems(browser=True)
    assert "1.57.0" in problem and services.PLAYWRIGHT_VERSION in problem and "playwright install" in problem
    del versions["argon2-cffi"]
    assert any("argon2-cffi is not installed" in text for text in services.package_problems(browser=False))


# --- the watchdog ------------------------------------------------------------------------

def _watchdog(mysql=("up",), redis=("up",), restart=True, **options):
    answers = {"mysql": list(mysql), "redis": list(redis)}
    restarted = []

    def ask(name):
        return lambda environ: answers[name].pop(0) if len(answers[name]) > 1 else answers[name][0]

    def restart_mysql(environ):
        restarted.append("mysql")
        return restart

    def restart_redis(environ):
        restarted.append("redis")
        return restart

    tripped = []
    dog = services.Watchdog({}, tripped.append, interval=0.01, mysql=ask("mysql"), redis=ask("redis"),
                            restart_mysql=restart_mysql, restart_redis=restart_redis, say=lambda text: None, **options)
    return dog, restarted, tripped


def test_a_server_that_answers_is_left_alone():
    dog, restarted, _ = _watchdog()
    assert dog.check() == "" and restarted == [] and dog.events == []


def test_a_mysql_that_is_gone_is_restarted_at_once():
    dog, restarted, _ = _watchdog(mysql=("down",))
    found = dog.check()
    assert restarted == ["mysql"] and "MySQL stopped answering and was restarted" in found and not dog.failed


def test_a_mysql_that_only_holds_its_port_must_miss_several_answers_before_it_is_restarted():
    dog, restarted, _ = _watchdog(mysql=("hung",))
    for _ in range(services.HUNG_AFTER - 1):
        assert dog.check() == ""
    assert restarted == []
    assert "restarted" in dog.check() and restarted == ["mysql"]
    dog, restarted, _ = _watchdog(mysql=("hung", "hung", "up", "hung", "hung"))
    assert [dog.check() for _ in range(5)] == [""] * 5 and restarted == []     # a good answer starts the count over


def test_redis_is_restarted_too_and_a_failed_restart_is_remembered():
    dog, restarted, _ = _watchdog(redis=("down",), restart=False)
    found = dog.check()
    assert restarted == ["redis"] and "Redis stopped answering and was NOT restarted" in found and dog.failed


def test_the_loop_trips_once_and_tells_the_runner():
    dog, restarted, tripped = _watchdog(mysql=("up", "up", "down"))
    dog.start()
    assert dog.tripped.wait(5)
    dog.stop()
    assert restarted == ["mysql"] and len(tripped) == 1 and "MySQL" in dog.what


# --- a phase under the watchdog ---------------------------------------------------------

class _Fake:
    """`scripts/test.py` with the watchdog and the pytest run replaced: attempt N does what `attempts[N]` says."""

    def __init__(self, monkeypatch, tmp_path, attempts, passed_per_attempt):
        self.commands, self.dogs, self.skips = [], [], []
        self.progress = tmp_path / "progress.jsonl"
        self.progress.write_text(json.dumps({"kind": "header"}) + "\n")
        self.attempts, self.passed_per_attempt = list(attempts), list(passed_per_attempt)
        fake = self

        class Dog:
            def __init__(self, environ, on_trip, **options):
                self.on_trip, self.events, self.what, self.failed = on_trip, [], "", False
                self.tripped = threading.Event()
                fake.dogs.append(self)

            def start(self):
                pass

            def stop(self):
                pass

        class Services:
            Watchdog = Dog

        monkeypatch.setattr(runner, "_services", lambda: Services)
        monkeypatch.setattr(runner, "run_phase", self.run_phase)

    def run_phase(self, command, environment, on_child=None):
        self.commands.append(command)
        if "--skip-passed" in command:      # the list is rewritten for each attempt: keep what this one was given
            with open(command[command.index("--skip-passed") + 1], encoding="utf-8") as skip:
                self.skips.append(json.load(skip))
        dog = self.dogs[-1]
        outcome = self.attempts.pop(0)
        with open(self.progress, "a", encoding="utf-8") as log:
            log.writelines(json.dumps({"kind": "test", "nodeid": nodeid, "outcome": "passed"}) + "\n"
                           for nodeid in self.passed_per_attempt.pop(0))
        if outcome == "restarted":
            dog.events.append("MySQL stopped answering and was restarted at 01:02:03Z")
            dog.what = dog.events[-1]
            dog.tripped.set()
            return -15, {}
        if outcome == "unrestartable":
            dog.events.append("MySQL stopped answering and was NOT restarted at 01:02:03Z")
            dog.what, dog.failed = dog.events[-1], True
            dog.tripped.set()
            return -15, {}
        return 0, {"passed": 4}


def test_a_phase_cut_short_by_a_restart_goes_on_without_the_tests_that_passed(monkeypatch, tmp_path):
    fake = _Fake(monkeypatch, tmp_path, ["restarted", "restarted", "done"],
                 [["a::t1", "a::t2"], ["a::t3"], ["a::t4"]])
    said = []
    code, counts, restarts = runner.run_watched("non-browser", "not browser", 2, ["--changed", "main"], {},
                                                str(fake.progress), say=said.append)
    assert code == 0 and len(restarts) == 2
    assert counts == {"passed": 4 + 2 + 1}                      # this attempt's, plus the two that passed in each stopped one
    assert fake.commands[0][-2:] == ["--changed", "main"]
    assert all("--changed" in command and command.count("--skip-passed") == 1 for command in fake.commands[1:])
    assert fake.skips == [["a::t1", "a::t2"], ["a::t1", "a::t2", "a::t3"]]
    assert len(said) == 2 and "restart 1 of at most" in said[0] and "without the 2 tests" in said[0]


def test_an_earlier_skip_list_is_kept_when_the_phase_starts_again(monkeypatch, tmp_path):
    earlier = tmp_path / "earlier.skip"
    earlier.write_text(json.dumps(["old::t0"]))
    fake = _Fake(monkeypatch, tmp_path, ["restarted", "done"], [["a::t1"], []])
    runner.run_watched("p", "m", 1, ["--skip-passed", str(earlier)], {}, str(fake.progress), say=lambda text: None)
    assert fake.skips == [["old::t0"], ["a::t1", "old::t0"]]


def test_a_phase_that_finished_before_the_stop_is_not_run_again(monkeypatch, tmp_path):
    fake = _Fake(monkeypatch, tmp_path, ["done"], [[]])
    original = fake.run_phase

    def finishing(command, environment, on_child=None):
        result = original(command, environment, on_child)
        fake.dogs[-1].tripped.set()                             # the stop came after the summary line was printed
        return result

    monkeypatch.setattr(runner, "run_phase", finishing)
    code, counts, _ = runner.run_watched("p", "m", 1, [], {}, str(fake.progress), say=lambda text: None)
    assert (code, counts) == (0, {"passed": 4}) and len(fake.commands) == 1


def test_a_server_that_cannot_be_restarted_ends_the_run_with_the_reason(monkeypatch, tmp_path):
    _Fake(monkeypatch, tmp_path, ["unrestartable"], [[]])
    with pytest.raises(runner.ServicesLost, match="NOT restarted"):
        runner.run_watched("p", "m", 1, [], {}, str(tmp_path / "progress.jsonl"), say=lambda text: None)


def test_a_phase_that_keeps_being_cut_short_is_given_up(monkeypatch, tmp_path):
    attempts = ["restarted"] * (runner.MAX_RESTARTS + 2)
    _Fake(monkeypatch, tmp_path, attempts, [[] for _ in attempts])
    with pytest.raises(runner.ServicesLost, match="something keeps killing them"):
        runner.run_watched("p", "m", 1, [], {}, str(tmp_path / "progress.jsonl"), say=lambda text: None)


# --- the command ---------------------------------------------------------------------------

def test_a_run_whose_preflight_fails_is_not_started_and_says_why(monkeypatch, capsys):
    class Broken:
        @staticmethod
        def preflight(environ, browser=True, say=print):
            return ["MySQL does not answer and could not be started"]

    ran = []
    monkeypatch.setattr(runner, "_services", lambda: Broken)
    monkeypatch.setattr(runner, "run_phase", lambda *args, **kwargs: ran.append(args) or (0, {}))
    assert runner.main(["quick"]) == 2 and ran == []
    out = capsys.readouterr().out
    assert "PREFLIGHT FAILED: MySQL does not answer" in out and "Nothing was run" in out


def test_the_env_command_only_checks(monkeypatch, capsys):
    ran = []
    monkeypatch.setattr(runner, "run_phase", lambda *args, **kwargs: ran.append(args) or (0, {}))
    assert runner.main(["env"]) == 0 and ran == []
    assert "in order" in capsys.readouterr().out


def test_no_preflight_runs_without_the_check_and_without_the_watchdog(monkeypatch):
    class Never:
        @staticmethod
        def preflight(*args, **kwargs):
            pytest.fail("checked")

        Watchdog = None

    seen = []
    monkeypatch.setattr(runner, "_services", lambda: Never)
    monkeypatch.setattr(runner, "run_phase", lambda command, environment, **options: seen.append(command) or (0, {"passed": 1}))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    assert runner.main(["quick", "--no-preflight"]) == 0 and len(seen) == 1


def test_the_closing_line_counts_the_restarts(monkeypatch, capsys):
    monkeypatch.setattr(runner, "run_watched", lambda *args, **kwargs: (0, {"passed": 3}, ["MySQL stopped answering and was restarted at 01:02:03Z"]))
    monkeypatch.setattr(runner, "save_run", lambda *args: None)
    assert runner.main(["quick"]) == 0
    last = capsys.readouterr().out.strip().splitlines()[-1]
    assert last.startswith("scripts/test.py quick: PASS") and "1 service restart(s): MySQL stopped answering" in last


def test_a_lost_service_stops_the_run_with_a_resume_hint(monkeypatch, capsys):
    def lost(*args, **kwargs):
        raise runner.ServicesLost("MySQL stopped answering and was NOT restarted")

    monkeypatch.setattr(runner, "run_watched", lost)
    assert runner.main(["quick"]) == 3
    out = capsys.readouterr().out
    assert "FAIL" in out and "NOT restarted" in out and "--resume" in out
