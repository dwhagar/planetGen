"""OPS.288: a full disk holds the run with a plain reason instead of ending it as "its process was gone" (`queue/diskgate.py`). The measuring and the controller's timing are tested with the drives made up; the hold of
a live run and its carry on use the real ledger and Redis, as `test_held_runs.py` does."""

import collections

import pytest

from planetgen.queue import diskgate, held, jobview, planrows, pressure, runrecord
from tests.test_held_runs import _job_state, _run_row, rows, run  # noqa: F401 -- fixtures
from tests.test_queue_page import conn, control_config, web_app   # noqa: F401 -- the fixtures

GB = 1024 ** 3
Usage = collections.namedtuple("Usage", "total used free")


class Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


def test_the_reserve_is_a_small_share_of_the_drive_with_a_floor_and_a_ceiling():
    assert diskgate.reserve_bytes(70 * GB) == diskgate.RESERVE_MIN_BYTES           # 0.5 percent of 70 GB is 358 MB, so the floor
    assert diskgate.reserve_bytes(1000 * GB) == int(0.005 * 1000 * GB)             # a big drive keeps more
    assert diskgate.reserve_bytes(2 * GB) == int(0.1 * 2 * GB)                     # a small one keeps no more than a tenth


@pytest.mark.parametrize("free_gb, level", [(100, pressure.CALM), (1.5, pressure.CALM), (0.9, pressure.PRESSED), (0.3, pressure.STARVED), (0, pressure.STARVED)])
def test_free_space_is_pressed_under_two_reserves_and_starved_under_one(free_gb, level):
    assert diskgate.level_of(int(free_gb * GB), 70 * GB) == level


def test_the_worst_drive_is_the_reading_and_one_drive_is_counted_once(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    figures = iter([Usage(70 * GB, 0, 40 * GB)])
    reading = diskgate.measure([("the database", str(tmp_path / "a")), ("the state folder", str(tmp_path / "b"))], usage=lambda path: next(figures))
    assert reading.where == "the database" and reading.level == pressure.CALM, "both places are on one drive: it is measured once"


def test_a_place_that_is_not_made_yet_is_measured_at_its_nearest_folder(tmp_path):
    reading = diskgate.measure([("the state folder", str(tmp_path / "not" / "made" / "yet"))], usage=lambda path: Usage(70 * GB, 70 * GB - GB // 4, GB // 4))
    assert reading.level == pressure.STARVED and reading.path == str(tmp_path) and "the state folder" in reading.words()
    assert diskgate.measure([]).level == pressure.CALM, "nothing to measure is calm"


def test_the_controller_holds_admission_when_a_drive_is_pressed_and_lets_go_when_it_is_calm():
    clock, free, said = Clock(), [GB // 2], []
    controller = diskgate.Controller(None, look=lambda: diskgate.Reading("the database", "/var/lib/mysql/", free[0], 70 * GB, diskgate.RESERVE_MIN_BYTES, diskgate.level_of(free[0], 70 * GB)), clock=clock,
                                     drain=lambda conn, run_id, reason: False)
    controller.attach(type("Supervisor", (), {"act": staticmethod(lambda what, reason="": said.append((what, reason)))})())
    controller.turn()
    assert controller.holding() is None, "not yet: the hold needs the drive short for a while"
    clock.now += pressure.HOLD_AFTER_SECONDS
    controller.turn()
    assert "the disk is full" in controller.holding() and "free of" in controller.holding() and said[-1][0] == "stopped admitting bulk and upkeep jobs"
    free[0] = 40 * GB
    controller.turn()
    clock.now += pressure.RELEASE_AFTER_SECONDS
    controller.turn()
    assert controller.holding() is None and said[-1][0] == "admitting again"


def _starved(connection, clock):
    free = int(0.1 * GB)
    look = lambda: diskgate.Reading("the database", "/var/lib/mysql/", free, 70 * GB, diskgate.RESERVE_MIN_BYTES, pressure.STARVED)   # noqa: E731
    controller = diskgate.Controller(connection, look=look, clock=clock)
    controller.turn()
    clock.now += pressure.HOLD_AFTER_SECONDS
    controller.turn()
    return controller


def test_a_starved_drive_holds_the_live_run_with_a_plain_reason_and_only_a_person_carries_it_on(conn, run, rows):   # noqa: F811
    clock = Clock()
    controller = _starved(conn, clock)
    assert runrecord.run_state(conn, run) == "held" and runrecord.active_run(conn, runrecord.GENERATION) is None
    assert runrecord.read_run(conn, run)["error"] == diskgate.DISK_REASON == "the disk is full: free space, then Retry held or finish-run"
    assert _run_row(conn, run)["state"] == "held" and _run_row(conn, run)["held_reason"] == diskgate.DISK_REASON
    assert runrecord._s(conn.hget(runrecord.run_key(run), held.HELD_FOR)) == diskgate.HELD_FOR_DISK
    assert held.resume_after_update(conn) == [] and runrecord.run_state(conn, run) == "held", "the end of an update does not carry a run held for the disk on"
    controller.turn()
    assert controller.held_runs == [run], "a run is held once"
    assert held.retry_held(conn, run, who="a person") == {"requeued": 3}, "Retry held carries it on after space is freed"
    assert runrecord.run_state(conn, run) == "running"


def test_a_run_ended_with_its_process_gone_says_the_disk_is_full_when_it_is(monkeypatch):
    assert runrecord.gone_text() == runrecord.GONE_TEXT or "the disk is full" in runrecord.gone_text()
    monkeypatch.setattr(diskgate, "full_text", lambda: None)
    assert runrecord.gone_text() == runrecord.GONE_TEXT and "process was gone" in runrecord.gone_text()
    monkeypatch.setattr(diskgate, "full_text", lambda: f"{diskgate.DISK_REASON} (50 MB free)")
    assert runrecord.gone_text().startswith("the disk is full: free space, then Retry held or finish-run") and "process was gone" not in runrecord.gone_text()


def test_the_text_for_a_full_disk_is_only_given_when_a_drive_is_starved(monkeypatch):
    monkeypatch.setattr(diskgate, "places", lambda control=None, connection=None: [("the state folder", "/")])
    monkeypatch.setattr(diskgate.shutil, "disk_usage", lambda path: Usage(70 * GB, 70 * GB, 0))
    assert diskgate.full_text().startswith(diskgate.DISK_REASON)
    monkeypatch.setattr(diskgate.shutil, "disk_usage", lambda path: Usage(70 * GB, GB, 69 * GB))
    assert diskgate.full_text() is None
