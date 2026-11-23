"""OPS.293: a follower that waits for units asks for a supervisor, and stops with the reason when none comes up (it waited for ever with no word: the supervisor exits when idle)."""

import pytest

from planetgen.queue import supervisor, units, wrapper


class Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


def _pump(monkeypatch, alive):
    asked = []
    monkeypatch.setattr(wrapper, "SUPERVISOR_WAIT_SECONDS", 60.0)        # the conftest turns the wait off for the other tests
    monkeypatch.setattr(wrapper, "ask_supervisor", lambda conn, url: asked.append(url))
    monkeypatch.setattr(wrapper.supervisor, "is_alive", lambda conn: alive())
    clock = Clock()
    return wrapper.supervisor_pump(object(), "redis://x", clock), clock, asked


def test_a_follower_asks_for_a_supervisor_at_every_call_and_says_nothing_while_one_is_alive(monkeypatch):
    pump, clock, asked = _pump(monkeypatch, lambda: True)
    for _ in range(5):
        clock.now += 600
        pump()
    assert asked == ["redis://x"] * 5


def test_a_follower_with_no_supervisor_stops_after_the_wait_with_the_reason_the_start_failed(monkeypatch):
    pump, clock, asked = _pump(monkeypatch, lambda: False)
    monkeypatch.setattr(supervisor, "_last_failure", "OSError: [Errno 12] Cannot allocate memory")
    pump()
    clock.now += wrapper.SUPERVISOR_WAIT_SECONDS - 1
    pump()                                  # a supervisor takes seconds to come up: not yet
    clock.now += 2
    with pytest.raises(units.RunStopped) as stopped:
        pump()
    text = str(stopped.value)
    assert "no supervisor has been alive for 1m 1s" in text and "Cannot allocate memory" in text and "planetgen.cli.supervisor" in text
    assert len(asked) == 3


def test_a_supervisor_that_comes_up_resets_the_wait(monkeypatch):
    state = {"alive": False}
    pump, clock, _asked = _pump(monkeypatch, lambda: state["alive"])
    pump()
    clock.now += 50
    state["alive"] = True
    pump()
    state["alive"] = False
    clock.now += 50
    pump()                                  # the wait starts again from here
    clock.now += 50
    pump()


def test_ensure_started_remembers_why_it_could_not_start_a_supervisor(monkeypatch):
    class Redis:
        def exists(self, key):
            return 0

        def set(self, key, value, nx=False):
            return True

        def delete(self, key):
            return 1

        def get(self, key):
            return None

    def refuse(argv):
        raise OSError("fork refused")

    assert supervisor.ensure_started(Redis(), start=refuse) is False
    assert "fork refused" in supervisor.last_failure()
    assert supervisor.ensure_started(Redis(), start=lambda argv: None) is True
    assert supervisor.last_failure() is None
