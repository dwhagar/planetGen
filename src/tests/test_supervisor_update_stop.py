"""OPS.298: a supervisor stops for an update it is not part of, and stays for the update's own commands (`check-math`, `warm`)."""

from planetgen.admin import updateguard
from planetgen.cli import supervisor as cli


def test_a_supervisor_outside_the_update_stops_when_one_runs(monkeypatch):
    monkeypatch.setattr(updateguard, "is_updating", lambda: True)
    monkeypatch.setattr(updateguard, "inside_update", lambda: False)
    assert cli.should_stop([]) is True


def test_a_supervisor_the_update_started_stays(monkeypatch):
    monkeypatch.setattr(updateguard, "is_updating", lambda: True)
    monkeypatch.setattr(updateguard, "inside_update", lambda: True)
    assert cli.should_stop([]) is False


def test_a_signal_stops_it_in_any_case(monkeypatch):
    monkeypatch.setattr(updateguard, "is_updating", lambda: False)
    monkeypatch.setattr(updateguard, "inside_update", lambda: False)
    assert cli.should_stop([True]) is True
    assert cli.should_stop([]) is False
