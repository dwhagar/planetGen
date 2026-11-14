"""OPS.194: a unit whose counter stands still but whose beat shows work is busy, never quiet, and nothing is stopped (`planetgen/queue/beatjudge.py`; the evidence is OPS.193)."""

from planetgen.queue import beatjudge, unitwatch
from planetgen.queue.beatjudge import Beat, History

NOW = 10_000.0
HISTORY = History(longest_gap=20.0, expected_seconds=600.0)        # the quiet limit is 300 s


def _beat(**fields):
    base = dict(beat=NOW - 3.0, started=NOW - 1000.0, done=5.0, estimate=100.0, moved=NOW - 900.0)
    base.update(fields)
    return Beat(**base)


def test_a_still_counter_with_fresh_evidence_is_busy_and_without_it_quiet():
    quiet = beatjudge.judge(_beat(), HISTORY, NOW)
    assert quiet.state == beatjudge.QUIET and quiet.actionable
    busy = beatjudge.judge(_beat(active=NOW - 8.0, sql="INSERT INTO scatter_masses"), HISTORY, NOW)
    assert busy.state == beatjudge.BUSY and not busy.actionable
    assert "in SQL: INSERT INTO scatter_masses" in busy.reason and "900 s" in busy.reason
    stale = beatjudge.judge(_beat(active=NOW - 120.0), HISTORY, NOW)
    assert stale.state == beatjudge.QUIET, "work seen minutes ago is no proof of work now"


def test_a_busy_unit_is_a_flag_only_past_the_learned_seconds_and_is_never_acted_on():
    inside = beatjudge.judge(_beat(started=NOW - 500.0, moved=NOW - 400.0, active=NOW - 5.0), HISTORY, NOW)
    assert inside.state == beatjudge.BUSY and inside.severity == ""
    past = beatjudge.judge(_beat(active=NOW - 5.0), HISTORY, NOW)
    assert past.state == beatjudge.BUSY and past.severity == beatjudge.FLAG and "learned" in past.reason
    unknown = beatjudge.judge(_beat(active=NOW - 5.0), History(longest_gap=20.0), NOW)
    assert unknown.severity == "", "nothing is learned for the command, so nothing is flagged"
    for verdict in (inside, past, unknown):
        assert not verdict.actionable and not beatjudge.agreed(verdict, verdict)


def test_a_counter_that_moved_recently_does_not_need_the_evidence_and_a_wait_is_still_a_wait():
    assert beatjudge.judge(_beat(done=20.0, started=NOW - 100.0, moved=NOW - 30.0, active=NOW - 2.0), HISTORY, NOW).state == beatjudge.WORKING
    assert beatjudge.judge(_beat(waiting_on="gate", active=NOW - 2.0), HISTORY, NOW).state == beatjudge.WAITING
    assert beatjudge.judge(_beat(beat=NOW - 60.0, active=NOW - 61.0), HISTORY, NOW).state == beatjudge.DEAD, "no beat is dead whatever the evidence said"


def test_the_page_names_busy_and_the_published_view_carries_the_statement():
    from planetgen.web import queue_page
    assert "busy" in queue_page.UNIT_STATE_LABELS and "busy" in beatjudge.STATES
    beat = _beat(active=NOW - 5.0, sql="OPTIMIZE TABLE x", sql_since=NOW - 90.0)
    shown = unitwatch.view(beatjudge.judge(beat, HISTORY, NOW), beat, NOW)
    assert shown["state"] == "busy" and shown["sql"] == "OPTIMIZE TABLE x" and shown["sql_seconds"] == 90.0
