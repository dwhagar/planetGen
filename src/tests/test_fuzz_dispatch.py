"""OPS.236: the dispatch rule (`planetgen/queue/dispatch.py`) over seeds and backlogs of the simulator (`tests/enginesim.py`, OPS.227): ten seeds at 30 and 300 units and two at 3,000, each with a
visitor's jobs arriving through the stage, a claim of 7,000 MB, a restarted unit that runs alone and a cut that needs the bulk lane to itself.

Every run checks in each turn that the claims admitted never come to more than the host less its reserve (but for one claim that is larger and runs alone), that a token is held by one job, that
the job that runs alone runs alone, that the kernel kills only when a job went past its claim or something outside the engine took the memory, and at the end that every job finished once.
The harder runs also have a decider that fails in one turn of ten (a failed turn admits nothing), claims that are too small, and something outside the engine holding 9.5 GB of the host for
a while, so the kernel kills jobs and each runs again alone on a claim grown by half.
"""

import pytest

from tests import enginesim

SITE_SECONDS = 1.5
SITE_SECONDS_WITH_FAILED_TURNS = 5.0
"""float: The 99th percentile of the wait of a visitor's job that does not come while something outside holds the memory (one that does waits for the memory): a turn at a time (one second) when the decider answers, a few more when it fails in one turn of ten."""

START_AFTER_SECONDS = 120.0
"""float: How long after every job of the bulk and maintenance lanes queued before it had finished the big claim, the restarted unit and the cut may start: later arrivals do not hold them back."""

PRESSURE = {30: [(60, 160, 9500.0)], 300: [(400, 700, 9500.0), (1500, 1800, 9500.0)], 3000: [(3000, 3300, 9500.0), (12000, 12300, 9500.0)]}
"""dict: Per backlog, when something outside the engine holds 9.5 GB of the host (from, to, MB; seconds of the simulated clock)."""


def _check(measures, broken, site):
    assert not broken, broken[:5]
    assert measures["wait_p99_site_calm"] <= site, measures
    for name in ("big-claim", "restarted-unit", "cut-next-stage"):
        assert measures[f"late_{name}"] <= START_AFTER_SECONDS, (name, measures)


@pytest.mark.parametrize("units", [30, 300])
@pytest.mark.parametrize("seed", range(10))
def test_in_every_seed_the_site_goes_at_once_and_the_big_claim_the_lone_unit_and_the_cut_all_start(seed, units):
    measures, broken = enginesim.run(seed, units, named=True, big=True)
    _check(measures, broken, SITE_SECONDS)


@pytest.mark.parametrize("units", [30, 300])
@pytest.mark.parametrize("seed", range(10))
def test_with_kills_small_claims_and_a_decider_that_fails_every_seed_still_keeps_its_budget_and_starts_everything(seed, units):
    measures, broken = enginesim.run(seed, units, named=True, big=True, crash=0.1, small=0.1, pressure=PRESSURE[units])
    _check(measures, broken, SITE_SECONDS_WITH_FAILED_TURNS)
    assert measures["crashed_turns"] > 0


def test_the_memory_pressure_kills_jobs_and_each_runs_again_once_alone():
    kills = 0
    for seed in range(3):
        measures, broken = enginesim.run(seed, 300, named=True, big=True, crash=0.1, small=0.1, pressure=PRESSURE[300])
        _check(measures, broken, SITE_SECONDS_WITH_FAILED_TURNS)
        kills += measures["oom_kills"]
    assert kills > 0, "the pressure windows kill nothing, so the retry path is not run"


@pytest.mark.parametrize("seed", [0, 1])
def test_a_backlog_of_three_thousand_units_is_served_the_same_way(seed):
    measures, broken = enginesim.run(seed, 3000, named=True, big=True, crash=0.1, small=0.1, pressure=PRESSURE[3000])
    _check(measures, broken, SITE_SECONDS_WITH_FAILED_TURNS)


def _admit_everything(waiting, host, held, now, used, window, places):
    return enginesim.dispatch.Turn([ticket.name for ticket in waiting], {})


def test_the_invariants_catch_a_rule_that_admits_every_job_and_one_that_never_admits_the_big_claim():
    """The sweep is only worth what its checks catch: a rule that lets everything go breaks the budget, and one that leaves a claim out leaves it unfinished."""
    _measures, broken = enginesim.run(0, 30, named=True, big=True, rule=_admit_everything)
    assert any("over the budget" in text or "killed" in text for text in broken), broken[:3]

    def never_big(waiting, host, held, now, used, window, places):
        return enginesim.dispatch.dispatch([ticket for ticket in waiting if ticket.name != "big-claim"], host, held, now, used, window, places)

    _measures, broken = enginesim.run(0, 30, named=True, big=True, rule=never_big, horizon=3000.0)
    assert "big-claim finished 0 times" in broken
