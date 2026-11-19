"""OPS.246: the kill matrix over a chain. One seeded galaxy run (the sectors, the charting that draws the sky, the settling of the paths) and the warm-up of the object index after it, with units killed at drawn commits all along the chain and the next attempts carrying on, ends as an undisturbed chain ends: the same galaxy, stars, cells and spheres (`test_engine_kill_backfill._digest`) and an index that holds every object once.

The gate runs a small fixed number of seeds. The end-of-version hunt runs it in full: `PLANETGEN_ENGINE_CHAIN_SEEDS=200` seeds and `PLANETGEN_ENGINE_CHAIN_KILLS=40` kills at most to a seed (the default is 2 seeds and 4 kills). Replay a seed with `pytest -k seed_3`.

`TRANSITIONS` names each transition of the engine and the test that kills it; the last test checks every one of them is a test that exists.
"""

import os
import random

import pytest

from planetgen.db import objectindex, store
from tests import kill_harness
from tests.test_engine_kill_backfill import _digest, _generate, _wipe
from tests.test_engine_kill_object_index import _kill_at_the_writes, _warm
from tests.test_sky_sliced_draw import PARTS, toy_tiers  # noqa: F401 (the fixture)

pytestmark = [pytest.mark.fresh_database, pytest.mark.slow]

SEEDS = ["control", *range(int(os.environ.get("PLANETGEN_ENGINE_CHAIN_SEEDS", "2")))]
MOST_KILLS = int(os.environ.get("PLANETGEN_ENGINE_CHAIN_KILLS", "4"))

TRANSITIONS = {
    "claim and start of a unit": ("test_engine_kill_processes.py", "test_a_job_killed_before_it_took_its_attempt_is_handed_out_again_and_finished_once"),
    "commit of a unit's piece": ("test_engine_kill_chain.py", "test_a_chain_run_with_killed_units_ends_as_an_undisturbed_chain_ends"),
    "end of a unit (its attempt written)": ("test_engine_kill_processes.py", "test_a_supervisor_killed_after_its_job_wrote_its_end_has_nothing_to_recover"),
    "hold of a unit that fails for good": ("test_engine_random_runs.py", "test_a_run_drawn_by_a_seed_ends_done_with_every_piece_once_and_every_stage_in_order"),
    "retry of a held run": ("test_engine_random_runs.py", "test_a_run_drawn_by_a_seed_ends_done_with_every_piece_once_and_every_stage_in_order"),
    "cancel and finish of a run": ("test_engine_random_runs.py", "test_a_run_drawn_by_a_seed_ends_done_with_every_piece_once_and_every_stage_in_order"),
    "end of a stage": ("test_engine_random_runs.py", "test_a_run_drawn_by_a_seed_ends_done_with_every_piece_once_and_every_stage_in_order"),
    "write of a window": ("test_engine_kill_launcher_windows.py", "test_every_commit_of_a_window_has_a_kill_that_leaves_the_stage_written_whole"),
    "launcher hand-over (a supervisor killed)": ("test_engine_kill_processes.py", "test_a_supervisor_killed_with_a_job_running_leaves_a_job_that_stops_at_its_piece_and_the_next_supervisor_finishes_it"),
    "template killed": ("test_engine_kill_processes.py", "test_a_template_and_its_job_both_killed_have_the_job_handed_out_again_once_and_finished"),
    "job killed by the kernel (claim grown)": ("test_engine_kill_processes.py", "test_a_job_the_kernel_kills_runs_again_alone_on_a_claim_grown_by_one_and_a_half_and_finishes_once"),
}
"""dict: Each named transition of the engine, and the (file, test) that kills it."""


def _chain(config):
    _generate(config)
    _warm(config)


def _holds_every_object_once(config):
    conn = store.get_connection(config)
    try:
        report = objectindex.check(conn)
        return all(not (kind["missing"] or kind["extra"] or kind["different"]) for kind in report.values())
    finally:
        conn.close()


@pytest.fixture
def undisturbed(mysql_config, toy_tiers):  # noqa: F811
    import _pytest.monkeypatch
    counter = kill_harness.Plan()
    patch = _pytest.monkeypatch.MonkeyPatch()
    try:
        _kill_at_the_writes(patch, counter)
        _chain(mysql_config)
    finally:
        patch.undo()
    return _digest(mysql_config), counter.seen


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_a_chain_run_with_killed_units_ends_as_an_undisturbed_chain_ends(mysql_config, toy_tiers, undisturbed, monkeypatch, seed):  # noqa: F811
    expected, commits = undisturbed
    assert commits > 0 and expected["bright_stars"][0] > 0 and _holds_every_object_once(mysql_config), "the undisturbed chain made a galaxy and an index to compare with"
    rng = random.Random(str(seed))
    kills = {} if seed == "control" else {number: rng.choice(("before", "after")) for number in rng.sample(range(commits), min(commits, rng.randint(1, MOST_KILLS)))}
    _wipe(mysql_config)
    plan = kill_harness.Plan(kills)
    _kill_at_the_writes(monkeypatch, plan)
    _chain(mysql_config)
    assert plan.fired or seed == "control", f"seed {seed}: no kill happened (commits {commits}, plan {kills})"
    assert _digest(mysql_config) == expected, f"seed {seed}: kills {plan.fired} changed the galaxy"
    assert _holds_every_object_once(mysql_config), f"seed {seed}: kills {plan.fired} left the index different from the tables"


def test_every_named_transition_of_the_engine_is_killed_by_a_test_that_exists():
    here = os.path.dirname(__file__)
    for name, (file, test) in TRANSITIONS.items():
        with open(os.path.join(here, file), encoding="utf-8") as handle:
            assert f"def {test}(" in handle.read(), f"{name}: {file} has no test {test}"
