"""
OPS.196, second part: no command goes silent. Each command of the clears, the drops and the reset (the last family, OPS.200 and OPS.206) is run on a small galaxy under the beat a worker gives a unit
(`unitbeat.UnitBeat`), and the longest stretch its count stood still (`longest_gap`, OPS.135: the pieces done, the second-level count, a phase or a statement's label moving it) is held to a budget.
A command that has put its work back into one loop of its own, or one piece that runs for minutes with no count, shows as a gap here before it shows as a stuck bar on the Queue page.
The registry guard is `test_commands_have_pieces.py`; the budget here is for a galaxy of a few sectors, so a gap over it is a command that waits for something that is not its work.
"""

import pytest

from planetgen.queue import commands, redisqueue, unitbeat, wrapper
from tests.test_api import _place_sector
from tests.test_galaxy_gen import _seed_skeleton

pytestmark = pytest.mark.fresh_database

BUDGET_SECONDS = 5.0
"""float: The longest a command of the list may stand without a count, a phase or a statement label moving, on a small galaxy."""

CASES = [
    ("clear-bright-stars", {"clear": "all"}),
    ("empty-phenomenon-scatter", {}),
    ("drop-phenomenon-indexes", {}),
    ("empty-nebula-field", {}),
    ("clear-scatter", {"part": "phenomena"}),
    ("clear-scatter", {"part": "nebula-field"}),
    ("clean-nebula-field", {"apply": False}),
    ("reset-galaxy", {}),
]
"""list: `(command, its arguments beside the database)`: what this test runs."""


@pytest.mark.parametrize("kind, extra", CASES, ids=[f"{kind}-{'-'.join(map(str, extra.values()))}".rstrip("-") for kind, extra in CASES])
def test_a_command_of_the_last_family_never_stands_without_a_count_past_the_budget(mysql_config, redis_server, kind, extra):
    _seed_skeleton(mysql_config, shells=[(0, 0), (1, 1)])
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    args = dict(extra, database=wrapper.database_args(mysql_config))
    with unitbeat.UnitBeat(redisqueue.connect(redis_server), "silence", 0, kind) as report:
        commands.run_alone(kind, args, report)
        gap = report.longest_gap()
    assert gap < BUDGET_SECONDS, f"{kind} stood {gap:.1f} s with no count, phase or statement moving (budget {BUDGET_SECONDS} s)"


def test_the_silence_meter_sees_a_unit_that_does_not_move(redis_server):
    """The meter is not a rubber stamp: a count that stands still for 30 s of a clock the test moves is a gap of 30 s."""
    now = [0.0]
    beat = unitbeat.UnitBeat(redisqueue.connect(redis_server), "silence", 0, "meter", clock=lambda: now[0])
    beat(1, 2)
    now[0] = 30.0
    beat(2, 2)
    assert beat.longest_gap() == 30.0
