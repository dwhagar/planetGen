# tests/test_backfill_edges.py

"""
OPS.242, the narrower edges: the backfill of a sector may start when that sector is saved (and the sectors of its list that the run saves are saved), while others are still being saved (`tests/backfill_edge_support.py`).
"""

import pytest

from planetgen import tuning
from tests.backfill_edge_support import prove
from tests.test_bright_star_scatter import LIGHT_STEP_MASSES

pytestmark = [pytest.mark.db, pytest.mark.usefixtures("large_spiral_fallback")]


@pytest.fixture(autouse=True)
def light_steps(monkeypatch):
    """The four-step walk of the scatter tests, so the toy galaxy's backfill stays quick."""
    monkeypatch.setattr(tuning, "BRIGHT_STAR_BACKFILL_RING_MASSES_SOL", LIGHT_STEP_MASSES)


def test_the_backfill_by_address_stores_what_plan_order_stores(mysql_config, redis_server):
    verdict = prove(mysql_config)
    print(verdict.line(), flush=True)
    assert verdict.passed, verdict.line() + "".join(f"\n    {line[:400]}" for line in [*verdict.undeclared, *verdict.differences[:20]])


def test_a_backfill_released_before_the_sector_it_lists_is_saved_stores_other_bytes(mysql_config, redis_server):
    """The control: a unit released at once draws stars into a sector of the run before it is saved, which is not what plan order stores."""
    verdict = prove(mysql_config, scope="none")
    print(verdict.line(), flush=True)
    assert verdict.differences, verdict.line()
