# tests/test_address_edges.py

"""
OPS.242, the narrower edges: link and settle of a sector may start when the sectors they read are saved, while others are still being saved (`tests/address_edge_support.py`).
"""

import pytest

from tests.address_edge_support import prove

pytestmark = [pytest.mark.db, pytest.mark.usefixtures("large_spiral_fallback")]


def test_link_and_settle_by_address_store_what_plan_order_stores(mysql_config, redis_server):
    verdict = prove(mysql_config)
    print(verdict.line(), flush=True)
    assert verdict.passed, verdict.line() + "".join(f"\n    {line[:400]}" for line in [*verdict.undeclared, *verdict.differences[:20]])


def test_a_scope_of_the_sector_alone_is_not_enough(mysql_config, redis_server):
    """The control: a unit released when its own sector is saved, before the sectors around it, stores other bytes than plan order, so the toy sees what the neighbors are for."""
    verdict = prove(mysql_config, scope="self")
    print(verdict.line(), flush=True)
    assert verdict.differences, verdict.line()
