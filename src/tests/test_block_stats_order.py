"""
TEST.201: the block totals (`galaxy_block_stats`) do not depend on the order of the rows they add up (docs/design/architecture/numeric-precision-audit.md
section 5). The same sectors are added in several random orders, by `block_stats.accumulate` and through `rebuild` on a database whose sectors were inserted
in each order, and the totals are required to be the same: the counts exactly, the real sums within `REL_BOUND`, the relative error of a 64-bit sum of this many
terms (each added sum is exact to 2**-53 and a block adds at most a few thousand). A failure means a total depends on the order the sectors came in.
"""

import random

import pytest

from planetgen.db import block_stats, store
from tests.test_block_stats import _addresses

REL_BOUND = 1e-12
"""float: How far apart two orders' real totals may be, as a share of the total: a few thousand 64-bit additions are good to about 1e-13."""

ORDERS = 6
"""int: How many random orders are tried beside the sorted one."""

REAL_COLUMNS = ("expected_systems", "age_sum", "luminosity_sol", "mass_sum_kg")


def _figures(addresses, seed=3):
    """Each sector's figures, drawn with the full 53 bits so that an order-dependent sum shows."""
    rng = random.Random(seed)
    return {address: dict(systems=rng.randrange(0, 9), expected=rng.uniform(0.5, 6.0), stars=rng.randrange(1, 40), mean_age=rng.uniform(0.01, 12.0),
                          luminosity=rng.uniform(1.0, 5e4) * 10 ** rng.uniform(-3, 3), habitable=rng.randrange(0, 3), mass=rng.uniform(1e30, 1e34))
            for address in addresses}


def _orders(addresses):
    rng = random.Random(11)
    yield list(addresses)
    for _ in range(ORDERS):
        shuffled = list(addresses)
        rng.shuffle(shuffled)
        yield shuffled


def _close(first, second):
    return abs(first - second) <= REL_BOUND * max(abs(first), abs(second), 1e-300)


def _accumulated(order, figures):
    cells = {}
    for address in order:
        f = figures[address]
        block_stats.accumulate(cells, address, f["systems"], f["expected"], f["stars"], f["mean_age"], f["luminosity"], f["habitable"], f["mass"])
    return cells


def test_accumulating_the_same_sectors_in_any_order_gives_the_same_block_sums():
    addresses = _addresses(count=70)
    figures = _figures(addresses)
    first = None
    for order in _orders(addresses):
        cells = _accumulated(order, figures)
        if first is None:
            first = cells
            continue
        assert cells.keys() == first.keys()
        for block, sums in cells.items():
            reference = first[block]
            assert [sums[0], sums[1], sums[3], sums[6]] == [reference[0], reference[1], reference[3], reference[6]], block
            for index in (2, 4, 5, 7):
                assert _close(sums[index], reference[index]), (block, index, sums[index], reference[index])


def test_a_block_with_many_sectors_of_very_different_size_adds_up_in_any_order():
    """One block that holds every sector, the luminosities spanning six orders of magnitude: the worst case for a plain float sum."""
    rng = random.Random(5)
    luminosities = [10 ** rng.uniform(-3, 6) for _ in range(3000)]
    reference = sum(sorted(luminosities))
    for _ in range(ORDERS):
        shuffled = list(luminosities)
        rng.shuffle(shuffled)
        assert _close(sum(shuffled), reference), "a plain 64-bit sum of 3000 terms of six decades is order dependent beyond the bound"


@pytest.mark.db
def test_rebuild_stores_the_same_totals_whatever_order_the_sectors_were_inserted_in(mysql_config):
    addresses = _addresses(count=70)
    figures = _figures(addresses)
    conn = store.get_connection(mysql_config)
    stored = []
    try:
        for order in _orders(addresses):
            with conn:
                conn.execute("DELETE FROM sector_stats")
                conn.execute("DELETE FROM sectors")
                for n, (shell, band, wedge) in enumerate(order):
                    f = figures[(shell, band, wedge)]
                    conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES (?, 4000, ?, ?, ?)",
                                 (f"Sector {n}", shell, band, wedge))
                    conn.execute(
                        "INSERT INTO sector_stats (shell_index, band_index, wedge_index, expected_systems, actual_systems, actual_stars, mean_age_gy,"
                        " total_luminosity_sol, bright_level_sol, habitable_worlds, total_mass_kg) VALUES (?, ?, ?, ?, ?, ?, ?, ?, 0, ?, ?)",
                        (shell, band, wedge, f["expected"], f["systems"], f["stars"], f["mean_age"], f["luminosity"], f["habitable"], f["mass"]))
                block_stats.rebuild(conn)
            rows = conn.execute("SELECT * FROM galaxy_block_stats ORDER BY m, shell_index, band_index, wedge_index").fetchall()
            stored.append([dict(row) for row in rows])
    finally:
        conn.close()
    reference = stored[0]
    assert reference
    for rows in stored[1:]:
        assert len(rows) == len(reference)
        for row, ref in zip(rows, reference):
            for column, value in ref.items():
                if column in REAL_COLUMNS:
                    assert _close(row[column], value), (column, row, ref)
                else:
                    assert row[column] == value, (column, row, ref)
