# tests/test_rescan_proofs.py

"""
The proof tests of the generation-rescan chain (see `tests/rescan_proof_support.py`): that what GEN.1476 keeps, GEN.1479 skips and GEN.1477 shares is, to the last bit, what
the pass over the whole galaxy gave. The harness is tried here on the code as it stands (the existing incremental refresh against the rebuild; the settle of the same sectors
twice); each item's test is wired when its PR is up and tells its owner (GEN.1476 Lane 12, GEN.1479 Lane 12, GEN.1477 Lane 13).
"""

import pytest

from planetgen.db import edits, sector_paths, store
from planetgen.generation.run_plan import rescore_habitability_unit
from tests import rescan_proof_support as proof
from tests.fuzz_support import deterministic_entropy
from tests.test_sector_paths_db import _save

pytestmark = pytest.mark.db


# --- the harness itself -------------------------------------------------------------------------------------------------------

def test_a_float_that_differs_in_the_last_place_or_in_sign_is_a_different_row():
    one = [(proof.bits(0.1 + 0.2), 3)]
    assert proof.first_difference(one, [(proof.bits(0.3), 3)], ["x", "n"]) == f"row 0 differs in ['x']: {one[0]!r} against {[(proof.bits(0.3), 3)][0]!r}"
    assert proof.first_difference([(proof.bits(0.0),)], [(proof.bits(-0.0),)]) is not None
    assert proof.first_difference([(proof.bits(1.5), None)], [(proof.bits(1.5), None)]) is None


def test_a_missing_row_is_named():
    message = proof.first_difference([(1,), (2,)], [(1,)])
    assert message.startswith("2 rows against 1") and "(2,)" in message


def test_the_touch_log_counts_rows_a_statement_processes_even_when_it_changes_nothing(mysql_config):
    galaxy = proof.small_galaxy(mysql_config)
    with proof.TouchLog(mysql_config, ["star_systems"]) as log:
        with proof.connection(mysql_config) as conn:
            conn.execute("SELECT * FROM star_systems").fetchall()
            assert log.touched() == []                                   # a read touches nothing
            first = conn.execute("SELECT id, name FROM star_systems ORDER BY id LIMIT 1").fetchone()
            conn.execute("UPDATE star_systems SET name = ? WHERE id = ?", (first["name"], first["id"]))      # the value it already holds
            conn.commit()
        assert log.touched() == [("star_systems", "U", first["id"])]
        log.clear()
        assert log.touched() == []
    assert galaxy["sector_ids"]


# --- GEN.1476: the running block totals equal a rebuild, byte for byte -----------------------------------------------------------

def test_the_block_lights_equal_a_rebuild_after_a_sector_is_deleted_and_made_again(mysql_config):
    """Wired on the refresh the stored lights have now (a saved sector takes its stars out of its blocks, a deleted one puts them back); GEN.1476's scatter units and
    sector edits call `assert_totals_equal_a_rebuild` in the same three steps."""
    galaxy = proof.small_galaxy(mysql_config)
    assert proof.assert_totals_equal_a_rebuild(mysql_config, "the galaxy as built") > 0
    victim = galaxy["sector_ids"][1]
    with proof.connection(mysql_config) as conn:
        assert edits.delete_sector_with_contents(conn, victim) is not None
        conn.commit()
    proof.assert_totals_equal_a_rebuild(mysql_config, "after a sector is deleted")
    with deterministic_entropy(proof.SMALL_GALAXY_SEED + 1):
        _save(mysql_config, systems=2, address=galaxy["addresses"][1])
    proof.assert_totals_equal_a_rebuild(mysql_config, "after the sector is made again")


# --- GEN.1479: a job run again with no change touches no row ---------------------------------------------------------------------

@pytest.mark.xfail(strict=True, reason="GEN.1479 is not built: rescore-habitability recomputes and writes every scored body whether or not its formulas changed, so a second run "
                                       "touches all 34 scored bodies of the small galaxy (the touch log counts them). When GEN.1479's stamp lands this passes; remove this mark then (strict, so it cannot be forgotten).")
def test_a_rescore_run_again_with_nothing_changed_touches_no_row(mysql_config):
    proof.small_galaxy(mysql_config)
    with proof.connection(mysql_config) as conn:
        low, high = conn.execute("SELECT MIN(id) AS low, MAX(id) AS high FROM star_systems").fetchone().values()
    params = {"low": low, "high": high, "database": {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user, "database": mysql_config.database}}
    with proof.TouchLog(mysql_config, ["planets", "moons", "star_systems"]) as log:
        rescore_habitability_unit(params, lambda *a, **k: None)       # the first run may write what was never scored
        log.clear()
        rescore_habitability_unit(params, lambda *a, **k: None)
        assert log.touched() == []


# --- GEN.1477: settled paths do not depend on how the neighbourhood is read ---------------------------------------------------------

def test_the_paths_a_settle_saves_are_the_same_bytes_every_time_it_runs(mysql_config):
    """The harness on the current settle: the rows are as bytes and in an order that does not depend on row ids. GEN.1477's test settles the same sectors with the
    neighbourhood read per sector and shared across a chunk (the chunk size is the knob) and compares with `assert_same_rows`."""
    galaxy = proof.small_galaxy(mysql_config)
    with proof.connection(mysql_config) as conn:
        assert sector_paths.settle_sectors(conn, galaxy["sector_ids"], workers=1) > 0
        first = proof.settled_rows(conn, galaxy["sector_ids"])
        assert sector_paths.settle_sectors(conn, galaxy["sector_ids"], workers=1) > 0
        second = proof.settled_rows(conn, galaxy["sector_ids"])
    assert first
    proof.assert_same_rows(first, second, "settled twice")


def test_the_paths_a_settle_saves_do_not_depend_on_how_many_sectors_share_the_neighbourhood_read(mysql_config, monkeypatch):
    """GEN.1477: the neighbourhood read per sector (a batch of one, the old way) and shared by a batch of 100 saves the same bytes, and so does the stage's cell-key order."""
    galaxy = proof.small_galaxy(mysql_config)
    ids = galaxy["sector_ids"]
    rows = {}
    for batch in (1, 2, 100):
        monkeypatch.setattr(sector_paths, "SETTLE_BATCH", batch)
        with proof.connection(mysql_config) as conn:
            order = sector_paths.sectors_in_cell_order(conn, ids)
            assert sector_paths.settle_sectors(conn, order, expand=False, workers=1) > 0
            rows[batch] = proof.settled_rows(conn, ids)
    assert rows[1]
    proof.assert_same_rows(rows[1], rows[2], "a read for each sector against a read shared by two")
    proof.assert_same_rows(rows[1], rows[100], "a read for each sector against a read shared by the batch")


# --- GEN.1482: the link stage writes the same rows in any batch ---------------------------------------------------------------

def test_the_link_rows_of_the_small_galaxy_are_there_and_are_bytes(mysql_config):
    """The harness on the current stage: the rows the link stage wrote when the sectors were saved, in an order that does not depend on row ids. GEN.1482's test links the same
    sectors in batches of 2, 20 and 100 (and a galaxy grown by a second run) and compares each with `assert_same_rows`."""
    proof.small_galaxy(mysql_config)
    with proof.connection(mysql_config) as conn:
        rows = proof.link_rows(conn)
    assert rows["star_systems"] and all(isinstance(value, (bytes, int, str, type(None))) for row in rows["nearest_systems"] for value in row)
    assert rows["nearest_systems"], "the three neighbouring sectors list each other's systems as nearest"
