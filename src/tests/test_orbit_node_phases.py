# tests/test_orbit_node_phases.py

"""
GEN.313: the orbital update advances the nodes of a built system's orbit tree (`orbit_nodes`). A due node's mean anomaly turns by
360 times the elapsed time over its own period, its epoch becomes the update's end, a second run over the same span moves
nothing, a scatter row's nodes (`scatter_nodes`) are never written, and the update writes no `binary_mutual_*` column of a pair (GEN.640 slice a: the pair is its node).
"""

import pytest

from planetgen.db import orbit_trees, store
from tests.member_systems import binary, with_member

pytestmark = pytest.mark.fresh_database

SECONDS_PER_YEAR = 365.25 * 86400.0


def _saved(mysql_config, wide):
    system, cfg = binary(wide)
    with_member(system, cfg)
    conn = store.get_connection(mysql_config)
    system_id = store.insert_star_system(conn, system, cfg)
    conn.commit()
    return conn, system_id


@pytest.mark.parametrize("wide", [False, True])
def test_every_node_of_a_built_system_turns_by_its_own_period_and_takes_the_update_end_as_its_epoch(mysql_config, wide):
    conn, system_id = _saved(mysql_config, wide)
    try:
        before = orbit_trees.read_nodes(conn, system_id)
        assert len(before) >= 2
        clock = store.orbit_clock(conn, 3.25 * max(node.period_years for node in before))     # long enough that every node is due
        counts = store.advance_orbital_phases(conn, clock)
        assert counts["orbit_nodes"] == len(before)
        after = orbit_trees.read_nodes(conn, system_id)
        for old, new in zip(before, after):
            turned = 360.0 * (clock.end_unix - old.epoch_unix) / (old.period_years * SECONDS_PER_YEAR)
            assert new.mean_anomaly_deg == pytest.approx((old.mean_anomaly_deg + turned) % 360.0, abs=1e-4)
            assert new.epoch_unix == pytest.approx(clock.end_unix)
            assert (new.separation_km, new.eccentricity, new.period_years, new.mass_solar) == (
                old.separation_km, old.eccentricity, old.period_years, old.mass_solar)
        assert store.advance_orbital_phases(conn, clock)["orbit_nodes"] == 0
        assert [n.mean_anomaly_deg for n in orbit_trees.read_nodes(conn, system_id)] == [n.mean_anomaly_deg for n in after]
    finally:
        conn.close()


def test_a_scatter_row_s_nodes_are_never_written_and_the_binary_keeps_its_own_phase(mysql_config):
    conn, system_id = _saved(mysql_config, False)
    try:
        conn.execute(
            "INSERT INTO scatter_nodes (source, source_id, node_no, mass_solar, separation_km, eccentricity, inclination_deg,"
            " ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, epoch_unix, period_years, kind)"
            " VALUES (1, 1, 0, 1.0, 1.0e6, 0.0, 0.0, 0.0, 0.0, 12.5, 100.0, 1.0, 'close')")
        conn.commit()
        columns = "binary_mutual_orbital_phase_deg AS p, binary_mutual_position_x_km AS x, binary_epoch_unix AS e"
        pair = orbit_trees.pair_places_km(conn, system_id)["node"]
        binary_before = dict(conn.execute(f"SELECT {columns} FROM star_systems WHERE id = ?", (system_id,)).fetchone())
        clock = store.orbit_clock(conn, pair.period_years * 3.25)
        store.advance_orbital_phases(conn, clock)
        row = conn.execute("SELECT mean_anomaly_deg AS m, epoch_unix AS e FROM scatter_nodes WHERE source_id = 1").fetchone()
        assert (row["m"], row["e"]) == (12.5, 100.0)
        # GEN.640 slice a: the pair's own columns are not written; its node moved by the turn from its epoch (the instant it was saved) to the update's end
        assert dict(conn.execute(f"SELECT {columns} FROM star_systems WHERE id = ?", (system_id,)).fetchone()) == binary_before
        moved = orbit_trees.pair_places_km(conn, system_id)["node"]
        turns = (clock.end_unix - pair.epoch_unix) / (pair.period_years * SECONDS_PER_YEAR)
        assert moved.mean_anomaly_deg == pytest.approx((pair.mean_anomaly_deg + turns * 360.0) % 360.0, abs=1e-4)
    finally:
        conn.close()
