"""PERF.284 and OPS.289: the `remnants` piece of the special phenomena unit builds the address index its block light reads need, and says where it is."""

from types import SimpleNamespace

from planetgen.db import store
from planetgen.generation import run_plan


def test_time_left_reads_the_pace_the_loop_kept():
    assert run_plan._time_left_text(0, 100, 5.0) is None
    assert run_plan._time_left_text(100, 100, 5.0) is None
    assert run_plan._time_left_text(25, 100, 60.0) == "about 3m left"
    assert run_plan._time_left_text(1, 3601, 1.0) == "about 1h left"


def test_the_remnants_piece_builds_the_address_index_when_the_scatter_took_it_off(mysql_config, monkeypatch):
    conn = store.get_connection(mysql_config)
    try:
        store.drop_phenomenon_scatter_index(conn, "idx_phenomenon_scatter_address")
        assert "idx_phenomenon_scatter_address" not in store.phenomenon_scatter_index_names(conn)
        said = []
        steps = []
        ctx = SimpleNamespace(db=conn, kept=lambda name, make: make(), step=lambda label, total: steps.append((label, total)),
                              phase=lambda text, **kw: said.append(text), commit_piece=lambda result: None, advance=lambda n: None)
        monkeypatch.setattr(store, "get_galaxy_shape", lambda c: SimpleNamespace())
        monkeypatch.setattr(store, "scatter_limits_of", lambda c, skeleton: None)
        monkeypatch.setattr(store, "backfilled_uncharted_addresses", lambda c: [(1, 1, 1), (1, 1, 2)])
        monkeypatch.setattr(store, "addresses_holding_below_cut_remnants", lambda c, mass: set())
        monkeypatch.setattr(run_plan.ringRemnants, "rows_for", lambda *a, **k: [])
        run_plan._write_special_remnants({"min_mass_solar": 8.0, "seed": 1}, ctx)
        assert "idx_phenomenon_scatter_address" in store.phenomenon_scatter_index_names(conn)
        assert any("address index" in text for text in said)
        assert steps == [("Neutron stars and black holes of the backfill rings", 2)]
    finally:
        conn.close()
