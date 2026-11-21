"""PERF.283: the forced-star tables are kept across processes in the sampling-table cache and read back whole."""

from planetgen.generation import forced_star


def test_a_table_is_kept_and_read_back_the_same(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_BAND_TABLE_CACHE", str(tmp_path))
    monkeypatch.setattr(forced_star.tuning, "FORCED_STAR_MASS_STEPS", 6)
    forced_star._table.cache_clear()
    built = forced_star._table(0.0)
    assert built and list(tmp_path.rglob("forced-star-*.json"))
    forced_star._table.cache_clear()
    monkeypatch.setattr(forced_star, "_build_table", lambda bucket: (_ for _ in ()).throw(AssertionError("rebuilt")))
    assert forced_star._table(0.0) == built
    forced_star._table.cache_clear()
