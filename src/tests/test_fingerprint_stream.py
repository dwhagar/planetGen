"""PERF.286: the fingerprint keeps a sector's rows at a time, not a galaxy's."""

from planetgen.db import fingerprint


def test_the_reader_lets_go_of_the_row_keys_of_a_sector_but_keeps_the_sector_and_config_keys():
    reader = fingerprint._Reader(None)
    reader._keys = {("planets", 1): "planets:a", ("stars", 2): "stars:b", ("sectors", 3): "sector:1 2 3", ("system_configs", 4): "config:abc"}
    reader.forget_rows()
    assert reader._keys == {("sectors", 3): "sector:1 2 3", ("system_configs", 4): "config:abc"}
