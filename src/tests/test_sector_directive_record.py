# tests/test_sector_directive_record.py

"""GEN.179: a sector keeps the directive it was drawn under and how the draw came out. `directives.record_of` makes the record from the directive,
the outcome and the sector kept; `run_sector.generate_directed_sector` puts it on the sector; `store.save_sector` writes it (`sector_directives`,
`sector_directive_minimums`) with the sector; the sector's API detail carries it as `directive`, its page shows it in the "Drawn under" row, and
`check-db counts` finds a record that contradicts itself."""

import argparse

import pytest

from planetgen.db import check, query, store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import directives, run_sector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.web import glance


class _Star:
    def __init__(self, type_, yerkes="V"):
        self.type, self.yerkes_class = type_, yerkes


class _System:
    def __init__(self, type_, hab=0):
        self.star, self.habitable_class_count = _Star(type_), hab


class _Entry:
    def __init__(self, type_, hab=0):
        self.star_system = _System(type_, hab)


class _Sector:
    def __init__(self, entries):
        self.entries = entries
        self.name = "Fake"
        self.directive_record = None


def _sector():
    return _Sector([_Entry("G2", hab=1), _Entry("K5", hab=2), _Entry("G0")])


def test_the_record_lists_each_minimum_with_what_was_wanted_and_what_the_sector_got():
    directive = directives.parse(["habitable>=2", "systems>=5", "type:G>=1"])
    record = directives.record_of(directive, directives.Outcome(directives.UNMET, 8, directives.shortfall(directive, _sector())), _sector())
    assert record.digest == directive.digest() and len(record.digest) == 16
    assert record.status == "unmet" and record.attempts == 8
    assert record.minimums == [("habitable", 2, 3), ("systems", 5, 3), ("type:G", 1, 2)]


def test_a_directive_that_asks_for_nothing_has_no_record():
    assert directives.record_of(directives.Directive(), directives.Outcome(directives.MET_NATURALLY, 1, {}), _sector()) is None
    assert directives.record_of(directives.parse(["systems>=0"]), directives.Outcome(directives.MET_NATURALLY, 1, {}), _sector()) is None


def test_every_outcome_status_is_one_the_table_accepts():
    assert set(directives.STATUSES) == {"met_naturally", "met_after_k_attempts", "met_forced", "unmet"}


def test_a_directed_sector_carries_its_record_and_a_plain_one_has_none(monkeypatch):
    monkeypatch.setattr(run_sector, "generate_sector", lambda args, **kwargs: ("Fake", _sector()))
    plain = argparse.Namespace(directive=None)
    assert run_sector.generate_directed_sector(plain)[1].directive_record is None
    asked = argparse.Namespace(directive=["systems>=2"], directive_attempts=3)
    record = run_sector.generate_directed_sector(asked)[1].directive_record
    assert record.status == "met_naturally" and record.attempts == 1 and record.minimums == [("systems", 2, 3)]


@pytest.mark.parametrize("record, text", [
    (None, None),
    ({"status": "met_naturally", "attempts": 1, "minimums": [{"key": "habitable", "wanted": 2, "got": 3}]},
     "at least 2 habitable worlds (got 3): met on the first draw"),
    ({"status": "met_after_k_attempts", "attempts": 4, "minimums": [{"key": "type:G", "wanted": 1, "got": 2}]},
     "at least 1 G stars (got 2): met after 4 draws"),
    ({"status": "unmet", "attempts": 8, "minimums": [{"key": "systems", "wanted": 5, "got": 3}]},
     "at least 5 star systems (got 3): not met after 8 draws; the closest draw was kept"),
])
def test_the_sector_page_says_what_was_asked_and_how_the_draw_came_out(record, text):
    row = glance._directive_row(record)
    assert (row is None) if text is None else (row["label"] == "Drawn under" and row["html"] == text)


def _saved_sector(mysql_config, record):
    sector = SpaceSector("Directive Sector", edge_ly=20.0)
    config = SystemConfig()
    config.PLANETS = False
    config.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=config), position=(1.0, 1.0, 1.0), system_config=config)
    sector.directive_record = record
    return store.save_sector(sector, config=mysql_config)


def _record():
    return directives.Record("0123456789abcdef", directives.MET_AFTER_ATTEMPTS, 3, [("habitable", 2, 2), ("systems", 1, 1)])


def test_a_saved_sector_keeps_its_record_and_the_api_detail_gives_it(mysql_config):
    sector_id = _saved_sector(mysql_config, _record())
    conn = store.get_connection(mysql_config)
    try:
        assert query.sector_directive(conn, sector_id) == {
            "digest": "0123456789abcdef", "status": "met_after_k_attempts", "attempts": 3,
            "minimums": [{"key": "habitable", "wanted": 2, "got": 2}, {"key": "systems", "wanted": 1, "got": 1}]}
        assert query.sector_detail(conn, sector_id)["directive"]["status"] == "met_after_k_attempts"
    finally:
        conn.close()


def test_a_sector_drawn_with_no_directive_stores_nothing(mysql_config):
    sector_id = _saved_sector(mysql_config, None)
    conn = store.get_connection(mysql_config)
    try:
        assert query.sector_directive(conn, sector_id) is None
        assert query.sector_detail(conn, sector_id)["directive"] is None
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_directives").fetchone()["n"] == 0
    finally:
        conn.close()


def test_the_record_goes_with_the_sector(mysql_config):
    sector_id = _saved_sector(mysql_config, _record())
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM star_systems WHERE sector_id = ?", (sector_id,))
        conn.execute("DELETE FROM sectors WHERE id = ?", (sector_id,))
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_directives").fetchone()["n"] == 0
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_directive_minimums").fetchone()["n"] == 0
    finally:
        conn.close()


def test_the_counts_check_finds_a_record_that_contradicts_itself(mysql_config):
    sector_id = _saved_sector(mysql_config, _record())
    conn = store.get_connection(mysql_config)
    try:
        ctx = check._Context(conn, mysql_config, check.Scope())
        assert check._counts_piece(ctx, "directives") == []
        conn.execute("UPDATE sector_directives SET status = 'unmet' WHERE sector_id = ?", (sector_id,))
        found = check._counts_piece(ctx, "directives")
        assert len(found) == 1 and "contradicts itself" in found[0]["what"] and not found[0]["warning"]
        conn.execute("UPDATE sector_directives SET status = 'met_forced' WHERE sector_id = ?", (sector_id,))
        conn.execute("UPDATE sector_directive_minimums SET got = 0 WHERE minimum_key = 'habitable'")
        assert len(check._counts_piece(ctx, "directives")) == 1
    finally:
        conn.close()
