"""
GEN.46: no newly chosen star system name is longer than two words. A base
name is one word or two (`utils.split_long_word`); the uniqueness
decorations (`nameUniqueness.py`, applied in `_db.py`) may only be used
while the name stays within two words, and otherwise the name is drawn
again. Names already stored are left as they are (Boss, 2026-10-02).

The pure resolver is checked here without a database; the rest uses the
`mysql_config` fixture (skipped without a MySQL test server).
"""

import random

import pytest

from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.nameUniqueness import (
    GREEK_ROMAN_CAPACITY, MAX_SYSTEM_NAME_WORDS, resolve_greek_roman_collision, word_count,
)
from stellarObjects.names import DIMINUTIVE_PREFIXES, GREEK_LETTERS
from planetgen.galaxy.sector import SpaceSector
from stellarObjects.systemData import StarSystem


# ---------------------------------------------------------------------------
# The pure resolver
# ---------------------------------------------------------------------------

def test_a_one_word_base_gets_the_greek_tier_only():
    for count in range(1, len(GREEK_LETTERS)):
        new_name, rename = resolve_greek_roman_collision("Vor", count, max_words=MAX_SYSTEM_NAME_WORDS)
        assert new_name == f"{GREEK_LETTERS[count]} Vor"
        assert rename == (("Vor", "Alpha Vor") if count == 1 else None)
    for count in range(len(GREEK_LETTERS), GREEK_ROMAN_CAPACITY + 2):
        assert resolve_greek_roman_collision("Vor", count, max_words=MAX_SYSTEM_NAME_WORDS) == (None, None)


def test_a_two_word_base_gets_no_decoration():
    assert resolve_greek_roman_collision("Xy Zz", 0, max_words=2) == ("Xy Zz", None)
    for count in range(1, GREEK_ROMAN_CAPACITY + 2):
        assert resolve_greek_roman_collision("Xy Zz", count, max_words=2) == (None, None)


def test_without_a_limit_the_old_tiers_are_unchanged():
    """Sectors keep every tier (the rule is for star system names)."""
    assert resolve_greek_roman_collision("Xy Zz", 1) == ("Beta Xy Zz", ("Xy Zz", "Alpha Xy Zz"))
    assert resolve_greek_roman_collision("Vor", len(GREEK_LETTERS))[0] == "Alpha Vor II"


# ---------------------------------------------------------------------------
# Against the database
# ---------------------------------------------------------------------------

def _system(name):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M5V"
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    system.name = name
    return system


def _insert_system(conn, name):
    with conn:
        system = _system(name)
        _db.insert_star_system(conn, system, system.system_config)
    return system.name


def _all_names(conn, table):
    return [row["name"] for row in conn.execute(f"SELECT name FROM {table}").fetchall()]


@pytest.fixture
def fresh_names(monkeypatch):
    """Draws again from a known pool, one and two words."""
    pool = iter(f"Fresh{n}" if n % 2 else f"Fresh{n} Vale" for n in range(10_000))
    monkeypatch.setattr(_db, "_regenerate_star_name", lambda: next(pool))


def test_one_base_name_thirty_times(mysql_config, fresh_names):
    conn = _db.get_connection(mysql_config)
    try:
        names = [_insert_system(conn, "Vor") for _ in range(30)]
        stored = _all_names(conn, "star_systems")
    finally:
        conn.close()
    assert sorted(f"{letter} Vor" for letter in GREEK_LETTERS) == sorted(n for n in stored if n.endswith(" Vor"))
    assert all(n.startswith("Fresh") for n in names[len(GREEK_LETTERS):])
    assert all(word_count(n) <= 2 for n in stored)
    assert len({n.casefold() for n in stored}) == len(stored)


def test_a_two_word_base_colliding_draws_a_fresh_name_and_the_holder_keeps_its_own(mysql_config, fresh_names):
    conn = _db.get_connection(mysql_config)
    try:
        first = _insert_system(conn, "Xy Zz")
        second = _insert_system(conn, "Xy Zz")
        stored = _all_names(conn, "star_systems")
    finally:
        conn.close()
    assert first == "Xy Zz"
    assert second.startswith("Fresh")
    assert sorted(stored) == sorted([first, second])


def test_a_system_after_a_two_word_sector_draws_a_fresh_name(mysql_config, fresh_names):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Xy Zz"))
        name = _insert_system(conn, "Xy Zz")
    finally:
        conn.close()
    assert name.startswith("Fresh")


def test_a_sector_after_a_greek_decorated_system_takes_a_fresh_name(mysql_config, fresh_names):
    """Putting a diminutive on "Beta Vor" would make three words: the
    systems keep their names and the sector draws another."""
    conn = _db.get_connection(mysql_config)
    try:
        _insert_system(conn, "Vor")
        _insert_system(conn, "Vor")
        with conn:
            sector = SpaceSector(name="Vor")
            _db.insert_sector(conn, sector)
        systems = _all_names(conn, "star_systems")
    finally:
        conn.close()
    assert sorted(systems) == ["Alpha Vor", "Beta Vor"]
    assert sector.name != "Vor" and "Vor" not in sector.name.split(" ")


def test_a_sector_after_a_diminutive_system_leaves_it_alone(mysql_config, fresh_names):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Quessa"))
        assert _insert_system(conn, "Quessa") == f"{DIMINUTIVE_PREFIXES[0]} Quessa"
        for _ in range(3):
            with conn:
                _db.insert_sector(conn, SpaceSector(name="Quessa"))
        systems = _all_names(conn, "star_systems")
        sectors = _all_names(conn, "sectors")
    finally:
        conn.close()
    assert systems == [f"{DIMINUTIVE_PREFIXES[0]} Quessa"]
    assert sorted(sectors) == sorted(f"{letter} Quessa" for letter in GREEK_LETTERS[:4])


def test_a_hand_given_long_name_keeps_the_old_decorations(mysql_config, fresh_names):
    """The rule is for generated names: a longer name given by hand is
    saved as given and, on a collision, decorated as before rather than
    swapped for a random one."""
    conn = _db.get_connection(mysql_config)
    try:
        first = _insert_system(conn, "Lonely Shell Star")
        second = _insert_system(conn, "Lonely Shell Star")
        with conn:
            _db.insert_sector(conn, SpaceSector(name="Far Out Place"))
        third = _insert_system(conn, "Far Out Place")
        stored = _all_names(conn, "star_systems")
    finally:
        conn.close()
    assert first == "Lonely Shell Star"
    assert second == "Beta Lonely Shell Star"
    assert "Alpha Lonely Shell Star" in stored
    assert third == f"{DIMINUTIVE_PREFIXES[0]} Far Out Place"


def test_many_colliding_names_never_pass_two_words(mysql_config, fresh_names):
    """Sectors and systems drawn from a small pool, so they collide over
    and over at both levels: every system name stays within two words and
    every sector and system name is unique."""
    pool = ["Vor", "Quessa", "Tellow", "Xy Zz", "Ana Rel", "Mervane"]
    rng = random.Random(46)
    conn = _db.get_connection(mysql_config)
    try:
        for _ in range(120):
            name = rng.choice(pool)
            if rng.random() < 0.2:
                with conn:
                    _db.insert_sector(conn, SpaceSector(name=name))
            else:
                _insert_system(conn, name)
        systems = _all_names(conn, "star_systems")
        sectors = _all_names(conn, "sectors")
    finally:
        conn.close()
    too_long = [name for name in systems if word_count(name) > MAX_SYSTEM_NAME_WORDS]
    assert not too_long
    folded = [name.casefold() for name in systems + sectors]
    assert len(folded) == len(set(folded))
