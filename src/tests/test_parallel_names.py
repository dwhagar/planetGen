# tests/test_parallel_names.py

"""
TEST.37: names under parallel saves. Several writers (each its own
connection, as each generation worker has) saving sectors and systems
with the same base name at the same moment: every saved name stays
unique, same-level collisions get Greek decorations, the diminutive tier
fills without reusing a word, two population passes at once don't fight
over species names, and `_unique_species_name` falls back to numbered
names and finally gives up with "could not find a free species name".

Every database test takes the `mysql_config` fixture (see `conftest.py`)
and is skipped, not failed, when no MySQL test server is reachable.
"""

import threading

import pytest

from stellarObjects import _db, population
from stellarObjects.config import SystemConfig
from stellarObjects.names import DIMINUTIVE_PREFIXES, GREEK_LETTERS
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem

from tests.test_population import _civilized_system

WRITERS = 4


@pytest.fixture
def mysql_config(mysql_config):
    """The fresh database with its schema already applied, as a real run's
    database has before its workers start: several first connections to
    an empty database at once would all apply the schema together."""
    _db.get_connection(mysql_config).close()
    return mysql_config


def _at_once(count, work):
    """Runs `work(index)` on `count` threads released together; returns
    their results in index order and re-raises the first error."""
    barrier = threading.Barrier(count)
    results = [None] * count
    errors = []

    def run(index):
        try:
            barrier.wait(timeout=30)
            results[index] = work(index)
        except BaseException as exc:  # noqa: BLE001 -- re-raised below
            errors.append(exc)

    threads = [threading.Thread(target=run, args=(index,)) for index in range(count)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=300)
    if errors:
        raise errors[0]
    return results


def _system(name):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M5V"
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    system.name = name
    return system


def _save_in_sector(config, index, *system_names):
    """What a generation worker does: saves one sector holding systems
    with these names (`save_sector`, retried on a deadlock)."""
    built = SpaceSector(name=f"Holder{index}x")
    for offset, name in enumerate(system_names):
        system = _system(name)
        built.add_system(system, position=(offset * 2.0, 0.0, 0.0), system_config=system.system_config)
    return _db.save_sector(built, config=config)


def _names(config, table):
    conn = _db.get_connection(config)
    try:
        return [row["name"] for row in conn.execute(f"SELECT name FROM {table}").fetchall()]
    finally:
        conn.close()


def _all_unique_names(config):
    names = [name.casefold() for name in _names(config, "sectors") + _names(config, "star_systems")]
    assert len(names) == len(set(names)), sorted(names)


# ---------------------------------------------------------------------------
# Same level: sectors with sectors, systems with systems
# ---------------------------------------------------------------------------

def test_sectors_saved_at_once_with_one_name_get_greek_letters(mysql_config):
    ids = _at_once(WRITERS, lambda _i: _db.save_sector(SpaceSector(name="Corvane"), config=mysql_config))
    assert len(set(ids)) == WRITERS
    names = _names(mysql_config, "sectors")
    assert sorted(names) == sorted(f"{letter} Corvane" for letter in GREEK_LETTERS[:WRITERS])
    conn = _db.get_connection(mysql_config)
    try:
        registry = conn.execute(
            "SELECT occurrence_count, first_sector_id FROM sector_name_registry WHERE base_name = 'Corvane'"
        ).fetchone()
    finally:
        conn.close()
    assert registry["occurrence_count"] == WRITERS
    assert registry["first_sector_id"] in ids


def test_systems_saved_at_once_with_one_name_get_greek_letters(mysql_config):
    ids = _at_once(WRITERS, lambda i: _save_in_sector(mysql_config, i, "Halveth"))
    assert len(set(ids)) == WRITERS
    names = _names(mysql_config, "star_systems")
    assert sorted(names) == sorted(f"{letter} Halveth" for letter in GREEK_LETTERS[:WRITERS])


def test_sectors_whose_systems_share_names_save_at_once(mysql_config):
    """Each worker's sector holds systems named like every other worker's
    (and twice within itself): the batched system-name reservation still
    hands out unique names."""
    sectors = []
    for index in range(WRITERS):
        built = SpaceSector(name=f"Reach{index}")
        for offset, base in enumerate(("Ostra", "Ostra", "Pellin")):
            system = _system(base)
            built.add_system(system, position=(offset * 2.0, 0.0, 0.0), system_config=system.system_config)
        sectors.append(built)
    _at_once(WRITERS, lambda i: _db.save_sector(sectors[i], config=mysql_config))
    names = _names(mysql_config, "star_systems")
    assert len(names) == 3 * WRITERS
    assert sum(name.endswith("Ostra") for name in names) == 2 * WRITERS
    _all_unique_names(mysql_config)


# ---------------------------------------------------------------------------
# Across levels: the diminutive tier
# ---------------------------------------------------------------------------

def test_systems_saved_at_once_after_a_sector_fill_the_diminutive_tier(mysql_config):
    _db.save_sector(SpaceSector(name="Mervane"), config=mysql_config)
    count = 6
    _at_once(count, lambda i: _save_in_sector(mysql_config, i, "Mervane"))
    names = _names(mysql_config, "star_systems")
    assert "Mervane" not in names
    # The first holder's "Little" gave way to "Alpha" when the second
    # arrived (a Greek letter replaces a diminutive); every later system
    # carries its own diminutive word, none reused, in the list's order.
    assert "Alpha Mervane" in names
    decorated = [name.split(" ") for name in names if name != "Alpha Mervane"]
    assert sorted(words[0] for words in decorated) == sorted(DIMINUTIVE_PREFIXES[1:count])
    assert sorted(words[1] for words in decorated) == sorted(GREEK_LETTERS[1:count])
    assert "Mervane" in _names(mysql_config, "sectors")
    _all_unique_names(mysql_config)


def test_sectors_saved_at_once_after_a_system_rename_only_the_system(mysql_config):
    _save_in_sector(mysql_config, 0, "Quessa")
    _at_once(WRITERS, lambda _i: _db.save_sector(SpaceSector(name="Quessa"), config=mysql_config))
    sectors = [name for name in _names(mysql_config, "sectors") if name.endswith("Quessa")]
    assert sorted(sectors) == sorted(f"{letter} Quessa" for letter in GREEK_LETTERS[:WRITERS])
    [system_name] = _names(mysql_config, "star_systems")
    # One sector base name: the system steps down one diminutive per
    # sector holder that collided with it.
    assert system_name.endswith("Quessa") and system_name != "Quessa"
    _all_unique_names(mysql_config)


def test_a_full_diminutive_tier_draws_a_fresh_name(mysql_config, monkeypatch):
    _db.save_sector(SpaceSector(name="Tellow"), config=mysql_config)
    count = len(DIMINUTIVE_PREFIXES) + 2
    fresh = iter(f"Fresh{n}" for n in range(1000))
    monkeypatch.setattr(_db, "_regenerate_star_name", lambda: next(fresh))
    for start in range(0, count, WRITERS):
        size = min(WRITERS, count - start)
        _at_once(size, lambda i, start=start: _save_in_sector(mysql_config, start + i, "Tellow"))
    names = _names(mysql_config, "star_systems")
    assert len(names) == count
    assert sum(name.endswith("Tellow") for name in names) == len(DIMINUTIVE_PREFIXES)
    assert sum(name.startswith("Fresh") for name in names) == 2
    _all_unique_names(mysql_config)


# ---------------------------------------------------------------------------
# Species names
# ---------------------------------------------------------------------------

def test_two_population_passes_at_once_name_each_world_once(mysql_config):
    for _ in range(3):
        _civilized_system(mysql_config)

    def run(_index):
        conn = _db.get_connection(mysql_config)
        try:
            return population.run_pass(conn)
        finally:
            conn.close()

    first, second = _at_once(2, run)
    conn = _db.get_connection(mysql_config)
    try:
        species = conn.execute("SELECT name, homeworld_planet_id FROM species").fetchall()
        watermark = population._watermark(conn)
        top = conn.execute("SELECT MAX(id) AS id FROM planets").fetchone()["id"]
    finally:
        conn.close()
    assert species
    assert first["new_species"] + second["new_species"] == len(species)
    assert sorted((first["new_species"], second["new_species"]))[0] == 0
    assert len({row["homeworld_planet_id"] for row in species}) == len(species)
    assert len({row["name"] for row in species}) == len(species)
    assert watermark == top


def test_two_passes_drawing_the_same_names_still_save_unique_ones(mysql_config, monkeypatch):
    """Both passes' name generators give the same sequence: the second
    pass, run after the first, skips every taken name."""
    for _ in range(2):
        _civilized_system(mysql_config)
    local = threading.local()

    def same_sequence():
        local.n = getattr(local, "n", 0) + 1
        return f"Kevra{local.n}"

    monkeypatch.setattr(population, "new_species_name", same_sequence)
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("DELETE FROM species")
        _at_once(2, lambda _i: _pass_on_own_connection(mysql_config, rescan=True))
        names = [row["name"] for row in conn.execute("SELECT name FROM species").fetchall()]
    finally:
        conn.close()
    assert names and len(names) == len(set(names))


def _pass_on_own_connection(config, **kwargs):
    conn = _db.get_connection(config)
    try:
        return population.run_pass(conn, **kwargs)
    finally:
        conn.close()


class _AllTaken:
    """A connection on which every species name is already stored."""

    def __init__(self, free=None):
        self.free = free
        self.asked = []

    def execute(self, sql, params=()):
        self.asked.append(params[0])
        return _Row(None if params[0] == self.free else {"1": 1})


class _Row:
    def __init__(self, row):
        self.row = row

    def fetchone(self):
        return self.row


def test_a_species_name_taken_every_draw_falls_back_to_a_number(monkeypatch):
    monkeypatch.setattr(population, "new_species_name", lambda: "Vorn")
    conn = _AllTaken(free="Vorn 4")
    assert population._unique_species_name(conn, taken=set()) == "Vorn 4"
    # Drawn names already in `taken` never reach the database.
    conn = _AllTaken(free="Vorn 2")
    assert population._unique_species_name(conn, taken={"Vorn"}) == "Vorn 2"
    assert conn.asked == ["Vorn 2"]


def test_no_free_species_name_at_all_is_an_error(monkeypatch):
    monkeypatch.setattr(population, "new_species_name", lambda: "Vorn")
    with pytest.raises(RuntimeError, match="could not find a free species name"):
        population._unique_species_name(_AllTaken(), taken=set())
