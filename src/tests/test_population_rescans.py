# tests/test_population_rescans.py

"""
TEST.38: the population pass run again and again as the galaxy grows
(`population.scan_life_worlds`, `refresh_civilizations`,
`refresh_territories` and the `population_state` watermark): a rescan
after new fills adds only the new worlds, leaves earlier species as they
were, never scans below the watermark, gives the same species in any
batch size, and `--rescan` and `--territories-only` do what they say.

Every test takes the `mysql_config` fixture (see `conftest.py`) and is
skipped, not failed, when no MySQL test server is reachable.
"""

import pytest

from stellarObjects import _db, population

from tests.test_population import _civilized_system, _place, _plain_system, _set_age, galaxy  # noqa: F401


def _species(conn):
    return {row["homeworld_planet_id"]: dict(row) for row in conn.execute(
        "SELECT id, name, homeworld_planet_id, star_system_id, life_stage, build, climate, size,"
        " civilization_age_years, era, spacefaring FROM species").fetchall()}


def _top_planet(conn):
    return conn.execute("SELECT MAX(id) AS id FROM planets").fetchone()["id"]


def _without_ids_and_names(species):
    return {planet: {key: value for key, value in row.items() if key not in ("id", "name")}
            for planet, row in species.items()}


@pytest.fixture
def conn(mysql_config):
    connection = _db.get_connection(mysql_config)
    try:
        yield connection
    finally:
        connection.close()


def test_a_rescan_after_new_fills_adds_only_the_new_worlds(mysql_config, conn):
    _civilized_system(mysql_config)
    first = population.run_pass(conn)
    before = _species(conn)
    assert first["new_species"] == len(before) >= 1
    assert population._watermark(conn) == _top_planet(conn)

    new_system = _civilized_system(mysql_config)
    # Plain systems have no intelligent life, but now and then still
    # a multicellular world, so this one may add a species too.
    plain = _plain_system(mysql_config)
    second = population.run_pass(conn)
    after = _species(conn)
    added = {planet: row for planet, row in after.items() if planet not in before}
    assert second["new_species"] == len(added) >= 1
    assert new_system in {row["star_system_id"] for row in added.values()}
    assert {row["star_system_id"] for row in added.values()} <= {new_system, plain}
    # Earlier species keep their ids, names and traits.
    assert {planet: after[planet] for planet in before} == before
    assert population._watermark(conn) == _top_planet(conn)

    assert population.run_pass(conn)["new_species"] == 0
    assert _species(conn) == after


def test_worlds_below_the_watermark_are_not_scanned_again(mysql_config, conn):
    _civilized_system(mysql_config)
    population.run_pass(conn)
    species = _species(conn)
    victim = min(species)
    with conn:
        conn.execute("DELETE FROM species WHERE homeworld_planet_id = ?", (victim,))
    # An ordinary pass starts at the watermark, so the gap stays...
    assert population.run_pass(conn)["new_species"] == 0
    assert victim not in _species(conn)
    # ...until a full rescan.
    counts = population.run_pass(conn, rescan=True)
    assert counts["new_species"] == len(species)
    assert victim in _species(conn)


def test_a_rescan_rebuilds_the_same_species_traits(mysql_config, conn):
    """Traits and ages come from `random.Random(planet_id)`, so a full
    rescan gives every world the same species again (names are drawn
    fresh)."""
    for _ in range(2):
        _civilized_system(mysql_config)
    population.run_pass(conn)
    before = _species(conn)
    counts = population.run_pass(conn, rescan=True)
    assert counts["new_species"] == len(before)
    assert _without_ids_and_names(_species(conn)) == _without_ids_and_names(before)


@pytest.mark.parametrize("batch", [1, 2, 7])
def test_any_scan_batch_size_finds_the_same_worlds(mysql_config, conn, monkeypatch, batch):
    for _ in range(2):
        _civilized_system(mysql_config)
    _plain_system(mysql_config)
    population.run_pass(conn)
    whole = _without_ids_and_names(_species(conn))
    monkeypatch.setattr(population, "SCAN_BATCH", batch)
    population.run_pass(conn, rescan=True)
    assert _without_ids_and_names(_species(conn)) == whole
    assert population._watermark(conn) == _top_planet(conn)


def test_an_empty_database_scans_nothing_and_sets_no_watermark(conn):
    counts = population.run_pass(conn)
    assert counts == {"new_species": 0, "species": 0, "spacefaring": 0, "polities": 0, "owned_systems": 0}
    assert population._watermark(conn) == 0


def test_territories_only_skips_new_worlds(mysql_config, conn):
    _civilized_system(mysql_config)
    population.run_pass(conn)
    mark = population._watermark(conn)
    count = len(_species(conn))
    _civilized_system(mysql_config)
    assert population.run_pass(conn, territories_only=True)["new_species"] == 0
    assert len(_species(conn)) == count
    assert population._watermark(conn) == mark
    assert population.run_pass(conn)["new_species"] >= 1


def test_refresh_civilizations_follows_changed_ages(mysql_config, conn, galaxy):  # noqa: F811
    population.run_pass(conn)
    with conn:
        _set_age(conn, galaxy["a"], 1e6)
        _set_age(conn, galaxy["b"], 500)    # below interstellar: no polity
    population.refresh_civilizations(conn)
    polities = conn.execute("SELECT id, capital_system_id, reach_ly FROM polities").fetchall()
    assert [row["capital_system_id"] for row in polities] == [galaxy["a"]]
    assert polities[0]["reach_ly"] == pytest.approx(population.reach_ly(1e6))
    # Older still: the same polity, a longer reach.
    with conn:
        conn.execute("UPDATE species SET civilization_age_years = 2e6 WHERE star_system_id = ?"
                     " AND civilization_age_years = 1e6", (galaxy["a"],))
    population.refresh_civilizations(conn)
    again = conn.execute("SELECT id, reach_ly FROM polities").fetchall()
    assert [row["id"] for row in again] == [polities[0]["id"]]
    assert again[0]["reach_ly"] == pytest.approx(population.reach_ly(2e6))
    # A second refresh with nothing changed founds nothing new.
    population.refresh_civilizations(conn)
    assert conn.execute("SELECT COUNT(*) AS n FROM polities").fetchone()["n"] == 1


def test_territories_grow_to_systems_filled_after_the_first_pass(mysql_config, conn, galaxy):  # noqa: F811
    population.run_pass(conn)
    with conn:
        _set_age(conn, galaxy["a"], 1e6)
    population.run_pass(conn)
    polity = conn.execute("SELECT id FROM polities WHERE capital_system_id = ?", (galaxy["a"],)).fetchone()["id"]
    owned_before = {row["star_system_id"] for row in conn.execute(
        "SELECT star_system_id FROM system_owners WHERE polity_id = ?", (polity,)).fetchall()}

    late = _plain_system(mysql_config)
    with conn:
        _place(conn, galaxy["sector"], late, -1.6)
    counts = population.run_pass(conn)
    owners = {row["star_system_id"]: row["polity_id"] for row in conn.execute(
        "SELECT star_system_id, polity_id FROM system_owners").fetchall()}
    assert owners[late] == polity
    assert owned_before | {late} <= {system for system, owner in owners.items() if owner == polity}
    assert counts["owned_systems"] == len(owners)
    # Only the late system's worlds can be new (a plain system now and
    # then still has multicellular life).
    on_late = conn.execute("SELECT COUNT(*) AS n FROM species WHERE star_system_id = ?", (late,)).fetchone()["n"]
    assert counts["new_species"] == on_late
