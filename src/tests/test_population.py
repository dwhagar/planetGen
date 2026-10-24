# tests/test_population.py

"""Population and politics (schema v44, POP.1 to POP.4): species, civilization
ages and eras, polities and territories, the `population` CLI and the API.
See docs/design/habitability/population-and-politics.md."""

import random
import sys

import pytest

from planetgen.web.app import create_app
from planetgen.api.config import Config
from planetgen.db import store
from planetgen.population import model as population
from planetgen.population import tech
from planetgen.util import draw
from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics.units import ly_to_pc

from planetgen.cli import generate as generate_cli
from planetgen.generation import run_population
from tests.publicids import pid


# ---------------------------------------------------------------------------
# Pure model
# ---------------------------------------------------------------------------

def test_timeline_of_reads_stage_and_ages():
    timeline = population.timeline_of("technological_civilization", 4.6, 4.5)
    assert timeline.life_stage == "technological_civilization"
    assert timeline.system_age_years == pytest.approx(4.6e9)
    assert timeline.milestone_age_years == pytest.approx(4.5e9)
    assert timeline.window_years == pytest.approx(1e8)

    multicellular = population.timeline_of("multicellularity", 2.0, 0.9)
    assert multicellular.life_stage == "multicellularity"
    assert multicellular.milestone_age_years == pytest.approx(9e8)


@pytest.mark.parametrize("milestone", ["abiogenesis", "photosynthesis", "complex_cells", None])
def test_simple_life_is_not_a_life_world(milestone):
    assert population.timeline_of(milestone, 2.0, 1.0) is None


def test_a_timeline_needs_both_ages():
    assert population.timeline_of("technological_civilization", None, 4.5) is None
    assert population.timeline_of("technological_civilization", 4.6, None) is None


def test_timeline_of_reads_real_generated_planets():
    """A generated planet's stored structure reads back as a timeline."""
    for _ in range(200):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.INTELLIGENT_LIFE = True
        cfg.HABITABLE_WORLD = True
        system = StarSystem(system_config=cfg)
        for planet in system.planets:
            timeline = population.timeline_of(getattr(planet, "life_stage", None),
                                              getattr(planet, "life_system_age_gyr", None),
                                              getattr(planet, "life_milestone_age_gyr", None))
            if timeline is not None and timeline.life_stage == "technological_civilization":
                assert timeline.milestone_age_years <= timeline.system_age_years
                return
    pytest.fail("no generated timeline reached a civilization")


def test_civilizations_are_rare_unless_forced():
    tech = population.Timeline("technological_civilization", 4.6e9, 4.5e9)
    multicellular = population.Timeline("multicellularity", 4.6e9, 1e9)
    rng = random.Random(3)
    assert population.has_civilization(tech, True, rng)
    assert not population.has_civilization(multicellular, True, rng)
    draws = sum(population.has_civilization(tech, None, rng) for _ in range(20000))
    assert draws == pytest.approx(20000 * tuning.CIVILIZATION_CHANCE, abs=20)


def test_civilization_age_stays_inside_the_window():
    rng = random.Random(1)
    low = tuning.CIVILIZATION_MIN_AGE_YEARS
    for window in (1e3, 1e6, 1e9):
        ages = [population.civilization_age(window, rng) for _ in range(500)]
        assert min(ages) >= low and max(ages) <= window
    assert population.civilization_age(0.0, rng) == low
    assert population.civilization_age(50.0, rng) == low


@pytest.mark.parametrize("age, era, spacefaring", [
    (100, "Industrial", False),
    (300, "Interplanetary", False),
    (1999, "Interplanetary", False),
    (2000, "Interstellar", True),
    (60_000, "Established", True),
    (2e6, "Ancient", True),
    (5e8, "Elder", True),
])
def test_eras(age, era, spacefaring):
    assert population.era_for_age(age) == (era, spacefaring)


def test_reach_grows_with_age_up_to_the_cap():
    assert population.reach_ly(1000) == 0.0
    assert population.reach_ly(2000) == pytest.approx(tuning.TERRITORY_BASE_REACH_LY)
    assert population.reach_ly(8000) == pytest.approx(2 * tuning.TERRITORY_BASE_REACH_LY)
    assert population.reach_ly(1e9) == tuning.TERRITORY_REACH_CAP_LY
    assert population.reach_ly(1e9, cap_ly=500) == 500


def test_species_traits_follow_the_homeworld():
    rng = random.Random(2)
    assert population.species_traits(0.5, 200, rng)["build"] == "gracile"
    assert population.species_traits(0.5, 200, rng)["climate"] == "cold"
    assert population.species_traits(2.0, 350, rng)["build"] == "robust"
    assert population.species_traits(2.0, 350, rng)["climate"] == "hot"
    traits = population.species_traits(None, None, rng)
    assert traits["build"] == "medium" and traits["climate"] == "temperate"
    assert traits["size"] in ("small", "medium", "large")


def test_resolve_claims_strongest_claim_wins():
    ly = ly_to_pc(1.0)
    polities = [(1, (0.0, 0.0, 0.0), 10.0), (2, (12 * ly, 0.0, 0.0), 50.0)]
    systems = {
        10: (0.0, 0.0, 0.0),        # polity 1's capital
        11: (3 * ly, 0.0, 0.0),     # 1: 10/3, 2: 50/9 -> 2 wins
        12: (-2 * ly, 0.0, 0.0),    # 1: 10/2, 2: 50/14 -> 1 wins
        13: (-9 * ly, 0.0, 0.0),    # only polity 1 reaches (2 is 21 ly away, inside 50)... 50/21 > 10/9
        14: (200 * ly, 0.0, 0.0),   # nobody
    }
    owners = population.resolve_claims(polities, systems)
    assert owners[10][0] == 1 and owners[10][1] == 0.0
    assert owners[11][0] == 2
    assert owners[12][0] == 1
    assert owners[13][0] == 2
    assert 14 not in owners


def test_polity_color_and_government_are_stable():
    uid = bytes.fromhex("000010000000004000000002")
    assert population.polity_color(uid) == population.polity_color(uid)
    assert population.polity_color(uid).startswith("#") and len(population.polity_color(uid)) == 7
    assert population.government_for(uid) in tuning.GOVERNMENT_FORMS
    assert population.government_for(uid) == population.government_for(uid)
    assert population.polity_color(uid) != population.polity_color(bytes.fromhex("000010000000004000000003"))
    with pytest.raises(ValueError, match="no unique ID"):
        population.government_for(None)


# ---------------------------------------------------------------------------
# Database pass
# ---------------------------------------------------------------------------

def _gets_a_species(conn, system_id):
    """Whether the population pass will put a species on a planet of the saved system: a world that reached a
    technological civilization (`has_civilization`; intelligent life is forced on in `_civilized_system`) and whose
    limits can hold its species (POP.14: no species where undrinkable water rules it out)."""
    from planetgen.db import species_fit
    for row in conn.execute("SELECT id, life_stage, life_system_age_gyr, life_milestone_age_gyr FROM planets "
                            "WHERE star_system_id = ?", (system_id,)).fetchall():
        timeline = population.timeline_of(row["life_stage"], row["life_system_age_gyr"], row["life_milestone_age_gyr"])
        if timeline is None or timeline.life_stage != "technological_civilization":
            continue
        try:
            species_fit.check_homeworld(conn, row["id"])
        except species_fit.PlacementError:
            continue
        return True
    return False


def _civilized_system(mysql_config):
    """A saved system with a habitable planet that reached a technological civilization and will get a species --
    retried, since generation is random. A system the pass would leave without a species is deleted again (TEST.296:
    one in a few kept by the old test of the timeline alone, and a test counting the species failed)."""
    for _ in range(200):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.INTELLIGENT_LIFE = True
        cfg.HABITABLE_WORLD = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        if any(getattr(p, "life_stage", None) == "technological_civilization" for p in system.planets):
            system_id = store.save_system(system, cfg, config=mysql_config)
            conn = store.get_connection(mysql_config)
            try:
                if _gets_a_species(conn, system_id):
                    return system_id
                with conn:
                    conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
            finally:
                conn.close()
    pytest.fail("could not generate a system with a civilization")


def _plain_system(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M5V"
    cfg.INTELLIGENT_LIFE = False
    cfg.BINARY_SYSTEM = False
    return store.save_system(StarSystem(system_config=cfg), cfg, config=mysql_config)


def _place(conn, sector_id, system_id, offset_pc):
    conn.execute(
        "UPDATE star_systems SET sector_id = ?, position_x_mpc = ?, position_y_mpc = 0, position_z_mpc = 0 "
        "WHERE id = ?", (sector_id, offset_pc * 1000.0, system_id),
    )


def _set_age(conn, system_id, age):
    """Every civilization industrial, then one in `system_id` `age` years
    old -- so each capital system has exactly one candidate polity."""
    conn.execute("UPDATE species SET civilization_age_years = 100 "
                 "WHERE life_stage = 'technological_civilization' AND star_system_id = ?", (system_id,))
    conn.execute("UPDATE species SET civilization_age_years = ? WHERE star_system_id = ? "
                 "AND life_stage = 'technological_civilization' ORDER BY id LIMIT 1", (age, system_id))


@pytest.fixture
def galaxy(mysql_config):
    """Two civilized systems and three plain ones placed in one sector."""
    capital_a = _civilized_system(mysql_config)
    capital_b = _civilized_system(mysql_config)
    plain = [_plain_system(mysql_config) for _ in range(3)]
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            sector_id = conn.execute(
                "INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc, "
                "shell_index, band_index, wedge_index) VALUES ('Test Reach', 4000, 8000, 0, 0, 8000, 1, 0, 0)"
            ).lastrowid
            _place(conn, sector_id, capital_a, -1.5)
            _place(conn, sector_id, capital_b, 1.5)
            _place(conn, sector_id, plain[0], -1.4)   # next to A
            _place(conn, sector_id, plain[1], 1.4)    # next to B
            _place(conn, sector_id, plain[2], 0.0)    # in between
    finally:
        conn.close()
    return {"a": capital_a, "b": capital_b, "plain": plain, "sector": sector_id}


def test_pass_names_species_founds_polities_and_draws_territories(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        counts = population.run_pass(conn)
        assert counts["new_species"] >= 2
        species = conn.execute("SELECT * FROM species ORDER BY id").fetchall()
        names = [row["name"] for row in species]
        assert len(set(names)) == len(names)
        civilized = {row["star_system_id"] for row in species if row["life_stage"] == "technological_civilization"}
        assert {galaxy["a"], galaxy["b"]} <= civilized
        # GEN.80: only a technological civilization gets a species.
        assert all(row["life_stage"] == "technological_civilization" for row in species)
        assert all(row["civilization_age_years"] is not None for row in species)
        for row in species:
            if row["civilization_age_years"] is not None:
                assert row["era"] == population.era_for_age(row["civilization_age_years"])[0]

        # A second pass finds nothing new.
        assert population.run_pass(conn)["new_species"] == 0

        # Make both capitals' civilizations spacefaring, A much older.
        with conn:
            _set_age(conn, galaxy["a"], 1e6)
            _set_age(conn, galaxy["b"], 4000)
        counts = population.run_pass(conn)
        polities = {row["capital_system_id"]: row for row in conn.execute("SELECT * FROM polities").fetchall()}
        assert galaxy["a"] in polities and galaxy["b"] in polities
        assert polities[galaxy["a"]]["reach_ly"] == pytest.approx(population.reach_ly(1e6))
        owners = {row["star_system_id"]: row["polity_id"] for row in
                  conn.execute("SELECT * FROM system_owners").fetchall()}
        a, b = polities[galaxy["a"]]["id"], polities[galaxy["b"]]["id"]
        assert owners[galaxy["a"]] == a and owners[galaxy["b"]] == b
        assert owners[galaxy["plain"][0]] == a
        assert owners[galaxy["plain"][1]] == b        # 0.1 pc from B's capital: B's claim is stronger
        assert owners[galaxy["plain"][2]] == a        # midway: A's longer reach wins
        assert counts["owned_systems"] == len(owners)

        # Back below interstellar: the polity dissolves and its systems free up.
        with conn:
            _set_age(conn, galaxy["b"], 500)
        population.run_pass(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM polities WHERE id = ?", (b,)).fetchone()["n"] == 0
        assert conn.execute("SELECT polity_id FROM system_owners WHERE star_system_id = ?",
                            (galaxy["plain"][1],)).fetchone()["polity_id"] == a

        # Deleting a capital's system takes its species, polity and territory with it.
        with conn:
            conn.execute("DELETE FROM star_systems WHERE id = ?", (galaxy["a"],))
        assert conn.execute("SELECT COUNT(*) AS n FROM polities").fetchone()["n"] == 0
        assert conn.execute("SELECT COUNT(*) AS n FROM system_owners").fetchone()["n"] == 0

        # --rescan starts over.
        assert population.run_pass(conn, rescan=True)["new_species"] >= 1
    finally:
        conn.close()


def test_a_pass_removes_species_stored_without_a_civilization(mysql_config, galaxy):
    # Earlier passes named a species on every life world (GEN.80).
    conn = store.get_connection(mysql_config)
    try:
        population.run_pass(conn)
        # Any stored planet will do, but a planet is one species' homeworld at
        # most (`uq_species_homeworld`): the first one may be a civilized
        # system's, so clear its species before the leftover takes its place.
        planet = conn.execute("SELECT id, star_system_id FROM planets ORDER BY id LIMIT 1").fetchone()
        with conn:
            conn.execute("DELETE FROM species WHERE homeworld_planet_id = ?", (planet["id"],))
            conn.execute(
                "INSERT INTO species (name, homeworld_planet_id, star_system_id, life_chemical, life_stage, build, "
                "climate, size, civilization_age_years, era, spacefaring) "
                "VALUES ('Leftover', ?, ?, 'carbon', 'multicellularity', 'average', 'temperate', 'medium', "
                "NULL, NULL, 0)", (planet["id"], planet["star_system_id"]),
            )
        population.run_pass(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM species WHERE name = 'Leftover'").fetchone()["n"] == 0
        assert conn.execute("SELECT COUNT(*) AS n FROM species WHERE civilization_age_years IS NULL"
                            ).fetchone()["n"] == 0
    finally:
        conn.close()


def test_population_cli(mysql_config, galaxy, monkeypatch):
    argv = ["planetgen", "population", "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
            "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
            "--mysql-database", mysql_config.database]
    monkeypatch.setattr(sys, "argv", argv)
    generate_cli.main()
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM species").fetchone()["n"] >= 2
    finally:
        conn.close()

    monkeypatch.setattr(sys, "argv", argv + ["--rescan", "--territories-only"])
    with pytest.raises(SystemExit):
        generate_cli.main()




# ---------------------------------------------------------------------------
# API
# ---------------------------------------------------------------------------

@pytest.fixture
def client(mysql_config):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def test_api(mysql_config, galaxy, client):
    conn = store.get_connection(mysql_config)
    try:
        population.run_pass(conn)
        with conn:
            # Both capitals, not just A: B left alone keeps its drawn ages,
            # and two spacefaring species there found two polities with one
            # capital, the second owning no system (TEST.69).
            _set_age(conn, galaxy["a"], 1e6)
            _set_age(conn, galaxy["b"], 4000)
        population.run_pass(conn)
        homeworld = conn.execute("SELECT homeworld_planet_id FROM species WHERE star_system_id = ? "
                                 "AND spacefaring = 1", (galaxy["a"],)).fetchone()["homeworld_planet_id"]
    finally:
        conn.close()

    listing = client.get("/api/species").get_json()
    assert listing["total"] >= 2 and listing["items"]
    spacefaring = client.get("/api/species?spacefaring=true").get_json()
    assert spacefaring["total"] >= 1 and all(item["spacefaring"] for item in spacefaring["items"])
    assert client.get("/api/species?spacefaring=maybe").status_code == 400

    one = listing["items"][0]
    assert client.get(f"/api/species/{one['id']}").get_json()["name"] == one["name"]
    assert client.get("/api/species/999999999").status_code == 404
    assert client.get(f"/api/planets/{pid('planet', homeworld, mysql_config)}/species").get_json()["homeworld_planet_id"] == pid("planet", homeworld, mysql_config)

    polities = client.get("/api/polities").get_json()
    assert polities["total"] >= 1
    polity = client.get(f"/api/polities/{polities['items'][0]['id']}").get_json()
    assert polity["systems"] and polity["systems"][0]["distance_ly"] == 0.0
    assert client.get("/api/polities/999999999").status_code == 404

    owner = client.get(f"/api/systems/{pid('system', galaxy['a'])}/owner").get_json()["owner"]
    assert owner is not None and owner["distance_ly"] == 0.0
    territories = client.get("/api/territories").get_json()
    assert territories["points"] and territories["polities"]


def test_fills_skip_population_unless_asked(monkeypatch):
    """Off by default (Boss, 2026-10-01): a sector/galaxy run only runs the
    pass with --population."""
    import argparse
    calls = []
    monkeypatch.setattr(store, "get_connection", lambda *a, **k: calls.append(1))
    run_population.run_population_after(argparse.Namespace(population=False))
    assert calls == []
    parser, _commands = generate_cli.build_parser()
    assert not hasattr(parser.parse_args(["galaxy"]), "population")      # `galaxy` always runs it (PERF.123)
    assert parser.parse_args(["sector"]).population is False
    assert parser.parse_args(["sector", "--population"]).population is True


def test_population_status(mysql_config, client):
    empty = client.get("/api/population").get_json()
    assert empty == {"generated": False, "species": False, "polities": False, "territories": False}
    conn = store.get_connection(mysql_config)
    try:
        assert population.population_status(conn) == empty
        with conn:
            conn.execute("DROP TABLE system_owners")
        assert population.population_status(conn) == empty     # pre-v44 shape
    finally:
        conn.close()


def test_population_status_after_a_pass(mysql_config, galaxy, client):
    conn = store.get_connection(mysql_config)
    try:
        population.run_pass(conn)
        with conn:
            _set_age(conn, galaxy["a"], 1e6)
        population.run_pass(conn)
    finally:
        conn.close()
    assert client.get("/api/population").get_json() == {
        "generated": True, "species": True, "polities": True, "territories": True}


# --- The Species and Polities tables (UX.41) ---------------------------------------------

def test_species_and_polities_sort_filter_and_count(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        population.run_pass(conn)
        with conn:
            _set_age(conn, galaxy["a"], 1e6)
        population.run_pass(conn)

        names = [row["name"] for row in population.list_species(conn, sort="name")]
        assert names == sorted(names)
        assert [row["name"] for row in population.list_species(conn, sort="name", descending=True)] == names[::-1]
        for key in population.SPECIES_SORTS:
            assert len(population.list_species(conn, sort=key)) == len(names)

        everything = population.list_species(conn)
        eras = {row["era"] for row in everything}
        one = next(iter(eras))
        found = population.list_species(conn, eras=[one])
        assert found and {row["era"] for row in found} == {one}
        assert population.count_species(conn, eras=[one]) == len(found)
        assert population.count_species(conn, spacefaring=True) + population.count_species(
            conn, spacefaring=False) == len(everything)

        facets = population.species_facets(conn)
        assert sum(o["count"] for o in facets["spacefaring"]) == len(everything)
        assert {o["value"]: o["count"] for o in facets["era"]}[one] == len(found)
        # A menu ignores its own filter and the other narrows.
        narrowed = population.species_facets(conn, eras=[one])
        assert sum(o["count"] for o in narrowed["spacefaring"]) == len(found)
        assert narrowed["era"] == facets["era"]

        with pytest.raises(ValueError):
            population.list_species(conn, sort="mass")

        polities = population.list_polities(conn)
        if polities:
            polity_facets = population.polity_facets(conn)
            assert sum(o["count"] for o in polity_facets["government"]) == len(polities)
            government = polities[0]["government"]
            narrowed = population.list_polities(conn, governments=[government])
            assert narrowed and population.count_polities(conn, governments=[government]) == len(narrowed)
            for key in population.POLITY_SORTS:
                assert len(population.list_polities(conn, sort=key, descending=True)) == len(polities)
            detail = population.polity_detail(conn, polities[0]["id"], sort="name", descending=True)
            systems = [row["name"] for row in detail["systems"]]
            assert systems == sorted(systems, reverse=True)
        with pytest.raises(ValueError):
            population.list_polities(conn, sort="mass")
    finally:
        conn.close()


def test_the_pass_announces_each_step_for_its_bar(mysql_config, galaxy):
    """UX.83: the pass names its steps (a bar of its own), four, or one for territories only."""
    conn = store.get_connection(mysql_config)
    try:
        stages = []
        population.run_pass(conn, on_stage=stages.append)
        assert len(stages) == population.pass_stage_count() == 4
        only = []
        population.run_pass(conn, territories_only=True, on_stage=only.append)
        assert len(only) == population.pass_stage_count(True) == 1
    finally:
        conn.close()


def test_the_pass_draws_tech_levels_once_and_leaves_hand_edits_alone(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        population.run_pass(conn)
        rows = conn.execute("SELECT * FROM species ORDER BY id").fetchall()
        assert rows
        for row in rows:
            indices = {domain: row[f"{domain}_index"] for domain in tech.DOMAINS}
            uid = conn.execute("SELECT uid FROM planets WHERE id = ?", (row["homeworld_planet_id"],)).fetchone()["uid"]
            assert all(0 <= value <= tech.MAX_INDEX for value in indices.values())
            assert float(row["tech_level"]) == tech.tech_level(indices)
            assert row["tech_source"] == "generated"
            # The draw reads the stream of the homeworld, so it can be redone exactly.
            assert indices == tech.draw_indices(row["civilization_age_years"],
                                                draw.Stream(f"tech:{population.object_key(uid)}"))
        before = {row["id"]: row["tech_level"] for row in rows}

        # A second pass draws nothing new and changes nothing.
        population.run_pass(conn)
        assert {row["id"]: row["tech_level"] for row in conn.execute("SELECT * FROM species").fetchall()} == before

        first, second = rows[0]["id"], rows[1]["id"]
        with conn:
            # A species from before the columns (no indices) is filled; an edited one is not.
            conn.execute("UPDATE species SET energy_index = NULL, tech_level = NULL, tech_source = NULL "
                         "WHERE id = ?", (first,))
            conn.execute("UPDATE species SET energy_index = NULL, tech_level = NULL, tech_source = 'admin' "
                         "WHERE id = ?", (second,))
        population.run_pass(conn)
        filled = conn.execute("SELECT * FROM species WHERE id = ?", (first,)).fetchone()
        assert filled["energy_index"] is not None and filled["tech_source"] == "generated"
        assert float(filled["tech_level"]) == float(before[first])
        untouched = conn.execute("SELECT * FROM species WHERE id = ?", (second,)).fetchone()
        assert untouched["energy_index"] is None and untouched["tech_level"] is None

        detail = population.species_detail(conn, first)
        assert detail["tech_band"] == tech.band_name(detail["tech_level"])
        indices = {domain: detail[f"{domain}_index"] for domain in tech.DOMAINS}
        assert detail["limiting_domain"] == tech.limiting_domain(indices)
        assert detail["lopsided"] == tech.is_lopsided(indices)
        assert population.species_detail(conn, second)["tech_band"] is None
    finally:
        conn.close()


def test_species_search_by_tech_level_and_index(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        population.run_pass(conn)
        with conn:
            ids = [row["id"] for row in conn.execute("SELECT id FROM species ORDER BY id").fetchall()]
            for level, species_id in zip((1.0, 3.0, 5.0), ids):
                conn.execute("UPDATE species SET tech_level = ?, energy_index = ? WHERE id = ?",
                             (level, int(level), species_id))
        everything = population.list_species(conn, sort="tech")
        levels = [row["tech_level"] for row in everything if row["tech_level"] is not None]
        assert levels == sorted(levels)
        assert [row["tech_level"] for row in population.list_species(conn, sort="tech", descending=True)
                ][0] == max(levels)

        high = population.list_species(conn, tech_ranges={"tech_level": (2.0, None)})
        assert high and all(row["tech_level"] >= 2.0 for row in high)
        assert population.count_species(conn, tech_ranges={"tech_level": (2.0, None)}) == len(high)
        narrow = population.list_species(conn, tech_ranges={"tech_level": (2.0, 4.0), "energy_index": (3, 3)})
        assert [row["tech_level"] for row in narrow] == [3.0]
        assert population.species_facets(conn, tech_ranges={"tech_level": (2.0, None)})["spacefaring"] is not None
        with pytest.raises(ValueError):
            population.list_species(conn, tech_ranges={"tech_source": (1, 2)})
    finally:
        conn.close()


def test_api_species_tech_filters(mysql_config, galaxy, client):
    conn = store.get_connection(mysql_config)
    try:
        population.run_pass(conn)
        with conn:
            # Every species low first, then one at 6.0, so no species the pass made outranks it (TEST.299: the pass
            # can make one at 6.55, which led the tech-descending list).
            conn.execute("UPDATE species SET tech_level = 1.0, defense_index = 1")
            conn.execute("UPDATE species SET tech_level = 6.0, defense_index = 6 ORDER BY id LIMIT 1")
    finally:
        conn.close()
    listing = client.get("/api/species?sort=tech&order=desc").get_json()
    item = listing["items"][0]
    assert item["tech_level"] == 6.0 and item["defense_index"] == 6
    assert item["tech_band"] == tech.band_name(6.0) and item["lopsided"] in (True, False)
    assert item["limiting_domain"] in tech.DOMAINS[:3]
    high = client.get("/api/species?tech_min=5").get_json()
    assert high["total"] >= 1 and all(row["tech_level"] >= 5 for row in high["items"])
    one = client.get("/api/species?defense_min=6&defense_max=6").get_json()
    assert one["total"] == 1 and one["items"][0]["defense_index"] == 6
    assert client.get("/api/species?tech_min=high").status_code == 400
    assert client.get("/api/species?energy_max=9").status_code == 400
