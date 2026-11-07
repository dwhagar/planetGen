"""
GEN.1047: the wiring of the oddity table (`physics/oddities.py` `apply`, docs/design/systems/weird-worlds.md 4.2 and 4.5): the rows whose gate holds
are rolled on their own streams and their effects run at the step of the chain their lever belongs to; a forced oddity skips the roll but never the
gate; the ids are kept on the body, stored as text and read back; the page and the API say what is odd.
"""

import json
import types

import pytest

from planetgen.physics import oddities


def _table(*rows):
    document = {"kinds": ["planet", "moon", "giant", "rogue", "small body", "star system"]}
    with open(oddities.DATA_FILE, encoding="utf-8") as handle:
        shipped = json.load(handle)
    document.update({key: shipped[key] for key in ("levers", "quantities")})
    document["rows"] = list(rows)
    return oddities.parse(document)[0]


def _row(row_id="W91", lever="L3", input_="core_mass_fraction", rate=0.5, **changes):
    row = {
        "id": row_id, "name": "Test oddity", "kind": "planet",
        "requires": [{"quantity": "body_type", "equals": "t"}, {"quantity": "mass_earth", "lt": 10.0}],
        "rate": {"value": rate, "source": "a test", "label": "RECALLED", "tier": "T2", "tweak": "L0", "registry": "OD1"},
        "apply": {"lever": lever, "input": input_, "effect": f"test_effect_{row_id.lower()}", "params": {}},
        "limits": [{"quantity": "mass_earth", "lt": 10.0}],
        "source": "a test", "text": "A test cause.",
    }
    row.update(changes)
    return row


def _body(mass_kg=3.0e24, seed=1.0, **extra):
    body = types.SimpleNamespace(body_type="t", mass=mass_kg, radius=6000.0, distance=seed, eccentricity=0.0, star=None, is_moon=False,
                                 oddities=(), system_config=types.SimpleNamespace(FORCE_ODDITIES=None))
    body.__dict__.update(extra)
    return body


@pytest.fixture
def effects(monkeypatch):
    calls = []
    for row_id in ("W91", "W92"):
        monkeypatch.setitem(oddities.EFFECTS, f"test_effect_{row_id.lower()}", lambda body, row, calls=calls: calls.append((row.id, id(body))))
    return calls


def test_a_body_that_passes_the_gate_and_the_roll_is_given_the_row_and_the_effect_runs(effects):
    rows = _table(_row(rate=1.0))
    body = _body()
    assert oddities.apply(body, rows=rows, dial=1.0) == ["W91"]
    assert body.oddities == ("W91",) and effects == [("W91", id(body))]


def test_a_body_that_fails_the_gate_gets_nothing_even_at_a_certain_rate(effects):
    rows = _table(_row(rate=1.0))
    body = _body(mass_kg=9.0e25)         # 15 Earth masses, over the row's limit
    assert oddities.apply(body, rows=rows, dial=1.0) == [] and body.oddities == () and not effects


def test_the_dial_at_zero_turns_every_row_off_and_an_empty_table_changes_nothing(effects):
    body = _body()
    assert oddities.apply(body, rows=_table(_row(rate=1.0)), dial=0.0) == [] and body.oddities == ()
    assert oddities.apply(body, rows=()) == [] and not effects


def test_a_run_with_its_own_dial_uses_it_and_a_run_without_takes_the_tunings(effects, monkeypatch):
    rows = _table(_row(rate=1.0))
    off = _body()
    off.system_config.WEIRDNESS_DIAL = 0.0
    assert oddities.apply(off, rows=rows) == []
    monkeypatch.setattr("planetgen.tuning.WEIRDNESS_DIAL", 0.0)
    assert oddities.apply(_body(), rows=rows) == []
    monkeypatch.setattr("planetgen.tuning.WEIRDNESS_DIAL", 1.0)
    assert oddities.apply(_body(), rows=rows) == ["W91"]


def test_only_the_rows_of_the_levers_asked_for_run(effects):
    rows = _table(_row(rate=1.0, lever="L3"), _row("W92", rate=1.0, lever="L1", input_="obliquity_deg"))
    body = _body()
    assert oddities.apply(body, levers=("L1",), rows=rows, dial=1.0) == ["W92"]
    assert oddities.apply(body, levers=("L3",), rows=rows, dial=1.0) == ["W91"]
    assert body.oddities == ("W92", "W91")
    assert oddities.apply(body, rows=rows, dial=1.0) == [], "a row a body has is not given twice"


def test_a_forced_oddity_skips_the_roll_but_not_the_gate(effects):
    rows = _table(_row(rate=0.0))
    forced = _body()
    forced.system_config.FORCE_ODDITIES = ("W91",)
    assert oddities.apply(forced, rows=rows, dial=1.0) == ["W91"]
    barred = _body(mass_kg=9.0e25)
    barred.system_config.FORCE_ODDITIES = ("W91",)
    assert oddities.apply(barred, rows=rows, dial=1.0) == []


def test_a_row_without_a_registered_effect_is_refused_loudly(effects):
    rows = _table(_row("W93", rate=1.0))
    with pytest.raises(oddities.OddityError, match="no module registered"):
        oddities.apply(_body(), rows=rows, dial=1.0)


def test_the_roll_depends_on_the_body_and_not_on_the_other_rows_or_the_order_of_calls(effects):
    only = _table(_row(rate=0.5))
    both = _table(_row(rate=0.5), _row("W92", rate=0.5))
    results_only, results_both = [], []
    for n in range(60):
        a, b = _body(seed=1.0 + n), _body(seed=1.0 + n)
        results_only.append("W91" in oddities.apply(a, rows=only, dial=1.0))
        results_both.append("W91" in oddities.apply(b, rows=both, dial=1.0))
    assert results_only == results_both and 10 < sum(results_only) < 50


def test_the_body_key_needs_no_draw_and_separates_a_moon_from_a_planet():
    planet, moon = _body(), _body(is_moon=True, primary_mass_kg=6.0e24)
    assert oddities.body_key(planet) == oddities.body_key(_body()) and oddities.body_key(planet) != oddities.body_key(moon)


def test_the_quantities_name_only_what_the_body_has():
    found = oddities.quantities_of(_body(collision_count=2))
    assert found["body_type"] == "t" and found["collision_count"] == 2 and 0.4 < found["mass_earth"] < 0.6
    assert "obliquity_deg" not in found and "host_teff_k" not in found


def test_the_ids_are_stored_as_comma_separated_text_and_read_back():
    assert oddities.column_text(("W10", "G1")) == "W10,G1" and oddities.column_text(()) is None
    assert oddities.from_column("W10,G1") == ("W10", "G1") and oddities.from_column(None) == () and oddities.from_column("") == ()


def test_a_planet_keeps_its_oddities_through_to_dict_and_from_dict():
    from planetgen.generation.system import StarSystem
    from tests.test_serialization import make_system

    system = next(candidate for candidate in (make_system() for _ in range(20)) if any(p.body_type != "a" for p in candidate.planets))
    planet = next(p for p in system.planets if p.body_type != "a")
    assert planet.oddities == ()
    planet.oddities = ("W10",)
    reloaded = StarSystem.from_dict(json.loads(json.dumps(system.to_dict())))
    again = next(p for p in reloaded.planets if p.body_type != "a" and p.distance == planet.distance)
    assert again.oddities == ("W10",)


def test_the_page_says_what_is_odd_one_row_for_each_oddity():
    from planetgen.web import item_pages

    found = {"oddities": [{"id": "W10", "name": "Iron-rich core", "text": "A giant impact knocked its mantle off."}]}
    context = item_pages.with_oddities({"glance": [{"label": "Mass", "html": "1 kg"}]}, found)
    assert [row["label"] for row in context["glance"]] == ["Mass", "Odd"]
    assert "Iron-rich core. A giant impact" in str(context["glance"][1]["html"])
    plain = {"glance": []}
    assert item_pages.with_oddities(plain, {"oddities": []}) is plain


@pytest.mark.db
def test_the_ids_are_saved_with_a_planet_and_a_moon_and_the_api_names_them(mysql_config):
    from planetgen.api import routes
    from planetgen.db import store
    from tests.test_storage_checksum import _system

    system, config = _system()
    planets = [p for p in system.planets if p.body_type != "a"]
    assert planets, "the test system holds a planet"
    planets[0].oddities = ("W10",)
    if planets[0].moons:
        planets[0].moons[0].oddities = ("W10",)
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, config)
        conn.commit()
        row = conn.execute("SELECT id, oddities FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()
        assert [r["oddities"] for r in row].count("W10") == 1
        planet_id = next(r["id"] for r in row if r["oddities"])
        found = routes.oddities_of(conn, "planet", planet_id)
        assert [item["id"] for item in found] == ["W10"] and "mantle" in found[0]["text"]
        assert routes.oddities_of(conn, "star", 1) == []
        loaded = store.load_star_system(conn, system_id)
        assert loaded.planets[system.planets.index(planets[0])].oddities == ("W10",)
        if planets[0].moons:
            assert conn.execute("SELECT COUNT(*) AS n FROM moons WHERE oddities = 'W10'").fetchone()["n"] == 1
    finally:
        conn.close()


# --- W10 wired: the iron-rich core of a planet that had a recorded merger (GEN.1048) -----------------------------------------------------------

def _systems_with_a_merged_rocky_planet(force):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem

    found = []
    for _ in range(60):
        config = SystemConfig()
        config.STAR_TYPE, config.BINARY_SYSTEM, config.FORCE_ODDITIES = "G2V", False, force
        system = StarSystem(system_config=config)
        found.extend(p for p in system.planets if p.body_type == "t" and p.collision_count >= 1 and p.mass < 10 * 5.972e24)
        if len(found) >= 4:
            break
    assert found, "no rocky planet with a recorded merger came out in 60 systems"
    return found


def test_a_forced_iron_rich_core_keeps_the_mass_raises_the_iron_and_takes_the_water_with_the_mantle():
    for planet in _systems_with_a_merged_rocky_planet(("W10",)):
        assert planet.oddities == ("W10",)
        assert planet.water_mass_fraction == 0.0 and planet.ice_mass_fraction == 0.0
        assert 0.0 < planet.metal_mass_fraction < 0.85 and abs(planet.metal_mass_fraction + planet.rock_mass_fraction - 1.0) < 1e-9
        assert oddities.limit_failures(oddities.BY_ID["W10"], oddities.quantities_of(planet)) == []


def test_a_planet_with_no_recorded_merger_never_gets_the_iron_rich_core_even_when_it_is_forced():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem

    config = SystemConfig()
    config.STAR_TYPE, config.BINARY_SYSTEM, config.FORCE_ODDITIES = "G2V", False, ("W10",)
    for _ in range(10):
        for planet in StarSystem(system_config=config).planets:
            if planet.body_type == "t" and planet.collision_count == 0:
                assert "W10" not in planet.oddities


def test_without_forcing_the_iron_rich_core_is_rare_among_merged_planets():
    planets = _systems_with_a_merged_rocky_planet(None)
    assert sum("W10" in p.oddities for p in planets) <= len(planets) // 2 + 1
