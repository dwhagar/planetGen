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
    planet.oddities = ("W10",)      # (whatever it rolled for itself is replaced: the test is the round trip)
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
        first = row[system.planets.index(planets[0])]
        assert first["oddities"] == "W10"
        found = routes.oddities_of(conn, "planet", first["id"])
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
        assert "W10" in planet.oddities
        assert planet.water_mass_fraction == 0.0 and planet.ice_mass_fraction == 0.0
        assert 0.0 < planet.metal_mass_fraction < 0.85 and abs(planet.metal_mass_fraction + planet.rock_mass_fraction + (planet.envelope_mass_fraction or 0.0) - 1.0) < 1e-9      # (a hydrogen envelope that settles later, GEN.1323, takes its share of all the solids)
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


# --- GEN.1048: W2, the Venus-type slow retrograde rotation ------------------------------------------------------------------------

def _venus(tilt=170.0, age=4.5, moons=(), pressure=9.2e6, distance=0.7233):
    """A Venus: its mass, radius, air, orbit and the Sun's luminosity and mass, a given birth tilt and age."""
    star = types.SimpleNamespace(mass=1.98847e30, luminosity=3.828e26, age=age, temperature=5772.0)
    return _body(mass_kg=4.8675e24, radius=6051.8, distance=distance, eccentricity=0.007, star=star, atmospheric_pressure=pressure, moons=list(moons),
                 _birth_spin=(1.0e-5, tilt, ()), axial_tilt_deg=tilt, rotation_period_hours=24.0, period=distance ** 1.5)


def _force_w2(body):
    body.system_config.FORCE_ODDITIES = ("W2",)
    return oddities.apply(body, levers=("L1",))


def test_a_born_upside_down_venus_is_given_its_243_day_retrograde_spin():
    body = _venus()
    assert _force_w2(body) == ["W2"] and body.oddities == ("W2",)
    from planetgen.physics import venus_spin

    assert venus_spin.set_on(body) == 180.0 and abs(body.rotation_period_hours / 24.0 - 243.0) < 1.0


def test_a_planet_born_upright_ends_prograde_so_the_oddity_is_declined_and_not_recorded():
    body = _venus(tilt=20.0)
    assert _force_w2(body) == [] and body.oddities == () and body._venus_spin is None


def test_a_planet_whose_tide_has_not_had_the_age_is_declined():
    body = _venus(age=0.5)
    assert _force_w2(body) == [] and body.oddities == ()


def test_a_planet_with_a_big_moon_or_a_thin_air_fails_the_gate_even_when_forced():
    moon = types.SimpleNamespace(mass=7.3e22)
    assert _force_w2(_venus(moons=[moon])) == []
    assert _force_w2(_venus(pressure=1.0e5)) == []


def test_a_day_under_twenty_days_is_declined_by_the_limit_of_the_row():
    body = _venus(distance=0.2)           # the same tide balance far closer in: the retrograde day is days, not months
    assert _force_w2(body) == [] and body.oddities == ()


def test_dropping_a_lever_takes_its_rows_off_and_a_second_roll_gives_them_back():
    body = _venus()
    _force_w2(body)
    oddities.drop(body, ("L3",))
    assert body.oddities == ("W2",)
    oddities.drop(body, ("L1",))
    assert body.oddities == ()
    assert _force_w2(body) == ["W2"]


def test_generate_spin_puts_the_state_on_the_rotation_and_the_tilt():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.planet import Planet
    from planetgen.generation.star import Star
    from planetgen.physics import planets
    from planetgen.util import draw

    with draw.bound("w2-generate-spin"):
        config = SystemConfig()
        config.STAR_TYPE, config.BINARY_SYSTEM, config.FORCE_ODDITIES = "G2V", False, ("W2",)
        star = Star(config)
        planet = Planet(config, star, star.habitable_zone, 0.7, moon_count=0)
        star.mass, star.luminosity, star.age = 1.98847e30, 3.828e26, 4.5
        planet.mass, planet.radius, planet.eccentricity, planet.distance, planet.atmospheric_pressure = 4.8675e24, 6051.8, 0.007, 0.7233, 9.2e6
        planets.generate_spin(planet, star.mass, birth=(planet._birth_spin[0], 170.0, planet._birth_spin[2]))
        assert planet.oddities == ("W2",) and planet.axial_tilt_deg == 180.0 and abs(planet.rotation_period_hours / 24.0 - 243.0) < 1.0
        planet.atmospheric_pressure = 1.0e5       # made again with a thin air: the row's gate fails and the Venus state goes
        planets.generate_spin(planet, star.mass, birth=(planet._birth_spin[0], 170.0, planet._birth_spin[2]))
        assert planet.oddities == () and planet._venus_spin is None


def test_a_planet_tilted_far_over_with_no_moon_is_named_for_its_seasons_and_nothing_on_it_changes():
    body = _venus(tilt=95.0, pressure=1.0e5)
    before = dict(body.__dict__)
    body.system_config.FORCE_ODDITIES = None
    assert oddities.apply(body, levers=("L1",), dial=1.0) == ["W1"]
    assert body.oddities == ("W1",) and {k: v for k, v in body.__dict__.items() if k != "oddities"} == {k: v for k, v in before.items() if k != "oddities"}


def test_the_seasons_row_wants_a_tilt_of_sixty_to_one_hundred_twenty_degrees_and_no_big_moon():
    for tilt, moons in ((30.0, ()), (150.0, ()), (95.0, [types.SimpleNamespace(mass=7.3e22)])):
        body = _venus(tilt=tilt, pressure=1.0e5, moons=moons)
        assert oddities.apply(body, levers=("L1",), dial=1.0) == [], tilt


def test_a_born_tilted_venus_that_the_tides_turn_over_is_not_named_for_its_seasons():
    body = _venus(tilt=100.0)         # born on its side, thick air: the tide balance turns it upside down (180) and its seasons are the upright ones
    ids = oddities.apply(body, levers=("L1",), dial=1.0, rows=None)
    assert ids == ["W2"] and body.axial_tilt_deg == 180.0


# --- GEN.1058: G1 (a giant tipped by an impact) --------------------------------------------------------------------------------

def _tipped_giant(tilt, mass_earth=14.5):
    star = types.SimpleNamespace(mass=1.98847e30, luminosity=3.828e26, age=4.5, temperature=5772.0, radius=6.957e5)
    return _body(mass_kg=mass_earth * 5.972e24, radius=25362.0, distance=19.2, eccentricity=0.047, star=star, body_type="g", moons=[], axial_tilt_deg=tilt)


def test_an_ice_giant_the_merger_laid_on_its_side_is_named_and_nothing_on_it_changes():
    body = _tipped_giant(97.8)
    before = dict(body.__dict__)
    assert oddities.apply(body, levers=("L1",), dial=1.0) == ["G1"]
    assert body.oddities == ("G1",) and {k: v for k, v in body.__dict__.items() if k != "oddities"} == {k: v for k, v in before.items() if k != "oddities"}


def test_the_giant_tilt_row_wants_sixty_degrees_and_a_mass_of_five_to_thirty_earths():
    for tilt, mass in ((28.3, 17.1), (97.8, 3.0), (97.8, 95.0), (170.0, 14.5)):
        body = _tipped_giant(tilt, mass)
        assert oddities.apply(body, levers=("L1",), dial=1.0) == (["G1"] if (tilt, mass) == (170.0, 14.5) else []), (tilt, mass)


def test_a_generated_ice_giant_is_read_for_the_giant_tilt_row_when_the_spin_is_made():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.planet import Planet
    from planetgen.generation.star import Star
    from planetgen.physics import planets
    from planetgen.util import draw

    with draw.bound("g1-generate-spin"):
        config = SystemConfig()
        config.STAR_TYPE, config.BINARY_SYSTEM, config.FORCE_ODDITIES = "G2V", False, ("G1",)
        star = Star(config)
        planet = Planet(config, star, star.habitable_zone, 19.2, moon_count=0, mass=14.5 * 5.972e24)
        assert planet.body_type == "g"
        planets.generate_spin(planet, star.mass, birth=(planet._birth_spin[0], 98.0, planet._birth_spin[2]))
        assert planet.oddities == ("G1",) and 97.0 < planet.axial_tilt_deg < 99.0
        planets.generate_spin(planet, star.mass, birth=(planet._birth_spin[0], 20.0, planet._birth_spin[2]))
        assert planet.oddities == ()


# --- GEN.1049: W3 (high eccentricity), W4 (polar orbit) and W5 (resonant chain) ---------------------------------------------------

def _orbit_body(kind="t", eccentricity=0.0, teff=5800.0, distance=1.0, radius_km=6.957e5, mass=6.0e24):
    star = types.SimpleNamespace(mass=1.98847e30, luminosity=3.828e26, age=4.5, temperature=teff, radius=radius_km)
    return _body(mass_kg=mass, radius=6371.0, distance=distance, eccentricity=eccentricity, star=star, body_type=kind, moons=[], orbital_inclination_deg=1.0)


def test_an_orbit_of_eccentricity_half_or_more_is_named_and_nothing_on_it_moves():
    body = _orbit_body(eccentricity=0.6)
    body.system_config.FORCE_ODDITIES = None
    assert oddities.apply(body, levers=("L2",), dial=1.0) == ["W3"] and body.eccentricity == 0.6
    giant = _orbit_body(kind="g", eccentricity=0.93, mass=1.0e27)
    assert oddities.apply(giant, levers=("L2",), dial=1.0) == ["W31"]
    assert oddities.apply(_orbit_body(eccentricity=0.3), levers=("L2",), dial=1.0) == []


def test_a_close_giant_round_a_hot_star_can_be_turned_to_a_polar_orbit_in_the_peak_of_the_measured_angles():
    from planetgen.physics import orbit_oddities

    giant = _orbit_body(kind="g", teff=7500.0, distance=0.03, mass=1.0e27)
    giant.system_config.FORCE_ODDITIES = ("W43",)
    assert oddities.apply(giant, levers=("L2",)) == ["W43"]
    low, high = orbit_oddities.POLAR_PEAK_DEG
    assert low <= giant.orbital_inclination_deg <= high
    again = _orbit_body(kind="g", teff=7500.0, distance=0.03, mass=1.0e27)
    again.system_config.FORCE_ODDITIES = ("W43",)
    oddities.apply(again, levers=("L2",))
    assert again.orbital_inclination_deg == giant.orbital_inclination_deg


def test_the_host_bands_gate_the_polar_rows_and_a_far_giant_never_gets_one():
    for teff, forced, distance in ((5000.0, "W43", 0.03), (7500.0, "W41", 0.03), (5000.0, "W41", 2.0)):
        giant = _orbit_body(kind="g", teff=teff, distance=distance, mass=1.0e27)
        giant.system_config.FORCE_ODDITIES = (forced,)
        assert oddities.apply(giant, levers=("L2",), dial=0.0) == [], (teff, forced, distance)
    assert oddities.BY_ID["W41"].rate.value < oddities.BY_ID["W42"].rate.value < oddities.BY_ID["W43"].rate.value


def _chain(count=3):
    """A row of planets on the Laplace chain of Io, Europa and Ganymede periods (2:1, 2:1) round a Sun, with the masses of three 10 Earth-mass planets."""
    star_mass = 1.98847e30
    planets = []
    for index in range(count):
        period_days = 3.0 * 2.0 ** index * 1.005 ** index
        planet = _orbit_body(mass=6.0e25, distance=(period_days / 365.25) ** (2.0 / 3.0))
        planet.period = period_days / 365.25
        planets.append(planet)
    return planets, star_mass


def test_a_chain_of_three_planets_gives_each_a_forced_eccentricity_and_a_pair_gives_none():
    from planetgen.physics import resonant_chains

    planets, star_mass = _chain(3)
    assert resonant_chains.give_forced_eccentricities(planets, star_mass) == 3
    assert all(p.oddity_quantities["forced_eccentricity"] > 0.0 for p in planets)
    pair, star_mass = _chain(2)
    assert resonant_chains.give_forced_eccentricities(pair, star_mass) == 0 and "forced_eccentricity" not in pair[0].oddity_quantities


def test_a_forced_resonance_raises_the_eccentricity_keeps_it_through_a_second_roll_and_declines_when_there_is_nothing_to_add():
    from planetgen.physics import resonant_chains

    planets, star_mass = _chain(3)
    resonant_chains.give_forced_eccentricities(planets, star_mass)
    body = planets[1]
    body.system_config.FORCE_ODDITIES = ("W5",)
    assert oddities.apply(body, levers=("L2",)) == ["W5"]
    held = body.eccentricity
    assert held >= 0.001 and body._orbit_birth is None
    oddities.drop(body, ("L2",))
    assert oddities.apply(body, levers=("L2",)) == ["W5"] and body.eccentricity == held      # the same planet again: still held, not raised twice
    rich = planets[2]
    rich.eccentricity = 0.5
    rich.system_config.FORCE_ODDITIES = ("W5",)
    assert "W5" not in oddities.apply(rich, levers=("L2",))


def test_a_planet_held_on_a_forced_eccentricity_is_flexed_by_its_star():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.planet import Planet
    from planetgen.generation.star import Star
    from planetgen.physics import planets
    from planetgen.util import draw

    with draw.bound("w5try"):
        config = SystemConfig()
        config.STAR_TYPE, config.BINARY_SYSTEM, config.FORCE_ODDITIES = "M3V", False, ("W5",)
        star = Star(config)
        planet = Planet(config, star, star.habitable_zone, 0.03, moon_count=0)
        before = planet.tidal_flux_w_m2 or 0.0
        planet.oddity_quantities = {"forced_eccentricity": 0.02}
        assert planets.finish_orbit_oddities(planet) is True
        assert planet.oddities == ("W5",) and abs(planet.eccentricity - 0.02) < 1e-12 and planet.tidal_flux_w_m2 > before


# --- GEN.1050: W14 (a disintegrating planet) and W9 (a ring from a fallen moon) ------------------------------------------------------

def _hot_rock(mass_earth=0.05, age=4.5, teff=3800.0, radius_km=4.2e5, distance=0.0097):
    star = types.SimpleNamespace(mass=0.6 * 1.98847e30, luminosity=1.0e25, age=age, temperature=teff, radius=radius_km)
    period = (distance ** 3 / 0.6) ** 0.5
    return _body(mass_kg=mass_earth * 5.972e24, radius=3000.0, distance=distance, eccentricity=0.0, star=star, period=period, moons=[])


def test_a_small_rock_with_a_dayside_over_two_thousand_kelvin_and_a_day_long_orbit_is_named_for_its_wind():
    body = _hot_rock()
    body.system_config.FORCE_ODDITIES = None
    assert oddities.quantities_of(body)["dayside_temperature_k"] > 2000.0
    assert oddities.apply(body, levers=("L5",), dial=1.0) == ["W14"]


def test_a_planet_that_the_wind_would_have_used_up_or_that_is_too_cool_or_too_heavy_is_not_named():
    for body in (_hot_rock(teff=5000.0), _hot_rock(distance=0.5), _hot_rock(mass_earth=0.5), _hot_rock(age=4.5, mass_earth=0.02, teff=3800.0)):
        body.system_config.FORCE_ODDITIES = ("W14",)
        assert oddities.apply(body, levers=("L5",)) == [], oddities.quantities_of(body)


def test_a_planet_whose_ring_came_from_a_moon_the_tide_pulled_in_is_named_and_one_from_a_comet_is_not():
    for origin, expected in (("disrupted_moon", ["W9"]), ("captured_body", []), ("moon_erosion", []), (None, [])):
        body = _body(ring_origin=origin)
        assert oddities.apply(body, levers=("L6",), dial=1.0) == expected, origin
    giant = _body(body_type="g", ring_origin="disrupted_moon")
    assert oddities.apply(giant, levers=("L6",), dial=1.0) == ["W91"]


# --- GEN.1051: W19 (nightside collapse) and W20 (the other climate branch) ----------------------------------------------------------

def _planet_at(star_type, distance_au, seed, force=(), dial=None, mass_earth=1.0):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.planet import Planet
    from planetgen.generation.star import Star
    from planetgen.util import draw

    with draw.bound(seed):
        config = SystemConfig()
        config.STAR_TYPE, config.BINARY_SYSTEM, config.FORCE_ODDITIES = star_type, False, tuple(force)
        if dial is not None:
            config.WEIRDNESS_DIAL = dial
        star = Star(config)
        return Planet(config, star, star.habitable_zone, distance_au, mass=mass_earth * 5.972e24, moon_count=0)


def test_a_locked_planet_with_a_thin_air_loses_it_to_its_dark_side_and_is_left_airless():
    lost = _planet_at("M3V", 0.12, "w19-M3V-2", force=("W19",))
    assert lost.oddities == ("W19",) and lost.atmospheric_pressure == 0.0
    kept = _planet_at("M3V", 0.12, "w19-M3V-2", dial=0.0)
    assert kept.oddities == () and kept.atmospheric_pressure > 0.0


def test_a_planet_that_is_not_locked_never_loses_its_air_to_a_dark_side_even_when_forced():
    far = _planet_at("G2V", 1.0, "w19-far", force=("W19",))
    assert "W19" not in far.oddities


def test_the_other_branch_of_a_bistable_climate_is_a_distinct_climate_from_the_first_answer():
    # a pinned seed whose planet has a second branch; the generator's draws moved at epoch 140 (GEN.1363, GEN.1374) and again with GEN.993 (the core's light elements), so the seed that had one was replaced (and again at GEN.1646, when the other branch began to repeat the first solve's draws, so the world keeps its inventory and a seed that had a second branch only through a new inventory lost it; this one is 20: 50 K and 0.34 apart)
    other = _planet_at("G2V", 1.2, "w20-G2V-20", force=("W20",))
    first = _planet_at("G2V", 1.2, "w20-G2V-20", dial=0.0)
    assert "W20" in other.oddities and first.oddities == ()
    assert abs(other.surface_temperature - first.surface_temperature) >= 20.0 and abs(other.bond_albedo - first.bond_albedo) >= 0.1


def test_the_other_branch_is_the_same_world_with_the_same_inventory_and_draws():
    """GEN.1646: the second solve of W20 repeats every draw of the first (design constants, mixture seed, late accretion); only the start temperature differs."""
    from planetgen.physics import surface_air

    airs = []
    real = surface_air.solve

    def spy(setting, start_temperature_k):
        air = real(setting, start_temperature_k)
        airs.append(air)
        return air

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(surface_air, "solve", spy)
        planet = _planet_at("G2V", 1.2, "w20-G2V-20", force=("W20",))
    assert "W20" in planet.oddities and len(airs) >= 2
    first, second = airs[0], airs[1]
    assert first.entry_state is not None and second.entry_state == first.entry_state
    assert second.constants == first.constants and second.late == first.late
    assert second.mixture.inventory == first.mixture.inventory
    assert second.temperature_k != first.temperature_k      # a different branch of the same world


def test_a_world_with_no_second_branch_is_not_named_for_one():
    hot = _planet_at("G2V", 0.4, "w20-hot", force=("W20",))
    assert "W20" not in hot.oddities


# --- GEN.1523: W26 (the Hycean box of `hycean_candidate`) ------------------------------------------------------------------------

def test_a_sub_neptune_inside_the_hycean_box_is_named_and_one_that_is_too_light_for_its_size_is_not():
    from planetgen.physics import hycean_candidate

    inside = _planet_at("K2V", 0.3, "hy-2-0.3-0", mass_earth=2.0, force=("W26",))
    found = oddities.quantities_of(inside)
    assert "W26" in inside.oddities and set(inside.oddities) <= {"W26", "W1"}      # the draw may also tip the planet into the tilt row W1
    assert hycean_candidate.hycean_candidate(found["mass_earth"], found["radius_earth"], found["equilibrium_temperature_k"])
    puffy = _planet_at("K2V", 1.0, "hy-2-1.0-0", mass_earth=2.0, force=("W26",))
    found = oddities.quantities_of(puffy)       # 2.0 Earth radii at 1.4 g/cc: under the box's density
    assert found["density_g_cm3"] < 1.5 and "W26" not in puffy.oddities
    assert not hycean_candidate.hycean_candidate(found["mass_earth"], found["radius_earth"], found["equilibrium_temperature_k"])


def test_the_gate_of_the_hycean_row_reads_the_same_box_as_the_module():
    row = oddities.BY_ID["W26"]
    from planetgen.physics import hycean_candidate as box

    bounds = {(c["quantity"], "min" if "min" in c else "lt"): c.get("min", c.get("lt")) for c in row.requires if c["quantity"] != "body_type"}
    assert bounds[("mass_earth", "min")] == box.HYCEAN_MASS_EARTH[0] and abs(bounds[("mass_earth", "lt")] - box.HYCEAN_MASS_EARTH[1]) < 1e-6
    assert bounds[("radius_earth", "min")] == box.HYCEAN_RADIUS_EARTH[0] and abs(bounds[("radius_earth", "lt")] - box.HYCEAN_RADIUS_EARTH[1]) < 1e-6
    assert bounds[("equilibrium_temperature_k", "min")] == box.HYCEAN_TEMPERATURE_K[0] and abs(bounds[("equilibrium_temperature_k", "lt")] - box.HYCEAN_TEMPERATURE_K[1]) < 1e-6
    assert bounds[("density_g_cm3", "min")] == box.HYCEAN_DENSITY_G_CM3[0] and abs(bounds[("density_g_cm3", "lt")] - box.HYCEAN_DENSITY_G_CM3[1]) < 1e-6


def test_the_dial_at_zero_names_no_hycean_candidate_and_a_moon_has_no_equilibrium_temperature_quantity():
    off = _planet_at("K2V", 0.3, "hy-2-0.3-0", mass_earth=2.0, dial=0.0)
    assert off.oddities == ()
    moon = _body(star=types.SimpleNamespace(luminosity=3.8e26), is_moon=True)
    assert "equilibrium_temperature_k" not in oddities.quantities_of(moon)


# --- GEN.1527: W16 (the helium-rich envelope of `helium_fractionation`) ----------------------------------------------------------

def test_a_planet_that_lost_most_of_its_envelope_ends_richer_in_helium_than_the_star_and_one_that_lost_little_does_not():
    from planetgen.physics import bulk_draw

    stripped = _planet_at("K2V", 0.1, "he2-K2V-6-0.1-0", mass_earth=6.0)         # lost more than half of the envelope it started with
    host = bulk_draw.host_ratios(stripped.star)
    assert stripped._envelope_helium_y > host.y_envelope and stripped.envelope_helium_fraction > host.y_envelope * (1.0 - host.z)
    far = _planet_at("G2V", 0.04, "he-3-0.04-0", mass_earth=3.0)
    assert getattr(far, "_envelope_helium_y", None) is None                     # no envelope at all


def test_the_helium_row_names_the_planets_that_reached_the_enhanced_share_and_no_others():
    from planetgen.physics import helium_fractionation as hf

    row = oddities.BY_ID["W16"]
    assert any(c["quantity"] == "envelope_helium_y" and c.get("min") == hf.MIN_HELIUM_FOR_ENHANCED for c in row.requires)
    rich = _body(_envelope_helium_y=0.5)
    assert oddities.apply(rich, levers=("L8",), dial=1.0) == ["W16"]
    for y in (0.3, None):
        poor = _body(_envelope_helium_y=y)
        assert oddities.apply(poor, levers=("L8",), dial=1.0) == []


# --- GEN.1047 isolation: the dial at 0 changes no other draw ---------------------------------------------------------------------

def _bodies_of(system):
    return [body for planet in system.all_planets if planet.body_type != "a" for body in (planet, *getattr(planet, "moons", ()))]


def test_the_dial_at_zero_gives_the_same_system_as_the_dial_at_one_wherever_nothing_was_rolled():
    """Every oddity draws from a stream of its own, so a system in which no row rolled in is the same system whatever the dial; one with an oddity differs only where the oddity acts."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    compared = 0
    for number in range(12):
        made = {}
        for dial in (0.0, 1.0):
            with draw.bound(f"isolation:{number}"):
                config = SystemConfig()
                config.STAR_TYPE, config.WEIRDNESS_DIAL = "G2V", dial
                made[dial] = StarSystem(system_config=config)
        assert not any(getattr(body, "oddities", ()) for body in _bodies_of(made[0.0])), "the dial at 0 rolls nothing"
        if not any(getattr(body, "oddities", ()) for body in _bodies_of(made[1.0])):
            assert json.dumps(made[0.0].to_dict(), sort_keys=True, default=str) == json.dumps(made[1.0].to_dict(), sort_keys=True, default=str), number
            compared += 1
    assert compared >= 4, "too few systems without an oddity to compare"


# --- UX.285: the dial in the settings file and page -----------------------------------------------------------------------------

def test_the_settings_dial_is_used_when_the_run_names_none_and_the_file_names_one(monkeypatch):
    from planetgen.util import settings

    named = settings.build({"generation": {"weirdness_dial": 0.0}}, {})
    unnamed = settings.build({}, {})
    assert named.generation.weirdness_dial == 0.0 and unnamed.generation.weirdness_dial == 1.0

    body = _venus(tilt=95.0, pressure=1.0e5)
    body.system_config.FORCE_ODDITIES = None
    body.system_config.WEIRDNESS_DIAL = None
    monkeypatch.setattr(settings, "get_settings", lambda environ=None: named)
    assert oddities.site_dial() == 0.0 and oddities.apply(body, levers=("L1",)) == []         # the file turns the oddities off

    monkeypatch.setattr(settings, "get_settings", lambda environ=None: unnamed)
    assert oddities.site_dial() is None                                                       # nothing named: the tuning's dial
    assert oddities.apply(body, levers=("L1",)) == ["W1"]


def test_a_run_that_names_its_own_dial_ignores_the_settings_file(monkeypatch):
    from planetgen.util import settings

    named = settings.build({"generation": {"weirdness_dial": 0.0}}, {})
    monkeypatch.setattr(settings, "get_settings", lambda environ=None: named)
    body = _venus(tilt=95.0, pressure=1.0e5)
    body.system_config.FORCE_ODDITIES = None
    body.system_config.WEIRDNESS_DIAL = 1.0
    assert oddities.apply(body, levers=("L1",)) == ["W1"]


def test_a_settings_file_that_does_not_fit_the_model_leaves_the_tunings_dial(monkeypatch):
    from planetgen.util import settings

    def broken(environ=None):
        raise settings.SettingsError("bad file")

    monkeypatch.setattr(settings, "get_settings", broken)
    assert oddities.site_dial() is None


def test_the_dial_is_an_option_of_the_world_generation_group():
    from planetgen.util import settings_view

    assert "generation" in settings_view.GROUPS
    from planetgen.util import settings

    field = settings.Generation.model_fields["weirdness_dial"]
    assert field.json_schema_extra["x-group"] == "generation" and field.default == 1.0


def test_a_forced_eccentricity_is_held_where_the_neighbours_keep_their_spacing():
    """The resonance raises an eccentricity after `space_orbits` placed the orbits, so the cap keeps the mutual Hill spacing of `clearance_ledger.planet_exclusion` (found with GEN.1365's draws)."""
    from planetgen.physics import clearance_ledger, resonant_chains

    host = 3.2e29
    planets = [types.SimpleNamespace(distance=distance, eccentricity=e, mass=1.5e24) for distance, e in ((0.0266, 0.0289), (0.03179, 0.0), (0.0385, 0.0025))]
    held = resonant_chains._clear_of_neighbours(planets, [0.0, 0.0119, 0.0], host)
    assert all(0.0 <= value <= 0.0119 for value in held)
    assert held[1] < 0.0119       # the middle planet has no room for the whole of it (a system the generator drew)
    for inner, outer, e_in, e_out in zip(planets, planets[1:], held, held[1:]):
        kappa = clearance_ledger.spacing_kappa(inner.mass, outer.mass, host)
        assert outer.distance * (1.0 - e_out - kappa / 2.0) >= inner.distance * (1.0 + e_in + kappa / 2.0) - 1e-12
    far = [types.SimpleNamespace(distance=distance, eccentricity=0.0, mass=1.5e24) for distance in (0.5, 1.0, 2.0)]
    assert resonant_chains._clear_of_neighbours(far, [0.01, 0.02, 0.03], host) == [0.01, 0.02, 0.03]    # a wide gap changes nothing
