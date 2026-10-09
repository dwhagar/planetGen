"""
Admin editing (TODO ADM.1): the delete and regenerate endpoints for one
planet, moon, asteroid belt, phenomenon or sector (ADM.8), the class
(ADM.6) and star (ADM.7) changes (`src/planetgen/api/edits.py`), the in-place writer behind them
(`planetgen/db/edits.py`) and the object edits
(`planetgen/admin/edits.py`), against a real throwaway database.
"""

import pytest
from markupsafe import escape

from planetgen.util import draw
from planetgen.web.app import create_app
from planetgen.api.config import Config
from planetgen.db import edits as editStore, store
from planetgen.admin import auth as adminAuth, edits as adminEdits
from planetgen.generation import validation
from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem


@pytest.fixture
def client(mysql_config):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


@pytest.fixture
def admin(client, mysql_config):
    store.get_connection(mysql_config).close()
    _username, password = adminAuth.bootstrap_control_schema(mysql_config)
    assert client.post("/api/auth/login", json={"username": "admin", "password": password}).status_code == 200
    assert client.post("/api/auth/change-credentials", json={
        "current_password": password, "new_username": "boss", "new_password": "a-long-new-password-1",
    }).status_code == 200
    return client


def _saved_system(mysql_config, want=None):
    """A saved single-star sector system with planets, moons and a belt
    where possible, retried until `want(system)` holds (by default, until
    it has at least one planet)."""
    if want is None:
        def want(system):
            return any(p.body_type != 'a' for p in system.planets)
    for _ in range(300):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MOONS = True
        cfg.ASTEROID_BELT = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(cfg)
        if want(system):
            sector = SpaceSector("Edit Sector", edge_ly=10.0)
            sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
            sector_id = store.save_sector(sector, config=mysql_config)
            conn = store.get_connection(mysql_config)
            try:
                system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?",
                                         (sector_id,)).fetchone()["id"]
            finally:
                conn.close()
            return sector_id, system_id
    pytest.skip("no suitable system generated")


def _with_moons(system):
    return any(p.body_type != 'a' and p.moons for p in system.planets)


def _rows(mysql_config, sql, params=()):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute(sql, params).fetchall()
    finally:
        conn.close()


def _load(mysql_config, system_id):
    conn = store.get_connection(mysql_config)
    try:
        return store.load_star_system(conn, system_id)
    finally:
        conn.close()


def test_edits_need_an_admin(client, mysql_config):
    store.get_connection(mysql_config).close()
    assert client.post("/api/planets/1/regenerate").status_code == 401
    assert client.delete("/api/planets/1").status_code == 401
    assert client.delete("/api/phenomena/nebula/1").status_code == 401
    assert client.post("/api/sectors/1/regenerate").status_code == 401


def test_save_system_edits_round_trips_an_unchanged_system(mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_moons)
    before = _load(mysql_config, system_id)
    planet_ids = sorted(r["id"] for r in _rows(mysql_config, "SELECT id FROM planets WHERE star_system_id = ?",
                                                   (system_id,)))
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            editStore.save_system_edits(conn, system_id, store.load_star_system(conn, system_id))
    finally:
        conn.close()
    after = _load(mysql_config, system_id)
    assert sorted(r["id"] for r in _rows(mysql_config, "SELECT id FROM planets WHERE star_system_id = ?",
                                         (system_id,))) == planet_ids
    for old, new in zip(before.planets, after.planets):
        assert getattr(old, "name", None) == getattr(new, "name", None)
        assert getattr(old, "planet_class", None) == getattr(new, "planet_class", None)
        assert old.distance == pytest.approx(new.distance)


def test_regenerate_planet_keeps_its_row_and_name(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_moons)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type != 'a' and p.moons)
    other_ids = {r["id"] for r in _rows(mysql_config, "SELECT id FROM planets WHERE star_system_id = ? AND id <> ?",
                                        (system_id, planet.db_id))}
    response = admin.post(f"/api/planets/{planet.db_id}/regenerate")
    assert response.status_code == 200, response.get_json()
    body = response.get_json()
    assert body["summary"].startswith(f"Regenerated {planet.name}")
    row = _rows(mysql_config, "SELECT name FROM planets WHERE id = ?", (planet.db_id,))[0]
    assert row["name"] == planet.name
    assert {r["id"] for r in _rows(mysql_config, "SELECT id FROM planets WHERE star_system_id = ? AND id <> ?",
                                   (system_id, planet.db_id))} == other_ids
    for moon in _rows(mysql_config, "SELECT name FROM moons WHERE planet_id = ?", (planet.db_id,)):
        assert moon["name"].startswith(planet.name)


def test_delete_planet_takes_its_moons(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_moons)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type != 'a' and p.moons)
    response = admin.delete(f"/api/planets/{planet.db_id}")
    assert response.status_code == 200, response.get_json()
    assert _rows(mysql_config, "SELECT id FROM planets WHERE id = ?", (planet.db_id,)) == []
    assert _rows(mysql_config, "SELECT id FROM moons WHERE planet_id = ?", (planet.db_id,)) == []
    assert admin.delete(f"/api/planets/{planet.db_id}").status_code == 404


def test_regenerate_and_delete_a_moon(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_moons)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type != 'a' and p.moons)
    moon = planet.moons[0]
    response = admin.post(f"/api/moons/{moon.db_id}/regenerate")
    assert response.status_code in (200, 409), response.get_json()
    if response.status_code == 200:
        assert _rows(mysql_config, "SELECT name FROM moons WHERE id = ?", (moon.db_id,))[0]["name"] == moon.name
    assert admin.delete(f"/api/moons/{moon.db_id}").status_code == 200
    assert _rows(mysql_config, "SELECT id FROM moons WHERE id = ?", (moon.db_id,)) == []


def test_regenerate_and_delete_a_belt(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, lambda s: any(p.body_type == 'a' for p in s.planets))
    belt = next(p for p in _load(mysql_config, system_id).planets if p.body_type == 'a')
    response = admin.post(f"/api/belts/{belt.db_id}/regenerate")
    assert response.status_code == 200, response.get_json()
    row = _rows(mysql_config, "SELECT lower_limit_km FROM asteroid_belts WHERE id = ?", (belt.db_id,))[0]
    assert row["lower_limit_km"] == pytest.approx(belt.lower_limit * 1.495978707e8, rel=1e-6)
    assert admin.delete(f"/api/belts/{belt.db_id}").status_code == 200
    assert _rows(mysql_config, "SELECT id FROM asteroid_belts WHERE id = ?", (belt.db_id,)) == []


def test_delete_refuses_to_drop_facilities_unless_told(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type != 'a')
    response = admin.post("/api/facilities", json={
        "name": "Relay One", "kind": "station", "placement": "orbital", "host_type": "planet",
        "host_id": planet.db_id})
    assert response.status_code == 201, response.get_json()
    assert admin.delete(f"/api/planets/{planet.db_id}").status_code == 409
    assert admin.delete(f"/api/planets/{planet.db_id}", json={"drop_facilities": True}).status_code == 200


def test_regenerate_planet_keeps_facilities_on_the_planet(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type != 'a')
    assert admin.post("/api/facilities", json={
        "name": "Relay Two", "kind": "station", "placement": "orbital", "host_type": "planet",
        "host_id": planet.db_id}).status_code == 201
    assert admin.post(f"/api/planets/{planet.db_id}/regenerate").status_code == 200
    assert len(_rows(mysql_config, "SELECT id FROM facilities WHERE planet_id = ?", (planet.db_id,))) == 1


def test_regenerated_system_validates(mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_moons)
    system = _load(mysql_config, system_id)
    planet = next(p for p in system.planets if p.body_type != 'a')
    body, owner = adminEdits.find_body(system, "planet", planet.db_id)
    result = adminEdits.regenerate_planet(system, body, owner)
    assert result.removed == []
    spacing = [p for p in validation.check_star_system(system) if "too close" in p.message
               and "moon" not in p.message]
    assert spacing == []


def _saved_nebula(mysql_config):
    sector_id, _system_id = _saved_system(mysql_config)
    cfg = SystemConfig()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            nebula_id = store.insert_nebula(conn, Nebula(cfg, name="Test Veil"), sector_id=sector_id,
                                          placement={"center_x_pc": 1.0, "center_y_pc": 2.0, "center_z_pc": 3.0,
                                                     "galactic_radius_pc": 3.7})
    finally:
        conn.close()
    return sector_id, nebula_id


def test_regenerate_phenomenon_keeps_id_name_and_place(admin, mysql_config):
    _sector_id, nebula_id = _saved_nebula(mysql_config)
    before = _rows(mysql_config, "SELECT * FROM nebulae WHERE id = ?", (nebula_id,))[0]
    response = admin.post(f"/api/phenomena/nebula/{nebula_id}/regenerate")
    assert response.status_code == 200, response.get_json()
    rows = _rows(mysql_config, "SELECT * FROM nebulae")
    assert len(rows) == 1
    after = rows[0]
    for key in ("id", "name", "sector_id", "center_x_pc", "center_y_pc", "center_z_pc"):
        assert after[key] == before[key]


def test_delete_phenomenon(admin, mysql_config):
    _sector_id, nebula_id = _saved_nebula(mysql_config)
    assert admin.delete(f"/api/phenomena/nebula/{nebula_id}").status_code == 200
    assert _rows(mysql_config, "SELECT id FROM nebulae") == []
    assert admin.delete(f"/api/phenomena/nebula/{nebula_id}").status_code == 404
    assert admin.delete("/api/phenomena/teapot/1").status_code == 404


def test_delete_sector_with_contents(admin, mysql_config):
    sector_id, system_id = _saved_system(mysql_config)
    response = admin.delete(f"/api/sectors/{sector_id}/contents")
    assert response.status_code == 200, response.get_json()
    assert response.get_json()["systems"] == 1
    assert _rows(mysql_config, "SELECT id FROM star_systems WHERE id = ?", (system_id,)) == []
    assert _rows(mysql_config, "SELECT id FROM sectors WHERE id = ?", (sector_id,)) == []
    assert admin.delete(f"/api/sectors/{sector_id}/contents").status_code == 404


def test_regenerate_sector_off_the_grid_is_refused(admin, mysql_config):
    sector_id, system_id = _saved_system(mysql_config)
    assert admin.post(f"/api/sectors/{sector_id}/regenerate").status_code == 409
    assert _rows(mysql_config, "SELECT id FROM star_systems WHERE id = ?", (system_id,))


def test_regenerate_sector_is_queued(admin, mysql_config, monkeypatch):
    """A placed sector with a density plan is queued, not regenerated in
    the request: 202 with a job id, and nothing deleted yet."""
    from planetgen.queue import api_jobs
    sector_id, system_id = _saved_system(mysql_config)
    monkeypatch.setattr(editStore, "sector_address", lambda conn, sector: (1, 2, 3))
    monkeypatch.setattr(store, "get_galaxy_shape", lambda conn: object())
    queued = []
    monkeypatch.setattr(api_jobs, "submit", lambda function, *args: queued.append((function, args)) or "cd" * 8)
    response = admin.post(f"/api/sectors/{sector_id}/regenerate")
    assert response.status_code == 202, response.get_json()
    assert response.get_json()["job_id"] == "cd" * 8
    assert queued == [(api_jobs.regenerate_sector, (sector_id, queued[0][1][1]))]
    assert _rows(mysql_config, "SELECT id FROM star_systems WHERE id = ?", (system_id,))


# ---------------------------------------------------------------------
# The web pages' buttons
# ---------------------------------------------------------------------

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.web import csrf  # noqa: E402


@pytest.fixture
def web_app(mysql_config, monkeypatch):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    def no_http(*args, **kwargs):
        raise AssertionError("HTTP transport used inside a Flask request")
    monkeypatch.setattr(apiclient, "_http_transport", no_http)
    application = create_app(RealConfig)
    application.testing = True
    return application


def _web_admin(web_app, mysql_config):
    client = web_app.test_client()
    store.get_connection(mysql_config).close()
    _username, password = adminAuth.bootstrap_control_schema(mysql_config)
    assert client.post("/api/auth/login", json={"username": "admin", "password": password}).status_code == 200
    assert client.post("/api/auth/change-credentials", json={
        "current_password": password, "new_username": "boss", "new_password": "a-long-new-password-1",
    }).status_code == 200
    return client


def _token(app, client):
    nonce = "n" * 43
    client.set_cookie(csrf.COOKIE_NAME, nonce)
    session = client.get_cookie("pg_admin_session")
    with app.app_context():
        return csrf._sign(nonce, session.value if session else "")


def _edit(app, client, url, action, target, **extra):
    return client.post(url, data={csrf.FIELD_NAME: _token(app, client), "edit_action": action,
                                  "edit_target": target, **extra})


def test_edit_panel_shows_only_for_an_admin(web_app, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    assert 'id="admin-menu"' not in web_app.test_client().get(f"/system/{system_id}").get_data(as_text=True)
    client = _web_admin(web_app, mysql_config)
    html = client.get(f"/system/{system_id}").get_data(as_text=True)
    assert 'id="admin-menu"' in html
    assert f'value="system:{system_id}"' in html


def test_system_page_regenerates_a_planet(web_app, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type != 'a')
    # TEST.71: a drawn name can hold an apostrophe, which the page escapes;
    # give it one every time.
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE planets SET name = ? WHERE id = ?", ("Ilq'Ot", planet.db_id))
    finally:
        conn.close()
    # An asteroid belt's id comes from its own table and can equal the
    # planet's, so match the body type too.
    planet = next(p for p in _load(mysql_config, system_id).planets
                  if p.body_type != 'a' and p.db_id == planet.db_id)
    assert planet.name == "Ilq'Ot"
    client = _web_admin(web_app, mysql_config)
    response = _edit(web_app, client, f"/system/{system_id}", "regenerate", f"planet:{planet.db_id}")
    assert response.status_code == 303
    page = client.get(response.headers["Location"]).get_data(as_text=True)
    assert f"Regenerated {escape(planet.name)}" in page


def test_system_page_refuses_a_body_from_another_system(web_app, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    client = _web_admin(web_app, mysql_config)
    response = _edit(web_app, client, f"/system/{system_id}", "delete", "planet:999999")
    page = client.get(response.headers["Location"]).get_data(as_text=True)
    assert "isn&#39;t something on this page" in page or "isn't something on this page" in page


def test_deleting_the_system_goes_back_to_its_sector(web_app, mysql_config):
    sector_id, system_id = _saved_system(mysql_config)
    client = _web_admin(web_app, mysql_config)
    response = _edit(web_app, client, f"/system/{system_id}", "delete", f"system:{system_id}")
    assert response.status_code == 303
    assert response.headers["Location"].endswith(f"/sector/{sector_id}")
    assert "System deleted." in client.get(response.headers["Location"]).get_data(as_text=True)


def test_sector_page_deletes_the_sector(web_app, mysql_config):
    sector_id, _system_id = _saved_system(mysql_config)
    client = _web_admin(web_app, mysql_config)
    assert 'id="admin-menu"' in client.get(f"/sector/{sector_id}").get_data(as_text=True)
    response = _edit(web_app, client, f"/sector/{sector_id}", "delete", f"sector:{sector_id}")
    assert response.status_code == 303
    assert response.headers["Location"].endswith("/sectors")
    assert "Sector deleted" in client.get(response.headers["Location"]).get_data(as_text=True)


def test_phenomenon_page_regenerates_and_deletes(web_app, mysql_config):
    _sector_id, nebula_id = _saved_nebula(mysql_config)
    client = _web_admin(web_app, mysql_config)
    url = f"/phenomenon/nebula/{nebula_id}"
    assert 'id="admin-menu"' in client.get(url).get_data(as_text=True)
    response = _edit(web_app, client, url, "regenerate", f"nebula:{nebula_id}")
    assert "Regenerated Test Veil." in client.get(response.headers["Location"]).get_data(as_text=True)
    response = _edit(web_app, client, url, "delete", f"nebula:{nebula_id}")
    assert response.headers["Location"].endswith("/phenomena")


# ---------------------------------------------------------------------
# ADM.6: class changes; ADM.7: star changes
# ---------------------------------------------------------------------

def _first_planet(system):
    return next(p for p in system.planets if p.body_type != 'a')


def _with_planet(system):
    return any(p.body_type != 'a' for p in system.planets)


def _with_reclassable_planet(system):
    """A planet first, and another class that fits its mass (ADM.27): a
    rare mass can fit only its own class."""
    if not _with_planet(system):
        return False
    planet = _first_planet(system)
    return any(c != planet.planet_class and adminEdits.class_fits_mass(c, planet.mass)
               for c in tuning.PLANET_CLASSES)


def test_class_options_lists_recommended_classes(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_planet)
    system = _load(mysql_config, system_id)
    response = admin.get(f"/api/systems/{system_id}/class-options")
    assert response.status_code == 200
    body = response.get_json()
    planet = _first_planet(system)
    recommended = body["recommended"][f"planet:{planet.db_id}"]
    assert planet.planet_class not in recommended
    assert set(recommended) <= set(body["all"])


def test_recommended_class_change_saves_and_validates(admin, mysql_config):
    def has_option(system):
        return any(adminEdits.recommended_classes(system, p, system.planets)
                   for p in system.planets if p.body_type != 'a')
    _sector_id, system_id = _saved_system(mysql_config, has_option)
    system = _load(mysql_config, system_id)
    planet = next(p for p in system.planets
                  if p.body_type != 'a' and adminEdits.recommended_classes(system, p, system.planets))
    new_class = adminEdits.recommended_classes(system, planet, system.planets)[0]
    response = admin.post(f"/api/planets/{planet.db_id}/class", json={"class": new_class})
    assert response.status_code == 200, response.get_json()
    row = _rows(mysql_config, "SELECT name, planet_class FROM planets WHERE id = ?", (planet.db_id,))[0]
    assert row["planet_class"] == new_class
    assert row["name"] == planet.name


def test_unrecommended_class_needs_force(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_planet)
    system = _load(mysql_config, system_id)
    planet = _first_planet(system)
    recommended = adminEdits.recommended_classes(system, planet, system.planets)
    other = next((c for c in sorted(tuning.PLANET_CLASSES)
                  if c not in recommended and c != planet.planet_class
                  and adminEdits.class_fits_mass(c, planet.mass)), None)
    if other is None:
        pytest.skip("every class this planet's mass fits is recommended")
    assert admin.post(f"/api/planets/{planet.db_id}/class", json={"class": other}).status_code == 409
    assert admin.post(f"/api/planets/{planet.db_id}/class", json={"class": "?"}).status_code == 400
    response = admin.post(f"/api/planets/{planet.db_id}/class", json={"class": other, "force": True})
    assert response.status_code == 200, response.get_json()
    assert _rows(mysql_config, "SELECT planet_class FROM planets WHERE id = ?",
                 (planet.db_id,))[0]["planet_class"] == other


def _unsaved_system_with_planet():
    for _ in range(200):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MOONS = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(cfg)
        if any(p.body_type == 't' for p in system.planets):
            return system
    pytest.skip("no rocky planet generated")


def test_class_change_regenerates_surface_conditions_keeping_orbit_mass_and_name():
    """ADM.27: the body is re-generated as the new class -- composition,
    atmosphere, temperature, pressure, life -- at its orbit, mass and name."""
    # A fixed seed, and a rocky planet that has another class its mass
    # fits: the test no longer depends on what ran before it (TEST.97).
    draw.set_run_seed(97)
    for _ in range(50):
        system = _unsaved_system_with_planet()
        planet = next(p for p in system.planets if p.body_type == 't')
        new_class = next((c for c in sorted(tuning.PLANET_CLASSES)
                          if c != planet.planet_class and adminEdits.class_fits_mass(c, planet.mass)), None)
        if new_class is not None:
            break
    assert new_class is not None, "no rocky planet with a second class its mass fits in 50 systems"
    name, mass, distance = planet.name, planet.mass, planet.distance
    adminEdits.change_class(system, planet, system.planets, new_class, force=True)
    data = tuning.PLANET_CLASSES[new_class]
    assert planet.planet_class == new_class
    assert planet.name == name
    assert planet.mass == pytest.approx(mass, rel=1e-9)
    assert planet.distance == pytest.approx(distance, rel=1e-6)
    assert planet.composition == data["composition"]
    assert planet.atmosphere == (data["atmosphere"] or "None")
    low, high = data["radius_range"]
    assert low <= planet.radius <= high
    volume_m3 = 4 / 3 * 3.141592653589793 * (planet.radius * 1000) ** 3
    assert planet.density == pytest.approx(mass / volume_m3 / 1000, rel=1e-6)
    assert planet.surface_temperature is not None


def test_class_change_the_mass_cant_fit_is_refused_and_changes_nothing():
    """ADM.27: a gas giant class for a rocky world is refused, even forced,
    with a message, and the body is left as it was."""
    system = _unsaved_system_with_planet()
    planet = next(p for p in system.planets if p.body_type == 't')
    assert not adminEdits.class_fits_mass("J", planet.mass)
    before = dict(vars(planet))
    with pytest.raises(ValueError, match="Earth masses; class J needs"):
        adminEdits.change_class(system, planet, system.planets, "J", force=True)
    assert "J" not in adminEdits.recommended_classes(system, planet, system.planets)
    after = vars(planet)
    assert {k: v for k, v in after.items() if not isinstance(v, list)} == \
        {k: v for k, v in before.items() if not isinstance(v, list)}


def test_change_star_keeps_classes_and_saves_the_type(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_planet)
    before = _load(mysql_config, system_id)
    response = admin.post(f"/api/systems/{system_id}/star", json={"star_type": "K2V", "drop_facilities": True})
    assert response.status_code == 200, response.get_json()
    body = response.get_json()
    assert body["star_type"] == "K2V"
    after = _load(mysql_config, system_id)
    assert after.star.type.split()[0] == "K2V"
    assert after.star.name == before.star.name
    kept = {p.db_id: p.planet_class for p in after.planets if p.body_type != 'a'}
    for planet in before.planets:
        if planet.body_type != 'a' and planet.db_id in kept:
            assert kept[planet.db_id] == planet.planet_class
    assert not validation.check_star_system(after) or body["warnings"]


def test_change_star_refuses_bad_types(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    assert admin.post(f"/api/systems/{system_id}/star", json={"star_type": "G10V"}).status_code == 400
    assert admin.post(f"/api/systems/{system_id}/star", json={"star": "G2V"}).status_code == 400
    assert admin.post("/api/systems/999999/star", json={"star_type": "G2V"}).status_code == 404


def test_change_star_refuses_a_binary(mysql_config):
    for _ in range(200):
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = True
        system = StarSystem(cfg)
        if system.binary_type is not None:
            break
    else:
        pytest.skip("no binary generated")
    with pytest.raises(ValueError, match="single star"):
        adminEdits.change_star(system, "K2V")


def test_system_page_changes_a_class(web_app, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_reclassable_planet)
    planet = _first_planet(_load(mysql_config, system_id))
    client = _web_admin(web_app, mysql_config)
    html = client.get(f"/system/{system_id}").get_data(as_text=True)
    assert "Change class" in html and "Change star" in html
    other = next(c for c in sorted(tuning.PLANET_CLASSES)
                 if c != planet.planet_class and adminEdits.class_fits_mass(c, planet.mass))
    response = _edit(web_app, client, f"/system/{system_id}", "class", f"planet:{planet.db_id}",
                     planet_class=f"force:{other}")
    assert response.status_code == 303
    page = client.get(response.headers["Location"]).get_data(as_text=True)
    assert f"to class {other}" in page
    response = _edit(web_app, client, f"/system/{system_id}", "class", f"system:{system_id}", planet_class="J")
    page = client.get(response.headers["Location"]).get_data(as_text=True)
    assert "something on this page" in page


def test_system_page_changes_the_star(web_app, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    client = _web_admin(web_app, mysql_config)
    response = _edit(web_app, client, f"/system/{system_id}", "star", f"system:{system_id}", star_type="k2v",
                     drop_facilities="1")
    assert response.status_code == 303
    client.get(response.headers["Location"])
    assert _load(mysql_config, system_id).star.type.split()[0] == "K2V"
