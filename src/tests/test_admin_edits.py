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
from planetgen.physics import composition_text
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from tests.publicids import pid


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
    response = admin.post(f"/api/planets/{pid('planet', planet.db_id)}/regenerate")
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
    response = admin.delete(f"/api/planets/{pid('planet', planet.db_id)}")
    assert response.status_code == 200, response.get_json()
    assert _rows(mysql_config, "SELECT id FROM planets WHERE id = ?", (planet.db_id,)) == []
    assert _rows(mysql_config, "SELECT id FROM moons WHERE planet_id = ?", (planet.db_id,)) == []
    assert admin.delete(f"/api/planets/{pid('planet', planet.db_id)}").status_code == 404


def test_regenerate_and_delete_a_moon(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_moons)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type != 'a' and p.moons)
    moon = planet.moons[0]
    response = admin.post(f"/api/moons/{pid('moon', moon.db_id)}/regenerate")
    assert response.status_code in (200, 409), response.get_json()
    if response.status_code == 200:
        assert _rows(mysql_config, "SELECT name FROM moons WHERE id = ?", (moon.db_id,))[0]["name"] == moon.name
    assert admin.delete(f"/api/moons/{pid('moon', moon.db_id)}").status_code == 200
    assert _rows(mysql_config, "SELECT id FROM moons WHERE id = ?", (moon.db_id,)) == []


def test_regenerate_and_delete_a_belt(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, lambda s: any(p.body_type == 'a' for p in s.planets))
    belt = next(p for p in _load(mysql_config, system_id).planets if p.body_type == 'a')
    response = admin.post(f"/api/belts/{pid('belt', belt.db_id)}/regenerate")
    assert response.status_code == 200, response.get_json()
    row = _rows(mysql_config, "SELECT lower_limit_km FROM asteroid_belts WHERE id = ?", (belt.db_id,))[0]
    assert row["lower_limit_km"] == pytest.approx(belt.lower_limit * 1.495978707e8, rel=1e-6)
    assert admin.delete(f"/api/belts/{pid('belt', belt.db_id)}").status_code == 200
    assert _rows(mysql_config, "SELECT id FROM asteroid_belts WHERE id = ?", (belt.db_id,)) == []


def test_delete_refuses_to_drop_facilities_unless_told(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type != 'a')
    response = admin.post("/api/facilities", json={
        "name": "Relay One", "kind": "station", "placement": "orbital", "host_type": "planet",
        "host_id": pid("planet", planet.db_id)})
    assert response.status_code == 201, response.get_json()
    assert admin.delete(f"/api/planets/{pid('planet', planet.db_id)}").status_code == 409
    assert admin.delete(f"/api/planets/{pid('planet', planet.db_id)}", json={"drop_facilities": True}).status_code == 200


def test_regenerate_planet_keeps_facilities_on_the_planet(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    planet = next(p for p in _load(mysql_config, system_id).planets if p.body_type != 'a')
    assert admin.post("/api/facilities", json={
        "name": "Relay Two", "kind": "station", "placement": "orbital", "host_type": "planet",
        "host_id": pid("planet", planet.db_id)}).status_code == 201
    assert admin.post(f"/api/planets/{pid('planet', planet.db_id)}/regenerate").status_code == 200
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
            store.assign_uids(conn, phenomenon=("nebulae", nebula_id))
    finally:
        conn.close()
    return sector_id, nebula_id


def test_regenerate_phenomenon_keeps_id_name_and_place(admin, mysql_config):
    _sector_id, nebula_id = _saved_nebula(mysql_config)
    before = _rows(mysql_config, "SELECT * FROM nebulae WHERE id = ?", (nebula_id,))[0]
    response = admin.post(f"/api/phenomena/nebula/{pid('nebula', nebula_id, mysql_config)}/regenerate")
    assert response.status_code == 200, response.get_json()
    rows = _rows(mysql_config, "SELECT * FROM nebulae")
    assert len(rows) == 1
    after = rows[0]
    assert before["uid"] is not None      # GEN.175: regenerating kept a NULL uid before, not the old one
    for key in ("id", "uid", "name", "sector_id", "center_x_pc", "center_y_pc", "center_z_pc"):
        assert after[key] == before[key]


def test_delete_phenomenon(admin, mysql_config):
    _sector_id, nebula_id = _saved_nebula(mysql_config)
    printed = pid("nebula", nebula_id, mysql_config)
    assert admin.delete(f"/api/phenomena/nebula/{printed}").status_code == 200
    assert _rows(mysql_config, "SELECT id FROM nebulae") == []
    assert admin.delete(f"/api/phenomena/nebula/{printed}").status_code == 404
    assert admin.delete("/api/phenomena/teapot/1").status_code == 404


def test_delete_sector_with_contents(admin, mysql_config):
    sector_id, system_id = _saved_system(mysql_config)
    response = admin.delete(f"/api/sectors/{pid('sector', sector_id)}/contents")
    assert response.status_code == 200, response.get_json()
    assert response.get_json()["systems"] == 1
    assert _rows(mysql_config, "SELECT id FROM star_systems WHERE id = ?", (system_id,)) == []
    assert _rows(mysql_config, "SELECT id FROM sectors WHERE id = ?", (sector_id,)) == []
    assert admin.delete(f"/api/sectors/{pid('sector', sector_id)}/contents").status_code == 404


def test_regenerate_sector_off_the_grid_is_refused(admin, mysql_config):
    sector_id, system_id = _saved_system(mysql_config)
    assert admin.post(f"/api/sectors/{pid('sector', sector_id)}/regenerate").status_code == 409
    assert _rows(mysql_config, "SELECT id FROM star_systems WHERE id = ?", (system_id,))


def test_regenerate_sector_is_queued(admin, mysql_config, monkeypatch):
    """A placed sector with a density plan is queued, not regenerated in
    the request: 202 with a job id, and nothing deleted yet."""
    from planetgen.queue import api_jobs
    sector_id, system_id = _saved_system(mysql_config)
    monkeypatch.setattr(editStore, "sector_address", lambda conn, sector: (1, 2, 3))
    monkeypatch.setattr(store, "get_galaxy_shape", lambda conn: object())
    queued = []
    monkeypatch.setattr(api_jobs, "submit", lambda function, *args, **options: queued.append((function, args)) or "cd" * 8)
    response = admin.post(f"/api/sectors/{pid('sector', sector_id)}/regenerate")
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
    assert 'id="admin-menu"' not in web_app.test_client().get(f"/system/{pid('system', system_id)}").get_data(as_text=True)
    client = _web_admin(web_app, mysql_config)
    html = client.get(f"/system/{pid('system', system_id)}").get_data(as_text=True)
    assert 'id="admin-menu"' in html
    assert f'value="system:{pid("system", system_id)}"' in html


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
    # UX.173: a body's Admin menu is on its own page, which takes the POST; the system page only takes the system's.
    response = _edit(web_app, client, f"/object/planet:{pid('planet', planet.db_id)}", "regenerate",
                     f"planet:{pid('planet', planet.db_id)}")
    assert response.status_code == 303
    page = client.get(response.headers["Location"]).get_data(as_text=True)
    assert f"Regenerated {escape(planet.name)}" in page


def test_system_page_refuses_a_body_from_another_system(web_app, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    client = _web_admin(web_app, mysql_config)
    response = _edit(web_app, client, f"/system/{pid('system', system_id)}", "delete", "planet:999999")
    page = client.get(response.headers["Location"]).get_data(as_text=True)
    assert "isn&#39;t something on this page" in page or "isn't something on this page" in page


def test_deleting_the_system_goes_back_to_its_sector(web_app, mysql_config):
    sector_id, system_id = _saved_system(mysql_config)
    client = _web_admin(web_app, mysql_config)
    response = _edit(web_app, client, f"/system/{pid('system', system_id)}", "delete", f"system:{pid('system', system_id)}")
    assert response.status_code == 303
    assert response.headers["Location"].endswith(f"/sector/{pid('sector', sector_id)}")
    assert "System deleted." in client.get(response.headers["Location"]).get_data(as_text=True)


def test_sector_page_deletes_the_sector(web_app, mysql_config):
    sector_id, _system_id = _saved_system(mysql_config)
    client = _web_admin(web_app, mysql_config)
    assert 'id="admin-menu"' in client.get(f"/sector/{pid('sector', sector_id)}").get_data(as_text=True)
    response = _edit(web_app, client, f"/sector/{pid('sector', sector_id)}", "delete", f"sector:{pid('sector', sector_id)}")
    assert response.status_code == 303
    assert response.headers["Location"].endswith("/sectors")
    assert "Sector deleted" in client.get(response.headers["Location"]).get_data(as_text=True)


def test_phenomenon_page_regenerates_and_deletes(web_app, mysql_config):
    _sector_id, nebula_id = _saved_nebula(mysql_config)
    client = _web_admin(web_app, mysql_config)
    nebula_pid = pid("nebula", nebula_id, mysql_config)
    url = f"/phenomenon/nebula/{nebula_pid}"
    assert 'id="admin-menu"' in client.get(url).get_data(as_text=True)
    response = _edit(web_app, client, url, "regenerate", f"nebula:{nebula_pid}")
    assert "Regenerated Test Veil." in client.get(response.headers["Location"]).get_data(as_text=True)
    response = _edit(web_app, client, url, "delete", f"nebula:{nebula_pid}")
    assert response.headers["Location"].endswith("/phenomena")


# ---------------------------------------------------------------------
# ADM.6: class changes; ADM.7: star changes
# ---------------------------------------------------------------------

def _first_planet(system):
    return next(p for p in system.planets if p.body_type != 'a')


def _with_planet(system):
    return any(p.body_type != 'a' for p in system.planets)


def test_change_star_derives_the_classes_again_and_saves_the_type(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config, _with_planet)
    before = _load(mysql_config, system_id)
    response = admin.post(f"/api/systems/{pid('system', system_id)}/star", json={"star_type": "K2V", "drop_facilities": True})
    assert response.status_code == 200, response.get_json()
    body = response.get_json()
    assert body["star_type"] == "K2V"
    after = _load(mysql_config, system_id)
    assert after.star.type.split()[0] == "K2V"
    assert after.star.name == before.star.name
    for planet in after.planets:  # GEN.646: the class follows the new light (made again by the edit), so it is not pinned to the old one
        if planet.body_type != 'a':
            assert planet.planet_class in tuning.PLANET_CLASSES
    assert not validation.check_star_system(after) or body["warnings"]


def test_change_star_refuses_bad_types(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    assert admin.post(f"/api/systems/{pid('system', system_id)}/star", json={"star_type": "G10V"}).status_code == 400
    assert admin.post(f"/api/systems/{pid('system', system_id)}/star", json={"star": "G2V"}).status_code == 400
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


def test_there_is_no_class_change_on_a_body_page(web_app, mysql_config):
    """GEN.646: a body's class is derived from its conditions, so the Admin menu has no Change class."""
    _sector_id, system_id = _saved_system(mysql_config, _with_planet)
    planet = _first_planet(_load(mysql_config, system_id))
    client = _web_admin(web_app, mysql_config)
    page = client.get(f"/object/planet:{pid('planet', planet.db_id)}").get_data(as_text=True)
    assert "Change class" not in page
    assert "Regenerate" in page


def test_system_page_changes_the_star(web_app, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    client = _web_admin(web_app, mysql_config)
    response = _edit(web_app, client, f"/system/{pid('system', system_id)}", "star", f"system:{pid('system', system_id)}", star_type="k2v",
                     drop_facilities="1")
    assert response.status_code == 303
    client.get(response.headers["Location"])
    assert _load(mysql_config, system_id).star.type.split()[0] == "K2V"


def _planet_uids(mysql_config, system_id):
    return [row["uid"] for row in _rows(mysql_config, "SELECT uid FROM planets WHERE star_system_id = ? ORDER BY id",
                                       (system_id,))]


def test_a_body_added_after_a_delete_gets_its_own_uid(mysql_config):
    """GEN.173 and GEN.174: delete a planet, then add one: the new row is saved with an ID (it was NULL), and
    not one another planet already has (it ranked among the current rows and took a used one)."""
    import copy
    _sector_id, system_id = _saved_system(mysql_config, lambda s: sum(p.body_type != 'a' for p in s.planets) >= 3)
    conn = store.get_connection(mysql_config)
    try:
        system = store.load_star_system(conn, system_id)
        planets = [p for p in system.planets if p.body_type != 'a']
        system.planets.remove(planets[0])                       # delete one, not the last
        with conn:
            editStore.save_system_edits(conn, system_id, system)
        system = store.load_star_system(conn, system_id)
        added = copy.deepcopy([p for p in system.planets if p.body_type != 'a'][-1])
        added.db_id = None
        added.name = "Added Planet"
        for moon in added.moons:
            moon.db_id = None
        system.planets.append(added)
        with conn:
            editStore.save_system_edits(conn, system_id, system)
    finally:
        conn.close()
    uids = _planet_uids(mysql_config, system_id)
    assert all(uid is not None for uid in uids), uids
    assert len(set(uids)) == len(uids), uids
