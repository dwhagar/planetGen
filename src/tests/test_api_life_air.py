# tests/test_api_life_air.py

"""UX.237: `GET /api/planets/<id>/life-air` and `/api/moons/<id>/life-air`, the data of the Life and air chart, worked out on request from the stored body."""

import copy

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import biogenic_oxygen as bo
from planetgen.physics import life_air
from tests.publicids import pid
from tests.test_api import client  # noqa: F401


def _living_planet_system(mysql_config, monkeypatch):
    """Saves a system whose first airy rocky planet has Class 1a life that made methane, then oxygen; returns `(system_id, planet row id)`."""
    monkeypatch.setattr(life_air, "methane_flux_mol_yr", lambda body, life_class, supplies=None: 1.0e12)
    monkeypatch.setattr(bo, "profile_of", lambda body, *args: bo.Profile(1.0, 0.0, 0.3e9, 1.0, 1.0, 0.3))
    for _ in range(40):
        config = SystemConfig()
        config.PLANETS = True
        config.MAX_PLANETS = True
        config.BINARY_SYSTEM = False
        system = StarSystem(system_config=config)
        for planet in system.planets:
            if planet.body_type == "t" and (planet.atmospheric_pressure or 0) > 20000.0 and planet.greenhouse_factor is not None:
                bo.store_air(planet, {"N2": 70000.0, "CO2": 30000.0})
                system.star.age = 3.0
                planet.life_class, planet.life_class_variant, planet.life_stage = 1, "a", "photosynthesis"
                planet.life_age_photosynthesis_gyr, planet.life_age_abiogenesis_gyr = 0.5, 0.2
                assert bo.apply(planet, 2.0e19)
                system_id = store.save_system(system, config, config=mysql_config)
                conn = store.get_connection(mysql_config)
                try:
                    rows = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND life_class = 1", (system_id,)).fetchall()
                finally:
                    conn.close()
                if rows:
                    return system_id, rows[0]["id"]
    pytest.skip("no airy rocky planet in these draws")


def test_a_living_planets_air_history_has_the_steps_and_stages(client, mysql_config, monkeypatch):
    _system_id, planet_id = _living_planet_system(mysql_config, monkeypatch)
    response = client.get(f"/api/planets/{pid('planet', planet_id)}/life-air")
    assert response.status_code == 200, response.get_data(as_text=True)
    body = response.get_json()
    assert body["id"] == pid("planet", planet_id) and body["life_class_id"] == "1a"
    history = body["life_air"]
    assert history["age_gyr"] == pytest.approx(3.0) and history["producers_gyr"] == pytest.approx(0.5)
    assert history["rose_gyr"] == pytest.approx(0.8)
    gases = {gas["formula"]: gas["steps"] for gas in history["gases"]}
    assert set(gases) == {"O2", "CO2", "CH4"}
    assert gases["O2"][0][1] == 0.0 and gases["O2"][-1][1] > 0.0 and gases["CH4"][1][1] > 0.0 and gases["CH4"][-1][1] == pytest.approx(0.0, abs=1e-6)
    keys = [stage["key"] for stage in history["stages"]]
    assert keys[:2] == ["abiogenesis", "photosynthesis"] and all(stage["age_gyr"] <= 3.0 for stage in history["stages"])
    assert history["speculative"] is None


def test_a_planet_without_a_living_class_has_no_history(client, mysql_config):
    config = SystemConfig()
    config.PLANETS = True
    config.MAX_PLANETS = True
    system_id = store.save_system(StarSystem(system_config=config), config, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        planet_id = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND life_stage IS NULL LIMIT 1", (system_id,)).fetchone()
    finally:
        conn.close()
    if planet_id is None:
        pytest.skip("every planet drawn has a life stage")
    body = client.get(f"/api/planets/{pid('planet', planet_id['id'])}/life-air").get_json()
    assert body["life_class_id"] is None and body["life_air"] is None


def test_an_unknown_planet_or_moon_is_a_404(client):
    assert client.get(f"/api/planets/{pid('planet', 999999)}/life-air").status_code == 404
    assert client.get(f"/api/moons/{pid('moon', 999999)}/life-air").status_code == 404
