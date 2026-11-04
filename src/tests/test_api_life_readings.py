# tests/test_api_life_readings.py

"""UX.251: `GET /api/planets/<id>/life-readings` and `/api/moons/<id>/life-readings`, the three words for each kind of life, worked out on request."""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from tests.publicids import pid
from tests.test_api import client  # noqa: F401


def _save_with_life(mysql_config, scores):
    """Saves a system whose first rocky planet has Class 1a life and the kind `scores`; returns the planet row id."""
    for _ in range(40):
        config = SystemConfig()
        config.PLANETS = True
        config.MAX_PLANETS = True
        config.BINARY_SYSTEM = False
        system = StarSystem(system_config=config)
        for planet in system.planets:
            if planet.body_type == "t" and planet.surface_temperature is not None:
                planet.life_class, planet.life_class_variant, planet.life_stage = 1, "a", "multicellularity"
                planet.life_microbial, planet.life_flora, planet.life_fauna, planet.life_sentient = scores
                system_id = store.save_system(system, config, config=mysql_config)
                conn = store.get_connection(mysql_config)
                try:
                    rows = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND life_class = 1", (system_id,)).fetchall()
                finally:
                    conn.close()
                if rows:
                    return rows[0]["id"]
    pytest.skip("no rocky planet in these draws")


def test_a_living_planet_has_the_four_kinds_with_their_words(client, mysql_config):
    planet_id = _save_with_life(mysql_config, (9, 5, 3, 0))
    response = client.get(f"/api/planets/{pid('planet', planet_id)}/life-readings")
    assert response.status_code == 200, response.get_data(as_text=True)
    body = response.get_json()
    assert body["id"] == pid("planet", planet_id) and body["life_class_id"] == "1a"
    readings = body["readings"]
    assert [r["kind"] for r in readings] == ["microbial", "flora", "fauna", "sentient"]
    assert [r["score"] for r in readings] == [9, 5, 3, 0]
    assert [r["abundance"] for r in readings] == ["Dominant", "Frequent", "Occasional", "Absent"]
    assert all(r["trend"] in ("Growing", "Stable", "Declining") and r["security"] in ("Secure", "Marginal", "Precarious", "Dying out")
               for r in readings[:3])
    assert readings[3]["trend"] is None and readings[3]["security"] is None
    assert readings[0]["text"].startswith("Microbes: dominant, ") and readings[3]["text"] == "Sentient life: absent."


def test_an_unknown_planet_or_moon_is_a_404(client):
    assert client.get(f"/api/planets/{pid('planet', 999999)}/life-readings").status_code == 404
    assert client.get(f"/api/moons/{pid('moon', 999999)}/life-readings").status_code == 404
