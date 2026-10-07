# tests/test_phenomenon_composition.py

"""
DB.2: an asteroid field's and an interstellar comet's composition rows
(`asteroid_field_composition`, `interstellar_comet_composition`) are read
back: `GET /api/phenomena/<type>/<id>` returns them as `composition`, and
the phenomenon page shows them. The parent row's `composition_summary` is
overwritten after saving, so the page can only show the right text by
reading the rows.
"""

import pytest

from stellarObjects import _db
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import InterstellarComet
from tests.test_admin_edits import web_app  # noqa: F401
from tests.test_api import client  # noqa: F401

_FIELD_COMPOSITION = [("iron", "high"), ("nickel", "moderate"), ("platinum", "trace")]
_COMET_COMPOSITION = ["water ice", "carbon monoxide ice"]


def _save(mysql_config, kind):
    cfg = SystemConfig()
    if kind == "asteroid-field":
        obj, table = AsteroidField(cfg), "asteroid_fields"
        obj.composition = list(_FIELD_COMPOSITION)
    else:
        obj, table = InterstellarComet(cfg), "interstellar_comets"
        obj.composition = list(_COMET_COMPOSITION)
    phenomenon_id = _db.save_phenomenon(obj, cfg, kind, config=mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        conn.execute(f"UPDATE {table} SET composition_summary = 'stale summary' WHERE id = ?", (phenomenon_id,))
        conn.commit()
    finally:
        conn.close()
    return phenomenon_id


def test_api_returns_asteroid_field_composition_rows(client, mysql_config):  # noqa: F811
    field_id = _save(mysql_config, "asteroid-field")
    body = client.get(f"/api/phenomena/asteroid_field/{field_id}").get_json()
    assert body["composition"] == [{"component": c, "concentration": n} for c, n in _FIELD_COMPOSITION]


def test_api_returns_interstellar_comet_composition_rows(client, mysql_config):  # noqa: F811
    comet_id = _save(mysql_config, "comet")
    body = client.get(f"/api/phenomena/interstellar_comet/{comet_id}").get_json()
    assert body["composition"] == _COMET_COMPOSITION


@pytest.mark.parametrize("kind, page_type, expected", [
    ("asteroid-field", "asteroid_field",
     "high concentrations of iron, moderate concentrations of nickel, and trace amounts of platinum"),
    ("comet", "interstellar_comet", "water ice and carbon monoxide ice"),
])
def test_phenomenon_page_shows_the_composition_rows(web_app, mysql_config, kind, page_type, expected):  # noqa: F811
    phenomenon_id = _save(mysql_config, kind)
    page = web_app.test_client().get(f"/phenomenon/{page_type}/{phenomenon_id}").get_data(as_text=True)
    assert expected in page
    assert "stale summary" not in page
