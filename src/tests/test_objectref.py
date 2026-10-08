# tests/test_objectref.py

"""
NAV.7: the one reference form for every object (`<kind>:<id>`), its Python
and JavaScript helpers, `queryDb.resolve_object` and `GET /api/objects/<ref>`.
"""

import json
import os
import shutil
import subprocess

import pytest

from planetgen.db import store
from planetgen.galaxy import objectref
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula

# Imported fixtures (see test_bughunt_api_gaps.py for why this works).
from tests.test_api import (  # noqa: F401
    _place_sector, _save_wide_binary_with_moons, admin_client, client, default_admin_client,
    first_admin_password, seeded_sector,
)

_STATIC = os.path.join(os.path.dirname(__file__), "..", "html", "static")

CASES = [
    ("system:12", ("system", 12)),
    ("  planet:3 ", ("planet", 3)),
    ("12", ("system", 12)),
    ("black_hole:0", ("black_hole", 0)),
    ("interstellar_comet:9", ("interstellar_comet", 9)),
]
BAD = ["", "galaxy", "planet", "planet:", "planet:-1", "planet:x", "ship:3", "moon:1:2", "system:1.5"]


@pytest.mark.parametrize("raw, expected", CASES)
def test_parse_reads_every_kind(raw, expected):
    assert objectref.parse(raw) == expected


@pytest.mark.parametrize("raw", BAD)
def test_parse_refuses_what_is_not_a_reference(raw):
    with pytest.raises(ValueError):
        objectref.parse(raw)


def test_format_round_trips_every_kind():
    for kind in objectref.KINDS:
        assert objectref.parse(objectref.format(kind, 7)) == (kind, 7)
    with pytest.raises(ValueError):
        objectref.format("ship", 1)
    with pytest.raises(ValueError):
        objectref.format("moon", -1)


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_browser_copy_matches_python():
    script = (
        "const m = await import(process.argv[1]);"
        "const raws = JSON.parse(process.argv[2]);"
        "console.log(JSON.stringify({parsed: raws.map(r => m.parse(r)), kinds: m.KINDS,"
        " formatted: m.KINDS.map(k => m.format(k, 7))}));"
    )
    module_url = "file://" + os.path.abspath(os.path.join(_STATIC, "objectref.js"))
    raws = [raw for raw, _ in CASES] + BAD
    out = json.loads(subprocess.run(
        ["node", "--input-type=module", "-e", script, module_url, json.dumps(raws)],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert out["kinds"] == list(objectref.KINDS)
    assert out["formatted"] == [objectref.format(kind, 7) for kind in objectref.KINDS]
    expected = [{"kind": kind, "id": object_id} for _, (kind, object_id) in CASES] + [None] * len(BAD)
    assert out["parsed"] == expected


def _ids(mysql_config, sql, params=()):
    conn = store.get_connection(mysql_config)
    try:
        return [row["id"] for row in conn.execute(sql, params).fetchall()]
    finally:
        conn.close()


def test_resolve_a_system_and_its_sector(client, seeded_sector):
    config, sector_id, system_ids = seeded_sector
    body = client.get(f"/api/objects/system:{system_ids[0]}").get_json()
    assert body["ref"] == f"system:{system_ids[0]}" and body["kind"] == "system"
    assert [p["ref"] for p in body["parents"]] == ["galaxy", f"sector:{sector_id}"]
    assert body["siblings"] == [f"system:{system_ids[1]}"]
    assert body["positions"]["sector_ly"] is not None
    assert body["positions"]["system_km"] == [0.0, 0.0, 0.0]
    # A bare number is a system.
    assert client.get(f"/api/objects/{system_ids[0]}").get_json() == body

    sector = client.get(f"/api/objects/sector:{sector_id}").get_json()
    assert sector["parents"] == [{"ref": "galaxy", "kind": "galaxy", "name": "Galaxy"}]
    assert sector["positions"]["sector_ly"] == [0.0, 0.0, 0.0]


def test_a_placed_sector_gives_galaxy_positions(client, mysql_config):
    sector_id = _place_sector(mysql_config, "Placed", (100.0, 200.0, 5.0))
    body = client.get(f"/api/objects/sector:{sector_id}").get_json()
    assert body["positions"]["galaxy_pc"] == [100.0, 200.0, 5.0]
    system_id = _ids(mysql_config, "SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,))[0]
    system = client.get(f"/api/objects/system:{system_id}").get_json()
    assert system["positions"]["galaxy_pc"] == pytest.approx([100.0, 200.0, 5.0])


def test_resolve_bodies_down_to_a_moon(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    moon_id = _ids(mysql_config, "SELECT id FROM moons WHERE star_system_id = ? ORDER BY id", (system_id,))[0]
    body = client.get(f"/api/objects/moon:{moon_id}").get_json()
    kinds = [p["kind"] for p in body["parents"]]
    assert kinds[0] == "galaxy" and kinds[-2:] == ["system", "planet"]
    assert body["parents"][-2]["ref"] == f"system:{system_id}"
    planet = client.get("/api/objects/" + body["parents"][-1]["ref"]).get_json()
    assert planet["kind"] == "planet" and planet["positions"]["system_km"] is not None
    assert body["ref"] not in body["siblings"]
    assert len(body["positions"]["system_km"]) == 3

    for star_id in _ids(mysql_config, "SELECT id FROM stars WHERE star_system_id = ?", (system_id,)):
        star = client.get(f"/api/objects/star:{star_id}").get_json()
        assert star["parents"][-1]["ref"] == f"system:{system_id}"
        assert len(star["positions"]["system_km"]) == 3
    # The two stars of a wide pair are apart, so they are each other's siblings.
    assert len(star["siblings"]) == 1


def test_resolve_refusals_are_json_404s(client, seeded_sector):
    for ref in ("moon:999999999", "ship:1", "planet", "system:99999999", "nebula:99999999"):
        response = client.get(f"/api/objects/{ref}")
        assert response.status_code == 404, ref
        assert isinstance(response.get_json()["error"], str)


def test_resolve_a_phenomenon(client, mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        nebula_id = store.insert_nebula(conn, Nebula(SystemConfig()), sector_id=None, placement={
            "center_x_pc": 5.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 5.5,
        })
        conn.commit()
    finally:
        conn.close()
    body = client.get(f"/api/objects/nebula:{nebula_id}").get_json()
    assert body["kind"] == "nebula" and body["parents"][0]["ref"] == "galaxy"
    assert body["positions"] == {"galaxy_pc": [5.0, 1.0, 2.0], "sector_ly": None, "system_km": None}
