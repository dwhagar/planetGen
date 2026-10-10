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
from tests.publicids import pid
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

PUBLIC_CASES = [
    ("system:FE81000A2B-0000005-000", ("system", "FE81000A2B-0000005-000")),
    ("  moon:fe81000a2b-0000005-003 ", ("moon", "FE81000A2B-0000005-003")),
    ("FE81000A2B-0000005-000", ("system", "FE81000A2B-0000005-000")),
    ("sector:100000000", ("sector", "100000000")),
    ("planet:3", ("planet", "3")),
]
PUBLIC_BAD = BAD + ["12", "100000000", "system:", "sector:xyz", "system:1-2-3-4"]


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


@pytest.mark.parametrize("raw, expected", PUBLIC_CASES)
def test_parse_public_keeps_the_id_as_text(raw, expected):
    assert objectref.parse_public(raw) == expected


@pytest.mark.parametrize("raw", PUBLIC_BAD)
def test_parse_public_refuses_what_is_not_a_public_reference(raw):
    with pytest.raises(ValueError):
        objectref.parse_public(raw)


def test_format_public_round_trips_every_kind():
    for kind in objectref.KINDS:
        assert objectref.parse_public(objectref.format_public(kind, "fe81000a2b-0000005-003")) == (
            kind, "FE81000A2B-0000005-003")
    with pytest.raises(ValueError):
        objectref.format_public("ship", "1")
    with pytest.raises(ValueError):
        objectref.format_public("moon", "-1")


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_browser_copy_matches_python():
    script = (
        "const m = await import(process.argv[1]);"
        "const raws = JSON.parse(process.argv[2]);"
        "console.log(JSON.stringify({parsed: raws.map(r => m.parse(r)), kinds: m.KINDS,"
        " formatted: m.KINDS.map(k => m.format(k, 'fe81000a2b-0000005-003'))}));"
    )
    module_url = "file://" + os.path.abspath(os.path.join(_STATIC, "objectref.js"))
    raws = [raw for raw, _ in PUBLIC_CASES] + PUBLIC_BAD
    out = json.loads(subprocess.run(
        ["node", "--input-type=module", "-e", script, module_url, json.dumps(raws)],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert out["kinds"] == list(objectref.KINDS)
    assert out["formatted"] == [objectref.format_public(kind, "FE81000A2B-0000005-003") for kind in objectref.KINDS]
    expected = ([{"kind": kind, "id": object_id} for _, (kind, object_id) in PUBLIC_CASES]
                + [None] * len(PUBLIC_BAD))
    assert out["parsed"] == expected


def _ids(mysql_config, sql, params=()):
    conn = store.get_connection(mysql_config)
    try:
        return [row["id"] for row in conn.execute(sql, params).fetchall()]
    finally:
        conn.close()


def test_resolve_a_system_and_its_sector(client, seeded_sector):
    config, sector_id, system_ids = seeded_sector
    first, second = pid("system", system_ids[0]), pid("system", system_ids[1])
    body = client.get(f"/api/objects/system:{first}").get_json()
    assert body["ref"] == f"system:{first}" and body["kind"] == "system"
    assert [p["ref"] for p in body["parents"]] == ["galaxy", f"sector:{pid('sector', sector_id)}"]
    assert body["siblings"] == [f"system:{second}"]
    assert body["positions"]["sector_ly"] is not None
    assert body["positions"]["system_km"] == [0.0, 0.0, 0.0]
    # A bare three-part ID is a system's.
    assert client.get(f"/api/objects/{first}").get_json() == body

    sector = client.get(f"/api/objects/sector:{pid('sector', sector_id)}").get_json()
    assert sector["parents"] == [{"ref": "galaxy", "kind": "galaxy", "name": "Galaxy"}]
    assert sector["positions"]["sector_ly"] == [0.0, 0.0, 0.0]


def test_a_placed_sector_gives_galaxy_positions(client, mysql_config):
    sector_id = _place_sector(mysql_config, "Placed", (100.0, 200.0, 5.0))
    body = client.get(f"/api/objects/sector:{pid('sector', sector_id)}").get_json()
    assert body["positions"]["galaxy_pc"] == [100.0, 200.0, 5.0]
    system_id = _ids(mysql_config, "SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,))[0]
    system = client.get(f"/api/objects/system:{pid('system', system_id)}").get_json()
    assert system["positions"]["galaxy_pc"] == pytest.approx([100.0, 200.0, 5.0])


def test_resolve_bodies_down_to_a_moon(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    moon_id = _ids(mysql_config, "SELECT id FROM moons WHERE star_system_id = ? ORDER BY id", (system_id,))[0]
    body = client.get(f"/api/objects/moon:{pid('moon', moon_id)}").get_json()
    kinds = [p["kind"] for p in body["parents"]]
    assert kinds[0] == "galaxy" and kinds[-2:] == ["system", "planet"]
    assert body["parents"][-2]["ref"] == f"system:{pid('system', system_id)}"
    planet = client.get("/api/objects/" + body["parents"][-1]["ref"]).get_json()
    assert planet["kind"] == "planet" and planet["positions"]["system_km"] is not None
    assert body["ref"] not in body["siblings"]
    assert len(body["positions"]["system_km"]) == 3

    for star_id in _ids(mysql_config, "SELECT id FROM stars WHERE star_system_id = ?", (system_id,)):
        star = client.get(f"/api/objects/star:{pid('star', star_id)}").get_json()
        assert star["parents"][-1]["ref"] == f"system:{pid('system', system_id)}"
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
    body = client.get(f"/api/objects/nebula:{pid('nebula', nebula_id)}").get_json()
    assert body["kind"] == "nebula" and body["parents"][0]["ref"] == "galaxy"
    assert body["positions"] == {"galaxy_pc": [5.0, 1.0, 2.0], "sector_ly": None, "system_km": None}


# --- NAV.24: the keep-out radius of every kind of object ------------------------

def test_a_compact_object_at_the_galactic_center_keeps_out_by_its_own_size():
    from planetgen.galaxy import keepout
    from planetgen.physics import constants

    mass = 4e6 * constants.SOLAR_MASS_TO_KG
    at_core = keepout.compact_keep_out(mass, 0.0, 1.2e7)
    assert (at_core.radius_km, at_core.basis) == (1.2e7, "radius")
    far = keepout.compact_keep_out(10 * constants.SOLAR_MASS_TO_KG, 8000.0, 30.0)
    assert far.basis == "galactic_hill" and far.radius_km > 30.0
    # The Hill radius grows with the cube root of the mass and linearly with distance.
    assert keepout.galactic_hill_radius_km(8 * 1e30, 2000.0) == pytest.approx(
        2 * keepout.galactic_hill_radius_km(1e30, 2000.0))
    assert keepout.galactic_hill_radius_km(1e30, 4000.0) == pytest.approx(
        2 * keepout.galactic_hill_radius_km(1e30, 2000.0))


def test_keep_out_for_bodies_and_systems(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    moon_id = _ids(mysql_config, "SELECT id FROM moons WHERE star_system_id = ? ORDER BY id", (system_id,))[0]
    planet_id = _ids(mysql_config, "SELECT planet_id AS id FROM moons WHERE id = ?", (moon_id,))[0]
    conn = store.get_connection(mysql_config)
    try:
        hill = conn.execute("SELECT hill_radius_km FROM planets WHERE id = ?", (planet_id,)).fetchone()
        perimeters = [r["system_perimeter_km"] for r in conn.execute(
            "SELECT system_perimeter_km FROM stars WHERE star_system_id = ?", (system_id,)).fetchall()]
    finally:
        conn.close()

    planet = client.get(f"/api/objects/planet:{pid('planet', planet_id)}").get_json()["keep_out"]
    assert planet["basis"] == "hill" and planet["radius_km"] == pytest.approx(hill["hill_radius_km"])
    assert client.get(f"/api/objects/moon:{pid('moon', moon_id)}").get_json()["keep_out"]["radius_km"] > 0

    system = client.get(f"/api/objects/system:{pid('system', system_id)}").get_json()["keep_out"]
    assert system["basis"] == "perimeter" and system["radius_km"] >= max(perimeters)
    star_id = _ids(mysql_config, "SELECT id FROM stars WHERE star_system_id = ?", (system_id,))[0]
    assert client.get(f"/api/objects/star:{pid('star', star_id)}").get_json()["keep_out"] == system


def test_keep_out_of_a_cloud_is_a_pass_through_with_a_note(client, mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        nebula_id = store.insert_nebula(conn, Nebula(SystemConfig()), sector_id=None, placement={
            "center_x_pc": 5.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 5.5})
        conn.commit()
    finally:
        conn.close()
    keep_out = client.get(f"/api/objects/nebula:{pid('nebula', nebula_id)}").get_json()["keep_out"]
    assert keep_out["radius_km"] is None and keep_out["basis"] == "none" and "passes through" in keep_out["note"]


def test_keep_out_of_compact_objects_and_rogues_uses_the_galactic_hill_radius(client, mysql_config):
    from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
    from planetgen.generation.phenomena.rogue import RoguePlanet

    placement = {"center_x_pc": 8000.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 8000.0}
    conn = store.get_connection(mysql_config)
    try:
        ids = {
            "black_hole": store.insert_black_hole(conn, BlackHole(SystemConfig()), placement=placement),
            "neutron_star": store.insert_neutron_star(conn, NeutronStar(SystemConfig()), placement=placement),
            "rogue_planet": store.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), placement=placement),
        }
        conn.commit()
    finally:
        conn.close()
    for kind, object_id in ids.items():
        keep_out = client.get(f"/api/objects/{kind}:{pid(kind, object_id)}").get_json()["keep_out"]
        assert keep_out["basis"] in ("galactic_hill", "radius") and keep_out["radius_km"] > 0, kind
