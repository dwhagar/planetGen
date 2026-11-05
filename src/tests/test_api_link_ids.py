# tests/test_api_link_ids.py

"""API.40: every field of an API answer that feeds a link to a system or a body page carries an object ID, never a
row number. The audit builds one grid-placed sector (with a binary system, comets, moons and a rogue planet), the
galaxy containers and the object index, fetches every endpoint that answers with systems, bodies, sectors or
references to them, and walks each JSON answer: a whole number under `id`, `*_id`, `*_ids` or a kind key, or a
`<kind>:<number>` reference, fails the test unless the endpoint is listed in `ROW_NUMBERS_KEPT` with its reason.
Every GET endpoint of the API is either audited here or listed in `NOT_AUDITED` with the reason it carries no link."""

import math
import re

import pytest

from planetgen.api.config import Config
from planetgen.db import containers as containers_db, objectindex, store
from planetgen.galaxy import object_uid, spherical
from planetgen.galaxy.drill import drill_chain_of, format_drill_key
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import clusters, containers
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem
from planetgen.physics.units import ly_to_pc
from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
from planetgen.web.app import create_app
from tests.test_api import _rebuild_block_stats
from tests.test_galaxy_containers import EDGE_PC, REACH_PC, SEED, SHAPE

pytestmark = pytest.mark.fresh_database

KIND_KEYS = {"system", "star", "planet", "moon", "belt", "comet", "host", "parent", "around", "owner", "target", "body",
             "object", "nebula", "phenomenon", "facility", "star_system", "capital", "sector"}
"""set: Keys that name an object, so a whole number under one is a row number."""

COUNT_KEYS = {"counts", "by_kind", "facets", "totals"}
"""set: Keys whose children are counts or choices keyed by a kind name, not objects."""

REFERENCE = re.compile(r"^(?P<kind>[a-z_]+):\d+$")
PAGE_PATH = re.compile(r"/(?:system|planet|moon|star|belt|comet|object|sector)s?/\d+(?:\b|$)")

ROW_NUMBERS_KEPT = {
    "api.containers_list": {"containers[].id": "container row", "containers[].component_id": "galaxy component row"},
    "api.container_detail": {"id": "container row", "component_id": "galaxy component row"},
    "api.container_contents": {"container.component_id": "galaxy component row"},
    "api.containers_at": {"containers[].id": "container row", "containers[].component_id": "galaxy component row"},
}
"""dict: The row numbers an answer still carries, with what they are. A container is not a system or a body page and its
`id` is not a link to one (`/api/containers/<id>` takes the printed container ID the answer calls `printed_id`); they wait
for the decision whether to print them too."""

NOT_AUDITED = {
    "admin": "administration pages carry no object", "auth": "accounts carry no object", "backups": "files, no object",
    "images": "files, no object", "settings": "settings, no object", "workqueue": "work runs, no object",
    "api.databases": "database names", "api.health": "a status", "api.job_status_route": "a job number",
    "api.page_job_status_route": "a job number", "facility_types.facility_type": "a facility type table",
    "facility_types.list_facility_types": "a facility type table",
    "api.facility": "no facility in the fixture; its `id` has an `ID_PATHS` rule",
    "api.nav_chart": "a drawing of `api.nav`, which is audited",
    "api.nebula_shape_route": "no nebula in the fixture; routed by printed ID", "api.nebula_surroundings_route": "no nebula in the fixture",
    "api.nebula_volume_route": "no nebula in the fixture", "api.galaxy_nebula_volume_route": "no nebula in the fixture",
    "api.galaxy_bright_stars_in_cell": "scattered stars, not built systems", "api.galaxy_profile": "density profile numbers",
    "api.galaxy_uncharted_sector": "uncharted cell numbers", "api.uncharted_phenomenon": "scattered phenomenon, no built object",
    "api.uncharted_system": "scattered star, no built object",
    "population.species_detail": "no species in the fixture; a species number is not an object",
    "population.polity_detail": "no polity in the fixture; a polity number is not an object",
    "population.species_habitable_in_system": "no species in the fixture; its `unaided[].id` and `equipped[].id` have `ID_PATHS` rules",
    "population.species_habitable_in_sector": "no species in the fixture; its `unaided[].id` and `equipped[].id` have `ID_PATHS` rules",
}
"""dict: Endpoint names (or the part before the dot) the audit does not fetch, with why."""


def _walk(node, path, found, under_counts=False):
    if isinstance(node, dict):
        for key, value in node.items():
            here = f"{path}.{key}" if path else key
            counted = under_counts or key in COUNT_KEYS
            named = key == "id" or key.endswith("_id") or key in KIND_KEYS
            if isinstance(value, int) and not isinstance(value, bool) and named and not counted:
                found.add(here)
            elif isinstance(value, list) and (key.endswith("_ids") or key == "ids") and any(isinstance(v, int) for v in value):
                found.add(here)
            elif isinstance(value, str) and not counted:
                match = REFERENCE.match(value)
                if match and match["kind"] != "sector" and (key in ("ref", "from", "to") or key.endswith("_ref")):
                    found.add(here)
                elif key in ("url", "href", "link") and PAGE_PATH.search(value):
                    found.add(here)
            _walk(value, here, found, counted)
    elif isinstance(node, list):
        for item in node[:200]:
            _walk(item, path + "[]", found, under_counts)


def test_the_walker_flags_row_numbers_and_leaves_printed_ids_counts_and_choices():
    found = set()
    _walk({"id": 5, "sector_id": "5000080001", "items": [{"system_id": 7, "ref": "planet:12"}], "parents": [{"ref": "sector:5000080001"}],
           "counts": {"planet": 3}, "by_kind": {"moon": 2}, "facets": {"x": [{"value": "microbial:1"}]}, "url": "/system/9"},
          "", found)
    assert found == {"id", "items[].system_id", "items[].ref", "url"}


def test_every_get_endpoint_is_audited_or_has_a_reason_it_carries_no_link():
    app = create_app(Config)
    audited = set(AUDITED_ENDPOINTS)
    for rule in app.url_map.iter_rules():
        if "GET" not in rule.methods or not rule.rule.startswith("/api"):
            continue
        reason = NOT_AUDITED.get(rule.endpoint) or NOT_AUDITED.get(rule.endpoint.split(".")[0])
        assert rule.endpoint in audited or reason, f"{rule.endpoint} ({rule.rule}) is neither audited nor given a reason"
    assert not set(ROW_NUMBERS_KEPT) - audited


AUDITED_ENDPOINTS = (
    "api.containers_list", "api.container_detail", "api.container_contents", "api.containers_at", "api.facility_list",
    "api.facility_orbit", "api.galaxy_cell", "api.galaxy_changes_route", "api.galaxy_charted_regions_route",
    "api.galaxy_locate_route", "api.galaxy_made", "api.galaxy_nebula_cover_route", "api.galaxy_nodes_route",
    "api.galaxy_phenomena", "api.galaxy_sectors", "api.galaxy_shape", "api.galaxy_stage_route", "api.galaxy_stamp_route",
    "api.galaxy_superbubbles", "api.galaxy_window_route", "api.lookup_route", "api.moon_life_air", "api.moon_life_readings",
    "api.nav", "api.near", "api.object_children_route", "api.object_list", "api.object_ref_route", "api.object_sheet_route",
    "api.phenomena", "api.phenomenon", "api.phenomenon_contents", "api.planet_life_air", "api.planet_life_readings",
    "api.search", "api.sector_detail", "api.sector_facilities", "api.sectors", "api.system_description", "api.system_detail",
    "api.system_facilities", "api.system_scene", "api.system_text", "api.systems", "api.table_estimates_route",
    "api.uncharted_systems", "population.moon_habitable_species", "population.planet_habitable_species",
    "population.planet_species", "population.polity_list", "population.population_status",
    "population.species_list", "population.system_habitable_species", "population.system_owner", "population.territories",
)
"""tuple: The endpoints `test_every_endpoint_answer_carries_object_ids_only` fetches."""


def _build(mysql_config):
    center = (5.0, 5.0, 5.0)
    address = spherical.sector_address_at(center, ly_to_pc(DEFAULT_SECTOR_EDGE_LY))
    sector = SpaceSector("Audit Sector", edge_ly=DEFAULT_SECTOR_EDGE_LY)
    for _ in range(400):
        config = SystemConfig()
        config.PLANETS = True
        config.MAX_PLANETS = True
        config.BINARY_SYSTEM = True
        system = StarSystem(system_config=config)
        if len(system.stars) > 1 and system.comets and any(getattr(p, "moons", None) for p in system.planets):
            break
    sector.add_system(system, position=(0.5, 0.5, 0.5), system_config=config)
    bare = SystemConfig()
    bare.PLANETS = False
    bare.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=bare), position=(-1.0, 1.0, 0.0), system_config=bare)
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    position = {"center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
                "galactic_radius_pc": math.dist(center, (0, 0, 0)), "shell_index": address[0], "band_index": address[1],
                "wedge_index": address[2]}
    store.save_sector(sector, config=mysql_config, galaxy_position=position)
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, clusters.generate_clusters(SEED, count=3), kind="cluster")
        conn.commit()
        made = containers.build(SEED, SHAPE, EDGE_PC, REACH_PC, store.get_galaxy_components(conn, kind="cluster"))
        containers_db.replace(conn, made)
        conn.commit()
        objectindex.build(conn)
        conn.commit()
    finally:
        conn.close()
    _rebuild_block_stats(mysql_config)
    return address, made[0]


def test_every_endpoint_answer_carries_object_ids_only(mysql_config):
    address, first = _build(mysql_config)

    class AuditConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(AuditConfig)
    app.testing = True
    client = app.test_client()
    sector_id = client.get("/api/sectors").get_json()["items"][0]["id"]
    systems = client.get("/api/systems").get_json()["items"]
    for row in systems:
        detail = client.get(f"/api/systems/{row['id']}").get_json()
        if detail["planets"]:
            break
    system_id, planet = row["id"], detail["planets"][0]
    moon = next(m for p in detail["planets"] for m in p["moons"])
    rogue = client.get("/api/phenomena").get_json()["items"][0]
    chain = [format_drill_key(b) for b in drill_chain_of(*address)]
    container = object_uid.format_id(object_uid.from_bytes(first.uid))
    cx, cy, cz = first.center_pc
    other = systems[-1]["id"]
    urls = [
        f"/api/containers/{container}", f"/api/containers/{container}/contents", f"/api/containers/at?x={cx}&y={cy}&z={cz}",
        "/api/containers", "/api/galaxy/charted-regions", "/api/galaxy/nebula-cover", "/api/galaxy/stamp", "/api/galaxy/changes",
        "/api/galaxy/shape", "/api/table-estimates", "/api/facilities/orbit?host_type=planet&host_id=" + planet["id"],
        "/api/facilities", "/api/sectors", f"/api/sectors/{sector_id}", f"/api/sectors/{sector_id}/facilities", "/api/systems",
        f"/api/systems/{system_id}", f"/api/systems/{system_id}/scene", f"/api/systems/{system_id}/description",
        f"/api/systems/{system_id}/facilities", f"/api/systems/{system_id}/text", f"/api/systems/{system_id}/owner",
        f"/api/systems/{system_id}/habitable-species", "/api/phenomena", f"/api/phenomena/{rogue['type']}/{rogue['id']}",
        f"/api/phenomena/{rogue['type']}/{rogue['id']}/contents", f"/api/objects/system:{system_id}",
        f"/api/objects/planet:{planet['id']}", f"/api/objects/planet:{planet['id']}/children",
        f"/api/objects/planet:{planet['id']}/sheet", f"/api/objects/system:{system_id}/children",
        f"/api/near?from=system:{system_id}&distance=30", f"/api/near?from=planet:{planet['id']}&distance=30",
        "/api/near?point=5,5,5&distance=30", "/api/search?q=Audit", f"/api/search?q={detail['name']}",
        f"/api/lookup?text={detail['name']}", f"/api/lookup?text={system_id}", f"/api/galaxy/nodes?keys={','.join(chain)}",
        f"/api/galaxy/stage?at={chain[-2]}", "/api/galaxy/stage", "/api/galaxy/sectors", "/api/galaxy/made?since=0",
        "/api/galaxy/phenomena", "/api/galaxy/locate?q=Audit", f"/api/galaxy/locate?q={detail['name']}",
        "/api/galaxy/cell?x=5&y=5&z=5", f"/api/galaxy/window?at={chain[-3]}", "/api/galaxy/superbubbles",
        f"/api/nav?from=system:{system_id}&to=system:{other}",
        *[f"/api/object-lists/{kind}" for kind in ("system", "star", "planet", "moon", "belt", "comet", "rogue_planet")],
        f"/api/object-lists/planet?sector_id={sector_id}", "/api/uncharted-systems", f"/api/planets/{planet['id']}/life-air",
        f"/api/planets/{planet['id']}/life-readings", f"/api/planets/{planet['id']}/habitable-species",
        f"/api/planets/{planet['id']}/species", f"/api/moons/{moon['id']}/life-air", f"/api/moons/{moon['id']}/life-readings",
        f"/api/moons/{moon['id']}/habitable-species", "/api/polities", "/api/species", "/api/territories", "/api/population",
    ]
    adapter = app.url_map.bind("localhost")
    reached, leaks, empty = set(), {}, []
    for url in urls:
        path = url.split("?")[0]
        endpoint = adapter.match(path)[0]
        reached.add(endpoint)
        response = client.get(url)
        assert response.status_code in (200, 404, 400), (url, response.status_code)
        body = response.get_json()
        found = set()
        _walk(body, "", found)
        extra = found - set(ROW_NUMBERS_KEPT.get(endpoint, {}))
        if extra:
            leaks[url] = sorted(extra)
        if response.status_code == 200 and not body:
            empty.append(url)
    assert not leaks, f"row numbers where an object ID belongs: {leaks}"
    assert reached == set(AUDITED_ENDPOINTS), sorted(reached ^ set(AUDITED_ENDPOINTS))
    lists = [client.get(f"/api/object-lists/{kind}").get_json()["items"] for kind in ("system", "planet", "moon", "star")]
    assert all(lists), "the audit needs items to walk"
    scoped = client.get(f"/api/object-lists/planet?sector_id={sector_id}").get_json()
    assert scoped["items"] and {i["sector_id"] for i in scoped["items"]} == {sector_id}
    assert client.get("/api/object-lists/planet?sector_id=4").status_code == 200  # a number is a sector ID in hex
    assert client.get("/api/object-lists/planet?sector_id=zz").status_code == 400
    assert client.get(f"/api/planets/{planet['id']}/species").get_json()["error"].endswith(planet["id"])
