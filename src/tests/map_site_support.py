"""
The map pages served with no database (TEST.70): the Flask site with the
`apiclient` functions the Galaxy Map and Sector Map pages call (and
`tilecache`'s reads behind the Galaxy Map's tiles and drill-down stages)
answered from fixture data held here, on a local server for headless
Chromium. `tests/test_web_map_behaviour.py` drives it; any browser test
that needs a map without MySQL can use `map_site` the same way.

The fixture galaxy is a real density shape (`build_galaxy_shape`, the
numbers `tests/test_js_unit.py` uses) with a few generated sectors near
the core; the drill-down stages are worked out from those sectors with
the same block math `queryDb.galaxy_stage` uses (`drill_chain_of`), so
the map can be walked from the whole galaxy down to a sector. The Sector
Map's sector holds stars across the luminosity range, a binary, a
nebula, a neutron star and two rogue planets.
"""

import math
import threading

import pytest
from werkzeug.serving import make_server

from planetgen.web.app import create_app
from planetgen.api.config import Config

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.web.lib import tilecache  # noqa: E402
from planetgen.galaxy.density import build_galaxy_shape, shape_with_terms  # noqa: E402
from planetgen.galaxy.drill import DrillBlock, drill_chain_of, format_drill_key, parse_drill_key  # noqa: E402
from planetgen.galaxy.skeleton import expected_system_count_at_density_1  # noqa: E402
from planetgen.physics.units import pc_to_ly  # noqa: E402
from planetgen.api.limiter import PAGE_LIMITS_OFF  # noqa: E402

DB = "planetgen_map_fixture"
EDGE_PC = 4.0
OUTER_RING = 3749
SOLAR_LUMINOSITY_W = 3.828e26
SUN_RADIUS_KM = 695700.0

SECTOR_ID = 7
"""The sector `/sector/<id>` shows: ring 0, layer 0, slot 1."""

# (id, ring, layer, slot, name): generated sectors, all in the core so a
# walk down the drill-down from the whole galaxy reaches them.
SECTORS = [
    (SECTOR_ID, 0, 0, 1, "Fixture Prime"),
    (8, 0, 0, 2, "Fixture Second"),
    (9, 1, 0, 3, "Fixture Ring One"),
    (10, 2, 1, 5, "Fixture High"),
]


def _star(star_type, temperature_k, radius_solar, luminosity_solar):
    return {
        "star_type": star_type, "temperature_k": temperature_k,
        "radius_km": radius_solar * SUN_RADIUS_KM, "luminosity_w": luminosity_solar * SOLAR_LUMINOSITY_W,
    }


# (id, name, x, y, z in mpc from the sector's center, stars): a red dwarf
# at the dim end, a Sun, a giant, a 1000 L_sun supergiant and a binary.
SYSTEMS = [
    (701, "Dimmest Dwarf", -1200.0, 300.0, -200.0, [_star("M8V", 2600.0, 0.11, 0.0003)]),
    (702, "Middling Sun", 400.0, -500.0, 300.0, [_star("G2V", 5772.0, 1.0, 1.0)]),
    (703, "Orange Giant", 900.0, 800.0, -600.0, [_star("K1III", 4500.0, 15.0, 90.0)]),
    (704, "Blue Beacon", -300.0, -900.0, 700.0, [_star("B1Ib", 21000.0, 25.0, 1000.0)]),
    (705, "Twin Lamps", 100.0, 1100.0, 100.0,
     [_star("F5V", 6500.0, 1.4, 3.0), _star("M2V", 3500.0, 0.4, 0.03)]),
]

PHENOMENA = [
    {"id": 31, "type": "nebula", "name": "Fixture Veil", "descriptor": "emission", "radius_ly": 1.2,
     "distance_ly": 3.0, "offset_x_ly": 2.5, "offset_y_ly": -2.0, "offset_z_ly": 1.0},
    {"id": 41, "type": "neutron_star", "name": "Fixture Pulsar", "descriptor": "pulsar", "radius_ly": 0.0,
     "distance_ly": 2.0, "offset_x_ly": -2.5, "offset_y_ly": -2.5, "offset_z_ly": -1.0},
    {"id": 51, "type": "rogue_planet", "name": "Fixture Wanderer", "descriptor": "frozen", "radius_ly": 0.0,
     "distance_ly": 1.5, "offset_x_ly": 1.5, "offset_y_ly": 2.5, "offset_z_ly": -1.5, "class": "T"},
    {"id": 52, "type": "rogue_planet", "name": "Fixture Drifter", "descriptor": "frozen", "radius_ly": 0.0,
     "distance_ly": 2.5, "offset_x_ly": -1.0, "offset_y_ly": 3.0, "offset_z_ly": 2.0, "class": "T"},
]


STAR_FIELD_CENTER_PC = (3000.0, 0.0, 0.0)
"""Where the Galaxy Map's fixture stars sit: a field 40 pc across, well
away from the generated sectors near the core."""


def _star_field():
    """The Galaxy Map's fixture stars (a tile's `stars`, the
    `galaxy_bright_stars_in_box` shape), luminosities 1e-4 to 1e5 L_sun."""
    import random

    from planetgen.galaxy.geometry import sector_address_at

    rng = random.Random(87)
    stars = []
    for i in range(240):
        lum = 10 ** rng.uniform(-4.0, 5.0)
        x = STAR_FIELD_CENTER_PC[0] + rng.uniform(-20.0, 20.0)
        y = STAR_FIELD_CENTER_PC[1] + rng.uniform(-20.0, 20.0)
        z = STAR_FIELD_CENTER_PC[2] + rng.uniform(-6.0, 6.0)
        ring, layer, slot = sector_address_at((x, y, z), EDGE_PC)
        temperature = 2600.0 + 30000.0 * (math.log10(lum) + 4.0) / 9.0 * rng.uniform(0.6, 1.0)
        stars.append({
            "id": 9000 + i, "name": None,
            "x": x, "y": y, "z": z,
            "luminosity_sol": float("%.4g" % lum), "temperature_k": round(temperature),
            "radius_sol": float("%.3g" % (lum ** 0.35)), "star_type": "fixture",
            "ring_index": ring, "layer_index": layer, "ring_slot_index": slot,
            # Every star "generated" (in a filled sector) so its panel links to a system page.
            "system_id": SYSTEMS[i % len(SYSTEMS)][0],
        })
    return stars


STAR_FIELD = _star_field()

NEBULA_CLOUD = {
    "type": "nebula", "id": 801, "name": "Fixture Cloud", "descriptor": "dark", "class": "D",
    "radius_pc": 30.0, "x": STAR_FIELD_CENTER_PC[0] - 45.0, "y": STAR_FIELD_CENTER_PC[1], "z": STAR_FIELD_CENTER_PC[2],
}
"""The Galaxy Map's fixture nebula: beside the star field, big enough to draw
from its shape mesh once the view is on that block (MAP.103)."""


NEBULA_OVER_THE_FIELD = {
    "type": "nebula", "id": 802, "name": "Fixture Pall", "descriptor": "dark", "class": "D",
    "radius_pc": 14.0, "x": STAR_FIELD_CENTER_PC[0], "y": STAR_FIELD_CENTER_PC[1], "z": STAR_FIELD_CENTER_PC[2],
}
"""A second fixture nebula, right over the star field, so a close view of that
block (all unfilled space: the generated sectors are in the core) is inside it."""
NEBULAE = (NEBULA_CLOUD, NEBULA_OVER_THE_FIELD)


def nebula_shape_payload(nebula_id, lod):
    """What `GET /api/nebulae/<id>/shape` answers for the fixture nebula."""
    import random

    from planetgen.galaxy import nebula_shape

    known = {cloud["id"] for cloud in NEBULAE} | {p["id"] for p in PHENOMENA if p["type"] == "nebula"}
    if int(nebula_id) not in known:
        raise apiclient.NotFoundError(f"no such nebula: {nebula_id}")
    vertices, faces = nebula_shape.draw_shape(random.Random(int(nebula_id))).mesh(lod)
    cloud = next((c for c in NEBULAE if c["id"] == int(nebula_id)), NEBULA_CLOUD)
    return {"id": int(nebula_id), "radius_ly": pc_to_ly(cloud["radius_pc"]),
            "center_pc": [cloud["x"], cloud["y"], cloud["z"]], "lod": lod,
            "vertices": [[round(float(v), 5) for v in vertex] for vertex in vertices],
            "faces": [[int(i) for i in face] for face in faces]}


POINT_FIELD = [
    {"type": "black_hole", "id": 61, "name": "Fixture Maw", "descriptor": "accreting", "luminosity_sol": 0.02,
     "x": STAR_FIELD_CENTER_PC[0] - 12.0, "y": 9.0, "z": 0.0},
    {"type": "neutron_star", "id": 41, "name": "Fixture Pulsar", "descriptor": "millisecond",
     "luminosity_sol": 0.5, "x": STAR_FIELD_CENTER_PC[0] + 11.0, "y": 3.0, "z": 0.0},
    {"type": "quasar", "id": 71, "name": "Fixture Beacon", "descriptor": "radio-loud", "luminosity_sol": 3e12,
     "x": STAR_FIELD_CENTER_PC[0] + 4.0, "y": 14.0, "z": 0.0},
]
"""The Galaxy Map's fixture point phenomena (a tile's `points`, MAP.80),
among the fixture stars."""


def galaxy_shape():
    """The stored shape `GET /api/galaxy/shape` would return."""
    shape = shape_with_terms(build_galaxy_shape(2800.0, 350.0, 200.0, 1.0, 2, math.radians(15.0), 0.4))
    shape.update({
        "edge_pc": EDGE_PC, "outer_ring_index": OUTER_RING,
        "expected_system_count_at_density_1": expected_system_count_at_density_1(pc_to_ly(EDGE_PC)),
    })
    return shape


def _sector_center_pc(ring, layer, slot):
    from planetgen.galaxy.geometry import sector_position_pc

    return sector_position_pc(ring, layer, slot, EDGE_PC)


def sector_detail(sector_id=SECTOR_ID):
    """`GET /api/sectors/<id>` for the fixture sector."""
    _id, ring, layer, slot, name = next(s for s in SECTORS if s[0] == sector_id)
    cx, cy, cz = _sector_center_pc(ring, layer, slot)
    systems = []
    for system_id, system_name, x, y, z, stars in (SYSTEMS if sector_id == SECTOR_ID else []):
        systems.append({
            "id": system_id, "name": system_name, "is_binary": int(len(stars) > 1), "binary_type": None,
            "stars": stars, "quadrant": "+x+y+z", "location": f"{name} -- near the middle",
            "center_distance_ly": round(math.sqrt(x * x + y * y + z * z) / 306.6, 2),
            "position_x_mpc": x, "position_y_mpc": y, "position_z_mpc": z,
        })
    return {
        "id": sector_id, "name": name, "edge_mpc": EDGE_PC * 1000.0, "edge_ly": pc_to_ly(EDGE_PC),
        "ring_index": ring, "layer_index": layer, "ring_slot_index": slot, "placed": True,
        "center_x_pc": cx, "center_y_pc": cy, "center_z_pc": cz, "wiki_url": None,
        "systems": systems,
        "phenomena": [dict(p) for p in PHENOMENA] if sector_id == SECTOR_ID else [],
        "neighbors": [],
    }


UNCHARTED_CELL = (0, 0, 0)
"""A cell of the core no sector was generated in, where the fixture's scatters left a bright star, a black
hole and a nebula (MAP.162)."""

SECTOR_STATS = {
    SECTOR_ID: (900, 1200, 4.5, 3000.0),
    8: (40, 60, 0.5, 900.0),
    9: (0, 0, None, 0.0),
    10: (6, 8, 11.0, 4.0),
}
"""Each fixture sector's stored systems, stars, mean star age (Gy) and summed
luminosity (DB.14): a dense average one, a young bright one, one with no
stars and a sparse old one."""


def galaxy_stage(at=None):
    """`queryDb.galaxy_stage` over `SECTORS`, with no database."""
    chains = [(s, drill_chain_of(s[1], s[2], s[3])) for s in SECTORS]
    if not at:
        counts = {}
        for s, chain in chains:
            counts.setdefault(chain[0], []).append(s[0])
        return {"at": None, "child_m": 243, "children": _children(counts), "sectors": None}
    block = parse_drill_key(at)
    level = (243, 27, 3).index(block.m)
    inside = [(s, chain) for s, chain in chains if chain[level] == block]
    child_m = (27, 3, 1)[level]
    counts = {}
    for s, chain in inside:
        counts.setdefault(chain[level + 1], []).append(s[0])
    sectors = None
    if child_m == 1:
        sectors = [{"ring": s[1], "layer": s[2], "slot": s[3], "id": s[0], "name": s[4],
                    "system_count": len(SYSTEMS) if s[0] == SECTOR_ID else 0, "stats": _stats([s[0]])}
                   for s, _chain in inside]
    return {"at": format_drill_key(block), "child_m": child_m, "children": _children(counts), "sectors": sectors}


def _stats(ids):
    """`queryDb._stage_stats` over these fixture sectors."""
    rows = [SECTOR_STATS[i] for i in ids]
    stars = sum(r[1] for r in rows)
    return {"systems": sum(r[0] for r in rows), "expected_systems": 0.0, "stars": stars,
            "mean_age_gy": sum(r[2] * r[1] for r in rows if r[1]) / stars if stars else None,
            "luminosity_sol": sum(r[3] for r in rows)}


def _children(counts):
    return [{"ring": b.ring, "wedge": b.wedge, "slab": b.slab, "generated": len(ids), "stats": _stats(ids)}
            for b, ids in sorted(counts.items(), key=lambda item: (item[0].slab, item[0].ring, item[0].wedge))]


def system_scene_payload(system_id):
    """The 3D system view's scene (`GET /api/systems/<id>/scene`) for a fixture
    system: a close pair with two planets (one with a moon), a belt and two
    comets, built from the hand-made scene of `test_body_positions.py`."""
    from tests.test_body_positions import AU_KM, sample_scene

    name = next((n for sid, n, *_rest in SYSTEMS if sid == int(system_id)), None)
    if name is None:
        raise apiclient.NotFoundError(f"No such system: {system_id}")
    scene = sample_scene()
    colors = ["#ffd9a0", "#9fb8ff", "#5f8a72", "#d9a441", "#8a7a66", "#c9d6e0", "#b9d6e8", "#b9d6e8"]
    n = 0
    for index, star in enumerate(scene["stars"], start=1):
        star.update(id=index, kind="star", name=f"{name} {'AB'[index - 1]}", role="primary" if index == 1 else "secondary",
                    star_type="G2V", radius_km=695700.0 * (1.0 if index == 1 else 0.6), mass_kg=2e30,
                    temperature_k=5772.0, luminosity_w=3.8e26, color=colors[n])
        n += 1
    for index, planet in enumerate(scene["planets"], start=1):
        planet.update(id=index, kind="planet", name=f"{name} {'I' * index}", planet_class="F" if index == 1 else "J",
                      radius_km=6371.0 if index == 1 else 69911.0, mass_kg=6e24, color=colors[n],
                      habitable=False, inhabited=False)
        n += 1
        for moon_index, moon in enumerate(planet["moons"], start=1):
            moon.update(id=moon_index, kind="moon", name=f"{planet['name']} {'abc'[moon_index - 1]}", planet_class="D",
                        radius_km=1737.0, mass_kg=7e22, color=colors[n], habitable=False, inhabited=False)
            n += 1
    for index, comet in enumerate(scene["comets"], start=1):
        comet.update(id=index, kind="comet", name=f"{name} Comet {index}", radius_km=5.0, color="#b9d6e8", is_active=True)
    scene["belts"] = [{"ref": "belt:1", "id": 1, "kind": "belt", "name": None, "body_type": "a", "around": "barycenter",
                       "inner_km": 2.5 * AU_KM, "outer_km": 3.2 * AU_KM, "distance_km": 2.8 * AU_KM}]
    scene["system"] = {"id": int(system_id), "ref": f"system:{int(system_id)}", "name": name, "is_binary": True,
                       "binary_configuration": "close", "heliopause_au": 120.0}
    scene["epoch"] = None
    scene["epoch_unix"] = None
    return scene


class FixtureApi:
    """The `apiclient`/`tilecache` reads the map pages make."""

    def __init__(self):
        self.shape = galaxy_shape()

    def get_galaxy_sectors(self, db):
        out = []
        for sector_id, ring, layer, slot, name in SECTORS:
            cx, cy, cz = _sector_center_pc(ring, layer, slot)
            out.append({"id": sector_id, "name": name, "x": cx, "y": cy, "z": cz,
                        "galactic_radius_pc": math.hypot(cx, cy), "ring_index": ring, "layer_index": layer,
                        "ring_slot_index": slot, "system_count": len(SYSTEMS) if sector_id == SECTOR_ID else 0})
        return out

    def get_galaxy_shape(self, db):
        return self.shape

    def get_polities(self, db, limit=None, offset=None):
        return {"items": [], "total": 0, "limit": limit, "offset": offset}

    def auth_me(self, cookie_header):
        return None

    def get_galaxy_changes(self, db, since=None):
        return {"stamp": "00000000000000f1", "state": "s", "full": since is None, "tiles": []}

    def get_galaxy_tiles(self, db, tile_keys):
        from planetgen.galaxy.viewport import parse_tile_key, tile_bounds_pc

        tiles = {}
        for key in tile_keys:
            lo, hi = tile_bounds_pc(*parse_tile_key(key))
            stars = [star for star in STAR_FIELD
                     if all(lo[i] <= star[axis] < hi[i] for i, axis in enumerate("xyz"))]
            points = [point for point in POINT_FIELD
                      if all(lo[i] <= point[axis] < hi[i] for i, axis in enumerate("xyz"))]
            # Like the real tile: every cloud whose sphere reaches into the box.
            clouds = []
            for cloud in NEBULAE:
                nearest = [min(max(cloud[axis], lo[i]), hi[i]) for i, axis in enumerate("xyz")]
                if math.dist(nearest, [cloud[axis] for axis in "xyz"]) <= cloud["radius_pc"]:
                    clouds.append(cloud)
            tiles[key] = {"placed": [], "planned": [], "filled": None, "clouds": clouds, "stars": stars,
                          "generated": None, "points": points}
        return {"tiles": tiles, "edge_pc": EDGE_PC, "has_shape": True}

    def get_galaxy_stage(self, db, at=None):
        return galaxy_stage(at)

    def get_population_status(self, db):
        return dict(apiclient.POPULATION_NONE)

    def get_galaxy_locate(self, db, q):
        return []

    def get_nebula_shape(self, db, nebula_id, lod="low"):
        return nebula_shape_payload(nebula_id, lod)

    def get_sector(self, db, sector_id):
        if int(sector_id) not in [s[0] for s in SECTORS]:
            raise apiclient.NotFoundError(f"No such sector: {sector_id}")
        return sector_detail(int(sector_id))

    def get_sector_facilities(self, db, sector_id):
        return []

    def get_system_scene(self, db, system_id):
        return system_scene_payload(system_id)

    def get_bright_stars_in_cell(self, db, ring_index, layer_index, ring_slot_index):
        return []

    def get_uncharted_sector(self, db, ring_index, layer_index, ring_slot_index):
        """`GET /api/galaxy/uncharted`: `UNCHARTED_CELL` holds a bright star and a black hole and a
        nebula that no sector was made for; every other cell is empty (MAP.162)."""
        from planetgen.galaxy.geometry import provisional_sector_designation

        address = (ring_index, layer_index, ring_slot_index)
        cx, cy, cz = _sector_center_pc(*address)
        generated = next((s[0] for s in SECTORS if s[1:4] == address), None)
        cell = {"ring_index": ring_index, "layer_index": layer_index, "ring_slot_index": ring_slot_index,
                "designation": provisional_sector_designation(*address), "center_pc": [cx, cy, cz],
                "edge_pc": EDGE_PC, "sector_id": generated, "stars": [], "scattered": []}
        if address == UNCHARTED_CELL:
            cell["stars"] = [{"id": 1, "x": cx + 0.5, "y": cy + 0.4, "z": cz, "luminosity_sol": 9000.0,
                              "temperature_k": 22000.0, "radius_sol": 6.0, "star_type": "B2V"}]
            cell["scattered"] = [
                {"id": 1, "kind": "black-hole", "subtype": "stellar", "type": "black_hole",
                 "x": cx - 0.6, "y": cy, "z": cz + 0.3},
                {"id": 2, "kind": "planetary-nebula", "subtype": None, "type": "nebula",
                 "x": cx, "y": cy - 0.7, "z": cz - 0.2},
            ]
        return cell


_APICLIENT = ("get_galaxy_sectors", "get_galaxy_shape", "get_polities", "auth_me", "get_galaxy_locate",
              "get_sector", "get_sector_facilities", "get_bright_stars_in_cell", "get_population_status",
              "get_galaxy_changes", "get_nebula_shape", "get_system_scene", "get_uncharted_sector")
_TILECACHE = ("get_galaxy_changes", "get_galaxy_tiles", "get_galaxy_stage")


class _SiteConfig(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "map-fixture-secret"
    RATELIMIT_PAGES = PAGE_LIMITS_OFF


@pytest.fixture(scope="module")
def map_site(tmp_path_factory):
    """The base URL of a local server for the map pages, on fixture data."""
    api = FixtureApi()
    patch = pytest.MonkeyPatch()
    try:
        for name in _APICLIENT:
            patch.setattr(apiclient, name, getattr(api, name))
        for name in _TILECACHE:
            patch.setattr(tilecache, name, getattr(api, name))
        patch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
        patch.setenv("PLANETGEN_TILE_CACHE_DIR", str(tmp_path_factory.mktemp("tiles")))
        patch.delenv("PLANETGEN_TILE_CACHE_MAX_MB", raising=False)
        app = create_app(_SiteConfig)
        app.testing = True
        server = make_server("127.0.0.1", 0, app, threaded=True)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            yield f"http://127.0.0.1:{server.server_port}"
        finally:
            server.shutdown()
            thread.join(timeout=10)
    finally:
        patch.undo()


__all__ = ["map_site", "DrillBlock", "SECTOR_ID", "SYSTEMS", "PHENOMENA", "SECTORS", "UNCHARTED_CELL"]
