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

from api.app import create_app
from api.config import Config

import web  # noqa: F401 -- puts src/html/lib on sys.path
import apiclient  # noqa: E402
import tilecache  # noqa: E402
from stellarObjects.galaxyDensity import build_galaxy_shape  # noqa: E402
from stellarObjects.galaxyDrill import DrillBlock, drill_chain_of, format_drill_key, parse_drill_key  # noqa: E402
from stellarObjects.galaxySkeleton import expected_system_count_at_density_1  # noqa: E402
from stellarObjects.utils import pc_to_ly  # noqa: E402
from api.limiter import PAGE_LIMITS_OFF  # noqa: E402

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

    from stellarObjects.galaxyGeometry import sector_address_at

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


def galaxy_shape():
    """The stored shape `GET /api/galaxy/shape` would return."""
    shape = build_galaxy_shape(2800.0, 350.0, 200.0, 1.0, 2, math.radians(15.0), 0.4)._asdict()
    shape.update({
        "edge_pc": EDGE_PC, "outer_ring_index": OUTER_RING,
        "expected_system_count_at_density_1": expected_system_count_at_density_1(pc_to_ly(EDGE_PC)),
    })
    return shape


def _sector_center_pc(ring, layer, slot):
    from stellarObjects.galaxyGeometry import sector_position_pc

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


def galaxy_stage(at=None):
    """`queryDb.galaxy_stage` over `SECTORS`, with no database."""
    chains = [(s, drill_chain_of(s[1], s[2], s[3])) for s in SECTORS]
    if not at:
        counts = {}
        for _s, chain in chains:
            counts[chain[0]] = counts.get(chain[0], 0) + 1
        return {"at": None, "child_m": 243, "children": _children(counts), "sectors": None}
    block = parse_drill_key(at)
    level = (243, 27, 3).index(block.m)
    inside = [(s, chain) for s, chain in chains if chain[level] == block]
    child_m = (27, 3, 1)[level]
    counts = {}
    for _s, chain in inside:
        counts[chain[level + 1]] = counts.get(chain[level + 1], 0) + 1
    sectors = None
    if child_m == 1:
        sectors = [{"ring": s[1], "layer": s[2], "slot": s[3], "id": s[0], "name": s[4],
                    "system_count": len(SYSTEMS) if s[0] == SECTOR_ID else 0} for s, _chain in inside]
    return {"at": format_drill_key(block), "child_m": child_m, "children": _children(counts), "sectors": sectors}


def _children(counts):
    return [{"ring": b.ring, "wedge": b.wedge, "slab": b.slab, "generated": n}
            for b, n in sorted(counts.items(), key=lambda item: (item[0].slab, item[0].ring, item[0].wedge))]


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
        from stellarObjects.galaxyViewport import parse_tile_key, tile_bounds_pc

        tiles = {}
        for key in tile_keys:
            lo, hi = tile_bounds_pc(*parse_tile_key(key))
            stars = [star for star in STAR_FIELD
                     if all(lo[i] <= star[axis] < hi[i] for i, axis in enumerate("xyz"))]
            tiles[key] = {"placed": [], "planned": [], "filled": None, "clouds": [], "stars": stars,
                          "generated": None}
        return {"tiles": tiles, "edge_pc": EDGE_PC, "has_shape": True}

    def get_galaxy_stage(self, db, at=None):
        return galaxy_stage(at)

    def get_population_status(self, db):
        return dict(apiclient.POPULATION_NONE)

    def get_galaxy_locate(self, db, q):
        return []

    def get_sector(self, db, sector_id):
        if int(sector_id) not in [s[0] for s in SECTORS]:
            raise apiclient.NotFoundError(f"No such sector: {sector_id}")
        return sector_detail(int(sector_id))

    def get_sector_facilities(self, db, sector_id):
        return []

    def get_bright_stars_in_cell(self, db, ring_index, layer_index, ring_slot_index):
        return []


_APICLIENT = ("get_galaxy_sectors", "get_galaxy_shape", "get_polities", "auth_me", "get_galaxy_locate",
              "get_sector", "get_sector_facilities", "get_bright_stars_in_cell", "get_population_status",
              "get_galaxy_changes")
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


__all__ = ["map_site", "DrillBlock", "SECTOR_ID", "SYSTEMS", "PHENOMENA", "SECTORS"]
