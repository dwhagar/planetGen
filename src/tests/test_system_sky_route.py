"""VIEW.31 (System sky S13): the sky route `GET /api/systems/<id>/sky` and `GET /system/<id>/sky`.

* `db.sky_sprites.objects_in_sphere`: the nebulae and remnants of a sphere, against a brute force over the stored rows;
* `galaxy.sky_sprites.sprites_from`: the sprite records the browser reads (`html/static/skysprites.js`), on synthetic objects;
* `api.system_sky.build_system_sky`: a real database, the packed list and its `meta` (the sprites and the report);
* the API route and the web route (limit checks, 404, the disk cache by version key, the clear).
"""

import base64
import math
import os

import numpy as np
import pytest

from planetgen import tuning
from planetgen.api import ids, system_sky
from planetgen.api.config import Config
from planetgen.db import sky_sprites, store
from planetgen.galaxy import sky, sky_sprites as sprite_maker
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.web.app import create_app
from planetgen.web.lib import apiclient, skycache, skypack, tilecache
from tests.publicids import pid
from tests.test_sky_sphere import add_bright, add_system
from tests.toy_shape import toy_shape

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
HOME = (8128.0, 0.0, 0.0)
FLAT = lambda radius: np.full(np.shape(radius), 220.0)          # a flat rotation curve, so the stand-in velocities need no galaxy


@pytest.fixture
def conn(mysql_config):
    shape = toy_shape(disk_scale_length_pc=3000.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=400.0, bulge_amplitude=2.0, arm_count=2,
                      pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=2000.0, config=mysql_config,
                            cell_star_ceiling=1e12)
    connection = store.get_connection(mysql_config)
    try:
        yield connection
    finally:
        connection.close()


def add_nebula(conn, point, radius_ly, family_class="D", extinction=None, density=None):
    """A placed nebula at `point` (parsecs) with its radius (and, when given, its dimming and density) pinned; returns its id."""
    nebula = Nebula(SystemConfig(), name=f"Veil {point[0]:.0f}")
    nebula.radius_ly = radius_ly
    nebula_id = store.insert_nebula(conn, nebula, placement={"center_x_pc": point[0], "center_y_pc": point[1], "center_z_pc": point[2], "galactic_radius_pc": 1.0})
    pins = {"nebula_class": family_class, "nebula_type": tuning.NEBULA_CLASSES[family_class]["family"]}
    if extinction is not None:
        pins["extinction_av"] = extinction
    if density is not None:
        pins["density_cm3"] = density
    sets = ", ".join(f"{name} = ?" for name in pins)
    conn.execute(f"UPDATE nebulae SET {sets} WHERE id = ?", (*pins.values(), nebula_id))
    return nebula_id


def add_remnant(conn, point, radius_ly, age_years):
    remnant = SupernovaRemnant(SystemConfig(), name="Old Shell")
    remnant.radius_ly = radius_ly
    remnant.age_years = age_years
    remnant.morphology = "shell"
    return store.insert_supernova_remnant(conn, remnant, placement={"center_x_pc": point[0], "center_y_pc": point[1], "center_z_pc": point[2], "galactic_radius_pc": 1.0})


# --- The objects of a sphere ---------------------------------------------------------------------------------------------------------------------

def test_the_objects_of_a_sphere_are_the_placed_nebulae_and_remnants_inside_it_sorted_by_kind_and_id(conn):
    near = add_nebula(conn, (8128.0, 100.0, 0.0), 40.0)
    far = add_nebula(conn, (8128.0, 5000.0, 0.0), 40.0)
    shell = add_remnant(conn, (8128.0, 0.0, 200.0), 30.0, 5000.0)
    conn.commit()
    found = sky_sprites.objects_in_sphere(conn, HOME, 1000.0)
    assert list(found.kind) == [sky_sprites.NEBULA, sky_sprites.REMNANT]
    assert list(found.ident) == [near, shell] and far not in list(found.ident)
    assert np.allclose(found.position_pc, [(8128.0, 100.0, 0.0), (8128.0, 0.0, 200.0)])
    assert found.family[1] == "supernova-remnant" and found.family[0] == tuning.NEBULA_CLASSES["D"]["family"]
    assert found.age_years[1] == 5000.0 and math.isnan(found.age_years[0])
    both = sky_sprites.objects_in_sphere(conn, HOME, 6000.0)
    assert list(both.ident) == [near, far, shell]
    assert len(sky_sprites.objects_in_sphere(conn, (0.0, 0.0, 1000.0), 10.0).ident) == 0


# --- The sprite records --------------------------------------------------------------------------------------------------------------------------

def objects(rows):
    """A `SpriteObjects` from `(kind, family, point, radius_ly, density, extinction, age)` rows at rest with an epoch of the sky's time."""
    n = len(rows)
    return sky_sprites.SpriteObjects(
        np.array([r[0] for r in rows], dtype=object), np.arange(1, n + 1), np.array([r[1] for r in rows], dtype=object),
        np.array([r[2] for r in rows], dtype=float).reshape(-1, 3), np.full(n, 1.0e9), np.zeros((n, 3)),
        np.array([r[3] for r in rows], dtype=float), np.array([r[4] for r in rows], dtype=float), np.array([r[5] for r in rows], dtype=float),
        np.array([r[6] for r in rows], dtype=float))


OBSERVER = sky.Observer(HOME, 1.0e9, (0.0, 0.0, 0.0), -1)
NO_DUST = {"av_per_kpc_at_sun": 0.0}


def test_a_bright_nebula_is_a_record_with_the_fields_the_browser_reads():
    rows = [("nebula", "emission", (8128.0 + 500.0, 0.0, 0.0), 100.0, 300.0, 0.4, math.nan)]
    records, report = sprite_maker.sprites_from(OBSERVER, objects(rows), 1.0e9, "g", FLAT, NO_DUST)
    assert report["in"] == 1 and report["sent"] == 1
    (record,) = records
    assert set(record) == {"id", "kind", "family", "direction", "distanceLy", "angularRadiusRad", "inside", "extinctionMag", "surfaceBrightness", "visible", "seen"}
    assert record["id"] == "nebula:1" and record["family"] == "emission" and record["kind"] == "nebula"
    assert record["direction"] == pytest.approx([1.0, 0.0, 0.0], abs=1e-6)
    assert record["distanceLy"] == pytest.approx(pc_to_ly(500.0), rel=1e-3)
    assert record["angularRadiusRad"] == pytest.approx(math.asin(100.0 / pc_to_ly(500.0)), rel=1e-3)
    assert record["inside"] is False and record["visible"] is True and record["seen"] is True
    assert isinstance(record["surfaceBrightness"], float) and record["surfaceBrightness"] < 23.0
    assert record["extinctionMag"] == pytest.approx(0.4 * math.sqrt(1.0 - (0.0 / 100.0) ** 2), rel=0.1)    # the chord through the middle


def test_a_dark_cloud_that_dims_is_sent_without_a_surface_brightness_and_a_faint_clear_one_is_not():
    rows = [("nebula", "dark", (8128.0, 400.0, 0.0), 80.0, 100.0, 2.0, math.nan),
            ("nebula", "diffuse", (8128.0, -400.0, 0.0), 80.0, 1.0, 0.0, math.nan)]
    records, report = sprite_maker.sprites_from(OBSERVER, objects(rows), 1.0e9, "g", FLAT, NO_DUST)
    assert [r["id"] for r in records] == ["nebula:1"]
    assert records[0]["surfaceBrightness"] is None and records[0]["visible"] is False and records[0]["extinctionMag"] > 1.0
    assert report["too_faint"] == 1 and report["sent"] == 1


def test_a_sprite_too_small_to_see_is_left_out_and_counted():
    rows = [("nebula", "emission", (8128.0 + 9000.0, 0.0, 0.0), 5.0, 300.0, 0.4, math.nan)]
    records, report = sprite_maker.sprites_from(OBSERVER, objects(rows), 1.0e9, "g", FLAT, NO_DUST)
    assert records == [] and report["too_small"] == 1


def test_a_remnant_whose_explosion_has_not_been_seen_is_left_out():
    point = (8128.0 + 100.0, 0.0, 0.0)                                  # 326 light years away: its light takes 326 years
    young = ("supernova_remnant", "supernova-remnant", point, 20.0, 100.0, 0.0, 100.0)
    old = ("supernova_remnant", "supernova-remnant", (8128.0, 100.0, 0.0), 20.0, 100.0, 0.0, 9000.0)
    records, report = sprite_maker.sprites_from(OBSERVER, objects([young, old]), 1.0e9, "g", FLAT, NO_DUST)
    assert report["not_seen"] == 1 and [r["id"] for r in records] == ["supernova_remnant:2"] and records[0]["seen"] is True


def test_an_observer_inside_a_nebula_is_told_so():
    rows = [("nebula", "emission", (8128.0 + 1.0, 0.0, 0.0), 100.0, 300.0, 0.4, math.nan)]
    records, report = sprite_maker.sprites_from(OBSERVER, objects(rows), 1.0e9, "g", FLAT, NO_DUST)
    assert report["inside"] == 1 and records[0]["inside"] is True and records[0]["angularRadiusRad"] == pytest.approx(math.pi)


def test_an_object_with_no_stored_velocity_or_epoch_still_gets_a_place_and_the_same_one_each_time():
    n = objects([("nebula", "emission", (8128.0 + 500.0, 0.0, 0.0), 100.0, 300.0, 0.4, math.nan)])
    n = n._replace(velocity_kms=np.full((1, 3), math.nan), epoch_unix=np.full(1, math.nan))
    first, _ = sprite_maker.sprites_from(OBSERVER, n, 1.0e9, "g", FLAT, NO_DUST)
    second, _ = sprite_maker.sprites_from(OBSERVER, n, 1.0e9, "g", FLAT, NO_DUST)
    assert first == second and len(first) == 1


def test_an_empty_sphere_gives_no_records():
    assert sprite_maker.sprites_from(OBSERVER, sky_sprites.empty(), 1.0e9) == ([], {"in": 0, "sent": 0, "too_small": 0, "too_faint": 0, "not_seen": 0, "inside": 0})


# --- The sky of a system ---------------------------------------------------------------------------------------------------------------------------

def a_neighbourhood(conn):
    conn.execute("UPDATE galaxy_shape SET bright_star_min_luminosity_sol = 800.0, bright_star_seed = 11 WHERE id = 1")
    own, _ = add_system(conn, HOME, [(1.0, 5772.0)], velocity=(0.0, 0.0, 0.0), name="Own")
    near, near_stars = add_system(conn, (8138.0, 0.0, 0.0), [(1.0, 5772.0)], velocity=(0.0, 0.0, 0.0), name="Near")
    bright = add_bright(conn, (8128.0, 60.0, 0.0), 2000.0)
    nebula = add_nebula(conn, (8128.0, 0.0, 400.0), 90.0, "D", density=300.0)
    conn.commit()
    pid("system", own, conn)                                    # a row saved by hand has no printed ID until it is asked for
    pid("nebula", nebula, conn)
    return own, near_stars[0], bright, nebula


def test_a_system_sky_holds_the_stars_it_sees_the_nebulae_and_the_report(conn):
    own, near_star, bright, nebula = a_neighbourhood(conn)
    built = system_sky.build_system_sky(conn, own, 7.5, tile_stamp="stamp-one")
    back = skypack.unpack_sky(built["packed"])
    assert built["count"] == back["count"] >= 2 and built["limit"] == 7.5 and built["version_key"] == skypack.version_key("stamp-one")
    assert back["limit"] == pytest.approx(7.5) and back["epoch_unix"] == built["epoch_unix"] == 1.0e9
    kinds = {v: k for k, v in skypack.KIND_CODES.items()}
    assert sorted(kinds[k] for k in back["kind"]).count("charted") == 1 and "bright_star" in {kinds[k] for k in back["kind"]}
    meta = back["meta"]
    assert meta["system"] == built["system"] == pid("system", own) and meta["version_key"] == built["version_key"] and meta["limit"] == 7.5
    assert [n["id"] for n in meta["nebulae"]] == [pid("nebula", nebula)] and built["sprites"] == 1
    report = meta["report"]
    assert report["visible"] == built["count"] and report["excluded_own_system"] == 1 and report["limit"] == 7.5
    assert set(report) >= {"stars_in", "tiers", "built_to_pc", "complete_above_sol", "wanted_pc", "filling", "sprites", "unbuilt_velocity"}
    assert report["sprites"]["sent"] == 1 and isinstance(report["filling"], bool) and report["wanted_pc"] > 0
    assert built["report"] == report


def test_the_charted_stars_of_a_sky_are_named_by_their_place_in_the_list_for_the_page_to_link(conn):
    """VIEW.37: `meta.charted` maps the place of each charted star to its printed star ID and name; no other tier is in it."""
    own, near_star, bright, _nebula = a_neighbourhood(conn)
    store.assign_uids(conn, system_ids=[conn.execute("SELECT star_system_id FROM stars WHERE id = ?", (near_star,)).fetchone()["star_system_id"]])
    conn.commit()                                               # a row saved by hand has no printed ID until it is issued one
    back = skypack.unpack_sky(system_sky.build_system_sky(conn, own, 7.5, "s")["packed"])
    charted = back["meta"]["charted"]
    places = [i for i, kind in enumerate(back["kind"]) if kind == skypack.KIND_CODES["charted"]]
    assert sorted(charted) == sorted(str(back["index"][i]) for i in places) and len(charted) == 1
    printed, name = next(iter(charted.values()))
    assert printed == ids.printed(conn, "star", near_star) and name == conn.execute("SELECT name FROM stars WHERE id = ?", (near_star,)).fetchone()["name"]


def test_a_sky_is_the_same_bytes_every_time_and_a_dimmer_limit_has_fewer_stars(conn):
    own, *_ = a_neighbourhood(conn)
    first = system_sky.build_system_sky(conn, own, 7.5, "s")
    second = system_sky.build_system_sky(conn, own, 7.5, "s")
    assert first["packed"] == second["packed"]
    shallow = system_sky.build_system_sky(conn, own, 3.0, "s")
    assert shallow["count"] < first["count"] and shallow["limit"] == 3.0


def test_a_sphere_of_too_many_sectors_has_its_cell_levels_left_unlisted(conn):
    from planetgen.db import sky_sphere
    conn.execute("UPDATE galaxy_shape SET bright_star_min_luminosity_sol = 800.0, bright_star_seed = 11 WHERE id = 1")
    small = sky_sphere.cell_levels(conn, HOME, 20.0)
    big = sky_sphere.cell_levels(conn, HOME, 40000.0)
    assert not small.truncated and len(small.cells) > 0
    assert big.truncated and big.cells == () and big.complete_above_sol == big.scatter_level_sol == 800.0


def test_a_system_with_no_place_or_no_row_is_refused(conn):
    config_id = conn.execute("INSERT INTO system_configs (name) VALUES (NULL)").lastrowid
    unplaced = conn.execute("INSERT INTO star_systems (system_config_id, name) VALUES (?, 'Nowhere')", (config_id,)).lastrowid
    conn.commit()
    with pytest.raises(system_sky.SkyError):
        system_sky.build_system_sky(conn, unplaced)
    with pytest.raises(system_sky.SkyError):
        system_sky.build_system_sky(conn, 999999)


def test_the_limit_is_checked_and_defaults_to_the_catalogues():
    assert skypack.DEFAULT_LIMIT == sky.DEFAULT_LIMIT
    assert skypack.check_limit(None) == skypack.check_limit("") == 7.5 and skypack.check_limit("6.25") == 6.25 and skypack.check_limit(8) == 8.0
    for bad in ("x", "nan", "inf", "1.99", "9.01", "-3"):
        with pytest.raises(skypack.SkyLimitError):
            skypack.check_limit(bad)
    assert skypack.version_key("a") == skypack.version_key("a") != skypack.version_key("b") and len(skypack.version_key("a")) == 16


# --- The API route -----------------------------------------------------------------------------------------------------------------------------

@pytest.fixture
def client(mysql_config):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def test_the_api_route_answers_the_packed_sky_as_base64(conn, client):
    own, *_ = a_neighbourhood(conn)
    answer = client.get(f"/api/systems/{pid('system', own)}/sky?limit=6")
    assert answer.status_code == 200
    body = answer.get_json()
    raw = base64.b64decode(body["data"])
    assert body["format"] == "PGK1" and skypack.is_whole(raw) and body["limit"] == 6.0 and body["count"] == skypack.unpack_sky(raw)["count"]
    assert body["system"] == pid("system", own) and len(body["version_key"]) == 16 and body["report"]["limit"] == 6.0
    assert client.get(f"/api/systems/{pid('system', own)}/sky").get_json()["limit"] == 7.5


def test_the_api_route_refuses_a_limit_out_of_range_and_an_unknown_system(conn, client):
    own, *_ = a_neighbourhood(conn)
    assert client.get(f"/api/systems/{pid('system', own)}/sky?limit=12").status_code == 400
    assert client.get(f"/api/systems/{pid('system', own)}/sky?limit=abc").status_code == 400
    assert client.get("/api/systems/999999/sky").status_code == 404


def test_a_system_band_is_the_glow_image_of_its_cell_and_the_sky_names_it(conn, monkeypatch):
    own, *_ = a_neighbourhood(conn)
    monkeypatch.setattr(system_sky, "_bands", system_sky.sky_band.MemoryBandStore())
    band = system_sky.build_system_band(conn, own, 7.5)
    assert (band["width"], band["height"], band["format"]) == (256, 128, "rgb16f") and len(band["data"]) == 256 * 128 * 6 and len(band["key"]) == 32
    image = np.frombuffer(band["data"], dtype="<f2").reshape(128, 256, 3).astype(float)
    assert np.isfinite(image).all() and image.max() > 0.0 and (image[..., 0] >= image[..., 2]).all()          # warm: red over blue
    assert band["cell"] == list(system_sky.sky_band.cell_of(HOME))
    again = system_sky.build_system_band(conn, own, 7.5)
    assert again["data"] == band["data"] and again["key"] == band["key"]
    assert system_sky.build_system_band(conn, own, 5.0)["key"] != band["key"]
    meta = skypack.unpack_sky(system_sky.build_system_sky(conn, own, 7.5, "s")["packed"])["meta"]
    assert meta["band"] == {"key": band["key"], "width": 256, "height": 128, "format": "rgb16f", "href": f"/system/{pid('system', own)}/sky/band?limit=7.5"}


def test_the_band_route_answers_the_half_floats_as_base64_and_checks_the_limit(conn, client, monkeypatch):
    own, *_ = a_neighbourhood(conn)
    monkeypatch.setattr(system_sky, "_bands", system_sky.sky_band.MemoryBandStore())
    answer = client.get(f"/api/systems/{pid('system', own)}/sky/band?limit=6")
    assert answer.status_code == 200
    body = answer.get_json()
    assert len(base64.b64decode(body["data"])) == 256 * 128 * 6 and body["limit"] == 6.0 and body["format"] == "rgb16f" and "seconds" in body
    assert client.get(f"/api/systems/{pid('system', own)}/sky/band?limit=1").status_code == 400
    assert client.get("/api/systems/999999/sky/band").status_code == 404


# --- The disk cache and the web route --------------------------------------------------------------------------------------------------------

def packed_sky(count=40, seed=3):
    from tests.test_skypack import a_sky
    result, _ = a_sky(count, seed)
    return result


class FakeApi:
    """`apiclient.get_system_sky` over a packed sky, counting its calls; `stamp` is the version it builds for."""

    def __init__(self):
        self.calls = []
        self.stamp = "aaaa"

    def __call__(self, db, system_id, limit=None):
        self.calls.append((db, system_id, limit))
        key = skypack.version_key(self.stamp)
        meta = {"system": system_id, "version_key": key, "limit": limit, "nebulae": [{"id": "n"}], "report": {"filling": False, "visible": 1}}
        raw = skypack.pack_sky(packed_sky(), 1.0e9, limit, meta=meta)
        return {"system": system_id, "limit": limit, "version_key": key, "data": base64.b64encode(raw).decode("ascii")}


@pytest.fixture
def cache(tmp_path, monkeypatch):
    api = FakeApi()
    state = {"stamp": "aaaa"}
    monkeypatch.setattr(tilecache, "cache_dir", lambda: str(tmp_path))
    monkeypatch.setattr(tilecache, "current_stamp", lambda db, root=None: {"stamp": api.stamp, "generation": "gen1", "history": []})
    monkeypatch.setattr(tilecache, "_read_remembered", lambda db_dir: {"stamp": api.stamp})
    monkeypatch.setattr(tilecache, "PRUNE_PROBABILITY", 0.0)
    monkeypatch.setattr(apiclient, "get_system_sky", api)
    return api, tmp_path


def sky_files(root):
    return sorted(name for _, _, names in os.walk(root) for name in names if name.startswith(skycache.NAME_PREFIX))


def test_a_sky_is_built_once_and_then_served_from_the_disk(cache):
    api, root = cache
    raw, header = skycache.fetch_sky("db", "FE81-0000005-000", "7.5")
    assert header["source"] == "built" and header["count"] > 0 and header["sprites"] == 1 and header["version_key"] == skypack.version_key("aaaa")
    again, header2 = skycache.fetch_sky("db", "FE81-0000005-000", 7.5)
    assert again == raw and header2["source"] == "cache" and header2["report"] == {"filling": False, "visible": 1} and len(api.calls) == 1
    assert len(sky_files(root)) == 1


def test_a_sky_per_limit_and_a_new_version_key_replaces_the_old_copy(cache):
    api, root = cache
    skycache.fetch_sky("db", "FE81-0000005-000", 7.5)
    skycache.fetch_sky("db", "FE81-0000005-000", 6.0)
    assert len(sky_files(root)) == 2 and len(api.calls) == 2
    api.stamp = "bbbb"                                                 # the galaxy changed: the stored copies are no longer for it
    _, header = skycache.fetch_sky("db", "FE81-0000005-000", 7.5)
    assert header["source"] == "built" and header["version_key"] == skypack.version_key("bbbb")
    assert len(sky_files(root)) == 2                                     # the 7.5 copy was replaced, the 6.0 one is stale but kept until its own rebuild
    assert any(skypack.version_key("aaaa") in name for name in sky_files(root))


def test_a_torn_file_is_built_again_and_the_clear_removes_every_sky(cache):
    api, root = cache
    skycache.fetch_sky("db", "FE81-0000005-000", 7.5)
    (path,) = [os.path.join(d, n) for d, _, names in os.walk(root) for n in names if n.startswith(skycache.NAME_PREFIX)]
    with open(path, "r+b") as handle:
        handle.truncate(30)
    _, header = skycache.fetch_sky("db", "FE81-0000005-000", 7.5)
    assert header["source"] == "built" and len(api.calls) == 2
    assert skycache.clear() == 1 and sky_files(root) == []


def test_a_cache_that_is_off_still_answers(monkeypatch):
    api = FakeApi()
    monkeypatch.setattr(tilecache, "cache_dir", lambda: None)
    monkeypatch.setattr(apiclient, "get_system_sky", api)
    _, header = skycache.fetch_sky("db", "FE81-0000005-000", 7.5)
    _, again = skycache.fetch_sky("db", "FE81-0000005-000", 7.5)
    assert header["source"] == again["source"] == "built" and len(api.calls) == 2


def test_the_web_route_sends_the_packed_list_and_says_where_it_came_from(cache):
    api, root = cache

    class SiteConfig(Config):
        WEB_DATABASE = "db"
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "sky-secret"

    app = create_app(SiteConfig)
    app.testing = True
    web = app.test_client()
    first = web.get("/system/FE81-0000005-000/sky?limit=6")
    assert first.status_code == 200 and first.mimetype == skypack.MIMETYPE and skypack.is_whole(first.data)
    assert first.headers["X-Sky-Source"] == "built" and first.headers["X-Sky-Version"] == skypack.version_key("aaaa") and first.headers["Cache-Control"] == "no-store"
    assert web.get("/system/FE81-0000005-000/sky?limit=6").headers["X-Sky-Source"] == "cache"
    assert web.get("/system/FE81-0000005-000/sky?limit=12").status_code == 400
    assert web.get("/system/FE81-0000005-000/sky?limit=12").get_json()["error"]

    def missing(db, system_id, limit=None):
        raise apiclient.NotFoundError("no such system")

    def broken(db, system_id, limit=None):
        raise apiclient.ApiError("down")

    import pytest as _pytest
    with _pytest.MonkeyPatch.context() as patch:
        patch.setattr(apiclient, "get_system_sky", missing)
        assert web.get("/system/FE81-0000005-000/sky?limit=5").status_code == 404
        patch.setattr(apiclient, "get_system_sky", broken)
        assert web.get("/system/FE81-0000005-000/sky?limit=5").status_code == 502


class FakeBandApi:
    """`apiclient.get_system_sky_band` over a flat image, counting its calls."""

    KEY = "ab" * 16

    def __init__(self):
        self.calls = []

    def __call__(self, db, system_id, limit=None):
        self.calls.append((db, system_id, limit))
        data = np.full((128, 256, 3), 1.0, dtype="<f2").tobytes()
        return {"key": self.KEY, "limit": limit, "cell": [0, 0, 0], "width": 256, "height": 128, "format": "rgb16f", "seconds": 0.1,
                "data": base64.b64encode(data).decode("ascii")}


def band_files(root):
    return sorted(name for _, _, names in os.walk(root) for name in names if name.startswith(skycache.BAND_PREFIX))


def test_a_band_is_built_once_stored_by_its_key_and_served_from_the_disk(cache, monkeypatch):
    api, root = cache
    band_api = FakeBandApi()
    monkeypatch.setattr(apiclient, "get_system_sky_band", band_api)
    raw, header = skycache.fetch_band("db", "FE81-0000005-000", "7.5")
    assert header == {"key": band_api.KEY, "width": 256, "height": 128, "format": "rgb16f", "source": "built"} and len(raw) == 256 * 128 * 6
    assert band_files(root) == [f"b{band_api.KEY}.bin"]
    again, header2 = skycache.fetch_band("db", "FE81-0000005-000", 7.5, band_api.KEY)
    assert again == raw and header2["source"] == "cache" and len(band_api.calls) == 1
    assert skycache.clear() == 1 and band_files(root) == []                     # the clear removes the glow images with the skies


def test_a_band_key_that_is_not_a_hash_or_a_torn_file_is_not_trusted(cache, monkeypatch):
    api, root = cache
    band_api = FakeBandApi()
    monkeypatch.setattr(apiclient, "get_system_sky_band", band_api)
    with pytest.raises(ValueError):
        skycache.fetch_band("db", "FE81-0000005-000", 7.5, "../../etc/passwd")
    skycache.fetch_band("db", "FE81-0000005-000", 7.5)
    (path,) = [os.path.join(d, n) for d, _, names in os.walk(root) for n in names if n.startswith(skycache.BAND_PREFIX)]
    with open(path, "r+b") as handle:
        handle.truncate(100)
    _, header = skycache.fetch_band("db", "FE81-0000005-000", 7.5, band_api.KEY)
    assert header["source"] == "built" and len(band_api.calls) == 2


def test_the_web_route_sends_the_glow_image_with_its_size_in_headers(cache, monkeypatch):
    api, root = cache
    band_api = FakeBandApi()
    monkeypatch.setattr(apiclient, "get_system_sky_band", band_api)

    class SiteConfig(Config):
        WEB_DATABASE = "db"
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "sky-secret"

    app = create_app(SiteConfig)
    app.testing = True
    web = app.test_client()
    first = web.get("/system/FE81-0000005-000/sky/band?limit=6")
    assert first.status_code == 200 and len(first.data) == 256 * 128 * 6 and first.mimetype == "application/octet-stream"
    assert first.headers["X-Sky-Source"] == "built" and first.headers["X-Sky-Band-Key"] == band_api.KEY and first.headers["X-Sky-Band-Format"] == "rgb16f"
    assert (first.headers["X-Sky-Band-Width"], first.headers["X-Sky-Band-Height"]) == ("256", "128")
    assert web.get(f"/system/FE81-0000005-000/sky/band?limit=6&key={band_api.KEY}").headers["X-Sky-Source"] == "cache"
    assert web.get("/system/FE81-0000005-000/sky/band?limit=12").status_code == 400
    assert web.get("/system/FE81-0000005-000/sky/band?limit=6&key=zz").status_code == 400
