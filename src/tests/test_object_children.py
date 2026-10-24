# tests/test_object_children.py

"""
UX.169: `GET /api/objects/<ref>/children` (`db/children.py`) pages through what belongs to an object: a sector's
systems, a system's stars, a star's planets and belts then its comets, a planet's moons; the counts match the tables.
"""

from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from tests.publicids import pid
from tests.test_api import _save_wide_binary_with_moons, client, seeded_sector  # noqa: F401


def _query(mysql_config, sql, params=()):
    conn = store.get_connection(mysql_config)
    try:
        return [dict(row) for row in conn.execute(sql, params).fetchall()]
    finally:
        conn.close()


def _page(client, ref, **args):
    response = client.get(f"/api/objects/{ref}/children", query_string=args)
    assert response.status_code == 200, response.get_json()
    return response.get_json()


def _all(client, ref, limit):
    """Every item of an object's children by following `next`, and how many pages that took."""
    items, token, pages = [], None, 0
    while True:
        page = _page(client, ref, limit=limit, **({"after": token} if token else {}))
        items += page["items"]
        pages += 1
        token = page["next"]
        if token is None:
            return items, pages


def _ref(kind, row_id):
    return f"{kind}:{pid(kind, row_id)}"


def test_a_wide_pair_holds_its_stars_and_each_star_holds_its_planets_and_their_moons(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    stars = _query(mysql_config, "SELECT id FROM stars WHERE star_system_id = ? ORDER BY COALESCE(slot, 0), id", (system_id,))
    system = _page(client, _ref("system", system_id))
    assert system["kind"] == "system" and system["total"] == len(stars) and system["next"] is None
    assert [item["id"] for item in system["items"]] == [pid("star", star["id"]) for star in stars]
    assert [item["ref"] for item in system["items"]] == [_ref("star", star["id"]) for star in stars]
    assert all(item["kind"] == "star" and item["name"] and item["summary"]["star_type"] for item in system["items"])

    planets_total = moons_total = 0
    for star in stars:
        planets = _query(mysql_config, "SELECT id FROM planets WHERE parent_type = 'star' AND parent_id = ? "
                                       "ORDER BY orbital_index, id", (star["id"],))
        belts = _query(mysql_config, "SELECT id FROM asteroid_belts WHERE parent_type = 'star' AND parent_id = ?",
                       (star["id"],))
        comets = _query(mysql_config, "SELECT id FROM comets WHERE parent_type = 'star' AND parent_id = ?", (star["id"],))
        page = _page(client, _ref("star", star["id"]))
        assert page["total"] == len(planets) + len(belts) + len(comets)
        kinds = [item["kind"] for item in page["items"]]
        assert kinds.count("planet") == len(planets) and kinds.count("belt") == len(belts)
        assert kinds.count("comet") == len(comets)
        assert kinds == sorted(kinds, key=lambda kind: kind == "comet"), "comets come last"
        listed_planets = [item["id"] for item in page["items"] if item["kind"] == "planet"]
        assert listed_planets == [pid("planet", planet["id"]) for planet in planets]
        for item in page["items"]:
            if item["kind"] != "planet":
                continue
            moons = _query(mysql_config, "SELECT id FROM moons WHERE parent_type = 'planet' AND parent_id = ?"
                                         " ORDER BY orbital_index, id", (_row_of("planet", item["id"], mysql_config),))
            moons_total += len(moons)
            assert item["child_count"] == len(moons) and item["counts"].get("moon", 0) == len(moons)
            held = _page(client, item["ref"])
            assert [m["id"] for m in held["items"]] == [pid("moon", m["id"]) for m in moons] and held["total"] == len(moons)
            assert all(m["kind"] == "moon" and m["child_count"] == 0 for m in held["items"])
        planets_total += len(planets)
        star_item = next(item for item in system["items"] if item["id"] == pid("star", star["id"]))
        assert star_item["child_count"] == page["total"]
        assert star_item["counts"].get("planet", 0) == len(planets)
    assert moons_total >= 1, "the system was drawn with a moon"
    counts = _query(mysql_config, "SELECT (SELECT COUNT(*) FROM planets WHERE star_system_id = ?) AS planets, "
                                  "(SELECT COUNT(*) FROM moons WHERE star_system_id = ?) AS moons", (system_id, system_id))[0]
    assert (planets_total, moons_total) == (counts["planets"], counts["moons"])


def _row_of(kind, printed, mysql_config):
    from planetgen.api import ids

    conn = store.get_connection(mysql_config)
    try:
        return ids.row_id(conn, kind, printed)
    finally:
        conn.close()


def test_a_page_of_one_follows_next_to_the_same_list(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    star_id = _query(mysql_config, "SELECT id FROM stars WHERE star_system_id = ? ORDER BY COALESCE(slot, 0), id",
                     (system_id,))[0]["id"]
    ref = _ref("star", star_id)
    whole = _page(client, ref)["items"]
    paged, pages = _all(client, ref, 1)
    assert [item["ref"] for item in paged] == [item["ref"] for item in whole] and pages == max(len(whole), 1)
    assert _page(client, ref, limit=len(whole))["next"] is None


def _sector_with_a_wide_pair(mysql_config):
    """A sector holding one wide binary with a planet that has a moon (generation is random, so retried); returns
    `(sector_id, system_id)`."""
    for _ in range(60):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = True
        system = StarSystem(system_config=cfg)
        if any(p.body_type != "a" and p.moons for p in system.planets):
            sector = SpaceSector("Children Sector", edge_ly=10.0)
            sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
            sector_id = store.save_sector(sector, config=mysql_config)
            return sector_id, _query(mysql_config, "SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,))[0]["id"]
    raise AssertionError("could not generate a wide binary with a moon")


def test_a_system_counts_what_lies_below_it(client, mysql_config):
    sector_id, system_id = _sector_with_a_wide_pair(mysql_config)
    counts = _query(mysql_config, "SELECT (SELECT COUNT(*) FROM stars WHERE star_system_id = ?) AS stars, "
                                  "(SELECT COUNT(*) FROM planets WHERE star_system_id = ?) AS planets, "
                                  "(SELECT COUNT(*) FROM moons WHERE star_system_id = ?) AS moons, "
                                  "(SELECT COUNT(*) FROM asteroid_belts WHERE star_system_id = ?) AS belts, "
                                  "(SELECT COUNT(*) FROM comets WHERE star_system_id = ?) AS comets", (system_id,) * 5)[0]
    page = _page(client, f"sector:{pid('sector', sector_id)}")
    item = next(item for item in page["items"] if item["id"] == pid("system", system_id))
    assert item["counts"] == {kind: counts[plural] for kind, plural in (
        ("star", "stars"), ("planet", "planets"), ("moon", "moons"), ("belt", "belts"), ("comet", "comets"))
        if counts[plural]}
    assert item["child_count"] == counts["stars"] and item["summary"]["is_binary"] is True
    assert len(item["summary"]["star_types"]) == counts["stars"]


def test_a_sector_holds_its_systems(client, seeded_sector):
    config, sector_id, system_ids = seeded_sector
    page = _page(client, f"sector:{pid('sector', sector_id)}")
    assert page["total"] == len(system_ids) == len(page["items"])
    assert [item["id"] for item in page["items"]] == [pid("system", system_id) for system_id in sorted(system_ids)]
    assert all(item["kind"] == "system" and item["summary"]["star_types"] for item in page["items"])
    assert sorted(item["summary"]["star_types"][0].split()[0] for item in page["items"]) == ["G2V", "M5V"]
    assert all(item["counts"]["star"] == 1 and item["child_count"] == 1 for item in page["items"])


def test_a_sector_with_thousands_of_systems_pages_without_gaps(client, seeded_sector):
    config, sector_id, system_ids = seeded_sector
    rows = _query(config, "SELECT system_config_id FROM star_systems WHERE id = ?", (system_ids[0],))
    conn = store.get_connection(config)
    try:
        conn.execute("SET SESSION cte_max_recursion_depth = 5000")
        conn.execute(
            "INSERT INTO star_systems (sector_id, system_config_id, name) "
            "WITH RECURSIVE n(i) AS (SELECT 1 UNION ALL SELECT i + 1 FROM n WHERE i < 2500) "
            "SELECT ?, ?, CONCAT('Bulk ', i) FROM n", (sector_id, rows[0]["system_config_id"]))
        conn.commit()
    finally:
        conn.close()
    ref = f"sector:{pid('sector', sector_id)}"
    first = _page(client, ref, limit=200)
    assert first["total"] == len(system_ids) + 2500 and len(first["items"]) == 200 and first["next"]
    items, pages = _all(client, ref, 200)
    assert len(items) == first["total"] and pages == -(-first["total"] // 200)
    assert len({item["id"] for item in items}) == len(items), "no system twice"
    assert all(item["child_count"] in (0, 1) for item in items)


def test_a_moon_a_belt_and_a_comet_hold_nothing(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    moon = _query(mysql_config, "SELECT id FROM moons WHERE star_system_id = ? LIMIT 1", (system_id,))[0]["id"]
    page = _page(client, _ref("moon", moon))
    assert page["items"] == [] and page["total"] == 0 and page["next"] is None and page["kind"] == "moon"


def test_a_missing_object_is_a_404_and_a_bad_page_a_400(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    for ref in ("moon:FFFFFFFFFF-0000005-099", "ship:1", "planet"):
        assert client.get(f"/api/objects/{ref}/children").status_code == 404, ref
    ref = _ref("system", system_id)
    for args in ({"limit": 0}, {"limit": 201}, {"after": "nope"}, {"after": "9.0.0.1"}, {"after": "-1.0.0.1"}):
        response = client.get(f"/api/objects/{ref}/children", query_string=args)
        assert response.status_code == 400, args
        assert isinstance(response.get_json()["error"], str)


def test_a_sector_lists_its_phenomena_after_its_systems(client, seeded_sector):
    from planetgen.generation.phenomena.nebula import Nebula

    config, sector_id, system_ids = seeded_sector
    conn = store.get_connection(config)
    try:
        nebula_id = store.insert_nebula(conn, Nebula(SystemConfig()), sector_id=sector_id, placement={
            "center_x_pc": 5.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 5.5})
        conn.commit()
    finally:
        conn.close()
    printed = pid("nebula", nebula_id)   # a row inserted by hand gets its Object ID here
    ref = f"sector:{pid('sector', sector_id)}"
    items, _pages = _all(client, ref, 1)
    assert [item["kind"] for item in items] == ["system"] * len(system_ids) + ["nebula"]
    assert items[-1]["id"] == printed and items[-1]["name"] and items[-1]["child_count"] == 0
    assert _page(client, ref)["total"] == len(system_ids) + 1
