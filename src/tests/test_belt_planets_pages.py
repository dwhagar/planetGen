# tests/test_belt_planets_pages.py

"""
GEN.1511: the pages of a belt that made planets (GEN.1510): the belt's information page lists the planets it made, the designation of a body made a planet leads to that planet, and a member that stayed has its page as before.
"""

import pytest

from planetgen.db import store
from planetgen.describe import belt_info, belt_member_info, sheet as sheet_module
from planetgen.galaxy import object_uid
from planetgen.generation import irregular_moons
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


def _config():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M3V"
    cfg.PLANETS = True
    cfg.MOONS = False
    cfg.BINARY_SYSTEM = False
    return cfg


@pytest.fixture
def saved_belt(mysql_config):
    """The id of a stored belt that made planets, with the planets it made."""
    config = _config()
    for seed in range(1, 30):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=config)
        made = [body for body in system.planets if getattr(body, "belt_rank", None) is not None]
        if not made:
            continue
        saved = store.save_system(system, config, config=mysql_config)
        with store.get_connection(mysql_config) as conn:
            loaded = store.load_star_system(conn, saved)
            order = next(b for b in loaded.planets if getattr(b, "belt_rank", None) is not None).belt_order
            belt = list(irregular_moons.belts_of(loaded.planets).items())
            belt = dict(belt)[order]
            return belt.db_id, [p for p in loaded.planets if getattr(p, "belt_rank", None) is not None], mysql_config
    pytest.fail("no belt that made planets in 29 systems")


def test_the_belt_page_lists_the_planets_it_made(saved_belt):
    belt_id, made, config = saved_belt
    with store.get_connection(config) as conn:
        page = belt_info.belt_info(conn, belt_id)
    listed = page["planets"]
    assert listed and [row["rank"] for row in listed] == sorted(row["rank"] for row in listed)
    assert {row["rank"] for row in listed} <= {p.belt_rank for p in made}
    first = listed[0]
    assert first["ref"].startswith("planet:") and first["planet_class"] and first["mass"]["unit"] == "kg" and first["designation"].count(".") == 1
    assert not {row["rank"] for row in listed} & {body["rank"] for body in page["bodies"]}          # a body is a member or a planet, not both
    assert abs(page["totals"]["difference"]["value"]) <= 1e-6 * page["totals"]["belt_mass"]["value"]


def test_the_designation_of_a_body_made_a_planet_leads_to_that_planet(saved_belt):
    belt_id, _made, config = saved_belt
    with store.get_connection(config) as conn:
        listed = belt_info.belt_info(conn, belt_id)["planets"][0]
        answer = belt_member_info.member_info(conn, belt_id, listed["rank"])
    assert answer["promoted"]["ref"] == listed["ref"] and answer["designation"] == listed["designation"] and "member" not in answer


def test_a_member_that_stayed_still_has_its_page(saved_belt):
    belt_id, _made, config = saved_belt
    with store.get_connection(config) as conn:
        page = belt_info.belt_info(conn, belt_id)
        if not page["bodies"]:
            pytest.skip("the belt kept no member")
        answer = belt_member_info.member_info(conn, belt_id, page["bodies"][0]["rank"])
    assert "member" in answer and "promoted" not in answer


def test_a_planet_made_from_a_belt_names_that_belt_on_its_sheet(saved_belt):
    """GEN.1511: the planet's own page keeps the belt it came from as its parent, as the belt's page lists the planet."""
    belt_id, made, config = saved_belt
    with store.get_connection(config) as conn:
        planet = made[0]
        belt_uid = object_uid.from_bytes(conn.execute("SELECT uid FROM asteroid_belts WHERE id = ?", (belt_id,)).fetchone()["uid"])
        found = sheet_module.sheet(conn, "planet", planet.db_id)
        rows = [row for group in found["groups"] for row in group["rows"] if row["key"] == sheet_module.BELT_ORIGIN_KEY]
        kept = next(b for b in store.load_star_system(conn, conn.execute("SELECT star_system_id FROM planets WHERE id = ?", (planet.db_id,)).fetchone()["star_system_id"]).planets if getattr(b, "belt_rank", None) is None and b.body_type != "a")
        plain = sheet_module.sheet(conn, "planet", kept.db_id)
    assert len(rows) == 1 and object_uid.designation(belt_uid, planet.belt_rank) in rows[0]["text"]
    assert not [row for group in plain["groups"] for row in group["rows"] if row["key"] == sheet_module.BELT_ORIGIN_KEY]


def test_the_system_scene_names_the_belt_a_planet_was_made_from(saved_belt):
    """GEN.1511: the planet's own record in the system scene (the API's planet record) names its belt; a planet that was not made from a belt has no such key."""
    from planetgen.web.maps import systemscene
    belt_id, made, config = saved_belt
    with store.get_connection(config) as conn:
        system_id = conn.execute("SELECT star_system_id FROM asteroid_belts WHERE id = ?", (belt_id,)).fetchone()["star_system_id"]
        belt_uid = object_uid.from_bytes(conn.execute("SELECT uid FROM asteroid_belts WHERE id = ?", (belt_id,)).fetchone()["uid"])
        scene = systemscene.build_scene(conn, system_id)
    by_name = {planet["name"]: planet for planet in scene["planets"]}
    with store.get_connection(config) as conn:
        listed = {row["name"]: row for row in belt_info.belt_info(conn, belt_id)["planets"]}      # the planets of this belt (the fixture's `made` may hold another belt's)
    assert listed
    for name, row in listed.items():
        found = by_name[name]["made_from_belt"]
        assert found == {"ref": f"belt:{belt_id}", "rank": row["rank"], "designation": object_uid.designation(belt_uid, row["rank"])}
    assert all("made_from_belt" not in planet for planet in scene["planets"] if planet["name"] not in {p.name for p in made})
