# tests/test_web_db_fields.py

"""
The Database thread's display fields (schemas v36 and v37) on the web
pages: a black hole's Class and a rogue planet's Mass Class on the
phenomenon page, a runaway or hypervelocity star on system rows and the
system page, the sector's star and interstellar-debris figures, and many
rogue planets folded into one Contents row, plus the nebula or remnant
a system or phenomenon sits inside (schema v39).
"""

from api.app import create_app
from api.config import Config

import web  # noqa: F401 -- puts src/html/lib on sys.path
import fmt  # noqa: E402
from web import sector_page, system_pages  # noqa: E402


class _Config(Config):
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


def _fields(kind, detail):
    return dict(system_pages.phenomenon_fields(kind, detail))


def test_black_hole_class_and_rogue_mass_class_rows():
    assert _fields("black_hole", {"mass_class": "supermassive"})["Class"] == "Supermassive"
    assert _fields("rogue_planet", {"mass_bin": "brown-dwarf"})["Mass Class"] == "Brown dwarf"
    assert _fields("rogue_planet", {"mass_bin": "sub-neptune"})["Mass Class"] == "Sub-Neptune"


def test_runaway_text():
    assert fmt.runaway_text({"runaway_class": None}) is None
    assert fmt.runaway_text({"runaway_class": "runaway", "runaway_speed_kms": 84.4}) == \
        "Runaway star, 84 km/s"
    assert fmt.runaway_text({"runaway_class": "hypervelocity", "runaway_speed_kms": 1234.0}) == \
        "Hypervelocity star, 1,234 km/s"


def test_debris_html():
    assert sector_page.debris_html(0) is None
    assert str(sector_page.debris_html(7.2e12)) == \
        "About 7&times;10<sup>12</sup> interstellar comets and planetesimals (estimated)"
    assert str(sector_page.debris_html(9.7e12)).startswith("About 1&times;10<sup>13</sup>")
    assert str(sector_page.debris_html(250)).startswith("About 250 ")


def _rogue(i, distance):
    return {"type": "rogue_planet", "id": i, "name": f"Rogue {i}", "descriptor": "gas giant",
            "radius_ly": 0, "distance_ly": distance}


def test_many_rogue_planets_fold_into_one_contents_row():
    sector = {
        "systems": [{
            "id": 1, "name": "Sol", "quadrant": "+X+Y+Z", "location": "", "is_binary": 0, "binary_type": None,
            "position_x_mpc": None, "position_y_mpc": None, "position_z_mpc": None, "center_distance_ly": 1.0,
            "runaway_class": "runaway", "runaway_speed_kms": 50.0,
            "stars": [{"star_type": "G2V", "temperature_k": 5800, "radius_km": 7e5, "luminosity_w": 3.8e26}],
        }],
        "phenomena": [_rogue(3, 4.0), _rogue(2, 2.0), {**_rogue(9, 3.0), "type": "black_hole",
                                                         "name": "Hole", "descriptor": "quiescent"}],
    }
    with create_app(_Config).test_request_context("/sector/1"):
        rows, _map = sector_page._contents(sector)
        group = [row for row in rows if row.get("members")]
        assert len(group) == 1 and group[0]["name"] == "2 rogue planets"
        assert [m["name"] for m in group[0]["members"]] == ["Rogue 2", "Rogue 3"]
        assert group[0]["distance_ly"] == 2.0
        assert rows[0]["details"] == "G2V, Runaway star, 50 km/s"
        assert len(rows) == 3

        sector["phenomena"] = [_rogue(3, 4.0)]
        rows, _map = sector_page._contents(sector)
        assert not any(row.get("members") for row in rows)
        assert any(row["name"] == "Rogue 3" for row in rows)


def test_phenomenon_class_and_contents_rows():
    nebula = _fields("nebula", {"nebula_class": "D", "nebula_type": "emission", "dominant_species": "H II",
                                "density_cm3": 120.0, "temperature_k": 8000.0, "extinction_av": 0.0})
    assert nebula["Class"] == "D: Classical H II region"
    assert nebula["Dominant Species"] == "H II"
    assert nebula["Density"] == "120 particles/cm³"
    assert nebula["Gas Temperature"] == "8,000 K"
    assert "Extinction" not in nebula  # stored 0: unset
    remnant = _fields("supernova_remnant", {"remnant_class": "T"})
    assert remnant["Class"] == "T: Pulsar wind nebula (plerion)"
    field = _fields("asteroid_field", {"field_class": "C3", "composition_family": "carbonaceous"})
    assert field["Class"] == "C3" and field["Composition Family"] == "Carbonaceous"


def test_diffuse_nebulae_have_a_sector_map_color():
    import starmap

    assert "diffuse" in starmap._NEBULA_TYPE_COLORS and "diffuse" in starmap._NEBULA_TYPE_ALPHA


def test_inside_text_and_contents_details():
    assert fmt.inside_text({"inside": None}) is None
    assert fmt.inside_text({"inside": {"type": "nebula", "id": 4, "name": "Veil", "class": "E"}}) == "Inside Veil"
    sector = {
        "systems": [{
            "id": 1, "name": "Sol", "quadrant": "+X+Y+Z", "location": "", "is_binary": 0, "binary_type": None,
            "position_x_mpc": None, "position_y_mpc": None, "position_z_mpc": None, "center_distance_ly": 1.0,
            "runaway_class": None, "runaway_speed_kms": None,
            "inside": {"type": "nebula", "id": 4, "name": "Veil", "class": "E"},
            "stars": [{"star_type": "G2V", "temperature_k": 5800, "radius_km": 7e5, "luminosity_w": 3.8e26}],
        }],
        "phenomena": [],
    }
    with create_app(_Config).test_request_context("/sector/1"):
        rows, _map = sector_page._contents(sector)
    assert rows[0]["details"] == "G2V, Inside Veil"


def test_phenomenon_inside_links_its_cloud(monkeypatch):
    calls = []

    def fake_get(db, kind, cloud_id):
        calls.append((kind, cloud_id))
        return {"id": cloud_id, "name": "Crab"}

    monkeypatch.setattr(system_pages.apiclient, "get_phenomenon", fake_get)
    with create_app(_Config).test_request_context("/phenomenon/rogue_planet/1"):
        assert system_pages._phenomenon_inside({"inside_nebula_id": None, "inside_remnant_id": None}) is None
        inside = system_pages._phenomenon_inside({"inside_nebula_id": None, "inside_remnant_id": 7})
    assert calls == [("supernova_remnant", 7)]
    assert inside["name"] == "Crab" and inside["url"].endswith("/phenomenon/supernova_remnant/7")


def test_system_detail_names_the_containing_nebula(mysql_config):
    import queryDb
    from stellarObjects import _db
    from tests.test_db_persistence import _placed_nebula, _sector_with_one_system

    sector_id = _db.save_sector(_sector_with_one_system("Misty"), config=mysql_config, galaxy_position={
        "center_x_pc": 100.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 100.0,
    })
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            nebula_id = _placed_nebula(conn, sector_id, (101.0, 0.0, 0.0), radius_ly=100.0, nebula_class="E")
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
        inside = queryDb.system_detail(conn, system_id)["inside"]
        assert (inside["type"], inside["id"]) == ("nebula", nebula_id)
    finally:
        conn.close()
