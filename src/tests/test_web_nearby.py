# tests/test_web_nearby.py

"""
The What's nearby page (`/nearby`, `web/nearby_page.py`, NAV.44) against a
small generated galaxy through the in-process API transport.
"""

import math
import re

import pytest
from markupsafe import escape

from planetgen import tuning
from planetgen.galaxy import spherical
from planetgen.api.config import Config
from planetgen.db import store
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.generation import run_galaxy, run_plan
from planetgen.physics.units import ly_to_pc
from tests import first_fill
from planetgen.web.app import create_app
from tests.toy_shape import toy_shape
from tests.publicids import pid

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = toy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
ADDRESSES = [(2, 0, 0), (2, 0, 1)]


@pytest.fixture
def client(mysql_config, monkeypatch):
    first_fill.without_the_deepest_floor(monkeypatch)   # the four systems asked for, not the deep floor's thousands per sector (GEN.633; TEST.475)
    monkeypatch.setattr(run_plan, "_draw_sector_bands", lambda *_args, **_kwargs: (0, 0))   # nor the bright-star bands of a capped cell (minutes of band tables; TEST.475)
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8,
                            expected_system_count_at_density_1=2000.0, config=mysql_config)
    store.replace_galaxy_shells([(shell, min(1, spherical.band_range(shell)[1])) for shell in range(9)],
                                config=mysql_config)
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 4
    for address in ADDRESSES:
        run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)

    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    app = create_app(RealConfig)
    app.testing = True
    return app.test_client()


def _first_system(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        row = dict(conn.execute("SELECT id, name, sector_id FROM star_systems ORDER BY id LIMIT 1").fetchone())
        row["pid"] = pid("system", row["id"], conn)
        row["sector_pid"] = pid("sector", row["sector_id"], conn)
        return row
    finally:
        conn.close()


def test_the_nav_section_is_current_on_both_of_its_pages(client):
    """UX.146: the Nav button in the bar is current on Plot Course and on Nearby Search, and the menu marks which."""
    for path, entry in (("/nav", "Plot Course"), ("/nearby", "Nearby Search")):
        html = client.get(path).get_data(as_text=True)
        bar = re.search(r'<nav class="site-sections" aria-label="Main">.*?</nav>', html, re.S).group(0)
        assert re.search(r'class="site-nav-trigger"[^>]*aria-current="page"|aria-current="page"[^>]*class="site-nav-trigger"',
                         bar), path
        assert re.findall(r'<sl-menu-item data-href="[^"]*" aria-current="page">([^<]*)<', bar) == [entry], path


def test_the_page_walks_from_a_sector_to_a_system_to_a_distance(client, mysql_config):
    system = _first_system(mysql_config)
    html = client.get("/nearby").get_data(as_text=True)
    assert f'<sl-option value="{system["sector_pid"]}">' in html and "Or a point in space" in html
    html = client.get(f"/nearby?sector={system['sector_pid']}").get_data(as_text=True)
    assert f'value="system:{system["pid"]}"' in html
    html = client.get(f"/nearby?place=system:{system['pid']}").get_data(as_text=True)
    assert 'name="distance"' in html and 'name="kinds"' in html


def _add_unbuilt_bright_star(mysql_config, system_id):
    """One pre-placed bright star, not built, 1 pc from the system."""
    conn = store.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT galaxy_x_pc x, galaxy_y_pc y, galaxy_z_pc z FROM star_systems WHERE id = ?",
                           (system_id,)).fetchone()
        point = (row["x"] + 1.0, row["y"], row["z"])
        store.insert_bright_stars(conn, [(*spherical.sector_address_at(point, EDGE_PC),
                                          *(int(round(v * 1000)) for v in point), "young", "B", "III", 1e31, 1e7,
                                          15000.0, 1e28, 0.05, 0.1, 8.0, None, 1, None, None)])   # the last value: the rotation period (GEN.846), which a scattered star stores
        conn.commit()
    finally:
        conn.close()


def test_the_list_shows_what_is_near_with_links_and_the_ungenerated_count(client, mysql_config):
    system = _first_system(mysql_config)
    _add_unbuilt_bright_star(mysql_config, system["id"])
    response = client.get(f"/nearby?place=system:{system['pid']}&distance=12")
    html = response.get_data(as_text=True)
    assert response.status_code == 200
    assert "Within 12 pc of" in html and escape(system["name"]) in html
    assert re.search(r'href="/system/[0-9A-F-]+"', html)
    assert "1 object in range is not generated yet (bright stars: 1)" in html
    assert 'id="nearby-results"' in html


def test_a_crowded_sphere_shows_its_early_stop_count_and_a_huge_one_is_a_message(client, mysql_config, monkeypatch):
    from planetgen.db import near
    system = _first_system(mysql_config)
    url = f"/nearby?place=system:{system['pid']}&distance=50"
    assert "star systems in range" not in client.get(url).get_data(as_text=True)
    monkeypatch.setattr(near, "COUNT_CAP", 1)
    html = client.get(url.replace("distance=50", "distance=49")).get_data(as_text=True)  # another URL: the page may be served from the cache
    assert "1+ star systems in range" in html
    monkeypatch.setattr(near, "SPHERE_LIMIT", 1)
    response = client.get(url.replace("distance=50", "distance=48"))
    assert response.status_code == 200
    assert "holds more than 1 systems" in response.get_data(as_text=True)
    assert "shrink the distance" in response.get_data(as_text=True)


def test_a_kind_filter_and_a_point_place(client, mysql_config):
    html = client.get("/nearby?place=0,0,0&distance=50&kinds=planet").get_data(as_text=True)
    assert "Within 50 pc of (0.00, 0.00, 0.00) pc" in html
    assert "<td>Star system</td>" not in html


def test_a_bad_distance_is_a_message_not_an_error_page(client, mysql_config):
    system = _first_system(mysql_config)
    response = client.get(f"/nearby?place=system:{system['pid']}&distance=500")
    assert response.status_code == 200
    assert "the largest distance is 50 pc" in response.get_data(as_text=True)


def test_an_unknown_place_is_a_404(client, mysql_config):
    assert client.get("/nearby?place=nonsense&distance=5").status_code == 404


def test_system_sector_and_phenomenon_pages_link_here(client, mysql_config):
    system = _first_system(mysql_config)
    assert f"/nearby?place=system:{system['pid']}" in client.get(f"/system/{pid('system', system['id'])}").get_data(as_text=True)
    assert f"/nearby?place=sector:{system['sector_pid']}" in client.get(
        f"/sector/{pid('sector', system['sector_id'])}").get_data(as_text=True)
    assert 'href="/nearby"' in client.get("/").get_data(as_text=True)
