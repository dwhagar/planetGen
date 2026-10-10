# tests/test_web_nearby.py

"""
The What's nearby page (`/nearby`, `web/nearby_page.py`, NAV.44) against a
small generated galaxy through the in-process API transport.
"""

import math
import re

import pytest

from planetgen import tuning
from planetgen.api.config import Config
from planetgen.db import store
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.geometry import sector_position_pc
from planetgen.generation import run_galaxy
from planetgen.physics.units import ly_to_pc
from planetgen.web.app import create_app

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
ADDRESSES = [(2, 0, 0), (2, 0, 1)]


@pytest.fixture
def client(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=8,
                            expected_system_count_at_density_1=2000.0, config=mysql_config)
    store.replace_galaxy_layers([(1, 6), (0, 8), (-1, 6)], config=mysql_config)
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
        return conn.execute("SELECT id, name, sector_id FROM star_systems ORDER BY id LIMIT 1").fetchone()
    finally:
        conn.close()


def test_the_page_walks_from_a_sector_to_a_system_to_a_distance(client, mysql_config):
    system = _first_system(mysql_config)
    html = client.get("/nearby").get_data(as_text=True)
    assert f'<option value="{system["sector_id"]}">' in html and "Or a point in space" in html
    html = client.get(f"/nearby?sector={system['sector_id']}").get_data(as_text=True)
    assert f'value="system:{system["id"]}"' in html
    html = client.get(f"/nearby?place=system:{system['id']}").get_data(as_text=True)
    assert 'name="distance"' in html and 'name="kinds"' in html


def test_the_list_shows_what_is_near_with_links_and_the_ungenerated_count(client, mysql_config):
    system = _first_system(mysql_config)
    response = client.get(f"/nearby?place=system:{system['id']}&distance=12")
    html = response.get_data(as_text=True)
    assert response.status_code == 200
    assert "Within 12 pc of" in html and system["name"] in html
    assert re.search(r'href="/system/\d+"', html)
    assert "are uncharted" in html or "is uncharted" in html
    assert 'id="nearby-results"' in html


def test_a_kind_filter_and_a_point_place(client, mysql_config):
    html = client.get("/nearby?place=0,0,0&distance=50&kinds=planet").get_data(as_text=True)
    assert "Within 50 pc of (0.00, 0.00, 0.00) pc" in html
    assert "<td>Star system</td>" not in html


def test_a_bad_distance_is_a_message_not_an_error_page(client, mysql_config):
    system = _first_system(mysql_config)
    response = client.get(f"/nearby?place=system:{system['id']}&distance=500")
    assert response.status_code == 200
    assert "the largest distance is 50 pc" in response.get_data(as_text=True)


def test_an_unknown_place_is_a_404(client, mysql_config):
    assert client.get("/nearby?place=nonsense&distance=5").status_code == 404


def test_system_sector_and_phenomenon_pages_link_here(client, mysql_config):
    system = _first_system(mysql_config)
    assert f"/nearby?place=system:{system['id']}" in client.get(f"/system/{system['id']}").get_data(as_text=True)
    assert f"/nearby?place=sector:{system['sector_id']}" in client.get(
        f"/sector/{system['sector_id']}").get_data(as_text=True)
    assert 'href="/nearby"' in client.get("/").get_data(as_text=True)
