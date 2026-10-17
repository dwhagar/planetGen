# tests/test_web_sector_profile.py

"""
UX.113: the sector page shows the profile of its cell, and `/sector/at/<shell>/<band>/<wedge>` shows it for any address,
with the heavy-cell warning before anything is run. The first tests fake the API like `test_web_sector_nav.py`; the
last uses a real throwaway database.
"""

import math

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.web.lib import apiclient
from planetgen.web.profile_view import FAR_RUN_SECTORS, HEAVY_CELL_SECONDS, profile_view

# Imported fixtures (see test_web_sector_nav.py).
from planetgen.web import csrf
from tests.test_web_sector_nav import _csrf, _log_in, app, client, db_app, db_client, fake  # noqa: F401


def _population(share, age=(0.0, 1.0), fe_h=((1.0, -0.1, 0.2),), dispersion=20.0, rotation="disk"):
    return {"density": share, "share": share, "age_gy": list(age), "fe_h": [list(p) for p in fe_h],
            "dispersion_kms": dispersion, "rotation": rotation}


def _profile(unresolved=None, seconds=3.0, superbubble=None, charted_sector_id=None):
    names = ("young", "intermediate", "old", "bulge", "halo", "cluster", "satellite", "nuclear_cluster", "nuclear_disc")
    shares = {"young": 0.2, "intermediate": 0.3, "old": 0.45, "halo": 0.05}
    populations = {name: _population(shares.get(name, 0.0)) for name in names}
    return {
        "shell_index": 40, "band_index": 0, "wedge_index": 3, "center_pc": [160.0, 0.0, 0.0], "edge_pc": 4.0,
        "volume_pc3": 64.0, "region": "disk", "density": 0.25, "expected_systems": 12.5,
        "mix": {name: p["density"] for name, p in populations.items()}, "populations": populations,
        "gas": {"cold_fraction": 0.6, "ambient": "disk", "ambient_density_cm3": 0.4, "mass_solar": 3.2,
                "superbubble": superbubble},
        "dark_matter": {"density_solar_pc3": 0.01, "mass_solar": 0.64},
        "kind_factors": {}, "type_ia_share": 0.25, "sector_id": charted_sector_id, "in_galaxy": True, "unresolved": unresolved,
        "estimate": {"systems": 12, "stars": 16, "bytes": 3_000_000, "seconds": seconds, "measured": True},
    }


def _uncharted(sector_id=None):
    return {"shell_index": 40, "band_index": 0, "wedge_index": 3, "designation": "ABC", "center_pc": [160.0, 0.0, 0.0],
            "edge_pc": 4.0, "sector_id": sector_id, "stars": [{"id": 1}], "scattered": []}


@pytest.fixture
def profiles(fake, monkeypatch):
    state = {"profile": _profile(), "cell": _uncharted(), "calls": []}

    def get_galaxy_profile(db, shell, band, wedge):
        state["calls"].append((shell, band, wedge))
        if state["profile"] is None:
            raise apiclient.NotFoundError("no galaxy has been planned")
        return state["profile"]

    monkeypatch.setattr(apiclient, "get_galaxy_profile", get_galaxy_profile)
    monkeypatch.setattr(apiclient, "get_uncharted_sector", lambda db, shell, band, wedge: state["cell"])
    return state


def test_a_charted_sector_page_shows_the_profile_of_its_address(client, profiles):
    html = client.get("/sector/5").get_data(as_text=True)
    assert 'id="sector-profile"' in html and "Stellar populations here" in html
    assert profiles["calls"] == [(5, 1, 20)]  # the fake sector's address
    assert "Old disk" in html and "45.0%" in html and "−0.10 ± 0.20" in html
    assert "Molecular gas is 60% of the solar circle" in html and "Dark matter:" in html
    assert "25% are Type Ia" in html
    assert 'id="heavy-cell"' not in html  # the sector exists already: nothing left to warn about


def test_a_sector_page_without_a_profile_still_renders(client, profiles):
    profiles["profile"] = None
    response = client.get("/sector/5")
    assert response.status_code == 200 and 'id="sector-profile"' not in response.get_data(as_text=True)


def test_an_uncharted_address_shows_its_profile_and_no_warning_for_a_light_cell(client, profiles):
    html = client.get("/sector/at/40/0/3").get_data(as_text=True)
    assert "Uncharted sector ABC" in html and "Nothing has been generated at shell 40, band 0, wedge 3" in html
    assert 'id="sector-profile"' in html and "about 12.5 star systems expected" in html
    assert "1 waiting bright star" in html and 'id="heavy-cell"' not in html
    assert "Generate page" not in html  # only an admin is pointed at it


def test_a_heavy_cell_warns_before_it_runs(client, profiles):
    unresolved = {"ceiling_systems": 5000.0, "luminosity_cap_sol": 30.0, "expected_systems": 2.0e6,
                  "unresolved_systems": 1.99e6, "unresolved_mass_solar": 1.0e6}
    profiles["profile"] = _profile(unresolved=unresolved, seconds=7300.0)
    html = client.get("/sector/at/40/0/3").get_data(as_text=True)
    assert 'id="heavy-cell"' in html and "more than the 5,000 a cell resolves" in html
    assert "Filling it takes about 2 h 01 m" in html and "3.0 MB" in html
    profiles["profile"] = _profile(seconds=HEAVY_CELL_SECONDS + 1)
    html = client.get("/sector/at/40/0/3").get_data(as_text=True)
    assert 'id="heavy-cell"' in html and "Filling it takes about 2 m 01 s" in html


def test_a_slow_estimate_from_an_unmeasured_server_says_so():
    view = profile_view(_profile(seconds=HEAVY_CELL_SECONDS) | {"estimate": {"systems": 1, "stars": 1, "bytes": 10,
                                                                           "seconds": HEAVY_CELL_SECONDS, "measured": False}})
    assert view["heavy"]["slow"] and view["heavy"]["measured"] is False
    assert profile_view(_profile(seconds=HEAVY_CELL_SECONDS), charted=True)["heavy"] is None


def test_a_charted_address_goes_to_its_sector_and_a_bad_one_is_a_404(client, profiles, monkeypatch):
    profiles["cell"] = _uncharted(sector_id=5)
    response = client.get("/sector/at/40/0/3")
    assert response.status_code == 302 and response.headers["Location"].endswith("/sector/5")

    def refused(db, shell, band, wedge):
        raise apiclient.ApiError("that cell is too far out", status_code=400)
    monkeypatch.setattr(apiclient, "get_uncharted_sector", refused)
    assert client.get("/sector/at/40/0/3").status_code == 404


def test_an_admin_is_pointed_at_the_generate_page(client, profiles, fake):
    fake.admin = {"username": "admin", "must_change_credentials": False}
    from planetgen.api.authz import SESSION_COOKIE_NAME
    client.set_cookie(SESSION_COOKIE_NAME, "token-value")
    assert "Generate page" in client.get("/sector/at/40/0/3").get_data(as_text=True)


def test_the_page_of_an_uncharted_address_on_a_real_galaxy(db_client, mysql_config):
    shape = build_galaxy_shape(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
                               bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=float(tuning.DEFAULT_SECTOR_EDGE_PC), outer_shell_index=1500,
                            expected_system_count_at_density_1=600.0, config=mysql_config, galaxy_seed=bytes(range(16)))
    response = db_client.get("/sector/at/40/0/3")
    html = response.get_data(as_text=True)
    assert response.status_code == 200
    assert "Uncharted sector" in html and 'id="sector-profile"' in html and "Stellar populations here" in html
    assert db_client.get("/sector/at/40/0/999999999").status_code == 404  # not a cell of that band
    assert "Traceback" not in html


def test_a_far_neighborhood_warns_before_it_runs(app, client, fake, profiles):
    """UX.113: the size-and-time dialog says so when the run holds more sectors than a run is warned about."""
    _log_in(client, fake)
    post = {"action": "generate_neighborhood", csrf.FIELD_NAME: _csrf(app, client)}
    assert 'id="far-run"' not in client.post("/sector/5", data=post).get_data(as_text=True)
    fake.estimate = dict(fake.estimate, sectors=FAR_RUN_SECTORS + 1)
    html = client.post("/sector/5", data=post).get_data(as_text=True)
    assert 'id="far-run"' in html and "A very large run" in html
