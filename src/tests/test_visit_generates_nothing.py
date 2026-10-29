# tests/test_visit_generates_nothing.py

"""PERF.149 rule 2: a visit to an uncharted sector generates nothing and queues nothing. The profile, the cell, the uncharted
scene and the locate answers are read-only on a cell no sector was generated in, no run starts, and no module of the site or
the API calls the functions that generate a sector. Making a sector is a run's work (`regenerate_sector_cell`, the Generate
page); an admin's visit queuing a one-sector run is PERF.204."""

import pathlib

import pytest

from planetgen.generation import run_galaxy
from planetgen.queue import redisqueue, runrecord

# Imported fixtures (see test_api_thin_routes.py).
from tests.test_api import client  # noqa: F401
from tests.test_api_galaxy_profile import _plan, _tables

ADDRESS = (1500 // 4, 3, 17)
"""tuple: A cell of the planned galaxy no sector was generated in."""

VISITS = (
    "/api/galaxy/profile?shell=%d&band=%d&wedge=%d",
    "/api/galaxy/cell?shell=%d&band=%d&wedge=%d",
    "/api/galaxy/uncharted?shell=%d&band=%d&wedge=%d",
)


@pytest.mark.parametrize("visit", VISITS)
def test_a_visit_to_an_uncharted_sector_changes_nothing_and_starts_no_run(client, mysql_config, visit):
    _plan(mysql_config)
    before = _tables(mysql_config)
    response = client.get(visit % ADDRESS)
    assert response.status_code == 200, response.get_data(as_text=True)
    assert response.get_json().get("sector_id") is None
    assert _tables(mysql_config) == before
    conn = redisqueue.connect()
    assert runrecord.recent_runs(conn, limit=10) == []
    assert redisqueue.queue("default", conn).count == 0


def test_the_sector_entry_point_that_generated_on_a_visit_is_gone():
    assert not hasattr(run_galaxy, "ensure_sector_generated")


def test_no_module_of_the_site_or_the_api_generates_a_sector_in_a_request():
    """The functions that save a sector are called by runs and API jobs (`generation/`, `queue/`), never by a page or a
    route directly."""
    root = pathlib.Path(__file__).resolve().parents[1] / "planetgen"
    generating = ("generate_and_save_sector_at", "regenerate_sector_cell", "ensure_sector_generated", "_cell_to_generate")
    found = {f"{folder}/{path.name}": [word for word in generating if word in path.read_text()]
             for folder in ("web", "api") for path in sorted((root / folder).rglob("*.py"))}
    assert {name: words for name, words in found.items() if words} == {}
