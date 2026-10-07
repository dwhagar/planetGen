# tests/test_neighborhood_edge.py

"""
GEN.65 A generation run fails from the web UI but not from the CLI
(docs/TODO.md). Boss thinks it happened while generating a neighborhood
centred close to the galaxy's edge, and the error text is lost, so this
generates neighborhoods through both paths -- the Generate page's job
(the form turned into `planetgen galaxy` arguments by `build_job` and run
by the job runner, as the page runs it) and the CLI in this process --
centred at the galaxy's rim, just inside and just outside it, in its top
and bottom layers, at the core and in a sparse region.

Needs a MySQL server (the `mysql_config` fixture); skips without one.
"""

import time

import pytest

from planetgen import tuning
from planetgen.galaxy.geometry import sector_position_pc
from planetgen.galaxy.skeleton import build_layer_extents
from planetgen.db import store
from planetgen.web import generate_page, jobs
from tests.test_galaxy_gen import EDGE_PC, _SKELETON_SHAPE, _mysql_argv, _run_cli, _seed_skeleton

E_VALUE = 0.1
"""float: The outline's threshold (a sector at the rim expects 1 / this
much density): about 13 rings out and 5 layers up and down."""

RADIUS_PC = EDGE_PC
"""float: One sector edge: the centre and its nearest neighbours."""

CORE_RADIUS_PC = 0.1
"""float: At the core a sector holds hundreds of systems, so there the
neighborhood is the centre sector alone."""

EXTENTS, OUTER_RING, _CONFIRMED = build_layer_extents(_SKELETON_SHAPE, EDGE_PC, 1.0 / E_VALUE)
TOP_LAYER, TOP_RING = EXTENTS[0]
BOTTOM_LAYER, BOTTOM_RING = EXTENTS[-1]
RIM = dict(EXTENTS)[0]

CENTRES = {
    "core": {"center_by": "position", "center_x_pc": "0", "center_y_pc": "0", "center_z_pc": "0"},
    "sparse": {"center_by": "address", "center_ring": str(RIM * 2 // 3), "center_layer": "0", "center_slot": "0"},
    "rim": {"center_by": "address", "center_ring": str(RIM), "center_layer": "0", "center_slot": "0"},
    "just inside": {"center_by": "address", "center_ring": str(RIM - 1), "center_layer": "0", "center_slot": "1"},
    "just outside": {"center_by": "address", "center_ring": str(RIM + 1), "center_layer": "0", "center_slot": "2"},
    "top layer": {"center_by": "address", "center_ring": str(TOP_RING), "center_layer": str(TOP_LAYER),
                  "center_slot": "0"},
    "bottom layer": {"center_by": "address", "center_ring": str(BOTTOM_RING), "center_layer": str(BOTTOM_LAYER),
                     "center_slot": "0"},
}


@pytest.fixture
def galaxy(mysql_config, tmp_path, monkeypatch):
    """A planned toy galaxy with a real outline, and a jobs directory."""
    monkeypatch.setenv("PLANETGEN_JOBS_DIR", str(tmp_path / "jobs"))
    _seed_skeleton(mysql_config, outer_ring_index=OUTER_RING, e_value=E_VALUE, layers=EXTENTS)
    return mysql_config


def _sectors(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return {(row["ring_index"], row["layer_index"], row["ring_slot_index"])
                for row in conn.execute("SELECT ring_index, layer_index, ring_slot_index FROM sectors").fetchall()}
    finally:
        conn.close()


def _inside(address):
    ring, layer, _slot = address
    extents = dict(EXTENTS)
    return layer in extents and ring <= extents[layer]


def _radius(form):
    return CORE_RADIUS_PC if form is CENTRES["core"] else RADIUS_PC


def _web_job(mysql_config, form, radius_pc=RADIUS_PC):
    """Runs the Generate page's job for `form` and returns it finished."""
    form = {"mode": "center", "center_radius_pc": f"{radius_pc:g}", **form}
    kind, title, steps = generate_page.build_job("galaxy", form, mysql_config.database, edge_pc=EDGE_PC)
    env = jobs.mysql_env(mysql_config, mysql_config.database)
    job_id = jobs.start_job(kind, title, steps, env=env, database=mysql_config.database)
    deadline = time.time() + 300
    while time.time() < deadline:
        job = jobs.get_job(job_id)
        if job and job["finished"]:
            return job, jobs.log_tail(job_id, max_bytes=16 * 1024)
        time.sleep(0.2)
    raise AssertionError(f"the job didn't finish: {jobs.log_tail(job_id)}")


def _cli_argv(form):
    argv, _description = generate_page.center_argv(
        {"center_radius_pc": f"{_radius(form):g}", **form}, EDGE_PC)
    return argv


def _named(form):
    """The address the form names (the position turned into its cell)."""
    argv = _cli_argv(form)
    return int(argv[argv.index("--ring") + 1]), int(argv[argv.index("--layer") + 1]), int(argv[argv.index("--slot") + 1])


def _check(mysql_config, form):
    sectors = _sectors(mysql_config)
    named = _named(form)
    # Only the named centre may lie outside the outline (GEN.81).
    assert {address for address in sectors if not _inside(address)} <= {named}
    if _inside(named):
        assert named in sectors


@pytest.mark.parametrize("where", sorted(CENTRES))
def test_a_neighborhood_from_the_generate_page(galaxy, where):
    job, output = _web_job(galaxy, CENTRES[where], _radius(CENTRES[where]))
    assert job["status"] == "succeeded", output
    assert "Traceback" not in output
    _check(galaxy, CENTRES[where])


@pytest.mark.parametrize("where", sorted(CENTRES))
def test_the_same_neighborhood_from_the_cli(galaxy, where):
    _run_cli(_cli_argv(CENTRES[where]) + _mysql_argv(galaxy))
    _check(galaxy, CENTRES[where])


@pytest.mark.parametrize("where", ["sparse", "rim", "top layer"])
def test_around_a_filled_sector_from_both_paths(galaxy, where):
    argv, _description = generate_page.center_argv({"center_radius_pc": "0.1", **CENTRES[where]}, EDGE_PC)
    _run_cli(argv + _mysql_argv(galaxy))
    conn = store.get_connection(galaxy)
    try:
        sector_id = conn.execute("SELECT id FROM sectors ORDER BY id LIMIT 1").fetchone()["id"]
    finally:
        conn.close()
    job, output = _web_job(galaxy, {"center_by": "sector", "center_sector": str(sector_id)})
    assert job["status"] == "succeeded", output
    _run_cli(["--center-sector", str(sector_id), "--radius-pc", f"{2 * RADIUS_PC:g}"]
             + _mysql_argv(galaxy))
    _check(galaxy, CENTRES[where])


def test_the_centres_span_the_outline():
    assert RIM > 6 and TOP_LAYER > 0 > BOTTOM_LAYER
    assert not _inside((RIM + 1, 0, 0)) and _inside((RIM, 0, 0))
    assert sector_position_pc(RIM, 0, 0, EDGE_PC)[0] > 0
    assert tuning.DEFAULT_GENERATE_RADIUS_PC > RADIUS_PC

