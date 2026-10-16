# tests/test_orbits_during_update.py

"""
OPS.54: the orbital update (`planetgen.cli.orbits`, which the nightly timer, a Generate job's step and a queued job
all run) moves nothing while an update runs, and a cancel never leaves the stored positions half moved. Real
database, the whole command.
"""

import os
import signal
import subprocess
import sys

import pytest

from planetgen.admin import maintenance, updateguard
from planetgen.cli import orbits
from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.galaxy.sector import SpaceSector
from planetgen.physics.units import pc_to_ly
from tests.test_galactic_motion import EDGE_PC, _place
from tests.test_update_guard import begin


def argv(config):
    return ["planetgen.cli.orbits", "--mysql-host", config.host, "--mysql-port", str(config.port),
            "--mysql-user", config.user, "--mysql-password", config.password or "",
            "--mysql-database", config.database]


def run_orbits(config, monkeypatch):
    monkeypatch.setattr(sys, "argv", argv(config))
    try:
        orbits.main()
    except SystemExit as exc:
        return exc.code or 0
    return 0


@pytest.fixture
def galaxy(mysql_config, monkeypatch):
    """A database with one system whose clocks are back-dated, so a run has something due."""
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    sector = SpaceSector("Clockwork", edge_ly=pc_to_ly(EDGE_PC))
    sector.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    _sector, position = _place(5, [], "unused")
    store.save_sector(sector, config=mysql_config, galaxy_position=position)
    assert run_orbits(mysql_config, monkeypatch) == 0           # the first run only sets the clock
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE orbit_simulation_state SET last_updated_at = NOW() - INTERVAL 400 DAY")
            conn.execute("UPDATE star_systems SET epoch_unix = NULL, next_update_due = NULL")
    finally:
        conn.close()
    return mysql_config


def snapshot(config):
    """Everything an orbital update changes: the clock, and each system's clock and place."""
    conn = store.get_connection(config)
    try:
        clock = conn.execute("SELECT last_updated_at FROM orbit_simulation_state").fetchall()
        systems = conn.execute("SELECT id, sector_id, epoch_unix, next_update_due, galaxy_x_pc, galaxy_y_pc,"
                               " galaxy_z_pc FROM star_systems ORDER BY id").fetchall()
        stars = conn.execute("SELECT id, galactic_orbital_phase_deg FROM stars ORDER BY id").fetchall()
        return clock, systems, stars
    finally:
        conn.close()


def test_while_an_update_runs_the_orbital_update_moves_nothing_and_exits_75(galaxy, monkeypatch, capsys):
    before = snapshot(galaxy)
    assert begin(sources=[])
    capsys.readouterr()
    assert run_orbits(galaxy, monkeypatch) == updateguard.EX_TEMPFAIL == 75
    captured = capsys.readouterr()
    assert "did not run and moved nothing" in captured.err and "the site is being updated" in captured.err
    assert "Moved" not in captured.out
    assert snapshot(galaxy) == before           # no row, no epoch, no due time and no clock changed


def test_once_the_update_has_finished_the_same_run_moves_the_objects(galaxy, monkeypatch, capsys):
    before = snapshot(galaxy)
    assert begin(sources=[])
    assert run_orbits(galaxy, monkeypatch) == 75
    assert snapshot(galaxy) == before
    updateguard.end(os.getpid(), ok=True)
    capsys.readouterr()
    assert run_orbits(galaxy, monkeypatch) == 0
    assert "Moved 1 object(s) along their galactic orbits" in capsys.readouterr().out
    assert snapshot(galaxy) != before


def test_a_site_a_failed_update_left_closed_refuses_it_too_until_an_admin_opens_it(galaxy, monkeypatch):
    before = snapshot(galaxy)
    assert begin(sources=[])
    updateguard.end(os.getpid(), ok=False, step="step 5/9", write=lambda text: None)
    assert run_orbits(galaxy, monkeypatch) == 75 and snapshot(galaxy) == before
    maintenance.turn_off(by="boss")
    assert run_orbits(galaxy, monkeypatch) == 0 and snapshot(galaxy) != before


def test_an_admins_own_maintenance_does_not_stop_the_orbital_update(galaxy, monkeypatch):
    maintenance.turn_on("an hour", by="boss")
    assert run_orbits(galaxy, monkeypatch) == 0


def test_the_command_a_queued_job_or_timer_runs_exits_75_without_touching_the_database():
    assert begin(sources=[])
    result = subprocess.run([sys.executable, "-m", "planetgen.cli.orbits", "--mysql-host", "no-such-host.invalid"],
                            capture_output=True, text=True, env=os.environ.copy(), timeout=120)
    assert result.returncode == 75 and "did not run and moved nothing" in result.stderr
    assert "error: could not open the database" not in result.stderr       # it never tried


def test_a_cancel_during_the_run_finishes_the_update_so_nothing_is_half_moved(galaxy, monkeypatch, capsys):
    before = snapshot(galaxy)
    real = orbits.advance_comet_orbits
    seen = {}

    def cancelled_here(conn, clock, **kwargs):
        seen["registered"] = [item["kind"] for item in updateguard.MarkerSource().items()]
        os.kill(os.getpid(), signal.SIGTERM)        # what the update's cancel sends, between two stages
        return real(conn, clock, **kwargs)

    monkeypatch.setattr(orbits, "advance_comet_orbits", cancelled_here)
    capsys.readouterr()
    assert run_orbits(galaxy, monkeypatch) == 0
    out = capsys.readouterr().out
    assert "Cancel requested: finishing this orbital update first" in out
    assert "Moved 1 object(s) along their galactic orbits" in out and "the update finished first" in out
    after = snapshot(galaxy)
    assert after[0] != before[0] and after[1] != before[1]      # the clock and the system moved together
    assert seen["registered"] == ["orbits"]                      # the guard could see the run while it ran
    assert updateguard.MarkerSource().items() == []              # and forgets it afterwards


def test_a_cancel_before_anything_moved_stops_at_once_and_changes_nothing(galaxy, monkeypatch, capsys):
    before = snapshot(galaxy)
    real = orbits.get_orbit_update_elapsed_years

    def cancelled_early(conn):
        os.kill(os.getpid(), signal.SIGTERM)
        return real(conn)

    monkeypatch.setattr(orbits, "get_orbit_update_elapsed_years", cancelled_early)
    assert run_orbits(galaxy, monkeypatch) == 130
    assert "Cancelled before any position moved" in capsys.readouterr().err
    assert snapshot(galaxy) == before


def test_the_signal_handlers_are_put_back(galaxy, monkeypatch):
    term, interrupt = signal.getsignal(signal.SIGTERM), signal.getsignal(signal.SIGINT)
    assert run_orbits(galaxy, monkeypatch) == 0
    assert signal.getsignal(signal.SIGTERM) is term and signal.getsignal(signal.SIGINT) is interrupt
