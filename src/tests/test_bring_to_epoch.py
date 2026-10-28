# tests/test_bring_to_epoch.py

"""
GEN.764: `galaxy/bring_to_epoch.py` brings a stored galaxy to the running epoch in part. It refuses a galaxy whose plan or
scatter is behind, generates again (in address order) the galaxy-placed sectors whose fill is behind, leaves a sector the
admin audit log names as edited, makes the species of an older population epoch again, and says each step. The epochs
are a made-up history here, so the tests do not move when a real epoch is added.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy import bring_to_epoch
from planetgen.population import model as population
from planetgen.queue import api_jobs
from planetgen.util import generator_epoch

from tests.test_bright_star_scatter import _seed_galaxy
from tests.test_migration_stage_epochs import _make_species, _system


def _history(monkeypatch, *affects):
    """Epochs 1, 2, 3, ... where epoch 1 changes nothing and the others change the stages `affects` names."""
    entries = [{"epoch": 1, "affects": []}] + [{"epoch": n, "affects": list(stages)} for n, stages in enumerate(affects, 2)]
    monkeypatch.setattr(generator_epoch, "epoch_history", lambda: tuple(entries))


class _Rows:
    def __init__(self, rows):
        self.rows = rows

    def fetchall(self):
        return self.rows


class _AuditLog:
    """A control connection that holds the given `(action, target)` audit rows."""

    def __init__(self, *rows):
        self.rows = [{"action": action, "target": target} for action, target in rows if target is not None]

    def execute(self, _sql, _args=()):
        return _Rows(self.rows)


def _sector(conn, name, address, epoch):
    shell, band, wedge = address if address else (None, None, None)
    return conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index, generator_epoch) "
                        "VALUES (?, 1000, ?, ?, ?, ?)", (name, shell, band, wedge, epoch)).lastrowid


@pytest.fixture
def galaxy(mysql_config):
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    with conn:      # the made-up histories start at epoch 1, so the plan is saved under it
        conn.execute("UPDATE galaxy_shape SET generator_epoch = 1 WHERE id = 1")
    yield conn
    conn.close()


def test_a_galaxy_at_the_running_epoch_has_nothing_to_do(galaxy, monkeypatch):
    _history(monkeypatch, ["fill"])
    _sector(galaxy, "Now", (2, 0, 1), 2)
    found = bring_to_epoch.survey(galaxy)
    assert found["refill"] == [] and found["refuse"] == []
    assert bring_to_epoch.report_lines(found) == ["The galaxy is at epoch 2: nothing to bring up."]


def test_a_plan_behind_is_refused_with_the_epochs_that_changed_it(galaxy, monkeypatch):
    _history(monkeypatch, ["fill"], ["plan", "fill"], ["fill"])
    galaxy.execute("UPDATE galaxy_shape SET generator_epoch = 2 WHERE id = 1")
    _sector(galaxy, "Old", (2, 0, 1), 2)
    found = bring_to_epoch.survey(galaxy)
    assert len(found["refuse"]) == 1 and "epoch 3 changed the galaxy's shape" in found["refuse"][0]
    assert "full regeneration" in found["refuse"][0]
    lines = bring_to_epoch.report_lines(found)
    assert lines[-1] == "Nothing was changed."
    assert not any("to be generated again" in line for line in lines)


def test_a_scatter_behind_is_refused_and_says_how_to_draw_it_again(galaxy, monkeypatch):
    _history(monkeypatch, ["scatter"], ["fill"])
    store.record_scatter_epoch(galaxy, ("phenomena",))
    galaxy.execute("UPDATE scatter_epochs SET generator_epoch = 1")
    found = bring_to_epoch.survey(galaxy)
    assert len(found["refuse"]) == 1
    assert "phenomena scatter" in found["refuse"][0] and "epoch 2 changed it" in found["refuse"][0]
    assert "planetgen clear-scatter" in found["refuse"][0]


def test_the_sectors_behind_are_listed_in_address_order_and_the_edited_one_is_left(galaxy, monkeypatch):
    _history(monkeypatch, ["fill"], ["population"], ["fill"])         # epochs 2 and 4 change the fill
    system, cfg = _system()
    with galaxy:
        late = _sector(galaxy, "Late", (3, 0, 1), 1)
        early = _sector(galaxy, "Early", (2, 0, 5), 2)
        edited = _sector(galaxy, "Edited", (2, 0, 7), 1)
        current = _sector(galaxy, "Current", (2, 0, 9), 4)
        _sector(galaxy, "By hand", None, 1)
        system_id = store.insert_star_system(galaxy, system, cfg)
        galaxy.execute("UPDATE star_systems SET sector_id = ? WHERE id = ?", (edited, system_id))
    audit = _AuditLog(("system.update", f"system:{system_id}"), ("sector.regenerate", f"sector:{current}"),
                      ("backup.take", None), ("system.update", "system:99999999"))
    found = bring_to_epoch.survey(galaxy, audit)
    assert [item["id"] for item in found["refill"]] == [early, late]
    assert [item["id"] for item in found["skipped"]] == [edited]
    assert found["skipped"][0]["why"] == "edited by an admin (system.update)"
    assert found["handmade"] == 1
    lines = bring_to_epoch.report_lines(found)
    assert lines[0] == "Bringing the galaxy to epoch 4; stages behind: fill, population."
    assert lines[-1] == "2 sector(s) to redo, 1 left alone."
    assert any("left alone: edited by an admin" in line for line in lines)


def test_a_run_generates_the_sectors_again_one_by_one_and_goes_on_after_a_failure(galaxy, mysql_config, monkeypatch):
    _history(monkeypatch, ["fill"])
    with galaxy:
        first = _sector(galaxy, "First", (2, 0, 1), 1)
        second = _sector(galaxy, "Second", (2, 0, 2), 1)
        third = _sector(galaxy, "Third", (2, 0, 3), 1)
    done = []

    def fake(sector_id, config):
        done.append(sector_id)
        if sector_id == second:
            raise RuntimeError("no room")
        return {"sector_id": sector_id, "sector_name": f"New {sector_id}"}

    monkeypatch.setattr(api_jobs, "regenerate_sector", fake)
    monkeypatch.setattr(population, "run_population", lambda conn: {"new_species": 0, "species": 0})
    monkeypatch.setattr(population, "run_territories", lambda conn: {"polities": 0, "owned_systems": 0})
    said = []
    result = bring_to_epoch.run(mysql_config, None, say=said.append)
    assert done == [first, second, third]
    assert result["redone"] == 2 and result["failed"] == [(result["failed"][0][0], "no room")]
    assert any("FAILED: no room" in line for line in said)
    assert said[-1].startswith("Done: 2 sector(s) generated again, 0 left alone, 1 failed")


def test_a_dry_run_and_a_refusal_change_nothing(galaxy, mysql_config, monkeypatch):
    _history(monkeypatch, ["fill"])
    with galaxy:
        _sector(galaxy, "Old", (2, 0, 1), 1)
    monkeypatch.setattr(api_jobs, "regenerate_sector", lambda *_: pytest.fail("a dry run generated a sector"))
    said = []
    result = bring_to_epoch.run(mysql_config, None, dry_run=True, say=said.append)
    assert result["redone"] == 0 and not result["refused"]
    assert any("to be generated again" in line for line in said)
    _history(monkeypatch, ["plan"], ["fill"])
    galaxy.execute("UPDATE galaxy_shape SET generator_epoch = 1 WHERE id = 1")
    galaxy.commit()
    result = bring_to_epoch.run(mysql_config, None, say=said.append)
    assert result["refused"] and result["redone"] == 0


def test_the_species_of_an_older_population_epoch_are_made_again_but_not_the_edited_ones(galaxy, mysql_config, monkeypatch):
    _history(monkeypatch, ["population"], ["population"])
    system, cfg = _system()
    with galaxy:
        sector = _sector(galaxy, "Sector", (2, 0, 1), 3)
        system_id = store.insert_star_system(galaxy, system, cfg)
        galaxy.execute("UPDATE star_systems SET sector_id = ? WHERE id = ?", (sector, system_id))
        _make_species(galaxy, system_id)
        galaxy.execute("UPDATE species SET generator_epoch = 1")
        galaxy.execute("INSERT INTO population_state (id, scanned_planet_id) VALUES (1, 500) "
                       "ON DUPLICATE KEY UPDATE scanned_planet_id = 500")
    homeworld = galaxy.execute("SELECT homeworld_planet_id AS id FROM species").fetchone()["id"]
    seen = []
    monkeypatch.setattr(population, "run_population", lambda conn: seen.append(1) or {"new_species": 4, "species": 4})
    monkeypatch.setattr(population, "run_territories", lambda conn: {"polities": 0, "owned_systems": 0})
    said = []
    result = bring_to_epoch.run(mysql_config, None, say=said.append)
    assert result["species"] == 1 and seen == [1]
    galaxy.commit()        # a fresh look at what the run committed
    assert galaxy.execute("SELECT COUNT(*) AS n FROM species").fetchone()["n"] == 0
    assert galaxy.execute("SELECT scanned_planet_id FROM population_state").fetchone()["scanned_planet_id"] == homeworld - 1

    with galaxy:
        _make_species(galaxy, system_id, name="Edited")
        galaxy.execute("UPDATE species SET generator_epoch = 1, tech_source = 'admin'")
    galaxy.commit()
    found = bring_to_epoch.survey(galaxy)
    assert found["species"] == {"redo": 0, "kept": 1}
    assert bring_to_epoch._redo_species(galaxy, [1]) == 0


@pytest.mark.db
@pytest.mark.usefixtures("unit_density_plans")
def test_a_sector_made_under_an_older_epoch_is_really_generated_again(mysql_config, redis_server, monkeypatch):
    """The whole path with nothing replaced: the sector is deleted and made again from its address, records the running
    epoch, holds the same systems (the code did not change), and a second run finds nothing to do."""
    from planetgen.generation import run_galaxy
    from tests.test_galaxy_gen import _all_sectors, _seed_skeleton, _shells
    _seed_skeleton(mysql_config, shells=_shells(2, 1))
    first = run_galaxy.ensure_sector_generated(0, 0, 0, config=mysql_config, backfill=False, settle=False)
    assert first["created"]
    conn = store.get_connection(mysql_config)
    try:
        before = conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?", (first["sector_id"],)).fetchone()["n"]
        _history(monkeypatch, ["fill"])
        with conn:
            conn.execute("UPDATE galaxy_shape SET generator_epoch = 1 WHERE id = 1")
            conn.execute("UPDATE sectors SET generator_epoch = 1")
    finally:
        conn.close()
    said = []
    result = bring_to_epoch.run(mysql_config, None, say=said.append)
    assert result["redone"] == 1 and not result["failed"], said
    [sector] = _all_sectors(mysql_config)
    assert sector["id"] != first["sector_id"]
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT generator_epoch FROM sectors WHERE id = ?", (sector["id"],)).fetchone()["generator_epoch"] == 2
        assert conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?", (sector["id"],)).fetchone()["n"] == before
    finally:
        conn.close()
    again = []
    assert bring_to_epoch.run(mysql_config, None, say=again.append)["redone"] == 0
    assert again == ["The galaxy is at epoch 2: nothing to bring up."]


def test_the_command_passes_the_options_and_exits_1_on_a_refusal(monkeypatch):
    from planetgen.cli import generate
    parser, _commands = generate.build_parser()
    args = parser.parse_args(["bring-to-epoch", "--dry-run"])
    assert args.command == "bring-to-epoch" and args.dry_run
    assert generate._COMMAND_HANDLERS["bring-to-epoch"] is generate.run_bring_to_epoch
    calls = []
    answer = {"refused": False, "failed": []}
    monkeypatch.setattr(bring_to_epoch, "run", lambda config, control, dry_run: calls.append(dry_run) or answer)
    generate.run_bring_to_epoch(args)
    assert calls == [True]
    answer["refused"] = True
    with pytest.raises(SystemExit) as stopped:
        generate.run_bring_to_epoch(args)
    assert stopped.value.code == 1
    assert generate.DATA_CHANGES["bring-to-epoch"] == ("generated", "population")
