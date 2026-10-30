# tests/test_bring_to_epoch.py

"""
GEN.764: `galaxy/bring_to_epoch.py` brings a stored galaxy to the running epoch in part. It refuses a galaxy whose plan
is behind, draws a scatter that is behind again (GEN.845), generates again (in address order) the galaxy-placed sectors whose fill is behind, leaves a sector the
admin audit log names as edited, makes the species of an older population epoch again, and says each step. The epochs
are a made-up history here, so the tests do not move when a real epoch is added.
"""

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import bring_to_epoch
from planetgen.population import model as population
from planetgen.queue import api_jobs
from planetgen.util import generator_epoch

from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy, light_steps  # noqa: F401 -- a fixture
from tests.test_migration_stage_epochs import _make_species, _system
from tests.sector_support import make_sector


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


def test_a_scatter_behind_is_drawn_again_with_the_settings_it_was_drawn_with(galaxy, monkeypatch):
    _history(monkeypatch, ["scatter"], ["fill"])
    store.record_scatter_epoch(galaxy, ("phenomena", "nebula-field"))
    galaxy.execute("UPDATE scatter_epochs SET generator_epoch = 1")
    galaxy.execute("UPDATE galaxy_shape SET phenomenon_scatter_seed = 7, phenomenon_min_mass_solar = 4.0 WHERE id = 1")
    found = bring_to_epoch.survey(galaxy)
    assert found["refuse"] == [] and found["redraw"] == ["phenomena", "nebula-field"]
    assert found["settings"]["compact"] == 4.0
    lines = bring_to_epoch.report_lines(found)
    assert "phenomena scatter: made under epoch 1, to be cleared and drawn again." in lines
    assert not bring_to_epoch._nothing_to_do(found)


def test_the_scatters_are_drawn_in_the_order_they_depend_on_each_other():
    none = {"mass_limit": 8.0}
    assert bring_to_epoch.redraw_order({"phenomena"}, none) == ["phenomena"]
    assert bring_to_epoch.redraw_order({"nebula-field", "phenomena"}, none) == ["phenomena", "nebula-field"]
    assert bring_to_epoch.redraw_order({"luminosity"}, none) == ["luminosity"]
    assert bring_to_epoch.redraw_order({"band"}, none) == ["luminosity"], "the band is part of the bright stars"
    assert bring_to_epoch.redraw_order({"mass"}, none) == ["mass", "luminosity"], "the bright stars are split from the massive ones"
    assert bring_to_epoch.redraw_order({"phenomena", "mass", "band"}, none) == ["mass", "luminosity", "phenomena"]
    assert bring_to_epoch.redraw_order({"luminosity"}, {"mass_limit": None}) == ["mass", "luminosity"], \
        "with no mass limit stored the bright stars cannot be drawn alone"


def test_a_scatter_command_gets_the_stored_settings_and_a_setting_it_refuses_stops_the_run(galaxy, mysql_config, monkeypatch):
    settings = {"floor": 5000.0, "mass_limit": 8.0, "compact": 2.0}
    mass = bring_to_epoch.scatter_args(mysql_config, "mass", settings)
    assert mass.mass and mass.bright_star_min_luminosity == 5000.0 and mass.phenomenon_min_mass == 8.0
    assert bring_to_epoch.scatter_args(mysql_config, "luminosity", settings).bright_star_min_luminosity == 5000.0
    assert bring_to_epoch.scatter_args(mysql_config, "phenomena", settings).compact_min_mass == 2.0
    assert bring_to_epoch.scatter_args(mysql_config, "phenomena", {**settings, "compact": 8.0}).compact_min_mass is None
    assert bring_to_epoch.scatter_args(mysql_config, "luminosity", {"floor": None, "mass_limit": None, "compact": None}).luminosity
    assert bring_to_epoch.scatter_args(mysql_config, "mass", None, clear=True).mass
    with pytest.raises(ValueError):
        bring_to_epoch.scatter_args(mysql_config, "mass", {**settings, "mass_limit": 7.5})
    _history(monkeypatch, ["scatter"], ["fill"])
    store.record_scatter_epoch(galaxy, ("mass", "luminosity"))
    galaxy.execute("UPDATE scatter_epochs SET generator_epoch = 1")
    galaxy.execute("UPDATE galaxy_shape SET bright_star_min_luminosity_sol = 5000, bright_star_seed = 1,"
                   " bright_star_mass_limit_sol = 7.5 WHERE id = 1")
    found = bring_to_epoch.survey(galaxy)
    assert len(found["refuse"]) == 1 and "mass scatter cannot be drawn again" in found["refuse"][0]
    assert bring_to_epoch.report_lines(found)[-1] == "Nothing was changed."


def test_a_run_empties_the_sectors_then_clears_and_draws_the_scatters_then_generates_the_sectors(galaxy, mysql_config, monkeypatch):
    from planetgen.generation import run_plan
    _history(monkeypatch, ["scatter", "fill"])
    with galaxy:
        first = _sector(galaxy, "First", (2, 0, 1), 1)
        second = _sector(galaxy, "Second", (2, 0, 2), 1)
    store.record_scatter_epoch(galaxy, ("mass", "luminosity", "phenomena"))
    galaxy.execute("UPDATE scatter_epochs SET generator_epoch = 1")
    galaxy.execute("UPDATE galaxy_shape SET bright_star_min_luminosity_sol = 5000, bright_star_seed = 1, bright_star_mass_limit_sol = 8,"
                   " phenomenon_scatter_seed = 7, phenomenon_min_mass_solar = 8 WHERE id = 1")
    galaxy.commit()
    steps = []
    stored = {"mass": True, "luminosity": True, "phenomena": True}
    monkeypatch.setattr(api_jobs, "empty_sector", lambda sector_id, config: steps.append(("empty", sector_id)) or
                        {"address": (2, 0, sector_id), "around": [], "deleted": {}})
    monkeypatch.setattr(api_jobs, "fill_emptied_sector", lambda emptied, config: steps.append(("fill", emptied["address"][2])) or
                        {"sector_id": emptied["address"][2], "sector_name": "New"})
    monkeypatch.setattr(api_jobs, "regenerate_sector", lambda *_: pytest.fail("a sector was deleted and made in one step"))
    monkeypatch.setattr(store, "scatter_stored", lambda conn, kind: stored[kind])
    monkeypatch.setattr(run_plan, "run_clear_scatter", lambda args: steps.append(("clear", run_plan.stages.clear_kind(args))))
    monkeypatch.setattr(run_plan, "run_scatter", lambda args: steps.append(
        ("draw", run_plan.stages.scatter_kind(args), args.bright_star_min_luminosity, args.phenomenon_min_mass, args.compact_min_mass)))
    monkeypatch.setattr(population, "run_population", lambda conn: {"new_species": 0, "species": 0})
    monkeypatch.setattr(population, "run_territories", lambda conn: {"polities": 0, "owned_systems": 0})
    said = []
    result = bring_to_epoch.run(mysql_config, None, say=said.append)
    assert not result["failed"] and result["redone"] == 2, said
    assert steps == [("empty", first), ("empty", second),
                     ("clear", "luminosity"), ("clear", "mass"),
                     ("draw", "mass", 5000.0, 8.0, None), ("draw", "luminosity", 5000.0, None, None),
                     ("clear", "phenomena"), ("draw", "phenomena", tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL, None, None),
                     ("fill", first), ("fill", second)]


def test_a_scatter_that_fails_leaves_the_emptied_sectors_and_says_what_to_do_by_hand(galaxy, mysql_config, monkeypatch):
    from planetgen.generation import run_plan
    _history(monkeypatch, ["scatter", "fill"])
    with galaxy:
        _sector(galaxy, "First", (2, 0, 1), 1)
    store.record_scatter_epoch(galaxy, ("nebula-field",))
    galaxy.execute("UPDATE scatter_epochs SET generator_epoch = 1")
    galaxy.commit()
    monkeypatch.setattr(api_jobs, "empty_sector", lambda sector_id, config: {"address": (2, 0, 1), "around": [], "deleted": {}})
    monkeypatch.setattr(api_jobs, "fill_emptied_sector", lambda *_: pytest.fail("a sector was generated from a broken scatter"))
    monkeypatch.setattr(store, "scatter_stored", lambda conn, kind: False)

    def broken(args):
        raise RuntimeError("no room on the disk")

    monkeypatch.setattr(run_plan, "run_scatter", broken)
    said = []
    result = bring_to_epoch.run(mysql_config, None, say=said.append)
    assert result["failed"] == [("scatter", "no room on the disk")] and result["redone"] == 0
    assert any("planetgen scatter --KIND" in line and "planetgen sectors" in line for line in said)


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
    short = bring_to_epoch.report_lines(found, limit=1)         # the Generate page names a few and counts the rest
    assert any(line == "... and 1 more sector(s) to be generated again." for line in short)
    assert not any("more sector(s) left alone" in line for line in short)


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
    first = make_sector(0, 0, 0, config=mysql_config, backfill=False, settle=False)
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


def test_the_status_for_the_generate_page_changes_nothing_and_says_what_a_run_would_do(galaxy, mysql_config, monkeypatch):
    _history(monkeypatch, ["fill"])
    with galaxy:
        sector = _sector(galaxy, "Behind", (2, 0, 5), 1)
    seen = []
    monkeypatch.setattr(api_jobs, "regenerate_sector", lambda sector_id, config: seen.append(sector_id))
    found = bring_to_epoch.status(mysql_config)
    assert found["error"] is None and not found["nothing"] and found["current"] == 2
    assert [item["id"] for item in found["refill"]] == [sector]
    assert found["lines"][-1] == "1 sector(s) to redo, 0 left alone."
    assert seen == []
    with galaxy:
        galaxy.execute("UPDATE sectors SET generator_epoch = 2")
    again = bring_to_epoch.status(mysql_config)
    assert again["nothing"] and again["lines"] == ["The galaxy is at epoch 2: nothing to bring up."]


def test_every_scatter_behind_is_really_drawn_again_with_the_same_settings(galaxy, mysql_config, monkeypatch):
    """The whole path with nothing replaced: the scatters drawn under epoch 1 are cleared and drawn again under epoch 2 from the
    settings they were drawn with, and (the code did not change) hold the same rows."""
    from planetgen.generation import run_plan
    _history(monkeypatch, ["scatter"])
    for kind in ("--mass", "--luminosity", "--phenomena", "--nebula-field"):
        run_plan.run_scatter(_scatter_args(mysql_config, kind))

    def snapshot():
        conn = store.get_connection(mysql_config)
        try:
            return {"stars": conn.execute("SELECT COUNT(*) AS n, COALESCE(ROUND(SUM(position_x_mpc), 3), 0) AS s FROM bright_stars").fetchone(),
                    "phenomena": conn.execute("SELECT COUNT(*) AS n FROM phenomenon_scatter").fetchone()["n"],
                    "field": conn.execute("SELECT COUNT(*) AS n FROM nebula_field").fetchone()["n"],
                    "settings": bring_to_epoch.scatter_settings(conn), "epochs": store.scatter_epochs(conn)}
        finally:
            conn.close()

    before = snapshot()
    assert before["stars"]["n"] > 0 and before["epochs"] == {"mass": 2, "luminosity": 2, "phenomena": 2, "nebula-field": 2}
    galaxy.execute("UPDATE scatter_epochs SET generator_epoch = 1")
    galaxy.commit()
    said = []
    result = bring_to_epoch.run(mysql_config, None, say=said.append)
    assert not result["failed"] and not result["refused"], said
    after = snapshot()
    assert after["epochs"] == before["epochs"], said
    assert {key: after[key] for key in ("stars", "phenomena", "field", "settings")} == \
        {key: before[key] for key in ("stars", "phenomena", "field", "settings")}
    again = []
    assert bring_to_epoch.run(mysql_config, None, say=again.append)["redone"] == 0
    assert again == ["The galaxy is at epoch 2: nothing to bring up."]
