# tests/test_habitability_fill.py

"""
`planetgen rescore-habitability` (GEN.217): the stored habitability numbers of an existing galaxy are
re-derived from the stored conditions, and come out the same as the ones a new generation stored.
"""

import math

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import derived_stamp, store
from planetgen.generation import run_plan
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import habitability_world
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy

COLUMNS = habitability_world.SCORE_FIELDS + ("flux_variation",)
"""tuple: The columns the fill writes."""


SEARCH_SEEDS = 3000
"""int: How many seeds `_save_systems` tries for a system of each kind; a close binary with a scored rocky body is under
1 in 200 (GEN.651: the draws moved, and the 300 it once took no longer hold one)."""


def _save_systems(mysql_config, wanted=("close", "wide", None)):
    """Saved systems of each wanted binary type with a scored rocky body, by `binary_type`."""
    saved = {}
    for seed in range(SEARCH_SEEDS):
        if len(saved) == len(wanted):
            break
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            cfg = SystemConfig()
            system = StarSystem(cfg, galactic_center_dist_ly=26000.0)
        kind = system.binary_type
        if kind in wanted and kind not in saved and any(
                getattr(p, "body_type", None) == "t" and p.phi4 is not None for p in system.planets):
            saved[kind] = store.save_system(system, cfg, config=mysql_config)
    assert set(saved) == set(wanted), "could not generate each kind of system"
    return saved


def _bodies(conn, system_id):
    rows = {}
    for table in ("planets", "moons"):
        for row in conn.execute(f"SELECT * FROM {table} WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall():
            rows[(table, row["id"])] = {name: row[name] for name in COLUMNS} | {"body_type": row["body_type"],
                                                                              "scored": row["energy_flux_w_m2"]}
    return rows


def _wipe(conn, system_id):
    """Scores gone, and the system no longer says the running code scored it (a galaxy scored before a change to the rules)."""
    with conn:
        for table in ("planets", "moons"):
            conn.execute(f"UPDATE {table} SET {', '.join(f'{name} = NULL' for name in COLUMNS)} "
                         "WHERE star_system_id = ? AND body_type <> 'g'", (system_id,))
        conn.execute("DELETE FROM system_summary WHERE star_system_id = ?", (system_id,))


def _stamp(conn, system_id):
    row = conn.execute("SELECT habitability_stamp FROM system_summary WHERE star_system_id = ?", (system_id,)).fetchone()
    return None if row is None else row["habitability_stamp"]


def _same(a, b):
    if a is None or b is None:
        return a is b
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)


def test_a_system_comes_back_with_the_numbers_generation_stored(mysql_config):
    saved = _save_systems(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        for kind, system_id in saved.items():
            before = _bodies(conn, system_id)
            assert any(b["phi4"] is not None for b in before.values())
            _wipe(conn, system_id)
            with conn:
                updated = store.refresh_habitability(conn, system_id)
            after = _bodies(conn, system_id)
            assert updated == sum(1 for b in before.values() if b["scored"] is not None and b["body_type"] != "g")
            for key, row in before.items():
                if row["body_type"] == "g" or row["scored"] is None:
                    continue
                for name in COLUMNS:
                    assert _same(after[key][name], row[name]), (kind, key, name, after[key][name], row[name])
    finally:
        conn.close()


def test_a_double_star_system_gets_its_flux_swing_back_and_a_lone_star_gets_its_orbits(mysql_config):
    saved = _save_systems(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        for kind, system_id in saved.items():
            _wipe(conn, system_id)
            with conn:
                store.refresh_habitability(conn, system_id)
            swings = [b["flux_variation"] for b in _bodies(conn, system_id).values()
                      if b["body_type"] != "g" and b["scored"] is not None]
            assert swings and all(value is not None for value in swings)
            if kind is None:
                # GEN.727: the swing of an eccentric orbit, 0 only for a circle
                assert all(0.0 <= value <= 1.0 for value in swings)
                assert any(value > 0.0 for value in swings)
            else:
                assert any(value > 0.0 for value in swings), kind
    finally:
        conn.close()


def test_a_stale_air_rating_is_corrected_without_touching_anything_else(mysql_config):
    system_id = _save_systems(mysql_config, wanted=(None,))[None]
    conn = store.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT * FROM planets WHERE star_system_id = ? AND body_type = 't' "
                           "AND phi4 IS NOT NULL ORDER BY id LIMIT 1", (system_id,)).fetchone()
        other = {name: row[name] for name in row if name not in COLUMNS and not name.endswith("_at")}
        with conn:
            conn.execute("UPDATE planets SET phi4_air = 0.123, tier_air = 3, phi4_suitability = NULL "
                         "WHERE id = ?", (row["id"],))
            store.refresh_habitability(conn, system_id)
        fresh = conn.execute("SELECT * FROM planets WHERE id = ?", (row["id"],)).fetchone()
        assert _same(fresh["phi4_air"], row["phi4_air"]) and fresh["tier_air"] == row["tier_air"]
        assert fresh["phi4_suitability"] is not None
        assert {name: fresh[name] for name in other} == other
    finally:
        conn.close()


def _fill_args(mysql_config):
    import argparse

    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_fill_arguments(parser)
    parser.set_defaults(command='rescore-habitability')   # the sub-command parser's own default (validate_fill_args reads it)
    args = parser.parse_args(["--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
                              "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
                              "--mysql-database", mysql_config.database, "--workers", "1"])
    generate_cli.validate_fill_args(args, parser)
    return args


def test_the_command_re_scores_every_system_and_a_second_run_changes_nothing(mysql_config):
    saved = _save_systems(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        expected = {sid: _bodies(conn, sid) for sid in saved.values()}
        for system_id in saved.values():
            _wipe(conn, system_id)
    finally:
        conn.close()
    args = _fill_args(mysql_config)
    original = run_plan.HABITABILITY_BATCH_SYSTEMS
    run_plan.HABITABILITY_BATCH_SYSTEMS = 1
    try:
        first = run_plan.fill_habitability(args)
        assert first == sum(1 for rows in expected.values() for b in rows.values()
                            if b["scored"] is not None and b["body_type"] != "g")
        assert run_plan.fill_habitability(args) == 0       # GEN.1479: the systems hold the running code's stamp, so nothing is touched
        args.all = True
        assert run_plan.fill_habitability(args) == first   # --all takes every system, as the command used to
    finally:
        run_plan.HABITABILITY_BATCH_SYSTEMS = original
    conn = store.get_connection(mysql_config)
    try:
        for system_id, rows in expected.items():
            after = _bodies(conn, system_id)
            for key, row in rows.items():
                if row["body_type"] == "g" or row["scored"] is None:
                    continue
                assert all(_same(after[key][name], row[name]) for name in COLUMNS), key
    finally:
        conn.close()


def test_the_fill_takes_no_setting_of_a_scatter():
    import argparse

    parser = argparse.ArgumentParser(prefix_chars="-+")
    generate_cli.add_fill_arguments(parser)
    assert parser.parse_args(["--workers", "2"]).workers == 2
    with pytest.raises(SystemExit):
        parser.parse_args(["--phenomenon-min-mass", "12"])
    with pytest.raises(SystemExit):
        parser.parse_args(["--fill-map-luminosity"])


def test_a_generated_body_stores_the_score_of_its_final_conditions():
    """GEN.217: a lifeless moon used to keep the score it was first given, before its temperature settled, so the
    chip disagreed with the description (which is worked out from the final conditions)."""
    checked = 0
    for seed in range(40):    # GEN.1334: a planet under an envelope is scored for nothing, so the same seeds hold fewer scored bodies
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets + system.secondary_planets:
            for body in [planet] + list(getattr(planet, "moons", [])):
                if getattr(body, "body_type", None) != "t" or getattr(body, "phi4", None) is None:
                    continue
                again = habitability_world.compute(habitability_world.values_from_body(body),
                                                   habitability_world.compact_host_name(getattr(body, "star", None)))
                for name in habitability_world.SCORE_FIELDS:
                    assert _same(getattr(body, name), again[name]), (seed, body.is_moon, name)
                checked += 1
    assert checked > 100


def test_a_saved_system_holds_the_stamp_of_the_running_code(mysql_config):
    saved = _save_systems(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        current = derived_stamp.current("habitability")
        assert {_stamp(conn, system_id) for system_id in saved.values()} == {current}
        _wipe(conn, saved[None])
        assert _stamp(conn, saved[None]) is None
        with conn:
            store.refresh_habitability(conn, saved[None])
        assert _stamp(conn, saved[None]) == current
    finally:
        conn.close()


def test_a_rerun_takes_only_the_systems_that_are_not_current(mysql_config):
    """GEN.1479: of three saved systems, the one scored by other code (or never) is re-scored and the others are not touched."""
    saved = _save_systems(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        stale = saved["close"]
        with conn:
            conn.execute("UPDATE system_summary SET habitability_stamp = habitability_stamp + 1 WHERE star_system_id = ?", (stale,))
            conn.execute("UPDATE planets SET phi4 = -1 WHERE star_system_id = ? AND phi4 IS NOT NULL", (saved["wide"],))
        untouched = conn.execute("SELECT phi4 FROM planets WHERE star_system_id = ? ORDER BY id", (saved["wide"],)).fetchall()
    finally:
        conn.close()
    original = run_plan.HABITABILITY_BATCH_SYSTEMS
    run_plan.HABITABILITY_BATCH_SYSTEMS = 1
    try:
        done = run_plan.fill_habitability(_fill_args(mysql_config))
    finally:
        run_plan.HABITABILITY_BATCH_SYSTEMS = original
    conn = store.get_connection(mysql_config)
    try:
        current = derived_stamp.current("habitability")
        assert done == conn.execute("SELECT (SELECT COUNT(*) FROM planets WHERE star_system_id = ? AND energy_flux_w_m2 IS NOT NULL AND body_type <> 'g')"
                                    " + (SELECT COUNT(*) FROM moons WHERE star_system_id = ? AND energy_flux_w_m2 IS NOT NULL AND body_type <> 'g') AS n",
                                    (stale, stale)).fetchone()["n"]
        assert _stamp(conn, stale) == current
        assert conn.execute("SELECT phi4 FROM planets WHERE star_system_id = ? ORDER BY id", (saved["wide"],)).fetchall() == untouched
    finally:
        conn.close()


def test_a_changed_formula_is_a_different_stamp():
    assert derived_stamp.current("habitability") == derived_stamp.current("habitability")
    sources = derived_stamp.sources("habitability")
    assert "planetgen.physics.habitability_world" in sources and "planetgen.physics.life_classes" in sources   # a lazy import counts
    assert not any(not name.startswith(("planetgen.physics.", "planetgen.tuning")) for name in sources)
    before = derived_stamp.current("habitability")
    original = dict(derived_stamp.RULE_VERSION)
    try:
        derived_stamp.current.cache_clear()
        derived_stamp.RULE_VERSION["habitability"] += 1
        assert derived_stamp.current("habitability") != before
    finally:
        derived_stamp.RULE_VERSION.update(original)
        derived_stamp.current.cache_clear()
    assert derived_stamp.current("habitability") == before
