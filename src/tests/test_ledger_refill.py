"""GEN.1498: `planetgen refill-effect-ledger` writes again the light and ultraviolet rows of a stored system, and leaves the rest of the ledger alone."""

import argparse
import math

from planetgen.cli import generate as generate_cli
from planetgen.db import derived_stamp, store
from planetgen.generation import run_plan
from planetgen.physics import effect_ledger, light_ledger
from tests.test_habitability_fill import _save_systems

LIGHT_CODES = tuple(effect_ledger.EFFECT_CODES[name] for name in light_ledger.EFFECT_NAMES)


def _rows(conn, system_id, codes=LIGHT_CODES):
    """The ledger rows of the effects `codes`, by (body, effect, source)."""
    marks = ", ".join("?" for _ in codes)
    found = conn.execute(f"SELECT * FROM effect_ledger WHERE star_system_id = ? AND effect IN ({marks})", (system_id, *codes)).fetchall()
    return {(r["body_type"], r["body_id"], r["effect"], r["source_type"], r["source_id"]): r for r in found}


def _stamp(conn, system_id):
    row = conn.execute("SELECT ledger_stamp FROM system_summary WHERE star_system_id = ?", (system_id,)).fetchone()
    return None if row is None else row["ledger_stamp"]


def _args(mysql_config):
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_fill_arguments(parser)
    parser.set_defaults(command='refill-effect-ledger')
    args = parser.parse_args(["--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
                              "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database, "--workers", "1"])
    generate_cli.validate_fill_args(args, parser)
    return args


def test_a_saved_system_holds_the_ledger_stamp_of_the_running_code(mysql_config):
    saved = _save_systems(mysql_config, wanted=(None,))
    conn = store.get_connection(mysql_config)
    try:
        assert _stamp(conn, saved[None]) == derived_stamp.current("ledger")
    finally:
        conn.close()


def test_a_refill_gives_back_the_rows_generation_wrote_and_a_rerun_does_nothing(mysql_config):
    saved = _save_systems(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        before = {sid: _rows(conn, sid) for sid in saved.values()}
        other = {sid: _rows(conn, sid, tuple(c for c in effect_ledger.EFFECT_CODES.values() if c not in LIGHT_CODES)) for sid in saved.values()}
        with conn:
            for sid in saved.values():
                conn.execute("DELETE FROM effect_ledger WHERE star_system_id = ? AND effect IN (?, ?)", (sid, *LIGHT_CODES))
                conn.execute("UPDATE system_summary SET ledger_stamp = NULL WHERE star_system_id = ?", (sid,))
    finally:
        conn.close()
    args = _args(mysql_config)
    first = run_plan.fill_effect_ledger(args)
    assert first == sum(len(rows) for rows in before.values())
    assert run_plan.fill_effect_ledger(args) == 0          # every system holds the running code's stamp, so none is read
    conn = store.get_connection(mysql_config)
    try:
        for sid, rows in before.items():
            after = _rows(conn, sid)
            assert set(after) == set(rows)
            for key, row in rows.items():
                assert math.isclose(after[key]["share"], row["share"], rel_tol=1e-3, abs_tol=1e-9), key
            assert set(_rows(conn, sid, tuple(c for c in effect_ledger.EFFECT_CODES.values() if c not in LIGHT_CODES))) == set(other[sid])
            assert _stamp(conn, sid) == derived_stamp.current("ledger")
    finally:
        conn.close()


def test_a_system_whose_stamp_is_another_codes_is_refilled_and_the_others_are_not(mysql_config):
    saved = _save_systems(mysql_config)
    stale = saved["close"]
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE system_summary SET ledger_stamp = ledger_stamp + 1 WHERE star_system_id = ?", (stale,))
            conn.execute("DELETE FROM effect_ledger WHERE star_system_id = ? AND effect IN (?, ?)", (saved["wide"], *LIGHT_CODES))
        wide_before = _rows(conn, saved["wide"])
        stale_rows = _rows(conn, stale)
    finally:
        conn.close()
    assert run_plan.fill_effect_ledger(_args(mysql_config)) == len(stale_rows)
    conn = store.get_connection(mysql_config)
    try:
        assert _rows(conn, saved["wide"]) == wide_before == {}      # the current system was not read, so what was removed stays removed
        assert _stamp(conn, stale) == derived_stamp.current("ledger")
    finally:
        conn.close()
