# tests/test_bughunt_end_to_end.py

"""
Tier 1 bug-hunt coverage: true top-to-bottom pipeline tests -- generate a
system via the real `generate.py system` CLI entry point, confirm it
landed in the database with the values the CLI was asked for, confirm the
JSON API serves the same data back, and confirm the web page renders that
same data (name, star type) correctly. Every other bug-hunt file in this
pass tests one layer in isolation (generation-core, DB persistence, API,
web pages); this file is the one that proves those layers actually agree
with each other end to end, which is exactly what "everything ends up in
the database where it's supposed to, and the web interface is accessing
it properly" (the request this whole pass was scoped from) means taken
literally.
"""

from tests.bughunt_support import mysql_argv, run_cli
from tests.webpage_support import live_api, run_page  # noqa: F401


def test_system_generated_via_cli_is_correct_through_db_api_and_webpage(mysql_config, live_api):
    # 1. Generate via the real CLI entry point -- the same code path a
    # real user's `python generate.py system ...` invocation runs.
    system_name = "Bughunt E2E Test System"
    run_cli(
        "system",
        ["--star-type", "K2V", "--name", system_name, "-planets"] + mysql_argv(mysql_config),
    )

    # 2. Confirm it landed in the database with the right values.
    conn = _db_get_connection(mysql_config)
    try:
        row = conn.execute(
            "SELECT ss.id, s.star_type, s.name AS star_name FROM star_systems ss "
            "JOIN stars s ON s.star_system_id = ss.id WHERE ss.name = ?",
            (system_name,),
        ).fetchone()
    finally:
        conn.close()
    assert row is not None, f"system {system_name!r} was not found in star_systems after CLI generation"
    assert row["star_type"].startswith("K2V"), f"expected a K2V star, got {row['star_type']!r}"
    system_id = row["id"]

    # 3. Confirm the JSON API serves back the same system, with the same
    # star type -- the read path a real client/browser actually uses,
    # not just a direct SQL SELECT.
    import json
    import urllib.request

    with urllib.request.urlopen(f"{live_api}/systems/{system_id}?db={mysql_config.database}") as resp:
        api_body = json.loads(resp.read())
    assert api_body["name"] == system_name
    assert api_body["stars"][0]["star_type"].startswith("K2V")

    # 4. Confirm the web page renders the same name and star type --
    # the actual thing a human visiting the site sees.
    result = run_page(live_api, "system.py", query={"db": mysql_config.database, "id": str(system_id)})
    assert result.status_code == 200
    assert system_name in result.body
    assert "K2V" in result.body


def test_sector_generated_via_cli_is_correct_through_db_api_and_webpage(mysql_config, live_api):
    sector_name = "Bughunt E2E Sector"
    run_cli(
        "sector",
        ["--name", sector_name, "--num-systems", "2", "-planets"] + mysql_argv(mysql_config),
    )

    conn = _db_get_connection(mysql_config)
    try:
        sector_row = conn.execute("SELECT id FROM sectors WHERE name = ?", (sector_name,)).fetchone()
        assert sector_row is not None, f"sector {sector_name!r} was not found in sectors after CLI generation"
        sector_id = sector_row["id"]
        system_count = conn.execute(
            "SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?", (sector_id,)
        ).fetchone()["n"]
    finally:
        conn.close()
    assert system_count == 2, f"expected 2 systems in the sector, found {system_count}"

    import json
    import urllib.request

    with urllib.request.urlopen(f"{live_api}/sectors/{sector_id}?db={mysql_config.database}") as resp:
        api_body = json.loads(resp.read())
    assert api_body["name"] == sector_name
    assert len(api_body["systems"]) == 2

    result = run_page(live_api, "sector.py", query={"db": mysql_config.database, "id": str(sector_id)})
    assert result.status_code == 200
    assert sector_name in result.body


def _db_get_connection(mysql_config):
    from stellarObjects import _db
    return _db.get_connection(mysql_config)
