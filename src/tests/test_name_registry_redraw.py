"""
TEST.85: a name `reserve_system_names` redraws mid-pass (its base had no
decoration left) must not count as a holder of the registry row its new
spelling matches later in the same pass. It did: that row then had one
holder more than its count, the existing count came out -1, and the
sector save failed with "existing_count must be >= 0, got -1".

Takes the `mysql_config` fixture (see `conftest.py`): skipped, not
failed, when no MySQL test server is configured.
"""

from stellarObjects import _db


def test_a_redrawn_name_is_not_counted_against_a_later_row_in_the_same_pass(mysql_config, monkeypatch):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            # A two-word base already held: a second "Aaa Bbb" has no
            # decoration within two words (GEN.46), so it is redrawn.
            assert _db.reserve_system_names(conn, ["Aaa Bbb"]) == [("Aaa Bbb", "Aaa Bbb", None)]
        # The redraw lands on a name later in the same batch (keys are
        # handled in sorted order: "aaa bbb" before "zed").
        monkeypatch.setattr(_db, "_regenerate_star_name", lambda: "Zed")
        with conn:
            reserved = _db.reserve_system_names(conn, ["Aaa Bbb", "Zed"])
        # "Zed" is first held by the batch's second name; the redrawn one
        # arrives after it as the second holder, so the two stay apart.
        names = [r[0] for r in reserved]
        assert names[0] == "Beta Zed"
        assert len(set(names)) == 2
        counts = {row["base_name"]: row["occurrence_count"] for row in
                  conn.execute("SELECT base_name, occurrence_count FROM system_name_registry").fetchall()}
        # "Aaa Bbb" stays one high, as a drawn-again name always has.
        assert counts == {"Aaa Bbb": 2, "Zed": 2}
    finally:
        conn.close()
