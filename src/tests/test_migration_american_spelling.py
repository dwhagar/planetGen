# tests/test_migration_american_spelling.py

"""
DOC.103: migration 0211 gives the stored words that kept a British spelling the American one, with no alias: the hydrosphere `vapour` of a planet
or moon becomes `vapor`, a bright star's population `nuclear_disc` becomes `nuclear_disk`, a belt's origin `outer-disc` becomes `outer-disk`
(in its list of causes too), and the ground column `ground_oxidised_fraction` is renamed `ground_oxidized_fraction`. The CHECKs that list a
renamed value are replaced by ones that list the new word, and a row keeps every other value.
"""

import pytest

from planetgen.db import store

pytestmark = pytest.mark.fresh_database


def test_a_migrated_galaxy_has_no_check_left_not_enforced_by_the_migration(mysql_config):
    store.get_connection(mysql_config).close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons", "bright_stars"):
            rows = conn.execute("SELECT CONSTRAINT_NAME, ENFORCED FROM information_schema.TABLE_CONSTRAINTS WHERE TABLE_SCHEMA = DATABASE()"
                                " AND TABLE_NAME = ? AND CONSTRAINT_TYPE = 'CHECK'", (table,)).fetchall()
            assert rows and all(row["ENFORCED"] == "YES" for row in rows), table
    finally:
        conn.close()
