"""
MAP.172: the stored nebula aggregates (`nebula_block_stats`, `nebula_block_classes`, `nebula_block_state`;
`db/nebula_blocks.py`). A rebuild holds the `nebula_field` clouds and the stored nebulae the field does not hold, the
nebulae stored since are added by `catch_up` to exactly what a rebuild would give, a plan makes the totals with the
field, and `plan --fill-nebula-blocks` makes them again.
"""

import numpy as np
import pytest

from planetgen import tuning
from planetgen.db import nebula_blocks, query, store
from planetgen.galaxy import nebula_blocks as rules, nebula_field
from planetgen.generation import run_galaxy, run_plan
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.physics.units import ly_to_pc

from tests.test_nebula_field_table import _plan_args, _seed_galaxy

pytestmark = pytest.mark.db

EDGE_PC = float(tuning.DEFAULT_SECTOR_EDGE_PC)


def _clouds(seed, count, letters="MN"):
    """`count` field clouds `(class, radius in ly, (x, y, z) in pc)` in a box 600 pc across."""
    rng = np.random.default_rng(seed)
    where = np.array([800.0, 200.0, 50.0]) + rng.uniform(-300, 300, (count, 3))
    return [(str(rng.choice(list(letters))), float(rng.uniform(5, 45)), tuple(float(v) for v in w))
            for w in where]


def _plant_field(conn, clouds):
    rows = []
    for number, (letter, radius_ly, (x, y, z)) in enumerate(clouds):
        i, j, k = nebula_field.cell_index((x, y, z))
        rows.append((i, j, k, number, letter, radius_ly, x, y, z))
    store.insert_nebula_field(conn, rows)
    conn.commit()


def _store_nebula(conn, letter, radius_ly, center):
    nebula = Nebula(SystemConfig(), nebula_class=letter)
    nebula.radius_ly = radius_ly
    nebula_id = store.insert_nebula(conn, nebula, placement={
        "center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
        "galactic_radius_pc": float(np.linalg.norm(center))})
    conn.commit()
    return nebula_id


def _link_to_field(conn, nebula_id):
    """What a sector's save does for a stored field cloud (`_mark_nebula_field_stored`): the nebula gets its object ID
    (there is no ID issuer on a bare test connection, so one is made up) and the field row takes it."""
    conn.execute("UPDATE nebulae SET uid = ? WHERE id = ?", (nebula_id.to_bytes(12, "big"), nebula_id))
    conn.execute("UPDATE nebula_field f JOIN nebulae n ON n.center_x_pc = f.center_x_pc AND n.center_y_pc = f.center_y_pc"
                 " AND n.center_z_pc = f.center_z_pc SET f.is_stored = 1, f.uid = n.uid WHERE n.id = ?", (nebula_id,))
    conn.commit()


def _arrays(clouds):
    return (rules.class_codes([c[0] for c in clouds]), np.array([c[2][0] for c in clouds]),
            np.array([c[2][1] for c in clouds]), np.array([c[2][2] for c in clouds]),
            np.array([ly_to_pc(c[1]) for c in clouds]))


def _table(conn):
    return {(r["m"], r["shell_index"], r["band_index"], r["wedge_index"]): r for r in conn.execute(
        "SELECT * FROM nebula_block_stats").fetchall()}


def _same_totals(found, expected):
    """The rows of two tables agree in everything but the bounding radius, which a catch-up only keeps as an upper bound."""
    assert found.keys() == expected.keys()
    for key, row in found.items():
        want = expected[key]
        assert row["nebula_count"] == want["nebula_count"], key
        for column in ("volume_pc3", "volume_diffuse_pc3", "volume_emission_pc3", "volume_reflection_pc3",
                       "volume_planetary_pc3", "volume_dark_pc3", "volume_remnant_pc3", "cover_volume_pc3", "cover"):
            assert row[column] == pytest.approx(want[column], rel=1e-9, abs=1e-9), (key, column)
        assert row["dominant_class"] == want["dominant_class"], key
        for column in ("centroid_x_pc", "centroid_y_pc", "centroid_z_pc"):
            if want[column] is None:
                assert row[column] is None
            else:
                assert row[column] == pytest.approx(want[column], rel=1e-9, abs=1e-6), (key, column)
        if want["bound_radius_pc"] is not None:
            assert row["bound_radius_pc"] >= want["bound_radius_pc"] - 1e-6, key


def test_a_rebuild_holds_the_field_and_the_stored_nebulae_the_field_does_not_hold(mysql_config):
    field = _clouds(1, 40)
    others = [("C", 3.0, (805.0, 195.0, 52.0)), ("F", 8.0, (830.0, 210.0, 40.0)), ("H", 2.5, (770.0, 190.0, 60.0)),
              ("P", 1.5, (810.0, 205.0, 45.0)), ("L", 12.0, (790.0, 215.0, 55.0))]
    conn = store.get_connection(mysql_config)
    try:
        _plant_field(conn, field)
        for letter, radius_ly, center in others:
            _store_nebula(conn, letter, radius_ly, center)
        for letter, radius_ly, center in field[:3]:    # stored field clouds: the field already counts them
            _link_to_field(conn, _store_nebula(conn, letter, radius_ly, center))
        assert nebula_blocks.through_id(conn) is None
        written = nebula_blocks.rebuild(conn, EDGE_PC)
        conn.commit()
        assert written == len(_table(conn))
        assert nebula_blocks.through_id(conn) == conn.execute("SELECT MAX(id) AS n FROM nebulae").fetchone()["n"]
        expected = {}
        everything = field + others
        codes, x, y, z, radius = _arrays(everything)
        for m, level in rules.all_rows(x, y, z, radius, codes, EDGE_PC).items():
            for row in level.blocks:
                expected[(m, row["shell_index"], row["band_index"], row["wedge_index"])] = row
        found = _table(conn)
        assert found.keys() == expected.keys()
        for key, row in found.items():
            want = expected[key]
            assert row["nebula_count"] == want["nebula_count"], key
            assert row["cover"] == pytest.approx(want["cover"], rel=1e-9, abs=1e-12), key
            assert row["dominant_class"] == want["dominant_class"], key
        level_three = [row for key, row in found.items() if key[0] == 3]
        assert sum(row["nebula_count"] for row in level_three) == len(everything)
        classes = conn.execute("SELECT SUM(volume_pc3) AS v FROM nebula_block_classes WHERE m = 3").fetchone()["v"]
        assert classes == pytest.approx(sum(row["volume_pc3"] for row in level_three), rel=1e-9)
    finally:
        conn.close()


def test_a_catch_up_adds_what_a_rebuild_would_hold_and_leaves_the_field_clouds_out(mysql_config):
    field = _clouds(2, 30)
    first = [("C", 3.0, (805.0, 195.0, 52.0)), ("F", 8.0, (830.0, 210.0, 40.0))]
    later = [("G", 5.0, (801.0, 199.0, 51.0)), ("H", 2.5, (770.0, 190.0, 60.0)), ("M", 40.0, (815.0, 190.0, 49.0)),
             ("L", 12.0, (790.0, 215.0, 55.0))]
    conn = store.get_connection(mysql_config)
    try:
        _plant_field(conn, field)
        for letter, radius_ly, center in first:
            _store_nebula(conn, letter, radius_ly, center)
        nebula_blocks.rebuild(conn, EDGE_PC)
        conn.commit()
        assert nebula_blocks.catch_up(conn, EDGE_PC) == 0
        for letter, radius_ly, center in later:
            _store_nebula(conn, letter, radius_ly, center)
        for letter, radius_ly, center in field[3:5]:    # a field cloud stored after the build is already counted
            _link_to_field(conn, _store_nebula(conn, letter, radius_ly, center))
        assert nebula_blocks.catch_up(conn, EDGE_PC) == len(later)
        caught_up = _table(conn)
        assert nebula_blocks.catch_up(conn, EDGE_PC) == 0
        classes_caught_up = conn.execute("SELECT * FROM nebula_block_classes ORDER BY m, shell_index, band_index,"
                                         " wedge_index, nebula_class").fetchall()
        nebula_blocks.rebuild(conn, EDGE_PC)
        conn.commit()
        _same_totals(caught_up, _table(conn))
        rebuilt = conn.execute("SELECT * FROM nebula_block_classes ORDER BY m, shell_index, band_index, wedge_index,"
                               " nebula_class").fetchall()
        assert [(r["m"], r["shell_index"], r["band_index"], r["wedge_index"], r["nebula_class"]) for r in rebuilt] == [
            (r["m"], r["shell_index"], r["band_index"], r["wedge_index"], r["nebula_class"]) for r in classes_caught_up]
        assert [r["volume_pc3"] for r in classes_caught_up] == pytest.approx([r["volume_pc3"] for r in rebuilt])
    finally:
        conn.close()


def test_a_catch_up_before_any_build_does_nothing(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        _store_nebula(conn, "C", 3.0, (805.0, 195.0, 52.0))
        assert nebula_blocks.catch_up(conn, EDGE_PC) is None
        assert conn.execute("SELECT COUNT(*) AS n FROM nebula_block_stats").fetchone()["n"] == 0
    finally:
        conn.close()


def test_read_gives_the_shares_the_dominant_class_and_the_cover_of_a_block(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        _plant_field(conn, _clouds(3, 25))
        _store_nebula(conn, "C", 3.0, (805.0, 195.0, 52.0))
        nebula_blocks.rebuild(conn, EDGE_PC)
        conn.commit()
        rows = _table(conn)
        key = max((k for k in rows if k[0] == 27), key=lambda k: rows[k]["cover"])
        found = nebula_blocks.read(conn, 27, [key[1:], (999999, 0, 0)])
        assert list(found) == [key[1:]]
        block = found[key[1:]]
        assert block["cover"] == rows[key]["cover"] and 0 < block["cover"] <= 1
        assert block["count"] == rows[key]["nebula_count"]
        assert sum(block["shares"].values()) == pytest.approx(1.0) or block["count"] == 0
        assert block["dominant_class"] in ("M", "N", "C", None)
        assert nebula_blocks.read(conn, 27, []) == {}
    finally:
        conn.close()


def test_a_plan_makes_the_totals_with_the_field_and_fill_nebula_blocks_makes_them_again(mysql_config, monkeypatch):
    scales = dict(tuning.PHENOMENON_RATE_SCALE)
    scales["molecular-cloud"] = 40.0  # the default galaxy is a small spiral (GEN.603): enough clouds for a field of more than 20
    monkeypatch.setattr(tuning, "PHENOMENON_RATE_SCALE", scales)
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    conn = store.get_connection(mysql_config)
    try:
        field = conn.execute("SELECT COUNT(*) AS n FROM nebula_field").fetchone()["n"]
        assert field > 20
        made = _table(conn)
        level_three = [row for key, row in made.items() if key[0] == 3]
        assert sum(row["nebula_count"] for row in level_three) == field
        assert nebula_blocks.through_id(conn) is not None
        conn.execute("DELETE FROM nebula_block_stats")
        conn.execute("DELETE FROM nebula_block_classes")
        conn.commit()
    finally:
        conn.close()
    args = _plan_args(mysql_config, "--fill-nebula-blocks")
    assert args.fill_nebula_blocks
    assert run_plan.fill_nebula_blocks(args) == len(made)
    conn = store.get_connection(mysql_config)
    try:
        _same_totals(_table(conn), made)
    finally:
        conn.close()


def test_a_galaxy_run_adds_the_nebulae_it_stored_once_at_its_end(mysql_config):
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        nebula_blocks.rebuild(conn, query.galaxy_edge_pc(conn))
        conn.commit()
        _store_nebula(conn, "C", 3.0, (805.0, 195.0, 52.0))
    finally:
        conn.close()
    assert run_galaxy.nebula_blocks_after_run(_plan_args(mysql_config)) == 1
    conn = store.get_connection(mysql_config)
    try:
        assert sum(r["nebula_count"] for k, r in _table(conn).items() if k[0] == 3) == 1
    finally:
        conn.close()
    assert run_galaxy.nebula_blocks_after_run(_plan_args(mysql_config)) == 0
