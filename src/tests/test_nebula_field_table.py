"""
DB.24: the `nebula_field` table. `planetgen plan` writes every kept class M and N cloud of the galaxy's molecular
cloud field into it from the seed, the rows equal what `clouds_reaching` gives a sector, a sector's save marks the
clouds it stores, and a plan run over a galaxy that already holds some links them back.
"""

import math

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import nebula_field
from planetgen.galaxy.geometry import sector_position_pc
from planetgen.generation import run_galaxy, run_plan

from tests.test_nebula_field import SHAPE
from tests.test_phenomenon_scatter import E_VALUE, EDGE_PC, EXTENTS, _plan_args

pytestmark = pytest.mark.db

REACH_PC = EDGE_PC * math.sqrt(3) / 2


GALAXY_SEED = bytes(range(16))
"""bytes: The test galaxy's seed, fixed. A new database drew a random one, and about one in sixteen drew a field with
no class M or N cloud reaching any of the sampled sectors, so the tests that need one failed on those (TEST.147)."""


def _seed_galaxy(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=8,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config, galaxy_seed=GALAXY_SEED)
    store.replace_galaxy_layers(EXTENTS, config=mysql_config)


@pytest.fixture(autouse=True)
def plenty_of_clouds(monkeypatch):
    """The test galaxy is tiny: raise the cloud rate so its field holds clouds of both table classes (a Milky Way sized shape, so the gas is there)."""
    scales = dict(tuning.PHENOMENON_RATE_SCALE)
    scales["molecular-cloud"] = 20.0
    monkeypatch.setattr(tuning, "PHENOMENON_RATE_SCALE", scales)


def _seed(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return store.get_galaxy_seed(conn)
    finally:
        conn.close()


def _table(mysql_config, where="1 = 1"):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute(f"SELECT * FROM nebula_field WHERE {where} ORDER BY cell_i, cell_j, cell_k, n").fetchall()
    finally:
        conn.close()


def _centers(rows):
    return {(row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]): row for row in rows}


def _sectors():
    return [(ring, layer, slot) for layer, outer in EXTENTS for ring in (0, outer // 2, outer) for slot in (0,)]


def test_the_table_holds_the_kept_classes_m_and_n_and_the_rows_equal_the_fields_clouds(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    rows = _table(mysql_config)
    assert len(rows) > 20
    assert {row["nebula_class"] for row in rows} == set(nebula_field.TABLE_CLASSES)
    by_center = _centers(rows)
    checked = 0
    for address in _sectors():
        position = sector_position_pc(*address, EDGE_PC)
        for nebula, center in nebula_field.clouds_reaching(_seed(mysql_config), SHAPE, position, REACH_PC):
            if nebula.nebula_class not in nebula_field.TABLE_CLASSES:
                continue
            row = by_center[center]
            assert row["nebula_class"] == nebula.nebula_class
            assert row["radius_ly"] == nebula.radius_ly
            assert (row["cell_i"], row["cell_j"], row["cell_k"]) == nebula_field.cell_index(center)
            checked += 1
    assert checked
    assert all(row["uid"] is None and not row["is_stored"] for row in rows)


def test_a_saved_sector_marks_the_clouds_it_stores_with_their_object_ids(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = next(address for address in _sectors() if any(
        nebula.nebula_class in nebula_field.TABLE_CLASSES
        for nebula, _ in nebula_field.clouds_reaching(_seed(mysql_config), SHAPE, sector_position_pc(*address, EDGE_PC), REACH_PC)))
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    marked = _table(mysql_config, "is_stored = 1")
    assert marked
    conn = store.get_connection(mysql_config)
    try:
        stored = {(r["center_x_pc"], r["center_y_pc"], r["center_z_pc"]): bytes(r["uid"]) for r in conn.execute(
            "SELECT center_x_pc, center_y_pc, center_z_pc, uid FROM nebulae WHERE nebula_class IN ('M', 'N')"
            " AND center_x_pc IS NOT NULL").fetchall()}
    finally:
        conn.close()
    assert {center: bytes(row["uid"]) for center, row in _centers(marked).items()} == stored


def test_a_plan_run_links_the_clouds_a_galaxy_already_holds(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 1
    address = next(address for address in _sectors() if any(
        nebula.nebula_class in nebula_field.TABLE_CLASSES
        for nebula, _ in nebula_field.clouds_reaching(_seed(mysql_config), SHAPE, sector_position_pc(*address, EDGE_PC), REACH_PC)))
    run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    before = _table(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    assert _table(mysql_config) == before


# --- GEN.205: --redo-nebula-field ------------------------------------------------------------------------------------

def _redo_args(mysql_config, mode):
    args = _plan_args(mysql_config)
    args.redo_nebula_field = mode
    return args


def _nebulae(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return {row["id"]: row for row in conn.execute(
            "SELECT id, nebula_class, center_x_pc, center_y_pc, center_z_pc FROM nebulae").fetchall()}
    finally:
        conn.close()


def _galaxy_at_the_old_rates(mysql_config, monkeypatch):
    """A galaxy planned and filled while every field cloud was kept, then the current keep rates restored."""
    current = dict(tuning.NEBULA_FIELD_KEEP)
    monkeypatch.setattr(tuning, "NEBULA_FIELD_KEEP", {letter: 1.0 for letter in current})
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 1
    for address in _sectors()[:2]:
        run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    monkeypatch.setattr(tuning, "NEBULA_FIELD_KEEP", current)


def test_the_dry_run_counts_the_rejected_clouds_and_changes_nothing(mysql_config, monkeypatch):
    _galaxy_at_the_old_rates(mysql_config, monkeypatch)
    before, table = _nebulae(mysql_config), _table(mysql_config)
    assert before
    result = run_plan.redo_nebula_field(_redo_args(mysql_config, "dry-run"))
    assert result["deleted"] == 0 and sum(result["rejected"].values()) > 0
    assert sum(result["stored"].values()) >= sum(result["rejected"].values())
    assert _nebulae(mysql_config) == before and _table(mysql_config) == table


def test_apply_deletes_the_rejected_clouds_and_draws_the_table_again(mysql_config, monkeypatch):
    _galaxy_at_the_old_rates(mysql_config, monkeypatch)
    before = _nebulae(mysql_config)
    dry = run_plan.redo_nebula_field(_redo_args(mysql_config, "dry-run"))
    result = run_plan.redo_nebula_field(_redo_args(mysql_config, "apply"))
    assert result["deleted"] == sum(dry["rejected"].values()) > 0
    after = _nebulae(mysql_config)
    assert len(after) == len(before) - result["deleted"]
    assert set(after) <= set(before)
    # What stays is exactly what the current rule keeps; run again: nothing more to delete.
    again = run_plan.redo_nebula_field(_redo_args(mysql_config, "dry-run"))
    assert sum(again["rejected"].values()) == 0
    # The table was drawn again at the current rates and links the clouds that stayed.
    table = _table(mysql_config)
    assert table and all(row["nebula_class"] in nebula_field.TABLE_CLASSES for row in table)
    stored_centers = {(r["center_x_pc"], r["center_y_pc"], r["center_z_pc"]) for r in after.values()
                      if r["nebula_class"] in nebula_field.TABLE_CLASSES}
    assert {(r["center_x_pc"], r["center_y_pc"], r["center_z_pc"]) for r in table if r["is_stored"]} == stored_centers


# --- GEN.222: the log says what each slab of the field was given --------------------------------------------------

def test_each_slab_line_lists_its_nebula_classes_and_the_layers_they_reach(mysql_config, monkeypatch):
    """Boss: the nebula field step gave only a final total; each slab says its kinds of nebula and its layers."""
    from planetgen.util import log
    _seed_galaxy(mysql_config)
    messages = []
    monkeypatch.setattr(log, "normal", lambda message, *args, **kwargs: messages.append(message))
    run_plan.scatter_phenomena(_plan_args(mysql_config, "--workers", "1"))
    lines = [m for m in messages if m.startswith("Nebula field, slab ")]
    table = _table(mysql_config)
    assert table and lines
    top = max(abs(layer) for layer, _outer in EXTENTS)
    total = 0
    for line in lines:
        head, tail = line.split(" clouds in ")
        total += int(head.rsplit(": ", 1)[1].replace(",", ""))
        layers, kinds = tail.split(": ", 1)
        numbers = [int(word) for word in layers.replace("layers ", "").replace("layer ", "").split(" to ")]
        assert numbers == sorted(numbers) and all(-top <= number <= top for number in numbers)
        assert " (class M)" in kinds or " (class N)" in kinds
    assert total == len(table)
    for letter in {row["nebula_class"] for row in table}:
        assert any(f"{tuning.NEBULA_CLASSES[letter]['name']} (class {letter})" in line for line in lines), letter


def test_a_slab_reports_the_layers_its_clouds_reach_within_the_galaxys_layers():
    text = run_plan._layer_range_text
    assert text((2, 2)) == "layer 2" and text((-1, 3)) == "layers -1 to 3"
