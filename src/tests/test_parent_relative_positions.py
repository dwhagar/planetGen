# tests/test_parent_relative_positions.py

"""
TEST.241 (section 5 point 3 of docs/design/unified-system-barycentre.md): nothing inside a system stores an absolute sector
or galactic position. A table whose rows are tied to a system or a body (`point_columns.BODY_PARENT_COLUMNS`) and that
holds a point stores it in the `system` frame, parent-relative; the few tables that do not are listed with their reason
and way out in `point_columns.ABSOLUTE_WITH_A_PARENT`, and a listed table that stopped needing it is refused.
"""

import sqlalchemy as sa

from planetgen.db import models, point_columns


def _columns(table):
    return {column.name for column in table.columns}


def holds_a_point(table):
    return point_columns.holds_a_point(_columns(table))


def absolute_positions_inside_systems(metadata, layouts=None, allowed=None):
    """The tables with a parent link whose point is not stored in the `system` frame and that are not allowed to."""
    layouts = point_columns.LAYOUTS if layouts is None else layouts
    allowed = point_columns.ABSOLUTE_WITH_A_PARENT if allowed is None else allowed
    found = []
    for name, table in sorted(metadata.tables.items()):
        names = _columns(table)
        if not (names & set(point_columns.BODY_PARENT_COLUMNS)) or not point_columns.holds_a_point(names):
            continue
        if name in point_columns.NOT_BODIES:
            continue
        layout = layouts.get(name)
        if layout is None or layout.frame != "system":
            if name not in allowed:
                found.append(f"{name}: has {sorted(names & set(point_columns.BODY_PARENT_COLUMNS))} and a "
                             f"{layout.frame if layout else 'no layout'}-frame point")
    return found


def test_no_table_inside_a_system_stores_an_absolute_position():
    assert absolute_positions_inside_systems(models.metadata) == []


def test_a_table_given_an_absolute_position_inside_a_system_is_refused():
    metadata = sa.MetaData()
    sa.Table("rings", metadata, sa.Column("id", sa.Integer, primary_key=True), sa.Column("planet_id", sa.Integer),
             sa.Column("center_x_pc", sa.Float), sa.Column("center_y_pc", sa.Float), sa.Column("center_z_pc", sa.Float))
    layouts = {"rings": point_columns.PointLayout(
        table="rings", frame="galaxy", position=("center_x_pc", "center_y_pc", "center_z_pc"), position_unit="pc",
        velocity=("v_x", "v_y", "v_z"), mass="mass_kg", mass_unit="kg", epoch="epoch_unix", due=None)}
    [problem] = absolute_positions_inside_systems(metadata, layouts, {})
    assert problem.startswith("rings:") and "galaxy-frame" in problem
    layouts["rings"] = point_columns.PointLayout(
        table="rings", frame="system", position=("center_x_pc", "center_y_pc", "center_z_pc"), position_unit="pc",
        velocity=("v_x", "v_y", "v_z"), mass="mass_kg", mass_unit="kg", epoch="epoch_unix", due=None)
    assert absolute_positions_inside_systems(metadata, layouts, {}) == []
    assert absolute_positions_inside_systems(metadata, {}, {"rings": "why"}) == []


def test_an_allowed_table_still_needs_its_exception():
    for name, reason in point_columns.ABSOLUTE_WITH_A_PARENT.items():
        table = models.metadata.tables[name]
        assert _columns(table) & set(point_columns.BODY_PARENT_COLUMNS), f"{name} no longer has a parent link"
        assert holds_a_point(table) and reason
        layout = point_columns.LAYOUTS[name]
        assert layout.frame != "system", f"{name} is stored parent-relative now: drop its exception"


def test_the_bodies_of_a_system_are_stored_in_the_system_frame():
    for name in ("planets", "moons", "comets", "asteroid_belts", "stars", "scatter_members"):
        assert point_columns.LAYOUTS[name].frame == "system", name
