# tests/test_galaxy_family_columns_store.py

"""
GEN.379: the family and component numbers of `galaxy_shape` (DB.66) are written by `save_galaxy_shape` and read by
`get_galaxy_shape`. A bar or thick disk number that is NULL is the structure draw's own value (never a fixed number), a
number pins it in the shape and its structure, and saving a plan again replaces what an earlier one stored.
"""

import math

import pytest

from planetgen.db import store
from planetgen.galaxy import density
from planetgen.galaxy import shape_draw

BASE = dict(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=700.0, bulge_amplitude=3.11,
            arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
FAMILY = ("family", "bar_present", "bar_axis_ratio_b", "bar_axis_ratio_c", "bar_angle_rad", "thick_disk_amplitude",
          "spheroid_amplitude", "spheroid_sersic_n", "spheroid_effective_radius_pc", "spheroid_axis_ratio_b",
          "spheroid_axis_ratio_c", "spheroid_orientation_rad", "body_clump_share")


def _save(config, shape):
    store.save_galaxy_shape(shape, 20.0, 3, 1.0e4, config=config)


def _read(config):
    conn = store.get_connection(config)
    try:
        return store.get_galaxy_shape(conn).shape
    finally:
        conn.close()


def _row(config):
    conn = store.get_connection(config)
    try:
        return conn.execute("SELECT * FROM galaxy_shape WHERE id = 1").fetchone()
    finally:
        conn.close()


def test_a_plan_with_nothing_pinned_stores_nulls_and_reads_the_structure_draw(mysql_config):
    shape = density.build_galaxy_shape(**BASE)
    _save(mysql_config, shape)
    row = _row(mysql_config)
    assert row["family"] == "spiral"
    for column in ("bar_present", "bar_axis_ratio_b", "bar_axis_ratio_c", "bar_angle_rad", "thick_disk_amplitude"):
        assert row[column] is None, column
    assert (row["spheroid_amplitude"], row["spheroid_axis_ratio_b"], row["body_clump_share"]) == (0.0, 1.0, 0.0)
    read = _read(mysql_config)
    assert read.family == "spiral" and read.bar_present is None and read.thick_disk_amplitude is None
    assert read.spheroid_orientation_rad == (0.0, 0.0, 0.0)


def test_the_stored_numbers_come_back_and_pin_the_structure(mysql_config):
    shape = density.build_galaxy_shape(**BASE)._replace(
        family="lenticular", bar_present=False, bar_axis_ratio_b=0.5, bar_axis_ratio_c=0.25, bar_angle_rad=0.7,
        thick_disk_amplitude=0.06, spheroid_amplitude=0.0, spheroid_sersic_n=3.5, spheroid_effective_radius_pc=1800.0,
        spheroid_axis_ratio_b=0.8, spheroid_axis_ratio_c=0.6, spheroid_orientation_rad=(0.1, 0.2, 0.3),
        body_clump_share=0.25)
    _save(mysql_config, shape)
    read = _read(mysql_config)
    for name in FAMILY:
        assert getattr(read, name) == getattr(shape, name), name
    assert read.structure.bar_present is False
    assert (read.structure.bar_axis_ratio_y, read.structure.bar_axis_ratio_z) == (0.5, 0.25)
    assert read.structure.bar_angle_deg == pytest.approx(math.degrees(0.7))
    assert density.thick_disk_amplitude(read) == 0.06
    assert {c.name: c.amplitude for c in density.density_components(read)}["thick_disk"] == 0.06


def test_saving_a_plan_again_replaces_what_the_earlier_one_stored(mysql_config):
    _save(mysql_config, density.build_galaxy_shape(**BASE)._replace(
        family="irregular", bar_present=True, thick_disk_amplitude=0.0, body_clump_share=0.5))
    assert _row(mysql_config)["family"] == "irregular"
    _save(mysql_config, density.build_galaxy_shape(**BASE))
    row = _row(mysql_config)
    assert row["family"] == "spiral" and row["bar_present"] is None and row["thick_disk_amplitude"] is None
    assert row["body_clump_share"] == 0.0


def test_the_drawn_structure_is_what_a_null_reads_back_as(mysql_config):
    _save(mysql_config, density.build_galaxy_shape(**BASE))
    conn = store.get_connection(mysql_config)
    try:
        info = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    # the structure is the seed's draw: the NULL columns changed nothing in it
    assert info.shape.structure == shape_draw.structure_for(info.galaxy_seed, None, None)
