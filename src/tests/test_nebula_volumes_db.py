"""
MAP.365: the points that draw the nebulae of a Galaxy Map node, read from the database (`db/nebula_volumes.py`): the nebulae that reach
the block, only the lattice points that lie in it, the same bytes every time, nothing for a block no nebula reaches, a 400-style
`ValueError` for a bad key.
"""

import math

import pytest

from planetgen import tuning
from planetgen.db import nebula_volumes, query, store
from planetgen.galaxy import nebula_volume, structure_draw
from planetgen.galaxy.drill import DrillBlock, drill_chain_of, format_drill_key
from planetgen.galaxy.spherical import sector_address_at

from tests import nebula_density_stand_in as stand_in
from tests.test_bright_star_scatter import toy_shape
from tests.test_nebula_blocks_db import _store_nebula

pytestmark = pytest.mark.db

EDGE_PC = float(tuning.DEFAULT_SECTOR_EDGE_PC)
CENTRE = (805.0, 195.0, 52.0)
SEED = bytes([39]) * 16
SHAPE = toy_shape(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0, bulge_amplitude=2.0, arm_count=2,
                  pitch_angle_rad=math.radians(15), arm_amplitude=0.4, structure=structure_draw.draw_structure(SEED))


class _Model:
    """The pieces of `galaxy.nebula_density` the reader asks for, with the stand-in density."""
    FORM_LEVEL = stand_in.FORM_LEVEL
    MAX_DETAIL = stand_in.MAX_DETAIL
    density_at = staticmethod(lambda seed, shape, points, field, detail: stand_in.density(points, field, detail))


@pytest.fixture
def planned(mysql_config):
    """A galaxy with a plan (a seed and a shape) and one stored 20 pc nebula."""
    conn = store.get_connection(mysql_config)
    try:
        if query.get_galaxy_shape(conn, light=True) is None:
            store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=2000.0,
                                    config=mysql_config, cell_star_ceiling=1e12, galaxy_seed=SEED)
        _store_nebula(conn, "C", 65.0, CENTRE)
        conn.commit()
    finally:
        conn.close()
    return mysql_config


def _keys(m):
    chain = drill_chain_of(*sector_address_at(CENTRE, EDGE_PC))
    return format_drill_key(next(b for b in chain if b.m == m))


def test_a_block_draws_the_points_of_the_nebulae_in_it(planned):
    conn = store.get_connection(planned)
    try:
        key = _keys(27)
        first = nebula_volumes.group_volume(conn, key, _Model)
        again = nebula_volumes.group_volume(conn, key, _Model)
        assert first == again
        assert first[:4] == nebula_volume.MAGIC
        count = int.from_bytes(first[6:10], "little")
        assert count > 20
        far = nebula_volumes.group_volume(conn, "27.9.0.1", _Model)
        assert int.from_bytes(far[6:10], "little") == 0
        with pytest.raises(ValueError):
            nebula_volumes.group_volume(conn, "not-a-key", _Model)
    finally:
        conn.close()


def _points(packed):
    """`{(ix, iy, iz): strength}` of a packed volume."""
    import numpy as np

    n = int.from_bytes(packed[6:10], "little")
    body = packed[14:]
    ix, iy, iz = (np.frombuffer(body[4 * n * k:4 * n * (k + 1)], dtype="<i4").tolist() for k in range(3))
    strength = list(body[12 * n:13 * n])
    return dict(zip(zip(ix, iy, iz), strength))


def test_one_nebula_alone_is_drawn_with_the_points_its_node_gives_it(planned):
    """MAP.142: the volume of one nebula holds the whole of its cloud and, at the same level, every point its node draws for it at the same
    strength, so a nebula looks the same on the Galaxy Map, the Sector Map and its own page."""
    conn = store.get_connection(planned)
    try:
        nebula_id = conn.execute("SELECT id FROM nebulae ORDER BY id DESC LIMIT 1").fetchone()["id"]
        alone = nebula_volumes.single_nebula_volume(conn, nebula_id, _Model)
        assert alone == nebula_volumes.single_nebula_volume(conn, nebula_id, _Model)
        assert alone[:4] == nebula_volume.MAGIC
        whole = _points(alone)
        assert len(whole) > 20
        block = _points(nebula_volumes.group_volume(conn, _keys(27), _Model))
        assert block and alone[5] == nebula_volumes.group_volume(conn, _keys(27), _Model)[5]
        assert all(whole.get(cell) == strength for cell, strength in block.items())
        assert len(whole) >= len(block)
        with pytest.raises(ValueError):
            nebula_volumes.single_nebula_volume(conn, 99999999, _Model)
    finally:
        conn.close()


def test_a_nebula_in_a_galaxy_never_planned_is_header_only(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        if query.get_galaxy_shape(conn, light=True) is not None:
            pytest.skip("this database already holds a plan")
        nebula_id = _store_nebula(conn, "C", 65.0, CENTRE)
        packed = nebula_volumes.single_nebula_volume(conn, nebula_id, _Model)
        assert int.from_bytes(packed[6:10], "little") == 0
    finally:
        conn.close()
