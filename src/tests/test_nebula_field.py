"""
The galaxy's molecular cloud field (GEN.47, `nebulaField`): dark-family
clouds drawn per field cell from the galaxy seed, so every sector and
every run agrees where they are; more in the arms and near the plane; a
generated neighborhood finds them, and each is stored once however many
sectors it reaches.
"""

import math

import pytest

from planetgen.generation import run_galaxy
from planetgen.db import store as _db
from planetgen.galaxy import nebula_field
from planetgen import tuning
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.geometry import sector_address_at, sector_position_pc
from planetgen.physics.units import ly_to_pc

from tests.bughunt_support import forced_system_config

# Milky Way scale: disk scale length 2.8 kpc, the solar circle at 2.82 of them.
SHAPE = build_galaxy_shape(2800.0, 350.0, 200.0, 1.0, 2, math.radians(15.0), 0.4)
SOLAR_RADIUS_PC = 2.82 * SHAPE.disk_scale_length_pc
ARM_ANGLE = math.log(SOLAR_RADIUS_PC / SHAPE.disk_scale_length_pc) / math.tan(SHAPE.pitch_angle_rad)
SEED = bytes(range(16))
EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
REACH_PC = EDGE_PC * math.sqrt(3) / 2


def _on_ring(angle, z=0.0, radius=SOLAR_RADIUS_PC):
    return (radius * math.cos(angle), radius * math.sin(angle), z)


def _key(cloud):
    nebula, center = cloud
    return (center, nebula.nebula_class, nebula.radius_ly, nebula.name)


def test_gas_is_densest_on_the_arms_and_near_the_plane():
    arm = nebula_field.gas_factor(_on_ring(ARM_ANGLE), SHAPE)
    between = nebula_field.gas_factor(_on_ring(ARM_ANGLE + math.pi / 2), SHAPE)
    high = nebula_field.gas_factor(_on_ring(ARM_ANGLE, z=300.0), SHAPE)
    assert arm > 1.5 > 1.0 > between > 0.0
    assert high < 0.05 * arm
    assert nebula_field.gas_factor((10.0, 10.0, 0.0), SHAPE) <= tuning.NEBULA_FIELD_MAX_GAS_FACTOR


def test_only_dark_family_clouds_and_the_same_ones_every_time():
    clouds = []
    for step in range(40):
        clouds += nebula_field.cell_clouds(SEED, SHAPE, nebula_field.cell_index(_on_ring(ARM_ANGLE + step * 0.01)))
    assert clouds
    assert {nebula.nebula_class for nebula, _center in clouds} <= set(nebula_field.FIELD_CLASSES)
    again = []
    for step in reversed(range(40)):
        again += nebula_field.cell_clouds(SEED, SHAPE, nebula_field.cell_index(_on_ring(ARM_ANGLE + step * 0.01)))
    assert sorted(map(_key, clouds)) == sorted(map(_key, again))
    other = []
    for step in range(40):
        other += nebula_field.cell_clouds(bytes(16), SHAPE, nebula_field.cell_index(_on_ring(ARM_ANGLE + step * 0.01)))
    assert sorted(map(_key, clouds)) != sorted(map(_key, other))


def test_the_field_leaves_the_sector_stream_alone():
    import random
    random.seed(7)
    expected = [random.random() for _ in range(3)]
    random.seed(7)
    nebula_field.clouds_reaching(SEED, SHAPE, _on_ring(ARM_ANGLE), REACH_PC)
    assert [random.random() for _ in range(3)] == expected


def test_neighboring_sectors_see_the_same_cloud():
    """A cloud reaching two sectors is the same cloud from both."""
    for step in range(400):
        center = _on_ring(ARM_ANGLE + step * 0.002)
        clouds = nebula_field.clouds_reaching(SEED, SHAPE, center, REACH_PC)
        big = [cloud for cloud in clouds if ly_to_pc(cloud[0].radius_ly) > 3 * EDGE_PC]
        if big:
            break
    else:
        pytest.fail("no big cloud on 400 arm sectors")
    # One sector's width toward the cloud's center: still inside its reach.
    cloud_center = big[0][1]
    gap = math.dist(center, cloud_center)
    neighbor = tuple(c + (t - c) * min(1.0, EDGE_PC / gap) for c, t in zip(center, cloud_center))
    keys = {_key(cloud) for cloud in nebula_field.clouds_reaching(SEED, SHAPE, neighbor, REACH_PC)}
    assert _key(big[0]) in keys


def test_a_neighborhood_finds_clouds_at_realistic_rates():
    """Arm sectors at the solar circle sit inside a dark cloud far more
    often than sectors between the arms, and a few hundred arm sectors
    find some (the old per-sector roll gave each sector about 3e-4)."""
    def share_inside(angle):
        hits = sum(bool(nebula_field.clouds_reaching(SEED, SHAPE, _on_ring(angle + step * 0.002), REACH_PC))
                   for step in range(300))
        return hits / 300
    arm, between = share_inside(ARM_ANGLE), share_inside(ARM_ANGLE + math.pi / 2)
    assert 0.03 < arm < 0.5
    assert between < arm


def _generate_at(mysql_config, address):
    args = run_galaxy._default_generation_args(mysql_config)
    args.num_systems = 1
    args.density = None
    position = sector_position_pc(*address, EDGE_PC)
    with forced_system_config(PLANETS=False):
        return run_galaxy.generate_and_save_sector_at(args, address, position, EDGE_PC)


def test_a_cloud_reaching_several_sectors_is_stored_once(mysql_config, monkeypatch):
    """Generated sectors store each field cloud reaching them once, at its
    own galaxy-frame center, and every one of those sectors' systems
    inside it says so."""
    _db.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=3000,
                          expected_system_count_at_density_1=1.0, config=mysql_config, galaxy_seed=SEED)
    _db.replace_galaxy_layers([(layer, 3000) for layer in range(-2, 3)], config=mysql_config)
    monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, "molecular-cloud", 20.0)

    # Three sectors in a row on the arm's crest at the solar circle.
    ring, layer, slot = sector_address_at(_on_ring(ARM_ANGLE), EDGE_PC)
    addresses = [(ring, layer, slot + step) for step in range(3)]
    for address in addresses:
        _generate_at(mysql_config, address)

    expected = {}
    for address in addresses:
        center = sector_position_pc(*address, EDGE_PC)
        for nebula, cloud_center in nebula_field.clouds_reaching(SEED, SHAPE, center, REACH_PC):
            expected[cloud_center] = nebula
    assert expected, "a raised rate must put a cloud over three arm sectors"

    conn = _db.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT id, name, nebula_class, radius_ly, center_x_pc, center_y_pc, center_z_pc"
                            " FROM nebulae WHERE nebula_type = 'dark'").fetchall()
        stored = {(row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]): row for row in rows}
        assert len(stored) == len(rows) == len(expected)
        for center, nebula in expected.items():
            row = stored[center]
            assert row["nebula_class"] == nebula.nebula_class
            assert row["radius_ly"] == pytest.approx(nebula.radius_ly)
        whole = [address for address in addresses
                 if any(math.dist(center, sector_position_pc(*address, EDGE_PC)) + REACH_PC
                        < ly_to_pc(nebula.radius_ly) for center, nebula in expected.items())]
        for address in whole:
            sector_id = _db.get_sector_id_at(conn, *address)
            outside = conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?"
                                   " AND inside_nebula_id IS NULL", (sector_id,)).fetchone()["n"]
            assert outside == 0
    finally:
        conn.close()
