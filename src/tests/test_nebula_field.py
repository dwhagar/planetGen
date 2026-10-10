"""
The galaxy's molecular cloud field (GEN.47, `nebulaField`): dark-family
clouds drawn per field cell from the galaxy seed, so every sector and
every run agrees where they are; more in the arms and near the plane; a
generated neighborhood finds them, and each is stored once however many
sectors it reaches.
"""

import math

import pytest

from planetgen.util import draw
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
SOLAR_RADIUS_PC = tuning.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH * SHAPE.disk_scale_length_pc
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
    draw.set_run_seed(7)
    expected = [draw.random() for _ in range(3)]
    draw.set_run_seed(7)
    nebula_field.clouds_reaching(SEED, SHAPE, _on_ring(ARM_ANGLE), REACH_PC)
    assert [draw.random() for _ in range(3)] == expected


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
        # GEN.75: a cloud holds a system only inside its shape, not its sphere.
        from planetgen.db import query as queryDb

        shapes = [(center, ly_to_pc(nebula.radius_ly), queryDb.nebula_shape(conn, stored[center]["id"])[1])
                  for center, nebula in expected.items()]
        for address in addresses:
            sector_id = _db.get_sector_id_at(conn, *address)
            sector_center = sector_position_pc(*address, EDGE_PC)
            for system in conn.execute(
                    "SELECT position_x_mpc, position_y_mpc, position_z_mpc, inside_nebula_id FROM star_systems"
                    " WHERE sector_id = ?", (sector_id,)).fetchall():
                point = _db.local_to_galaxy_pc(sector_center, tuple(
                    system[f"position_{axis}_mpc"] / 1000.0 for axis in "xyz"))
                inside = any(shape.contains(tuple((point[i] - center[i]) / radius for i in range(3)))
                             for center, radius, shape in shapes)
                assert (system["inside_nebula_id"] is not None) == inside
    finally:
        conn.close()


def _some_clouds(count=12):
    found = []
    for step in range(600):
        for cloud in nebula_field.clouds_reaching(SEED, SHAPE, _on_ring(ARM_ANGLE + step * 0.004), REACH_PC):
            if all(cloud[1] != other[1] for other in found):
                found.append(cloud)
        if len(found) >= count:
            break
    return found


def test_the_centroid_lies_inside_the_clouds_bounding_sphere():
    for nebula, _center in _some_clouds(6):
        centroid = nebula.get_shape().interior_centroid()
        assert math.dist(centroid, (0.0, 0.0, 0.0)) <= 1.0
        assert centroid == nebula.get_shape().interior_centroid()


def test_a_cloud_is_born_in_the_sector_holding_its_centroid(monkeypatch):
    """GEN.176: the ID names the sector holding the centre of the space the cloud fills, rounded to 1 mpc."""
    from planetgen.galaxy import object_uid

    clouds = _some_clouds(8)
    assert clouds
    for nebula, center in clouds:
        value = object_uid.unpack(nebula_field.cloud_object_id(SEED, SHAPE, nebula, center, EDGE_PC))
        mpc = nebula_field.centroid_mpc(nebula, center)
        assert value.serial_kind == object_uid.SERIAL_FIELD and value.body == 0
        assert (value.ring, value.layer, value.slot) == sector_address_at(tuple(v / 1000.0 for v in mpc), EDGE_PC)


def test_field_cloud_ids_are_unique_and_do_not_depend_on_the_order_asked():
    clouds = _some_clouds(12)
    first = [nebula_field.cloud_object_id(SEED, SHAPE, nebula, center, EDGE_PC) for nebula, center in clouds]
    cache = {}
    backwards = [nebula_field.cloud_object_id(SEED, SHAPE, nebula, center, EDGE_PC, cache)
                 for nebula, center in reversed(clouds)]
    assert first == list(reversed(backwards))
    assert len(set(first)) == len(first)


def test_two_clouds_born_in_one_sector_take_ranks_in_cell_order(monkeypatch):
    """The rank counts the clouds born in the sector that come first by field cell and draw order."""
    from planetgen.galaxy import object_uid

    monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, "molecular-cloud", 30.0)   # crowd the cells
    clouds = _some_clouds(16)
    by_sector = {}
    for nebula, center in clouds:
        value = object_uid.unpack(nebula_field.cloud_object_id(SEED, SHAPE, nebula, center, EDGE_PC))
        by_sector.setdefault((value.ring, value.layer, value.slot), []).append(value.serial)
    assert all(len(set(serials)) == len(serials) for serials in by_sector.values())


def test_a_stored_cloud_has_the_same_id_whichever_sector_is_saved_first(mysql_config, monkeypatch):
    """GEN.176: the ID of a cloud reached by several sectors is worked out from the seed, so saving the sectors
    the other way round gives every cloud the ID it had."""
    from planetgen.galaxy import object_uid
    from tests.db_schema_support import scratch_database

    monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, "molecular-cloud", 20.0)
    ring, layer, slot = sector_address_at(_on_ring(ARM_ANGLE), EDGE_PC)
    addresses = [(ring, layer, slot + step) for step in range(3)]

    def ids_after(config, order):
        _db.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=3000,
                              expected_system_count_at_density_1=1.0, config=config, galaxy_seed=SEED)
        _db.replace_galaxy_layers([(lay, 3000) for lay in range(-2, 3)], config=config)
        for address in order:
            _generate_at(config, address)
        conn = _db.get_connection(config)
        try:
            return {(row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]): bytes(row["uid"])
                    for row in conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc, uid FROM nebulae"
                                            " WHERE nebula_type = 'dark'").fetchall()}
        finally:
            conn.close()

    forward = ids_after(mysql_config, addresses)
    assert forward, "a raised rate must put a cloud over three arm sectors"
    with scratch_database(mysql_config) as other:
        backward = ids_after(other, list(reversed(addresses)))
    assert forward == backward
    kinds = {object_uid.unpack(object_uid.from_bytes(uid)).serial_kind for uid in forward.values()}
    assert kinds == {object_uid.SERIAL_FIELD}


def test_a_centroid_is_rounded_to_a_milliparsec_before_its_sector_is_taken():
    """GEN.176's face case: noise below a milliparsec cannot move a cloud across a sector face."""
    from types import SimpleNamespace

    face = ly_to_pc(0.0) + EDGE_PC * 100.0       # a point on a layer face, parsecs
    def cloud_with(offset):
        shape = SimpleNamespace(interior_centroid=lambda: (0.0, 0.0, offset))
        return SimpleNamespace(get_shape=lambda: shape, radius_ly=1.0)

    centre = (SOLAR_RADIUS_PC, 0.0, face)
    quiet = cloud_with(0.0)
    noisy = cloud_with(1e-9)                      # a billionth of the radius: far below 1 mpc
    assert nebula_field.centroid_mpc(quiet, centre) == nebula_field.centroid_mpc(noisy, centre)
    assert nebula_field.birth_sector(quiet, centre, EDGE_PC) == nebula_field.birth_sector(noisy, centre, EDGE_PC)
