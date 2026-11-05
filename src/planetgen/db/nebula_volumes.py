# planetgen/db/nebula_volumes.py

"""
The points that draw the nebulae of one Galaxy Map node (MAP.365).

`group_volume(conn, key)` reads the nebulae whose centre lies in the block `key` (`m.shell.band.wedge`, the node key of the map's
ladder) and the nebulae around it, and has `galaxy.nebula_volume` sample their shape from the model's density
(`galaxy.nebula_density`). The same key always gives the same bytes for the same database contents; nothing is stored.
"""

from planetgen import tuning
from planetgen.db import query
from planetgen.physics import units
from planetgen.galaxy import nebula_blocks, nebula_volume
from planetgen.galaxy.drill import parse_node_key
from planetgen.galaxy.nodecut import block_box

MAX_NEBULAE = 600
"""int: Most nebulae of a block drawn (the biggest first); the rest are too small to show at the lattice the others need."""


def family_of(cloud):
    """The family (`nebula_blocks.FAMILIES`) of a cloud row from `query.galaxy_clouds_in_box`, from its class letter."""
    letter = cloud.get("class")
    if letter in tuning.NEBULA_CLASSES:
        return tuning.NEBULA_CLASSES[letter]["family"]
    return "supernova-remnant" if cloud.get("type") == "supernova_remnant" else "diffuse"


def inside_block(block, edge_pc):
    """`inside(points_pc) -> bool array` for `sample_group`: which (n, 3) points lie in `block`, by the cut the ladder draws
    (`nebula_blocks.locate`)."""
    def inside(points):
        shell, band, wedge = nebula_blocks.locate(points[:, 0], points[:, 1], points[:, 2], block.m, edge_pc)
        return (shell == block.shell) & (band == block.band) & (wedge == block.wedge)
    return inside


def block_nebulae(conn, block, edge_pc):
    """The nebulae whose sphere reaches `block`, each `{"key", "x", "y", "z", "radius_pc", "family"}`, the biggest first."""
    lo, hi = block_box(block, edge_pc)
    clouds = query.galaxy_clouds_in_box(conn, lo, hi, max_clouds=MAX_NEBULAE)
    return [{"key": f"{cloud['type']}:{cloud['id']}", "x": cloud["x"], "y": cloud["y"], "z": cloud["z"],
             "radius_pc": cloud["radius_pc"], "family": family_of(cloud)} for cloud in clouds]


def group_volume(conn, key, model=None):
    """
    The packed points (`nebula_volume.pack`) that draw the nebulae in block `key`: every lattice point in the block that a nebula
    holds, with its strength from the model's density.

    Args:
        conn: An open, read-only connection.
        key (str): A node key, `m.shell.band.wedge`.
        model (module, optional): Where the density is read: `density_at`, `FORM_LEVEL` and `MAX_DETAIL` (default `galaxy.nebula_density`).

    Returns:
        bytes: The packed volume; header-only for a block no nebula reaches or a galaxy that has never been planned.

    Raises:
        ValueError: On a malformed key.
    """
    block = parse_node_key(key)
    skeleton = query.get_galaxy_shape(conn)
    nebulae = [] if skeleton is None or skeleton.galaxy_seed is None else block_nebulae(conn, block, skeleton.edge_pc)
    if not nebulae:
        return _empty()
    if model is None:
        from planetgen.galaxy import nebula_density as model
    seed, shape = skeleton.galaxy_seed, skeleton.shape

    def density(points, field, detail):
        return model.density_at(seed, shape, points, field, detail)

    volume = nebula_volume.sample_group(
        nebulae, inside_block(block, skeleton.edge_pc), density, model.FORM_LEVEL, model.MAX_DETAIL, nebula_blocks.FAMILIES)
    return nebula_volume.pack(volume)


def single_nebula_volume(conn, nebula_id, model=None):
    """
    The packed points (`nebula_volume.pack`) that draw one nebula alone (MAP.142, the Sector Map and the nebula page): the same sampler and density as
    `group_volume`, over this nebula's sphere only and at the lattice level its own radius asks for, so a nebula looks the same in a node of the Galaxy
    Map and by itself.

    Args:
        conn: An open, read-only connection.
        nebula_id (int): `nebulae.id`.
        model (module, optional): As in `group_volume`.

    Returns:
        bytes: The packed volume; header-only for a nebula never placed or a galaxy that has never been planned.

    Raises:
        ValueError: If no such nebula exists.
    """
    row = conn.execute(
        "SELECT radius_ly, center_x_pc, center_y_pc, center_z_pc, nebula_class FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
    if row is None:
        raise ValueError(f"no such nebula: {nebula_id!r}")
    skeleton = query.get_galaxy_shape(conn)
    if skeleton is None or skeleton.galaxy_seed is None or row["center_x_pc"] is None or not row["radius_ly"]:
        return _empty()
    if model is None:
        from planetgen.galaxy import nebula_density as model
    seed, shape = skeleton.galaxy_seed, skeleton.shape

    def density(points, field, detail):
        return model.density_at(seed, shape, points, field, detail)

    cloud = {"class": row["nebula_class"], "type": "nebula"}
    nebula = {"key": f"nebula:{nebula_id}", "x": row["center_x_pc"], "y": row["center_y_pc"], "z": row["center_z_pc"],
              "radius_pc": units.ly_to_pc(row["radius_ly"]), "family": family_of(cloud)}
    volume = nebula_volume.sample_group(
        [nebula], None, density, model.FORM_LEVEL, model.MAX_DETAIL, nebula_blocks.FAMILIES)
    return nebula_volume.pack(volume)


def _empty():
    """The header-only volume of a nebula with nothing to draw."""
    return nebula_volume.pack(nebula_volume.sample_group([], None, None, 0.0, 0, nebula_blocks.FAMILIES))
