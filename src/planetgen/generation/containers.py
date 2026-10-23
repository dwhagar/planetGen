# planetgen/generation/containers.py

"""
The galaxy's containers (DB.85): each globular cluster and each superbubble as a row of `galaxy_containers`.

A container is a region of space that holds objects: a globular cluster inside its tidal radius, a superbubble inside
its cavity. Both come from the galaxy seed alone, so nothing here is fixed content: the clusters are the plan's
`galaxy_components` (`generation.clusters`), the bubbles are the superbubble field (`galaxy.superbubble_field`, GEN.254),
and `build` lists them with the object ID each is named by.

The object ID is the sector holding the container's centre (its birth sector), the `SERIAL_REGION` kind and the
container's rank among the containers born in that sector, counted in a fixed order (clusters by name and seed, then
bubbles by column and draw), so the same seed gives the same IDs on every run.
"""

import dataclasses

from planetgen.galaxy import object_uid as objectUid, spherical as galaxySpherical, superbubble_field

GLOBULAR_CLUSTER = "globular_cluster"
SUPERBUBBLE = "superbubble"
KINDS = (GLOBULAR_CLUSTER, SUPERBUBBLE)
"""tuple: The kinds of container (the table's CHECK), in the order their IDs are ranked."""

GLOBULAR_PROFILE = "king"
"""str: The density profile of the `cluster` components that are globular clusters (the young open clusters are not)."""


@dataclasses.dataclass(frozen=True)
class Container:
    """
    One row of `galaxy_containers`, before it has an `id`.

    Attributes:
        kind (str): One of `KINDS`.
        name (str or None): A cluster's name; `None` for a bubble.
        component_id (int or None): The `galaxy_components` row of a cluster.
        bubble (tuple or None): `(i, j, n)`, the column and draw of a bubble in the superbubble field.
        center_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        extent_pc (float): A cluster's tidal radius or a bubble's radius, parsecs.
        birth (tuple): `(shell, band, wedge)` of the sector holding the centre.
        uid (bytes): The container's object ID.
    """

    kind: str
    name: str | None
    component_id: int | None
    bubble: tuple | None
    center_pc: tuple
    extent_pc: float
    birth: tuple
    uid: bytes


def _is_globular(component):
    return component.kind == "cluster" and component.profile == GLOBULAR_PROFILE


def _wanted(galaxy_seed, shape, reach_pc, components):
    """`(sort key, kind, name, component id, bubble address, centre, extent)` of every container, unsorted."""
    wanted = []
    for component in components:
        if _is_globular(component):
            wanted.append(((0, component.name or "", component.seed), GLOBULAR_CLUSTER, component.name, component.id,
                           None, tuple(component.center_pc), component.tidal_radius_pc))
    if galaxy_seed is not None and shape is not None and reach_pc > 0.0:
        for bubble in superbubble_field.all_bubbles(galaxy_seed, shape, reach_pc):
            wanted.append(((1, bubble.index[0], bubble.index[1], bubble.n), SUPERBUBBLE, None, None,
                           (bubble.index[0], bubble.index[1], bubble.n), tuple(bubble.center_pc), bubble.radius_pc))
    return wanted


def build(galaxy_seed, shape, edge_pc, reach_pc, components):
    """
    Every globular cluster and superbubble of a galaxy as a container, in ID order.

    Args:
        galaxy_seed (bytes or int or None): The galaxy's seed, which the superbubble field draws from; no bubbles
            when `None`.
        shape (GalaxyShape): The galaxy's shape, for the bubble field's gas tracer.
        edge_pc (float): The sector edge, to find the sector holding each centre.
        reach_pc (float): How far from the axis a bubble's centre may lie: the charted galaxy's reach.
        components (iterable[galaxy.components.Component]): The stored components, each with its `id`; the globular
            clusters (a `cluster` with the `king` profile) become containers.

    Returns:
        list[Container]: Clusters first (by name, then seed), then bubbles (by column, then draw).
    """
    wanted = sorted(_wanted(galaxy_seed, shape, reach_pc, components), key=lambda item: item[0])
    ranks = {}
    out = []
    for _key, kind, name, component_id, bubble, center, extent in wanted:
        birth = galaxySpherical.sector_address_at(center, edge_pc)
        rank = ranks.get(birth, 0)
        ranks[birth] = rank + 1
        uid = objectUid.to_bytes(objectUid.pack(*birth, objectUid.SERIAL_REGION, rank))
        out.append(Container(kind, name, component_id, bubble, center, float(extent), birth, uid))
    return out
