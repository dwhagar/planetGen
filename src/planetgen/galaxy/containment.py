# planetgen/galaxy/containment.py

"""
Which cloud holds a point (GEN.357)
===================================

The geometry the store (`store.refresh_containment`, after a sector is saved) and sector generation
(`generation.cloud`, while a system is placed) share, so a system is inside the same cloud whichever
tests it. A container is a dict with `column` (`"inside_nebula_id"` or `"inside_remnant_id"`), `id`,
`center` and `radius_pc` in the galaxy frame, and for a nebula `shape`, a function returning its
`NebulaShape` (GEN.75). Pure geometry: no database, no I/O.
"""

from planetgen.util import detmath


def innermost_container(point_pc, containers, own_radius_pc=0.0, own=None):
    """
    The smallest container in `containers` whose sphere holds `point_pc`
    -- for a nebula (`own_radius_pc` > 0), only a container larger than it,
    and never itself (`own`, a `(column, id)` pair). A nebula holds a point
    only inside its shape (GEN.75, `planetgen.galaxy.nebula_shape`), which
    its container's `shape` function supplies; its sphere is the quick
    first test. A container without `shape` (a remnant) is a plain sphere.

    Returns:
        dict or None: The container, or `None` when the point is in open space.
    """
    best = None
    for container in containers:
        if own is not None and (container["column"], container["id"]) == own:
            continue
        if container["radius_pc"] <= own_radius_pc:
            continue
        if detmath.dist(point_pc, container["center"]) > container["radius_pc"]:
            continue
        if best is not None and container["radius_pc"] >= best["radius_pc"]:
            continue
        shape = container.get("shape")
        if shape is not None:
            radius = container["radius_pc"]
            offset = tuple((point_pc[i] - container["center"][i]) / radius for i in range(3))
            if not shape().contains(offset):
                continue
        best = container
    return best
