# planetgen/galaxy/containment.py

"""
Which cloud holds a point (GEN.357)
===================================

The geometry the store (`store.refresh_containment`, after a sector is saved) and sector generation
(`generation.cloud`, while a system is placed) share, so a system is inside the same cloud whichever
tests it. A container is a dict with `column` (`"inside_nebula_id"` or `"inside_remnant_id"`), `id`,
`center` and `radius_pc` in the galaxy frame, and for a nebula `inside`, a function of a point that says whether the nebula's
contour holds it (`contour_test`, GEN.1267). Pure geometry and the galaxy's density: no database, no I/O.
"""

from planetgen.util import detmath


def innermost_container(point_pc, containers, own_radius_pc=0.0, own=None):
    """
    The smallest container in `containers` whose sphere holds `point_pc`
    -- for a nebula (`own_radius_pc` > 0), only a container larger than it,
    and never itself (`own`, a `(column, id)` pair). A nebula holds a point
    only inside its contour, the density contour the map draws
    (`nebula_density.inside`, GEN.1267), which its container's `inside`
    function supplies; its sphere is the quick first test. A container
    without `inside` (a remnant) is a plain sphere.

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
        inside = container.get("inside")
        if inside is not None and not inside(point_pc):
            continue
        best = container
    return best


def contour_test(galaxy, nebula_class, center_pc, radius_pc):
    """
    The `inside` function of a nebula's container (GEN.1267): given a point (galaxy-frame parsecs), whether the nebula's contour holds it. The
    nebula is read as its class's family, its center and its radius, against the galaxy's seed and shape, so the generator and the stored-row
    refresh agree and a map draws the same outline (`nebula_density.inside`). The density is read only when asked, and a call costs a few
    milliseconds, so the caller tests the sphere first (`innermost_container` does).

    Args:
        galaxy (callable): `galaxy()` gives `(galaxy_seed, galaxy_shape)`, called on the first point that reaches the density (a skeleton can
            take a minute to build and most containers never see a point inside their sphere); or `None` for a database that has never been
            planned, which has no density, so the nebula is its sphere.
        nebula_class (str): The nebula's class letter, which says whether the gas or the young stars anchor it.
        center_pc (tuple): The nebula's center, parsecs.
        radius_pc (float): The nebula's radius, parsecs.
    """
    from planetgen.galaxy import nebula_density

    field = nebula_density.field_of_class(nebula_class)
    center = tuple(center_pc)

    def inside(point_pc):
        found = galaxy()
        if found is None:
            return True
        return bool(nebula_density.inside(found[0], found[1], point_pc, center, radius_pc, field))

    return inside
