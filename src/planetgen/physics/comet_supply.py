# planetgen/physics/comet_supply.py

"""
Where a star's own comets come from (GEN.1765, slice 2; docs/design/systems/comet-orbital-realism.md)
=====================================================================================================

Comets are the planetesimals of a star's outer disk that a giant has thrown about. Where the star has such a disk, two of the hand-set shapes of GEN.595 are read from it:

1. **The family shares.** The giant's scatter step (`physics/scatter.scatter_shares`, kept on the belt as `scatter_fates`) gives the shares of the disk in each population. The
   `scattered_disk` and the `centaurs` are the comets that come inward and are held by the giant: the Jupiter family. The `detached` population is the source of the cloud, whose
   comets are the Halley-type and the long-period ones, split between the two as the default weights split them (`tuning.COMET_PERIOD_CLASSES`). The shares of the star's outer disks
   are added, each weighted by the mass its disk had.
2. **The top of the nucleus draw.** No nucleus is larger than the largest body the star's outer disks hold (`largest_body_km`), so the top of
   `tuning.BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM` is cut to it, never under the draw's lower end.

A star with no outer disk that a giant scattered keeps the default weights and the default top. The chance that a star has comets at all and how many it has stay the hand-set
`SYSTEM_COMET_CHANCE` and `SYSTEM_COMET_COUNT_RANGE`: the census of GEN.1765 found no star in 398 with a scattered outer disk and a giant to scatter it, so a rule of "comets only
where a disk was scattered" would leave nearly every star without any.

Use:

    weights = family_weights(belts)             # {family: weight} or None
    top_km = nucleus_top_km(belts)              # km or None
"""

from planetgen import tuning

OUTER_DISK = "outer-disk"
"""str: The `origin` of a belt made from a star's outer disk (`physics/belt_origin`)."""

_INWARD = ("scattered_disk", "centaurs")
_CLOUD = ("detached",)


def _outer_disks(belts):
    return [belt for belt in belts if getattr(belt, "origin", None) == OUTER_DISK]


def family_weights(belts):
    """
    The weights of the three comet families for a star whose outer disks a giant scattered.

    Args:
        belts (iterable): The star's belts (`AsteroidBelt`); each outer-disk one with `scatter_fates` and a `primordial_mass_earth` counts.

    Returns:
        dict or None: `{family: weight}` over the keys of `tuning.COMET_PERIOD_CLASSES` when the disks put mass in the inward or the cloud populations; `None` when no outer disk was scattered (the
            caller keeps the default weights).
    """
    inward = cloud = 0.0
    for belt in _outer_disks(belts):
        fates = getattr(belt, "scatter_fates", None)
        mass = getattr(belt, "primordial_mass_earth", None)
        if not fates or not mass or not mass > 0.0:
            continue
        inward += mass * sum(fates.get(name, 0.0) for name in _INWARD)
        cloud += mass * sum(fates.get(name, 0.0) for name in _CLOUD)
    if not inward + cloud > 0.0:
        return None
    classes = tuning.COMET_PERIOD_CLASSES
    halley, long_period = classes["halley_type"]["weight"], classes["long_period"]["weight"]
    halley_share = halley / (halley + long_period)
    return {"jupiter_family": inward, "halley_type": cloud * halley_share, "long_period": cloud * (1.0 - halley_share)}


def nucleus_top_km(belts):
    """
    The largest nucleus the star's outer disks can give, km: their largest body, or `None` when the star has no outer disk (the draw keeps its own top).
    """
    sizes = [belt.largest_body_km for belt in _outer_disks(belts) if getattr(belt, "largest_body_km", None)]
    return max(sizes) if sizes else None
