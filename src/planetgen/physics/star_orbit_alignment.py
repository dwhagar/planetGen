"""How well the orbit of a star's companion lines up with the plane it is measured against (GEN.1500).

A pair of stars forms with a random relative orientation when it is wide (a core that fragments carries no memory of its
parent's spin), and is pulled toward the plane of the pair it orbits when it is close (the disks and the accretion torques
that fed it last long enough to line the angular momenta up). The observed alignment of the disks and orbits of young
multiples, of triples and of planets with a wide companion is complete inside a few tens of AU and gone by about a thousand
([R3333], recalled from the review, not read). So the angle a companion's orbit makes with its reference plane (the plane of the
pair it orbits, or for the first companion the plane of the star's planets) is drawn here from two populations: aligned, a
small prograde tilt (a Rayleigh law); and random, the isotropic sphere the draw used for every separation before. The share
of the aligned one falls with the separation, log-linearly from all of them to none.

The module holds the law and the geometry (`compose_plane`), and `generation/multistar.py` is the caller.
"""

import math

ALIGNED_FULL_AU = 30.0
"""float: Separation (AU) at and inside which every companion orbit is drawn aligned (stand-in AL1)."""

ALIGNED_NONE_AU = 1000.0
"""float: Separation (AU) at and beyond which no companion orbit is drawn aligned (stand-in AL1)."""

ALIGNED_SCALE_DEG = 5.0
"""float: The scale (degrees) of the Rayleigh law of the tilt of an aligned orbit (stand-in AL2)."""


def aligned_share(separation_au):
    """The share of companion orbits at `separation_au` that are drawn aligned: 1 at `ALIGNED_FULL_AU` and inside, 0 at `ALIGNED_NONE_AU` and beyond, linear in the logarithm of the separation between."""
    if separation_au <= ALIGNED_FULL_AU:
        return 1.0
    if separation_au >= ALIGNED_NONE_AU:
        return 0.0
    return math.log(ALIGNED_NONE_AU / separation_au) / math.log(ALIGNED_NONE_AU / ALIGNED_FULL_AU)


def draw_mutual_inclination(rng, separation_au):
    """The angle (radians, 0 to pi) between a companion's orbit and its reference plane, drawn from the law at `separation_au`. Draws three numbers from `rng`, always, so a caller's draw order does not depend on which population came up."""
    pick, size, _unused = rng.random(), rng.random(), rng.random()
    if pick < aligned_share(separation_au):
        return math.radians(ALIGNED_SCALE_DEG) * math.sqrt(-2.0 * math.log(1.0 - size))
    return math.acos(1.0 - 2.0 * size)


def compose_plane(reference_inclination, reference_node, mutual, azimuth):
    """`(inclination, node)` in radians of the plane that makes the angle `mutual` with the reference plane `(reference_inclination, reference_node)`, turned by `azimuth` about the reference plane's normal.

    The normal of a plane of inclination `i` and ascending node `n` is `(sin i sin n, -sin i cos n, cos i)`, the form in which the angle between two planes is `stability.mutual_inclination`.
    """
    si, ci = math.sin(reference_inclination), math.cos(reference_inclination)
    sn, cn = math.sin(reference_node), math.cos(reference_node)
    normal = (si * sn, -si * cn, ci)
    along = (cn, sn, 0.0)                                             # on the reference plane, toward its ascending node
    across = (-ci * sn, ci * cn, si)                                 # normal x along: the other axis of the reference plane
    sm, cm = math.sin(mutual), math.cos(mutual)
    sa, ca = math.sin(azimuth), math.cos(azimuth)
    x, y, z = (cm * normal[k] + sm * (ca * along[k] + sa * across[k]) for k in range(3))
    inclination = math.acos(max(-1.0, min(1.0, z)))
    node = math.atan2(x, -y) % (2.0 * math.pi)
    return inclination, node
