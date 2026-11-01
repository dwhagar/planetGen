# planetgen/generation/clusters.py

"""
A generated galaxy's globular clusters (GEN.159; docs/design/galaxy/globular-clusters.md sections 8 and 9, slice 4).

`generate_clusters` draws the whole cluster system from the galaxy seed alone: each cluster's mass (lognormal about
2e5 solar masses), half-mass radius and King concentration, whether its core has collapsed, its [Fe/H] (two groups, the
metal-poor halo clusters and the metal-rich ones of the bulge and disk), age, and its place and velocity in the galaxy
frame. It is the generator the Milky Way catalogue will replace row by row once Boss supplies one; until then these
rows stand in for the real 157. How many it draws follows the galaxy's halo mass (`cluster_count`, GEN.403). The constants are in `tuning` ("a generated galaxy's globular cluster system").
"""

import math

from planetgen import tuning
from planetgen.galaxy import components, enclosed_mass, potential, shape_draw
from planetgen.physics import king
from planetgen.util import draw
from planetgen.util import detmath

_GROUP_HALO = 0
_ATTEMPTS = 64


def _stream(galaxy_seed, part=None):
    """The stream every draw of the cluster system comes from: the seed, bytes or an int, under the system's own name (and
    `part`'s, for the count, which must not move the clusters themselves)."""
    text = galaxy_seed.hex() if isinstance(galaxy_seed, (bytes, bytearray)) else str(galaxy_seed)
    return draw.Stream(f"globular-clusters:{text}" if part is None else f"globular-clusters-{part}:{text}")


def _clipped_gauss(stream, mean, sigma, low, high):
    return min(max(stream.gauss(mean, sigma), low), high)


def _mass(stream):
    low, high = tuning.CLUSTER_MASS_RANGE_SOL
    for _ in range(_ATTEMPTS):
        mass = tuning.CLUSTER_MASS_PEAK_SOL * 10.0 ** stream.gauss(0.0, tuning.CLUSTER_MASS_SIGMA_DEX)
        if low <= mass <= high:
            return mass
    return min(max(mass, low), high)


def _direction(stream):
    """A unit vector, uniform on the sphere."""
    z = stream.uniform(-1.0, 1.0)
    phi = stream.uniform(0.0, 2.0 * math.pi)
    rho = math.sqrt(1.0 - z * z)
    return rho * math.cos(phi), rho * math.sin(phi), z


def _halo_distance(stream):
    """A distance from the centre for a density falling as r**-3.5 between the halo's limits: the number inside r grows
    as r**-0.5 away from the inner limit, which inverts exactly."""
    low, high = tuning.CLUSTER_HALO_RADIUS_RANGE_PC
    shrink = 1.0 - (low / high) ** 0.5
    return low * (1.0 - stream.random() * shrink) ** -2.0


def _bulge_distance(stream):
    low, _high = tuning.CLUSTER_HALO_RADIUS_RANGE_PC
    return max(tuning.CLUSTER_BULGE_MEDIAN_PC * 10.0 ** stream.gauss(0.0, 0.3), low / 3.0)


def _place(stream, group, scale, pot):
    """`(position, velocity)` in the galaxy frame (the position before `scale` is applied): a halo cluster on an isotropic orbit, a metal-rich one
    flattened towards the plane and rotating with the disk at the circular speed of the galaxy's own potential `pot` at its (scaled) radius (GEN.339)."""
    if group == _GROUP_HALO:
        dx, dy, dz = _direction(stream)
        distance = _halo_distance(stream)
        dispersion = tuning.CLUSTER_HALO_DISPERSION_KMS
        rotation = 0.0
    else:
        dx, dy, dz = _direction(stream)
        distance = _bulge_distance(stream)
        dz *= tuning.CLUSTER_BULGE_FLATTENING
        dispersion = tuning.CLUSTER_BULGE_DISPERSION_KMS
        rotation = None
    x, y, z = distance * dx, distance * dy, distance * dz
    cylindrical = detmath.hypot(x, y) or 1.0
    if rotation is None:
        rotation = potential.circular_speed_kms(pot, cylindrical * scale)
    velocity = (
        -rotation * y / cylindrical + stream.gauss(0.0, dispersion),
        rotation * x / cylindrical + stream.gauss(0.0, dispersion),
        stream.gauss(0.0, dispersion),
    )
    return (x, y, z), velocity


def _structure(stream, galactocentric_pc, scale=1.0):
    """`(half_mass_radius_pc, king_w0, tidal_radius_pc, core_collapsed)` of one cluster at this distance from the
    centre: the King concentration is lowered where the Jacobi radius would cut the cluster's edge off."""
    low, high = tuning.CLUSTER_HALF_MASS_RADIUS_RANGE_PC
    reference_pc = tuning.CLUSTER_RADIUS_REFERENCE_SCALE_LENGTHS * tuning.HALO_SCALE_LENGTH_PC * scale
    growth = (max(galactocentric_pc, 1000.0 * scale) / reference_pc) ** tuning.CLUSTER_RADIUS_GROWTH
    radius = tuning.CLUSTER_HALF_MASS_RADIUS_PC * 10.0 ** stream.gauss(
        0.0, tuning.CLUSTER_HALF_MASS_RADIUS_SIGMA_DEX) * growth
    radius = min(max(radius, low), high)
    collapsed = stream.random() < tuning.CLUSTER_CORE_COLLAPSED_FRACTION
    if collapsed:
        w0 = tuning.KING_W0_RANGE[1]
    else:
        w0 = king.nearest_w0(stream.gauss(tuning.CLUSTER_CONCENTRATION_MEAN, tuning.CLUSTER_CONCENTRATION_SIGMA))
    return radius, w0, collapsed


def _fe_h(stream, group):
    """The cluster's [Fe/H]: the halo group is the first of `CLUSTER_FE_H_GROUPS`, the bulge group the second."""
    _fraction, mean, sigma = tuning.CLUSTER_FE_H_GROUPS[0 if group == _GROUP_HALO else 1]
    low, high = tuning.METALLICITY_RANGE_FE_H
    return _clipped_gauss(stream, mean, sigma, low, high)


def mean_cluster_mass_solar():
    """The mean mass of one cluster, solar masses: the mean of the lognormal mass function (peak `CLUSTER_MASS_PEAK_SOL`,
    width `CLUSTER_MASS_SIGMA_DEX`) cut to `CLUSTER_MASS_RANGE_SOL`, from the truncated lognormal's closed form."""
    sigma = tuning.CLUSTER_MASS_SIGMA_DEX * math.log(10.0)
    low, high = (math.log(m / tuning.CLUSTER_MASS_PEAK_SOL) / sigma for m in tuning.CLUSTER_MASS_RANGE_SOL)

    def cdf(z):
        return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
    kept = cdf(high) - cdf(low)
    return tuning.CLUSTER_MASS_PEAK_SOL * math.exp(0.5 * sigma * sigma) * (cdf(high - sigma) - cdf(low - sigma)) / kept


def cluster_count(galaxy_seed, halo_mass_solar, family_factor=1.0):
    """
    How many globular clusters the galaxy holds (GEN.403): the system's mass over the mean cluster mass, the system's mass
    being `eta` times the halo's with a log-normal scatter (`tuning.CLUSTER_SYSTEM_FRACTION_RANGE`,
    `CLUSTER_SYSTEM_SCATTER_DEX`), rounded up or down at random so the mean count is exact for a halo that holds well under
    one cluster.

    Args:
        galaxy_seed (bytes or int): The galaxy's seed; `eta`, the scatter and the rounding come from it.
        halo_mass_solar (float or None): The galaxy's halo mass; `None` or `0` is no halo and no clusters.
        family_factor (float): A galaxy family's specific-frequency factor (1 for the first family; GEN.387 sets the
            others), multiplying the system's mass.

    Returns:
        int: The cluster count, `0` or more.
    """
    if not halo_mass_solar or halo_mass_solar <= 0.0:
        return 0
    stream = _stream(galaxy_seed, "count")
    low, high = (math.log10(v) for v in tuning.CLUSTER_SYSTEM_FRACTION_RANGE)
    eta = 10.0 ** stream.uniform(low, high)
    system = eta * halo_mass_solar * 10.0 ** stream.gauss(0.0, tuning.CLUSTER_SYSTEM_SCATTER_DEX) * family_factor
    expected = system / mean_cluster_mass_solar()
    whole = math.floor(expected)
    return whole + (1 if stream.random() < expected - whole else 0)


def generate_clusters(galaxy_seed, count, epoch_unix=0.0, scale=1.0, shape=None):
    """
    Draws a galaxy's globular clusters from its seed.

    Args:
        galaxy_seed (bytes or int): The galaxy's seed; the same seed gives the same clusters.
        count (int): How many clusters (`cluster_count` from the galaxy's halo mass).
        epoch_unix (float): When the positions and velocities hold.
        scale (float): How much larger than the Milky Way the galaxy is, by its disk scale length: every distance from
            the centre is multiplied by it, so a toy galaxy keeps its clusters inside itself. The clusters themselves
            (mass, size, concentration) do not scale.
        shape (GalaxyShape, optional): The galaxy's shape, whose potential and enclosed mass set a metal-rich cluster's rotation speed and every
            cluster's Jacobi radius (GEN.339); `None` takes the shape the seed draws with every lever on Auto (`shape_draw.default_shape`).

    Returns:
        list[Component]: `count` clusters named "Generated cluster 001" and on, in drawing order.
    """
    stream = _stream(galaxy_seed)
    seed_bytes = potential.seed_bytes(galaxy_seed)
    if shape is None:
        shape = shape_draw.default_shape(seed_bytes)
    pot = potential.for_galaxy(shape, seed_bytes)
    halo_fraction = tuning.CLUSTER_FE_H_GROUPS[0][0]
    clusters = []
    for index in range(count):
        group = _GROUP_HALO if stream.random() < halo_fraction else 1
        position, velocity = _place(stream, group, scale, pot)
        position = tuple(c * scale for c in position)
        distance = math.sqrt(sum(c * c for c in position))
        mass = _mass(stream)
        radius, w0, collapsed = _structure(stream, distance, scale)
        w0 = king.limit_w0(w0, radius, king.jacobi_radius_pc(mass, distance, enclosed_mass.enclosed_mass_solar(shape, seed_bytes, max(distance, 1.0))))
        clusters.append(components.Component(
            kind="cluster", name=f"Generated cluster {index + 1:03d}", source="generated",
            center_pc=position, velocity_kms=velocity, epoch_unix=epoch_unix, mass_solar=mass,
            half_mass_radius_pc=radius, profile="king", age_gy=stream.uniform(*tuning.CLUSTER_AGE_RANGE_GY),
            metallicity_fe_h=_fe_h(stream, group), seed=stream.getrandbits(63), king_w0=w0,
            tidal_radius_pc=king.tidal_radius_pc(w0, radius), core_collapsed=collapsed and w0 == tuning.KING_W0_RANGE[1]))
    return clusters
