# planetgen/generation/open_clusters.py

"""
A generated galaxy's young open clusters (GEN.244; docs/design/galaxy/intergalactic-and-halo-contents.md 5.14).

`generate_open_clusters` draws them from the galaxy seed and the galaxy's own disk, with no cluster typed in:

1. **Formed.** The disk's star formation history (`galaxy.disk_history`, GEN.469) says how much mass the disk formed in
   each ring and each step of the last Gy. The share of it that goes into bound clusters is the cluster formation
   efficiency of the ring's star formation rate density (`physics.open_cluster_evolution.cluster_formation_efficiency`),
   and the mass function (slope -2) turns that mass into a Poisson count of clusters and gives each its initial mass.
2. **Placed.** Each lies at the radius of its ring, at the azimuth and height stars of its age have there
   (`galaxy.density.population_densities`: the `young` population under 0.1 Gy, the `intermediate` one over it, so the
   arms hold the youngest), on the circular orbit with that population's dispersion, with the disk's [Fe/H] at its radius.
3. **Lost.** Stellar evolution and tidal disruption take mass from it at a rate that follows the mean density of its place
   (`physics.open_cluster_evolution.present_masses`); one that is gone, or lighter than the mass function's low end, is
   not a component. The central ring's clusters are drawn under the most massive cluster its gas can form
   (`physics.cluster_maximum_mass`, GEN.1195) and the passing clouds of the ring take mass from them besides
   (`physics.open_cluster_evolution.cloud_shock_mass_lost_sol`).
4. **Sized.** A half-mass radius from a half-mass density, grown with age and held inside a sixth of its Jacobi radius,
   the profile Plummer's.

The result is a `Component` of kind `cluster` and profile `plummer` (a globular cluster is a `king`), with the age
of the cluster; `galaxy.component_field` puts its stars in the population of its age's quarter-decade
(`tuning.OPEN_CLUSTER_POPULATIONS`). An OB association is the unbound form of the same stars and is already the disk's
young population; the clusters here are the bound part, which adds its density to the cell.
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import components, core_gas, disk_history, enclosed_mass, potential
from planetgen.galaxy import structure as galaxy_structure
from planetgen.galaxy import density as galaxy_density
from planetgen.physics import cluster_maximum_mass, component_profiles, constants, metallicity
from planetgen.physics import open_cluster_evolution as evolution
from planetgen.physics.units import pc_to_ly
from planetgen.util import draw
from planetgen.util.poisson import poisson_count


_ARM_ROUNDS = 24
_MAX_AGE_GY = tuning.OPEN_CLUSTER_AGE_RANGE_MYR[1] / 1000.0
_MIN_AGE_GY = tuning.OPEN_CLUSTER_AGE_RANGE_MYR[0] / 1000.0


def _stream(galaxy_seed):
    """The stream every draw of the young clusters comes from: the seed, bytes or an int, under the system's own name."""
    text = galaxy_seed.hex() if isinstance(galaxy_seed, (bytes, bytearray)) else str(galaxy_seed)
    return draw.Stream(f"open-clusters:{text}")


def _normal(stream, count):
    """Standard normal draws (Box-Muller) from `stream`."""
    u = stream.random_array(2 * count).reshape(2, count)
    return np.sqrt(-2.0 * np.log(np.maximum(u[0], 1.0e-300))) * np.cos(2.0 * math.pi * u[1])


def _disk_population(age_gy):
    """`young` or `intermediate` (the disk population stars of this age belong to), by the cluster's age."""
    return np.where(np.asarray(age_gy) < tuning.OPEN_CLUSTER_DISK_POPULATION_EDGE_GY, 0, 1)


def _disk_maximum_masses(found, radii, efficiency):
    """
    The most massive cluster each cell of the formation table can form, solar masses (`physics.cluster_maximum_mass`, GEN.1196):
    from the ring's gas (`found.gas_now`, the gas of the last Gy), the epicyclic frequency of the flat rotation curve at its
    radius (`tuning.CLUSTER_HALO_CIRCULAR_KMS`), the cell's bound share (its cluster formation efficiency) and the disc's gas
    dispersion (`tuning.DISC_GAS_DISPERSION_KMS`). A ring with no gas forms none above the least the mass function holds.
    """
    kappa = np.array([cluster_maximum_mass.epicyclic_frequency_per_myr(tuning.CLUSTER_HALO_CIRCULAR_KMS, float(r))
                      for r in radii])
    gas = np.asarray(found.gas_now, dtype=float)

    def cap(ring, bound):
        if gas[ring] <= 0.0:
            return 0.0
        return cluster_maximum_mass.maximum_cluster_mass_sol(
            float(gas[ring]), float(kappa[ring]), float(bound), tuning.DISC_GAS_DISPERSION_KMS)
    rings = np.arange(len(radii))[:, None]
    return np.vectorize(cap, otypes=[float])(rings, efficiency)


def _formation_table(shape):
    """
    `(ages_gy, step_gy, radii_pc, ring_width_pc, mass_formed, efficiency, maximum_mass)` for the steps of the last Gy and every
    ring of the disk's history: the mass formed in each cell of the table (solar masses), the cluster formation efficiency of
    its rate density and the most massive cluster its gas can form (`_disk_maximum_masses`).
    """
    found = disk_history.history_of(shape)
    steps = np.nonzero(found.age_gy <= _MAX_AGE_GY)[0]
    if len(steps) < 2:
        steps = np.arange(max(len(found.age_gy) - 2, 0), len(found.age_gy))
    step_gy = float(found.age_gy[0]) / (len(found.age_gy) - 1)
    radii = found.radius_scale_lengths * shape.disk_scale_length_pc
    ring_width = float(radii[1] - radii[0]) if len(radii) > 1 else float(radii[0])
    area_pc2 = 2.0 * math.pi * radii * ring_width
    rate = found.rate[:, steps]                                  # solar masses per Gy per pc2
    mass_formed = rate * area_pc2[:, None] * step_gy
    efficiency = evolution.cluster_formation_efficiency(rate * 1.0e-3)     # per year per kpc2 from per Gy per pc2
    return (found.age_gy[steps], step_gy, radii, ring_width, mass_formed, efficiency,
            _disk_maximum_masses(found, radii, efficiency))


def _arm_phase(radius_pc, shape):
    return (shape.spiral_reference_angle_rad
            + np.log(radius_pc / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad))


def _dispersion_kms(shape, age_gy):
    """Each cluster's one-axis velocity dispersion, km/s: the thin disc's age heating of the cluster's bin
    (`galaxy.structure.disc_dispersion_kms` is the three-axis total), isotropic."""
    bins = np.clip(np.searchsorted(np.asarray(tuning.OPEN_CLUSTER_AGE_EDGES_MYR), age_gy * 1000.0, side="right") - 1,
                   0, len(tuning.OPEN_CLUSTER_POPULATIONS) - 1)
    per_bin = np.array([galaxy_structure.disc_dispersion_kms(shape.structure, name) / constants.ISOTROPIC_TOTAL_OVER_ONE_AXIS
                        for name in tuning.OPEN_CLUSTER_POPULATIONS])
    return per_bin[bins]


def _disk_density(x, y, z, young, shape):
    """The `young` density where `young` is true, the `intermediate` density elsewhere, at arrays of points."""
    found = galaxy_density.population_densities_array(x, y, z, shape)
    return np.maximum(np.where(young, found["young"], found["intermediate"]), 0.0)


def _place(stream, radius_pc, age_gy, shape):
    """
    `(x, y, z)` arrays: each cluster at its ring's radius, at the height of stars of its age (the population's
    `exp(-|z| / h)` profile) and at an azimuth drawn by rejection against the density of stars of that age, which carries the
    arms (the youngest the most).
    """
    count = len(radius_pc)
    young = age_gy < tuning.OPEN_CLUSTER_DISK_POPULATION_EDGE_GY
    ratio = galaxy_density.population_height_ratios(shape)
    height = shape.disk_scale_height_pc * np.where(young, ratio["young"], ratio["intermediate"])
    u = np.clip(2.0 * stream.random_array(count) - 1.0, -1.0 + 1.0e-12, 1.0 - 1.0e-12)
    z = -np.sign(u) * height * np.log1p(-np.abs(u))      # the inverse of the profile's two-sided distribution
    phase = _arm_phase(radius_pc, shape)
    arms = max(int(shape.arm_count), 1)
    phi = stream.random_array(count) * 2.0 * math.pi
    pending = np.arange(count)
    for _ in range(_ARM_ROUNDS):
        if not len(pending):
            break
        r, zz, yy = radius_pc[pending], z[pending], young[pending]
        here = _disk_density(r * np.cos(phi[pending]), r * np.sin(phi[pending]), zz, yy, shape)
        best = np.zeros(len(pending))
        for offset in (0.0, math.pi / arms, math.pi / (2.0 * arms)):
            angle = phase[pending] + offset
            best = np.maximum(best, _disk_density(r * np.cos(angle), r * np.sin(angle), zz, yy, shape))
        keep = stream.random_array(len(pending)) * best <= here
        keep |= best <= 0.0
        pending = pending[~keep]
        if len(pending):
            phi[pending] = stream.random_array(len(pending)) * 2.0 * math.pi
    return radius_pc * np.cos(phi), radius_pc * np.sin(phi), z


def _fe_h(stream, radius_pc, age_gy):
    """Each cluster's [Fe/H]: the disk's mean at its radius for stars of its age, with that population's scatter."""
    distance_ly = pc_to_ly(radius_pc)
    pick = _disk_population(age_gy)
    out = np.zeros(len(radius_pc))
    for index, name in enumerate(("young", "intermediate")):
        row = tuning.METALLICITY_BY_POPULATION[name]
        gradient = metallicity.disc_gradient_fe_h(distance_ly) if row["gradient"] else 0.0
        out = np.where(pick == index, row["offset"] + gradient + row["sigma"] * _normal(stream, len(radius_pc)), out)
    low, high = tuning.METALLICITY_RANGE_FE_H
    return np.clip(out, low, high)


def ambient_density_msun_pc3(x, y, z, shape, expected_at_density_1, edge_pc):
    """
    The mean density of stars where the clusters lie, solar masses per cubic parsec: the model's stellar density
    (relative to its calibration point, `galaxy.density.total_density_array`) times the systems a cell holds at
    density 1 over the cell's volume, times the mean system mass. Gas is left out, so the density is lower than the
    stars-and-gas one Lamers et al. 2005 use and the disruption slightly slower.
    """
    relative = galaxy_density.total_density_array(x, y, z, shape)
    return relative * expected_at_density_1 / edge_pc ** 3 * galaxy_density.mean_system_mass_solar(shape)


def _ring_maximum_mass_sol(pot, core):
    """
    The most massive cluster the central ring forms, solar masses (`physics.cluster_maximum_mass`, GEN.1195): from the
    ring's gas surface density and the epicyclic frequency of the nucleus's circular speed at the ring's radius, and the
    share of the new stars that are bound. No ring forms none (0); a ring of no gas data is held to the mass function's top.
    """
    if core.rate_msun_yr <= 0.0 or core.gas_surface_density_msun_pc2 <= 0.0:
        return tuning.OPEN_CLUSTER_MASS_RANGE_SOL[1]
    speed = float(potential.circular_speed_kms(pot, core.radius_pc))
    kappa = cluster_maximum_mass.epicyclic_frequency_per_myr(speed, core.radius_pc)
    return cluster_maximum_mass.maximum_cluster_mass_sol(core.gas_surface_density_msun_pc2, kappa, core.bound_share,
                                                         tuning.CORE_RING_GAS_DISPERSION_KMS)


def _core_ring(stream, shape, core, pot, max_mass_sol=None):
    """
    The clusters the central ring forms (`galaxy.core_gas`): `(age_gy, x, y, z, fe_h, speed_kms, dispersion_kms)` arrays
    of those formed over `tuning.CORE_STAR_FORMATION_SPAN_GY` at the ring's rate (the share of it that goes into bound
    clusters, the mass function's mean mass for the count), at the ring's radius, width and height, with the nuclear
    disc's [Fe/H] and dispersion, on the circular orbit of the nucleus (the galaxy's potential `pot`, GEN.338, bulge and black hole included). Their initial masses are drawn with the others,
    under `max_mass_sol`, the most massive cluster the ring can form (the count is of the mass function held to it).
    """
    span_gy = min(tuning.CORE_STAR_FORMATION_SPAN_GY, _MAX_AGE_GY)
    mean_count = core.rate_msun_yr * core.bound_share * span_gy * 1.0e9 / evolution.mean_initial_mass_sol(max_mass_sol)
    count = poisson_count(mean_count, stream) if mean_count > 0.0 else 0
    age = np.maximum(stream.random_array(count) * span_gy, 1.0e-9)
    radius = np.abs(core.radius_pc + core.sigma_pc * _normal(stream, count))
    phi = stream.random_array(count) * 2.0 * math.pi
    z = 2.0 * core.height_pc * np.arctanh(np.clip(2.0 * stream.random_array(count) - 1.0, -1.0 + 1.0e-12, 1.0 - 1.0e-12))
    row = tuning.METALLICITY_BY_POPULATION["nuclear_disc"]
    fe_h = np.clip(row["offset"] + row["sigma"] * _normal(stream, count), *tuning.METALLICITY_RANGE_FE_H)
    speed = potential.circular_speed_kms(pot, radius)
    dispersion = np.full(count, tuning.STELLAR_POPULATION_DISPERSION_KMS["nuclear_disc"])
    return age, radius * np.cos(phi), radius * np.sin(phi), z, fe_h, speed, dispersion


def _disk_formed(stream, shape, pot):
    """
    `(age_gy, initial_mass, (x, y, z), fe_h, speed_kms, dispersion_kms)` of the clusters the disk's star formation history
    forms (see the module docstring, steps 1 and 2): empty arrays for a disk that formed none.
    """
    empty = (np.zeros(0), np.zeros(0), (np.zeros(0), np.zeros(0), np.zeros(0)), np.zeros(0), np.zeros(0), np.zeros(0))
    ages_gy, step_gy, radii, ring_width, mass_formed, efficiency, maximum = _formation_table(shape)
    weight = mass_formed * efficiency / evolution.mean_initial_mass_sol(maximum)       # clusters formed per cell
    total = float(weight.sum())
    if total <= 0.0:
        return empty
    floor = None
    if total > tuning.OPEN_CLUSTER_MAX_COUNT:      # a program limit: keep the heavier ones (`tuning.OPEN_CLUSTER_MAX_COUNT`)
        floor = evolution.low_mass_keeping_share(tuning.OPEN_CLUSTER_MAX_COUNT / total)
        weight = weight * (tuning.OPEN_CLUSTER_MAX_COUNT / total)
        total = float(tuning.OPEN_CLUSTER_MAX_COUNT)
    count = poisson_count(total, stream)
    if count == 0:
        return empty
    cumulative = np.cumsum(weight.reshape(-1))
    cumulative /= cumulative[-1]
    cell = np.minimum(np.searchsorted(cumulative, stream.random_array(count), side="right"), cumulative.size - 1)
    ring, step = np.divmod(cell, weight.shape[1])
    age = np.minimum(np.maximum(ages_gy[step] + (stream.random_array(count) - 0.5) * step_gy, 0.0), _MAX_AGE_GY)
    age = np.maximum(age, 1.0e-9)       # a cluster younger than `_MIN_AGE_GY` is still embedded (GEN.241's nebulae), not drawn
    radius = radii[ring] + (stream.random_array(count) - 0.5) * ring_width
    radius = np.maximum(radius, 0.5 * ring_width)
    initial = evolution.initial_masses(stream.random_array(count), floor, maximum[ring, step])
    place = _place(stream, radius, age, shape)
    fe_h = _fe_h(stream, radius, age)
    return (age, initial, place, fe_h, potential.circular_speed_kms(pot, radius), _dispersion_kms(shape, age))


def generate_open_clusters(galaxy_seed, shape, edge_pc, expected_at_density_1, epoch_unix=0.0):
    """
    Draws a galaxy's young open clusters from its seed and shape.

    Args:
        galaxy_seed (bytes or int): The galaxy's seed; the same seed and shape give the same clusters.
        shape (GalaxyShape): The galaxy's shape (its disk and its star formation history).
        edge_pc (float): The sector edge, parsecs.
        expected_at_density_1 (float): The systems a cell holds at `relative_density` 1.
        epoch_unix (float): When the positions and velocities hold.

    Returns:
        list[Component]: The clusters alive today, named "Open cluster 00001" and on, oldest draw first.
    """
    stream = _stream(galaxy_seed)
    c_env_myr = stream.uniform(*tuning.OPEN_CLUSTER_DISRUPTION_CENV_MYR)
    seed_bytes = potential.seed_bytes(galaxy_seed)
    pot = potential.for_galaxy(shape, seed_bytes)
    disk = _disk_formed(stream, shape, pot)
    core = core_gas.core_gas(shape)
    ring_cap = _ring_maximum_mass_sol(pot, core)
    ring = _core_ring(stream, shape, core, pot, ring_cap)
    ring_initial = evolution.initial_masses(stream.random_array(len(ring[0])), high_mass_sol=ring_cap)
    age = np.concatenate((disk[0], ring[0]))
    initial = np.concatenate((disk[1], ring_initial))
    x, y, z = (np.concatenate((disk[2][k], ring[1 + k])) for k in range(3))
    fe_h = np.concatenate((disk[3], ring[4]))
    speed = np.concatenate((disk[4], ring[5]))
    dispersion = np.concatenate((disk[5], ring[6]))
    count = len(age)
    if count == 0:
        return []
    rho = ambient_density_msun_pc3(x, y, z, shape, expected_at_density_1, edge_pc)
    t0 = evolution.disruption_t0_yr(rho, c_env_myr)
    metals = tuning.OPEN_CLUSTER_SOLAR_Z * 10.0 ** fe_h
    mass = evolution.present_masses(initial, age * 1.0e9, t0, metals)
    density_z = _normal(stream, count)
    expansion = stream.random_array(count)
    in_ring = np.arange(count) >= len(disk[0])
    mass = np.where(in_ring, np.maximum(mass - evolution.cloud_shock_mass_lost_sol(
        evolution.birth_half_mass_radii_pc(initial, density_z), evolution.expansion_factors(expansion), age,
        tuning.CORE_RING_CLOUD_SURFACE_DENSITY_MSUN_PC2, tuning.CORE_RING_CLOUD_DENSITY_MSUN_PC3), 0.0), mass)
    galactocentric = np.sqrt(x * x + y * y + z * z)
    size = evolution.half_mass_radii_pc(initial, age, density_z, expansion)
    size = np.minimum(size, evolution.jacobi_radii_pc(mass, galactocentric, enclosed_mass.enclosed_mass_solar(shape, seed_bytes, np.maximum(galactocentric, 1.0))) / tuning.SATELLITE_TIDAL_OVER_HALF_MASS)
    kept = np.nonzero((mass >= tuning.OPEN_CLUSTER_MASS_RANGE_SOL[0]) & (size > 0.0) & (age >= _MIN_AGE_GY))[0]
    kicks = np.stack([_normal(stream, count), _normal(stream, count), _normal(stream, count)]) * dispersion
    cylindrical = np.maximum(np.hypot(x, y), 1.0)
    velocity = (-speed * y / cylindrical + kicks[0], speed * x / cylindrical + kicks[1], kicks[2])
    seeds = [stream.getrandbits(63) for _ in range(len(kept))]
    found = []
    for number, i in enumerate(kept):
        found.append(components.Component(
            kind="cluster", name=f"Open cluster {number + 1:05d}", source="generated",
            center_pc=(float(x[i]), float(y[i]), float(z[i])),
            velocity_kms=(float(velocity[0][i]), float(velocity[1][i]), float(velocity[2][i])),
            epoch_unix=epoch_unix, mass_solar=float(mass[i]), half_mass_radius_pc=float(size[i]), profile="plummer",
            age_gy=float(age[i]), metallicity_fe_h=float(fe_h[i]), seed=seeds[number],
            tidal_radius_pc=component_profiles.cut_radius_pc("plummer", float(size[i]))))
    return found
