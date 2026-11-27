# planetgen/generation/cluster_blue_stragglers.py

"""
The blue stragglers of the globular clusters (GEN.326)
======================================================

A blue straggler is a main-sequence star brighter and bluer than its cluster's turn-off: the product of a merger or a mass transfer, a
star of 1 to 1.7 solar masses in a cluster whose stars are 0.85 solar masses and under. The number in a cluster's core goes as the mass in
the core to the 0.38 power [R2397] (0.40 in [R5367]), almost independent of how often the stars meet at a fixed core mass, so the count
here is `N = N_1e4 (M_core / 1e4 solar masses)^0.38` with the normalization `tuning.CLUSTER_BLUE_STRAGGLERS_AT_1E4_CORE_MASS`, a stand-in (the
sources give the slope only). The core mass is the King model's mass inside its scale radius (`physics/king.py`), the cluster's own mass and
half-mass radius being the component's.

Each blue straggler is a row of the components' pass (`component_pass.scatter_blue_stragglers`, written beside the stars of `scatter_stars`), in the core: a position by the inverse of the King
profile's mass table inside the core, a mass from the initial mass function between 1 and 1.7 solar masses, and an age uniform up to the
time a blue straggler stays one (`tuning.BLUE_STRAGGLER_LIFETIME_GY`, 1.5 Gy) or the main-sequence life of its mass. The cluster's own stars
are drawn as before from their own stream; the blue stragglers come from a stream of their own, so adding them changes none of them. They are placed
whatever the tiers keep: a blue straggler is a few solar luminosities, under every tier's floor, and the sector fill draws a cluster's ordinary stars only. A
star outside a cluster is never one, a young open cluster and a satellite draw none, and a cluster in a sector the pass skips (`skip_addresses`)
keeps none there.
"""

import math

from planetgen import tuning
from planetgen.galaxy import component_field
from planetgen.generation import bright_stars
from planetgen.physics import constants, imf, king
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw
from planetgen.util.poisson import poisson_count

REFERENCE_CORE_MASS_SOLAR = 1.0e4
"""float: The core mass the normalization is given at, solar masses: derived, the 1e4 in the name of `tuning.CLUSTER_BLUE_STRAGGLERS_AT_1E4_CORE_MASS` (stand-in row BN2 of registry section 6.331); a unit of that law, not a second number."""


def has_blue_stragglers(component):
    """Whether `component` is a globular cluster (a King profile, not a young open cluster or a satellite)."""
    return component.kind == "cluster" and component.profile == "king"


def core_mass_solar(component):
    """The mass inside the King scale radius of a globular cluster, solar masses (its mass times the share of it the King model holds there)."""
    radius_pc = king.scale_radius_pc(component.king_w0, component.half_mass_radius_pc)
    return component.mass_solar * king.mass_within_pc(component.king_w0, component.half_mass_radius_pc, radius_pc)


def expected_count(component):
    """The mean number of blue stragglers of a cluster, 0 for a component that is not a globular cluster."""
    if not has_blue_stragglers(component):
        return 0.0
    return (tuning.CLUSTER_BLUE_STRAGGLERS_AT_1E4_CORE_MASS
            * (core_mass_solar(component) / REFERENCE_CORE_MASS_SOLAR) ** tuning.CLUSTER_BLUE_STRAGGLER_CORE_MASS_EXPONENT)


def _star(component, population, rng):
    """`(mass, effective age, params)` of one blue straggler of `component`: `star_params` of a star of the drawn mass at the drawn age, or
    `None` when that is not a main-sequence star (the age is under the time it takes to reach the main sequence)."""
    low, high = tuning.BLUE_STRAGGLER_MASS_RANGE_SOL
    mass = imf.sample(low, high, rng)
    fe_h = evolution.birth_fe_h(rng, population, component_field.fe_h_peaks_of(component))
    age = rng.uniform(0.0, min(tuning.BLUE_STRAGGLER_LIFETIME_GY, evolution.main_sequence_lifetime_gy(mass, fe_h)))
    state = evolution.evolve_star(mass, age, rng, fe_h=fe_h, population=population)
    if state is None or state["yerkes_class"] not in ("V", "IV"):
        return None         # a merger product that is still contracting, or a giant, is no blue straggler of this draw
    return mass, age, evolution.star_params(mass, age, state)


def blue_straggler_rows(component, edge_pc, seed, skip_addresses, positions_of, address_of):
    """
    The rows (`store.BRIGHT_STAR_COLUMNS` order) of the blue stragglers of one cluster: Poisson in `expected_count`, each in the core, from the
    stream `{seed}:objects:component:cluster:{name}:blue-stragglers`. `positions_of(component, uniform)` and `address_of(point, edge_pc)`
    are `component_pass.positions_mpc` and `component_pass._address`.
    """
    if not has_blue_stragglers(component):
        return
    rng = draw.Stream(f"{seed}:objects:component:{component.kind}:{component.key}:blue-stragglers")
    count = poisson_count(expected_count(component), rng)
    if count <= 0:
        return
    radius_pc = king.scale_radius_pc(component.king_w0, component.half_mass_radius_pc)
    in_core = king.mass_within_pc(component.king_w0, component.half_mass_radius_pc, radius_pc)
    population = component.kind
    uniform = rng.random_array(3 * count).reshape(3, count)
    uniform[0] *= in_core       # the inverse transform of the profile's mass table, over the core's share only
    points = positions_of(component, uniform)
    for index in range(count):
        point = tuple(int(value) for value in points[:, index])
        address = address_of(point, edge_pc)
        found = _star(component, population, rng)
        if found is None or address in skip_addresses:
            continue
        yield bright_stars._row(*address, point, population, found[2], rng, component.id)
