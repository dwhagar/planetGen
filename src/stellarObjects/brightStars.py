"""
Bright-star pre-placement: every star at least
`program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL` bright is drawn and
placed galaxy-wide right after `generate.py plan`, into `bright_stars`,
while every sector stays unfilled. A sector's fill later builds a full
system around each of its pre-placed stars (`fill_context`, used by
`generate.generate_sector`) and draws the rest of its systems from dimmer
stars only, so its expected total is unchanged.

The scatter works ring by ring, never sector by sector (a full galaxy
has billions of cells): each ring's expected count per stellar
population is averaged over angle bins, drawn as a Poisson count, and
each star then lands in a qualifying slot of a bin picked in proportion
to that population's density.
"""

import math
import random

from . import program_constants
from .galaxyDensity import population_densities, predicted_star_count
from .galaxyGeometry import (
    galaxy_to_local_pc,
    layer_bounds_pc,
    layer_center_z_pc,
    ring_bounds_pc,
    ring_radius_pc,
    ring_sector_count,
    sector_address_at,
    sector_position_pc,
)
from .spaceSector import _sample_poisson_count
from .stellarPopulation import bright_star_fraction, pick_population, sample_bright_stars
from .utils import pc_to_ly

POPULATIONS = ("young", "intermediate", "old", "bulge")
"""tuple: The stellar populations `galaxyDensity.population_densities`
splits a position's density into, each scattered on its own."""

ANGLE_BINS = 32
"""int: How many angle bins a ring's density is averaged over (fewer when
the ring has fewer slots). The spiral arms' density varies smoothly over
a bin this size."""

SLOT_REDRAWS = 8
"""int: How many angles a star tries before giving up when it keeps
landing in a slot too sparse to ever be filled (only at a qualifying
edge, where a bin straddles it)."""

MPC_PER_PC = 1000


def _densities(position_pc, shape):
    """`population_densities`, with a negative term (a shape with a
    negative bulge amplitude) counted as none."""
    return {population: max(density, 0.0) for population, density in population_densities(position_pc, shape).items()}


def _qualifies(position_pc, shape, expected_at_density_1):
    return predicted_star_count(position_pc, shape, expected_at_density_1) >= 1.0


def _ring_bins(ring_index, layer_index, shape, expected_at_density_1, edge_pc):
    """Per angle bin, the population densities at the ring's centerline
    (zeroed where the bin's own center doesn't qualify)."""
    slots = ring_sector_count(ring_index)
    count = min(slots, ANGLE_BINS)
    radius = ring_radius_pc(ring_index, edge_pc)
    z = layer_center_z_pc(layer_index, edge_pc)
    bins = []
    for k in range(count):
        theta = (k + 0.5) * 2 * math.pi / count
        point = (radius * math.cos(theta), radius * math.sin(theta), z)
        if _qualifies(point, shape, expected_at_density_1):
            bins.append(_densities(point, shape))
        else:
            bins.append(None)
    return slots, bins


def _place_one(rng, weights, ring_index, layer_index, slots, shape, expected_at_density_1, edge_pc):
    """A uniform point in a qualifying slot of a bin picked by `weights`,
    as `(slot, (x, y, z))` in whole milliparsecs (as stored), or `None` if
    every try landed in a slot too sparse to be filled. A point whose
    rounding carries it over the cell's edge is redrawn too."""
    total = sum(weights)
    bin_width = 2 * math.pi / len(weights)
    slot_width = 2 * math.pi / slots
    r_low, r_high = ring_bounds_pc(ring_index, edge_pc)
    z_low, z_high = layer_bounds_pc(layer_index, edge_pc)
    for _ in range(SLOT_REDRAWS):
        pick = rng.random() * total
        for k, weight in enumerate(weights):
            if pick < weight:
                break
            pick -= weight
        theta = (k + rng.random()) * bin_width
        slot = min(int(theta / slot_width), slots - 1)
        if not _qualifies(sector_position_pc(ring_index, layer_index, slot, edge_pc), shape, expected_at_density_1):
            continue
        theta = (slot + rng.random()) * slot_width
        # Uniform over the cell's area, which grows with radius.
        radius = math.sqrt(rng.uniform(r_low * r_low, r_high * r_high))
        point = (radius * math.cos(theta), radius * math.sin(theta), rng.uniform(z_low, z_high))
        stored = tuple(round(value * MPC_PER_PC) for value in point)
        if sector_address_at(tuple(value / MPC_PER_PC for value in stored), edge_pc) == (ring_index, layer_index, slot):
            return slot, stored
    return None


def scatter(shape, extents, edge_pc, expected_at_density_1, min_luminosity_sol, seed,
            skip_addresses=None, on_layer=None):
    """
    Draws and places every bright star in the galaxy's outline.

    Args:
        shape (GalaxyShape): The galaxy's shape.
        extents (list): `(layer_index, outer_ring_index)` per layer, as
            `galaxySkeleton.build_layer_extents` returns.
        edge_pc (float): The sector edge.
        expected_at_density_1 (float): Systems per sector at density 1.
        min_luminosity_sol (float): The threshold, in solar luminosities.
        seed (int): The scatter's seed; the same seed and outline give the
            same stars.
        skip_addresses (set, optional): `(ring, layer, slot)` cells to leave
            out (sectors already filled).
        on_layer (callable, optional): Called as `on_layer(done, total)`
            after each layer.

    Yields:
        tuple: One row per star, in `_db.BRIGHT_STAR_COLUMNS` order.
    """
    rng = random.Random(seed)
    skip_addresses = skip_addresses or set()
    fractions = {population: bright_star_fraction(min_luminosity_sol, population) for population in POPULATIONS}
    for done, (layer_index, outer_ring) in enumerate(extents, start=1):
        placed = {population: [] for population in POPULATIONS}
        for ring_index in range(outer_ring + 1):
            slots, bins = _ring_bins(ring_index, layer_index, shape, expected_at_density_1, edge_pc)
            slots_per_bin = slots / len(bins)
            for population in POPULATIONS:
                weights = [densities[population] if densities else 0.0 for densities in bins]
                mean = expected_at_density_1 * slots_per_bin * sum(weights) * fractions[population]
                if mean <= 0.0:
                    continue
                for _ in range(_sample_poisson_count(mean, rng=rng)):
                    spot = _place_one(rng, weights, ring_index, layer_index, slots, shape,
                                      expected_at_density_1, edge_pc)
                    if spot is None or (ring_index, layer_index, spot[0]) in skip_addresses:
                        continue
                    placed[population].append((ring_index, spot[0], spot[1]))
        for population, spots in placed.items():
            if not spots:
                continue
            stars = sample_bright_stars(len(spots), min_luminosity_sol, population, rng)
            for (ring_index, slot, (x, y, z)), params in zip(spots, stars):
                yield (
                    ring_index, layer_index, slot,
                    x, y, z,
                    population, params["type"], params["yerkes_class"], params["mass_kg"],
                    params["radius_km"], params["temperature_k"], params["luminosity_w"],
                    params["age_gy"], params["lifespan_gy"], params["initial_mass_sol"],
                    params["phase_end_age_gy"], rng.getrandbits(63),
                )
        if on_layer is not None:
            on_layer(done, len(extents))


def star_params(row):
    """A `bright_stars` row back as the `stellarEvolution.star_params`
    dict `Star.from_params` takes."""
    return {
        "type": row["star_type"], "yerkes_class": row["yerkes_class"], "mass_kg": row["mass_kg"],
        "radius_km": row["radius_km"], "temperature_k": row["temperature_k"],
        "luminosity_w": row["luminosity_w"], "age_gy": row["age_gy"], "lifespan_gy": row["lifespan_gy"],
        "initial_mass_sol": row["initial_mass_sol"], "phase_end_age_gy": row["phase_end_age_gy"],
    }


def local_position_ly(row, center_pc):
    """A `bright_stars` row's stored point as a sector-local `(x, y, z)`
    in light-years (the frame `SpaceSector` places systems in)."""
    point = (row["position_x_mpc"] / MPC_PER_PC, row["position_y_mpc"] / MPC_PER_PC,
             row["position_z_mpc"] / MPC_PER_PC)
    return tuple(pc_to_ly(value) for value in galaxy_to_local_pc(center_pc, point))


class FillContext:
    """
    What a galaxy-placed sector's fill needs beyond its arguments: the
    stellar population mix at its center, and its pre-placed bright stars
    with the threshold they were scattered at.

    Attributes:
        center_pc (tuple): The sector's center, galaxy-frame parsecs.
        densities (dict): `population_densities` at the center.
        bright_rows (list): The sector's unfilled `bright_stars` rows.
        min_luminosity_sol (float or None): The scatter's threshold;
            `None` when no scatter ran (no star is capped then).
    """

    def __init__(self, center_pc, shape, bright_rows=(), min_luminosity_sol=None):
        self.center_pc = center_pc
        self.densities = _densities(center_pc, shape)
        self.bright_rows = list(bright_rows)
        self.min_luminosity_sol = min_luminosity_sol

    def bright_share(self):
        """The share of this position's stars at or above the threshold
        (0 when no scatter ran): the part of the sector's expected count
        its pre-placed stars already stand for."""
        if self.min_luminosity_sol is None:
            return 0.0
        total = sum(self.densities.values())
        if total <= 0.0:
            return 0.0
        return sum(density * bright_star_fraction(self.min_luminosity_sol, population)
                   for population, density in self.densities.items()) / total

    def apply(self, system_config, rng=random):
        """Gives one of the sector's own (dim) systems its population and,
        after a scatter, the luminosity cap."""
        system_config.POPULATION = pick_population(self.densities, rng)
        if self.min_luminosity_sol is not None:
            system_config.MAX_STAR_LUMINOSITY_SOL = self.min_luminosity_sol
        return system_config
