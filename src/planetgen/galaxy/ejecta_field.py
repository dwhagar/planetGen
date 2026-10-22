# planetgen/galaxy/ejecta_field.py

"""
The galaxy's ejecta (GEN.238)
=============================

A neutron star born with a kick above the galaxy's escape speed, and a black hole ejected, leave it and stream out for
ever: a steady outflow that fills the whole sphere, with a density that falls as `1 / r^2`. The code reads them as one
population per kind, `ejecta`, apart from the bound stars and their remnants: the unbound share of the kind's galaxy
total `N` (a share of the neutron stars that the galaxy's escape speed sets, `galaxy.natal_kick`, about 2% of the black holes) leaves at the speed `v` (300 to 500 km/s
above escape) and has been leaving for `T` (10 Gyr), so

    n(r) = N / (4 pi r^2 v T)    for r < v T, and none beyond

(docs/design/intergalactic-and-halo-contents.md 5.6). Against the halo's bound neutron stars the ejecta take over
near 50 kpc: 1.3e-7 of a 4 pc cell at 50 kpc, 3.3e-8 at 100 kpc. They are not pre-placed (about 1e7 rows inside the
262 kpc reach would be needed); `SectorProfile.ejecta` holds the mean count a cell draws at fill.

Everything here is a pure function of the galaxy's shape, its sector edge and its seed, and holds no state.
"""

import functools
import math

from planetgen import tuning
from planetgen.galaxy import remnant_density
from planetgen.galaxy import density as galaxyDensity
from planetgen.galaxy import seed as galaxySeed
from planetgen.physics import constants as physical_constants
from planetgen.util import draw

KINDS = ("neutron-star", "black-hole")
"""tuple: The kinds that have an ejecta population."""

_KM_PER_PC = physical_constants.PARSEC_M / physical_constants.KM_M


@functools.lru_cache(maxsize=64)
def galaxy_systems(shape, expected_at_density_1, edge_pc):
    """How many systems the galaxy holds: its components' and its halo's integral of the density (unnormalized units,
    as `galaxy.density.component_masses`), times the normalization and the systems a cube of `edge_pc` holds at density
    1. The clusters' and satellites' own systems (a percent) are not added."""
    visible = galaxyDensity.model_terms(shape)["visible_mass"]
    return expected_at_density_1 * shape.k_norm * (visible + galaxyDensity.halo_mass(shape)) / edge_pc ** 3


def unbound_fraction(kind, shape=None):
    """The share of the galaxy's `kind` that is ejected: a black hole's constant, a neutron star's the galaxy's own
    (`Structure.neutron_star_unbound_share`, `galaxy.natal_kick`; the Milky Way structure's without a shape)."""
    if kind == "black-hole":
        return tuning.BLACK_HOLE_UNBOUND_FRACTION
    from planetgen.galaxy import structure
    return (structure.REFERENCE_STRUCTURE if shape is None else shape.structure).neutron_star_unbound_share


def unbound_count(kind, shape, expected_at_density_1, edge_pc):
    """The ejected `kind` the galaxy holds: the unbound share of what it formed. What the galaxy's metallicity makes of
    each kind (`remnant_density.formed_per_star`, the bound ones and the ejected together), so bound plus ejected is the
    whole (GEN.563 and GEN.554)."""
    formed = remnant_density.formed_per_star(shape, kind)
    unbound = unbound_fraction(kind, shape)
    return formed * galaxy_systems(shape, expected_at_density_1, edge_pc) * unbound


def mean_speed_kms():
    """The speed the population's density is taken at, km/s: the middle of `EJECTA_SPEED_RANGE_KMS`."""
    return sum(tuning.EJECTA_SPEED_RANGE_KMS) / 2.0


def outer_radius_pc():
    """Where the population ends, `v T`, parsecs."""
    return mean_speed_kms() * tuning.EJECTA_AGE_GY * 1.0e9 * physical_constants.SECONDS_PER_YEAR / _KM_PER_PC


def density_pc3(kind, radius_pc, shape, expected_at_density_1, edge_pc):
    """The ejecta of `kind` per cubic parsec `radius_pc` from the galactic centre."""
    outer = outer_radius_pc()
    if radius_pc >= outer or radius_pc <= 0.0:
        return 0.0
    return unbound_count(kind, shape, expected_at_density_1, edge_pc) / (
        4.0 * math.pi * radius_pc ** 2 * outer)


def expected_in_cell(position_pc, volume_pc3, shape, expected_at_density_1, edge_pc):
    """`{kind: mean}`: how many ejecta of each kind a cell of `volume_pc3` centred on `position_pc` holds on average,
    the mean of its Poisson draw, at most `tuning.EJECTA_MAX_PER_CELL` a kind."""
    radius = math.sqrt(sum(c * c for c in position_pc))
    return {kind: min(density_pc3(kind, radius, shape, expected_at_density_1, edge_pc) * volume_pc3,
                      tuning.EJECTA_MAX_PER_CELL)
            for kind in KINDS}


def draw_velocity_kms(rng=draw):
    """The velocity of one ejected object, galactic axes, km/s: a speed uniform in `EJECTA_SPEED_RANGE_KMS` along a
    random direction (`rng` supplies `random`)."""
    z = 2.0 * rng.random() - 1.0
    phi = 2.0 * math.pi * rng.random()
    ring = math.sqrt(max(0.0, 1.0 - z * z))
    low, high = tuning.EJECTA_SPEED_RANGE_KMS
    speed = low + (high - low) * rng.random()
    return (speed * ring * math.cos(phi), speed * ring * math.sin(phi), speed * z)
