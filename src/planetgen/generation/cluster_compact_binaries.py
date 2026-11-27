# planetgen/generation/cluster_compact_binaries.py

"""
Compact binaries of the globular clusters (GEN.330; docs/design/systems/multi-star-systems.md 13.4)
===================================================================================================

A black hole, a neutron star or a millisecond pulsar of the field has the companion odds of `multistar.MEMBER_STATE_BY_PLACED`. In a globular cluster the odds follow the
cluster: the dense core swaps, hardens and ejects, so a black hole is likelier to hold a star, a millisecond pulsar is recycled by a companion it keeps (and more often ends
beside a neutron star), and a black hole's black hole is formed by the density. This module gives those odds for one cluster (`tables`) and the widest pair the cluster keeps
(`axis_limit_au`), and the scatter's companion draw (`scatter_companions.phenomenon_companion`, `phenomenon_tree`) reads them for a row that carries a `cluster_id`. A row of
the field is drawn exactly as before.

**The odds** (stand-ins 6.340 CB1 to CB5, `tuning.CLUSTER_*`):

* a black hole: a star `CLUSTER_BLACK_HOLE_STAR_ODDS` (10 percent), a black hole `black_hole_pair_odds` (density-driven, at most `CLUSTER_BLACK_HOLE_PAIR_MAX_ODDS`), a neutron star as the field's;
* a millisecond pulsar: a companion of any kind `CLUSTER_MILLISECOND_PULSAR_COMPANION_ODDS` (80 percent) with a neutron star `CLUSTER_MILLISECOND_PULSAR_NEUTRON_STAR_ODDS`
  and the other states in the field's proportions;
* every other class (the young neutron stars, the white dwarfs, the stars) as the field's.

**The soft pairs** are cut by the hard-soft rule of the cluster's planets (GEN.162:a, `cluster_unbinding`): a pair whose circular speed is under the relative speed of the cluster's
stars (the virial dispersion of its mass and half-mass radius, `virial_dispersion_km_s`, `relative_speed_km_s`) is soft and has been ionised, so the tree is drawn with
`max_axis_au` at the orbit where the two speeds are equal (`hard_soft_axis_au`). A pair inside it is hard and survives.

Not drawn from the cluster: the low-mass X-ray binary share of a black hole's star companions (`multistar.draw_low_mass_xray_pair`) stays the field's share of the star odds, so
a cluster black hole's tenfold star odds make tenfold the X-ray binaries per black hole, the over-representation the clusters show.
"""

import functools

from planetgen import tuning
from planetgen.generation import mass_ratio, multistar
from planetgen.physics import cluster_unbinding

def black_hole_pair_odds(mass_solar, half_mass_radius_pc):
    """
    The share of a cluster's (of this mass and half-mass radius) black holes that have a black hole for a companion (GEN.330, CB4 and CB5): `MAX x / (1 + x)` of the half-mass density
    over `tuning.CLUSTER_BLACK_HOLE_PAIR_DENSITY_REFERENCE_SOLAR_PC3`, so it follows the density and saturates at `tuning.CLUSTER_BLACK_HOLE_PAIR_MAX_ODDS`.
    """
    density = cluster_unbinding.half_mass_density(mass_solar, half_mass_radius_pc)
    x = density / tuning.CLUSTER_BLACK_HOLE_PAIR_DENSITY_REFERENCE_SOLAR_PC3
    return tuning.CLUSTER_BLACK_HOLE_PAIR_MAX_ODDS * x / (1.0 + x)


def _recycled_pulsar_table():
    """The cluster's millisecond pulsar table: the field's, its non-neutron-star states rescaled so the total is the cluster's and the neutron star's share is the cluster's."""
    field = multistar.MEMBER_STATE_BY_PLACED["recycled-pulsar"]
    fixed = tuning.CLUSTER_MILLISECOND_PULSAR_NEUTRON_STAR_ODDS
    others = {state: row for state, row in field.items() if state != "neutron-star"}
    scale = (tuning.CLUSTER_MILLISECOND_PULSAR_COMPANION_ODDS - fixed) / sum(row[0] for row in others.values())
    table = {}
    for state, (default, low, high) in field.items():
        if state == "neutron-star":
            table[state] = (fixed, low, high)
        else:
            table[state] = (default * scale, low * scale, high * scale)
    return table


@functools.lru_cache(maxsize=4096)
def tables_of(mass_solar, half_mass_radius_pc):
    """`tables` for a cluster's mass and radius (the cache key): the field's dict with the cluster's classes replaced."""
    tables = dict(multistar.MEMBER_STATE_BY_PLACED)
    hole = dict(tables["black-hole"])
    star = hole["star"]
    hole["star"] = (tuning.CLUSTER_BLACK_HOLE_STAR_ODDS, star[1], star[2])
    pair = hole["black-hole"]
    hole["black-hole"] = (black_hole_pair_odds(mass_solar, half_mass_radius_pc), pair[1], pair[2])
    tables["black-hole"] = hole
    tables["recycled-pulsar"] = _recycled_pulsar_table()
    return tables


def tables(cluster):
    """
    `multistar.MEMBER_STATE_BY_PLACED` for a placed object in `cluster` (a `galaxy.components.Component`): the field's table with the black hole's and the millisecond pulsar's
    classes replaced by the cluster's. `cluster` `None` is the field's table itself.
    """
    if cluster is None:
        return multistar.MEMBER_STATE_BY_PLACED
    return tables_of(float(cluster.mass_solar), float(cluster.half_mass_radius_pc))


def axis_limit_au(cluster, placed_mass_solar):
    """
    The widest axis a pair of a placed object of `placed_mass_solar` keeps in `cluster` (the hard-soft rule of GEN.162:a applied to the pair): the orbit whose circular speed about
    the placed object and a companion of the mean mass ratio equals the relative speed of the cluster's stars (`cluster_unbinding.relative_speed_km_s` of the virial dispersion), `G M / V^2`.
    A wider pair is soft against the cluster's encounters and has been ionised (a pair inside it is hard and survives, hardened by them).
    """
    speed = cluster_unbinding.relative_speed_km_s(cluster_unbinding.virial_dispersion_km_s(cluster.mass_solar, cluster.half_mass_radius_pc))
    total = placed_mass_solar * (1.0 + mass_ratio.mean_ratio(placed_mass_solar))
    return cluster_unbinding.hard_soft_axis_au(total, speed)
