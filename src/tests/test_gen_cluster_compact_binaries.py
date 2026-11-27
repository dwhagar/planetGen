"""
GEN.330: the compact binaries of a globular cluster (`planetgen.generation.cluster_compact_binaries`, docs/design/systems/multi-star-systems.md 13.4). A black hole, neutron star or millisecond
pulsar row that carries a `cluster_id` draws its companion at the cluster's odds (a black hole with a star 10 percent, a millisecond pulsar with a companion 80 percent with more neutron stars,
a black hole with a black hole from the cluster's density), its pair is cut at the hard-soft orbit of the cluster, and the rows per system follow the cluster's odds; a row of the field is
drawn exactly as before.
"""

import dataclasses
import math
import random
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.generation import cluster_compact_binaries as cluster_binaries
from planetgen.generation import component_pass, multistar, phenomenon_scatter, scatter_companions
from planetgen.galaxy import components as galaxyComponents
from planetgen.physics import cluster_unbinding
from tests.test_gen_scatter_companions import PHEN, _hit_seeds, _phenomenon_row, _seeds

TYPICAL = SimpleNamespace(mass_solar=2.0e5, half_mass_radius_pc=5.0, age_gy=11.0)
COMPACT = SimpleNamespace(mass_solar=1.0e6, half_mass_radius_pc=3.0, age_gy=12.0)
SPARSE = SimpleNamespace(mass_solar=2.0e4, half_mass_radius_pc=20.0, age_gy=12.0)
DRAWS = 40000


def _share(hits, total):
    """`(share, four-sigma half width)` of `hits` in `total`."""
    share = hits / total
    return share, 4.0 * math.sqrt(max(share * (1.0 - share), 1.0e-12) / total)


def _within(measured, expected, total):
    return abs(measured - expected) <= 4.0 * math.sqrt(expected * (1.0 - expected) / total)


# --- The odds -----------------------------------------------------------------------------------------------------


def test_no_cluster_is_the_field_table_itself():
    assert cluster_binaries.tables(None) is multistar.MEMBER_STATE_BY_PLACED


def test_a_clusters_black_holes_hold_a_star_at_the_clusters_rate_and_every_other_class_stays_the_fields():
    table = cluster_binaries.tables(TYPICAL)
    assert table["black-hole"]["star"][0] == tuning.CLUSTER_BLACK_HOLE_STAR_ODDS == 0.10
    assert table["black-hole"]["neutron-star"] == multistar.MEMBER_STATE_BY_PLACED["black-hole"]["neutron-star"]
    for name in set(multistar.MEMBER_STATE_BY_PLACED) - {"black-hole", "recycled-pulsar"}:
        assert table[name] == multistar.MEMBER_STATE_BY_PLACED[name]


def test_a_clusters_millisecond_pulsars_have_a_companion_at_the_clusters_rate_with_more_neutron_stars():
    field = multistar.MEMBER_STATE_BY_PLACED["recycled-pulsar"]
    table = cluster_binaries.tables(TYPICAL)["recycled-pulsar"]
    assert sum(row[0] for row in table.values()) == pytest.approx(tuning.CLUSTER_MILLISECOND_PULSAR_COMPANION_ODDS)
    assert table["neutron-star"][0] == tuning.CLUSTER_MILLISECOND_PULSAR_NEUTRON_STAR_ODDS > field["neutron-star"][0]
    others = [state for state in field if state != "neutron-star"]
    ratios = {table[state][0] / field[state][0] for state in others}
    assert max(ratios) == pytest.approx(min(ratios))                   # the other states keep the field's proportions


def test_the_black_hole_pairs_follow_the_density_and_saturate():
    sparse, typical, compact = (cluster_binaries.tables(c)["black-hole"]["black-hole"][0] for c in (SPARSE, TYPICAL, COMPACT))
    field = multistar.MEMBER_STATE_BY_PLACED["black-hole"]["black-hole"][0]
    assert sparse < typical < compact < tuning.CLUSTER_BLACK_HOLE_PAIR_MAX_ODDS
    assert sparse < field < compact             # the field's 0.2 percent sits between the sparse and the dense clusters
    huge = SimpleNamespace(mass_solar=1.0e9, half_mass_radius_pc=1.0)
    assert cluster_binaries.black_hole_pair_odds(huge.mass_solar, huge.half_mass_radius_pc) == pytest.approx(tuning.CLUSTER_BLACK_HOLE_PAIR_MAX_ODDS, rel=1e-3)


# --- The screen ---------------------------------------------------------------------------------------------------


def test_a_cluster_black_hole_has_a_star_at_ten_percent_and_the_field_one_at_one():
    seeds = _seeds(DRAWS, start=301)
    for cluster, expected in ((TYPICAL, 0.10), (None, multistar.MEMBER_STATE_BY_PLACED["black-hole"]["star"][0])):
        stars = sum(1 for s in seeds if (scatter_companions.phenomenon_companion("black-hole", "stellar", s, cluster) or (0, 0, ""))[2] == "star")
        assert _within(stars / DRAWS, expected, DRAWS), (cluster, stars / DRAWS)


def test_a_cluster_millisecond_pulsar_has_a_companion_at_eighty_percent_and_a_neutron_star_at_two_and_a_half(monkeypatch):
    monkeypatch.setattr(phenomenon_scatter, "NEUTRON_STAR_RECYCLED_SHARE", 1.0)       # every neutron star a millisecond pulsar, so the draw is of the pulsar table alone
    seeds = _seeds(DRAWS, start=307)
    hits = [hit for s in seeds if (hit := scatter_companions.phenomenon_companion("neutron-star", None, s, TYPICAL)) is not None]
    assert all(hit[1] for hit in hits)
    assert _within(len(hits) / DRAWS, tuning.CLUSTER_MILLISECOND_PULSAR_COMPANION_ODDS, DRAWS)
    assert _within(sum(1 for hit in hits if hit[2] == "neutron-star") / DRAWS, tuning.CLUSTER_MILLISECOND_PULSAR_NEUTRON_STAR_ODDS, DRAWS)
    field = sum(1 for s in seeds if scatter_companions.phenomenon_companion("neutron-star", None, s) is not None)
    assert _within(field / DRAWS, sum(row[0] for row in multistar.MEMBER_STATE_BY_PLACED["recycled-pulsar"].values()), DRAWS)


def test_a_cluster_black_hole_has_a_black_hole_at_the_density_driven_rate():
    seeds = _seeds(DRAWS, start=311)
    for cluster in (SPARSE, COMPACT):
        expected = cluster_binaries.tables(cluster)["black-hole"]["black-hole"][0]
        pairs = sum(1 for s in seeds if (scatter_companions.phenomenon_companion("black-hole", "stellar", s, cluster) or (0, 0, ""))[2] == "black-hole")
        assert _within(pairs / DRAWS, expected, DRAWS), (cluster, pairs / DRAWS, expected)


def test_a_row_of_the_field_is_drawn_exactly_as_before():
    seeds = _seeds(3000, start=313)
    for kind, subtype in (("black-hole", "stellar"), ("neutron-star", None), ("black-hole", "stellar-pair"), ("neutron-star", "pair"), ("planetary-nebula", None)):
        field = [scatter_companions.phenomenon_companion(kind, subtype, s) for s in seeds]
        assert field == [scatter_companions.phenomenon_companion(kind, subtype, s, None) for s in seeds]
    row = _phenomenon_row("black-hole", "stellar", _hit_seeds("black-hole", "stellar", 1, state="star", start=317)[0], mass=16.0)
    assert scatter_companions.phenomenon_tree(row) == scatter_companions.phenomenon_tree(row, None)


def test_a_cluster_pair_row_is_always_a_pair_and_its_black_hole_share_follows_the_cluster():
    seeds = _seeds(DRAWS, start=319)
    hits = [scatter_companions.phenomenon_companion("black-hole", "stellar-pair", s, COMPACT) for s in seeds]
    assert all(hit is not None for hit in hits)
    table = cluster_binaries.tables(COMPACT)["black-hole"]
    expected = table["black-hole"][0] / (table["black-hole"][0] + table["neutron-star"][0])
    assert _within(sum(1 for hit in hits if hit[2] == "black-hole") / DRAWS, expected, DRAWS)


# --- The soft pairs ----------------------------------------------------------------------------------------------


def test_the_hard_soft_orbit_is_where_the_pairs_circular_speed_is_the_clusters_relative_speed():
    speed = cluster_unbinding.relative_speed_km_s(cluster_unbinding.virial_dispersion_km_s(TYPICAL.mass_solar, TYPICAL.half_mass_radius_pc))
    for mass in (1.4, 10.0, 30.0):
        limit = cluster_binaries.axis_limit_au(TYPICAL, mass)
        total = mass * (1.0 + cluster_binaries.mass_ratio.mean_ratio(mass))
        assert not cluster_unbinding.soft(limit * 0.99, total, speed) and cluster_unbinding.soft(limit * 1.01, total, speed)
    assert cluster_binaries.axis_limit_au(COMPACT, 10.0) < cluster_binaries.axis_limit_au(TYPICAL, 10.0) < cluster_binaries.axis_limit_au(SPARSE, 10.0)
    assert cluster_binaries.axis_limit_au(TYPICAL, 10.0) < cluster_binaries.axis_limit_au(TYPICAL, 30.0)


def _widest_axis(tree):
    return max(node.separation_au for node in tree.nodes)


def test_a_cluster_pair_is_never_wider_than_the_cluster_keeps_and_a_field_pair_can_be():
    limit = cluster_binaries.axis_limit_au(TYPICAL, 16.0)
    seeds = _hit_seeds("black-hole", "stellar", 150, state="star", start=7)
    rows = [_phenomenon_row("black-hole", "stellar", seed, mass=16.0) for seed in seeds]
    inside = [_widest_axis(tree) for row in rows if (tree := scatter_companions.phenomenon_tree(row, TYPICAL)) is not None]
    outside = [_widest_axis(tree) for row in rows if (tree := scatter_companions.phenomenon_tree(row)) is not None]
    assert inside and max(inside) <= limit * (1.0 + 1e-9)
    assert max(outside) > limit


def test_a_cluster_system_is_as_old_as_the_cluster_and_a_field_one_is_not():
    seed = _hit_seeds("black-hole", "stellar", 1, state="star", start=331)[0]
    row = _phenomenon_row("black-hole", "stellar", seed, mass=16.0)
    rng = scatter_companions.draw.Stream(f"{seed}:members")
    rng.uniform(*tuning.COMPACT_REMNANT_AGE_RANGE_GY)
    by_hand = multistar.draw_member_tree(multistar.Placed("black-hole", 16.0), TYPICAL.age_gy * 1.0e9, rng, companion="star",
                                         max_axis_au=cluster_binaries.axis_limit_au(TYPICAL, 16.0))
    assert scatter_companions.phenomenon_tree(row, TYPICAL) == (by_hand if by_hand.nodes else None)


# --- The rows per system -----------------------------------------------------------------------------------------


def test_a_cluster_has_the_rows_per_system_its_own_odds_give():
    for kind, subtype in (("black-hole", "stellar"), ("neutron-star", None)):
        ratio = phenomenon_scatter.cluster_rate_ratio(kind, subtype, TYPICAL)
        assert ratio == pytest.approx(phenomenon_scatter.row_density_divisor(kind, subtype) / phenomenon_scatter.row_density_divisor(kind, subtype, TYPICAL))
        assert 0.97 < ratio <= 1.0 + 1e-9               # a few percent of the class are the companions of a row's placed object
    pair = phenomenon_scatter.cluster_rate_ratio("black-hole", "stellar-pair", COMPACT)
    assert pair == pytest.approx(phenomenon_scatter.pair_odds("black-hole", "stellar-pair", COMPACT) / phenomenon_scatter.pair_odds("black-hole", "stellar-pair")
                                 * phenomenon_scatter.row_density_divisor("black-hole", "stellar-pair") / phenomenon_scatter.row_density_divisor("black-hole", "stellar-pair", COMPACT))
    assert pair > 1.0 > phenomenon_scatter.cluster_rate_ratio("black-hole", "stellar-pair", SPARSE)      # a dense cluster has more black hole pairs than the field, a sparse one fewer
    for kind, subtype in (("planetary-nebula", None), ("supernova-remnant", None), ("black-hole", "intermediate")):
        assert phenomenon_scatter.cluster_rate_ratio(kind, subtype, TYPICAL) == 1.0


def test_the_components_pass_gives_a_cluster_its_ratio_and_a_satellite_the_fields():
    cluster = dataclasses.replace(_component("cluster"), mass_solar=COMPACT.mass_solar, half_mass_radius_pc=COMPACT.half_mass_radius_pc)
    satellite = _component("satellite")
    rate = 1.0e-6
    ratio = phenomenon_scatter.cluster_rate_ratio("black-hole", "stellar-pair", cluster)
    assert ratio > 1.0
    assert component_pass._class_mean(cluster, "black-hole", "stellar-pair", rate, None) == pytest.approx(component_pass.systems_of(cluster) * rate * ratio)
    assert component_pass._class_mean(satellite, "black-hole", "stellar-pair", rate, None) == pytest.approx(component_pass.systems_of(satellite) * rate)


def _component(kind):
    return galaxyComponents.Component(kind=kind, name=f"{kind} 1", source="generated", center_pc=(0.0, 0.0, 0.0), velocity_kms=(0.0, 0.0, 0.0), epoch_unix=0.0,
                                      mass_solar=2.0e5, half_mass_radius_pc=5.0, profile="king", age_gy=11.0, metallicity_fe_h=-1.5, seed=7, king_w0=7,
                                      tidal_radius_pc=60.0, id=1)


# --- The writers --------------------------------------------------------------------------------------------------


def test_the_scatter_writer_reads_the_rows_cluster_and_draws_its_companion_at_the_clusters_odds(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        cluster = dataclasses.replace(_component("cluster"), mass_solar=TYPICAL.mass_solar, half_mass_radius_pc=TYPICAL.half_mass_radius_pc, age_gy=TYPICAL.age_gy)
        store.replace_galaxy_components(conn, [cluster], kind="cluster")
        conn.commit()
        stored = store.get_galaxy_components(conn, kind="cluster")[0]
        seeds = _seeds(400, start=401)
        field_rows = [_phenomenon_row("black-hole", "stellar", s, mass=16.0) for s in seeds]
        cluster_rows = [tuple(stored.id if name == "cluster_id" else value for name, value in zip(store.PHENOMENON_SCATTER_COLUMNS, row)) for row in field_rows]
        by_field = {index for index, _ in scatter_companions.phenomenon_trees(field_rows, conn)}
        by_cluster = {index for index, _ in scatter_companions.phenomenon_trees(cluster_rows, conn)}
        assert by_field == {index for index, _ in scatter_companions.phenomenon_trees(field_rows)}      # a field row ignores the connection
        assert len(by_cluster) > 3 * len(by_field)                                                       # 10 percent of the cluster's black holes hold a star, 1 percent of the field's
        assert by_cluster == {index for index, _ in scatter_companions.phenomenon_trees(cluster_rows, conn)}
        assert {index for index, _ in scatter_companions.phenomenon_trees(cluster_rows)} == by_field    # without the connection the cluster is unknown
        assert PHEN["cluster_id"] == store.PHENOMENON_SCATTER_COLUMNS.index("cluster_id")
    finally:
        conn.close()
