# tests/test_black_hole_density.py

"""
GEN.563: the stellar-mass black holes follow the stars' metallicity (docs/design/black-hole-density.md). A galaxy draws
the share of its massive stars that fail to explode at solar and at zero metallicity (cut at the neutron star budget),
the share falls as the metals rise, a place's black holes per star are the populations' share of it, and every scatter
class draws from a stream of its own.
"""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import black_hole_density, density, ejecta_field, profile, structure_draw
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.generation import object_first, phenomenon_reach as reach, phenomenon_scatter as scatter
from planetgen.physics import black_hole_formation as formation
from planetgen.physics.units import ly_to_pc
from tests import test_gen_scatter_tiers as tiers

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
E_VALUE = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
SUN = (8128.0, 0.0, 0.0)
OLD_DENSITY_PC3 = 7.0e-5
"""The constant GEN.563 replaced: 1e8 black holes in the Galaxy at the reference density."""


def _shape(solar, free):
    return SHAPE._replace(structure=SHAPE.structure._replace(black_hole_solar_share=solar, black_hole_metal_free_share=free))


# --- The formation laws -----------------------------------------------------------------------------------------


def test_the_massive_stars_are_six_tenths_of_a_percent_of_the_stars_born():
    assert formation.imf_massive_share() == pytest.approx(0.00637, rel=0.01)


def test_a_young_window_has_made_four_fifths_of_its_remnants_and_an_older_one_all():
    assert formation.remnant_formed_share(0.0, 0.1) == pytest.approx(0.79, abs=0.03)
    for window in ((0.1, 3.0), (3.0, 10.0), (10.0, 13.0)):
        assert formation.remnant_formed_share(*window) == 1.0
    assert formation.remnant_formed_share(0.0, 0.0) == 0.0


def test_the_failed_share_is_linear_in_the_metals_and_held_at_the_solar_value_above_solar():
    f = formation.failed_share
    assert f(0.0, 0.1, 0.3) == pytest.approx(0.1) and f(1.0, 0.1, 0.3) == pytest.approx(0.1)
    assert f(-3.5, 0.1, 0.3) == pytest.approx(0.3, abs=1e-3)
    assert f(math.log10(0.5), 0.1, 0.3) == pytest.approx(0.2)           # half the solar metals: halfway
    assert f(np.array([-1.0, 0.0, 0.5]), 0.1, 0.3).tolist() == pytest.approx([0.28, 0.1, 0.1])


def test_the_neutron_star_budget_caps_the_share_between_12_and_25_percent():
    assert formation.share_cap(0.4) == pytest.approx(0.121, abs=0.002)
    assert formation.share_cap(0.2) == pytest.approx(0.247, abs=0.002)
    assert formation.share_cap(0.2) > formation.share_cap(0.3) > formation.share_cap(0.4)


def test_the_draws_keep_the_shares_in_order_and_inside_the_cap():
    for solar_u in (0.0, 0.3, 1.0):
        for free_u in (0.0, 0.5, 1.0):
            solar, free = formation.draw_failed_shares(solar_u, free_u, 0.15)
            assert 0.0 <= solar <= 0.15 and solar <= free <= 0.15 + 1e-12


def test_a_seed_draws_its_shares_in_range_and_the_same_every_time():
    """T7: over a thousand seeds."""
    low = tuning.BLACK_HOLE_FAILED_SHARE_METAL_FREE_RANGE
    caps = []
    for i in range(1000):
        seed = i.to_bytes(16, "big")
        solar, free = structure_draw.draw_black_hole_shares(seed)
        cap = formation.share_cap(ejecta_field.unbound_fraction("neutron-star", seed))
        caps.append(cap)
        assert 0.0 <= solar <= min(0.45, cap) + 1e-12
        assert solar <= free <= max(cap, solar) + 1e-12
        assert free <= max(low[1], solar) + 1e-12
        assert (solar, free) == structure_draw.draw_black_hole_shares(seed)
    assert 0.121 <= min(caps) and max(caps) <= 0.247 + 1e-9
    assert len({structure_draw.draw_black_hole_shares(i.to_bytes(16, "big")) for i in range(50)}) == 50


def test_the_structure_carries_the_seeds_shares():
    seed = bytes(range(16))
    structure = structure_draw.draw_structure(seed)
    assert (structure.black_hole_solar_share, structure.black_hole_metal_free_share) == structure_draw.draw_black_hole_shares(seed)


# --- A place's black holes -----------------------------------------------------------------------------------------


def test_the_failed_share_is_monotone_in_the_metals_and_the_halo_makes_at_least_the_disks():
    """T3: for any draw, a metal-poor place makes at least as many as a metal-rich one."""
    for solar, free in ((0.0, 0.15), (0.2, 0.2), (0.05, 0.3), (0.12, 0.12)):
        shape = _shape(solar, free)
        halo = black_hole_density.place_factor((15000.0, 0.0, 15000.0), shape)
        disk = black_hole_density.place_factor(SUN, shape)
        assert halo >= disk
        low = formation.failed_share(-3.5, solar, free)
        assert low >= formation.failed_share(0.0, solar, free)
    assert black_hole_density.place_factor(SUN, _shape(0.0, 0.0)) == 0.0


def test_the_arrays_and_the_one_point_form_agree_and_the_profile_reads_the_same_factor():
    points = [SUN, (3000.0, 0.0, 0.0), (15000.0, 0.0, 15000.0), (0.0, 0.0, 500.0), (400.0, 100.0, 20.0)]
    x, y, z = (np.array([p[axis] for p in points]) for axis in range(3))
    weights = black_hole_density.weights_array(x, y, z, SHAPE)
    total = density.total_density_array(x, y, z, SHAPE)
    for index, point in enumerate(points):
        assert weights[index] / total[index] == pytest.approx(black_hole_density.place_factor(point, SHAPE), rel=1e-6)
        found = profile.profile_at(point, SHAPE, EDGE_PC, E_VALUE)
        assert found.kind_factors["black-hole"] == pytest.approx(black_hole_density.place_factor(point, SHAPE), rel=0.02)


def test_the_old_constant_is_inside_the_new_range_and_the_new_rates_are_reported(capsys):
    """T4: 7e-5 per pc3 over the reference density is a per-star 5.0e-4; the new solar-circle rate is a multiple."""
    old = OLD_DENSITY_PC3 / tuning.REFERENCE_STELLAR_DENSITY_PC3
    lowest = black_hole_density.base_rate_per_star() * black_hole_density.place_factor(SUN, _shape(0.0, 0.15))
    highest = black_hole_density.base_rate_per_star() * black_hole_density.place_factor(SUN, _shape(0.247, 0.247))
    assert lowest < old < highest
    middle = black_hole_density.base_rate_per_star() * black_hole_density.place_factor(SUN, SHAPE)
    print(f"black holes per star at the solar circle: {middle:.2e} ({middle / old:.2f} times the old {old:.1e})")


def test_the_intermediate_mass_holes_keep_their_rate_whatever_the_metals():
    assert scatter.class_rate_per_star("black-hole", "intermediate") == pytest.approx(1.5e-9)
    assert scatter.class_rate_per_star("black-hole", "stellar") == pytest.approx(black_hole_density.base_rate_per_star())


def test_the_below_cut_draw_scales_with_the_place_factor():
    centre = (4000.0, 0.0, 0.0)
    cells = [(shell, band, 0) for shell in range(2, 6) for band in range(-2, 3)]

    def total(shape):
        return sum(count for cell in cells for kind, subtype, _range, count in
                   scatter.below_cut_draws(cell, centre, shape, 5000.0, 20.0, 7)[0] if kind == "black-hole")
    poor, rich = total(_shape(0.0, 0.05)), total(_shape(0.15, 0.15))
    assert rich > 1.5 * poor
    assert total(_shape(0.0, 0.0)) == 0


# --- The galaxy's totals ------------------------------------------------------------------------------------------


def test_bound_plus_ejected_black_holes_are_what_the_galaxy_formed():
    """T5."""
    systems = ejecta_field.galaxy_systems(SHAPE, E_VALUE, EDGE_PC)
    formed = black_hole_density.formed_per_star(SHAPE) * systems
    bound = black_hole_density.base_rate_per_star() * black_hole_density.mean_factor(SHAPE) * systems
    ejected = ejecta_field.unbound_count("black-hole", SHAPE, E_VALUE, EDGE_PC)
    assert ejected == pytest.approx(tuning.BLACK_HOLE_UNBOUND_FRACTION * formed)
    assert bound + ejected == pytest.approx(formed)


def test_neutron_stars_and_black_holes_formed_fit_among_the_massive_stars_in_every_galaxy():
    """T9: over a seed sweep, and the black holes stay within 0.5 to 3.1 times the old count."""
    old = OLD_DENSITY_PC3 / tuning.REFERENCE_STELLAR_DENSITY_PC3
    budget = formation.imf_massive_share()
    ratios = []
    for i in range(40):
        seed = (i * 7919).to_bytes(16, "big")
        structure = structure_draw.draw_structure(seed)
        shape = SHAPE._replace(structure=structure)
        unbound = ejecta_field.unbound_fraction("neutron-star", seed)
        neutron = tuning.phenomenon_rate_per_star("neutron-star") * (1.0 + unbound)
        holes = black_hole_density.formed_per_star(shape) * (1.0 - budget)       # per star born
        assert neutron + holes <= budget * (1.0 + 1e-9), (i, neutron, holes, budget)
        ratios.append(black_hole_density.base_rate_per_star() * black_hole_density.mean_factor(shape) / old)
    assert min(ratios) > 0.3 and max(ratios) < 3.4


# --- The scatter ---------------------------------------------------------------------------------------------------


def test_each_class_draws_from_a_stream_of_its_own_so_the_black_holes_move_no_other_row(monkeypatch):
    """T2: more black holes change none of the neutron stars', planetary nebulae's or remnants' rows."""
    def rows(shell=4):
        return list(scatter.scatter_shell(tiers.SHAPE, shell, 3, tiers.EDGE_PC, 4000.0, 11, min_mass_solar=1.0))
    before = rows()
    base = black_hole_density.base_rate_per_star()
    monkeypatch.setattr(black_hole_density, "base_rate_per_star", lambda: 3.0 * base)
    after = rows()

    def others(found):
        return [row for row in found if row[3] != "black-hole"]
    assert others(before) == others(after) and others(before)
    holes_before = sum(1 for row in before if row[3] == "black-hole" and row[4] == "stellar")
    holes_after = sum(1 for row in after if row[3] == "black-hole" and row[4] == "stellar")
    assert holes_after > 1.5 * holes_before > 0


def test_the_reach_bound_holds_for_the_lowest_and_the_highest_draw():
    """T6: no weight exceeds its bound whatever the galaxy's shares."""
    rng = np.random.default_rng(3)
    for solar, free in ((0.0, 0.15), (0.247, 0.247), (0.05, 0.32)):
        shape = tiers.SHAPE._replace(structure=tiers.SHAPE.structure._replace(
            black_hole_solar_share=solar, black_hole_metal_free_share=free))
        for shell, high_band in [(5, -1), (60, 10), (800, 60)]:
            groups = object_first.ReachGroups.for_shell(shell, high_band, tiers.EDGE_PC)
            bound = reach.class_bounds(shape, groups, "black-hole", "stellar")
            for index in range(len(groups)):
                radius = np.cbrt(rng.uniform(groups.r_in[index] ** 3, groups.r_out[index] ** 3, 200))
                cos_colat = rng.uniform(groups.cos_low[index], groups.cos_high[index], 200)
                lon = rng.uniform(0.0, 2.0 * math.pi, 200)
                plane = radius * np.sqrt(1.0 - cos_colat ** 2)
                weight = reach._weights("black-hole", "stellar", plane * np.cos(lon), plane * np.sin(lon),
                                        radius * cos_colat, shape)
                assert np.all(weight <= bound[index] * (1.0 + 1e-9)), (solar, free, shell, index)
            rows = list(reach.scatter_shell(shape, shell, high_band, tiers.EDGE_PC, 4000.0, 5, min_mass_solar=1.0,
                                            check=True))
            assert rows
