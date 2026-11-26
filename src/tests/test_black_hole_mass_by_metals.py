# tests/test_black_hole_mass_by_metals.py

"""GEN.554: a stellar-mass black hole's mass follows the metals (stand-in BHM1, L2). The registered law (uniform over 5 to 20 solar
masses) is the anchor, the holes of the disk around the sun; [R1466]'s delayed engine over the stars that fail sets how much the span
grows as the metals fall, held at 50 solar masses ([R5365]). Counts, the scatter's cut and the neutron stars do not move."""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import remnant_density
from planetgen.physics import fryer_remnants, remnant_mass
from planetgen.util import draw
from tests.test_phenomenon_scatter import SHAPE


def test_the_fits_are_the_papers_equations_at_the_sun_and_at_a_tenth_of_it():
    # Eq. 5 at 24 solar masses, solar metals: 1.1 + 0.2 e^(13/4) - 3 e^(-0.8).
    assert fryer_remnants.delayed_baryonic_mass_sol(24.0, 1.0) == pytest.approx(1.1 + 0.2 * math.exp(3.25) - 3.0 * math.exp(-0.8))
    # Eq. 7 at 40 solar masses: the lesser of 33.35 + 6 x 6 and 40 - (52 - 18.35) at the sun, and with a tenth of the metals.
    assert fryer_remnants.delayed_baryonic_mass_sol(40.0, 1.0) == pytest.approx(40.0 - (1.3 * 40.0 - 18.35))
    assert fryer_remnants.delayed_baryonic_mass_sol(40.0, 0.1) == pytest.approx(40.0 - math.sqrt(0.1) * (1.3 * 40.0 - 18.35))
    # Eq. 9 takes over from the wind-stripped Eq. 7 above 50: 1.8 + log10(100 - 89) at 100 solar masses.
    assert fryer_remnants.delayed_baryonic_mass_sol(100.0, 1.0) == pytest.approx(1.8 + math.log10(11.0))
    assert fryer_remnants.delayed_black_hole_mass_sol(40.0, 1.0) == pytest.approx(0.9 * 6.35)


def test_the_fits_join_across_30_solar_masses_and_the_holes_grow_as_the_metals_fall():
    for z in (1.0, 0.3, 0.05):
        assert fryer_remnants.delayed_baryonic_mass_sol(29.999, z) == pytest.approx(fryer_remnants.delayed_baryonic_mass_sol(30.001, z), rel=1e-2)   # the paper's two fits meet within a percent
    for mass in (24.0, 35.0, 45.0, 70.0):
        assert fryer_remnants.delayed_baryonic_mass_sol(mass, 0.1) > fryer_remnants.delayed_baryonic_mass_sol(mass, 1.0)


def test_the_trend_is_one_for_the_disk_around_the_sun_and_rises_as_the_metals_fall():
    fe_h = np.array([0.4, 0.0, -0.2, -0.5, -1.0, -2.0, -3.5])
    trend = fryer_remnants.mass_trend(fe_h)
    assert np.all(np.diff(trend) >= 0.0)
    assert trend[0] == trend[1]                                   # the fit stops at the solar metals
    assert trend[0] < 1.0 < trend[3]
    assert fryer_remnants.mass_trend(-0.5) == pytest.approx(float(trend[3]))
    assert 5.0 < trend[-1] < 9.0                                  # 6 times the disk's at no metals


def test_the_stretch_holds_the_span_to_the_top_and_leaves_the_anchor_alone_at_the_disks_metals():
    low, high = tuning.BLACK_HOLE_MASS_RANGE_SOLAR
    stretch = remnant_mass.black_hole_stretch(np.array([0.4, 0.0, -0.2, -3.5]))
    assert stretch[0] == stretch[1] < 1.0
    assert 1.0 < stretch[2] < 3.0
    assert low + (high - low) * stretch[3] == pytest.approx(tuning.BLACK_HOLE_MASS_TOP_SOLAR)
    assert remnant_mass.stretched_mass_sol(low, 7.0) == low        # the lower edge stays
    assert remnant_mass.stretched_mass_sol(high, 1.0) == high


def test_a_hole_drawn_at_a_metal_poor_place_keeps_its_rank_in_the_law_and_none_below_the_edge():
    stretch = remnant_mass.black_hole_stretch(-1.0)
    anchor = remnant_mass.draw_mass_sol("black-hole", rng=draw.Stream("stretch"))
    stretched = remnant_mass.draw_mass_sol("black-hole", rng=draw.Stream("stretch"), stretch=stretch)
    assert stretched == pytest.approx(5.0 + (anchor - 5.0) * stretch)
    cut = remnant_mass.draw_mass_sol("black-hole", (14.0, 20.0), draw.Stream("stretch"), stretch=stretch)
    assert 5.0 + 9.0 * stretch <= cut <= 5.0 + 15.0 * stretch      # the cut is in the anchor law's masses
    masses = [remnant_mass.draw_mass_sol("black-hole", rng=draw.Stream(f"m{i}"), stretch=stretch) for i in range(400)]
    assert 5.0 <= min(masses) and max(masses) <= tuning.BLACK_HOLE_MASS_TOP_SOLAR
    assert max(masses) > 40.0                                       # a tenth of the metals reaches the heavy holes
    plain = [remnant_mass.draw_mass_sol("black-hole", rng=draw.Stream(f"m{i}")) for i in range(400)]
    assert plain == [remnant_mass.draw_mass_sol("black-hole", rng=draw.Stream(f"m{i}"), stretch=None) for i in range(400)]


def test_the_neutron_stars_and_the_share_above_the_cut_do_not_follow_the_metals():
    assert remnant_mass.share_above("black-hole", 14.0) == pytest.approx(0.4)
    assert remnant_mass.draw_mass_sol("neutron-star", rng=draw.Stream("ns")) == remnant_mass.draw_mass_sol("neutron-star", rng=draw.Stream("ns"), stretch=None)


def test_the_stretch_of_a_place_is_the_densities_weighed_mean_of_its_populations():
    inner = remnant_density.place_stretch((1000.0, 0.0, 0.0), SHAPE)
    outer = remnant_density.place_stretch((9000.0, 0.0, 0.0), SHAPE)
    assert outer > inner > 0.0                                     # the disk's [Fe/H] falls outward, the holes grow
    array = remnant_density.stretch_array(np.array([1000.0, 9000.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]), SHAPE)
    assert array.tolist() == pytest.approx([inner, outer])
    stretches = remnant_density.population_stretches(30000.0)
    assert stretches["old"] > stretches["young"]                   # the old disk is the poorer in metals


def test_a_black_hole_of_a_star_system_takes_the_stretch_of_its_own_metals():
    from planetgen.generation import multistar

    def holes(fe_h):
        found = [multistar._remnant("black-hole", 8.0e9, 0.0, draw.Stream(f"member{fe_h}:{i}"), fe_h=fe_h) for i in range(300)]
        return np.array([leaf.mass for leaf in found if leaf is not None])
    rich, poor = holes(0.0), holes(-1.0)
    assert len(rich) > 100 and len(poor) > 100
    assert 5.0 <= rich.min() and rich.max() <= 5.0 + 15.0 * remnant_mass.black_hole_stretch(0.0) + 1e-9
    assert poor.max() > 30.0 and poor.mean() > 1.8 * rich.mean()       # three times the span, a quarter of it above the edge
    assert poor.max() <= tuning.BLACK_HOLE_MASS_TOP_SOLAR
