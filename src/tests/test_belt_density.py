# tests/test_belt_density.py

"""
GEN.1093: the density factor of a belt (`planetgen.physics.belt_density`): the tail normalised to the mass left after the listed bodies, mass and count conservation, the
bins of the size table adding to the tail, the registry units, a main belt like the Solar System's as a test target, and the inputs refused.
"""

import math

import pytest

from planetgen.physics import belt_density as bd
from planetgen.physics import constants, unit_registry as units

MAIN_MASS = 0.6e-3 * constants.EARTH_MASS_TO_KG
# the largest asteroids by diameter (km) and mass (kg): test targets, not inputs
CERES_ET_AL = ((939, 9.38e20), (525, 2.59e20), (512, 2.04e20), (434, 8.7e19), (350, 3.0e19), (300, 2.0e19), (270, 1.5e19), (250, 1.2e19), (230, 1.0e19))
LISTED_D = tuple(d for d, _ in CERES_ET_AL)
LISTED_M = tuple(m for _, m in CERES_ET_AL)


def _density(**kw):
    return bd.belt_density(2.2, 3.4, MAIN_MASS, bd.MAIN_BELT, bd.BULK_DENSITY_KG_M3, **kw)


def test_the_mass_per_square_au_is_the_belt_mass_over_its_annulus():
    result = _density()
    assert result.area_au2 == pytest.approx(math.pi * (3.4 ** 2 - 2.2 ** 2))
    assert result.mass_per_au2_kg == pytest.approx(MAIN_MASS / result.area_au2)


def test_listed_bodies_and_the_tail_add_to_the_belt_mass():
    result = _density(listed_diameters_km=LISTED_D, listed_masses_kg=LISTED_M)
    assert result.listed_mass_kg == pytest.approx(sum(LISTED_M))
    assert result.listed_mass_kg + result.tail_mass_kg == pytest.approx(MAIN_MASS, rel=1e-12)
    assert result.tail_top_km == min(LISTED_D)


def test_without_a_ledger_the_tail_is_the_whole_belt():
    result = _density()
    assert result.listed_mass_kg == 0.0 and result.tail_mass_kg == pytest.approx(MAIN_MASS, rel=1e-12)
    assert result.tail_top_km == bd.MAIN_BELT.top_km
    assert _density(largest_km=500.0).tail_top_km == 500.0


def test_the_count_above_a_diameter_falls_with_the_diameter_and_the_listed_bodies_are_exact():
    result = _density(listed_diameters_km=LISTED_D, listed_masses_kg=LISTED_M)
    objects = [c.objects for c in result.counts]
    assert objects == sorted(objects, reverse=True)
    above_100 = result.counts[2]
    assert above_100.listed == len(LISTED_D) and above_100.objects > above_100.listed
    # above the smallest listed body there is nothing but the listed bodies
    high = bd.belt_density(2.2, 3.4, MAIN_MASS, listed_diameters_km=LISTED_D, listed_masses_kg=LISTED_M, diameters_km=(300.0, 600.0))
    assert [c.objects for c in high.counts] == [5.0, 1.0] and all(c.per_au2 == c.objects / high.area_au2 for c in high.counts)


def test_the_solar_systems_belt_is_a_main_belt_within_a_few_times_its_counts():
    # Ivezic et al. 2001: 530,000 objects above 1 km; some two hundred above 100 km (test targets); the nine biggest are listed
    result = _density(listed_diameters_km=LISTED_D, listed_masses_kg=LISTED_M)
    one, ten, hundred = result.counts
    assert 530_000 / 5 < one.objects < 530_000 * 40
    assert 50 < hundred.objects < 600
    assert ten.objects > hundred.objects * 100


def test_more_mass_is_more_objects_in_proportion_and_the_shape_does_not_change():
    light = bd.belt_density(2.2, 3.4, MAIN_MASS)
    heavy = bd.belt_density(2.2, 3.4, 10.0 * MAIN_MASS)
    for a, b in zip(light.counts, heavy.counts):
        assert b.objects == pytest.approx(10.0 * a.objects, rel=1e-9)
    assert heavy.mass_per_au2_kg == pytest.approx(10.0 * light.mass_per_au2_kg)


def test_the_density_per_square_au_is_the_count_over_the_annulus():
    narrow = bd.belt_density(2.2, 3.4, MAIN_MASS)
    wide = bd.belt_density(1.2, 4.4, MAIN_MASS)
    assert narrow.counts[0].objects == pytest.approx(wide.counts[0].objects)
    assert narrow.counts[0].per_au2 / wide.counts[0].per_au2 == pytest.approx(wide.area_au2 / narrow.area_au2)


def test_a_diameter_under_the_floor_is_answered_at_the_floor_and_flagged():
    below = bd.belt_density(2.2, 3.4, MAIN_MASS, diameters_km=(0.1, bd.MAIN_BELT.floor_km))
    assert below.counts[0].below_floor and not below.counts[1].below_floor
    assert below.counts[0].objects == pytest.approx(below.counts[1].objects)


def test_a_belt_the_listed_bodies_fill_has_no_tail():
    full = bd.belt_density(2.2, 3.4, sum(LISTED_M), listed_diameters_km=LISTED_D, listed_masses_kg=LISTED_M)
    assert full.tail_mass_kg == 0.0 and [c.objects for c in full.counts] == [float(len(LISTED_D))] * 2 + [float(len(LISTED_D))]
    over = bd.belt_density(2.2, 3.4, 0.5 * sum(LISTED_M), listed_diameters_km=LISTED_D, listed_masses_kg=LISTED_M)
    assert over.tail_mass_kg == 0.0


def test_the_bins_of_the_size_table_add_to_the_tail():
    edges = (0.4, 1.0, 5.0, 10.0, 40.0, 100.0, 230.0)
    kw = dict(belt_mass_kg=MAIN_MASS, listed_diameters_km=LISTED_D, listed_masses_kg=LISTED_M)
    bins = [bd.tail_bin(lo, hi, **kw) for lo, hi in zip(edges, edges[1:])]
    tail = bd.belt_density(2.2, 3.4, MAIN_MASS, listed_diameters_km=LISTED_D, listed_masses_kg=LISTED_M)
    assert sum(mass for _, mass in bins) == pytest.approx(tail.tail_mass_kg, rel=1e-9)
    total = sum(count for count, _ in bins)
    above_floor = bd.belt_density(2.2, 3.4, MAIN_MASS, listed_diameters_km=LISTED_D, listed_masses_kg=LISTED_M, diameters_km=(0.4,)).counts[0]
    assert total == pytest.approx(above_floor.objects - above_floor.listed, rel=1e-9)
    assert all(count > 0.0 and mass > 0.0 for count, mass in bins)


def test_a_stretch_with_q_of_three_has_the_same_mass_per_decade_of_size():
    law = bd.SizeLaw("three", (1.0, 10.0, 100.0), (3.0, 3.0))
    a = bd.tail_bin(1.0, 10.0, 1e20, law)[1]
    b = bd.tail_bin(10.0, 100.0, 1e20, law)[1]
    assert a == pytest.approx(b, rel=1e-9)
    assert a + b == pytest.approx(1e20, rel=1e-9)


def test_a_law_with_q_under_three_puts_its_mass_at_the_top_and_over_three_at_the_bottom():
    shallow = bd.SizeLaw("shallow", (1.0, 100.0), (2.0,))
    steep = bd.SizeLaw("steep", (1.0, 100.0), (4.0,))
    assert bd.tail_bin(10.0, 100.0, 1e20, shallow)[1] > bd.tail_bin(1.0, 10.0, 1e20, shallow)[1]
    assert bd.tail_bin(1.0, 10.0, 1e20, steep)[1] > bd.tail_bin(10.0, 100.0, 1e20, steep)[1]


def test_the_kuiper_laws_are_a_steep_bright_side_and_a_shallow_faint_side_and_hold_their_mass():
    for law in (bd.KUIPER_HOT, bd.KUIPER_COLD):
        result = bd.belt_density(39.0, 48.0, 0.02 * constants.EARTH_MASS_TO_KG, law, 1000.0)
        assert result.tail_mass_kg == pytest.approx(0.02 * constants.EARTH_MASS_TO_KG, rel=1e-12)
        assert result.counts[0].objects > result.counts[1].objects > result.counts[2].objects > 0.0


def test_the_density_units_are_in_the_registry_and_convert():
    assert units.unit(bd.OBJECTS_UNIT).kind == "surface_number_density" and units.unit(bd.MASS_UNIT).kind == "column_mass"
    assert units.convert(1.0, "au-2", "km-2") == pytest.approx(1.0 / (constants.AU_M / 1e3) ** 2)
    assert units.convert(1.0, "kg_au2", "kg_m2") == pytest.approx(1.0 / constants.AU_M ** 2)
    assert units.convert(constants.EARTH_MASS_TO_KG, "kg_au2", "Mearth_au2") == pytest.approx(1.0)
    assert units.KINDS["surface_number_density"].dimension == units.dimension(length=-2)


def test_the_diameter_of_a_mass():
    assert bd.diameter_of(9.38e20, 2160.0) == pytest.approx(939.0, rel=0.02)
    assert bd.diameter_of(1e19) == pytest.approx((6e19 / (math.pi * 2500.0)) ** (1 / 3) / 1e3)


def test_inputs_the_model_cannot_use_are_refused():
    with pytest.raises(bd.BeltDensityError):
        bd.belt_density(3.0, 2.0, MAIN_MASS)
    with pytest.raises(bd.BeltDensityError):
        bd.belt_density(2.0, 3.0, MAIN_MASS, listed_diameters_km=(100.0,), listed_masses_kg=())
    with pytest.raises(bd.BeltDensityError):
        bd.belt_density(2.0, 3.0, MAIN_MASS, diameters_km=(0.0,))
    with pytest.raises(bd.BeltDensityError):
        bd.belt_density(2.0, 3.0, MAIN_MASS, density_kg_m3=0.0)
    with pytest.raises(bd.BeltDensityError):
        bd.SizeLaw("bad", (5.0, 1.0), (2.0,))
    with pytest.raises(bd.BeltDensityError):
        bd.SizeLaw("bad", (1.0, 5.0), (2.0, 3.0))
    with pytest.raises(bd.BeltDensityError):
        bd.tail_bin(5.0, 1.0, MAIN_MASS)
