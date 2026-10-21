# tests/test_shape_draw.py

"""GEN.401: a spiral galaxy's shape (size, height, bulge share, bar length, arms) drawn from its seed by the laws of
`galaxy/shape_draw.py`: the same seed gives the same shape, and over many seeds the draws follow the stated laws."""

import argparse
import math
import statistics

import pytest

from planetgen import tuning
from planetgen.cli.generate import add_plan_arguments
from planetgen.galaxy import density, mass_chain, shape_draw


def _seed(number):
    return number.to_bytes(16, "big")


SEEDS = range(1, 2001)


@pytest.fixture(scope="module")
def draws():
    return [shape_draw.draw_galaxy(_seed(number)) for number in SEEDS]


def _spread(values):
    return statistics.mean(values), statistics.stdev(values)


def test_the_same_seed_gives_the_same_shape_and_two_seeds_differ():
    assert shape_draw.draw_galaxy(_seed(5)) == shape_draw.draw_galaxy(_seed(5))
    assert shape_draw.draw_galaxy(_seed(5)) != shape_draw.draw_galaxy(_seed(6))
    assert shape_draw.draw_galaxy(None) == shape_draw.draw_galaxy(mass_chain.FALLBACK_SEED)


def test_only_the_spiral_family_and_the_six_options_exist():
    with pytest.raises(ValueError, match="elliptical"):
        shape_draw.draw_galaxy(_seed(1), family="elliptical")
    with pytest.raises(ValueError, match="not a shape option"):
        shape_draw.draw_galaxy(_seed(1), bar_angle_deg=27.0)


def test_the_disk_length_follows_the_size_mass_relation_of_shen_2003(draws):
    """`ln(R50 / median)` over its mass-dependent width is a standard normal; and 3.9 kpc at 5e10 solar masses."""
    assert shape_draw.size_median_kpc(5e10) == pytest.approx(3.9, abs=0.05)
    assert shape_draw.size_sigma_ln(1e8) == pytest.approx(0.47, abs=0.01)
    assert shape_draw.size_sigma_ln(1e12) == pytest.approx(0.34, abs=0.01)
    scores = []
    for d in draws:
        mass = d.chain.stellar_mass_solar
        r50 = d.disk_scale_length_pc * shape_draw.EXPONENTIAL_HALF_LIGHT_OVER_SCALE / 1000.0
        scores.append(math.log(r50 / shape_draw.size_median_kpc(mass)) / shape_draw.size_sigma_ln(mass))
    mean, spread = _spread(scores)
    assert abs(mean) < 4.0 / math.sqrt(len(scores))
    assert spread == pytest.approx(1.0, rel=0.08)


def test_heavier_galaxies_have_longer_disks(draws):
    ordered = sorted(draws, key=lambda d: d.chain.stellar_mass_solar)
    third = len(ordered) // 3
    light, heavy = ordered[:third], ordered[-third:]
    assert (statistics.median(d.disk_scale_length_pc for d in heavy)
            > 1.3 * statistics.median(d.disk_scale_length_pc for d in light))


def test_the_height_is_the_length_over_kregels_ratio(draws):
    mean, spread = _spread([math.log10(d.disk_scale_length_pc / d.disk_scale_height_pc) for d in draws])
    assert mean == pytest.approx(math.log10(shape_draw.KREGEL_Q), abs=4.0 * 0.13 / math.sqrt(len(draws)))
    assert spread == pytest.approx(shape_draw.KREGEL_SCATTER_DEX, rel=0.08)


def test_the_late_share_falls_with_mass_as_kelvins_table_does():
    shares = [shape_draw.late_share(10.0 ** m) for m in (8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 12.0)]
    assert shares[0] == shares[1] == pytest.approx(0.28 / 0.46)  # held below the first bound
    assert shares[-1] == shares[-2] == pytest.approx(0.06 / 0.28)  # held above the last
    assert all(a >= b for a, b in zip(shares, shares[1:]))
    assert shape_draw.late_share(10.0 ** 9.75) == pytest.approx(0.5 * (0.36 / 0.64 + 0.32 / 0.74))


def test_the_late_stage_is_drawn_in_the_tables_share(draws):
    for low, high in ((9.4, 9.9), (10.2, 10.6), (10.6, 11.4)):
        chosen = [d for d in draws if low <= math.log10(d.chain.stellar_mass_solar) < high]
        if len(chosen) < 100:
            continue
        expected = statistics.mean(shape_draw.late_share(d.chain.stellar_mass_solar) for d in chosen)
        late = sum(d.stage == "late" for d in chosen) / len(chosen)
        assert abs(late - expected) < 4.0 * math.sqrt(expected * (1.0 - expected) / len(chosen)), (low, high)


def test_the_bulge_share_follows_graham_and_worleys_medians_by_stage(draws):
    for stage, median, sigma in (("early", shape_draw.LOG_BD_EARLY, shape_draw.LOG_BD_EARLY_SIGMA),
                                 ("late", shape_draw.LOG_BD_LATE, shape_draw.LOG_BD_LATE_SIGMA)):
        logs = [math.log10(d.bulge_to_total / (1.0 - d.bulge_to_total)) for d in draws if d.stage == stage]
        mean, spread = _spread(logs)
        assert mean == pytest.approx(median, abs=4.0 * sigma / math.sqrt(len(logs)))
        assert spread == pytest.approx(sigma, rel=0.1)
    early = statistics.median(d.bulge_to_total for d in draws if d.stage == "early")
    late = statistics.median(d.bulge_to_total for d in draws if d.stage == "late")
    assert early == pytest.approx(0.24, abs=0.03) and late == pytest.approx(0.06, abs=0.02)


def test_the_bulge_amplitude_gives_the_drawn_bulge_to_total_in_the_models_own_masses():
    for number in (3, 11, 40, 77):
        d = shape_draw.draw_galaxy(_seed(number))
        shape = density.build_galaxy_shape(
            d.disk_scale_length_pc, d.disk_scale_height_pc, d.bulge_scale_radius_pc, d.bulge_amplitude,
            d.arm_count, math.radians(d.pitch_angle_deg), 0.4, structure=d.structure)
        masses = density.component_masses(shape)
        total = masses["bulge"] + masses["thin_disk"] + masses["thick_disk"]
        assert masses["bulge"] / total == pytest.approx(d.bulge_to_total, rel=0.02)


def test_the_bar_is_between_a_fifth_and_three_disk_scale_lengths(draws):
    for d in draws:
        ratio = d.bar_length_pc / d.disk_scale_length_pc
        assert shape_draw.BAR_LENGTH_OVER_H_RANGE[0] <= ratio <= shape_draw.BAR_LENGTH_OVER_H_RANGE[1]
        assert d.bulge_scale_radius_pc == pytest.approx(d.bar_length_pc / shape_draw.BAR_LENGTH_PER_SCALE_RADIUS)
    early = statistics.mean(d.bar_length_pc / d.disk_scale_length_pc for d in draws if d.stage == "early")
    late = statistics.mean(d.bar_length_pc / d.disk_scale_length_pc for d in draws if d.stage == "late")
    assert early > 1.5 * late  # the 0.2 floor lifts the late mean from 0.57


def test_the_arm_counts_follow_harts_shares(draws):
    for arms, percent in shape_draw.ARM_COUNT_SHARES:
        share = sum(d.arm_count == arms for d in draws) / len(draws)
        expected = percent / 100.0
        assert abs(share - expected) < 4.0 * math.sqrt(expected * (1.0 - expected) / len(draws)) + 0.005, arms


def test_the_pitch_angle_falls_with_mass_inside_harts_range(draws):
    assert shape_draw.pitch_median_deg(10.0 ** 9.6) == pytest.approx(26.0)
    assert shape_draw.pitch_median_deg(10.0 ** 11.1) == pytest.approx(13.8)
    low, high = shape_draw.PITCH_RANGE_DEG
    assert all(low <= d.pitch_angle_deg <= high for d in draws)
    ordered = sorted(draws, key=lambda d: d.chain.stellar_mass_solar)
    third = len(ordered) // 3
    assert (statistics.mean(d.pitch_angle_deg for d in ordered[:third])
            > statistics.mean(d.pitch_angle_deg for d in ordered[-third:]) + 2.0)


def test_an_option_given_replaces_its_draw_and_the_bulge_follows_the_disk_used():
    base = shape_draw.draw_galaxy(_seed(9))
    longer = shape_draw.draw_galaxy(_seed(9), disk_scale_length_pc=2.0 * base.disk_scale_length_pc)
    assert longer.disk_scale_length_pc == 2.0 * base.disk_scale_length_pc
    assert longer.bulge_to_total == base.bulge_to_total and longer.arm_count == base.arm_count
    assert longer.bulge_scale_radius_pc == pytest.approx(2.0 * base.bulge_scale_radius_pc)  # the bar scales with the disk
    assert longer.disk_scale_height_pc == pytest.approx(2.0 * base.disk_scale_height_pc)  # so does the height
    fixed = shape_draw.draw_galaxy(_seed(9), arm_count=4, pitch_angle_deg=12.0)
    assert (fixed.arm_count, fixed.pitch_angle_deg) == (4, 12.0)
    assert shape_draw.draw_galaxy(_seed(9), arm_count=None) == base


# --- The Structure dial (GEN.604): flatness, bulge share and bulge size are ratios ---------------------------------

def _drawn_flatness(d):
    return d.disk_scale_length_pc / d.disk_scale_height_pc


def test_flatness_alone_sets_the_height_as_a_ratio_of_the_length_and_nothing_else_of_the_disk():
    base = shape_draw.draw_galaxy(_seed(9))
    flat = shape_draw.draw_galaxy(_seed(9), flatness=12.0)
    assert flat.disk_scale_height_pc == pytest.approx(flat.disk_scale_length_pc / 12.0)
    assert flat.disk_scale_length_pc == base.disk_scale_length_pc
    assert (flat.bulge_to_total, flat.bulge_scale_radius_pc, flat.arm_count, flat.pitch_angle_deg, flat.stage) == (
        base.bulge_to_total, base.bulge_scale_radius_pc, base.arm_count, base.pitch_angle_deg, base.stage)
    # the same share of a different disk is a different bulge density: the amplitude is recomputed
    assert flat.bulge_amplitude != base.bulge_amplitude


def test_a_pinned_flatness_keeps_its_ratio_when_the_length_moves():
    for length in (500.0, 2600.0, 9000.0):
        d = shape_draw.draw_galaxy(_seed(9), disk_scale_length_pc=length, flatness=8.0)
        assert _drawn_flatness(d) == pytest.approx(8.0)


def test_the_drawn_flatness_residual_is_the_same_at_a_new_length():
    base = shape_draw.draw_galaxy(_seed(9))
    longer = shape_draw.draw_galaxy(_seed(9), disk_scale_length_pc=3.0 * base.disk_scale_length_pc)
    assert _drawn_flatness(longer) == pytest.approx(_drawn_flatness(base))


def test_bulge_share_alone_sets_the_share_and_solves_the_amplitude_for_it():
    base = shape_draw.draw_galaxy(_seed(9))
    for share in (0.05, 0.4, 0.8):
        d = shape_draw.draw_galaxy(_seed(9), bulge_to_total=share)
        assert d.bulge_to_total == share
        assert (d.disk_scale_length_pc, d.disk_scale_height_pc, d.bulge_scale_radius_pc) == (
            base.disk_scale_length_pc, base.disk_scale_height_pc, base.bulge_scale_radius_pc)
        shape = density.build_galaxy_shape(
            d.disk_scale_length_pc, d.disk_scale_height_pc, d.bulge_scale_radius_pc, d.bulge_amplitude,
            d.arm_count, math.radians(d.pitch_angle_deg), 0.4, structure=d.structure)
        masses = density.component_masses(shape)
        total = masses["bulge"] + masses["thin_disk"] + masses["thick_disk"]
        assert masses["bulge"] / total == pytest.approx(share, rel=0.02)


def test_a_bulge_share_of_zero_is_a_disc_with_no_bulge_density():
    d = shape_draw.draw_galaxy(_seed(9), bulge_to_total=0.0)
    assert d.bulge_amplitude == 0.0
    shape = density.build_galaxy_shape(
        d.disk_scale_length_pc, d.disk_scale_height_pc, d.bulge_scale_radius_pc, d.bulge_amplitude,
        d.arm_count, math.radians(d.pitch_angle_deg), 0.4, structure=d.structure)
    assert density.component_masses(shape)["bulge"] == 0.0


def test_bulge_size_alone_sets_the_radius_as_a_ratio_of_the_length_and_the_bar_with_it():
    base = shape_draw.draw_galaxy(_seed(9))
    for ratio in (0.05, 0.3, 1.0):
        d = shape_draw.draw_galaxy(_seed(9), bulge_size=ratio)
        assert d.bulge_scale_radius_pc == pytest.approx(ratio * d.disk_scale_length_pc)
        assert d.bar_length_pc == pytest.approx(shape_draw.BAR_LENGTH_PER_SCALE_RADIUS * d.bulge_scale_radius_pc)
        assert (d.bulge_to_total, d.disk_scale_height_pc) == (base.bulge_to_total, base.disk_scale_height_pc)
    longer = shape_draw.draw_galaxy(_seed(9), disk_scale_length_pc=2.0 * base.disk_scale_length_pc, bulge_size=0.3)
    assert longer.bulge_scale_radius_pc == pytest.approx(0.6 * base.disk_scale_length_pc)


def test_the_retired_absolute_structure_options_are_not_options():
    for name in ("disk_scale_height_pc", "bulge_scale_radius_pc", "bulge_amplitude"):
        with pytest.raises(ValueError):
            shape_draw.draw_galaxy(_seed(9), **{name: 10.0})


def test_a_plan_with_no_shape_options_draws_every_one_from_its_seed():
    parser = argparse.ArgumentParser()
    add_plan_arguments(parser)
    args = parser.parse_args([])
    for option in ("disk_scale_length", "flatness", "bulge_share", "bulge_size",
                   "arm_count", "pitch_angle_deg", "shape", "nuclear_amplitude", "arm_density", "interarm_density"):
        assert getattr(args, option) is None, option
    shape, drawn = shape_draw.plan_shape(args, _seed(21))
    assert shape.disk_scale_length_pc == drawn.disk_scale_length_pc
    assert shape.arm_count == drawn.arm_count and shape.k_norm > 0
    other, _ = shape_draw.plan_shape(args, _seed(22))
    assert other.disk_scale_length_pc != shape.disk_scale_length_pc
    explicit = parser.parse_args(["--disk-scale-length", "2000pc", "--arm-count", "3"])
    shape, _ = shape_draw.plan_shape(explicit, _seed(21))
    assert shape.disk_scale_length_pc == 2000.0 and shape.arm_count == 3


def test_no_milky_way_shape_is_left_in_the_tuning():
    for name in ("GALAXY_DISK_SCALE_LENGTH_PC_DEFAULT", "DEFAULT_DISK_SCALE_LENGTH_PC", "DEFAULT_BULGE_AMPLITUDE"):
        assert not hasattr(tuning, name)
