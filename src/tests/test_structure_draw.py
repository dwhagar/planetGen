# tests/test_structure_draw.py

"""GEN.532 and GEN.533: the thick disk, the bar, the stellar halo and the populations' motions are drawn from the galaxy
seed (`galaxy.structure_draw`, `galaxy.structure`), not stored, and the density model reads them."""

import math
import re
from pathlib import Path

import pytest

from planetgen.galaxy import density, mass_chain, shape_draw, star_formation, structure as galaxyStructure
from planetgen.galaxy import structure_draw as sd


def _seed(number):
    return number.to_bytes(16, "big")


SEEDS = [_seed(n) for n in range(1, 1501)]


def _structures():
    return [(mass_chain.draw_mass_chain(s).stellar_mass_solar, sd.draw_structure(s)) for s in SEEDS]


STRUCTURES = _structures()


def test_the_same_seed_gives_the_same_structure_and_two_seeds_differ():
    assert sd.draw_structure(_seed(7)) == sd.draw_structure(_seed(7))
    assert sd.draw_structure(_seed(7)) != sd.draw_structure(_seed(8))
    assert sd.draw_structure(None) == sd.draw_structure(mass_chain.FALLBACK_SEED)


def test_the_shape_and_a_loaded_shape_have_the_same_structure():
    """No column holds the structure, so the load draws it again from the seed and must get what the plan used."""
    for number in (3, 11, 40):
        assert shape_draw.draw_galaxy(_seed(number)).structure == shape_draw.structure_for(_seed(number))


def test_the_circular_speed_follows_the_baryonic_tully_fisher_relation():
    assert sd.circular_speed_kms(5.0e10) == pytest.approx((5.0e10 / 47.0) ** 0.25)
    for mass, structure in STRUCTURES:
        assert abs(math.log10(structure.vc_kms / sd.circular_speed_kms(mass))) <= sd.SCATTER_LIMIT * sd.BTFR_V_SCATTER_DEX + 1e-12


def test_the_thick_disk_mass_ratio_sits_between_the_two_fits_and_inside_the_bound_of_the_hot_component():
    inside = 0
    for mass, structure in STRUCTURES:
        low, high = sorted(sd.thick_mass_ratio_bounds(structure.vc_kms))
        log_ratio = math.log10(structure.thick_mass_ratio)
        assert math.log10(low) - sd.SCATTER_LIMIT * sd.COMERON_MASS_RATIO_SCATTER_DEX - 1e-9 <= log_ratio
        assert log_ratio <= math.log10(high) + sd.SCATTER_LIMIT * sd.COMERON_MASS_RATIO_SCATTER_DEX + 1e-9
        inside += low <= structure.thick_mass_ratio <= high
    # The centre is spread evenly (in the log) between the fits, so about half the draws stay between them
    # once the 0.27 dex scatter is added.
    assert 0.25 < inside / len(STRUCTURES) < 0.75
    # The bound of Comeron et al. 2021 holds for the mean relation of either fit, with a margin of the hot
    # component's own 0.13 dex (galaxy by galaxy 39 percent of the real galaxies exceed it).
    for log_mass in (9.48, 10.0, 10.7, 11.3):
        mass = 10.0 ** log_mass
        vc = sd.circular_speed_kms(mass)
        for ratio in sd.thick_mass_ratio_bounds(vc):
            assert ratio <= sd.hot_component_bound(mass) * 10.0 ** sd.HOT_COMPONENT_SCATTER_DEX


def test_the_thick_disk_lengths_and_heights_follow_their_fits():
    for _mass, structure in STRUCTURES:
        median = sd.thick_length_ratio_median(structure.vc_kms)
        assert abs(math.log10(structure.thick_length_ratio / median)) <= sd.SCATTER_LIMIT * sd.COMERON_LENGTH_RATIO_SCATTER_DEX + 1e-9
        margin = sd.SCATTER_LIMIT * sd.THICK_HEIGHT_RATIO_SCATTER_DEX
        assert sd.THICK_HEIGHT_RATIO_RANGE[0] * 10.0 ** -margin <= structure.thick_height_ratio
        assert structure.thick_height_ratio <= sd.THICK_HEIGHT_RATIO_RANGE[1] * 10.0 ** margin
    heights = sorted(structure.thick_height_ratio for _, structure in STRUCTURES)
    assert 2.5 < heights[len(heights) // 2] < 4.2


def test_the_stellar_halo_follows_the_placeholder_law_with_the_dragonfly_width():
    deviations = []
    for mass, structure in STRUCTURES:
        assert sd.HALO_RATIO_RANGE[0] <= structure.halo_stellar_ratio <= sd.HALO_RATIO_RANGE[1]
        deviations.append(math.log10(structure.halo_stellar_ratio / sd.halo_ratio_median(mass)))
    mean = sum(deviations) / len(deviations)
    spread = math.sqrt(sum((d - mean) ** 2 for d in deviations) / (len(deviations) - 1))
    # Cut at the 0.1 and 30 percent shares and 2 sigma, so a little under the 1.01 dex.
    assert 0.6 < spread < sd.HALO_RATIO_SCATTER_DEX
    assert sd.halo_ratio_median(5.0e10) == pytest.approx(0.05)
    assert sd.halo_ratio_median(2.0e11) == pytest.approx(0.05 * 4.0 ** 0.79)


def test_the_bar_fraction_follows_erwin_2018_by_mass_and_every_parameter_is_inside_its_range():
    assert sd.bar_fraction(10.0 ** 9.48) == pytest.approx(0.69, abs=0.01)
    assert sd.bar_fraction(10.0 ** 10.5) == pytest.approx(0.58, abs=0.01)
    assert sd.bar_fraction(10.0 ** 11.3) == pytest.approx(0.21, abs=0.01)
    expected = sum(sd.bar_fraction(mass) for mass, _ in STRUCTURES)
    barred = sum(structure.bar_present for _, structure in STRUCTURES)
    sigma = math.sqrt(sum(sd.bar_fraction(mass) * (1.0 - sd.bar_fraction(mass)) for mass, _ in STRUCTURES))
    assert abs(barred - expected) <= 4.0 * sigma
    for _mass, structure in STRUCTURES:
        assert 0.0 <= structure.bar_angle_deg < 180.0
        assert 0.0 < structure.bar_axis_ratio_z <= structure.bar_axis_ratio_y <= 1.0
        if structure.bar_present:
            assert 1.0 - sd.BAR_ELLIPTICITY_RANGE[1] <= structure.bar_axis_ratio_y <= 1.0 - sd.BAR_ELLIPTICITY_RANGE[0]
            low = min(sd.THIN_BAR_HEIGHT_RANGE[0], sd.BOXY_PEANUT_HEIGHT_RANGE[0])
            high = max(sd.THIN_BAR_HEIGHT_RANGE[1], sd.BOXY_PEANUT_HEIGHT_RANGE[1])
            assert low <= structure.bar_axis_ratio_z <= high
        else:
            assert sd.UNBARRED_AXIS_RATIO_Y_RANGE[0] <= structure.bar_axis_ratio_y <= sd.UNBARRED_AXIS_RATIO_Y_RANGE[1]
            assert sd.UNBARRED_AXIS_RATIO_Z_RANGE[0] <= structure.bar_axis_ratio_z <= sd.UNBARRED_AXIS_RATIO_Z_RANGE[1]


def test_a_barred_bulge_is_a_bar_and_an_unbarred_one_is_rounder():
    barred = [s for _, s in STRUCTURES if s.bar_present]
    round_ = [s for _, s in STRUCTURES if not s.bar_present]
    assert len(barred) > 100 and len(round_) > 100
    assert sum(s.bar_axis_ratio_y for s in barred) / len(barred) < sum(s.bar_axis_ratio_y for s in round_) / len(round_) - 0.2
    # Peanuts take over with mass (Erwin and Debattista 2017): half the bars at 10^10.4.
    assert sd.boxy_peanut_fraction(10.0 ** 10.37) == pytest.approx(0.5, abs=0.01)
    assert sd.boxy_peanut_fraction(10.0 ** 9.48) < 0.05 < 0.95 < sd.boxy_peanut_fraction(10.0 ** 11.3)
    # The bar's angle is uniform: a quarter of the seeds in each quarter of 0 to 180 degrees, within sampling error.
    for low in (0.0, 45.0, 90.0, 135.0):
        count = sum(low <= s.bar_angle_deg < low + 45.0 for _, s in STRUCTURES)
        assert abs(count - len(STRUCTURES) / 4.0) <= 4.0 * math.sqrt(len(STRUCTURES) * 0.25 * 0.75)


def test_the_young_discs_height_for_a_milky_way_twin_is_inside_the_sourced_band():
    """GEN.1182 and GEN.1188: young OB stars, HII regions, Cepheids and open clusters of the Milky Way have a scale height of
    about 50 pc, 36 to 66 by tracer (Bobylev and Bajkova 2015 and 2016, Urquhart 2011, Bonatto 2006, Joshi 2016, Hao 2021,
    Moreira 2025), against a 300 pc old thin disc: young over old 0.12 to 0.23. A Milky Way check, not a law for every galaxy
    (six low-mass edge-on discs show 0.40 to 0.59, Seth 2005)."""
    reference = galaxyStructure.REFERENCE_STRUCTURE
    young = galaxyStructure.population_height_ratios(reference)["young"]
    assert 0.12 <= young <= 0.23
    assert 35.0 <= 300.0 * young <= 70.0
    low, high = 5.0, 9.0
    assert low <= sd.GAS_DISPERSION_KMS[0] <= high, "the CO line widths of GEN.536 span 5 to 9 km/s"


def test_the_populations_heights_and_dispersions_come_from_one_law():
    reference = galaxyStructure.REFERENCE_STRUCTURE
    ratios = galaxyStructure.population_height_ratios(reference)
    sigma = {name: galaxyStructure.disc_vertical_dispersion_kms(reference, name) for name in ratios}
    assert ratios["old"] == 1.0
    for name in ratios:
        assert ratios[name] == pytest.approx((sigma[name] / sigma["old"]) ** 2)
    assert 0.0 < ratios["young"] < ratios["intermediate"] < 1.0
    # The Milky Way check of the gap plan (4.3): the old disc's vertical dispersion at the Sun is within 35 percent of
    # the 21.1 km/s of Sharma and Hayden 2020.
    assert abs(sigma["old"] / 21.1 - 1.0) < 0.35
    totals = [galaxyStructure.disc_dispersion_kms(reference, name) for name in ratios]
    assert totals == sorted(totals)
    for _mass, structure in STRUCTURES[:200]:
        assert 0.0 < galaxyStructure.population_height_ratios(structure)["young"] < 1.0
        assert structure.disc_central_dispersion_kms > structure.disc_gas_dispersion_kms


def test_the_halo_and_bulge_dispersions_follow_the_circular_speed():
    reference = galaxyStructure.REFERENCE_STRUCTURE
    assert galaxyStructure.halo_dispersion_kms(reference) == pytest.approx(math.sqrt(3.0 / 3.5) * 220.0)
    assert galaxyStructure.bulge_dispersion_kms(reference, 0.0) == pytest.approx(math.sqrt(3.0) * 110.0)
    assert galaxyStructure.bulge_dispersion_kms(reference, 0.5) == pytest.approx(math.sqrt(3.0) * 220.0)
    scatter = [structure.bulge_dispersion_scatter_dex for _, structure in STRUCTURES]
    assert max(abs(value) for value in scatter) <= sd.SCATTER_LIMIT * sd.BULGE_DISPERSION_SCATTER_DEX + 1e-12
    mean = sum(scatter) / len(scatter)
    spread = math.sqrt(sum((v - mean) ** 2 for v in scatter) / (len(scatter) - 1))
    assert spread == pytest.approx(sd.BULGE_DISPERSION_SCATTER_DEX, rel=0.1)


def test_the_density_model_reads_the_drawn_thick_disk_and_bulge():
    for number in (2, 5, 9, 13):
        drawn = shape_draw.draw_galaxy(_seed(number))
        shape = density.build_galaxy_shape(
            drawn.disk_scale_length_pc, drawn.disk_scale_height_pc, drawn.bulge_scale_radius_pc, drawn.bulge_amplitude,
            drawn.arm_count, math.radians(drawn.pitch_angle_deg), 0.4, structure=drawn.structure)
        masses = density.integrated_masses(shape, steps=48)
        assert masses["thick_disk"] / masses["thin_disk"] == pytest.approx(drawn.structure.thick_mass_ratio, rel=0.03)
        assert density.bulge_to_total(shape) == pytest.approx(drawn.bulge_to_total, rel=1e-9)
        terms = density.model_terms(shape)
        assert terms["thick_disk_scale_length_pc"] == pytest.approx(
            drawn.structure.thick_length_ratio * shape.disk_scale_length_pc)
        assert terms["thick_disk_scale_height_pc"] == pytest.approx(
            drawn.structure.thick_height_ratio * shape.disk_scale_height_pc)
        assert terms["bulge_scale_y_pc"] == pytest.approx(drawn.structure.bar_axis_ratio_y * shape.bulge_scale_radius_pc)


def test_an_unbarred_galaxy_has_no_bar_gap():
    barred = next(s for _, s in STRUCTURES if s.bar_present)
    unbarred = next(s for _, s in STRUCTURES if not s.bar_present)
    kwargs = dict(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
                  bulge_amplitude=3.0, arm_count=2, pitch_angle_rad=0.26, arm_amplitude=0.4)
    with_bar = density.build_galaxy_shape(**kwargs, structure=barred)
    without = density.build_galaxy_shape(**kwargs, structure=unbarred)
    inside = 0.5 * star_formation.gap_radius_pc(with_bar)
    assert star_formation.recent_factor(inside, with_bar) < 0.5 * star_formation.recent_factor(inside, without)


def test_the_stored_shape_row_gains_nothing():
    """The structure is rebuilt from the seed: no column of `galaxy_shape` holds any field of it. The component columns
    DB.66 added (the bar, the thick disk, the spheroid of a family) are NULL or at their default for a spiral, which
    means "as the structure draw gives them"; only a family that sets them stores them (docs/design/galaxy/galaxy-types.md
    section 10)."""
    schema = (Path(__file__).resolve().parents[1] / "planetgen" / "db" / "schema.sql").read_text()
    table = re.search(r"CREATE TABLE IF NOT EXISTS galaxy_shape \((.*?)\n\)", schema, re.S)
    assert table is not None
    overrides = ("bar_present", "bar_axis_ratio_b", "bar_axis_ratio_c", "bar_angle_rad", "thick_disk_amplitude")
    for name in galaxyStructure.Structure._fields:
        if name not in overrides:
            assert not re.search(rf"\b{name}\b", table.group(1)), name
    for name in overrides:
        assert re.search(rf"\b{name}\s+(TINYINT|DOUBLE),", table.group(1)), f"{name} must stay nullable with no default (DB.66)"
