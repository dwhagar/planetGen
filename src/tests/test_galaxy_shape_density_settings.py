"""
The galaxy shape density settings (ADM.49): the arm density, the inter-arm
density, the core density and the bulge density each change the density where
they act, and only there; the settings reach the stored shape from the
command line and the Generate page.
"""

import argparse
import math

import pytest

from planetgen import tuning
from planetgen.cli.generate import add_plan_arguments
from planetgen.db import store
from planetgen.galaxy import density
from planetgen.web import generate_page
from tests.toy_shape import toy_shape

BASE = dict(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
            bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15.0))


def _shape(arm_density=1.4, interarm_density=0.6, bulge=3.11, nuclear=0.0):
    amplitude, level = density.arm_terms(arm_density, interarm_density)
    return toy_shape(**{**BASE, "bulge_amplitude": bulge}, arm_amplitude=amplitude,
                     arm_level=level, nuclear_amplitude=nuclear)


def _on_arm(shape, radius_pc=12000.0):
    """The arm crest's azimuth at `radius_pc`, found as the inter-arm minimum's opposite."""
    theta = density._interarm_angle(radius_pc, shape) + math.pi / shape.arm_count
    return (radius_pc * math.cos(theta), radius_pc * math.sin(theta), 0.0)


def _between_arms(shape, radius_pc=12000.0):
    theta = density._interarm_angle(radius_pc, shape)
    return (radius_pc * math.cos(theta), radius_pc * math.sin(theta), 0.0)


def test_the_usual_densities_are_the_old_arm_contrast():
    amplitude, level = density.arm_terms(1.4, 0.6)
    assert amplitude == pytest.approx(0.4) and level == pytest.approx(1.0)
    assert density.arm_terms(2.0, 0.0) == (1.0, 1.0)


@pytest.mark.parametrize("arm, interarm", [(0.0, 0.0), (1.0, 2.0), (1.0, -0.1), (-1.0, 0.0)])
def test_impossible_arm_densities_are_refused(arm, interarm):
    with pytest.raises(ValueError):
        density.arm_terms(arm, interarm)


def test_a_denser_arm_changes_the_density_on_the_arm_and_not_between_arms_much():
    base, dense = _shape(), _shape(arm_density=2.4)
    arm_before = density.relative_density(_on_arm(base), base)
    arm_after = density.relative_density(_on_arm(dense), dense)
    assert arm_after > arm_before * 1.2
    # Normalised at the inter-arm point, so the trough stays where it was.
    assert density.relative_density(_between_arms(dense), dense) == pytest.approx(
        density.relative_density(_between_arms(base), base), rel=0.5)


def test_a_denser_inter_arm_space_raises_the_trough_relative_to_the_crest():
    base, filled = _shape(), _shape(interarm_density=1.2)
    contrast_before = density.relative_density(_on_arm(base), base) / density.relative_density(_between_arms(base), base)
    contrast_after = density.relative_density(_on_arm(filled), filled) / density.relative_density(_between_arms(filled), filled)
    assert contrast_after < contrast_before
    assert density.relative_density(_between_arms(filled), filled) > 0.0


def test_the_crest_and_the_trough_of_the_thin_disk_are_the_two_densities():
    shape = _shape(arm_density=2.0, interarm_density=0.5)
    crest = density._arm_cosine(*_on_arm(shape)[:2], shape)
    trough = density._arm_cosine(*_between_arms(shape)[:2], shape)
    assert shape.arm_level * (1 + shape.arm_amplitude * crest) == pytest.approx(2.0)
    assert shape.arm_level * (1 + shape.arm_amplitude * trough) == pytest.approx(0.5)


def test_the_nuclear_structures_raise_the_density_at_the_centre_and_not_far_away():
    plain, nuclear = _shape(), _shape(nuclear=1.0)
    assert density.relative_density((0.0, 0.0, 0.0), nuclear) > 10 * density.relative_density((0.0, 0.0, 0.0), plain)
    far = (6000.0, 0.0, 0.0)
    assert density.relative_density(far, nuclear) == pytest.approx(density.relative_density(far, plain), rel=0.01)
    # Past the nuclear disc's few scale lengths only the bulge is left.
    out = (8 * density.model_terms(nuclear)["nuclear_disc_length_pc"], 0.0, 0.0)
    assert density.relative_density(out, nuclear) < 1.05 * density.relative_density(out, plain)


def test_a_denser_bulge_changes_the_density_in_the_bulge():
    plain, dense = _shape(), _shape(bulge=8.0)
    inside = (800.0, 300.0, 100.0)
    assert density.relative_density(inside, dense) > 1.5 * density.relative_density(inside, plain)


def test_the_bounds_still_bound_with_the_nuclear_structures_and_a_dense_arm():
    from planetgen.galaxy import skeleton
    shape = _shape(arm_density=3.0, interarm_density=0.2, nuclear=2.0)
    for radius in (0.0, 300.0, 3000.0, 9000.0):
        for z in (0.0, 200.0, 1000.0):
            bound = skeleton.bound_relative_density_at(shape, radius, z)
            for step in range(24):
                theta = step * math.pi / 12
                point = (radius * math.cos(theta), radius * math.sin(theta), z)
                assert density.relative_density(point, shape) <= bound * (1 + 1e-9) + density.halo_density(point, shape)


def test_the_populations_still_sum_to_the_density_with_the_nuclear_structures():
    shape = _shape(arm_density=2.0, interarm_density=0.3, nuclear=1.0)
    for point in ((0.0, 0.0, 0.0), (50.0, 20.0, 5.0), _on_arm(shape), _between_arms(shape)):
        parts = density.population_densities(point, shape)
        assert sum(parts.values()) == pytest.approx(density.relative_density(point, shape), rel=1e-9)


def test_the_command_line_and_the_page_offer_the_same_defaults():
    parser = argparse.ArgumentParser()
    add_plan_arguments(parser)
    args = parser.parse_args([])
    assert (args.arm_density, args.interarm_density, args.nuclear_amplitude) == (None, None, None)
    fixed = {row.name: row.auto_value for row in generate_page.PLAN_FIELDS}
    assert fixed["arm_density"] == 1.4 and fixed["interarm_density"] == 0.6 and fixed["nuclear_amplitude"] == 1.0


def test_the_page_turns_the_settings_into_plan_options():
    argv = generate_page.plan_argv({"arm_density": "2.2", "interarm_density": "0.4", "nuclear_amplitude": "2",
                                    "bulge_share": "0.5"})
    for flag, value in (("--arm-density", "2.2"), ("--interarm-density", "0.4"), ("--nuclear-amplitude", "2.0"),
                        ("--bulge-share", "0.5")):
        assert argv[argv.index(flag) + 1] == value


def test_a_plan_stores_the_settings_and_the_stored_shape_gives_the_same_density(mysql_config):
    from tests.test_galaxy_gen import _build_real_skeleton
    _build_real_skeleton(mysql_config, ["--arm-density", "2.0", "--interarm-density", "0.5", "--nuclear-amplitude", "0.5",
                                        "--bulge-share", "0.4"])
    conn = store.get_connection(mysql_config)
    try:
        stored = store.get_galaxy_shape(conn).shape
    finally:
        conn.close()
    assert stored.arm_level == pytest.approx(1.25) and stored.arm_amplitude == pytest.approx(0.6)
    assert stored.nuclear_amplitude == 0.5
    masses = density.component_masses(stored)
    assert masses["bulge"] / (masses["bulge"] + masses["thin_disk"] + masses["thick_disk"]) == pytest.approx(0.4, rel=0.02)
    terms = density.shape_with_terms(stored)
    full = density.model_terms(stored._replace(nuclear_amplitude=1.0))
    assert terms["nuclear_cluster_mass"] == pytest.approx(0.5 * full["nuclear_cluster_mass"])
    assert terms["nuclear_disc_mass"] == pytest.approx(tuning.NUCLEAR_DISC_CLUSTER_MASS_RATIO * terms["nuclear_cluster_mass"])
