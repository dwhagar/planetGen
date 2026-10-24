# tests/test_galaxy_census.py

"""
GEN.607: the Census dial. `--richness` is the systems each solar mass of stars makes (the inverse of the mean system
mass): Auto is what the generator derives for the galaxy's populations (GEN.545), a pinned value sets the galaxy's
number of systems and the mean mass of its systems, and the plan stores it with the other pins.
"""

import math
import types

import pytest

from planetgen.galaxy import density, dials, shape_draw, system_mass
from tests.test_shape_pins import SEED, _load, _save

pytestmark = pytest.mark.fresh_database


def _args(**given):
    return types.SimpleNamespace(**{**{row.dest: None for row in dials.LEVERS}, **given})


def _shape(**given):
    return shape_draw.plan_shape(_args(**given), SEED)[0]


def test_the_richness_is_an_advanced_row_of_the_census_dial_with_a_range_and_a_source():
    row = next(row for row in dials.LEVERS if row.name == "richness")
    assert row.dial == "census" and row.advanced and row.flag == "--richness"
    assert (row.minimum, row.maximum) == (1.0, 3.0)
    assert row.source.strip() and row.status in dials.STATUSES
    assert row.problem(0.5) is not None and row.problem(3.5) is not None and row.problem(2.0) is None


def test_auto_is_the_mean_system_mass_the_generator_derives():
    shape = _shape()
    assert shape.system_mass_scale == 1.0
    derived = density.mean_system_mass_solar(shape)
    assert 1.0 < 1.0 / derived < 3.0     # inside the band the dial offers


def test_a_pinned_richness_sets_the_mean_system_mass_and_the_number_of_systems():
    auto, pinned = _shape(), _shape(richness=1.5)
    assert density.mean_system_mass_solar(pinned) == pytest.approx(1.0 / 1.5, rel=1e-9)
    stellar = pinned.structure.stellar_mass_solar
    # k_norm holds the galaxy's systems (stellar mass over the mean system mass) in the same unit as before
    ratio = (1.5 * density.mean_system_mass_solar(auto))
    assert pinned.k_norm / auto.k_norm == pytest.approx(ratio, rel=1e-9)
    assert stellar * 1.5 == pytest.approx(stellar / density.mean_system_mass_solar(pinned), rel=1e-9)


def test_every_population_mean_follows_the_pinned_richness_by_one_factor():
    auto, pinned = _shape(), _shape(richness=2.5)
    factor = pinned.system_mass_scale
    for name, mean in system_mass.population_means_solar(auto).items():
        assert system_mass.population_means_solar(pinned)[name] == pytest.approx(mean / factor, rel=1e-12)


def test_the_richness_alone_changes_nothing_but_the_systems_count():
    auto, pinned = _shape(), _shape(richness=1.2)
    for field in ("disk_scale_length_pc", "disk_scale_height_pc", "bulge_amplitude", "arm_count", "masses"):
        assert getattr(auto, field) == getattr(pinned, field)
    assert pinned.k_norm < auto.k_norm      # fewer systems for the same mass


def test_the_plan_stores_the_richness_and_the_loaded_shape_has_the_same_mean_system_mass(mysql_config):
    saved, _drawn = _save(mysql_config, richness=1.7)
    info, shape_json, _central = _load(mysql_config)
    assert info.pins == {"richness": 1.7}
    assert density.mean_system_mass_solar(info.shape) == pytest.approx(1.0 / 1.7, rel=1e-6)
    assert info.shape.system_mass_scale == pytest.approx(saved.system_mass_scale, rel=1e-6)
    assert "system_mass_scale" not in shape_json


def test_an_unpinned_galaxy_loads_with_the_derived_mean_system_mass(mysql_config):
    saved, _drawn = _save(mysql_config)
    info, _shape_json, _central = _load(mysql_config)
    assert info.pins == {} and info.shape.system_mass_scale == 1.0
    assert density.mean_system_mass_solar(info.shape) == pytest.approx(density.mean_system_mass_solar(saved), rel=1e-6)


def test_the_page_turns_the_richness_into_its_flag():
    from planetgen.web import generate_page
    argv = generate_page.plan_argv({"richness": "1.8"})
    assert argv[argv.index("--richness") + 1] == "1.8"
    assert math.isfinite(float(argv[argv.index("--richness") + 1]))
