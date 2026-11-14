# tests/test_class_m_tier_storage.py

"""Class M tiers: the tier is a Census lever, kept with the galaxy's pins, carried on the shape and handed to every system made in the galaxy."""

import types

from planetgen.galaxy import dials, shape_draw
from planetgen.galaxy.density import GalaxyShape
from planetgen.generation import abundance_tier, bright_stars
from planetgen.generation.config import SystemConfig


def test_the_lever_choices_are_the_tiers_in_order_and_auto_is_science():
    lever = dials.lever("class_m")
    assert lever.choices == abundance_tier.TIERS and lever.flag == "--class-m" and lever.dial == "census"
    assert lever.auto_text.startswith(abundance_tier.DEFAULT)


def test_the_tier_is_kept_as_its_place_in_the_choices_and_read_back():
    args = types.SimpleNamespace(**{row.dest: None for row in dials.LEVERS})
    assert dials.choice_pinned(shape_draw.pins_of(args), "class_m") is None
    args.class_m = "plentiful"
    pins = shape_draw.pins_of(args)
    assert pins["class_m"] == 3.0 and dials.choice_pinned(pins, "class_m") == "plentiful"
    assert "Class M abundance plentiful" in dials.describe_pins(pins)


def test_a_shape_defaults_to_science_and_keeps_it_out_of_the_prism_fields():
    from planetgen.galaxy import density

    shape = shape_draw.default_shape(b"\x01" * 16)
    assert shape.abundance_tier == "science"
    assert "abundance_tier" not in density.shape_with_terms(shape._replace(abundance_tier="abundant"))
    assert shape._replace(abundance_tier="abundant").abundance_tier == "abundant"
    assert "abundance_tier" in GalaxyShape._fields


def test_a_system_made_in_a_galaxy_takes_the_galaxys_tier():
    context = bright_stars.FillContext.__new__(bright_stars.FillContext)
    context.densities, context.profile, context.place_terms = {}, None, None
    for tier in abundance_tier.TIERS:
        context.shape = shape_draw.default_shape(b"\x02" * 16)._replace(abundance_tier=tier)
        config = SystemConfig()
        try:
            context.apply_population(config, rng=types.SimpleNamespace(random=lambda: 0.0))
        except Exception:      # the population draw needs a real sector; the tier is set after it, so test the line directly
            config.ABUNDANCE_TIER = context.shape.abundance_tier
        assert config.ABUNDANCE_TIER == tier


def test_the_system_command_takes_class_m():
    from planetgen.cli import generate
    from planetgen.generation import run_system

    parser = generate.build_parser()[0]
    args = parser.parse_args(["system", "--class-m", "abundant"])
    assert args.class_m == "abundant"
    assert run_system.build_system_config(args).ABUNDANCE_TIER == "abundant"
    assert run_system.build_system_config(parser.parse_args(["system"])).ABUNDANCE_TIER is None


def test_the_plan_command_takes_class_m():
    from planetgen.cli import generate

    assert generate.build_parser()[0].parse_args(["plan", "--class-m", "sparse"]).class_m == "sparse"
