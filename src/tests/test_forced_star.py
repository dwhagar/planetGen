"""
Forced star types are stars of the model (GEN.941, `generation/forced_star.py`).
"""
import pytest

from planetgen import tuning
from planetgen.generation import forced_star
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import constants, post_main_sequence


def make_star(star_type, age=None):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    if age is not None:
        cfg.AGE = age
    return Star(cfg)


@pytest.mark.parametrize("star_type", ["G2V", "M5V", "K2III", "F5IV", "A0V", "B2V", "M3IAB", "G5VII"])
def test_a_type_the_model_makes_is_one_of_its_stars(star_type):
    for _ in range(10):
        star = make_star(star_type)
        assert star.type.startswith(star_type)
        assert star.initial_mass_sol is not None
        assert star.age <= star.lifespan * 1.001
        assert star.age <= tuning.UNIVERSE_AGE_GY * 1.001


def test_a_giant_is_a_star_that_has_left_its_main_sequence():
    for _ in range(20):
        star = make_star("K2III")
        line = post_main_sequence.timeline(star.initial_mass_sol, star.metallicity_fe_h)
        assert star.age >= line.main_sequence_gy * 0.999


def test_a_sun_like_type_is_a_sun_like_star_of_the_model():
    sun_like = [make_star("G2V") for _ in range(20)]
    for star in sun_like:
        assert 0.6 <= star.mass / constants.SOLAR_MASS_TO_KG <= 2.0
        assert 5400 <= star.temperature <= 6000


def test_a_type_the_model_never_makes_is_its_nearest_star_at_the_temperature_asked():
    letter, subclass = "B", 6
    assert (letter, subclass, "III") not in forced_star._table(0.0)
    star = make_star("B6III")
    assert star.type.startswith("B6III")
    assert star.yerkes_class == "III"
    assert star.temperature == forced_star.target_temperature_k(letter, subclass)


def test_the_aliases_are_the_models_classes_under_their_own_names():
    hypergiant = make_star("O5IA+")
    assert hypergiant.yerkes_class == "IA+"
    assert hypergiant.type.startswith("O5IA+")
    white_dwarf = make_star("A0D")
    assert white_dwarf.yerkes_class == "D"
    assert white_dwarf.lifespan == float("inf")


def test_a_subdwarf_is_an_old_metal_poor_star():
    for _ in range(10):
        star = make_star("K5VI")
        assert star.age >= tuning.SUBDWARF_MIN_AGE_GY
        assert tuning.SUBDWARF_FE_H_RANGE[0] <= star.metallicity_fe_h <= tuning.SUBDWARF_FE_H_RANGE[1]


def test_the_age_bias_picks_the_older_or_younger_half():
    young = [make_star("G2V", "young").age for _ in range(30)]
    old = [make_star("G2V", "old").age for _ in range(30)]
    assert sum(old) / len(old) > sum(young) / len(young)


def test_the_class_pages_read_the_model():
    sun = forced_star.letter_summary("G")
    assert 0.8 <= sun.mass_sol[0] and sun.mass_sol[1] <= 1.2
    assert 0 < sun.share < 1
    giants = forced_star.class_summary("III")
    assert giants.luminosity_sol[0] < giants.luminosity_sol[1]
    with pytest.raises(ValueError):
        forced_star.class_summary("XX")
