# tests/test_protection.py

"""
GEN.208: the habitability paragraph (`planetgen.physics.protection`) says where a body falls in each
of the six categories and what gear it needs, from one table of limits
(`planetgen.physics.protection_limits`) that the score's bands are built from.
"""

import random

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import habitability as hab
from planetgen.physics import habitability_world as hw
from planetgen.physics import protection, protection_limits
from planetgen.physics.habitability import World
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy
from tests.test_habitability import EXAMPLES


def _bold(text):
    return f"**{text}**"


def _world(**changes):
    values = dict(pressure_kpa=101.3, temperature_c=15.0, gases={"N2": 0.78, "O2": 0.21, "CO2": 0.0004}, ph=8.1,
                  relative_humidity=0.5, life_stage="biosphere")
    values.update(changes)
    return World(**values)


# --- the one table ----------------------------------------------------------------------

def test_the_score_bands_are_built_from_the_limit_table():
    limit = protection_limits.limit
    assert hab.PRESSURE_BAND == (limit("pressure_low_outer_kpa"), limit("armstrong_kpa"), limit("mask_min_kpa"),
                                 limit("pressure_no_gear_kpa"), limit("pressure_blue_max_kpa"),
                                 limit("pressure_narcosis_kpa"), limit("pressure_hyperbaric_max_kpa"),
                                 limit("pressure_outer_kpa"))
    assert hab.TEMPERATURE_BAND[1:7] == (limit("cold_suit_c"), limit("cold_clothing_c"), limit("comfort_low_c"),
                                         limit("comfort_high_c"), limit("heat_garment_c"), limit("heat_max_c"))
    assert hab.DOSE_BAND[4:7] == (limit("dose_blue_sv_yr"), limit("dose_green_sv_yr"), limit("dose_yellow_sv_yr"))
    assert hab.WET_BULB_CRIT_C == limit("wet_bulb_c")
    # The mask limit is derived from the Armstrong limit and the inspired oxygen the table also holds.
    assert hab.MASK_MIN_KPA == pytest.approx(limit("mask_min_kpa"))


def test_where_the_notes_and_the_built_table_differ_the_more_permissive_value_is_in_use():
    pairs = {name: (built, notes, value) for name, built, notes, value in protection_limits.chosen_pairs()}
    assert pairs == {
        "mask_min_kpa": (14.3, 20.0, 14.3), "pressure_no_gear_kpa": (50.0, 70.0, 50.0),
        "pressure_hyperbaric_max_kpa": (10000.0, 5000.0, 10000.0), "inspired_o2_min_kpa": (8.0, 19.5, 8.0),
        "cold_suit_c": (-50.0, -40.0, -50.0), "wet_bulb_c": (31.0, 30.0, 31.0)}
    for name, (built, notes, value) in pairs.items():
        row = protection_limits.LIMITS[name]
        pick = min if row.looser == "lower" else max
        assert value == pick(built, notes)


# --- the gear ladder --------------------------------------------------------------------

@pytest.mark.parametrize("pressure,gear", [
    (101.0, None), (30.0, "a breathing mask"), (10.0, "a counterpressure or pressurised suit"),
    (3.0, "a full gas-pressurised suit"), (2_000.0, "a hyperbaric pressure suit"),
    (50_000.0, "a rigid armoured pressure suit")])
def test_the_pressure_gear_follows_the_ladder(pressure, gear):
    assert protection.pressure_need(_world(pressure_kpa=pressure))[1] == gear


@pytest.mark.parametrize("temperature,gear", [
    (15.0, None), (-10.0, "warm clothing"), (-30.0, "an insulated suit"), (-80.0, "a sealed suit with active heating"),
    (40.0, "cooling garments"), (80.0, "an actively cooled sealed suit"),
    (200.0, "a heat-shielded suit with a liquid cooling garment")])
def test_the_temperature_gear_follows_the_ladder(temperature, gear):
    assert protection.temperature_need(_world(temperature_c=temperature, relative_humidity=0.3))[1] == gear


def test_a_hot_wet_day_needs_cooling_garments_and_says_why():
    rank, gear, why = protection.temperature_need(_world(temperature_c=38.0, relative_humidity=0.9))
    assert gear == "cooling garments" and "wet bulb" in why
    assert protection.temperature_need(_world(temperature_c=25.0, relative_humidity=0.9))[1] is None


def test_the_air_gear_is_a_mask_a_scrubber_or_a_sealed_suit_with_its_own_air():
    assert protection.air_need(_world())[1] is None
    assert protection.air_need(_world(pressure_kpa=30.0, gases={"O2": 0.2}))[1] == "a breathing mask with supplied oxygen"
    assert protection.air_need(_world(gases={"O2": 0.21, "CO2": 0.03}))[1] == "a mask with a scrubber"
    assert protection.air_need(_world(gases={"O2": 0.21, "CO2": 0.6}))[1] == "a sealed suit with its own air"


def test_radiation_needs_a_shelter_not_a_suit():
    assert protection.radiation_need(_world(surface_dose_sv_yr=0.002))[1] is None
    assert protection.radiation_need(_world(surface_dose_sv_yr=0.05))[1] is None
    assert protection.radiation_need(_world(surface_dose_sv_yr=0.5))[1] == "a shielded shelter for long stays"
    rank, gear, why = protection.radiation_need(_world(surface_dose_sv_yr=50.0))
    assert gear == "a hardened, shielded habitat" and "past the Yellow band" in why and "no suit shields" in why


def test_no_liquid_water_means_it_has_to_be_brought():
    assert "brought or made from ice" in protection.water_need(_world(solvent=None))[1]
    assert protection.water_need(_world())[1] is None
    assert "treatment" in protection.water_need(_world(ph=4.5))[1]


# --- the paragraph ----------------------------------------------------------------------

@pytest.mark.parametrize("name", list(EXAMPLES))
def test_every_example_world_gets_one_paragraph_naming_each_category_its_colour_and_value(name):
    world = EXAMPLES[name][0]
    scores = hab.scores(world)
    text = protection.paragraph(world, _bold)
    assert "\n" not in text
    assert text.startswith(f"**PHI-6 habitability** is {scores['PHI-4']:.2f}, and a visitor needs ")
    assert "PHI-4" not in text
    for domain in hab.DOMAINS:
        assert f"**{domain}: {scores[domain][1]}** ({scores[domain][0]:.2f})." in text, domain
    assert text.index("**Pressure:") < text.index("**Temperature:") < text.index("**Air:") < text.index(
        "**Water:") < text.index("**Radiation:") < text.index("**Suitability:")


def test_the_lines_of_an_airless_frozen_world_read_like_the_examples_boss_gave():
    """GEN.216: Boss's own example lines (2026-10-15) are the wording."""
    world = _world(pressure_kpa=0.0, temperature_c=-227.0, gases={}, solvent=None, life_stage="none",
                   water_source=None, surface_dose_sv_yr=0.742)
    text = protection.paragraph(world, _bold)
    for line in (
            "**PHI-6 habitability** is 0.00, and a visitor needs a sealed suit.",
            "**Pressure: Red** (0.00). The air pressure is 0 kPa, so a pressure suit is needed.",
            "**Temperature: Red** (0.00). The air temperature is -227 C, so an environment suit is needed.",
            "**Air: Red** (0.00). There is no atmosphere (0.0 kPa).",
            "**Water: Red** (0.00). There is no liquid water to support life.",
            "**Radiation: Yellow** (0.45). The radiation dose is a mean of 742 mSv a year, inside the Yellow band "
            "(up to 1 Sv a year)."):
        assert line in text, line
    assert text.index("Suitability: Red") > text.index("**Radiation:")
    assert "There is no food or water readily available, and only limited shelter." in text


def test_a_habitable_world_says_it_needs_no_equipment_in_the_first_line():
    text = protection.paragraph(_world(), _bold)
    assert text.startswith("**PHI-6 habitability** is 1.00, and a visitor needs no equipment. **Pressure: Blue** (1.00).")
    assert "**Suitability: Blue** (1.00). Water, food and shelter are all readily available." in text
    assert "Also considered" not in text


def test_the_suitability_line_names_what_is_and_is_not_readily_available():
    text = protection.paragraph(_world(life_stage="none", water_source="liquid"), _bold)
    assert "There is no food readily available, but there is water and shelter available." in text
    dry = protection.paragraph(_world(life_stage="none", solvent=None, water_source=None), _bold)
    assert "There is no food or water readily available" in dry
    assert "Suitability: Red**" in dry


def test_the_gear_of_a_cold_world_is_named_in_short_words_and_the_first_line_says_a_sealed_suit():
    text = protection.paragraph(_world(temperature_c=-80.0), _bold)
    assert "The air temperature is -80 C, so an environment suit is needed." in text
    assert "and a visitor needs a sealed suit." in text


def test_the_gear_of_a_hot_high_pressure_world_names_each_suit():
    text = protection.paragraph(_world(pressure_kpa=2_000.0, temperature_c=80.0), _bold)
    assert "The air pressure is 2,000 kPa, so a hyperbaric suit is needed." in text
    assert "The air temperature is 80 C, so an environment suit is needed." in text
    assert "and a visitor needs a sealed suit." in text


def test_a_planet_of_a_compact_host_says_its_dose_is_rated_lethal():
    world = _world(surface_dose_sv_yr=1000.0)
    text = protection.paragraph(world, _bold, compact_host="pulsar")
    assert "A planet of a pulsar is rated as a lethal dose." in text and "so a hardened habitat is needed" in text


def test_a_sample_of_generated_planets_and_moons_each_have_one_consistent_paragraph():
    seen = 0
    for seed in range(8):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            config = SystemConfig()
            config.MARKDOWN = True
            system = StarSystem(config, galactic_center_dist_ly=26000.0)
        planets = [planet for planet in system.planets if hasattr(planet, "moons")]
        bodies = planets + [moon for planet in planets for moon in planet.moons]
        for body in bodies:
            if body.body_type == "g":
                assert body._habitability_paragraphs() == []
                continue
            paragraphs = body._habitability_paragraphs()
            assert len(paragraphs) == 1
            values = hw.values_from_body(body)
            world = hw.world_from(values, lethal_host=None)
            scores = hab.phi4(world)
            for domain in hab.DOMAINS:
                assert f"**{domain}: {scores[domain][1]}** (" in paragraphs[0]
            seen += 1
    assert seen > 10


def test_the_paragraph_lengths_stay_short_enough_to_read():
    rng = random.Random(208)
    for _ in range(200):
        world = _world(pressure_kpa=10 ** rng.uniform(-2, 5), temperature_c=rng.uniform(-200, 400),
                       surface_dose_sv_yr=10 ** rng.uniform(-3, 3), life_stage=rng.choice(hab.LIFE_STAGES))
        assert len(protection.paragraph(world, _bold)) < 1100


def test_the_air_line_of_a_world_with_no_atmosphere_says_so_and_is_red():
    """GEN.215: Air is Red and the text says there is no atmosphere, never that the gases are fine."""
    world = _world(pressure_kpa=0.0, gases={})
    text = protection.paragraph(world, _bold)
    assert "**Air: Red** (0.00). There is no atmosphere (0.0 kPa). **Water:" in text
    assert "within their long-term limits" not in text
