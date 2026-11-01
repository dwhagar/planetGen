# tests/test_protection.py

"""
GEN.208: the habitability paragraph (`planetgen.physics.protection`) says where a body falls in each
of the six categories and what gear it needs, from one table of limits
(`planetgen.physics.protection_limits`) that the score's bands are built from.
"""

import random

import pytest

from planetgen import describe
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


def _paragraph(world, compact_host=None):
    """The rating written as one paragraph, as the pages' ratings table says it (the text the assertions below read)."""
    result = protection.rating(world, compact_host)
    gear = f" {result['gear']}" if result["gear"] else ""
    parts = [f"{_bold('PHI-6 habitability')} is {result['overall']:.2f}.{gear}"]
    for row in result["rows"]:
        parts.append(f"{_bold(row['domain'] + ': ' + row['tier'])} ({row['score']:.2f}). {row['text']}")
    return " ".join(parts)


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


def _chosen_pairs():
    """The limits where the built table and the research notes gave two different values, as `(name, built, notes, value in use)`."""
    out = []
    for row in protection_limits.LIMIT_ROWS:
        if row.built is None or row.notes is None or row.built == row.notes:
            continue
        out.append((row.name, row.built, row.notes, row.value))
    return out


def test_where_the_notes_and_the_built_table_differ_the_more_permissive_value_is_in_use():
    pairs = {name: (built, notes, value) for name, built, notes, value in _chosen_pairs()}
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
    (15.0, None), (-10.0, "warm clothing"), (-30.0, "an insulated suit"), (-80.0, "an environment suit with active heating"),
    (40.0, "cooling garments"), (80.0, "an actively cooled environment suit"),
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
    assert protection.air_need(_world(gases={"O2": 0.21, "CO2": 0.6}))[1] == "an environment suit with its own air"


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
    text = _paragraph(world)
    assert "\n" not in text
    assert text.startswith(f"**PHI-6 habitability** is {scores['PHI-4']:.2f}")
    assert "PHI-4" not in text
    for domain in hab.DOMAINS:
        assert f"**{domain}: {scores[domain][1]}** ({scores[domain][0]:.2f})." in text, domain
    assert text.index("**Pressure:") < text.index("**Temperature:") < text.index("**Air:") < text.index(
        "**Water:") < text.index("**Radiation:") < text.index("**Suitability:")


def test_the_lines_of_an_airless_frozen_world_read_like_the_examples_boss_gave():
    """GEN.216: Boss's own example lines (2026-10-15) are the wording."""
    world = _world(pressure_kpa=0.0, temperature_c=-227.0, gases={}, solvent=None, life_stage="none",
                   water_source=None, surface_dose_sv_yr=0.742)
    text = _paragraph(world)
    for line in (
            "**PHI-6 habitability** is 0.00. Required: hardened pressure suit + dust-tight seals. Recommended: hardened habitat.",
            "**Pressure: Red** (0.00). The air pressure is 0 kPa. Required: hardened pressure suit. "
            "**Temperature: Red** (0.00). The air temperature is -227 C. Required: hardened pressure suit. "
            "**Air: Red** (0.00). There is no atmosphere (0 kPa). Required: hardened pressure suit + dust-tight seals.",
            "**Water: Red** (0.00). There is no liquid water to support life.",
            "**Radiation: Yellow** (0.45). The radiation dose is a mean of 742 mSv a year, inside the Yellow band "
            "(up to 1 Sv a year). Required: hardened pressure suit. Recommended: hardened habitat."):
        assert line in text, line
    assert text.index("Suitability: Red") > text.index("**Radiation:")
    assert "There is no food or water readily available, and only limited shelter." in text


def test_a_habitable_world_leaves_both_equipment_lines_off():
    text = _paragraph(_world())
    assert text.startswith("**PHI-6 habitability** is 1.00. **Pressure: Blue** (1.00).")
    assert "Required" not in text and "Recommended" not in text
    assert "**Suitability: Blue** (1.00). Water, food and shelter are all readily available." in text
    assert "Also considered" not in text


def test_the_suitability_line_names_what_is_and_is_not_readily_available():
    text = _paragraph(_world(life_stage="none", water_source="liquid"))
    assert "There is no food readily available, but there is water and shelter available." in text
    dry = _paragraph(_world(life_stage="none", solvent=None, water_source=None))
    assert "There is no food or water readily available" in dry
    assert "Suitability: Red**" in dry


def test_suitability_does_not_call_water_available_when_the_water_row_asks_for_treatment():
    """GEN.691: the water rung is reachability, the Water row is chemistry; reachable water that is undrinkable is
    said to be reachable and the treatment is repeated, never "there is water available"."""
    heavy = _world(ph=0.5, life_stage="none", water_source="liquid")
    rows = {row["domain"]: row["text"] for row in protection.rating(heavy)["rows"]}
    assert "undrinkable without heavy treatment" in rows["Water"]
    assert rows["Suitability"].startswith("There is no food readily available, but there is shelter available, and the water "
                                          "is reachable but it is undrinkable without heavy treatment."), rows["Suitability"]
    fed = _world(ph=0.5, life_stage="biosphere", water_source="liquid")
    rows = {row["domain"]: row["text"] for row in protection.rating(fed)["rows"]}
    assert rows["Suitability"].startswith("Food and shelter are readily available, and the water is reachable but it is "
                                          "undrinkable without heavy treatment."), rows["Suitability"]
    mild = _world(ph=2.0, life_stage="none", water_source="liquid")
    rows = {row["domain"]: row["text"] for row in protection.rating(mild)["rows"]}
    assert "the water is reachable but it needs treatment before drinking" in rows["Suitability"]


def test_suitability_never_says_water_is_available_against_a_water_row_that_wants_treatment():
    for ph in (0.5, 2.0, 4.0, 6.0, 8.1, 10.0, 12.0, 13.5):
        for life in ("none", "microbial", "biosphere"):
            for source in ("liquid", "ice", "hydrated", None):
                world = _world(ph=ph, life_stage=life, water_source=source)
                text = {row["domain"]: row["text"] for row in protection.rating(world)["rows"]}["Suitability"]
                treatment = protection.water_need(world)[1]
                if treatment and world.solvent == "water" and "reachable" not in text:
                    # no plain "water ... available" without the treatment note: only a rung that is not good says so
                    assert "water and" not in text and "there is water" not in text and "Water, food" not in text, (ph, life, source, text)
                if "the water is reachable" in text:
                    assert treatment in text, (ph, life, source, text)


def test_the_gear_of_a_cold_world_is_named_in_short_words_and_the_first_line_names_the_heated_suit():
    text = _paragraph(_world(temperature_c=-80.0))
    assert "**Temperature: Red** (0.23). The air temperature is -80 C. Required: heated sealed suit + self-contained air. **Air:" in text
    assert "Required: heated sealed suit + self-contained air." in text


def test_the_gear_of_a_hot_high_pressure_world_names_each_suit():
    text = _paragraph(_world(pressure_kpa=2_000.0, temperature_c=80.0))
    assert "The air pressure is 2,000 kPa. Required: heated diving suit + mixed-gas rig. **Temperature:" in text
    assert "Required: heated diving suit + mixed-gas rig." in text


def test_a_planet_of_a_compact_host_says_its_dose_is_rated_lethal():
    world = _world(surface_dose_sv_yr=1000.0)
    text = _paragraph(world, compact_host="pulsar")
    assert "A planet of a pulsar is rated as a lethal dose." in text and "Required: hardened suit + hardened habitat. A planet of a pulsar" in text


def test_a_sample_of_generated_planets_and_moons_each_carry_the_ratings_block():
    seen = 0
    for seed in range(8):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        document = describe.describe(system)
        by_name = {section["title"]: section for section in describe.model.sections(document)}
        planets = [planet for planet in system.planets if hasattr(planet, "moons")]
        bodies = planets + [moon for planet in planets for moon in planet.moons]
        for body in bodies:
            blocks = [b for b in by_name[body.name]["blocks"] if b["type"] == "ratings"]
            if body.body_type == "g":
                assert blocks == []
                continue
            values = hw.values_from_body(body)
            world = hw.world_from(values, lethal_host=None)
            scores = hab.phi4(world)
            assert len(blocks) == 1
            for row in blocks[0]["rows"]:
                assert row["tier"] == scores[row["domain"]][1]
            seen += 1
    assert seen > 10


def test_the_paragraph_lengths_stay_short_enough_to_read():
    rng = random.Random(208)
    for _ in range(200):
        world = _world(pressure_kpa=10 ** rng.uniform(-2, 5), temperature_c=rng.uniform(-200, 400),
                       surface_dose_sv_yr=10 ** rng.uniform(-3, 3), life_stage=rng.choice(hab.LIFE_STAGES))
        assert len(_paragraph(world)) < 1100


def test_the_air_line_of_a_world_with_no_atmosphere_says_so_and_is_red():
    """GEN.215: Air is Red and the text says there is no atmosphere, never that the gases are fine."""
    world = _world(pressure_kpa=0.0, gases={})
    text = _paragraph(world)
    assert "**Air: Red** (0.00). There is no atmosphere (0 kPa). Required: " in text and "dust-tight seals" in text
    assert "within their long-term limits" not in text


# --- ported from the retired habitability_explain (GEN.218) ---------------------------------------------------------

def _earth():
    return World(pressure_kpa=101.0, temperature_c=15.0, gases={"O2": 0.21, "CO2": 0.0004}, gravity_ms2=9.8,
                 life_stage="biosphere")


def test_an_earth_like_world_is_one_paragraph_of_six_blue_categories_with_their_values():
    text = _paragraph(_earth())
    assert text.startswith("**PHI-6 habitability** is 1.00. **Pressure: Blue**")
    for domain in ("Pressure", "Temperature", "Air", "Water", "Radiation", "Suitability"):
        assert f"**{domain}: Blue** (1.00)." in text
    assert "The air pressure is 101 kPa, so no gear is needed." in text
    assert "The air temperature is 15 C, so no gear is needed." in text


def test_a_gas_past_its_long_term_limit_asks_for_a_scrubber_and_is_named_once_in_the_one_gas_table():
    world = World(pressure_kpa=101.0, temperature_c=15.0, gases={"O2": 0.21, "CO2": 0.03}, gravity_ms2=9.8)
    text = _paragraph(world)
    assert "Required: CO2 scrubber. **Water:" in text and "**Air: Yellow**" in text and "carbon dioxide" in text.lower()
    assert not hasattr(protection, "GAS_NAMES")
    assert protection._gas("H2S") == "hydrogen sulfide"


def _worlds(count=400):
    rng = random.Random(91)
    for _ in range(count):
        yield World(
            pressure_kpa=10 ** rng.uniform(-1, 5), temperature_c=rng.uniform(-150, 500),
            gases={"O2": rng.choice([0.0, 0.05, 0.21, 0.4, 0.9]), "CO2": rng.choice([0.0, 0.0004, 0.01, 0.05, 0.5]),
                   "CO": rng.choice([0.0, 1e-6, 1e-3]), "H2S": rng.choice([0.0, 1e-6, 1e-3])},
            gravity_ms2=rng.uniform(1, 25), relative_humidity=rng.uniform(0.05, 0.99), solvent="water",
            ph=rng.choice([None, 2.0, 7.0, 10.0, 13.0]), water_activity=rng.choice([1.0, 0.8, 0.5]),
            chaotropicity_kj_kg=rng.choice([0.0, 50.0, 90.0]), surface_dose_sv_yr=10 ** rng.uniform(-3, 3),
            life_stage=rng.choice(hab.LIFE_STAGES), water_source=rng.choice(["liquid", "ice", "hydrated", "vapour", None]),
            life_on_surface=rng.random() < 0.8)


def test_the_rows_colours_and_scores_always_equal_the_stored_ones():
    for world in _worlds():
        stored = hab.phi4(world)
        result = protection.rating(world)
        assert result["overall"] == stored["PHI-4"]
        for row in result["rows"]:
            assert (row["score"], row["tier"]) == stored[row["domain"]]


def _system_document():
    for seed in range(1, 60):
        random.seed(seed)
        document = describe.describe(StarSystem(SystemConfig()))
        if "PHI-6 habitability" in describe.render(document, "text"):
            return document
    raise AssertionError("no rocky body in 60 systems")


def test_a_generated_planet_or_moon_description_carries_the_ratings_in_each_format():
    document = _system_document()
    for fmt, header in (("text", "Category     Rating"), ("markdown", "| Category | Rating | Score | Notes |")):
        rendered = describe.render(document, fmt)
        assert "PHI-6 habitability is" in rendered and header in rendered
        for domain in ("Pressure", "Temperature", "Air", "Water", "Radiation", "Suitability"):
            assert domain in rendered


def test_there_is_one_describer_of_the_six_categories():
    """`protection.rating` words the six categories. `physics/habitability_explain.py` came back with GEN.445 to name the
    composition behind the Radiation term only (`radiation_terms`), so it may exist, but it may not grow a second describer
    of the categories: its public functions are the radiation ones."""
    from planetgen.physics import habitability_explain
    public = {name for name, value in vars(habitability_explain).items()
              if not name.startswith("_") and getattr(value, "__module__", None) == habitability_explain.__name__}
    assert public == {"RadiationTerms", "air_shares", "radiation_terms"}
    assert not hasattr(protection, "radiation_terms")


def test_a_small_dose_prints_in_the_ladders_scientific_form_not_as_a_python_exponent():
    text = protection.radiation_need(_world(surface_dose_sv_yr=2.0e-7))[2]
    assert "e-" not in text and "10" in text


def _rows(world):
    return {row["domain"]: row for row in protection.rating(world)["rows"]}


def test_an_item_sits_under_the_rows_that_ask_for_it_not_under_a_blue_one():
    """GEN.702: a cold world's own air is for the cold (Air is Blue), a hardened pressure suit is for the pressure and the
    radiation, and a suit that seals keeps out the cold and the thin air as well."""
    cold = _rows(_world(temperature_c=-80.0))
    assert "self-contained air" in cold["Temperature"]["text"] and "Required" not in cold["Air"]["text"]
    mars = _rows(World(pressure_kpa=0.6, temperature_c=-60.0, gases={"CO2": 0.95}, surface_dose_sv_yr=0.236,
                       dose_gcr_sv_yr=0.236, dose_sep_sv_yr=0.0, dose_ground_sv_yr=0.001, gravity_ms2=3.7))
    for domain in ("Pressure", "Temperature", "Radiation"):
        assert "hardened pressure suit" in mars[domain]["text"], domain
    assert "hardened habitat" in mars["Radiation"]["text"] and "hardened habitat" not in mars["Pressure"]["text"]


def test_radiation_need_names_the_ultraviolet_when_it_is_the_worse_term():
    """GEN.445: the description agrees with the score; a Blue dose under unfiltered ultraviolet is not 'inside the Blue band'."""
    rank, gear, why = protection.radiation_need(_world(surface_dose_sv_yr=0.002, uv_index=1000.0))
    assert gear == "a hardened, shielded habitat" and "ultraviolet" in why and "Blue" not in why
    rank, gear, why = protection.radiation_need(_world(surface_dose_sv_yr=0.002, uv_index=100.0))
    assert gear == "a shielded shelter for long stays" and "ultraviolet" in why
    assert "ultraviolet" not in protection.radiation_need(_world(surface_dose_sv_yr=0.5, uv_index=1.0))[2]
