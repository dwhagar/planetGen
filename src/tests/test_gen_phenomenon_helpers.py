"""
Phenomenon class helpers (TODO TEST.36).

Direct tests of the module-level class and designation helpers behind
the phenomena -- nebula (`nebulaData`), supernova remnant
(`supernovaRemnantData`), black hole (`compactRemnant`), rogue planet and
interstellar comet (`roguePlanetData`), star-bound comet (`cometData`)
and asteroid field (`asteroidFieldData`) -- plus every
`get_table_properties` (Star, Planet, BinaryStarProxy, WideBinaryPair,
BlackHole, NeutronStar): keys, string values, and the conditional wobble
rows.
"""
import math
from types import SimpleNamespace

import pytest

from stellarObjects import asteroidFieldData, cometData, nebulaData, physical_constants, program_constants
from stellarObjects import roguePlanetData, supernovaRemnantData
from stellarObjects.compactRemnant import BlackHole, NeutronStar, infer_black_hole_mass_class
from stellarObjects.config import SystemConfig
from stellarObjects.nebulaData import NEBULA_CLASS_LETTERS, REMNANT_CLASS_LETTERS, Nebula
from stellarObjects.planetData import Planet
from stellarObjects.roguePlanetData import RoguePlanet
from stellarObjects.starData import Star
from stellarObjects.systemData import StarSystem

from tests.fuzz_support import deterministic_entropy

PC = program_constants
CLASSES = PC.NEBULA_CLASSES
JUPITER = physical_constants.JUPITER_MASS_TO_KG
EARTH = physical_constants.EARTH_MASS_TO_KG


@pytest.fixture(autouse=True)
def _seeded():
    with deterministic_entropy(3601):
        yield


def cfg(**overrides):
    config = SystemConfig()
    for key, value in overrides.items():
        setattr(config, key, value)
    return config


# --- nebula ---------------------------------------------------------------

def test_nebula_and_remnant_letters_partition_the_class_table():
    assert set(NEBULA_CLASS_LETTERS).isdisjoint(REMNANT_CLASS_LETTERS)
    assert set(NEBULA_CLASS_LETTERS) | set(REMNANT_CLASS_LETTERS) == set(CLASSES)
    assert all(CLASSES[letter]["family"] == "supernova-remnant" for letter in REMNANT_CLASS_LETTERS)


def test_draw_in_range_is_uniform_from_zero_and_log_uniform_otherwise():
    for _ in range(100):
        assert 0.0 <= nebulaData._draw_in_range(0.0, 0.1) <= 0.1
        assert 10.0 <= nebulaData._draw_in_range(10.0, 1e6) <= 1e6


def test_mid_of_range_is_arithmetic_from_zero_and_geometric_otherwise():
    assert nebulaData._mid_of_range(0.0, 0.3) == pytest.approx(0.15)
    assert nebulaData._mid_of_range(10.0, 1000.0) == pytest.approx(100.0)
    assert nebulaData._mid_of_range(5.0, 5.0) == pytest.approx(5.0)


@pytest.mark.parametrize("letter", list(CLASSES))
def test_draw_and_typical_contents_stay_in_the_class_ranges(letter):
    data = CLASSES[letter]
    for contents in (nebulaData.draw_class_contents(letter), nebulaData.typical_class_contents(letter)):
        species, density, temperature, extinction = contents
        assert species == data["species"]
        assert data["density_range_cm3"][0] <= density <= data["density_range_cm3"][1]
        assert data["temperature_range_k"][0] <= temperature <= data["temperature_range_k"][1]
        assert data["extinction_range_av"][0] <= extinction <= data["extinction_range_av"][1]


def test_typical_contents_are_fixed():
    assert nebulaData.typical_class_contents("D") == nebulaData.typical_class_contents("D")


def test_choose_weighted_class_picks_from_the_given_letters():
    assert nebulaData.choose_weighted_class(("Q",), "test") == "Q"
    for _ in range(50):
        assert nebulaData.choose_weighted_class(("A", "N"), "test") in ("A", "N")


def test_choose_weighted_class_rejects_no_letters():
    with pytest.raises(IndexError):
        nebulaData.choose_weighted_class((), "test")


@pytest.mark.parametrize("family", list(PC.NEBULA_FAMILIES))
def test_inferred_nebula_class_stays_in_its_family(family):
    for radius in (0.01, 0.2, 1.0, 15.0, 120.0, 1e6):
        letter = nebulaData.infer_nebula_class(family, radius)
        assert CLASSES[letter]["family"] == family


def test_inferred_nebula_class_prefers_the_most_common_fitting_class():
    # 0.3 ly fits H (0.1-0.5, freq 1), K (0.1-2, 1) and L (0.2-3, 1) but
    # not J (0.5-3, freq 2); 1.0 ly fits J, which is the most common.
    assert nebulaData.infer_nebula_class("planetary", 1.0) == "J"
    assert nebulaData.infer_nebula_class("planetary", 0.3) in ("H", "K", "L")


def test_inferred_nebula_class_falls_back_to_the_most_common_in_family():
    # 1e6 ly fits nothing; N (freq 3) is the dark family's most common.
    assert nebulaData.infer_nebula_class("dark", 1e6) == "N"


def test_inferred_nebula_class_for_an_unknown_family_uses_any_nebula_class():
    letter = nebulaData.infer_nebula_class("no-such-family", 1.0)
    assert letter in NEBULA_CLASS_LETTERS


def test_nebula_constructor_errors():
    with pytest.raises(ValueError, match=r"nebula_class must be one of \[.*\], got 'R'"):
        Nebula(cfg(), nebula_class="R")
    with pytest.raises(ValueError, match=r"nebula_type must be one of \[.*\], got 'supernova-remnant'"):
        Nebula(cfg(), nebula_type="supernova-remnant")


@pytest.mark.parametrize("letter", NEBULA_CLASS_LETTERS)
def test_nebula_class_name_matches_its_letter(letter):
    nebula = Nebula(cfg(), nebula_class=letter)
    assert nebula.class_name == CLASSES[letter]["name"]
    assert nebula.nebula_type == CLASSES[letter]["family"]
    low, high = CLASSES[letter]["radius_range_ly"]
    assert low <= nebula.radius_ly <= high


# --- supernova remnant ----------------------------------------------------

@pytest.mark.parametrize("kind", ["neutron_star", "black_hole", None])
def test_a_type_ia_remnant_is_always_w(kind):
    assert supernovaRemnantData.remnant_classes_for("Type Ia", kind) == ("W",)


def test_core_collapse_remnant_classes_follow_the_core():
    ns = supernovaRemnantData.remnant_classes_for("Type II", "neutron_star")
    bh = supernovaRemnantData.remnant_classes_for("Type II", "black_hole")
    bare = supernovaRemnantData.remnant_classes_for("Type II", None)
    assert {"T", "U"} <= set(ns)
    assert not {"T", "U"} & set(bh) and not {"T", "U"} & set(bare)
    for classes in (ns, bh, bare):
        assert "W" not in classes
        assert set(classes) <= set(REMNANT_CLASS_LETTERS)


def test_unknown_core_kind_has_no_remnant_class():
    assert supernovaRemnantData.remnant_classes_for("Type II", "quark_star") == ()


@pytest.mark.parametrize("morphology, progenitor, age, expected", [
    ("plerion", "Type Ia", 100.0, "W"),
    ("plerion", "Type II", 100.0, "T"),
    ("composite", "Type Ib", 5000.0, "U"),
    ("shell", "Type II", 1000.0, "R"),
    ("shell", "Type II", 3000.0, "R"),
    ("shell", "Type II", 50000.0, "V"),
    ("shell", "Type II", 10000.0, "S"),
    ("shell", "Type II", 10.0, "S"),
    ("shell", "Type II", 1e9, "S"),
])
def test_infer_remnant_class(morphology, progenitor, age, expected):
    assert supernovaRemnantData.infer_remnant_class(morphology, progenitor, age) == expected


# --- black hole -----------------------------------------------------------

def test_black_hole_mass_class_boundaries():
    smbh_low = PC.BLACK_HOLE_SUPERMASSIVE_MASS_RANGE_SOLAR[0]
    imbh_low = PC.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR[0]
    assert infer_black_hole_mass_class(smbh_low) == "supermassive"
    assert infer_black_hole_mass_class(smbh_low * 0.999) == "intermediate"
    assert infer_black_hole_mass_class(imbh_low) == "intermediate"
    assert infer_black_hole_mass_class(imbh_low * 0.999) == "stellar"
    assert infer_black_hole_mass_class(0.0) == "stellar"


@pytest.mark.parametrize("mass_class", ["stellar", "intermediate", "Supermassive", ""])
def test_black_hole_rejects_an_unsupported_mass_class(mass_class):
    with pytest.raises(ValueError, match=rf"mass_class must be None or 'supermassive', got {mass_class!r}"):
        BlackHole(cfg(), mass_class=mass_class)


def test_black_hole_mass_class_label():
    labels = {"stellar": "Stellar-Mass", "intermediate": "Intermediate-Mass", "supermassive": "Supermassive"}
    hole = BlackHole(cfg())
    for mass_class, label in labels.items():
        hole.mass_class = mass_class
        assert hole.mass_class_label == label
    hole.mass_class = "bogus"
    with pytest.raises(KeyError):
        hole.mass_class_label


# --- rogue planet and interstellar comet ----------------------------------

def test_infer_rogue_mass_bin_boundaries():
    bins = PC.ROGUE_PLANET_MASS_BINS
    for name, (low, high, _rate) in bins.items():
        assert roguePlanetData.infer_rogue_mass_bin(low * EARTH * 1.0001) == name
        assert roguePlanetData.infer_rogue_mass_bin(high * EARTH * 0.9999) == name
    brown_low = PC.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER[0] * JUPITER
    assert roguePlanetData.infer_rogue_mass_bin(brown_low) == "brown-dwarf"
    assert roguePlanetData.infer_rogue_mass_bin(brown_low * 0.999) == "jupiter"
    assert roguePlanetData.infer_rogue_mass_bin(0.0) == "terrestrial"


@pytest.mark.parametrize("planet_type", ["t", "g"])
def test_rogue_planet_classes_are_rogue_flagged_and_typed(planet_type):
    classes = roguePlanetData.rogue_planet_classes(planet_type)
    assert classes
    for code in classes:
        assert PC.PLANET_CLASSES[code]["r"] and PC.PLANET_CLASSES[code]["type"] == planet_type


def test_rogue_planet_classes_for_an_unknown_type_are_empty():
    assert roguePlanetData.rogue_planet_classes("x") == []
    assert roguePlanetData.rogue_planet_class_candidates("x", 6000.0, EARTH) == []
    assert roguePlanetData.choose_rogue_planet_class("x", 6000.0, EARTH) is None
    assert roguePlanetData.default_rogue_planet_class("x", 6000.0, EARTH) is None


def test_rogue_candidates_fit_radius_and_mass():
    candidates = roguePlanetData.rogue_planet_class_candidates("t", 6000.0, EARTH)
    assert candidates
    for code in candidates:
        low, high = PC.PLANET_CLASSES[code]["radius_range"]
        assert low <= 6000.0 <= high


@pytest.mark.parametrize("radius_km", [1.0, 1e7])
def test_rogue_candidates_fall_back_to_the_nearest_class(radius_km):
    candidates = roguePlanetData.rogue_planet_class_candidates("t", radius_km, EARTH)
    assert len(candidates) == 1 and candidates[0] in roguePlanetData.rogue_planet_classes("t")


def test_choose_and_default_rogue_class():
    assert roguePlanetData.choose_rogue_planet_class("g", 70000.0, 50 * JUPITER, mass_bin="brown-dwarf") is None
    assert roguePlanetData.default_rogue_planet_class("g", 70000.0, 50 * JUPITER, mass_bin="brown-dwarf") is None
    candidates = roguePlanetData.rogue_planet_class_candidates("g", 30000.0, 1e26)
    for _ in range(20):
        assert roguePlanetData.choose_rogue_planet_class("g", 30000.0, 1e26) in candidates
    default = roguePlanetData.default_rogue_planet_class("g", 30000.0, 1e26)
    assert default in candidates
    assert default == roguePlanetData.default_rogue_planet_class("g", 30000.0, 1e26)
    best = max(PC.PLANET_CLASS_PROBABILITIES.get(code, 0.0) for code in candidates)
    assert PC.PLANET_CLASS_PROBABILITIES.get(default, 0.0) == best


@pytest.mark.parametrize("mass_bin", ["huge", "Jupiter", ""])
def test_rogue_planet_rejects_an_unknown_mass_bin(mass_bin):
    with pytest.raises(ValueError, match=rf"mass_bin must be one of .*, got {mass_bin!r}"):
        RoguePlanet(cfg(), mass_bin=mass_bin)


@pytest.mark.parametrize("composition, expected", [
    ([], "unknown composition"),
    (["water ice"], "water ice"),
    (["water ice", "dust"], "water ice and dust"),
    (["a", "b", "c"], "a, b, and c"),
])
def test_comet_composition_summary(composition, expected):
    assert roguePlanetData.format_comet_composition_summary(composition) == expected


def test_interstellar_comet_designation():
    assert roguePlanetData.interstellar_comet_designation("FE81", 3) == "I/FE81-3"
    assert roguePlanetData.interstellar_comet_designation("", 3) == "I/3"
    assert roguePlanetData.interstellar_comet_designation(None, 12) == "I/12"


# --- star-bound comet -----------------------------------------------------

def comet_stub(orbit_type, period):
    return SimpleNamespace(orbit_type=orbit_type, orbital_period_years=period)


@pytest.mark.parametrize("orbit_type, period, prefix", [
    ("elliptical", 5.0, "P"),
    ("elliptical", cometData.PERIODIC_COMET_MAX_PERIOD_YEARS - 1e-9, "P"),
    ("elliptical", cometData.PERIODIC_COMET_MAX_PERIOD_YEARS, "C"),
    ("elliptical", 1e5, "C"),
    ("elliptical", None, "C"),
    ("parabolic", None, "C"),
    ("parabolic", 5.0, "C"),
])
def test_comet_designation_prefix(orbit_type, period, prefix):
    assert cometData.comet_designation("Sol", 2, comet_stub(orbit_type, period)) == f"{prefix}/Sol-2"


@pytest.mark.parametrize("name, expected", [
    ("P/Old-3", "P/New-3"),
    ("C/Old A-12", "C/New A-12"),
    ("C/Older-1", None),
    ("P/Old-x", None),
    ("X/Old-1", None),
    ("P-Old-1", None),
    ("P/Old", None),
    ("", None),
    (None, None),
])
def test_rename_comet_designation(name, expected):
    assert cometData.rename_comet_designation(name, "Old", "New") == expected


def test_activity_chance_falls_with_perihelion_and_is_bounded():
    threshold = PC.COMET_ACTIVITY_PERIHELION_THRESHOLD_AU
    chances = [cometData._activity_chance(q) for q in (0.0, threshold / 2, threshold, threshold * 10)]
    assert chances == sorted(chances, reverse=True)
    assert chances[0] == PC.COMET_ACTIVITY_MAX_CHANCE
    assert chances[2] == chances[3] == pytest.approx(PC.COMET_ACTIVITY_MIN_CHANCE)
    assert cometData._activity_chance(-1.0) == PC.COMET_ACTIVITY_MAX_CHANCE


# --- asteroid field -------------------------------------------------------

def test_asteroid_field_size_digit_clamps():
    low, high = PC.ASTEROID_FIELD_SIZE_DIGIT_RANGE
    au = asteroidFieldData.AU_PER_LY
    assert asteroidFieldData.asteroid_field_size_digit(1e-9) == low
    assert asteroidFieldData.asteroid_field_size_digit(1e9) == high
    assert asteroidFieldData.asteroid_field_size_digit(1000 / au) == 3
    digits = [asteroidFieldData.asteroid_field_size_digit(r) for r in (1e-5, 1e-3, 1e-2, 0.1, 1.0, 10.0)]
    assert digits == sorted(digits)


@pytest.mark.parametrize("family", list(PC.ASTEROID_FIELD_COMPOSITIONS))
def test_composition_family_round_trips_through_its_letters(family):
    for letter in PC.ASTEROID_FIELD_COMPOSITIONS[family]["letters"].values():
        assert asteroidFieldData.composition_family_for_letter(letter) == family


@pytest.mark.parametrize("letter", ["Z", "", "a", "C3"])
def test_composition_family_for_an_unknown_letter_raises(letter):
    with pytest.raises(ValueError, match=rf"unknown asteroid field class letter {letter!r}"):
        asteroidFieldData.composition_family_for_letter(letter)


def test_asteroid_field_class_rejects_unknown_family_or_density():
    with pytest.raises(KeyError):
        asteroidFieldData.asteroid_field_class("plasma", "dense", 1.0)
    with pytest.raises(KeyError):
        asteroidFieldData.asteroid_field_class("icy", "very dense", 1.0)


def test_asteroid_field_designation():
    assert asteroidFieldData.asteroid_field_designation("C3", "FE81000A2B", 1) == "AF C3-FE81000A2B-01"
    assert asteroidFieldData.asteroid_field_designation("C3", None, 7) == "AF C3-07"
    assert asteroidFieldData.asteroid_field_designation("C3", "", 123) == "AF C3-123"


# --- get_table_properties -------------------------------------------------

def assert_all_strings(properties, allow_none=()):
    for key, value in properties.items():
        if key in allow_none and value is None:
            continue
        assert isinstance(value, str) and value, key


def test_star_table_properties():
    star = Star(cfg(STAR_TYPE="G2V"))
    props = star.get_table_properties()
    assert set(props) == {"type", "radius", "mass", "temp", "lum", "hab", "orbit", "loc"}
    assert_all_strings(props)
    assert props["loc"] == star.name and props["temp"].endswith(" K") and props["hab"].startswith("Between ")
    star.reflex_offset_x = 1e-4
    assert "pulled by its own planets" in star.get_table_properties()["wobble"]


def test_planet_table_properties():
    star = Star(cfg(STAR_TYPE="G2V"))
    planet = Planet(star.system_config, star, star.habitable_zone, star.habitable_zone[1] * 3,
                    planet_class="J", moon_count=0)
    props = planet.get_table_properties()
    assert set(props) == {"class", "distance", "period", "speed", "radius", "gravity"}
    assert props["class"] == "J" and props["gravity"].endswith(" g")
    assert_all_strings(props)
    planet.reflex_offset_y = 1e-6
    assert "pulled by its own moons" in planet.get_table_properties()["moon_wobble"]
    planet.planet_class = None
    assert planet.get_table_properties()["class"] is None


def test_close_binary_table_properties():
    system = StarSystem(cfg(STAR_TYPE="G2V", PLANETS=False, BINARY_SYSTEM=True, WIDE_BINARY=False))
    props = system.star.get_table_properties()
    assert {"type", "mass", "lum", "hab", "separation", "mutual_orbit", "orbit", "loc"} <= set(props)
    assert_all_strings(props)
    assert "per orbit" in props["mutual_orbit"]


def test_wide_binary_table_properties():
    system = StarSystem(cfg(STAR_TYPE="G2V", PLANETS=False, BINARY_SYSTEM=True, WIDE_BINARY=True))
    pair = system.wide_binary
    props = pair.get_table_properties()
    assert set(props) == {"separation", "eccentricity", "mutual_orbit", "wobble", "primary_limit", "secondary_limit"}
    assert_all_strings(props)
    assert float(props["eccentricity"]) == pytest.approx(pair.eccentricity, abs=5e-4)
    assert props["primary_limit"].startswith(f"{pair.primary.name}'s planetary system is stable out to ")
    assert "from the barycenter" in props["wobble"]


@pytest.mark.parametrize("mass_class", [None, "supermassive"])
def test_black_hole_table_properties(mass_class):
    hole = BlackHole(cfg(), mass_class=mass_class)
    props = hole.get_table_properties()
    assert set(props) == {"type", "mass", "event_horizon", "spin", "disk", "orbit", "loc"}
    assert_all_strings(props)
    assert props["type"] == f"{hole.mass_class_label} Black Hole"
    assert props["disk"] in ("Active accretion disk detected", "None detected")
    assert props["spin"].endswith("(dimensionless)")
    if mass_class == "supermassive":
        assert props["orbit"] == "None (the galaxy's center)"
    hole.reflex_offset_z = 1e-5
    assert "wobble" in hole.get_table_properties()


def test_neutron_star_table_properties():
    star = NeutronStar(cfg())
    props = star.get_table_properties()
    assert set(props) == {"type", "mass", "radius", "spin_period", "magnetic_field", "surface_temp", "orbit", "loc"}
    assert_all_strings(props)
    assert props["magnetic_field"].endswith(" G") and props["surface_temp"].endswith(" K")
    assert math.isfinite(float(props["magnetic_field"][:-2]))
    star.reflex_offset_x = 1e-5
    assert "wobble" in star.get_table_properties()
