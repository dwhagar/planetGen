"""
Star and evolution helpers (TODO TEST.35).

Direct tests of `Star.calculate_heliosphere`, the population-model star
path (`Star._generate_from_population_model` / `stellar_evolution.evolve_star`
/ `star_params`), the age windows (`age_window_gy`,
`population_age_range_gy`, `sample_star_age_gy`), the main-sequence radius,
luminosity and temperature helpers, the white dwarf radius, the Yerkes
class helpers, and the raise messages around them. Behavior only:
reference values live in the physics reference gate.
"""
import math
import random

import pytest

from stellarObjects import starData, stellarPopulation
from planetgen.physics import constants, stellar_evolution
from planetgen import tuning
from stellarObjects.config import SystemConfig
from stellarObjects.starData import Star
from planetgen.physics.stellar_evolution import YERKES_CLASS_NAMES

from tests.fuzz_support import deterministic_entropy

PC = tuning


@pytest.fixture(autouse=True)
def _seeded():
    with deterministic_entropy(3501):
        yield


def make_star(star_type=None, **overrides):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return Star(cfg)


# --- calculate_heliosphere ------------------------------------------------

@pytest.mark.parametrize("star_type", ["G2V", "O5V", "B3V", "M5V", "K0III", "B0IA", "G2VII", "F5IV"])
def test_heliosphere_is_a_positive_finite_au_value(star_type):
    star = make_star(star_type)
    radius_au = star.calculate_heliosphere()
    assert isinstance(radius_au, float)
    assert math.isfinite(radius_au) and radius_au > 0


def test_heliosphere_delegates_to_the_static_helper():
    star = make_star("K1V")
    assert star.calculate_heliosphere() == Star._calculate_heliosphere_radius_static(
        star.mass, star.luminosity, star.radius, star.type, star.yerkes_class)


def test_heliosphere_rejects_zero_luminosity_cleanly():
    """The cool-dwarf wind scales with a negative power of luminosity, so
    zero is out of domain: a clean ValueError, not a ZeroDivisionError."""
    with pytest.raises(ValueError, match="_calculate_heliosphere_radius_static: input out of range"):
        Star._calculate_heliosphere_radius_static(constants.SOLAR_MASS_TO_KG, 0.0, 696000.0, "G2V", "V")


def test_heliosphere_grows_with_a_stronger_wind():
    """Within the evolved tier, more luminosity means more mass loss."""
    mass, radius = 5 * constants.SOLAR_MASS_TO_KG, 3e7
    sizes = [Star._calculate_heliosphere_radius_static(mass, lum * constants.SOLAR_LUMINOSITY, radius,
                                                       "K0III", "III") for lum in (10, 100, 1000)]
    assert sizes == sorted(sizes) and sizes[0] < sizes[-1]


def test_hot_dwarf_wind_beats_a_cool_dwarf_of_the_same_size():
    """The O/B radiation-driven tier applies only to spectral class O or B."""
    args = (2 * constants.SOLAR_MASS_TO_KG, 20 * constants.SOLAR_LUMINOSITY, 1.5e6)
    hot = Star._calculate_heliosphere_radius_static(*args, "B5V", "V")
    cool = Star._calculate_heliosphere_radius_static(*args, "A5V", "V")
    assert hot > cool


def test_hypergiant_wind_is_faster_than_a_supergiants():
    args = (30 * constants.SOLAR_MASS_TO_KG, 3e5 * constants.SOLAR_LUMINOSITY, 5e8)
    assert (Star._calculate_heliosphere_radius_static(*args, "B0IA", "0")
            > Star._calculate_heliosphere_radius_static(*args, "B0IA", "IA"))


def test_empty_star_type_falls_to_the_cool_dwarf_tier():
    args = (constants.SOLAR_MASS_TO_KG, constants.SOLAR_LUMINOSITY, 696000.0)
    assert (Star._calculate_heliosphere_radius_static(*args, "", "V")
            == Star._calculate_heliosphere_radius_static(*args, "G2V", "V"))


@pytest.mark.parametrize("arg_index", [0, 1, 2])
@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_heliosphere_rejects_non_finite_inputs(arg_index, bad):
    args = [constants.SOLAR_MASS_TO_KG, constants.SOLAR_LUMINOSITY, 696000.0]
    args[arg_index] = bad
    with pytest.raises(ValueError, match="must be a finite number"):
        Star._calculate_heliosphere_radius_static(*args, "G2V", "V")


def test_heliosphere_rejects_a_zero_radius():
    with pytest.raises(ValueError, match="_calculate_heliosphere_radius_static: input out of range"):
        Star._calculate_heliosphere_radius_static(constants.SOLAR_MASS_TO_KG,
                                                  constants.SOLAR_LUMINOSITY, 0.0, "G2V", "V")


# --- population-model star path -------------------------------------------

def test_random_star_comes_from_the_population_model():
    for _ in range(10):
        star = make_star()
        assert star.initial_mass_sol is not None
        assert PC.IMF_BREAKS_SOL[0] <= star.initial_mass_sol <= PC.IMF_BREAKS_SOL[-1]
        assert star.yerkes_class in YERKES_CLASS_NAMES
        assert YERKES_CLASS_NAMES[star.yerkes_class] in star.type
        assert star.age <= star.phase_end_age_gy
        assert star.mass > 0 and star.radius > 0 and star.luminosity > 0


def test_population_model_honours_a_given_mass_and_age():
    star = make_star()
    star._generate_from_population_model(initial_mass_sol=1.0, age_gy=4.6)
    assert star.initial_mass_sol == 1.0
    assert star.yerkes_class == "V"
    assert star._model_age_lifespan[0] == 4.6
    assert star.type.startswith("G")


def test_population_model_clamps_a_mass_below_the_imf_floor():
    star = make_star()
    star._generate_from_population_model(initial_mass_sol=0.001, age_gy=1.0)
    assert star.initial_mass_sol == PC.IMF_BREAKS_SOL[0]


def test_population_model_redraws_a_collapsed_companion():
    """A given mass that has already collapsed (30 Msun at 5 Gy) falls
    back to a fresh living star rather than failing."""
    assert stellar_evolution.evolve_star(30.0, 5.0) is None
    star = make_star()
    star._generate_from_population_model(initial_mass_sol=30.0, age_gy=5.0)
    assert star.yerkes_class not in ("BH", "NS")
    assert star.initial_mass_sol is not None


def test_population_model_respects_a_luminosity_cap():
    cap = 2.0
    for _ in range(10):
        star = make_star(MAX_STAR_LUMINOSITY_SOL=cap)
        assert star.yerkes_class == "VII" or star.luminosity / constants.SOLAR_LUMINOSITY < cap


def test_star_params_shape_and_units():
    state = stellar_evolution.evolve_star(1.0, 4.6)
    params = stellar_evolution.star_params(1.0, 4.6, state)
    assert set(params) == {"type", "yerkes_class", "mass_kg", "radius_km", "temperature_k", "luminosity_w",
                           "age_gy", "lifespan_gy", "initial_mass_sol", "phase_end_age_gy"}
    assert params["temperature_k"] % 100 == 0
    assert params["mass_kg"] == pytest.approx(constants.SOLAR_MASS_TO_KG)
    assert params["type"].endswith("Main Sequence Star")


def test_white_dwarf_params_keep_the_mass_radius_relation():
    state = stellar_evolution.evolve_star(2.0, 9.0)
    assert state["yerkes_class"] == "VII"
    params = stellar_evolution.star_params(2.0, 9.0, state)
    assert params["radius_km"] == stellar_evolution.white_dwarf_radius_km(state["mass_sol"])
    assert math.isinf(params["lifespan_gy"]) and math.isinf(params["phase_end_age_gy"])


@pytest.mark.parametrize("mass, ms_fraction, yerkes", [
    (1.0, 0.1, "V"), (1.0, 1.04, "IV"), (1.0, 1.15, "III"), (1.0, 1.25, "VII"),
    (6.0, 1.15, "II"), (6.0, 1.25, "VII"),
])
def test_evolve_star_phases_follow_age(mass, ms_fraction, yerkes):
    """Phase boundaries are multiples of the main-sequence lifetime."""
    age = stellar_evolution.main_sequence_lifetime_gy(mass) * ms_fraction
    rng = random.Random(7)
    state = stellar_evolution.evolve_star(mass, age, rng)
    assert state["yerkes_class"] == yerkes
    assert age < state["phase_end_gy"]


def test_evolve_star_supergiant_luminosity_honours_a_minimum():
    mass = 20.0
    t_ms = stellar_evolution.main_sequence_lifetime_gy(mass)
    rng = random.Random(11)
    growth_low, growth_high = PC.SUPERGIANT_LUMINOSITY_GROWTH_RANGE
    floor = stellar_evolution.main_sequence_luminosity_sol(mass) * (growth_low + growth_high) / 2
    for _ in range(20):
        state = stellar_evolution.evolve_star(mass, t_ms * 1.05, rng, min_luminosity_sol=floor)
        assert state["luminosity_sol"] >= floor
        assert state["yerkes_class"] in ("IB", "IAB", "IA", "0")


# --- age windows ----------------------------------------------------------

def test_age_window_unbiased_keeps_the_whole_range():
    assert stellar_evolution.age_window_gy(2.0, 8.0) == (2.0, 8.0)
    assert stellar_evolution.age_window_gy(2.0, 8.0, "anything else") == (2.0, 8.0)


def test_age_window_young_keeps_the_first_part():
    low, high = stellar_evolution.age_window_gy(2.0, 8.0, "young")
    assert low == 2.0
    assert high == pytest.approx(2.0 + 6.0 * PC.YOUNG_STAR_AGE_LIFESPAN_RATIO)


def test_age_window_old_drops_the_first_part():
    low, high = stellar_evolution.age_window_gy(2.0, 8.0, "old")
    assert high == 8.0
    assert low == pytest.approx(2.0 + 6.0 * PC.OLD_STAR_AGE_LIFESPAN_RATIO)


def test_age_window_of_an_empty_range_stays_empty():
    for bias in (None, "young", "old"):
        assert stellar_evolution.age_window_gy(3.0, 3.0, bias) == (3.0, 3.0)


def test_population_age_ranges():
    assert stellar_evolution.population_age_range_gy() == PC.STAR_FORMATION_AGE_RANGE_GY
    for name, window in PC.STELLAR_POPULATION_AGE_RANGES_GY.items():
        assert stellar_evolution.population_age_range_gy(name) == window


@pytest.mark.parametrize("population", ["halo", "", "Young"])
def test_unknown_population_raises(population):
    with pytest.raises(ValueError, match=rf"unknown stellar population {population!r}; expected one of "):
        stellar_evolution.population_age_range_gy(population)


@pytest.mark.parametrize("population", [None, "young", "intermediate", "old", "bulge"])
@pytest.mark.parametrize("bias", [None, "young", "old"])
def test_sampled_ages_stay_in_their_window(population, bias):
    rng = random.Random(3)
    low, high = stellar_evolution.age_window_gy(*stellar_evolution.population_age_range_gy(population), bias)
    for _ in range(50):
        assert low <= stellar_evolution.sample_star_age_gy(bias, rng, population) <= high


# --- radius, luminosity and temperature helpers ---------------------------

MASSES = [0.08, 0.2, 0.5, 0.9, 0.999, 1.0, 1.5, 3.0, 10.0, 40.0, 150.0]


def test_main_sequence_radius_is_increasing_and_solar_at_one():
    radii = [stellar_evolution.main_sequence_radius_sol(m) for m in MASSES]
    assert radii == sorted(radii)
    assert stellar_evolution.main_sequence_radius_sol(1.0) == 1.0
    assert stellar_evolution.main_sequence_radius_sol(1 - 1e-9) == pytest.approx(1.0, rel=1e-6)


def test_main_sequence_luminosity_is_increasing_and_positive():
    lums = [stellar_evolution.main_sequence_luminosity_sol(m) for m in MASSES]
    assert all(lum > 0 for lum in lums)
    assert lums == sorted(lums)


def test_main_sequence_lifetime_falls_with_mass():
    lifetimes = [stellar_evolution.main_sequence_lifetime_gy(m) for m in MASSES]
    assert lifetimes == sorted(lifetimes, reverse=True)
    assert stellar_evolution.main_sequence_lifetime_gy(1.0) == PC.SOLAR_MS_LIFESPAN_GY


def test_effective_temperature_is_solar_for_the_sun_and_monotonic():
    assert stellar_evolution.effective_temperature_k(1.0, 1.0) == PC.SUN_EFFECTIVE_TEMPERATURE_K
    assert stellar_evolution.effective_temperature_k(10.0, 1.0) > PC.SUN_EFFECTIVE_TEMPERATURE_K
    assert stellar_evolution.effective_temperature_k(1.0, 10.0) < PC.SUN_EFFECTIVE_TEMPERATURE_K
    assert stellar_evolution.effective_temperature_k(0.0, 1.0) == 0.0


def test_effective_temperature_rejects_a_zero_radius():
    with pytest.raises(ZeroDivisionError):
        stellar_evolution.effective_temperature_k(1.0, 0.0)


def test_spectral_letter_clamps_outside_the_table():
    ranges = constants.TEMP_RANGES
    assert stellar_evolution.spectral_letter_and_subclass(ranges["O"][1] * 10) == ("O", 0)
    assert stellar_evolution.spectral_letter_and_subclass(1.0) == ("M", constants.SUBCLASS_MAX_VALUE)


def test_spectral_subclass_runs_hot_to_cool_within_a_letter():
    low, high = constants.TEMP_RANGES["G"]
    hot = stellar_evolution.spectral_letter_and_subclass(high - 1)
    cool = stellar_evolution.spectral_letter_and_subclass(low + 1)
    assert hot[0] == cool[0] == "G"
    assert hot[1] == 0 and cool[1] == constants.SUBCLASS_MAX_VALUE


# --- white dwarf radius ---------------------------------------------------

def test_white_dwarf_radius_shrinks_with_mass():
    masses = [0.17, 0.4, 0.6, 1.0, 1.38]
    radii = [stellar_evolution.white_dwarf_radius_km(m) for m in masses]
    assert all(r > 0 for r in radii)
    assert radii == sorted(radii, reverse=True)
    assert stellar_evolution.white_dwarf_radius_km(1.0) == constants.WHITE_DWARF_BASE_RADIUS_KM


def test_white_dwarf_radius_is_earth_sized_across_the_mass_range():
    for mass in PC.WD_MASS_RANGE_SOL:
        assert 1000 < stellar_evolution.white_dwarf_radius_km(mass) < 20000


def test_white_dwarf_radius_rejects_zero_mass():
    with pytest.raises(ZeroDivisionError):
        stellar_evolution.white_dwarf_radius_km(0.0)


# --- Yerkes class ---------------------------------------------------------

def test_supergiant_yerkes_class_below_every_threshold_is_ib():
    lowest = min(PC.SUPERGIANT_YERKES_THRESHOLDS_SOL.values())
    assert stellar_evolution._supergiant_yerkes_class(lowest * 0.999) == "IB"
    assert stellar_evolution._supergiant_yerkes_class(0.0) == "IB"


def test_supergiant_yerkes_class_at_each_threshold():
    for name, threshold in PC.SUPERGIANT_YERKES_THRESHOLDS_SOL.items():
        assert stellar_evolution._supergiant_yerkes_class(threshold) == name


def test_supergiant_yerkes_class_brightens_monotonically():
    order = ["IB", *sorted(PC.SUPERGIANT_YERKES_THRESHOLDS_SOL, key=PC.SUPERGIANT_YERKES_THRESHOLDS_SOL.get)]
    seen = [stellar_evolution._supergiant_yerkes_class(10 ** (k / 4)) for k in range(0, 40)]
    ranks = [order.index(cls) for cls in seen]
    assert ranks == sorted(ranks)
    assert seen[-1] == "0"


def test_yerkes_names_cover_every_parseable_class_and_alias_d():
    for yerkes in ("IA+", "IAB", "VII", "III", "IA", "IB", "II", "IV", "VI", "0", "V", "D"):
        assert yerkes in YERKES_CLASS_NAMES
    assert YERKES_CLASS_NAMES["D"] == YERKES_CLASS_NAMES["VII"] == "White Dwarf"


@pytest.mark.parametrize("star_type, yerkes", [("G2V", "V"), ("K0III", "III"), ("B0IA", "IA"), ("G2VII", "VII"),
                                                ("F5IV", "IV"), ("g2v", "V")])
def test_forced_type_sets_the_yerkes_class_and_name(star_type, yerkes):
    star = make_star(star_type)
    assert star.yerkes_class == yerkes
    assert YERKES_CLASS_NAMES[yerkes] in star.type


@pytest.mark.parametrize("star_type", ["DA", "G2", "X2V", "G10V", "G2Vjunk", "G2 V", "2GV"])
def test_bad_star_type_raises(star_type):
    with pytest.raises(ValueError, match=r"Invalid star type format. Expected format is e.g., G2V."):
        make_star(star_type)


# --- other raise messages -------------------------------------------------

@pytest.mark.parametrize("low, high", [(1.0, 1.0), (5.0, 2.0), (200.0, None), (None, 0.01)])
def test_imf_empty_range_raises(low, high):
    with pytest.raises(ValueError, match=r"empty IMF range \["):
        stellar_evolution.sample_imf_mass_sol(low, high)


def test_imf_draw_stays_inside_its_truncation():
    rng = random.Random(5)
    for _ in range(200):
        assert 0.4 <= stellar_evolution.sample_imf_mass_sol(0.4, 0.6, rng) <= 0.6


def test_no_living_star_raises_after_max_redraws(monkeypatch):
    monkeypatch.setattr(stellar_evolution, "evolve_star", lambda *args, **kwargs: None)
    with pytest.raises(ValueError, match=f"no living star in {PC.STAR_MODEL_MAX_REDRAWS} draws"):
        stellar_evolution.sample_living_star(rng=random.Random(1))


def test_evolved_mass_sampler_raises_for_a_range_too_light_to_have_evolved():
    with pytest.raises(ValueError, match=r"Could not sample a mass in \[0.1, 0.2\] Msun"):
        starData._sample_evolved_star_mass_sol(0.1, 0.2)


def test_evolved_mass_sampler_accepts_only_masses_that_have_evolved():
    for _ in range(50):
        mass = starData._sample_evolved_star_mass_sol(0.5, 3.0)
        assert 0.5 <= mass <= 3.0
        assert stellar_evolution.main_sequence_lifetime_gy(mass) <= PC.UNIVERSE_AGE_GY


@pytest.mark.parametrize("threshold", [math.nan, math.inf, "10", None])
def test_bright_threshold_must_be_finite(threshold):
    with pytest.raises(ValueError, match="luminosity threshold must be a finite number"):
        stellarPopulation._check_threshold(threshold)


def test_bright_threshold_must_clear_the_brightest_white_dwarf():
    too_low = PC.WD_LUMINOSITY_RANGE_SOL[1] / 2
    with pytest.raises(ValueError, match="must be at least the brightest white dwarf"):
        stellarPopulation._check_threshold(too_low)


def test_bright_band_must_not_be_empty():
    with pytest.raises(ValueError, match=r"empty luminosity band \[100, 50\) Lsun"):
        stellarPopulation.sample_bright_stars(1, 100.0, max_luminosity_sol=50.0)
