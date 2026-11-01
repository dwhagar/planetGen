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

from planetgen.generation import star_population
from planetgen.physics import constants, imf, main_sequence, stellar_evolution
from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
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
                           "age_gy", "lifespan_gy", "initial_mass_sol", "phase_end_age_gy", "metallicity_fe_h"}
    assert params["temperature_k"] % 100 == 0
    assert params["mass_kg"] == pytest.approx(constants.SOLAR_MASS_TO_KG, rel=1e-3)      # MIST: the Sun has lost 1e-4 by now
    assert params["type"].endswith("Main Sequence Star")


def test_white_dwarf_params_keep_the_radius_of_the_cooling_sequences():
    state = stellar_evolution.evolve_star(2.0, 9.0)
    assert state["yerkes_class"] == "VII"
    params = stellar_evolution.star_params(2.0, 9.0, state)
    assert params["radius_km"] == state["radius_km"]
    assert math.isinf(params["lifespan_gy"]) and math.isinf(params["phase_end_age_gy"])


@pytest.mark.parametrize("mass, ms_fraction, yerkes", [
    (1.0, 0.1, "V"), (1.0, 1.04, "IV"), (1.0, 1.09, "III"), (1.0, 1.25, "VII"),
    (6.0, 1.05, "II"), (6.0, 1.2, "VII"),
])
def test_evolve_star_phases_follow_age(mass, ms_fraction, yerkes):
    """Phase boundaries are multiples of the main-sequence lifetime: a solar-mass star leaves the subgiant phase at 1.05
    of it and its giant branch ends at 1.12 (GEN.502, [R2055]); a 6 solar mass star is a bright giant
    burning helium in its core until 1.13 of it (the MIST track, GEN.935) and a white dwarf after."""
    age = stellar_evolution.main_sequence_lifetime_gy(mass) * ms_fraction
    rng = random.Random(7)
    state = stellar_evolution.evolve_star(mass, age, rng, fe_h=0.0)
    assert state["yerkes_class"] == yerkes
    assert age < state["phase_end_gy"]


def test_evolve_star_supergiant_luminosity_follows_the_track_from_the_end_of_the_main_sequence():
    mass = 20.0
    t_ms = stellar_evolution.main_sequence_lifetime_gy(mass)
    rng = random.Random(11)
    floor = main_sequence.end_state(mass, 0.0).luminosity_sol * 0.9
    for _ in range(20):
        state = stellar_evolution.evolve_star(mass, t_ms * 1.01, rng, min_luminosity_sol=floor, fe_h=0.0)
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


@pytest.mark.parametrize("population", ["mystery", "", "Young"])
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


def test_main_sequence_radius_is_increasing_and_near_solar_for_the_sun_at_its_age():
    radii = [main_sequence.state(m, 0.3 * stellar_evolution.main_sequence_lifetime_gy(m)).radius_sol for m in MASSES]
    assert radii == sorted(radii)
    assert main_sequence.state(1.0, 4.57).radius_sol == pytest.approx(1.0, rel=0.1)


def test_main_sequence_luminosity_is_increasing_and_positive():
    lums = [main_sequence.zams_state(m).luminosity_sol for m in MASSES]
    assert all(lum > 0 for lum in lums)
    assert lums == sorted(lums)


def test_a_main_sequence_star_brightens_with_age_and_ends_on_its_end_row():
    for mass in (0.5, 1.0, 3.0, 10.0):
        t_ms = stellar_evolution.main_sequence_lifetime_gy(mass)
        ages = [t_ms * share for share in (0.25, 0.5, 0.75, 0.999)]
        lums = [stellar_evolution.evolve_star(mass, age, fe_h=0.0)["luminosity_sol"] for age in ages]
        assert lums == sorted(lums), mass
        assert lums[-1] == pytest.approx(main_sequence.end_state(mass, 0.0).luminosity_sol, rel=0.15)


def test_the_age_a_main_sequence_star_reaches_a_brightness_is_where_its_row_does():
    mass, t_ms = 2.0, stellar_evolution.main_sequence_lifetime_gy(2.0)
    first, last = (main_sequence.state(mass, age).luminosity_sol for age in (0.0, t_ms))
    target = (first * last) ** 0.5
    age = stellar_evolution.main_sequence_bright_from_gy(mass, target, t_ms)
    assert 0.0 < age < t_ms
    assert main_sequence.state(mass, age).luminosity_sol == pytest.approx(target, rel=1e-3)
    assert stellar_evolution.main_sequence_bright_from_gy(mass, first / 2, t_ms) == 0.0
    assert stellar_evolution.main_sequence_bright_from_gy(mass, last * 2, t_ms) is None


def test_main_sequence_lifetime_falls_with_mass():
    lifetimes = [stellar_evolution.main_sequence_lifetime_gy(m) for m in MASSES]
    assert lifetimes == sorted(lifetimes, reverse=True)
    assert stellar_evolution.main_sequence_lifetime_gy(1.0) == pytest.approx(10.0, rel=0.05)  # GEN.467: 10.2 Gy


def test_effective_temperature_is_solar_for_the_sun_and_monotonic():
    assert stellar_evolution.effective_temperature_k(1.0, 1.0) == PC.SUN_EFFECTIVE_TEMPERATURE_K
    assert stellar_evolution.effective_temperature_k(10.0, 1.0) > PC.SUN_EFFECTIVE_TEMPERATURE_K
    assert stellar_evolution.effective_temperature_k(1.0, 10.0) < PC.SUN_EFFECTIVE_TEMPERATURE_K
    assert stellar_evolution.effective_temperature_k(0.0, 1.0) == 0.0


def test_effective_temperature_rejects_a_zero_radius():
    with pytest.raises(ZeroDivisionError):
        stellar_evolution.effective_temperature_k(1.0, 0.0)


def test_spectral_letter_clamps_outside_the_table():
    assert stellar_evolution.spectral_letter_and_subclass(1.0e6) == ("O", 0)
    assert stellar_evolution.spectral_letter_and_subclass(1.0) == ("M", constants.SUBCLASS_MAX_VALUE)


@pytest.mark.parametrize("temperature, expected", [
    (44900, ("O", 3)), (33300, ("O", 9)), (20600, ("B", 2)), (9700, ("A", 0)), (6550, ("F", 5)), (5770, ("G", 2)),
    (5270, ("K", 0)), (3930, ("K", 9)), (3850, ("M", 0)), (3060, ("M", 5)), (2380, ("M", 9)),
])
def test_spectral_types_are_the_mamajek_sequence(temperature, expected):
    """GEN.943: the Sun is G2, K0 is 5,270 K and the K/M edge is near 3,900 K (the types of the Pecaut and Mamajek table)."""
    assert stellar_evolution.spectral_letter_and_subclass(temperature) == expected


def test_the_k_m_boundary_is_near_3900_k():
    assert stellar_evolution.spectral_letter_and_subclass(3900)[0] == "K"
    assert stellar_evolution.spectral_letter_and_subclass(3880)[0] == "M"


def test_spectral_subclass_runs_hot_to_cool_within_a_letter():
    previous = (-1, -1)
    for kelvin in range(5980, 5330, -10):
        letter, subclass = stellar_evolution.spectral_letter_and_subclass(kelvin)
        assert letter == "G"
        assert subclass >= previous[1] or previous[0] != letter
        previous = (letter, subclass)
    assert stellar_evolution.spectral_letter_and_subclass(5980) == ("G", 0)
    assert stellar_evolution.spectral_letter_and_subclass(5330) == ("G", 9)


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
        imf.sample(low, high)


def test_imf_draw_stays_inside_its_truncation():
    rng = random.Random(5)
    for _ in range(200):
        assert 0.4 <= imf.sample(0.4, 0.6, rng) <= 0.6


def test_no_living_star_raises_after_max_redraws(monkeypatch):
    monkeypatch.setattr(stellar_evolution, "evolve_star", lambda *args, **kwargs: None)
    with pytest.raises(ValueError, match=f"no living star in {PC.STAR_MODEL_MAX_REDRAWS} draws"):
        stellar_evolution.sample_living_star(rng=random.Random(1))


@pytest.mark.parametrize("threshold", [math.nan, math.inf, "10", None])
def test_bright_threshold_must_be_finite(threshold):
    with pytest.raises(ValueError, match="luminosity threshold must be a finite number"):
        star_population.check_threshold(threshold)


def test_bright_threshold_must_clear_the_brightest_white_dwarf():
    too_low = PC.WD_LUMINOSITY_RANGE_SOL[1] / 2
    with pytest.raises(ValueError, match="must be at least the brightest white dwarf"):
        star_population.check_threshold(too_low)


def test_bright_band_must_not_be_empty():
    with pytest.raises(ValueError, match=r"empty luminosity band \[100, 50\) Lsun"):
        star_population.sample_bright_stars(1, 100.0, max_luminosity_sol=50.0)


# --- GEN.200: the map luminosity of a mass ---------------------------------

@pytest.mark.parametrize("mass_sol, expected_sol", [
    (0.1, 0.23 * 0.1 ** 2.3),
    (0.43, 0.43 ** 4),
    (1.0, 1.0),
    (2.0, 1.5 * 2.0 ** 3.5),
    (20.0, 1.5 * 20.0 ** 3.5),  # the linear band starts where the 1.5 M^3.5 band ends
    (4.0e6, 1.5 * 20.0 ** 2.5 * 4.0e6),
])
def test_map_luminosity_follows_the_tables_bands(mass_sol, expected_sol):
    assert stellar_evolution.map_luminosity_sol(mass_sol) == pytest.approx(expected_sol, rel=1e-12)


def test_map_luminosity_matches_the_figures_in_the_item():
    assert stellar_evolution.map_luminosity_sol(20.0) == pytest.approx(53_700, rel=0.01)
    assert stellar_evolution.map_luminosity_sol(4.0e6) == pytest.approx(1.0e10, rel=0.1)
    assert 1.5 * 20.0 ** 2.5 == pytest.approx(2683, rel=1e-3)


def test_map_luminosity_rises_with_mass_over_the_whole_range():
    masses = [10 ** (exponent / 50) for exponent in range(-100, 400)]  # 0.01 to 1e8 Msun
    values = [stellar_evolution.map_luminosity_sol(m) for m in masses]
    assert all(later > earlier for earlier, later in zip(values, values[1:]))


@pytest.mark.parametrize("join_sol, step", [(0.43, 0.036), (2.0, 0.062), (20.0, 0.0)])
def test_map_luminosity_steps_only_a_little_at_the_joins(join_sol, step):
    below = stellar_evolution.map_luminosity_sol(join_sol * (1 - 1e-9))
    above = stellar_evolution.map_luminosity_sol(join_sol)
    assert above / below - 1 == pytest.approx(step, abs=0.005)


def test_map_luminosity_is_within_a_factor_of_the_star_generators_main_sequence():
    # The map's table is its own (tuning.py); the generators read MIST (GEN.652). Over a main sequence the two keep
    # within a factor of four of the star's luminosity at its middle age.
    for mass in (0.2, 0.5, 1.0, 3.0, 10.0):
        middle = main_sequence.state(mass, 0.5 * stellar_evolution.main_sequence_lifetime_gy(mass)).luminosity_sol
        assert 0.25 < stellar_evolution.map_luminosity_sol(mass) / middle < 4.0, mass


# --- planetary-nebula central stars (GEN.149) ------------------------------

@pytest.mark.parametrize("star_type", list(tuning.PLANETARY_NEBULA_CENTRAL_STAR_TYPES))
def test_hot_white_dwarf_central_stars_are_the_models_hot_white_dwarfs(star_type):
    """A central star is a white dwarf the star model makes at that temperature (GEN.941): a remnant of 0.5 to 1.3 Msun
    whose luminosity and radius follow from its cooling and mass, as hot as its letter says."""
    for _ in range(20):
        star = make_star(star_type)
        mass_sol = star.mass / constants.SOLAR_MASS_TO_KG
        assert star.yerkes_class == "VII"
        assert 0.5 <= mass_sol <= 1.3
        assert star.temperature >= constants.COOL_WHITE_DWARF_BELOW_K
        assert star.lifespan == float("inf")


def test_hot_white_dwarf_temperature_follows_the_subclass():
    assert make_star("O3VII").temperature > make_star("O9VII").temperature >= make_star("B0VII").temperature
    assert make_star("B0VII").temperature >= constants.TEMP_RANGES["B"][1] * 0.9


def test_cool_white_dwarfs_keep_their_own_mass_and_temperature_ranges():
    star = make_star("M2VII")
    assert star.mass / constants.SOLAR_MASS_TO_KG <= 1.3
    assert star.temperature <= constants.TEMP_RANGES["M"][1]


def test_a_phenomenon_is_drawn_as_bright_as_a_star_of_its_mass():
    """MAP.208: a 10 Msun black hole about as bright as a 10 Msun star, a little dimmer; a neutron star equal."""
    star = stellar_evolution.map_luminosity_sol(10.0)
    black_hole = stellar_evolution.phenomenon_map_luminosity_sol("black_hole", 10.0)
    assert black_hole == pytest.approx(star * tuning.BLACK_HOLE_MAP_BRIGHTNESS)
    assert black_hole < star
    assert stellar_evolution.phenomenon_map_luminosity_sol("neutron_star", 1.4) == pytest.approx(
        stellar_evolution.map_luminosity_sol(1.4))
    assert stellar_evolution.phenomenon_map_luminosity_sol("quasar", 4.0e6) == pytest.approx(
        stellar_evolution.map_luminosity_sol(4.0e6))


@pytest.mark.parametrize("mass", [None, 0, -1.0])
def test_a_massless_marker_is_drawn_as_a_one_solar_mass_star(mass):
    assert stellar_evolution.phenomenon_map_luminosity_sol("nebula", mass) == pytest.approx(
        stellar_evolution.map_luminosity_sol(tuning.MASSLESS_MARKER_MASS_SOL))
    assert tuning.MASSLESS_MARKER_MASS_SOL == 1.0
    assert stellar_evolution.phenomenon_map_luminosity_sol("black_hole", mass) < stellar_evolution.phenomenon_map_luminosity_sol("quasar", mass)
