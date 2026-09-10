"""
Exhaustive star-type matrix regression tests.

Builds a `Star` (no planets) for every (spectral class x subclass x Yerkes
class) combination the generator accepts and checks age/lifespan, the
gravitational sphere of influence (Hill/galactic-tidal radius), and the
heliosphere against physically-grounded sanity bounds.

This is the test that originally caught two real bugs: the hypergiant wind
constant was ~9 orders of magnitude too high (heliosphere radii of ~1e11 AU,
bigger than the galaxy) and the O/B dwarf wind constant was ~3-4 orders of
magnitude too low. It also caught giants/supergiants/subgiants/hypergiants
drawing their age from their *current temperature letter's* main-sequence
lifespan table (e.g. a red giant reported as hundreds of billions of years
old) instead of their own generated mass. See stellarObjects/starData.py and
physical_constants.py for the fixes.

Run with: pytest tests/test_star_matrix.py
"""
import math

import pytest

from stellarObjects.config import SystemConfig
from stellarObjects.starData import Star, _sample_evolved_star_mass_sol
from stellarObjects import physical_constants as pc
from stellarObjects import program_constants as prog_c

SPECTRAL_CLASSES = ["O", "B", "A", "F", "G", "K", "M"]
SUBCLASSES = list(range(10))
YERKES_CLASSES = ["0", "IA+", "IA", "IAB", "IB", "II", "III", "IV", "V", "VI", "VII"]

# A handful of independent trials per combo, since luminosity/mass/age are
# randomized within physically-valid ranges for that type.
TRIALS = 5

ALL_STAR_TYPES = [
    f"{spec}{sub}{yk}"
    for spec in SPECTRAL_CLASSES
    for sub in SUBCLASSES
    for yk in YERKES_CLASSES
]


def make_star(star_type):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    return Star(cfg)


@pytest.mark.parametrize("star_type", ALL_STAR_TYPES)
def test_star_generates_without_error(star_type):
    for _ in range(TRIALS):
        make_star(star_type)


@pytest.mark.parametrize("star_type", ALL_STAR_TYPES)
def test_age_never_exceeds_lifespan(star_type):
    for _ in range(TRIALS):
        s = make_star(star_type)
        if s.lifespan == float('inf'):
            continue
        assert s.age <= s.lifespan * 1.001, (
            f"{star_type}: age {s.age} Gy exceeds lifespan {s.lifespan} Gy"
        )


@pytest.mark.parametrize("star_type", [f"M{sub}V" for sub in SUBCLASSES] + [f"K{sub}V" for sub in SUBCLASSES])
def test_long_lived_main_sequence_age_never_exceeds_universe_age(star_type):
    """
    K and M main-sequence dwarfs can have theoretical lifespans of tens to
    thousands of billions of years (STAR_EVOLUTION["M"]["lifespan_gy"] runs
    up to 5500 Gy) -- real astrophysics, since no red dwarf has ever
    actually died of old age. Before UNIVERSE_AGE_GY was introduced, age was
    drawn as up to 90% of that lifespan with no cosmological ceiling,
    producing stars "918 billion years old" -- older than the ~13.8 Gy
    universe itself, even though age <= lifespan still held. This is that
    regression test.
    """
    for _ in range(TRIALS):
        s = make_star(star_type)
        assert s.age <= prog_c.UNIVERSE_AGE_GY * 1.001, (
            f"{star_type}: age {s.age} Gy exceeds the age of the universe "
            f"({prog_c.UNIVERSE_AGE_GY} Gy)"
        )

        cfg = SystemConfig()
        cfg.STAR_TYPE = star_type
        cfg.AGE = "old"
        s_old = Star(cfg)
        assert s_old.age <= prog_c.UNIVERSE_AGE_GY * 1.001, (
            f"{star_type} (AGE='old'): age {s_old.age} Gy exceeds the age of "
            f"the universe ({prog_c.UNIVERSE_AGE_GY} Gy)"
        )


@pytest.mark.parametrize("star_type", ALL_STAR_TYPES)
def test_age_and_lifespan_are_positive_and_finite_or_inf(star_type):
    for _ in range(TRIALS):
        s = make_star(star_type)
        assert s.age > 0
        assert s.lifespan > 0
        assert s.lifespan == float('inf') or math.isfinite(s.lifespan)
        assert math.isfinite(s.age)


@pytest.mark.parametrize("yerkes_class", ["VII", "D"])
def test_white_dwarfs_have_infinite_lifespan(yerkes_class):
    for spec in SPECTRAL_CLASSES:
        s = make_star(f"{spec}5{yerkes_class}")
        assert s.lifespan == float('inf'), (
            f"{spec}5{yerkes_class}: expected infinite lifespan for a white dwarf, got {s.lifespan}"
        )
        assert prog_c.WHITE_DWARF_MIN_AGE_GY <= s.age <= prog_c.WHITE_DWARF_MAX_AGE_GY


def test_subdwarf_age_never_exceeds_universe_age():
    """
    Regression test for "subdwarf (Yerkes VI) age modeling" (see
    CHANGELOG.md and docs/TODO.md's now-resolved "Future ideas" entry):
    before Yerkes VI got its own age-generation path, it was routed through
    the generic evolved-star branch, which derives age/lifespan from the
    star's own generated mass -- but class VI's *entire* allowed mass range
    (0.1-0.8 Msun) implies a main-sequence lifespan longer than the age of
    the universe on its own, and generate_star deliberately samples a
    subdwarf's mass uniformly across that whole range (unlike every other
    evolved class, see _sample_evolved_star_mass_sol's docstring for why).
    That combination meant a generated subdwarf's age came out equal to its
    own (pre-Big-Bang) implied main-sequence lifespan, routinely hundreds of
    billions of years old. Every subdwarf's age must now stay within the
    age of the universe regardless.
    """
    for spec in SPECTRAL_CLASSES:
        for _ in range(TRIALS):
            s = make_star(f"{spec}5VI")
            assert s.age <= prog_c.UNIVERSE_AGE_GY * 1.001, (
                f"{spec}5VI: age {s.age} Gy exceeds the age of the universe "
                f"({prog_c.UNIVERSE_AGE_GY} Gy)"
            )
            assert prog_c.SUBDWARF_MIN_AGE_GY * 0.999 <= s.age <= prog_c.SUBDWARF_MAX_AGE_GY * 1.001
            assert math.isfinite(s.lifespan)
            assert s.age <= s.lifespan * 1.001


@pytest.mark.parametrize("yerkes_class", ["II", "III", "IV", "IB", "IAB", "IA", "IA+", "0"])
def test_evolved_stars_are_at_least_as_old_as_their_own_main_sequence_lifespan(yerkes_class):
    """
    A giant/supergiant/subgiant/bright-giant/hypergiant must already have
    finished its main-sequence phase to be observed as that class -- its age
    can never be less than the main-sequence lifespan implied by its own
    generated mass (see Star._calculate_initial_star_age_and_lifespan).
    """
    for spec in SPECTRAL_CLASSES:
        s = make_star(f"{spec}5{yerkes_class}")
        mass_sol = s.mass / pc.SOLAR_MASS_TO_KG
        ms_lifespan = 10.0 * mass_sol ** -2.5
        assert s.age >= ms_lifespan * 0.999, (
            f"{spec}5{yerkes_class}: age {s.age} Gy is younger than its own "
            f"implied main-sequence lifespan {ms_lifespan} Gy"
        )


@pytest.mark.parametrize("yerkes_class", ["II", "III", "IV", "IB", "IAB", "IA", "IA+", "0"])
def test_evolved_star_mass_never_implies_a_pre_big_bang_star(yerkes_class):
    """
    Regression test for "Evolved-star mass sampling can imply a
    pre-Big-Bang star" (see TODO.md's "Future ideas" / CHANGELOG.md): an
    evolved-class star's own generated mass must never imply a
    main-sequence lifespan longer than the age of the universe -- such a
    progenitor couldn't have finished its main-sequence phase yet, so a
    star with that mass could never actually be observed as an evolved
    class today. Class III (Giant)'s allowed mass range (0.8-8 Msun, see
    physical_constants.YERKES_MASS_CONSTRAINTS) is the one that actually
    brushes the ~0.88 Msun cutoff, so more trials run there to exercise the
    reject-and-resample loop (Star._sample_evolved_star_mass_sol).
    """
    trials = 200 if yerkes_class == "III" else 20
    for spec in SPECTRAL_CLASSES:
        for _ in range(trials):
            s = make_star(f"{spec}5{yerkes_class}")
            mass_sol = s.mass / pc.SOLAR_MASS_TO_KG
            ms_lifespan = 10.0 * mass_sol ** -2.5
            assert ms_lifespan <= prog_c.UNIVERSE_AGE_GY * 1.001, (
                f"{spec}5{yerkes_class}: mass {mass_sol} Msun implies a "
                f"main-sequence lifespan of {ms_lifespan} Gy, longer than the "
                f"age of the universe ({prog_c.UNIVERSE_AGE_GY} Gy) -- this "
                f"progenitor couldn't have finished its main-sequence phase yet."
            )


def test_sample_evolved_star_mass_sol_rejects_pre_big_bang_masses():
    """
    Unit-level regression test for Star._sample_evolved_star_mass_sol
    itself: every mass it returns for class III (Giant)'s allowed range
    (0.8-8 Msun -- the one class where the ~0.88 Msun cutoff actually falls
    inside the range) must have an implied main-sequence lifespan within
    the age of the universe, even though plain uniform sampling over that
    same range would sometimes land below the cutoff.
    """
    min_mass_sol, max_mass_sol = pc.YERKES_MASS_CONSTRAINTS["III"]
    for _ in range(500):
        mass_sol = _sample_evolved_star_mass_sol(min_mass_sol, max_mass_sol)
        assert min_mass_sol <= mass_sol <= max_mass_sol
        ms_lifespan = prog_c.SOLAR_MS_LIFESPAN_GY * mass_sol ** prog_c.MS_LIFESPAN_MASS_EXPONENT
        assert ms_lifespan <= prog_c.UNIVERSE_AGE_GY * 1.001, (
            f"sampled mass {mass_sol} Msun implies a main-sequence lifespan of "
            f"{ms_lifespan} Gy, longer than the age of the universe"
        )


def test_sample_evolved_star_mass_sol_raises_when_range_is_entirely_invalid():
    """
    A mass range whose *every* value implies a main-sequence lifespan
    longer than the age of the universe (e.g. class VI's own 0.1-0.8 Msun
    range, which is exactly why Star.generate_star excludes VI from this
    check) must fail loudly rather than looping forever or silently
    returning a self-contradictory mass.
    """
    with pytest.raises(ValueError):
        _sample_evolved_star_mass_sol(*pc.YERKES_MASS_CONSTRAINTS["VI"])


@pytest.mark.parametrize("star_type", ALL_STAR_TYPES)
def test_heliosphere_and_perimeter_are_positive_and_finite(star_type):
    for _ in range(TRIALS):
        s = make_star(star_type)
        assert math.isfinite(s.heliosphere_radius)
        assert math.isfinite(s.system_perimeter)
        assert s.heliosphere_radius > 0
        assert s.system_perimeter > 0


def test_heliosphere_never_reaches_absurd_galactic_scale():
    """
    Regression guard for the hypergiant wind-constant bug: heliosphere radii
    should stay well under galactic scale (the Milky Way is ~100,000 ly
    across) for every star type, including the most extreme hypergiants.
    """
    ten_thousand_ly_in_au = 10_000 / pc.AU_TO_LY
    for spec in SPECTRAL_CLASSES:
        for _ in range(TRIALS):
            s = make_star(f"{spec}5" + "0")  # Yerkes hypergiant class
            assert s.heliosphere_radius < ten_thousand_ly_in_au, (
                f"{spec}50: heliosphere {s.heliosphere_radius} AU is absurdly large"
            )


def test_sun_like_star_heliosphere_matches_voyager_measurement():
    """
    A G2V star's heliosphere should land in the right ballpark of the Sun's
    actual, Voyager-measured heliopause distance (~120 AU).
    """
    values = [make_star("G2V").heliosphere_radius for _ in range(20)]
    median = sorted(values)[len(values) // 2]
    assert 40 <= median <= 400, f"G2V median heliosphere {median} AU is not solar-like"


def test_gravitational_sphere_of_influence_scales_with_mass_cube_root():
    """
    Regression guard for the Hill/galactic-tidal sphere ("sphere of
    influence") formula: system_perimeter / mass^(1/3) must be a constant
    (the formula is a pure function of mass and galactic position, both
    fixed), and for a solar-mass star it should land within the literature
    range for the Sun's own galactic tidal radius (~100,000-200,000 AU,
    i.e. the commonly-cited ~1.5-2 ly edge of the Oort cloud).
    """
    ratios = []
    for star_type in ["G2V", "M5V", "K3III", "B2V", "A0VII"]:
        s = make_star(star_type)
        mass_sol = s.mass / pc.SOLAR_MASS_TO_KG
        ratios.append(s.system_perimeter / mass_sol ** (1 / 3))

    for r in ratios[1:]:
        assert r == pytest.approx(ratios[0], rel=1e-6)

    assert 80_000 <= ratios[0] <= 250_000


@pytest.mark.parametrize("star_type", ALL_STAR_TYPES)
def test_galactic_orbit_is_positive_and_finite(star_type):
    for _ in range(TRIALS):
        s = make_star(star_type)
        assert math.isfinite(s.galactic_orbital_speed_kms)
        assert math.isfinite(s.galactic_orbital_period_gy)
        assert s.galactic_orbital_speed_kms > 0
        assert s.galactic_orbital_period_gy > 0


def test_galactic_orbit_is_independent_of_stellar_mass():
    """
    Unlike `system_perimeter` (a Hill sphere, scaling with the star's own
    mass), circular orbital speed/period around the galactic center depends
    only on distance from the galactic center -- every star type generated
    without an explicit `galactic_center_dist_ly` falls back to the same
    fixed `physical_constants.GALACTIC_CENTER_DISTANCE_LY`, so its
    galactic-orbit values should be identical (not just proportional)
    across wildly different masses.
    """
    speeds = []
    periods = []
    for star_type in ["G2V", "M5V", "K3III", "B2V", "A0VII"]:
        s = make_star(star_type)
        speeds.append(s.galactic_orbital_speed_kms)
        periods.append(s.galactic_orbital_period_gy)

    for speed, period in zip(speeds[1:], periods[1:]):
        assert speed == pytest.approx(speeds[0], rel=1e-9)
        assert period == pytest.approx(periods[0], rel=1e-9)


def test_sun_like_star_galactic_orbit_matches_real_measurements():
    """
    Falling back to the fixed `GALACTIC_CENTER_DISTANCE_LY` (Sol's own
    galactocentric distance), a generated star's circular orbital
    speed/period should land near the Sun's own measured values: ~220-240
    km/s circular velocity and a ~225-250 million-year "galactic year".
    """
    s = make_star("G2V")
    assert 150 <= s.galactic_orbital_speed_kms <= 260
    assert 0.15 <= s.galactic_orbital_period_gy <= 0.30


@pytest.mark.parametrize("star_type", ALL_STAR_TYPES)
def test_galactic_position_change_interval_is_positive_and_finite(star_type):
    """
    No generated star falls back to the `float('inf')` degenerate case
    (that only happens at galactic_center_dist_ly <= 0, structurally
    unreachable through this fallback path -- see
    `utils.calculate_galactic_orbit`'s own `distance_ly <= 0` guard).
    """
    for _ in range(TRIALS):
        s = make_star(star_type)
        assert math.isfinite(s.galactic_position_change_interval_hours)
        assert s.galactic_position_change_interval_hours > 0


def test_galactic_position_change_interval_matches_manual_formula():
    from stellarObjects.utils import position_change_interval_hours

    s = make_star("G2V")
    expected = position_change_interval_hours(s.radius, s.galactic_orbital_speed_kms)
    assert s.galactic_position_change_interval_hours == pytest.approx(expected, rel=1e-9)


def test_sun_like_star_galactic_position_change_interval_is_a_couple_hours():
    """
    Sanity check against the hand-computed order of magnitude: a Sol-like
    star's own diameter (~1.39 million km) divided by its ~206 km/s
    galactic orbital speed comes out to roughly 1.5-2.5 hours.
    """
    s = make_star("G2V")
    assert 0.5 <= s.galactic_position_change_interval_hours <= 5.0


_MAIN_SEQUENCE_ONLY_NOTE_FRAGMENTS = [
    note["evolutionary_constraint_notes"]
    for note in prog_c.STAR_EVOLUTION.values()
    if "evolutionary_constraint_notes" in note
]


@pytest.mark.parametrize("yerkes_class", ["VII", "III", "IB", "0"])
def test_evolved_star_narrative_never_quotes_main_sequence_notes(yerkes_class):
    """
    Regression guard for the Star.to_paragraph_list bug where a giant,
    supergiant, or white dwarf's age sentence quoted a main-sequence-only
    evolutionary note belonging to its current-temperature letter (e.g. a
    white dwarf's age exceeding the "main-sequence lifespan" quoted for it).
    """
    for spec in SPECTRAL_CLASSES:
        s = make_star(f"{spec}5{yerkes_class}")
        age_sentence = s.to_paragraph_list()[1]
        for fragment in _MAIN_SEQUENCE_ONLY_NOTE_FRAGMENTS:
            assert fragment not in age_sentence, (
                f"{spec}5{yerkes_class}: narrative wrongly includes a main-sequence-only note"
            )
