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
from stellarObjects.starData import Star
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
