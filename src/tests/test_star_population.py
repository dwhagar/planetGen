"""
The star population model (`planetgen.physics.stellar_evolution`, used by
`Star.generate_star` for every random star).

The shares are checked against the real census the model is built to
reproduce (see docs/design and the "Star population model" section of
program_constants): most stars are M dwarfs, about 5% are white dwarfs,
giants are a fraction of a percent and supergiants vanishingly rare.

Run with: pytest src/tests/test_star_population.py
"""
import collections
import math
import random

import pytest

from planetgen.physics import constants as phys_c, imf, main_sequence
from planetgen import tuning as prog_c
from planetgen.physics import alpha_tracks, stellar_abundances, stellar_evolution as se
from planetgen.generation import star_population as sp
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import STAR_TYPE_PATTERN, Star
from planetgen.generation.system import StarSystem
from planetgen.util import draw

SAMPLE_SIZE = 200_000


def _category(state, rng_letter):
    yerkes = state["yerkes_class"]
    if yerkes == "V":
        return f"{rng_letter} dwarf"
    if yerkes == "IV":
        return "subgiant"
    if yerkes in ("III", "II"):
        return "giant"
    if yerkes == "VII":
        return "white dwarf"
    if yerkes in ("PMS", "protostar"):
        return "pre-main-sequence"
    return "supergiant"


@pytest.fixture(scope="module")
def population():
    rng = random.Random(20260930)
    counts = collections.Counter()
    states = []
    for _ in range(SAMPLE_SIZE):
        mass, age, state = se.sample_living_star(rng=rng)
        letter, _ = se.spectral_letter_and_subclass(state["temperature_k"])
        counts[_category(state, letter)] += 1
        states.append((mass, age, state))
    return counts, states


# (category, low %, high %): the model's expected share with room for
# sampling noise and later tuning of the letter boundaries.
EXPECTED_SHARES = [
    ("M dwarf", 68.0, 80.0),
    # GEN.897 and GEN.943: the K/M edge is 3,890 to 3,900 K (M0V), so K is the observed 12% of the main-sequence stars, 8.6% of this sample,
    # which also holds the pre-main-sequence stars, white dwarfs and giants (K was 11 to 18% with the edge at 3700 K, where it ran to 19%)
    ("K dwarf", 7.0, 12.0),
    ("G dwarf", 2.0, 6.0),
    ("F dwarf", 1.2, 3.0),
    # GEN.652: on the MIST main sequence a star is A-type (7500 to 10000 K) between about 1.75 and 2.4 solar masses and
    # cools out of it before its main sequence ends, narrower than the piecewise law's 1.4 to 2.4, so A and B are rarer.
    ("A dwarf", 0.1, 0.4),
    ("B dwarf", 0.03, 0.15),
    ("white dwarf", 4.0, 7.0),
    ("subgiant", 0.08, 0.3),  # GEN.502: a light star's subgiant phase is 5% of its main sequence, not 10%
    ("giant", 0.12, 0.5),
    # GEN.1174: a star is pre-main-sequence until its radius has settled (BHAC15, [R0254]): 0.19 Gy at 0.2 solar masses, 0.40 at 0.1, 0.05 at 0.8, no longer
    # the 0.1 to 0.7 Gy of MIST's nuclear arrival (GEN.529, kept as `nuclear_arrival_yr`). The share is the part of the stars younger than that, 1.8 to 1.9% of this sample.
    ("pre-main-sequence", 1.2, 6.0),
]


@pytest.mark.parametrize("category, low, high", EXPECTED_SHARES)
def test_star_shares_match_the_census(population, category, low, high):
    counts, _ = population
    share = 100 * counts[category] / SAMPLE_SIZE
    assert low <= share <= high, f"{category}: {share:.3f}% (want {low}-{high}%)"


def test_supergiants_are_vanishingly_rare(population):
    counts, _ = population
    # Real share ~1e-6; allow a few in 200,000 for noise.
    assert counts["supergiant"] <= 5


def test_main_sequence_luminosity_follows_mass(population):
    _, states = population
    for mass, age, state in states:
        if state["yerkes_class"] == "V":
            row = main_sequence.state(mass, age, alpha_tracks.track_fe_h(state["metallicity_fe_h"], age))
            assert state["luminosity_sol"] == pytest.approx(row.luminosity_sol)


def test_every_state_is_in_its_own_phase(population):
    _, states = population
    for mass, age, state in states:
        t_ms = se.main_sequence_lifetime_gy(mass, alpha_tracks.track_fe_h(state["metallicity_fe_h"], age))
        if state["yerkes_class"] in ("V", "PMS", "protostar"):
            assert age < t_ms
        else:
            assert age >= t_ms
        assert age < state["phase_end_gy"]


def test_white_dwarfs_are_earth_sized_and_dim(population):
    _, states = population
    white_dwarfs = [state for _, _, state in states if state["yerkes_class"] == "VII"]
    assert white_dwarfs
    low, high = prog_c.WD_MASS_RANGE_SOL
    # The heaviest the initial-final mass relation makes (1.33 solar masses, [R0961]) is smaller than the
    # zero-temperature radius at the top of the range (2,700 km); a cold star is never smaller than that.
    smallest = stellar_abundances.white_dwarf_radius_solar(high) * phys_c.SOLAR_RADIUS_M / 1000.0
    for state in white_dwarfs:
        assert low <= state["mass_sol"] <= high
        assert smallest < state["radius_km"] < 20000          # a white dwarf just formed is inflated: 15,700 km at 93,000 K and 0.62 solar masses (GEN.742)
        assert state["luminosity_sol"] <= prog_c.WD_LUMINOSITY_RANGE_SOL[1]


def test_imf_sampler_respects_truncation():
    rng = random.Random(1)
    masses = [imf.sample(min_mass_sol=1.4, rng=rng) for _ in range(5000)]
    assert min(masses) >= 1.4
    assert max(masses) <= prog_c.IMF_BREAKS_SOL[-1]
    # Kroupa: about half of all stars are below 0.3 Msun.
    all_masses = [imf.sample(rng=rng) for _ in range(20000)]
    below = sum(m < 0.3 for m in all_masses) / len(all_masses)
    assert 0.4 < below < 0.65


def test_collapsed_massive_stars_return_none():
    assert se.evolve_star(20.0, 5.0) is None


@pytest.mark.parametrize("letter_temperature", [(35000, "O"), (5772, "G"), (3000, "M"), (100000, "O"), (1000, "M")])
def test_spectral_letter_from_temperature(letter_temperature):
    temperature, letter = letter_temperature
    assert se.spectral_letter_and_subclass(temperature)[0] == letter


def test_random_stars_are_self_consistent():
    for _ in range(300):
        star = Star(SystemConfig())
        assert (STAR_TYPE_PATTERN.fullmatch(star.type.split()[0])
                or star.yerkes_class in ("PMS", "protostar"))  # GEN.242: a contracting star has no numeral
        assert star.initial_mass_sol is not None
        assert 0 <= star.age < star.phase_end_age_gy
        if star.yerkes_class == "V":
            mass_sol = star.mass / phys_c.SOLAR_MASS_TO_KG
            luminosity_sol = star.luminosity / phys_c.SOLAR_LUMINOSITY
            row = main_sequence.state(star.initial_mass_sol, star.age, alpha_tracks.track_fe_h(star.metallicity_fe_h, star.age))
            assert luminosity_sol == pytest.approx(row.luminosity_sol, rel=1e-6)


def test_large_stars_are_massive_and_alive():
    cfg = SystemConfig()
    cfg.LARGE_STAR = True
    for _ in range(200):
        star = Star(cfg)
        assert star.initial_mass_sol >= prog_c.LARGE_STAR_MIN_MASS_SOL
        assert star.yerkes_class != "VII"


def test_binary_companions_share_the_primarys_age_and_are_dimmer():
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = True
    cfg.PLANETS = False
    for _ in range(150):
        system = StarSystem(system_config=cfg)
        primary, secondary = system.primary_star, system.secondary_star
        assert secondary.mass <= primary.mass
        if primary.yerkes_class == "V" and secondary.yerkes_class == "V":
            assert secondary.luminosity <= primary.luminosity * (1 + 1e-9)
        assert math.isclose(primary.age, secondary.age, rel_tol=1e-9)


def test_planets_respect_their_stars_age_and_history():
    for _ in range(150):      # a whole system with moons, the air and the heat costs about 4 s, so 400 of them passed 900 s
        system = StarSystem(system_config=SystemConfig())
        hosts = [(system.star, system.planets)]
        if system.binary_type == "wide":
            hosts.append((system.secondary_star, system.secondary_planets))
        for star, bodies in hosts:
            engulfed_au = system._engulfment_radius_au(star)
            for body in bodies:
                inner_edge = body.lower_limit if body.body_type == "a" else body.distance
                assert inner_edge >= engulfed_au
                if star.age < prog_c.PLANET_MIN_STAR_AGE_GY:
                    assert body.body_type == "a"
                if star.age < prog_c.LIFE_MIN_STAR_AGE_GY and body.body_type != "a":
                    # the class letter is derived from the conditions and says nothing of life (GEN.915); a star this young has none
                    assert body.life_class is None
                    assert all(moon.life_class is None for moon in body.moons)


def test_a_requested_habitable_world_makes_a_star_too_young_or_dead_for_life_rarer_not_impossible():
    """A request only raises the odds (GEN.1130): a star younger than life needs, or a white dwarf, is kept one time in `REQUEST_ODDS`, so the share
    of them falls to a fraction of what the plain draw gives, and is not forced to zero."""
    def unfit_share(requested):
        draw.set_run_seed(5)
        cfg = SystemConfig()
        cfg.HABITABLE_WORLD = requested
        stars = [Star(cfg) for _ in range(600)]
        return sum(star.age < prog_c.LIFE_MIN_STAR_AGE_GY or star.yerkes_class == "VII" for star in stars) / len(stars)

    plain, requested = unfit_share(False), unfit_share(True)
    assert plain > 0.0
    assert requested <= 0.5 * plain


def test_from_params_rebuilds_a_star_without_rolling():
    rng = random.Random(7)
    mass, age, state = se.sample_living_star(rng=rng)
    params = se.star_params(mass, age, state)
    star = Star.from_params(params, SystemConfig(), name="Beacon")
    assert star.name == "Beacon"
    assert (star.type, star.mass, star.radius, star.temperature, star.luminosity, star.age) == (
        params["type"], params["mass_kg"], params["radius_km"], params["temperature_k"],
        params["luminosity_w"], params["age_gy"])
    assert star.habitable_zone and star.system_perimeter > 0

    system = StarSystem(system_config=SystemConfig(), primary_star_params=params)
    assert system.primary_star.mass == params["mass_kg"] or getattr(system, "secondary_star", system.primary_star).mass == params["mass_kg"]


@pytest.mark.parametrize("population, floor", [("old", 1000.0), ("halo", 50.0)])
def test_a_light_stars_cell_bound_holds_and_is_tight_enough_to_draw(population, floor):
    """GEN.502: the phase fractions of a light star move with its mass, so a cell's bound takes them over its own
    masses; a bound over every light star accepted one draw in a hundred."""
    table = sp._bright_table(floor, population)
    window = table["window"]
    cumulative = table["cumulative"]
    rng = random.Random(3)
    accepted = 0.0
    previous = 0.0
    for index, (low, high, alpha, bound) in enumerate(table["cells"]):
        share = (cumulative[index] - previous) / cumulative[-1]
        previous = cumulative[index]
        measures = []
        for _ in range(6):
            mass = imf.sample_power_law(low, high, alpha, rng)
            measures.append(sum(max(0.0, min(end, window[1]) - max(start, window[0])) * chance
                                for start, end, chance in sp._bright_windows(mass, floor)))
        assert max(measures) <= bound * (1 + 1e-9), (low, high)
        accepted += share * sum(measures) / len(measures) / bound
    assert accepted > 0.02, accepted


@pytest.mark.parametrize("min_luminosity_sol", [50.0, 2500.0, 5.0e4])
@pytest.mark.parametrize("population", ["old", "young"])
def test_a_bright_table_bound_holds_for_every_star_of_its_fe_h_bin(population, min_luminosity_sol):
    """GEN.651: a star draws from the table of its own [Fe/H] bin, so a cell's bound has to be at least the measure of
    every mass in the cell at every composition in the bin; this is what `BRIGHT_TABLE_BOUND_MARGIN` is calibrated by."""
    rng = random.Random(651)
    step = prog_c.BRIGHT_TABLE_FE_H_STEP_DEX
    for fe_bin in sp.FE_H_BINS:
        table = sp._bright_table(min_luminosity_sol, population, None, sp._law(population, None), fe_bin)
        for a, b, alpha, bound in table["cells"]:
            for _ in range(8):
                mass = a + (b - a) * rng.random()
                fe_h = fe_bin + (rng.random() - 0.5) * step
                measure = sp._bright_measure(mass, min_luminosity_sol, table["law"], fe_h)
                assert measure <= bound, f"bin {fe_bin}, mass {mass}, [Fe/H] {fe_h}: {measure} over bound {bound}"


@pytest.mark.parametrize("fe_h", [-0.2, 0.0, 0.3])
def test_a_bright_star_is_drawn_at_any_composition_that_can_make_one(fe_h):
    rng = random.Random(7)
    kept = 0
    for _ in range(5):
        params = sp.sample_bright_stars(1, 2500.0, "old", rng, fe_h_source=lambda index: fe_h)[0]
        assert params["luminosity_w"] >= 2500.0 * phys_c.SOLAR_LUMINOSITY
        kept += params["metallicity_fe_h"] == fe_h
    assert kept == 5


def test_a_metal_rich_place_has_no_giant_that_bright_in_its_window_and_its_metal_poor_tail_makes_them_all():
    """GEN.651: a place's bright stars are its prior times each composition's bright share, so a metal-rich place with a
    metal-poor tail gets its 6000 Lsun giants from the tail (the MIST thermal pulses of a metal-poor star reach higher than a
    metal-rich one's: none of the latter in an old window does)."""
    rng = random.Random(8)
    fe_hs = [0.1] * 3 + [-1.4]
    seen = iter(range(10 ** 6))
    drawn = [sp.sample_bright_stars(1, 6000.0, "old", rng, fe_h_source=lambda index: fe_hs[next(seen) % 4])[0]
             for _ in range(5)]
    assert all(params["metallicity_fe_h"] == -1.4 for params in drawn)


def test_a_populations_counts_are_worked_out_at_its_own_composition_not_the_suns():
    """GEN.661: a metal-poor population's old stars keep more of their giants bright, so its share of bright stars is not
    the solar composition's."""
    assert sp.reference_fe_h("halo") < -1.0 and sp.reference_fe_h("young") == pytest.approx(0.0, abs=0.1)
    assert sp.reference_fe_h("cluster") < -0.5 and sp.reference_fe_h(None) == sp.reference_fe_h("intermediate")
    # GEN.935: on the MIST tracks the halo's own-composition gain in bright stars is 1.09 (it was over 1.2 on the hand-set giants)
    for population, least in (("old", 1.02), ("halo", 1.05)):
        solar = sp._bright_fraction(500.0, population, None, sp._law(population, None), 0.0)       # the share at the Sun's composition
        own = sp.bright_star_fraction(500.0, population)
        assert own > least * solar, (population, own, solar)
    assert sp.massive_star_fraction(8.0, "young") > 0.0


def test_a_bright_share_is_worked_out_without_a_sampling_tables_bounds(monkeypatch):
    """The share of stars that are bright needs the measure of a few stars a cell and no bound: the version-key battery
    asks it of 23 populations, and building a table's bounds for each took minutes."""
    def refuse(*_args, **_kwargs):
        raise AssertionError("a bright share read a measure bound")
    monkeypatch.setattr(sp, "_bright_measure_bound", refuse)
    sp._bright_fraction.cache_clear()
    share = sp.bright_star_fraction(2500.0, "old")
    assert 0.0 < share < 1.0
    assert sp.bright_star_fraction(2500.0, "old", mass_range=(8.0, None)) <= share
