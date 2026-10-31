# tests/test_wolf_rayet.py

"""
GEN.243: the Wolf-Rayet state (`physics/wolf_rayet.py`): the metallicity threshold, the phase at the end of a very
massive star's supergiant window, the subtypes, `evolve_star` and `star_params` for such a star, the bright pass that
draws them, the Milky Way count and the companions of the stars.
"""

import math
import random
import re

import pytest

from planetgen import tuning
from planetgen.generation import star_population
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.datasets import isochrones
from planetgen.physics import alpha_tracks, imf, post_main_sequence
from planetgen.physics import stellar_evolution as evolution
from planetgen.physics import wolf_rayet
from planetgen.util import draw
from tests import render_support


def _window(mass_sol):
    t_ms = evolution.main_sequence_lifetime_gy(mass_sol)
    return t_ms, post_main_sequence.timeline(mass_sol, 0.0).end_gy


def test_the_threshold_rises_as_metallicity_falls():
    """25 solar masses at the Sun's metallicity, about 35 at half and 56 at a fifth of it (the note's figures)."""
    assert wolf_rayet.threshold_mass_sol() == pytest.approx(25.0)
    assert wolf_rayet.threshold_mass_sol(math.log10(0.5)) == pytest.approx(35.0, abs=1.0)
    assert wolf_rayet.threshold_mass_sol(math.log10(0.2)) == pytest.approx(56.0, abs=1.0)
    assert wolf_rayet.threshold_mass_sol(0.3) < 25.0


def test_the_phase_is_the_last_third_of_a_myr_and_never_more_than_half_the_supergiant_window():
    assert wolf_rayet.phase_gy(2.0e-3) == pytest.approx(tuning.WR_PHASE_MYR / 1000.0)
    assert wolf_rayet.phase_gy(0.4e-3) == pytest.approx(0.2e-3)
    for mass in (25.0, 40.0, 60.0, 100.0, 150.0):
        t_ms, t_end = _window(mass)
        start = wolf_rayet.start_gy(mass, t_ms, t_end)
        assert t_ms < start < t_end
        assert t_end - start <= tuning.WR_PHASE_MYR / 1000.0 + 1.0e-12
    t_ms, t_end = _window(20.0)
    assert wolf_rayet.start_gy(20.0, t_ms, t_end) is None


def test_a_star_is_wolf_rayet_only_in_the_last_stretch_and_only_above_the_threshold():
    rng = random.Random(1)
    t_ms, t_end = _window(60.0)
    start = wolf_rayet.start_gy(60.0, t_ms, t_end)
    assert evolution.evolve_star(60.0, 0.5 * (t_ms + start), rng, fe_h=0.0)["yerkes_class"] in ("IA", "0", "IAB", "IB")
    state = evolution.evolve_star(60.0, 0.5 * (start + t_end), rng, fe_h=0.0)
    assert state["yerkes_class"] == wolf_rayet.CLASS and state["wr_subtype"] in ("WN", "WC", "WO")
    assert evolution.evolve_star(60.0, t_end + 1.0e-6, rng, fe_h=0.0) is None  # collapsed
    low_t_ms, low_end = _window(30.0)
    # 30 solar masses is a Wolf-Rayet star at the Sun's metallicity and keeps its supergiant at a fifth of it.
    start_30 = wolf_rayet.start_gy(30.0, low_t_ms, low_end)
    assert evolution.evolve_star(30.0, low_end - 1.0e-5, rng, fe_h=0.0)["yerkes_class"] == wolf_rayet.CLASS
    assert start_30 < low_end - 1.0e-5
    assert evolution.evolve_star(30.0, low_end - 1.0e-5, rng, fe_h=math.log10(0.2))["yerkes_class"] != wolf_rayet.CLASS


def test_the_state_follows_the_ranges_of_the_note():
    rng = random.Random(2)
    t_ms, t_end = _window(60.0)
    start = wolf_rayet.start_gy(60.0, t_ms, t_end)
    for _ in range(300):
        age = rng.uniform(start, t_end)
        state = evolution.evolve_star(60.0, age, rng, fe_h=0.0)
        low, high = tuning.WR_LUMINOSITY_RANGE_SOL
        assert low <= state["luminosity_sol"] <= high
        low_k = min(r[0] for r in tuning.WR_TEMPERATURE_RANGE_K.values())
        high_k = max(r[1] for r in tuning.WR_TEMPERATURE_RANGE_K.values())
        assert low_k <= state["temperature_k"] <= high_k
        # GEN.939/940: the present mass is the track's, and the luminosity follows it by the helium-star relation.
        assert state["mass_sol"] == pytest.approx(post_main_sequence.timeline(60.0, 0.0).state(age).mass_sol)
        intercept, slope, sigma = isochrones.helium_star_relation()
        assert abs(math.log10(state["luminosity_sol"]) - (intercept + slope * math.log10(state["mass_sol"]))) <= 6 * sigma
        assert state["phase_end_gy"] == t_end


def test_the_subtype_follows_the_phase_and_wo_needs_a_heavy_star():
    assert [wolf_rayet.subtype(60.0, p) for p in (0.0, 0.59, 0.61, 0.89, 0.91, 0.99)] == [
        "WN", "WN", "WC", "WC", "WO", "WO"]
    assert wolf_rayet.subtype(30.0, 0.95) == "WC"


def test_the_text_names_the_subtype_and_a_subclass_in_its_range():
    rng = random.Random(3)
    t_ms, t_end = _window(100.0)
    start = wolf_rayet.start_gy(100.0, t_ms, t_end)
    seen = set()
    for _ in range(200):
        age = rng.uniform(start, t_end)
        state = evolution.evolve_star(100.0, age, rng, fe_h=0.0)
        params = evolution.star_params(100.0, age, state)
        match = re.fullmatch(r"(W[NCO])([0-9]) [A-Za-z-]+ Wolf-Rayet Star", params["type"])
        assert match, params["type"]
        low, high = sorted(tuning.WR_SUBCLASS_RANGE[match.group(1)])
        assert low <= int(match.group(2)) <= high
        assert params["yerkes_class"] == "W"
        seen.add(match.group(1))
    assert seen == {"WN", "WC", "WO"}


def test_a_wolf_rayet_star_builds_a_system_and_rebuilds_from_its_record():
    draw.set_run_seed(5)
    t_ms, t_end = _window(60.0)
    age = 0.5 * (wolf_rayet.start_gy(60.0, t_ms, t_end) + t_end)
    params = evolution.star_params(60.0, age, evolution.evolve_star(60.0, age, fe_h=0.0))
    config = SystemConfig()
    star = Star.from_params(params, config, name="Wr")
    assert star.yerkes_class == "W" and star.type.startswith("W")
    assert Star.from_dict(star.to_dict(), config).type == star.type
    assert "Wolf-Rayet" in star.type
    config.PLANETS = True
    for seed in range(20):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=SystemConfig(), primary_star_params=params)
        assert not [p for p in system.planets if p.body_type != "a"]  # too young for planets


@pytest.mark.parametrize("subtype, age_fraction, low, high", [("WN", 0.2, 0.44, 0.60), ("WC", 0.75, 0.82, 0.96)])
def test_the_companions_are_those_of_the_dsilva_survey(subtype, age_fraction, low, high):
    """52% of WN and 90% of WC stars have a companion (`tuning.WR_COMPANION_CHANCE`)."""
    assert wolf_rayet.companion_chance(f"{subtype}6 Blue Wolf-Rayet Star") == tuning.WR_COMPANION_CHANCE[subtype]
    assert wolf_rayet.companion_chance("G2V Yellow Main Sequence Star") is None
    t_ms, t_end = _window(60.0)
    start = wolf_rayet.start_gy(60.0, t_ms, t_end)
    age = start + age_fraction * (t_end - start)
    found = trials = 0
    for seed in range(400):
        draw.set_run_seed(seed)
        params = evolution.star_params(60.0, age, evolution.evolve_star(60.0, age, fe_h=0.0))
        if not params["type"].startswith(subtype) and subtype == "WN":
            continue
        config = SystemConfig()
        config.PLANETS = False
        system = StarSystem(system_config=config, primary_star_params=params)
        trials += 1
        found += system.binary_type is not None
    assert trials > 100 and low < found / trials < high, (found, trials)


def test_the_bright_pass_draws_wolf_rayet_stars_from_young_populations_and_none_from_old_ones(monkeypatch):
    draw.set_run_seed(3)
    stars = star_population.sample_bright_stars(2000, 1.0e5, population="young")
    share = sum(star["yerkes_class"] == "W" for star in stars) / len(stars)
    assert 0.01 < share < 0.15
    with pytest.raises(ValueError):
        star_population.sample_bright_stars(5, 1.0e5, population="old")
    # A floor above their luminosities leaves none. Few tries are kept there (about 1 in 500 once the metal-poor stars read their
    # alpha-enhanced tracks, GEN.1024), so the cap on the tries a star may take is raised for this one draw.
    monkeypatch.setattr(tuning, "STAR_MODEL_MAX_REDRAWS", 20000)
    stars = star_population.sample_bright_stars(300, tuning.WR_LUMINOSITY_RANGE_SOL[1] * 1.01, population="young")
    assert all(star["yerkes_class"] != "W" for star in stars)


def _stars_born_at_or_above(threshold_sol):
    """Stars born at or above `threshold_sol` per solar mass formed, by the repository's IMF."""
    breaks, slopes = tuning.IMF_BREAKS_SOL, tuning.IMF_SLOPES
    factors = [1.0]
    for i in range(1, len(slopes)):
        factors.append(factors[-1] * breaks[i] ** (slopes[i] - slopes[i - 1]))
    number = mass = 0.0
    for i, alpha in enumerate(slopes):
        low, high = breaks[i], breaks[i + 1]
        mass += factors[i] * imf.power_law_integral(low, high, alpha - 1.0)
        low = max(low, threshold_sol)
        if low < high:
            number += factors[i] * imf.power_law_integral(low, high, alpha)
    return number / mass


def test_the_milky_way_holds_about_twelve_hundred_wolf_rayet_stars():
    """The note's check: a star formation rate of 1.65 solar masses a year with this IMF, the threshold and the phase
    give the 1,200 +- 100 known (and expected) in the Galaxy, within the 300 the test allows."""
    rate_per_year = 1.65 * _stars_born_at_or_above(wolf_rayet.threshold_mass_sol())
    count = rate_per_year * tuning.WR_PHASE_MYR * 1.0e6
    assert 900.0 < count < 1500.0, count


def test_none_is_born_old_enough_to_be_one_in_an_old_population():
    """Every star that reaches the phase is younger than the longest-lived threshold star."""
    t_ms, t_end = _window(wolf_rayet.threshold_mass_sol())
    assert t_end * 1000.0 < tuning.NEBULA_CLASS_MAX_AGE_MYR.get("X", 10.0)
    assert t_end < min(low for low, _ in tuning.STELLAR_POPULATION_AGE_RANGES_GY.values() if low >= 3.0)


def _wolf_rayet_star(seed, mass_sol=60.0):
    draw.set_run_seed(seed)
    t_ms, t_end = _window(mass_sol)
    start = wolf_rayet.start_gy(mass_sol, t_ms, t_end)
    age = start + draw.random() * (t_end - start)
    return Star.from_params(evolution.star_params(mass_sol, age, evolution.evolve_star(mass_sol, age, fe_h=0.0)),
                            SystemConfig(), name="Wr")


def test_the_ring_class_is_a_wolf_rayet_only_emission_class_of_two_to_thirty_five_parsecs():
    data = tuning.NEBULA_CLASSES["X"]
    assert data["family"] == "emission" and data["host_only"]
    low, high = data["radius_range_ly"]
    assert low == pytest.approx(2.0 * 3.2616, rel=0.03) and high == pytest.approx(35.0 * 3.2616, rel=0.03)
    assert ("W", 0, 9, ("X",), 0.4) in tuning.NEBULA_HOST_RULES
    assert tuning.NEBULA_CLASS_MAX_AGE_MYR["X"] == 10.0
    # The longest a star that reaches the phase lives (the threshold mass at the Sun's metallicity) is under the limit.
    assert _window(wolf_rayet.threshold_mass_sol())[1] * 1000.0 < tuning.NEBULA_CLASS_MAX_AGE_MYR["X"]


def test_no_nebula_drawn_without_a_host_is_a_ring():
    from planetgen.generation.phenomena import nebula as nebula_module
    from planetgen.generation import phenomena_plausibility
    assert "X" not in nebula_module.FREE_CLASS_LETTERS and "X" in nebula_module.NEBULA_CLASS_LETTERS
    assert "X" not in phenomena_plausibility._NEBULA_CLASS_SHARES
    assert sum(phenomena_plausibility._NEBULA_CLASS_SHARES.values()) == pytest.approx(1.0)
    draw.set_run_seed(1)
    for _ in range(300):
        assert nebula_module.Nebula(SystemConfig(), nebula_type="emission").nebula_class != "X"


def test_a_wolf_rayet_star_lights_a_ring_in_forty_percent_and_the_ring_names_it():
    from planetgen.generation import nebula_contents
    rings = []
    trials = 500
    for seed in range(trials):
        star = _wolf_rayet_star(seed)
        nebula = nebula_contents.hosted_nebula_of(star)
        if nebula is not None:
            rings.append((star, nebula))
    assert 0.32 < len(rings) / trials < 0.48, len(rings)
    low, high = tuning.NEBULA_CLASSES["X"]["radius_range_ly"]
    for star, nebula in rings:
        assert nebula.nebula_class == "X" and nebula.nebula_type == "emission"
        assert nebula.hosts == f"{star.type.split(' ', 1)[0]} x1"
        assert low <= nebula.radius_ly <= high and nebula.q_total is None
        assert nebula.age_myr == pytest.approx(star.age * 1000.0)
        assert "wolf-rayet ring nebula" in render_support.facts_text(render_support.phenomenon_section(nebula))


def test_a_wolf_rayet_star_older_than_the_limit_lights_nothing_and_neither_does_a_remnant_class_star():
    import types
    from planetgen.generation import nebula_contents
    old = types.SimpleNamespace(type="WN6 Blue Wolf-Rayet Star", age=0.02, mass=1.0e32, temperature=60000.0)
    assert nebula_contents.hosted_nebula_of(old) is None
    ordinary = types.SimpleNamespace(type="G2V Yellow Main Sequence Star", age=0.004, mass=1.0e30, temperature=5800.0)
    assert nebula_contents.hosted_nebula_of(ordinary) is None


def test_the_class_tables_hold_the_wolf_rayet_class():
    from planetgen.generation import star_labels
    from planetgen.generation import forced_star
    found = forced_star.class_summary("W")
    assert found.luminosity_sol == tuning.WR_LUMINOSITY_RANGE_SOL
    low, high = found.mass_sol
    for initial in (tuning.WR_MIN_INITIAL_MASS_SOL, 40.0, 60.0, 100.0, tuning.IMF_BREAKS_SOL[-1]):
        line = post_main_sequence.timeline(initial, 0.0)
        assert low <= line.state(line.end_gy - 1.0e-9).mass_sol <= high, initial
    assert star_labels.star_label("WN6", "W") == ("Wolf-Rayet star", "Wolf-Rayet stars")


@pytest.mark.parametrize("fe_h, initial", [(0.0, 30.0), (0.0, 60.0), (0.0, 150.0), (-0.7, 60.0), (-0.7, 150.0), (0.3, 40.0)])
def test_the_luminosity_of_a_wolf_rayet_star_follows_its_present_mass_by_the_mist_helium_star_relation(fe_h, initial):
    """GEN.939: a stripped star's light is set by its mass (the relation of MIST's own Wolf-Rayet rows), not drawn over the class's
    range whatever star it came from; the present mass is the track's (GEN.940)."""
    rng = random.Random(7)
    track = alpha_tracks.track_segments(fe_h)[0][2]          # GEN.1024: a young star reads the track of the thin branch
    line = post_main_sequence.timeline(initial, track)
    start = wolf_rayet.start_gy(initial, line.main_sequence_gy, line.end_gy, track)
    intercept, slope, sigma = isochrones.helium_star_relation()
    assert slope > 0.5 and 0.0 < sigma < 0.3
    for age in [start + (line.end_gy - start) * share for share in (0.0, 0.3, 0.7, 0.99)]:
        state = evolution.evolve_star(initial, age, rng, fe_h=fe_h)
        assert state["yerkes_class"] == "W"
        assert state["mass_sol"] == pytest.approx(line.state(age).mass_sol)
        assert abs(math.log10(state["luminosity_sol"]) - (intercept + slope * math.log10(state["mass_sol"]))) <= 6 * sigma


def test_a_heavier_wolf_rayet_star_is_brighter_and_the_chance_to_clear_a_floor_follows_it():
    assert wolf_rayet.bright_chance(30.0, 1.0e5) == pytest.approx(1.0)
    assert wolf_rayet.bright_chance(10.0, 1.0e6) < wolf_rayet.bright_chance(20.0, 1.0e6) < wolf_rayet.bright_chance(40.0, 1.0e6)
    assert wolf_rayet.bright_chance(40.0, tuning.WR_LUMINOSITY_RANGE_SOL[1] * 1.01) == 0.0
    rng = random.Random(1)
    light = [wolf_rayet.state(60.0, 12.0, 0.5, 0.0, 1.0, rng)["luminosity_sol"] for _ in range(200)]
    heavy = [wolf_rayet.state(60.0, 35.0, 0.5, 0.0, 1.0, rng)["luminosity_sol"] for _ in range(200)]
    assert min(heavy) > max(light) * 0.8 and sum(heavy) > 2.0 * sum(light)
