"""
GEN.935: a star after its main sequence reads its luminosity, temperature and radius from the MIST tracks
(`physics/post_main_sequence.py`, `datasets/evolved_tracks.py`) instead of drawing them, and the figures of the 2026-10-29
star class audit (reports/star-class-audit-2026-10-29.md, sections 1 and 6) hold: the first ascent of 1 and 2 solar mass
stars (MIST 28 Lsun at 4,507 K and 8.7 Rsun; 164 Lsun at 4,383 K and 22 Rsun) and its tip ([R4024: Table 2]:
1,650 Lsun at 4,550 K for Z = 0.0001, 2,410 Lsun at 3,640 K and 123 Rsun for Z = 0.004), core helium burning at 3, 5 and 8
solar masses (MIST 113 to 144 Lsun, 1,507 Lsun, 6,477 Lsun at 4,097 K), and no red supergiant above the Humphreys-Davidson limit
([R1013]: log L about 5.5, about 320,000 Lsun).
"""

import math
import random

import numpy as np
import pytest

from planetgen import tuning
from planetgen.generation import star_population as sp
from planetgen.physics import constants, imf, main_sequence, post_main_sequence as pm, remnant_mass
from planetgen.physics import stellar_evolution as se


def _at(mass, window, share, fe_h=0.0):
    """The state `share` of the way through a window of a star."""
    line = pm.timeline(mass, fe_h)
    w = line.window(window)
    return line.state(w.start_gy + share * (w.end_gy - w.start_gy))


def _radius_sol(state):
    return state.radius_sol


@pytest.mark.parametrize("mass,fe_h", [(0.9, -1.5), (1.0, 0.0), (1.8, 0.0), (2.2, 0.0), (3.0, 0.0), (5.0, -0.5), (7.9, 0.0),
                                        (8.0, 0.0), (15.0, 0.0), (40.0, 0.3), (120.0, 0.0)])
def test_the_windows_follow_each_other_from_the_end_of_the_main_sequence_to_the_end_of_the_track(mass, fe_h):
    line = pm.timeline(mass, fe_h)
    assert line.windows[0].start_gy == pytest.approx(se.main_sequence_lifetime_gy(mass, fe_h))
    for before, after in zip(line.windows, line.windows[1:]):
        assert before.end_gy == pytest.approx(after.start_gy) and before.f_end == pytest.approx(after.f_start, rel=1e-6, abs=1e-9)
    assert all(w.end_gy >= w.start_gy and w.f_end >= w.f_start for w in line.windows)
    names = [w.name for w in line.windows]
    expected = ["supergiant"] if mass >= tuning.SUPERGIANT_MIN_MASS_SOL else (
        ["subgiant", "giant", "horizontal-branch", "asymptotic-giant"] if mass < tuning.HELIUM_FLASH_MASS_SOL else
        ["subgiant", "giant", "core-helium-burning", "asymptotic-giant"])
    assert names[:len(expected)] == expected or (mass >= 7.0 and names[:3] == expected[:3])
    assert 1.03 < se.remnant_age_fraction(mass, fe_h) < 1.4


def test_core_helium_burning_giants_have_the_luminosities_of_the_mist_tracks():
    three = _at(3.0, "core-helium-burning", 0.5)
    assert 100.0 < three.luminosity_sol < 160.0 and 4800.0 < three.temperature_k < 5200.0
    assert 12.0 < three.radius_sol < 18.0
    five = _at(5.0, "core-helium-burning", 0.5)
    assert 1200.0 < five.luminosity_sol < 1900.0 and 40.0 < five.radius_sol < 80.0
    eight = _at(8.0, "supergiant", 0.5)
    assert 5500.0 < eight.luminosity_sol < 8500.0 and 3800.0 < eight.temperature_k < 4500.0 and 120.0 < eight.radius_sol < 200.0


def test_the_hertzsprung_gap_is_a_small_part_of_the_life_of_a_star_of_two_and_a_half_solar_masses_or_more():
    for mass in (2.5, 3.0, 5.0, 7.0):
        line = pm.timeline(mass, 0.0)
        subgiant = line.window("subgiant")
        share = (subgiant.end_gy - subgiant.start_gy) / line.main_sequence_gy
        assert 0.001 < share < 0.01, mass


def test_the_first_ascent_of_a_one_and_a_two_solar_mass_star_is_a_line_in_luminosity_and_temperature():
    for mass, luminosity, temperature, radius in ((1.0, 28.0, 4507.0, 8.7), (2.0, 164.0, 4383.0, 22.0)):
        line = pm.timeline(mass, 0.0)
        giant = line.window("giant")
        found = None
        for share in np.linspace(0.0, 1.0, 4001)[:-1]:
            state = line.state(giant.start_gy + share * (giant.end_gy - giant.start_gy))
            if state.luminosity_sol >= luminosity:
                found = state
                break
        assert found.temperature_k == pytest.approx(temperature, rel=0.03) and found.radius_sol == pytest.approx(radius, rel=0.06)
    # The temperature falls as the star climbs.
    line = pm.timeline(1.0, 0.0)
    giant = line.window("giant")
    cool = [line.state(giant.start_gy + share * (giant.end_gy - giant.start_gy)) for share in np.linspace(0.0, 0.999, 60)]
    temperatures = [state.temperature_k for state in cool]
    assert temperatures[0] > temperatures[-1] + 1000.0
    assert all(b <= a + 1.0 for a, b in zip(temperatures, temperatures[1:]))


def test_the_tip_of_the_giant_branch_gets_larger_and_cooler_toward_solar_metallicity():
    """[R4024: Table 2] (1 solar mass): Z = 0.0001 gives 1,650 Lsun at 4,550 K (65 Rsun), Z = 0.004 gives 2,410 Lsun at
    3,640 K (123 Rsun). The code drew at most 54 and 65 Rsun before; the MIST tracks, at the matching [Fe/H], are within 12%."""
    poor = pm.timeline(1.0, math.log10(1.0e-4 / tuning.MS_LIFETIME_REFERENCE_Z))
    rich = pm.timeline(1.0, math.log10(4.0e-3 / tuning.MS_LIFETIME_REFERENCE_Z))
    tips = []
    for line in (poor, rich):
        giant = line.window("giant")
        tips.append(line.state(giant.end_gy - 1.0e-9))
    assert tips[0].radius_sol == pytest.approx(65.0, rel=0.12) and tips[0].temperature_k == pytest.approx(4550.0, rel=0.05)
    assert tips[1].radius_sol == pytest.approx(123.0, rel=0.12) and tips[1].temperature_k == pytest.approx(3640.0, rel=0.05)
    assert tips[1].radius_sol > 1.7 * tips[0].radius_sol


def test_no_red_supergiant_is_brighter_than_the_humphreys_davidson_limit_and_none_is_smaller_than_it_was_leaving_the_main_sequence():
    rng = random.Random(3)
    red_above = smaller = total = 0
    for _ in range(3000):
        mass = imf.sample(8.0, 150.0, rng)
        line = pm.timeline(mass, 0.0)
        # (Past the first percent of the window: a star contracts a little at the very end of its hydrogen burning, the hook.)
        age = line.main_sequence_gy + rng.uniform(0.01, 1.0) * (line.end_gy - line.main_sequence_gy)
        state = se.evolve_star(mass, age, rng, fe_h=0.0)
        if state is None or state["yerkes_class"] == "W":
            continue
        total += 1
        luminosity, temperature = state["luminosity_sol"], state["temperature_k"]
        red_above += temperature < 4500.0 and luminosity > 3.2e5
        radius = math.sqrt(luminosity) * (tuning.SUN_EFFECTIVE_TEMPERATURE_K / temperature) ** 2
        smaller += mass < 45.0 and radius < 0.9 * main_sequence.end_state(mass, 0.0).radius_sol
    assert total > 2000 and red_above == 0 and smaller == 0


def test_a_supergiant_cools_from_the_hot_end_of_the_main_sequence_to_a_red_one_below_the_limit_and_a_hot_one_above():
    """A star that has left the main sequence is cooler than it was; a 15 solar mass star ends as a red supergiant, a 60 solar mass
    one never reaches the red."""
    line = pm.timeline(15.0, 0.0)
    window = line.window("supergiant")
    early, late = line.state(window.start_gy + 1.0e-6), line.state(window.end_gy - 1.0e-6)
    assert early.temperature_k > 15000.0 > 4500.0 > late.temperature_k and late.radius_sol > 300.0
    hot = pm.timeline(60.0, 0.0)
    states = [hot.state(hot.windows[0].start_gy + share * (hot.end_gy - hot.windows[0].start_gy)) for share in np.linspace(0, 0.999, 80)]
    assert min(state.temperature_k for state in states) > 7000.0


def test_stars_inside_the_instability_strip_are_flagged_as_classical_cepheids_only_when_luminous_and_of_the_right_mass():
    rng = random.Random(8)
    line = pm.timeline(5.0, 0.0)
    flagged = clear = 0
    for age in np.linspace(line.main_sequence_gy, line.end_gy * 0.999999, 600):
        state = se.evolve_star(5.0, float(age), rng, fe_h=0.0)
        in_strip = tuning.CEPHEID_TEMPERATURE_RANGE_K[0] <= state["temperature_k"] <= tuning.CEPHEID_TEMPERATURE_RANGE_K[1]
        if state.get("phase") == "classical-cepheid":
            flagged += 1
            assert in_strip and state["luminosity_sol"] >= tuning.CEPHEID_MIN_LUMINOSITY_SOL
            params = se.star_params(5.0, float(age), state)
            assert "Classical Cepheid Star" in params["type"]
        elif in_strip and state["luminosity_sol"] >= tuning.CEPHEID_MIN_LUMINOSITY_SOL and state["yerkes_class"] != "IV":
            clear += 1
    assert flagged > 10 and clear == 0
    # The yellow supergiants of 20 and 40 solar masses are supergiants, not Cepheids.
    for mass in (20.0, 40.0):
        line = pm.timeline(mass, 0.0)
        for age in np.linspace(line.main_sequence_gy, line.end_gy * 0.999, 200):
            state = se.evolve_star(mass, float(age), rng, fe_h=0.0)
            assert state is None or state.get("phase") != "classical-cepheid"


def test_the_bright_pieces_are_exactly_the_ages_the_star_is_that_bright():
    rng = random.Random(5)
    for mass, floor in ((0.95, 300.0), (1.5, 100.0), (3.0, 150.0), (5.0, 1400.0), (12.0, 30000.0), (30.0, 3.0e5)):
        line = pm.timeline(mass, 0.0)
        pieces = line.bright_pieces(floor, skip=("horizontal-branch",))
        assert pieces, mass
        for _ in range(200):
            age = rng.uniform(line.main_sequence_gy, line.end_gy)
            if line.window_at(age).name == "horizontal-branch":
                continue
            inside = any(start <= age <= end for _, start, end in pieces)
            luminosity = line.state(age).luminosity_sol
            if inside:
                assert luminosity >= floor * (1.0 - 1.0e-3), (mass, age)
            else:
                assert luminosity <= floor * (1.0 + 1.0e-3), (mass, age)


def test_the_bright_star_tables_mirror_the_tracks_for_helium_burning_giants_and_red_supergiants():
    """A scatter that asks for stars of at least 1,000 Lsun draws core-helium-burning giants of 5 solar masses and up (the old
    model had a giant luminosity of 1,000 to 10,000 Lsun for them whatever the age), and one for 100,000 Lsun draws red
    supergiants, every one at least that bright."""
    types = [star["type"] for star in sp.sample_bright_stars(300, 1000.0, "intermediate", random.Random(2))]
    assert any("Core Helium Burning Giant" in kind for kind in types)
    young = sp.sample_bright_stars(300, 1.0e5, "young", random.Random(2))
    assert any("Supergiant" in star["type"] and "Red" in star["type"] for star in young)
    assert all(star["luminosity_w"] >= 1.0e5 * constants.SOLAR_LUMINOSITY * 0.999 for star in young)


def test_a_thirty_five_solar_mass_star_is_a_yellow_hypergiant_for_part_of_its_track():
    line = pm.timeline(35.0, 0.0)
    window = line.windows[0]
    temperatures = [line.state(float(age)).temperature_k for age in np.linspace(window.start_gy, window.end_gy * 0.999, 400)]
    assert any(5000.0 <= t <= 9000.0 for t in temperatures)


def test_the_phase_end_and_lifespan_say_when_the_star_leaves_its_window_and_becomes_a_remnant():
    rng = random.Random(4)
    line = pm.timeline(3.0, 0.0)
    helium = line.window("core-helium-burning")
    state = se.evolve_star(3.0, (helium.start_gy + helium.end_gy) / 2.0, rng, fe_h=0.0)
    assert state["phase"] == "core-helium-burning" and state["yerkes_class"] == "III"
    assert state["phase_end_gy"] == pytest.approx(helium.end_gy) and state["lifespan_gy"] == pytest.approx(line.end_gy)
    bright = se.evolve_star(5.0, pm.timeline(5.0, 0.0).window("core-helium-burning").start_gy * 1.01, rng, fe_h=0.0)
    assert bright["yerkes_class"] == "II"
    assert se.evolve_star(3.0, line.end_gy * 1.001, rng, fe_h=0.0)["yerkes_class"] == "VII"
    assert se.evolve_star(15.0, pm.timeline(15.0, 0.0).end_gy * 1.001, rng, fe_h=0.0) is None


def test_the_youngest_progenitor_is_the_lightest_star_whose_track_has_ended_by_that_age():
    for age_yr in (1.0e8, 5.0e8, 3.0e9, 1.3e10):
        mass = se.youngest_progenitor_sol(age_yr)
        assert se.remnant_formation_time_gy(mass) * 1.0e9 == pytest.approx(age_yr, rel=1e-6)
        assert se.remnant_formation_time_gy(mass * 0.97) * 1.0e9 > age_yr


def _cluster(age_gy, low, high, count=6000):
    """The stars of a coeval cluster of masses `low` to `high` (a fine grid, not a draw), those still shining: `(mass, state)`."""
    rng = random.Random(3)
    found = []
    for mass in np.exp(np.linspace(math.log(low), math.log(high), count)):
        state = se.evolve_star(float(mass), age_gy, rng, fe_h=0.0)
        if state is not None and state["yerkes_class"] != "VII":
            found.append((float(mass), state))
    return found


def test_a_one_gigayear_cluster_has_a_red_clump_and_a_ten_gigayear_one_a_horizontal_branch_and_a_cool_bright_asymptotic_branch():
    """The synthetic cluster check of the diagram (GEN.741 folded into GEN.935): the helium-burning stars of an intermediate-age
    cluster sit in one clump, a globular cluster has its horizontal branch below 51 Lsun and an asymptotic branch that is cool and
    runs up past 1,000 Lsun, and no star is hotter than its main sequence allows while on the way."""
    clump = [s for _, s in _cluster(1.0, 1.6, 2.6) if s.get("phase") == "core-helium-burning"]
    assert len(clump) > 100
    assert all(40.0 <= s["luminosity_sol"] <= 125.0 and 4500.0 <= s["temperature_k"] <= 5100.0 for s in clump)
    old = _cluster(10.0, 0.8, 1.1)
    branch = [s for _, s in old if s.get("phase") == "horizontal-branch"]
    asymptotic = [s for _, s in old if s.get("phase") == "asymptotic-giant"]
    assert branch and all(35.0 <= s["luminosity_sol"] <= 50.0 for s in branch)
    assert asymptotic and max(s["luminosity_sol"] for s in asymptotic) > 1000.0
    assert all(s["temperature_k"] < 4700.0 for s in asymptotic)
    young = [s for _, s in _cluster(0.1, 3.0, 6.0) if s.get("phase") in ("core-helium-burning", "classical-cepheid")]
    assert young and all(1000.0 <= s["luminosity_sol"] <= 3000.0 for s in young)
    assert all(5000.0 <= s["temperature_k"] <= 7000.0 for s in young if s["phase"] == "classical-cepheid")


@pytest.mark.parametrize("mass, age_share", [(3.0, 1.0150), (1.0, 1.08)])
def test_the_planets_of_a_giant_read_its_new_luminosity_and_radius(mass, age_share):
    """A giant's system is built from the state of the track: its luminosity and temperature put the habitable zone where the
    track's figures say (well out from the main-sequence one), and no planet is inside the radius the star swallows."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.star import Star
    from planetgen.generation.system import StarSystem
    from planetgen.physics.orbits import calculate_habitable_zone
    t_ms = se.main_sequence_lifetime_gy(mass)
    state = se.evolve_star(mass, t_ms * age_share, random.Random(1), fe_h=0.0)
    assert state["yerkes_class"] in ("III", "II", "IV")
    params = se.star_params(mass, t_ms * age_share, state)
    star = Star.from_params(params, SystemConfig(), name="Giant")
    assert star.luminosity == params["luminosity_w"] and star.radius == params["radius_km"]
    main = se.evolve_star(mass, t_ms * 0.5, random.Random(1), fe_h=0.0)
    dwarf = Star.from_params(se.star_params(mass, t_ms * 0.5, main), SystemConfig(), name="Dwarf")
    assert star.habitable_zone == calculate_habitable_zone(star.luminosity, star.temperature)
    assert star.habitable_zone[0] > 1.5 * dwarf.habitable_zone[0]
    for seed in range(5):
        system = StarSystem(system_config=SystemConfig(), primary_star_params=params)
        engulfed_au = system._engulfment_radius_au(system.primary_star)
        for body in system.planets:
            assert (body.lower_limit if body.body_type == "a" else body.distance) >= engulfed_au


# ---- GEN.940: the present mass ----

def _present_masses(mass, per_window=40):
    line = pm.timeline(mass, 0.0)
    found = []
    for window in line.windows:
        for share in np.linspace(0.0, 0.9999999, per_window):
            found.append(line.state(window.start_gy + share * (window.end_gy - window.start_gy)).mass_sol)
    return found


@pytest.mark.parametrize("mass", [0.85, 0.9, 1.0, 1.5, 2.0, 3.0, 5.0, 6.0, 7.0, 7.9, 8.0, 15.0, 40.0])
def test_a_star_never_gains_mass_after_its_main_sequence(mass):
    masses = _present_masses(mass)
    assert all(later <= earlier * (1.0 + 1.0e-4) for earlier, later in zip(masses, masses[1:])), mass
    assert masses[0] <= mass * 1.0001 and masses[-1] > 0.0


@pytest.mark.parametrize("mass", [0.85, 0.9, 1.0, 1.5, 2.0, 3.0, 5.0, 6.0, 7.0, 7.9])
def test_a_star_of_up_to_eight_solar_masses_ends_its_last_giant_phase_at_the_mass_of_the_white_dwarf_it_becomes(mass):
    dwarf = remnant_mass.white_dwarf_mass_sol(mass)
    assert _present_masses(mass)[-1] == pytest.approx(dwarf, rel=0.1, abs=0.05)
    assert min(_present_masses(mass)) >= dwarf * 0.9999


def test_the_red_giant_branch_loses_the_mass_of_a_reimers_wind_at_the_fitted_efficiency_not_mists():
    """MIST's eta_R = 0.1 takes 0.046 solar masses of a one solar mass star by the tip; McDonald and Zijlstra's 0.477 takes 4.77
    times as much, so the horizontal-branch star of an old cluster is 0.55 to 0.7 solar masses (globular clusters' turn-offs
    are 0.8 to 0.9)."""
    line = pm.timeline(1.0, 0.0)
    tip = line.window("giant")
    at_tip = line.state(tip.end_gy - 1.0e-9)
    mist_loss = 1.0 - line.track.star_mass[line.track.index_of(605)]
    scale = tuning.RED_GIANT_BRANCH_REIMERS_ETA / tuning.MIST_REIMERS_ETA
    assert 1.0 - at_tip.mass_sol == pytest.approx(scale * mist_loss, rel=0.02)
    for turn_off in (0.85, 0.9):
        branch = pm.timeline(turn_off, 0.0).window("horizontal-branch")
        mass = pm.timeline(turn_off, 0.0).state((branch.start_gy + branch.end_gy) / 2.0).mass_sol
        assert 0.55 <= mass <= 0.7, (turn_off, mass)


def test_a_star_lighter_than_the_lightest_track_keeps_its_own_mass_on_the_giant_branch():
    line = pm.timeline(0.85, 0.0)
    assert line.state(line.windows[0].start_gy).mass_sol == pytest.approx(0.85, rel=0.01)


def test_massive_stars_keep_the_winds_of_the_track_and_lose_a_tenth_to_two_fifths_of_their_mass():
    rng = random.Random(2)
    for mass, low, high in ((20.0, 0.60, 0.95), (40.0, 0.5, 0.9)):
        line = pm.timeline(mass, 0.0)
        end = line.state(line.end_gy - 1.0e-9).mass_sol / mass
        assert low <= end <= high, (mass, end)
        on_sequence = se.evolve_star(mass, line.main_sequence_gy * 0.9, rng, fe_h=0.0)
        assert on_sequence["mass_sol"] < mass and on_sequence["mass_sol"] > 0.8 * mass
    sun_like = se.evolve_star(1.0, 4.6, rng, fe_h=0.0)
    assert sun_like["mass_sol"] == pytest.approx(1.0, rel=2.0e-3)


def test_the_mass_of_a_star_in_its_record_is_the_present_one():
    rng = random.Random(3)
    mass, age = 1.0, pm.timeline(1.0, 0.0).windows[2].start_gy + 0.01
    state = se.evolve_star(mass, age, rng, fe_h=0.0)
    params = se.star_params(mass, age, state)
    assert params["mass_kg"] == pytest.approx(state["mass_sol"] * constants.SOLAR_MASS_TO_KG)
    assert params["mass_kg"] < 0.9 * mass * constants.SOLAR_MASS_TO_KG and params["initial_mass_sol"] == mass
