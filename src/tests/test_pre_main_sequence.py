# tests/test_pre_main_sequence.py

"""
GEN.242: the pre-main-sequence state (`physics/pre_main_sequence.py`): the derived track and its contraction time, the
stage a star of a mass and age is in, `evolve_star` and `star_params` for a young star, the disk descriptor on the
`Star` and in `stars`, the companions of a protostar (Perseus' multiplicity) and the words of the star page.
"""

import math
import statistics

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.describe import build as describe_build
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, main_sequence, pre_main_sequence as pms, stellar_abundances as sa
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw

MYR = 1.0e6


@pytest.mark.parametrize("mass, mist_myr", [
    (0.2, (708.0, 794.0)), (0.3, (398.0, 447.0)), (0.5, (141.0, 158.0)), (1.0, (40.0, 45.0)), (1.5, (14.0, 16.0)),
    (2.0, (8.9, 10.0)), (3.0, (3.2, 3.6)), (5.0, (0.9, 1.0)),
])
def test_the_nuclear_arrival_is_the_age_mist_leaves_the_pre_main_sequence(mass, mist_myr):
    """GEN.947: the age at which the MIST isochrones (solar metallicity) hold no contracting star of this mass; the isochrone ages
    are 0.05 dex apart, so the two ends are the ages either side of the arrival and the model is within 0.06 dex of them."""
    low, high = mist_myr
    assert abs(math.log10(pms.nuclear_arrival_yr(mass) / MYR / math.sqrt(low * high))) < 0.06


@pytest.mark.parametrize("mass, bhac15_myr", [
    (0.1, 554.0), (0.08, 700.0),
])
def test_below_the_isochrones_the_nuclear_arrival_is_the_models_own(mass, bhac15_myr):
    """GEN.529: where MIST's rows do not reach the mass the nuclear arrival is the model's own time, within 0.5 dex of the BHAC15 tracks' (the age
    at which the radius is within 2% of its minimum, [R0254], solar metallicity)."""
    assert pms.nuclear_arrival_yr(mass) == pytest.approx(pms._track_own_yr(mass, 0.0))
    assert abs(math.log10(pms.nuclear_arrival_yr(mass) / MYR / bhac15_myr)) < 0.5


@pytest.mark.parametrize("mass, teff_1_myr, radius_first_rsun", [
    (0.1, 2935.0, 1.033), (0.3, 3428.0, 2.215), (0.5, 3803.0, 2.458), (1.0, 4377.0, 3.096), (1.4, 4658.0, 3.621),
])
def test_the_hayashi_temperature_and_birthline_radius_are_the_bhac15_ones(mass, teff_1_myr, radius_first_rsun):
    """The effective temperature at 1 Myr and the radius at the table's first age (0.5 Myr) of the BHAC15 tracks."""
    assert pms.hayashi_temperature_k(mass) == pytest.approx(teff_1_myr, rel=0.08)
    assert pms.birthline_radius_m(mass) / constants.SOLAR_RADIUS_M == pytest.approx(radius_first_rsun, rel=0.3)


def test_the_longest_contraction_shortcut_covers_the_lightest_star():
    lightest = sa.hydrogen_burning_limit_sol(sa.PROTO_Z)       # the least a star is at solar metallicity (GEN.936)
    assert pms.contraction_time_yr(lightest) < tuning.PRE_MAIN_SEQUENCE_LONGEST_YR
    assert pms.contraction_time_yr(lightest) > 0.8 * tuning.PRE_MAIN_SEQUENCE_LONGEST_YR


def test_a_heavier_star_contracts_faster_and_the_old_law_is_gone():
    times = [pms.contraction_time_yr(m) for m in (0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 7.0)]
    assert times == sorted(times, reverse=True)
    assert pms.contraction_time_yr(0.1) < 1.0e9  # the old `31 Myr * m^-2.5` law gave 10 Gyr here


@pytest.mark.parametrize("mass", [0.1, 0.3, 1.0, 3.0, 7.0])
def test_the_track_dims_and_heats_towards_the_main_sequence_without_a_jump(mass):
    end = pms.contraction_time_yr(mass)
    ages = [end * fraction for fraction in (0.001, 0.01, 0.1, 0.3, 0.6, 0.9, 0.999)]
    states = [pms._track_state(mass, age) for age in ages]
    luminosities = [s["luminosity_sol"] for s in states]
    temperatures = [s["temperature_k"] for s in states]
    assert luminosities == sorted(luminosities, reverse=True)
    assert temperatures == sorted(temperatures)
    last = states[-1]
    arrival = main_sequence.zams_state(mass)
    assert last["luminosity_sol"] == pytest.approx(arrival.luminosity_sol, rel=0.02)
    assert last["temperature_k"] == pytest.approx(arrival.temperature_k, rel=0.03)


def test_the_stage_follows_the_mass_and_the_age():
    assert pms.stage(1.0, 0.2 * MYR) == pms.PROTOSTAR
    assert pms.stage(1.0, 2.0 * MYR) == pms.PRE_MAIN_SEQUENCE
    assert pms.stage(1.0, 11.0 * MYR) == pms.PRE_MAIN_SEQUENCE  # GEN.529: no 10 Myr cap, the track takes 37 Myr
    assert pms.stage(1.0, 45.0 * MYR) is None
    assert pms.stage(0.2, 100.0 * MYR) == pms.PRE_MAIN_SEQUENCE  # a red dwarf is still shrinking at 100 Myr
    assert pms.stage(0.2, 150.0 * MYR) == pms.PRE_MAIN_SEQUENCE   # GEN.1174: the 0.2 solar mass radius settles at 188 Myr (BHAC15), not at MIST's 700 to 800
    assert pms.stage(0.2, 400.0 * MYR) is None
    assert pms.stage(0.1, 350.0 * MYR) == pms.PRE_MAIN_SEQUENCE
    assert pms.stage(0.1, 450.0 * MYR) is None
    assert pms.stage(20.0, 0.1 * MYR) is None  # 8 solar masses or more have none
    assert pms.stage(3.0, 1.0 * MYR) == pms.PRE_MAIN_SEQUENCE
    assert pms.stage(3.0, 1.0 * MYR) == pms.PRE_MAIN_SEQUENCE  # GEN.1174: the radius time of 3 solar masses is 1.2 Myr (MIST's nuclear arrival 3.2 to 3.6)
    assert pms.stage(3.0, 4.0 * MYR) is None  # contracted already
    assert pms.stage(5.0, 0.2 * MYR) == pms.PROTOSTAR


def test_evolve_star_gives_a_young_star_its_state_and_leaves_older_ones_alone():
    young = evolution.evolve_star(1.0, 0.0003, fe_h=0.0)
    assert young["yerkes_class"] == "protostar" and young["disk_class"] == "I" and young["radius_km"] is None
    assert young["lifespan_gy"] > 5.0  # the lifespan of the star it will be
    assert evolution.evolve_star(1.0, 0.005, fe_h=0.0)["yerkes_class"] == "PMS"
    for mass, age in ((1.0, 0.5), (1.0, 0.05), (0.2, 0.3), (10.0, 0.0001)):
        assert evolution.evolve_star(mass, age, fe_h=0.0)["yerkes_class"] == "V"
    for mass, age in ((1.0, 0.011), (0.2, 0.1), (0.2, 0.15)):
        assert evolution.evolve_star(mass, age, fe_h=0.0)["yerkes_class"] == "PMS"


def test_the_type_text_names_the_stage_without_a_numeral():
    protostar = evolution.star_params(1.0, 0.0003, evolution.evolve_star(1.0, 0.0003, fe_h=0.0))
    assert protostar["type"].endswith("Class I Protostar") and protostar["yerkes_class"] == "protostar"
    young = evolution.star_params(1.0, 0.005, evolution.evolve_star(1.0, 0.005, fe_h=0.0))
    assert young["type"].endswith("Pre-Main-Sequence Star") and young["yerkes_class"] == "PMS"
    assert young["phase_end_age_gy"] == pytest.approx(pms.contraction_time_yr(1.0) / 1e9)
    assert evolution.star_params(1.0, 0.5, evolution.evolve_star(1.0, 0.5, fe_h=0.0))["type"].startswith("G")


def test_the_disk_class_follows_the_age():
    assert [pms.disk_class(age * MYR) for age in (0.1, 0.3, 1.0, 2.9, 3.1, 8.0)] == ["0", "I", "II", "II", "III", "III"]


def test_a_disk_has_the_ranges_of_the_note_and_an_inner_edge_inside_it():
    for seed in range(200):
        stream = draw.Stream(f"disk-test:{seed}")
        disk = pms.draw_disk(1.0, 1.0 * MYR, 1.5, stream)
        assert disk["disk_kind"] == "class-ii"
        low, high = tuning.DISK_MASS_FRACTION_RANGE
        assert low <= disk["disk_mass_solar"] <= high
        assert tuning.DISK_OUTER_RADIUS_RANGE_AU[0] <= disk["disk_outer_au"] <= tuning.DISK_OUTER_RADIUS_RANGE_AU[1]
        assert disk["disk_inner_au"] < disk["disk_outer_au"]
    assert pms.draw_disk(1.0, 5.0 * MYR, 1.0) is None  # Class III: the disk is gone
    assert pms.draw_disk(1.0, 0.1 * MYR, 1.5, draw.Stream("a"))["disk_kind"] == "class-0"


def test_the_inner_edge_is_where_a_blackbody_grain_sublimates():
    """For a solar luminosity and 1500 K that is 0.034 AU (the blackbody radius `sqrt(L / 16 pi sigma T^4)`)."""
    assert pms.sublimation_radius_au(1.0) == pytest.approx(0.034, rel=0.03)
    assert pms.sublimation_radius_au(100.0) == pytest.approx(10 * pms.sublimation_radius_au(1.0))


def test_a_protostar_star_holds_a_disk_and_a_main_sequence_one_holds_no_young_one():
    draw.set_run_seed(5)
    config = SystemConfig()
    protostar = Star.from_params(evolution.star_params(1.0, 0.0003, evolution.evolve_star(1.0, 0.0003)), config,
                                 name="Young")
    assert protostar.yerkes_class == "protostar" and protostar.disk_kind == "class-i"
    assert protostar.disk_inner_au < protostar.disk_outer_au
    assert 0.01 <= protostar.disk_mass_solar / (protostar.mass / 1.98840987e30) <= 0.3
    rebuilt = Star.from_dict(protostar.to_dict(), config)
    assert (rebuilt.disk_kind, rebuilt.disk_mass_solar, rebuilt.disk_outer_au) == (
        protostar.disk_kind, protostar.disk_mass_solar, protostar.disk_outer_au)
    again = Star.from_params(evolution.star_params(1.0, 0.0003, evolution.evolve_star(1.0, 0.0003)), config, name="Young")
    assert again.disk_outer_au == protostar.disk_outer_au  # a stream of its own: the same star, the same disk
    old = Star.from_params(evolution.star_params(1.0, 4.6, evolution.evolve_star(1.0, 4.6)), config, name="Old")
    assert old.disk_kind in (None, "debris") and old.disk_mass_solar is None    # at most the cold belt of GEN.787, which has no mass


def test_the_star_page_says_the_star_is_contracting_and_shows_the_disk():
    draw.set_run_seed(5)
    config = SystemConfig()
    star = Star.from_params(evolution.star_params(1.0, 0.0003, evolution.evolve_star(1.0, 0.0003)), config, name="Young")
    rows = {row["key"]: row for row in describe_build._star_facts(star)["rows"]}
    assert "disk" in rows and rows["disk"]["label"] == "Disk"
    assert "contracting towards the main sequence" in str(describe_build._star_age_sentence(star))


def test_a_protostar_has_the_companions_of_the_perseus_survey():
    """[R4379]: 57% of Class 0 and 23% of Class I protostars have a companion (the field's is 26 to 45%)."""
    assert pms.companion_chance(1.0, 0.1 * MYR) == 0.57
    assert pms.companion_chance(1.0, 0.3 * MYR) == 0.23
    assert pms.companion_chance(1.0, 2.0 * MYR) is None  # a pre-main-sequence star keeps the field's rate
    assert pms.companion_chance(1.0, 1.0e9) is None
    shares = {}
    for label, age_gy in (("0", 0.0001), ("I", 0.0003)):
        found = 0
        trials = 400
        for seed in range(trials):
            draw.set_run_seed(seed)
            config = SystemConfig()
            config.PLANETS = False
            params = evolution.star_params(1.0, age_gy, evolution.evolve_star(1.0, age_gy))
            system = StarSystem(system_config=config, primary_star_params=params)
            found += system.binary_type is not None
        shares[label] = found / trials
    assert 0.47 < shares["0"] < 0.67 and 0.15 < shares["I"] < 0.31, shares


def test_a_protostars_system_holds_no_planets():
    for seed in range(20):
        draw.set_run_seed(seed)
        params = evolution.star_params(1.0, 0.0003, evolution.evolve_star(1.0, 0.0003))
        system = StarSystem(system_config=SystemConfig(), primary_star_params=params)
        assert [p for p in system.planets if p.body_type != "a"] == []


def test_the_disk_is_stored_with_the_star_and_comes_back(mysql_config):
    draw.set_run_seed(9)
    config = SystemConfig()
    config.PLANETS = False
    config.BINARY_SYSTEM = False
    params = evolution.star_params(1.0, 0.0003, evolution.evolve_star(1.0, 0.0003))
    system = StarSystem(system_config=config, primary_star_params=params)
    conn = store.get_connection(mysql_config)
    try:
        system_id = store.insert_star_system(conn, system, config)
        conn.commit()
        row = conn.execute("SELECT yerkes_class, disk_kind, disk_mass_solar, disk_inner_au, disk_outer_au "
                           "FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()
        assert row["yerkes_class"] == "protostar" and row["disk_kind"] == "class-i"
        assert row["disk_outer_au"] == system.star.disk_outer_au
        loaded = store.load_star_system(conn, system_id)
        assert loaded.star.disk_kind == "class-i" and loaded.star.disk_mass_solar == system.star.disk_mass_solar
    finally:
        conn.close()


def test_a_population_draw_is_rarely_young():
    """Under 10 Myr a thousandth of the disk's stars are young; with the track taking a red dwarf 0.1 to 0.7 Gyr
    (GEN.529) about 3.6% are still contracting, nearly all of them red dwarfs; 8 solar masses and over have no state."""
    import random
    rng = random.Random(7)
    draws = [evolution.sample_living_star(rng=rng) for _ in range(40000)]
    young = [draw_[2]["yerkes_class"] in ("PMS", "protostar") for draw_ in draws]
    assert 0.015 < statistics.mean(young) < 0.06
    masses = [draw_[2]["mass_sol"] for draw_, is_young in zip(draws, young) if is_young]
    assert statistics.median(masses) < 0.5


@pytest.mark.parametrize("mass, bhac15_myr", [
    (0.08, 830.5), (0.1, 402.4), (0.15, 240.1), (0.2, 187.5), (0.3, 130.2), (0.5, 88.2), (0.8, 50.3),
])
def test_the_contraction_time_is_the_age_the_bhac15_radius_settles(mass, bhac15_myr):
    """GEN.1174: the age at which the radius of the BHAC15 track ([R0254], solar metallicity) is within 5 percent of its value at 5 Gyr."""
    assert pms.contraction_time_yr(mass) / MYR == pytest.approx(bhac15_myr, rel=1e-6)
    assert pms.radius_settled_yr(mass) == pytest.approx(pms.contraction_time_yr(mass))


def test_the_radius_settled_time_falls_with_the_mass_without_a_step():
    masses = [0.07 + 0.005 * i for i in range(190)]      # to 1.0; past it the model's own time (MIST rows) has small bumps
    times = [pms.radius_settled_yr(m) for m in masses]
    assert all(later <= earlier * (1 + 1e-12) for earlier, later in zip(times, times[1:]))
    assert max(math.log(a / b) for a, b in zip(times, times[1:])) < 0.3      # the steepest is the BHAC15 fall from 830 to 507 Myr between 0.08 and 0.09


def test_the_track_ends_on_the_main_sequence_when_the_radius_settles():
    for mass in (0.1, 0.3, 1.0):
        end = pms._track_state(mass, pms.contraction_time_yr(mass))
        arrival = main_sequence.zams_state(mass, 0.0)
        assert end["luminosity_sol"] == pytest.approx(arrival.luminosity_sol, rel=0.02)
        assert end["radius_sol"] == pytest.approx(arrival.radius_sol, rel=0.02)


def test_the_nuclear_arrival_is_kept_beside_the_radius_time():
    state = pms.state(0.2, 100.0 * MYR)
    assert state["phase_end_gy"] == pytest.approx(pms.contraction_time_yr(0.2) / 1e9)
    assert state["nuclear_arrival_gy"] == pytest.approx(pms.nuclear_arrival_yr(0.2) / 1e9)
    assert state["nuclear_arrival_gy"] > 3.0 * state["phase_end_gy"]      # 0.73 Gyr against 0.19 at 0.2 solar masses
