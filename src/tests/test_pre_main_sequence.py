# tests/test_pre_main_sequence.py

"""
GEN.242: the pre-main-sequence state (`physics/pre_main_sequence.py`): the derived track and its contraction time, the
stage a star of a mass and age is in, `evolve_star` and `star_params` for a young star, the disk descriptor on the
`Star` and in `stars`, the companions of a protostar (Perseus' multiplicity) and the words of the star page.
"""

import statistics

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.describe import build as describe_build
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import pre_main_sequence as pms
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw

MYR = 1.0e6


@pytest.mark.parametrize("mass, low_myr, high_myr", [
    (0.1, 500.0, 2000.0), (0.3, 100.0, 600.0), (0.5, 50.0, 250.0), (1.0, 20.0, 80.0), (2.0, 4.0, 15.0), (5.0, 0.3, 3.0),
])
def test_the_contraction_time_is_near_what_published_tracks_give(mass, low_myr, high_myr):
    """The anchors are remembered, not tabulated: a solar mass reaches the main sequence in some 40 Myr and 0.1 solar
    masses in about 1 Gyr, the old `31 Myr * m^-2.5` law gave 10 Gyr there."""
    assert low_myr < pms.contraction_time_yr(mass) / MYR < high_myr


def test_a_heavier_star_contracts_faster_and_the_old_law_is_gone():
    times = [pms.contraction_time_yr(m) for m in (0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 7.0)]
    assert times == sorted(times, reverse=True)
    assert pms.contraction_time_yr(0.1) < 2.0e9


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
    assert last["luminosity_sol"] == pytest.approx(evolution.main_sequence_luminosity_sol(mass), rel=0.02)
    assert last["temperature_k"] == pytest.approx(
        evolution.effective_temperature_k(evolution.main_sequence_luminosity_sol(mass),
                                          evolution.main_sequence_radius_sol(mass)), rel=0.03)


def test_the_stage_follows_the_mass_and_the_age():
    assert pms.stage(1.0, 0.2 * MYR) == pms.PROTOSTAR
    assert pms.stage(1.0, 2.0 * MYR) == pms.PRE_MAIN_SEQUENCE
    assert pms.stage(1.0, 11.0 * MYR) is None  # no state is drawn at 10 Myr or more
    assert pms.stage(0.2, 9.0 * MYR) == pms.PRE_MAIN_SEQUENCE
    assert pms.stage(0.2, 12.0 * MYR) is None
    assert pms.stage(20.0, 0.1 * MYR) is None  # 8 solar masses or more have none
    assert pms.stage(5.0, 0.7 * MYR) == pms.PRE_MAIN_SEQUENCE
    assert pms.stage(5.0, 0.9 * MYR) is None  # contracted already: the model takes 0.87 Myr
    assert pms.stage(5.0, 0.2 * MYR) == pms.PROTOSTAR


def test_evolve_star_gives_a_young_star_its_state_and_leaves_older_ones_alone():
    young = evolution.evolve_star(1.0, 0.0003)
    assert young["yerkes_class"] == "protostar" and young["disk_class"] == "I" and young["radius_km"] is None
    assert young["lifespan_gy"] > 5.0  # the lifespan of the star it will be
    assert evolution.evolve_star(1.0, 0.005)["yerkes_class"] == "PMS"
    for mass, age in ((1.0, 0.5), (1.0, 0.011), (0.2, 0.05), (10.0, 0.0001)):
        assert evolution.evolve_star(mass, age)["yerkes_class"] == "V"


def test_the_type_text_names_the_stage_without_a_numeral():
    protostar = evolution.star_params(1.0, 0.0003, evolution.evolve_star(1.0, 0.0003))
    assert protostar["type"].endswith("Class I Protostar") and protostar["yerkes_class"] == "protostar"
    young = evolution.star_params(1.0, 0.005, evolution.evolve_star(1.0, 0.005))
    assert young["type"].endswith("Pre-Main-Sequence Star") and young["yerkes_class"] == "PMS"
    assert young["phase_end_age_gy"] == pytest.approx(tuning.PRE_MAIN_SEQUENCE_MAX_AGE_YR / 1e9)
    assert evolution.star_params(1.0, 0.5, evolution.evolve_star(1.0, 0.5))["type"].startswith("G2V")


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


def test_a_protostar_star_holds_a_disk_and_a_main_sequence_one_does_not():
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
    assert old.disk_kind is None and old.disk_mass_solar is None


def test_the_star_page_says_the_star_is_contracting_and_shows_the_disk():
    draw.set_run_seed(5)
    config = SystemConfig()
    star = Star.from_params(evolution.star_params(1.0, 0.0003, evolution.evolve_star(1.0, 0.0003)), config, name="Young")
    rows = {row["key"]: row for row in describe_build._star_facts(star)["rows"]}
    assert "disk" in rows and rows["disk"]["label"] == "Disk"
    assert "contracting towards the main sequence" in str(describe_build._star_age_sentence(star))
    assert "contracting" in " ".join(star.to_paragraph_list())
    assert "Disk" in star.get_table_properties()["disk"].capitalize() or "disk" in star.get_table_properties()["disk"]


def test_a_protostar_has_the_companions_of_the_perseus_survey():
    """Tobin et al. 2016: 57% of Class 0 and 23% of Class I protostars have a companion (the field's is 26 to 45%)."""
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
    """About a thousandth of the disk's stars are under 10 Myr; the share of those drawn as contracting is the
    lighter ones' (8 solar masses and over have no state)."""
    import random
    rng = random.Random(7)
    young = [evolution.sample_living_star(rng=rng)[2]["yerkes_class"] in ("PMS", "protostar") for _ in range(40000)]
    assert 0.0002 < statistics.mean(young) < 0.003
