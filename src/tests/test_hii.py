"""
H II regions from their stars' photons (GEN.150): `physics.hii` and the
nebulae built on it.

Run with: pytest tests/test_hii.py
"""
import math

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.physics import hii, ms_lifetime
from planetgen.util import draw


def test_stromgren_radii_match_the_design_table():
    """Design note, n = 100 per cc: O3V 17 ly, O5V 12, O7V 8.8, O9V 5.6, B0V 3.5, B1V 1.9, B2V 0.8."""
    expected = {"O3V": 17.0, "O5V": 12.0, "O7V": 8.8, "O9V": 5.6, "B0V": 3.5, "B1V": 1.9, "B2V": 0.8}
    for kind, anchor_ly in expected.items():
        q = 10.0 ** next(log_q for name, _m, log_q in hii.TYPE_ANCHORS if name == kind)
        radius_ly = hii.stromgren_radius_pc(q, 100.0) * hii.LY_PER_PC
        assert radius_ly == pytest.approx(anchor_ly, rel=0.1), kind


def test_photon_rate_falls_with_spectral_type_and_stops_at_b5():
    rates = [hii.log_q_of_type("O", s) for s in (3, 5, 7, 9)] + [hii.log_q_of_type("B", s) for s in (0, 1, 2, 3, 5)]
    assert rates == sorted(rates, reverse=True)
    assert hii.log_q_of_type("B", 6) is None
    assert hii.log_q_of_type("A", 0) is None
    assert hii.log_q_of_mass(hii.MIN_IONIZING_MASS_SOL - 0.1) is None
    assert hii.log_q_of_mass(150.0) == hii.log_q_of_type("O", 3)


def test_the_front_keeps_growing_and_thins_the_gas():
    q = 1e49
    r_young, n_young = hii.region(q, 100.0, 0.1)
    r_old, n_old = hii.region(q, 100.0, 3.0)
    assert r_old > r_young
    assert n_old < n_young <= 100.0
    assert r_old == pytest.approx(54.0, rel=0.1)
    assert hii.region(0.0, 100.0, 3.0) == (0.0, 100.0)


def test_cluster_photons_per_solar_mass_follow_the_imf():
    draw.set_run_seed(3)
    mass = 1e5
    q = sum(10.0 ** hii.log_q_of_mass(m) for m in hii.cluster_hosts(mass, 0.03))  # before any of them has left the main sequence
    assert q / mass == pytest.approx(5.7e46, rel=0.4)


def test_old_clusters_lose_their_hottest_stars():
    draw.set_run_seed(4)
    young = hii.cluster_hosts(1e5, 0.2)
    old = hii.cluster_hosts(1e5, 10.0)
    assert max(young) > 40.0
    assert max(old) < ms_lifetime.mass_for_lifetime_gy(0.010) * 1.001  # the lifetime fit (GEN.467): about 18 Msun
    assert len(old) < len(young)


def test_host_summary_lists_the_hottest_first():
    counts = hii.host_counts([60.0, 41.0, 8.5, 8.6, 17.0])
    assert hii.host_summary(counts) == "O3V x1, O5V x1, O9V x1, B2V x2"


@pytest.mark.parametrize("letter", ["C", "D", "E", "G"])
def test_an_ionized_nebula_is_sized_by_its_stars(letter):
    draw.set_run_seed(5)
    for _ in range(8):
        nebula = Nebula(SystemConfig(), nebula_class=letter)
        data = tuning.NEBULA_CLASSES[letter]
        assert nebula.hosts and nebula.age_myr <= tuning.NEBULA_CLASS_MAX_AGE_MYR[letter]
        assert nebula.q_total and nebula.q_total > 0.0
        assert data["radius_range_ly"][0] <= nebula.radius_ly <= data["radius_range_ly"][1]
        assert data["density_range_cm3"][0] <= nebula.density_cm3 <= data["density_range_cm3"][1]


def test_a_given_host_and_age_set_the_nebula():
    draw.set_run_seed(6)
    nebula = Nebula(SystemConfig(), nebula_class="C", host=("O", 7), age_myr=0.1)
    assert nebula.age_myr == 0.1
    assert nebula.hosts == "O7V x1"
    assert math.log10(nebula.q_total) == pytest.approx(hii.log_q_of_type("O", 7))


def test_a_reflection_nebula_has_hosts_and_no_photons():
    draw.set_run_seed(7)
    nebula = Nebula(SystemConfig(), nebula_class="F")
    assert nebula.hosts and nebula.q_total is None


def test_clusters_add_stars_to_classes_d_and_e_only():
    draw.set_run_seed(8)
    cluster = Nebula(SystemConfig(), nebula_class="D", host=("O", 5), age_myr=3.0)
    alone = Nebula(SystemConfig(), nebula_class="C", host=("O", 5), age_myr=3.0)
    assert cluster.hosts.count(",") > 0 and "," not in alone.hosts
    assert cluster.q_total > alone.q_total


def test_a_dark_nebula_has_no_hosts():
    draw.set_run_seed(9)
    nebula = Nebula(SystemConfig(), nebula_class="M")
    assert nebula.hosts is None and nebula.age_myr is None and nebula.q_total is None


def test_the_hosted_rules_respect_the_age_table():
    assert set(tuning.NEBULA_CLASS_MAX_AGE_MYR) >= {c for rule in tuning.NEBULA_HOST_RULES for c in rule[3]}
    assert all(0.0 < rule[4] <= 0.5 for rule in tuning.NEBULA_HOST_RULES)
