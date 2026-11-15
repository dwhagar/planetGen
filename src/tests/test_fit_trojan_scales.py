# tests/test_fit_trojan_scales.py

"""
The fit of the trojan swarm scales TJ9 to TJ11 (GEN.1643, scripts/fit_trojan_scales.py): the Rayleigh equivalent of a distribution, the libration integral of the Jupiter paper, and the published numbers
the constants of `trojan_swarm` and `trojans` rest on (GEN.1647: the tadpole limit is an amplitude, 78 degrees, not the 23.9 degrees of the nearest approach).
"""

import importlib.util
import math
import os

import pytest

from planetgen.physics import trojan_swarm, trojans

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
_spec = importlib.util.spec_from_file_location("fit_trojan_scales", os.path.join(ROOT, "scripts", "fit_trojan_scales.py"))
fit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fit)


def test_a_rayleigh_law_is_its_own_equivalent():
    mean, median = fit.moments(fit.generalised(1.0, 0.05, 2.0), 0.0, 0.5)
    assert mean / fit.RAYLEIGH_MEAN == pytest.approx(0.05, rel=1e-3)
    assert median / fit.RAYLEIGH_MEDIAN == pytest.approx(0.05, rel=1e-3)


def test_the_libration_integral_has_its_minimum_at_the_point_and_its_edge_at_the_separatrix():
    assert fit.libration_function(60.0) == pytest.approx(fit.C_MIN)
    assert fit.libration_function(180.0) == pytest.approx(fit.C_MAX)
    assert fit.libration_amplitude_deg(fit.C_MIN) == 0.0
    amplitudes = [fit.libration_amplitude_deg(1.5 + 0.05 * k) for k in range(1, 21)]
    assert amplitudes == sorted(amplitudes)


def test_the_generator_keeps_the_nearest_approach_apart_from_the_amplitude_limit():
    # 23.9 degrees is where the tadpole ends nearest the secondary; the amplitude (half the range of the angle) at the same orbit is 78 degrees
    assert fit.closest_approach_deg() == pytest.approx(trojans.TADPOLE_NEAREST_APPROACH_DEG, abs=0.05)
    assert fit.separatrix_amplitude_deg() == pytest.approx(trojans.TADPOLE_LIMIT_DEG, abs=0.05)


def test_the_closed_form_mean_of_a_cut_rayleigh_is_what_the_generator_draws():
    n = 4000
    for sigma, limit in ((10.0, 23.9), (13.0, 78.0)):
        drawn = sum(trojan_swarm.truncated_rayleigh((i + 0.5) / n, sigma, limit) for i in range(n)) / n
        assert fit.truncated_rayleigh_mean(sigma, limit) == pytest.approx(drawn, rel=2e-3)


def test_the_jupiter_trojans_reach_past_the_nearest_approach_but_not_the_cut():
    shares = [fit.share_beyond(params, trojans.TADPOLE_NEAREST_APPROACH_DEG) for params in fit.JUPITER.values()]
    assert all(0.1 < s < 0.3 for s in shares)  # about a fifth of the real tadpoles librate wider than 23.9 degrees
    assert all(fit.share_beyond(params, trojans.TADPOLE_LIMIT_DEG) < 1e-6 for params in fit.JUPITER.values())  # none past the edge of the zone


def test_the_constants_sit_on_the_fit():
    got = fit.fit(fit.JUPITER)
    both = got["L4+L5"]
    assert both["libration"]["eq_mean"] == pytest.approx(trojan_swarm.LIBRATION_SIGMA_DEG, abs=0.5)  # TJ9
    assert both["ecc"]["eq_mean"] == pytest.approx(trojan_swarm.ECC_SIGMA, abs=0.005)  # TJ10
    assert both["inc"]["eq_mean"] == pytest.approx(trojan_swarm.INCLINATION_SIGMA_DEG, abs=0.5)  # TJ11
    drawn = fit.truncated_rayleigh_mean(trojan_swarm.LIBRATION_SIGMA_DEG, trojans.TADPOLE_LIMIT_DEG)
    assert drawn == pytest.approx(both["libration"]["mean"], abs=0.6)  # the cut draw has the Jupiter mean amplitude
    for cloud in ("L4", "L5"):
        for key in ("libration", "ecc", "inc"):
            row = got[cloud][key]
            assert abs(row["eq_mean"] - row["eq_median"]) < 0.12 * row["eq_mean"]  # the shapes are close to Rayleigh


def test_neptune_agrees_with_jupiter_on_the_libration_and_the_eccentricity():
    both = fit.fit(fit.JUPITER)["L4+L5"]
    assert fit.NEPTUNE["libration_table4"]["libration"] <= both["libration"]["eq_mean"] <= fit.NEPTUNE["single"]["libration"]
    assert both["ecc"]["eq_mean"] == pytest.approx(fit.NEPTUNE["single"]["ecc"], rel=0.35)
    assert fit.NEPTUNE["ranges"]["inc"][0] <= fit.NEPTUNE["single"]["inc"] <= fit.NEPTUNE["ranges"]["inc"][1]


def test_the_swarm_scales_are_the_ones_the_registry_names():
    assert (trojan_swarm.LIBRATION_SIGMA_DEG, trojan_swarm.ECC_SIGMA, trojan_swarm.INCLINATION_SIGMA_DEG) == (13.0, 0.06, 12.0)


def test_the_report_prints_every_cloud_and_the_weighted_pair():
    text = "\n".join(fit.report())
    for needle in ("L4+L5", "L5", "Neptune", "tadpole zone"):
        assert needle in text
