# tests/test_mass_radius_sample.py

"""
GEN.505: the mass-radius sample selection model (`planetgen.physics.mass_radius_sample`): the fitted error model and its refusals, K
and the transit depth of a planet, the chance a planet passes the cut (mass error under 25 percent, radius error under 8 percent),
one planet measured, the agreement of the quadrature with many measurements and with the archive's pass fractions, and the
registration as the selection model the tuning harness reads.
"""

import json
import math
import os
import statistics
import subprocess
import sys

import pytest

from planetgen.generation import tuning_harness as th
from planetgen.physics import mass_radius_sample as mr
from planetgen.physics.elements import DataError
from planetgen.util import draw

PLANET = {"mass_earth": 5.0, "radius_earth": 1.9, "period_d": 4.0, "host_mass_solar": 0.3, "host_radius_solar": 0.3}


def written(tmp_path, mutate=None):
    with open(mr.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    if mutate:
        mutate(document)
    path = tmp_path / "model.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return str(path)


def refused(tmp_path, match, mutate):
    with pytest.raises(DataError, match=match):
        mr.load(written(tmp_path, mutate))


# --- the model -------------------------------------------------------------------------------------------------------------


def test_the_model_holds_the_cut_and_two_fitted_subsamples_of_the_archive():
    assert mr.MODEL.mass_cut == 0.25 and mr.MODEL.radius_cut == 0.08
    assert set(mr.MODEL.subsamples) == {"small", "m-dwarf"}
    small, mdwarf = mr.MODEL.subsamples["small"], mr.MODEL.subsamples["m-dwarf"]
    assert small.n > 500 and mdwarf.n > 80 and mdwarf.n < small.n
    assert 0.4 < small.pass_fraction < 0.65 and 0.6 < mdwarf.pass_fraction < 0.9
    for sub in (small, mdwarf):
        assert 0.5 < sub.sigma_k.slope < 1.2 and 0.5 < sub.sigma_depth.slope < 1.2, "the error grows about with the signal"
        assert sub.sigma_k.scatter > 0 and sub.sigma_depth.scatter > 0 and -1 < sub.correlation < 1
    assert "pscomppars" in mr.MODEL.source and mr.MODEL.retrieved


def test_the_loader_refuses_a_broken_model(tmp_path):
    refused(tmp_path, "out of range", lambda d: d["cut"].update(mass_error_share=1.5))
    refused(tmp_path, "not positive", lambda d: d["subsamples"]["small"]["ln_sigma_k"].update(scatter=0))
    refused(tmp_path, "outside", lambda d: d["subsamples"]["small"].update(residual_correlation=1.0))
    refused(tmp_path, "not a sample", lambda d: d["subsamples"]["small"].update(passing=10 ** 6))
    refused(tmp_path, "no sub-samples", lambda d: d.update(subsamples={}))
    refused(tmp_path, "not a finite number", lambda d: d.update(k_per_earth_mass_ms="x"))
    with pytest.raises(DataError, match="cannot read"):
        mr.load(str(tmp_path / "missing.json"))
    with pytest.raises(DataError, match="none of"):
        mr.pass_probability(PLANET, "giants")


def test_k_and_the_transit_depth_of_a_planet_follow_the_stored_quantities():
    # Earth around the Sun: 0.0895 m/s and (R_earth/R_sun)^2 = 84 ppm
    earth = {"mass_earth": 1.0, "radius_earth": 1.0, "period_d": 365.25, "host_mass_solar": 1.0, "host_radius_solar": 1.0}
    assert mr.rv_amplitude_ms(earth) == pytest.approx(0.08946, rel=1e-3)
    assert mr.transit_depth_ppm(earth) == pytest.approx(83.9, rel=0.01)
    assert mr.rv_amplitude_ms(dict(earth, period_d=365.25 / 8)) == pytest.approx(2 * mr.rv_amplitude_ms(earth), rel=1e-9)
    assert mr.transit_depth_ppm(dict(earth, radius_earth=2.0)) == pytest.approx(4 * mr.transit_depth_ppm(earth), rel=1e-9)


def test_the_bivariate_quadrature_matches_the_closed_forms():
    assert mr.both_below(0.3, 0.7, 0.0) == pytest.approx(mr._phi(0.3) * mr._phi(0.7), abs=1e-4)
    assert mr.both_below(0.5, 40.0, 0.6) == pytest.approx(mr._phi(0.5), abs=1e-4)
    assert mr.both_below(-40.0, 0.5, 0.6) == 0.0
    # P(Z1 < 0, Z2 < 0) = 1/4 + asin(rho)/(2 pi)
    assert mr.both_below(0.0, 0.0, 0.5) == pytest.approx(0.25 + math.asin(0.5) / (2 * math.pi), abs=1e-3)
    with pytest.raises(DataError, match="inside"):
        mr.both_below(0.0, 0.0, 1.0)


# --- one planet ------------------------------------------------------------------------------------------------------------


def test_a_cooler_host_and_a_bigger_signal_are_more_likely_to_make_the_sample():
    assert mr.pass_probability(PLANET, "m-dwarf") > mr.pass_probability(PLANET, "small")
    small_planet = dict(PLANET, mass_earth=0.5, radius_earth=0.8)
    big_planet = dict(PLANET, mass_earth=20.0, radius_earth=3.5)
    assert mr.pass_probability(big_planet) > mr.pass_probability(PLANET) > mr.pass_probability(small_planet)
    for item in (PLANET, small_planet, big_planet):
        assert 0.0 < mr.pass_probability(item) < 1.0


def test_an_observation_is_four_draws_repeats_on_a_seed_and_flags_the_cut():
    with draw.bound(draw.Stream("mr:a")):
        first = mr.observe(PLANET, "m-dwarf")
    with draw.bound(draw.Stream("mr:a")):
        again = mr.observe(PLANET, "m-dwarf")
    assert first == again
    assert first.passes == (first.mass_error_share < 0.25 and first.radius_error_share < 0.08)
    assert first.density_g_cm3 > 0 and first.mass_earth > 0 and first.radius_earth > 0
    true_density = PLANET["mass_earth"] * mr.M_EARTH_G / (4.0 / 3.0 * math.pi * (PLANET["radius_earth"] * mr.R_EARTH_CM) ** 3)
    assert 0.3 * true_density < first.density_g_cm3 < 3.0 * true_density


def test_many_observations_pass_at_the_rate_the_quadrature_gives():
    n = 4000
    passed = 0
    for i in range(n):
        with draw.bound(draw.Stream(f"mr:mc{i}")):
            passed += mr.observe(PLANET, "small").passes
    assert passed / n == pytest.approx(mr.pass_probability(PLANET, "small"), abs=0.03)


# --- the archive -----------------------------------------------------------------------------------------------------------


def test_the_model_applied_to_the_archives_planets_recovers_the_pass_fraction_it_was_fitted_to():
    """The Done test: the cut applied to the Gaussian summary of the archive's planets gives the published pass fraction."""
    for name, tolerance in (("small", 0.06), ("m-dwarf", 0.06)):
        sub = mr.MODEL.subsamples[name]
        assert mr.population_pass_fraction(name) == pytest.approx(sub.pass_fraction, abs=tolerance), name


def test_the_fit_script_reads_a_saved_query_and_reproduces_the_file(tmp_path):
    header = "pl_name,hostname,pl_orbper,pl_rade,pl_radeerr1,pl_radeerr2,pl_bmasse,pl_bmasseerr1,pl_bmasseerr2,pl_bmassprov,st_teff,st_mass,st_rad,tran_flag"
    lines = [header]
    for i in range(40):
        mass, radius = 1.0 + 0.3 * i, 1.0 + 0.04 * i
        lines.append(f"p{i},h{i},{3.0 + i},{radius},{0.05 * radius * (1 + i % 5)},-{0.05 * radius * (1 + i % 5)},{mass},"
                     f"{0.1 * mass * (1 + i % 7)},-{0.1 * mass * (1 + i % 7)},Mass,{3500 + 50 * i},0.5,0.5,1")
    lines.append("skip,h,5,1,0.1,-0.1,1,0.1,-0.1,Msini,3500,0.5,0.5,0")
    csv_path = tmp_path / "query.csv"
    csv_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
    try:
        import fit_mass_radius_sample as script
    finally:
        sys.path.pop(0)
    document = script.build(csv_path.read_text(encoding="utf-8"), "2026-01-01")
    assert document["subsamples"]["small"]["n"] == 40 and 0 < document["subsamples"]["small"]["passing"] < 40
    assert 0 < document["subsamples"]["m-dwarf"]["n"] < 40
    assert document["retrieved"] == "2026-01-01" and json.loads(json.dumps(document)) == document


# --- the harness -----------------------------------------------------------------------------------------------------------


def target(**overrides):
    raw = {"id": "mr-mean", "row": "tab:r", "statistic": "mean measured density of the precisely measured planets", "kind": "mean",
           "value": 5.0, "sigma": 1.0, "n": 50, "source": "[R2774]", "label": "READ", "selection_model": mr.MODEL_ID,
           "field": "measured_density_g_cm3", "validity": {"radius_earth": [None, 4.0]}}
    raw.update(overrides)
    return th.TargetTable({"version": 1, "targets": [raw]}).get("mr-mean")


def test_the_model_is_registered_under_its_id_and_keeps_the_items_that_pass():
    assert th.SELECTION_MODELS[mr.MODEL_ID] is mr.precise_sample
    assert mr.precise_sample({"passes": True}) == 1.0 and mr.precise_sample({"passes": False}) == 0.0
    assert mr.validity("small") == {"radius_earth": (None, 4.0)}
    assert mr.validity("m-dwarf")["host_teff_k"] == (None, 4000.0)


def test_a_target_naming_the_model_imports_it_on_first_use():
    code = ("from planetgen.generation import tuning_harness as t; "
            "assert t.selection_model('mass-radius-sample') is not None; print('ok')")
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert result.stdout.strip() == "ok"


def test_the_harness_keeps_the_planets_that_pass_in_proportion_to_the_chance():
    t = target()
    sampler = mr.observed_sampler(lambda i: dict(PLANET), "small")
    items = th.simulate(t, sampler, 1500)
    assert len(items) / 1500 == pytest.approx(mr.pass_probability(PLANET, "small"), abs=0.04)
    assert all(item["passes"] and item["mass_error_share"] < 0.25 and item["radius_error_share"] < 0.08 for item in items)
    values = th.values_of(t, items)
    assert statistics.fmean(values) == pytest.approx(
        PLANET["mass_earth"] * mr.M_EARTH_G / (4.0 / 3.0 * math.pi * (PLANET["radius_earth"] * mr.R_EARTH_CM) ** 3), rel=0.15)
    assert th.simulate(t, mr.observed_sampler(lambda i: None), 5) == []
