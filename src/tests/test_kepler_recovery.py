# tests/test_kepler_recovery.py

"""
GEN.504: the Kepler transit-plus-recovery selection model (`planetgen.physics.kepler_recovery`): the fitted model and its refusals, the
orbit and transit geometry, the noise model, the recovery curve, the chance Kepler finds a planet, the agreement with the injections
the model was fitted to, the fit script on a synthetic injection set, and the registration as the selection model the tuning harness
reads.
"""

import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

from planetgen.generation import tuning_harness as th
from planetgen.physics import kepler_recovery as kr
from planetgen.physics.elements import DataError
from planetgen.util import draw

SUN_LIKE = {"period_d": 10.0, "radius_earth": 2.0, "host_radius_solar": 1.0, "host_mass_solar": 1.0, "kepmag": 14.0}
EARTH = {"period_d": 365.25, "radius_earth": 1.0, "host_radius_solar": 1.0, "host_mass_solar": 1.0, "kepmag": 12.0}


def written(tmp_path, mutate=None):
    with open(kr.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    if mutate:
        mutate(document)
    path = tmp_path / "model.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return str(path)


def refused(tmp_path, match, mutate):
    with pytest.raises(DataError, match=match):
        kr.load(written(tmp_path, mutate))


# --- the model -------------------------------------------------------------------------------------------------------------


def test_the_model_holds_the_recovery_curve_the_noise_lines_and_the_injections_it_was_fitted_to():
    m = kr.MODEL
    assert m.injections > 140000 and 40000 < m.recovered < 50000
    assert m.recovery.c0 < 1 and 15 < m.recovery.a < 50 and 2.5 < m.recovery.n0 < 3.2
    assert [line.duration_h for line in m.noise] == [1.5, 3.0, 6.0, 12.0, 15.0]
    assert all(line.slope_kepmag < 0 and line.slope_kepmag2 > 0 and line.scatter > 0.2 for line in m.noise)
    assert len(m.kepmag_quantiles) == 21 and m.kepmag_quantiles[0] < 6 and m.kepmag_quantiles[-1] > 17
    assert m.period_range_d[0] < 1 and m.period_range_d[1] > 400
    assert "DR25" in m.source and m.retrieved
    assert kr.validity()["period_d"] == m.period_range_d


def test_the_loader_refuses_a_broken_model(tmp_path):
    refused(tmp_path, "out of range", lambda d: d["recovery"].update(c0=1.5))
    refused(tmp_path, "not a finite number", lambda d: d["recovery"].update(a="x"))
    refused(tmp_path, "not positive", lambda d: d["limb"].update(c0=0))
    refused(tmp_path, "not in order", lambda d: d["noise"]["durations"].reverse())
    refused(tmp_path, "at least two", lambda d: d["noise"].update(durations=d["noise"]["durations"][:1]))
    refused(tmp_path, "not sorted", lambda d: d["noise"].update(kepmag_quantiles=[3.0, 2.0, 1.0]))
    refused(tmp_path, "must be positive", lambda d: d.update(span_days=0))
    refused(tmp_path, "validity ranges", lambda d: d.update(validity={}))
    with pytest.raises(DataError, match="cannot read"):
        kr.load(str(tmp_path / "missing.json"))
    with pytest.raises(DataError, match="no kepmag"):
        kr.cdpp({k: v for k, v in SUN_LIKE.items() if k != "kepmag"}, 3.0)


# --- geometry --------------------------------------------------------------------------------------------------------------


def test_the_gamma_cdf_matches_the_closed_forms_and_scipy():
    from scipy import special
    for x in (0.1, 1.0, 5.0):
        assert kr.gamma_cdf(x, 1.0) == pytest.approx(1.0 - math.exp(-x), rel=1e-10)
        assert kr.gamma_cdf(x, 0.5) == pytest.approx(math.erf(math.sqrt(x)), rel=1e-9)
    for x, a in ((5.0, 26.7), (8.0, 26.7), (20.0, 26.7), (50.0, 26.7), (140.0, 100.0)):
        assert kr.gamma_cdf(x, a) == pytest.approx(special.gammainc(a, x), rel=1e-9, abs=1e-300)
    assert kr.gamma_cdf(0.0, 3.0) == 0.0 and kr.gamma_cdf(-1.0, 3.0) == 0.0


def test_the_earth_has_the_semi_major_axis_and_the_duration_it_has():
    assert kr.scaled_semi_major_axis(EARTH) == pytest.approx(215.0, rel=0.01)
    assert kr.transit_duration_h(EARTH, 0.0) == pytest.approx(13.0, rel=0.02)
    assert kr.transit_duration_h(EARTH, 0.9) < kr.transit_duration_h(EARTH, 0.0)
    assert kr.transit_duration_h(EARTH, 1.2) == 0.0, "a planet whose chord misses the star has no transit"
    assert kr.radius_ratio(EARTH) == pytest.approx(0.00916, rel=0.01)
    assert kr.transit_count(EARTH) == pytest.approx(kr.MODEL.transit_count_factor * kr.MODEL.span_days / 365.25)
    assert kr.transit_duration_h(dict(EARTH, host_mass_solar=1e-6, period_d=0.5), 0.0) == 0.0, "an orbit inside the star"


def test_the_limb_darkening_factor_falls_with_the_impact_parameter():
    assert kr.depth_factor(0.0) > kr.depth_factor(0.5) > kr.depth_factor(0.9) > 0.5
    assert kr.depth_factor(0.0) == pytest.approx(kr.MODEL.limb[0])


# --- noise -----------------------------------------------------------------------------------------------------------------


def test_a_fainter_or_noisier_host_has_more_noise_and_a_longer_duration_less():
    assert kr.cdpp(dict(SUN_LIKE, kepmag=15.0), 6.0) > kr.cdpp(SUN_LIKE, 6.0) > kr.cdpp(dict(SUN_LIKE, kepmag=11.0), 6.0)
    assert kr.cdpp(dict(SUN_LIKE, cdpp_z=1.0), 6.0) > kr.cdpp(SUN_LIKE, 6.0) > kr.cdpp(dict(SUN_LIKE, cdpp_z=-1.0), 6.0)
    assert kr.cdpp(SUN_LIKE, 1.5) > kr.cdpp(SUN_LIKE, 6.0) > kr.cdpp(SUN_LIKE, 15.0)
    assert kr.cdpp(SUN_LIKE, 0.2) == kr.cdpp(SUN_LIKE, 1.5) and kr.cdpp(SUN_LIKE, 40.0) == kr.cdpp(SUN_LIKE, 15.0)
    assert 1e-5 < kr.cdpp(SUN_LIKE, 3.0) < 1e-3, "tens to hundreds of ppm for a 14th magnitude star"


def test_the_noise_draw_gives_a_magnitude_in_the_targets_range_and_repeats_on_a_seed():
    planet = {k: v for k, v in SUN_LIKE.items() if k != "kepmag"}
    with draw.bound(draw.Stream("kr:a")):
        first = kr.with_noise(planet)
    with draw.bound(draw.Stream("kr:a")):
        assert kr.with_noise(planet) == first
    assert kr.MODEL.kepmag_quantiles[0] <= first["kepmag"] <= kr.MODEL.kepmag_quantiles[-1] and "cdpp_z" in first
    given = dict(planet, kepmag=10.0, cdpp_z=0.5)
    assert kr.with_noise(given) == given, "a given magnitude and draw are kept"
    medians = []
    for i in range(500):
        with draw.bound(draw.Stream(f"kr:m{i}")):
            medians.append(kr.with_noise(planet)["kepmag"])
    assert 14.0 < sorted(medians)[250] < 15.2


# --- recovery --------------------------------------------------------------------------------------------------------------


def test_the_recovery_curve_rises_with_the_signal_and_needs_three_transits_and_a_short_duration():
    curve = [kr.recovery_probability(mes, 10.0, 6.0) for mes in (3, 6, 7.1, 9, 12, 30)]
    assert curve == sorted(curve) and curve[0] < 0.01 and 0.1 < curve[2] < 0.35 and curve[-1] == pytest.approx(kr.MODEL.recovery.c0, abs=0.03)
    assert kr.recovery_probability(30.0, 1.0, 6.0) < 1e-3 and kr.recovery_probability(30.0, 6.0, 6.0) > 0.85
    assert kr.recovery_probability(30.0, 10.0, 200.0) < 0.01
    assert kr.recovery_probability(0.0, 10.0, 6.0) == 0.0 and kr.recovery_probability(30.0, 10.0, 0.0) == 0.0


def test_the_curve_agrees_with_the_injection_calibration_it_was_fitted_to():
    """The Done test: the stored calibration (observed against predicted recovery by MES) agrees bin by bin."""
    for row in kr.MODEL.calibration:
        if row["n"] > 2000:
            assert row["predicted"] == pytest.approx(row["observed"], abs=0.02), row
    assert kr.MODEL.predicted_recovered == pytest.approx(kr.MODEL.recovered, rel=0.005)


def test_a_host_with_a_modelled_noise_level_recovers_the_injections_within_a_few_percent():
    check = kr.MODEL.generated_host
    assert check["rows"] > 140000
    assert check["predicted"] == pytest.approx(check["observed"], rel=0.08)


# --- one planet ------------------------------------------------------------------------------------------------------------


def test_the_chance_is_the_geometry_times_the_recovery_and_is_a_probability():
    p = kr.detection_probability(SUN_LIKE)
    assert 0 < p <= 1.0 / kr.scaled_semi_major_axis(SUN_LIKE)
    for b in (0.0, 0.5, 0.95):
        assert 0 <= kr.recovery_at(SUN_LIKE, b) <= kr.MODEL.recovery.c0
    assert kr.recovery_at(SUN_LIKE, 1.0) == 0.0 and kr.recovery_at(SUN_LIKE, -0.1) == 0.0, "outside the fit"
    assert kr.detection_probability(dict(SUN_LIKE, impact_parameter=0.3)) == kr.recovery_at(SUN_LIKE, 0.3)
    assert kr.detection_probability(dict(SUN_LIKE, impact_parameter=1.4)) == 0.0
    assert kr.detection_probability(dict(SUN_LIKE, host_mass_solar=1e-6, period_d=0.5)) == 0.0


def test_a_shorter_period_a_bigger_planet_and_a_brighter_star_are_easier_to_find():
    base = kr.detection_probability(SUN_LIKE)
    assert kr.detection_probability(dict(SUN_LIKE, period_d=3.0)) > base > kr.detection_probability(dict(SUN_LIKE, period_d=30.0))
    assert kr.detection_probability(dict(SUN_LIKE, radius_earth=4.0)) >= 0.99 * base > kr.detection_probability(dict(SUN_LIKE, radius_earth=1.0))
    assert kr.detection_probability(dict(SUN_LIKE, kepmag=11.0)) > base > kr.detection_probability(dict(SUN_LIKE, kepmag=16.0))
    assert kr.detection_probability(dict(SUN_LIKE, radius_earth=1.0, period_d=200.0)) < 0.002
    assert kr.detection_probability(dict(SUN_LIKE, cdpp_z=-1.0)) > base > kr.detection_probability(dict(SUN_LIKE, cdpp_z=1.0))


def test_the_expected_number_found_is_the_sum_of_the_chances():
    items = [dict(SUN_LIKE, period_d=p) for p in (2.0, 5.0, 12.0, 40.0)]
    assert kr.expected_found(items) == pytest.approx(sum(kr.detection_probability(i) for i in items))


# --- the fit script --------------------------------------------------------------------------------------------------------


def synthetic(n_stars=300, n_inj=4000):
    rng = np.random.default_rng(7)
    kp = rng.uniform(10.0, 16.0, n_stars)
    radius = rng.uniform(0.7, 1.5, n_stars)
    lines = ["kepid,kepmag,teff,radius,mass,dutycycle,dataspan,rrmscdpp01p5,rrmscdpp03p0,rrmscdpp06p0,rrmscdpp12p0,rrmscdpp15p0"]
    cdpp = []
    for i in range(n_stars):
        base = math.exp(8.5 - 0.9 * kp[i] + 0.05 * kp[i] ** 2) * radius[i] ** 0.2
        row = [base * f for f in (1.0, 0.75, 0.55, 0.42, 0.38)]
        cdpp.append(row)
        lines.append(f"{1000 + i},{kp[i]:.3f},5700,{radius[i]:.3f},1.0,0.88,1458.931," + ",".join(f"{v:.3f}" for v in row))
    names = ["KIC_ID", "i_period", "N_Transit", "i_ror", "i_b", "i_dur", "Expected_MES", "Recovered"]
    out = ["\\name = synthetic", "| " + " | ".join(names) + " |", "| " + " | ".join(["float"] * len(names)) + " |",
           "| " + " | ".join(["x"] * len(names)) + " |", "| " + " | ".join(["null"] * len(names)) + " |"]
    for _ in range(n_inj):
        s = int(rng.integers(n_stars))
        period = math.exp(rng.uniform(math.log(0.5), math.log(500.0)))
        ror = math.exp(rng.uniform(math.log(0.005), math.log(0.06)))
        b = rng.uniform(0.0, 1.0)
        dur = max(0.5, 13.0 * (period / 365.25) ** (1 / 3) * math.sqrt(max(1e-3, 1 - b * b)) * radius[s] ** 0.5)
        count = max(0.0, 0.96 * 1458.931 / period)
        c = cdpp[s][2] * 1e-6
        mes = (1.17 - 0.09 * b * b - 0.51 * b ** 4) * math.sqrt(count) * ror * ror / c
        p = 0.94 * (1 - math.exp(-((mes / 8.0) ** 4))) if count >= 3 else 0.0
        out.append(f"{1000 + s} {period:.5f} {count:.4f} {ror:.5f} {b:.4f} {dur:.4f} {mes:.4f} {int(rng.random() < p)}")
    return "\n".join(out) + "\n", "\n".join(lines) + "\n"


def script():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
    try:
        import fit_kepler_recovery
    finally:
        sys.path.pop(0)
    return fit_kepler_recovery


def test_the_fit_script_recovers_the_planted_limb_law_and_a_rising_recovery_curve():
    fit = script()
    injections, stars = synthetic()
    document = fit.build(injections, stars, "2026-01-01")
    assert document["injections"]["rows"] == 4000 and document["injections"]["used"] > 3000
    # the planted effective depth is 1.17 - 0.09 b^2 - 0.51 b^4 for the 6 h CDPP; the script reads the CDPP at the injection's duration
    assert 0.5 < document["limb"]["c0"] < 2.0 and document["limb"]["c4"] < 0.5
    assert document["transit_count_factor"] == pytest.approx(0.96 * 1458.931 / 1458.931, rel=0.02)
    assert 0 < document["recovery"]["c0"] <= 1 and document["recovery"]["a"] > 0
    assert document["noise"]["stars"] == 300 and len(document["noise"]["durations"]) == 5
    assert document["noise"]["durations"][0]["slope_kepmag"] == pytest.approx(-0.9, abs=0.4)
    low = [r for r in document["calibration"] if r["mes_max"] <= 5 and r["n"] > 50]
    assert all(r["observed"] < 0.1 for r in low)
    assert document["retrieved"] == "2026-01-01" and json.loads(json.dumps(document)) == document
    assert fit.close(document, json.loads(json.dumps(document))) and not fit.close(document, dict(document, span_days=1.0))


# --- the harness -----------------------------------------------------------------------------------------------------------


def target(**overrides):
    raw = {"id": "kepler-mean", "row": "tab:r", "statistic": "mean radius of the planets Kepler would have found", "kind": "mean",
           "value": 2.5, "sigma": 0.5, "n": 1000, "source": "Kepler DR25 injections", "label": "DATA",
           "selection_model": kr.MODEL_ID, "field": "radius_earth", "validity": {"period_d": [0.5, 500.0]}}
    raw.update(overrides)
    return th.TargetTable({"version": 1, "targets": [raw]}).get("kepler-mean")


def test_the_model_is_registered_under_its_id_and_imported_on_first_use():
    assert th.SELECTION_MODELS[kr.MODEL_ID] is kr.kepler_recovery
    code = ("from planetgen.generation import tuning_harness as t; "
            "assert t.selection_model('kepler-recovery') is not None; print('ok')")
    assert subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout.strip() == "ok"


def test_the_harness_keeps_a_planet_with_its_probability_and_skips_a_planet_the_row_did_not_draw():
    t = target()
    base = {k: v for k, v in SUN_LIKE.items() if k != "kepmag"}
    items = th.simulate(t, kr.observed_sampler(lambda i: dict(base, period_d=2.0 + 0.01 * i)), 3000)
    assert 0 < len(items) < 3000 and all("kepmag" in item and "cdpp_z" in item for item in items)
    assert th.simulate(t, kr.observed_sampler(lambda i: None), 5) == []
    total = 0.0
    for i in range(3000):
        with draw.bound(draw.Stream(f"tuning:{t.id}:{i}")):
            total += kr.detection_probability(kr.with_noise(dict(base, period_d=2.0 + 0.01 * i)))
    assert len(items) / 3000 == pytest.approx(total / 3000, abs=0.02)
