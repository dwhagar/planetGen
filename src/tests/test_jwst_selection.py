# tests/test_jwst_selection.py

"""
GEN.506: the JWST atmosphere target-selection model (`planetgen.physics.jwst_selection`): the fitted model and its refusals, the
transmission and emission metrics of [R2315], the chance a planet has a published JWST spectrum, host magnitudes for a
generated planet, the fit script recovering a planted selection, and the registration as the selection model the harness reads.
"""

import json
import math
import os
import random
import statistics
import subprocess
import sys

import pytest

from planetgen.generation import tuning_harness as th
from planetgen.physics import jwst_selection as js
from planetgen.physics.elements import DataError
from planetgen.util import draw

GIANT = {"radius_earth": 14.0, "mass_earth": 89.0, "teq_k": 1100.0, "host_radius_solar": 0.9, "host_teff_k": 5400.0, "j_mag": 10.2,
         "k_mag": 9.8}
SUB_NEPTUNE = {"radius_earth": 2.74, "mass_earth": 8.17, "teq_k": 596.0, "host_radius_solar": 0.216, "host_teff_k": 3250.0, "j_mag": 9.75,
               "k_mag": 8.78}
EARTH = {"radius_earth": 1.0, "mass_earth": 1.0, "teq_k": 255.0, "host_radius_solar": 1.0, "host_teff_k": 5772.0, "j_mag": 12.0, "k_mag": 11.8}


def script():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
    try:
        import fit_jwst_selection
    finally:
        sys.path.pop(0)
    return fit_jwst_selection


def written(tmp_path, mutate=None):
    with open(js.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    if mutate:
        mutate(document)
    path = tmp_path / "model.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return str(path)


def refused(tmp_path, match, mutate):
    with pytest.raises(DataError, match=match):
        js.load(written(tmp_path, mutate))


# --- the model -------------------------------------------------------------------------------------------------------------


def test_the_model_is_the_fit_to_the_archives_transiting_planets_and_reproduces_the_count_in_every_radius_class():
    model = js.MODEL
    assert model.planets > 4000 and 50 < model.observed < 150 and model.auc > 0.9
    assert model.tsm_slope > 0 and model.esm_slope > 0, "brighter signals are chosen more often"
    assert model.class_edges[-1] == math.inf and len(model.class_offsets) == len(model.class_edges) == 4
    assert model.class_offsets[0] == 0.0 and all(o < 0 for o in model.class_offsets[1:])
    for entry in model.calibration:
        assert entry["predicted"] == pytest.approx(entry["observed"], abs=0.05)
    assert sum(e["observed"] for e in model.calibration) == model.observed
    assert sum(e["planets"] for e in model.calibration) == model.planets
    assert "pscomppars" in model.source and model.retrieved


def test_the_loader_refuses_a_broken_model(tmp_path):
    refused(tmp_path, "do not match", lambda d: d["coefficients"].update(class_offsets=[0.0]))
    refused(tmp_path, "not contiguous", lambda d: d["radius_classes_earth"][1].__setitem__(0, 2.0))
    refused(tmp_path, "out of range", lambda d: d["host_magnitudes"].update(sd_j=0))
    refused(tmp_path, "not a sample", lambda d: d.update(observed=0))
    refused(tmp_path, "does not cover", lambda d: d.update(calibration=d["calibration"][:2]))
    refused(tmp_path, "not a finite number", lambda d: d["coefficients"].update(log10_tsm="steep"))
    with pytest.raises(DataError, match="cannot read"):
        js.load(str(tmp_path / "missing.json"))


# --- the metrics -----------------------------------------------------------------------------------------------------------


def test_the_transmission_metric_follows_equation_1_with_the_scale_factor_of_the_radius_class():
    by_hand = 1.26 * 2.74 ** 3 * 596.0 / (8.17 * 0.216 ** 2) * 10 ** (-9.75 / 5)
    assert js.tsm(SUB_NEPTUNE) == pytest.approx(by_hand, rel=1e-12)
    below = js.tsm(dict(SUB_NEPTUNE, radius_earth=1.4, mass_earth=3.0))
    above = js.tsm(dict(SUB_NEPTUNE, radius_earth=1.6, mass_earth=3.0))
    assert below / above == pytest.approx((0.190 / 1.26) * (1.4 / 1.6) ** 3, rel=1e-12), "the scale factor changes at 1.5 Earth radii"
    assert js.tsm(dict(SUB_NEPTUNE, j_mag=SUB_NEPTUNE["j_mag"] + 5)) == pytest.approx(js.tsm(SUB_NEPTUNE) / 10.0)


def test_the_emission_metric_follows_equation_4_and_scales_with_the_dayside_blackbody():
    # GJ 1132 b, the benchmark of [R2315], has an ESM of 7.5 by construction (their section 3.2)
    gj1132b = {"radius_earth": 1.13, "teq_k": 529.0, "host_radius_solar": 0.207, "host_teff_k": 3270.0, "k_mag": 8.32}
    assert js.esm(gj1132b) == pytest.approx(7.5, rel=0.25)
    assert js.esm(dict(gj1132b, k_mag=gj1132b["k_mag"] + 5)) == pytest.approx(js.esm(gj1132b) / 10.0)
    assert js.esm(dict(gj1132b, teq_k=700.0)) > js.esm(gj1132b)


# --- one planet ------------------------------------------------------------------------------------------------------------


def test_a_bright_hot_large_scale_height_planet_is_more_likely_to_have_a_spectrum_than_a_cold_faint_one():
    assert js.probability(SUB_NEPTUNE) > 0.2 and js.probability(GIANT) > 0.1
    assert js.probability(EARTH) < 1e-3
    fainter = dict(SUB_NEPTUNE, j_mag=SUB_NEPTUNE["j_mag"] + 2, k_mag=SUB_NEPTUNE["k_mag"] + 2)
    assert js.probability(fainter) < js.probability(SUB_NEPTUNE)
    cooler = dict(SUB_NEPTUNE, teq_k=300.0)
    assert js.probability(cooler) < js.probability(SUB_NEPTUNE)
    for item in (GIANT, SUB_NEPTUNE, EARTH):
        assert 0.0 < js.probability(item) < 1.0 and js.jwst_target(item) == js.probability(item)


def test_the_radius_classes_are_those_of_the_fit():
    assert [js.radius_class(r) for r in (0.5, 1.49, 1.5, 3.9, 4.0, 9.9, 10.0, 20.0)] == [0, 0, 1, 1, 2, 2, 3, 3]


# --- host magnitudes -------------------------------------------------------------------------------------------------------


def test_a_generated_planet_gets_a_host_magnitude_from_what_the_archives_hosts_show_and_repeats_on_a_seed():
    base = {k: v for k, v in SUB_NEPTUNE.items() if k not in ("j_mag", "k_mag")}
    with draw.bound(draw.Stream("jwst:a")):
        first = js.with_magnitudes(base)
    with draw.bound(draw.Stream("jwst:a")):
        again = js.with_magnitudes(base)
    assert first == again and "j_mag" in first and first["radius_earth"] == base["radius_earth"]
    js_values, colours = [], []
    for i in range(1500):
        with draw.bound(draw.Stream(f"jwst:m{i}")):
            item = js.with_magnitudes(base)
        js_values.append(item["j_mag"])
        colours.append(item["j_mag"] - item["k_mag"])
    assert min(js_values) >= js.MODEL.min_j and max(js_values) <= js.MODEL.max_j
    assert statistics.fmean(js_values) == pytest.approx(js.MODEL.mean_j, abs=0.25)
    expected = js.MODEL.colour_intercept + js.MODEL.colour_slope * math.log10(base["host_teff_k"])
    assert statistics.fmean(colours) == pytest.approx(expected, abs=0.01)
    assert statistics.stdev(colours) == pytest.approx(js.MODEL.colour_scatter, rel=0.15)


# --- the fit ---------------------------------------------------------------------------------------------------------------


def planted(seed=11, count=4000):
    """A saved archive query of a planet population with a known selection: logit P = -6 + 2 log10 TSM + 0.5 log10 ESM + class offset."""
    fit = script()
    rng = random.Random(seed)
    planets = ["pl_name,pl_rade,pl_bmasse,pl_eqt,st_teff,st_rad,sy_jmag,sy_kmag,tran_flag"]
    spectra = ["pl_name,spec_type,facility"]
    offsets = (0.0, -1.0, -1.5, -2.0)
    for i in range(count):
        radius = math.exp(rng.uniform(math.log(0.6), math.log(18.0)))
        teff = rng.uniform(3000.0, 6500.0)
        host_radius = 0.2 + 0.25 * (teff - 3000.0) / 1000.0
        teq = rng.uniform(300.0, 1800.0)
        mass = 1.436 * radius ** 1.7
        j_mag = rng.uniform(7.0, 14.0)
        k_mag = j_mag - 0.9
        item = fit.tsm(radius, mass, teq, host_radius, j_mag), fit.esm(radius, teq, teff, host_radius, k_mag)
        klass = next(i2 for i2, (low, high) in enumerate(fit.RADIUS_BINS) if low <= radius < high)
        eta = -2.0 + 2.0 * math.log10(item[0]) + 0.5 * math.log10(item[1]) + offsets[klass] - 3.0
        planets.append(f"p{i},{radius},{mass},{teq},{teff},{host_radius},{j_mag},{k_mag},1")
        if rng.random() < 1.0 / (1.0 + math.exp(-eta)):
            spectra.append(f"p{i},Transmission,NASA 6.5m James Webb Space Telescope (JWST) Satellite Mission")
        else:
            spectra.append(f"p{i},Transmission,Hubble Space Telescope satellite")
    return "\n".join(planets) + "\n", "\n".join(spectra) + "\n", offsets


def test_the_fit_script_recovers_a_planted_selection():
    """The Done test: the fit applied to a planet population with a known selection returns the planted coefficients."""
    planets, spectra, offsets = planted()
    document = script().build(planets, spectra, "2026-01-01")
    coefficients = document["coefficients"]
    assert 100 < document["observed"] < document["planets"] / 2
    assert coefficients["log10_tsm"] == pytest.approx(2.0, abs=0.35)
    assert coefficients["log10_esm"] == pytest.approx(0.5, abs=0.35)
    for fitted, planted_offset in zip(coefficients["class_offsets"][1:], offsets[1:]):
        assert fitted == pytest.approx(planted_offset, abs=0.6)
    for entry in document["calibration"]:
        assert entry["predicted"] == pytest.approx(entry["observed"], abs=0.05)
    assert 0.5 < document["auc_logistic"] < 1.0
    assert json.loads(json.dumps(document)) == document


def test_the_fit_script_ignores_non_transiting_planets_incomplete_rows_and_direct_images():
    fit = script()
    planets = ("pl_name,pl_rade,pl_bmasse,pl_eqt,st_teff,st_rad,sy_jmag,sy_kmag,tran_flag\n"
               "a,2,5,600,3500,0.3,10,9,1\nb,2,5,600,3500,0.3,10,9,0\nc,2,,600,3500,0.3,10,9,1\n")
    spectra = ("pl_name,spec_type,facility\na,Direct Imaging,NASA 6.5m James Webb Space Telescope (JWST) Satellite Mission\n"
               "a,Transmission,NASA 6.5m James Webb Space Telescope (JWST) Satellite Mission\n")
    rows = fit.sample(planets, spectra)
    assert len(rows) == 1 and rows[0][0] is True
    rows = fit.sample(planets, "pl_name,spec_type,facility\na,Direct Imaging,NASA 6.5m James Webb Space Telescope (JWST) Satellite Mission\n")
    assert rows[0][0] is False


# --- the harness -----------------------------------------------------------------------------------------------------------


def target(**overrides):
    raw = {"id": "jwst-mean", "row": "tab:r", "statistic": "mean radius of the planets with a JWST spectrum", "kind": "mean",
           "value": 8.0, "sigma": 2.0, "n": 78, "source": "NASA Exoplanet Archive spectra table", "label": "DATA",
           "selection_model": js.MODEL_ID, "field": "radius_earth", "validity": {"j_mag": [3.9, 16.3]}}
    raw.update(overrides)
    return th.TargetTable({"version": 1, "targets": [raw]}).get("jwst-mean")


def test_the_model_is_registered_under_its_id_and_imported_on_first_use():
    assert th.SELECTION_MODELS[js.MODEL_ID] is js.jwst_target
    code = ("from planetgen.generation import tuning_harness as t; "
            "assert t.selection_model('jwst-target-selection') is not None; print('ok')")
    assert subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout.strip() == "ok"


def test_the_harness_keeps_a_planet_with_its_probability_and_skips_a_planet_the_row_did_not_draw():
    t = target()
    base = {k: v for k, v in SUB_NEPTUNE.items() if k not in ("j_mag", "k_mag")}
    items = th.simulate(t, js.observed_sampler(lambda i: dict(base)), 1500)
    assert 0 < len(items) < 1500 and all("j_mag" in item and "k_mag" in item for item in items)
    assert th.simulate(t, js.observed_sampler(lambda i: None), 5) == []
    # the same streams give the same magnitudes, so the expected share is the mean probability over them
    total = 0.0
    for i in range(1500):
        with draw.bound(draw.Stream(f"tuning:{t.id}:{i}")):
            total += js.probability(js.with_magnitudes(base))
    assert len(items) / 1500 == pytest.approx(total / 1500, abs=0.02)
