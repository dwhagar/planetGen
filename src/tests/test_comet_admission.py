# tests/test_comet_admission.py

"""
GEN.507: the comet admission selection model (`planetgen.physics.comet_admission`): the fitted model and its refusals, the chance a
comet has a published composition by group and perihelion distance, the heliocentric distance it is measured at (truncated at the
perihelion), the sampler for the tuning harness, the agreement with the database's counts, and the registration as the selection
model the tuning harness reads.
"""

import json
import math
import os
import statistics
import subprocess
import sys

import pytest

from planetgen.generation import tuning_harness as th
from planetgen.physics import comet_admission as ca
from planetgen.physics.elements import DataError
from planetgen.util import draw

COMET = {"period_class": "jupiter_family", "perihelion_distance_au": 1.3}


def written(tmp_path, mutate=None):
    with open(ca.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    if mutate:
        mutate(document)
    path = tmp_path / "model.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return str(path)


def refused(tmp_path, match, mutate):
    with pytest.raises(DataError, match=match):
        ca.load(written(tmp_path, mutate))


# --- the model -------------------------------------------------------------------------------------------------------------


def test_the_model_holds_nine_cells_and_two_distance_lines_fitted_to_the_database():
    assert len(ca.MODEL.cells) == 9 and {c.group for c in ca.MODEL.cells} == set(ca.GROUPS)
    assert ca.MODEL.known_comets == sum(c.comets for c in ca.MODEL.cells) > 2000
    assert ca.MODEL.measured_comets == sum(c.measured for c in ca.MODEL.cells) == 93
    assert set(ca.MODEL.distance) == {"parent", "radical"}
    for line in ca.MODEL.distance.values():
        assert 0.2 < line.slope < 1.0 and line.scatter > 0 and line.rows > 100 and line.comets > 30
    assert "[R3782]" in ca.MODEL.source and ca.MODEL.retrieved
    assert ca.validity() == {"perihelion_distance_au": (0.3, 8.0)}


def test_the_loader_refuses_a_broken_model(tmp_path):
    refused(tmp_path, "not a group and a count", lambda d: d["admission"][0].update(measured=10 ** 6))
    refused(tmp_path, "not a group and a count", lambda d: d["admission"][0].update(group="trojan"))
    refused(tmp_path, "no perihelion range", lambda d: d["admission"][0].update(q_max_au=0.1))
    refused(tmp_path, "not contiguous", lambda d: d["admission"][1].update(q_min_au=1.2))
    refused(tmp_path, "missing", lambda d: d.update(admission=[c for c in d["admission"] if c["group"] != "halley_type"]))
    refused(tmp_path, "no radical", lambda d: d["distance"].pop("radical"))
    refused(tmp_path, "not positive", lambda d: d["distance"]["parent"].update(scatter=0))
    refused(tmp_path, "not a finite number", lambda d: d["distance"]["parent"].update(slope="x"))
    with pytest.raises(DataError, match="cannot read"):
        ca.load(str(tmp_path / "missing.json"))
    with pytest.raises(DataError, match="none of"):
        ca.measurement_distance_au(COMET, "metals")
    with pytest.raises(DataError, match="none of"):
        ca.admission_probability({"period_class": "trojan", "perihelion_distance_au": 1.0})


# --- admission -------------------------------------------------------------------------------------------------------------


def test_the_probability_is_the_jeffreys_share_of_the_cell():
    cell = ca.cell_of(COMET)
    assert (cell.group, cell.q_min_au, cell.q_max_au) == ("jupiter_family", 1.0, 2.0)
    assert ca.admission_probability(COMET) == pytest.approx((cell.measured + 0.5) / (cell.comets + 1.0))
    assert ca.admission_probability(COMET) > ca.admission_probability(dict(COMET, perihelion_distance_au=4.0))


def test_a_comet_without_a_period_class_is_long_period_and_an_empty_cell_is_not_zero():
    parabolic = {"period_class": None, "perihelion_distance_au": 1.5}
    assert ca.group_of(parabolic) == "long_period"
    assert ca.admission_probability(parabolic) == ca.admission_probability(dict(parabolic, period_class="long_period"))
    empty = ca.cell_of({"period_class": "halley_type", "perihelion_distance_au": 5.0})
    assert empty.measured == 0 and ca.admission_probability({"period_class": "halley_type", "perihelion_distance_au": 5.0}) > 0


def test_a_perihelion_outside_the_fit_takes_the_nearest_bin():
    assert ca.cell_of(dict(COMET, perihelion_distance_au=0.05)).q_min_au == 0.3
    assert ca.cell_of(dict(COMET, perihelion_distance_au=30.0)).q_max_au == 8.0
    assert ca.cell_of(dict(COMET, perihelion_distance_au=1.0)).q_min_au == 1.0


def test_the_expected_number_measured_among_the_known_comets_is_the_database_count():
    items = []
    for cell in ca.MODEL.cells:
        items += [{"period_class": cell.group, "perihelion_distance_au": 0.5 * (cell.q_min_au + cell.q_max_au)}] * cell.comets
    assert ca.expected_measured(items) == pytest.approx(ca.MODEL.measured_comets, abs=0.06 * ca.MODEL.measured_comets)


# --- distance --------------------------------------------------------------------------------------------------------------


def test_the_distance_is_never_inside_the_perihelion_and_repeats_on_a_seed():
    for q in (0.4, 1.0, 2.5, 6.0):
        for i in range(200):
            with draw.bound(draw.Stream(f"ca:{q}:{i}")):
                assert ca.measurement_distance_au({"period_class": None, "perihelion_distance_au": q}, "radical") >= q * (1 - 1e-9)
    with draw.bound(draw.Stream("ca:a")):
        first = ca.measurement_distance_au(COMET)
    with draw.bound(draw.Stream("ca:a")):
        assert ca.measurement_distance_au(COMET) == first


def test_the_median_distance_follows_the_fitted_line_where_the_truncation_does_not_bite():
    q = 0.5
    line = ca.MODEL.distance["parent"]
    values = []
    for i in range(1500):
        with draw.bound(draw.Stream(f"ca:med{i}")):
            values.append(ca.measurement_distance_au({"period_class": "jupiter_family", "perihelion_distance_au": q}))
    expected = math.exp(line.intercept + line.slope * math.log(q))
    assert statistics.median(values) == pytest.approx(expected, rel=0.06)
    assert statistics.median(values) > q


def test_radicals_are_measured_closer_in_than_parents_at_large_perihelion_only_through_their_own_line():
    parent, radical = ca.MODEL.distance["parent"], ca.MODEL.distance["radical"]
    assert radical.slope < parent.slope and radical.scatter > parent.scatter


# --- the fit script --------------------------------------------------------------------------------------------------------


def test_the_fit_script_joins_a_saved_database_to_the_orbits_and_reproduces_the_file(tmp_path):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
    try:
        import fit_comet_admission as script
    finally:
        sys.path.pop(0)
    orbits = {"fields": ["pdes", "class", "q"],
              "data": [[f"{n}", "JFc", str(1.0 + 0.1 * n)] for n in range(1, 30)] + [["C/2020 F3", "COM", "0.44"], ["C/2019 Y4", "HYP", "5.0"]]}
    for n in range(1, 30):
        orbits["data"][n - 1][0] = f"{n}P"
    header = "type,number,desig,X,Rh(au)"
    rows = [header] + [f"P,{n},,{'CN' if n % 2 else 'H2O'},{1.2 + 0.15 * n}" for n in range(1, 25)]
    rows += [f"C,,C/2020 F3,H2O,0.8"]
    document = script.build("\n".join(rows) + "\n", json.dumps(orbits), "2026-01-01")
    jf = [c for c in document["admission"] if c["group"] == "jupiter_family"]
    assert sum(c["comets"] for c in jf) == 29 and sum(c["measured"] for c in jf) == 24
    assert document["distance"]["parent"]["rows"] == 13 and document["distance"]["radical"]["rows"] == 12
    assert document["retrieved"] == "2026-01-01" and json.loads(json.dumps(document)) == document


# --- the harness -----------------------------------------------------------------------------------------------------------


def target(**overrides):
    raw = {"id": "comet-mean", "row": "tab:r", "statistic": "mean heliocentric distance of the measured comets", "kind": "mean",
           "value": 2.0, "sigma": 0.3, "n": 223, "source": "[R3782] database", "label": "DATA",
           "selection_model": ca.MODEL_ID, "field": "rh_au", "validity": {"perihelion_distance_au": [0.3, 8.0]}}
    raw.update(overrides)
    return th.TargetTable({"version": 1, "targets": [raw]}).get("comet-mean")


def test_the_model_is_registered_under_its_id_and_imported_on_first_use():
    assert th.SELECTION_MODELS[ca.MODEL_ID] is ca.comet_admission
    code = ("from planetgen.generation import tuning_harness as t; "
            "assert t.selection_model('comet-admission') is not None; print('ok')")
    assert subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout.strip() == "ok"


def test_the_harness_keeps_a_comet_with_its_probability_and_gives_it_a_distance():
    t = target()
    base = {"period_class": "jupiter_family", "perihelion_distance_au": 1.4}
    items = th.simulate(t, ca.observed_sampler(lambda i: dict(base)), 3000)
    assert 0 < len(items) < 3000 and all(item["rh_au"] >= 1.4 * (1 - 1e-9) for item in items)
    assert len(items) / 3000 == pytest.approx(ca.admission_probability(base), abs=0.04)
    assert th.simulate(t, ca.observed_sampler(lambda i: None), 5) == []
