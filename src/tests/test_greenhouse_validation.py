# tests/test_greenhouse_validation.py

"""GEN.1854: the greenhouse table is measured live on the validation airs kept in the repository (`tests/fixtures/greenhouse_validation.json.gz`: airs of the generator the radiative-convective model solved, seeds 5000 to 5060 of 11-26,
and the hydrogen envelope and giant airs), so that its accuracy on today's airs cannot drift unseen again: the table of GEN.1796 had, and read 1.75 percent (clear cool airs) and 2.5 percent (all) on them while its summary said 1 and 2.
The summary test of `test_greenhouse_table.py` reads what the file states; this one reads what the table does."""

import gzip
import json
import os

import numpy as np
import pytest

from planetgen.physics import greenhouse_table as gt

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "greenhouse_validation.json.gz")
HYDROGEN_RICH_SHARE = 0.1
"""float: The hydrogen share from which an air is stated apart (`fit_greenhouse_table.HYDROGEN_RICH_SHARE`)."""
COOL_K = 300.0
CLEAR_MEDIAN = 0.0105
"""float: The median relative error the clear cool airs may reach (the file before the refit of GEN.1854 read 0.0175 on them)."""


@pytest.fixture(scope="module")
def errors():
    with gzip.open(FIXTURE, "rt", encoding="utf-8") as handle:
        data = json.load(handle)

    def read(rows):
        model = np.array([row["model_k"] for row in rows])
        table = np.array([gt.warming(row["partials"], row["gravity"], row["flux"], row["aloft_share"], row["aloft_pressure_pa"]).surface_temperature_k for row in rows])
        return np.abs(table / model - 1.0), np.abs(table - model), model

    rows = [row for row in data["validation"] if row["partials"].get("H2", 0.0) / sum(row["partials"].values()) < HYDROGEN_RICH_SHARE]
    relative, kelvin, model = read(rows)
    classes = np.array([row["class"] for row in rows])
    envelope_relative, _kelvin, _model = read(data["envelope"])
    return {"relative": relative, "kelvin": kelvin, "model": model, "classes": classes, "envelope": envelope_relative}


def test_the_table_reads_the_airs_of_the_generator_within_what_the_file_states(errors):
    assert len(errors["relative"]) >= 1500 and len(errors["envelope"]) >= 300
    cool = errors["model"] < COOL_K
    assert cool.sum() >= 500
    assert np.percentile(errors["relative"][cool], 50) < 0.02
    assert np.percentile(errors["relative"], 50) < 0.02
    clear = cool & (errors["classes"] == "clear")
    assert clear.sum() >= 200
    assert np.percentile(errors["relative"][clear], 50) < CLEAR_MEDIAN
    assert np.percentile(errors["relative"][clear], 90) < 0.10
    assert np.percentile(errors["kelvin"][clear], 99) < 50.0
    for name in ("haze", "other_infrared_gases"):
        chosen = cool & (errors["classes"] == name)
        assert np.percentile(errors["relative"][chosen], 50) < 0.06, name
        assert np.percentile(errors["relative"][chosen], 90) < 0.30, name


def test_the_table_reads_the_hydrogen_envelopes_and_giants_within_what_the_file_states(errors):
    assert np.percentile(errors["envelope"], 50) < 0.05
    assert np.percentile(errors["envelope"], 90) < 0.20
    assert np.percentile(errors["envelope"], 99) < 0.40
