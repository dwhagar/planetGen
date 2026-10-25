"""
The stand-in registry of the atmosphere mixture draw (GEN.675, docs/design/research/gap-filling-plan.md section 6.7): every row
names a constant that exists and a test that exists, and every stand-in constant of the draw has a row.
"""

import importlib
import pathlib
import re

import pytest

PLAN = pathlib.Path(__file__).resolve().parents[2] / "docs" / "design" / "research" / "gap-filling-plan.md"
TESTS = pathlib.Path(__file__).resolve().parent

STANDINS = (
    "mixture_evolved.WATER_RESIDUAL", "mixture_evolved.SULFUR_SINK", "mixture_evolved.WATER_PHOTOLYSED",
    "mixture_evolved.OCEAN_WATER_WT_PERCENT", "mixture_evolved.CARBON_PER_OCEAN_WATER", "mixture_evolved.OXYGEN_SINK",
    "mixture_evolved.TITAN_HUMIDITY", "mixture_evolved.TITAN_HYDROGEN", "mixture_evolved.FROST_ICE",
    "mixture_lava.ALKALI_ACTIVITY_COEFFICIENT", "mixture_lava.COMPOSITION_SCATTER", "mixture_coma.DOMINANCE",
    "mixture_rocky.INVENTORY", "mixture_rocky.MAGMATIC_TEMPERATURE_K", "mixture_rocky.QUENCH_TEMPERATURE_K",
    "mixture_rocky.OXIDATION_TIME_GY", "mixture_trace.ABUNDANCES", "mixture_escape.EXOBASE_RATIO",
    "mixture_escape.XUV_EXPONENT", "mixture_escape.ESCAPE_WIDTH", "mixture_escape.LOST_BELOW", "mixture_escape.LOST_TOTAL",
    "mixture_evolved.TITAN_MAX_METHANE", "mixture_oxygen.INERT_FLOOR", "mixture_subject.PHOTOSPHERE_BAR",
    "mixture_subject.DAYSIDE_OVER_EQUILIBRIUM",
    "surface_pressure.PROCESSING_ENERGY_J_PER_KG", "surface_pressure.DEGAS_EFFICIENCY_H", "surface_pressure.DEGAS_EFFICIENCY",
    "surface_pressure.WEATHERING_ALPHA", "surface_pressure.WEATHERING_T_E_K", "surface_pressure.OUTGASSING_EXPONENT",
    "surface_pressure.HEAT_FLUX_EARTH_W_M2",
)


def _rows():
    text = PLAN.read_text()
    block = text.split("<!-- atmosphere-standins:begin -->")[1].split("<!-- atmosphere-standins:end -->")[0]
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")] for line in block.strip().splitlines()]
    return [row for row in rows[2:] if len(row) == 10]


def test_the_registry_has_rows():
    assert len(_rows()) >= len(STANDINS)


@pytest.mark.parametrize("row", _rows(), ids=lambda row: row[0])
def test_a_row_names_a_constant_and_a_test_that_exist(row):
    fill_id, kind, tier, tweak, evidence, constant, knobs, validity, falsifier, replaced_by = row
    assert kind in ("D", "B", "T", "P")
    assert evidence and knobs and validity and replaced_by
    module, name = constant.split(".")
    assert hasattr(importlib.import_module(f"planetgen.physics.{module}"), name), constant
    if not falsifier.startswith("none"):
        file, test = falsifier.split("::")
        assert re.search(rf"^def {re.escape(test)}\b", (TESTS / f"{file}.py").read_text(), re.M), falsifier
    else:
        assert len(falsifier) > len("none: ")


def test_every_stand_in_constant_has_a_row():
    named = {row[5] for row in _rows()}
    assert set(STANDINS) <= named
    assert len({row[0] for row in _rows()}) == len(_rows())
