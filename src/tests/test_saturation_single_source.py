# tests/test_saturation_single_source.py

"""One vapour-pressure function (GEN.632): no module outside `physics/saturation.py` holds a vapour-pressure constant, and the two
callers that did (the nitrogen frost of an old rogue-planet surface and the hydrosphere's boiling point) read the saturation curves."""

import pathlib

import pytest

from planetgen.physics import hydrosphere, saturation

SOURCE = pathlib.Path(saturation.__file__).resolve().parents[1]
OLD_NAMES = ("nitrogen_vapor_pressure_pa", "NITROGEN_TRIPLE_POINT_K", "NITROGEN_TRIPLE_POINT_PA", "NITROGEN_SUBLIMATION_J_MOL",
             "WATER_VAPORIZATION_J_MOL", "WATER_BOILING_K", "WATER_CRITICAL_K")


def test_the_old_vapour_pressure_names_are_gone_from_the_source():
    for path in SOURCE.rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        for name in OLD_NAMES:
            assert name not in text, f"{name} came back in {path}"
    for path in SOURCE.rglob("*.py"):
        if path.name != "saturation.py":
            assert "WATER_TRIPLE_POINT_PA =" not in path.read_text(encoding="utf-8"), path


def test_nitrogen_frost_pressure_is_the_saturation_curve():
    assert 0.5 < saturation.pressure_pa("N2", 40.0) < 20      # Pluto's air: a trace of nitrogen over its frost
    assert saturation.pressure_pa("N2", 20.0) < 1e-6


def test_water_boils_where_the_saturation_curve_says():
    assert hydrosphere.water_boiling_k(101325.0) == pytest.approx(373.15, abs=0.1)
    for temperature in (300.0, 450.0, 600.0):
        pressure = saturation.pressure_pa("H2O", temperature)
        assert hydrosphere.water_boiling_k(pressure) == pytest.approx(temperature, abs=1e-4)
    assert hydrosphere.water_boiling_k(1e9) == pytest.approx(647.1, abs=0.1)
    assert hydrosphere.water_boiling_k(0.0) == hydrosphere.WATER_MELTING_K
