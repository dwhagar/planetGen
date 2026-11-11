"""Each physics module that other modules reach through a long chain imports on its own in a fresh interpreter (GEN.1401: `hydrosphere` and `bulk_draw` failed with a circular import)."""

import subprocess
import sys

import pytest

MODULES = ("hydrosphere", "bulk_draw", "icy_tide", "planet_type", "orbit_elements", "interior_heat", "envelope_capture", "surface_composition", "aquifers")


@pytest.mark.parametrize("module", MODULES)
def test_a_physics_module_imports_in_a_fresh_interpreter(module):
    done = subprocess.run([sys.executable, "-c", f"import planetgen.physics.{module}"], capture_output=True, text=True, timeout=120)
    assert done.returncode == 0, done.stderr[-600:]


def test_the_melting_point_of_water_has_one_definition():
    from planetgen.physics import aquifers, hydrosphere, icy_tide
    assert hydrosphere.WATER_MELTING_K is aquifers.WATER_MELTING_K
    assert icy_tide.flux_from_tide.__defaults__[-1] == aquifers.WATER_MELTING_K
