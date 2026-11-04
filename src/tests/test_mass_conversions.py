# tests/test_mass_conversions.py

"""GEN.934: a mass is converted with the unit registry (`unit_registry.converter`), not with a hand multiply by a mass constant."""

import ast
import pathlib

import pytest

import planetgen
from planetgen.physics import constants, unit_registry

ROOT = pathlib.Path(planetgen.__file__).parent
MASS_CONSTANTS = {"SOLAR_MASS_TO_KG", "EARTH_MASS_TO_KG", "JUPITER_MASS_TO_KG"}
SKIPPED = {pathlib.PurePath("physics/constants.py"), pathlib.PurePath("physics/unit_registry.py")}

REMAINING = {
    "generation/multistar.py": 1,
    "generation/open_clusters.py": 1,
    "generation/place_terms.py": 1,
    "physics/accretion.py": 5,
    "physics/body_status.py": 1,
    "physics/mathcheck.py": 1,
    "physics/planet_survival.py": 1,
    "physics/planets.py": 1,
    "physics/rogue_surface.py": 1,
    "physics/rogue_volatiles.py": 1,
    "physics/venus_spin.py": 1,
}
"""dict: The multiplies and divides by a mass constant that stay by hand, by file. Each sits inside a longer product (`G * m * C`, `3 * m * C`), where a converter call would
regroup the multiplication and change the last bit of a generated number, or has the constant on the left of a divide (`C / m`), which is no unit conversion. A new site is a
failure (use `unit_registry.converter`); a site that is converted or removed must come off this list."""


def _sites():
    found = {}
    for path in sorted(ROOT.rglob("*.py")):
        relative = path.relative_to(ROOT)
        if relative in SKIPPED:
            continue
        count = 0
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Mult, ast.Div)):
                count += sum(1 for operand in (node.left, node.right) if isinstance(operand, ast.Attribute) and operand.attr in MASS_CONSTANTS)
        if count:
            found[relative.as_posix()] = count
    return found


def test_no_new_hand_written_mass_conversion():
    found = _sites()
    new = {k: v for k, v in found.items() if v > REMAINING.get(k, 0)}
    assert not new, f"use unit_registry.converter instead of a multiply or divide by a mass constant: {new}"


def test_the_list_of_kept_sites_has_no_stale_entry():
    found = _sites()
    stale = {k: v for k, v in REMAINING.items() if found.get(k, 0) < v}
    assert not stale, f"these sites are converted or gone: take them off REMAINING: {stale}"


@pytest.mark.parametrize("unit, constant", [("Msun", constants.SOLAR_MASS_TO_KG), ("Mearth", constants.EARTH_MASS_TO_KG), ("Mjup", constants.JUPITER_MASS_TO_KG)])
def test_the_converters_are_bit_for_bit_the_hand_multiply_and_divide(unit, constant):
    to_kg, from_kg = unit_registry.converter(unit, "kg"), unit_registry.converter("kg", unit)
    for x in (1.0, 0.075, 3.141592653589793, 1.0e-9, 123456.789, 2.5e30):
        assert to_kg(x) == x * constant
        assert from_kg(x) == x / constant
