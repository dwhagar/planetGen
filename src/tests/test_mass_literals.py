"""GEN.1152: the mass of Jupiter, the Sun and Earth is written once, in `physics/constants.py`; the modules that had their own rounded copy read it."""

import ast
import pathlib

from planetgen import tuning
from planetgen.physics import constants, giant_interior, moon_composition, ring_origin, stellar_abundances

OLD_LITERALS = (1.89813e27, 1.989e30, 317.83)
"""tuple: The rounded copies the sites had: Jupiter's mass in kg (moon_composition, stellar_abundances), the Sun's in kg (ring_origin), Jupiter in Earth masses (giant_interior)."""


def test_every_site_reads_the_one_value():
    assert moon_composition.JUPITER_MASS_KG == constants.JUPITER_MASS_TO_KG
    assert ring_origin.SATURN["star_mass_kg"] == constants.SOLAR_MASS_TO_KG
    assert giant_interior.EARTH_MASS_JUPITER == constants.JUPITER_MASS_EARTH
    assert stellar_abundances.JUPITER_MASS_SOL == constants.JUPITER_MASS_TO_KG / constants.SOLAR_MASS_TO_KG
    assert tuning.ROGUE_PLANET_MASS_BIN_EDGES_EARTH["saturn"][1] == constants.JUPITER_MASS_EARTH == tuning.ROGUE_PLANET_MASS_BIN_EDGES_EARTH["jupiter"][0]


def test_jupiter_in_earth_masses_is_the_ratio_of_the_two_conversion_factors():
    assert constants.JUPITER_MASS_EARTH == constants.JUPITER_MASS_TO_KG / constants.EARTH_MASS_TO_KG
    assert abs(constants.JUPITER_MASS_EARTH - 317.83) < 0.01


def test_no_module_of_the_physics_writes_the_rounded_copy_again():
    root = pathlib.Path(__file__).resolve().parents[1] / "planetgen"
    found = []
    for path in sorted(root.rglob("*.py")):
        if path.name == "constants.py":
            continue
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Constant) and isinstance(node.value, float) and any(abs(node.value - old) <= 1e-9 * old for old in OLD_LITERALS):
                found.append(f"{path.relative_to(root)}:{node.lineno} {node.value}")
    assert not found, f"use the value in physics/constants.py: {found}"
